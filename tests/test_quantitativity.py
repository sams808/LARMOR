"""larmor.quantitativity (Qt-free): TopSpin ct1t2.txt reading, the sibling-T1
search, the title parser, the recovery / flip arithmetic at its limits and
the Check texts -- on inline fixtures copied from the real files, synthetic
sample folders written to tmp_path, and real NMRFAM anchors that skip
cleanly elsewhere (conftest.require)."""
import dataclasses
import json
from pathlib import Path

import pytest

from conftest import _D, LAW_CA_11B, require
from larmor import fithealth, quantitativity as Q, satrec
from larmor.io import bruker, scan
from larmor.recipe import Param, SiteModel

P5BI812 = _D / "2026-05/04272026_P5-Bi8-12_SS_ALP"
CA12F = _D / "2026-01/01222026_SR31648_2Ca-12F_SS_ALP"
BASE2CA = _D / "2026-01/01202026_SR31649_Base2Ca_SS_ALP"

# Base1Ca/23/pdata/1/ct1t2.txt (legacy layout), verbatim
CT1T2_LEGACY = """Dataset :
 C:/Users/samso/Desktop/WSU_work/NMR/NMRFAM/DATA/01202026_SR31649_Base1Ca_SS_ALP/23/pdata/1
AREA fit :
 I[t]=I[0]+P*exp(-t/T1)

7 points for Integral 1,  Integral Region from 27.382 to 5.846 ppm
Results     Comp. 1

I[0]  =    9.932e-01
P     =   -5.560e-01
T1    =       4.643s
SD    =    1.207e-02

    tau    ppm     integral    intensity

    0.000n    14.092     eliminated     eliminated
  250.000m    14.092     eliminated     eliminated
  500.000m    14.092   2.5015e+10     eliminated
    1.000s    14.092   2.6209e+10     eliminated
   32.000s    14.092   4.9895e+10     eliminated



7 points for Integral 2,  Integral Region from 5.846 to -7.618 ppm
Results     Comp. 1

I[0]  =    9.957e-01
P     =   -7.106e-01
T1    =       3.896s
SD    =    1.350e-02

    tau    ppm     integral    intensity

    0.000n     1.114     eliminated   5.1332e+08
   32.000s     1.072   6.4429e+10     eliminated


"""

# 2025-12/pos8A2O-20F/3/pdata/1/ct1t2.txt (SIMFIT, four components), verbatim
CT1T2_SIMFIT_MULTI = """SIMFIT RESULTS
==============

Dataset : C:/Users/samso/Desktop/WSU_work/NMR/NMRFAM/DATA/pos8A2O-20F/3/pdata/1/ct1t2.txt

AREA fit : Saturation-Recovery(T1) : I[t]=I[0](1-exp(-t/T1))

12 points for Integral 1,  Integral Region from -17.702 to -176.036 ppm

Converged after 1212 iterations!

Results     Comp. 1       Comp. 2       Comp. 3       Comp. 4

I[0]  =    1.324e-01     6.084e-01     1.283e+01    -1.257e+01
T1    =      60.764m        1.759s        4.539s        4.466s

RSS   =    3.079e-03
SD    =    1.602e-02

Point       Tau          Expt          Calc       Difference

    1       0.000n     4.577e-02     0.000e+00   -4.577e-02
   12      25.600s     1.000e+00     9.958e-01   -4.191e-03

============================================================

"""

# 2026-07/07062026_SR31648_Peyton-S2_SS_ALP/1901 (INTENSITY fit on a peak point)
CT1T2_PEAK = """Dataset :
 /data/soudani_data/07062026_SR31648_Peyton-S2_SS_ALP/1901/pdata/1
INTENSITY fit :
 I[t]=I[0]+P*exp(-t/T1)

4 points for Peak 1,  Peak Point at -68.561 ppm
Results     Comp. 1

I[0]  =    1.149e+00
P     =   -1.085e+00
T1    =       8.094s
SD    =    3.399e-03

    tau    ppm     integral    intensity

    0.000n   -68.561  -3.4032e+10     eliminated
"""


def _simfit_single(t1_tok: str) -> str:
    """The SIMFIT fixture with a single-component T1 line (Kogarkoite style)."""
    import re

    return re.sub(r"^T1\s*=.*$", f"T1    =     {t1_tok}", CT1T2_SIMFIT_MULTI,
                  flags=re.MULTILINE)


# ------------------------------------------------------------- fixtures
def _jcamp(path: Path, **params):
    """A minimal TopSpin JCAMP parameter file; a list value is written in the
    '##$D= (0..63)' array layout nmrglue parses."""
    lines = ["##TITLE= Parameter file, TopSpin 4.1", "##JCAMPDX= 5.0",
             "##DATATYPE= Parameter Values", "##ORIGIN= Bruker BioSpin GmbH",
             "##OWNER= nmr"]
    for k, v in params.items():
        if isinstance(v, (list, tuple)):
            lines.append(f"##${k}= (0..{len(v) - 1})")
            lines.append(" ".join(str(x) for x in v))
        else:
            lines.append(f"##${k}= {v}")
    lines.append("##END=")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _expno(root: Path, expno, nuc="11B", pulprog="zg", d1=14.0, p1=0.425,
           plw1=100.0, ns=256, td=7988, swh=100000, title="", vdlist=None,
           ct1t2=None, raw=None):
    e = root / str(expno)
    _jcamp(e / "acqus", NUC1=f"<{nuc}>", PULPROG=f"<{pulprog}>", BF1=192.4,
           SFO1=192.43, O1=0, SW_h=swh, TD=td, NS=ns, PARMODE=0,
           PROBHD="<16_Solenoid (PMAS16)>", D=[0.2, d1, 0, 0], P=[0, p1, 0, 0],
           PLW=[0, plw1, 0, 0])
    pd = e / "pdata" / "1"
    pd.mkdir(parents=True, exist_ok=True)
    (pd / "title").write_text(title, encoding="utf-8")
    if raw is None:
        raw = "ser" if "satrec" in pulprog else "fid"
    (e / raw).write_bytes(b"")
    if raw == "fid":
        (pd / "1r").write_bytes(b"")
    if vdlist is not None:
        (e / "vdlist").write_text("\n".join(str(v) for v in vdlist) + "\n")
    if ct1t2 is not None:
        (pd / "ct1t2.txt").write_text(ct1t2, encoding="utf-8")
    return e


VDLIST = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256]
CT1T2_DIVERGED = CT1T2_LEGACY.replace("4.643s", "3500037.796s")


def _sample(tmp_path: Path) -> Path:
    """The synthetic sample of the design: two 23Na, two 11B, three 31P and
    three 27Al EXPNOs with and without TopSpin fits."""
    root = tmp_path / "01202026_SR31649_Base1Ca_SS_ALP"
    _expno(root, 12, "23Na", "satrect1", d1=0.5, title="23Na Satrec 993m",
           vdlist=VDLIST, ct1t2=CT1T2_DIVERGED)
    _expno(root, 13, "23Na", "zg", d1=7.5, p1=0.5, plw1=150, ns=128)
    _expno(root, 23, "11B", "satrect1", d1=0.05, title="11B zg - satrec 4.64 3.74s",
           vdlist=VDLIST, ct1t2=CT1T2_LEGACY)
    _expno(root, 24, "11B", "zg", title="11B with short tip angle")
    _expno(root, 3101, "31P", "satrect1", d1=20, p1=1.25, plw1=207, ns=1,
           vdlist=[0.25, 1, 4, 16, 64, 256, 1024, 4096])
    _expno(root, 3102, "31P", "zg", d1=300, p1=1.25, plw1=207, ns=14, td=9590,
           title="31P\nP1(90)=3.750; 30 deg tip (pi/6) = 1.25")
    _expno(root, 3112, "31P", "satrect1", d1=20, p1=1.25, plw1=207, ns=1,
           vdlist=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048],
           ct1t2=CT1T2_LEGACY.replace("27.382 to 5.846", "44.082 to -59.276")
           .replace("4.643s", "83.295s"))
    _expno(root, 2701, "27Al", "satrect1", d1=0.05, p1=0.375, plw1=207, ns=2,
           vdlist=VDLIST[:11],
           ct1t2=CT1T2_LEGACY.replace("27.382 to 5.846", "100.592 to -20.595")
           .replace("4.643s", "644.856m"))
    _expno(root, 2703, "27Al", "satrect1", d1=0.05, p1=0.375, plw1=207, ns=2,
           vdlist=VDLIST[:11])
    _expno(root, 2704, "27Al", "zg", d1=25, p1=0.375, plw1=207, ns=512, td=9590,
           title="27Al\nP1(90)=3.125; 11deg tip = 0.375")
    # the stray operator export 2Ca-12F carries next to its EXPNOs
    (root / "01202026_SR31649_Base1Ca_SS_ALP_24_3_norm-to-area.txt").write_text("x")
    return root


def _sites(*ppms, labels=None):
    out = []
    for k, p in enumerate(ppms):
        out.append(SiteModel(model="gauss_lor",
                             label=(labels[k] if labels else f"S{k}"),
                             params={"isotropic_chemical_shift_ppm": Param(p),
                                     "shift_fwhm_ppm": Param(5.0),
                                     "amplitude": Param(1.0),
                                     "gl": Param(0.5, vary=False)}))
    return out


def _base1ca_t1() -> Q.T1Source:
    return Q.T1Source(expno="/x/23", kind="ct1t2", regions=[
        satrec.T1Region(1, 27.382, 5.846, 4.643, (4.643,), 1.207e-2, 7),
        satrec.T1Region(2, 5.846, -7.618, 3.896, (3.896,), 1.350e-2, 7)],
        vdlist_max_s=256.0, plausible=True, title_notes_s=[4.64, 3.74])


def _acq(**over) -> Q.Acquisition:
    base = dict(expno="/x/24", nucleus="11B", pulprog="zg", kind="Single pulse",
                ns=256, d1_s=14.0, aq_s=0.04, p1_us=0.425, plw1_w=100.0,
                probhd="16_Solenoid (PMAS16)", title="11B with short tip angle",
                p90_us_title=None, flip_deg_title=None, t1_multiple_claimed=None)
    base.update(over)
    return Q.Acquisition(**base)


def _facts(acq=None, t1="base1ca", spin=1.5, status="ok", sibling="/x/23"):
    if t1 == "base1ca":
        t1 = _base1ca_t1()
    return Q.AcqFacts(acquisition=acq if acq is not None else _acq(), t1=t1,
                      sibling_expno=sibling, spin=spin, folder="/x",
                      t1_status=status if t1 is not None else
                      (status if status != "ok" else "none"))


# -------------------------------------------------------------- ct1t2
def test_read_ct1t2_legacy_two_regions(tmp_path):
    pd = tmp_path / "pdata" / "1"
    pd.mkdir(parents=True)
    (pd / "ct1t2.txt").write_text(CT1T2_LEGACY, encoding="utf-8")
    regions = satrec.read_ct1t2(pd)
    assert len(regions) == 2
    r1, r2 = regions
    assert (r1.index, r1.hi_ppm, r1.lo_ppm, r1.t1_s, r1.sd, r1.npoints) == \
        (1, 27.382, 5.846, 4.643, 1.207e-2, 7)
    assert (r2.index, r2.hi_ppm, r2.lo_ppm, r2.t1_s, r2.sd, r2.npoints) == \
        (2, 5.846, -7.618, 3.896, 1.350e-2, 7)
    assert r1.components_s == (4.643,)
    assert satrec.read_ct1t2(pd / "ct1t2.txt") == regions     # the file itself
    src = Q.T1Source(expno="23", kind="ct1t2", regions=regions)
    assert src.t1_for_ppm(15.0) == (4.643, r1, True)
    assert src.t1_for_ppm(0.0) == (3.896, r2, True)
    assert src.t1_for_ppm(40.0) == (4.643, r1, False)          # nearest region
    assert src.t1_for_ppm(-30.0) == (3.896, r2, False)
    assert Q.T1Source(expno="", kind="ct1t2").t1_for_ppm(1.0) == (None, None, False)


def test_read_ct1t2_simfit_multicomponent_peak_and_units(tmp_path):
    f = tmp_path / "ct1t2.txt"
    f.write_text(CT1T2_SIMFIT_MULTI, encoding="utf-8")
    (r,) = satrec.read_ct1t2(f)
    assert r.components_s == pytest.approx((0.060764, 1.759, 4.539, 4.466))
    assert r.t1_s == 4.539                                  # the longest component
    assert (r.hi_ppm, r.lo_ppm, r.npoints, r.sd) == (-17.702, -176.036, 12, 1.602e-2)
    f.write_text(_simfit_single("2.126s"), encoding="utf-8")
    assert satrec.read_ct1t2(f)[0].t1_s == 2.126
    f.write_text(_simfit_single("515.046m"), encoding="utf-8")
    assert satrec.read_ct1t2(f)[0].t1_s == pytest.approx(0.515046)
    f.write_text(CT1T2_PEAK, encoding="utf-8")
    (pk,) = satrec.read_ct1t2(f)
    assert pk.hi_ppm == pk.lo_ppm == -68.561 and pk.t1_s == 8.094 and pk.npoints == 4
    # a fit that never produced a result, a missing file, garbage: [] and no raise
    f.write_text("Dataset :\n /data/x/1/pdata/1\n", encoding="utf-8")
    assert satrec.read_ct1t2(f) == []
    assert satrec.read_ct1t2(tmp_path / "nowhere") == []
    f.write_bytes(bytes(range(256)) * 4)
    assert satrec.read_ct1t2(f) == []
    assert satrec._seconds("250m") == 0.25
    assert satrec._seconds("100u") == pytest.approx(1e-4)
    assert satrec._seconds("4s") == 4.0
    assert satrec._seconds("0.000n") == 0.0
    assert satrec._seconds("3") == 3.0
    with pytest.raises(ValueError):
        satrec._seconds("abc")


def test_ct1t2_diverged_fit_is_implausible(tmp_path):
    e = _expno(tmp_path, 12, "23Na", "satrect1", vdlist=VDLIST, ct1t2=CT1T2_DIVERGED)
    src = Q.read_t1_source(e)
    assert src.plausible is False and "did not converge" in src.note
    assert "3.5e+06" in src.note and "256 s" in src.note
    assert src.vdlist_max_s == 256.0 and src.kind == "ct1t2"
    ok = _expno(tmp_path, 13, "23Na", "satrect1", vdlist=VDLIST,
                ct1t2=CT1T2_LEGACY.replace("4.643s", "53s"))
    src = Q.read_t1_source(ok)
    assert src.plausible is True and src.regions[0].t1_s == 53.0
    novd = _expno(tmp_path, 14, "23Na", "satrect1", ct1t2=CT1T2_LEGACY)
    src = Q.read_t1_source(novd)
    assert src.plausible is True and src.vdlist_max_s is None
    empty = _expno(tmp_path, 15, "23Na", "satrect1", vdlist=VDLIST)
    src = Q.read_t1_source(empty)
    assert src.regions == [] and "no T1 result" in src.note
    assert Q.read_t1_source(tmp_path / "99") is None


# -------------------------------------------------------------- titles
@pytest.mark.parametrize("title,p90,flip,mult", [
    ("P1(90)=3.750; 30 deg tip (pi/6) = 1.25", 3.75, 30.0, None),
    ("P1(90)=3.125; 11deg tip = 0.375", 3.125, 11.0, None),
    ("27Al zg - spectrum with 11 degree tip angle and 5*T1", None, 11.0, 5.0),
    ("23Na zg spectrum using 5xT1 and 18 degree pulse", None, 18.0, 5.0),
    ("11B zg - 5*T1 with 18 degree tip angle", None, 18.0, 5.0),
    ("11B with short tip angle", None, None, None),
    ("VT 25 degrees", None, None, None),
    ("VT ambient - no flow", None, None, None),
    ("", None, None, None),
])
def test_parse_title_phrases_from_the_tree(title, p90, flip, mult):
    got = Q.parse_title(title)
    assert (got["p90_us"], got["flip_deg"], got["t1_multiple"]) == (p90, flip, mult)


@pytest.mark.parametrize("title,notes", [
    ("11B zg - satrec 4.64 3.74s", [4.64, 3.74]),
    ("T1 approximately 53s", [53.0]),
    ("Longest T1 at 27.6 ms", [0.0276]),
    ("23Na Satrec 993m", [0.993]),
    ("T1 10.6 29.7s", [10.6, 29.7]),
    ("T1 approximately xx ms; uncertainty high", []),
    ("19F Hahn echo with 1 rotor period", []),
])
def test_parse_title_t1_notes(title, notes):
    assert Q.parse_title(title)["t1_notes_s"] == pytest.approx(notes)


# ---------------------------------------------------------- arithmetic
def test_flip_limit_effective_angle_and_steady_state_recovery():
    assert Q.flip_limit_deg(1.5) == 15.0
    assert Q.flip_limit_deg(2.5) == 10.0
    assert Q.flip_limit_deg(3.5) == pytest.approx(7.5)
    assert Q.flip_limit_deg(0.5) is None and Q.flip_limit_deg(None) is None
    assert Q.flip_limit_deg(1.0) is None                    # integer spin
    assert Q.effective_flip_deg(11, 2.5) == pytest.approx(33.0)
    assert Q.effective_flip_deg(30, 0.5) == 30.0
    assert Q.effective_flip_deg(100, 2.5) == 180.0          # capped
    assert Q.spin_text(1.5) == "3/2" and Q.spin_text(0.5) == "1/2"
    T1 = 4.643
    assert Q.steady_state_recovery(3 * T1, T1) == pytest.approx(1 - 2.718281828 ** -3, abs=1e-6)
    assert Q.steady_state_recovery(3 * T1, T1) == pytest.approx(0.9502, abs=1e-4)
    assert Q.steady_state_recovery(3 * T1, T1, 90, None) == \
        pytest.approx(Q.steady_state_recovery(3 * T1, T1))
    assert Q.steady_state_recovery(3 * T1, T1, 0.01, 1.5) > 0.9999
    rs = [Q.steady_state_recovery(t, T1) for t in (1.0, 2.0, 5.0, 20.0)]
    assert rs == sorted(rs)                                  # monotonic in recycle
    # Base1Ca/24 vs /23 at the assumed 90 degrees
    assert Q.steady_state_recovery(14.04, 4.643) == pytest.approx(0.951, abs=1e-3)
    assert Q.steady_state_recovery(14.04, 3.896) == pytest.approx(0.973, abs=1e-3)
    # the 2026-05 31P series: 3.6 T1 at a 30-degree pulse recovers 99.6 %,
    # while the saturation-recovery form would flag it at 97.3 %
    assert Q.steady_state_recovery(300.05, 83.295, 30, 0.5) == pytest.approx(0.996, abs=1e-3)
    assert Q.steady_state_recovery(300.05, 83.295, 30, 0.5) >= Q.RECOVERY_MIN
    assert Q.steady_state_recovery(300.05, 83.295) == pytest.approx(0.973, abs=1e-3)
    assert Q.steady_state_recovery(300.05, 83.295) < Q.RECOVERY_MIN
    # 2Ca-12F Hahn echo, two 19F regions
    assert Q.steady_state_recovery(60.0, 28.424) == pytest.approx(0.879, abs=1e-3)
    assert Q.steady_state_recovery(60.0, 12.685) == pytest.approx(0.991, abs=1e-3)
    # the group's 5*T1 rule at 90 degrees clears the 99 % threshold
    assert Q.steady_state_recovery(5 * T1, T1) == pytest.approx(0.9933, abs=1e-4)
    assert Q.steady_state_recovery(5 * T1, T1) >= Q.RECOVERY_MIN
    assert Q.steady_state_recovery(1.0, 0.0) == 1.0          # no T1: no judgement


# --------------------------------------------------------------- check
def test_check_levels_and_texts_synthetic():
    facts = _facts()
    chk = Q.check(facts, _sites(15.0, 0.0))
    assert [round(r.recovery, 3) for r in chk.recoveries] == [0.951, 0.973]
    assert [r.exact for r in chk.recoveries] == [True, True]
    assert chk.recovery_min < Q.RECOVERY_MIN
    assert chk.recovery_spread == pytest.approx(0.0214, abs=2e-3)
    assert chk.flip_assumed_90 and chk.flip_deg is None and chk.flip_source == ""
    assert chk.flip_effective_deg == 90.0 and chk.excitation_judged
    assert chk.recovery_text() == "D1 = 3.0–3.6 T1 → 95–97 % (90° assumed)"
    assert chk.excitation_text() is None
    assert chk.excitation_unknown_text() == "flip angle unknown (I = 3/2)"
    assert "90° pulse" in chk.excitation_unknown_detail()
    d = chk.recovery_detail()
    for bit in ("14.04 s", "EXPNO 23", "27.4 … 5.8 ppm", "4.64 s", "3.0 T1", "95.1 %",
                "biased by up to 2.2 %", "90° assumed", "note in EXPNO 23: T1 4.64 3.74 s"):
        assert bit in d, (bit, d)
    assert chk.recovery_ok() is False and chk.t1_status == "ok"
    # a typed 90-degree pulse turns the short pulse into its steady state
    chk = Q.check(facts, _sites(15.0, 0.0), {"p90_us": 3.5})
    assert chk.flip_deg == pytest.approx(10.93, abs=0.01) and chk.flip_source == "user"
    assert chk.flip_effective_deg == pytest.approx(21.86, abs=0.02)
    assert chk.recovery_min >= Q.RECOVERY_MIN and not chk.flip_assumed_90
    assert chk.recovery_text().startswith("D1 = 3.0–3.6 T1 at 10.9° → ")
    assert "(90° assumed)" not in chk.recovery_text()
    assert any(line.startswith("recycle 14 s = 3.0–3.6 T1 → ") for line in chk.passing_lines())
    assert "flip 10.9° within the linear regime (≤ 15°)" in chk.passing_lines()
    # a typed flip angle wins over the typed pulse
    chk = Q.check(facts, _sites(15.0), {"p90_us": 3.5, "flip_deg": 8.0})
    assert chk.flip_deg == 8.0 and chk.flip_source == "user"
    # the title's P1(90) with a longer P1: over the 15-degree limit
    acq = _acq(title="P1(90)=2.5", p90_us_title=2.5, p1_us=0.5)
    chk = Q.check(_facts(acq), _sites(15.0))
    assert chk.flip_deg == 18.0 and chk.flip_source == "P1(90)"
    assert chk.excitation_over()
    assert chk.excitation_text() == "flip 18° > 15° limit (I = 3/2)"
    assert "just above the limit (20 %)" in chk.excitation_detail()
    assert "Edén Eq. 38" in chk.excitation_detail()
    # the title's stated tip alone
    chk = Q.check(_facts(_acq(flip_deg_title=11.0)), _sites(15.0))
    assert chk.flip_deg == 11.0 and chk.flip_source == "title" and not chk.excitation_over()
    # spin-1/2: excitation not judged, any angle quantitative
    chk = Q.check(_facts(_acq(nucleus="31P"), spin=0.5), _sites(15.0))
    assert chk.excitation_judged is False and chk.flip_limit_deg is None
    assert "spin-½: any flip angle is quantitative" in chk.passing_lines()
    # an echo program: 90 degrees by construction, excitation not judged
    chk = Q.check(_facts(_acq(pulprog="hahnecho.nmrfam", kind="Hahn echo")), _sites(15.0))
    assert chk.excitation_judged is False and chk.flip_source == "echo(90)"
    assert chk.flip_deg == 90.0 and not chk.flip_assumed_90
    assert chk.recovery_text() == "D1 = 3.0 T1 at 90° → 95 %"
    assert "echo/CPMG: excitation uniformity not judged" in chk.passing_lines()
    # a per-site T1 override beats the TopSpin region; a scalar T1 fills the rest
    chk = Q.check(facts, _sites(15.0, 0.0, labels=["A", "B"]),
                  {"t1_by_site": {"A": 20.0}, "t1_s": 2.0})
    a, b = chk.recoveries
    assert a.t1_s == 20.0 and a.source == "per_site" and a.exact
    assert b.t1_s == 2.0 and b.source == "user"
    assert "typed / per-site T1" in chk.recovery_detail()
    # backgrounds carry no delta_iso meaning: skipped
    bg = SiteModel(model="spectrum", label="bg", params={
        "isotropic_chemical_shift_ppm": Param(0.0), "amplitude": Param(1.0)})
    chk = Q.check(facts, _sites(15.0) + [bg])
    assert [r.label for r in chk.recoveries] == ["S0"]
    # site dicts (the recipe dict of the desktop) are accepted too
    chk = Q.check(facts, [s.__dict__ | {"params": {k: {"value": p.value} for k, p in s.params.items()}}
                          for s in _sites(1.0)])
    assert chk.recoveries[0].t1_s == 3.896
    # no T1 anywhere: the unknown wording; no sibling: 'none'
    chk = Q.check(_facts(t1=None, status="missing"), _sites(15.0))
    assert chk.recoveries == [] and chk.recovery_min is None
    assert chk.t1_status == "missing"
    assert chk.recovery_text() is None
    assert chk.recovery_unknown_text() == "recycle 14 s — T1 unknown"
    assert "ct1t2" in chk.recovery_unknown_detail()
    chk = Q.check(_facts(t1=None, status="none", sibling=None), _sites(15.0))
    assert chk.t1_status == "none" and "no T1 EXPNO" in chk.recovery_unknown_detail()
    chk = Q.check(_facts(t1=None, status="implausible"), _sites(15.0))
    assert chk.t1_status == "implausible" and "did not converge" in chk.recovery_unknown_detail()
    # every chip text stays within the strip's bound, whatever the labels
    long = _sites(*range(6), labels=["a-very-long-site-label-%d" % k for k in range(6)])
    chk = Q.check(facts, long)
    for t in (chk.recovery_text(), chk.excitation_unknown_text()):
        assert len(t) <= fithealth.TEXT_LIMIT
    # plain scalars only (Health.__eq__ compares the Check)
    json.dumps(dataclasses.asdict(chk))
    assert Q.check(facts, _sites(15.0)) == Q.check(facts, _sites(15.0))


def test_facts_for_resolves_any_bruker_path_fails_soft_and_picks_the_fitted_sibling(tmp_path):
    root = _sample(tmp_path)
    f = Q.facts_for(root / "24")
    assert f is not None and f == Q.facts_for(root / "24" / "pdata" / "1" / "1r")
    a = f.acquisition
    assert (a.d1_s, a.p1_us, a.plw1_w, a.ns, a.nucleus, a.kind) == \
        (14.0, 0.425, 100.0, 256, "11B", "Single pulse")
    assert a.aq_s == pytest.approx(0.03994) and a.recycle_s == pytest.approx(14.03994)
    assert a.probhd == "16_Solenoid (PMAS16)" and a.pulprog == "zg"
    assert f.spin == 1.5 and f.t1_status == "ok"
    assert f.t1.expno.endswith("23") and f.t1.name == "23" and f.sibling_expno.endswith("23")
    assert [r.t1_s for r in f.t1.regions] == [4.643, 3.896]
    assert f.t1.title_notes_s == [4.64, 3.74] and f.t1.vdlist_max_s == 256.0
    assert f.folder == str(root)
    # the nearest sibling has no ct1t2: the fitted one wins, the nearest is named
    f = Q.facts_for(root / "3102")
    assert f.t1.expno.endswith("3112") and f.sibling_expno.endswith("3101")
    assert f.acquisition.p90_us_title == 3.75 and f.acquisition.flip_deg_title == 30.0
    assert f.spin == 0.5
    f = Q.facts_for(root / "2704")
    assert f.t1.expno.endswith("2701") and f.sibling_expno.endswith("2703")
    assert f.t1.regions[0].t1_s == pytest.approx(0.644856)
    # a diverged TopSpin fit is not a T1
    f = Q.facts_for(root / "13")
    assert f.t1 is None and f.t1_status == "implausible" and f.sibling_expno.endswith("12")
    assert "did not converge" in f.t1_note
    chk = Q.check(f, _sites(-10.0))
    assert chk.recoveries == [] and chk.t1_status == "implausible"
    # not a Bruker path at all
    assert Q.facts_for("src") is None
    assert Q.facts_for("") is None and Q.facts_for(None) is None
    (tmp_path / "x.fxmla").write_text("<dmfit/>")
    assert Q.facts_for(tmp_path / "x.fxmla") is None
    # a sample without any T1 EXPNO of the nucleus
    lone = tmp_path / "lone"
    _expno(lone, 5, "29Si", "zg", d1=60)
    _expno(lone, 6, "27Al", "satrect1", vdlist=VDLIST, ct1t2=CT1T2_LEGACY)
    f = Q.facts_for(lone / "5")
    assert f.t1 is None and f.sibling_expno is None and f.t1_status == "none"
    assert "no T1 EXPNO of 29Si" in f.t1_note
    chk = Q.check(f, _sites(-100.0))
    assert chk.recovery_unknown_text() == "recycle 60 s — T1 unknown"
    # the readers behind it fail soft
    assert bruker.read_acqus_meta(tmp_path / "nowhere") == {}
    assert Q.read_acquisition(tmp_path / "nowhere") is None
    assert scan.find_sibling_t1(tmp_path / "nowhere" / "1", "11B") == []


def test_real_base1ca_24_and_p5bi8_12_3102_and_2ca12f_6():
    path = require(LAW_CA_11B[1])
    f = Q.facts_for(path)
    a = f.acquisition
    assert (a.d1_s, a.ns, a.p1_us, a.plw1_w) == (14.0, 256, 0.425, 100.0)
    assert a.probhd == "16_Solenoid (PMAS16)" and a.kind == "Single pulse"
    assert a.flip_deg_title is None and a.p90_us_title is None     # 'short tip angle'
    assert f.t1.expno.endswith("23") and f.t1_status == "ok"
    assert [(r.hi_ppm, r.lo_ppm, r.t1_s) for r in f.t1.regions] == pytest.approx(
        [(27.382, 5.846, 4.643), (5.846, -7.618, 3.896)], abs=1e-3)
    assert f.t1.title_notes_s == [4.64, 3.74]
    chk = Q.check(f, _sites(15.0, 1.0))
    assert chk.recovery_min == pytest.approx(0.951, abs=1e-3)
    assert chk.recovery_min < Q.RECOVERY_MIN and chk.flip_assumed_90
    assert chk.recovery_text() == "D1 = 3.0–3.6 T1 → 95–97 % (90° assumed)"
    assert chk.excitation_judged and chk.flip_deg is None

    f = Q.facts_for(require(P5BI812 / "3102"))
    a = f.acquisition
    assert (a.d1_s, a.p1_us, a.plw1_w, a.p90_us_title) == (300.0, 1.25, 207.0, 3.75)
    chk = Q.check(f, _sites(0.0))
    assert chk.flip_deg == 30.0 and chk.flip_source == "P1(90)"
    assert f.spin == 0.5 and chk.excitation_judged is False
    assert f.t1.expno.endswith("3112") and f.sibling_expno.endswith("3101")
    assert f.t1.regions[0].t1_s == pytest.approx(83.295)
    (r,) = chk.recoveries
    assert r.ratio == pytest.approx(3.6, abs=0.01)
    assert r.recovery == pytest.approx(0.996, abs=1e-3) and chk.recovery_ok()

    f = Q.facts_for(require(P5BI812 / "2704"))
    chk = Q.check(f, _sites(55.0))
    assert chk.flip_deg == pytest.approx(10.8) and chk.flip_limit_deg == 10.0
    assert chk.excitation_text() == "flip 10.8° > 10° limit (I = 5/2)"
    assert "just above the limit (8 %)" in chk.excitation_detail()
    assert f.t1.expno.endswith("2701") and chk.recovery_ok()

    f = Q.facts_for(require(CA12F / "6"))
    a = f.acquisition
    assert a.kind == "Hahn echo" and a.d1_s == 60.0
    chk = Q.check(f, _sites(-200.0, -120.0))
    assert chk.flip_source == "echo(90)"
    assert [(r.hi_ppm, r.lo_ppm, r.t1_s) for r in f.t1.regions] == pytest.approx(
        [(-87.487, -162.517, 12.685), (-162.517, -236.681, 28.424)], abs=1e-3)
    assert f.t1.title_notes_s == [53.0]
    assert chk.recovery_min == pytest.approx(0.879, abs=1e-3)
    assert chk.recovery_text().startswith("D1 = 2.1–4.7 T1 at 90° → 88–99 %")
    assert "note in EXPNO 1: T1 53 s" in chk.recovery_detail()

    f = Q.facts_for(require(BASE2CA / "13"))
    assert f.t1 is None and f.t1_status == "implausible"
    assert f.sibling_expno.endswith("12") and "did not converge" in f.t1_note
