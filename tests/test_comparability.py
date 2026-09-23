"""Comparability check (larmor.comparability): reading acqus / procs /
auditp, the majority rule and its levels, the TopSpin -> LARMOR ops mapping,
the reprocess path through loader.apply_processing -- on synthetic JCAMP
fixtures, plus the real five-glass 11B series and the 27Al 2702 fid/1r pair
that pin the phase and TDeff conventions."""
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from conftest import BRUKER_1R, BRUKER_FID, LAW_CA_11B, require
from larmor import comparability as C
from larmor import processing as proc
from larmor.io import bruker


# ---------------------------------------------------------------- fixtures
def _jcamp(path: Path, params: dict, arrays: dict | None = None):
    lines = ["##TITLE= Parameter file, TopSpin 4.1", "##JCAMPDX= 5.0",
             "##DATATYPE= Parameter Values", "##ORIGIN= Bruker BioSpin GmbH",
             "##OWNER= nmr"]
    for k, v in params.items():
        lines.append(f"##${k}= {v}")
    for k, vals in (arrays or {}).items():
        lines.append(f"##${k}= (0..{len(vals) - 1})")
        lines.append(" ".join(str(x) for x in vals))
    lines.append("##END=")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _arr(i, v, n=64):
    out = [0] * n
    out[i] = v
    return out


def _fake_expno(tmp_path, sample, expno="24", *, lb=0, tdeff=1024, si=32768, wdw=1,
                gb=0, ssb=0, ph_mod=1, phc0=-160.0, phc1=29.9, absg=0, bc_mod=0,
                fcor=1, sf=192.431528, sfo1=192.430693, bf1=192.430693, nuc="11B",
                pulprog="zg", ns=256, rg=194.07, td=7988, sw_h=100000, d1=12.5,
                p1=0.425, plw1=100.0, pl1=None, probe="Z1 3.2mm",
                date=1_768_800_000, with_1r=True, with_fid=False, with_procs=True,
                audit=None, title="11B with short tip angle\nRotor SR31649"):
    """A synthetic EXPNO (acqus + pdata/1/procs + title [+ 1r stub, fid stub,
    auditp.txt]); returns the pdata/1/1r path (the Explorer's batch path) or
    the EXPNO folder when ``with_1r`` is False. read_params never reads the
    data files, so 4-byte stubs suffice."""
    e = tmp_path / sample / str(expno)
    arrays = {"D": _arr(1, d1), "P": _arr(1, p1)}
    if plw1 is not None:
        arrays["PLW"] = _arr(1, plw1)
    if pl1 is not None:
        arrays["PL"] = _arr(1, pl1)
    _jcamp(e / "acqus", {"NUC1": f"<{nuc}>", "PULPROG": f"<{pulprog}>", "NS": ns,
                         "TD": td, "RG": rg, "SW_h": sw_h, "SFO1": sfo1, "BF1": bf1,
                         "O1": (sfo1 - bf1) * 1e6, "DATE": date,
                         "PROBHD": f"<{probe}>", "MASR": 20000}, arrays)
    pd = e / "pdata" / "1"
    pd.mkdir(parents=True, exist_ok=True)
    if with_procs:
        _jcamp(pd / "procs", {"WDW": wdw, "LB": lb, "GB": gb, "SSB": ssb,
                              "TDeff": tdeff, "SI": si, "FCOR": fcor,
                              "PH_mod": ph_mod, "PHC0": phc0, "PHC1": phc1,
                              "ABSG": absg, "BC_mod": bc_mod, "FT_mod": 6,
                              "PKNL": "yes", "SF": sf, "OFFSET": 100.0,
                              "SW_p": sw_h, "SR": (sf - bf1) * 1e6})
    (pd / "title").write_text(title, encoding="utf-8")
    if with_1r:
        (pd / "1r").write_bytes(b"\0\0\0\0")
    if with_fid:
        (e / "fid").write_bytes(b"\0\0\0\0")
    if audit is not None:
        (pd / "auditp.txt").write_text(audit, encoding="utf-8")
    return pd / "1r" if with_1r else e


#: 04272026_P5-Bi8-12_SS_ALP/3104/pdata/1/auditp.txt, as TopSpin wrote it
AUDIT_3104 = """##TITLE= Audit trail, TopSpin 4.5.0
##JCAMPDX= 5.01
##ORIGIN= Bruker BioSpin GmbH & Co. KG
##OWNER= samso
$$ C:/Users/x/DATA/2026-5/04272026_P5-Bi8-12_SS_ALP/3104/pdata/1/auditp.txt
##AUDIT TRAIL=  $$ (NUMBER, WHEN, WHO, WHERE, PROCESS, VERSION, WHAT)
(   1,<2026-05-02 01:29:39.981 -0500>,<nmr>,<biocwk-01407l.ad.wisc.edu>,<go4>,<TopSpin 3.6.2>,
      <created by zg
\tstarted at 2026-05-02 00:19:36.943 -0500,
\tPOWCHK enabled, PULCHK disabled,
       configuration hash MD5:
       92 D7 9E 89 28 7C 5A 14 42 79 B3 51 EB F8 B1 F7
       data hash MD5: 9590
       69 BD 39 94 1C DD B6 2E 90 9A 1F 6F 2E 19 6B 0D>)
(   2,<2026-05-26T09:34:30.911-0700>,<samso>,<Ordi_Sami>,<proc1d>,<TopSpin 4.5.0>,
      <Start of raw data processing
       em LB = 0 TDeff = 384 SI = 32K
       fid data hash MD5: 9590
       69 BD 39 94 1C DD B6 2E 90 9A 1F 6F 2E 19 6B 0D
       data hash MD5: 32K
       EE 6F 85 C6 1B 26 0D 18 CA E4 70 53 6E F7 CD 67>)
(   3,<2026-05-26T09:34:31.079-0700>,<samso>,<Ordi_Sami>,<proc1d>,<TopSpin 4.5.0>,
      <ft FT_mod = 6 PKNL = 1 SI = 32K
       data hash MD5: 32K
       9D C9 62 5B 3D 93 16 A3 18 F5 39 0D 80 26 B3 45>)
(   4,<2026-05-26T09:34:31.298-0700>,<samso>,<Ordi_Sami>,<proc1d>,<TopSpin 4.5.0>,
      <apk
       data hash MD5: 32K
       CD 74 6F 26 9F 7D 6E B9 BE F3 80 67 C6 C1 36 93>)
(   5,<2026-05-26T09:34:31.557-0700>,<samso>,<Ordi_Sami>,<proc1d>,<TopSpin 4.5.0>,
      <abs n ABSG = 5
       data hash MD5: 32K
       62 68 DC EF 5B 00 76 8A 4C 5B E1 7E 03 84 67 C4>)
##END=

$$ hash MD5
$$ 4B E7 60 55 3A 9D F8 1A 3C 8F DB 2F 70 DE C6 D2
"""

#: 2702/pdata/1/auditp.txt: the efp form, an acquisition record and a user
#: comment (the 3.6.2 -> 4.5.0 hand-over), plus an earlier superseded chain
AUDIT_2702 = """##TITLE= Audit trail, TopSpin 3.6.2
##JCAMPDX= 5.01
##ORIGIN= Bruker BioSpin GmbH
##OWNER= nmr
$$ /data/soudani_data/04272026_P1-Bi1-12_SS_ALP/2702/pdata/1/auditp.txt
##AUDIT TRAIL=  $$ (NUMBER, WHEN, WHO, WHERE, PROCESS, VERSION, WHAT)
(   1,<2026-04-29 15:59:19.361 -0500>,<nmr>,<biocwk-01407l.ad.wisc.edu>,<go4>,<TopSpin 3.6.2>,
      <acquisition in progress>)
(   2,<2026-04-29 16:00:00.000 -0500>,<nmr>,<biocwk-01407l.ad.wisc.edu>,<proc1d>,<TopSpin 3.6.2>,
      <Start of raw data processing
       efp LB = 0 FT_mod = 6 PKNL = 1 PHC0 = 269.3149 PHC1 = 36.40312 TDeff = 256 SI = 32K
       data hash MD5: 32K
       00 31 1D BE 64 C7 B8 97 1A C9 24 83 49 74 0F 2B>)
(   3,<2026-04-29 16:01:00.000 -0500>,<samso>,<Ordi_Sami>,<cprserver>,<TopSpin 4.5.0>,
      <changed parameter F1 TDeff to 512>)
(   4,<2026-04-29 16:02:17.249 -0500>,<nmr>,<biocwk-01407l.ad.wisc.edu>,<proc1d>,<TopSpin 3.6.2>,
      <Start of raw data processing
       efp LB = 0 FT_mod = 6 PKNL = 1 PHC0 = 269.3149 PHC1 = 36.40312 TDeff = 512 SI = 32K
       data hash MD5: 32K
       FA 31 1D BE 64 C7 B8 97 1A C9 24 83 49 74 0F 2B>)
$$ C:/Users/x/DATA/2026-5/04272026_P1-Bi1-12_SS_ALP/2702/pdata/1/auditp.txt
(   5,<2026-05-26T08:52:29.057-0700>,<samso>,<Ordi_Sami>,<audit>,<TopSpin 4.5.0>,
      <user comment:
       TITLE=27Al failed, MAS stopped
       P1(90)=3.125; 11deg tip = 0.375
       Rotor RS2427418
       MASR 20kHz
       ALPTITLE_END
       data hash MD5: 32K
       FA 31 1D BE 64 C7 B8 97 1A C9 24 83 49 74 0F 2B>)
##END=
"""


def _rule(key):
    return next(r for r in C.KEY_RULES if r.key == key)


# ---------------------------------------------------------------- reading
def test_read_params_flat_keys_from_synthetic_jcamp(tmp_path):
    p1r = _fake_expno(tmp_path, "01192026_SR31649_Base0Ca_SS_ALP")
    p = C.read_params(p1r)
    assert p is not None
    assert p.acq["D1"] == 12.5 and p.acq["P1"] == 0.425 and p.acq["PLW1"] == 100.0
    assert p.acq["PLW1_unit"] == "W"
    assert p.acq["NS"] == 256 and p.acq["TD"] == 7988 and p.acq["PULPROG"] == "zg"
    assert p.acq["PROBHD"] == "Z1 3.2mm" and p.acq["NUC1"] == "11B"
    assert p.acq["SFO1"] == pytest.approx(192.430693)
    assert p.proc["WDW"] == 1 and p.proc["LB"] == 0 and p.proc["TDeff"] == 1024
    assert p.proc["SI"] == 32768 and p.proc["PH_mod"] == 1
    assert p.proc["PHC0"] == pytest.approx(-160.0) and p.proc["PHC1"] == pytest.approx(29.9)
    assert p.proc["SF"] == pytest.approx(192.431528) and p.proc["PKNL"] is True
    assert p.procno == 1 and p.sample == "01192026_SR31649_Base0Ca_SS_ALP"
    assert Path(p.expno) == tmp_path / p.sample / "24"
    assert p.title.startswith("11B with short tip angle")
    assert p.has_fid is False and p.commands == []
    # an EXPNO folder, a pdata folder and a fid path all resolve to the same
    assert C.read_params(tmp_path / p.sample / "24").expno == p.expno
    assert C.read_params(p1r.parent).expno == p.expno
    with_fid = _fake_expno(tmp_path, "S2", with_fid=True)
    q = C.read_params(with_fid)
    assert q.has_fid is True
    assert C.read_params(Path(q.expno) / "fid").expno == q.expno
    # TDeff 0 (or >= TD) means the whole fid
    r = C.read_params(_fake_expno(tmp_path, "S3", tdeff=0))
    assert _rule("TDeff").display(r.proc["TDeff"], r) == "all"
    assert _rule("TDeff").display(1024, p) == "1024"
    # an older acqus carries PL (dB) and no PLW (W)
    s = C.read_params(_fake_expno(tmp_path, "S4", plw1=None, pl1=-3.0))
    assert s.acq["PLW1"] == -3.0 and s.acq["PLW1_unit"] == "dB"
    assert _rule("PLW1").display(s.acq["PLW1"], s).endswith("dB")
    assert _rule("PLW1").display(p.acq["PLW1"], p) == "100 W"
    # no procs yet (freshly acquired): proc rows empty, still a Bruker spectrum
    t = C.read_params(_fake_expno(tmp_path, "S5", with_procs=False, with_1r=False,
                                  with_fid=True))
    assert t is not None and t.proc == {} and t.has_fid


def test_non_bruker_inputs_degrade_to_none(tmp_path):
    csvp = tmp_path / "a.csv"
    csvp.write_text("# nucleus = 11B\n0 1\n1 2\n", encoding="utf-8")
    rec = tmp_path / "a.recipe.json"
    rec.write_text(json.dumps({"sample": "x"}), encoding="utf-8")
    bare = tmp_path / "folder"
    bare.mkdir()
    assert C.read_params(csvp) is None
    assert C.read_params(rec) is None
    assert C.read_params(bare) is None
    assert C.read_params(tmp_path / "missing" / "1r") is None
    cmp = C.compare([None, None])
    assert cmp.reports == [] and cmp.level == "none" and cmp.summary() == ""
    assert C.compare([]).level == "none"
    assert C.compare([None]).level == "none"
    a = _fake_expno(tmp_path, "A")
    b = _fake_expno(tmp_path, "B", lb=100)
    before = bruker.snapshot(tmp_path)
    pa, pb = C.read_params(a), C.read_params(b)
    cmp = C.compare([pa, None, pb], ["A", "csv", "B"])
    bruker.verify_untouched(tmp_path, before)          # nothing written
    assert cmp.n_bruker == 2 and cmp.level == "check"
    lb = cmp.report("LB")
    assert lb.display == ["0", "—", "100"]
    assert lb.deviants == [2] and cmp.signature(1) == "" and cmp.details(1) == []
    assert C.caveat_note(cmp, 1) == ""
    assert all(1 not in r.deviants for r in cmp.reports)
    # one Bruker row only -> nothing to compare
    assert C.compare([pa, None]).level == "none"


# ---------------------------------------------------------------- the majority rule
def test_compare_flags_shape_deviants_against_the_majority(tmp_path):
    lbs = [0, 0, 100, 100, 100]
    nss = [64, 128, 256, 256, 256]
    ps = [C.read_params(_fake_expno(tmp_path, f"S{k}", lb=lb, ns=ns, phc0=-160.0 - k))
          for k, (lb, ns) in enumerate(zip(lbs, nss))]
    cmp = C.compare(ps, [f"S{k}" for k in range(5)])
    lb = cmp.report("LB")
    assert lb.majority == 100 and lb.deviants == [0, 1] and lb.level == "check"
    assert "0 / 100 Hz" in lb.message and "(EM)" in lb.message
    assert not lb.tie and not lb.varies
    assert cmp.report("PHC0").level == "info" and cmp.report("NS").level == "info"
    assert cmp.report("PHC0").varies and cmp.report("NS").deviants == [0, 1]
    assert {r.key for r in cmp.flagged()} == {"LB"}
    assert cmp.signature(0) == "LB 0 Hz" and cmp.signature(2) == ""
    assert cmp.details(0) == ["LB 0 Hz — series majority 100 Hz"]
    s = cmp.summary()
    assert s.startswith("⚠") and "LB 0 / 100 Hz" in s and "processed differently" in s
    assert "NS 64 / 128 / 256" in s          # an info key the summary still names
    assert "PHC0" not in s
    assert cmp.any_flagged and cmp.level == "check"
    assert C.caveat_note(cmp, 0) == ("acquired/processed differently from the series "
                                     "majority: LB 0 Hz (majority 100 Hz)")
    assert C.caveat_note(cmp, 2) == ""


def test_varies_tie_and_two_row_rules(tmp_path):
    d1s = [12.5, 14, 18, 29, 36]
    ps = [C.read_params(_fake_expno(tmp_path, f"D{k}", d1=d)) for k, d in enumerate(d1s)]
    cmp = C.compare(ps)
    d1 = cmp.report("D1")
    assert d1.varies and d1.deviants == [] and d1.flagged
    assert "D1 varies 12.5–36 s" in cmp.summary()
    assert all(cmp.signature(k) == "" for k in range(5))     # no per-spectrum flag
    # 5 % relative: 12.5 and 12.9 are the same recycle delay
    same = C.compare([C.read_params(_fake_expno(tmp_path, "E0", d1=12.5)),
                      C.read_params(_fake_expno(tmp_path, "E1", d1=12.9))])
    assert not same.report("D1").differs
    # an exact 2-vs-2 tie takes the first spectrum's value and says so
    tie = C.compare([C.read_params(_fake_expno(tmp_path, f"T{k}", lb=lb))
                     for k, lb in enumerate([100, 100, 0, 0])])
    lb = tie.report("LB")
    assert lb.majority == 100 and lb.tie and lb.deviants == [2, 3]
    assert "tie" in lb.message
    # two rows: the second is the deviant
    two = C.compare([C.read_params(_fake_expno(tmp_path, "U0", lb=0)),
                     C.read_params(_fake_expno(tmp_path, "U1", lb=100))])
    assert two.report("LB").deviants == [1] and two.signature(1) == "LB 100 Hz"
    # three identical members: alike
    ok = C.compare([C.read_params(_fake_expno(tmp_path, f"I{k}")) for k in range(3)])
    assert ok.flagged() == [] and ok.level == "ok" and not ok.any_flagged
    s = ok.summary()
    assert s.startswith("✓") and "3 spectra" in s
    for bit in ("EM 0 Hz", "TDeff 1024", "SI 32k", "NS 256", "D1 12.5 s"):
        assert bit in s, s


def test_bad_level_and_sf_ppm_rule(tmp_path):
    b = C.read_params(_fake_expno(tmp_path, "B11", nuc="11B"))
    al = C.read_params(_fake_expno(tmp_path, "Al27", nuc="27Al", sfo1=156.28, bf1=156.28,
                                   sf=156.2816))
    cmp = C.compare([b, al])
    assert cmp.level == "bad"
    nuc = cmp.report("NUC1")
    assert nuc.level == "bad" and nuc.deviants == [1] and "mixed nuclei" in nuc.message
    assert "11B" in nuc.message and "27Al" in nuc.message
    sfo = cmp.report("SFO1")
    assert sfo.level == "bad" and sfo.deviants == [1]
    assert "Larmor frequencies" in sfo.message
    assert cmp.summary().startswith("⚠") and "shared model" in cmp.summary()
    # the same field: SFO1 differs by less than 5 % -> not flagged, info
    close = C.compare([C.read_params(_fake_expno(tmp_path, "F0", sfo1=192.430693)),
                       C.read_params(_fake_expno(tmp_path, "F1", sfo1=192.5))])
    assert close.report("SFO1").level == "info" and not close.report("SFO1").differs
    # referencing: 2.5 ppm apart is flagged in ppm, 0.001 ppm is not
    sf = C.compare([C.read_params(_fake_expno(tmp_path, "R0", sf=242.7925, sfo1=242.792156,
                                              bf1=242.792156, nuc="31P")),
                    C.read_params(_fake_expno(tmp_path, "R1", sf=242.7919, sfo1=242.792156,
                                              bf1=242.792156, nuc="31P"))])
    r = sf.report("SF")
    assert r.level == "check" and r.deviants == [1]
    assert "ppm" in r.message and "Referencing audit" in r.message
    assert "referenced differently" in sf.summary()
    fine = C.compare([C.read_params(_fake_expno(tmp_path, "Q0", sf=242.7925)),
                      C.read_params(_fake_expno(tmp_path, "Q1", sf=242.7925 * (1 + 1e-9)))])
    assert not fine.report("SF").differs
    # GB only for GM, SSB only for SINE/QSINE; named codes in the messages
    em_gm = C.compare([C.read_params(_fake_expno(tmp_path, "W0", wdw=1, lb=50)),
                       C.read_params(_fake_expno(tmp_path, "W1", wdw=2, lb=-20, gb=0.1))])
    assert em_gm.report("WDW").message.startswith("WDW EM / GM")
    assert em_gm.report("GB").display == ["—", "0.1"]
    assert em_gm.report("SSB") is None
    q = C.compare([C.read_params(_fake_expno(tmp_path, "X0", wdw=4, ssb=2)),
                   C.read_params(_fake_expno(tmp_path, "X1", wdw=4, ssb=3))])
    assert q.report("SSB").deviants == [1] and q.report("GB") is None
    ph = C.compare([C.read_params(_fake_expno(tmp_path, "P0", ph_mod=1)),
                    C.read_params(_fake_expno(tmp_path, "P1", ph_mod=2))])
    assert ph.report("PH_mod").message.startswith("PH_mod pk / mc")


# ---------------------------------------------------------------- auditp
def test_parse_auditp_extracts_processing_commands(tmp_path):
    f = tmp_path / "auditp.txt"
    f.write_text(AUDIT_3104, encoding="utf-8")
    assert C.parse_auditp(f) == ["em LB = 0 TDeff = 384 SI = 32K",
                                 "ft FT_mod = 6 PKNL = 1 SI = 32K", "apk",
                                 "abs n ABSG = 5"]
    assert C.parse_auditp(tmp_path) == C.parse_auditp(f)      # the pdata folder
    f.write_text(AUDIT_2702, encoding="utf-8")
    # the LAST 'Start of raw data processing' chain; go4 / cprserver / audit excluded
    assert C.parse_auditp(f) == ["efp LB = 0 FT_mod = 6 PKNL = 1 PHC0 = 269.3149 "
                                 "PHC1 = 36.40312 TDeff = 512 SI = 32K"]
    assert C.parse_auditp(tmp_path / "nowhere.txt") == []
    f.write_text(AUDIT_3104[:len(AUDIT_3104) // 2], encoding="utf-8")
    assert isinstance(C.parse_auditp(f), list)                # truncated: no raise
    (tmp_path / "junk.txt").write_bytes(b"\xff\xfe\x00garbage")
    assert C.parse_auditp(tmp_path / "junk.txt") == []
    # the commands row rides along as an info row of the comparison
    a = C.read_params(_fake_expno(tmp_path, "A", audit=AUDIT_3104))
    b = C.read_params(_fake_expno(tmp_path, "B", audit=AUDIT_2702))
    assert a.commands[0].startswith("em LB = 0") and b.commands[0].startswith("efp")
    cmp = C.compare([a, b])
    row = cmp.report("commands")
    assert row.level == "info" and row.differs and row not in cmp.flagged()
    assert "apk" in row.display[0] and "efp" in row.display[1]
    assert "commands" in cmp.to_text(False) and "commands" not in cmp.to_text(True)


# ---------------------------------------------------------------- ops mapping
def test_ops_for_maps_topspin_codes_and_template(tmp_path):
    from larmor import qcpmg

    p = C.read_params(_fake_expno(tmp_path, "S", lb=100, tdeff=1024, td=7988))
    ops = C.ops_for(p)
    names = [o["op"] for o in ops]
    assert names == ["tdeff", "em", "zf", "ft", "phase"]
    assert ops[0] == {"op": "tdeff", "points": 512}        # TopSpin real points, halved
    assert ops[1] == {"op": "em", "lb_hz": 100.0}
    assert ops[2] == {"op": "zf", "si": 32768}
    off = qcpmg.carrier_ppm({"larmor_MHz": 192.430693, "sf_MHz": 192.431528})[0]
    assert ops[3]["offset_ppm"] == pytest.approx(off) and off == pytest.approx(-4.34, abs=0.01)
    assert ops[4] == {"op": "phase", "p0": 160.0, "p1": pytest.approx(29.9),
                      "pivot_frac": 1.0}
    assert proc.chain_start_domain(ops) == "time"
    assert all(o["op"] in proc.OPS for o in ops)
    # LB 0 -> no window step; TDeff 0 -> no tdeff step; FCOR 1 -> no fcor
    q = C.read_params(_fake_expno(tmp_path, "Q", lb=0, tdeff=0))
    assert [o["op"] for o in C.ops_for(q)] == ["zf", "ft", "phase"]
    f = C.read_params(_fake_expno(tmp_path, "F", fcor=0.5))
    assert {"op": "fcor", "factor": 0.5} in C.ops_for(f)
    # window codes
    gm = C.read_params(_fake_expno(tmp_path, "G", wdw=2, lb=-20, gb=0.1))
    assert {"op": "gm", "lb_hz": -20.0, "gb": 0.1} in C.ops_for(gm)
    qs = C.read_params(_fake_expno(tmp_path, "QS", wdw=4, ssb=2))
    assert {"op": "sine", "ssb": 2.0, "power": 2} in C.ops_for(qs)
    sn = C.read_params(_fake_expno(tmp_path, "SN", wdw=3, ssb=0))
    assert {"op": "sine", "ssb": 2.0, "power": 1} in C.ops_for(sn)
    tr = C.read_params(_fake_expno(tmp_path, "TR", wdw=9, lb=0))
    assert {"op": "traf", "lb_hz": 10.0} in C.ops_for(tr)
    with pytest.raises(ValueError, match="WDW 5"):
        C.ops_for(C.read_params(_fake_expno(tmp_path, "TP", wdw=5)))
    # phase modes
    mc = C.read_params(_fake_expno(tmp_path, "MC", ph_mod=2))
    assert C.ops_for(mc)[-1] == {"op": "magnitude"}
    no = C.read_params(_fake_expno(tmp_path, "NO", ph_mod=0))
    assert C.ops_for(no)[-1]["op"] == "ft"
    assert C.ops_for(p, phase="auto")[-1] == {"op": "autophase"}
    assert C.ops_for(p, phase="none")[-1]["op"] == "ft"
    with pytest.raises(ValueError):
        C.ops_for(p, phase="sideways")
    # a template overrides only the shape keys; phase and offset stay the spectrum's
    t = {"wdw": 0, "lb_hz": 0.0, "gb": 0.0, "ssb": 0.0, "tdeff": 1024, "si": 32768,
         "fcor": 1.0}
    ops_t = C.ops_for(p, t)
    assert [o["op"] for o in ops_t] == ["tdeff", "zf", "ft", "phase"]
    assert ops_t[-1]["p0"] == 160.0 and ops_t[2]["offset_ppm"] == pytest.approx(off)
    assert C.ops_for(gm, t)[-1]["p0"] == 160.0        # its own PHC0, not p's
    assert C.template_from(p) == {"wdw": 1, "lb_hz": 100.0, "gb": 0.0, "ssb": 0.0,
                                  "tdeff": 1024, "si": 32768, "fcor": 1.0}
    same = C.compare([C.read_params(_fake_expno(tmp_path, f"I{k}")) for k in range(3)])
    assert C.majority_template(same) == {"wdw": 1, "lb_hz": 0.0, "gb": 0.0, "ssb": 0.0,
                                         "tdeff": 1024, "si": 32768, "fcor": 1.0}
    d = C.describe_template(C.majority_template(same))
    assert "TDeff 1024" in d and "EM 0 Hz" in d and "SI 32k" in d and "own PHC0" in d
    assert "no window" in C.describe_template(t) and "autophase" in C.describe_template(t, "auto")
    assert "no phase" in C.describe_template({**t, "tdeff": 0}, "none")
    assert "TDeff all" in C.describe_template({**t, "tdeff": 0})


def test_to_text_to_csv_and_without(tmp_path):
    lbs = [0, 0, 100, 100, 100]
    labels = [f"S{k}" for k in range(5)]
    ps = [C.read_params(_fake_expno(tmp_path, lab, lb=lb, phc0=-160.0 - k))
          for k, (lab, lb) in enumerate(zip(labels, lbs))]
    cmp = C.compare(ps, labels)
    short = cmp.to_text(True)
    lines = [ln for ln in short.splitlines() if ln and not ln.startswith("*")]
    assert len(lines) == 2 and lines[1].split()[1] == "LB"      # header + LB
    assert "PHC0" not in short and "* differs" in short
    assert "LB (Hz)" in short and "0 *" in short
    full = cmp.to_text(False)
    assert "NS" in full and "PHC0" in full          # (no auditp here: no commands row)
    rows = cmp.rows()
    assert rows[0] == ["group", "key", "label", "level", "majority", *labels]
    lb_row = next(r for r in rows if r[1] == "LB")
    assert lb_row == ["processing", "LB", "LB (Hz)", "check", "100", "0", "0", "100",
                      "100", "100"]
    out = tmp_path / "cmp.csv"
    cmp.to_csv(out)
    with open(out, newline="", encoding="utf-8") as f:
        got = list(csv.reader(f))
    assert got == [[str(c) for c in r] for r in rows]
    less = cmp.without(C.PROC_SHAPE_KEYS)
    assert less.report("LB") is None and less.report("TDeff") is None
    assert less.report("D1") is not None and less.report("PHC0") is not None
    assert less.level == "ok" and less.signature(0) == ""
    assert cmp.report("LB") is not None                          # the original intact
    assert "no differences" in less.to_text(True)


def test_reprocess_requires_a_fid_and_replays_through_loader(tmp_path, monkeypatch):
    from larmor import loader

    p = C.read_params(_fake_expno(tmp_path, "S", lb=100))
    ops = C.ops_for(p)
    with pytest.raises(ValueError, match="no raw fid"):
        C.reprocess(p, ops)
    q = C.read_params(_fake_expno(tmp_path, "T", lb=100, with_fid=True))
    seen = {}

    def fake_apply(recipe, ppm, amp, source_path=None):
        seen["recipe"], seen["source"] = recipe, source_path
        return np.linspace(-10, 10, 8), np.ones(8), ["replayed 5 processing step(s)"]

    monkeypatch.setattr(loader, "apply_processing", fake_apply)
    ppm, amp, notes = C.reprocess(q, ops)
    assert seen["source"] == q.expno and bruker.is_expno(seen["source"])
    assert seen["recipe"].processing == ops and seen["recipe"].processing_from_raw
    assert seen["recipe"].nucleus == "11B"
    assert seen["recipe"].larmor_frequency_MHz == pytest.approx(192.430693)
    assert ppm.size == 8 and amp.size == 8 and notes == ["replayed 5 processing step(s)"]


# ---------------------------------------------------------------- real data
def test_law_ca_11b_series_real_data():
    require(LAW_CA_11B[0])
    paths = [e / "pdata" / "1" / "1r" for e in LAW_CA_11B]
    before = bruker.snapshot(LAW_CA_11B[0])
    ps = [C.read_params(p) for p in paths]
    bruker.verify_untouched(LAW_CA_11B[0], before)
    assert all(p is not None and p.has_fid for p in ps)
    assert ps[0].commands and ps[0].commands[0].startswith("efp")
    cmp = C.compare(ps, ["Base0Ca", "Base1Ca", "Base2Ca", "Base3Ca", "Base4Ca"])
    assert cmp.level == "check"
    lb = cmp.report("LB")
    assert lb.values == [0, 0, 100, 100, 100] and lb.majority == 100
    assert lb.deviants == [0, 1]
    assert cmp.report("WDW").display == ["EM"] * 5 and not cmp.report("WDW").differs
    assert not cmp.report("TDeff").differs and cmp.report("TDeff").values == [1024] * 5
    d1 = cmp.report("D1")
    assert d1.varies and d1.values == [12.5, 14, 18, 29, 36]
    for key in ("NS", "P1", "PLW1", "RG", "SF", "PULPROG", "SI", "TD", "SW_h", "PROBHD"):
        assert not cmp.report(key).flagged, key
    assert cmp.report("PHC0").level == "info" and cmp.report("PHC0").differs
    assert {r.key for r in cmp.flagged()} == {"LB", "D1"}
    assert cmp.signature(0) == "LB 0 Hz" and cmp.signature(4) == ""
    s = cmp.summary()
    assert "LB 0 / 100 Hz (EM)" in s and "D1 varies 12.5–36 s" in s
    assert C.majority_template(cmp)["lb_hz"] == 100.0
    assert C.majority_template(cmp)["tdeff"] == 1024


def test_reprocess_2702_from_fid_matches_the_1r():
    """The two conventions the whole reprocess path rests on, pinned on the
    27Al 2702 fid/1r pair: TopSpin's stored phase replays as
    phase(p0=-PHC0, p1=+PHC1, pivot 1.0) -- r = +0.9994 against the 1r (the
    opposite sign gives r = -0.83) -- and TDeff counts real points, so the
    tdeff step takes TDeff // 2 complex points: unhalved, the peak lands
    0.9 ppm off the 1r's instead of 0.4 ppm (the FWHM alone does not
    discriminate: 31.3 ppm either way, 1r 31.5 ppm)."""
    require(BRUKER_FID)
    require(BRUKER_1R)
    p = C.read_params(BRUKER_1R)
    assert p.has_fid and p.proc["TDeff"] == 512 and p.proc["PHC0"] == pytest.approx(269.3149)
    ops = C.ops_for(p)
    assert ops[0] == {"op": "tdeff", "points": 256}
    assert ops[-1] == {"op": "phase", "p0": pytest.approx(-269.3149),
                       "p1": pytest.approx(36.40312), "pivot_frac": 1.0}
    ref = bruker.read(BRUKER_1R)
    x1, y1 = np.asarray(ref.axes[0].values), np.asarray(ref.data, float)
    peak1 = float(x1[np.argmax(y1)])
    sel = (x1 > peak1 - 60) & (x1 < peak1 + 60)

    def fwhm(x, y):
        i = int(np.argmax(y))
        half = y[i] / 2
        lo, hi = i, i
        while lo > 0 and y[lo] > half:
            lo -= 1
        while hi < y.size - 1 and y[hi] > half:
            hi += 1
        return abs(float(x[hi] - x[lo]))

    def against_1r(ops_):
        ppm, amp, notes = C.reprocess(p, ops_)
        assert amp.size == 32768 and notes and "replayed" in notes[0]
        yi = np.interp(x1, ppm, amp)
        r = float(np.corrcoef(yi[sel], y1[sel])[0, 1])
        return r, float(ppm[np.argmax(amp)]), fwhm(ppm, amp)

    r, peak, w = against_1r(ops)
    assert r > 0.99, r                                        # measured 0.9994
    assert abs(peak - peak1) < 0.5, (peak, peak1)             # measured 0.43 ppm
    assert abs(w - fwhm(x1, y1)) < 0.1 * fwhm(x1, y1)
    r_auto, peak_auto, _ = against_1r(C.ops_for(p, phase="auto"))
    assert r_auto > 0.98 and abs(peak_auto - peak1) < 0.5      # the fallback
    wrong = [dict(o, p0=-o["p0"], p1=-o["p1"]) if o["op"] == "phase" else o for o in ops]
    assert against_1r(wrong)[0] < 0                           # measured -0.83
    unhalved = [dict(o, points=512) if o["op"] == "tdeff" else o for o in ops]
    assert abs(against_1r(unhalved)[1] - peak1) > 0.5         # measured 0.92 ppm
