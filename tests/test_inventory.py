"""Session inventory (larmor.inventory) and the sample identity it rests on
(larmor.io.scan): the folder-name rule on every real folder shape, the
extended acqus reader, roles and the two-tier production pick on a synthetic
month, every flag kind, the SR join, the outputs, the CLI -- plus the real
NMRFAM / MagLab sessions as anchors behind conftest.require()."""
from pathlib import Path

import pytest

from conftest import BRUKER_1R, DATA_ROOT, LAW_CA_11B, MAGLAB_35CL, require
from larmor import inventory as I
from larmor import referencing as R
from larmor.io import scan

SF_H, BF1_H = 599.771479328723, 599.772
BF1_AL, BF1_P, BF1_NA = 156.281744, 242.792156, 158.6
T0 = 1_777_400_000


# ---------------------------------------------------------------- builders
def _jcamp(path: Path, **params):
    lines = ["##TITLE= Parameter file, TopSpin 4.1", "##JCAMPDX= 5.0",
             "##DATATYPE= Parameter Values", "##ORIGIN= Bruker BioSpin GmbH",
             "##OWNER= nmr"]
    for k, v in params.items():
        lines.append(f"##${k}= {v}")
    lines.append("##END=")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _expno(root, sample, expno, nuc, *, ns=1, d1=1.0, pulprog="zg", title="",
           has_1r=True, twod=False, popt=False, bf1=BF1_AL, sf=None, date=0):
    """One synthetic EXPNO: acqus (NS, DATE, the D array over two lines),
    pdata/1/title, an empty 1r (or acqu2s + ser + 2rr for a 2D), procs when
    ``sf`` is given, popt.array on request."""
    e = Path(root) / sample / str(expno)
    d = [0.2, float(d1)] + [0.0] * 62
    d_text = "(0..63)\n" + " ".join(f"{v:g}" for v in d[:32]) + "\n" \
        + " ".join(f"{v:g}" for v in d[32:])
    _jcamp(e / "acqus", NUC1=f"<{nuc}>", BF1=bf1, SFO1=bf1, O1=0,
           DATE=date or (T0 + int(expno)), PROBHD="<Z1 3.2mm>",
           PULPROG=f"<{pulprog}>", PARMODE=1 if twod else 0, NS=ns, D=d_text)
    pd = e / "pdata" / "1"
    pd.mkdir(parents=True, exist_ok=True)
    (pd / "title").write_text(title, encoding="utf-8")
    if sf is not None:
        _jcamp(pd / "procs", SF=sf, SI=1024, SR=(sf - bf1) * 1e6)
    if twod:
        _jcamp(e / "acqu2s", NUC1=f"<{nuc}>")
        (e / "ser").write_bytes(b"\0" * 8)
        (pd / "2rr").write_bytes(b"\0" * 8)
    elif has_1r:
        (pd / "1r").write_bytes(b"\0" * 8)
    (e / "fid").write_bytes(b"\0" * 8)
    if popt:
        (e / "popt.array").write_text("popt", encoding="utf-8")
    return e


A = "04272026_P5-Bi8-12_SS_ALP"
B = "04272026_P1-Bi1-12_SS_ALP"
C = "04272026_RS40339_P5-Bi1-12"
D = "05082026_P5-Bi8-12_SS_ALP"
E = "01192026_SR31649_Base0Ca_SS_ALP"


def _month(tmp_path) -> Path:
    """The synthetic 2026-05: folders A-E plus a TopSpin ~TEMP."""
    root = tmp_path / "2026-05"
    _expno(root, A, 2701, "27Al", ns=2, pulprog="satrect1", twod=True, title="27Al satrec")
    _expno(root, A, 2702, "27Al", ns=512, d1=5, title="27Al")
    _expno(root, A, 2704, "27Al", ns=512, d1=25, title="27Al")
    _expno(root, A, 2799, "27Al", ns=24, d1=3, title="27Al")
    for k in (3102, 3103, 3104):
        _expno(root, A, k, "31P", ns=14, d1=300, bf1=BF1_P, title="31P")
    _expno(root, A, 3112, "31P", ns=2, pulprog="satrect1", twod=True, bf1=BF1_P,
           title="31P satrec")
    _expno(root, A, 3113, "31P", ns=28, d1=800, bf1=BF1_P, title="31P")
    _expno(root, A, 3114, "31P", ns=140, d1=300, bf1=BF1_P,
           title="31P\nP1(90)=3.750\n\nRotor RS2427418\nMASR 20 kHz\nALP")
    _expno(root, B, 2702, "27Al", ns=56, title="27Al failed, MAS stopped")
    _expno(root, B, 2703, "27Al", ns=40, title="27Al")
    _expno(root, B, 2704, "27Al", ns=512, title="27Al")
    _expno(root, C, 1, "1H", ns=1, bf1=BF1_H, title="1H spectrum for later referencing")
    _expno(root, D, 2702, "27Al", ns=512, has_1r=False, title="27Al")
    _expno(root, D, 2799, "27Al", ns=28, title="27Al")
    _expno(root, D, 11, "27Al", ns=96, pulprog="mp3qdfsz", popt=True, title="27Al 3Q DFS Z")
    _expno(root, D, 3102, "31P", ns=64, bf1=BF1_P, title="31P")
    _expno(root, E, 31, "27Al", ns=64, title="27Al zg power check")
    _expno(root, E, 33, "27Al", ns=512, title="27Al zg power check")
    _expno(root, E, 36, "27Al", ns=4, popt=True, title="27Al zg power check")
    _expno(root, "~TEMP", 1, "27Al", ns=1, title="")
    (root / "txt").mkdir()
    return root


def _maglab(tmp_path) -> Path:
    """An EXPNO-per-sample set: '35Cl_2025-12/<n>' with 'MM/DD/YYYY\\nSample X'."""
    root = tmp_path / "DATA" / "35Cl_2025-12"
    spec = [(1, "LAW3CL0CA", 40960), (2, "LAW3CL4CA", 40960), (3, "LAW3CL4CA", 307200),
            (10, "CS-3Cl", 1024), (15, "CS-3Cl", 307200)]
    for n, name, ns in spec:
        _expno(tmp_path / "DATA", "35Cl_2025-12", n, "35Cl", ns=ns, d1=0.1,
               pulprog="qcpmg_dec.ih", bf1=83.3, popt=True,
               title=f"12/09/2025\nSample {name}\n\nRotor-synchronized CPMG\n300k scans")
    return root


def _pick(inv, folder, key, nuc):
    return inv.cell((folder, key), nuc)


def _row(inv, folder, expno):
    return next(r for r in inv.rows if r.folder == folder and r.expno == expno)


# ---------------------------------------------------------------- identity
@pytest.mark.parametrize("folder,key,date,rotor,suffix", [
    ("01192026_SR31649_Base0Ca_SS_ALP", "Base0Ca", "01192026", "SR31649", "_SS_ALP"),
    ("04272026_P5-Bi8-12_SS_ALP", "P5-Bi8-12", "04272026", "", "_SS_ALP"),
    ("05082026_P5-Bi8-12_SS_ALP", "P5-Bi8-12", "05082026", "", "_SS_ALP"),
    ("04272026_RS40339_P5-Bi1-12", "P5-Bi1-12", "04272026", "RS40339", ""),
    ("06052026_RS40175_LAW3Cl2Ca_SS_ALP", "LAW3Cl2Ca", "06052026", "RS40175", "_SS_ALP"),
    ("07062026_SR31649_2Ca12F_SS_ALP", "2Ca12F", "07062026", "SR31649", "_SS_ALP"),
    ("08272026_SR31649_2Ca-12F_ALP_SS", "2Ca-12F", "08272026", "SR31649", "_ALP_SS"),
    ("031172026_SR31649_Cryolite_SS_ALP", "Cryolite", "031172026", "SR31649", "_SS_ALP"),
    ("07062026_SR31648_NBF-SS_SS_ALP", "NBF-SS", "07062026", "SR31648", "_SS_ALP"),
    ("08312026_SR31602LS_11B27Al23Na19F_setup_ALP", "11B27Al23Na19F_setup", "08312026",
     "SR31602LS", "_ALP"),
    ("08272026_SR31649_Na75.12F_enrich_ALP_SS", "Na75.12F_enrich", "08272026", "SR31649",
     "_ALP_SS"),
    ("08272026_test", "test", "08272026", "", ""),
    ("07022026_Phoenix_1p6mm_setup", "Phoenix_1p6mm_setup", "07022026", "", ""),
    ("01232026_SR31648_Na3Al3Si3O8F8_regular_SS_ALP", "Na3Al3Si3O8F8_regular", "01232026",
     "SR31648", "_SS_ALP"),
    ("Cryolite", "Cryolite", "", "", ""),
    ("35Cl_2025-12", "35Cl_2025-12", "", "", ""),
])
def test_name_parts_strips_date_rotor_and_one_suffix(folder, key, date, rotor, suffix):
    n = scan.name_parts(folder)
    assert (n.key, n.date, n.rotor, n.suffix) == (key, date, rotor, suffix)
    assert n.folder == folder and n.from_folder == bool(date or rotor or suffix)
    assert scan.name_parts(f"C:/x/2026-05/{folder}").key == key   # a full path too
    assert scan.sample_key(folder) == key


def test_sample_name_layers_folder_month_title_sample_and_title(tmp_path):
    # (a) the folder convention -- the title is never consulted
    n = scan.sample_name(tmp_path / "2026-01" / E / "24",
                         "11B with short tip angle\n\nRotor SR31649")
    assert (n.key, n.source, n.title_first) == ("Base0Ca", "folder", "11B with short tip angle")
    assert n.rotor == "SR31649" and n.folder == E
    n = scan.sample_name(tmp_path / "2026-03b" / "03232026_P5-Bi8-12_SS_ALP" / "3104",
                         "ALP Restart 2026-03-30\n31P")
    assert (n.key, n.source) == ("P5-Bi8-12", "folder")
    # (b) a plain folder directly under a session month
    n = scan.sample_name(tmp_path / "2025-12" / "neg8A2O-0F" / "33", "27Al zg - spectrum")
    assert (n.key, n.source) == ("neg8A2O-0F", "month-folder")
    # (c) the MagLab "Sample X" title line
    n = scan.sample_name(tmp_path / "DATA" / "35Cl_2025-12" / "1", "12/09/2025\nSample LAW3CL0CA")
    assert (n.key, n.source, n.title_first) == ("LAW3CL0CA", "title-sample", "12/09/2025")
    assert scan.sample_name(tmp_path / "DATA" / "set" / "2", "x\nSample: CS-3Cl").key == "CS-3Cl"
    assert scan.sample_name(tmp_path / "DATA" / "set" / "3", "sample = NS3").key == "NS3"
    # (d) the title's first line (Tutorial 7 quotes '08/10/2026'), else the folder
    n = scan.sample_name(tmp_path / "DATA" / "81Br_2026-08" / "30", "08/10/2026\n81 Br")
    assert (n.key, n.source) == ("08/10/2026", "title")
    n = scan.sample_name(tmp_path / "foreign" / "bar" / "1", "")
    assert (n.key, n.source, n.title_first) == ("bar", "fallback", "")


def test_sample_label_and_disambiguate(tmp_path):
    root = _month(tmp_path)
    p = str(root / A / "3114" / "pdata" / "1" / "1r")
    assert scan.sample_label(p, {}) == "P5-Bi8-12"              # through sample_name
    assert scan.sample_label(p, {"sample": "31P", "nucleus": "31P"}) == "P5-Bi8-12"
    assert scan.sample_label(p, {"sample": "my glass"}) == "my glass"
    assert scan.sample_label("C:/v/x.fid/fid", {}) == "x"
    assert scan.sample_label(f"C:/nowhere/2026-05/{A}/3114/pdata/1/1r", {}) == "P5-Bi8-12"
    assert scan.sample_label("C:/nowhere/g0_raw.csv", {}) == "g0_raw"
    assert scan.sample_label("C:/nowhere/P1_31P.fxml", {}) == "P1_31P"
    a, d = f"C:/x/2026-05/{A}/3114/pdata/1/1r", f"C:/x/2026-05/{D}/3102/pdata/1/1r"
    assert scan.disambiguate(["P5-Bi8-12", "P5-Bi8-12"], [a, d]) == \
        ["P5-Bi8-12 (04272026)", "P5-Bi8-12 (05082026)"]
    same = [f"C:/x/2026-05/{A}/2702/pdata/1/1r", f"C:/x/2026-05/{A}/2704/pdata/1/1r"]
    assert scan.disambiguate(["P5-Bi8-12", "P5-Bi8-12"], same) == \
        ["P5-Bi8-12 · 2702", "P5-Bi8-12 · 2704"]
    # the same day twice: the folder name tells them apart
    c = f"C:/x/2026-05/{C}/1/pdata/1/1r"
    b1 = "C:/x/2026-05/04272026_P5-Bi1-12_SS_ALP/2702/pdata/1/1r"
    assert scan.disambiguate(["P5-Bi1-12", "P5-Bi1-12"], [c, b1]) == \
        [f"P5-Bi1-12 ({C})", "P5-Bi1-12 (04272026_P5-Bi1-12_SS_ALP)"]
    assert scan.disambiguate(["x", "y"], [a, d]) == ["x", "y"]
    assert scan.disambiguate(["x", "x"], [a]) == ["x", "x"]       # length mismatch: untouched


def test_read_experiment_carries_ns_d1_date_procs_popt(tmp_path):
    e = _expno(tmp_path / "m", "s", 3114, "31P", ns=140, d1=300, date=1777694973,
               popt=True, title="31P\nnote line")
    (e / "pdata" / "2").mkdir()
    (e / "pdata" / "notaproc").mkdir()
    info = scan.read_experiment(e)
    assert info.ns == 140 and info.date == 1777694973.0
    assert info.d1_s == pytest.approx(300.0)
    assert info.date_iso.startswith("2026-") and info.n_procs == 2 and info.has_popt
    assert info.title == "31P" and info.title_full == "31P\nnote line"
    bare = tmp_path / "m" / "s" / "1"
    _jcamp(bare / "acqus", NUC1="<27Al>", PULPROG="<zg>")
    b = scan.read_experiment(bare)
    assert (b.ns, b.d1_s, b.date, b.n_procs, b.has_popt, b.title_full) == (0, 0.0, 0.0, 0, False, "")
    assert b.date_iso == ""


def test_inventory_root_layouts(tmp_path):
    root = _month(tmp_path)
    for p in (root, root / A, root / A / "3114", root / A / "3114" / "pdata" / "1",
              root / A / "3114" / "pdata" / "1" / "1r", root / "~TEMP" / "1"):
        assert I.inventory_root(p) == root, p
    mag = _maglab(tmp_path)
    assert I.inventory_root(mag) == mag                       # an EXPNO-per-sample set
    assert I.inventory_root(mag / "1") == mag
    assert I.inventory_root(mag / "1" / "pdata" / "1" / "1r") == mag
    plain = tmp_path / "2025-12" / "Cryolite"
    _expno(tmp_path / "2025-12", "Cryolite", 16, "23Na", title="23Na")
    assert I.inventory_root(plain) == tmp_path / "2025-12"   # plain folder under a month


# ---------------------------------------------------------------- roles
@pytest.mark.parametrize("title,expected", [
    ("27Al zg power check", "setup"), ("23Na 1D for power optimization", "setup"),
    ("11B  non-quantitatiev", "setup"),
    ("19F with no fluorine in sample to check empty rotor signal", "setup"),
    ("27Al zg - background spectrum for 27Al", "setup"),
    ("1H spectrum for later referencing", "setup"),
    ("27Al failed, MAS stopped", "failed"), ("zg 29Si trash", "failed"),
    ("31P", ""), ("11B with short tip angle", ""),
    ("Rotor RS40175 with LAW3Cl3Ca for MAS QCPMG", ""),
    ("19F Hahn echo\n\nSR is set for IUPAC referencing from adamantane", ""),
    ("31P\nrun aborted after 3 h", "failed"),
])
def test_classify_title_words(title, expected):
    assert I.classify_title(title) == expected


@pytest.mark.parametrize("ns_frac,pick_2799", [(I.DEFAULT_NS_FRAC, False), (0.01, True)])
def test_roles_and_two_tier_production_pick_on_synthetic_month(tmp_path, ns_frac, pick_2799):
    inv = I.build(_month(tmp_path), ns_frac=ns_frac, sr_audit=False)
    assert inv.root == tmp_path / "2026-05"
    assert not any(r.folder.startswith("~") for r in inv.rows)
    assert len(inv.rows) == 21 and inv.n_folders() == 5
    role = {(r.folder, r.expno): r.role for r in inv.rows}
    assert role[(A, 2701)] == "arrayed" and role[(A, 3112)] == "arrayed"
    if pick_2799:
        assert role[(A, 2799)] == "production" and role[(A, 2704)] == "candidate"
    else:
        assert role[(A, 2799)] == "short" and "5 %" in _row(inv, A, 2799).reasons[0]
        assert role[(A, 2704)] == "production" and role[(A, 2702)] == "candidate"
        prod = _row(inv, A, 2704)
        assert any("highest EXPNO with pdata/1/1r is 2799" in s and "picked 2704" in s
                   for s in prod.reasons)
        assert all(role[(A, k)] == "short" for k in (3102, 3103, 3104, 3113))
        assert role[(B, 2703)] == "short"
    assert role[(A, 3114)] == "production"
    assert role[(B, 2702)] == "failed" and "failed" in _row(inv, B, 2702).reasons[0]
    assert role[(B, 2704)] == "production"
    assert role[(C, 1)] == "reference" and not _row(inv, C, 1).pick
    assert role[(D, 2702)] == "unprocessed"
    assert role[(D, 11)] == "setup"
    assert "2D pulse program acquired as 1D" in _row(inv, D, 11).reasons
    assert role[(D, 2799)] == "production"
    assert any("2702" in f and "no pdata/1/1r" in f for f in _row(inv, D, 2799).flags)
    assert role[(D, 3102)] == "production"
    # every row of block E is a 'power check': the pick falls to tier 2
    assert role[(E, 33)] == "production"
    assert any("power check" in s for s in _row(inv, E, 33).reasons)
    assert any("no unflagged spectrum" in s for s in _row(inv, E, 33).reasons)
    assert role[(E, 31)] == "setup"
    if not pick_2799:                       # 64 of 512 is short at the default fraction only
        assert any("12 %" in s for s in _row(inv, E, 31).reasons)
    assert role[(E, 36)] == "setup" and any("popt" in s for s in _row(inv, E, 36).reasons)
    assert set(inv.counts()) <= set(I.ROLES)


def test_grid_labels_picks_and_manual_override(tmp_path):
    inv = I.build(_month(tmp_path), sr_audit=False)
    assert inv.nuclei() == ["27Al", "31P", "1H"]
    labels = inv.labels()
    assert labels[(A, "P5-Bi8-12")] == "P5-Bi8-12 (04272026)"
    assert labels[(D, "P5-Bi8-12")] == "P5-Bi8-12 (05082026)"
    assert labels[(B, "P1-Bi1-12")] == "P1-Bi1-12" and labels[(C, "P5-Bi1-12")] == "P5-Bi1-12"
    assert inv.samples() == [(E, "Base0Ca"), (B, "P1-Bi1-12"), (C, "P5-Bi1-12"),
                             (A, "P5-Bi8-12"), (D, "P5-Bi8-12")]
    picks = inv.picks("31P")
    assert [r.openable for r in picks] == [str(tmp_path / "2026-05" / A / "3114" / "pdata" / "1" / "1r"),
                                           str(tmp_path / "2026-05" / D / "3102" / "pdata" / "1" / "1r")]
    g = inv.grid()
    assert g[((A, "P5-Bi8-12"), "1H")] is None and g[((A, "P5-Bi8-12"), "27Al")].expno == 2704
    assert inv.cell((C, "P5-Bi1-12"), "1H").role == "reference"   # shown, never picked
    assert len(inv.picks()) == 6 and inv.counts()["production"] == 6
    # manual override: exclusive per block, recorded as a reason
    inv.set_pick((A, "P5-Bi8-12"), "27Al", 2702)
    block = inv.candidates((A, "P5-Bi8-12"), "27Al")
    assert [r.expno for r in block if r.pick] == [2702]
    assert _row(inv, A, 2702).role == "production" and "chosen by user" in _row(inv, A, 2702).reasons
    assert _row(inv, A, 2704).role == "candidate" and not _row(inv, A, 2704).pick
    assert inv.grid()[((A, "P5-Bi8-12"), "27Al")].expno == 2702
    inv.set_pick((A, "P5-Bi8-12"), "27Al", None)                # clears the block
    assert inv.grid()[((A, "P5-Bi8-12"), "27Al")] is None
    with pytest.raises(ValueError):
        inv.set_pick((A, "P5-Bi8-12"), "27Al", 9999)


def test_title_flags_rotor_sibling_nucleus_kind_duplicate(tmp_path):
    root = tmp_path / "2026-06"
    L2, L3 = "06052026_RS40175_LAW3Cl2Ca_SS_ALP", "06052026_RS40175_LAW3Cl3Ca_SS_ALP"
    t = "Rotor RS40175 with LAW3Cl3Ca for MAS QCPMG\nUsing short pulses"
    _expno(root, L2, 2, "35Cl", ns=48000, pulprog="qcpmg", bf1=107.85, title=t)
    _expno(root, L3, 2, "35Cl", ns=48000, pulprog="qcpmg", bf1=107.85, title=t)
    _expno(root, E, 42, "27Al", ns=1536, pulprog="mp3qdfsz", twod=True,
           title="27Al zg power check\n\nRotor SR31648 - packed with cryolite\nMASR 35.714 kHz")
    _expno(root, E, 33, "27Al", ns=512, title="27Al zg\n\nRotor SR31649")
    _expno(root, "01192026_SR31648_Cryolite_SS_ALP", 16, "23Na", ns=8, bf1=BF1_NA,
           title="1H spectrum for later referencing")
    _expno(root, "01192026_SR31648_Cryolite_SS_ALP", 17, "23Na", ns=8, bf1=BF1_NA,
           title="23Na 1D for power optimization")
    _expno(root, "01192026_SR31648_Cryolite_SS_ALP", 1, "1H", ns=1, bf1=BF1_H,
           title="1H for referencing")
    _expno(root, A, 3114, "31P", ns=140, bf1=BF1_P, title="31P\n\nRotor RS2427418")
    _expno(root, "06052026_An_SS_ALP", 1, "35Cl", ns=8, bf1=107.85, title="An 35Cl")
    _expno(root, "06062026_RS40175_LAW3Cl2Ca_SS_ALP", 2, "35Cl", ns=48000, pulprog="qcpmg",
           bf1=107.85, title="35Cl QCPMG")
    inv = I.build(root, sr_audit=False)
    f = {(r.folder, r.expno): r.flags for r in inv.rows}
    assert any("LAW3Cl3Ca" in s and "another sample" in s for s in f[(L2, 2)])
    assert not any("another sample" in s for s in f[(L3, 2)])             # its own key
    assert any("SR31648" in s and "SR31649" in s for s in f[(E, 42)])       # rotor
    assert not any("another sample" in s for s in f[(E, 42)])             # line 3 only
    assert any("zg" in s and "mp3qdfsz" in s for s in f[(E, 42)])          # kind
    assert _row(inv, E, 42).role == "2D"
    assert f[(E, 33)] == []
    assert f[(A, 3114)] == []                                             # no folder rotor
    cry = "01192026_SR31648_Cryolite_SS_ALP"
    assert any("starts with 1H" in s and "23Na" in s for s in f[(cry, 16)])
    assert not any("starts with" in s for s in f[(cry, 17)])
    assert not any("another sample" in s for s in f[(cry, 16)])          # 'An' is 2 letters
    # two folders share the key LAW3Cl2Ca: both picks carry the duplicate flag
    dup = [r for r in inv.rows if r.sample == "LAW3Cl2Ca" and r.pick]
    assert len(dup) == 2 and all(any("shared with folder" in s for s in r.flags) for r in dup)
    assert inv.labels()[(L2, "LAW3Cl2Ca")] == "LAW3Cl2Ca (06052026)"
    assert inv.n_flags() == 4          # L2/2, E/42, Cryolite/16, 06062026/2


def test_maglab_set_samples_from_titles_and_block_picks(tmp_path):
    inv = I.build(_maglab(tmp_path), sr_audit=False)
    assert inv.root.name == "35Cl_2025-12" and inv.n_folders() == 1
    labels = set(inv.labels().values())
    assert labels == {"LAW3CL0CA", "LAW3CL4CA", "CS-3Cl"}
    cs = _pick(inv, "35Cl_2025-12", "CS-3Cl", "35Cl")
    assert cs.expno == 15 and _row(inv, "35Cl_2025-12", 10).role == "short"
    assert _pick(inv, "35Cl_2025-12", "LAW3CL4CA", "35Cl").expno == 3
    assert _row(inv, "35Cl_2025-12", 2).role == "short"
    # popt.array in EVERY EXPNO of the set is boilerplate, not a demotion
    assert not any("popt" in s for r in inv.rows for s in r.reasons)
    assert all(r.role in ("production", "short") for r in inv.rows)


def test_sr_join_uses_the_audit_verdicts(tmp_path):
    root = tmp_path / "2026-05"
    _expno(root, C, 1, "1H", ns=1, bf1=BF1_H, sf=SF_H, title="1H spectrum for later referencing")
    al_ok = R.expected_sf_MHz(SF_H, "27Al")
    _expno(root, B, 2701, "27Al", ns=512, sf=al_ok, title="27Al")
    _expno(root, B, 2702, "27Al", ns=512, sf=BF1_AL, title="27Al forgot xiref")
    _expno(root, B, 3102, "31P", ns=140, bf1=BF1_P, has_1r=False, title="31P")
    _expno(root, "06052026_LAW", 5, "35Cl", ns=8, bf1=107.85, sf=107.8503, title="35Cl")
    inv = I.build(root)
    v = {(r.folder, r.expno): r.sr_verdict for r in inv.rows}
    assert v[(B, 2701)] == "ok" and v[(B, 2702)] == "unreferenced"
    assert v[(B, 3102)] == "unprocessed" and v[("06052026_LAW", 5)] == "no reference"
    assert v[(C, 1)] == "1H reference"
    assert "xiref" in _row(inv, B, 2702).sr_note
    off = I.build(root, sr_audit=False)
    assert all(r.sr_verdict == "" for r in off.rows)


def test_outputs_csv_picks_grid_text(tmp_path):
    inv = I.build(_month(tmp_path), sr_audit=False)
    out = I.to_csv(inv, tmp_path / "inv.csv")
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == ("path,session,folder,sample,label,expno,nucleus,ndim,kind,pulse_program,"
                        "ns,d1_s,has_1r,has_2rr,n_procs,has_popt,date,sr_verdict,sr_note,role,"
                        "reasons,flags,pick,title")
    assert len(lines) == 1 + len(inv.rows)
    assert any(",production," in ln and ",3114," in ln for ln in lines)
    txt = I.picks_text(inv, "31P")
    body = [ln for ln in txt.splitlines() if not ln.startswith("#")]
    assert txt.startswith("# sample\tnucleus\tEXPNO\tpath")
    assert [ln.split("\t")[:3] for ln in body] == [["P5-Bi8-12 (04272026)", "31P", "3114"],
                                                   ["P5-Bi8-12 (05082026)", "31P", "3102"]]
    assert all(ln.split("\t")[3].endswith("1r") for ln in body)
    grid = I.grid_text(inv)
    rows = grid.splitlines()
    assert rows[0].split() == ["sample", "27Al", "31P", "1H"]
    a_line = next(ln for ln in rows if ln.startswith("P5-Bi8-12 (04272026)"))
    assert a_line.split()[-3:] == ["2704*", "3114*", "-"]      # 3114: duplicate-name flag
    b_line = next(ln for ln in rows if ln.startswith("P1-Bi1-12"))
    assert b_line.split()[-3:] == ["2704", "-", "-"]           # clean: demotions are lower EXPNOs
    assert "1 ref" in next(ln for ln in rows if ln.startswith("P5-Bi1-12"))
    short = I.text_report(inv)
    assert "2799" in short and "failed" in short and f"{B}/2704" not in short
    assert f"{B}/2704" in I.text_report(inv, all_rows=True)


def test_cli_inventory(tmp_path, capsys):
    from larmor import cli

    root = _month(tmp_path)
    out = tmp_path / "inv.csv"
    rc = cli.main(["inventory", str(root / A / "2702"), "--csv", str(out), "--no-sr"])
    assert rc == 0
    text = capsys.readouterr().out
    assert "P5-Bi8-12" in text and "3114" in text and "summary:" in text
    assert "5 sample folder" in text and "6 production pick" in text
    assert out.exists() and (tmp_path / "inv_picks.txt").exists()
    rc = cli.main(["inventory", str(root), "--picks", "--nucleus", "31P", "--no-sr"])
    assert rc == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2 and all(ln.endswith("1r") for ln in lines)
    assert lines[0].endswith(str(Path(A) / "3114" / "pdata" / "1" / "1r"))
    empty = tmp_path / "empty"
    empty.mkdir()
    assert cli.main(["inventory", str(empty), "--no-sr"]) == 1
    assert "give the month folder" in capsys.readouterr().err
    rc = cli.main(["inventory", str(root), "--ns-frac", "0.01", "--all", "--no-sr"])
    assert rc == 0
    text = capsys.readouterr().out
    assert f"{A}/2799" in text and "production" in text


# ---------------------------------------------------------------- real sessions
def test_real_2026_05_picks_failed_title_and_unprocessed_flag():
    require(BRUKER_1R)
    inv = I.build(BRUKER_1R.parents[4], sr_audit=False)
    assert inv.root.name == "2026-05"
    assert not any(r.folder.startswith("~") for r in inv.rows)
    assert inv.n_folders() == 13 and len(inv.rows) == 135
    P5, P1, P5b = "04272026_P5-Bi8-12_SS_ALP", "04272026_P1-Bi1-12_SS_ALP", "05082026_P5-Bi8-12_SS_ALP"
    assert _pick(inv, P5, "P5-Bi8-12", "27Al").expno == 2704
    assert _row(inv, P5, 2799).role == "short"
    assert _pick(inv, P5, "P5-Bi8-12", "31P").expno == 3114
    assert _pick(inv, P1, "P1-Bi1-12", "27Al").expno == 2704
    assert _row(inv, P1, 2702).role == "failed" and _row(inv, P1, 2703).role == "short"
    assert _pick(inv, P1, "P1-Bi1-12", "31P").expno == 3107
    p = _pick(inv, P5b, "P5-Bi8-12", "27Al")
    assert p.expno == 2799
    assert any("2702" in s and "2703" in s and "no pdata/1/1r" in s for s in p.flags)
    labels = inv.labels()
    assert labels[(P5, "P5-Bi8-12")] == "P5-Bi8-12 (04272026)"
    assert labels[(P5b, "P5-Bi8-12")] == "P5-Bi8-12 (05082026)"
    assert inv.nuclei() == ["27Al", "31P", "1H"]
    assert inv.cell(("04272026_RS40339_P5-Bi1-12", "P5-Bi1-12"), "1H").role == "reference"
    assert len(inv.picks()) == 24 and inv.n_flags() == 6


def test_real_2026_01_five_names_rotor_flag_and_power_check_block():
    require(LAW_CA_11B[0])
    inv = I.build(LAW_CA_11B[0].parents[1], sr_audit=False)
    labels = inv.labels()
    names = set(labels.values())
    assert {"Base0Ca", "Base1Ca", "Base2Ca", "Base3Ca", "Base4Ca"} <= names
    five = [r for r in inv.picks("11B") if labels[r.sample_id].startswith("Base")]
    assert [labels[r.sample_id] for r in five] == ["Base0Ca", "Base1Ca", "Base2Ca", "Base3Ca", "Base4Ca"]
    assert [r.expno for r in five] == [24] * 5
    assert [r.openable for r in five] == [str(e / "pdata" / "1" / "1r") for e in LAW_CA_11B]
    b0 = "01192026_SR31649_Base0Ca_SS_ALP"
    assert _row(inv, b0, 21).role == "setup" and _row(inv, b0, 22).role == "short"
    assert _row(inv, b0, 23).role == "arrayed"
    al = _pick(inv, b0, "Base0Ca", "27Al")
    assert al.expno == 33 and any("power check" in s for s in al.reasons)
    assert any("SR31648" in s and "SR31649" in s for s in al.flags)
    b1 = "01202026_SR31649_Base1Ca_SS_ALP"
    assert any("SR31648" in s for s in _row(inv, b1, 24).flags)
    r42 = _row(inv, b0, 42)
    assert r42.role == "2D" and any("mp3qdfsz" in s for s in r42.flags)
    assert len(inv.rows) == 133 and inv.n_folders() == 16


def test_real_2026_06_35cl_sibling_title_flag():
    root = require(DATA_ROOT / "Desktop/WSU_work/NMR/NMRFAM/DATA/2026-06_35Cl")
    inv = I.build(root, sr_audit=False)
    assert inv.n_folders() == 5 and len(inv.rows) == 13
    l2 = "06052026_RS40175_LAW3Cl2Ca_SS_ALP"
    for k in (2, 3):
        assert any("LAW3Cl3Ca" in s and "another sample" in s for s in _row(inv, l2, k).flags)
    assert _pick(inv, l2, "LAW3Cl2Ca", "35Cl").expno == 3 and _row(inv, l2, 2).role == "short"
    l3 = "06052026_RS40175_LAW3Cl3Ca_SS_ALP"
    assert not any("another sample" in s for s in _row(inv, l3, 2).flags)
    assert inv.cell((l3, "LAW3Cl3Ca"), "1H").role == "reference"
    assert inv.n_flags() == 2


def test_real_maglab_35cl_samples_from_title():
    require(MAGLAB_35CL / "1" / "acqus")
    inv = I.build(MAGLAB_35CL, sr_audit=False)
    assert inv.root == MAGLAB_35CL and inv.n_folders() == 1
    names = set(inv.labels().values())
    assert {"LAW3CL0CA", "LAW3CL4CA", "CS-3Cl"} <= names and len(names) == 14
    f = "35Cl_2025-12"
    assert _pick(inv, f, "LAW3CL4CA", "35Cl").expno == 3 and _row(inv, f, 2).role == "short"
    assert _pick(inv, f, "CS-3Cl", "35Cl").expno == 15 and _row(inv, f, 10).role == "short"
    assert all(r.role in ("production", "short", "candidate") for r in inv.rows)
    assert len(inv.picks()) == 14 and len(inv.rows) == 16
