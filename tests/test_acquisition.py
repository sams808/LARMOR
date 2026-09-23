"""larmor.acquisition: the flat acquisition block of an EXPNO, the title
rules, the Experimental sentences that never state what the files do not
support, Table S1 with its varying columns, the per-nucleus paragraph and
the referencing evidence -- on synthetic JCAMP EXPNOs and on the real 3102,
2702 and MagLab 81Br anchors (read-only, snapshot/verify_untouched)."""
import json
from pathlib import Path

import numpy as np
import pytest

from conftest import BRUKER_1R, EXPNO_3102, LAW_CA_11B, MAGLAB_81BR, require
from larmor import acquisition as A


# ----------------------------------------------------------- synthetic EXPNOs
def jcamp(path: Path, **params):
    """A TopSpin-style JCAMP parameter file; a list value becomes a
    ``(0..n-1)`` array block the way acqus writes D / P / PLW."""
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


def make_expno(root: Path, sample: str, expno, *, nuc="27Al", bf1=156.281744,
               sf=156.2816, ns=512, d1=5.0, p1=0.375, plw1=100.0, lb=100.0,
               wdw=1, tdeff=2048, si=1024, ph_mod=1, absg=5, masr=22000, date=1777400000,
               title="P1(90)=3.125; 11deg tip = 0.375\nRotor RS43009\nMASR 22kHz",
               probe="Z1 3.2mm", pulprog="zg", uxnmr: dict | None = None,
               booking_khz=None, procnos=(1,), scale=1.0, with_1r=True,
               td=4096, sw=100000, with_fid=False) -> Path:
    """One synthetic EXPNO ``root/sample/expno`` with acqus, one pdata per
    procno (procs + title + a Gaussian 1r scaled by ``scale * procno``),
    optional uxnmr.info / experiment_addenda.xml sidecars and, with
    ``with_fid``, a decaying complex fid nmrglue can read."""
    e = root / sample / str(expno)
    jcamp(e / "acqus", NUC1=f"<{nuc}>", BF1=bf1, SFO1=bf1, O1=0, SW_h=sw, TD=td,
          NS=ns, DS=4, RG=100, DATE=date, PULPROG=f"<{pulprog}>", MASR=masr,
          PROBHD=f"<{probe}>", INSTRUM="<spect>", BYTORDA=0, DTYPA=0, GRPDLY=0,
          D=[0.2, d1, 0], P=[0, p1, 0], PLW=[0, plw1, 0])
    if with_fid:
        t = np.arange(td // 2) / sw
        sig = 1e5 * np.exp(-t * 200.0) * np.exp(2j * np.pi * 1000.0 * t)
        raw = np.empty(td, dtype="<i4")
        raw[0::2], raw[1::2] = sig.real.astype("<i4"), sig.imag.astype("<i4")
        (e / "fid").write_bytes(raw.tobytes())
    for procno in procnos:
        pd = e / "pdata" / str(procno)
        if sf is not None:
            jcamp(pd / "procs", SF=sf, SI=si, OFFSET=200.0, SW_p=sw, BYTORDP=0, DTYPP=0,
                  NC_proc=0, WDW=wdw, LB=lb, GB=0, SSB=0, TDeff=tdeff, PH_mod=ph_mod,
                  ABSG=absg, FT_mod=6, PHC0=0, PHC1=0)
        pd.mkdir(parents=True, exist_ok=True)
        (pd / "title").write_text(title, encoding="utf-8")
        if with_1r and sf is not None:
            x = np.arange(si)
            y = (scale * procno * 1000 * np.exp(-0.5 * ((x - si / 2) / 20) ** 2)).astype("<i4")
            (pd / "1r").write_bytes(y.tobytes())
    if uxnmr:
        (e / "uxnmr.info").write_text(
            "CONFIGURATION INFORMATION\n=========================\n\n"
            + "".join(f"{k:<13}: {v}\n" for k, v in uxnmr.items()), encoding="latin-1")
    if booking_khz is not None:
        (e / "experiment_addenda.xml").write_text(
            "<experiment_addenda><magic_angle_spinning>true</magic_angle_spinning>"
            f"<magic_angle_spinning_rate>{booking_khz}</magic_angle_spinning_rate>"
            "</experiment_addenda>", encoding="utf-8")
    return e


UX600 = {"Release": "TopSpin 3.6.2", "System": "Avance III HD 600 NMR spectrometer",
         "1H-frequency": "599.77 MHz"}


def series_11b(root: Path) -> list[Path]:
    """Five 11B EXPNO 24 like the Tutorial-4 set: identical except D1
    12.5/14/18/29/36 s, LB 0/0/100/100/100 Hz and the rotor (SR31649 on the
    first, SR31648 on the rest)."""
    out = []
    for k, (name, d1, lb, rotor) in enumerate((
            ("01192026_SR31649_Base0Ca_SS_ALP", 12.5, 0, "SR31649"),
            ("01202026_SR31649_Base1Ca_SS_ALP", 14.0, 0, "SR31648"),
            ("01202026_SR31649_Base2Ca_SS_ALP", 18.0, 100, "SR31648"),
            ("01202026_SR31648_Base3Ca_SS_ALP", 29.0, 100, "SR31648"),
            ("01202026_SR31649_Base4Ca_SS_ALP", 36.0, 100, "SR31648"))):
        out.append(make_expno(
            root, name, 24, nuc="11B", bf1=192.430693, sf=192.431528, ns=256, d1=d1,
            p1=0.425, plw1=100.0, lb=lb, tdeff=1024, si=32768, masr=4200,
            date=1768921303 + 3600 * k, uxnmr=UX600, booking_khz=35.714,
            title=f"11B with short tip angle\n\nRotor {rotor}\nMASR 35.714 kHz\n"
                  "VT ambient - no flow\nALP"))
    return out


def _untouched(path: Path):
    from larmor.io import bruker

    class _Ctx:
        def __enter__(self):
            self.before = bruker.snapshot(path)
            return self

        def __exit__(self, *a):
            bruker.verify_untouched(path, self.before)
    return _Ctx()


# ------------------------------------------------------------- real anchors
def test_read_block_real_3102_matches_topspin():
    e = require(EXPNO_3102)
    with _untouched(e):
        b = A.read_block(e)
    json.dumps(b)                                    # every value a JSON scalar
    assert set(b) == set(A.BLOCK_KEYS)
    assert b["nucleus"] == "31P" and b["ns"] == 14 and b["ds"] == 0
    assert b["d1_s"] == 300.0 and b["p1_us"] == 1.25 and b["plw1_W"] == 207.0
    assert b["probe"] == "SPRB600511_7297 (MAS)" and b["pulprog"] == "zg"
    assert b["p90_us"] == 3.75 and b["flip_deg"] == 30.0 and b["flip_source"] == "title"
    assert b["rotor"] == "RS2427418" and b["vt_note"] and b["decoupling_note"] == ""
    assert b["spectrometer"] == "Avance III HD 600" and b["topspin"] == "3.6.2"
    assert b["magnet_1h_MHz"] == 599.77 and b["b0_T"] == pytest.approx(14.09, abs=0.01)
    assert (b["mas_acqus_Hz"], b["mas_title_Hz"], b["mas_booking_Hz"]) == (4200, 20000, 22000)
    # N9's rule: three disagreeing sources keep the highest, flagged
    assert b["spin_rate_Hz"] == 22000.0 and b["mas_uncertain"] is True
    assert b["mas_source"] == "highest"
    assert b["sf_MHz"] == pytest.approx(242.79194514) and b["sr_hz"] == pytest.approx(-210.86, abs=0.1)
    assert b["td"] == 9590 and b["sw_Hz"] == 100000 and b["tdeff"] == 384 and b["si"] == 32768
    assert b["lb_hz"] == 0 and b["wdw"] == "EM" and b["ph_mod"] == "pk" and b["absg"] == 5
    assert b["date_utc"].startswith("2026-05-02T04:09:33")
    assert b["data_file"] == "pdata/1/1r" and b["procno"] == 1
    assert b["sample"] == "P5-Bi8-12" and b["sample_folder"] == "04272026_P5-Bi8-12_SS_ALP"
    # the same block from the 1r file and from the pdata folder
    assert A.read_block(e / "pdata" / "1" / "1r") == b == A.read_block(e / "pdata" / "1")


def test_static_maglab_81br_block():
    e = require(MAGLAB_81BR / "30")
    with _untouched(e):
        b = A.read_block(e)
    assert b["nucleus"] == "81Br" and b["ns"] == 512000
    assert b["p1_us"] == 2.8 and b["plw1_W"] == 63.0 and "NHMFL 3.2mm" in b["probe"]
    assert b["pulprog"] == "WCPMG.jk" and b["spectrometer"] == "Avance Neo 800"
    assert b["magnet_1h_MHz"] == 799.74 and b["b0_T"] == pytest.approx(18.78, abs=0.01)
    assert b["mas_acqus_Hz"] == 14000 and b["mas_uncertain"] is True
    assert b["mas_booking_Hz"] is None and b["wdw"] == "none" and b["ph_mod"] == "mc"
    assert b["flip_source"] == "" and b["flip_deg"] is None
    # the recipe confirmed static: the paragraph says so and never prints the stale MASR
    txt = A.paragraph([b], spin_rates={b["expno_path"]: (0.0, False)})
    assert "static" in txt and "14 kHz" not in txt and "14000" not in txt
    # unconfirmed: a bracket naming the stale source and the static program
    unc = A.paragraph([b])
    assert "[confirm MAS rate: acqus 14 kHz; static pulse program WCPMG.jk]" in unc


def test_real_2702_title_flip_and_series_variation():
    e = require(BRUKER_1R)
    b = A.read_block(e)
    assert b["flip_deg"] == 11.0 and b["flip_source"] == "title" and b["p90_us"] == 3.125
    assert b["ns"] == 56 and b["p1_us"] == 0.375 and b["data_file"] == "pdata/1/1r"
    paths = [require(p) for p in LAW_CA_11B]
    blocks = [A.read_block(p) for p in paths]
    t = A.table(blocks)
    assert t.varying == {"d1_s", "lb_hz", "rotor"}
    assert t.constants["ns"] == 256 and t.constants["p1_us"] == 0.425
    assert t.constants["plw1_W"] == 100.0 and t.constants["spin_rate_Hz"] == 35714.0
    para = A.paragraph(blocks)
    assert "recycle delays of 12.5–36 s (Table S1)" in para
    assert "LB = 0–100 Hz (Table S1)" in para and "SR31649 / SR31648" in para
    assert "short tip angle" in para and "°" not in para.split("[")[0].split("pulse")[1][:40]


# ------------------------------------------------------------- synthetic
def test_read_block_synthetic_expno_without_sidecars(tmp_path):
    from larmor.referencing import xi_ratio

    e = make_expno(tmp_path / "2026-05", "04272026_RS43009_S1_SS_ALP", 2701)
    b = A.read_block(e)
    assert b["spectrometer"] == "" and b["topspin"] == "4.1"
    assert b["magnet_1h_MHz"] == pytest.approx(156.281744 / xi_ratio("27Al"), abs=0.01)
    assert b["magnet_1h_MHz"] == pytest.approx(599.77, abs=0.01)
    assert b["mas_booking_Hz"] is None
    assert b["flip_deg"] == 11.0 and b["flip_source"] == "title"
    assert b["spin_rate_Hz"] == 22000.0 and b["mas_uncertain"] is False   # acqus + title agree
    assert b["lb_hz"] == 100.0 and b["tdeff"] == 2048 and b["ph_mod"] == "pk"
    assert b["rotor"] == "RS43009" and b["sample"] == "S1"
    assert b["sr_hz"] == pytest.approx((156.2816 - 156.281744) * 1e6)
    # no acqus -> {}; a truncated procs never raises
    assert A.read_block(tmp_path / "nothing") == {}
    (e / "pdata" / "1" / "procs").write_text("##TITLE= broken\n##$SF= ", encoding="utf-8")
    b2 = A.read_block(e)
    assert b2 and b2["nucleus"] == "27Al"
    # uxnmr.info wins over BF1/Xi and gives the spectrometer name
    e2 = make_expno(tmp_path / "2026-05", "04272026_RS43009_S2_SS_ALP", 2702, uxnmr=UX600,
                    booking_khz=22)
    b3 = A.read_block(e2)
    assert b3["spectrometer"] == "Avance III HD 600" and b3["magnet_1h_MHz"] == 599.77
    assert b3["topspin"] == "4.1" and b3["mas_booking_Hz"] == 22000.0
    assert b3["mas_source"] == "all"


def test_flip_angle_and_title_rules():
    assert A.flip_angle(1.25, 3.75, "P1(90)=3.750; 30 deg tip (pi/6) = 1.25") == (30.0, "title")
    v, src = A.flip_angle(0.375, 3.125, "P1(90)=3.125")
    assert v == pytest.approx(10.8) and src == "computed"
    assert A.flip_angle(0.425, None, "11B with short tip angle") == (None, "")
    assert A.parse_title("11B with short tip angle")["tip_words"] == "short tip angle"
    t = "27Al zg - spectrum with 11 degree tip angle and 5*T1 with 19F decoupling"
    assert A.flip_angle(1.0, None, t) == (11.0, "title")
    assert A.parse_title(t)["decoupling_note"] == "with 19F decoupling"
    assert A.parse_title("27Al zg without 19F decoupling")["decoupling_note"] == "without 19F decoupling"
    assert A.parse_title("27Al zg")["decoupling_note"] == ""
    assert A.flip_angle(1.0, None, "23Na zg spectrum using 5xT1 and 30 degree pulse") == (30.0, "title")
    assert A.flip_angle(1.0, None, "180 deg pulse") == (None, "")
    assert A.parse_title("VT ambient - no flow")["vt_note"] == "ambient temperature, no VT gas flow"
    assert A.parse_title("VT at ambient, no flow")["vt_note"] == "ambient temperature, no VT gas flow"
    assert A.parse_title("31P zg")["vt_note"] == ""
    assert A.parse_title("Rotor SR31649 with Base0Ca")["rotor"] == "SR31649"
    assert A.parse_title("MASR 20 kHz")["mas_title_Hz"] == 20000.0


def _block_3102_like(**over) -> dict:
    b = {k: None for k in A.BLOCK_KEYS}
    b.update({
        "expno_path": "X:/DATA/2026-05/04272026_P5-Bi8-12_SS_ALP/3102", "expno": "3102",
        "procno": 1, "data_file": "pdata/1/1r", "sample": "P5-Bi8-12",
        "sample_folder": "04272026_P5-Bi8-12_SS_ALP", "title": "31P", "nucleus": "31P",
        "bf1_MHz": 242.792156, "sfo1_MHz": 242.792156, "sf_MHz": 242.79194514,
        "sr_hz": -210.86, "magnet_1h_MHz": 599.77, "b0_T": 14.0866,
        "spectrometer": "Avance III HD 600", "topspin": "3.6.2",
        "probe": "SPRB600511_7297 (MAS)", "pulprog": "zg", "ns": 14, "ds": 0, "rg": 194.07,
        "d1_s": 300.0, "aq_s": 0.04795, "p1_us": 1.25, "plw1_W": 207.0, "p90_us": 3.75,
        "flip_deg": 30.0, "flip_source": "title", "tip_words": "", "rotor": "RS2427418",
        "vt_note": "ambient temperature, no VT gas flow", "decoupling_note": "",
        "sw_Hz": 100000.0, "td": 9590, "date_utc": "2026-05-02T04:09:33Z",
        "mas_acqus_Hz": 4200.0, "mas_title_Hz": 20000.0, "mas_booking_Hz": 22000.0,
        "spin_rate_Hz": 22000.0, "mas_uncertain": True, "mas_source": "highest",
        "wdw": "EM", "lb_hz": 0.0, "gb": 0.0, "ssb": 0.0, "si": 32768, "tdeff": 384,
        "ph_mod": "pk", "phc0": 75.75, "phc1": -53.6, "absg": 5})
    b.update(over)
    return b


def test_sentences_3102_like_never_state_what_the_data_do_not_support():
    b = _block_3102_like()
    txt = " ".join(A.sentences(b, spin_rate_Hz=20000, mas_uncertain=False, sr_hz=-210.86,
                               larmor_ops=[{"op": "baseline", "order": 3}]))
    for must in ("242.79 MHz", "14.1 T", "599.77 MHz", "Avance III HD 600", "SPRB600511_7297",
                 "20.0 kHz", "1.25 µs", "30°", "3.75 µs", "207 W", "300 s", "14 transients",
                 "zg", "exponential", "TDeff", "384", "32768", "baseline (order 3)",
                 "ambient temperature", "SR = −210.86 Hz", "[state the reference standard]",
                 "TopSpin 3.6.2"):
        assert must in txt, must
    for never in ("decoupl", "4200", "adamantane", " K"):
        assert never not in txt, never
    # evidence given: the adamantane clause and no bracket at all
    ev = " ".join(A.sentences(b, spin_rate_Hz=20000, mas_uncertain=False, sr_hz=-210.86,
                              referencing_text="referenced indirectly (IUPAC Ξ) to the ¹H "
                                               "resonance of adamantane at 1.82 ppm"))
    assert "adamantane" in ev and "1.82 ppm" in ev and "[" not in ev
    # unconfirmed: the bracketed candidate list, never a silent guess
    unc = " ".join(A.sentences(b))
    assert "[confirm MAS rate: acqus 4.2 kHz, title 20 kHz, booking 22 kHz]" in unc
    assert "spinning rate was" not in unc
    # the recipe's own mas_rate block feeds the bracket and its confirmation
    conf = " ".join(A.sentences(b, spin_rate_Hz=20000, mas_uncertain=True,
                                mas_rate={"acqus_Hz": 4200.0, "title_Hz": 20000.0,
                                          "booking_Hz": 22000.0, "confirmed": "2026-09-23T10:00"}))
    assert "20.0 kHz" in conf and "[confirm" not in conf
    # static and confirmed: 'static', no MASR anywhere
    st = " ".join(A.sentences(b, spin_rate_Hz=0.0, mas_uncertain=False))
    assert "static" in st and "4200" not in st and "kHz" not in st.split(".")[1]
    # no flip evidence: the pulse length with the title's words, no degrees
    plain = " ".join(A.sentences(_block_3102_like(flip_deg=None, flip_source="",
                                                  p90_us=None, tip_words="short tip angle"),
                                 spin_rate_Hz=20000, mas_uncertain=False))
    assert "1.25 µs pulse (short tip angle) at 207 W" in plain and "°" not in plain.split("[")[0]
    comp = " ".join(A.sentences(_block_3102_like(flip_deg=10.8, flip_source="computed",
                                                  p90_us=3.125, p1_us=0.375),
                                 spin_rate_Hz=20000, mas_uncertain=False))
    assert "≈10.8° flip angle, from the title's 90° pulse of 3.125 µs" in comp
    # decoupling only from the title
    dec = " ".join(A.sentences(_block_3102_like(decoupling_note="with 19F decoupling"),
                               spin_rate_Hz=20000, mas_uncertain=False))
    assert "with 19F decoupling" in dec
    # no spectrometer name: the field alone, no 'Bruker'
    nos = " ".join(A.sentences(_block_3102_like(spectrometer=""), spin_rate_Hz=20000,
                               mas_uncertain=False))
    assert "Bruker" not in nos and "operating at B₀ = 14.1 T" in nos
    assert A.sentences({}) == []


def _five_11b_blocks() -> list[dict]:
    out = []
    for k, (name, d1, lb, rotor) in enumerate((
            ("Base0Ca", 12.5, 0.0, "SR31649"), ("Base1Ca", 14.0, 0.0, "SR31648"),
            ("Base2Ca", 18.0, 100.0, "SR31648"), ("Base3Ca", 29.0, 100.0, "SR31648"),
            ("Base4Ca", 36.0, 100.0, "SR31648"))):
        out.append(_block_3102_like(
            expno_path=f"X:/DATA/2026-01/{name}/24", expno="24", sample=name,
            sample_folder=f"0120202_SR31649_{name}_SS_ALP", nucleus="11B",
            bf1_MHz=192.430693, sfo1_MHz=192.430693, sf_MHz=192.431528, sr_hz=835.31,
            probe="16_Solenoid (PMAS16)", ns=256, d1_s=d1, p1_us=0.425, plw1_W=100.0,
            p90_us=None, flip_deg=None, flip_source="", tip_words="short tip angle",
            rotor=rotor, td=7988, date_utc=f"2026-01-2{k}T10:00:00Z",
            mas_acqus_Hz=4200.0, mas_title_Hz=35714.0, mas_booking_Hz=35714.0,
            spin_rate_Hz=35714.0, mas_uncertain=False, mas_source="title+booking",
            lb_hz=lb, tdeff=1024, absg=0))
    return out


def test_table_flags_varying_columns_ranges_and_latex_is_ascii():
    blocks = _five_11b_blocks()
    t = A.table(blocks)
    assert t.varying == {"d1_s", "lb_hz", "rotor"}
    assert not ({"sample", "expno", "date_utc", "expno_path"} & t.varying)
    assert t.constants["ns"] == 256 and t.constants["p1_us"] == 0.425
    assert t.varying_headers() == ["Rotor", "D1 (s)", "LB (Hz)"]
    csv_txt = A.table_csv(t)
    lines = csv_txt.strip().splitlines()
    assert len(lines) == 7 and "D1 (s)" in lines[0] and lines[-1].startswith("# varies")
    head = lines[0].split(",")
    last = lines[-1].split(",")
    assert last[head.index("D1 (s)")] == "yes" and last[head.index("NS")] == ""
    assert "12.5" in csv_txt and "SR31649" in csv_txt
    tex = A.table_latex(t, caption="Acquisition parameters.", compact=True)
    assert max(ord(c) for c in tex) < 0x80
    assert r"\begin{tabular}{llllll}" in tex          # sample, EXPNO, nucleus + 3 varying
    assert "Sample & EXPNO & Nucleus & Rotor & D1 (s) & LB (Hz)" in tex
    assert "NS = 256" in tex and "P1 (us) = 0.425" in tex and r"\textbf{12.5}" in tex
    assert "Columns in bold vary across the series" in tex and r"\toprule" in tex
    assert tex.count("\\\\\n") >= 7                      # every row ends with \\
    full = A.table_latex(t, compact=False)
    assert "Sample & EXPNO & Nucleus & $\\nu_0$ (MHz)" in full and "NS &" in full
    assert max(ord(c) for c in full) < 0x80
    md = A.table_markdown(t)
    assert "**D1 (s)**" in md and "**LB (Hz)**" in md and "| NS |" in md and "**12.5**" in md
    para = A.paragraph(blocks)
    assert "12.5–36 s" in para and "LB = 0–100 Hz" in para and "Table S1" in para
    assert "recycle delay of 12.5 s" not in para and "recycle delays of 12.5–36 s (Table S1)" in para
    assert "rotors SR31649 / SR31648 (Table S1)" in para
    assert "spectra were acquired" in para and "35.7 kHz" in para
    # a single block: no ranges, no Table S1
    one = A.paragraph(blocks[:1])
    assert "Table S1" not in one and "12.5 s" in one and "spectrum was acquired" in one
    assert A.paragraph([]) == "" and A.table([]).rows == []


def test_paragraph_groups_by_nucleus_and_uses_confirmed_rates():
    al = _block_3102_like(expno_path="X:/DATA/2026-05/S/2701", expno="2701", nucleus="27Al",
                          bf1_MHz=156.281744, sfo1_MHz=156.281744, sf_MHz=156.281609,
                          mas_acqus_Hz=4200.0, mas_title_Hz=20000.0, mas_booking_Hz=22000.0,
                          spin_rate_Hz=22000.0, mas_uncertain=True)
    p = _block_3102_like()
    txt = A.paragraph([al, p], spin_rates={al["expno_path"]: (22000.0, False)})
    paras = txt.split("\n\n")
    assert len(paras) == 2
    assert "27Al" in paras[0] and "156.28 MHz" in paras[0]
    assert "31P" in paras[1] and "242.79 MHz" in paras[1]
    assert paras[0].count("Avance III HD 600") == 1 and paras[1].count("Avance III HD 600") == 1
    # the confirmed rate wins over the block's own uncertain flag
    assert "22.0 kHz" in paras[0] and "[confirm" not in paras[0]
    # the 31P block stays unconfirmed: a bracket
    assert "[confirm MAS rate: acqus 4.2 kHz, title 20 kHz, booking 22 kHz]" in paras[1]
    # a rate that varies across a nucleus group is a range
    al2 = dict(al, expno_path="X:/DATA/2026-05/S/2703", expno="2703")
    var = A.paragraph([al, al2], spin_rates={al["expno_path"]: (20000.0, False),
                                             al2["expno_path"]: (22000.0, False)})
    assert "20.0–22.0 kHz (Table S1)" in var


def test_block_for_recipe_and_referencing_evidence(tmp_path, monkeypatch):
    from larmor import referencing as R

    block = _block_3102_like()
    monkeypatch.setattr(A, "read_block", lambda *a, **k: (_ for _ in ()).throw(RuntimeError))
    assert A.block_for_recipe({"acquisition": block}) == block          # no disk access
    assert A.block_for_recipe({"source_kind": "csv", "source_path": "x.csv"}) == {}
    assert A.block_for_recipe({}) == {}
    monkeypatch.undo()
    e = make_expno(tmp_path / "2026-05", "04272026_RS43009_S1_SS_ALP", 2701)
    b = A.block_for_recipe({"source_kind": "bruker", "source_path": str(e)})
    assert b["nucleus"] == "27Al" and b["expno_path"] == str(e)

    # a synthetic month with a referenced 1H spectrum: the 27Al whose SR is
    # the Xi prediction gets the adamantane clause, the unreferenced one none
    root = tmp_path / "month" / "2026-05"
    sf_h = 599.771479328723
    make_expno(root, "04272026_RS40339_P5", 1, nuc="1H", bf1=599.772, sf=sf_h,
               title="1H spectrum for later referencing", date=1777400000)
    ok = make_expno(root, "04272026_P1-Bi1-12", 2701, nuc="27Al", bf1=156.281744,
                    sf=R.expected_sf_MHz(sf_h, "27Al"), date=1777403600)
    bad = make_expno(root, "04272026_P1-Bi1-12", 2702, nuc="27Al", bf1=156.281744,
                     sf=156.281744, date=1777407200)
    A._audit_rows.cache_clear()
    ev = A.referencing_evidence([A.read_block(ok)])
    assert ev and "adamantane" in ev and "1.82" in ev and "Ξ" in ev
    assert A.referencing_evidence([A.read_block(bad)]) is None
    assert A.referencing_evidence([A.read_block(ok), A.read_block(bad)]) is None
    # a month without any referenced 1H
    root2 = tmp_path / "other" / "2026-06"
    lone = make_expno(root2, "S", 5, nuc="27Al", bf1=156.281744,
                      sf=R.expected_sf_MHz(sf_h, "27Al"))
    assert A.referencing_evidence([A.read_block(lone)]) is None
    # the recipe's own audit block wins, with its adamantane shift
    note = ("SR 0.00 → -135.31 Hz from the ¹H reference P5/1 (adamantane at 1.85 ppm, "
            "Ξ indirect)")
    ev2 = A.referencing_evidence([A.read_block(bad)],
                                 recipe={"provenance": {"referencing": {"note": note}}})
    assert ev2 and "1.85 ppm" in ev2
    assert A.referencing_evidence([]) is None
    # LARMOR steps in prose, unknown op as the Processing-steps dialog prints it
    d = A.describe_processing([{"op": "em", "lb_hz": 100}, {"op": "twopoint_bg"},
                               {"op": "sr", "sr_hz": 120}, {"op": "foo", "a": 1}])
    assert "exponential apodization (LB 100 Hz)" in d and "two-point linear baseline" in d
    assert "re-referencing (SR 120 Hz)" in d and "foo (a=1)" in d and " and foo" in d
    assert A.describe_processing([]) == "" and A.describe_processing(None) == ""
    assert set(A.OP_PHRASES) >= {"em", "gm", "sine", "traf", "tdeff", "zf", "ft", "phase",
                                 "autophase", "baseline", "iterbaseline", "flat_baseline",
                                 "twopoint_bg", "subtract_avg", "sr", "extract",
                                 "wurst_correct", "hilbert", "ift", "magnitude", "scale",
                                 "normalize", "fcor", "lp", "swap_echo", "echo_apodize"}


def test_summary_lines_and_folder_expnos(tmp_path):
    b = _block_3102_like()
    lines = A.summary_lines(b, {"larmor": "0.13.0", "git_commit": "e3c400dfd297",
                                "mrsimulator": "1.0.0", "lmfit": "1.3.4"}, 22000.0, True)
    txt = "\n".join(lines)
    for must in ("NS 14", "D1 300 s", "30°", "acqus 4200 Hz", "title 20 kHz", "booking 22 kHz",
                 "(confirm)", "EM LB 0 Hz", "TDeff 384", "SI 32768", "ABSG 5",
                 "Avance III HD 600 (599.77 MHz, 14.1 T)", "fitted with LARMOR 0.13.0",
                 "commit e3c400d", "mrsimulator 1.0.0", "lmfit 1.3.4"):
        assert must in txt, must
    assert A.summary_lines({}) == []
    assert "fitted with" not in "\n".join(A.summary_lines(b))
    # a sample folder gives its EXPNOs, a month folder the whole session
    month = tmp_path / "2026-05"
    e1 = make_expno(month, "04272026_S1_SS_ALP", 1)
    e2 = make_expno(month, "04272026_S1_SS_ALP", 12)
    e3 = make_expno(month, "04272026_S2_SS_ALP", 3)
    assert A.folder_expnos(month / "04272026_S1_SS_ALP") == [str(e1), str(e2)]
    assert set(A.folder_expnos(month)) == {str(e1), str(e2), str(e3)}
    assert A.folder_expnos(e1) == [str(e1)] and A.folder_expnos(tmp_path / "none") == []
