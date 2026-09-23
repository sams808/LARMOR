import pytest

from larmor.io import bruker

from conftest import EXPNO_1903, require


def test_read_expno_1903_readonly():
    path = require(EXPNO_1903)
    before = bruker.snapshot(path)
    exp = bruker.read_expno(path, verify=True)  # raises on any modification
    bruker.verify_untouched(path, before)

    assert exp.nucleus == "19F"
    assert exp.sfo1_MHz == pytest.approx(564.27052307)
    assert exp.pulse_program == "hahnecho.nmrfam"
    assert exp.fid is not None and exp.fid.size > 0
    if (EXPNO_1903 / "pdata" / "1" / "1r").exists():
        assert exp.processed is not None
        assert exp.processed_ppm is not None
        assert exp.processed_ppm.size == exp.processed.size
    else:  # the user reorganized this dataset; raw-fid path must still work
        assert exp.processed is None


def test_metadata_conflict_detected():
    """acqus MASR (4200 Hz) disagrees with the title (35.714 kHz) in this EXPNO."""
    exp = bruker.read_expno(require(EXPNO_1903))
    assert any("MAS rate" in c for c in exp.conflicts)


def test_is_expno():
    assert bruker.is_expno(require(EXPNO_1903))
    assert not bruker.is_expno(EXPNO_1903 / "pdata")


class TestResolveMas:
    """_resolve_mas: static (rate 0) is a first-class outcome. A recorded
    MASR=0 used to fall through to the 35714 Hz fallback, so a static 81Br
    WCPMG glass dataset was silently simulated as fast MAS."""

    def test_positive_masr_is_trusted(self):
        assert bruker._resolve_mas({"MASR": 20000}, "") == (20000.0, False)

    def test_title_rate_alone_is_trusted(self):
        rate, unc = bruker._resolve_mas({}, "sample X MAS 12.5 kHz")
        assert rate == pytest.approx(12500.0) and unc is False

    def test_disagreeing_rates_flagged(self):
        rate, unc = bruker._resolve_mas({"MASR": 4200}, "MAS 35.714 kHz")
        assert rate == pytest.approx(35714.0) and unc is True

    def test_masr_zero_alone_is_static_but_flagged(self):
        assert bruker._resolve_mas({"MASR": 0}, "") == (0.0, True)

    def test_masr_zero_plus_static_title_is_certain(self):
        assert bruker._resolve_mas({"MASR": 0}, "81Br glass, static") == (0.0, False)

    def test_masr_zero_plus_wcpmg_pulprog_is_certain(self):
        acq = {"MASR": 0, "PULPROG": "<wcpmg.jp>"}
        assert bruker._resolve_mas(acq, "") == (0.0, False)

    def test_static_title_alone_is_certain(self):
        # operator-typed "static" is trusted exactly like operator-typed kHz
        assert bruker._resolve_mas({}, "PBi glass static 81Br") == (0.0, False)
        assert bruker._resolve_mas({}, "verre statique") == (0.0, False)

    def test_static_pulprog_alone_is_flagged(self):
        assert bruker._resolve_mas({"PULPROG": "wurst_echo"}, "") == (0.0, True)

    def test_qcpmg_alone_is_NOT_static_evidence(self):
        # MAS-QCPMG is a real technique; qcpmg must keep the old fallback
        rate, unc = bruker._resolve_mas({"PULPROG": "qcpmg"}, "")
        assert rate == bruker.MAS_FALLBACK_HZ and unc is True

    def test_positive_masr_under_static_pulprog_is_flagged(self):
        # the real MagLab 81Br WCPMG set: a static probe, MASR left at the
        # previous session's 14 kHz -- the rate still wins (a hint does not
        # prove static) but it must not load silently as MAS
        rate, unc = bruker._resolve_mas(
            {"MASR": 14000, "PULPROG": "<WCPMG.jk>"}, "81 Br\n0 Ca Br")
        assert rate == pytest.approx(14000.0) and unc is True

    def test_masr_zero_against_title_rate_is_flagged(self):
        # controller says not spinning, operator typed 20 kHz: a disagreement
        rate, unc = bruker._resolve_mas({"MASR": 0}, "MAS 20 kHz")
        assert rate == pytest.approx(20000.0) and unc is True

    def test_static_word_does_not_match_inside_words(self):
        # "electrostatic" must not read as a static experiment
        rate, unc = bruker._resolve_mas({}, "electrostatic sample holder")
        assert rate == bruker.MAS_FALLBACK_HZ and unc is True

    def test_nothing_found_keeps_the_flagged_fallback(self):
        assert bruker._resolve_mas({}, "") == (bruker.MAS_FALLBACK_HZ, True)


def _jcamp(path, **params):
    lines = ["##TITLE= Parameter file, TopSpin 4.1", "##JCAMPDX= 5.0",
             "##DATATYPE= Parameter Values", "##ORIGIN= Bruker BioSpin GmbH"]
    for k, v in params.items():
        if isinstance(v, (list, tuple)):
            lines.append(f"##${k}= (0..{len(v) - 1})")
            lines.append(" ".join(str(x) for x in v))
        else:
            lines.append(f"##${k}= {v}")
    lines.append("##END=")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_meta_1d_reads_ns_pulse_power_and_aq(tmp_path):
    """The additive acquisition keys the quantitativity chips read (Base1Ca/24
    values); every pre-existing key is unchanged."""
    # nmrglue delivers the <...> strings stripped, PROBHD sometimes not
    acqus = {"NUC1": "11B", "SFO1": 192.43, "SW_h": 100000, "TD": 7988,
             "NS": 256, "D": [0.2, 14.0], "P": [0, 0.425], "PLW": [0, 100],
             "PROBHD": "<16_Solenoid (PMAS16)>", "PULPROG": "zg"}
    m = bruker._meta_1d(acqus, "title", tmp_path)
    assert m["ns"] == 256 and m["d1_s"] == 14.0
    assert m["aq_s"] == pytest.approx(0.03994)
    assert m["p_us"][1] == 0.425 and m["plw_w"][1] == 100.0
    assert m["probhd"] == "16_Solenoid (PMAS16)"
    assert m["d"] == [0.2, 14.0] and m["td"] == 7988 and m["sw_Hz"] == 100000.0
    assert m["nucleus"] == "11B" and m["pulse_program"] == "zg" and m["title"] == "title"
    assert m["expno"] == str(tmp_path) and "sf_MHz" not in m     # no procs here
    bare = bruker._meta_1d({"NUC1": "<27Al>", "SFO1": 130.3}, "", tmp_path)
    assert bare["d1_s"] is None and bare["aq_s"] is None and bare["ns"] == 0
    assert bare["p_us"] == [] and bare["plw_w"] == [] and bare["probhd"] == ""
    # the fid-free reader gives the same keys from a folder
    e = tmp_path / "24"
    _jcamp(e / "acqus", NUC1="<11B>", SFO1=192.43, SW_h=100000, TD=7988, NS=256,
           D=[0.2, 14.0, 0, 0], P=[0, 0.425, 0, 0], PLW=[0, 100, 0, 0],
           PROBHD="<16_Solenoid (PMAS16)>", PULPROG="<zg>")
    (e / "pdata" / "1").mkdir(parents=True)
    (e / "pdata" / "1" / "title").write_text("11B with short tip angle\n")
    r = bruker.read_acqus_meta(e)
    assert r["ns"] == 256 and r["d1_s"] == 14.0 and r["aq_s"] == pytest.approx(0.03994)
    assert r["p_us"][1] == 0.425 and r["plw_w"][1] == 100.0
    assert r["probhd"] == "16_Solenoid (PMAS16)" and r["nucleus"] == "11B"
    assert r["title"].startswith("11B with short tip angle")
    assert set(m) == set(r)
    assert bruker.read_acqus_meta(tmp_path / "nothing") == {}
