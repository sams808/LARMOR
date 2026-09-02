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
