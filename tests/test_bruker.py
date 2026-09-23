import pytest

from larmor.io import bruker

from conftest import BRUKER_1R, EXPNO_1903, LAW_CA_11B, require


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

    # -- the third source: the NMRFAM booking sidecar (larmor.masrate) --

    def test_booking_agreeing_with_title_outvotes_acqus(self):
        # the 405 title == booking EXPNOs: acqus 4200 is a leftover, no badge
        assert bruker._resolve_mas({"MASR": 4200}, "MASR 35.714 kHz",
                                   booking_Hz=35714.0) == (35714.0, False)

    def test_three_way_disagreement_is_flagged_highest_wins(self):
        # 2026-05: booked 22, titled 20, acqus 4200 -> 22000 flagged
        assert bruker._resolve_mas({"MASR": 4200}, "MASR 20 kHz",
                                   booking_Hz=22000.0) == (22000.0, True)
        # 2026-06_35Cl: booked 10, titled 20, acqus 5000 -> 20000 flagged
        assert bruker._resolve_mas({"MASR": 5000}, "MASR 20 kHz",
                                   booking_Hz=10000.0) == (20000.0, True)

    def test_booking_alone_and_masr_zero_against_booking(self):
        assert bruker._resolve_mas({}, "", booking_Hz=22000.0) == (22000.0, False)
        assert bruker._resolve_mas({"MASR": 0}, "", booking_Hz=20000.0) == (20000.0, True)
        # acqus and booking agree, no title: certain
        assert bruker._resolve_mas({"MASR": 20000}, "", booking_Hz=20000.0) == (20000.0, False)

    def test_title_typos_and_zero(self):
        # 'MASR 35741 kHz' (1900/1901): re-read as 35741 Hz, booking corroborates
        assert bruker._resolve_mas({"MASR": 4200}, "MASR 35741 kHz",
                                   booking_Hz=35714.0) == (35714.0, False)
        # 'MASR 35.714 Hz' (the 1103 series, no booking): 35714 flagged, not 4200
        assert bruker._resolve_mas({"MASR": 4200}, "MASR 35.714 Hz") == (35714.0, True)
        # a typo alone never becomes 35.7 MHz -- and never loads silently
        rate, unc = bruker._resolve_mas({}, "MASR 35714 kHz")
        assert rate == 35714.0 and unc is True
        # '0 kHz' is the operator declaring static (Phoenix setup EXPNOs)
        assert bruker._resolve_mas({}, "MASR 0 kHz") == (0.0, False)
        assert bruker._resolve_mas({"MASR": 4200}, "MASR 0 kHz") == (4200.0, True)
        # implausible both ways: dropped, the rest is flagged
        assert bruker._resolve_mas({"MASR": 4200}, "MASR 900 kHz") == (4200.0, True)


class TestMasSidecar:
    """_meta_1d records the three sources; _conflicts words them."""

    @staticmethod
    def _expno(tmp_path, booking, title="Rotor RS2427418\nMASR 20 kHz"):
        from test_referencing import _expno

        e = _expno(tmp_path / "2026-05", "04272026_P5-Bi8-12_SS_ALP", 3104,
                   "31P", 242.792156, 242.792156, "1777476884", title=title)
        if booking is not None:
            (e / "experiment_addenda.xml").write_text(
                "<experiment_addenda>\n<magic_angle_spinning>true"
                "</magic_angle_spinning>\n<magic_angle_spinning_rate>"
                f"{booking}</magic_angle_spinning_rate>\n</experiment_addenda>\n",
                encoding="utf-8")
        return e

    def test_meta_1d_and_conflict_wordings(self, tmp_path):
        acq = {"MASR": 4200, "NUC1": "<31P>"}
        title = "Rotor RS2427418\nMASR 20 kHz"
        e = self._expno(tmp_path, 22, title)
        meta = bruker._meta_1d(acq, title, e)
        assert meta["spin_rate_Hz"] == 22000.0 and meta["mas_uncertain"] is True
        assert meta["mas_source"] == "highest" and meta["rotor"] == "RS2427418"
        assert meta["mas_sources"] == {
            "acqus_Hz": 4200.0, "title_Hz": 20000.0, "title_note": "",
            "title_repaired": False, "title_raw": "MASR 20 kHz",
            "booking_Hz": 22000.0, "booking_flag": True}
        assert meta["masr_Hz"] == 4200.0            # the raw reading is kept
        (line,) = bruker._conflicts(meta, title)
        assert line.startswith("MAS rate: acqus says 4200 Hz, the title says "
                               "20000 Hz and the booking sidecar says 22000 Hz")
        assert line.endswith("-- confirm before fitting")

        # acqus outvoted: today's prefix, then the sidecar sentence
        e = self._expno(tmp_path / "b", "35.714", "MASR 35.714 kHz")
        meta = bruker._meta_1d(acq, "MASR 35.714 kHz", e)
        assert meta["spin_rate_Hz"] == pytest.approx(35714.0)
        assert meta["mas_uncertain"] is False and meta["mas_source"] == "title+booking"
        (line,) = bruker._conflicts(meta, "")
        assert line.startswith("MAS rate: acqus says 4200 Hz but the title says 35714 Hz")
        assert "the booking sidecar agrees -- acqus outvoted, using 35714 Hz" in line

        # two sources, no sidecar: the pre-sidecar sentence verbatim
        e = self._expno(tmp_path / "c", None, "MASR 35.714 kHz")
        meta = bruker._meta_1d(acq, "MASR 35.714 kHz", e)
        assert meta["mas_sources"]["booking_Hz"] is None
        assert meta["mas_sources"]["booking_flag"] is None
        assert bruker._conflicts(meta, "") == [
            "MAS rate: acqus says 4200 Hz but the title says 35714 Hz -- "
            "confirm before fitting"]

        # acqus vs booking, no title rate
        e = self._expno(tmp_path / "d", 10, "35Cl zg")
        meta = bruker._meta_1d({"MASR": 5000, "NUC1": "<35Cl>"}, "35Cl zg", e)
        assert bruker._conflicts(meta, "") == [
            "MAS rate: acqus says 5000 Hz but the booking sidecar says 10000 Hz"
            " -- confirm before fitting"]

        # title outvoted by acqus + booking
        e = self._expno(tmp_path / "e", 20, "MASR 22 kHz")
        meta = bruker._meta_1d({"MASR": 20000, "NUC1": "<31P>"}, "MASR 22 kHz", e)
        assert meta["spin_rate_Hz"] == 20000.0 and meta["mas_uncertain"] is False
        (line,) = bruker._conflicts(meta, "")
        assert line == ("MAS rate: the title says 22000 Hz but acqus says 20000 Hz"
                        " and the booking sidecar agrees -- title outvoted, "
                        "using 20000 Hz")

        # a typo note is its own entry; all agreeing -> nothing
        e = self._expno(tmp_path / "f", "35.714", "MASR 35741 kHz")
        meta = bruker._meta_1d(acq, "MASR 35741 kHz", e)
        assert meta["spin_rate_Hz"] == pytest.approx(35714.0)
        lines = bruker._conflicts(meta, "")
        assert any("outvoted" in ln for ln in lines)
        assert any('title says "MASR 35741 kHz", read as 35741 Hz' in ln
                   for ln in lines)
        e = self._expno(tmp_path / "g", 20, "MASR 20 kHz")
        meta = bruker._meta_1d({"MASR": 20000, "NUC1": "<31P>"}, "MASR 20 kHz", e)
        assert meta["mas_source"] == "all" and bruker._conflicts(meta, "") == []

    def test_summary_shows_booking_line_only_when_present(self):
        kw = dict(path="X", nucleus="27Al", sfo1_MHz=156.28, pulse_program="zg",
                  td=1024, sw_Hz=100000.0, masr_Hz=4200.0, title="t\nMASR 20kHz",
                  fid=None, processed=None, processed_ppm=None)
        plain = bruker.BrukerExperiment(**kw).summary.splitlines()
        assert len(plain) == 5 and plain[3] == "MASR (acqus): 4200.0 Hz"
        assert plain[4] == "title: t"
        with_side = bruker.BrukerExperiment(
            **kw, mas_sources={"booking_Hz": 22000.0}).summary.splitlines()
        assert len(with_side) == 6
        assert with_side[3] == "MASR (acqus): 4200.0 Hz"
        assert with_side[4] == "MAS (booking sidecar): 22000 Hz"
        assert with_side[5] == "title: t"


def test_read_2702_three_way_and_1903_outvoted():
    """Real anchors: 2702 (acqus 4200 / title 20kHz / booking 22) is flagged
    three-way; 1903 and the tutorial-4 LAW EXPNOs (title == booking 35.714,
    acqus 4200) are settled with acqus outvoted."""
    data = bruker.read(require(BRUKER_1R))
    m = data.meta
    assert m["spin_rate_Hz"] == 22000.0 and m["mas_uncertain"] is True
    assert m["mas_source"] == "highest" and m["rotor"] == "RS2427418"
    assert (m["mas_sources"]["acqus_Hz"], m["mas_sources"]["title_Hz"],
            m["mas_sources"]["booking_Hz"]) == (4200.0, 20000.0, 22000.0)
    mas_lines = [w for w in data.warnings if w.startswith("MAS rate:")]
    assert len(mas_lines) == 1 and "booking sidecar says 22000 Hz" in mas_lines[0]

    exp = bruker.read_expno(require(EXPNO_1903))
    assert exp.mas_sources["booking_Hz"] == pytest.approx(35714.0)
    assert exp.mas_sources["title_Hz"] == pytest.approx(35714.0)
    meta = bruker._meta_1d({"MASR": 4200, "NUC1": "<19F>"}, exp.title, EXPNO_1903)
    assert meta["spin_rate_Hz"] == pytest.approx(35714.0)
    assert meta["mas_uncertain"] is False and meta["mas_source"] == "title+booking"
    assert any("MAS rate" in c and "outvoted" in c for c in exp.conflicts)
    assert "MAS (booking sidecar): 35714 Hz" in exp.summary

    law = bruker.read_expno(require(LAW_CA_11B[0]))
    assert law.conflicts[0].startswith(
        "MAS rate: acqus says 4200 Hz but the title says 35714 Hz")
    assert "outvoted" in law.conflicts[0]
    assert law.mas_sources["booking_Hz"] == pytest.approx(35714.0)
