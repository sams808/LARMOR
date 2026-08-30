"""Smoke tests for the guided QCPMG dialog (offscreen).

These pin the WIRING, not the physics (that is tests/test_qcpmg.py): that the
six stages exist, the interactive echo-top pick stays in sync with its
spinbox, the matched-LB button writes the field it advertises, both spectra
can be overlaid at once, and the result carries its full provenance.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _synthetic_qcpmg(tmp_path, period=128, n_echoes=24, top=64, t2_echoes=6.0):
    """A Bruker-less stand-in: the dialog only needs a readable time-domain
    source, so drive it through the same public API the loader produces."""
    t = (np.arange(period) - top)
    echo = np.exp(-np.abs(t) / 9.0) * np.exp(2j * np.pi * 0.06 * t)
    fid = np.concatenate([echo * np.exp(-k / t2_echoes) for k in range(n_echoes)])
    return fid.astype(complex)


def _dialog_with(qapp, fid, period=128, top=64, sw=50000.0):
    from larmor.desktop.qcpmg_dialog import QcpmgDialog

    d = QcpmgDialog(None, None)
    d._loading = True
    d.fid = fid
    d.meta = {"sw_Hz": sw, "larmor_MHz": 78.0, "nucleus": "35Cl",
              "title": "synthetic", "expno": "1"}
    d._carrier, d._referenced = 0.0, True
    d._period_src = "manual"
    d.period.setValue(period)
    d.periodHz.setValue(sw / period)
    d.nEch.setMaximum(fid.size // period)
    d.nEch.setValue(fid.size // period)
    d.top.setMaximum(period - 1)
    d.top.setValue(top)
    d._loading = False
    d._recompute()
    for w in (d.btnCsv, d.btnSend):
        w.setEnabled(True)
    return d


def test_dialog_constructs_with_no_source(qapp):
    from larmor.desktop.qcpmg_dialog import QcpmgDialog

    d = QcpmgDialog(None, None)
    assert d.tabs.count() == 6
    titles = [d.tabs.tabText(i) for i in range(6)]
    assert titles[0].startswith("1") and titles[-1].startswith("6")
    # nothing to send before anything is loaded
    assert not d.btnSend.isEnabled() and not d.btnCsv.isEnabled()
    d.close()


def test_echo_top_line_and_spinbox_stay_in_sync(qapp, tmp_path):
    """The interactive pick replaces ssNake's 'type the number' step -- the
    draggable line and the numeric field are one value, and neither may
    recurse into an endless recompute."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.topLine.setValue(70.0)
    d._on_top_line()
    assert d.top.value() == 70
    d.top.setValue(55)
    assert int(round(float(d.topLine.value()))) == 55
    # a drag past the end of the echo is clamped, not wrapped
    d.topLine.setValue(9999.0)
    d._on_top_line()
    assert d.top.value() == d.period.value() - 1
    d.close()


def test_matched_lb_button_writes_the_lb_field(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    assert d.t2 is not None and d.t2.ok
    assert f"{d.t2.lb_Hz:,.1f}" in d.btnMatched.text()   # the value is on the button
    d.lb.setValue(0.0)
    d._use_matched()
    assert d.lb.value() == pytest.approx(round(d.t2.lb_Hz, 1), rel=1e-6)
    d.close()


def test_matched_lb_button_disabled_when_the_fit_fails(qapp, tmp_path):
    """A failed T2 must not offer a matched filter derived from nothing."""
    d = _dialog_with(qapp, np.zeros(128 * 12, complex))
    assert d.t2 is not None and not d.t2.ok
    assert not d.btnMatched.isEnabled()
    assert "did not converge" in d.lbl3.text()
    d.close()


def test_sum_echo_and_spikelets_can_be_shown_together(qapp, tmp_path):
    """ssNake's 'Superimposed with spikelets' figure: the old dialog had a
    mutually-exclusive mode combo, so this was structurally impossible."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.showSum.setChecked(True)
    d.showSpk.setChecked(True)
    assert len(d.p_spec.getPlotItem().listDataItems()) == 2
    assert "⊕" in d.lbl5.text()
    d.showSpk.setChecked(False)
    assert len(d.p_spec.getPlotItem().listDataItems()) == 1
    d.close()


def test_clicking_a_decay_point_excludes_it_and_refits(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    assert d.keep is not None and d.keep.all()
    before = d.t2.T2_s

    class _Pt:
        def data(self):
            return 3
    d._decay_clicked(None, [_Pt()])
    assert not d.keep[3]                       # excluded
    assert d.t2 is not None                    # and refitted
    assert d.t2.T2_s != before or d.t2.n_used < d.keep.size
    d.close()


def test_send_carries_full_provenance(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    got = {}
    d.accepted_1d.connect(lambda p, a, m: got.update(ppm=p, amp=a, meta=m))
    d._send()
    m = got["meta"]
    assert got["ppm"].size == got["amp"].size > 0
    for key in ("qcpmg_period_pts", "qcpmg_echo_top", "qcpmg_n_echoes",
                "qcpmg_lb_Hz", "qcpmg_gb_Hz", "qcpmg_p0_deg", "qcpmg_p1_deg",
                "qcpmg_carrier_ppm", "qcpmg_T2_s", "qcpmg_matched_lb_Hz"):
        assert key in m, key
    # a QCPMG sum echo is not spinning-sideband modulated: never inherit a
    # stale MAS rate from acqus
    assert m["spin_rate_Hz"] == 0.0 and m["mas_uncertain"] is False
    d.close()


def test_period_in_points_and_hz_stay_consistent(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path), period=128, sw=50000.0)
    d.periodHz.setValue(250.0)                 # 50000/250 = 200 points
    assert d.period.value() == 200
    d.period.setValue(100)                     # -> 500 Hz
    assert d.periodHz.value() == pytest.approx(500.0, rel=1e-6)
    d.close()


def test_send_flushes_a_pending_debounce(qapp, tmp_path):
    """Send inside the 150 ms debounce window used to ship the STALE spectrum
    while the meta recorded the NEW settings -- provenance that lies about how
    the data was processed."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    stale = d._spec.copy()                     # spectrum at LB = 0
    d.lb.setValue(800.0)                       # queues a recompute...
    assert d._timer.isActive()                 # ...which has NOT run yet
    got = {}
    d.accepted_1d.connect(lambda p, a, m: got.update(amp=a, meta=m))
    d._send()                                  # must flush, then emit
    assert not d._timer.isActive()
    assert got["meta"]["qcpmg_lb_Hz"] == 800.0
    assert not np.array_equal(got["amp"], stale)
    d.close()


def test_failed_second_load_keeps_nothing_sendable(qapp, tmp_path, monkeypatch):
    """A second _load that reads fine but cannot be processed (e.g. an aborted
    3-point fid) used to leave the OLD spectrum sendable under the NEW file's
    name; now nothing is sendable and the state is fully reset."""
    from larmor.io import bruker

    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    assert d.btnSend.isEnabled()

    class _Tiny:
        domain, ndim = "time", 1
        data = np.ones(3, complex)
        meta = {"sw_Hz": 50000.0, "title": "aborted", "expno": "2"}
    monkeypatch.setattr(bruker, "read", lambda _s: _Tiny())
    d._load(str(tmp_path / "2" / "fid"))
    assert "cannot process" in d.res.text()
    assert not d.btnSend.isEnabled() and not d.btnCsv.isEnabled()
    assert d._spec is None and d._spec_raw is None
    got = {}
    d.accepted_1d.connect(lambda p, a, m: got.update(meta=m))
    d._send()
    assert not got                             # nothing emitted
    d._copy_csv()                              # must not raise either
    d.close()


def test_period_change_resets_echo_bookkeeping(qapp, tmp_path):
    """Exclusions and the echo count are indices INTO one particular split;
    after a period change they would point at different stretches of data, so
    they are reset (excluded={3} at 128 pts is not echo 3 at 64 pts)."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))

    class _Pt:
        def data(self):
            return 3
    d._decay_clicked(None, [_Pt()])
    assert d.excluded == {3}
    d.dropFirst.setValue(1)
    d.period.setValue(64)
    assert d.excluded == set()
    assert d.nEch.value() == d.fid.size // 64
    assert d.dropFirst.value() == 0
    d.close()


def test_spinboxes_do_not_react_per_keystroke(qapp):
    """Typing '293' digit by digit must not act on the intermediate '29' --
    _on_period_pts would clamp the echo top to 28 and the recompute would
    silently run from the wrong point."""
    from PySide6.QtWidgets import QDoubleSpinBox, QSpinBox

    from larmor.desktop.qcpmg_dialog import QcpmgDialog
    d = QcpmgDialog(None, None)
    boxes = d.findChildren(QSpinBox) + d.findChildren(QDoubleSpinBox)
    assert boxes
    assert all(not b.keyboardTracking() for b in boxes)
    d.close()


def test_period_hz_field_always_shows_a_realisable_spacing(qapp, tmp_path):
    """Typing 392 Hz quantises to the same 128-point period; the field must
    snap back to the 390.625 Hz that period actually realises, not keep
    showing a spacing the data cannot have."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path), period=128, sw=50000.0)
    d.periodHz.setValue(392.0)
    assert d.period.value() == 128
    # the field shows 2 decimals, so the snap lands on 390.62
    assert d.periodHz.value() == pytest.approx(390.625, abs=0.011)
    d.close()


def test_copy_csv_includes_both_time_axes(qapp, tmp_path):
    """The CSV must never say a bare 'T2' -- the ssNake pseudo-axis value and
    the physical one differ by the echo length and are both in the wild."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d._copy_csv()
    text = QApplication.clipboard().text()
    assert "echo,t_s,intensity,kept" in text
    assert "T2_physical_s=" in text and "T2_ssnake_D1_s=" in text
    assert "LB_physical_Hz=" in text and "LB_ssnake_D1_Hz=" in text
    assert "period_pts=" in text and "echo_top=" in text
    d.close()


def test_dialog_fits_the_screen_and_can_be_shrunk(qapp, tmp_path):
    """The fixed 1180x820 spilled past a laptop screen and the layout minimum
    made shrinking impossible -- the dialog must open inside the available
    screen and accept being resized well below its default."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.show()
    qapp.processEvents()
    avail = d.screen().availableGeometry()
    assert d.width() <= avail.width() and d.height() <= avail.height()
    assert d.isSizeGripEnabled()
    d.resize(760, 540)
    qapp.processEvents()
    assert d.width() <= 800 and d.height() <= 620
    d.close()


def test_dialog_shrinks_far_below_the_old_floor(qapp, tmp_path):
    """Stages sit in scroll areas now, so no layout minimum can stop the user
    from making the window as small as they like (content scrolls instead)."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.show()
    qapp.processEvents()
    d.resize(520, 400)
    qapp.processEvents()
    assert d.width() <= 560 and d.height() <= 440
    d.close()


def test_dialog_remembers_the_size_the_user_chose(qapp, tmp_path, monkeypatch):
    """The size/position the dialog is closed at comes back on the next open
    (gated on LARMOR_NO_SESSION so tests do not pollute real settings)."""
    from PySide6.QtCore import QSettings

    s = QSettings("LARMOR", "app")
    old = s.value("qcpmgDialogGeometry")
    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    try:
        d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
        d.show(); qapp.processEvents()
        d.resize(804, 604); qapp.processEvents()
        d.done(0)
        d2 = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
        d2.show(); qapp.processEvents()
        # the platform may re-clamp by a few frame pixels; the point is that
        # the user's chosen size comes back, not the built-in default
        assert abs(d2.width() - 804) <= 10 and abs(d2.height() - 604) <= 10
        d2.done(0)
    finally:
        if old is None:
            s.remove("qcpmgDialogGeometry")
        else:
            s.setValue("qcpmgDialogGeometry", old)


def test_fields_dialog_table_is_visible_at_default_size(qapp):
    """The infinite-field dialog's per-field table was squeezed to ~1.5 rows
    by the plot's fixed minimum; it must get real height from the start."""
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    d = QcpmgFieldsDialog(None, "35Cl")
    d.show()
    qapp.processEvents()
    assert d.table.height() >= 140
    d.close()


def test_figure_package_builds_and_exports_three_formats(qapp, tmp_path,
                                                         monkeypatch):
    """'Export figure package' writes one publication-layout figure of the
    whole workflow as .png + .svg + .pdf next to wherever the user points."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    fig = d._package_figure()
    assert len(fig.axes) == 4                     # train / decay / spectrum / band
    base = tmp_path / "pkg" ; (tmp_path / "pkg").mkdir()
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getSaveFileName",
        staticmethod(lambda *a, **k: (str(base / "qcpmg_fig"), "")))
    d._export_package()
    for ext in ("png", "svg", "pdf"):
        f = base / f"qcpmg_fig.{ext}"
        assert f.exists() and f.stat().st_size > 5_000, ext
    assert "figure package written" in d.res.text()
    d.close()


def test_large_p1_gets_an_honest_warning(qapp, tmp_path):
    """A top BETWEEN samples legitimately needs |p1| up to ~180 deg (it is
    the exact compensation of a sub-dwell shift); beyond ~200 deg the period
    or top is genuinely off and the stage-5 label must say so (the symptom
    the user sees is dips flanking the line)."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d._phased = True
    d.p1.setValue(150.0)                  # legitimate sub-dwell compensation
    d._redraw_spec()
    assert "period or top is off" not in d.lbl5.text()
    d.p1.setValue(450.0)                  # beyond any sub-dwell explanation
    d._redraw_spec()
    assert "period or top is off" in d.lbl5.text()
    d.p1.setValue(0.0)
    d._redraw_spec()
    assert "period or top is off" not in d.lbl5.text()
    d.close()


def test_save_as_dataset_roundtrips_through_the_reader(qapp, tmp_path,
                                                       monkeypatch):
    """Stage 5 'Save as dataset…' writes a LARMOR .csv that loads back with
    the same axis, data and nucleus."""
    from larmor.io import spectra

    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    out = tmp_path / "qcpmg_1_sumecho.csv"
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getSaveFileName",
        staticmethod(lambda *a, **k: (str(out), "")))
    d._save_dataset()
    assert out.exists()
    ppm, amp, meta = spectra.read_csv(out)
    assert meta.get("nucleus") == "35Cl"
    assert ppm.size == d._ppm.size
    assert np.allclose(np.sort(ppm), np.sort(d._ppm))
    d.close()


def test_send_to_infinite_field_accumulates_across_sessions(qapp, tmp_path):
    """Stage 6 '→ infinite-field δiso…': the shared dialog keeps rows from
    one processing session to the next, which is the whole point — process
    field 1, send; process field 2, send; Compute."""
    from larmor.desktop import qcpmg_fields_dialog as qfd

    qfd._shared = None                       # isolate from other tests
    try:
        d1 = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
        d1._autophase()
        d1._send_infinite()
        assert "sent to infinite-field" in d1.res.text()
        shared = qfd._shared
        assert shared is not None
        rows_after_first = shared.table.rowCount()
        d1.close()

        d2 = _dialog_with(qapp, _synthetic_qcpmg(tmp_path, t2_echoes=4.0))
        d2.meta["larmor_MHz"] = 160.0        # "the other field"
        d2._autophase()
        d2._send_infinite()
        assert qfd._shared is shared          # same instance, accumulated
        assert shared.table.rowCount() == rows_after_first + 1
        # both dataset rows are supervisable (spectrum attached)
        ds_rows = [r for r in range(shared.table.rowCount())
                   if shared._row_ds_id(r) is not None]
        assert len(ds_rows) == 2
        shared.close()
    finally:
        qfd._shared = None


def test_load_verifies_a_guessed_period_against_the_data(qapp, tmp_path,
                                                         monkeypatch):
    """The 81Br lesson: CNST7 unset, so the period was ASSUMED
    rotor-synchronised from MASR — and was wrong by a factor of ~2.7. A
    guessed (not recorded) period is now verified against the data's own
    echo-repeat correlation and replaced when the data disagrees."""
    from larmor.io import bruker

    per = 250
    t = np.arange(per) - 125
    echo = np.exp(-np.abs(t) / 20.0) * np.exp(2j * np.pi * 0.03 * t)
    rng = np.random.default_rng(5)
    fid = np.concatenate([echo * np.exp(-k / 6.0) for k in range(24)])
    fid = fid + rng.normal(0, 0.01, fid.size) + 1j * rng.normal(0, 0.01, fid.size)

    class _Fake:
        domain, ndim = "time", 1
        data = fid
        # MASR would suggest sw/masr = 100 pts -- wrong by 2.5x
        meta = {"sw_Hz": 50000.0, "larmor_MHz": 78.0, "nucleus": "35Cl",
                "masr_Hz": 500.0, "title": "t", "expno": "9"}
    monkeypatch.setattr(bruker, "read", lambda _s: _Fake())
    from larmor.desktop.qcpmg_dialog import QcpmgDialog
    d = QcpmgDialog(None, None)
    d._load(str(tmp_path / "9" / "fid"))
    assert d.period.value() == per
    assert d._period_src == "echo correlation"
    assert "MEASURED from the data" in d.lbl1.text()
    d.close()


def test_magnitude_mode_needs_no_phase_and_unlocks_stage6(qapp, tmp_path):
    """'magnitude (mc)': the phase controls go dead (they mean nothing), the
    trace becomes |spectrum| with its floor removed, and stage 6 measures
    without demanding a phase first."""
    from larmor import qcpmg

    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d._phased = False
    d._measure()
    assert d._cg is None and "phase the spectrum first" in d.lbl6.text()

    d.magMode.setChecked(True)
    assert not d.p0.isEnabled() and not d.btnAutophase.isEnabled()
    assert not d.p2.isEnabled()
    assert "magnitude" in d.showSum.text()      # never call it "absorption"
    assert np.allclose(d._spec, qcpmg.magnitude_spectrum(d._spec_raw))
    assert "magnitude (mc)" in d.lbl5.text()
    d._measure()
    assert d._cg is not None                      # measured without phasing
    assert "magnitude spectrum" in d.lbl6.text()

    # a phase change must not alter a magnitude trace
    before = d._spec.copy()
    d.p0.setValue(55.0)
    d._rephase()
    assert np.allclose(d._spec, before)

    d.magMode.setChecked(False)
    assert d.p0.isEnabled() and d.btnAutophase.isEnabled() and d.p2.isEnabled()
    assert "absorption" in d.showSum.text()
    d.close()


def test_magnitude_mode_is_recorded_and_reset(qapp, tmp_path, monkeypatch):
    """The mode travels with the result (provenance + CSV) and does not leak
    into the next dataset."""
    from PySide6.QtWidgets import QApplication, QMessageBox
    from larmor.io import bruker

    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.magMode.setChecked(True)
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: QMessageBox.Yes))
    got = {}
    d.accepted_1d.connect(lambda p, a, m: got.update(meta=m))
    d._send()
    assert got["meta"]["qcpmg_magnitude"] is True
    d._copy_csv()
    assert "spectrum_mode=magnitude(mc)" in QApplication.clipboard().text()

    class _Other:
        domain, ndim = "time", 1
        data = _synthetic_qcpmg(tmp_path)
        meta = {"sw_Hz": 50000.0, "larmor_MHz": 78.0, "nucleus": "35Cl",
                "title": "next", "expno": "2"}
    monkeypatch.setattr(bruker, "read", lambda _s: _Other())
    d._load(str(tmp_path / "2" / "fid"))
    assert not d.magMode.isChecked()               # a fresh sample, fresh mode
    assert d.p0.isEnabled()
    d.close()


def test_second_order_phase_is_offered_and_applied(qapp, tmp_path):
    """Stage 5 exposes p2 and applies it; Autophase fills it only when the
    spectrum needs it (a swept-pulse chirp), leaving it 0 otherwise."""
    from larmor import qcpmg

    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.p2.setValue(120.0)
    d._rephase()
    assert np.allclose(d._spec, qcpmg.phase_spectrum(
        d._spec_raw, d.p0.value(), d.p1.value(), p2_deg=120.0).real)
    assert "p2 120" in d.lbl5.text()
    assert d._phased                             # p2 alone counts as phased

    d.p2.setValue(0.0)
    d._autophase()                               # an ordinary echo
    assert d.p2.value() == 0.0
    d.close()


def test_infinite_field_handoff_records_the_mode(qapp, tmp_path):
    """A magnitude δcg and an absorption δcg are different observables; the
    multi-field table must know which is which and say so."""
    from larmor.desktop import qcpmg_fields_dialog as qfd

    qfd._shared = None
    try:
        d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
        d._autophase()
        d._send_infinite()                       # absorption
        d.magMode.setChecked(True)
        d.meta["larmor_MHz"] = 160.0
        d._send_infinite()                       # magnitude
        shared = qfd._shared
        modes = [v.get("magnitude") for v in shared._ds.values()]
        assert modes == [False, True]
        assert "mix magnitude" in shared.wresult.text()
        shared.close()
        d.close()
    finally:
        qfd._shared = None


def test_failed_recompute_leaves_no_stale_trace(qapp, tmp_path):
    """After a failure the plots and stage labels must not keep showing the
    last good spectrum while the controls describe something else."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    assert d._spec is not None
    d.period.setValue(d.fid.size * 2)            # a period longer than the train
    d._flush()
    assert d._spec is None and d._spec_raw is None
    assert not d.btnSend.isEnabled()
    assert d.lbl5.text() == "" and d.lbl6.text() == ""
    assert not d.p_spec.getPlotItem().listDataItems()
    d.magMode.setChecked(True)                   # must not redraw dead data
    assert d._spec is None
    d.close()
