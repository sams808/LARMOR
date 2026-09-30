"""The pybaselines dialog (Process > Baseline > pybaselines): it builds
offscreen for every method, params() carries the method key and every
parameter, the preview on the shipped 27Al spectrum gives a baseline of the
right length, an invalid value shows in the status label without raising,
the log slider and its spin box stay in step -- and the recorded step
survives the processing panel's next Apply (the panel carries it), both on
the panel alone and through MainWindow.apply_pybaseline."""
import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from larmor import pybaseline as pyb  # noqa: E402  (Qt-free)

REPO = Path(__file__).resolve().parents[1]
EXPNO_3616 = REPO / "examples" / "pCABS2-4" / "3616"     # 27Al zg, 2048 points


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def al27():
    from larmor.loader import load_any

    ppm, amp, _recipe, _meta, _warnings = load_any(EXPNO_3616)
    return np.asarray(ppm, float), np.asarray(amp, float)


def _dialog(ppm, amp, method=None):
    from larmor.desktop.pybaseline_dialog import PybaselineDialog

    return PybaselineDialog(None, ppm, amp, method=method)


@pytest.mark.parametrize("key", list(pyb.METHODS))
def test_dialog_builds_for_every_method_with_complete_params(qapp, al27, key):
    ppm, amp = al27
    dlg = _dialog(ppm, amp, key)
    try:
        assert dlg.method() == key and dlg.method_label() == pyb.METHODS[key].label
        p = dlg.params()
        spec = pyb.METHODS[key]
        assert p["method"] == key
        assert set(p) == {"method", *(q.key for q in spec.params)}
        # it opens on the length-scaled defaults, and what it records is
        # valid for the op
        assert p == {"method": key, **pyb.default_params(key, amp.size)}
        pyb.validate_params(key, {k: v for k, v in p.items() if k != "method"}, amp.size)
        assert spec.description in dlg.description.text()
        assert dlg._factor == 1                    # 2048 points: previewed in full
        dlg._preview()
        assert dlg.c_base.yData is not None and len(dlg.c_base.yData) == amp.size
        assert np.isfinite(dlg.c_base.yData).all()
        assert len(dlg.c_corr.yData) == amp.size and len(dlg.c_raw.xData) == amp.size
        assert np.allclose(dlg.c_corr.yData, amp - dlg.c_base.yData)
        assert spec.label in dlg.status.text() and "failed" not in dlg.status.text()
        assert dlg.btn_apply.isEnabled()
    finally:
        dlg.close()


def test_preview_on_3616_and_the_decimated_64k_preview(qapp, al27):
    ppm, amp = al27
    dlg = _dialog(ppm, amp)
    try:
        assert dlg.method() == "arpls"              # the default method
        dlg._preview()
        base = np.asarray(dlg.c_base.yData)
        assert base.shape == amp.shape
        assert np.abs(base).max() < 0.5 * np.abs(amp).max()   # a baseline, not the peaks
        # a spectrum's worth of the recorded step
        assert dlg.params() == {"method": "arpls", "lam": 1e5, "diff_order": 2,
                                "max_iter": 50}
    finally:
        dlg.close()
    # a 64k-point spectrum is previewed on every 8th point with rescaled
    # parameters; the recorded parameters stay the full-length ones
    x = np.linspace(500.0, -500.0, 65536)
    y = np.interp(x, ppm, amp)
    dlg2 = _dialog(x, y, "rolling_ball")
    try:
        assert dlg2._factor == 8 and dlg2._py.size == 8192
        dlg2._preview()
        assert len(dlg2.c_base.yData) == 8192 and len(dlg2.c_raw.xData) == 8192
        assert "every 8th point (8192 of 65536)" in dlg2.status.text()
        assert "Apply uses all points" in dlg2.status.text()
        assert dlg2.params()["half_window"] == pyb.suggest_half_window(65536) == 2048
    finally:
        dlg2.close()


def test_invalid_value_shows_in_the_status_label_without_raising(qapp, al27):
    ppm, amp = al27
    dlg = _dialog(ppm, amp, "snip")
    try:
        w = dlg._widgets["max_half_window"]
        assert w.maximum() == amp.size // 2          # the spin box itself clamps
        w.setMaximum(10 ** 6)                        # so hand it a bad value
        w.setValue(5000)
        dlg._timer.stop()
        dlg._preview()                               # no exception
        assert "SNIP: max half-window (points) must be between 1 and 1024 (got 5000)" \
            in dlg.status.text()
        assert not dlg.btn_apply.isEnabled()
        assert dlg.c_base.yData is None or len(dlg.c_base.yData) == 0   # no stale curve
        w.setValue(64)
        dlg._preview()
        assert dlg.btn_apply.isEnabled() and "SNIP" in dlg.status.text()
        assert len(dlg.c_base.yData) == amp.size
        # a failure of pybaselines itself is reported the same way
        from larmor.desktop import pybaseline_dialog as mod

        real = mod.pyb.compute

        def boom(*a, **k):
            raise RuntimeError("solver blew up")
        mod.pyb.compute = boom
        try:
            dlg._preview()
        finally:
            mod.pyb.compute = real
        assert "SNIP failed: solver blew up" in dlg.status.text()
        assert not dlg.btn_apply.isEnabled()
    finally:
        dlg.close()


def test_log_slider_spinbox_reset_and_method_switch(qapp, al27):
    from PySide6.QtGui import QValidator

    ppm, amp = al27
    dlg = _dialog(ppm, amp, "arpls")
    try:
        sb, sl = dlg._widgets["lam"], dlg._sliders["lam"]
        assert sb.value() == 1e5 and sl.value() == 50
        assert sl.minimum() == 20 and sl.maximum() == 120     # 1e2 .. 1e12, ticks per decade
        sl.setValue(70)
        assert sb.value() == pytest.approx(1e7) and dlg.params()["lam"] == pytest.approx(1e7)
        sb.setValue(3e4)
        assert sl.value() == round(np.log10(3e4) * 10)
        sb.stepBy(1)
        assert sb.value() == pytest.approx(6e4)
        assert sb.textFromValue(1e5) == "1e+05" and sb.textFromValue(250.0) == "250"
        assert sb.valueFromText("2e6") == 2e6
        assert sb.validate("1e5", 3)[0] == QValidator.State.Acceptable
        assert sb.validate("1e", 2)[0] == QValidator.State.Intermediate
        assert sb.validate("abc", 3)[0] == QValidator.State.Invalid
        dlg._widgets["diff_order"].setValue(3)
        dlg.reset_defaults()
        assert dlg.params() == {"method": "arpls", "lam": 1e5, "diff_order": 2,
                                "max_iter": 50}
        assert sl.value() == 50
        # switching the method rebuilds the form on that method's defaults
        dlg.set_method("penalized_poly")
        assert dlg.method() == "penalized_poly"
        assert set(dlg._widgets) == {"poly_order", "cost_function", "max_iter"}
        assert dlg.params()["cost_function"] == "asymmetric_truncated_quadratic"
        cf = dlg._widgets["cost_function"]
        cf.setCurrentIndex(cf.findData("symmetric_huber"))
        assert dlg.params()["cost_function"] == "symmetric_huber"
        assert pyb.METHODS["penalized_poly"].description in dlg.description.text()
        dlg._timer.stop()
        dlg._preview()
        assert "Penalized poly" in dlg.status.text()
        # a family heading made current by code (the popup never offers one)
        # bounces back to the method, form and edits intact
        assert dlg.combo.itemData(0) is None
        row = dlg.combo.currentIndex()
        dlg.combo.setCurrentIndex(0)
        assert dlg.method() == "penalized_poly" and dlg.combo.currentIndex() == row
        assert dlg.params()["cost_function"] == "symmetric_huber"
        # a bool parameter records a Python bool
        dlg.set_method("snip")
        dlg._widgets["decreasing"].setChecked(True)
        assert dlg.params()["decreasing"] is True
    finally:
        dlg.close()


def test_apply_accepts_and_cancel_rejects(qapp, al27):
    from PySide6.QtWidgets import QDialog

    ppm, amp = al27
    dlg = _dialog(ppm, amp)
    dlg._preview()
    dlg.btn_apply.click()
    assert dlg.result() == QDialog.Accepted
    dlg2 = _dialog(ppm, amp)
    dlg2.bb.rejected.emit()
    assert dlg2.result() == QDialog.Rejected


def test_panel_carries_a_recorded_pybaseline_step_through_its_next_apply(qapp):
    """The panel emits the ABSOLUTE chain from its widgets: a step it has no
    widget for survives only because sync_from_ops carries it and _emit
    re-appends it after the widgets' chain."""
    from larmor.desktop.panels import ProcessingPanel

    p = ProcessingPanel()
    step = {"op": "pybaseline", "method": "arpls", **pyb.default_params("arpls", 2048)}
    chain = [{"op": "phase", "p0": 12.0, "p1": 0.0}, step]
    assert p.sync_from_ops(chain, False) is True
    got = []
    p.apply_requested.connect(lambda ops, raw: got.append((ops, raw)))
    p.p0v.setValue(20.0)
    p.btnApply.click()
    ops, raw = got[-1]
    assert raw is False
    assert [o["op"] for o in ops] == ["phase", "pybaseline"]
    assert ops[0]["p0"] == 20.0 and ops[1] == step


def test_main_window_action_appends_the_step_and_the_panel_keeps_it(qapp, monkeypatch):
    """apply_pybaseline appends to the recorded chain (an earlier phase is not
    discarded), records it in the recipe and syncs the panel, so the panel's
    next Apply keeps the step and the spectrum. Needs a MainWindow: the
    append and the sync live in the action itself."""
    from PySide6.QtWidgets import QDialog, QMessageBox

    from larmor.desktop.app import MainWindow
    from larmor.desktop.pybaseline_dialog import PybaselineDialog

    popups = []
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: popups.append(a)))
    monkeypatch.setattr(PybaselineDialog, "exec", lambda self: QDialog.Accepted)
    win = MainWindow()
    try:
        win.load_source(str(EXPNO_3616), keep_fit=False)
        qapp.processEvents()
        n = win.exp_amp.size
        win.apply_processing([{"op": "phase", "p0": 5.0, "p1": 0.0}], False)
        qapp.processEvents()
        phased = win.exp_amp.copy()
        win.apply_pybaseline("airpls")
        qapp.processEvents()
        ops = win.recipe["processing"]
        assert [o["op"] for o in ops] == ["phase", "pybaseline"]
        assert ops[1] == {"op": "pybaseline", "method": "airpls",
                          **pyb.default_params("airpls", n)}
        assert ops[1]["lam"] == pyb.suggest_lam(n, 2, 1e6)
        assert not np.array_equal(phased, win.exp_amp)
        corrected = win.exp_amp.copy()
        assert win.proc_panel._carried_ops == [ops[1]]
        # the panel's next Apply (widgets untouched) keeps the step and the spectrum
        win.proc_panel.btnApply.click()
        qapp.processEvents()
        assert [o["op"] for o in win.recipe["processing"]] == ["phase", "pybaseline"]
        assert np.allclose(win.exp_amp, corrected)
        # and a second, different correction stacks on top of the first
        win.apply_pybaseline("mor")
        qapp.processEvents()
        names = [o["op"] for o in win.recipe["processing"]]
        assert names == ["phase", "pybaseline", "pybaseline"]
        assert win.recipe["processing"][2]["method"] == "mor"
        assert "applied" in win.statusBar().currentMessage()
        # the action passes a QAction's checked flag through as "no method"
        win.apply_pybaseline(False)
        qapp.processEvents()
        assert win.recipe["processing"][3]["method"] == "arpls"
        assert not popups
    finally:
        win.close()
