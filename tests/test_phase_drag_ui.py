"""F2 -- TopSpin drag-to-phase: offscreen tests for the mode.

Three layers, nothing skips and no real data is needed: the AnchoredViewBox
routing and the SpectrumView gesture alone (duck-typed MouseDragEvents), the
ProcessingPanel's set_phase / sync_phase_from alone, and MainWindow flows on a
synthetic CSV spectrum (one undo entry per drag, drag == typing, live off,
pivot re-expression, mutual exclusion, the FID-display gate). Kept out of
tests/test_desktop.py to shrink the shared-file footprint. The window is
never shown, so assertions use isHidden() (explicit state), never isVisible().
"""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

import pyqtgraph as pg  # noqa: E402
from PySide6.QtCore import QEvent, Qt  # noqa: E402
from PySide6.QtGui import QKeyEvent, QKeySequence  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from larmor.phasedrag import wrap_p0  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")  # never inherit a real session
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


class _Drag:
    """Duck-typed pyqtgraph MouseDragEvent (scene coordinates = pixels)."""

    def __init__(self, button=Qt.LeftButton, down=(100.0, 100.0),
                 pos=(100.0, 100.0), last=None, start=False, finish=False,
                 modifiers=Qt.NoModifier):
        self._button = button
        self._down = pg.Point(*down)
        self._pos = pg.Point(*pos)
        self._last = pg.Point(*(last if last is not None else down))
        self._start, self._finish = start, finish
        self._mods = modifiers
        self._accepted = False

    def scenePos(self):
        return self._pos

    def lastScenePos(self):
        return self._last

    def buttonDownScenePos(self, *_):
        return self._down

    def button(self):
        return self._button

    def modifiers(self):
        return self._mods

    def isStart(self):
        return self._start

    def isFinish(self):
        return self._finish

    def accept(self):
        self._accepted = True

    def isAccepted(self):
        return self._accepted


def _synthetic_xy():
    """Two Lorentzian absorption lines (a real dispersion part under Hilbert)."""
    x = np.linspace(-50.0, 50.0, 2048)
    amp = 1.0 / (1 + (x / 1.5) ** 2) + 0.4 / (1 + ((x - 20.0) / 2.5) ** 2)
    return x, amp


def _load(win, qapp, tmp_path):
    """Write the synthetic spectrum as a LARMOR CSV and open it for real, so
    source_path, the recipe (larmor_frequency_MHz, processing=[]) and the
    view are exactly what a user's load leaves behind."""
    path = tmp_path / "phase.csv"
    x, amp = _synthetic_xy()
    with open(path, "w", encoding="utf-8") as f:
        f.write("# LARMOR spectrum\n# nucleus=31P\n# larmor_MHz=162.0\n"
                "# spin_rate_Hz=10000.0\n")
        f.write("ppm,intensity\n")
        for xi, yi in zip(x, amp):
            f.write(f"{xi},{yi}\n")
    win.load_source(str(path), keep_fit=False)
    qapp.processEvents()
    assert win.central_stack.currentWidget() is win.view
    assert win.exp_ppm.size == x.size
    return win.exp_ppm.copy(), win.exp_amp.copy()


def _expected(win, base_x, base_amp, p0, p1, frac):
    """What typing p0/p1 into the panel produces: the panel's chain
    [hilbert, phase] applied to the unprocessed base, sorted like
    apply_processing sorts."""
    from larmor import processing as proc

    s = proc.apply(proc.from_processed(base_x, base_amp,
                                       win.recipe["larmor_frequency_MHz"]),
                   [{"op": "hilbert"},
                    {"op": "phase", "p0": p0, "p1": p1, "pivot_frac": frac}])
    return s.y.real[np.argsort(s.x_ppm)]


def _process_action(win, prefix):
    """The Process-menu action (any depth: Phase ▸, Baseline ▸ …) whose label
    starts with ``prefix``."""
    from PySide6.QtWidgets import QMenu

    def walk(menu):
        for a in menu.actions():
            if a.text().startswith(prefix):
                return a
            sub = next((m for m in menu.findChildren(QMenu)
                        if m.menuAction() == a), None)
            if sub is not None:
                hit = walk(sub)
                if hit is not None:
                    return hit
        return None

    for top in win.menuBar().actions():
        if top.text() == "&Process":
            hit = walk(top.menu())
            if hit is not None:
                return hit
    raise AssertionError(f"no Process action starting with {prefix!r}")


# ============================================================== SpectrumView
def test_viewbox_routes_only_left_canvas_drags_without_ctrl_to_the_handler(
        qapp, monkeypatch):
    from larmor.desktop.plot import AnchoredViewBox, SpectrumView

    view = SpectrumView()
    vb = view.getPlotItem().getViewBox()
    assert isinstance(vb, AnchoredViewBox)
    base_calls = []
    monkeypatch.setattr(pg.ViewBox, "mouseDragEvent",
                        lambda self, ev, axis=None: base_calls.append(axis))

    # mode off: every drag pans/zooms as before
    assert vb.phase_drag_handler is None
    vb.mouseDragEvent(_Drag())
    assert base_calls == [None]

    view.set_phase_drag_mode(True)
    assert vb.phase_drag_handler is not None
    assert view.phase_drag_active()
    assert vb.phase_drag_takes(_Drag())
    assert not vb.phase_drag_takes(_Drag(modifiers=Qt.ControlModifier))
    assert not vb.phase_drag_takes(_Drag(), axis=0)          # axis-item drag
    assert not vb.phase_drag_takes(_Drag(button=Qt.RightButton))
    assert not vb.phase_drag_takes(_Drag(button=Qt.MiddleButton))
    ev = _Drag(start=True)
    vb.mouseDragEvent(ev)
    assert base_calls == [None] and ev.isAccepted()          # base NOT called
    vb.mouseDragEvent(_Drag(button=Qt.MiddleButton))
    assert base_calls == [None, None]                        # middle still pans
    assert view.cursor().shape() == Qt.SizeAllCursor

    view.set_phase_drag_mode(False)
    assert vb.phase_drag_handler is None
    assert not view.phase_drag_active()
    assert view.cursor().shape() == Qt.ArrowCursor
    vb.mouseDragEvent(_Drag())
    assert len(base_calls) == 3


def test_view_emits_cumulative_degrees_with_axis_lock_and_release(qapp):
    from larmor.desktop.plot import SpectrumView

    view = SpectrumView()
    view.set_phase_drag_mode(True)
    moves, released = [], []
    view.phase_dragged.connect(lambda a, b: moves.append((a, b)))
    view.phase_drag_released.connect(lambda: released.append(1))

    view._phase_drag(_Drag(down=(100, 100), last=(100, 100), pos=(140, 102),
                           start=True))
    assert moves[-1] == (10.0, 0.0)                         # 40 px right = +10 deg
    view._phase_drag(_Drag(down=(100, 100), last=(140, 102), pos=(140, 60)))
    assert moves[-1] == (10.0, 0.0)                         # locked to p0
    view._phase_drag(_Drag(down=(100, 100), last=(140, 60), pos=(140, 60),
                           finish=True))
    assert released == [1]

    # a new gesture: up = +p1, Shift = fine, applied per event
    view._phase_drag(_Drag(down=(100, 100), last=(100, 100), pos=(101, 50),
                           start=True))
    assert moves[-1] == (0.0, 50.0)
    view._phase_drag(_Drag(down=(100, 100), last=(101, 50), pos=(101, 40),
                           modifiers=Qt.ShiftModifier))
    assert moves[-1] == pytest.approx((0.0, 51.0))
    view._phase_drag(_Drag(down=(100, 100), last=(101, 40), pos=(101, 40),
                           finish=True))
    assert released == [1, 1]

    # below the lock threshold nothing is emitted, the release still is
    n = len(moves)
    view._phase_drag(_Drag(down=(100, 100), last=(100, 100), pos=(102, 100),
                           start=True))
    view._phase_drag(_Drag(down=(100, 100), last=(102, 100), pos=(102, 100),
                           finish=True))
    assert len(moves) == n and released == [1, 1, 1]


def test_pivot_survives_dock_hide_while_phasing_and_reports_moves(qapp):
    from larmor.desktop.plot import SpectrumView

    view = SpectrumView()
    view.set_phase_drag_mode(True)                 # armed before any data
    assert view._pivot is None                     # nothing to pin yet
    x, amp = _synthetic_xy()
    view.set_experiment(x, amp)                    # data arrives in the mode
    assert view._pivot is not None
    view.set_phase_drag_mode(False)
    view.show_phase_pivot(False)
    assert view._pivot is None
    view.set_experiment(x, amp)
    assert view._pivot is None                     # mode off: the dock governs

    view.set_phase_drag_mode(True)
    assert view._pivot is not None                 # created without any dock
    view.show_phase_pivot(False)                   # the dock's hide: ignored
    assert view._pivot is not None
    f0 = view.phase_pivot_frac()
    got = []
    view.pivot_moved.connect(got.append)
    view._pivot.setValue(20.0)
    view._pivot.sigPositionChangeFinished.emit(view._pivot)
    assert len(got) == 1
    assert got[0] == view.phase_pivot_frac() and got[0] != f0
    view.set_phase_drag_mode(False)
    view.show_phase_pivot(False)
    assert view._pivot is None


def test_fid_display_suspends_the_gesture_and_the_spectrum_restores_it(qapp):
    """F6 handoff: no phasing while the FID is displayed -- a left drag pans
    the FID, and the ppm axis brings the gesture back."""
    from larmor.desktop.plot import SpectrumView

    view = SpectrumView()
    x, amp = _synthetic_xy()
    view.set_experiment(x, amp)
    view.set_phase_drag_mode(True)
    vb = view.getPlotItem().getViewBox()
    assert vb.phase_drag_handler is not None
    t = np.linspace(0.0, 10.0, 256)
    view.set_fid(t, np.exp(-t), "FID (real)")
    assert view.domain == "time"
    assert vb.phase_drag_handler is None           # pans, does not phase
    assert not vb.phase_drag_takes(_Drag())
    assert view.cursor().shape() == Qt.ArrowCursor
    assert view.phase_drag_active()                # the mode itself is kept
    view.set_experiment(x, amp)
    assert view.domain == "freq"
    assert vb.phase_drag_handler is not None
    assert view.cursor().shape() == Qt.SizeAllCursor


# =========================================================== ProcessingPanel
def test_panel_set_phase_wraps_clamps_applies_once_and_sync_is_silent(qapp):
    from larmor.desktop.panels import ProcessingPanel

    p = ProcessingPanel()
    got = []
    p.apply_requested.connect(lambda ops, raw: got.append((ops, raw)))
    p.chkHilbert.setChecked(True)

    p.set_phase(190.0, 800.0)
    assert p.p0v.value() == -170.0 and p.p0.value() == -170
    assert p.p1v.value() == 720.0                  # clamped by the spin range
    assert len(got) == 1
    assert got[0][0] == [{"op": "hilbert"},
                         {"op": "phase", "p0": -170.0, "p1": 720.0}]
    assert not p._live_timer.isActive()
    QTest.qWait(200)
    qapp.processEvents()
    assert len(got) == 1                           # the debounce did not fire too

    p.set_phase(10.0, 0.0, apply=False)
    assert len(got) == 1 and p.p0v.value() == 10.0
    assert not p._live_timer.isActive()

    p.sync_phase_from([{"op": "hilbert"},
                       {"op": "phase", "p0": 33.0, "p1": -40.0, "pivot_frac": 0.5}])
    assert (p.p0v.value(), p.p1v.value()) == (33.0, -40.0)
    p.sync_phase_from(None)
    assert p.phase_values() == (0.0, 0.0)
    p.sync_phase_from([{"op": "baseline", "order": 3}])
    assert p.phase_values() == (0.0, 0.0)
    assert len(got) == 1                           # silent throughout

    modes = []
    p.phase_drag_mode.connect(modes.append)
    p.btnDrag.setChecked(True)
    assert modes == [True]
    assert p.btnDrag.isCheckable()

    p.p0v.setValue(100.0)
    p._nudge_p0(90.0)
    assert p.p0v.value() == -170.0                 # refactor onto wrap_p0 kept the rule
    assert p.phase_values() == (p.p0v.value(), p.p1v.value())

    # a fractional phase survives an integer crossing (the slider's echo used
    # to clobber 35.96 to 35.0) while the slider still drives the spin box
    p.set_phase(35.96, 0.25, apply=False)
    assert p.phase_values() == (35.96, 0.25)
    assert (p.p0.value(), p.p1.value()) == (35, 0)
    p.p0.setValue(40)
    assert p.p0v.value() == 40.0
    p._live_timer.stop()


# ================================================================ MainWindow
def test_drag_is_one_undo_step_equals_typing_and_resyncs_the_panel(
        qapp, win, tmp_path):
    x, before = _load(win, qapp, tmp_path)
    pp = win.proc_panel
    n0 = len(win.undo_stack)
    assert not pp.chkHilbert.isChecked()

    pp.btnDrag.setChecked(True)
    assert pp.chkHilbert.isChecked()               # pdata: Hilbert first ticked
    assert "Hilbert" in win.statusBar().currentMessage()
    assert not win.proc_dock.isHidden()
    assert win.view._pivot is not None
    assert win.view.getPlotItem().getViewBox().phase_drag_handler is not None

    for d in (20.0, 40.0, 60.0):                   # three moves of one gesture
        win.view.phase_dragged.emit(d, 0.0)
    win.view.phase_drag_released.emit()
    frac = win.view.phase_pivot_frac()
    assert pp.phase_values() == (60.0, 0.0)
    assert np.allclose(win.exp_amp, _expected(win, x, before, 60.0, 0.0, frac),
                       rtol=1e-9, atol=1e-9)       # drag == typing
    assert len(win.undo_stack) == n0 + 1           # ONE snapshot for three moves
    assert win.recipe["processing"] == [
        {"op": "hilbert"},
        {"op": "phase", "p0": 60.0, "p1": 0.0, "pivot_frac": frac}]
    assert "p0 +60.00°" in win.statusBar().currentMessage()

    win.undo()
    assert np.array_equal(win.exp_amp, before)
    assert pp.p0v.value() == 0.0                   # the controls follow the recipe
    assert not [o for o in win.recipe.get("processing") or []
                if o["op"] == "phase"]

    win.redo()
    assert np.allclose(win.exp_amp, _expected(win, x, before, 60.0, 0.0, frac))
    assert pp.p0v.value() == 60.0

    # a second gesture is a second undo step
    win.view.phase_dragged.emit(0.0, 30.0)
    win.view.phase_drag_released.emit()
    assert pp.phase_values() == (60.0, 30.0)
    assert len(win.undo_stack) == n0 + 2

    pp.btnDrag.setChecked(False)
    assert win.view.getPlotItem().getViewBox().phase_drag_handler is None
    assert win.view._pivot is not None             # the panel is still open
    assert "phase kept" in win.statusBar().currentMessage()


def test_live_off_moves_numbers_only_and_applies_once_on_release(
        qapp, win, tmp_path):
    x, before = _load(win, qapp, tmp_path)
    pp = win.proc_panel
    pp.chkLive.setChecked(False)
    pp.btnDrag.setChecked(True)
    n0 = len(win.undo_stack)

    win.view.phase_dragged.emit(30.0, 0.0)
    assert np.array_equal(win.exp_amp, before)     # numbers only
    assert pp.p0v.value() == 30.0
    assert len(win.undo_stack) == n0 + 1

    win.view.phase_drag_released.emit()            # the one deferred apply
    frac = win.view.phase_pivot_frac()
    assert np.allclose(win.exp_amp, _expected(win, x, before, 30.0, 0.0, frac))
    assert len(win.undo_stack) == n0 + 1
    assert win.recipe["processing"][-1]["p0"] == 30.0


def test_pivot_move_while_phasing_keeps_the_spectrum_and_reexpresses_p0(
        qapp, win, tmp_path):
    _load(win, qapp, tmp_path)
    pp = win.proc_panel
    pp.btnDrag.setChecked(True)
    pp.set_phase(0.0, 180.0)                       # a first-order twist to expose
    amp0 = win.exp_amp.copy()
    n0 = len(win.undo_stack)
    f_old = win.view.phase_pivot_frac()

    win.view._pivot.setValue(20.0)
    win.view._pivot.sigPositionChangeFinished.emit(win.view._pivot)
    f_new = win.view.phase_pivot_frac()
    assert f_new != f_old
    assert pp.p0v.value() == pytest.approx(wrap_p0(180.0 * (f_new - f_old)),
                                           abs=0.011)   # 2-decimal spin box
    assert pp.p1v.value() == 180.0
    assert np.allclose(win.exp_amp, amp0, atol=3e-4 * amp0.max())   # unchanged
    step = win.recipe["processing"][-1]
    assert step["op"] == "phase" and step["pivot_frac"] == f_new
    assert step["p0"] == pp.p0v.value()
    assert len(win.undo_stack) == n0                # no undo entry for a pivot move

    # with the mode off the pivot move changes nothing in the panel
    pp.btnDrag.setChecked(False)
    p0_kept = pp.p0v.value()
    amp1 = win.exp_amp.copy()
    win.view._pivot.setValue(-10.0)
    win.view._pivot.sigPositionChangeFinished.emit(win.view._pivot)
    assert pp.p0v.value() == p0_kept
    assert np.array_equal(win.exp_amp, amp1)


def test_mode_excludes_pickers_menu_toggles_and_esc_exits(qapp, win, tmp_path):
    _load(win, qapp, tmp_path)
    pp = win.proc_panel
    vb = win.view.getPlotItem().getViewBox()
    name = next(iter(win._model_actions))
    win._set_add_mode(name)
    pp.btnTpPick.setChecked(True)

    pp.btnDrag.setChecked(True)
    assert not pp.btnTpPick.isChecked()            # the picker is released
    assert win.view._add_mode is None              # so is add mode
    assert vb.phase_drag_handler is not None

    pp.btnBlPick.setChecked(True)                  # a picker releases the mode
    assert not pp.btnDrag.isChecked()
    assert vb.phase_drag_handler is None
    assert win.view._pivot is not None             # the dock is not hidden
    pp.btnBlPick.setChecked(False)

    pp.btnDrag.setChecked(True)
    win.keyPressEvent(QKeyEvent(QEvent.KeyPress, Qt.Key_Escape, Qt.NoModifier))
    assert not pp.btnDrag.isChecked()              # Esc leaves the mode

    win.start_phase_drag()                         # the menu entry toggles
    assert pp.btnDrag.isChecked()
    win.start_phase_drag()
    assert not pp.btnDrag.isChecked()

    act = _process_action(win, "&Drag to phase")
    assert act.shortcut() == QKeySequence("Ctrl+P")
    assert not act.isCheckable()

    # on the 2D view the gesture is inert
    pp.btnDrag.setChecked(True)
    n0, p0 = len(win.undo_stack), pp.p0v.value()
    win.central_stack.setCurrentWidget(win.view2d)
    win.view.phase_dragged.emit(10.0, 0.0)
    win.view.phase_drag_released.emit()
    assert len(win.undo_stack) == n0 and pp.p0v.value() == p0
    win.central_stack.setCurrentWidget(win.view)


def test_fid_view_in_the_workbench_gates_the_gesture(qapp, win, tmp_path):
    """Ctrl+T while phasing: the FID is shown, drags pan and emits are
    ignored; back on the spectrum the gesture works again."""
    _load(win, qapp, tmp_path)
    pp = win.proc_panel
    vb = win.view.getPlotItem().getViewBox()
    pp.btnDrag.setChecked(True)
    assert vb.phase_drag_handler is not None

    win._toggle_time_domain(True)
    qapp.processEvents()
    assert win.view.domain == "time"
    assert vb.phase_drag_handler is None
    n0, p0 = len(win.undo_stack), pp.p0v.value()
    win.view.phase_dragged.emit(10.0, 0.0)
    win.view.phase_drag_released.emit()
    assert len(win.undo_stack) == n0 and pp.p0v.value() == p0
    assert pp.btnDrag.isChecked()                  # the mode is sticky

    win._toggle_time_domain(False)
    qapp.processEvents()
    assert win.view.domain == "freq"
    assert vb.phase_drag_handler is not None
    win.view.phase_dragged.emit(10.0, 0.0)
    win.view.phase_drag_released.emit()
    assert pp.p0v.value() == p0 + 10.0
    assert len(win.undo_stack) == n0 + 1


# ===================================================================== docs
def test_manuals_describe_drag_to_phase():
    from larmor.desktop.help_dialog import help_path

    guide = help_path("spectra-1d").read_text(encoding="utf-8")
    assert "Drag to phase" in guide and "Ctrl+P" in guide
    ref = help_path("processing-reference").read_text(encoding="utf-8")
    assert "Drag to phase" in ref and "Ctrl+P" in ref
    assert "adjustable step size" not in ref
