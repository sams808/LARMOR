"""View > Y axis: raw intensity / normalise to maximum / to area / to the
area of a region.

Sam (2026-09): "Make it possible to switch the plotting area to normalized
Y: to the max, or to the area (total, or a user-selected X region)."

ONE display factor, computed from the active spectrum (larmor.display
.y_factor), multiplies everything drawn for it -- experiment, model,
components, residual, animation frames, paddles -- and the inverse maps the
interactive handles back to raw amplitudes, so the recipe, the fit and every
export stay raw. Overlays are normalised by their OWN maximum / area under
the same mode, on top of their per-overlay scale."""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtCore import QPointF, Qt  # noqa: E402
from PySide6.QtWidgets import QApplication, QDialog  # noqa: E402

from larmor import display  # noqa: E402
from larmor.io import spectra  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def view(qapp):
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    v.resize(900, 600)
    v.show()
    qapp.processEvents()
    yield v
    v.close()


@pytest.fixture
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _data(n=2048):
    x = np.linspace(200.0, -100.0, n)
    y = 5.0e6 * np.exp(-((x - 60.0) / 15.0) ** 2) + 8.0e5 * np.exp(-((x + 20.0) / 8.0) ** 2)
    return x, y


def _model(x):
    lo, hi = float(np.min(x)), float(np.max(x))
    xk = np.linspace(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo), 3001)
    return xk, 4.6e6 * np.exp(-((xk - 60.0) / 15.0) ** 2)


def _label(view):
    return view.getPlotItem().getAxis("left").labelText


# ------------------------------------------------------------- bare view
def test_max_mode_scales_every_drawn_item_by_one_factor(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    xk, yk = _model(x)
    view.set_model(xk, yk, [yk * 0.6, yk * 0.4], ["a", "b"], set(), x, y)
    view.set_paddles([(0, 60.0, 4.0e6, 30.0, True)])
    got = []
    view.paddle_moved.connect(lambda *a: got.append(a))
    assert _label(view) == "intensity"

    view.set_y_mode("max")
    f = 1.0 / float(y.max())
    assert view.y_scale() == pytest.approx(f)
    assert view.y_mode() == ("max", None)
    assert float(view._exp.yData.max()) == pytest.approx(1.0)
    assert np.allclose(view._exp.yData, y * f)
    # model and components: the same factor (masked to the data)
    mask = view._model_mask(xk)
    assert np.allclose(view._model.yData, yk[mask] * f)
    assert np.allclose(view._components[0].yData, 0.6 * yk[mask] * f)
    assert np.allclose(view._components[1].yData, 0.4 * yk[mask] * f)
    # residual strip and its zero line
    yi = np.interp(x, xk[mask], yk[mask])
    offset = -0.10 * float(y.max())
    assert np.allclose(view._resid.yData, ((y - yi) + offset) * f)
    assert view._resid_zero.value() == pytest.approx(offset * f)
    # the paddle sits at the display amplitude and reports raw
    pad = view._paddles[0]
    assert pad._amp == pytest.approx(4.0e6 * f)
    pad._handle_dragged("top", QPointF(62.0, 0.5))
    assert got[-1][0] == 0 and got[-1][1] == pytest.approx(62.0)
    assert got[-1][2] == pytest.approx(0.5 / f)              # raw amplitude
    # animation frames follow too
    view.start_fit_animation()
    view.set_fit_frame(xk, yk, 3, rms=0.1)
    assert float(view._anim_main.yData.max()) == pytest.approx(float(yk.max()) * f)
    view.stop_fit_animation()
    assert _label(view) == "intensity (normalised to max)"

    # back to raw: everything as stored
    view.set_y_mode("raw")
    assert view.y_scale() == 1.0
    assert np.array_equal(view._exp.yData, y)
    assert np.allclose(view._model.yData, yk[mask])
    assert view._paddles[0]._amp == pytest.approx(4.0e6)
    assert _label(view) == "intensity"


def test_area_and_region_modes(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    view.set_y_mode("area")
    assert display.area(x, view._exp.yData) == pytest.approx(1.0)
    assert _label(view) == "intensity (normalised to area)"
    view.set_y_mode("region", (-100.0, 0.0))                  # any order
    assert view.y_mode() == ("region", (0.0, -100.0))
    assert display.area(x, view._exp.yData, region=(0.0, -100.0)) == pytest.approx(1.0)
    assert display.area(x, view._exp.yData) > 1.0
    assert _label(view) == "intensity (normalised to area 0…−100 ppm)"
    # a region mode without a region is raw
    view.set_y_mode("region", None)
    assert view.y_mode() == ("raw", None) and view.y_scale() == 1.0
    # a new spectrum recomputes the factor for the current mode
    view.set_y_mode("max")
    x2 = np.linspace(50.0, -50.0, 512)
    y2 = 3.0 * np.exp(-(x2 / 5.0) ** 2)
    view.set_experiment(x2, y2)
    assert float(view._exp.yData.max()) == pytest.approx(1.0)
    assert view.y_scale() == pytest.approx(1.0 / float(y2.max()))


def test_zoom_follows_and_click_to_add_maps_back_to_raw(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    view.zoom_full()
    view.setYRange(0.0, 2.0e6, padding=0)
    view.set_y_mode("max")
    f = 1.0 / float(y.max())
    _, (y0, y1) = view.getPlotItem().getViewBox().viewRange()
    assert (y0, y1) == pytest.approx((0.0, 2.0e6 * f))
    got = []
    view.add_requested.connect(lambda px, amp: got.append((px, amp)))
    view.set_add_mode("gauss_lor")
    vb = view.getPlotItem().getViewBox()
    scene = vb.mapViewToScene(QPointF(60.0, 0.2))             # inside the 0...0.4 view

    class _Click:
        def button(self):
            return Qt.LeftButton

        def scenePos(self):
            return scene

        def accept(self):
            pass

    view._on_click(_Click())
    assert got and got[0][0] == pytest.approx(60.0, abs=0.5)
    assert got[0][1] == pytest.approx(0.2 / f, rel=1e-3)      # raw amplitude
    # manual baseline anchors: placed at display y, read back raw
    view.set_baseline_mode(True)
    view._add_baseline_anchor(100.0, 0.2)
    view._add_baseline_anchor(-50.0, 0.4)
    pts = view.baseline_anchors()
    assert pts[0][1] == pytest.approx(0.4 / f) and pts[1][1] == pytest.approx(0.2 / f)
    base = view.baseline_curve(x)
    assert base.max() == pytest.approx(0.4 / f) and base.min() == pytest.approx(0.2 / f)
    assert float(view._bl_curve.yData.max()) == pytest.approx(0.4)   # drawn scaled
    view.set_y_mode("raw")                                      # anchors follow
    pts = view.baseline_anchors()
    assert pts[0][1] == pytest.approx(0.4 / f) and pts[1][1] == pytest.approx(0.2 / f)


def test_a_new_active_spectrum_carries_the_display_unit_handles(qapp, view):
    """Group E (2026-09-24 review, severity high): set_experiment() assigned
    self._y_scale directly instead of going through the rescale
    _apply_y_scale performs, so a baseline anchor placed before a load stayed
    put on screen while the factor under it changed -- and
    baseline_anchors(), which divides the anchor's display y by the CURRENT
    factor, then returned a RAW level the user never clicked.
    apply_manual_baseline() subtracts exactly that from the data."""
    x = np.linspace(200.0, -100.0, 1024)
    wings = 2.0 * np.ones_like(x)
    yA = wings + 100.0 * np.exp(-((x - 60.0) / 5.0) ** 2)
    view.set_experiment(x, yA)
    view.set_y_mode("max")
    fA = view.y_scale()
    assert fA == pytest.approx(1.0 / float(yA.max()))
    view.set_paddles([(0, 60.0, 100.0, 5.0, True)])
    view.set_baseline_mode(True)
    view._add_baseline_anchor(-80.0, 2.0 * fA)          # on the drawn wings
    view._add_baseline_anchor(80.0, 2.0 * fA)
    assert [p[1] for p in view.baseline_anchors()] == pytest.approx([2.0, 2.0])

    # a new active spectrum: the SAME wings, twice the peak -- so the factor
    # halves while the raw level under the anchors is still 2.0
    yB = wings + 200.0 * np.exp(-((x - 60.0) / 5.0) ** 2)
    view.set_experiment(x, yB)
    fB = view.y_scale()
    assert fB == pytest.approx(1.0 / float(yB.max()))
    assert fB < 0.6 * fA
    # the raw level read back out of the handles is the one that was placed
    assert [p[1] for p in view.baseline_anchors()] == pytest.approx([2.0, 2.0])
    assert float(view.baseline_curve(x).max()) == pytest.approx(2.0)
    # ... and on screen they still sit on the drawn wings
    assert [float(t.pos().y()) for t in view._bl_anchors] == \
        pytest.approx([2.0 * fB, 2.0 * fB])
    assert float(view._bl_curve.yData.max()) == pytest.approx(2.0 * fB)
    # the paddle is a display-unit handle too: it follows and still reports raw
    got = []
    view.paddle_moved.connect(lambda *a: got.append(a))
    assert view._paddles[0]._amp == pytest.approx(100.0 * fB)
    view._paddles[0]._handle_dragged("top", QPointF(60.0, 50.0 * fB))
    assert got[-1][2] == pytest.approx(50.0)                 # raw amplitude
    view.set_baseline_mode(False)


# ------------------------------------------------------------- MainWindow
def _active(win):
    ppm = np.linspace(80.0, -20.0, 600)
    amp = 100.0 * np.exp(-((ppm - 30.0) / 6.0) ** 2)
    win._display_1d(ppm, amp, "27Al", 130.3, 20000.0, "active", "active")
    return ppm, amp


def _site(amp_value):
    return {"model": "gauss_lor", "label": "A", "params": {
        "isotropic_chemical_shift_ppm": {"value": 30.0, "stderr": None, "vary": True,
                                         "min": None, "max": None, "expr": None},
        "shift_fwhm_ppm": {"value": 10.0, "stderr": None, "vary": True,
                           "min": 0.1, "max": None, "expr": None},
        "gl": {"value": 0.5, "stderr": None, "vary": False, "min": 0.0,
               "max": 1.0, "expr": None},
        "amplitude": {"value": amp_value, "stderr": None, "vary": True,
                      "min": 0.0, "max": None, "expr": None}}}


def test_menu_normalises_the_plot_overlays_by_their_own_max_and_keeps_raw(win, tmp_path):
    ppm, amp = _active(win)
    for name, scale in (("ref1.csv", 5.0), ("ref2.csv", 900.0)):
        p = spectra.write_csv(tmp_path / name, ppm, scale * np.exp(-((ppm - 10.0) / 4.0) ** 2),
                              {"nucleus": "27Al"})
        assert win.add_overlay_path(str(p))
    win.recipe["sites"] = [_site(80.0)]
    win._update_paddles()
    xk = np.linspace(-40.0, 100.0, 1401)
    yk = 80.0 * np.exp(-((xk - 30.0) / 6.0) ** 2)
    win._first_sim = False
    win._sim_done(xk, yk, [yk])
    assert win._y_axis_actions["raw"].isChecked()

    win._y_axis_actions["max"].trigger()
    f = 1.0 / float(amp.max())
    assert win.view.y_mode() == ("max", None) and win.view.y_scale() == pytest.approx(f)
    assert win._y_axis_actions["max"].isChecked()
    assert float(win.view._exp.yData.max()) == pytest.approx(1.0)
    assert float(win.view._model.yData.max()) == pytest.approx(80.0 * f)
    assert "display only" in win.statusBar().currentMessage()
    # each overlay reaches 1.0 by its OWN maximum, whatever its intensity ...
    for it in win.view._overlay_items:
        assert float(it.yData.max()) == pytest.approx(1.0)
    # ... times its own scale, plus offsets as fractions of the displayed span
    win.overlay_set_scale(1, 2.0)
    assert float(win.view._overlay_items[1].yData.max()) == pytest.approx(2.0)
    win.overlay_set_yoff(0, 0.5)
    span_disp = float((amp.max() - amp.min()) * f)
    assert float(win.view._overlay_items[0].yData.max()) == pytest.approx(1.0 + 0.5 * span_disp)
    assert "×2" in win.datasets_panel._row_widgets[1]["lab"].text()
    # the paddle drag writes the RAW amplitude into the recipe
    pad = win.view._paddles[0]
    assert pad._amp == pytest.approx(80.0 * f)
    pad._handle_dragged("top", QPointF(31.0, 0.5))
    prm = win.recipe["sites"][0]["params"]
    assert prm["amplitude"]["value"] == pytest.approx(0.5 / f)
    assert prm["isotropic_chemical_shift_ppm"]["value"] == pytest.approx(31.0)
    pad._drag_finished()
    # the data and the stored overlays are untouched
    assert np.array_equal(win.exp_amp, amp)
    assert float(win._overlays[1]["amp"].max()) == pytest.approx(900.0, rel=1e-3)
    # autoscale Y works in display units
    win.autoscale_y()
    _, (y0, y1) = win.view.getPlotItem().getViewBox().viewRange()
    assert 1.0 < y1 < 1.5 and y0 < 0.0
    assert win.view.getPlotItem().getAxis("left").labelText == "intensity (normalised to max)"

    win._y_axis_actions["raw"].trigger()
    assert win.view.y_scale() == 1.0
    assert np.array_equal(win.view._exp.yData, amp)
    assert float(win.view._overlay_items[1].yData.max()) == pytest.approx(2.0 * 900.0, rel=1e-3)


def test_region_mode_asks_for_a_range_and_a_cancel_keeps_the_mode(win, monkeypatch):
    from larmor.desktop.ynorm_dialog import RegionDialog

    ppm, amp = _active(win)
    win.view.set_zones([[60.0, 0.0]])
    seen = {}

    def fake_exec(self):
        seen["zones_enabled"] = self.btnZones.isEnabled()
        self.use_zones()
        seen["from_zones"] = self.region()
        self.set_region((50.0, 10.0))
        return QDialog.Accepted

    monkeypatch.setattr(RegionDialog, "exec", fake_exec)
    win._y_axis_actions["region"].trigger()
    assert seen == {"zones_enabled": True, "from_zones": (60.0, 0.0)}
    assert win.view.y_mode() == ("region", (50.0, 10.0))
    assert display.area(ppm, win.view._exp.yData, region=(50.0, 10.0)) == pytest.approx(1.0)
    assert win._y_axis_actions["region"].isChecked()
    assert "50…10 ppm" in win.view.getPlotItem().getAxis("left").labelText

    monkeypatch.setattr(RegionDialog, "exec", lambda self: QDialog.Rejected)
    win._y_axis_actions["max"].trigger()
    win._y_axis_actions["region"].trigger()                   # cancelled
    assert win.view.y_mode() == ("max", None)
    assert win._y_axis_actions["max"].isChecked()
    assert not win._y_axis_actions["region"].isChecked()


def test_region_dialog_shortcuts_and_validation(qapp):
    from larmor.desktop.ynorm_dialog import RegionDialog

    dlg = RegionDialog(None, None, zones=[[30.0, 10.0], [80.0, 50.0]],
                       view_range=(120.0, -20.0))
    assert dlg.region() == (120.0, -20.0)                    # the view by default
    dlg.use_zones()
    assert dlg.region() == (80.0, 10.0)                      # the zones' outer bounds
    dlg.hi.setValue(5.0)
    dlg.lo.setValue(5.0)
    dlg.accept()
    assert dlg.result() != QDialog.Accepted                  # an empty region stays open
    dlg.lo.setValue(15.0)
    assert dlg.region() == (15.0, 5.0)                       # any order
    dlg.accept()
    assert dlg.result() == QDialog.Accepted
    bare = RegionDialog(None, (100.0, -20.0))
    assert not bare.btnZones.isEnabled() and not bare.btnView.isEnabled()
    assert bare.region() == (100.0, -20.0)


def test_mode_is_remembered_in_qsettings(win, monkeypatch):
    from larmor.desktop import mw_chrome, mw_menus

    class FakeSettings:
        store = {}

        def __init__(self, *args):
            pass

        def setValue(self, key, value):
            FakeSettings.store[key] = value

        def value(self, key, default=None, type=None):
            return FakeSettings.store.get(key, default)

    monkeypatch.setattr(mw_chrome, "QSettings", FakeSettings)
    monkeypatch.setattr(mw_menus, "QSettings", FakeSettings)
    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    _active(win)
    win._set_y_mode("max")
    assert FakeSettings.store["yAxisMode"] == "max"
    assert win._saved_y_mode() == ("max", None)
    win.view.set_y_mode("region", (100.0, -20.0))
    win._set_y_mode("area")
    assert FakeSettings.store == {"yAxisMode": "area", "yAxisRegion": ""}
    FakeSettings.store.update({"yAxisMode": "region", "yAxisRegion": "100,-20"})
    assert win._saved_y_mode() == ("region", (100.0, -20.0))
    FakeSettings.store.update({"yAxisMode": "region", "yAxisRegion": ""})
    assert win._saved_y_mode() == ("raw", None)
    FakeSettings.store.update({"yAxisMode": "bogus"})
    assert win._saved_y_mode() == ("raw", None)
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    FakeSettings.store.update({"yAxisMode": "max", "yAxisRegion": ""})
    assert win._saved_y_mode() == ("raw", None)               # tests never inherit a choice
