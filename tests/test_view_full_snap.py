"""The Full view and the snap-back.

Sam (2026-09): "When clicking View all on the plot window, make sure it
resets to the good size to view the full spectrum. Also if the user goes too
far and exceeds the window limits, go back to the full view."

* zoom_full (View > Zoom > Full spectrum, the sidebar's Full) sets the x
  range to the data extent + FULL_MARGIN_FRAC of its span and the y range to
  the data's min...max with room for the residual strip -- from the arrays
  (larmor.display.full_extents), never pyqtgraph's auto-range, which the
  model items would widen; the same on the FID's ms axis.
* A REQUESTED x range that no longer intersects the data, or is wider than
  SNAP_WIDTH_FACTOR data spans, is replaced by the Full view: every mouse
  pan / wheel zoom and every programmatic setXRange funnels through
  AnchoredViewBox.setRange, where SpectrumView._snap_guard inspects the
  request before pyqtgraph clamps it to the limits. A zoom inside the data
  is kept.
* ViewBox.setLimits keeps the view within +-LIMIT_X_SPANS data spans in x.
"""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtCore import QPointF  # noqa: E402
from PySide6.QtWidgets import QApplication, QMenu, QToolBar  # noqa: E402

from larmor import display  # noqa: E402

EXAMPLE_27AL = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "examples", "pCABS2-4", "3616")


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


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _xr(view):
    (x0, x1), _ = view.getPlotItem().getViewBox().viewRange()
    return min(x0, x1), max(x0, x1)


def _yr(view):
    _, (y0, y1) = view.getPlotItem().getViewBox().viewRange()
    return min(y0, y1), max(y0, y1)


def _data(n=4096):
    x = np.linspace(640.0, -580.0, n)                  # high -> low ppm
    y = 2.0e6 * np.exp(-((x - 60.0) / 30.0) ** 2) - 3.0e4    # a negative floor
    return x, y


def _wide_model(x):
    lo, hi = float(np.min(x)), float(np.max(x))
    span = hi - lo
    xk = np.linspace(lo - 0.4 * span, hi + 0.4 * span, 2718)
    return xk, 1.8e6 * np.exp(-((xk - 60.0) / 30.0) ** 2)


# ------------------------------------------------------------- the Full view
def test_zoom_full_is_the_data_extent_not_the_auto_range(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    xk, yk = _wide_model(x)
    view.set_model(xk, yk, [yk * 0.5], ["a"], set(), x, y)
    view.set_paddles([(0, 60.0, 1.8e6, 30.0, True), (1, 1400.0, 0.3e6, 30.0, False)])
    view.setXRange(20.0, 100.0, padding=0)                 # the user zoomed in
    view.zoom_full()
    qapp.processEvents()
    m = display.FULL_MARGIN_FRAC * 1220.0
    lo, hi = _xr(view)
    assert abs(lo - (-580.0)) <= m + 1e-6 and abs(hi - 640.0) <= m + 1e-6
    assert lo < -580.0 < 640.0 < hi                        # data fully inside
    (_, _), (ylo, yhi) = display.full_extents(x, y)
    assert _yr(view) == pytest.approx((ylo, yhi))
    assert yhi > float(y.max()) and ylo < -display.RESID_ROOM_FRAC * float(y.max())
    # neither the wide model nor the far paddle widened the view
    assert hi < 700.0
    vb = view.getPlotItem().getViewBox()
    assert not any(vb.autoRangeEnabled())                 # the view holds


def test_zoom_full_frames_the_visible_overlays(qapp, view):
    """Group D (2026-09-24 review): zoom_full built its range from
    _display_arrays(), the active trace alone, so a compared spectrum over a
    wider ppm range -- or lifted by the Datasets dock's stack offset -- was
    outside the view the Full command had just set. Before the batch
    zoom_full was enableAutoRange(), which framed the overlay curves."""
    x, y = _data()                                      # 640 ... -580 ppm
    view.set_experiment(x, y)
    ox = x + 900.0                                      # 1540 ... 320
    oy = y + 3.0e6                                      # stacked above
    view.set_overlays([(ox, oy, "#ff0000", "ref")])
    view.set_paddles([(0, 6000.0, 1.8e6, 30.0, False)])  # a far sideband copy
    view.zoom_full()
    qapp.processEvents()
    lo, hi = _xr(view)
    ylo, yhi = _yr(view)
    assert lo < -580.0 and hi > 1540.0                  # the overlay is inside
    assert hi < 2000.0                                  # the far paddle is not
    assert yhi > float(oy.max()) and ylo < float(y.min())
    assert ylo < -display.RESID_ROOM_FRAC * float(y.max())   # the resid strip
    assert not any(view.getPlotItem().getViewBox().autoRangeEnabled())

    # the snap-back target frames exactly the same thing
    view.setXRange(20.0, 100.0, padding=0)
    view.setXRange(-9000.0, -8000.0, padding=0)         # entirely off the data
    qapp.processEvents()
    assert _xr(view) == pytest.approx((lo, hi))
    assert _yr(view) == pytest.approx((ylo, yhi))

    # View > Overlays (hidden) takes them back out of Full
    view.set_overlays_hidden(True)
    view.zoom_full()
    qapp.processEvents()
    assert _xr(view)[1] < 700.0
    assert _yr(view) == pytest.approx(display.full_extents(x, y)[1])
    view.set_overlays_hidden(False)
    view.set_overlays([])
    view.zoom_full()
    qapp.processEvents()
    assert _xr(view)[1] < 700.0


def test_zoom_full_on_the_fid_axis_and_back(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    t = np.linspace(0.0, 50.0, 512)
    fid = 7.0e5 * np.exp(-t / 12.0)
    view.set_fid(t, fid, "FID (real)")
    view.setXRange(5.0, 6.0, padding=0)
    view.zoom_full()
    qapp.processEvents()
    assert view.domain == "time"
    assert _xr(view) == pytest.approx((-1.0, 51.0))       # 0...50 ms + 2 %
    (_, _), (ylo, yhi) = display.full_extents(t, fid)
    assert _yr(view) == pytest.approx((ylo, yhi))
    assert view.getPlotItem().getAxis("left").labelText == "intensity"
    view.set_experiment(x, y)                               # the spectrum is back
    view.zoom_full()
    lo, hi = _xr(view)
    m = display.FULL_MARGIN_FRAC * 1220.0
    assert lo == pytest.approx(-580.0 - m) and hi == pytest.approx(640.0 + m)


def test_zoom_full_without_data_falls_back_to_auto_range(qapp, view):
    view.zoom_full()                                        # no crash, no data
    assert any(view.getPlotItem().getViewBox().autoRangeEnabled())


# ------------------------------------------------------------- snap-back
def test_programmatic_ranges_outside_the_data_snap_back(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    view.zoom_full()
    full = _xr(view), _yr(view)
    for req in ((1000.0, 1200.0), (-2000.0, -1000.0), (700.0, 650.0)):
        view.setXRange(*req, padding=0)                     # entirely off the data
        qapp.processEvents()
        assert (_xr(view), _yr(view)) == full
    view.setXRange(-3000.0, 3000.0, padding=0)              # 6000 > 3 x 1220
    assert (_xr(view), _yr(view)) == full


def test_a_zoom_inside_the_data_is_kept(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    view.zoom_full()
    view.setXRange(100.0, 20.0, padding=0)
    assert _xr(view) == pytest.approx((20.0, 100.0))
    view.setYRange(-1e5, 5e5, padding=0)
    assert _yr(view) == pytest.approx((-1e5, 5e5))
    view.setXRange(600.0, 900.0, padding=0)                 # half outside: kept
    assert _xr(view) == pytest.approx((600.0, 900.0))
    # a later model or paddle changes nothing
    xk, yk = _wide_model(x)
    view.set_model(xk, yk, [yk], ["a"], set(), x, y)
    view.set_paddles([(0, 60.0, 1.8e6, 30.0, True)])
    qapp.processEvents()
    assert _xr(view) == pytest.approx((600.0, 900.0))
    assert _yr(view) == pytest.approx((-1e5, 5e5))


def test_mouse_pan_and_wheel_paths_snap_back_too(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    view.zoom_full()
    full = _xr(view)
    vb = view.getPlotItem().getViewBox()
    view.setXRange(100.0, 20.0, padding=0)
    vb.translateBy(x=2000.0)                                # dragged far off the data
    assert _xr(view) == pytest.approx(full)
    view.setXRange(100.0, 20.0, padding=0)
    vb.scaleBy((100.0, 1.0), QPointF(60.0, 0.0))            # wheel-zoomed way out
    assert _xr(view) == pytest.approx(full)
    view.setXRange(100.0, 20.0, padding=0)
    vb.scaleBy((0.5, 1.0), QPointF(60.0, 0.0))              # zoomed in: kept
    assert _xr(view) == pytest.approx((40.0, 80.0))


def test_limits_keep_the_mouse_within_one_data_span(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    vb = view.getPlotItem().getViewBox()
    assert vb.state["limits"]["xLimits"] == pytest.approx([-580.0 - 1220.0, 640.0 + 1220.0])
    ylo, yhi = vb.state["limits"]["yLimits"]
    assert ylo < float(y.min()) and yhi > float(y.max())
    # a request that still shows the data but reaches past the envelope is
    # shifted back inside it
    view.setXRange(-1900.0, -500.0, padding=0)
    lo, hi = _xr(view)
    assert lo >= -1800.0 - 1e-6 and hi > -580.0
    # the envelope follows the data: a new spectrum, new limits
    x2 = np.linspace(100.0, -100.0, 512)
    view.set_experiment(x2, np.exp(-(x2 / 10.0) ** 2))
    assert vb.state["limits"]["xLimits"] == pytest.approx([-300.0, 300.0])
    # and the FID axis: the ms envelope, then the ppm one again
    t = np.linspace(0.0, 50.0, 256)
    view.set_fid(t, np.exp(-t / 10.0), "FID (real)")
    assert vb.state["limits"]["xLimits"] == pytest.approx([-50.0, 100.0])
    view.set_experiment(x2, np.exp(-(x2 / 10.0) ** 2))
    assert vb.state["limits"]["xLimits"] == pytest.approx([-300.0, 300.0])


def test_overlays_and_far_paddles_extend_what_counts_as_data(qapp, view):
    x, y = _data()
    view.set_experiment(x, y)
    view.zoom_full()
    # an overlay far from the active spectrum: looking at it is not 'outside'
    view.set_overlays([(x + 1500.0, y, "#ff0000", "ref")])
    view.setXRange(1000.0, 1200.0, padding=0)
    assert _xr(view) == pytest.approx((1000.0, 1200.0))
    view.set_overlays([])
    view.setXRange(1000.0, 1200.0, padding=0)                # gone: off the data
    assert _xr(view)[1] < 700.0
    # a linked sideband paddle outside the data can be reached
    view.set_paddles([(1, 1100.0, 0.3e6, 30.0, False)])
    view.setXRange(1000.0, 1200.0, padding=0)
    assert _xr(view) == pytest.approx((1000.0, 1200.0))


# ------------------------------------------------------------- MainWindow
def test_full_button_and_menu_reset_to_the_data(qapp, win):
    win.load_source(EXAMPLE_27AL, keep_fit=False)
    qapp.processEvents()
    assert win.exp_ppm.size == 2048
    (xl, xh), (yl, yh) = display.full_extents(win.exp_ppm, win.exp_amp)
    win.view.setXRange(120.0, -20.0, padding=0)
    win.view.setYRange(0.0, 0.5 * float(win.exp_amp.max()), padding=0)
    sidebar = [tb for tb in win.findChildren(QToolBar) if tb.objectName() == "sidebar"][0]
    full = [a for a in sidebar.actions() if a.text() == "Full"][0]
    full.trigger()
    qapp.processEvents()
    assert _xr(win.view) == pytest.approx((xl, xh))
    assert _yr(win.view) == pytest.approx((yl, yh))
    # a wide simulation landing afterwards moves nothing
    xk, yk = _wide_model(win.exp_ppm)
    win._first_sim = False
    win.recipe["sites"] = [{"model": "gauss_lor", "label": "A", "params": {
        "isotropic_chemical_shift_ppm": {"value": 55.0, "stderr": None, "vary": True,
                                         "min": None, "max": None, "expr": None},
        "shift_fwhm_ppm": {"value": 20.0, "stderr": None, "vary": True,
                           "min": 0.1, "max": None, "expr": None},
        "gl": {"value": 0.5, "stderr": None, "vary": False, "min": 0.0,
               "max": 1.0, "expr": None},
        "amplitude": {"value": float(win.exp_amp.max()), "stderr": None,
                      "vary": True, "min": 0.0, "max": None, "expr": None}}}]
    win._sim_done(xk, yk, [yk])
    qapp.processEvents()
    assert _xr(win.view) == pytest.approx((xl, xh))
    # View > Zoom > Full spectrum does the same after a zoom
    win.view.setXRange(80.0, 40.0, padding=0)
    menus = [m for m in win.menuBar().findChildren(QMenu) if m.title() == "&Zoom"]
    full_spectrum = [a for a in menus[0].actions() if a.text() == "&Full spectrum"]
    full_spectrum[0].trigger()
    assert _xr(win.view) == pytest.approx((xl, xh))
    # panning off the spectrum comes back to it
    win.view.setXRange(xh + 100.0, xh + 300.0, padding=0)
    assert _xr(win.view) == pytest.approx((xl, xh))
