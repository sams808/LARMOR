"""The spectrum view never moves because of a model, simulation or fit.

Sam's report (2026-09): "sometimes when importing a spectrum the plot window
drifts out, same sometimes when fitting" -- the x axis reached 2500, 5000,
7000 ppm on data spanning 640 ... -580 ppm. Two mechanisms, both pyqtgraph
auto-range (enabled by zoom_full at load until the user zooms):

* the model is simulated on the Czjzek KERNEL axis (engine.make_context: at
  least 150 kHz wide, 1.25x the data span otherwise), and set_model drew it
  in full, so every simulation stretched the auto range to the kernel axis;
* the fit-animation label was pinned to the current view corner each frame
  and counted in the auto-range bounds, so each frame grew the range by the
  label's extent plus padding -- a feedback loop that ran away over a fit.

Now every model-side item is added with ignoreBounds=True and drawn only
across the experiment (+ a 2 % margin); the view range is the data's at load
and the user's afterwards.
"""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtWidgets import QApplication  # noqa: E402

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


def _settle(qapp):
    """Let pyqtgraph's auto-range finish its own first pass (it refines the
    pixel padding on the next event cycle) before taking a reference."""
    for _ in range(3):
        qapp.processEvents()


def _xr(view):
    (x0, x1), _ = view.getPlotItem().getViewBox().viewRange()
    return (round(min(x0, x1), 3), round(max(x0, x1), 3))


def _yr(view):
    _, (y0, y1) = view.getPlotItem().getViewBox().viewRange()
    return (round(y0, 6), round(y1, 6))


def _wide_model(x_exp):
    """A model on an axis 40 % wider than the data on each side (a kernel
    axis is at least that wide for a narrow 27Al spectrum)."""
    lo, hi = float(np.min(x_exp)), float(np.max(x_exp))
    span = hi - lo
    xk = np.linspace(lo - 0.4 * span, hi + 0.4 * span, 2718)
    yk = np.exp(-((xk - 60.0) / 30.0) ** 2)
    return xk, yk


def test_set_model_on_a_wide_axis_leaves_the_view_alone(qapp, view):
    x = np.linspace(-580.0, 640.0, 32768)[::-1]
    y = np.exp(-((x - 60.0) / 30.0) ** 2)
    view.set_experiment(x, y)
    view.getPlotItem().enableAutoRange()
    _settle(qapp)
    before_x, before_y = _xr(view), _yr(view)
    assert before_x[0] < -580.0 < 640.0 < before_x[1]      # the data, padded

    xk, yk = _wide_model(x)
    for _ in range(5):
        view.set_model(xk, yk, [yk * 0.5, yk * 0.4], ["a", "b"], set(), x, y)
        qapp.processEvents()
        assert _xr(view) == before_x
        assert _yr(view) == before_y
    # the model is drawn only across the data (+ a small margin)
    mx = view._model.xData
    m = view.MODEL_MARGIN_FRAC * 1220.0
    assert mx.min() >= -580.0 - m - 1e-9 and mx.max() <= 640.0 + m + 1e-9
    assert len(mx) < len(xk)
    for comp in view._components:
        assert comp.xData.min() >= -580.0 - m - 1e-9
        assert comp.xData.max() <= 640.0 + m + 1e-9
    # the residual is on the experiment axis
    assert view._resid.xData.min() >= -580.0 and view._resid.xData.max() <= 640.0


def test_fit_animation_frames_do_not_drift_the_view(qapp, view):
    x = np.linspace(-580.0, 640.0, 4096)[::-1]
    y = np.exp(-((x - 60.0) / 30.0) ** 2)
    view.set_experiment(x, y)
    view.getPlotItem().enableAutoRange()
    _settle(qapp)
    before = _xr(view), _yr(view)
    xk, yk = _wide_model(x)
    view.start_fit_animation()
    for i in range(12):                        # the runaway grew every frame
        view.set_fit_frame(xk, yk * (1 + 0.05 * i), i, rms=0.1 / (i + 1))
        qapp.processEvents()
        assert (_xr(view), _yr(view)) == before
    view.stop_fit_animation()
    qapp.processEvents()
    assert (_xr(view), _yr(view)) == before


def test_paddles_and_markers_outside_the_data_do_not_pull_the_range(qapp, view):
    x = np.linspace(-100.0, 200.0, 2048)[::-1]
    y = np.exp(-((x - 60.0) / 30.0) ** 2)
    view.set_experiment(x, y)
    view.getPlotItem().enableAutoRange()
    _settle(qapp)
    before = _xr(view), _yr(view)
    # a linked sideband copy two rotor periods away sits far outside the data
    view.set_paddles([(0, 60.0, 1.0, 30.0, True), (1, 60.0 + 400.0, 0.3, 30.0, False)])
    view.set_markers([(0, 60.0, True), (1, 460.0, False)])
    qapp.processEvents()
    assert (_xr(view), _yr(view)) == before


def test_user_zoom_and_fit_result_keep_the_window(qapp, win):
    """MainWindow level: load the bundled 27Al example, zoom to a window, add
    a line, fit it (the real engine, gauss_lor) and hand the result to
    _fit_done -- the view range is untouched by the simulation and the fit."""
    from larmor import fit as fitmod
    from larmor.recipe import Recipe

    win.load_source(EXAMPLE_27AL, keep_fit=False)
    _settle(qapp)
    assert win.exp_ppm.size == 2048
    loaded = _xr(win.view), _yr(win.view)

    # a wide model landing while the load range is still the auto range
    xk, yk = _wide_model(win.exp_ppm)
    win._first_sim = False                    # not the one-time y autoscale
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
    assert (_xr(win.view), _yr(win.view)) == loaded

    # the user zooms in: that window is theirs until they change it
    win.view.setXRange(120.0, -20.0, padding=0)
    win.view.setYRange(-0.2 * float(win.exp_amp.max()),
                       1.1 * float(win.exp_amp.max()), padding=0)
    qapp.processEvents()
    zoomed = _xr(win.view), _yr(win.view)
    assert zoomed != loaded

    result = fitmod.fit(Recipe.from_dict(win.recipe), win.exp_ppm, win.exp_amp,
                        window_ppm=(120.0, -20.0))
    win._fit_done(result)
    qapp.processEvents()
    assert (_xr(win.view), _yr(win.view)) == zoomed
    # and a later simulation on the kernel-wide axis
    win._sim_done(xk, yk, [yk])
    qapp.processEvents()
    assert (_xr(win.view), _yr(win.view)) == zoomed
    # the drawn model stays within the data
    assert win.view._model.xData.max() <= float(win.exp_ppm.max()) + 0.02 * 767.0 + 1e-6
