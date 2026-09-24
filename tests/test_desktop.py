"""Offscreen smoke tests for the desktop app (no window is shown)."""
import json
import os

import numpy as np
import pytest

from conftest import (BRUKER_1R, BRUKER_2RR_DQSQ, BRUKER_2RR_MQMAS,
                      BRUKER_2RR_PSEUDO, BRUKER_FID, BRUKER_SER, CAALGLASS,
                      require)

pyside = pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402


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


def test_open_any_type_never_rejects(qapp, win):
    """Every data type opens with a basic display: 1D → workbench, 2D → the
    contour view, raw fid/ser → a preview. None is rejected."""
    def is_1d():
        return win.central_stack.currentWidget() is win.view

    # 1D processed spectrum → fit workbench
    win.load_source(str(require(BRUKER_1R)), keep_fit=False)
    qapp.processEvents()
    assert is_1d() and win.exp_ppm.size > 0 and win.recipe["nucleus"] == "27Al"

    # real MQMAS 2rr → the 2D contour view (was previously a rejection)
    win.load_source(str(require(BRUKER_2RR_MQMAS)), keep_fit=False)
    qapp.processEvents()
    assert not is_1d()
    assert win.view2d.data is not None and win.view2d.data.z.ndim == 2

    # pseudo-2D relaxation 2rr → the 2D view (the exact case the user hit)
    win.load_source(str(require(BRUKER_2RR_PSEUDO)), keep_fit=False)
    qapp.processEvents()
    assert not is_1d()

    # raw fid → 1D magnitude preview on the workbench
    win.load_source(str(require(BRUKER_FID)), keep_fit=False)
    qapp.processEvents()
    assert is_1d() and win.exp_ppm.size > 0

    # raw ser → 2D preview
    win.load_source(str(require(BRUKER_SER)), keep_fit=False)
    qapp.processEvents()
    assert not is_1d()

    # a DQ/SQ 2D correlation (indirect dim = double quantum) → the 2D view
    win.load_source(str(require(BRUKER_2RR_DQSQ)), keep_fit=False)
    qapp.processEvents()
    assert not is_1d()
    assert win.view2d.data is not None and win.view2d.data.z.ndim == 2

    # pull a 1D projection out of the 2D → back to the workbench
    win.view2d._emit_projection("skyline")
    qapp.processEvents()
    assert is_1d() and win.exp_ppm.size > 0


def test_overlay_cockpit(qapp, win):
    """Overlays draw behind the active spectrum, honour visibility/stack offset,
    and 'make active' promotes an overlay to the fit target."""
    x = np.linspace(-100, 100, 200)
    win._display_1d(x, np.exp(-(x ** 2) / 200), "27Al", 100.0, None, "A", "srcA")
    qapp.processEvents()
    win._add_overlay("B", x, 0.5 * np.exp(-((x - 20) ** 2) / 200), "srcB")
    win._add_overlay("C", x, 0.3 * np.exp(-((x + 30) ** 2) / 200), "srcC")
    qapp.processEvents()
    assert len(win._overlays) == 2
    assert len(win.view._overlay_items) == 2

    win.datasets_panel.offset.setValue(0.4); qapp.processEvents()
    assert len(win.view._overlay_items) == 2      # still drawn, now stacked
    win.overlay_visibility(0, False); qapp.processEvents()
    assert len(win.view._overlay_items) == 1      # hidden one dropped
    win.overlay_remove(0); qapp.processEvents()
    assert len(win._overlays) == 1


def test_workspaces_switch_close(qapp, win):
    """Opening a 2D and extracting a trace create separate workspaces; switching
    restores their state, and closing frees them."""
    from larmor import twod

    x = np.linspace(-100, 100, 200)
    win._display_1d(x, np.exp(-(x / 20) ** 2), "27Al", 156.0, None, "A", "sA")
    qapp.processEvents()
    win._model_actions["gauss_lor"].setChecked(True)
    win.add_site_at(0.0, 1.0); qapp.processEvents(); win._sync_active()
    assert len(win.workspaces) == 1 and win.workspaces[0]["has_fit"]

    f2 = np.linspace(-50, 50, 80); f1 = np.linspace(-30, 30, 40)
    Z = (np.exp(-((f2[None, :] - 10) / 6) ** 2)
         * np.exp(-((f1[:, None] - 5) / 6) ** 2))
    win._show_2d(twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="1H",
                             larmor_MHz=500.0), "HMQC", "2D")
    qapp.processEvents()
    win.view2d.set_projection_1d("f2", f2, np.exp(-((f2 - 10) / 6) ** 2))
    assert len(win.workspaces) == 2 and win.active_ws == 1

    win.view2d._emit_projection("skyline"); qapp.processEvents()   # trace → new ws
    assert len(win.workspaces) == 3 and win.central_stack.currentWidget() is win.view

    win.back_to_2d(); qapp.processEvents()
    assert win.active_ws == 1 and win.central_stack.currentWidget() is win.view2d
    assert win.view2d._hmqc["f2"] is not None            # 2D state restored

    win.switch_workspace(0); qapp.processEvents()
    assert len(win.recipe["sites"]) == 1                 # the fit came back

    win.close_workspace(1); qapp.processEvents()
    assert len(win.workspaces) == 2


def test_hmqc_explorer_pick_and_back_nav(qapp, win, tmp_path):
    """Arming an HMQC projection pick routes the next Explorer activation to the
    overlay (colored), and Back-to-2D restores the map with overlays intact."""
    from larmor import twod
    from larmor.io import spectra

    f2 = np.linspace(-50, 50, 120); f1 = np.linspace(-30, 30, 60)
    Z = (np.exp(-((f2[None, :] - 10) / 4) ** 2)
         * np.exp(-((f1[:, None] - 5) / 4) ** 2))
    win._show_2d(twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="1H",
                             larmor_MHz=500.0), "HMQC", "2D")
    qapp.processEvents()
    csv = tmp_path / "n15.csv"
    spectra.write_csv(csv, f1, np.exp(-((f1 - 5) / 4) ** 2), {"nucleus": "15N"})

    win.load_projection_1d("f1")
    assert win._proj_pick_axis == "f1"
    win._explorer_open(str(csv))              # simulate the Explorer click
    qapp.processEvents()
    assert win._proj_pick_axis is None
    assert win.view2d._hmqc["f1"] is not None
    assert win.view2d._hmqc_color["f1"] == win.PROJ_COLOR["f1"]

    win.view2d._emit_uncorrelated("f1"); qapp.processEvents()
    assert win.central_stack.currentWidget() is win.view      # went to workbench
    win.back_to_2d(); qapp.processEvents()
    assert win.central_stack.currentWidget() is win.view2d
    assert win.view2d._hmqc["f1"] is not None                 # state preserved


def test_hmqc_uncorrelated_features(qapp):
    """Overlaying a 1D on an HMQC projection and subtracting the (scaled)
    projection isolates the features that do NOT correlate."""
    from larmor import twod
    from larmor.desktop.twod_view import Contour2DView

    f2 = np.linspace(-50, 50, 200); f1 = np.linspace(-30, 30, 80)
    Z = (np.exp(-((f2[None, :] - 10) / 3) ** 2)
         * np.exp(-((f1[:, None] - 5) / 3) ** 2))          # one cross-peak @ F2=10
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="1H", larmor_MHz=500.0)
    v = Contour2DView(); v.set_data(d.normalized(), "HMQC"); qapp.processEvents()

    # a 1D with the correlated peak (@10) AND an uncorrelated one (@-20)
    oned = np.exp(-((f2 - 10) / 3) ** 2) + 0.7 * np.exp(-((f2 + 20) / 3) ** 2)
    v.set_projection_1d("f2", f2, oned); qapp.processEvents()
    assert v.btnHmqc.isChecked()

    got = {}
    v.slice_to_fit.connect(lambda p, a, lab: got.update(ppm=p, amp=a))
    v._emit_uncorrelated("f2"); qapp.processEvents()
    ppm, amp = got["ppm"], got["amp"]
    at10 = amp[int(np.argmin(np.abs(ppm - 10)))]
    atm20 = amp[int(np.argmin(np.abs(ppm + 20)))]
    assert abs(at10) < 0.15            # correlated peak removed
    assert atm20 > 0.5                 # uncorrelated peak retained


def test_twod_fit_wiring(qapp, win):
    """A displayed 2D gets a recipe, click-to-add places 2D sites, the fitted
    overlay renders, and run_fit routes to the 2D path (rejecting 1D-only models)."""
    from larmor import twod

    f2 = np.linspace(-60, 60, 48); f1 = np.linspace(-30, 30, 24)
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1,
                    z=np.abs(np.random.RandomState(3).randn(24, 48)),
                    nucleus="27Al", larmor_MHz=195.5)
    win._show_2d(d, "syn", "2D"); qapp.processEvents()
    assert win.central_stack.currentWidget() is win.view2d
    assert win.recipe is not None and win.recipe["nucleus"] == "27Al"
    assert win._data2d_fittable

    win._model_actions["czjzek"].setChecked(True)
    win.add_site_2d(-10.0, 5.0); qapp.processEvents()
    assert len(win.recipe["sites"]) == 1 and win.recipe["sites"][0]["model"] == "czjzek"

    win.view2d.set_model(d.z, f2, f1); qapp.processEvents()
    assert win.view2d._model is not None

    # a pseudo-2D (relaxation) array is not fittable
    dp = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=d.z, nucleus="27Al")
    dp.notes = ["pseudo-2D (arrayed)"]
    win._show_2d(dp, "relax", "pseudo-2D"); qapp.processEvents()
    assert not win._data2d_fittable


def test_calibrate_and_measure(qapp):
    """Calibrate shifts the axis so a picked peak reads the target; the 2D
    measure readout reports Δ in ppm and Hz."""
    from larmor import twod
    from larmor.desktop.twod_view import Contour2DView

    # 2D axis calibration
    f2 = np.linspace(-100, 100, 64); f1 = np.linspace(-40, 40, 16)
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=np.random.RandomState(2).randn(16, 64),
                    larmor_MHz=100.0)
    v = Contour2DView(); v.set_data(d.normalized(), "syn"); qapp.processEvents()
    v._shift_axes(5.0, -3.0); qapp.processEvents()
    assert v.data.f2_ppm[0] == pytest.approx(f2[0] + 5.0)
    assert v.data.f1_ppm[0] == pytest.approx(f1[0] - 3.0)

    v.btnMeasure.setChecked(True); qapp.processEvents()
    assert len(v._mtargets) == 2
    v._mtargets[0].setPos((-50.0, 10.0))
    v._mtargets[1].setPos((-80.0, 4.0)); qapp.processEvents()
    assert "ΔF2 30.00 ppm" in v.cursor.text()
    assert "3000 Hz" in v.cursor.text()          # 30 ppm × 100 MHz


def test_twod_phasing_and_contours(qapp):
    """The 2D contour view: phase_1d matches Data2D.phased, and the pick →
    phase-traces → apply flow mutates the committed data and returns to the map."""
    from larmor import twod
    from larmor.desktop.twod_view import Contour2DView

    # phase_1d must equal the row-by-row result of Data2D.phased (so the live
    # single-row preview is faithful to what Apply will do)
    f2 = np.linspace(-100, 100, 256)
    f1 = np.linspace(-50, 50, 32)
    z = np.random.RandomState(1).randn(32, 256)
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=z)
    dp = d.phased("f2", 33.0, 110.0, pivot_ppm=10.0)
    row = twod.phase_1d(z[7], f2, 33.0, 110.0, 10.0)
    assert np.allclose(row, dp.z[7], atol=1e-9)
    assert d.phased("f2", 0.0, 0.0) is d          # no-op returns self

    v = Contour2DView()
    v.set_data(d.normalized(), "synthetic")
    qapp.processEvents()
    for s in ("both", "negative", "positive"):
        v.sign.setCurrentText(s); qapp.processEvents()
    v.btnPhase.setChecked(True)
    v._pick_axis = "f2"; v._pivot = 10.0; v._picks = [7]
    v._enter_phasing(); qapp.processEvents()
    assert v.stack.currentWidget() is v.phase_glw and len(v._pref) == 1
    v._nudge_p0(90.0); v.p1v.setValue(30.0); qapp.processEvents()
    before = v._committed.z.copy()
    v._apply_phase(); qapp.processEvents()
    assert not np.allclose(before, v._committed.z)
    assert v.stack.currentWidget() is v.glw       # returned to the contour map


@pytest.mark.slow
def test_load_fit_quantify_undo(qapp, win):
    require(CAALGLASS)
    win.load_source(str(CAALGLASS))
    qapp.processEvents()
    assert win.recipe is not None and len(win.recipe["sites"]) == 5
    assert win.exp_ppm.size == 8192

    # from_dict must never mutate the window's recipe (regression:
    # shallow-copy bug found by this very scenario)
    from larmor.recipe import Recipe

    before = json.dumps(win.recipe)
    Recipe.from_dict(win.recipe)
    assert json.dumps(win.recipe) == before

    from larmor import fit as fitmod

    result = fitmod.fit(Recipe.from_dict(win.recipe),
                        win.exp_ppm, win.exp_amp, window_ppm=(150.0, -80.0))
    win._fit_done(result)
    qapp.processEvents()
    assert result.rmsd < 0.01
    assert win.qtable.rowCount() == 5

    # undo restores the site count
    n = len(win.recipe["sites"])
    win.snapshot()
    win.recipe["sites"].pop()
    win.on_structure_changed()
    win.undo()
    assert len(win.recipe["sites"]) == n


def test_add_site_and_paddles(qapp, win):
    require(CAALGLASS)
    win.load_source(str(CAALGLASS))
    qapp.processEvents()
    n = len(win.recipe["sites"])
    win._model_actions["gauss_lor"].setChecked(True)
    win.add_site_at(42.0, 123.0)
    assert len(win.recipe["sites"]) == n + 1
    new = win.recipe["sites"][-1]
    assert new["model"] == "gauss_lor"
    assert new["params"]["isotropic_chemical_shift_ppm"]["value"] == 42.0
    assert new["params"]["amplitude"]["value"] == 123.0
    # a dmfit-style paddle exists for it
    assert any(p.index == n for p in win.view._paddles)

    # paddle drag updates position, amplitude AND width in the recipe
    win.on_paddle_moved(n, 55.5, 200.0, 8.0)
    prm = win.recipe["sites"][n]["params"]
    assert prm["isotropic_chemical_shift_ppm"]["value"] == 55.5
    assert prm["amplitude"]["value"] == 200.0
    assert prm["shift_fwhm_ppm"]["value"] == 8.0
    win.on_paddle_released(n)

    # the fit-parameters table shows one row per line
    assert win.lines_table.table.rowCount() == n + 1


def test_add_mode_does_not_wipe_the_plot_menu(qapp):
    """pyqtgraph's setMenuEnabled(False) DESTROYS the ViewBoxMenu object and
    setMenuEnabled(True) builds a brand-new default one -- set_add_mode used
    that to suppress right-click while placing lines, which silently wiped
    "Export figure..."/"Send to Plotting studio" (added once at construction)
    the first time add-mode was ever exited. Must survive every cycle."""
    from larmor.desktop.plot import SpectrumView

    view = SpectrumView()
    vb = view.getPlotItem().getViewBox()
    before = [a.text() for a in vb.menu.actions()]
    assert "Export figure…" in before and "Send to Plotting studio" in before

    view.set_add_mode("gauss_lor")
    assert vb.menu is None                       # pyqtgraph tears it down here
    view.set_add_mode(None)
    after = [a.text() for a in vb.menu.actions()]
    assert after == before                        # fully restored, no duplicates

    # a second cycle must not duplicate or drop anything either
    view.set_add_mode("czjzek")
    view.set_add_mode(None)
    assert [a.text() for a in vb.menu.actions()] == before
    view.close()


def test_cofit_dialog_plots_both_datasets(qapp):
    """The co-fit dialog draws one plot per dataset: a 1D exp-vs-fit overlay and
    a 2D experiment/model contour overlay (not just the text report)."""
    from types import SimpleNamespace

    import pyqtgraph as pg

    from larmor import twod
    from larmor.desktop.cofit_dialog import CofitDialog
    from larmor.recipe import Recipe

    x = np.linspace(-50, 120, 400)
    amp = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    base = {"kind": "1d", "label": "MAS", "ppm": x, "amp": amp,
            "nucleus": "27Al", "larmor": 195.5}
    f2 = np.linspace(-40, 110, 60); f1 = np.linspace(-30, 90, 50)
    Z = np.exp(-0.5 * (((f2[None, :] - 55) / 12) ** 2
                       + ((f1[:, None] - 60) / 10) ** 2))
    d2 = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="27Al", larmor_MHz=195.5)
    d2ds = {"kind": "2d", "label": "MQMAS", "data2d": d2,
            "nucleus": "27Al", "larmor": 195.5}

    recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 195.5, "sites": [
        {"model": "czjzek", "label": "AlIV", "params": {
            "isotropic_chemical_shift_ppm": {"value": 60.0, "vary": True,
                                             "min": 0, "max": 120},
            "sigma_Cq_MHz": {"value": 2.0, "vary": True, "min": 0.2, "max": 8},
            "shift_fwhm_ppm": {"value": 12.0, "vary": True, "min": 1, "max": 30},
            "line_fwhm_ppm": {"value": 4.0, "vary": True, "min": 0},
            "amplitude": {"value": 1.0, "vary": True, "min": 0}}}]}

    dlg = CofitDialog(None, recipe, base)
    dlg.datasets = [base, d2ds]
    dlg._result = SimpleNamespace(
        recipes=[Recipe.from_dict(recipe)], rmsd=[0.031, 0.048],
        per_dataset=[{"kind": "1d", "x": x, "y_fit": 1.05 * amp},
                     {"kind": "2d", "f2": f2, "f1": f1, "z_fit": Z * 0.98,
                      "per_site": [Z * 0.98]}])
    dlg._plot_result()
    n_plots = sum(1 for i in range(dlg._plot_v.count())
                  if isinstance(dlg._plot_v.itemAt(i).widget(), pg.PlotWidget))
    assert n_plots == 2
    dlg.close()


def test_overlay_1d_on_2d_superposition(qapp, win):
    """Overlay 1D ▸ Current 1D → F2/F1 superposes the current spectrum on the
    map's projection; Clear removes it."""
    from larmor import twod

    f2 = np.linspace(-40, 110, 60); f1 = np.linspace(-30, 90, 50)
    Z = np.exp(-0.5 * (((f2[None, :] - 55) / 12) ** 2
                       + ((f1[:, None] - 60) / 10) ** 2))
    d2 = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="27Al", larmor_MHz=195.5)
    win._show_2d(d2, "MQMAS", "2D"); qapp.processEvents()

    x = np.linspace(-50, 120, 400); amp = np.exp(-0.5 * ((x - 55) / 9) ** 2)
    win.exp_ppm, win.exp_amp = x, amp
    win.overlay_1d_on_2d("f2", "current")
    assert win.view2d._hmqc["f2"] is not None
    win.overlay_1d_on_2d("f1", "current")
    assert win.view2d._hmqc["f1"] is not None
    win.view2d.clear_projection_1d()
    assert win.view2d._hmqc["f2"] is None and win.view2d._hmqc["f1"] is None

    # no current 1D → no crash, nothing overlaid
    win.exp_ppm, win.exp_amp = np.array([]), np.array([])
    win.overlay_1d_on_2d("f2", "current")
    assert win.view2d._hmqc["f2"] is None


def test_twod_axis_orientation_defaults_and_flip(qapp):
    """The 2D map defaults to the standard NMR/dmfit convention: F2 high-ppm
    LEFT (invertX) and F1 high-ppm TOP (NOT invertY). Both are user-flippable
    (view only) and the projections follow the contour."""
    from larmor.desktop.twod_view import Contour2DView

    dv = Contour2DView()
    assert dv._flip == {"f2": True, "f1": False}
    vb = dv.p_main.getViewBox()
    assert vb.xInverted() is True and vb.yInverted() is False
    # projections share the contour's directions
    assert dv.p_top.getViewBox().xInverted() is True
    assert dv.p_left.getViewBox().yInverted() is False
    # menu checkmarks match the state
    assert dv.actFlipF2.isChecked() and not dv.actFlipF1.isChecked()

    # flipping F1 inverts the contour AND its projection together
    dv.toggle_axis_flip("f1", True)
    assert dv.p_main.getViewBox().yInverted() is True
    assert dv.p_left.getViewBox().yInverted() is True
    dv.toggle_axis_flip("f1", False)
    assert dv.p_main.getViewBox().yInverted() is False
    dv.deleteLater()


def test_cofit_dialog_editable_param_grid_and_preview(qapp):
    """The co-fit dialog exposes an editable parameter grid (value + Fix + the
    MQMAS F1 reference) and a Preview that simulates the current values."""
    from PySide6.QtCore import Qt

    from larmor import twod
    from larmor.desktop.cofit_dialog import CofitDialog

    x = np.linspace(-50, 120, 300); amp = np.exp(-0.5 * ((x - 60) / 8) ** 2) * 1e6
    base = {"kind": "1d", "label": "MAS", "ppm": x, "amp": amp,
            "nucleus": "27Al", "larmor": 130.3}
    f2 = np.linspace(-40, 110, 48); f1 = np.linspace(-30, 90, 40)
    Z = np.exp(-0.5 * (((f2[None, :] - 50) / 12) ** 2 + ((f1[:, None] - 65) / 12) ** 2))
    d2 = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="27Al", larmor_MHz=130.3)
    d2ds = {"kind": "2d", "label": "2rr", "data2d": d2,
            "nucleus": "27Al", "larmor": 130.3}
    recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
              "mqmas_f1_ref_ppm": 0.0, "mqmas_f1_ref_vary": True, "sites": [
                  {"model": "czjzek", "label": "AlIV", "params": {
                      "isotropic_chemical_shift_ppm": {"value": 60.0, "vary": True,
                                                       "min": 0, "max": 120},
                      "sigma_Cq_MHz": {"value": 1.6, "vary": True, "min": 0.2,
                                       "max": 8},
                      "shift_fwhm_ppm": {"value": 13.0, "vary": True, "min": 1,
                                         "max": 30},
                      "line_fwhm_ppm": {"value": 4.0, "vary": True, "min": 0},
                      "amplitude": {"value": 1e6, "vary": True, "min": 0}}}]}

    dlg = CofitDialog(None, recipe, base)
    dlg.datasets = [base, d2ds]
    dlg._build_param_rows()
    # a row per site-param plus a trailing MQMAS F1-reference row
    assert dlg._row_map[-1] == ("f1ref",)
    assert dlg.params_table.rowCount() == 6

    # edit δiso, fix σ, fix F1 ref at 12 ppm
    dlg.params_table.item(0, 2).setText("62")
    dlg.params_table.item(1, 3).setCheckState(Qt.Checked)
    lr = dlg.params_table.rowCount() - 1
    dlg.params_table.item(lr, 2).setText("12")
    dlg.params_table.item(lr, 3).setCheckState(Qt.Checked)
    dlg._apply_param_edits()
    sp = dlg.base_recipe["sites"][0]["params"]
    assert sp["isotropic_chemical_shift_ppm"]["value"] == 62.0
    assert sp["sigma_Cq_MHz"]["vary"] is False
    assert dlg.base_recipe["mqmas_f1_ref_ppm"] == 12.0
    assert dlg.base_recipe["mqmas_f1_ref_vary"] is False

    # Preview simulates (no fit) and draws both overlays
    import pyqtgraph as pg
    dlg._preview()
    n = sum(1 for i in range(dlg._plot_v.count())
            if isinstance(dlg._plot_v.itemAt(i).widget(), pg.PlotWidget))
    assert n == 2
    dlg.close()


def test_twod_display_modes_and_per_site_colours(qapp):
    """The 2D view offers contour/density/filled/values display modes and draws
    a fitted model's per-site components in the shared site colours."""
    from larmor import twod
    from larmor.desktop.plot import site_color
    from larmor.desktop.twod_view import Contour2DView

    f2 = np.linspace(-40, 110, 60); f1 = np.linspace(-30, 90, 50)
    Z = np.exp(-0.5 * (((f2[None, :] - 55) / 12) ** 2 + ((f1[:, None] - 65) / 12) ** 2))
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="27Al", larmor_MHz=195.5)
    dv = Contour2DView(); dv.set_data(d, "syn")
    assert [dv.disp.itemText(i) for i in range(dv.disp.count())] == \
        ["contour", "density", "filled", "contour+values"]
    for mode in ("contour", "density", "filled", "contour+values"):
        dv.disp.setCurrentText(mode)          # each mode redraws without error
    # per-site overlay keeps the components for site-coloured drawing
    s0 = np.exp(-0.5 * (((f2[None, :] - 55) / 8) ** 2 + ((f1[:, None] - 65) / 8) ** 2))
    s1 = 0.4 * np.exp(-0.5 * (((f2[None, :] - 20) / 8) ** 2 + ((f1[:, None] - 30) / 8) ** 2))
    dv.set_model(s0 + s1, f2, f1, per_site=[s0, s1])
    assert dv._model_sites is not None and len(dv._model_sites) == 2
    assert site_color(0) != site_color(1)
    dv.deleteLater()


def _cofit_czjzek_recipe():
    return {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
            "mqmas_f1_ref_ppm": 0.0, "mqmas_f1_ref_vary": True, "sites": [
                {"model": "czjzek", "label": "AlIV", "params": {
                    "isotropic_chemical_shift_ppm": {"value": 60.0, "vary": True,
                                                     "min": 0, "max": 120},
                    "sigma_Cq_MHz": {"value": 1.6, "vary": True, "min": 0.3,
                                     "max": 4},
                    "shift_fwhm_ppm": {"value": 13.0, "vary": True, "min": 1,
                                       "max": 30},
                    "line_fwhm_ppm": {"value": 4.0, "vary": True, "min": 0},
                    "amplitude": {"value": 1e6, "vary": True, "min": 0}}}]}


def _cofit_setup(win, r):
    """Put win into co-fit with a 1D and a self-consistent 2D, two recipes."""
    import json

    from larmor import twod
    from larmor.recipe import Recipe

    x = np.linspace(-50, 120, 400); amp = np.exp(-0.5 * ((x - 60) / 11) ** 2) * 1e6
    win.exp_ppm, win.exp_amp = x, amp; win.recipe = r
    kd = twod.Data2D(f2_ppm=np.linspace(-40, 110, 72),
                     f1_ppm=np.linspace(-30, 90, 60), z=np.zeros((60, 72)),
                     nucleus="27Al", larmor_MHz=130.3)
    k = twod._kernel_for(Recipe.from_dict(r), kd)
    Z, _ = twod.simulate_2d(Recipe.from_dict(r), k)
    d2 = twod.Data2D(f2_ppm=k.f2_ppm, f1_ppm=k.f1_ppm, z=Z / Z.max(),
                     nucleus="27Al", larmor_MHz=130.3)
    win._cofit = {"d1": (x, amp, "MAS"), "d2": (d2, "2rr"),
                  "home": win._cofit_home(),
                  "r1": json.loads(json.dumps(r)), "r2": json.loads(json.dumps(r)),
                  "tie": set(win._default_tie(r))}
    win.central_stack.setCurrentWidget(win.cofit_page)
    win._cofit_tie_rebuild(); win._cofit_rebuild_tables()
    return x, amp, d2


def test_cofit_two_recipes_tables_tie_and_fit(qapp, win):
    """Co-fit page has a parameter table per dataset (two independent recipes);
    the tie bar only offers the model's real params; untied params fit
    independently; Preview overlays each panel; Run fits both."""
    r = _cofit_czjzek_recipe()
    _cofit_setup(win, r)
    assert win._cofit_active()
    # tie bar shows only czjzek params that influence the lineshape (no Cq/eta/…)
    tie_keys = set(win._cofit_tie.keys())
    assert tie_keys == {"isotropic_chemical_shift_ppm", "sigma_Cq_MHz",
                        "shift_fwhm_ppm", "line_fwhm_ppm"}
    # two independent tables bound to r1 and r2
    assert win.cofit_table1d._recipe is win._cofit["r1"]
    assert win.cofit_table2d._recipe is win._cofit["r2"]

    # untie δiso, then a 1D edit must NOT propagate to the 2D recipe
    win._cofit_tie_toggled("isotropic_chemical_shift_ppm", False)
    win._cofit["r1"]["sites"][0]["params"][
        "isotropic_chemical_shift_ppm"]["value"] = 58.0
    win._cofit_on_edit(1)
    assert win._cofit["r2"]["sites"][0]["params"][
        "isotropic_chemical_shift_ppm"]["value"] == 60.0     # decorrelated

    # a TIED edit does propagate
    win._cofit["r1"]["sites"][0]["params"]["sigma_Cq_MHz"]["value"] = 2.2
    win._cofit_on_edit(1)
    assert win._cofit["r2"]["sites"][0]["params"]["sigma_Cq_MHz"]["value"] == 2.2

    win._cofit_simulate_now()
    assert win.cofit_view1d._model is not None
    assert win.cofit_view2d._model_sites is not None

    win.run_cofit_fit()
    assert "2D" in win.cofit_rmsd.text()
    win.close_cofit()
    assert win.central_stack.currentWidget() is not win.cofit_page


def test_cofit_split_stays_even_with_empty_2d(qapp, win):
    """The co-fit page splits 1D | 2D in half side by side — the 2D panel must
    keep its half even before a 2D map is added (it used to collapse to zero)."""
    from PySide6.QtCore import Qt

    win.resize(1400, 800); win.show(); qapp.processEvents()
    assert win.cofit_split.orientation() == Qt.Horizontal
    assert win.cofit_split.childrenCollapsible() is False

    x = np.linspace(-50, 120, 300); amp = np.exp(-0.5 * ((x - 60) / 11) ** 2) * 1e6
    win.exp_ppm, win.exp_amp = x, amp
    win._cofit = {"d1": (x, amp, "MAS"), "d2": None, "r1": None, "r2": None,
                  "tie": set()}
    win.central_stack.setCurrentWidget(win.cofit_page)
    win._cofit_split_even(); win._cofit_refresh_panels()
    for _ in range(4):
        qapp.processEvents()
    sizes = win.cofit_split.sizes(); total = sum(sizes) or 1
    assert sizes[1] > 0.4 * total, f"2D half collapsed: {sizes}"


def test_cofit_pauses_and_resumes_on_view_switch(qapp, win):
    """Switching the central view while co-fitting must leave co-fit cleanly and
    Decomposition ▸ Co-fit resumes it with both datasets still loaded."""
    r = _cofit_czjzek_recipe()
    _cofit_setup(win, r)

    win.central_stack.setCurrentWidget(win.view)      # user switches away
    assert win._cofit is None                          # no zombie state
    assert win._cofit_last is not None and win._cofit_last["d2"] is not None

    win.open_cofit()                                   # resume
    assert win.central_stack.currentWidget() is win.cofit_page
    assert win._cofit["d1"] is not None and win._cofit["d2"] is not None
    assert win._cofit.get("r1") and win._cofit.get("r2")

    win.close_cofit()
    assert win._cofit_last is None                     # deliberate close forgets it


def test_cofit_reseeds_on_new_dataset(qapp, win):
    """A stashed co-fit must NOT be restored over a *different* dataset: opening
    Co-fit while a new spectrum is loaded seeds fresh from that spectrum, rather
    than resurrecting the previous (e.g. example) co-fit."""
    r = _cofit_czjzek_recipe()
    _cofit_setup(win, r)                               # co-fit on dataset A
    win.source_path = "sampleA"; win._cofit["home"] = win._cofit_home()

    win.central_stack.setCurrentWidget(win.view)       # switch away -> stashed
    assert win._cofit is None and win._cofit_last is not None

    # user now opens a different 1D dataset (new "workspace")
    xb = np.linspace(-50, 120, 350)
    ampb = np.exp(-0.5 * ((xb - 30) / 9) ** 2) * 5e5
    win.exp_ppm, win.exp_amp = xb, ampb
    win.source_path = "sampleB"
    win.recipe = {**r, "sample": "B"}

    win._cofit_add_dataset = lambda *a, **k: None      # skip the "add 2nd" dialog
    win.open_cofit()
    assert win.central_stack.currentWidget() is win.cofit_page
    d1 = win._cofit["d1"]
    assert d1 is not None and d1[0].shape == xb.shape       # it's dataset B
    assert np.allclose(d1[1], ampb)                          # not the stash
    assert win._cofit["d2"] is None                          # fresh: no old 2D


def test_2d_reverse_tracks_model_overlay(qapp):
    """Reversing (or transposing) the 2D data axis mirrors the fitted model
    overlay too, and a later Compute (set_model) re-lands in the same
    orientation — so the co-fit model stays registered with the experiment."""
    from larmor import twod
    from larmor.desktop.twod_view import Contour2DView

    f2 = np.linspace(-20, 80, 40); f1 = np.linspace(-10, 60, 30)
    z = np.random.RandomState(0).rand(30, 40)
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=z, nucleus="27Al", larmor_MHz=130.3)
    v = Contour2DView(); v.set_data(d, "exp")
    mz = np.outer(np.hanning(30), np.hanning(40))
    v.set_model(mz, f2, f1, per_site=[mz])
    assert np.allclose(v._model[0], mz)

    v._op("rev_f2")                                    # user reverses F2
    assert v._model_ops == ["rev_f2"]
    assert np.allclose(v._model[0], np.flip(mz, 1))    # model followed the data
    assert np.allclose(v._model_sites[0], np.flip(mz, 1))

    v.set_model(mz, f2, f1, per_site=[mz])             # "Compute" re-simulates
    assert np.allclose(v._model[0], np.flip(mz, 1))    # still reversed (ops replayed)

    v._op("rev_f2")                                    # reverse again = undo
    assert np.allclose(v._model[0], mz)
    v.close()


def test_shift_recipe_positions_moves_only_positions():
    """The re-reference helper shifts absolute ppm positions, not widths."""
    from larmor.desktop.app import MainWindow
    r = {"sites": [
        {"params": {"isotropic_chemical_shift_ppm": {"value": 10.0},
                    "sigma_Cq_MHz": {"value": 2.0}}},
        {"params": {"shift_ppm": {"value": -5.0}}}]}
    MainWindow._shift_recipe_positions(r, 3.0)
    assert r["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] == 13.0
    assert r["sites"][0]["params"]["sigma_Cq_MHz"]["value"] == 2.0        # width kept
    assert r["sites"][1]["params"]["shift_ppm"]["value"] == -2.0


def test_calibrate_moves_model_with_axis(qapp, win, monkeypatch):
    """1D calibrate re-references the axis AND moves the fitted sites the same
    way, so the model stays on the peaks (was: only the axis shifted)."""
    from PySide6.QtWidgets import QInputDialog
    POS = "isotropic_chemical_shift_ppm"
    x = np.linspace(-40, 120, 400)
    win.exp_ppm = x; win.exp_amp = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    win.recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3, "sites": [
        {"model": "gauss_lor", "label": "A", "params": {
            POS: {"value": 60.0, "vary": True},
            "fwhm_ppm": {"value": 8.0, "vary": True},
            "gl": {"value": 0.0, "vary": True},
            "amplitude": {"value": 1.0, "vary": True}}}]}
    win.view.set_experiment(x, win.exp_amp)
    monkeypatch.setattr(QInputDialog, "getDouble",
                        staticmethod(lambda *a, **k: (70.0, True)))
    win.on_calibrate_picked(60.0)                    # "this 60-peak is really 70"
    assert abs(win.recipe["sites"][0]["params"][POS]["value"] - 70.0) < 1e-6


def test_2d_calibrate_moves_model_overlay(qapp):
    """2D re-referencing (_shift_axes) moves the fitted model overlay with the
    axes, like 1D calibrate."""
    from larmor import twod
    from larmor.desktop.twod_view import Contour2DView
    f2 = np.linspace(-20, 80, 40); f1 = np.linspace(-10, 60, 30)
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=np.random.RandomState(0).rand(30, 40),
                    nucleus="27Al", larmor_MHz=130.3)
    v = Contour2DView(); v.set_data(d, "exp")
    v.set_model(np.outer(np.hanning(30), np.hanning(40)), f2, f1)
    v._shift_axes(5.0, -3.0)
    assert np.allclose(v._model[1], f2 + 5.0)        # model F2 followed
    assert np.allclose(v._model[2], f1 - 3.0)        # model F1 followed
    v.close()


def test_fit_progress_bar_ticks(qapp, win):
    """The fit progress bar advances on lmfit iterations (iter_cb wired through
    the fit functions and the worker's progress signal)."""
    from larmor.desktop.app import FitWorker
    from larmor.engine import make_context, simulate_site
    from larmor.recipe import Recipe, SiteModel, Param

    x = np.linspace(-20, 120, 300)
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="gauss_lor", label="p", params={
            "isotropic_chemical_shift_ppm": Param(60, min=0, max=120),
            "shift_fwhm_ppm": Param(8, min=1, max=40),
            "gl": Param(0.5, min=0, max=1),
            "amplitude": Param(1e6, min=0)})])
    ctx = make_context(r, exp_ppm=x)
    y = np.sum([simulate_site(s, ctx) for s in r.sites], axis=0)
    ticks = []
    fw = FitWorker(r.to_dict(), x, y * 1.05, (120, -20))
    fw.progress.connect(lambda it, rms: ticks.append((it, rms)))
    fw.run()                                          # synchronous
    assert ticks and ticks[-1][0] >= 1
    # the bar runs as a busy/indeterminate marquee while fitting (no fixed % that
    # can stick near the end); the iteration count shows in its text
    win._progress_start("t")
    assert (win.progress.minimum(), win.progress.maximum()) == (0, 0)
    win._progress_tick(20, 0.01)
    assert "iter 20" in win.progress.format()
    win._progress_end(True)
    assert (win.progress.minimum(), win.progress.maximum()) == (0, 100)
    assert win.progress.value() == 100


def test_twopoint_background_subtracts_the_line_through_two_picks(qapp, win):
    """The 2-point background tool subtracts the straight line through the two
    picked points (a tilted flat baseline), leaves the peak, and is undoable."""
    x = np.linspace(-50, 50, 600)
    peak = 100.0 * np.exp(-0.5 * (x / 5.0) ** 2)
    tilt = 0.5 * x + 20.0                              # a sloped background
    win.exp_ppm = x; win.exp_amp = peak + tilt
    win.recipe = {"nucleus": "31P", "larmor_frequency_MHz": 162.0, "sites": [
        {"model": "gauss_lor", "label": "A", "params": {
            "isotropic_chemical_shift_ppm": {"value": 0.0, "vary": True},
            "shift_fwhm_ppm": {"value": 5.0, "vary": True},
            "gl": {"value": 1.0, "vary": False},
            "amplitude": {"value": 100.0, "vary": True}}}]}
    win.view.set_experiment(x, win.exp_amp)
    before = win.exp_amp.copy()
    # pick two baseline points, one each side, sitting on the tilt
    win._twopoint_mode(True)
    for xp in (-45.0, 45.0):
        win.view._add_baseline_anchor(xp, 0.5 * xp + 20.0)
    win._twopoint_mode(False)                          # turning off applies it
    edge = np.concatenate([win.exp_amp[:50], win.exp_amp[-50:]])
    assert abs(float(np.mean(edge))) < 1.0             # background flattened to ~0
    assert win.exp_amp.max() == pytest.approx(100.0, abs=0.5)   # peak preserved
    assert not win.proc_panel.btnTpPick.isChecked()    # pick mode released
    win.undo()
    assert np.allclose(win.exp_amp, before)            # undoable


def test_report_summary_never_forces_window_width(qapp, win):
    """A long fit summary must not widen the Report dock (and thus force the
    whole window past the screen): the label wraps and is horizontally
    shrinkable, so it contributes ~0 to the minimum width regardless of text."""
    from PySide6.QtWidgets import QSizePolicy

    assert win.results_summary.wordWrap()
    assert (win.results_summary.sizePolicy().horizontalPolicy()
            == QSizePolicy.Ignored)
    win.results_summary.setText(
        "RMSD 0.01  ·  χ² 1e6  ·  " + "unidentifiable (see Correlations)  ·  " * 12)
    # the label's minimum width stays near the widest word, not the sentence
    assert win.results_summary.minimumSizeHint().width() < 200


def test_theme_menu_has_a_hidden_aesthetic_submenu(win):
    """The 4 just-for-fun styles live in their own "More styles…" submenu
    (+ "Normal" to clear the override) -- separate from the 10 normal
    presets' own action group, so the everyday Theme list is unaffected."""
    from larmor.desktop import theme

    assert win._aesthetic_group is not None
    labels = {a.text() for a in win._aesthetic_group.actions()}
    assert labels == {"Normal"} | set(theme.aesthetic_names())
    assert win._theme_group is not None
    assert {a.text() for a in win._theme_group.actions()} == set(theme.names())


def test_set_aesthetic_override_applies_live_and_normal_restores_it(win):
    """Choosing a hidden style writes the setting AND applies it immediately
    (no restart needed); choosing "Normal" restores whichever normal theme
    was active before, also live."""
    from PySide6.QtCore import QSettings
    from larmor.desktop import theme

    settings = QSettings("LARMOR", "app")
    original_override = settings.value("appearanceOverride", "")
    original_theme = settings.value("theme", theme.DEFAULT)
    try:
        win._set_theme(original_theme)                 # known starting point
        win._set_aesthetic_override("Vaporwave")
        assert settings.value("appearanceOverride", "") == "Vaporwave"
        assert theme.active().name == "Vaporwave"       # applied live

        win._set_aesthetic_override("")                 # back to Normal
        assert settings.value("appearanceOverride", "") == ""
        assert theme.active().name == original_theme    # restored live
    finally:
        settings.setValue("appearanceOverride", original_override)
        win._set_theme(original_theme)


def _write_csv_spectrum(path, centres=((0.0, 3.0, 1.0), (15.0, 5.0, 0.45))):
    """A tiny CSV spectrum with the metadata header io/spectra reads."""
    ppm = np.linspace(-60.0, 100.0, 512)
    amp = np.zeros_like(ppm)
    for c, w, a in centres:
        amp += a * np.exp(-4 * np.log(2) * ((ppm - c) / w) ** 2)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# LARMOR spectrum\n# nucleus=11B\n"
                "# larmor_MHz=160.46\n# spin_rate_Hz=20000.0\n")
        f.write("ppm,intensity\n")
        for x, y in zip(ppm, amp):
            f.write(f"{x},{y}\n")
    return ppm, amp


def test_add_background_spectrum_attaches_a_real_trace(win, qapp, tmp_path,
                                                       monkeypatch):
    """Decomposition ▸ Add background spectrum on a recipe that already has
    lines. load_any returns (ppm, amp, recipe, ...); unpacking it as
    (recipe, ppm, amp, ...) put the recipe DICT into `amp`, so np.asarray(...,
    float) raised out of the slot and nothing was ever added."""
    from PySide6.QtWidgets import QFileDialog
    from larmor import models as model_registry
    from larmor.engine import make_context, simulate_site
    from larmor.recipe import Recipe

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert win.exp_amp.size > 0

    for centre in (0.0, 15.0):                 # two ordinary lines first
        m = model_registry.get("gauss_lor")
        params = {pd.name: {"value": pd.default, "stderr": None,
                            "vary": pd.vary, "min": pd.min, "max": pd.max,
                            "expr": None} for pd in m.params}
        params["isotropic_chemical_shift_ppm"]["value"] = centre
        win.recipe["sites"].append(
            {"model": "gauss_lor", "label": f"L{centre:g}", "params": params})
    win.on_structure_changed()
    assert len(win.recipe["sites"]) == 2

    bg = tmp_path / "bg.csv"
    bg_ppm, _ = _write_csv_spectrum(bg, centres=((-45.0, 4.0, 1.0),))
    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(bg), "")))
    win.add_background_spectrum()

    sites = win.recipe["sites"]
    assert len(sites) == 3, "the background spectrum was not added"
    s = sites[-1]
    assert s["model"] == "spectrum"
    assert len(s["ref"]["ppm"]) == bg_ppm.size
    assert len(s["ref"]["amp"]) == bg_ppm.size
    assert max(abs(v) for v in s["ref"]["amp"]) == pytest.approx(1.0)

    # and it renders something, both as-is and shifted
    r = Recipe.from_dict(win.recipe)
    ctx = make_context(r, exp_ppm=win.exp_ppm)
    assert np.any(simulate_site(r.sites[-1], ctx) != 0.0)
    r.sites[-1].params["shift_ppm"].value = 20.0
    assert np.any(simulate_site(r.sites[-1], ctx) != 0.0)


def test_spectrum_component_amplitude_may_be_negative(win, qapp, tmp_path,
                                                      monkeypatch):
    """A shifted copy of the same spectrum is how a satellite/sideband manifold
    is cancelled, so the spectrum component's amplitude must not be clamped at
    zero the way an ordinary lineshape's is."""
    from PySide6.QtWidgets import QFileDialog

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()

    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(data), "")))
    win.add_background_spectrum()
    assert win.recipe["sites"][-1]["params"]["amplitude"]["min"] is None


def test_add_current_spectrum_line_uses_the_processed_on_screen_trace(
        win, qapp, tmp_path, monkeypatch):
    """Decomposition ▸ Add a copy of this spectrum: the reference trace must be
    the PROCESSED data as displayed (re-reading the file would give the raw
    one, which is not what has to cancel), the shift must default to the MAS
    sideband spacing, and it must be held fixed — a copy of the same data with
    a free amplitude AND a free shift is degenerate with the whole model."""
    from PySide6.QtWidgets import QDialog

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()

    # a marker that exists only in the on-screen trace, never in the file
    win.exp_amp = win.exp_amp.copy()
    win.exp_amp[10] = 99.0

    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Accepted)
    win.add_current_spectrum_line()

    s = win.recipe["sites"][-1]
    assert s["model"] == "spectrum"
    ref_amp = np.asarray(s["ref"]["amp"], float)
    assert ref_amp.size == win.exp_ppm.size
    assert int(np.argmax(np.abs(ref_amp))) == 10, "not the on-screen trace"
    assert abs(ref_amp).max() == pytest.approx(1.0)          # unit peak

    nur_ppm = (win.recipe["spin_rate_Hz"]
               / win.recipe["larmor_frequency_MHz"])
    assert s["params"]["shift_ppm"]["value"] == pytest.approx(nur_ppm, rel=1e-3)
    assert s["params"]["shift_ppm"]["vary"] is False         # held by default
    assert s["params"]["amplitude"]["vary"] is True


def test_add_current_spectrum_line_cancelled_adds_nothing(win, qapp, tmp_path,
                                                          monkeypatch):
    from PySide6.QtWidgets import QDialog

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    before = len(win.recipe["sites"])

    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Rejected)
    win.add_current_spectrum_line()
    assert len(win.recipe["sites"]) == before


def test_axis_unit_khz_display(win, qapp, tmp_path, monkeypatch):
    """View > Axis unit: kHz relabels the bottom axis with ROUND kHz ticks
    (display value = ppm * SFO / 1000) while every plotted coordinate stays
    in ppm; switching back restores the ppm labelling."""
    from larmor.desktop.plot import ScaledAxis, axis_factor

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()

    ax = win.view.getPlotItem().getAxis("bottom")
    assert isinstance(ax, ScaledAxis)
    sfo = win.recipe["larmor_frequency_MHz"]        # 160.46

    win._set_axis_unit("kHz")
    f = axis_factor("kHz", sfo)
    assert ax._factor == pytest.approx(f)
    # ticks must be round in the DISPLAY unit, not in ppm
    vals = ax.tickValues(-60.0, 100.0, 800)
    major = vals[0][1]
    assert major, "no major ticks"
    for v in major:
        disp = v * f
        assert abs(disp - round(disp, 6)) < 1e-9 or \
            abs(disp * 10 - round(disp * 10)) < 1e-6, disp
    # tick STRINGS speak kHz: each real (round-in-kHz) tick renders exactly
    spacing = vals[0][0]
    strs = ax.tickStrings(major, 1.0, spacing)
    for v, txt in zip(major, strs):
        assert float(txt) == pytest.approx(v * f, abs=1e-6), (v, txt)
    # the cursor formatter leads with the unit and keeps ppm visible
    assert "kHz" in win._format_x(100.0) and "ppm" in win._format_x(100.0)

    win._set_axis_unit("ppm")
    assert ax._factor == 1.0
    assert win._format_x(100.0) == "100.00 ppm"


def test_integrals_dialog_reports_fwhm_khz_when_sfo_known(win, qapp, tmp_path):
    from larmor.desktop.integrate_dialog import IntegralsDialog

    data = tmp_path / "sample.csv"
    ppm, amp = _write_csv_spectrum(data)
    dlg = IntegralsDialog(win, ppm, amp, sfo_MHz=160.46)
    try:
        heads = [dlg.table.horizontalHeaderItem(j).text()
                 for j in range(dlg.table.columnCount())]
        assert heads[-1] == "FWHM (kHz)"
        dlg._recompute()
        assert dlg.table.rowCount() >= 1
        fwhm_ppm = float(dlg.table.item(0, 4).text())
        fwhm_khz = float(dlg.table.item(0, 5).text())
        assert fwhm_khz == pytest.approx(fwhm_ppm * 160.46 / 1000.0, abs=2e-3)
        # without an SFO the column is absent (plain solution-style use)
        dlg2 = IntegralsDialog(win, ppm, amp)
        assert dlg2.table.columnCount() == 5
        dlg2.close()
    finally:
        dlg.close()


def test_undo_snapshots_share_reference_traces(win, qapp, tmp_path, monkeypatch):
    """C1: a spectrum component's "ref" trace is shared across undo states
    (it is immutable by contract), not re-serialised on every action — and
    undo/redo still round-trips parameter edits correctly."""
    from PySide6.QtWidgets import QFileDialog

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(data), "")))
    win.add_background_spectrum()
    ref = win.recipe["sites"][-1]["ref"]

    win.snapshot()
    snap = win.undo_stack[-1]
    assert snap["recipe"]["sites"][-1]["ref"] is ref      # shared, not copied
    # but params are deep-copied: editing the live recipe leaves undo intact
    win.recipe["sites"][-1]["params"]["amplitude"]["value"] = 123.0
    assert snap["recipe"]["sites"][-1]["params"]["amplitude"]["value"] != 123.0

    win.undo(); qapp.processEvents()
    assert win.recipe["sites"][-1]["params"]["amplitude"]["value"] != 123.0
    win.redo(); qapp.processEvents()
    assert win.recipe["sites"][-1]["params"]["amplitude"]["value"] == 123.0


def test_session_is_a_file_not_the_registry(win, qapp, tmp_path, monkeypatch):
    """C2: the session persists to LOCALAPPDATA/LARMOR/session.json (atomic
    write), never to QSettings; the write path is gated on LARMOR_NO_SESSION
    (it was unguarded, so the test suite itself polluted the registry)."""
    fake_home = tmp_path / "appdata"
    monkeypatch.setenv("LOCALAPPDATA", str(fake_home))
    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()

    # under LARMOR_NO_SESSION (test default) nothing is written at all
    win._flush_session()
    assert not (fake_home / "LARMOR" / "session.json").exists()

    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    win._flush_session()
    f = fake_home / "LARMOR" / "session.json"
    assert f.exists()
    d = json.loads(f.read_text(encoding="utf-8"))
    assert d["source"] == str(data)
    assert d["recipe"]["nucleus"] == "11B"
    # and the debounce path schedules rather than writing synchronously
    f.unlink()
    win._persist_session()
    assert not f.exists() and win._session_timer.isActive()
    win._session_timer.stop()


def test_vocs_dialog_stitches_to_workbench(win, qapp, tmp_path, monkeypatch):
    """A3: three offset sub-spectra -> the VOCS dialog stitches them and the
    workbench receives the combined pattern with its provenance."""
    from PySide6.QtWidgets import QFileDialog

    from larmor.desktop.vocs_dialog import VocsDialog
    from larmor.io import spectra

    def gauss(x, c, w, a):
        return a * np.exp(-4 * np.log(2) * ((x - c) / w) ** 2)

    paths = []
    for i, centre in enumerate((-400.0, 0.0, 400.0)):
        x = np.linspace(centre - 320.0, centre + 320.0, 641)
        y = gauss(x, -300.0, 150.0, 1.0) + gauss(x, 250.0, 120.0, 0.7)
        f = tmp_path / f"off{i}.csv"
        spectra.write_csv(f, x, y, {"nucleus": "81Br",
                                    "larmor_MHz": 216.0})
        paths.append(str(f))

    dlg = VocsDialog(win)
    got = {}
    dlg.applied.connect(lambda p, a, n, s: got.update(
        ppm=np.asarray(p), amp=np.asarray(a), notes=n, sources=s))
    monkeypatch.setattr(QFileDialog, "getOpenFileNames",
                        staticmethod(lambda *a, **k: (paths, "")))
    dlg._add()
    assert dlg.btnApply.isEnabled()
    dlg._apply()
    assert got["ppm"].size > 100 and len(got["sources"]) == 3
    assert got["ppm"][0] < -600 and got["ppm"][-1] > 600

    win._vocs_to_workbench(got["ppm"], got["amp"], got["notes"],
                           got["sources"])
    assert win.exp_ppm.size == got["ppm"].size
    assert win.recipe["provenance"]["vocs_sources"] == got["sources"]
    assert win.recipe["mas_uncertain"] is True
    dlg.close()


def test_wurst_correct_divides_profile_and_records_provenance(
        win, qapp, tmp_path, monkeypatch):
    """A4: Process > WURST excitation profile divides the on-screen spectrum
    by the computed sweep weighting and records the parameters."""
    from PySide6.QtWidgets import QDialog

    from larmor.processing import wurst_profile

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    before = win.exp_amp.copy()
    sfo = win.recipe["larmor_frequency_MHz"]

    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Accepted)
    win.open_wurst_correct()

    prov = win.recipe["provenance"]["wurst_correct"]
    w = wurst_profile(win.exp_ppm, sfo, prov["centre_ppm"],
                      prov["sweep_kHz"], n=prov["n"], floor=prov["floor"])
    assert np.allclose(win.exp_amp, before / w, rtol=1e-12)
    # the profile is non-trivial: flat mid-sweep, clamped at the floor edges
    assert w.max() == pytest.approx(1.0, abs=1e-6)
    assert w.min() == pytest.approx(prov["floor"])
    win.undo(); qapp.processEvents()                 # and it is undoable
    assert np.allclose(win.exp_amp, before)


def test_kernel_warm_worker_spin_gate_and_dedup(win, qapp, tmp_path,
                                                monkeypatch):
    """B6: loading a quadrupolar 1D dataset kicks off ONE background kernel
    pre-build; spin-1/2 nuclei and repeat loads of the same dataset do not.
    (The build itself is exercised via the engine tests; here the worker is
    stubbed so the test stays fast.)"""
    import larmor.desktop.mw_files as fmod      # the module _warm_kernel reads

    started = []

    class FakeWorker:
        def __init__(self, nucleus, lar, spin, ppm):
            started.append(nucleus)
        def start(self):
            pass
        def isRunning(self):
            return False

    monkeypatch.setattr(fmod, "KernelWarmWorker", FakeWorker)
    monkeypatch.delenv("LARMOR_NO_KERNEL_WARM", raising=False)

    data = tmp_path / "sample.csv"
    _write_csv_spectrum(data)                       # 11B: spin 3/2
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert started == ["11B"]

    win._warm_kernel()                              # same dataset: deduped
    assert started == ["11B"]

    # a spin-1/2 dataset still constructs the worker (the spin gate lives in
    # the worker's run() so the lookup cost is off the GUI thread) -- what we
    # assert here is that a NEW dataset key warms again
    win.recipe["nucleus"] = "29Si"
    win.recipe["larmor_frequency_MHz"] = 99.3
    win._warmed_key = None
    win._warm_kernel()
    assert started == ["11B", "29Si"]


def test_staticct_dialog_reads_a_simulated_pattern(win, qapp, monkeypatch):
    """A8 end-to-end: markers placed on a simulated static 81Br pattern's
    true features produce a reading near the truth, and seeding creates a
    quad_ct site carrying it."""
    from larmor import engine
    from larmor.desktop.staticct_dialog import StaticCtDialog
    from larmor.recipe import Param, Recipe, SiteModel
    from larmor.staticct import ct_static_features

    cq, eta, diso, lar = 30.0, 0.5, -300.0, 216.0
    site = SiteModel(model="quad_ct", label="s", params={
        "isotropic_chemical_shift_ppm": Param(diso), "Cq_MHz": Param(cq),
        "eta": Param(eta), "shift_fwhm_ppm": Param(10.0),
        "amplitude": Param(1.0)})
    r = Recipe(nucleus="81Br", larmor_frequency_MHz=lar, spin_rate_Hz=0.0,
               sites=[site])
    x = np.linspace(-5200, 3000, 8001)
    gx, y, _ = engine.simulate(r, exp_ppm=x)

    win._display_1d(gx, np.clip(y, 0, None), "81Br", lar, 0.0, "sim", "sim")
    qapp.processEvents()

    dlg = StaticCtDialog(win, win.exp_ppm, win.exp_amp, "81Br", 1.5, lar)
    dlg.seed_site.connect(win._staticct_seed)   # what open_staticct wires
    f = ct_static_features(cq, eta, 1.5, lar)
    dists = [min(abs(e - h) for h in f.horns) for e in f.edges]
    dlg.horn_a.setValue(f.horns[0] + diso)
    dlg.horn_b.setValue(f.horns[1] + diso)
    dlg.edge.setValue(f.edges[int(np.argmax(dists))] + diso)
    qapp.processEvents()
    assert dlg.reading is not None and dlg.reading.ok
    assert dlg.reading.Cq_MHz == pytest.approx(cq, rel=0.05)
    assert dlg.reading.eta == pytest.approx(eta, abs=0.05)

    n_before = len(win.recipe["sites"])
    dlg._seed()
    sites = win.recipe["sites"]
    assert len(sites) == n_before + 1
    got = sites[-1]["params"]
    assert got["Cq_MHz"]["value"] == pytest.approx(cq, rel=0.05)
    assert got["eta"]["value"] == pytest.approx(eta, abs=0.05)
    dlg.close()


def test_wurst_correction_is_undoable(win, qapp, monkeypatch):
    """Process > WURST excitation profile rewrites exp_amp, so its undo
    snapshot must carry the axis: Ctrl+Z restores the spectrum, not only the
    recipe (Tutorial 7 says the step is undoable)."""
    from PySide6.QtWidgets import QDialog

    win._display_1d(np.linspace(-3000, 1000, 2001), np.ones(2001), "81Br",
                    216.0, 0.0, "sim", "sim")
    qapp.processEvents()
    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Accepted)
    before = win.exp_amp.copy()
    win.open_wurst_correct()
    assert not np.allclose(win.exp_amp, before)
    assert "wurst_correct" in (win.recipe.get("provenance") or {})
    win.undo()
    assert np.allclose(win.exp_amp, before)


# ------------------------------------------------ spinning-sideband offer (F3)
def _write_csv_manifold(path, rate_hz=20000.0, header_rate_hz=None, ratio=0.3,
                        n=2, fwhm=4.0):
    """A CSV 11B manifold at 160.46 MHz: a centreband at 0 ppm and ±k·νrot
    Gaussians of height ratio**k on a −400…400 ppm axis (4096 points); the
    header carries ``header_rate_hz`` (default: the true rate)."""
    lar = 160.46
    ppm = np.linspace(-400.0, 400.0, 4096)
    amp = np.zeros_like(ppm)
    spacing = rate_hz / lar
    for k in range(-n, n + 1):
        amp += ratio ** abs(k) * np.exp(
            -4 * np.log(2) * ((ppm - k * spacing) / fwhm) ** 2)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# LARMOR spectrum\n# nucleus=11B\n"
                f"# larmor_MHz={lar}\n# spin_rate_Hz={header_rate_hz or rate_hz}\n")
        f.write("ppm,intensity\n")
        for x, y in zip(ppm, amp):
            f.write(f"{x},{y}\n")
    return ppm, amp


def _banner_shown(win) -> bool:
    b = getattr(win, "ssb_banner", None)
    return b is not None and not b.isHidden()


def test_sideband_offer_appears_and_adds_a_linked_manifold_in_one_undo_step(
        win, qapp, tmp_path):
    """Loading a spectrum that repeats at ±νrot shows the banner with its
    comb guides; [Add linked manifold] builds one centreband plus one linked
    line per matched order -- position = parent ± k·νrot as a constraint the
    table renders as 'A+124.6', every shape parameter tied, amplitude free
    -- valid for the fit translator, undone in ONE step, and never offered
    again on a model that already carries a manifold."""
    from larmor import cellparse
    from larmor.fit import _make_params
    from larmor.recipe import Recipe

    data = tmp_path / "manifold.csv"
    _write_csv_manifold(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert _banner_shown(win)
    det = win._ssb_detection
    assert det.ok and det.nu_rot_Hz == pytest.approx(20000.0, rel=3e-3)
    assert len(win.ssb_banner.guides) == 5                # centre + 4 orders
    assert "⚠" not in win.ssb_banner.label.text()
    assert "20 000 Hz" in win.ssb_banner.label.text()
    assert win.recipe["sites"] == []

    win.ssb_banner.btnManifold.click()
    sites = win.recipe["sites"]
    assert len(sites) == 5
    parent, bands = sites[0], sites[1:]
    assert parent["model"] == "gauss_lor" and "sb" not in parent["label"]
    assert parent["params"]["isotropic_chemical_shift_ppm"]["value"] == \
        pytest.approx(0.0, abs=0.1)
    assert parent["params"]["shift_fwhm_ppm"]["value"] == pytest.approx(4.0, abs=0.4)
    assert parent["params"]["amplitude"]["value"] == pytest.approx(1.0, abs=0.05)
    spacing = 20000.0 / 160.46
    seen = []
    for s in bands:
        assert s["model"] == "gauss_lor"
        k = int(s["label"][-4:-2])                        # '...+1sb' -> +1
        seen.append(k)
        p = s["params"]
        off = k * det.spacing_ppm
        sign = "+" if off >= 0 else "-"
        assert p["isotropic_chemical_shift_ppm"]["expr"] == \
            f"s0.isotropic_chemical_shift_ppm {sign} {abs(off):.6g}"
        shown = cellparse.format_link(p["isotropic_chemical_shift_ppm"]["expr"],
                                      "isotropic_chemical_shift_ppm")
        assert shown.startswith("A+") or shown.startswith("A-")
        assert p["isotropic_chemical_shift_ppm"]["value"] == \
            pytest.approx(k * spacing, abs=0.5)
        assert p["shift_fwhm_ppm"]["expr"] == "s0.shift_fwhm_ppm"
        assert p["gl"]["expr"] == "s0.gl"
        assert p["amplitude"]["expr"] is None and p["amplitude"]["vary"] is True
        assert p["amplitude"]["value"] == pytest.approx(0.3 ** abs(k), rel=0.15)
        assert s["label"].endswith(f"{k:+d}sb")
    assert sorted(seen) == [-2, -1, 1, 2]
    _make_params(Recipe.from_dict(win.recipe))            # constraints valid
    assert win.recipe["spin_rate_Hz"] == 20000.0          # certain: untouched
    assert not _banner_shown(win) and win.ssb_banner.guides == []
    assert "linked sideband" in win.statusBar().currentMessage()

    win.undo()
    assert win.recipe["sites"] == []
    win.redo()
    assert len(win.recipe["sites"]) == 5
    win._maybe_offer_sidebands()
    assert not _banner_shown(win)          # linked positions ARE a manifold


def test_sideband_offer_copy_action_adds_held_copies_seeded_from_the_teeth(
        win, qapp, tmp_path):
    data = tmp_path / "manifold.csv"
    ppm, amp = _write_csv_manifold(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert _banner_shown(win)
    win.ssb_banner.btnCopy.click()
    sites = win.recipe["sites"]
    assert len(sites) == 2 and all(s["model"] == "spectrum" for s in sites)
    spacing = 20000.0 / 160.46
    shifts = sorted(s["params"]["shift_ppm"]["value"] for s in sites)
    assert shifts[0] == pytest.approx(-spacing, rel=1e-3)
    assert shifts[1] == pytest.approx(spacing, rel=1e-3)
    for s in sites:
        p = s["params"]
        assert p["shift_ppm"]["vary"] is False
        assert p["amplitude"]["vary"] is True and p["amplitude"]["min"] is None
        assert p["amplitude"]["value"] == pytest.approx(0.3 * amp.max(), rel=0.15)
        assert max(abs(v) for v in s["ref"]["amp"]) == pytest.approx(1.0)
        assert len(s["ref"]["ppm"]) == win.exp_ppm.size
    assert sorted(s["label"] for s in sites) == ["copy+1sb", "copy-1sb"]
    assert "shifted cop" in win.statusBar().currentMessage()
    win.undo()
    assert win.recipe["sites"] == []


def test_sideband_offer_warns_on_a_wrong_rate_and_writes_it_only_when_uncertain(
        win, qapp, tmp_path):
    from PySide6.QtCore import QSettings

    # (a) a certain rate 1.5 % off: the banner warns, the rate is left alone
    data = tmp_path / "manifold.csv"
    ppm, amp = _write_csv_manifold(data, rate_hz=20000.0, header_rate_hz=20300.0)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert _banner_shown(win)
    assert "⚠" in win.ssb_banner.label.text()
    assert "20 300" in win.ssb_banner.label.text()
    assert win._ssb_detection.nu_rot_Hz == pytest.approx(20000.0, rel=3e-3)
    win.ssb_banner.btnManifold.click()
    assert win.recipe["spin_rate_Hz"] == 20300.0
    assert "νrot set" not in win.statusBar().currentMessage()

    # (b) the same rate flagged uncertain: the write is due, undo restores it
    win._display_1d(ppm, amp, "11B", 160.46, 20300.0, "m", "m")
    qapp.processEvents()
    win.recipe["mas_uncertain"] = True
    win._update_mas_label()
    assert not win.mas_label.isHidden()
    win._maybe_offer_sidebands()
    assert _banner_shown(win)
    win.ssb_banner.btnManifold.click()
    assert win.recipe["spin_rate_Hz"] == pytest.approx(20000.0, rel=3e-3)
    assert win.recipe["mas_uncertain"] is False
    assert win.mas_label.isHidden()
    assert "νrot set to 20 000 Hz" in win.statusBar().currentMessage()
    assert "20000" in win.exp_label.text()
    win.undo()
    assert win.recipe["spin_rate_Hz"] == 20300.0
    assert "20300" in win.exp_label.text()
    assert win.recipe["sites"] == []

    # (c) the drop-down `sidebands` model line, seeded from the teeth
    win._display_1d(ppm, amp, "11B", 160.46, 20300.0, "m", "m")
    qapp.processEvents()
    win.recipe["mas_uncertain"] = True
    win._maybe_offer_sidebands()
    assert _banner_shown(win)
    win.ssb_banner.actModelLine.trigger()
    sites = win.recipe["sites"]
    assert len(sites) == 1 and sites[0]["model"] == "sidebands"
    p = sites[0]["params"]
    assert p["n_ssb"]["value"] == 2 and p["n_ssb"]["vary"] is False
    assert 0.2 < p["ssb_ratio"]["value"] < 0.4
    assert p["isotropic_chemical_shift_ppm"]["value"] == pytest.approx(0.0, abs=0.5)
    assert p["shift_fwhm_ppm"]["value"] == pytest.approx(4.0, abs=0.4)
    assert win.recipe["spin_rate_Hz"] == pytest.approx(20000.0, rel=3e-3)
    assert "`sidebands` line" in win.statusBar().currentMessage()
    assert QSettings("LARMOR", "app").value("ssbAutoOffer", True, type=bool) \
        == win.actSsbOffer.isChecked()


def test_sideband_offer_is_silent_for_plain_static_and_toggled_off_data(
        win, qapp, tmp_path, monkeypatch):
    from PySide6.QtCore import QSettings

    # (a) the default fixture: two lines, a 20 kHz header whose spacing
    # (124.6 ppm) exceeds half the 160 ppm axis -> nothing, and the menu
    # action says why
    plain = tmp_path / "plain.csv"
    _write_csv_spectrum(plain)
    win.load_source(str(plain), keep_fit=False)
    qapp.processEvents()
    assert not _banner_shown(win)
    win.detect_sidebands()
    assert not _banner_shown(win)
    assert win.statusBar().currentMessage().startswith("no ±νrot repeat")

    # (b) a confirmed-static recipe never runs the detector, even on data
    # that repeats
    ppm, amp = _write_csv_manifold(tmp_path / "m.csv")
    calls = []
    real = win._run_sideband_detection
    monkeypatch.setattr(win, "_run_sideband_detection",
                        lambda **kw: calls.append(kw) or real(**kw))
    win._display_1d(ppm, amp, "11B", 160.46, 0.0, "s", "s")
    qapp.processEvents()
    assert not _banner_shown(win) and calls == []
    # ... but a rate flagged uncertain is a question the data may answer
    win.recipe["mas_uncertain"] = True
    win._maybe_offer_sidebands()
    assert calls == [{"scan": True}] and _banner_shown(win)
    assert win._ssb_detection.scanned
    monkeypatch.undo()

    # (c) the View toggle, persisted in QSettings
    settings = QSettings("LARMOR", "app")
    before = settings.value("ssbAutoOffer", True, type=bool)
    try:
        win.actSsbOffer.setChecked(False)
        assert settings.value("ssbAutoOffer", True, type=bool) is False
        assert not _banner_shown(win)
        data = tmp_path / "manifold.csv"
        _write_csv_manifold(data)
        win.load_source(str(data), keep_fit=False)
        qapp.processEvents()
        assert not _banner_shown(win)
        win.actSsbOffer.setChecked(True)               # no reload needed
        assert settings.value("ssbAutoOffer", True, type=bool) is True
        assert _banner_shown(win)
    finally:
        settings.setValue("ssbAutoOffer", bool(before))
        win.actSsbOffer.setChecked(bool(before))


def test_sideband_offer_escape_dismiss_and_workspace_switch(win, qapp, tmp_path):
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QKeySequence
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QMenu

    data = tmp_path / "manifold.csv"
    _write_csv_manifold(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert _banner_shown(win)
    QTest.keyClick(win, Qt.Key_Escape)
    assert not _banner_shown(win) and win.ssb_banner.guides == []
    assert not any(a.isChecked() for a in win._model_actions.values())
    assert win._ssb_dismissed is not None
    # the same trace, re-processed to the same values: still dismissed
    win.apply_processing([{"op": "scale", "factor": 1.0}], False)
    qapp.processEvents()
    assert not _banner_shown(win)
    # another document in a new workspace, then back: the dismissal survives
    manifold_ws = win.active_ws
    plain = tmp_path / "plain.csv"
    _write_csv_spectrum(plain)
    win._ws_mode = "new"
    win.load_source(str(plain), keep_fit=False)
    qapp.processEvents()
    assert not _banner_shown(win)
    win.switch_workspace(manifold_ws)
    qapp.processEvents()
    assert win.exp_ppm.size == 4096 and not _banner_shown(win)
    # the menu action forgets the dismissal and re-offers
    win.detect_sidebands()
    assert _banner_shown(win)
    assert "Return adds the linked manifold" in win.statusBar().currentMessage()
    # a real rate change is a new key: Esc, then re-evaluate at a new rate
    QTest.keyClick(win, Qt.Key_Escape)
    assert not _banner_shown(win)
    win.recipe["spin_rate_Hz"] = 20150.0                   # 0.75 % off, in window
    win._maybe_offer_sidebands()
    assert _banner_shown(win) and "⚠" in win.ssb_banner.label.text()

    # menu presence and the shortcut
    dec = next(m for m in win.menuBar().findChildren(QMenu)
               if m.title() == "&Decomposition")
    act = next(a for a in dec.actions()
               if a.text().startswith("&Detect spinning sidebands"))
    assert act.shortcut() == QKeySequence("Ctrl+Shift+D")
    assert act is win.actDetectSsb
    view = next(m for m in win.menuBar().findChildren(QMenu)
                if m.title() == "&View")
    toggle = next(a for a in view.actions()
                  if a.text() == "&Offer spinning-sideband detection on load")
    assert toggle.isCheckable() and toggle is win.actSsbOffer


# ------------------------------------------------ MAS rate: three sources
def _mas_block(session, uncertain=True, **over):
    """The 2702-shaped provenance block: acqus 4200 / title 20 kHz /
    booking 22 kHz, three-way, on rotor RS2427418 in the 2026-05 session."""
    b = {"acqus_Hz": 4200.0, "title_Hz": 20000.0, "title_note": "",
         "title_raw": "MASR 20kHz", "booking_Hz": 22000.0, "booking_flag": True,
         "source": "highest", "uncertain": uncertain, "session": session,
         "rotor": "RS2427418", "nucleus": "31P", "confirmed": None,
         "confirmed_note": ""}
    b.update(over)
    return b


def test_experiment_dialog_lists_sources_and_use_buttons(win, qapp, tmp_path):
    from larmor.desktop.dialogs import ExperimentDialog

    ppm = np.linspace(-200, 200, 512)
    win._display_1d(ppm, np.exp(-(ppm / 5) ** 2), "31P", 242.79, 22000.0, "t", "e")
    win.recipe["provenance"] = {"mas_rate": _mas_block(str(tmp_path / "2026-05"))}
    win.recipe["mas_uncertain"] = True
    dlg = ExperimentDialog(win, win.recipe)
    assert set(dlg.source_rows) == {"acqus", "title", "booking", "measured"}
    from PySide6.QtWidgets import QLabel
    texts = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert any("4 200 Hz" in t for t in texts)
    assert any("20 000 Hz" in t for t in texts)
    assert any("22 000 Hz" in t for t in texts)
    assert any("experiment_addenda.xml" in t for t in texts)
    assert any('"MASR 20kHz"' in t for t in texts)
    assert any(t.startswith("Sources disagree") for t in texts)
    assert dlg.source_rows["measured"].text() == "Measure"
    assert not dlg.source_rows["measured"].isEnabled()    # no detector passed
    assert dlg.mas.value() == 22000.0
    dlg.source_rows["title"].click()
    assert dlg.mas.value() == 20000.0
    dlg.source_rows["acqus"].click()
    assert dlg.mas.value() == 4200.0
    dlg.source_rows["booking"].click()
    assert dlg.mas.value() == 22000.0
    assert dlg.remember_mas.isChecked() and not dlg.remember_mas.isHidden()
    assert "2026-05 · rotor RS2427418 · 31P" in dlg.remember_mas.text()
    assert dlg.forget_requested is False
    dlg.source_rows["title"].click()
    dlg._accept()
    assert win.recipe["spin_rate_Hz"] == 20000.0
    dlg.close()

    # a certain block: the checkbox is offered but unchecked; a confirmed
    # block shows the Forget button and its verdict
    win.recipe["provenance"]["mas_rate"] = _mas_block(
        str(tmp_path / "2026-05"), uncertain=False, source="title+booking",
        booking_Hz=20000.0)
    dlg = ExperimentDialog(win, win.recipe)
    assert not dlg.remember_mas.isChecked()
    texts = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert any("Title and booking sidecar agree; acqus MASR (4 200 Hz) is "
               "outvoted." == t for t in texts)
    assert not hasattr(dlg, "btnForget")
    dlg.close()
    win.recipe["provenance"]["mas_rate"] = _mas_block(
        str(tmp_path / "2026-05"), uncertain=False, source="confirmed",
        confirmed="2026-09-22T10:00:00")
    dlg = ExperimentDialog(win, win.recipe)
    texts = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert any(t == "Confirmed on 2026-09-22 for 2026-05 · rotor RS2427418 · 31P."
               for t in texts)
    dlg.btnForget.click()
    assert dlg.forget_requested is True and not dlg.remember_mas.isChecked()
    dlg.close()


def test_experiment_dialog_without_block_is_unchanged(win, qapp, tmp_path,
                                                      monkeypatch):
    from PySide6.QtWidgets import QGroupBox
    from larmor.desktop.dialogs import ExperimentDialog

    data = tmp_path / "plain.csv"
    _write_csv_manifold(data, rate_hz=20000.0)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert not (win.recipe.get("provenance") or {}).get("mas_rate")
    dlg = ExperimentDialog(win, win.recipe)
    assert dlg.source_rows == {}
    assert dlg.findChildren(QGroupBox) == []
    assert dlg.remember_mas.isHidden() and not dlg.remember_mas.isChecked()
    assert dlg.forget_requested is False
    assert dlg.nucleus.text() == "11B"
    assert dlg.larmor.value() == pytest.approx(160.46)
    assert dlg.mas.value() == 20000.0
    dlg.mas.setValue(12500.0)
    dlg._accept()
    assert win.recipe["spin_rate_Hz"] == 12500.0
    dlg.close()
    # edit_experiment on a CSV recipe writes nothing to the store
    from larmor import masrate
    store = tmp_path / "store.jsonl"
    monkeypatch.setenv("LARMOR_MAS_LOG", str(store))
    assert win._commit_mas_choice(dlg) == ""
    assert not store.exists()
    assert masrate.log_path() == store


def test_experiment_dialog_shows_the_acquisition_summary(qapp):
    """N4: a recipe with an acquisition block gets a read-only summary row
    under SR (acqus facts, the flip angle with its source, the MAS sources
    and verdict); a recipe without one looks as before."""
    from larmor.desktop.dialogs import ExperimentDialog

    rec = {"nucleus": "31P", "larmor_frequency_MHz": 242.79, "spin_rate_Hz": 22000.0,
           "mas_uncertain": True, "sr_hz": -210.86,
           "acquisition": {"pulprog": "zg", "ns": 14, "d1_s": 300.0, "p1_us": 1.25,
                           "flip_deg": 30.0, "flip_source": "title", "plw1_W": 207.0,
                           "probe": "SPRB600511_7297 (MAS)", "mas_acqus_Hz": 4200.0,
                           "mas_title_Hz": 20000.0, "mas_booking_Hz": 22000.0,
                           "spin_rate_Hz": 22000.0, "mas_uncertain": True, "procno": 1,
                           "wdw": "EM", "lb_hz": 0.0, "tdeff": 384, "si": 32768,
                           "ph_mod": "pk", "absg": 5}}
    dlg = ExperimentDialog(None, rec)
    txt = dlg.acq_label.text()
    for must in ("NS 14", "D1 300 s", "30°", "acqus 4200 Hz", "title 20 kHz", "booking 22 kHz",
                 "(confirm)", "SPRB600511_7297", "TDeff 384"):
        assert must in txt, must
    assert "fitted with" not in txt
    dlg.close()
    dlg2 = ExperimentDialog(None, {**rec, "software": {"larmor": "0.13.0", "mrsimulator": "1.0.0"}})
    assert "fitted with LARMOR 0.13.0" in dlg2.acq_label.text()
    dlg2.close()
    plain = ExperimentDialog(None, {"nucleus": "11B", "larmor_frequency_MHz": 160.0,
                                    "spin_rate_Hz": 0.0})
    assert not hasattr(plain, "acq_label")
    plain.close()


def test_edit_experiment_remembers_and_clears_badge(win, qapp, monkeypatch, tmp_path):
    from larmor import masrate
    from larmor.desktop.dialogs import ExperimentDialog

    store = tmp_path / "mas_confirmations.jsonl"
    monkeypatch.setenv("LARMOR_MAS_LOG", str(store))
    ppm = np.linspace(-200, 200, 512)
    win._display_1d(ppm, np.exp(-(ppm / 5) ** 2), "31P", 242.79, 22000.0, "t", "e")
    session = str(tmp_path / "2026-05")
    win.recipe["provenance"] = {"mas_rate": _mas_block(session)}
    win.recipe["mas_uncertain"] = True
    win._update_exp_label()                  # strip tooltip + badge
    assert not win.mas_label.isHidden()
    assert "acqus 4 200 Hz · title 20 000 Hz · booking 22 000 Hz disagree" \
        in win.mas_label.toolTip()
    assert "highest of disagreeing sources" in win.exp_label.toolTip()

    def fake_exec(dlg_self):             # the user kept 22 000 Hz and pressed OK
        assert dlg_self.remember_mas.isChecked()
        dlg_self.mas.setValue(22000.0)
        dlg_self._accept()
        return 1

    monkeypatch.setattr(ExperimentDialog, "exec", fake_exec)
    win.edit_experiment()
    key = (session, "RS2427418", "31P")
    ev = masrate.MasEvidence.from_dict(_mas_block(session))
    assert masrate.lookup(key, ev)["rate_Hz"] == 22000.0
    assert win.recipe["mas_uncertain"] is False
    assert win.mas_label.isHidden()
    block = win.recipe["provenance"]["mas_rate"]
    assert block["source"] == "confirmed" and block["confirmed"]
    msg = win.statusBar().currentMessage()
    assert "νrot 22 000 Hz confirmed" in msg
    assert "remembered for 2026-05 · rotor RS2427418 · 31P" in msg
    assert "confirmed on" in win.exp_label.toolTip()

    # a second visit with [Forget]: the store gets a reversal
    def fake_exec_forget(dlg_self):
        assert hasattr(dlg_self, "btnForget")
        dlg_self.btnForget.click()
        dlg_self._accept()
        return 1

    monkeypatch.setattr(ExperimentDialog, "exec", fake_exec_forget)
    win.edit_experiment()
    assert masrate.lookup(key, ev) is None
    assert block["confirmed"] is None and block["source"] == "highest"
    assert "forgotten" in win.statusBar().currentMessage()
    assert win.recipe["mas_uncertain"] is False
    lines = store.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2 and json.loads(lines[1])["action"] == "forget"

    # OK with Remember unchecked: nothing written, badge still cleared
    def fake_exec_plain(dlg_self):
        dlg_self.remember_mas.setChecked(False)
        dlg_self._accept()
        return 1

    win.recipe["mas_uncertain"] = True
    monkeypatch.setattr(ExperimentDialog, "exec", fake_exec_plain)
    win.edit_experiment()
    assert len(store.read_text(encoding="utf-8").splitlines()) == 2
    assert win.recipe["mas_uncertain"] is False and win.mas_label.isHidden()
    assert "remembered" not in win.statusBar().currentMessage()


def test_experiment_dialog_measure_uses_sideband_detector(win, qapp, tmp_path):
    from larmor.desktop.dialogs import ExperimentDialog

    ppm, amp = _write_csv_manifold(tmp_path / "m.csv", rate_hz=20000.0,
                                   header_rate_hz=20300.0)
    win._display_1d(ppm, amp, "11B", 160.46, 20300.0, "m", "m")
    qapp.processEvents()
    win.recipe["provenance"] = {"mas_rate": _mas_block(
        str(tmp_path / "2026-05"), title_Hz=20300.0, booking_Hz=None,
        source="highest", nucleus="11B")}
    dlg = ExperimentDialog(win, win.recipe,
                           measure=lambda: win._run_sideband_detection(scan=True))
    assert set(dlg.source_rows) == {"acqus", "title", "measured"}
    btn = dlg.source_rows["measured"]
    assert btn.isEnabled() and btn.text() == "Measure"
    det = dlg._measure()
    assert det is not None and det.ok
    assert dlg.measured_label.text().startswith("Spinning sidebands: νrot 20 0")
    assert btn.text() == "Use"
    btn.click()
    assert dlg.mas.value() == pytest.approx(20000.0, rel=3e-3)
    dlg.close()


def test_datasets_compare_button_and_overlay_marker(qapp, win, monkeypatch):
    """Datasets dock: 'Compare acquisition…' wakes once an overlay with a
    source exists, an overlay acquired / processed unlike the compared set's
    majority gets a ⚠ marker on its detail line, the comparison is cached
    by the tuple of sources (a stack-offset change does not re-read), and a
    set without Bruker files says so in the status bar."""
    from pathlib import Path
    from conftest import fake_spectrum_params
    from larmor import comparability
    from larmor.desktop import comparability_dialog, datasets

    assert not win.datasets_panel.btnCompare.isEnabled()
    fakes = {"srcA": fake_spectrum_params("srcA", lb=0),
             "srcB": fake_spectrum_params("srcB", lb=100)}
    monkeypatch.setattr(comparability, "read_params",
                        lambda p, procno=None: fakes.get(Path(str(p)).name))
    real_compare = comparability.compare
    calls = {"n": 0}

    def counting_compare(params, labels=None):
        calls["n"] += 1
        return real_compare(params, labels)

    monkeypatch.setattr(comparability, "compare", counting_compare)
    x = np.linspace(-100, 100, 200)
    win._display_1d(x, np.exp(-(x ** 2) / 200), "11B", 192.43, None, "A", "srcA")
    win.source_path = "srcA"          # _display_1d leaves it to load_source
    win._refresh_overlays()
    qapp.processEvents()
    assert not win.datasets_panel.btnCompare.isEnabled()      # no overlay yet
    win._add_overlay("B", x, 0.5 * np.exp(-((x - 20) ** 2) / 200), "srcB")
    qapp.processEvents()
    assert win.datasets_panel.btnCompare.isEnabled()
    assert win._overlays[0]["comparability"] == "LB 100 Hz"
    detail = datasets.DatasetsPanel._detail_of(win._overlays[0])
    assert "⚠" in detail and "LB 100 Hz" in detail
    n = calls["n"]
    assert n >= 1
    win.datasets_panel.offset.setValue(0.4)                    # a second _refresh_overlays
    qapp.processEvents()
    assert calls["n"] == n                                     # cached by the source tuple
    win.overlay_set_color(0, "#123456")
    assert calls["n"] == n
    seen = []
    monkeypatch.setattr(comparability_dialog.ComparabilityDialog, "exec",
                        lambda self: seen.append(self._cmp) or 0)
    win.datasets_panel.compare_requested.emit()
    qapp.processEvents()
    assert len(seen) == 1 and seen[0].report("LB").deviants == [1]
    assert seen[0].labels[1] == "B"
    # nothing Bruker among the compared spectra: a status-bar line, no dialog
    monkeypatch.setattr(comparability, "read_params", lambda p, procno=None: None)
    win._add_overlay("C", x, 0.3 * np.exp(-((x + 30) ** 2) / 200), "srcC")
    qapp.processEvents()
    assert win._overlays[0]["comparability"] == ""              # recomputed: sources changed
    win.compare_overlays()
    assert "no Bruker acquisition files" in win.statusBar().currentMessage()
    assert len(seen) == 1


def test_fit_worker_stop_and_cancel_return_within_a_second_mid_evaluation(qapp, monkeypatch):
    """WA: Stop / Cancel during a LONG residual evaluation (six slow sites,
    0.48 s per evaluation) come back within a second -- the cancel hook the
    worker registers is checked between site renders; the latency is one
    render to notice plus the final full-grid render of the kept point (six
    renders again), so the sleeps stay short. Stop returns a result
    (parameters kept), Cancel reports its mode (the app reverts)."""
    import dataclasses
    import time

    from larmor.desktop.app import FitWorker
    from larmor.engine import make_context, simulate_site
    from larmor.models import base as mbase
    from larmor.recipe import Param, Recipe, SiteModel

    model = mbase.REGISTRY["gauss_lor"]
    orig = model.render

    def slow(values, ctx):
        time.sleep(0.08)
        return orig(values, ctx)

    monkeypatch.setitem(mbase.REGISTRY, "gauss_lor",
                        dataclasses.replace(model, render=slow))
    x = np.linspace(-20.0, 120.0, 300)
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="gauss_lor", label=f"p{i}", params={
            "isotropic_chemical_shift_ppm": Param(20.0 + 15.0 * i, min=-50, max=150),
            "shift_fwhm_ppm": Param(6.0, min=1, max=40),
            "gl": Param(0.5, vary=False),
            "amplitude": Param(1e6, min=0)}) for i in range(6)])
    ctx = make_context(r, exp_ppm=x)
    y = 1.3 * np.sum([simulate_site(s, ctx) for s in r.sites], axis=0)

    for mode in ("stop", "cancel"):
        fw = FitWorker(r.to_dict(), x, y, (120.0, -20.0))
        got, ticks = {}, []
        fw.done.connect(lambda res, m, got=got: got.update(res=res, mode=m))
        fw.failed.connect(lambda msg, got=got: got.update(failed=msg))
        fw.progress.connect(lambda it, rms, ticks=ticks: ticks.append(time.perf_counter()))
        fw.start()
        t0 = time.perf_counter()
        while len(ticks) < 2 and time.perf_counter() - t0 < 60:
            qapp.processEvents()
        assert len(ticks) >= 2
        # 0.2 s into an evaluation (three of its six renders still to come)
        while time.perf_counter() - ticks[-1] < 0.2:
            qapp.processEvents()
        t_req = time.perf_counter()
        fw.request_stop(mode)
        while fw.isRunning() and time.perf_counter() - t_req < 30:
            qapp.processEvents()
        latency = time.perf_counter() - t_req
        fw.wait()
        qapp.processEvents()
        assert latency < 1.0, latency
        assert "failed" not in got, got
        assert got.get("mode") == mode
        assert got.get("res") is not None and got["res"].lmfit_result.aborted


def test_linked_sidebands_follow_a_change_of_spin_rate(qapp, win, monkeypatch):
    """WA: sidebands added at 20 kHz move to +-22000/nu0 ppm when nu_rot
    becomes 22 kHz -- through the Experiment dialog, the detector's Use and
    a direct refresh -- constraint, value and drawn paddle alike."""
    from PySide6.QtWidgets import QDialog

    from larmor import sidebands as sb
    from larmor.desktop import dialogs as _dialogs

    win.load_source(str(require(CAALGLASS)), keep_fit=False)
    qapp.processEvents()
    lar = float(win.recipe["larmor_frequency_MHz"])
    win.recipe["spin_rate_Hz"] = 20000.0
    win.recipe["mas_uncertain"] = False
    win.recipe["sites"] = []
    win._model_actions["gauss_lor"].setChecked(True)
    win.add_site_at(60.0, 1000.0)
    qapp.processEvents()

    # Add spinning sidebands... (its dialog accepted with the defaults: one
    # forward, one backward, linked)
    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Accepted)
    win.add_sidebands()
    qapp.processEvents()
    assert len(win.recipe["sites"]) == 3
    d20 = 20000.0 / lar
    copies = {s["sideband"]["k"]: s for s in win.recipe["sites"][1:]}
    assert set(copies) == {1, -1}
    for k, s in copies.items():
        p = s["params"]["isotropic_chemical_shift_ppm"]
        assert s["sideband"]["parent"] == 0
        assert p["expr"] == sb.linked_position_expr(0, k, d20)
        assert p["value"] == pytest.approx(60.0 + k * d20)

    def drawn():
        return {p.index: p._pos for p in win.view._paddles}

    assert drawn()[1] == pytest.approx(60.0 + d20)

    # 1) the Experiment dialog: OK with nu_rot typed as 22 kHz
    def fake_exec(self):
        self.mas.setValue(22000.0)
        self._accept()
        return QDialog.Accepted

    monkeypatch.setattr(_dialogs.ExperimentDialog, "exec", fake_exec)
    win.edit_experiment()
    qapp.processEvents()
    d22 = 22000.0 / lar
    assert win.recipe["spin_rate_Hz"] == 22000.0
    for k, s in copies.items():
        p = s["params"]["isotropic_chemical_shift_ppm"]
        assert p["expr"] == sb.linked_position_expr(0, k, d22)
        assert p["value"] == pytest.approx(60.0 + k * d22)
    pos = drawn()
    assert pos[1] == pytest.approx(60.0 + d22) and pos[2] == pytest.approx(60.0 - d22)
    # the table shows the new offset in the letter form
    from larmor import cellparse
    assert cellparse.format_link(copies[1]["params"]["isotropic_chemical_shift_ppm"]["expr"],
                                 "isotropic_chemical_shift_ppm") == f"A+{d22:.6g}"
    # the simulated components sit on the new comb
    from larmor import engine
    from larmor.recipe import Recipe
    x, total, per_site = engine.simulate(Recipe.from_dict(win.recipe),
                                         exp_ppm=win.exp_ppm)
    for i, k in ((1, 1), (2, -1)):           # within the axis step (~0.5 ppm)
        assert x[int(np.argmax(per_site[i]))] == pytest.approx(60.0 + k * d22, abs=1.0)
    # Ctrl+Z brings the 20 kHz comb back together with the rate
    win.undo()
    assert win.recipe["spin_rate_Hz"] == 20000.0
    assert win.recipe["sites"][1]["params"]["isotropic_chemical_shift_ppm"]["expr"] == \
        sb.linked_position_expr(0, 1, d20)

    # 2) the detector's rate write (an uncertain recorded rate answered by
    #    the data) moves them too
    class _Det:
        nu_rot_Hz = 24000.0
    win.recipe["mas_uncertain"] = True
    assert win._ssb_set_rate_if_due(_Det()) is True
    d24 = 24000.0 / lar
    assert win.recipe["sites"][2]["params"]["isotropic_chemical_shift_ppm"]["expr"] == \
        sb.linked_position_expr(0, -1, d24)

    # 3) a recipe without the marker (made before it existed) keeps its
    #    constant offsets
    win.recipe["sites"][1].pop("sideband")
    win.recipe["spin_rate_Hz"] = 26000.0
    assert win._refresh_sideband_exprs() == 1          # only the marked copy
    assert win.recipe["sites"][1]["params"]["isotropic_chemical_shift_ppm"]["expr"] == \
        sb.linked_position_expr(0, 1, d24)


def test_baseline_iterative_dialog_constructs_previews_and_its_buttons_work(
        qapp, win, monkeypatch):
    """The 'Baseline iterative' dialog (Yon et al. 2020) crashed on open for
    months (PlotItem has no getPlotItem) and no test caught it: construct it
    on real data, check the live preview, the parameter change path, Apply
    (accept) and Cancel (reject), and the params it hands back -- and the
    workbench action behind it (which compared against ``dlg.Accepted`` on
    the INSTANCE, an AttributeError in PySide6 6.x, so Apply crashed after
    the dialog closed; found by this test)."""
    from PySide6.QtWidgets import QDialog, QDialogButtonBox

    from larmor.desktop.baseline_dialog import BaselineDialog

    win.load_source(str(require(BRUKER_1R)), keep_fit=False)
    qapp.processEvents()
    dlg = BaselineDialog(win, win.exp_ppm, win.exp_amp)
    qapp.processEvents()
    # the preview ran on construction: three curves on the two plots
    assert dlg.c_raw.xData is not None and len(dlg.c_raw.xData) == win.exp_ppm.size
    assert dlg.c_base.yData is not None and np.isfinite(dlg.c_base.yData).all()
    assert dlg.c_corr.yData is not None and len(dlg.c_corr.yData) == win.exp_ppm.size
    assert "iterations" in dlg.status.text()
    assert dlg.params() == {"dead_time_pts": 0, "smoothness": 1.0,
                            "threshold_factor": 1.0}
    # a spinbox edit re-previews (debounced) with the new parameters
    dlg.sp_smooth.setValue(4.0)
    dlg.sp_thr.setValue(0.6)
    dlg._timer.stop()
    dlg._preview()
    assert dlg.params()["smoothness"] == 4.0 and dlg.params()["threshold_factor"] == 0.6
    assert np.isfinite(dlg.c_base.yData).all()
    # the buttons: Apply accepts, Cancel rejects
    bb = dlg.findChild(QDialogButtonBox)
    bb.button(QDialogButtonBox.Apply).click()
    assert dlg.result() == QDialog.Accepted
    dlg2 = BaselineDialog(win, win.exp_ppm, win.exp_amp)
    bb2 = dlg2.findChild(QDialogButtonBox)
    bb2.button(QDialogButtonBox.Cancel).click()
    assert dlg2.result() == QDialog.Rejected
    # and the workbench action applies the step as recorded processing
    before = win.exp_amp.copy()
    monkeypatch.setattr(BaselineDialog, "exec", lambda self: QDialog.Accepted)
    win.apply_iterbaseline()
    qapp.processEvents()
    ops = win.recipe.get("processing") or []
    assert any(op.get("op") == "iterbaseline" for op in ops)
    assert not np.array_equal(before, win.exp_amp)
