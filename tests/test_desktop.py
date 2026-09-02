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
    import larmor.desktop.app as appmod

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
