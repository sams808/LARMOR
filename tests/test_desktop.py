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
def win(qapp):
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
