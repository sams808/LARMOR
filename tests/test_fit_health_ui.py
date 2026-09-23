"""MainWindow integration of the fit-health strip: fed by every fit (Auto Fit
included), re-evaluated on edits, reset on document changes, restored by a
workspace switch, click-throughs, the F7 menu and the Panels toggle.

Synthetic data throughout (the 300-point gauss_lor recipe of
test_desktop.py::test_fit_progress_bar_ticks); one real-data test skips
cleanly when CaAlGlass.fxmla is absent. _fit_done / _sim_done are called
directly (as test_desktop.py does) instead of spinning worker threads.
"""
import os

import numpy as np
import pytest

from conftest import CAALGLASS, require

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402

STALE = " · edited since fit"


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")   # never inherit a real session
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _recipe():
    from larmor.recipe import Param, Recipe, SiteModel

    return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="gauss_lor", label="p", params={
            "isotropic_chemical_shift_ppm": Param(60, min=0, max=120),
            "shift_fwhm_ppm": Param(8, min=1, max=40),
            "gl": Param(0.5, min=0, max=1),
            "amplitude": Param(1e6, min=0)})])


def _synthetic():
    from larmor.engine import make_context, simulate_site

    x = np.linspace(-20, 120, 300)
    r = _recipe()
    ctx = make_context(r, exp_ppm=x)
    y = np.sum([simulate_site(s, ctx) for s in r.sites], axis=0) * 1.05
    y = y + np.random.RandomState(0).normal(0, 0.005 * y.max(), y.size)
    return x, y, r


def _fitted(win, qapp, title="t", src="src"):
    """Synthetic spectrum on the workbench, fitted, fed through _fit_done."""
    from larmor import fit as fitmod
    from larmor.recipe import Recipe

    x, y, r = _synthetic()
    win._display_1d(x, y, "27Al", 130.3, None, title, src)
    win.recipe["sites"] = r.to_dict()["sites"]
    result = fitmod.fit(Recipe.from_dict(win.recipe), x, y, window_ppm=(120, -20))
    win._fit_done(result)
    qapp.processEvents()
    return x, y, result


def test_fit_done_feeds_strip_report_header_chi2_and_status(qapp, win):
    x, y, result = _fitted(win, qapp)
    strip = win.health_strip
    assert strip.isVisibleTo(win)
    assert strip.pill.text().startswith(("✓ Fit", "⚠ Fit"))
    assert win._health.fitted and not win._health.stale
    assert win._health is win._health_fit
    assert win.results_summary.text().startswith("RMSD ")
    assert win.lines_table.chi2.text() == win._health.chi_text()
    assert strip.rmsd_chip.text() == win._health.chi_text()
    assert win._last_lmfit is result.lmfit_result
    assert "fit done  ·  " in win.statusBar().currentMessage()
    assert win.actHealth.isEnabled()
    assert win.qtable.rowCount() == 1                # quantify still ran


def test_edit_after_fit_goes_stale_then_compute_keeps_the_verdict(qapp, win):
    x, y, result = _fitted(win, qapp)
    sentence = win.results_summary.text()
    p = win.recipe["sites"][0]["params"]["isotropic_chemical_shift_ppm"]
    orig = p["value"]
    p["value"] = orig + 15
    win._sim_done(x, 0.7 * y, [0.7 * y])            # the debounced sim landing
    t = win.health_strip.pill.text()
    assert t.endswith(STALE) and t.startswith(("⚠ Model", "✗ Model", "Model"))
    assert win._health.stale and win._health is not win._health_fit
    assert win.results_summary.text() == sentence   # the header documents the fit
    p["value"] = orig                               # back to the fitted values
    win._sim_done(np.asarray(result.x_ppm), np.asarray(result.y_fit), result.per_site)
    assert win.health_strip.pill.text().startswith(("✓ Fit", "⚠ Fit"))
    assert win._health is win._health_fit           # signature equal → early return


def test_paddle_drag_only_marks_stale_and_release_does_the_full_pass(qapp, win):
    x, y, result = _fitted(win, qapp)
    chips_before = [c.text() for c in win.health_strip.chips]
    win._paddle_live = True
    win.recipe["sites"][0]["params"]["amplitude"]["value"] *= 0.5
    win._sim_done(x, 0.5 * y, [0.5 * y])
    assert win.health_strip.pill.text().endswith(STALE)
    assert [c.text() for c in win.health_strip.chips] == chips_before   # no new chip yet
    assert not win._health.stale                    # only the pill hint moved
    win.on_paddle_released(0)
    assert win._health.stale
    texts = [c.text() for c in win.health_strip.chips]
    assert any(t.startswith("residual ") and t.endswith("× noise") for t in texts), texts


def test_unphysical_edit_shows_red_live_flag_and_focus_click_through(qapp, win):
    from larmor.desktop import theme

    x, y, result = _fitted(win, qapp)
    win.recipe["sites"][0]["params"]["gl"]["value"] = 1.4
    win._sim_done(x, y, [y])
    strip = win.health_strip
    assert strip.pill.text() == "✗ Model: not physical" + STALE
    assert "Gauss/Lorentz mix = 1.4 outside [0, 1]" in strip.pill.toolTip()
    assert strip.colours("bad")[0] == theme.active().model
    assert [c.text() for c in strip.chips][0] == "unphysical ×1"
    win._health_focus_param(0, "gl")
    t = win.lines_table.table
    assert t.currentRow() == 0
    assert t.currentColumn() == 2 + win.lines_table._used_keys.index("gl")
    assert win.lines_table.select_param(0, "no_such_param") is False
    assert win.lines_table.select_param(7, "gl") is False


def test_strip_hides_off_the_1d_view_and_resets_on_new_document(qapp, win):
    x, y, result = _fitted(win, qapp)
    win.central_stack.setCurrentWidget(win.view2d)
    assert not win.health_strip.isVisibleTo(win)
    win.central_stack.setCurrentWidget(win.view)
    assert win.health_strip.isVisibleTo(win)
    win._display_1d(x, y, "27Al", 130.3, None, "other", "src2")
    assert win._health is None and win._health_fit is None
    assert win._last_lmfit is None
    assert not win.health_strip.isVisibleTo(win)     # no sites yet
    win.show_correlations()
    assert win.statusBar().currentMessage() == "run a fit first to see correlations"
    win.show_fit_health()
    assert win.statusBar().currentMessage() == "run a fit first (F5)"
    assert not win.actHealth.isEnabled()


def test_workspace_switch_restores_the_verdict_and_covariance(qapp, win):
    x, y, result = _fitted(win, qapp)
    h0, lm0 = win._health, win._last_lmfit
    win._ws_mode = "new"
    x2 = np.linspace(-30, 130, 320)
    win._display_1d(x2, np.exp(-0.5 * ((x2 - 50) / 10) ** 2), "27Al", 130.3,
                    None, "second", "src2")
    assert win.active_ws == 1 and win._health is None
    assert "no fit yet" in win.health_strip.pill.text()
    win.switch_workspace(0)
    qapp.processEvents()
    assert win._health is h0 and win._last_lmfit is lm0
    assert win.health_strip.pill.text() == h0.pill_text()
    assert not win.health_strip.pill.text().endswith(STALE)
    # the re-simulation the switch requests finds equal signatures
    win._sim_done(np.asarray(result.x_ppm), np.asarray(result.y_fit), result.per_site)
    assert win._health is h0


def test_health_menu_reuses_the_decomposition_actions_and_f7_is_bound(qapp, win):
    assert win.actHealth.shortcut().toString() == "F7"
    assert not win.actHealth.isEnabled()             # no verdict yet
    _fitted(win, qapp)
    extra = [win.actCorr, win.actErrors, win.actMC, win.actChi2]
    m = win.health_strip.details_menu(extra)
    acts = m.actions()
    for a in extra:
        assert any(x is a for x in acts), a.text()
    texts = [a.text() for a in acts]
    assert "Show fit report" in texts and "What the flags mean…" in texts
    assert win.actCorr.text().startswith("Parameter correlations")
    assert win.actMC.text().startswith("Monte-&Carlo")
    assert win.actChi2.text().startswith("χ² map")


def test_residual_chip_turns_the_residual_view_on(qapp, win):
    win.actResid.setChecked(False)
    win.view.show_residual = False
    win._health_show_residual()
    assert win.actResid.isChecked() and win.view.show_residual is True
    win._health_show_residual()                      # never toggles it off
    assert win.actResid.isChecked() and win.view.show_residual is True
    win.health_strip.show_residual.emit()
    assert win.actResid.isChecked() and win.view.show_residual is True


def test_panels_toggle_hides_the_strip_and_persists(qapp, win):
    from PySide6.QtCore import QSettings

    s = QSettings("LARMOR", "app")
    saved = s.value("fitHealthStrip")
    try:
        _fitted(win, qapp)
        assert win.health_strip.isVisibleTo(win)
        win.actHealthStrip.setChecked(False)
        assert not win.health_strip.isVisibleTo(win)
        assert s.value("fitHealthStrip", type=bool) is False
        win.actHealthStrip.setChecked(True)
        assert win.health_strip.isVisibleTo(win)
        assert s.value("fitHealthStrip", type=bool) is True
    finally:
        if saved is None:
            s.remove("fitHealthStrip")
        else:
            s.setValue("fitHealthStrip", saved)


def test_theme_switch_restyles_the_strip(qapp, win):
    from larmor.desktop import theme

    prev = theme.active().name
    _fitted(win, qapp)
    text = win.health_strip.pill.text()
    try:
        win._apply_theme_live("Nord")
        assert win.health_strip.colours("ok")[0] == theme.get("Nord").measure
        assert win.health_strip.pill.text() == text
    finally:
        win._apply_theme_live(prev)


def test_auto_fit_result_feeds_the_strip(qapp, win, monkeypatch):
    from PySide6.QtWidgets import QInputDialog

    from larmor import autofit

    x, y, result = _fitted(win, qapp)
    win._health_reset()                              # forget the plain fit's verdict
    assert win._health is None and win._last_lmfit is None
    monkeypatch.setattr(QInputDialog, "getInt",
                        staticmethod(lambda *a, **k: (2, True)))

    def stub(r, ppm, amp, window_ppm=None, n_starts=12, **kw):
        return autofit.AutoFitResult(recipe=r, best_rmsd=result.rmsd,
                                     trials=[result.rmsd], n_improved=0,
                                     result=result)

    monkeypatch.setattr(autofit, "auto_fit", stub)
    win.run_auto_fit()
    assert win.health_strip.pill.text().startswith(("✓ Fit", "⚠ Fit"))
    assert win._last_lmfit is result.lmfit_result
    assert win.lines_table.chi2.text() == win._health.chi_text()
    assert win.actHealth.isEnabled()


@pytest.mark.slow
def test_real_glass_fit_verdict(qapp, win):
    require(CAALGLASS)
    from larmor import fit as fitmod
    from larmor.recipe import Recipe

    win.load_source(str(CAALGLASS), keep_fit=True)
    qapp.processEvents()
    assert win._health is None                       # a loaded fit has no covariance
    result = fitmod.fit(Recipe.from_dict(win.recipe), win.exp_ppm, win.exp_amp,
                        window_ppm=(150.0, -80.0))
    win._fit_done(result)
    qapp.processEvents()
    h = win._health
    assert h is not None and h.level in ("ok", "check", "bad")
    # the Czjzek model lives on its own 16k-point axis: the residual is judged
    # after interpolation onto the data (the legacy path silently skipped it)
    assert h.noise_ratio is not None and h.runs_z is not None
    texts = [f.text for f in h.flags]
    assert [c.text() for c in win.health_strip.chips] == \
        [t if len(t) <= 48 else t[:47] + "…" for t in texts]
    assert win.qtable.rowCount() == 5                # quantify still ran
    # three overlapping Czjzek sites on one 1D lineshape: the degenerate pair
    # is the documented outcome (docs/tutorials/01, section 4)
    assert "degenerate" in h.kinds() and h.level == "bad"
    assert win._health_fit is h and win._last_lmfit is result.lmfit_result
