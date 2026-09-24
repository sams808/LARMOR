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

from conftest import CAALGLASS, LAW_CA_11B, require

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


def _two_site_fitted(win, qapp):
    """A two-line 11B synthetic (15 ppm / 1 ppm) on the workbench, fitted,
    fed through _fit_done -- the shape a BO3/BO4 tagging needs."""
    from larmor import fit as fitmod
    from larmor.engine import make_context, simulate_site
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 40, 600)
    truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="B3", params={
            "isotropic_chemical_shift_ppm": Param(15.0, min=0, max=30),
            "shift_fwhm_ppm": Param(8.0, min=1, max=40),
            "gl": Param(1.0, vary=False), "amplitude": Param(100.0, min=0)}),
        SiteModel(model="gauss_lor", label="B4", params={
            "isotropic_chemical_shift_ppm": Param(1.0, min=-10, max=10),
            "shift_fwhm_ppm": Param(4.0, min=1, max=40),
            "gl": Param(1.0, vary=False), "amplitude": Param(60.0, min=0)})])
    ctx = make_context(truth, exp_ppm=x)
    y = np.sum([simulate_site(s, ctx) for s in truth.sites], axis=0)
    y = y + np.random.RandomState(1).normal(0, 0.5, y.size)
    win._display_1d(x, y, "11B", 160.0, None, "t", "src")
    win.recipe["sites"] = truth.to_dict()["sites"]
    result = fitmod.fit(Recipe.from_dict(win.recipe), x, y, window_ppm=(40, -20))
    win._fit_done(result)
    qapp.processEvents()
    return x, y, result


def test_report_family_rows_appear_after_tagging_and_stay_absent_when_untagged(qapp, win):
    """N3: the Report (F6) adds bold Σ family rows and the N4 row only once
    lines are tagged, on the covariance basis right after a fit (the
    tooltip says so), falls back to the flagged independent basis after a
    value edit, and Copy CSV / Copy methods carry the block."""
    _two_site_fitted(win, qapp)
    assert win.qtable.rowCount() == 2
    assert win._last_quant["families"] == [] and win._last_quant["ratios"] == []
    win.recipe["sites"][0]["family"] = "BO3"
    win.recipe["sites"][1]["family"] = "BO4"
    win.on_family_changed()
    assert win.qtable.rowCount() == 2 + 2 + 1
    assert win.qtable.item(2, 0).text().startswith("Σ BO3  (A)")
    assert win.qtable.item(3, 0).text().startswith("Σ BO4  (B)")
    assert win.qtable.item(4, 0).text() == "N4 = BO4/(BO3+BO4)"
    assert "±" in win.qtable.item(4, 3).text() and "±" in win.qtable.item(2, 3).text()
    assert win.qtable.item(2, 0).font().bold() and not win.qtable.item(0, 0).font().bold()
    assert win._last_quant["family_basis"] == "covariance"
    assert "covariance" in win.qtable.item(4, 0).toolTip()
    n4 = win._last_quant["ratios"][0]
    assert n4["defined"] and 0.0 < n4["value"] < 1.0 and n4["err"] > 0
    assert "Report (F6) updated" in win.statusBar().currentMessage()
    # the per-site rows are the same numbers as before the tagging
    assert win.qtable.item(0, 3).text() and win.qtable.item(1, 3).text()
    # Copy CSV: the site table, a blank line, the family block, the basis line
    win.copy_csv()
    text = QApplication.clipboard().text()
    assert "\n\nfamily,lines,integral,fraction_pct,fraction_err_pct\n" in text
    assert "\nBO3,A," in text and "\nBO4,B," in text
    assert "ratio,description,value,err\nN4,BO4/(BO3+BO4)," in text
    assert text.rstrip().endswith("# family/ratio errors: covariance")
    win.copy_methods()
    assert "structural families (BO3, BO4)" in QApplication.clipboard().text()
    assert "covariance between line amplitudes" in QApplication.clipboard().text()
    # an amplitude edited after the fit: the covariance is stale -> flagged
    win.recipe["sites"][0]["params"]["amplitude"]["value"] *= 1.01
    win.on_family_changed()
    assert win.qtable.rowCount() == 5
    assert win._last_quant["family_basis"] == "independent"
    assert "independent" in win.qtable.item(4, 0).toolTip()
    win.copy_methods()
    assert "treated as independent" in QApplication.clipboard().text()
    # untagging both lines removes the block again
    win.recipe["sites"][0].pop("family")
    win.recipe["sites"][1].pop("family")
    win.on_family_changed()
    assert win.qtable.rowCount() == 2
    win.copy_csv()
    assert "family" not in QApplication.clipboard().text()


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
    assert any(t.startswith("residual ") and "× noise" in t for t in texts), texts


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
    # is the documented outcome (docs/tutorials/01, section 4) -- amber, a
    # thing to check, since every value is physical
    assert "degenerate" in h.kinds() and h.level == "check"
    assert "physical" not in h.kinds()
    assert win._health_fit is h and win._last_lmfit is result.lmfit_result


# ------------------------------------------------------ quantitativity chips
def _chips(win) -> dict:
    return {c.property("flag_kind"): c for c in win.health_strip.chips}


def test_tail_chip_widens_the_window_and_requantifies(qapp, win):
    x, y, result = _fitted(win, qapp)
    sigma = 8.0 / 2.3548
    assert win.qtable.columnCount() == 5
    assert win.qtable.horizontalHeaderItem(4).text() == "outside window (%)"
    # quantify integrates the VIEW; offscreen the auto-range has not been
    # applied yet, so set the fit window explicitly: no tail chip there
    from larmor.quantify import TAIL_LIMIT_PCT

    win.view.setXRange(-20.0, 120.0, padding=0)
    win.run_quantify(show=False)
    win._health_from_result(result)
    assert "tail" not in win._health.kinds() and win._health.tail_checked
    # the fitted mix is part Lorentzian: its wings leave ~1.7 % outside
    assert 0.0 <= float(win.qtable.item(0, 4).text()) < TAIL_LIMIT_PCT
    win.view.setXRange(57.0, 63.0, padding=0)
    win.run_quantify(show=False)
    win._health_from_result(result)
    qapp.processEvents()
    chip = _chips(win)["tail"]
    assert chip.text().startswith("tail outside window: p ")
    assert float(win.qtable.item(0, 4).text()) > 2
    assert "⚠ tail outside window" in win.results_summary.text()
    chip.click()
    qapp.processEvents()
    # the fitted mix is half Lorentzian, whose 99.5 % window (~±250 ppm) is
    # clipped to the acquired -20 … 120 ppm axis: the view spans the whole
    # spectrum, ~1.7 % stays beyond it (under the 2 % limit) and is named
    (x0, x1), _ = win.view.getPlotItem().getViewBox().viewRange()
    assert max(x0, x1) >= 60 + 2.8 * sigma - 0.6 and min(x0, x1) <= 60 - 2.8 * sigma + 0.6
    assert max(x0, x1) == pytest.approx(120.0, abs=0.5)
    assert min(x0, x1) == pytest.approx(-20.0, abs=0.5)
    outside = float(win.qtable.item(0, 4).text())
    assert 0.0 <= outside < TAIL_LIMIT_PCT
    assert "tail" not in win._health.kinds()
    assert win._health is win._health_fit and not win._health.stale
    msg = win.statusBar().currentMessage()
    assert "window widened to 120 … -20 ppm" in msg
    assert "of p lies beyond the acquired spectrum" in msg
    assert "tails inside the window" in win.health_strip.pill.toolTip()
    # a pure Gaussian (gl fixed) is contained: the widened window leaves ≤ 0.5 %
    win.recipe["sites"][0]["params"]["gl"]["value"] = 1.0
    win.recipe["sites"][0]["params"]["gl"]["vary"] = False
    win.view.setXRange(57.0, 63.0, padding=0)
    win.run_quantify(show=False)
    assert float(win.qtable.item(0, 4).text()) > 30
    win._health_widen_window()
    (x0, x1), _ = win.view.getPlotItem().getViewBox().viewRange()
    assert max(x0, x1) == pytest.approx(60 + 2.807 * sigma, abs=0.6)
    assert min(x0, x1) == pytest.approx(60 - 2.807 * sigma, abs=0.6)
    assert float(win.qtable.item(0, 4).text()) <= 0.5
    assert "beyond the acquired spectrum" not in win.statusBar().currentMessage()
    win.copy_csv()
    head = QApplication.clipboard().text().splitlines()[0]
    assert head.endswith(",tail_outside_pct")


def test_acquisition_chips_from_stubbed_facts_and_click_throughs(qapp, win, monkeypatch):
    from larmor import quantitativity as Q
    from larmor.desktop import dialogs, satrec_dialog
    from larmor.satrec import T1Region

    acq = Q.Acquisition(expno="/x/2704", nucleus="27Al", pulprog="zg", kind="Single pulse",
                        ns=512, d1_s=1.0, aq_s=0.02, p1_us=0.375, plw1_w=207.0,
                        probhd="SPRB600511_7297 (MAS)", title="P1(90)=3.125; 11deg tip",
                        p90_us_title=3.125, flip_deg_title=11.0, t1_multiple_claimed=None)
    t1 = Q.T1Source(expno="/x/2701", kind="ct1t2",
                    regions=[T1Region(1, 120.0, -40.0, 0.557, (0.557,), None, 10)])
    facts = Q.AcqFacts(acquisition=acq, t1=t1, sibling_expno="/x/2701", spin=2.5,
                       folder="/x", t1_status="ok")
    monkeypatch.setattr(Q, "facts_for", lambda src: facts)
    _fitted(win, qapp)
    chips = _chips(win)
    assert chips["recovery"].text().startswith("D1 = 1.8 T1 at 10.8° → ")
    assert chips["excitation"].text() == "flip 10.8° > 10° limit (I = 5/2)"
    assert "just above the limit (8 %)" in chips["excitation"].toolTip()
    assert "⚠ D1 = 1.8 T1" in win.results_summary.text()
    assert win._acq_cache == {"src": facts}
    # the recycle chip opens the Relaxation tool on the sibling EXPNO
    opened = []
    monkeypatch.setattr(satrec_dialog.SatrecDialog, "__init__",
                        lambda self, parent, expno: opened.append((parent, expno)))
    monkeypatch.setattr(satrec_dialog.SatrecDialog, "exec", lambda self: 0)
    win.health_strip.open_relaxation.emit()
    assert opened and opened[-1][0] is win and str(opened[-1][1]).endswith("2701")
    # the F7 menu offers the per-site measurement and the tool on that EXPNO
    win.health_strip.show_details = lambda extra: setattr(win, "_f7_extra", list(extra))
    win.show_fit_health()
    texts = [a.text() for a in win._f7_extra]
    assert "Measure T1 per site from EXPNO 2701 (uses this fit)…" in texts
    assert "Open T1 measurement (EXPNO 2701)…" in texts
    # the flip chip opens Experiment parameters on the 90° pulse; a typed
    # flip angle re-judges the last fit without a refit
    seen = {}

    def fake_exec(dlg):
        seen["focus"] = dlg.p90.hasFocus() or dlg.focusWidget() is dlg.p90
        seen["p90"] = dlg.p90.value()
        seen["rows"] = dlg.t1 is not None
        dlg.flip.setValue(8.0)
        dlg._accept()
        return 1

    monkeypatch.setattr(dialogs.ExperimentDialog, "exec", fake_exec)
    win.health_strip.enter_flip.emit()
    qapp.processEvents()
    assert seen["p90"] == 3.125 and seen["rows"]           # prefilled from P1(90)=
    assert win.recipe["provenance"]["quantitativity"]["flip_deg"] == 8.0
    chips = _chips(win)
    assert "excitation" not in chips
    assert chips["recovery"].text().startswith("D1 = 1.8 T1 at 8° → ")
    assert not chips["recovery"].text().endswith("(90° assumed)")
    assert win._health is win._health_fit and not win._health.stale
    win._health_reset()
    assert win._acq_cache == {}


def test_real_base1ca_recovery_chip(qapp, win):
    path = require(LAW_CA_11B[1])
    from larmor import fit as fitmod
    from larmor.recipe import Param, Recipe, SiteModel

    # the 1r itself: the EXPNO has two procnos and load_source would ask
    win.load_source(str(path / "pdata" / "1" / "1r"), keep_fit=False)
    qapp.processEvents()
    assert win.recipe.get("source_path", "").endswith("24")
    sites = [SiteModel(model="gauss_lor", label=lab, params={
        "isotropic_chemical_shift_ppm": Param(pos, min=pos - 6, max=pos + 6),
        "shift_fwhm_ppm": Param(8.0, min=2, max=30), "gl": Param(0.5, vary=False),
        "amplitude": Param(float(win.exp_amp.max()), min=0)})
        for lab, pos in (("BO3", 15.0), ("BO4", 1.0))]
    win.recipe["sites"] = Recipe(nucleus="11B", larmor_frequency_MHz=192.43,
                                 sites=sites).to_dict()["sites"]
    result = fitmod.fit(Recipe.from_dict(win.recipe), win.exp_ppm, win.exp_amp,
                        window_ppm=(40.0, -20.0))
    win._fit_done(result)
    qapp.processEvents()
    h = win._health
    (rec,) = [f for f in h.flags if f.kind == "recovery"]
    assert rec.level == "check" and rec.text.endswith("(90° assumed)")
    assert rec.text.startswith("D1 = 3.0–3.6 T1 → 95–97 %")
    assert "EXPNO 23" in rec.detail and "4.64" in rec.detail
    # the unknown flip angle is no chip: a tooltip line names the fact
    assert "excitation" not in h.kinds()
    assert any(u.startswith("flip angle unknown (I = 3/2)") for u in h.unchecked)
    assert "flip angle unknown (I = 3/2)" in win.health_strip.pill.toolTip()
    assert h.acquisition.facts.t1.expno.endswith("23")
    assert h.acquisition.facts.acquisition.d1_s == 14.0
    # the 90° pulse typed once clears both chips without a refit
    win.recipe["provenance"] = {"quantitativity": {"p90_us": 3.5}}
    win._health_requantify()
    h = win._health
    assert not ({"recovery", "excitation"} & {f.kind for f in h.flags if f.level == "check"})
    assert "excitation" not in h.kinds()
    passing = h.passing()
    assert any(p.startswith("recycle 14 s = 3.0–3.6 T1 → ") for p in passing), passing
    assert "flip 10.9° within the linear regime (≤ 15°)" in passing
