"""The error tools on a fit with ONE line, and what they must say.

Sami could not "calculate error on a simple fit of a czjzek for 23Na". Run
on a real single-line 23Na Czjzek fit, every tool returned -- but the χ²
profile dialog sat frozen for the ~15 s its process pool took to start (the
whole scan was 2 s of work), the profile's best value disagreed with the
table because the scan found a lower χ² than the fit had reached and said
nothing, and the Report gave the single line a population of 100 ± 18 %.
"""
import time
import warnings

import numpy as np
import pytest

from larmor import autofit, parallel, quantify
from larmor.engine import make_context, simulate_site
from larmor.recipe import Param, Recipe, SiteModel


def _one_line(noise=0.01, seed=0):
    x = np.linspace(-40.0, 60.0, 400)
    truth = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(12.0, min=-40, max=60),
            "shift_fwhm_ppm": Param(6.0, min=0.5, max=40),
            "gl": Param(0.5, min=0, max=1, vary=False),
            "amplitude": Param(100.0, min=0)})])
    ctx = make_context(truth, exp_ppm=x)
    y = np.sum([simulate_site(s, ctx) for s in truth.sites], axis=0)
    y = y + np.random.default_rng(seed).normal(0.0, noise * 100.0, x.size)
    return truth, x, y


# ------------------------------------------------------------- populations
def test_single_line_population_has_no_error():
    from larmor import fit as fitmod

    truth, x, y = _one_line()
    r = Recipe.from_dict(truth.to_dict())
    fitmod.fit(r, x, y)
    q = quantify.quantify(r, window_ppm=(60.0, -40.0))
    assert len(q["rows"]) == 1
    assert q["rows"][0]["fraction_pct"] == 100.0
    assert q["rows"][0]["fraction_err_pct"] == 0.0
    assert q["rows"][0]["integral_err"] is not None and q["rows"][0]["integral_err"] > 0


def test_fraction_error_propagates_the_normalisation():
    """f_i = I_i / ΣI: the derivative goes through the total. Checked against
    a numerical derivative and against the old σ_i / T shortcut, which is
    only right when the OTHER integrals carry no error."""
    ints, errs = [60.0, 40.0], [3.0, 2.0]

    def frac(i, vals):
        return 100.0 * vals[i] / sum(vals)

    for i in range(2):
        analytic = quantify.fraction_err_pct(ints, errs, i)
        h = 1e-6
        var = 0.0
        for j in range(2):
            up = list(ints); up[j] += h
            dn = list(ints); dn[j] -= h
            d = (frac(i, up) - frac(i, dn)) / (2 * h)      # percentage points per unit
            var += (d * errs[j]) ** 2
        assert analytic == pytest.approx(np.sqrt(var), rel=1e-6)
        # the old shortcut (σ_i / T) is the same only when the other line has
        # no error; with both uncertain it understates one and overstates none
        shortcut = 100.0 * errs[i] / sum(ints)
        assert analytic != pytest.approx(shortcut, rel=1e-3)
    # one line: 100 % exactly
    assert quantify.fraction_err_pct([60.0], [3.0], 0) == 0.0
    # an unknown error for the line itself: unknown; unknown errors of the
    # OTHER lines contribute nothing
    assert quantify.fraction_err_pct([60.0, 40.0], [None, 2.0], 0) is None
    only_own = quantify.fraction_err_pct([60.0, 40.0], [3.0, None], 0)
    assert only_own == pytest.approx(100.0 * (40.0 / 100.0 ** 2) * 3.0)


# ------------------------------------------------------------- χ² profile
def test_profile_says_when_the_fit_was_not_at_the_minimum():
    from larmor import fit as fitmod

    truth, x, y = _one_line()
    r = Recipe.from_dict(truth.to_dict())
    fitmod.fit(r, x, y)
    at_min = autofit.error_profile(r, x, y, site=0, param="shift_fwhm_ppm",
                                   n_points=9, span=3.0, parallel=False)
    assert at_min.fit_at_minimum and at_min.chi2_fit is not None
    assert not any("had not converged" in n for n in at_min.notes)

    # the same recipe pushed off its minimum: the scan refits the rest and
    # lands lower than the recipe's own χ² -- it must say so
    off = Recipe.from_dict(r.to_dict())
    off.sites[0].params["isotropic_chemical_shift_ppm"].value += 3.0
    off.sites[0].params["isotropic_chemical_shift_ppm"].stderr = 0.2
    prof = autofit.error_profile(off, x, y, site=0, param="shift_fwhm_ppm",
                                 n_points=9, span=3.0, parallel=False)
    assert prof.chi2_fit > prof.level68
    assert prof.fit_at_minimum is False
    assert any("had not converged" in n for n in prof.notes)
    assert prof.chi2_min < prof.chi2_fit


def test_profile_flags_an_interval_narrower_than_its_step():
    from larmor import fit as fitmod

    truth, x, y = _one_line()
    r = Recipe.from_dict(truth.to_dict())
    fitmod.fit(r, x, y)
    # a span of 30 error bars: the 1σ crossings fall inside the first step
    wide = autofit.error_profile(r, x, y, site=0, param="isotropic_chemical_shift_ppm",
                                 n_points=9, span=30.0, parallel=False)
    assert any("narrower than the scan step" in n for n in wide.notes)
    fine = autofit.error_profile(r, x, y, site=0, param="isotropic_chemical_shift_ppm",
                                 n_points=9, span=3.0, parallel=False)
    assert not any("narrower than the scan step" in n for n in fine.notes)


# ------------------------------------------------------------- the solver
def test_solver_scaling_reaches_the_deeper_minimum_of_the_shipped_27al_fit():
    """The shipped 27Al recipe (three Czjzek sites) sat at RMSD 0.0468, a
    shallow minimum the unscaled trust region could not leave because the
    amplitude (~4e6) dominated its step test; with x_scale="jac" the plain
    fit reaches 0.0458 -- the minimum a 12-start Auto fit also finds -- and
    a refit from there stays put."""
    from pathlib import Path

    from larmor import fit as fitmod
    from larmor import loader

    assert fitmod.X_SCALE == "jac"
    path = Path(__file__).resolve().parents[1] / "examples" / "pCABS2-4_27Al.recipe.json"
    ppm, amp, rec, meta, warns = loader.load_any(str(path))
    win = tuple(rec.get("fit_window_ppm") or (150.0, -80.0))
    r = Recipe.from_dict(rec)
    res = fitmod.fit(r, ppm, amp, window_ppm=win)
    assert res.rmsd < 0.0462
    assert not res.lmfit_result.aborted
    again = Recipe.from_dict(r.to_dict())
    res2 = fitmod.fit(again, ppm, amp, window_ppm=win)
    assert res2.rmsd == pytest.approx(res.rmsd, abs=2e-5)
    pos = again.sites[0].params["isotropic_chemical_shift_ppm"].value
    assert pos == pytest.approx(r.sites[0].params["isotropic_chemical_shift_ppm"].value,
                                abs=0.2)


# ------------------------------------------------------------- the pool
def test_auto_parallel_times_the_first_item_and_stays_sequential_when_short(monkeypatch):
    calls = []

    def fn(item):
        calls.append(item)
        return item * 2

    seen = []
    res, ran = autofit._map_items(fn, list(range(12)), parallel="auto",
                                  on_result=lambda i, r: seen.append((i, r)))
    assert res == [i * 2 for i in range(12)] and ran is False
    assert seen == [(i, i * 2) for i in range(12)]        # in order, first included

    # a slow first item with many to go: the pool is chosen (patched to a
    # sequential stand-in so the test never spawns processes)
    used = {}

    def fake_map(f, items, **kw):
        used["use"] = kw.get("use_processes")
        out = []
        for k, it in enumerate(items):
            r = f(it)
            out.append(r)
            if kw.get("on_result"):
                kw["on_result"](k, r)
        return out

    monkeypatch.setattr(parallel, "parallel_map", fake_map)
    monkeypatch.setattr(autofit, "POOL_BREAKEVEN_S", 0.05)

    def slow(item):
        time.sleep(0.02)
        return item

    res2, ran2 = autofit._map_items(slow, list(range(12)), parallel="auto")
    assert res2 == list(range(12))
    assert used["use"] is True and ran2 is True
    monkeypatch.setattr(autofit, "POOL_BREAKEVEN_S", 1e9)
    res3, ran3 = autofit._map_items(slow, list(range(12)), parallel="auto")
    assert res3 == list(range(12)) and used["use"] is False and ran3 is False


def test_parallel_map_falls_back_to_sequential_when_the_pool_breaks(monkeypatch):
    from concurrent.futures.process import BrokenProcessPool

    class BrokenPool:
        def submit(self, fn, item):
            raise BrokenProcessPool("a worker died")

    monkeypatch.setattr(parallel, "shared_pool", lambda n=None: BrokenPool())
    shut = []
    monkeypatch.setattr(parallel, "shutdown_shared_pool", lambda: shut.append(1))
    beats = []
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        res = parallel.parallel_map(lambda v: v + 1, list(range(20)),
                                    use_processes=True, heartbeat=lambda: beats.append(1))
    assert res == list(range(1, 21))
    assert shut == [1]
    assert any("process pool failed" in str(x.message) for x in w)
    assert len(beats) >= 20                 # the sequential path pulses too


def test_heartbeat_pulses_during_a_real_pool_start():
    """Through a real (tiny) pool: the heartbeat fires while nothing has
    finished yet, so a dialog can pump its event loop."""
    beats = []
    res = parallel.parallel_map(_double, list(range(10)), max_workers=2,
                                use_processes=True, heartbeat=lambda: beats.append(1))
    assert res == [2 * i for i in range(10)]
    assert len(beats) >= 1


def _double(v):
    return 2 * v


# ------------------------------------------------------------- Monte-Carlo
def test_monte_carlo_names_parameters_the_refit_moved():
    from larmor import fit as fitmod

    truth, x, y = _one_line()
    r = Recipe.from_dict(truth.to_dict())
    fitmod.fit(r, x, y)
    ok = autofit.monte_carlo_errors(Recipe.from_dict(r.to_dict()), x, y,
                                    window_ppm=(60.0, -40.0), n_trials=10,
                                    parallel=False)
    assert ok.moved == [] and "had not converged" not in ok.summary
    off = Recipe.from_dict(r.to_dict())
    off.sites[0].params["isotropic_chemical_shift_ppm"].value += 3.0
    mc = autofit.monte_carlo_errors(off, x, y, window_ppm=(60.0, -40.0),
                                    n_trials=10, parallel=False)
    assert "s0.isotropic_chemical_shift_ppm" in mc.moved
    assert "had not converged" in mc.summary


# ------------------------------------------------------------- the dialog
pytest.importorskip("PySide6")


@pytest.fixture(scope="module")
def qapp():
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


def test_errors_dialog_scans_with_progress_and_notes(qapp):
    from larmor import fit as fitmod
    from larmor.desktop.tool_dialogs import ErrorsDialog

    truth, x, y = _one_line()
    r = Recipe.from_dict(truth.to_dict())
    fitmod.fit(r, x, y)
    dlg = ErrorsDialog(None, r.to_dict(), x, y, (60.0, -40.0))
    try:
        assert dlg.site.count() == 1
        assert [dlg.param.itemText(i) for i in range(dlg.param.count())] == [
            "isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "amplitude"]
        assert dlg.btnStop.isEnabled() is False
        dlg.points.setValue(7)
        dlg.param.setCurrentText("shift_fwhm_ppm")
        dlg._run()
        assert dlg.prof is not None and dlg.prof.fit_at_minimum
        assert dlg.prog.value() == 7 and dlg.prog.maximum() == 7
        assert "s0.shift_fwhm_ppm" in dlg.res.text()
        assert dlg.btnRun.isEnabled() and not dlg.btnStop.isEnabled()
        # a stop before the third point (the Stop button sets the flag while
        # the scan runs): the dialog says so instead of calling it a failure
        def stop_now(*a, **k):
            dlg._stop = True
            raise RuntimeError("chi-square profile failed: too few valid points")

        import larmor.autofit as af
        real = af.error_profile
        af.error_profile = stop_now
        try:
            dlg._run()
        finally:
            af.error_profile = real
        assert dlg.res.text().startswith("stopped")
    finally:
        dlg.close()
