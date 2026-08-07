"""Regressions for Monte-Carlo error estimation (autofit.monte_carlo_errors).

The MC method (dmfit ▸ Errors ▸ Monte Carlo) must: recover the injected noise
level, return an unbiased mean per parameter, give a positive, sensible σ, fire
progress, and honour an early stop. If this drifts, quoted MC errors are wrong.
"""
import numpy as np
import pytest

from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine, autofit


def _synthetic(noise=2.0, seed=1):
    x = np.linspace(-20, 60, 1500)                 # ascending, as the app stores
    truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                   sites=[SiteModel(model="gauss_lor", label="A", params={
                       "isotropic_chemical_shift_ppm": Param(15.0),
                       "shift_fwhm_ppm": Param(6.0, min=0.1),
                       "amplitude": Param(100.0, min=0.0),
                       "gl": Param(1.0, vary=False)})])
    _, model, _ = engine.simulate(truth, exp_ppm=x)
    data = model + np.random.default_rng(seed).normal(0.0, noise, x.size)
    return truth, x, data


def test_recovers_noise_and_unbiased_errors():
    truth, x, data = _synthetic(noise=2.0)
    fired = []
    res = autofit.monte_carlo_errors(
        Recipe.from_dict(truth.to_dict()), x, data, window_ppm=(-10, 40),
        n_trials=40, seed=3, progress=lambda k, n: fired.append(k))
    assert res.n_ok == 40
    assert len(fired) == 40
    assert res.noise == pytest.approx(2.0, rel=0.25)      # recovers injected σ
    by = {p.param: p for p in res.params}
    # mean is unbiased (close to the best fit) and σ is positive & small for a
    # strong, well-determined line
    for name, best in (("isotropic_chemical_shift_ppm", 15.0),
                       ("shift_fwhm_ppm", 6.0), ("amplitude", 100.0)):
        p = by[name]
        assert p.std > 0
        assert abs(p.mean - best) < max(5 * p.std, 0.5)
    assert by["isotropic_chemical_shift_ppm"].std < 0.5    # ppm, tightly fit
    assert np.isfinite(by["amplitude"].pct)


def test_manual_noise_and_early_stop():
    truth, x, data = _synthetic()
    # a should_stop that trips after 5 trials
    state = {"k": 0}

    def stop():
        state["k"] += 1
        return state["k"] > 5

    res = autofit.monte_carlo_errors(
        Recipe.from_dict(truth.to_dict()), x, data, window_ppm=(-10, 40),
        n_trials=100, seed=0, noise=3.0, should_stop=stop)
    assert res.noise == 3.0                # honoured the manual noise
    assert res.n_ok <= 6                   # stopped early
    assert res.trials == 100


def test_kernel_model_does_not_crash_on_grid_mismatch():
    """A Czjzek site simulates on the KERNEL's own fixed grid, never on exp_ppm
    (engine.needs_kernel) -- monte_carlo_errors must interpolate that grid onto
    exp_ppm before comparing to exp_amp, or every kernel-model MC run crashes
    with a shape mismatch (was: IndexError, 2048 vs len(exp_ppm))."""
    x = np.linspace(-40, 120, 700)                 # deliberately NOT the kernel size
    truth = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="czjzek", label="Al", params={
            "isotropic_chemical_shift_ppm": Param(60.0),
            "sigma_Cq_MHz": Param(3.0, min=0.05),
            "shift_fwhm_ppm": Param(8.0, min=0.1),
            "line_fwhm_ppm": Param(1.0, min=0.0),
            "amplitude": Param(100.0, min=0.0)})])
    kx, model, _ = engine.simulate(truth)           # the kernel's own grid
    y = np.interp(x, np.sort(kx), model[np.argsort(kx)])
    data = y + np.random.default_rng(2).normal(0.0, y.max() * 0.01, x.size)

    res = autofit.monte_carlo_errors(
        Recipe.from_dict(truth.to_dict()), x, data, window_ppm=(120, -40),
        n_trials=12, seed=1)
    assert res.n_ok == 12
    amp = next(p for p in res.params if p.param == "amplitude")
    assert amp.std > 0 and np.isfinite(amp.std)
    assert abs(amp.mean - 100.0) < max(10 * amp.std, 20.0)


def test_parallel_matches_sequential_for_the_same_seed():
    """parallel=True must be a pure execution-strategy change, not a
    different answer: the synthetic noise draws are generated up front in a
    single sequential pass over `rng` regardless of how the trials
    themselves get scheduled, so the SAME seed gives BIT-IDENTICAL results
    whether trials run one at a time or across a process pool."""
    truth, x, data = _synthetic()
    seq = autofit.monte_carlo_errors(
        Recipe.from_dict(truth.to_dict()), x, data, window_ppm=(-10, 40),
        n_trials=12, seed=7, parallel=False)
    par = autofit.monte_carlo_errors(
        Recipe.from_dict(truth.to_dict()), x, data, window_ppm=(-10, 40),
        n_trials=12, seed=7, parallel=True, max_workers=2)
    assert par.n_ok == seq.n_ok == 12
    by_seq = {p.label: p for p in seq.params}
    by_par = {p.label: p for p in par.params}
    assert set(by_seq) == set(by_par)
    for label, ps in by_seq.items():
        pp = by_par[label]
        assert pp.mean == pytest.approx(ps.mean, rel=1e-9)
        assert pp.std == pytest.approx(ps.std, rel=1e-9)


def test_report_lists_every_free_parameter():
    truth, x, data = _synthetic()
    res = autofit.monte_carlo_errors(
        Recipe.from_dict(truth.to_dict()), x, data, window_ppm=(-10, 40),
        n_trials=15, seed=1)
    text = res.report()
    assert "Monte-Carlo errors" in text
    for name in ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "amplitude"):
        assert f"s0.{name}" in text
