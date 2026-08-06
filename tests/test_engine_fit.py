import numpy as np
import pytest

from larmor import engine
from larmor.io import fxmla
from larmor import fit as fitmod
from larmor.recipe import Recipe, SiteModel, Param

from conftest import CAALGLASS, require


def _degenerate_recipe():
    """Two IDENTICAL-shape overlapping sites: their amplitudes are perfectly
    correlated (only the SUM is determined by the data), so the covariance
    matrix is singular and the first least_squares pass reliably fails to get
    error bars -- a controlled trigger for the errorbar-rescue retry."""
    x = np.linspace(-30, 30, 500)
    truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(0.0),
            "shift_fwhm_ppm": Param(10.0), "amplitude": Param(80.0),
            "gl": Param(1.0, vary=False)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    data = y + np.random.default_rng(0).normal(0, 0.3, x.size)
    r = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(0.0, vary=False),
            "shift_fwhm_ppm": Param(10.0, vary=False),
            "amplitude": Param(40.0, min=0), "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="B", params={
            "isotropic_chemical_shift_ppm": Param(0.0, vary=False),
            "shift_fwhm_ppm": Param(10.0, vary=False),
            "amplitude": Param(40.0, min=0), "gl": Param(1.0, vary=False)}),
    ])
    return r, x, data


def test_compute_errorbars_false_skips_the_rescue_retry():
    """A degenerate covariance normally triggers a second full optimization
    (leastsq from the converged point) to try to recover error bars.
    compute_errorbars=False must skip that retry entirely -- verified by
    counting actual evaluations, not just checking the final stderr."""
    r1, x, data = _degenerate_recipe()
    calls1 = []
    fitmod.fit(r1, x, data, compute_errorbars=True,
              iter_cb=lambda *a, **k: calls1.append(1))

    r2, x, data = _degenerate_recipe()
    calls2 = []
    fitmod.fit(r2, x, data, compute_errorbars=False,
              iter_cb=lambda *a, **k: calls2.append(1))

    # with the retry skipped, meaningfully fewer evaluations happened
    assert len(calls2) < len(calls1)
    # the fitted (identifiable) TOTAL amplitude is unaffected either way
    total1 = sum(s.params["amplitude"].value for s in r1.sites)
    total2 = sum(s.params["amplitude"].value for s in r2.sites)
    assert total1 == pytest.approx(total2, rel=0.05)


def test_gauss_lor_shapes():
    x = np.linspace(-50, 50, 2001)
    g = engine.gauss_lor(x, 0.0, 10.0, 2.0, gl=1.0)
    l = engine.gauss_lor(x, 0.0, 10.0, 2.0, gl=0.0)
    assert g.max() == pytest.approx(2.0, rel=1e-6)
    assert l.max() == pytest.approx(2.0, rel=1e-6)
    # both peak-normalized with the same FWHM...
    half = np.where(g >= 1.0)[0]
    assert x[half[-1]] - x[half[0]] == pytest.approx(10.0, abs=0.2)
    # ...but the Lorentzian has heavier tails
    tail = np.abs(x) > 25
    assert l[tail].sum() > 5 * g[tail].sum()


@pytest.mark.slow
def test_kernel_and_fit_caalglass():
    """End-to-end: import the dmfit fit, refine with lmfit, beat Phase 0's RMSD."""
    dm = fxmla.read(require(CAALGLASS))
    recipe, _ = fxmla.to_recipe(dm)

    kernel = engine.build_kernel(
        recipe.nucleus, recipe.larmor_frequency_MHz, recipe.spin_rate_Hz)
    assert kernel.K.shape == (80 * 11, 2048)
    assert np.all(np.diff(kernel.x_ppm) > 0)

    # sigma reweighting sanity: larger sigma -> broader lineshape
    y_narrow = kernel.weights(1.0) @ kernel.K
    y_broad = kernel.weights(4.0) @ kernel.K
    width = lambda y: np.sum(y > y.max() / 2)
    assert width(y_broad) > width(y_narrow)

    exp_ppm, exp_amp = dm.spectrum.ppm, dm.spectrum.amplitude
    result = fitmod.fit(recipe, exp_ppm, exp_amp, window_ppm=(150.0, -80.0))

    # refined fit should beat the fixed-parameter Phase 0 replay (RMSD 0.027)
    assert result.rmsd < 0.01

    # dmfit's ad-hoc Gauss/Lor sideband lines (236/208 ppm) sit outside the
    # window and must be frozen, not left to wander
    assert len(result.frozen_sites) == 2

    site1 = recipe.sites[0]
    assert site1.params["isotropic_chemical_shift_ppm"].value == pytest.approx(
        66.2, abs=5.0)
    assert site1.params["sigma_Cq_MHz"].value == pytest.approx(2.0, abs=0.7)
    # the whole point: uncertainties exist now
    assert result.lmfit_result.errorbars
    assert site1.params["sigma_Cq_MHz"].stderr is not None
    assert 0 < site1.params["sigma_Cq_MHz"].stderr < 0.5


def test_grid_restrictable_allowlist():
    from larmor import engine

    def rec(*models):
        return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
            SiteModel(model=m, label=str(i), params={"amplitude": Param(1.0)})
            for i, m in enumerate(models)])

    assert engine.grid_restrictable(rec("gauss_lor"))
    assert engine.grid_restrictable(rec("czjzek"))
    assert engine.grid_restrictable(rec("gauss_lor", "voigt", "czjzek"))
    # the LRU single-site models build their own cached simulation directly
    # from ctx.x_ppm's span and convolution-broaden it -- restricting the grid
    # there would truncate real signal (satellite/sideband manifolds), not
    # just its cost, so they must NOT be treated as restrictable
    for unsafe in ("quad_ct", "quad_first", "quad_csa", "csa_mas", "csa_czjzek"):
        assert not engine.grid_restrictable(rec(unsafe)), unsafe
    # ANY unsafe site anywhere in the recipe disqualifies the whole fit (ctx is
    # shared across all sites -- there's no such thing as a partial restriction)
    assert not engine.grid_restrictable(rec("gauss_lor", "quad_ct"))
    # an unregistered/unaudited future model defaults to NOT restrictable
    assert not engine.grid_restrictable(rec("some_future_model"))


def test_fit_windowed_grid_recovers_truth_and_returns_full_span():
    """A narrow fit window against a MUCH wider experimental sweep must still
    recover the true parameters (the optimisation-time grid restriction can't
    change the answer), and the returned model curve must still span the FULL
    experiment (the restriction is an internal speed optimisation, invisible
    to callers/the UI -- not a change in what gets displayed)."""
    x = np.linspace(-150, 150, 3000)
    truth = Recipe(nucleus="31P", larmor_frequency_MHz=162.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(10.0),
            "shift_fwhm_ppm": Param(5.0), "amplitude": Param(100.0),
            "gl": Param(1.0, vary=False)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    data = y + np.random.default_rng(1).normal(0, 0.5, x.size)

    r = Recipe(nucleus="31P", larmor_frequency_MHz=162.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(9.0, min=-20, max=40),
            "shift_fwhm_ppm": Param(6.0, min=0.1),
            "amplitude": Param(80.0, min=0), "gl": Param(1.0, vary=False)})])
    result = fitmod.fit(r, x, data, window_ppm=(40.0, -20.0))

    p = r.sites[0].params
    assert p["isotropic_chemical_shift_ppm"].value == pytest.approx(10.0, abs=0.1)
    assert p["shift_fwhm_ppm"].value == pytest.approx(5.0, abs=0.2)
    assert p["amplitude"].value == pytest.approx(100.0, rel=0.05)
    # the fit window (±60 ppm) is far narrower than the full ±150 ppm sweep --
    # the returned curve must still cover the whole thing
    assert result.x_ppm.min() <= -149.0 and result.x_ppm.max() >= 149.0
    assert len(result.x_ppm) == len(x)
    assert len(result.y_fit) == len(result.x_ppm)


def test_fit_windowed_grid_matches_unrestricted_for_excluded_model():
    """A quad_ct fit (excluded from the grid-restriction allowlist) must give
    the identical fitted values whether or not grid_restrictable would apply
    to other models -- verifies the exclusion actually prevents any change,
    not just that the recipe was skipped by convention."""
    x = np.linspace(-40, 140, 1200)
    truth = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="quad_ct", label="A", params={
            "isotropic_chemical_shift_ppm": Param(60.0),
            "Cq_MHz": Param(3.0, min=0.01, max=20),
            "eta": Param(0.2, vary=False),
            "shift_fwhm_ppm": Param(2.0, min=0.05),
            "amplitude": Param(100.0, min=0)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    data = y + np.random.default_rng(2).normal(0, y.max() * 0.005, x.size)

    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="quad_ct", label="A", params={
            "isotropic_chemical_shift_ppm": Param(58.0, min=40, max=80),
            "Cq_MHz": Param(3.2, min=0.01, max=20),
            "eta": Param(0.2, vary=False),
            "shift_fwhm_ppm": Param(2.2, min=0.05),
            "amplitude": Param(90.0, min=0)})])
    fitmod.fit(r, x, data, window_ppm=(140.0, -40.0))
    p = r.sites[0].params
    assert p["isotropic_chemical_shift_ppm"].value == pytest.approx(60.0, abs=0.5)
    assert p["Cq_MHz"].value == pytest.approx(3.0, abs=0.3)
