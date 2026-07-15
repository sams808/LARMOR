import numpy as np
import pytest

from larmor import fit as fitmod
from larmor.io import fxmla
from larmor.recipe import Param, Recipe, SiteModel

from conftest import CAALGLASS, require


def _two_site_recipe():
    return Recipe(
        nucleus="27Al", larmor_frequency_MHz=195.483, spin_rate_Hz=33296.0,
        sites=[
            SiteModel(model="gauss_lor", label="a", params={
                "isotropic_chemical_shift_ppm": Param(10.0),
                "shift_fwhm_ppm": Param(5.0, min=0.1),
                "amplitude": Param(100.0, min=0.0),
                "gl": Param(1.0, vary=False),
            }),
            SiteModel(model="gauss_lor", label="b", params={
                "isotropic_chemical_shift_ppm": Param(-10.0),
                "shift_fwhm_ppm": Param(5.0, min=0.1),
                "amplitude": Param(50.0, min=0.0),
                "gl": Param(1.0, vary=False),
            }),
        ],
    )


def test_translate_expr():
    r = _two_site_recipe()
    assert fitmod.translate_expr("0.5 * s0.amplitude", r) == "0.5 * s0_amp"
    assert fitmod.translate_expr("s1.shift_fwhm_ppm + 2", r) == "s1_fwhm + 2"
    assert fitmod.translate_expr(
        "s0.amplitude / s1.amplitude", r) == "s0_amp / s1_amp"
    # non-reference text passes through untouched
    assert fitmod.translate_expr("2 * sin(1)", r) == "2 * sin(1)"


def test_translate_expr_errors():
    r = _two_site_recipe()
    with pytest.raises(fitmod.ConstraintError, match="site s9 does not exist"):
        fitmod.translate_expr("s9.amplitude", r)
    with pytest.raises(fitmod.ConstraintError, match="no parameter 'sigma_Cq_MHz'"):
        fitmod.translate_expr("s0.sigma_Cq_MHz", r)  # gauss_lor has no sigma


def test_make_params_resolves_links():
    r = _two_site_recipe()
    r.sites[1].params["amplitude"].expr = "0.5 * s0.amplitude"
    params = fitmod._make_params(r)
    assert params["s1_amp"].value == pytest.approx(50.0)
    assert not params["s1_amp"].vary
    # invalid expression fails at build time, not inside the minimizer
    r.sites[0].params["amplitude"].expr = "nonsense("
    with pytest.raises(fitmod.ConstraintError):
        fitmod._make_params(r)


def test_constrained_fit_synthetic():
    """Fit synthetic two-Gaussian data with a locked 2:1 amplitude ratio."""
    from larmor.engine import gauss_lor

    x = np.linspace(-60, 60, 4000)
    truth = gauss_lor(x, 12.0, 8.0, 200.0, 1.0) + gauss_lor(x, -8.0, 8.0, 100.0, 1.0)
    rng = np.random.default_rng(0)
    y = truth + rng.normal(0, 1.0, x.size)

    r = _two_site_recipe()
    r.sites[1].params["amplitude"].expr = "0.5 * s0.amplitude"
    r.sites[1].params["shift_fwhm_ppm"].expr = "s0.shift_fwhm_ppm"

    result = fitmod.fit(r, x, y)
    a0 = r.sites[0].params["amplitude"]
    a1 = r.sites[1].params["amplitude"]
    assert a0.value == pytest.approx(200.0, rel=0.02)
    assert a1.value == pytest.approx(0.5 * a0.value, rel=1e-9)   # ratio is exact
    assert a1.stderr == pytest.approx(0.5 * a0.stderr, rel=1e-6)  # error propagates
    assert r.sites[1].params["shift_fwhm_ppm"].value == pytest.approx(
        r.sites[0].params["shift_fwhm_ppm"].value)
    # expr survives the fit and the JSON round trip
    assert r.sites[1].params["amplitude"].expr == "0.5 * s0.amplitude"
    back = Recipe.from_dict(r.to_dict())
    assert back.sites[1].params["amplitude"].expr == "0.5 * s0.amplitude"


@pytest.mark.slow
def test_constrained_fit_caalglass():
    """The ssNake-style use case on real data: lock amp ratio + share width.

    Ratio 0.29 is near the freely-refined optimum, so the constrained fit
    stays healthy: no parameters at bounds, full error bars.
    """
    dm = fxmla.read(require(CAALGLASS))
    recipe, _ = fxmla.to_recipe(dm)
    recipe.sites[1].params["amplitude"].expr = "0.29 * s0.amplitude"
    recipe.sites[1].params["shift_fwhm_ppm"].expr = "s0.shift_fwhm_ppm"

    result = fitmod.fit(recipe, dm.spectrum.ppm, dm.spectrum.amplitude,
                        window_ppm=(150.0, -80.0))
    assert result.lmfit_result.errorbars
    a0 = recipe.sites[0].params["amplitude"]
    a1 = recipe.sites[1].params["amplitude"]
    assert a1.value == pytest.approx(0.29 * a0.value, rel=1e-9)
    assert result.rmsd < 0.01  # constrained, so worse than free (0.0025) but still good


@pytest.mark.slow
def test_bad_constraint_is_diagnosed():
    """A ratio that fights the data drives parameters to bounds -- LARMOR
    must say so instead of silently returning a fit without uncertainties."""
    dm = fxmla.read(require(CAALGLASS))
    recipe, _ = fxmla.to_recipe(dm)
    recipe.sites[1].params["amplitude"].expr = "0.5 * s0.amplitude"  # too big
    recipe.sites[1].params["shift_fwhm_ppm"].expr = "s0.shift_fwhm_ppm"

    result = fitmod.fit(recipe, dm.spectrum.ppm, dm.spectrum.amplitude,
                        window_ppm=(150.0, -80.0))
    assert result.at_bounds, "expected at-bound parameters to be reported"
    assert any("at a bound" in n for n in recipe.notes)
    # covariance retry with pinned boundary params should recover error bars
    assert result.lmfit_result.errorbars
