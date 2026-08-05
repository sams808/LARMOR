"""χ² surface over a parameter pair: the minimum sits at the fitted values."""
import numpy as np

from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine, chi2map


def _recipe_and_data():
    x = np.linspace(-20, 60, 400)
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(15.0, min=-10, max=40),
            "shift_fwhm_ppm": Param(6.0, min=0.1),
            "amplitude": Param(100.0, min=0),
            "gl": Param(1.0, vary=False)})])
    _, m, _ = engine.simulate(rec, exp_ppm=x)
    return rec, x, m


def test_varying_params_excludes_fixed_and_gl():
    rec, _, _ = _recipe_and_data()
    params = chi2map.varying_params(rec.to_dict())
    names = {pn for _, pn, _ in params}
    assert "isotropic_chemical_shift_ppm" in names
    assert "gl" not in names                         # fixed / mix excluded


def test_chi2_minimum_is_at_the_fitted_values():
    rec, x, m = _recipe_and_data()
    A, B, Z, (a0, b0) = chi2map.chi2_surface(
        rec.to_dict(), x, m, (40.0, -20.0),
        axis_a=(0, "isotropic_chemical_shift_ppm"),
        axis_b=(0, "shift_fwhm_ppm"), n=11)
    ib, ia = np.unravel_index(np.argmin(Z), Z.shape)
    # the χ² minimum is at (or adjacent to) the centre = the fitted values
    assert abs(A[ia] - a0) <= (A[1] - A[0]) + 1e-9
    assert abs(B[ib] - b0) <= (B[1] - B[0]) + 1e-9
    assert Z.min() < Z.max()                         # a real surface
