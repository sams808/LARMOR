"""A recipe that omits a parameter's bounds must still inherit the model's
physical bounds at fit time — so a released width can never go negative."""
import numpy as np

from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine, fit as fitmod


def _recipe(width_min):
    return Recipe(nucleus="31P", larmor_frequency_MHz=162.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(5.0, vary=False),
            "shift_fwhm_ppm": Param(6.0, min=width_min),   # None = no recipe bound
            "amplitude": Param(100.0, min=0.0),
            "gl": Param(1.0, vary=False)})])


def test_make_params_falls_back_to_model_min():
    # recipe width has min=None → lmfit gets the model's physical floor (0.1)
    p = fitmod._make_params(_recipe(None))
    name = fitmod._lmfit_name(0, _recipe(None).sites[0], "shift_fwhm_ppm")
    assert p[name].min == 0.1
    # an explicit recipe bound is honoured over the model default
    p2 = fitmod._make_params(_recipe(2.0))
    assert p2[name].min == 2.0


def test_released_width_cannot_go_negative():
    r = _recipe(None)
    x = np.linspace(-60, 60, 1200)
    _, m, _ = engine.simulate(r, exp_ppm=x)
    data = m + np.random.default_rng(0).normal(0, m.max() * 0.02, x.size)
    r.sites[0].params["shift_fwhm_ppm"].vary = True    # release the width
    fitmod.fit(r, x, data, window_ppm=(60, -60))
    assert r.sites[0].params["shift_fwhm_ppm"].value >= 0.1   # stayed physical
