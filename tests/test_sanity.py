"""Physical sanity checks on a fitted model (flags unphysical parameters)."""
import numpy as np

from larmor.recipe import Recipe, SiteModel, Param
from larmor import sanity


def _rec(**over):
    p = {"isotropic_chemical_shift_ppm": Param(15.0), "shift_fwhm_ppm": Param(6.0),
         "amplitude": Param(100.0), "gl": Param(0.5, vary=False)}
    p.update(over)
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0,
                  sites=[SiteModel(model="gauss_lor", label="A", params=p)])


def test_clean_fit_has_no_warnings():
    assert sanity.check_recipe(_rec(), window=(40.0, -10.0)) == []


def test_negative_amplitude_flagged():
    w = sanity.check_recipe(_rec(amplitude=Param(-5.0)))
    assert any("negative amplitude" in x["message"] for x in w)


def test_negative_width_flagged():
    w = sanity.check_recipe(_rec(shift_fwhm_ppm=Param(-1.0)))
    assert any(x["param"] == "shift_fwhm_ppm" for x in w)


def test_gl_out_of_range_flagged():
    w = sanity.check_recipe(_rec(gl=Param(1.4, vary=False)))
    assert any("mix" in x["message"] for x in w)


def test_eta_out_of_range_flagged():
    from larmor.recipe import SiteModel as SM
    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.0, sites=[
        SM(model="quad_ct", label="Q", params={
            "isotropic_chemical_shift_ppm": Param(60.0),
            "Cq_MHz": Param(3.0), "eta": Param(1.5),
            "line_fwhm_ppm": Param(5.0), "amplitude": Param(100.0)})])
    w = sanity.check_recipe(rec)
    assert any("η" in x["message"] for x in w)


def test_site_outside_window_flagged():
    w = sanity.check_recipe(_rec(), window=(40.0, 20.0))    # 15 ppm is below 20
    assert any("outside the fit window" in x["message"] for x in w)
    assert sanity.summarize(w).startswith("⚠ physical check")


def test_non_finite_flagged():
    w = sanity.check_recipe(_rec(shift_fwhm_ppm=Param(np.nan)))
    assert any("not finite" in x["message"] for x in w)
