"""Reusable constraint sets: capture links/bounds/fixes and re-apply them."""
from larmor.recipe import Recipe, SiteModel, Param
from larmor import constraint_library as clib


def _rec():
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(15.0, min=0, max=30),
            "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0),
            "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="B", params={
            "isotropic_chemical_shift_ppm": Param(2.0),
            "shift_fwhm_ppm": Param(3.0),
            "amplitude": Param(50.0, expr="0.5 * s0.amplitude"),
            "gl": Param(1.0, vary=False)})])


def test_capture_finds_link_bound_and_fix():
    cset = clib.capture(_rec())
    s0 = cset["sites"][0]
    assert s0["isotropic_chemical_shift_ppm"]["min"] == 0
    assert s0["isotropic_chemical_shift_ppm"]["max"] == 30
    assert s0["gl"]["vary"] is False
    assert cset["sites"][1]["amplitude"]["expr"] == "0.5 * s0.amplitude"


def test_apply_restores_constraints_on_fresh_recipe():
    cset = clib.capture(_rec())
    fresh = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(15.0),
            "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0),
            "gl": Param(1.0)}),
        SiteModel(model="gauss_lor", label="B", params={
            "isotropic_chemical_shift_ppm": Param(2.0),
            "shift_fwhm_ppm": Param(3.0), "amplitude": Param(50.0),
            "gl": Param(1.0)})])
    applied = clib.apply(fresh, cset)
    assert "s0.gl" in applied
    assert fresh.sites[0].params["isotropic_chemical_shift_ppm"].min == 0
    assert fresh.sites[0].params["gl"].vary is False
    assert fresh.sites[1].params["amplitude"].expr == "0.5 * s0.amplitude"


def test_apply_skips_missing_sites_and_params():
    cset = clib.capture(_rec())
    one = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(15.0),
            "amplitude": Param(100.0)})])
    applied = clib.apply(one, cset)          # site B absent, gl absent → skipped
    assert all(a.startswith("s0.") for a in applied)


def test_describe():
    assert "link" in clib.describe(clib.capture(_rec()))
    assert clib.describe({"sites": []}) == "empty"
