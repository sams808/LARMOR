"""Copy-ready outputs: LaTeX results table + methods sentence."""
from larmor.recipe import Recipe, SiteModel, Param
from larmor import methods, quantify
import numpy as np


def _czjzek_recipe():
    return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="czjzek", label="Al(IV)", params={
            "isotropic_chemical_shift_ppm": Param(60.0),
            "sigma_Cq_MHz": Param(3.0), "line_fwhm_ppm": Param(8.0),
            "amplitude": Param(100.0)}),
        SiteModel(model="czjzek", label="Al(VI)", params={
            "isotropic_chemical_shift_ppm": Param(5.0),
            "sigma_Cq_MHz": Param(2.0), "line_fwhm_ppm": Param(6.0),
            "amplitude": Param(30.0)})])


def test_latex_table_has_sites_and_populations():
    rec = _czjzek_recipe()
    q = quantify.quantify(rec, (120.0, -40.0))
    tex = methods.latex_table(rec.to_dict(), q)
    assert r"\begin{tabular}" in tex and r"\toprule" in tex
    assert "Al(IV)" in tex and "Al(VI)" in tex
    assert "δiso (ppm)" in tex
    # a population column with numbers
    assert "pop." in tex


def test_methods_sentence_names_nucleus_model_and_errors():
    s = methods.methods_sentence(_czjzek_recipe().to_dict())
    assert "27Al" in s and "130.3 MHz" in s
    assert "Czjzek" in s and "2 sites" in s
    assert "covariance" in s
    mc = methods.methods_sentence(_czjzek_recipe().to_dict(), error_method="montecarlo")
    assert "Monte-Carlo" in mc


def test_latex_handles_gauss_lor_without_quad_columns():
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="B(3)", params={
            "isotropic_chemical_shift_ppm": Param(15.0),
            "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0),
            "gl": Param(1.0, vary=False)})])
    tex = methods.latex_table(rec.to_dict())
    assert "C_Q" not in tex          # no quadrupolar columns for a gl-only model
    assert "FWHM (ppm)" in tex
