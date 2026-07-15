import numpy as np
import pytest

from larmor.io import fxmla
from larmor.recipe import Recipe

from conftest import CAALGLASS, CAALGLASS_MQ, require


def test_parse_caalglass():
    dm = fxmla.read(require(CAALGLASS))
    assert dm.fit_mode == "Fit 1D"
    assert dm.version == "20110208"
    assert len(dm.dimensions) == 1
    dim = dm.dimensions[0]
    assert dim.nucleus == "27Al"
    assert dim.frequency_MHz == pytest.approx(195.483)
    assert len(dim.lines) == 9

    models = [ln.model_name for ln in dim.lines]
    assert models.count("CzSimple") == 3
    assert models.count("Gaus/Lor") == 2
    assert models.count("ss band") == 4

    cz1 = dim.lines[0]
    assert cz1.params["pos"].value == pytest.approx(66.17629762)
    assert cz1.params["pos"].unit == "ppm"
    assert cz1.params["pos"].fix_flag is True
    assert cz1.params["sCZ_CQ"].value == pytest.approx(4548.650849)
    assert cz1.params["CQ"].value == pytest.approx(2 * cz1.params["sCZ_CQ"].value)


def test_parse_embedded_spectrum():
    dm = fxmla.read(require(CAALGLASS))
    spec = dm.spectrum
    assert spec is not None
    assert not dm.is_2d
    assert spec.amplitude.size == 8192
    # tallest point of the experimental spectrum sits at 58.1 ppm (Phase 0)
    peak_ppm = spec.ppm[np.argmax(spec.amplitude)]
    assert peak_ppm == pytest.approx(58.1, abs=0.5)


def test_parse_mqmas_file():
    dm = fxmla.read(require(CAALGLASS_MQ))
    assert dm.fit_mode == "MQMAS"
    assert dm.is_2d
    dim = dm.dimensions[0]
    assert dim.nucleus == "27Al"
    assert len(dim.lines) == 2
    assert all(ln.model_name == "CzSimple" for ln in dim.lines)
    assert dim.lines[0].params["CQ"].value == pytest.approx(8163.637295)


def test_to_recipe_sigma_convention():
    dm = fxmla.read(require(CAALGLASS))
    recipe, warnings = fxmla.to_recipe(dm)

    czjzek = [s for s in recipe.sites if s.model == "czjzek"]
    assert len(czjzek) == 3
    # THE Phase 0 conversion: sigma[MHz] = sCZ_CQ[kHz] / 2 / 1000
    assert czjzek[0].params["sigma_Cq_MHz"].value == pytest.approx(
        4548.650849 / 2000.0)
    assert czjzek[0].params["isotropic_chemical_shift_ppm"].value == pytest.approx(
        66.176, abs=0.001)

    gl = [s for s in recipe.sites if s.model == "gauss_lor"]
    assert len(gl) == 2
    # ss bands are skipped but reported, never silently dropped
    assert sum("ss band" in w for w in warnings) == 4

    assert recipe.nucleus == "27Al"
    assert recipe.spin_rate_Hz == pytest.approx(33296.15741)
    assert len(recipe.source_sha256) == 64


def test_recipe_roundtrip(tmp_path):
    dm = fxmla.read(require(CAALGLASS))
    recipe, _ = fxmla.to_recipe(dm)
    path = tmp_path / "r.json"
    recipe.save(path)
    back = Recipe.load(path)
    assert back.nucleus == recipe.nucleus
    assert len(back.sites) == len(recipe.sites)
    assert back.sites[0].params["sigma_Cq_MHz"].value == pytest.approx(
        recipe.sites[0].params["sigma_Cq_MHz"].value)
