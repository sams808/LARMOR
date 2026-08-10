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


def test_read_new_format_csdm_fxmla_with_invalid_xml(tmp_path):
    """dmfit #20230120 writes .fxmla files whose <ExpData><Data> holds a CSDM
    JSON payload embedding the raw Bruker metadata -- audit-trail text full
    of bare '<...>' tokens, so the FILE is not well-formed XML and a strict
    parse raises. The lenient fallback must still recover the (clean)
    <FitParameters> island AND the base64 float32 spectrum from the JSON.
    Found on a real fit dmfit itself wrote (27Al zg, 2026-08-10) that
    LARMOR previously could not open at all."""
    import base64
    import json

    amp = np.array([0.0, 1.0, 3.0, 7.0, 3.0, 1.0, 0.0, 0.0], dtype="<f4")
    payload = {"csdm": {
        "version": "1.0", "description": "27Al zg power check",
        "dimensions": [{"type": "linear", "count": 8,
                        "quantity_name": "Frequency",
                        "increment": "-10 Hz",
                        "coordinates_offset": "100 Hz",
                        "origin_offset": "156.28 MHz"}],
        "dependent_variables": [{"type": "internal", "encoding": "base64",
                                 "numeric_type": "float32",
                                 "quantity_type": "scalar",
                                 "components": [
                                     base64.b64encode(amp.tobytes()).decode()]}],
        "application": {"bruker_topspin_parameters": {
            # the raw '<' tokens that make the whole file invalid XML
            "auditp.txt": "(1,<2026-01-19 14:42:48>,<user@host.edu>,<go4>)"}},
    }}
    text = (
        "<dmfit>\n<FitParameters>\n"
        "<DMFitVersion>dmfit #20230120</DMFitVersion>\n"
        "<FitModeAsc>Fit 1D</FitModeAsc>\n"
        "<Dimension>F2<nucleus>27Al</nucleus><frequency>156.28</frequency>"
        "<spinrate>35714</spinrate>\n"
        "<line><ModelName>CzSimple</ModelName><Group>"
        '<amp Fix="*">129.38</amp><pos Fix="*" Unit="ppm">64.58</pos>'
        '<dCS Fix="*" Unit="ppm">9.55</dCS><lb>-3.86</lb>'
        '<CQ>4756.57</CQ><sCZ_CQ Fix="*">2378.29</sCZ_CQ><d>5</d>'
        "</Group></line>\n"
        "</Dimension>\n</FitParameters>\n"
        "<ExpData>\n<Data>\n" + json.dumps(payload) + "\n</Data>\n</ExpData>\n"
        "</dmfit>\n")
    p = tmp_path / "newformat.fxmla"
    p.write_text(text, encoding="utf-8")

    # the whole point: a strict XML parse of this file fails
    import xml.etree.ElementTree as ET
    with pytest.raises(ET.ParseError):
        ET.fromstring(text)

    dm = fxmla.read(p)
    assert dm.version == "dmfit #20230120"
    assert len(dm.dimensions) == 1
    dim = dm.dimensions[0]
    assert dim.nucleus == "27Al" and dim.frequency_MHz == pytest.approx(156.28)
    assert dim.lines[0].params["sCZ_CQ"].value == pytest.approx(2378.29)

    assert dm.spectrum is not None and not dm.is_2d
    assert dm.spectrum.amplitude == pytest.approx(amp)
    # ppm axis: (X0 + i*dX)/Sf, Hz over MHz
    assert dm.spectrum.ppm[0] == pytest.approx(100.0 / 156.28)
    assert dm.spectrum.ppm[1] == pytest.approx(90.0 / 156.28)
    assert dm.comment == "27Al zg power check"

    recipe, warnings = fxmla.to_recipe(dm)
    s = recipe.sites[0]
    assert s.model == "czjzek"
    assert s.params["sigma_Cq_MHz"].value == pytest.approx(2378.29 / 2000.0)
    assert s.params["isotropic_chemical_shift_ppm"].value == pytest.approx(64.58)
