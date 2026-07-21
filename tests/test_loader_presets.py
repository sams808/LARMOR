"""Shared loader and recipe-embedded processing presets."""
import json

import numpy as np
import pytest

from larmor import loader
from larmor.recipe import Param, Recipe, SiteModel

from conftest import CAALGLASS, CAALGLASS_MQ, EXPNO_1901, require


def test_load_fxmla():
    ppm, amp, recipe, meta, warnings = loader.load_any(require(CAALGLASS))
    assert ppm.size == 8192
    assert np.all(np.diff(ppm) > 0)
    assert recipe["nucleus"] == "27Al"
    assert "dmfit" in meta


def test_load_2d_fxmla_is_refused_with_guidance():
    with pytest.raises(ValueError, match="2D .MQMAS. dmfit file"):
        loader.load_any(require(CAALGLASS_MQ))


def test_load_unknown_source(tmp_path):
    p = tmp_path / "x.xyz"
    p.write_text("hello")
    with pytest.raises(ValueError, match="unrecognized source"):
        loader.load_any(p)


def test_recipe_roundtrip_keeps_processing(tmp_path):
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=195.5,
               processing=[{"op": "baseline", "order": 2},
                           {"op": "sr", "sr_hz": 120.0}],
               processing_from_raw=False)
    p = tmp_path / "r.json"
    r.save(p)
    back = Recipe.load(p)
    assert back.processing == r.processing
    assert back.processing_from_raw is False


def test_from_dict_tolerates_unknown_fields():
    """A recipe from a NEWER LARMOR must still open, with a note."""
    d = Recipe(nucleus="27Al").to_dict()
    d["some_future_field"] = {"a": 1}
    r = Recipe.from_dict(d)
    assert r.nucleus == "27Al"
    assert any("unknown recipe fields" in n for n in r.notes)


def test_apply_processing_on_pdata_arrays():
    x = np.linspace(-50, 50, 800)
    from larmor.engine import gauss_lor

    y = gauss_lor(x, 0.0, 5.0, 10.0, 1.0) + 3.0        # constant offset
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=195.5,
               processing=[{"op": "baseline", "order": 1}])
    ppm, amp, notes = loader.apply_processing(r, x, y)
    assert abs(amp[0]) < 0.3          # offset removed
    assert amp.max() == pytest.approx(10.0, rel=0.05)
    assert any("replayed" in n for n in notes)


def test_apply_processing_needs_raw_reports_clearly():
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=195.5,
               processing=[{"op": "em", "lb_hz": 50}, {"op": "ft"}],
               processing_from_raw=True)
    x = np.linspace(-10, 10, 100)
    with pytest.raises(ValueError, match="not a Bruker EXPNO"):
        loader.apply_processing(r, x, np.zeros(100), source_path=None)


def test_reopened_recipe_replays_its_processing(tmp_path):
    """The reproducibility contract: a saved recipe re-derives the exact
    spectrum it was fitted against, from the untouched source file."""
    src = require(CAALGLASS)
    ppm0, amp0, rd, _, _ = loader.load_any(src)

    r = Recipe.from_dict(rd)
    r.source_path = str(src)
    r.processing = [{"op": "scale", "factor": 2.0}]
    path = tmp_path / "with_proc.recipe.json"
    r.save(path)

    ppm1, amp1, rd1, meta, warnings = loader.load_any(path)
    assert amp1.max() == pytest.approx(2.0 * amp0.max(), rel=1e-9)
    assert any("replayed" in w for w in warnings)
    assert rd1["processing"] == [{"op": "scale", "factor": 2.0}]

    # and without replay the raw source data comes back
    ppm2, amp2, _, _, _ = loader.load_any(path, replay=False)
    assert amp2.max() == pytest.approx(amp0.max(), rel=1e-9)


def test_broken_processing_warns_but_still_loads(tmp_path):
    src = require(CAALGLASS)
    _, _, rd, _, _ = loader.load_any(src)
    r = Recipe.from_dict(rd)
    r.source_path = str(src)
    r.processing = [{"op": "not_a_real_op"}]
    path = tmp_path / "broken.recipe.json"
    r.save(path)
    ppm, amp, _, _, warnings = loader.load_any(path)
    assert ppm.size > 0                        # data still usable
    assert any("processing replay failed" in w for w in warnings)
