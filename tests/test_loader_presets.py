"""Shared loader and recipe-embedded processing presets."""
import numpy as np
import pytest

from larmor import loader
from larmor.recipe import Recipe

from conftest import BRUKER_1R, CAALGLASS, CAALGLASS_MQ, MAGLAB_35CL, require


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


def test_bruker_sample_is_the_folder_key():
    """A Bruker source is named by its sample folder (date and operator tokens
    off), not by the title's pulse note; the title's first line and the raw
    folder survive in the provenance."""
    _, _, rec, _, _ = loader.load_any(require(BRUKER_1R))
    assert rec["sample"] == "P1-Bi1-12"
    assert rec["provenance"]["title"] == "27Al failed, MAS stopped"
    assert rec["provenance"]["sample_folder"] == "04272026_P1-Bi1-12_SS_ALP"
    assert rec["nucleus"] == "27Al" and rec["larmor_frequency_MHz"] > 100


def test_maglab_sample_comes_from_the_title_sample_line():
    """EXPNO-per-sample layout: the folder carries no sample, the title's
    second line reads 'Sample LAW3CL0CA'."""
    _, _, rec, _, _ = loader.load_any(require(MAGLAB_35CL / "1" / "pdata" / "1" / "1r"))
    assert rec["sample"] == "LAW3CL0CA"
    assert rec["provenance"]["sample_folder"] == "35Cl_2025-12"


def test_saved_recipe_keeps_its_own_sample(tmp_path):
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=156.28, sample="my glass")
    p = tmp_path / "r.recipe.json"
    r.save(p)
    assert Recipe.load(p).sample == "my glass"


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


def test_apply_processing_reapodize_pipeline_replays_from_the_processed_arrays():
    """A re-apodized 1r/CSV records [hilbert, ift, em, ft]: it starts in the
    frequency domain, so it must replay from the processed arrays -- the old
    'any time-domain op means raw fid' rule refused it for CSV sources."""
    x = np.linspace(-50, 50, 1024)
    g = 0.5
    y = g ** 2 / ((x - 12.3) ** 2 + g ** 2)
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=100.0,
               processing=[{"op": "hilbert"}, {"op": "ift"},
                           {"op": "em", "lb_hz": 50}, {"op": "ft"}])
    ppm, amp, notes = loader.apply_processing(r, x, y, source_path=None)
    assert np.allclose(ppm, x, atol=1e-9)
    assert amp.size == y.size
    assert abs(ppm[int(np.argmax(amp))] - 12.3) < 0.5
    assert amp.max() < y.max()                         # broadened, not replaced
    assert any("replayed 4" in n for n in notes)


def test_time_domain_ops_constant_matches_the_ops_that_refuse_frequency_data():
    """TIME_DOMAIN_OPS is exactly the set of ops that raise on frequency-domain
    input; every other op accepts a frequency-domain spectrum."""
    from larmor import processing as proc

    x = np.linspace(-10, 10, 64)
    y = np.exp(-(x / 3.0) ** 2) + 0j
    minimal = {"tdeff": {"points": 1}, "shift_fid": {"points": 0},
               "swap_echo": {"point": 1}, "lp": {"n_predict": 0}}
    for name in proc.TIME_DOMAIN_OPS:
        s = proc.from_processed(x, y.copy(), 100.0, sw_Hz=1000.0)
        with pytest.raises(ValueError):
            proc.OPS[name](s, **minimal.get(name, {}))
    accepts_freq = ("phase", "hilbert", "magnitude", "sr", "scale", "offset",
                    "real", "imag", "conj", "subtract_avg", "flat_baseline",
                    "normalize", "ift", "baseline", "scale_sw")
    for name in accepts_freq:
        assert name not in proc.TIME_DOMAIN_OPS
        s = proc.from_processed(x, y.copy(), 100.0, sw_Hz=1000.0)
        proc.OPS[name](s)                               # must not raise


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
