import numpy as np
import pytest

from larmor import models, processing
from larmor.models.base import SimContext
from larmor.quantify import quantify
from larmor.recipe import Param, Recipe, SiteModel


def test_registry_contents():
    names = {m["name"] for m in models.describe_all()}
    assert {"gauss_lor", "czjzek", "quad_ct", "csa_mas"} <= names
    # every model exposes the two params the fit engine relies on
    for m in models.REGISTRY.values():
        assert "isotropic_chemical_shift_ppm" in m.param_names
        assert "amplitude" in m.param_names
        assert m.key_of("amplitude") == "amp"


def test_unknown_model_message():
    with pytest.raises(ValueError, match="unknown site model"):
        models.get("nope")


@pytest.mark.slow
def test_quad_ct_physics():
    """Second-order QIS must shift the CT peak BELOW delta_iso, more for
    larger Cq -- the physics dmfit users check first."""
    ctx = SimContext("27Al", 195.483, 20000.0, np.linspace(-100, 150, 2048))
    peaks = {}
    for cq in (2.0, 6.0):
        v = {"isotropic_chemical_shift_ppm": 50.0, "Cq_MHz": cq, "eta": 0.2,
             "shift_fwhm_ppm": 1.0, "amplitude": 1.0}
        y = models.get("quad_ct").render(v, ctx)
        peaks[cq] = ctx.x_ppm[np.argmax(y)]
    assert peaks[2.0] < 50.0
    assert peaks[6.0] < peaks[2.0]  # bigger Cq -> bigger QIS


def test_quantify_fractions():
    x = np.linspace(-100, 100, 2000)
    recipe = Recipe(nucleus="27Al", larmor_frequency_MHz=195.5, spin_rate_Hz=0,
                    sites=[
        SiteModel(model="gauss_lor", label="a", params={
            "isotropic_chemical_shift_ppm": Param(20.0),
            "shift_fwhm_ppm": Param(8.0), "amplitude": Param(2.0, stderr=0.1),
            "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="b", params={
            "isotropic_chemical_shift_ppm": Param(-20.0),
            "shift_fwhm_ppm": Param(8.0), "amplitude": Param(1.0, stderr=0.1),
            "gl": Param(1.0, vary=False)}),
    ])
    q = quantify(recipe, window_ppm=(100.0, -100.0))
    fracs = [r["fraction_pct"] for r in q["rows"]]
    # equal widths, amplitude ratio 2:1 -> integral fractions 66.7 / 33.3
    assert fracs[0] == pytest.approx(66.67, abs=0.5)
    assert fracs[1] == pytest.approx(33.33, abs=0.5)
    assert q["rows"][0]["fraction_err_pct"] is not None


def test_processing_pipeline_synthetic():
    """em -> zf -> ft -> phase roundtrip on a synthetic fid."""
    sw, n, freq_hz = 10000.0, 1024, 1234.0
    t = np.arange(n) / sw
    fid = np.exp(2j * np.pi * freq_hz * t) * np.exp(-t * 30.0)
    s = processing.Spectrum1D(x_ppm=None, y=fid, sfo1_MHz=100.0, sw_Hz=sw,
                              domain="time")
    s = processing.apply(s, [{"op": "em", "lb_hz": 20},
                             {"op": "zf", "factor": 2},
                             {"op": "ft"}])
    assert s.domain == "freq" and s.y.size == 2048
    # the peak must sit at freq_hz / sfo1 = 12.34 ppm... axis sign: fftshifted
    peak_ppm = s.x_ppm[np.argmax(np.abs(s.y))]
    assert abs(abs(peak_ppm) - 12.34) < 0.2

    # deliberately dephase, then autophase must restore a positive peak
    processing.op_phase(s, p0=90.0)
    assert s.y.real.max() < 0.9 * np.abs(s.y).max()
    processing.op_autophase(s)
    assert s.y.real.max() > 0.95 * np.abs(s.y).max()


def test_processing_baseline():
    x = np.linspace(-50, 50, 3000)
    from larmor.models.analytic import gauss_lor

    signal = gauss_lor(x, 0.0, 4.0, 10.0, 1.0)
    drift = 0.5 + 0.02 * x + 0.001 * x ** 2
    s = processing.from_processed(x, signal + drift, 100.0)
    s = processing.apply(s, [{"op": "baseline", "order": 2}])
    edges = np.concatenate([s.y.real[:200], s.y.real[-200:]])
    assert np.abs(edges).max() < 0.15   # drift removed at the signal-free edges
    assert s.y.real.max() == pytest.approx(10.0, rel=0.05)


def test_ops_validation():
    s = processing.from_processed(np.linspace(-1, 1, 10), np.zeros(10), 100.0)
    with pytest.raises(ValueError, match="unknown processing op"):
        processing.apply(s, [{"op": "sorcery"}])
    with pytest.raises(ValueError, match="time-domain"):
        processing.apply(s, [{"op": "em", "lb_hz": 10}])
