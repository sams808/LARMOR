import numpy as np
import pytest

from larmor import models, processing
from larmor.models.base import SimContext
from larmor.quantify import quantify
from larmor.recipe import Param, Recipe, SiteModel


def test_registry_contents():
    names = {m["name"] for m in models.describe_all()}
    assert {"gauss_lor", "czjzek", "quad_ct", "csa_mas"} <= names
    # every analytic model exposes the two params the fit engine relies on;
    # "spectrum" is a data-backed background component (no analytic centre)
    for m in models.REGISTRY.values():
        assert "amplitude" in m.param_names
        assert m.key_of("amplitude") == "amp"
        if m.name not in ("spectrum", "function"):
            assert "isotropic_chemical_shift_ppm" in m.param_names


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


def test_quantify_untagged_recipe_is_byte_identical_and_site_integrals_matches_rows():
    """N3: an untagged recipe's rows are exactly the pre-family result (same
    keys, 66.67/33.33, same note) with empty family lists; site_integrals()
    is the same trapezoid the rows carry; families are sums of AREAS."""
    from larmor.quantify import site_integrals

    def two(fam_a="", fam_b="", wa=8.0, wb=8.0):
        return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0, sites=[
            SiteModel(model="gauss_lor", label="a", family=fam_a, params={
                "isotropic_chemical_shift_ppm": Param(20.0),
                "shift_fwhm_ppm": Param(wa), "amplitude": Param(2.0, stderr=0.1),
                "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="b", family=fam_b, params={
                "isotropic_chemical_shift_ppm": Param(-20.0),
                "shift_fwhm_ppm": Param(wb), "amplitude": Param(1.0, stderr=0.1),
                "gl": Param(1.0, vary=False)})])

    q = quantify(two(), window_ppm=(100.0, -100.0))
    assert set(q["rows"][0]) == {"site", "label", "model", "position_ppm",
                                 "position_err", "integral", "integral_err",
                                 "tail_outside_pct", "fraction_pct",
                                 "fraction_err_pct"}
    assert [r["fraction_pct"] for r in q["rows"]] == pytest.approx([66.67, 33.33], abs=0.5)
    assert q["note"].startswith("fraction errors are first-order")
    assert q["families"] == [] and q["ratios"] == [] and q["untagged"] == [0, 1]
    ints, win = site_integrals(two(), (100.0, -100.0))
    assert win == (100.0, -100.0)
    assert list(ints) == pytest.approx([r["integral"] for r in q["rows"]], rel=1e-12)
    # tagging changes nothing in the rows -- only adds the family block
    qt = quantify(two("BO3", "BO4"), window_ppm=(100.0, -100.0))
    assert qt["rows"] == q["rows"]
    assert [f["family"] for f in qt["families"]] == ["BO3", "BO4"]
    assert qt["family_basis"] == "independent"
    # equal peak heights, widths 8 and 4 ppm: family % follow the AREAS
    qa = quantify(Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        two("BO3", "BO4", 8.0, 4.0).sites[0],
        SiteModel(model="gauss_lor", label="b", family="BO4", params={
            "isotropic_chemical_shift_ppm": Param(-20.0),
            "shift_fwhm_ppm": Param(4.0), "amplitude": Param(2.0),
            "gl": Param(1.0, vary=False)})]), window_ppm=(100.0, -100.0))
    fams = {f["family"]: f["fraction_pct"] for f in qa["families"]}
    assert fams["BO3"] == pytest.approx(66.7, abs=0.3)
    assert fams["BO4"] == pytest.approx(33.3, abs=0.3)
    n4 = next(r for r in qa["ratios"] if r["name"] == "N4")
    assert n4["value"] == pytest.approx(0.333, abs=0.003)


def test_quantify_nan_stderr_reported_as_no_error_not_as_nan():
    """An ill-conditioned covariance (amplitude at/near a bound, degenerate
    with another free parameter) makes lmfit report stderr as NaN rather than
    None -- and NaN is truthy, so a naive `if amp.stderr` would treat it as a
    real error and the Report table would literally print "± nan"."""
    recipe = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="a", params={
            "isotropic_chemical_shift_ppm": Param(15.0, stderr=float("nan")),
            "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0,
                                                             stderr=float("nan")),
            "gl": Param(1.0, vary=False)})])
    row = quantify(recipe, window_ppm=(40.0, -20.0))["rows"][0]
    assert row["integral_err"] is None
    assert row["fraction_err_pct"] is None
    assert row["position_err"] is None


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


def test_quantify_tail_outside_pct_and_widen():
    """The share of a line outside the integration window, and the window
    that contains 99.5 % of every line (the fit-health 'widen' click)."""
    import math

    from larmor.quantify import TAIL_COVERAGE, window_containing_tails

    def gl(pos, label, model="gauss_lor"):
        return SiteModel(model=model, label=label, params={
            "isotropic_chemical_shift_ppm": Param(pos), "shift_fwhm_ppm": Param(8.0),
            "amplitude": Param(1.0), "gl": Param(1.0, vary=False)})

    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[gl(60.0, "p")])
    sigma = 8.0 / 2.3548
    q = quantify(rec, window_ppm=(120.0, -20.0))
    assert q["rows"][0]["tail_outside_pct"] < 0.01
    assert q["axis_ppm"] == [300.0, -300.0]
    assert "lower bound" in q["note"]
    # half the line: the window edges are interpolated on the cumulative area
    assert quantify(rec, window_ppm=(60.0, -20.0))["rows"][0]["tail_outside_pct"] == \
        pytest.approx(50.0, abs=0.5)
    narrow = quantify(rec, window_ppm=(63.0, 57.0))["rows"][0]["tail_outside_pct"]
    assert narrow == pytest.approx(100.0 * (1 - math.erf(3.0 / (sigma * math.sqrt(2)))), abs=2)
    hi, lo, clipped = window_containing_tails(rec, TAIL_COVERAGE)
    assert hi == pytest.approx(60 + 2.807 * sigma, abs=0.5)
    assert lo == pytest.approx(60 - 2.807 * sigma, abs=0.5)
    assert clipped == {}
    assert quantify(rec, window_ppm=(hi, lo))["rows"][0]["tail_outside_pct"] <= 0.5 + 1e-6
    # clipped to an acquired axis that ends at 50 ppm: the share beyond it is named
    hi2, lo2, clipped = window_containing_tails(rec, exp_ppm=np.linspace(-20, 50, 300))
    assert hi2 == 50.0 and lo2 == pytest.approx(lo, abs=1e-9)
    assert list(clipped) == ["p"] and clipped["p"] > 0.1
    # one tail per row; a background row has none
    two = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3,
                 sites=[gl(60.0, "p"), gl(-40.0, "q"),
                        SiteModel(model="function", label="bg", func="a + 0*x", params={
                            "isotropic_chemical_shift_ppm": Param(0.0),
                            "amplitude": Param(1.0), "a": Param(0.1)})])
    rows = quantify(two, window_ppm=(70.0, -60.0))["rows"]
    assert rows[0]["tail_outside_pct"] > 0.05 and rows[1]["tail_outside_pct"] < 0.01
    assert rows[2]["tail_outside_pct"] is None
    hi3, lo3, _ = window_containing_tails(two)
    assert hi3 == pytest.approx(hi, abs=0.5) and lo3 == pytest.approx(-40 - 2.807 * sigma, abs=0.5)
    # the fractions of test_quantify_fractions are untouched by the new keys
    q = quantify(Recipe(nucleus="27Al", larmor_frequency_MHz=195.5, spin_rate_Hz=0, sites=[
        SiteModel(model="gauss_lor", label="a", params={
            "isotropic_chemical_shift_ppm": Param(20.0), "shift_fwhm_ppm": Param(8.0),
            "amplitude": Param(2.0, stderr=0.1), "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="b", params={
            "isotropic_chemical_shift_ppm": Param(-20.0), "shift_fwhm_ppm": Param(8.0),
            "amplitude": Param(1.0, stderr=0.1), "gl": Param(1.0, vary=False)})]),
        window_ppm=(100.0, -100.0))
    assert [r["fraction_pct"] for r in q["rows"]] == pytest.approx([66.67, 33.33], abs=0.5)
