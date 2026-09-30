"""pybaselines-based baseline correction (larmor.pybaseline): the registry,
compute(), the length-scaling rules, the decimated preview, the processing
op and its replay through a recipe's processing chain."""
import inspect
import json
import sys

import numpy as np
import pytest

from larmor import processing as P
from larmor import pybaseline as pyb


def synthetic(n=2048, seed=0, noise=0.3):
    """Two Gaussians (σ 2 and 4 ppm) on a curved baseline, descending axis."""
    x = np.linspace(200.0, -100.0, n)
    peaks = (100 * np.exp(-0.5 * ((x - 60) / 2.0) ** 2)
             + 60 * np.exp(-0.5 * ((x - 10) / 4.0) ** 2))
    base = 8.0 + 0.05 * (x - 50) + 0.0008 * (x - 50) ** 2
    y = peaks + base + np.random.RandomState(seed).normal(0, noise, n)
    return x, y, base


def rms(a):
    return float(np.sqrt(np.mean(np.square(np.asarray(a, float)))))


# ------------------------------------------------------------- registry
def test_registry_is_the_thirteen_methods_under_pybaselines_names():
    assert list(pyb.METHODS) == [
        "arpls", "asls", "airpls", "iarpls", "drpls", "aspls",
        "modpoly", "imodpoly", "penalized_poly",
        "snip", "rolling_ball", "mor", "mormol"]
    from pybaselines import Baseline

    for key, spec in pyb.METHODS.items():
        assert spec.key == key and spec.label and spec.description
        assert spec.family in pyb.FAMILIES
        keys = [p.key for p in spec.params]
        assert keys and len(keys) == len(set(keys))
        # every exposed parameter is one pybaselines takes, under its own name
        sig = inspect.signature(getattr(Baseline, key))
        assert set(keys) <= set(sig.parameters), (key, set(keys) - set(sig.parameters))
        for p in spec.params:
            assert p.kind in ("float", "int", "bool", "choice") and p.tooltip
            if p.kind == "choice":
                assert p.default in p.choices
            if p.log:
                assert p.kind == "float" and p.min is not None and p.min > 0
            if p.scale == "window":
                assert p.kind == "int" and p.max is None
            if p.scale == "lam":
                assert p.log and p.min == pyb.LAM_MIN and p.max == pyb.LAM_MAX
        assert spec.param("nope") is None
    assert {s.family for s in pyb.METHODS.values()} == set(pyb.FAMILIES)
    assert [s.key for s in pyb.METHODS.values() if s.family == "whittaker"] == \
        ["arpls", "asls", "airpls", "iarpls", "drpls", "aspls"]
    assert "10.5281/zenodo.5608581" in pyb.CITATION


@pytest.mark.parametrize("key", list(pyb.METHODS))
def test_every_method_reduces_the_baseline_error(key):
    """With its length-scaled defaults each method recovers the curved
    baseline under two Gaussians: measured ratios 0.00-0.15 of the error of
    doing nothing (asPLS and ModPoly the worst); the bound is 0.25."""
    x, y, base = synthetic()
    p = pyb.default_params(key, y.size)
    est = pyb.compute(x, y, key, **p)
    assert est.shape == y.shape and est.dtype == float and np.isfinite(est).all()
    assert rms(est - base) < 0.25 * rms(y - base)
    # a descending ppm axis (NMR) and the ascending one give the same baseline
    est_up = pyb.compute(x[::-1], y[::-1], key, **p)
    assert np.allclose(est_up[::-1], est, atol=1e-6 * float(np.abs(est).max()))


def test_suggest_lam_rule_and_window_default():
    # the registry λ is written for 2048 points; ×16 per doubling at order 2
    assert pyb.suggest_lam(2048) == 1e5
    assert pyb.suggest_lam(2048, lam_ref=1e6) == 1e6
    assert pyb.suggest_lam(4096) == 1e6           # 1.6e6 rounds to the decade
    assert pyb.suggest_lam(16384) == 1e9
    assert pyb.suggest_lam(65536) == 1e11
    assert pyb.suggest_lam(2048, diff_order=1) == 1e5
    assert pyb.suggest_lam(8) == pyb.LAM_MIN
    assert pyb.suggest_lam(10 ** 7) == pyb.LAM_MAX
    vals = [pyb.suggest_lam(n) for n in (256, 1024, 2048, 4096, 32768, 65536)]
    assert vals == sorted(vals)
    assert pyb.suggest_half_window(2048) == 64
    assert pyb.suggest_half_window(10) == 2


def test_default_params_follow_the_spectrum_length():
    assert pyb.default_params("arpls") == {"lam": 1e5, "diff_order": 2, "max_iter": 50}
    assert pyb.default_params("arpls", 16384)["lam"] == 1e9
    d = pyb.default_params("asls", 2048)
    assert d["lam"] == 1e6 and d["p"] == 0.01
    assert pyb.default_params("snip", 4096) == {
        "max_half_window": 128, "decreasing": False, "smooth_half_window": 0,
        "filter_order": 2}
    assert pyb.default_params("mor", 16) == {"half_window": 2}
    assert pyb.default_params("rolling_ball", 65536)["half_window"] == 2048
    assert pyb.resolved_max(pyb.METHODS["mor"].param("half_window"), 512) == 256
    assert pyb.resolved_max(pyb.METHODS["arpls"].param("lam"), 512) == pyb.LAM_MAX


def test_preview_params_and_decimate():
    x = np.linspace(100.0, -100.0, 65536)
    y = np.zeros(65536)
    xd, yd, f = pyb.decimate(x, y, 8192)
    assert f == 8 and yd.size == 8192 and xd.size == 8192
    assert xd[0] == x[0] and xd[1] == x[8]
    assert pyb.decimate(x[:2048], y[:2048], 8192)[2] == 1
    assert pyb.decimate(x[:8193], y[:8193], 8192)[2] == 2
    assert pyb.decimate(None, y, 8192)[0] is None
    p = {"lam": 1e11, "diff_order": 2, "max_iter": 50}
    q = pyb.preview_params("arpls", p, 8)
    assert q["lam"] == pytest.approx(1e11 / 8 ** 4) and q["max_iter"] == 50
    assert p["lam"] == 1e11                                    # the input is not touched
    assert pyb.preview_params("arpls", {"lam": 1e11, "diff_order": 1, "max_iter": 50},
                              8)["lam"] == pytest.approx(1e11 / 64)
    s = pyb.preview_params("snip", {"max_half_window": 2048, "decreasing": False,
                                    "smooth_half_window": 0, "filter_order": 2}, 8)
    assert s["max_half_window"] == 256 and s["smooth_half_window"] == 0
    assert s["filter_order"] == 2 and s["decreasing"] is False
    assert pyb.preview_params("mor", {"half_window": 3}, 8)["half_window"] == 1
    assert pyb.preview_params("mor", {"half_window": 64}, 1) == {"half_window": 64}


@pytest.mark.parametrize("key", ["arpls", "mor"])
def test_decimated_preview_tracks_the_full_length_baseline(key):
    """The dialog previews a 64k spectrum on every 8th point with λ ÷ 8⁴
    (arPLS) or the window ÷ 8 (Mor): measured rms differences to the full
    result 0.03 and 0.30 on a baseline spanning ~26 units."""
    x, y, base = synthetic(65536)
    p = pyb.default_params(key, y.size)
    full = pyb.compute(x, y, key, **p)
    xd, yd, f = pyb.decimate(x, y, 8192)
    prev = pyb.compute(xd, yd, key, **pyb.preview_params(key, p, f))
    assert f == 8 and prev.size == 8192
    assert rms(prev - full[::f]) < 0.02 * float(np.ptp(base))


def test_compute_rejects_bad_input_readably():
    x, y, _ = synthetic(512)
    with pytest.raises(ValueError, match="unknown pybaselines method 'nope'"):
        pyb.compute(x, y, "nope")
    with pytest.raises(ValueError, match=r"arPLS: λ \(smoothness\) must be between "
                                         r"100 and 1e\+12 \(got 1e\+20\)"):
        pyb.compute(x, y, "arpls", lam=1e20)
    with pytest.raises(ValueError, match="arPLS has no parameter 'bogus'"):
        pyb.compute(x, y, "arpls", bogus=1)
    with pytest.raises(ValueError, match="whole number"):
        pyb.compute(x, y, "arpls", diff_order=1.5)
    with pytest.raises(ValueError, match="must be one of 2, 4, 6, 8"):
        pyb.compute(x, y, "snip", filter_order=3)
    with pytest.raises(ValueError, match="true or false"):
        pyb.compute(x, y, "snip", decreasing="maybe")
    with pytest.raises(ValueError, match=r"Mor: half-window \(points\) must be "
                                         r"between 1 and 256 \(got 300\)"):
        pyb.compute(x, y, "mor", half_window=300)
    with pytest.raises(ValueError, match="must be a number"):
        pyb.compute(x, y, "asls", p="small")
    with pytest.raises(ValueError, match="at least 8 points"):
        pyb.compute(x[:4], y[:4], "arpls")
    with pytest.raises(ValueError, match="x has 100 points, y has 512"):
        pyb.compute(x[:100], y, "arpls")
    # what a hand-edited recipe may hold: a number as text, an integral float
    v = pyb.validate_params("arpls", {"lam": "1e5", "diff_order": 2.0}, 512)
    assert v == {"lam": 1e5, "diff_order": 2, "max_iter": 50}
    assert type(v["diff_order"]) is int and type(v["lam"]) is float
    assert pyb.validate_params("snip", {"filter_order": "4"}, 512)["filter_order"] == 4
    assert pyb.validate_params("snip", {"decreasing": 1}, 512)["decreasing"] is True


def test_compute_bridges_nan_and_inf_points():
    x, y, base = synthetic(1024)
    p = pyb.default_params("arpls", 1024)
    ref = pyb.compute(x, y, "arpls", **p)
    yn = y.copy()
    yn[100:110] = np.nan
    yn[500] = np.inf
    est = pyb.compute(x, yn, "arpls", **p)
    assert est.shape == y.shape and np.isfinite(est).all()
    assert rms(est - ref) < 0.05 * rms(y - base)
    # subtracting keeps exactly the broken points broken
    corr = yn - est
    assert np.isnan(corr[100:110]).all() and np.isinf(corr[500])
    rest = np.delete(corr, list(range(100, 110)) + [500])
    assert np.isfinite(rest).all()
    with pytest.raises(ValueError, match="finite"):
        pyb.compute(x, np.full(1024, np.nan), "arpls")
    # a broken x axis falls back to the point index (the polynomials fit in x)
    xb = x.copy()
    xb[3] = np.nan
    assert np.isfinite(pyb.compute(xb, y, "modpoly", poly_order=3)).all()
    assert np.isfinite(pyb.compute(None, y, "modpoly", poly_order=3)).all()


def test_missing_pybaselines_gives_the_install_hint(monkeypatch):
    monkeypatch.setitem(sys.modules, "pybaselines", None)
    assert pyb.available() is False
    x, y, _ = synthetic(256)
    with pytest.raises(ImportError, match="pip install pybaselines"):
        pyb.compute(x, y, "arpls")
    monkeypatch.undo()
    assert pyb.available() is True


# --------------------------------------------------------- processing op
def test_op_is_registered_as_a_frequency_domain_op():
    assert P.OPS["pybaseline"] is P.op_pybaseline
    assert "pybaseline" not in P.TIME_DOMAIN_OPS
    assert "pybaseline" not in P.DOMAIN_AGNOSTIC_OPS
    assert P.chain_start_domain([{"op": "pybaseline", "method": "arpls"}]) == "freq"
    assert P.chain_start_domain([{"op": "pybaseline", "method": "arpls"},
                                 {"op": "em", "lb_hz": 10}]) == "freq"


def test_op_subtracts_the_baseline_from_the_real_part_only():
    x, y, base = synthetic(1024)
    imag = np.linspace(0.0, 1.0, 1024)
    p = pyb.default_params("arpls", 1024)
    out = P.apply(P.from_processed(x, y + 1j * imag, 100.0),
                  [{"op": "pybaseline", "method": "arpls", **p}])
    est = pyb.compute(x, y, "arpls", **p)
    assert np.allclose(out.y.real, y - est) and np.allclose(out.y.imag, imag)
    assert np.array_equal(out.x_ppm, x) and out.domain == "freq"
    assert rms(out.y.real - (y - base)) < 0.25 * rms(y - base)
    # a step recording only the method takes the length-scaled defaults
    out2 = P.apply(P.from_processed(x, y, 100.0),
                   [{"op": "pybaseline", "method": "mor"}])
    assert np.allclose(out2.y.real,
                       y - pyb.compute(x, y, "mor", **pyb.default_params("mor", 1024)))
    # the default method is arPLS
    out3 = P.apply(P.from_processed(x, y, 100.0), [{"op": "pybaseline"}])
    assert np.allclose(out3.y.real, y - est)
    s_time = P.Spectrum1D(x_ppm=None, y=np.zeros(64, complex), sfo1_MHz=100.0,
                          sw_Hz=1e4, domain="time")
    with pytest.raises(ValueError, match="frequency-domain"):
        P.op_pybaseline(s_time)
    with pytest.raises(ValueError, match="arPLS has no parameter 'order'"):
        P.apply(P.from_processed(x, y, 100.0),
                [{"op": "pybaseline", "method": "arpls", "order": 3}])


def test_recipe_chain_holding_the_op_round_trips_through_the_loader():
    """What a saved fit holds: the step in recipe.processing, through JSON,
    replayed by loader.apply_processing exactly as processing.apply does."""
    from larmor.loader import apply_processing
    from larmor.recipe import Recipe

    x, y, base = synthetic(2048)
    step = {"op": "pybaseline", "method": "airpls",
            **pyb.default_params("airpls", 2048)}
    r = Recipe(sample="s", source_kind="csv", source_path="x.csv", nucleus="27Al",
               larmor_frequency_MHz=104.26,
               processing=[{"op": "phase", "p0": 0.0, "p1": 0.0}, step])
    r2 = Recipe.from_dict(json.loads(json.dumps(r.to_dict())))
    assert r2.processing[1] == step
    ppm, amp, notes = apply_processing(r2, x, y)
    assert notes == ["replayed 2 processing step(s) from the recipe"]
    assert np.all(np.diff(ppm) > 0)                    # the loader sorts ascending
    ref = P.apply(P.from_processed(x, y, 104.26), r2.processing)
    order = np.argsort(ref.x_ppm)
    assert np.allclose(amp, ref.y.real[order]) and np.allclose(ppm, ref.x_ppm[order])
    assert rms(amp - (y - base)[order]) < 0.25 * rms(y - base)


def test_methods_sentence_names_the_method():
    from larmor import acquisition as A

    assert "pybaseline" in A.OP_PHRASES
    assert A.describe_processing([{"op": "pybaseline", "method": "arpls",
                                   "lam": 1e5}]) == "arPLS baseline (pybaselines)"
    assert A.describe_processing([{"op": "pybaseline", "method": "rolling_ball",
                                   "half_window": 64}]) == \
        "Rolling ball baseline (pybaselines)"
    assert A.describe_processing([{"op": "pybaseline", "method": "zzz"}]) == \
        "zzz baseline (pybaselines)"
    assert A.describe_processing([{"op": "em", "lb_hz": 50},
                                  {"op": "pybaseline", "method": "snip"}]) == \
        "exponential apodization (LB 50 Hz) and SNIP baseline (pybaselines)"
