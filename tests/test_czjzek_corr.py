"""Czjzek with a correlated (delta_iso, C_Q) shift (`czjzek_corr`):
delta_iso = pos + slope (C_Q - <C_Q>) on the d = 5 kernel weights.

Guards: the pure shift-sum primitive (area, centroid and variance identities);
slope = 0 reproduces czjzek through render() AND engine.simulate; the centre
of gravity is slope-invariant (the pivot at <C_Q>); the variance ordering that
IS sign-safe (Var(-s) > Var(+s), because Cov(delta_2, C_Q) < 0); seeding
parity with czjzek; the table column; and slope recovery from a synthetic
spectrum.
"""
import os

import numpy as np
import pytest

from larmor import engine, models
from larmor.models.base import SimContext
from larmor.models.quadrupolar import _shift_sum
from larmor.recipe import Param, Recipe, SiteModel

CTX_AL = SimContext("27Al", 195.483, 20000.0, np.linspace(-150, 200, 2048))


def _cog(x, y):
    return float((x * y).sum() / y.sum())


def _var(x, y):
    c = _cog(x, y)
    return float(((x - c) ** 2 * y).sum() / y.sum())


def test_shift_sum_helper_conserves_area_and_moves_centroid():
    x = np.linspace(-200.0, 200.0, 4001)
    rng = np.random.default_rng(1)
    centres = rng.uniform(-40.0, 40.0, 20)
    widths = rng.uniform(3.0, 8.0, 20)
    amps = rng.uniform(0.2, 1.0, 20)
    Y = np.array([a * np.exp(-0.5 * ((x - c) / w) ** 2)
                  for a, c, w in zip(amps, centres, widths)])
    shifts = rng.uniform(-15.0, 15.0, 20)
    y0 = Y.sum(axis=0)
    y = _shift_sum(x, Y, shifts)
    area_q = np.trapezoid(Y, x, axis=1)
    # linear interpolation onto a translated grid conserves each row's sum
    assert np.trapezoid(y, x) == pytest.approx(np.trapezoid(y0, x), rel=1e-9)
    # ... and moves the centroid by exactly the area-weighted mean shift
    mean_shift = float((area_q * shifts).sum() / area_q.sum())
    assert _cog(x, y) - _cog(x, y0) == pytest.approx(mean_shift, abs=1e-6)
    # the variance grows by Var_a(s) + 2 Cov_a(cog_q, s): exact for a mixture
    # of translated components (interpolation adds < h^2/4 per row)
    cog_q = np.array([_cog(x, row) for row in Y])
    wa = area_q / area_q.sum()
    var_s = float((wa * shifts ** 2).sum() - mean_shift ** 2)
    cov = float((wa * cog_q * shifts).sum() - (wa * cog_q).sum() * mean_shift)
    assert _var(x, y) - _var(x, y0) == pytest.approx(var_s + 2.0 * cov,
                                                     rel=1e-3, abs=0.01)
    # zero shifts: the plain row sum
    assert np.allclose(_shift_sum(x, Y, np.zeros(20)), y0, atol=1e-12)
    # rows below tol * max weight are skipped
    wts = np.ones(20)
    wts[3] = 1e-9
    keep = np.ones(20, dtype=bool)
    keep[3] = False
    expect = _shift_sum(x, Y[keep], shifts[keep])
    assert np.allclose(_shift_sum(x, Y, shifts, weights=wts), expect, atol=1e-12)


def _recipe(model, values):
    params = {k: Param(v) for k, v in values.items()}
    return Recipe(nucleus="27Al", larmor_frequency_MHz=195.483,
                  spin_rate_Hz=20000.0,
                  sites=[SiteModel(model=model, label="s", params=params)])


@pytest.mark.slow
def test_slope_zero_equals_czjzek():
    v = {"isotropic_chemical_shift_ppm": 60.0, "sigma_Cq_MHz": 1.5,
         "shift_fwhm_ppm": 5.0, "line_fwhm_ppm": 0.0, "amplitude": 1.0}
    y_c = models.get("czjzek").render(dict(v), CTX_AL)
    y_0 = models.get("czjzek_corr").render(
        {**v, "shift_slope_ppm_per_MHz": 0.0}, CTX_AL)
    assert np.isfinite(y_0).all()
    assert np.allclose(y_0, y_c, rtol=1e-5, atol=1e-6)
    # through engine.simulate: both models take the kernel axis (needs_kernel)
    x = np.linspace(-150.0, 200.0, 2048)
    gx1, y1, _ = engine.simulate(_recipe("czjzek", v), exp_ppm=x)
    gx2, y2, _ = engine.simulate(
        _recipe("czjzek_corr", {**v, "shift_slope_ppm_per_MHz": 0.0}), exp_ppm=x)
    assert np.array_equal(gx1, gx2)
    assert np.allclose(y1, y2, rtol=1e-5, atol=1e-6)


@pytest.mark.slow
def test_centre_of_gravity_is_slope_invariant_and_variance_ordering_is_sign_robust():
    base = {"isotropic_chemical_shift_ppm": 60.0, "sigma_Cq_MHz": 1.5,
            "shift_fwhm_ppm": 2.0, "line_fwhm_ppm": 0.0, "amplitude": 1.0}
    x = CTX_AL.x_ppm

    def y(s):
        return models.get("czjzek_corr").render(
            {**base, "shift_slope_ppm_per_MHz": s}, CTX_AL)

    ys = {s: y(s) for s in (-3.0, -1.0, 0.0, 1.0)}
    for s, yy in ys.items():
        assert np.isfinite(yy).all(), s
        assert yy.max() == pytest.approx(1.0, abs=0.01), s
    # the pivot at <C_Q> makes the centre of gravity slope-invariant
    c0 = _cog(x, ys[0.0])
    assert _cog(x, ys[-1.0]) == pytest.approx(c0, abs=0.1)
    assert _cog(x, ys[1.0]) == pytest.approx(c0, abs=0.1)
    # Cov(delta_2, C_Q) < 0 (high-C_Q rows carry the most negative second-
    # order shift), so Var(-s) - Var(+s) = -4 s Cov > 0 for every s > 0 -- the
    # ordering that holds wherever the variance minimum sits (a positive
    # slope first NARROWS the pattern, so Var(+s) > Var(0) is not safe)
    assert _var(x, ys[-1.0]) > _var(x, ys[1.0])
    assert _var(x, ys[-3.0]) > _var(x, ys[0.0])


def test_seeds_and_table_column():
    from larmor import constraints_util, estimate, fit

    # start_values parity with czjzek on a synthetic static 81Br pattern
    quad = Recipe(nucleus="81Br", larmor_frequency_MHz=216.0, spin_rate_Hz=0.0,
                  sites=[SiteModel(model="quad_ct", label="q", params={
                      "isotropic_chemical_shift_ppm": Param(0.0),
                      "Cq_MHz": Param(30.0), "eta": Param(0.6),
                      "shift_fwhm_ppm": Param(150.0), "amplitude": Param(1.0)})])
    x = np.linspace(-6000.0, 6000.0, 4001)
    _, yq, _ = engine.simulate(quad, exp_ppm=x)
    ref = estimate.start_values("czjzek", x, yq, "81Br", 216.0)
    assert ref.get("sigma_Cq_MHz", 0.0) > 0.0
    assert estimate.start_values("czjzek_corr", x, yq, "81Br", 216.0) == ref
    assert estimate.start_values("czjzek_d", x, yq, "81Br", 216.0) == ref
    # the four partitions + the kernel axis rule
    assert "czjzek_corr" in fit._SIMULATED_MODELS
    assert "czjzek_corr" in engine._GRID_RESTRICTABLE
    assert "czjzek_corr" in constraints_util._NOT_PEAK_FWHM_MODELS
    assert estimate._WIDTH_KEY["czjzek_corr"] == ("sigma_Cq_MHz", True)
    assert engine.needs_kernel(Recipe(
        nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
            SiteModel(model="czjzek_corr", label="c",
                      params={"amplitude": Param(1.0)})]))
    # the slope gets its own, labelled table column
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.table import LinesTable

    m = models.get("czjzek_corr")
    params = {p.name: {"value": p.default, "vary": p.vary,
                       "min": p.min, "max": p.max} for p in m.params}
    rec = {"nucleus": "27Al", "larmor_frequency_MHz": 195.483,
           "spin_rate_Hz": 20000.0,
           "sites": [{"model": "czjzek_corr", "label": "c", "params": params}]}
    t = LinesTable()
    t.rebuild(rec, set())
    assert "shift_slope_ppm_per_MHz" in t._used_keys
    col = 2 + t._used_keys.index("shift_slope_ppm_per_MHz")
    assert t.table.cellWidget(0, col) is not None
    assert "ppm/MHz" in t.table.horizontalHeaderItem(col).text()


@pytest.mark.slow
def test_fit_recovers_slope_from_synthetic():
    from larmor import fit as fitmod

    truth = {"isotropic_chemical_shift_ppm": 60.0, "sigma_Cq_MHz": 2.0,
             "shift_slope_ppm_per_MHz": 1.0, "shift_fwhm_ppm": 4.0,
             "line_fwhm_ppm": 0.0, "amplitude": 1.0}

    def rec(vals):
        params = models.get("czjzek_corr").defaults()
        for k, val in vals.items():
            params[k].value = val
        params["line_fwhm_ppm"].vary = False
        return Recipe(nucleus="27Al", larmor_frequency_MHz=195.483,
                      spin_rate_Hz=20000.0,
                      sites=[SiteModel(model="czjzek_corr", label="c",
                                       params=params)])

    x = np.linspace(-150.0, 200.0, 2048)
    gx, ytrue, _ = engine.simulate(rec(truth), exp_ppm=x)
    y = np.interp(x, gx, ytrue)
    y += np.random.default_rng(0).normal(0.0, 0.01, x.size)
    start = rec({**truth, "shift_slope_ppm_per_MHz": 0.0})
    res = fitmod.fit(start, x, y)
    p = res.recipe.sites[0].params
    assert p["shift_slope_ppm_per_MHz"].value == pytest.approx(1.0, abs=0.3)
    assert p["sigma_Cq_MHz"].value == pytest.approx(2.0, rel=0.15)
    assert not res.at_bounds
