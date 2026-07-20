"""Guided-T1 zone integration + outlier-robust fitting, and fit bounds."""
import numpy as np
import pytest

from larmor import series
from larmor import fit as fitmod
from larmor.engine import gauss_lor
from larmor.recipe import Param, Recipe, SiteModel

from conftest import BRUKER_SER, require


# ---------------------------------------------------------------- zones
def test_integrate_zones_picks_the_right_peak():
    x = np.linspace(-50, 50, 2000)
    # two peaks with DIFFERENT build-ups: fast at +20, slow at -20
    delays = np.array([0.01, 0.05, 0.2, 0.5, 1.0, 3.0, 10.0, 30.0])
    slices = np.empty((delays.size, x.size))
    for k, t in enumerate(delays):
        a_fast = 1 - np.exp(-t / 0.3)
        a_slow = 1 - np.exp(-t / 6.0)
        slices[k] = a_fast * gauss_lor(x, 20.0, 4.0, 100.0, 1.0) + \
            a_slow * gauss_lor(x, -20.0, 4.0, 80.0, 1.0)
    zones = [(26.0, 14.0), (-14.0, -26.0)]      # around each peak
    ig = series.integrate_zones(x, slices, zones)
    assert ig.shape == (2, delays.size)
    # each zone's build-up recovers its own T1
    r_fast = series.fit_buildup(delays, ig[0])
    r_slow = series.fit_buildup(delays, ig[1])
    assert r_fast["tau"] == pytest.approx(0.3, rel=0.2)
    assert r_slow["tau"] == pytest.approx(6.0, rel=0.2)


def test_fit_buildup_excludes_outliers():
    delays = np.array([0.01, 0.05, 0.2, 0.5, 1.0, 3.0, 10.0, 30.0])
    clean = 1 - np.exp(-delays / 2.0)
    spoiled = clean.copy()
    spoiled[3] = 5.0             # a wild outlier at index 3

    bad = series.fit_buildup(delays, spoiled)             # keeps everything
    keep = np.ones(delays.size, bool); keep[3] = False
    good = series.fit_buildup(delays, spoiled, keep=keep)  # drops the outlier
    assert abs(good["tau"] - 2.0) < abs(bad["tau"] - 2.0)
    assert good["tau"] == pytest.approx(2.0, rel=0.15)


def test_fit_buildup_needs_enough_points():
    delays = np.array([0.1, 1.0, 10.0])
    keep = np.array([True, False, False])
    with pytest.raises(ValueError, match="at least"):
        series.fit_buildup(delays, np.array([0.1, 0.5, 0.9]), keep=keep)


@pytest.mark.slow
def test_zone_integration_on_real_ser():
    """The whole guided pipeline on the user's real saturation-recovery ser:
    process slices, integrate a zone on the relaxed-slice peak, fit T1 > 0."""
    from larmor import satrec

    expno = require(BRUKER_SER).parent                 # .../32 (owns the ser)
    x, slices = satrec.process_slices(expno, lb_hz=200.0)
    delays, _ = series.read_delays(expno)
    n = min(len(delays), slices.shape[0])
    delays, slices = delays[:n], slices[:n]
    peak = float(x[int(np.argmax(np.abs(slices[-1])))])
    ig = series.integrate_zones(x, slices, [(peak + 30, peak - 30)])
    r = series.fit_buildup(delays, ig[0])
    assert r["tau"] > 0 and np.isfinite(r["tau"])


# ---------------------------------------------------------------- bounds
def _line(pos0, **kw):
    return Recipe(nucleus="27Al", larmor_frequency_MHz=195.5,
                  sites=[SiteModel(model="gauss_lor", params={
                      "isotropic_chemical_shift_ppm": Param(pos0, **kw),
                      "shift_fwhm_ppm": Param(6.0, min=1),
                      "amplitude": Param(80.0, min=0),
                      "gl": Param(1.0, vary=False)})])


def test_bounds_are_enforced_and_flagged():
    x = np.linspace(-50, 50, 1500)
    y = gauss_lor(x, 10.0, 6.0, 100.0, 1.0)
    r = _line(0.0, min=-5.0, max=5.0)          # truth (10) is OUTSIDE the box
    res = fitmod.fit(r, x, y, window_ppm=(50, -50))
    pos = r.sites[0].params["isotropic_chemical_shift_ppm"]
    assert -5.0 - 1e-6 <= pos.value <= 5.0 + 1e-6      # stayed in the box
    assert pos.value == pytest.approx(5.0, abs=0.05)   # pinned at the max
    assert "s0.isotropic_chemical_shift_ppm" in res.at_bounds


def test_no_bounds_recovers_truth():
    x = np.linspace(-50, 50, 1500)
    y = gauss_lor(x, 10.0, 6.0, 100.0, 1.0)
    r = _line(0.0)                              # unconstrained
    fitmod.fit(r, x, y, window_ppm=(50, -50))
    assert r.sites[0].params["isotropic_chemical_shift_ppm"].value == \
        pytest.approx(10.0, abs=0.5)


def test_bounds_survive_recipe_roundtrip(tmp_path):
    r = _line(2.0, min=-1.0, max=8.0)
    p = tmp_path / "r.json"
    r.save(p)
    back = Recipe.load(p)
    q = back.sites[0].params["isotropic_chemical_shift_ppm"]
    assert q.min == -1.0 and q.max == 8.0
