"""Auto Fit (multi-start) and Errors Analysis (chi-square profiles)."""
import numpy as np
import pytest

from larmor import autofit
from larmor.engine import gauss_lor
from larmor.recipe import Param, Recipe, SiteModel


def _synthetic():
    x = np.linspace(-60, 60, 2000)
    y = (gauss_lor(x, 12.0, 8.0, 100.0, 1.0)
         + gauss_lor(x, -14.0, 6.0, 60.0, 1.0))
    rng = np.random.default_rng(3)
    return x, y + rng.normal(0, 0.4, x.size)


def _recipe(pos_a=8.0, pos_b=-8.0):
    return Recipe(nucleus="27Al", larmor_frequency_MHz=195.5, spin_rate_Hz=0,
                  sites=[
        SiteModel(model="gauss_lor", label="a", params={
            "isotropic_chemical_shift_ppm": Param(pos_a, min=-40, max=40),
            "shift_fwhm_ppm": Param(6.0, min=1.0, max=25.0),
            "amplitude": Param(80.0, min=0.0, max=400.0),
            "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="b", params={
            "isotropic_chemical_shift_ppm": Param(pos_b, min=-40, max=40),
            "shift_fwhm_ppm": Param(6.0, min=1.0, max=25.0),
            "amplitude": Param(50.0, min=0.0, max=400.0),
            "gl": Param(1.0, vary=False)}),
    ])


def test_auto_fit_recovers_truth():
    x, y = _synthetic()
    r = _recipe()
    res = autofit.auto_fit(r, x, y, n_starts=8, seed=1)
    assert res.best_rmsd < 0.02
    # trials are sorted best-first and the winner is the reported best
    assert res.trials[0] == pytest.approx(res.best_rmsd)
    assert len(res.trials) >= 2
    a = r.sites[0].params
    assert a["isotropic_chemical_shift_ppm"].value == pytest.approx(12.0, abs=1.0)
    assert any("auto fit" in n for n in r.notes)


def test_auto_fit_escapes_a_local_minimum():
    """Start both lines on top of each other: a plain fit gets stuck, the
    multi-start must do at least as well and usually better."""
    from larmor import fit as fitmod

    x, y = _synthetic()
    plain = _recipe(pos_a=0.0, pos_b=1.0)
    plain_res = fitmod.fit(plain, x, y)

    auto = _recipe(pos_a=0.0, pos_b=1.0)
    res = autofit.auto_fit(auto, x, y, n_starts=14, spread=0.6, seed=5)
    assert res.best_rmsd <= plain_res.rmsd + 1e-9


def test_auto_fit_progress_callback():
    x, y = _synthetic()
    seen = []
    autofit.auto_fit(_recipe(), x, y, n_starts=3,
                     progress=lambda i, n, b: seen.append((i, n)))
    assert seen[-1][0] == seen[-1][1] == 4      # 3 restarts + the initial fit


def test_error_profile_brackets_1sigma():
    from larmor import fit as fitmod

    x, y = _synthetic()
    r = _recipe()
    fitmod.fit(r, x, y)
    prof = autofit.error_profile(r, x, y, site=0,
                                 param="isotropic_chemical_shift_ppm",
                                 n_points=11, span=3.0)
    assert prof.chi2.min() == prof.chi2_min
    # the profile minimum sits at the fitted value
    assert prof.best_value == pytest.approx(12.0, abs=1.0)
    lo, hi = prof.ci68
    assert lo is not None and hi is not None
    assert lo < prof.best_value < hi
    # a well-determined position: 1sigma interval is tight
    assert (hi - lo) < 4.0
    # and 2sigma is wider than 1sigma
    lo95, hi95 = prof.ci95
    if lo95 is not None and hi95 is not None:
        assert (hi95 - lo95) >= (hi - lo)
    assert "isotropic_chemical_shift_ppm" in prof.summary


def test_error_profile_is_parabolic_near_minimum():
    """chi2 must rise on BOTH sides of the minimum -- the signature the
    covariance assumes and that this tool verifies."""
    from larmor import fit as fitmod

    x, y = _synthetic()
    r = _recipe()
    fitmod.fit(r, x, y)
    prof = autofit.error_profile(r, x, y, site=0, param="amplitude",
                                 n_points=9, span=2.5)
    imin = int(np.argmin(prof.chi2))
    assert 0 < imin < len(prof.chi2) - 1
    assert prof.chi2[0] > prof.chi2_min
    assert prof.chi2[-1] > prof.chi2_min
