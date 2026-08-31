"""Starting values measured from the spectrum, and the engine limits that
made a broad-line fit impossible before them."""
import warnings

import numpy as np
import pytest

from larmor import engine, estimate
from larmor.recipe import Param, Recipe, SiteModel

warnings.filterwarnings("ignore")


def _quad(nucleus, mhz, cq, eta=0.6, width=50.0, pos=0.0):
    return Recipe(nucleus=nucleus, larmor_frequency_MHz=mhz, spin_rate_Hz=0.0,
                  sites=[SiteModel(model="quad_ct", label="q", params={
                      "isotropic_chemical_shift_ppm": Param(pos),
                      "Cq_MHz": Param(cq), "eta": Param(eta),
                      "shift_fwhm_ppm": Param(width), "amplitude": Param(1.0)})])


def test_kernel_window_follows_the_data():
    """THE bug behind 'my model is a flat line': the Czjzek kernel was built
    on a FIXED 150 kHz window, which is only 696 ppm at 216 MHz. A site
    outside it rendered as exactly zero."""
    x_wide = np.linspace(-6000.0, 6000.0, 4001)
    sw, ref = engine.kernel_window(x_wide, 216.0)
    assert sw > engine.KERNEL_MIN_SW_HZ
    assert sw / 216.0 >= 12000.0                      # covers the axis in ppm
    # a narrow axis keeps the historical window exactly (no silent change to
    # every fit that already worked)
    sw2, ref2 = engine.kernel_window(np.linspace(-100.0, 200.0, 1000), 130.3)
    assert (sw2, ref2) == (engine.KERNEL_MIN_SW_HZ, 30.0)


def test_a_wide_czjzek_site_is_not_rendered_as_zero():
    """A real 81Br site at 617 ppm fell outside the kernel's own -318..378 ppm
    window and came back all zeros -- the flat red line the user saw."""
    x = np.linspace(-6000.0, 6000.0, 3001)
    rec = Recipe(nucleus="81Br", larmor_frequency_MHz=216.0, spin_rate_Hz=0.0,
                 sites=[SiteModel(model="ext_czjzek", label="A", params={
                     "isotropic_chemical_shift_ppm": Param(617.0),
                     "Cq_MHz": Param(5.0), "eta": Param(0.2), "eps": Param(0.3),
                     "shift_fwhm_ppm": Param(50.0), "line_fwhm_ppm": Param(10.0),
                     "amplitude": Param(1.0)})])
    _, y, _ = engine.simulate(rec, exp_ppm=x)
    assert y.max() > 0.0
    # and the width must keep responding to Cq instead of saturating at the
    # kernel's edge (it stopped at ~694 ppm however large Cq got)
    widths = []
    for cq in (10.0, 20.0, 30.0):
        rec.sites[0].params["Cq_MHz"] = Param(cq)
        gx, yy, _ = engine.simulate(rec, exp_ppm=x)
        widths.append(estimate.band_width_ppm(gx, yy)[1])
    assert widths[1] > 1.3 * widths[0] and widths[2] > 1.2 * widths[1]


def test_cq_ceiling_reaches_the_heavy_halides():
    """40 MHz was below what a real 81Br glass needs (~34), so the ceiling
    sat on top of the answer."""
    from larmor.models.base import REGISTRY
    for name in ("ext_czjzek", "quad_ct", "quad_csa", "quad_first"):
        cq = next(p for p in REGISTRY[name].params if p.name == "Cq_MHz")
        assert cq.max >= 100.0, name


def test_diff_step_is_coarse_for_simulated_models_only():
    """A grid-simulated model's Jacobian is quantisation noise at scipy's
    ~1e-8 step; a closed-form one is smooth and wants the default."""
    from larmor.fit import diff_step_for
    assert diff_step_for(_quad("81Br", 216.0, 20.0)) == pytest.approx(1e-3)
    analytic = Recipe(nucleus="1H", larmor_frequency_MHz=400.0, sites=[
        SiteModel(model="gauss_lor", label="g", params={
            "isotropic_chemical_shift_ppm": Param(0.0),
            "shift_fwhm_ppm": Param(1.0), "gl": Param(0.5),
            "amplitude": Param(1.0)})])
    assert diff_step_for(analytic) is None


def test_band_width_measures_the_band_not_a_horn():
    x = np.linspace(-200.0, 200.0, 4001)
    y = np.exp(-((x - 30.0) / 20.0) ** 2)
    centre, fwhm = estimate.band_width_ppm(x, y)
    assert centre == pytest.approx(30.0, abs=1.0)
    assert fwhm == pytest.approx(2 * 20.0 * np.sqrt(np.log(2)), rel=0.05)
    # a low fraction reports the BREADTH, which is what a powder pattern needs
    assert estimate.band_width_ppm(x, y, frac=0.10)[1] > fwhm


def test_start_values_land_in_the_right_ballpark():
    """The point of seeding from data: a broad 81Br pattern must not start
    from a 5 MHz needle. Within a factor of ~2 is enough for the fit."""
    x = np.linspace(-6000.0, 6000.0, 4001)
    _, y, _ = engine.simulate(_quad("81Br", 216.0, 30.0, width=150.0), exp_ppm=x)
    got = estimate.start_values("quad_ct", x, y, "81Br", 216.0)
    assert 15.0 < got["Cq_MHz"] < 60.0            # the truth is 30
    # an analytic model just takes the observed width
    g = estimate.start_values("gauss_lor", x, y, "81Br", 216.0)
    assert g["shift_fwhm_ppm"] == pytest.approx(
        estimate.band_width_ppm(x, y)[1], rel=1e-6)
    # a flat/empty trace asks for nothing rather than inventing a number
    assert estimate.start_values("quad_ct", x, np.zeros_like(x), "81Br", 216.0) == {}
