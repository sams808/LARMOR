"""Physics tests for the extended model catalogue and fit zones."""
import numpy as np
import pytest

from larmor import models
from larmor import fit as fitmod
from larmor.models.base import SimContext
from larmor.recipe import Param, Recipe, SiteModel


CTX_AL = SimContext("27Al", 195.483, 20000.0, np.linspace(-150, 200, 2048))


def test_registry_has_seven_models():
    names = {m["name"] for m in models.describe_all()}
    assert {"gauss_lor", "czjzek", "ext_czjzek", "quad_ct", "quad_first",
            "quad_csa", "csa_mas"} <= names


@pytest.mark.slow
def test_ext_czjzek_physics():
    v = {"isotropic_chemical_shift_ppm": 30.0, "Cq_MHz": 6.0, "eta": 0.2,
         "eps": 0.2, "shift_fwhm_ppm": 2.0, "amplitude": 1.0}
    y1 = models.get("ext_czjzek").render(v, CTX_AL)
    y2 = models.get("ext_czjzek").render({**v, "eps": 1.0}, CTX_AL)
    assert np.isfinite(y1).all() and y1.max() == pytest.approx(1.0, abs=0.01)
    # QIS pulls the peak below delta_iso
    assert CTX_AL.x_ppm[y1.argmax()] < 30.0
    width = lambda y: int((y > y.max() / 2).sum())
    assert width(y2) > width(y1)   # more disorder -> broader


@pytest.mark.slow
def test_quad_first_satellites():
    ctx = SimContext("23Na", 132.29, 12500.0, np.linspace(-800, 800, 4096))
    v = {"isotropic_chemical_shift_ppm": 0.0, "Cq_MHz": 1.5, "eta": 0.1,
         "shift_fwhm_ppm": 1.0, "amplitude": 1.0}
    y = models.get("quad_first").render(v, ctx)
    far = np.abs(ctx.x_ppm) > 150
    assert y[far].sum() / y.sum() > 0.05   # satellite manifold present


@pytest.mark.slow
def test_quad_csa_differs_from_quad_ct():
    vqc = {"isotropic_chemical_shift_ppm": 30.0, "Cq_MHz": 4.0, "eta_q": 0.2,
           "zeta_ppm": 80.0, "eta_cs": 0.4, "shift_fwhm_ppm": 2.0,
           "amplitude": 1.0}
    vct = {"isotropic_chemical_shift_ppm": 30.0, "Cq_MHz": 4.0, "eta": 0.2,
           "shift_fwhm_ppm": 2.0, "amplitude": 1.0}
    yqc = models.get("quad_csa").render(vqc, CTX_AL)
    yct = models.get("quad_ct").render(vct, CTX_AL)
    assert float(np.abs(yqc - yct).max()) > 0.03


def _two_gauss_recipe():
    return Recipe(
        nucleus="27Al", larmor_frequency_MHz=195.5, spin_rate_Hz=0,
        sites=[
            SiteModel(model="gauss_lor", label="a", params={
                "isotropic_chemical_shift_ppm": Param(15.0),
                "shift_fwhm_ppm": Param(4.0, min=0.1),
                "amplitude": Param(100.0, min=0.0),
                "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="b", params={
                "isotropic_chemical_shift_ppm": Param(-15.0),
                "shift_fwhm_ppm": Param(4.0, min=0.1),
                "amplitude": Param(50.0, min=0.0),
                "gl": Param(1.0, vary=False)}),
        ],
    )


def test_fit_zones_union():
    """With zones set, data OUTSIDE the zones must not influence the fit."""
    from larmor.models.analytic import gauss_lor

    x = np.linspace(-60, 60, 3000)
    truth = gauss_lor(x, 15.0, 4.0, 100.0, 1.0) + gauss_lor(x, -15.0, 4.0, 50.0, 1.0)
    # a huge artifact far outside both zones
    artifact = gauss_lor(x, 45.0, 3.0, 500.0, 1.0)
    y = truth + artifact

    r = _two_gauss_recipe()
    r.fit_zones = [[25.0, 5.0], [-5.0, -25.0]]   # around each real peak only
    result = fitmod.fit(r, x, y)
    a0 = r.sites[0].params["amplitude"].value
    a1 = r.sites[1].params["amplitude"].value
    assert a0 == pytest.approx(100.0, rel=0.05)   # artifact ignored
    assert a1 == pytest.approx(50.0, rel=0.05)
    # zones survive the JSON round trip
    back = Recipe.from_dict(r.to_dict())
    assert back.fit_zones == [[25.0, 5.0], [-5.0, -25.0]]


def test_position_offset_link_hz_equivalent():
    """The Hz->ppm conversion used by the position-link dialog."""
    r = _two_gauss_recipe()
    larmor = r.larmor_frequency_MHz
    offset_hz = 1955.0
    expected_ppm = offset_hz / larmor    # 10 ppm at 195.5 MHz
    r.sites[1].params["isotropic_chemical_shift_ppm"].expr = \
        f"s0.isotropic_chemical_shift_ppm + {expected_ppm:.6g}"
    params = fitmod._make_params(r)
    assert params["s1_pos"].value == pytest.approx(15.0 + 10.0, abs=1e-6)


def test_spectrum_background_component_fits():
    """An external measured spectrum can be added as a fit component whose
    amplitude (and ppm shift) are optimized alongside the analytic lines."""
    x = np.linspace(-100, 100, 1024)
    bg = np.exp(-((x + 30) / 25) ** 2)            # unit-peak background
    peak = 1.0 / (1.0 + ((x - 10) / 2.0) ** 2)
    data = 3.0 * bg + 5.0 * peak

    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=100.0)
    rec.sites = [
        SiteModel(model="spectrum", label="bg",
                  ref={"ppm": x.tolist(), "amp": bg.tolist()},
                  params={"amplitude": Param(1.0, min=0.0),
                          "shift_ppm": Param(0.0, min=-10, max=10)}),
        SiteModel(model="gauss_lor", label="pk", params={
            "isotropic_chemical_shift_ppm": Param(10.0),
            "shift_fwhm_ppm": Param(4.0, min=0.1),
            "amplitude": Param(1.0, min=0.0),
            "gl": Param(0.5, min=0.0, max=1.0)}),
    ]
    res = fitmod.fit(rec, x, data, window_ppm=(100.0, -100.0))
    assert res.recipe.sites[0].params["amplitude"].value == pytest.approx(3.0, abs=0.05)
    assert res.recipe.sites[1].params["amplitude"].value == pytest.approx(5.0, abs=0.05)
    assert res.rmsd < 1e-3
