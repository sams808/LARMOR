"""Regressions for the dmfit 'Amorphous' model (Gaussian Cq/eta distribution).

The Amorphous model fits BO3 in 11B (and other well-defined quadrupolar sites
with modest Gaussian disorder). These guard: (a) the narrow-distribution limit
collapses to the exact second-order CT shift, so Cq/delta_iso are trustworthy;
(b) a broad Cq distribution really broadens/shifts the line; (c) the dmfit .fxml
import maps every shape parameter (with the kHz->MHz Cq conversion).
"""
import numpy as np
import pytest


def _render(model, params, nucleus="11B", larmor=160.4606, mas=20000.0,
            x=None):
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import engine
    if x is None:
        x = np.linspace(-60, 60, 5000)
    site = SiteModel(model=model, label="s",
                     params={k: Param(v) for k, v in params.items()})
    r = Recipe(nucleus=nucleus, larmor_frequency_MHz=larmor, spin_rate_Hz=mas,
               sites=[site])
    _, _, per = engine.simulate(r, exp_ppm=x)
    return x, np.clip(per[0], 0.0, None)


@pytest.mark.slow
def test_amorphous_narrow_limit_matches_second_order_shift():
    """A near-delta Cq distribution must reproduce the analytic CT centroid
    delta_iso + delta_2 -- i.e. the fitted Cq/delta_iso are physically exact."""
    from larmor.convert import ct_second_order_shift_ppm, pq_from_cq_eta
    cq, eta, pos = 2.561, 0.2, 0.0
    x, y = _render("amorphous", dict(
        isotropic_chemical_shift_ppm=pos, Cq_MHz=cq, eta=eta,
        Cq_fwhm_MHz=0.0, eta_fwhm=0.0, shift_fwhm_ppm=0.0,
        line_fwhm_ppm=0.0, gl=0.0, amplitude=1.0))
    centroid = float((x * y).sum() / y.sum())
    d2 = ct_second_order_shift_ppm(pq_from_cq_eta(cq, eta), 1.5, 160.4606)
    assert centroid == pytest.approx(pos + d2, abs=0.25)   # ppm


@pytest.mark.slow
def test_amorphous_cq_distribution_broadens_and_shifts():
    """Widening the Cq distribution must broaden the line and pull its centre
    of gravity to lower frequency (more second-order quadrupolar shift)."""
    common = dict(isotropic_chemical_shift_ppm=18.0, Cq_MHz=2.6, eta=0.15,
                  eta_fwhm=0.0, shift_fwhm_ppm=0.0, line_fwhm_ppm=0.0,
                  gl=0.0, amplitude=1.0)
    x, y_narrow = _render("amorphous", {**common, "Cq_fwhm_MHz": 0.05})
    _, y_broad = _render("amorphous", {**common, "Cq_fwhm_MHz": 0.9})

    def cog(y):
        return (x * y).sum() / y.sum()

    def variance(y):                       # robust width (the raw CT horns make
        c = cog(y)                         # half-max FWHM non-monotonic)
        return (y * (x - c) ** 2).sum() / y.sum()

    assert variance(y_broad) > variance(y_narrow)
    assert cog(y_broad) < cog(y_narrow)                    # more neg. shift


def test_amorphous_registered_with_expected_params():
    """The model is in the registry (so it appears in the app) with the dmfit
    parameter set."""
    from larmor.models import base
    m = base.get("amorphous")
    assert m.needs_quadrupolar
    names = set(m.param_names)
    assert {"Cq_MHz", "eta", "Cq_fwhm_MHz", "eta_fwhm", "shift_fwhm_ppm",
            "line_fwhm_ppm", "gl", "amplitude",
            "isotropic_chemical_shift_ppm"} <= names


def test_fxml_import_maps_amorphous(tmp_path):
    """dmfit .fxml 'Amorphous' lines import with the kHz->MHz Cq conversion and
    all shape parameters (flat .fxml layout, Fit=* attributes)."""
    from larmor.io import fxmla
    xml = """<?xml version="1.0" encoding="utf-8" ?>
<NMRFit><FitParameters><FitModeAsc>Fit 1D</FitModeAsc>
<Dimension>F2<nucleus>11B</nucleus><frequency>160.46</frequency>
<spinrate>20000</spinrate>
<line><ModelName>Amorphous</ModelName>
<amp Fit="*">353.65</amp><pos Fit="*" Unit="ppm">18.39</pos><gl>0</gl>
<dCS Unit="ppm">3.23</dCS><em_au>50</em_au><lb Unit="ppm">0.46</lb>
<CQ>2561.14</CQ><etaQ>0.2</etaQ><FWHM_CQ>350</FWHM_CQ></line>
</Dimension></FitParameters></NMRFit>"""
    p = tmp_path / "fit.fxml"
    p.write_text(xml, encoding="utf-8")
    dm = fxmla.read(p)
    recipe, warnings = fxmla.to_recipe(dm)
    assert len(recipe.sites) == 1
    s = recipe.sites[0]
    assert s.model == "amorphous"
    assert s.params["Cq_MHz"].value == pytest.approx(2.56114, abs=1e-4)   # kHz->MHz
    assert s.params["Cq_fwhm_MHz"].value == pytest.approx(0.350, abs=1e-4)
    assert s.params["eta"].value == pytest.approx(0.2)
    assert s.params["shift_fwhm_ppm"].value == pytest.approx(3.23)
    # dmfit's Amorphous amp (area) is converted to a peak amp on import, so it is
    # no longer the raw 353.65, and a note records the conversion.
    assert s.params["amplitude"].value != pytest.approx(353.65, abs=1.0)
    assert s.params["amplitude"].value > 0
    assert any("area" in n.lower() and "peak" in n.lower() for n in recipe.notes)


def test_lorentz_convolve_never_grows_past_the_signal_length():
    """Regression: np.convolve(a, v, mode="same") returns length
    max(len(a), len(v)), NOT len(a) -- if the Lorentzian kernel (sized off
    the FWHM, previously unbounded) ever came out longer than the signal,
    the "broadened" array silently grew past the signal's length, which
    crashed the caller's final np.interp(ctx.x_ppm, kernel.x_ppm, y, ...)
    with "fp and xp are not of the same length". Found fitting a real 11B
    glass dataset (LARMOR validation, 2026-08): an lmfit errorbar-rescue
    retry step pushed a poorly-determined FWHM parameter far outside typical
    values while probing the Jacobian."""
    from larmor.models.quadrupolar import _lorentz_convolve

    y = np.zeros(50)
    for fwhm_pts in (0.0, 1.0, 40.0, 400.0, 4000.0):
        out = _lorentz_convolve(y, fwhm_pts)
        assert len(out) == len(y)
    # a single-point signal is the tightest edge case (half must clamp to 0)
    assert len(_lorentz_convolve(np.zeros(1), 100.0)) == 1
    assert len(_lorentz_convolve(np.zeros(0), 100.0)) == 0


def test_amorphous_renders_without_crashing_at_an_extreme_lorentzian_fwhm():
    """Integration-level guard for the same bug: the full amorphous render
    path (kernel reweight -> broaden -> interp back onto ctx.x_ppm) must not
    crash even when line_fwhm_ppm is driven far outside sane bounds relative
    to a short experimental grid -- exactly the scenario that crashed a real
    batch fit's errorbar-rescue retry."""
    x, y = _render("amorphous", dict(
        isotropic_chemical_shift_ppm=0.0, Cq_MHz=2.6, eta=0.2,
        Cq_fwhm_MHz=0.3, eta_fwhm=0.1, shift_fwhm_ppm=3.0,
        line_fwhm_ppm=500.0, gl=0.0, amplitude=100.0),
        x=np.linspace(-60, 60, 50))          # a short grid, like the real case
    assert len(y) == 50
    assert np.all(np.isfinite(y))
