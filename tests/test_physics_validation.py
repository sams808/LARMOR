"""Physics-validation regressions backing docs/validation.md.

These guard the claims that (a) the mrsimulator-simulated central-transition
lineshape agrees with LARMOR's analytic 2nd-order quadrupolar shift formula, and
(b) the Czjzek convention relations are exact. If either drifts, a fit's δiso or
C_Q could be quietly wrong — so this must stay green.
"""
import numpy as np
import pytest


def test_quad_ct_centroid_matches_analytic_second_order_shift():
    """The simulated CT centroid at δiso=0 must equal the analytic δ₂
    (convert.ct_second_order_shift_ppm) — the formula used for MQMAS/QCPMG δiso."""
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import engine
    from larmor.convert import ct_second_order_shift_ppm, pq_from_cq_eta

    lar = 130.3
    for cq, eta in ((4.0, 0.3), (2.5, 0.0), (3.0, 0.8)):
        site = SiteModel(model="quad_ct", label="q", params={
            "isotropic_chemical_shift_ppm": Param(0.0), "Cq_MHz": Param(cq),
            "eta": Param(eta), "shift_fwhm_ppm": Param(0.3),
            "amplitude": Param(1.0)})
        r = Recipe(nucleus="27Al", larmor_frequency_MHz=lar, sites=[site])
        x = np.linspace(-160, 40, 4000)
        _, _, per = engine.simulate(r, exp_ppm=x)
        c = np.clip(per[0], 0, None)
        centroid = float((x * c).sum() / c.sum())
        analytic = ct_second_order_shift_ppm(pq_from_cq_eta(cq, eta), 2.5, lar)
        assert centroid == pytest.approx(analytic, abs=0.5)   # ppm


def test_czjzek_convention_relations_exact():
    """dmfit sCZ_CQ = 2σ (mode of |Cq|) and √⟨PQ²⟩ = √5·σ (Eq. 7)."""
    from larmor import czjzek_dist as cd
    for sigma in (1.0, 1.8, 3.0):
        cq = cd.suggested_cq_axis(sigma, 1200)
        mode = cq[int(np.argmax(cd.marginal_cq(sigma, cq)))]
        assert 1.7 * sigma < mode < 2.1 * sigma           # ≈ 2σ (dmfit sCZ_CQ)
        assert cd.rms_pq(sigma) == pytest.approx(np.sqrt(5) * sigma)


def test_dipolar_and_efg_constants():
    """Spot-check the dipolar and EFG conversion constants against known values."""
    from larmor.convert import dipolar_Hz, cq_from_efg
    # 1H-1H at 1.5 Å ≈ 35.6 kHz
    assert dipolar_Hz(42.577, 42.577, 1.5) == pytest.approx(35.6e3, rel=0.02)
    # EFG constant: Cq[MHz] = 234.9647 * Q[barn] * Vzz[a.u.]
    assert cq_from_efg(1.0, 1.0) == pytest.approx(234.9647, rel=1e-4)
