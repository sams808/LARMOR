"""The Czjzek probability distribution p(C_Q, η) implied by a fitted σ.

Every function here takes LARMOR's STORED width ``sigma_Cq_MHz`` -- the σ of
mrsimulator's ``CzjzekDistribution`` (the standard deviation of each of the
five Gaussian tensor components, in C_Q units).  The width that appears in
Czjzek's own formula (Czjzek 1981; d'Espinose de Lacaillerie 2008, Eq. 6) is
twice that:

    σ_Cz = 2·σ        (= dmfit's sCZ_CQ box; mrsimulator: ``sigma_ = 2*sigma``)

so, in the stored σ, the joint density of the general (dimensionality d)
Czjzek distribution reads

    p_d(C_Q, η) ∝ C_Q^(d−1) · η · (1 − η²/9) · exp[ −C_Q²(1 + η²/3) / (2 σ_Cz²) ]

with d = 5 the fully isotropic case (Le Caër & Brand's Gaussian Isotropic
Model, LARMOR's ``czjzek``) and 1 ≤ d < 5 dmfit CzSimple's empirical ``<d>``
family (LARMOR's ``czjzek_d``).  Because r² = C_Q²(1 + η²/3) = P_Q² enters
only through the exponential, P_Q follows a chi law with d degrees of freedom
and scale σ_Cz, which gives the closed-form invariants used below:

    √⟨P_Q²⟩ = √d · σ_Cz = 2√d · σ          (d = 5: 2√5·σ ≈ 4.47σ)
    mode of P_Q = σ_Cz · √(d − 1)          (d = 5: 4σ)
    mode of the |C_Q| marginal ≈ 1.87 σ_Cz ≈ 3.73σ  (d = 5; ≈ dmfit's CQ = 4σ)

The kernel weights the fit engine uses (``engine.CzjzekKernel.weights``) are
mrsimulator's analytic d = 5 density on the same grid; ``czjzek_weights`` here
reproduces them to floating-point precision at d = 5 (pinned by
tests/test_physics_validation.py), so the read-outs in the P(C_Q) dialog, the
batch table and the Methods sentence describe exactly the distribution the
fit used.  Qt-free.
"""
from __future__ import annotations

import numpy as np


def _sigma_cz(sigma_MHz: float) -> float:
    """Czjzek-paper width σ_Cz = 2σ from LARMOR's stored σ (floored)."""
    return 2.0 * max(float(sigma_MHz), 1e-6)


def czjzek_pdf(sigma_MHz: float, cq: np.ndarray, eta: np.ndarray,
               d: float = 5.0) -> np.ndarray:
    """Normalised joint PDF p_d(C_Q, η) on a grid; ``cq``, ``eta`` are 1-D axes.
    Returns a 2-D array of shape (len(eta), len(cq)) (row = η, col = C_Q)."""
    s_cz = _sigma_cz(sigma_MHz)
    d = float(d)
    CQ, ETA = np.meshgrid(np.asarray(cq, float), np.asarray(eta, float))
    p = (CQ ** (d - 1.0) * ETA * (1.0 - ETA ** 2 / 9.0)
         * np.exp(-CQ ** 2 * (1.0 + ETA ** 2 / 3.0) / (2.0 * s_cz ** 2)))
    p[p < 0] = 0.0
    s = p.sum()
    return p / s if s > 0 else p


def marginal_cq(sigma_MHz: float, cq: np.ndarray, d: float = 5.0) -> np.ndarray:
    """Marginal P(C_Q) = ∫ p_d(C_Q, η) dη (η from 0 to 1), normalised to sum 1.
    For d = 5 it peaks near C_Q ≈ 3.7σ (1.87 σ_Cz)."""
    eta = np.linspace(0.0, 1.0, 101)
    p = czjzek_pdf(sigma_MHz, cq, eta, d)
    m = p.sum(axis=0)
    s = m.sum()
    return m / s if s > 0 else m


def rms_pq(sigma_MHz: float, d: float = 5.0) -> float:
    """√⟨P_Q²⟩ of the Czjzek distribution: P_Q is chi-distributed with d degrees
    of freedom and scale σ_Cz = 2σ, so ⟨P_Q²⟩ = d·σ_Cz² and √⟨P_Q²⟩ = 2√d·σ
    (d = 5: 2√5·σ; d'Espinose 2008 Eq. 7 in the σ_Cz convention)."""
    return float(2.0 * np.sqrt(float(d)) * float(sigma_MHz))


def mode_pq(sigma_MHz: float, d: float = 5.0) -> float:
    """Mode of the quadrupolar product P_Q (the chi variable): σ_Cz·√(d − 1) =
    2σ·√(d − 1); 4σ at d = 5, 0 at d = 1."""
    return float(2.0 * float(sigma_MHz) * np.sqrt(max(float(d) - 1.0, 0.0)))


def suggested_cq_axis(sigma_MHz: float, n: int = 300) -> np.ndarray:
    """A C_Q axis covering the distribution (0 to 10σ = 5σ_Cz; only 4.5e-5 of
    the d = 5 mass lies beyond it -- 21 % lies beyond 5σ)."""
    return np.linspace(0.0, max(10.0 * float(sigma_MHz), 1.0), n)


def czjzek_weights(sigma_MHz: float, d: float, cq_grid: np.ndarray,
                   eta_grid: np.ndarray) -> np.ndarray:
    """Row weights for a (C_Q, η) kernel basis (``engine.CzjzekKernel.K``) from
    the general-d Czjzek density, in the kernel's own row order
    (``np.meshgrid(cq_grid, eta_grid, indexing='xy').ravel()`` -- η-major).

    Reproduces mrsimulator's analytic d = 5 weights exactly, including its
    halving of the η = 1 boundary row (``analytical_distributions.py``), so a
    ``czjzek_d`` site at d = 5 renders identically to ``czjzek`` on the same
    kernel.  Normalised to sum 1 (all-zero input returns zeros)."""
    cq_grid = np.asarray(cq_grid, float)
    eta_grid = np.asarray(eta_grid, float)
    s_cz = _sigma_cz(sigma_MHz)
    d = float(d)
    CQ, ETA = np.meshgrid(cq_grid, eta_grid, indexing="xy")
    p = (CQ ** (d - 1.0) * ETA * (1.0 - ETA ** 2 / 9.0)
         * np.exp(-CQ ** 2 * (1.0 + ETA ** 2 / 3.0) / (2.0 * s_cz ** 2)))
    p[ETA == 1.0] /= 2.0                      # the η = 1 boundary, as mrsimulator
    w = p.ravel()
    s = w.sum()
    return w / s if s > 0 else w
