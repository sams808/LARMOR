"""The Czjzek probability distribution p(C_Q, η) implied by a fitted σ.

For a standard (d = 5) Czjzek distribution the joint PDF is (d'Espinose de
Lacaillerie et al. 2008, Eq. 6; here in C_Q rather than ν_Q units, with σ =
`sigma_Cq_MHz`):

    p(C_Q, η) ∝ C_Q⁴ · η · (1 − η²/9) · exp[ −C_Q²(1 + η²/3) / (2σ²) ]

C_Q and η are coupled (the η prefactor peaks near η ≈ 0.6), and the single width
σ sets the whole spread. Plotting it turns the fitted σ into the physical
distribution it stands for — the marginal P(C_Q) peaks at C_Q ≈ 2σ.
"""
from __future__ import annotations

import numpy as np


def czjzek_pdf(sigma_MHz: float, cq: np.ndarray, eta: np.ndarray) -> np.ndarray:
    """Un-normalised joint PDF p(C_Q, η) on a grid; cq, eta are 1-D axes.
    Returns a 2-D array of shape (len(eta), len(cq)) (row = η, col = C_Q)."""
    sigma = max(float(sigma_MHz), 1e-6)
    CQ, ETA = np.meshgrid(np.asarray(cq, float), np.asarray(eta, float))
    p = (CQ ** 4 * ETA * (1.0 - ETA ** 2 / 9.0)
         * np.exp(-CQ ** 2 * (1.0 + ETA ** 2 / 3.0) / (2.0 * sigma ** 2)))
    p[p < 0] = 0.0
    s = p.sum()
    return p / s if s > 0 else p


def marginal_cq(sigma_MHz: float, cq: np.ndarray) -> np.ndarray:
    """Marginal P(C_Q) = ∫ p(C_Q, η) dη (η from 0 to 1). Peaks near C_Q = 2σ."""
    eta = np.linspace(0.0, 1.0, 101)
    p = czjzek_pdf(sigma_MHz, cq, eta)
    m = p.sum(axis=0)
    s = m.sum()
    return m / s if s > 0 else m


def rms_pq(sigma_MHz: float) -> float:
    """√⟨P_Q²⟩ for a Czjzek distribution: ⟨C_Q²(1+η²/3)⟩ = 5σ² (Eq. 7), so
    √⟨P_Q²⟩ = √5·σ."""
    return float(np.sqrt(5.0) * float(sigma_MHz))


def suggested_cq_axis(sigma_MHz: float, n: int = 300) -> np.ndarray:
    """A C_Q axis covering the bulk of the distribution (0 to ~5σ)."""
    return np.linspace(0.0, max(5.0 * float(sigma_MHz), 1.0), n)
