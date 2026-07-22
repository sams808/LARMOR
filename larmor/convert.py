"""NMR conversion calculators (ssNake Utilities parity): chemical shift,
quadrupole coupling, and dipolar distance. Pure, tested physics.

Units: gyromagnetic ratios in MHz/T (mrsimulator convention), Cq/PQ in MHz,
distances in angstrom, fields in MHz (¹H) or Tesla.
"""
from __future__ import annotations

import numpy as np

MU0_4PI = 1.0e-7                     # T·m/A
HBAR = 1.054571817e-34              # J·s


# ---------------------------------------------------------------- chemical shift
def ppm_to_Hz(ppm: float, sfo_MHz: float) -> float:
    return ppm * sfo_MHz            # sfo in MHz, ppm·MHz = Hz


def Hz_to_ppm(Hz: float, sfo_MHz: float) -> float:
    return Hz / sfo_MHz if sfo_MHz else 0.0


# ---------------------------------------------------------------- quadrupole
def pq_from_cq_eta(cq_MHz: float, eta: float) -> float:
    """Quadrupolar product PQ = SOQE = Cq·√(1 + η²/3) (MHz)."""
    return cq_MHz * np.sqrt(1.0 + eta ** 2 / 3.0)


def cq_from_pq_eta(pq_MHz: float, eta: float) -> float:
    return pq_MHz / np.sqrt(1.0 + eta ** 2 / 3.0)


def nu_q(cq_MHz: float, spin: float) -> float:
    """First-order quadrupolar frequency νQ = 3Cq / [2I(2I−1)] (MHz)."""
    I = spin
    denom = 2.0 * I * (2.0 * I - 1.0)
    return 3.0 * cq_MHz / denom if denom else 0.0


def cq_from_nu_q(nu_q_MHz: float, spin: float) -> float:
    I = spin
    return nu_q_MHz * 2.0 * I * (2.0 * I - 1.0) / 3.0


def ct_second_order_shift_ppm(pq_MHz: float, spin: float, larmor_MHz: float
                              ) -> float:
    """Isotropic second-order quadrupolar shift of the central transition (ppm),
    δ₂ = −(3/40)·[I(I+1)−3/4]/[I²(2I−1)²]·(PQ/ν₀)²·1e6. Always negative."""
    I = spin
    coeff = (I * (I + 1.0) - 0.75) / (I ** 2 * (2.0 * I - 1.0) ** 2)
    return -(3.0 / 40.0) * coeff * (pq_MHz / larmor_MHz) ** 2 * 1e6


def cq_from_efg(vzz_au: float, quad_moment_barn: float) -> float:
    """Cq (MHz) from an EFG Vzz (atomic units) and the nuclear Q (barn):
    Cq = e·Q·Vzz/h ≈ 234.9647 · Q(barn) · Vzz(a.u.)  [MHz]."""
    return 234.9647 * quad_moment_barn * vzz_au


# ---------------------------------------------------------------- dipolar
def _gamma_rad(g_MHz_per_T: float) -> float:
    return g_MHz_per_T * 2.0 * np.pi * 1.0e6      # rad·s⁻¹·T⁻¹


def dipolar_Hz(g1_MHz_T: float, g2_MHz_T: float, r_angstrom: float) -> float:
    """Dipolar coupling constant b/2π (Hz) between two nuclei r Å apart:
    d = (μ0/4π)·γ1·γ2·ħ / r³, in Hz."""
    if r_angstrom <= 0:
        return 0.0
    g1, g2 = _gamma_rad(g1_MHz_T), _gamma_rad(g2_MHz_T)
    r = r_angstrom * 1.0e-10
    d_rad = MU0_4PI * g1 * g2 * HBAR / r ** 3      # rad/s
    return d_rad / (2.0 * np.pi)


def distance_from_dipolar(g1_MHz_T: float, g2_MHz_T: float, d_Hz: float
                          ) -> float:
    """Internuclear distance (Å) from a dipolar coupling |d| (Hz)."""
    if d_Hz == 0:
        return float("inf")
    g1, g2 = _gamma_rad(g1_MHz_T), _gamma_rad(g2_MHz_T)
    r3 = MU0_4PI * abs(g1 * g2) * HBAR / (2.0 * np.pi * abs(d_Hz))
    return (r3 ** (1.0 / 3.0)) / 1.0e-10
