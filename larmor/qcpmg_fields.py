"""Infinite-field extrapolation of the isotropic chemical shift from the
central-transition centre of gravity measured at two (or more) magnetic fields.

For a half-integer quadrupolar nucleus the observed centre of gravity of the
central-transition MAS spectrum carries a second-order quadrupolar shift that
scales as 1/ν0². Measuring δcg at several fields and extrapolating to
1/ν0² → 0 removes it, giving the true isotropic chemical shift δiso and the
quadrupolar coupling C_Q (Sandland et al. 2004, Eq. 1; Baasner et al. 2014,
Fig. 6; Freude & Haase 1993; Schmidt et al. 2000).

    δcg = δiso − (10⁶/40) · C_Q²(3+η²) / [ν0² · I²(2I−1)²] · (I(I+1) − 3/4)      (1)

so a plot of δcg (ppm) vs 1/ν0² (MHz⁻²) is a straight line whose intercept is
δiso and whose slope gives C_Q (with an assumed η, conventionally 0.7).

**Central-transition assumption / selective vs non-selective pulses.** Equation
(1) is the shift of the *central transition* centre of gravity. In the large-C_Q
limit (C_Q ≳ 1.5 MHz for the systems here) only the ½ ↔ −½ transition is excited
even by a non-selective (hard) pulse — the satellites are too broad — so the
measured centroid is the CT centroid at both fields regardless of pulse
selectivity (Baasner et al. 2014). It is therefore valid to combine a
non-selective field with a CT-selective field as long as both are in that limit
and the centroid is taken over the CT band only. `per_field_selective` records
the excitation for provenance; it does not change Eq. (1) in this limit.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _spin_factor(spin: float) -> float:
    """[I(I+1) − 3/4] / [I²(2I−1)²] for the CT second-order shift."""
    I = float(spin)
    return (I * (I + 1) - 0.75) / (I ** 2 * (2 * I - 1) ** 2)


def cq_from_slope(slope_ppm_MHz2: float, spin: float, eta: float = 0.7) -> float:
    """Invert Eq. (1): C_Q (MHz) from the slope of δcg vs 1/ν0² (ppm·MHz²).

    slope = −(10⁶/40)·C_Q²(3+η²)·[spin factor]  →  C_Q = √(−slope / A),
    A = (10⁶/40)·(3+η²)·[spin factor].
    """
    A = (1.0e6 / 40.0) * (3.0 + eta ** 2) * _spin_factor(spin)
    val = -slope_ppm_MHz2 / A if A else 0.0
    return float(np.sqrt(val)) if val > 0 else 0.0


def dcg_at_field(delta_iso_ppm: float, cq_MHz: float, larmor_MHz: float,
                 spin: float, eta: float = 0.7) -> float:
    """Forward Eq. (1): predicted δcg (ppm) at a given Larmor frequency."""
    A = (1.0e6 / 40.0) * (3.0 + eta ** 2) * _spin_factor(spin)
    return delta_iso_ppm - A * cq_MHz ** 2 / larmor_MHz ** 2


@dataclass
class FieldPoint:
    larmor_MHz: float                 # observe (Larmor) frequency of the nucleus
    dcg_ppm: float                    # measured CT centre of gravity
    dcg_err_ppm: float = 0.0
    ct_selective: bool | None = None  # provenance only (see module docstring)
    label: str = ""


@dataclass
class InfiniteFieldResult:
    delta_iso_ppm: float
    delta_iso_err_ppm: float
    cq_MHz: float
    cq_err_MHz: float
    pq_MHz: float
    eta: float
    spin: float
    slope: float                      # ppm·MHz² (δcg vs 1/ν0²)
    intercept: float                  # == delta_iso_ppm
    points: list[FieldPoint] = field(default_factory=list)

    def line(self, inv_nu2: np.ndarray) -> np.ndarray:
        """δcg on the fit line for given 1/ν0² values (for plotting)."""
        return self.intercept + self.slope * np.asarray(inv_nu2, float)


def infinite_field_diso(points: list[FieldPoint], spin: float,
                        eta: float = 0.7) -> InfiniteFieldResult:
    """Fit δcg = δiso + slope·(1/ν0²) across fields and return δiso, C_Q, P_Q.

    Needs ≥ 2 fields. With exactly 2 the line is exact (errors from the δcg
    uncertainties are propagated); with > 2 a weighted least-squares line is fit.
    """
    if len(points) < 2:
        raise ValueError("need the centre of gravity at at least two fields")
    x = np.array([1.0 / p.larmor_MHz ** 2 for p in points])   # 1/ν0²  (MHz⁻²)
    y = np.array([p.dcg_ppm for p in points])
    err = np.array([p.dcg_err_ppm or 1.0 for p in points])
    w = 1.0 / err ** 2

    # weighted linear fit y = a + b x  (a = δiso, b = slope)
    sw = w.sum()
    sx = (w * x).sum(); sy = (w * y).sum()
    sxx = (w * x * x).sum(); sxy = (w * x * y).sum()
    denom = sw * sxx - sx * sx
    if denom == 0:
        raise ValueError("the two fields are too close to extrapolate")
    b = (sw * sxy - sx * sy) / denom
    a = (sy - b * sx) / sw
    # parameter variances from the weighted fit
    var_a = sxx / denom
    var_b = sw / denom

    cq = cq_from_slope(b, spin, eta)
    # dC_Q/db = C_Q / (2b) → σ_Cq = |C_Q/(2b)|·σ_b
    cq_err = abs(cq / (2.0 * b)) * np.sqrt(var_b) if b else 0.0
    pq = cq * np.sqrt(1.0 + eta ** 2 / 3.0)
    return InfiniteFieldResult(
        delta_iso_ppm=a, delta_iso_err_ppm=float(np.sqrt(var_a)),
        cq_MHz=cq, cq_err_MHz=cq_err, pq_MHz=pq, eta=eta, spin=spin,
        slope=b, intercept=a, points=list(points))


@dataclass
class WidthSplit:
    wq_lo_ppm: float          # quadrupolar width at the LOWER field (ppm)
    wq_hi_ppm: float          # quadrupolar width at the higher field (ppm)
    wcsd_ppm: float           # chemical-shift-distribution width (field-independent)
    ok: bool                  # False when the split is unphysical (see note)
    note: str = ""


def two_field_widths(nu1_MHz: float, fwhm1_ppm: float,
                     nu2_MHz: float, fwhm2_ppm: float) -> WidthSplit:
    """Separate the CT linewidth into a quadrupolar part W_q (∝ 1/ν0², broader at
    low field) and a chemical-shift-distribution part W_csd (field-independent in
    ppm) from the FWHM measured at two fields (Sandland et al. 2004, Eq. 2).

    In ppm, W_q ∝ 1/ν0² and W_csd is constant, so with the lower field = 1:
        FWHM1² = W_q1² + W_csd²
        FWHM2² = W_q1²·(ν1/ν2)⁴ + W_csd²
    → W_q1² = (FWHM1² − FWHM2²)/(1 − (ν1/ν2)⁴),  W_csd² = FWHM1² − W_q1².
    (Sandland writes it in Hz; the ppm form here is equivalent.)
    """
    # order so field 1 is the lower field
    if nu1_MHz > nu2_MHz:
        nu1_MHz, fwhm1_ppm, nu2_MHz, fwhm2_ppm = nu2_MHz, fwhm2_ppm, nu1_MHz, fwhm1_ppm
    r = (nu1_MHz / nu2_MHz) ** 4
    if abs(1.0 - r) < 1e-9:
        return WidthSplit(0, 0, 0, False, "the two fields are too close")
    wq1_sq = (fwhm1_ppm ** 2 - fwhm2_ppm ** 2) / (1.0 - r)
    wcsd_sq = fwhm1_ppm ** 2 - wq1_sq
    ok = wq1_sq >= 0 and wcsd_sq >= 0
    note = ("" if ok else
            "widths do not separate: the higher-field line is broader than the "
            "quadrupolar model allows — check the FWHM values or the CT band.")
    wq1 = float(np.sqrt(max(wq1_sq, 0.0)))
    wcsd = float(np.sqrt(max(wcsd_sq, 0.0)))
    wq2 = wq1 * (nu1_MHz / nu2_MHz) ** 2         # W_q at the higher field
    return WidthSplit(wq1, wq2, wcsd, ok, note)


def centre_of_gravity(ppm: np.ndarray, amp: np.ndarray,
                      lo_ppm: float | None = None,
                      hi_ppm: float | None = None) -> float:
    """Intensity-weighted centre of gravity of a spectrum (over an optional
    ppm window — restrict it to the central-transition band)."""
    ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
    a = np.clip(amp, 0.0, None)
    if lo_ppm is not None and hi_ppm is not None:
        m = (ppm >= min(lo_ppm, hi_ppm)) & (ppm <= max(lo_ppm, hi_ppm))
        ppm, a = ppm[m], a[m]
    s = a.sum()
    return float((ppm * a).sum() / s) if s > 0 else float("nan")
