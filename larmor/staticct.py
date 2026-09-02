"""Read (C_Q, η, δiso) off a STATIC second-order CT powder pattern — a
measurement from the pattern's singular features, not a fit.

For a static central transition the second-order quadrupolar frequency is

    ν(θ, φ) = −(ν_Q²/(6 ν₀)) · [I(I+1) − ¾] · (A cos⁴θ + B cos²θ + C)

with ν_Q = 3 C_Q / (2I(2I−1)) and the Freude–Haase coefficients

    A(φ) = −27/8 + (9/4) η cos2φ − (3/8) η² cos²2φ
    B(φ) = +30/8 − (1/2) η²   − 2 η cos2φ + (3/4) η² cos²2φ
    C(φ) = −3/8 + (1/3) η²    − (1/4) η cos2φ − (3/8) η² cos²2φ

Two independent checks pin this expression in the tests: its analytic powder
average reproduces ``convert.ct_second_order_shift_ppm`` exactly, and its
predicted features land on an mrsimulator static simulation.

The features themselves are found NUMERICALLY, as the peaks of the ideal
powder histogram of ν(θ, φ). The classical closed-form feature lists
enumerate critical points on the φ = 0°/90° branches only — but the
DOMINANT horn of the static CT pattern comes from an interior-φ saddle
(measured: the tallest simulated maximum sat at +604 ppm for an ⁸¹Br-like
case while the branch enumeration put nothing there), so a histogram of the
exact frequency surface is both simpler and correct.

Feature positions scale as C_Q²/ν₀², so their RATIOS depend only on (η, I):
``read_cq_eta`` matches the measured horn/edge ratios against a per-spin
reference table (cached), then C_Q comes from the horn separation and δiso
from the absolute positions. That is the classic singularity reading,
defensible where a single-field lineshape fit is not.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

#: reference field for the cached per-spin ratio table; positions rescale as
#: (Cq / ref_Cq)² · (ref_nu0 / nu0)²
_REF_CQ = 10.0
_REF_NU0 = 100.0


def _abc(eta: float, cos2phi):
    e = float(eta)
    c = np.asarray(cos2phi, float)
    a = -27.0 / 8.0 + (9.0 / 4.0) * e * c - (3.0 / 8.0) * (e * c) ** 2
    b = (30.0 / 8.0 - 0.5 * e * e - 2.0 * e * c
         + (3.0 / 4.0) * (e * c) ** 2)
    cc = (-3.0 / 8.0 + e * e / 3.0 - 0.25 * e * c
          - (3.0 / 8.0) * (e * c) ** 2)
    return a, b, cc


def _prefactor_ppm(Cq_MHz: float, spin: float, larmor_MHz: float) -> float:
    """−(ν_Q²/(6 ν₀))·[I(I+1)−¾], expressed directly in ppm."""
    I = float(spin)
    nu_q_Hz = 3.0 * Cq_MHz * 1e6 / (2.0 * I * (2.0 * I - 1.0))
    F = I * (I + 1.0) - 0.75
    return -(nu_q_Hz ** 2) / (6.0 * larmor_MHz * 1e6) * F / larmor_MHz


def shift_ppm(theta, phi, Cq_MHz: float, eta: float, spin: float,
              larmor_MHz: float):
    """Second-order CT shift (ppm, relative to δiso) at orientation (θ, φ)."""
    ct2 = np.cos(np.asarray(theta, float)) ** 2
    a, b, c = _abc(eta, np.cos(2.0 * np.asarray(phi, float)))
    return _prefactor_ppm(Cq_MHz, spin, larmor_MHz) * ((a * ct2 + b) * ct2 + c)


@dataclass
class StaticFeatures:
    #: divergent maxima of the ideal powder lineshape, ppm rel. δiso,
    #: sorted by POSITION (up to 4 exist; the two tallest are kept)
    horns: tuple
    #: (min, max) ppm rel. δiso — the pattern's support limits
    edges: tuple


def _powder_histogram(eta: float, spin: float, *, n_ct: int, n_phi: int,
                      bins: int):
    """(bin centres, counts) of the ideal static-CT lineshape at the
    reference (Cq, ν₀). cosθ sampled uniformly = area-faithful powder."""
    ct = np.linspace(0.0, 1.0, n_ct)
    phi = np.linspace(0.0, np.pi / 2.0, n_phi)      # symmetry quadrant
    a, b, c = _abc(eta, np.cos(2.0 * phi))
    ct2 = ct[:, None] ** 2
    nu = _prefactor_ppm(_REF_CQ, spin, _REF_NU0) * (
        (a[None, :] * ct2 + b[None, :]) * ct2 + c[None, :])
    counts, edges = np.histogram(nu.ravel(), bins=bins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, counts, (float(nu.min()), float(nu.max()))


def _features_ref(eta: float, spin: float, *, n_ct: int = 601,
                  n_phi: int = 301, bins: int = 900) -> StaticFeatures:
    centres, counts, support = _powder_histogram(
        eta, spin, n_ct=n_ct, n_phi=n_phi, bins=bins)
    # light smoothing so bin noise does not fake maxima
    k = np.array([1.0, 2.0, 3.0, 2.0, 1.0])
    sm = np.convolve(counts.astype(float), k / k.sum(), mode="same")
    local = ((sm[1:-1] >= sm[:-2]) & (sm[1:-1] > sm[2:])
             & (sm[1:-1] > 0.15 * sm.max()))
    idx = np.where(local)[0] + 1
    # merge runs of adjacent bins into one feature at the weighted centre
    horns: list[tuple[float, float]] = []
    if idx.size:
        splits = np.where(np.diff(idx) > 3)[0]
        for grp in np.split(idx, splits + 1):
            w = sm[grp]
            horns.append((float(np.sum(centres[grp] * w) / np.sum(w)),
                          float(w.max())))
    horns.sort(key=lambda t: -t[1])                 # tallest first
    keep = sorted(h for h, _w in horns[:2])
    return StaticFeatures(horns=tuple(keep), edges=support)


def ct_static_features(Cq_MHz: float, eta: float, spin: float,
                       larmor_MHz: float) -> StaticFeatures:
    """Horn and edge positions of the static CT pattern (ppm rel. δiso)."""
    f = _features_ref(float(eta), float(spin))
    scale = ((Cq_MHz / _REF_CQ) ** 2) * ((_REF_NU0 / larmor_MHz) ** 2)
    return StaticFeatures(horns=tuple(h * scale for h in f.horns),
                          edges=(f.edges[0] * scale, f.edges[1] * scale))


@lru_cache(maxsize=8)
def _eta_table(spin: float, n_eta: int = 201) -> tuple:
    """η → reference features at (_REF_CQ, _REF_NU0); ratios are what the
    inversion matches, and they depend only on (η, spin)."""
    rows = []
    for eta in np.linspace(0.0, 1.0, n_eta):
        f = _features_ref(float(eta), float(spin),
                          n_ct=401, n_phi=201, bins=700)
        rows.append((float(eta), f))
    return tuple(rows)


@dataclass
class StaticReading:
    Cq_MHz: float
    eta: float
    delta_iso_ppm: float
    Pq_MHz: float
    #: worst mismatch between the given positions and the solution's (ppm)
    resid_ppm: float
    ok: bool
    message: str = ""


def read_cq_eta(horn_a_ppm: float, horn_b_ppm: float, edge_ppm: float,
                spin: float, larmor_MHz: float,
                cq_max_MHz: float = 120.0) -> StaticReading:
    """Invert the feature map from the two horns and one INFORMATIVE edge.

    The low-frequency limit of the pattern often coincides with a horn (it
    does for the whole η range at spin 3/2), in which case it adds nothing —
    give the opposite, step-like limit instead; both assignments are tried
    and degenerate ones (edge within 5 % of a horn) are skipped. Matching
    runs on position RATIOS (translation- and scale-invariant), then C_Q
    comes from the horn separation and δiso from the absolute positions.
    """
    if larmor_MHz <= 0 or spin < 1.0:
        return StaticReading(0, 0, 0, 0, np.inf, False,
                             "needs a quadrupolar nucleus and a field")
    h_lo, h_hi = sorted((float(horn_a_ppm), float(horn_b_ppm)))
    edge = float(edge_ppm)
    sep_obs = h_hi - h_lo
    if sep_obs <= 0:
        return StaticReading(0, 0, 0, 0, np.inf, False,
                             "the two horn positions coincide")
    r_obs = (edge - h_lo) / sep_obs
    if min(abs(r_obs), abs(r_obs - 1.0)) < 0.05:
        return StaticReading(0, 0, 0, 0, np.inf, False,
                             "that edge sits on a horn and adds nothing — "
                             "mark the opposite (step) limit of the pattern")

    best = None
    for eta, f in _eta_table(float(spin)):
        if len(f.horns) < 2:
            continue
        sep_ref = f.horns[1] - f.horns[0]
        if sep_ref <= 0:
            continue
        for ref_edge in f.edges:
            r_ref = (ref_edge - f.horns[0]) / sep_ref
            if min(abs(r_ref), abs(r_ref - 1.0)) < 0.05:
                continue                            # degenerate assignment
            miss = abs(r_ref - r_obs)
            if best is None or miss < best[0]:
                best = (miss, eta, ref_edge, f)
    if best is None:
        return StaticReading(0, 0, 0, 0, np.inf, False,
                             "no η reproduces this feature pattern")
    _, eta, ref_edge, f = best
    sep_ref = f.horns[1] - f.horns[0]
    # positions scale as Cq²/ν₀² about the reference point
    cq = _REF_CQ * float(np.sqrt(sep_obs / sep_ref)) * (larmor_MHz / _REF_NU0)
    if not (0 < cq <= cq_max_MHz):
        return StaticReading(cq, eta, 0, 0, np.inf, False,
                             f"C_Q {cq:.1f} MHz is outside 0–{cq_max_MHz:g}")
    sol = ct_static_features(cq, eta, spin, larmor_MHz)
    if len(sol.horns) < 2:
        return StaticReading(cq, eta, 0, 0, np.inf, False,
                             "solution lost a horn — feature positions "
                             "are inconsistent")
    ref_pos = [sol.horns[0], sol.horns[1],
               sol.edges[0] if ref_edge == f.edges[0] else sol.edges[1]]
    obs_pos = [h_lo, h_hi, edge]
    diso = float(np.mean([o - r for o, r in zip(obs_pos, ref_pos)]))
    resid = float(max(abs(o - (r + diso))
                      for o, r in zip(obs_pos, ref_pos)))
    pq = cq * float(np.sqrt(1.0 + eta * eta / 3.0))
    return StaticReading(Cq_MHz=cq, eta=float(eta), delta_iso_ppm=diso,
                         Pq_MHz=pq, resid_ppm=resid, ok=True)
