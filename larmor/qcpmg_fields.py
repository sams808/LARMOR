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
and the centroid is taken over the CT band only. ``FieldPoint.ct_selective``
records the operator's DECLARATION of the excitation for provenance (True /
False / None = not declared); it is not read from the data and does not
change Eq. (1) in this limit.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = [
    "DEFAULT_ETA", "FieldPoint", "InfiniteFieldResult", "WidthSplit",
    "cq_from_slope", "dcg_at_field", "weighted_line", "infinite_field_diso",
    "two_field_widths", "centre_of_gravity",
]


#: conventional assumed eta for the two-field extrapolation (Stebbins &
#: Du 2002). The dialogs seed their spinboxes from THIS constant, so the
#: call sites can no longer drift apart; the user's spinbox still
#: overrides per analysis.
DEFAULT_ETA = 0.7


def _spin_factor(spin: float) -> float:
    """[I(I+1) − 3/4] / [I²(2I−1)²] for the CT second-order shift."""
    I = float(spin)
    return (I * (I + 1) - 0.75) / (I ** 2 * (2 * I - 1) ** 2)


def cq_from_slope(slope_ppm_MHz2: float, spin: float, eta: float = DEFAULT_ETA) -> float:
    """Invert Eq. (1): C_Q (MHz) from the slope of δcg vs 1/ν0² (ppm·MHz²).

    slope = −(10⁶/40)·C_Q²(3+η²)·[spin factor]  →  C_Q = √(−slope / A),
    A = (10⁶/40)·(3+η²)·[spin factor].
    """
    A = (1.0e6 / 40.0) * (3.0 + eta ** 2) * _spin_factor(spin)
    val = -slope_ppm_MHz2 / A if A else 0.0
    return float(np.sqrt(val)) if val > 0 else 0.0


def dcg_at_field(delta_iso_ppm: float, cq_MHz: float, larmor_MHz: float,
                 spin: float, eta: float = DEFAULT_ETA) -> float:
    """Forward Eq. (1): predicted δcg (ppm) at a given Larmor frequency."""
    A = (1.0e6 / 40.0) * (3.0 + eta ** 2) * _spin_factor(spin)
    return delta_iso_ppm - A * cq_MHz ** 2 / larmor_MHz ** 2


#: smallest dcg uncertainty the fit will weight with (ppm). A jitter sigma of
#: 1e-14 ppm (a window that is the whole axis) would otherwise give a weight
#: of 1e28 and cancel the normal-equation denominator to exactly 0.0. Both
#: dialogs and the fitter apply this SAME floor, so a displayed "0.1" is the
#: value that was used.
ERR_FLOOR_PPM = 0.1

#: lever-arm ratio (x.max - x.min) / x.max in 1/nu0^2 below which two fields
#: cannot be extrapolated (5 % in 1/nu0^2 = Larmor frequencies within ~2.5 %)
MIN_LEVER_ARM = 0.05
#: below this the extrapolation is valid but poorly conditioned (a warning)
SHORT_LEVER_ARM = 0.20


@dataclass
class FieldPoint:
    larmor_MHz: float                 # observe (Larmor) frequency of the nucleus
    dcg_ppm: float                    # measured CT centre of gravity
    dcg_err_ppm: float = 0.0          # 0 / NaN / negative = no uncertainty known
    ct_selective: bool | None = None  # operator's declaration; provenance only
    label: str = ""
    # ---- provenance (reported, never fitted) ------------------------------
    magnitude: bool | None = None     # δcg measured on |spectrum| (mc)? None = unknown
    window: tuple | None = None       # (lo, hi) ppm the CG was integrated over
    window_mode: str = ""             # minima | manual | manifold | centreband
    source: str = ""                  # file / EXPNO / 'workspace'
    rotor_Hz: float = 0.0             # MAS rate of the acquisition (0 = static/unknown)
    flags: tuple = ()                 # quality flags of the measurement ('! ...')
    cg_sequence: str = ""             # 'CG(w, 1.5w, 2w, 3w) = ...' convergence record
    drift_ppm: float = float("nan")   # |CG(2w) - CG(w)|
    lb_Hz: float | None = None        # processing line broadening
    sf_MHz: float | None = None       # referenced spectrometer frequency (procs SF)
    sr_hz: float | None = None        # spectral reference SF - BF1 (0 = unreferenced)
    referenced: bool | None = None    # the ppm axis carries a reference
    ref_dev_ppm: float | None = None  # SF vs a same-session 1H reference (ppm)

    @property
    def has_err(self) -> bool:
        """True when dcg_err_ppm is a usable (finite, positive) uncertainty."""
        e = float(self.dcg_err_ppm)
        return bool(np.isfinite(e) and e > 0.0)

    @classmethod
    def from_measurement(cls, larmor_MHz: float, meas, *, magnitude=None,
                         source: str = "", rotor_Hz: float = 0.0,
                         ct_selective=None, label: str = "",
                         dcg_ppm: float | None = None,
                         dcg_err_ppm: float | None = None,
                         meta: dict | None = None,
                         ref_sf_h_MHz: float | None = None,
                         nucleus: str = "") -> "FieldPoint":
        """A point from a :class:`larmor.qcpmg.CgMeasurement`: the error is
        max(jitter sigma, convergence drift) floored at ERR_FLOOR_PPM, and
        window, mode, flags and the convergence sequence travel with it.
        ``dcg_ppm`` / ``dcg_err_ppm`` override the measured values when the
        user edited the cells. ``meta`` (a saved dataset's header or a
        Bruker meta) supplies LB, SF, SR and the referenced flag;
        ``ref_sf_h_MHz`` (the session's 1H reference SF) lets the SF be
        checked against ``nucleus``'s expected frequency."""
        err = meas.sigma_ppm if dcg_err_ppm is None else dcg_err_ppm
        err = max(float(err), ERR_FLOOR_PPM) if np.isfinite(err) else 0.0
        conv = getattr(meas, "convergence", None)
        pt = cls(float(larmor_MHz),
                 float(meas.cg_ppm if dcg_ppm is None else dcg_ppm), err,
                 ct_selective, label, magnitude=magnitude,
                 window=tuple(meas.window), window_mode=meas.mode,
                 source=source, rotor_Hz=float(rotor_Hz or 0.0),
                 flags=tuple(meas.flags),
                 cg_sequence=conv.sequence() if conv is not None else "",
                 drift_ppm=float(meas.drift_ppm))
        pt.fill_referencing(meta, ref_sf_h_MHz, nucleus)
        return pt

    def fill_referencing(self, meta: dict | None, ref_sf_h_MHz: float | None = None,
                         nucleus: str = "") -> None:
        """LB, SF, SR and the referenced flag from ``meta``; the deviation
        of SF from a session 1H reference when one is given."""
        meta = meta or {}
        for key, attr in (("lb_Hz", "lb_Hz"), ("qcpmg_lb_Hz", "lb_Hz"),
                          ("sf_MHz", "sf_MHz"), ("sr_hz", "sr_hz")):
            v = meta.get(key)
            if v is not None and getattr(self, attr) is None:
                try:
                    setattr(self, attr, float(v))
                except (TypeError, ValueError):
                    pass
        if meta.get("referenced") is not None:
            self.referenced = bool(meta["referenced"])
        elif self.sr_hz is not None:
            from larmor.referencing import UNREFERENCED_HZ
            self.referenced = abs(self.sr_hz) > UNREFERENCED_HZ
        if ref_sf_h_MHz and self.sf_MHz and nucleus:
            self.ref_dev_ppm = reference_deviation_ppm(self.sf_MHz, nucleus, ref_sf_h_MHz)

    def provenance(self) -> str:
        """One line: 'window -206.8 … -35.4 ppm (minima) · magnitude · LB 75 Hz
        · SR +3982.9 Hz · MAS 16000 Hz · LAW0Ca-3Cl_850_MHz.csv'."""
        bits = []
        if self.window is not None:
            lo, hi = min(self.window), max(self.window)
            bits.append(f"window {lo:.1f} ... {hi:.1f} ppm"
                        + (f" ({self.window_mode})" if self.window_mode else ""))
        if self.magnitude is not None:
            bits.append("magnitude (mc)" if self.magnitude else "absorption")
        if self.lb_Hz is not None:
            bits.append(f"LB {self.lb_Hz:.0f} Hz")
        if self.sr_hz is not None:
            bits.append(f"SR {self.sr_hz:+.1f} Hz")
        if self.rotor_Hz:
            bits.append(f"MAS {self.rotor_Hz:.0f} Hz")
        if self.source:
            bits.append(self.source)
        return " · ".join(bits)


#: |SF - expected| beyond this (ppm) against the session's 1H reference is a
#: referencing error (-0.043 ppm on every referenced NMRFAM EXPNO passes,
#: the SF = BF1 EXPNO's -2.49 ppm fails)
REF_TOLERANCE_PPM = 0.5


def reference_deviation_ppm(sf_MHz: float, nucleus: str, sf_h_MHz: float) -> float:
    """(SF - expected)/expected in ppm, where expected is the frequency the
    nucleus should have on the same magnet given the 1H reference SF
    (larmor.referencing.expected_sf_MHz, IUPAC Xi ratios)."""
    from larmor.referencing import expected_sf_MHz
    exp = expected_sf_MHz(float(sf_h_MHz), nucleus)
    return float((float(sf_MHz) - exp) / exp * 1e6)


def referencing_checks(points) -> list[str]:
    """Report lines for the per-field referencing: a hard warning for an
    unreferenced axis (|SR| below UNREFERENCED_HZ) or an SF that deviates
    from the session 1H reference by more than REF_TOLERANCE_PPM; an
    informational line when SR is known but no 1H reference was given."""
    from larmor.referencing import UNREFERENCED_HZ
    out = []
    for p in points:
        if p.sr_hz is None and p.referenced is None:
            continue
        unref = (p.referenced is False) or (p.sr_hz is not None
                                            and abs(p.sr_hz) <= UNREFERENCED_HZ)
        if unref:
            out.append(f"! {p.larmor_MHz:.4f} MHz: SR = "
                       f"{(p.sr_hz if p.sr_hz is not None else 0.0):+.1f} Hz -- "
                       "unreferenced ppm axis; a rigid offset enters delta_iso "
                       "with the field's lever (x-1.1 low / x+2.1 high field)")
        elif p.ref_dev_ppm is not None:
            if abs(p.ref_dev_ppm) > REF_TOLERANCE_PPM:
                out.append(f"! {p.larmor_MHz:.4f} MHz: SF deviates "
                           f"{p.ref_dev_ppm:+.2f} ppm from the session 1H reference")
            else:
                out.append(f"referencing {p.larmor_MHz:.4f} MHz: SF within "
                           f"{p.ref_dev_ppm:+.3f} ppm of the session 1H reference")
        else:
            out.append(f"referencing {p.larmor_MHz:.4f} MHz: SR {p.sr_hz:+.1f} Hz "
                       "(not checked against a 1H reference -- enter the session's "
                       "1H SF to check)")
    return out




@dataclass
class InfiniteFieldResult:
    delta_iso_ppm: float
    delta_iso_err_ppm: float          # NaN when nothing could be propagated
    cq_MHz: float                     # 0.0 when the slope is not significantly < 0
    cq_err_MHz: float                 # 0.0 in that case (see cq_upper_2sigma_MHz)
    pq_MHz: float
    eta: float
    spin: float
    slope: float                      # ppm·MHz² (δcg vs 1/ν0²)
    intercept: float                  # == delta_iso_ppm
    points: list[FieldPoint] = field(default_factory=list)
    #: σ(P_Q) = σ(C_Q)·sqrt(1+η²/3) = P_Q·σ_b/(2|b|): the η-FREE uncertainty
    pq_err_MHz: float = 0.0
    #: C_Q at η = 0 and η = 1 for the same slope (+7.9 % / -6.6 % about η 0.7):
    #: the systematic the assumed η adds, kept OUT of the +-
    cq_eta_range: tuple = (0.0, 0.0)
    slope_err: float = float("nan")   # σ_b (ppm·MHz²)
    #: 2-σ upper bound on C_Q when the slope is not significantly negative
    cq_upper_2sigma_MHz: float = 0.0
    #: non-empty when C_Q could not be quoted as a value (positive or
    #: non-significant slope) or when the uncertainties were not propagated
    note: str = ""
    #: non-fatal quality warnings (short lever arm, poorly constrained)
    warning: str = ""
    lever_arm: float = 0.0            # (x.max - x.min) / x.max in 1/ν0²
    weighted: bool = True             # False: equal weights (a sigma was missing)
    # ---- goodness of fit (n >= 3 has redundancy; n == 2 is exact) ----------
    residuals_ppm: tuple = ()         # y - a - b x per point
    chi2: float = 0.0                 # sum w (y - a - b x)^2
    dof: int = 0                      # n - 2
    chi2_red: float = float("nan")    # chi2 / dof (NaN when dof == 0)
    p_value: float = float("nan")     # survival probability of chi2 at dof
    #: +- scaled by sqrt(chi2_red) when the a-priori errors under-predict the
    #: scatter (never scaled down); equal to the a-priori values otherwise
    delta_iso_err_scaled_ppm: float = float("nan")
    cq_err_scaled_MHz: float = float("nan")
    #: True when no input sigma was usable and the +- come from the scatter
    errors_from_scatter: bool = False

    def line(self, inv_nu2: np.ndarray) -> np.ndarray:
        """δcg on the fit line for given 1/ν0² values (for plotting)."""
        return self.intercept + self.slope * np.asarray(inv_nu2, float)

    @property
    def cq_is_bound(self) -> bool:
        """True when C_Q is reported as an upper bound, not a value."""
        return self.cq_upper_2sigma_MHz > 0.0 and self.cq_MHz == 0.0

    @property
    def scaled(self) -> bool:
        """True when the scaled +- differ from the a-priori ones (chi2_red > 1
        with real input errors)."""
        return (np.isfinite(self.delta_iso_err_scaled_ppm)
                and np.isfinite(self.delta_iso_err_ppm)
                and self.delta_iso_err_scaled_ppm > self.delta_iso_err_ppm * (1 + 1e-9))

    @property
    def misfit(self) -> bool:
        """True when the scatter about the line is improbable (p < 0.01)."""
        return bool(np.isfinite(self.p_value) and self.p_value < 0.01)


def _quad_A(spin: float, eta: float) -> float:
    """A in slope = −A·C_Q² (ppm·MHz² per MHz²)."""
    return (1.0e6 / 40.0) * (3.0 + eta ** 2) * _spin_factor(spin)


def weighted_line(x, y, err) -> tuple[float, float, np.ndarray]:
    """Weighted least-squares line y = slope·x + intercept with weights
    1/err², in closed form. Returns ``(slope, intercept, cov)`` with ``cov``
    the 2×2 parameter covariance ``[[var(slope), cov], [cov, var(intercept)]]``
    propagated from ``err`` (unscaled -- the caller decides whether to
    inflate it by the reduced χ²). Raises ValueError('degenerate abscissae')
    when every x is the same.

    The DFT shielding calibration (larmor.shiftcal) calls it; the QCPMG
    extrapolation (:func:`infinite_field_diso`) carries the same algebra
    inline because it adds the σ policy (equal weights when a σ is missing,
    the ERR_FLOOR_PPM floor) and the χ² scaling around it. The denominator
    is the sum over pairs Σ w_i w_j (x_i − x_j)², which is free of the
    sw·sxx − sx² cancellation; tests/test_shiftcal.py pins that the two
    agree to rounding.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    err = np.asarray(err, float)
    w = 1.0 / err ** 2

    # weighted linear fit y = a + b x  (a = intercept, b = slope)
    sw = w.sum()
    sx = (w * x).sum(); sy = (w * y).sum()
    sxx = (w * x * x).sum(); sxy = (w * x * y).sum()
    # the denominator as a sum over pairs (algebraically sw*sxx - sx*sx)
    denom = 0.0
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            denom += w[i] * w[j] * (x[i] - x[j]) ** 2
    if not denom > 0.0:
        raise ValueError("degenerate abscissae")
    b = (sw * sxy - sx * sy) / denom
    a = (sy - b * sx) / sw
    # parameter variances from the weighted fit
    var_a = sxx / denom
    var_b = sw / denom
    cov_ab = -sx / denom
    return float(b), float(a), np.array([[var_b, cov_ab], [cov_ab, var_a]])


def infinite_field_diso(points: list[FieldPoint], spin: float,
                        eta: float = DEFAULT_ETA) -> InfiniteFieldResult:
    """Fit δcg = δiso + slope·(1/ν0²) across fields and return δiso, C_Q, P_Q.

    Needs ≥ 2 fields. With exactly 2 the line is exact (errors from the δcg
    uncertainties are propagated); with > 2 a weighted least-squares line is
    fit.

    Uncertainty policy. A point whose ``dcg_err_ppm`` is 0, negative or NaN
    has NO known uncertainty. The fit is then UNWEIGHTED (equal weights, said
    so in ``note``): with two points nothing can be propagated and the ± are
    NaN; with three or more the ± come from the scatter about the line
    (ordinary least squares scaled by χ²/(n−2)). Nothing is ever invented
    for a missing σ. Valid σ are floored at ``ERR_FLOOR_PPM``.

    C_Q policy. The slope is converted through C_Q² = −b/A BEFORE the square
    root, with σ(C_Q²) = σ_b/A. When C_Q² ≥ 2σ the value is quoted as usual;
    when the slope is within ±2σ of zero C_Q is reported as a 2-σ upper
    bound (``cq_upper_2sigma_MHz``) with ``cq_MHz = cq_err_MHz = 0.0`` and a
    ``note``; when the slope is significantly POSITIVE Eq. (1) cannot have
    produced it and the note says the two δcg are not comparable.
    """
    if len(points) < 2:
        raise ValueError("need the centre of gravity at at least two fields")
    if any(not (float(p.larmor_MHz) > 0.0) for p in points):
        raise ValueError(
            "Larmor frequency must be > 0 MHz (the observe frequency of THIS "
            "nucleus, not the 1H frequency)")
    n = len(points)
    x = np.array([1.0 / float(p.larmor_MHz) ** 2 for p in points])  # 1/ν0² (MHz⁻²)
    y = np.array([float(p.dcg_ppm) for p in points])
    err = np.array([float(p.dcg_err_ppm) for p in points])
    valid = np.isfinite(err) & (err > 0.0)

    notes: list[str] = []
    warnings: list[str] = []
    weighted = bool(valid.all())
    if weighted:
        w = 1.0 / np.maximum(err, ERR_FLOOR_PPM) ** 2
    else:
        w = np.ones(n)
        rows = ", ".join(str(i + 1) for i in np.where(~valid)[0])
        notes.append("uncertainties not propagated: equal weights (row "
                     f"{rows} has no/invalid sigma)")

    xr = float((x.max() - x.min()) / x.max())
    if xr < MIN_LEVER_ARM:
        raise ValueError(
            "fields differ by less than 5 % in 1/nu0^2 (Larmor frequencies "
            f"within ~2.5 %; lever arm {xr:.4f}) -- cannot extrapolate to "
            "infinite field")

    # weighted linear fit y = a + b x  (a = δiso, b = slope); the denominator
    # as a sum over pairs is free of the sw*sxx - sx*sx cancellation
    sw = w.sum()
    sx = (w * x).sum(); sy = (w * y).sum()
    sxx = (w * x * x).sum(); sxy = (w * x * y).sum()
    denom = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            denom += w[i] * w[j] * (x[i] - x[j]) ** 2
    if not denom > 1e-12 * sw * sxx:
        raise ValueError(
            f"the fields are too close in 1/nu0^2 to extrapolate (lever arm {xr:.4f})")
    b = (sw * sxy - sx * sy) / denom
    a = (sy - b * sx) / sw
    resid = y - a - b * x
    chi2 = float((w * resid ** 2).sum())
    dof = n - 2
    # parameter variances of the weighted fit (unit-variance weights)
    var_a = max(sxx / denom, 0.0)
    var_b = max(sw / denom, 0.0)
    chi2_red = chi2 / dof if dof > 0 else float("nan")
    p_value = float("nan")
    if dof > 0 and weighted:
        from scipy.stats import chi2 as _chi2_dist
        p_value = float(_chi2_dist.sf(chi2, dof))
    errors_from_scatter = False
    scale = 1.0
    if weighted:
        sig_a, sig_b = float(np.sqrt(var_a)), float(np.sqrt(var_b))
        if dof > 0 and chi2_red > 1.0:
            scale = float(np.sqrt(chi2_red))          # never scaled down
    elif dof > 0:
        s2 = chi2 / dof                       # OLS residual variance
        sig_a, sig_b = float(np.sqrt(var_a * s2)), float(np.sqrt(var_b * s2))
        notes.append("no input errors given: +- from scatter about the line")
        errors_from_scatter = True
    else:
        sig_a = sig_b = float("nan")          # exact line, nothing to estimate

    A = _quad_A(spin, eta)
    cq2 = -b / A
    sig_cq2 = sig_b / A if np.isfinite(sig_b) else float("nan")
    cq = cq_err = cq_upper = 0.0
    if np.isfinite(sig_cq2) and cq2 < 2.0 * sig_cq2:
        cq_upper = float(np.sqrt(max(cq2 + 2.0 * sig_cq2, 0.0)))
        if cq2 > -2.0 * sig_cq2:
            notes.append(
                f"slope not significantly negative: C_Q <= {cq_upper:.3f} MHz "
                "(2-sigma upper bound) -- check that both dcg are the same "
                "observable (mode, window, referencing)")
        else:
            notes.append(
                "slope significantly POSITIVE: Eq. 1 cannot produce this; the "
                "two dcg are not comparable (mode, window or referencing differ)")
    elif cq2 <= 0.0:
        # no σ to judge significance (two points, no uncertainties) and the
        # slope is not negative: Eq. (1) cannot have produced it
        notes.append(
            "slope is positive or zero: Eq. 1 cannot produce this; the two "
            "dcg are not comparable (mode, window or referencing differ)")
    else:
        cq = float(np.sqrt(cq2))
        # dC_Q/d(C_Q²) = 1/(2 C_Q)  →  σ_Cq = σ(C_Q²)/(2 C_Q) = |C_Q/(2b)|·σ_b
        cq_err = float(sig_cq2 / (2.0 * cq)) if np.isfinite(sig_cq2) else float("nan")
    eta_factor = float(np.sqrt(1.0 + eta ** 2 / 3.0))
    pq = cq * eta_factor
    pq_err = cq_err * eta_factor if np.isfinite(cq_err) else float("nan")
    cq_eta_range = (cq_from_slope(float(b), spin, 0.0), cq_from_slope(float(b), spin, 1.0))

    if xr < SHORT_LEVER_ARM:
        warnings.append(
            f"short lever arm ({xr:.2f} of 1/nu0^2): the intercept is an "
            "extrapolation over more than 5x the spanned range")
    if valid.any() and np.isfinite(sig_a) and sig_a > 10.0 * float(err[valid].max()):
        warnings.append(
            f"poorly constrained: sigma(delta_iso) = {sig_a:.1f} ppm is more "
            "than 10x the largest dcg uncertainty")
    sep = float(y.max() - y.min())
    for p in points:
        d = float(getattr(p, "drift_ppm", float("nan")))
        if np.isfinite(d) and sep > 0 and d > 0.10 * sep:
            warnings.append(
                f"CG drift {d:.1f} ppm at {p.larmor_MHz:.3f} MHz is {100 * d / sep:.0f} % "
                f"of the {sep:.1f} ppm separation between the fields -- the window "
                "cuts the pattern")

    if dof > 0 and weighted and np.isfinite(p_value) and p_value < 0.01:
        warnings.append(
            f"scatter about the line is improbable (chi2/dof = {chi2:.3g}/{dof}, "
            f"p = {p_value:.2g}): an outlier, a window that caught a sideband or "
            "a mode mismatch -- the scaled +- are quoted alongside")
    return InfiniteFieldResult(
        delta_iso_ppm=float(a), delta_iso_err_ppm=sig_a,
        cq_MHz=cq, cq_err_MHz=cq_err, pq_MHz=pq, eta=eta, spin=spin,
        slope=float(b), intercept=float(a), points=list(points),
        pq_err_MHz=pq_err, cq_eta_range=cq_eta_range,
        slope_err=sig_b, cq_upper_2sigma_MHz=cq_upper,
        note="; ".join(notes), warning="; ".join(warnings),
        lever_arm=xr, weighted=weighted,
        residuals_ppm=tuple(float(v) for v in resid), chi2=chi2, dof=dof,
        chi2_red=chi2_red, p_value=p_value,
        delta_iso_err_scaled_ppm=sig_a * scale if np.isfinite(sig_a) else float("nan"),
        cq_err_scaled_MHz=cq_err * scale if np.isfinite(cq_err) else float("nan"),
        errors_from_scatter=errors_from_scatter)


#: plausible static-field range for a solid-state NMR magnet (tesla): 4 T
#: (~170 MHz 1H) to 36 T (the 1.5 GHz series-connected hybrid)
B0_PLAUSIBLE_T = (4.0, 36.0)


def implied_B0_T(nu_MHz: float, nucleus: str) -> float | None:
    """Static field (T) that ``nu_MHz`` implies for ``nucleus``, or None for a
    blank/unknown nucleus. Uses |γ| (mrsimulator's gyromagnetic_ratio is
    SIGNED: 29Si is −8.4655 MHz/T)."""
    if not nucleus or not (float(nu_MHz or 0.0) > 0.0):
        return None
    try:
        from mrsimulator.spin_system.isotope import ISOTOPE_DATA
        d = ISOTOPE_DATA.get(str(nucleus))
    except Exception:                                         # noqa: BLE001
        return None
    if not d:
        return None
    gamma = abs(float(d.get("gyromagnetic_ratio", 0.0) or 0.0))
    return float(nu_MHz) / gamma if gamma > 0 else None


def field_plausibility_warning(nu_MHz: float, nucleus: str) -> str:
    """'' when ``nu_MHz`` puts ``nucleus`` in a plausible magnet, else a
    one-line warning naming the implied field -- the 1H frequency (800 /
    1100 MHz) typed for 35Cl would imply 200+ T and a C_Q ~10x too large."""
    b0 = implied_B0_T(nu_MHz, nucleus)
    if b0 is None:
        return ""
    lo, hi = B0_PLAUSIBLE_T
    if lo <= b0 <= hi:
        return ""
    hint = (" -- did you enter the 1H frequency?" if b0 > hi
            else " -- is this the observe frequency in MHz (not GHz)?")
    return (f"nu0 = {float(nu_MHz):g} MHz would put {nucleus} at {b0:.1f} T"
            + hint)


#: Gaussian FWHM / sigma
GAUSS_FWHM_PER_SIGMA = 2.3548200450309493


@dataclass
class WidthSplit:
    wq_lo_ppm: float          # quadrupolar width at the LOWEST field (ppm)
    wq_hi_ppm: float          # quadrupolar width at the highest field (ppm)
    #: FIELD-INDEPENDENT width (ppm): everything constant in ppm -- the
    #: distribution of isotropic shifts AND the chemical-shift anisotropy --
    #: so an upper bound on shift disorder unless the CSA is known small.
    #: (The name keeps the token 'csd' for API stability.)
    wcsd_ppm: float
    ok: bool                  # False when the split is unphysical (see note)
    note: str = ""
    n_fields: int = 2
    fields_MHz: tuple = ()
    chi2: float = float("nan")          # sum w (FWHM_i^2 - model)^2
    residuals_ppm2: tuple = ()          # FWHM_i^2 - model_i (ppm^2)
    weighted: bool = False              # 1/sigma^2 weights were used
    wq_err_ppm: float = float("nan")    # at the lowest field
    wcsd_err_ppm: float = float("nan")
    gate: str = ""                      # non-empty in the spurious regime
    kind: str = "fwhm"                  # 'fwhm' (Sandland Eq. 2) | 'second_moment'

    @property
    def label(self) -> str:
        return ("W_csd" if self.kind == "fwhm" else
                "W_csd (Gaussian-equivalent FWHM = 2.3548 sigma_csd)")


#: |FWHM_hi/FWHM_lo / (nu_lo/nu_hi)^2 - 1| below this: the widths scale as a
#: pure quadrupolar pattern and W_csd is not resolved
WIDTH_GATE_PURE_QUAD = 0.15
#: W_csd below this fraction of W_q(low field) is in the regime where the
#: quadrature model returns spurious values (10-17 ppm on a pure CT pattern)
WIDTH_GATE_CSD_FRACTION = 0.30


def _width_split(points, nu_ref, lb_Hz, kind: str) -> WidthSplit:
    """Shared solver for the FWHM and second-moment splits.

    ``points``: (nu_MHz, width_ppm[, sigma_width_ppm]); the model is
    width_i^2 = Wq_ref^2 (nu_ref/nu_i)^4 + Wcsd^2, linear in the two
    squares. n == 2 is solved exactly (the Sandland closed form); n > 2 by
    least squares (1/(2 W sigma)^2 weights when every point has a sigma),
    with scipy's nnls when the unconstrained solution goes negative. A
    processing line broadening ``lb_Hz`` (one value or one per point) is
    removed in quadrature, in ppm, at each field first.
    """
    pts = []
    for i, p in enumerate(points):
        nu, w = float(p[0]), float(p[1])
        sig = float(p[2]) if len(p) > 2 and p[2] is not None else float("nan")
        lb = 0.0
        if lb_Hz is not None:
            lb = float(lb_Hz[i] if isinstance(lb_Hz, (list, tuple)) else lb_Hz)
        pts.append((nu, w, sig, lb))
    if len(pts) < 2:
        return WidthSplit(0, 0, 0, False, "need the width at two or more fields",
                          n_fields=len(pts), kind=kind)
    pts.sort(key=lambda t: t[0])
    nus = np.array([t[0] for t in pts]); ws = np.array([t[1] for t in pts])
    sigs = np.array([t[2] for t in pts]); lbs = np.array([t[3] for t in pts])
    if not np.all(nus > 0) or not np.all(np.isfinite(ws)):
        return WidthSplit(0, 0, 0, False, "widths or fields are not finite",
                          n_fields=len(pts), fields_MHz=tuple(nus), kind=kind)
    nu_ref = float(nu_ref) if nu_ref else float(nus[0])
    w2 = ws ** 2 - (lbs / nus) ** 2                    # LB removed in quadrature
    x = (nu_ref / nus) ** 4
    if np.ptp(x) < 1e-9 * max(1.0, float(np.max(x))):
        return WidthSplit(0, 0, 0, False, "the fields are too close",
                          n_fields=len(pts), fields_MHz=tuple(nus), kind=kind)
    weighted = bool(np.all(np.isfinite(sigs) & (sigs > 0)))
    wts = 1.0 / (2.0 * ws * sigs) ** 2 if weighted else np.ones(ws.size)
    A = np.column_stack([x, np.ones_like(x)])
    sw = np.sqrt(wts)
    sol, *_ = np.linalg.lstsq(A * sw[:, None], w2 * sw, rcond=None)
    a, c = float(sol[0]), float(sol[1])                # Wq_ref^2, Wcsd^2
    ok = a >= 0.0 and c >= 0.0
    if not ok and ws.size > 2:
        from scipy.optimize import nnls
        sol_nn, _ = nnls(A * sw[:, None], w2 * sw)
        a, c = float(sol_nn[0]), float(sol_nn[1])
    a, c = max(a, 0.0), max(c, 0.0)
    model = a * x + c
    resid = w2 - model
    chi2 = float((wts * resid ** 2).sum())
    # parameter errors from the (weighted) normal equations
    try:
        cov = np.linalg.inv((A * wts[:, None]).T @ A)
        if not weighted and ws.size > 2:
            cov = cov * chi2 / (ws.size - 2)
        var_a, var_c = max(float(cov[0, 0]), 0.0), max(float(cov[1, 1]), 0.0)
    except np.linalg.LinAlgError:
        var_a = var_c = float("nan")
    wq = float(np.sqrt(a)); wcsd = float(np.sqrt(c))
    wq_err = float(np.sqrt(var_a) / (2 * wq)) if wq > 0 and (weighted or ws.size > 2) else float("nan")
    wcsd_err = float(np.sqrt(var_c) / (2 * wcsd)) if wcsd > 0 and (weighted or ws.size > 2) else float("nan")
    wq_hi = wq * (nu_ref / float(nus[-1])) ** 2
    note = "" if ok else (
        "widths do not separate: the higher-field line is broader than the "
        "quadrupolar model allows — check the FWHM values or the CT band.")
    gate = ""
    ratio = float(ws[-1] / ws[0]) if ws[0] > 0 else float("nan")
    pure = (float(nus[0]) / float(nus[-1])) ** 2
    if np.isfinite(ratio) and abs(ratio / pure - 1.0) < WIDTH_GATE_PURE_QUAD:
        gate = (f"widths scale as a pure quadrupolar pattern (FWHM ratio {ratio:.3f} "
                f"vs (nu_lo/nu_hi)^2 = {pure:.3f}): W_csd is not resolved")
    elif ok and wq > 0 and wcsd < WIDTH_GATE_CSD_FRACTION * wq:
        gate = (f"W_csd = {wcsd:.1f} ppm is below 0.3 W_q ({wq:.1f} ppm): in this "
                "regime the quadrature model returns spurious W_csd (10-17 ppm on "
                "a pure CT pattern) -- treat it as an upper bound")
    return WidthSplit(wq, wq_hi, wcsd, ok, note, n_fields=int(ws.size),
                      fields_MHz=tuple(float(v) for v in nus), chi2=chi2,
                      residuals_ppm2=tuple(float(v) for v in resid),
                      weighted=weighted, wq_err_ppm=wq_err, wcsd_err_ppm=wcsd_err,
                      gate=gate, kind=kind)


def multi_field_widths(points, nu_ref: float | None = None,
                       lb_Hz=None) -> WidthSplit:
    """Separate the CT linewidth into a quadrupolar part W_q (∝ 1/ν0²,
    broader at low field) and a FIELD-INDEPENDENT part W_csd (constant in
    ppm: the distribution of isotropic shifts AND the chemical-shift
    anisotropy) from the FWHM measured at two or more fields (Sandland et
    al. 2004, Eq. 2, generalised to N fields):

        FWHM_i² = W_q,ref² (ν_ref/ν_i)⁴ + W_csd²

    ``points``: (nu_MHz, fwhm_ppm[, sigma_fwhm_ppm]). Eq. 2 is exact only for
    Gaussian-convolved lines: a pure second-order CT pattern returns a
    non-zero W_csd (10-17 ppm static) and a CSA span comparable to W_q
    inflates both -- the ``gate`` field says when the split is in that
    regime. Prefer :func:`second_moment_split` where the whole band is
    measurable.
    """
    return _width_split(points, nu_ref, lb_Hz, "fwhm")


def second_moment_split(points, nu_ref: float | None = None) -> WidthSplit:
    """The same separation on the SECOND MOMENT (intensity-weighted standard
    deviation, :func:`larmor.qcpmg.second_moment_ppm`) instead of the FWHM:

        σ_i² = σ_q,ref² (ν_ref/ν_i)⁴ + σ_csd²

    Variances add exactly under convolution, so this split needs no
    Gaussian assumption (injected 0 / 0.3 / 30 / 60 ppm shift distributions
    come back as 0.1 / 0.3 / 30.0 / 60.0 on static CT patterns). It is
    window-sensitive -- the window must hold the whole CT band and exclude
    spinning sidebands, which are fixed in Hz. The result is reported as the
    Gaussian-EQUIVALENT FWHM = 2.3548 σ, labelled as such.
    ``points``: (nu_MHz, sigma_ppm[, err_ppm]).
    """
    pts = [(p[0], GAUSS_FWHM_PER_SIGMA * float(p[1]),
            (GAUSS_FWHM_PER_SIGMA * float(p[2]) if len(p) > 2 and p[2] is not None
             else None)) for p in points]
    return _width_split(pts, nu_ref, None, "second_moment")


def two_field_widths(nu1_MHz: float, fwhm1_ppm: float,
                     nu2_MHz: float, fwhm2_ppm: float) -> WidthSplit:
    """Two-field closed form of :func:`multi_field_widths` (Sandland Eq. 2):
    with the lower field = 1,
        FWHM1² = W_q1² + W_csd²
        FWHM2² = W_q1²·(ν1/ν2)⁴ + W_csd²
    → W_q1² = (FWHM1² − FWHM2²)/(1 − (ν1/ν2)⁴),  W_csd² = FWHM1² − W_q1².
    Field order does not matter; an unphysical split is flagged ok=False.
    """
    return multi_field_widths([(nu1_MHz, fwhm1_ppm), (nu2_MHz, fwhm2_ppm)])


def centre_of_gravity(ppm: np.ndarray, amp: np.ndarray,
                      lo_ppm: float | None = None,
                      hi_ppm: float | None = None) -> float:
    """Intensity-weighted centre of gravity of a spectrum over an optional
    ppm window (restrict it to the central-transition band).

    A thin wrapper on :func:`larmor.qcpmg.centre_of_gravity` -- the SIGNED
    estimator the dataset rows and the batch grid use -- so every route into
    the fit shares one definition. The former ``np.clip(amp, 0)`` turned
    zero-mean noise into a positive pedestal that pulled the CG towards the
    window centre (bias -3.7 ppm on a Gaussian at S/N 20 over a wide window
    against -0.1 ppm signed). With no window the first-minima window of
    :func:`larmor.qcpmg.cg_window` is used.
    """
    from larmor import qcpmg
    window = None if lo_ppm is None or hi_ppm is None else (float(lo_ppm), float(hi_ppm))
    return float(qcpmg.centre_of_gravity(np.asarray(ppm, float),
                                         np.asarray(amp, float), window)[0])


def point_provenance(p: FieldPoint) -> dict:
    """A JSON-able record of one point for the figure spec sidecar."""
    return {"larmor_MHz": p.larmor_MHz, "dcg_ppm": p.dcg_ppm,
            "dcg_err_ppm": p.dcg_err_ppm,
            "window": list(p.window) if p.window is not None else None,
            "window_mode": p.window_mode, "magnitude": p.magnitude,
            "source": p.source, "rotor_Hz": p.rotor_Hz, "lb_Hz": p.lb_Hz,
            "sf_MHz": p.sf_MHz, "sr_hz": p.sr_hz, "referenced": p.referenced,
            "ct_selective": p.ct_selective, "flags": list(p.flags),
            "cg_sequence": p.cg_sequence}


_FIELD_TOKEN_RE = re.compile(r"[_\-\s]+\d+(?:p\d+)?_?(?:MHz|GHz|T)$", re.IGNORECASE)
_SUMECHO_SUFFIX_RE = re.compile(r"\s*(?:·|\||-)?\s*QCPMG sum echo.*$")
_DATE_ONLY_RE = re.compile(r"\s*\d{1,4}[/.-]\d{1,2}[/.-]\d{2,4}\s*")


def sample_label(meta: dict | None, path: str | Path | None) -> str:
    """A sample name for a dataset: the header's sample line with the
    ' · QCPMG sum echo (...)' suffix and a leading date token removed
    (qcpmg.sample_name), else the file stem with its trailing field token
    removed ('LAW0Ca-3Cl_850_MHz' -> 'LAW0Ca-3Cl')."""
    from larmor.qcpmg import sample_name

    raw = str((meta or {}).get("sample", "") or "")
    raw = _SUMECHO_SUFFIX_RE.sub("", raw).strip()
    name = sample_name(raw) if raw else ""
    if not name and (meta or {}).get("title"):
        name = sample_name(str(meta["title"]))          # a Bruker title
    if _DATE_ONLY_RE.fullmatch(name or ""):
        name = ""                     # a title that is only a date names nothing
    if not name and path:
        name = _FIELD_TOKEN_RE.sub("", Path(path).stem).strip()
    return name


def spectrum_mode_from_meta(meta: dict) -> bool | None:
    """True (magnitude) / False (absorption) / None (unknown) from a spectrum's
    metadata: the ``spectrum_mode`` header a saved dataset carries, else the
    legacy ', magnitude)' suffix of its sample line, else procs ``ph_mod``
    (2 = mc) of a Bruker 1r."""
    mode = str(meta.get("spectrum_mode", "") or "").lower()
    if mode.startswith("magnitude"):
        return True
    if mode.startswith("absorption"):
        return False
    sample = str(meta.get("sample", "") or "")
    if "QCPMG sum echo" in sample:
        return "magnitude" in sample.lower()
    if meta.get("ph_mod") is not None:
        try:
            return int(meta["ph_mod"]) == 2
        except (TypeError, ValueError):
            return None
    return None


def read_field_spectrum(path: str) -> dict:
    """Load one field's spectrum for either multi-field dialog.

    Returns ``{"ppm", "amp", "larmor", "nucleus", "magnitude", "source",
    "meta", "seed"}``. A LARMOR .csv (the sum-echo dataset *Save as
    dataset...* writes) supplies the Larmor frequency, nucleus and processing
    mode from its header; a Bruker 1r supplies them from acqus/procs (mode
    from PH_mod) and is checked for a spikelet comb, whose window is then
    seeded from the envelope (:func:`larmor.qcpmg.seed_window`).
    """
    from larmor import qcpmg
    from larmor.io import bruker

    try:
        ref = bruker.resolve(path)
    except (ValueError, FileNotFoundError):
        ref = None
    if ref is not None:
        d = bruker.read(path)
        if d.ndim != 1 or d.domain != "freq":
            raise ValueError("not a processed 1D spectrum")
        ppm = np.asarray(d.axes[0].values, float)
        amp = np.asarray(d.data, float)
        meta = dict(d.meta)
        larmor = float(meta.get("larmor_MHz", 0.0) or 0.0)
        nucleus = str(meta.get("nucleus", "") or "").strip()
        source = f"Bruker 1r {ref.expno}/pdata/{ref.procno}"
    else:
        from larmor.io import spectra
        from larmor.loader import load_any

        ppm, amp, recipe, _summary, _warn = load_any(path)
        ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
        larmor = float(recipe.get("larmor_frequency_MHz") or 0.0)
        nucleus = str(recipe.get("nucleus") or "").strip()
        meta = {}
        if str(recipe.get("source_kind", "")) == "csv":
            _, _, meta = spectra.read_csv(path)
        source = f"dataset {Path(path).name}"
    if not larmor:
        raise ValueError("no Larmor frequency in the file — save it from "
                         "the QCPMG dialog, which records one")
    rotor, rotor_note = rotor_rate_of(meta, ppm, amp, larmor)
    if ref is not None and meta.get("sf_MHz") is not None and "referenced" not in meta:
        from larmor.referencing import UNREFERENCED_HZ
        meta["referenced"] = abs(float(meta.get("sr_hz", 0.0) or 0.0)) > UNREFERENCED_HZ
    return {"ppm": ppm, "amp": amp, "larmor": larmor, "nucleus": nucleus,
            "magnitude": spectrum_mode_from_meta(meta), "source": source,
            "meta": meta, "seed": qcpmg.seed_window(ppm, amp, meta),
            "rotor_Hz": rotor, "rotor_note": rotor_note}


def rotor_rate_of(meta: dict, ppm=None, amp=None, larmor_MHz: float = 0.0
                  ) -> tuple[float, str]:
    """(rotor_Hz, note): the MAS rate a spectrum was acquired at -- a saved
    dataset's ``qcpmg_rotor_Hz`` header (its ``spin_rate_Hz`` is 0 by design,
    see QcpmgDialog.dataset_meta), else a Bruker meta's resolved
    ``spin_rate_Hz``. When the trace shows a sideband repeat
    (larmor.sidebands.detect) that disagrees with the recorded rate by more
    than 5 % the note says so; ``mas_uncertain`` is passed on as a note."""
    rate = 0.0
    if meta.get("qcpmg_rotor_Hz") is not None:
        rate = float(meta.get("qcpmg_rotor_Hz") or 0.0)
    elif meta.get("spin_rate_Hz") is not None and "larmor_MHz" in meta:
        rate = float(meta.get("spin_rate_Hz") or 0.0)       # Bruker meta
    notes = []
    if rate > 0 and meta.get("mas_uncertain"):
        notes.append(f"rotor rate {rate:.0f} Hz is flagged uncertain "
                     "(acqus MASR and the title disagree)")
    if rate > 0 and ppm is not None and amp is not None and larmor_MHz:
        try:
            from larmor import sidebands
            det = sidebands.detect(ppm, amp, larmor_MHz, rate)
            if det.ok and abs(det.nu_rot_Hz - rate) > 0.05 * rate:
                notes.append(f"sideband spacing measured {det.nu_rot_Hz:.0f} Hz "
                             f"vs recorded {rate:.0f} Hz")
        except Exception:                                     # noqa: BLE001
            pass
    return rate, "; ".join(notes)


def fit_samples(rows, spin: float, eta: float = DEFAULT_ETA
                ) -> dict[str, "InfiniteFieldResult | str"]:
    """Extrapolate SEVERAL samples at once.

    ``rows`` is an iterable of ``(sample, FieldPoint)``. Returns
    ``{sample: InfiniteFieldResult}``, with a plain string in place of the
    result for any sample that could not be fitted (fewer than two fields,
    two fields too close) -- one bad sample must never sink the batch.
    Insertion order is preserved so a report reads in the order entered.
    """
    grouped: dict[str, list[FieldPoint]] = {}
    for sample, pt in rows:
        grouped.setdefault(str(sample or ""), []).append(pt)
    out: dict[str, InfiniteFieldResult | str] = {}
    for sample, pts in grouped.items():
        try:
            out[sample] = infinite_field_diso(pts, spin=spin, eta=eta)
        except Exception as exc:                              # noqa: BLE001
            out[sample] = str(exc)
    return out


def _pm(value: float, digits: int) -> str:
    """'+- 1.23' or '+- --' for a NaN (nothing could be propagated)."""
    return f"+- {value:.{digits}f}" if np.isfinite(value) else "+- --"


def _pm_scaled(value: float, scaled: float, digits: int, unit: str) -> str:
    """'+- 1.36 ppm (a priori) / +- 8.82 (scaled by sqrt(chi2/dof))' when the
    scatter exceeds the input errors, else the plain '+- 1.36 ppm'."""
    base = f"{_pm(value, digits)} {unit}"
    if np.isfinite(scaled) and np.isfinite(value) and scaled > value * (1 + 1e-9):
        return f"{base} (a priori) / +- {scaled:.{digits}f} (scaled by sqrt(chi2/dof))"
    return base


def fmt_result_lines(res: InfiniteFieldResult) -> list[str]:
    """The fitted-number lines of one sample, shared by the report and the
    dialogs so the two can never disagree on how a bound is worded."""
    out = [f"delta_iso = {res.delta_iso_ppm:.2f} "
           f"{_pm_scaled(res.delta_iso_err_ppm, res.delta_iso_err_scaled_ppm, 2, 'ppm')}"]
    if res.cq_is_bound or (res.cq_MHz == 0.0 and res.note):
        up = res.cq_upper_2sigma_MHz
        bound = f"<= {up:.3f} MHz (2-sigma upper bound)" if up > 0 else "not determined"
        pq_up = up * float(np.sqrt(1.0 + res.eta ** 2 / 3.0))
        out.append(f"P_Q {'<= ' + format(pq_up, '.3f') + ' MHz (2-sigma upper bound, eta-independent)' if up > 0 else 'not determined'}")
        out.append(f"C_Q {bound}   (eta = {res.eta:g} assumed)")
    else:
        eta_factor = float(np.sqrt(1.0 + res.eta ** 2 / 3.0))
        pq_scaled = (res.cq_err_scaled_MHz * eta_factor
                     if np.isfinite(res.cq_err_scaled_MHz) else float("nan"))
        out.append(f"P_Q = {res.pq_MHz:.3f} "
                   f"{_pm_scaled(res.pq_err_MHz, pq_scaled, 3, 'MHz')} (eta-independent)")
        lo, hi = sorted(res.cq_eta_range)
        out.append(f"C_Q = {res.cq_MHz:.3f} "
                   f"{_pm_scaled(res.cq_err_MHz, res.cq_err_scaled_MHz, 3, 'MHz')}   "
                   f"(eta = {res.eta:g} assumed; {lo:.3f}-{hi:.3f} over eta 0-1)")
    out.append(f"slope = {res.slope:.6g} {_pm(res.slope_err, 0)} ppm.MHz^2"
               if np.isfinite(res.slope_err) else
               f"slope = {res.slope:.6g} +- -- ppm.MHz^2")
    if res.dof > 0:
        p_txt = f" (p = {res.p_value:.2g})" if np.isfinite(res.p_value) else ""
        out.append(f"chi2/dof = {res.chi2:.3g}/{res.dof}{p_txt}"
                   + ("   ! improbable scatter" if res.misfit else "")
                   + ("   (no input errors: +- from the scatter)"
                      if res.errors_from_scatter else ""))
    else:
        out.append("chi2/dof = exact (2 points, no redundancy)")
    return out


def mode_label(magnitude) -> str:
    """'magnitude' / 'absorption' / '?' for a point's processing mode."""
    return "?" if magnitude is None else ("magnitude" if magnitude else "absorption")


def mixed_modes(points) -> bool:
    """True when the points whose mode is KNOWN mix magnitude and absorption
    -- two different observables that must not share one line."""
    known = {bool(p.magnitude) for p in points if p.magnitude is not None}
    return len(known) > 1


def report_text(results: dict, spin: float, eta: float, nucleus: str = "",
                widths: dict | None = None) -> str:
    """A plain-text report of an infinite-field extrapolation, for the lab
    book or the SI: every input point, every fitted number with its
    uncertainty, and the assumptions that were made."""
    lines = [
        "QCPMG infinite-field extrapolation of the isotropic chemical shift",
        "=" * 66,
        f"nucleus            : {nucleus or '(unspecified)'}",
        f"spin I             : {spin:g}",
        f"eta (ASSUMED)      : {eta:g}   <- not determined by this method",
        "model              : dcg = diso + slope/nu0^2   (Sandland 2004 Eq. 1)",
        "",
    ]
    ok = [(s, r) for s, r in results.items()
          if isinstance(r, InfiniteFieldResult)]
    bad = [(s, r) for s, r in results.items() if not isinstance(r, InfiniteFieldResult)]

    for sample, res in ok:
        lines.append(f"--- {sample or '(unnamed)'} " + "-" * max(0, 60 - len(sample)))
        lines.append("    nu0 (MHz)      dcg (ppm)   +- err    resid   mode        "
                     "CT-selective (declared)")
        for i, p in enumerate(res.points):
            sel = "?" if p.ct_selective is None else ("yes" if p.ct_selective else "no")
            err = f"{p.dcg_err_ppm:7.2f}" if p.has_err else f"{'n/a':>7s}"
            resid = (f"{res.residuals_ppm[i]:7.2f}" if i < len(res.residuals_ppm)
                     else f"{'':7s}")
            lines.append(f"    {p.larmor_MHz:10.4f}  {p.dcg_ppm:11.2f}  "
                         f"{err}  {resid}  {mode_label(p.magnitude):10s}  {sel}")
            prov = p.provenance()
            if prov:
                lines.append(f"                 {prov}")
            for fl in p.flags:
                lines.append(f"                 {fl}")
            if p.cg_sequence:
                lines.append(f"                 {p.cg_sequence}")
        for ln in referencing_checks(res.points):
            lines.append(f"    {ln}")
        if mixed_modes(res.points):
            lines.append("    ! NOT COMPARABLE: mixed magnitude/absorption dcg -- "
                         "reprocess both fields the same way")
        elif res.points and all(p.magnitude for p in res.points):
            lines.append("    all points measured on magnitude (mc) spectra -- dcg is the "
                         "|spectrum| centroid, not the absorption one;")
            lines.append("    rectified noise pulls it toward the window centre, growing "
                         "with window width and 1/(S/N); compare with the")
            lines.append("    phased spectrum on the same window before quoting P_Q to "
                         "better than a few %")
        for ln in fmt_result_lines(res):
            key, _, rest = ln.partition(" ")
            lines.append(f"    {key:14s} {rest}")
        if res.note:
            for part in res.note.split("; "):
                lines.append(f"    note           : {part}")
        if res.warning:
            for part in res.warning.split("; "):
                lines.append(f"    ! {part}")
        w = (widths or {}).get(sample)
        if w is not None and getattr(w, "ok", False):
            fields = ", ".join(f"{f:.1f}" for f in w.fields_MHz)
            what = ("FWHM, Sandland Eq. 2" if w.kind == "fwhm"
                    else "second moment, Gaussian-equivalent FWHM = 2.3548 sigma")
            lines.append(f"    width split over {w.n_fields} fields ({fields} MHz; {what}"
                         + ("; 1/sigma^2 weights" if w.weighted else "") + ")")
            lines.append(f"    W_q            = {w.wq_lo_ppm:8.1f} ppm (low field)"
                         f" / {w.wq_hi_ppm:.1f} ppm (high field)"
                         + (f"   +- {w.wq_err_ppm:.1f} (low)" if np.isfinite(w.wq_err_ppm) else ""))
            lines.append(f"    W_csd          = {w.wcsd_ppm:8.1f} ppm "
                         f"(field-independent: shift distribution + CSA)"
                         + (f"   +- {w.wcsd_err_ppm:.1f}" if np.isfinite(w.wcsd_err_ppm) else ""))
            if w.n_fields > 2 and np.isfinite(w.chi2):
                lines.append(f"    width chi2     = {w.chi2:.3g} over {w.n_fields - 2} dof")
            if w.gate:
                lines.append(f"    ! {w.gate}")
        elif w is not None and getattr(w, "note", ""):
            lines.append(f"    width split    : {w.note}")
        lines.append("")

    if ok:
        lines.append("Summary")
        lines.append("-" * 66)
        lines.append("sample                         diso (ppm)        "
                     "P_Q (MHz)          C_Q (MHz)")
        for sample, res in ok:
            if res.cq_is_bound or (res.cq_MHz == 0.0 and res.note):
                up = res.cq_upper_2sigma_MHz
                pq = (f"<= {up * float(np.sqrt(1.0 + res.eta ** 2 / 3.0)):.3f} (2-sigma)"
                      if up > 0 else "not determined")
                cq = f"<= {up:.3f} (2-sigma)" if up > 0 else "not determined"
            else:
                pq = f"{res.pq_MHz:6.3f} {_pm(res.pq_err_MHz, 3)}"
                cq = f"{res.cq_MHz:6.3f} {_pm(res.cq_err_MHz, 3)}"
            lines.append(f"{(sample or '(unnamed)')[:28]:28s}  "
                         f"{res.delta_iso_ppm:7.2f} {_pm(res.delta_iso_err_ppm, 2):<9s} "
                         f"{pq:18s} {cq}")
        lines.append("")
    for sample, why in bad:
        lines.append(f"NOT FITTED  {sample or '(unnamed)'}: {why}")
    if bad:
        lines.append("")
    lines.append("CT-selective is the operator's declaration, recorded for")
    lines.append("provenance; it does not enter the fit.")
    lines.append("delta_iso is the intercept at 1/nu0^2 -> 0. P_Q = C_Q sqrt(1+eta^2/3)")
    lines.append("follows from the slope WITHOUT eta (C_Q^2 (3+eta^2) = 3 P_Q^2), so quote")
    lines.append("P_Q +- when eta is unknown; C_Q needs the assumed eta and its range over")
    lines.append("eta = 0-1 is a systematic kept out of the +-. A slope that is not")
    lines.append("significantly negative gives only an upper bound on C_Q; a")
    lines.append("significantly positive slope means the dcg values are not the")
    lines.append("same observable (mode, window or referencing differ).")
    return "\n".join(lines)
