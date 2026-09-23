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

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


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
                         dcg_err_ppm: float | None = None) -> "FieldPoint":
        """A point from a :class:`larmor.qcpmg.CgMeasurement`: the error is
        max(jitter sigma, convergence drift) floored at ERR_FLOOR_PPM, and
        window, mode, flags and the convergence sequence travel with it.
        ``dcg_ppm`` / ``dcg_err_ppm`` override the measured values when the
        user edited the cells."""
        err = meas.sigma_ppm if dcg_err_ppm is None else dcg_err_ppm
        err = max(float(err), ERR_FLOOR_PPM) if np.isfinite(err) else 0.0
        conv = getattr(meas, "convergence", None)
        return cls(float(larmor_MHz),
                   float(meas.cg_ppm if dcg_ppm is None else dcg_ppm), err,
                   ct_selective, label, magnitude=magnitude,
                   window=tuple(meas.window), window_mode=meas.mode,
                   source=source, rotor_Hz=float(rotor_Hz or 0.0),
                   flags=tuple(meas.flags),
                   cg_sequence=conv.sequence() if conv is not None else "",
                   drift_ppm=float(meas.drift_ppm))


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

    def line(self, inv_nu2: np.ndarray) -> np.ndarray:
        """δcg on the fit line for given 1/ν0² values (for plotting)."""
        return self.intercept + self.slope * np.asarray(inv_nu2, float)

    @property
    def cq_is_bound(self) -> bool:
        """True when C_Q is reported as an upper bound, not a value."""
        return self.cq_upper_2sigma_MHz > 0.0 and self.cq_MHz == 0.0


def _quad_A(spin: float, eta: float) -> float:
    """A in slope = −A·C_Q² (ppm·MHz² per MHz²)."""
    return (1.0e6 / 40.0) * (3.0 + eta ** 2) * _spin_factor(spin)


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
    if weighted:
        sig_a, sig_b = float(np.sqrt(var_a)), float(np.sqrt(var_b))
    elif dof > 0:
        s2 = chi2 / dof                       # OLS residual variance
        sig_a, sig_b = float(np.sqrt(var_a * s2)), float(np.sqrt(var_b * s2))
        notes.append("no input errors given: +- from scatter about the line")
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
    pq = cq * float(np.sqrt(1.0 + eta ** 2 / 3.0))

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

    return InfiniteFieldResult(
        delta_iso_ppm=float(a), delta_iso_err_ppm=sig_a,
        cq_MHz=cq, cq_err_MHz=cq_err, pq_MHz=pq, eta=eta, spin=spin,
        slope=float(b), intercept=float(a), points=list(points),
        slope_err=sig_b, cq_upper_2sigma_MHz=cq_upper,
        note="; ".join(notes), warning="; ".join(warnings),
        lever_arm=xr, weighted=weighted)


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


def fmt_result_lines(res: InfiniteFieldResult) -> list[str]:
    """The fitted-number lines of one sample, shared by the report and the
    dialogs so the two can never disagree on how a bound is worded."""
    out = [f"delta_iso = {res.delta_iso_ppm:.2f} {_pm(res.delta_iso_err_ppm, 2)} ppm"]
    if res.cq_is_bound or (res.cq_MHz == 0.0 and res.note):
        up = res.cq_upper_2sigma_MHz
        bound = f"<= {up:.3f} MHz (2-sigma upper bound)" if up > 0 else "not determined"
        out.append(f"C_Q {bound}   (eta = {res.eta:g} assumed)")
        pq_up = up * float(np.sqrt(1.0 + res.eta ** 2 / 3.0))
        out.append(f"P_Q {'<= ' + format(pq_up, '.3f') + ' MHz (2-sigma upper bound)' if up > 0 else 'not determined'}")
    else:
        out.append(f"C_Q = {res.cq_MHz:.3f} {_pm(res.cq_err_MHz, 3)} MHz   "
                   f"(eta = {res.eta:g} assumed)")
        out.append(f"P_Q = {res.pq_MHz:.3f} MHz")
    out.append(f"slope = {res.slope:.6g} {_pm(res.slope_err, 0)} ppm.MHz^2"
               if np.isfinite(res.slope_err) else
               f"slope = {res.slope:.6g} +- -- ppm.MHz^2")
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
        lines.append("    nu0 (MHz)      dcg (ppm)   +- err   mode        "
                     "CT-selective (declared)")
        for p in res.points:
            sel = "?" if p.ct_selective is None else ("yes" if p.ct_selective else "no")
            err = f"{p.dcg_err_ppm:7.2f}" if p.has_err else f"{'n/a':>7s}"
            lines.append(f"    {p.larmor_MHz:10.4f}  {p.dcg_ppm:11.2f}  "
                         f"{err}   {mode_label(p.magnitude):10s}  {sel}")
            for fl in p.flags:
                lines.append(f"                 {fl}")
            if p.cg_sequence:
                lines.append(f"                 {p.cg_sequence}")
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
            lines.append(f"    W_q            = {w.wq_lo_ppm:8.1f} ppm (low field)"
                         f" / {w.wq_hi_ppm:.1f} ppm (high field)")
            lines.append(f"    W_csd          = {w.wcsd_ppm:8.1f} ppm "
                         f"(field-independent)")
        elif w is not None and getattr(w, "note", ""):
            lines.append(f"    width split    : {w.note}")
        lines.append("")

    if ok:
        lines.append("Summary")
        lines.append("-" * 66)
        lines.append("sample                         diso (ppm)      C_Q (MHz)")
        for sample, res in ok:
            if res.cq_is_bound or (res.cq_MHz == 0.0 and res.note):
                cq = (f"<= {res.cq_upper_2sigma_MHz:.3f} (2-sigma)"
                      if res.cq_upper_2sigma_MHz > 0 else "not determined")
            else:
                cq = f"{res.cq_MHz:6.3f} {_pm(res.cq_err_MHz, 3)}"
            lines.append(f"{(sample or '(unnamed)')[:28]:28s}  "
                         f"{res.delta_iso_ppm:7.2f} {_pm(res.delta_iso_err_ppm, 2):<9s} "
                         f"{cq}")
        lines.append("")
    for sample, why in bad:
        lines.append(f"NOT FITTED  {sample or '(unnamed)'}: {why}")
    if bad:
        lines.append("")
    lines.append("CT-selective is the operator's declaration, recorded for")
    lines.append("provenance; it does not enter the fit.")
    lines.append("delta_iso is the intercept at 1/nu0^2 -> 0; C_Q follows from the")
    lines.append("slope with the assumed eta, so its accuracy is limited by that")
    lines.append("assumption. Quote P_Q when eta is unknown. A slope that is not")
    lines.append("significantly negative gives only an upper bound on C_Q; a")
    lines.append("significantly positive slope means the dcg values are not the")
    lines.append("same observable (mode, window or referencing differ).")
    return "\n".join(lines)
