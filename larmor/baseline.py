"""Iterative baseline correction for 1D spectra (Yon et al. 2020).

An adaptation of the core algorithm of the *MY Baseline Corrector*:

    M. Yon, F. Fayon, D. Massiot, V. Sarou-Kanian, "Iterative baseline
    correction algorithm for dead time truncated one-dimensional solid-state
    MAS NMR spectra", Solid State Nucl. Magn. Reson. 110, 101699 (2020).
    doi:10.1016/j.ssnmr.2020.101699 ;  https://github.com/maximeYon/Baseline_Corrector

It removes the rolling baseline produced by the truncation of the first FID
points by the receiver **dead time** in pulse-acquire (single-pulse) MAS NMR --
the case where a polynomial baseline fails. Two ideas do the work:

1. **Histogram filter** -- the baseline/noise points are selected automatically
   at each iteration by thresholding on the intensity histogram: the noise forms
   the dominant peak of the histogram (its mode), and everything up to the width
   of that peak above the mode is treated as baseline; taller points are peaks
   and are excluded. No manual baseline regions.
2. **Dead-time (time-domain) restriction** -- the estimated baseline is
   transformed to the time domain and truncated to the first points (set by the
   dead time), because a dead-time truncation artefact only produces *broad*,
   slowly varying baseline structure. This keeps the correction from eating real
   (narrow) peaks.

The baseline is built up **iteratively**: fit -> restrict -> subtract -> repeat
until the correction is negligible.

This is a from-scratch NumPy/SciPy port of the published method (not the MATLAB
GUI); cite the paper above when you use it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BaselineResult:
    baseline: np.ndarray      # the estimated baseline (same length as y)
    corrected: np.ndarray     # y - baseline
    n_iter: int
    converged: bool


def _histogram_band(resid: np.ndarray, factor: float = 1.0) -> tuple[float, float]:
    """Intensity band of the noise/baseline (Yon et al., histogram filter).

    The noise is the tallest peak of the intensity histogram; the band spans one
    noise-peak-width either side of its mode. Points inside the band are
    baseline; taller points are (positive) peaks and lower points are negative
    excursions/artefacts -- both are excluded (the published SymThresh, used here
    by default because a single upper cut would sweep negative spikes into the
    baseline and blow up the spline)."""
    n = resid.size
    bins = max(int(n / 10), 16)
    hist, edges = np.histogram(resid, bins=bins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    # light smoothing of the histogram (MATLAB smooth(hi, 8))
    k = np.ones(8) / 8.0
    hist = np.convolve(hist.astype(float), k, mode="same")
    mode_idx = int(np.argmax(hist))
    max_h = hist[mode_idx]
    noise_width = int(np.count_nonzero(hist >= max_h / 2.0))   # bins at >= half max
    w = int(round(noise_width * factor))
    hi = float(centres[min(mode_idx + w, bins - 1)])
    lo = float(centres[max(mode_idx - w, 0)])
    return lo, hi


def _smooth_spline(x: np.ndarray, y: np.ndarray, mask: np.ndarray,
                   lam: float) -> np.ndarray:
    """Penalised smoothing spline through the selected (baseline) points,
    evaluated on the full axis. `lam` is the smoothing strength (larger =
    stiffer/smoother); x is normalised to [0, 1] so `lam` is length-independent."""
    from scipy.interpolate import make_smoothing_spline

    xs = x[mask]
    ys = y[mask]
    if xs.size < 6:                       # too few points -> flat (median) baseline
        return np.full_like(x, float(np.median(ys)) if ys.size else 0.0)
    xn = (x - x[0]) / (x[-1] - x[0])
    xsn = xn[mask]
    # make_smoothing_spline needs strictly increasing, unique x
    xsn, uidx = np.unique(xsn, return_index=True)
    ys = ys[uidx]
    # normalise y so `lam` (the 2nd-derivative penalty) is intensity-independent:
    # the spline penalty scales with y^2, so without this a 1e6-intensity spectrum
    # and a unit one would need very different lam.
    yscale = float(np.std(ys)) or 1.0
    try:
        spl = make_smoothing_spline(xsn, ys / yscale, lam=lam)
        return np.asarray(spl(xn), float) * yscale
    except Exception:
        # fall back to a robust low-order polynomial if the spline solver fails
        coef = np.polyfit(xsn, ys / yscale, 3)
        return np.polyval(coef, xn) * yscale


def _fid_restrict(baseline: np.ndarray, keep_pts: int) -> np.ndarray:
    """Restrict a baseline estimate to broad components by keeping only the
    first `keep_pts` time-domain points (the dead-time window), per Yon et al.

    Mirrors the MATLAB: ifftshift(ifft(fftshift(.))) -> keep the central window
    around t=0 -> fftshift(fft(fftshift(.))) -> real."""
    n = baseline.size
    fid = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(baseline)))
    half = max(int(keep_pts) // 2, 1)
    c = n // 2
    keep = np.zeros(n, dtype=complex)
    lo, hi = max(c - half, 0), min(c + half, n)
    keep[lo:hi] = fid[lo:hi]
    out = np.fft.fftshift(np.fft.fft(np.fft.fftshift(keep)))
    return np.asarray(out.real, float)


def iterative_baseline(y: np.ndarray, *, x: np.ndarray | None = None,
                       dead_time_pts: int = 0, smoothness: float = 1.0,
                       threshold_factor: float = 1.0, max_iter: int = 60,
                       tol: float = 1e-6) -> BaselineResult:
    """Estimate and remove a rolling baseline (Yon et al. 2020).

    Parameters
    ----------
    y : the real spectrum (1-D).
    x : optional abscissa (only its spacing matters); defaults to point index.
    dead_time_pts : if > 0, restrict the baseline to this many time-domain points
        (~ 2 * dead_time / dwell). 0 disables the dead-time restriction (then the
        method is a plain iterative histogram-thresholded smoothing baseline).
    smoothness : smoothing-spline strength (larger = stiffer/smoother baseline;
        the default is stiff, appropriate for a slowly rolling baseline).
    threshold_factor : scales the histogram threshold (dmfit's AdvanceFilter/100).
    max_iter : hard cap on iterations.
    tol : iterate until the per-iteration correction, normalised by the peak,
        drops below `tol` OR the incremental correction falls below the noise
        floor (nothing baseline-like left to remove).
    """
    y = np.asarray(y, float)
    n = y.size
    idx = np.arange(n, dtype=float) if x is None else np.asarray(x, float)
    peak = float(np.max(np.abs(y))) or 1.0

    resid = y.copy()
    baseline_total = np.zeros(n)
    converged = False
    noise = None
    it = 0
    for it in range(1, max_iter + 1):
        lo, hi = _histogram_band(resid, factor=threshold_factor)
        mask = (resid > lo) & (resid < hi)
        if np.count_nonzero(mask) < 6:
            break
        if noise is None:                       # robust noise from the 1st pass
            sel = resid[mask]
            noise = float(np.median(np.abs(sel - np.median(sel)))) * 1.4826 or 1.0
        est = _smooth_spline(idx, resid, mask, lam=smoothness)
        if dead_time_pts and dead_time_pts > 1:
            est = _fid_restrict(est, dead_time_pts)
        resid = resid - est
        baseline_total += est
        crit = float(np.sum(np.abs(est))) / (n * peak)
        # stop once the correction is negligible vs the peak OR below the noise
        if crit < tol or float(np.sqrt(np.mean(est ** 2))) < 0.1 * noise:
            converged = True
            break
    return BaselineResult(baseline=baseline_total, corrected=y - baseline_total,
                          n_iter=it, converged=converged)
