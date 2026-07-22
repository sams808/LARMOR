"""Region measurements: integral, FWHM, centre of mass, S/N — the numbers you
read off a spectrum without fitting (dmfit/TopSpin integration, ssNake FWHM /
Centre of Mass / Integrals). Pure and tested."""
from __future__ import annotations

import numpy as np

_trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def _region_mask(x: np.ndarray, region) -> np.ndarray:
    if region is None:
        return np.ones(np.shape(x), bool)
    lo, hi = min(region), max(region)
    return (x >= lo) & (x <= hi)


def integrate(x: np.ndarray, y: np.ndarray, region=None) -> float:
    """∫ y dx over a ppm region (trapezoidal; sign follows the x order)."""
    m = _region_mask(x, region)
    if m.sum() < 2:
        return 0.0
    return float(abs(_trap(np.asarray(y)[m], np.asarray(x)[m])))


def centre_of_mass(x: np.ndarray, y: np.ndarray, region=None) -> float:
    """Intensity-weighted centroid (ppm) over a region."""
    m = _region_mask(x, region)
    xs, ys = np.asarray(x)[m], np.asarray(y)[m]
    w = ys - min(ys.min(), 0.0)
    s = w.sum()
    return float((xs * w).sum() / s) if s else float(xs.mean() if xs.size else 0.0)


def fwhm(x: np.ndarray, y: np.ndarray, region=None) -> float:
    """Full width at half maximum (ppm) of the tallest peak in a region, by
    linear interpolation of the half-max crossings."""
    m = _region_mask(x, region)
    xs, ys = np.asarray(x, float)[m], np.asarray(y, float)[m]
    if xs.size < 3:
        return 0.0
    o = np.argsort(xs); xs, ys = xs[o], ys[o]
    base = ys.min()
    i0 = int(np.argmax(ys))
    half = base + 0.5 * (ys[i0] - base)

    def cross(seq_i, direction):
        i = i0
        while 0 <= i + direction < ys.size and ys[i] > half:
            i += direction
        if i == i0 or not (0 <= i < ys.size):
            return xs[i0]
        j = i - direction                       # last point above half
        if ys[j] == ys[i]:
            return xs[i]
        t = (half - ys[j]) / (ys[i] - ys[j])
        return xs[j] + t * (xs[i] - xs[j])

    left = cross(i0, -1)
    right = cross(i0, +1)
    return float(abs(right - left))


def snr(x: np.ndarray, y: np.ndarray, signal_region=None, noise_region=None
        ) -> float:
    """Peak signal in signal_region ÷ RMS of noise_region (edges if unset)."""
    y = np.asarray(y, float)
    ms = _region_mask(x, signal_region)
    if noise_region is None:
        edge = max(5, y.size // 10)
        noise = np.concatenate([y[:edge], y[-edge:]])
    else:
        noise = y[_region_mask(x, noise_region)]
    rms = float(np.std(noise - np.median(noise))) or 1e-12
    sig = float(np.max(np.abs(y[ms] - np.median(noise)))) if ms.any() else 0.0
    return sig / rms


def integrate_regions(x, y, regions) -> list[dict]:
    """Per-region {range, integral, percent, centre, fwhm}. Percents sum to 100
    over the supplied regions."""
    rows = [{"range": (max(r), min(r)),
             "integral": integrate(x, y, r),
             "centre": centre_of_mass(x, y, r),
             "fwhm": fwhm(x, y, r)} for r in regions]
    total = sum(row["integral"] for row in rows) or 1.0
    for row in rows:
        row["percent"] = 100.0 * row["integral"] / total
    return rows
