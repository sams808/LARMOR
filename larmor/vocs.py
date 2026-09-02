"""Frequency-stepped (VOCS) acquisition: stitch sub-spectra into one pattern.

A pattern wider than the probe/pulse bandwidth is acquired as a series of
EXPNOs at stepped transmitter offsets (Variable Offset Cumulative Spectroscopy;
Massiot et al. 1995). Each processed sub-spectrum carries a correctly
referenced ppm axis (the Bruker SF/SR referencing is absolute, so different
O1 settings land on one shared ppm scale), which makes stitching a pure
resampling-and-combining problem:

- **skyline**: max across pieces at each point -- the standard VOCS
  projection; immune to how many pieces cover a point, sensitive to a noisy
  piece's baseline.
- **average**: mean over the pieces that COVER each point (coverage-weighted,
  so overlap regions are not double-counted); needs the pieces to be on a
  consistent intensity scale (same NS/RG or normalised).

Edges: excitation rolls off toward each piece's edges, so ``trim_frac`` drops
that fraction of every piece's span from BOTH its ends before combining
(default 10 %). Set 0 to keep everything.

Qt-free; the desktop wraps this in Tools > Stitch frequency-stepped (VOCS).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class StitchResult:
    ppm: np.ndarray
    amp: np.ndarray
    #: per output point, how many pieces contributed (after trimming)
    coverage: np.ndarray
    #: human-readable record for recipe notes / provenance
    notes: list = field(default_factory=list)


def stitch(pieces: list[tuple[np.ndarray, np.ndarray]], *,
           method: str = "skyline", trim_frac: float = 0.10,
           npts: int | None = None,
           labels: list[str] | None = None) -> StitchResult:
    """Combine ``pieces`` = [(ppm, amp), ...] into one spectrum.

    The output grid spans the union of the (trimmed) pieces at the FINEST
    spacing any piece provides (or ``npts`` when given), so no piece is
    downsampled below its own resolution.
    """
    if not pieces:
        raise ValueError("no spectra to stitch")
    if method not in ("skyline", "average"):
        raise ValueError(f"unknown stitch method {method!r}")
    trim_frac = float(np.clip(trim_frac, 0.0, 0.45))

    prepared: list[tuple[np.ndarray, np.ndarray]] = []
    finest = np.inf
    lo_all, hi_all = np.inf, -np.inf
    for ppm, amp in pieces:
        x = np.asarray(ppm, float).ravel()
        y = np.asarray(amp, float).ravel()
        if x.size < 2 or y.size != x.size:
            raise ValueError("each piece needs matching 1D ppm/amp arrays")
        order = np.argsort(x)
        x, y = x[order], y[order]
        span = x[-1] - x[0]
        cut = trim_frac * span
        keep = (x >= x[0] + cut) & (x <= x[-1] - cut)
        x, y = x[keep], y[keep]
        if x.size < 2:
            raise ValueError("a piece vanished entirely under trim_frac "
                             f"{trim_frac:g} -- lower it")
        finest = min(finest, float(np.median(np.diff(x))))
        lo_all, hi_all = min(lo_all, x[0]), max(hi_all, x[-1])
        prepared.append((x, y))

    n = int(npts) if npts else int(round((hi_all - lo_all) / finest)) + 1
    n = int(np.clip(n, 2, 262144))
    grid = np.linspace(lo_all, hi_all, n)

    stack = np.full((len(prepared), n), np.nan)
    for i, (x, y) in enumerate(prepared):
        inside = (grid >= x[0]) & (grid <= x[-1])
        stack[i, inside] = np.interp(grid[inside], x, y)

    coverage = np.sum(~np.isnan(stack), axis=0)
    covered = coverage > 0
    out = np.zeros(n)
    with np.errstate(invalid="ignore"):
        if method == "skyline":
            out[covered] = np.nanmax(stack[:, covered], axis=0)
        else:
            out[covered] = np.nanmean(stack[:, covered], axis=0)

    names = labels or [f"piece {i + 1}" for i in range(len(prepared))]
    # covered runs: rising edges (counting a covered grid START as one)
    rising = int(np.sum(np.diff(covered.astype(int)) == 1)) + int(covered[0])
    gaps = rising - 1
    notes = [f"VOCS stitch ({method}) of {len(prepared)} sub-spectra, "
             f"trim {trim_frac:.0%} per edge: " + ", ".join(names)]
    if gaps > 0:
        notes.append(f"WARNING: {gaps} uncovered gap(s) between pieces -- "
                     "the offsets do not tile the pattern; those points are "
                     "zero, not measured")
    return StitchResult(ppm=grid, amp=out, coverage=coverage, notes=notes)
