"""Multi-experiment correlation decomposition (roadmap Priority 9).

ARCHITECTURE ONLY — this engine is intentionally NOT wired into the app yet
(see ROADMAP.md). It generalizes the HMQC "1D − projection" difference to any
set of shared-nucleus observables (1D spectra, 2D projections, REDOR ΔS/S…):

  * every dataset contributes an *observable* on a shared ppm axis;
  * `align` puts them on a common grid;
  * `intersection` returns the feature present in ALL of them (correlated);
  * `difference` returns what is specific to one target after removing scaled
    contributions of the others (un-correlated).

The two-dataset case (a 1D and an HMQC projection) reproduces the existing
HMQC un-correlated-feature extraction exactly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Observable:
    """A 1D trace on a shared axis, from any experiment."""

    label: str
    ppm: np.ndarray
    amp: np.ndarray
    kind: str = "1d"                # "1d" | "projection" | "redor" | ...

    def on(self, grid: np.ndarray) -> np.ndarray:
        o = np.argsort(self.ppm)
        return np.interp(grid, np.asarray(self.ppm)[o], np.asarray(self.amp)[o],
                         left=0.0, right=0.0)


def common_grid(observables: list[Observable], n: int | None = None,
                mode: str = "union") -> np.ndarray:
    """A shared ascending ppm grid spanning the observables (union or
    intersection of their ranges)."""
    los = [float(np.min(o.ppm)) for o in observables]
    his = [float(np.max(o.ppm)) for o in observables]
    if mode == "intersection":
        lo, hi = max(los), min(his)
    else:
        lo, hi = min(los), max(his)
    if n is None:
        n = max(len(o.ppm) for o in observables)
    return np.linspace(lo, hi, int(n))


def scale_to(target: np.ndarray, source: np.ndarray, region_mask=None) -> float:
    """Least-squares non-negative scale s minimising ||target − s·source||²
    over an optional mask (default: everywhere)."""
    t, s = np.asarray(target, float), np.asarray(source, float)
    if region_mask is not None:
        t, s = t[region_mask], s[region_mask]
    denom = float(s @ s)
    return max(0.0, float(t @ s) / denom) if denom > 0 else 0.0


def align(observables: list[Observable], grid: np.ndarray | None = None
          ) -> tuple[np.ndarray, list[np.ndarray]]:
    grid = common_grid(observables) if grid is None else grid
    return grid, [o.on(grid) for o in observables]


def intersection(observables: list[Observable], grid: np.ndarray | None = None,
                 ) -> tuple[np.ndarray, np.ndarray]:
    """The feature present in ALL observables (correlated): each is peak-scaled
    to unit height, then the pointwise minimum is taken. Returns (grid, amp)."""
    grid, amps = align(observables, grid)
    norm = [a / (np.abs(a).max() or 1.0) for a in amps]
    return grid, np.minimum.reduce(norm) if norm else grid * 0.0


def difference(target: Observable, references: list[Observable],
               grid: np.ndarray | None = None, region=None,
               ) -> tuple[np.ndarray, np.ndarray]:
    """What is specific to ``target`` after removing least-squares-scaled
    references (un-correlated). Returns (grid, amp)."""
    grid = common_grid([target, *references]) if grid is None else grid
    t = target.on(grid)
    mask = None
    if region is not None:
        lo, hi = min(region), max(region)
        mask = (grid >= lo) & (grid <= hi)
    out = t.copy()
    for ref in references:
        r = ref.on(grid)
        out = out - scale_to(t, r, mask) * r
    return grid, out
