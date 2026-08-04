"""Batch fit: one shared model fitted to many 1D spectra at once.

The everyday case is a whole sample series measured the same way: you want ONE
set of lineshape/position parameters that describes every spectrum, with only
the **amplitudes** free per spectrum (relative populations change, the sites do
not). That is a co-fit with *everything except amplitude* shared — and it is far
more robust than fitting each spectrum alone, because the shared parameters are
constrained by all the data together.

An optional second stage **releases** chosen parameters: they come off the
shared tie and are allowed to drift by a small ±fraction around the shared value,
independently per spectrum (a "relaxation" — e.g. let δ_iso wander ±5 % while
widths stay locked).

Built on larmor.multifit.fit_cofit; 1D only. Qt-free and testable.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from larmor.multifit import fit_cofit
from larmor.recipe import Recipe


def all_but_amplitude(recipes: list[Recipe]) -> tuple[str, ...]:
    """Every parameter name present in the model except amplitude (the default
    shared set: identical shapes/positions, free amplitudes)."""
    names: set[str] = set()
    for r in recipes:
        for s in r.sites:
            names.update(s.params.keys())
    names.discard("amplitude")
    return tuple(sorted(names))


@dataclass
class BatchFitResult:
    recipes: list[Recipe]                 # per-spectrum fitted recipes (final)
    labels: list[str]
    rmsd: list[float]
    per_dataset: list[dict]               # {"x", "y_fit"} per spectrum
    shared: tuple[str, ...]               # parameters still shared at the end
    released: tuple[str, ...]             # parameters loosened per spectrum
    release_frac: float = 0.0
    lmfit_result: object = None
    warnings: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        rel = (f" · released {', '.join(self.released)} (±{self.release_frac:.0%})"
               if self.released else "")
        return (f"batch fit: {len(self.recipes)} spectra, "
                f"{len(self.shared)} shared parameters{rel} · "
                f"mean RMSD {np.mean(self.rmsd):.4f}")


def _relax_bounds(value: float, frac: float, cur_min, cur_max):
    """A small ±frac window around `value`, clipped to any physical bounds."""
    d = abs(value) * frac
    if d < 1e-12:
        d = frac                       # value ≈ 0 → treat frac as an absolute floor
    lo, hi = value - d, value + d
    if cur_min is not None:
        lo = max(lo, cur_min)
    if cur_max is not None:
        hi = min(hi, cur_max)
    return lo, hi


def batch_fit(entries: list[tuple], *, share: tuple[str, ...] | None = None,
              release: tuple[str, ...] = (), release_frac: float = 0.1,
              iter_cb=None) -> BatchFitResult:
    """Fit one shared model to several 1D spectra (amplitudes free per spectrum).

    `entries` is a list of ``(recipe, ppm, amp, window)`` — every recipe must
    have the same sites/models (apply one model to all). `share` defaults to
    *all but amplitude*. If `release` is given, those parameters are then let
    loose by ±`release_frac` around their shared value, per spectrum.
    """
    if len(entries) < 2:
        raise ValueError("batch fit needs at least two spectra")
    recipes = [e[0] for e in entries]
    windows = [e[3] for e in entries]
    if share is None:
        share = all_but_amplitude(recipes)

    def _conv(recs):
        return [(recs[k], (np.asarray(entries[k][1], float),
                           np.asarray(entries[k][2], float)))
                for k in range(len(entries))]

    # ---- stage 1: shared shape, free amplitudes ----
    res = fit_cofit(_conv(recipes), share=share, windows=windows, iter_cb=iter_cb)
    final_share, released = share, ()

    # ---- stage 2 (optional): release chosen parameters, slightly, per spectrum
    released = tuple(p for p in release if p in share)
    if released:
        share2 = tuple(p for p in share if p not in released)
        master = res.recipes[0]
        for i, site in enumerate(master.sites):
            for p in released:
                mp = site.params.get(p)
                if mp is None or mp.expr:
                    continue
                v = float(mp.value)
                for rec in res.recipes:
                    pp = rec.sites[i].params.get(p)
                    if pp is None or pp.expr:
                        continue
                    pp.value = v
                    pp.min, pp.max = _relax_bounds(v, release_frac, pp.min, pp.max)
                    pp.vary = True
        res = fit_cofit(_conv(res.recipes), share=share2, windows=windows,
                        iter_cb=iter_cb)
        final_share = share2

    labels = [(r.sample or f"spectrum {k + 1}")
              for k, r in enumerate(res.recipes)]
    return BatchFitResult(
        recipes=res.recipes, labels=labels, rmsd=res.rmsd,
        per_dataset=res.per_dataset, shared=final_share, released=released,
        release_frac=release_frac if released else 0.0,
        lmfit_result=res.lmfit_result)


def shared_table(result: BatchFitResult) -> list[dict]:
    """The fitted values for a report: shared parameters (once, from the master)
    and the per-spectrum amplitudes + any released parameters."""
    rows = []
    master = result.recipes[0]
    for i, site in enumerate(master.sites):
        for pn, p in site.params.items():
            if pn in result.shared:
                rows.append({"scope": "shared", "site": f"s{i}",
                             "label": site.label or site.model, "param": pn,
                             "value": p.value, "stderr": p.stderr})
    for k, rec in enumerate(result.recipes):
        for i, site in enumerate(rec.sites):
            for pn, p in site.params.items():
                if pn == "amplitude" or pn in result.released:
                    rows.append({"scope": result.labels[k], "site": f"s{i}",
                                 "label": site.label or site.model, "param": pn,
                                 "value": p.value, "stderr": p.stderr})
    return rows
