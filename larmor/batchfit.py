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


def free_amplitudes(recipes: list[Recipe]) -> None:
    """Make every site's amplitude free and allow it to reach zero, in place.

    A batch fit refines amplitudes per spectrum, so a recipe that *locked* an
    amplitude (``vary=False``) or bounded it above zero (``min>0``) must be
    overridden — otherwise a component can neither adapt per spectrum nor vanish
    where it is absent. Linked amplitudes (``expr``) are left alone (they follow
    their master)."""
    for r in recipes:
        for s in r.sites:
            amp = s.params.get("amplitude")
            if amp is None or getattr(amp, "expr", None):
                continue
            amp.vary = True
            if amp.min is None or amp.min > 0:
                amp.min = 0.0


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
              iter_cb=None, tol=None) -> BatchFitResult:
    """Apply one recipe to several 1D spectra, fitting **only the amplitudes**
    (per spectrum), unless parameters are explicitly released.

    ``entries`` is a list of ``(recipe, ppm, amp, window)`` — every recipe carries
    the same model. The loaded recipe is treated as the answer for lineshape: every
    parameter is **held fixed at its recipe value except the amplitude**, which is
    always free per spectrum (and may reach zero). Parameters named in ``release``
    are additionally freed and fit **per spectrum**, allowed to drift by
    ±``release_frac`` around their recipe value. (``share`` is accepted for
    backward compatibility and ignored.)
    """
    if len(entries) < 2:
        raise ValueError("batch fit needs at least two spectra")
    recipes = [e[0] for e in entries]
    windows = [e[3] for e in entries]
    released = tuple(release)

    # amplitudes: always free per spectrum, allowed to reach zero (a line can be
    # absent in some spectra) — overriding any recipe lock/positive lower bound
    free_amplitudes(recipes)
    # every OTHER parameter is FIXED at the recipe value, except those the user
    # released (which are fit per spectrum within ±release_frac of the recipe value)
    for r in recipes:
        for s in r.sites:
            for pn, p in s.params.items():
                if pn == "amplitude" or getattr(p, "expr", None):
                    continue                     # amplitude free; links follow master
                if pn in released:
                    p.vary = True
                    p.min, p.max = _relax_bounds(float(p.value), release_frac,
                                                 p.min, p.max)
                else:
                    p.vary = False               # held at the recipe value

    def _conv(recs):
        return [(recs[k], (np.asarray(entries[k][1], float),
                           np.asarray(entries[k][2], float)))
                for k in range(len(entries))]

    # nothing is tied across spectra: amplitudes (and any released params) are
    # optimised independently per spectrum; the fixed shape is the recipe's
    res = fit_cofit(_conv(recipes), share=(), windows=windows,
                    iter_cb=iter_cb, tol=tol)

    fixed = tuple(p for p in all_but_amplitude(recipes) if p not in released)
    labels = [(r.sample or f"spectrum {k + 1}")
              for k, r in enumerate(res.recipes)]
    return BatchFitResult(
        recipes=res.recipes, labels=labels, rmsd=res.rmsd,
        per_dataset=res.per_dataset, shared=fixed, released=released,
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
