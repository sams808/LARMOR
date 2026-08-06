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
    # error analysis (filled by batch_error_analysis; covariance is the default
    # from the fit itself). error_detail maps a method name -> per-spectrum
    # {(site, param): ParamError}.
    error_method: str = "covariance"
    error_detail: dict = field(default_factory=dict)

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


# --------------------------------------------------------------------------
# Per-spectrum error analysis. The covariance stderr comes free with the fit;
# a stronger estimate (Monte-Carlo parametric bootstrap, or a χ² profile — the
# same tools as the single-fit dialogs, in larmor.autofit) can be computed on
# top and exported. Each method's result is stored per spectrum so the user can
# switch which one they export without losing the others.

@dataclass
class ParamError:
    site: int
    param: str
    label: str                            # "s0.amplitude"
    value: float
    stderr: float | None
    ci68: tuple = (None, None)            # (lo, hi) from a χ² profile, else Nones
    pct: float | None = None              # |stderr / value| * 100


def _free_params(rec: Recipe):
    """(site, param) of every free, non-linked parameter of a fitted recipe —
    for a batch these are the amplitudes plus any released shape parameters."""
    return [(i, pn) for i, s in enumerate(rec.sites)
            for pn, p in s.params.items()
            if p.vary and not getattr(p, "expr", None)]


def _snapshot_covariance(result: BatchFitResult) -> list[dict]:
    """Capture the least-squares covariance stderr currently on the recipes,
    so it survives a later Monte-Carlo/profile pass overwriting Param.stderr."""
    out = []
    for rec in result.recipes:
        d = {}
        for i, pn in _free_params(rec):
            p = rec.sites[i].params[pn]
            se = p.stderr
            pct = abs(se / p.value) * 100.0 if se and p.value else None
            d[(i, pn)] = ParamError(i, pn, f"s{i}.{pn}", float(p.value), se,
                                    (None, None), pct)
        out.append(d)
    return out


def batch_error_analysis(result: BatchFitResult, data: list[tuple], *,
                         method: str = "montecarlo", n_trials: int = 200,
                         seed: int = 0, n_points: int = 15, span: float = 3.0,
                         progress=None, should_stop=None) -> BatchFitResult:
    """Estimate per-spectrum parameter errors on a finished batch fit.

    ``data`` is ``[(ppm, amp, window), ...]`` aligned with ``result.recipes``.
    ``method``:

    * ``"covariance"`` — the least-squares covariance stderr already on the fit
      (no refit; just recorded as the export method).
    * ``"montecarlo"`` — parametric bootstrap: refit ``n_trials`` synthetic
      noisy copies per spectrum (``autofit.monte_carlo_errors``). Captures
      correlations and non-linearity.
    * ``"profile"`` — χ² profile per free parameter: scan it, refit the rest,
      read a real 1σ interval off the χ² curve (``autofit.error_profile``).

    The chosen method's error is written back into each recipe's ``Param.stderr``
    (so saved fits carry it) and stored in ``result.error_detail[method]``.
    ``progress(k, n, j, tot)`` = spectrum k of n, sub-step j of tot.
    """
    from larmor import autofit

    method = {"mc": "montecarlo", "monte-carlo": "montecarlo",
              "errors": "profile", "error": "profile",
              "chi2": "profile"}.get(method.lower(), method.lower())

    # preserve the covariance errors before any overwrite, so it stays exportable
    if "covariance" not in result.error_detail:
        result.error_detail["covariance"] = _snapshot_covariance(result)

    if method == "covariance":
        result.error_method = "covariance"
        return result

    n = len(result.recipes)
    detail: list[dict] = []
    for k, rec in enumerate(result.recipes):
        if should_stop is not None and should_stop():
            break
        ppm, amp, window = data[k]
        rows: dict = {}
        if method == "montecarlo":
            def pcb(j, tot, _k=k):
                if progress:
                    progress(_k, n, j, tot)
            mc = autofit.monte_carlo_errors(
                rec, ppm, amp, window_ppm=window, n_trials=n_trials, seed=seed,
                progress=pcb, should_stop=should_stop)
            for mp in mc.params:
                rec.sites[mp.site].params[mp.param].stderr = mp.std
                rows[(mp.site, mp.param)] = ParamError(
                    mp.site, mp.param, mp.label, mp.best, mp.std, (None, None),
                    mp.pct)
        elif method == "profile":
            free = _free_params(rec)
            for j, (i, pn) in enumerate(free):
                if should_stop is not None and should_stop():
                    break
                if progress:
                    progress(k, n, j, len(free))
                try:
                    ep = autofit.error_profile(rec, ppm, amp, window_ppm=window,
                                               site=i, param=pn,
                                               n_points=n_points, span=span)
                except Exception:
                    continue
                lo, hi = ep.ci68
                se = (hi - lo) / 2.0 if (lo is not None and hi is not None) else None
                rec.sites[i].params[pn].stderr = se
                pct = abs(se / ep.best_value) * 100.0 if se and ep.best_value else None
                rows[(i, pn)] = ParamError(i, pn, f"s{i}.{pn}", ep.best_value,
                                           se, ep.ci68, pct)
        else:
            raise ValueError(f"unknown error method: {method!r}")
        detail.append(rows)
        if progress:
            progress(k + 1, n, 0, 1)

    result.error_detail[method] = detail
    result.error_method = method
    return result


def error_table(result: BatchFitResult, method: str | None = None) -> list[dict]:
    """Per-parameter rows for a CSV export carrying the selected error method.

    Shared (held) parameters appear once with no error; the per-spectrum free
    parameters (amplitudes + released) carry value, stderr, %-error, and — for a
    χ² profile — the 1σ interval, using ``method``'s stored errors (falling back
    to the covariance stderr currently on the recipes)."""
    method = (method or result.error_method or "covariance")
    detail = result.error_detail.get(method)
    if detail is None:
        detail = (_snapshot_covariance(result) if method == "covariance"
                  else [{} for _ in result.recipes])

    rows: list[dict] = []
    master = result.recipes[0]
    for i, site in enumerate(master.sites):        # shared / held params, once
        for pn, p in site.params.items():
            if pn in result.shared:
                rows.append({"scope": "shared", "site": f"s{i}",
                             "label": site.label or site.model, "param": pn,
                             "value": p.value, "stderr": None, "sigma_pct": None,
                             "ci68_lo": None, "ci68_hi": None,
                             "error_method": method})
    for k, rec in enumerate(result.recipes):       # per-spectrum free params
        d = detail[k] if k < len(detail) else {}
        for i, site in enumerate(rec.sites):
            for pn, p in site.params.items():
                if pn == "amplitude" or pn in result.released:
                    pe = d.get((i, pn))
                    stderr = pe.stderr if pe else p.stderr
                    lo, hi = pe.ci68 if (pe and pe.ci68) else (None, None)
                    pct = (pe.pct if pe else
                           (abs(stderr / p.value) * 100.0
                            if stderr and p.value else None))
                    rows.append({"scope": result.labels[k], "site": f"s{i}",
                                 "label": site.label or site.model, "param": pn,
                                 "value": p.value, "stderr": stderr,
                                 "sigma_pct": pct, "ci68_lo": lo, "ci68_hi": hi,
                                 "error_method": method})
    return rows


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
