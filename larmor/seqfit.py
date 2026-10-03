"""Sequential (forward–backward) fitting of a spectral series.

The everyday case is a composition/temperature series whose lineshape parameters
evolve **smoothly** from one end-member to the other. Instead of tying every
spectrum to one shared model (see :mod:`larmor.batchfit`), you fit one spectrum,
then **carry its fitted parameters forward** as the starting point for the next,
fit that, and so on to the far end-member — then optionally sweep **back** to
smooth the trajectory. Each spectrum keeps its own fit; the neighbour only warm-
starts it, so parameters track the physics of the series without being forced
equal.

The members need not share one model. Lines are paired across spectra by their
**label** (:mod:`larmor.components`): a component the user removed on one
spectrum or added on another is simply absent there — nothing is carried onto
it, and the trajectory smoother skips it. Members listed in ``fixed`` are
**kept**: the sweep starts their neighbours from them but never refits them.

Qt-free and testable; the series mode of the main window (``mw_series``) and
the CLI ``larmor seqfit`` both call :func:`run_sequential`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from larmor import fit as fitmod
from larmor.batchfit import all_but_amplitude
from larmor.components import component_map, pair_sites
from larmor.recipe import Recipe


def seed_from(dst: Recipe, src: Recipe, params=None) -> list:
    """Copy fitted parameter **values** from ``src`` into ``dst`` on the lines
    that pair (:func:`larmor.components.pair_sites`: by label, by index only
    when a label cannot identify its line), clipped into ``dst``'s own bounds.
    ``dst`` keeps its bounds, vary flags, links and its own nucleus/Larmor.
    ``params=None`` → every parameter; ``()`` → nothing. Returns the pairing
    (per ``dst`` line the ``src`` index or None)."""
    pairs = pair_sites(dst.sites, src.sites)
    for i, site in enumerate(dst.sites):
        j = pairs[i]
        if j is None:
            continue
        src_site = src.sites[j]
        for pn, p in site.params.items():
            if p.expr:                          # linked — follows its master
                continue
            if params is not None and pn not in params:
                continue
            sp = src_site.params.get(pn)
            if sp is None:
                continue
            v = float(sp.value)
            if p.min is not None and np.isfinite(p.min):
                v = max(v, p.min)
            if p.max is not None and np.isfinite(p.max):
                v = min(v, p.max)
            p.value = v
    return pairs


def _rmsd(recipe: Recipe, ppm, amp, window) -> float:
    from larmor import engine
    x, total, _ = engine.simulate(recipe, exp_ppm=np.asarray(ppm, float))
    yi = np.interp(ppm, x, total)
    d = amp - yi
    if window:
        lo, hi = min(window), max(window)
        sel = (ppm >= lo) & (ppm <= hi)
        if sel.any():
            d = d[sel]
    return float(np.sqrt(np.mean(np.asarray(d, float) ** 2)))


def _smooth(vals: np.ndarray, window: int) -> np.ndarray:
    """Odd-window moving average with edge padding (a light trajectory smoother)."""
    vals = np.asarray(vals, float)
    w = int(window)
    if w < 3 or vals.size < 3:
        return vals
    w = min(w, vals.size)
    if w % 2 == 0:
        w += 1
    pad = w // 2
    ext = np.pad(vals, pad, mode="edge")
    return np.convolve(ext, np.ones(w) / w, mode="valid")[:vals.size]


def smooth_trajectories(recipes, params, order, window: int) -> None:
    """Smooth each parameter's value across the series (in ``order``) and write it
    back — used between passes so the next sweep starts from a smoother guess.
    Components are paired by label (:func:`larmor.components.component_map`);
    one that only some members carry is smoothed over those members, and one
    present in fewer than three is left alone."""
    if window < 3 or len(order) < 3:
        return
    comps = component_map([recipes[k] for k in order])
    for c in comps:
        where = [(j, k) for j, k in enumerate(order) if c["index"][j] is not None]
        if len(where) < 3:
            continue
        for pn in params:
            vals = []
            for j, k in where:
                p = recipes[k].sites[c["index"][j]].params.get(pn)
                if p is None:
                    vals = None
                    break
                vals.append(p.value)
            if vals is None:
                continue
            sm = _smooth(np.array(vals, float), window)
            for (j, k), v in zip(where, sm):
                p = recipes[k].sites[c["index"][j]].params.get(pn)
                if p is not None and not p.expr:
                    p.value = float(v)


@dataclass
class SeqFitResult:
    recipes: list                    # per-spectrum fitted recipes (final)
    labels: list
    rmsd: list                       # final RMSD per spectrum (series order)
    per_dataset: list                # {"x","y_fit"} per spectrum
    history: list                    # per-pass {"pass","direction","rmsd","mean"}
    passes: int
    propagated: tuple                # parameters carried between spectra
    warnings: list = field(default_factory=list)
    fixed: tuple = ()                # members kept as seeds, never refitted

    @property
    def summary(self) -> str:
        trend = ""
        if len(self.history) >= 2:
            trend = (f" · mean RMSD {self.history[0]['mean']:.4g} → "
                     f"{self.history[-1]['mean']:.4g}")
        kept = f", {len(self.fixed)} kept" if self.fixed else ""
        return (f"sequential fit: {len(self.recipes)} spectra{kept}, {self.passes} pass"
                f"{'es' if self.passes != 1 else ''}{trend}")


def run_sequential(entries, *, passes: int = 2, start: str = "first",
                   propagate=None, smooth: int = 0, tol=None,
                   progress=None, should_stop=None, fixed=None) -> SeqFitResult:
    """Fit a series by warm-starting each spectrum from its fitted neighbour.

    ``entries`` = list of ``(recipe, ppm, amp, window)`` in **series order** (each
    recipe already carries the model to fit). Passes alternate direction
    (forward, backward, …); ``start`` picks the first direction. ``propagate``
    defaults to *all but amplitude* (positions/widths/quadrupolar carry; each
    amplitude is re-fit fresh); an empty tuple carries nothing, so every
    spectrum simply refits from its own model. ``smooth`` (window ≥ 3) smooths
    the parameter trajectories between passes. ``fixed`` lists the series
    indices to keep: their recipes are never refitted or seeded, but they seed
    their neighbours like any fitted member (their RMSD is still measured).
    ``progress(pass, k, rmsd)`` fires after each spectrum; ``should_stop()``
    aborts between spectra.
    """
    if len(entries) < 2:
        raise ValueError("sequential fit needs at least two spectra")
    recipes = [e[0] for e in entries]
    ppms = [np.asarray(e[1], float) for e in entries]
    amps = [np.asarray(e[2], float) for e in entries]
    windows = [e[3] for e in entries]
    n = len(entries)
    fixed = frozenset(int(k) for k in (fixed or ()) if 0 <= int(k) < n)
    if len(fixed) >= n:
        raise ValueError("every spectrum is kept — release at least one to fit")
    if propagate is None:
        propagate = all_but_amplitude(recipes)
    propagate = tuple(propagate)
    base = list(range(n)) if start == "first" else list(range(n - 1, -1, -1))

    # an iter_cb that aborts the individual lmfit run the moment a stop is asked
    # for — so Stop is responsive DURING a fit, not only between spectra (a single
    # slow fit otherwise makes the sweep feel unstoppable)
    def _abort_cb(params, it, resid, *a, **k):
        return True if (should_stop is not None and should_stop()) else None

    history = []
    stopped = False
    for p in range(max(1, passes)):
        order = base if p % 2 == 0 else base[::-1]
        prev = None
        rmsds = [float("nan")] * n
        for k in order:
            if should_stop is not None and should_stop():
                stopped = True
                break
            if k in fixed:
                # kept: measured as it is, and the seed of the next one
                r = _rmsd(recipes[k], ppms[k], amps[k], windows[k])
                rmsds[k] = r
                if progress is not None:
                    progress(p, k, r)
                prev = k
                continue
            if prev is not None:
                seed_from(recipes[k], recipes[prev], propagate)
            fitmod.fit(recipes[k], ppms[k], amps[k],
                       window_ppm=windows[k], tol=tol, iter_cb=_abort_cb)
            r = _rmsd(recipes[k], ppms[k], amps[k], windows[k])
            rmsds[k] = r
            if progress is not None:
                progress(p, k, r)
            prev = k
            if should_stop is not None and should_stop():
                stopped = True
                break
        finite = [v for v in rmsds if np.isfinite(v)]
        history.append({"pass": p, "direction": "→" if order[0] < order[-1] else "←",
                        "rmsd": rmsds,
                        "mean": float(np.mean(finite)) if finite else float("nan")})
        if stopped:
            break
        if smooth and p < passes - 1:
            smooth_trajectories(recipes, propagate, base, smooth)

    from larmor import engine
    per, final_rmsd = [], []
    for k in range(n):
        x, total, _ = engine.simulate(recipes[k], exp_ppm=ppms[k])
        per.append({"x": x, "y_fit": total})
        final_rmsd.append(_rmsd(recipes[k], ppms[k], amps[k], windows[k]))
    labels = [(r.sample or f"spectrum {k + 1}") for k, r in enumerate(recipes)]
    return SeqFitResult(recipes=recipes, labels=labels, rmsd=final_rmsd,
                        per_dataset=per, history=history, passes=max(1, passes),
                        propagated=propagate, fixed=tuple(sorted(fixed)))
