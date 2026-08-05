"""Sequential (forward–backward) fitting of a spectral series.

The everyday case is a composition/temperature series whose lineshape parameters
evolve **smoothly** from one end-member to the other. Instead of tying every
spectrum to one shared model (see :mod:`larmor.batchfit`), you fit one spectrum,
then **carry its fitted parameters forward** as the starting point for the next,
fit that, and so on to the far end-member — then optionally sweep **back** to
smooth the trajectory. Each spectrum keeps its own fit; the neighbour only warm-
starts it, so parameters track the physics of the series without being forced
equal.

Qt-free and testable; the desktop dialog drives it interactively and with an
auto "N passes forward–backward" mode plus optional trajectory smoothing.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from larmor import fit as fitmod
from larmor.batchfit import all_but_amplitude
from larmor.recipe import Recipe


def seed_from(dst: Recipe, src: Recipe, params=None) -> None:
    """Copy fitted parameter **values** from ``src`` into ``dst`` (matched by site
    index and name), clipped into ``dst``'s own bounds. ``dst`` keeps its bounds,
    vary flags, links and its own nucleus/Larmor. ``params=None`` → all params."""
    for i, site in enumerate(dst.sites):
        if i >= len(src.sites):
            break
        src_site = src.sites[i]
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
    back — used between passes so the next sweep starts from a smoother guess."""
    if window < 3 or len(order) < 3:
        return
    nsite = len(recipes[order[0]].sites)
    for i in range(nsite):
        for pn in params:
            try:
                vals = np.array([recipes[k].sites[i].params[pn].value
                                 for k in order], float)
            except (KeyError, IndexError):
                continue
            sm = _smooth(vals, window)
            for j, k in enumerate(order):
                p = recipes[k].sites[i].params.get(pn)
                if p is not None and not p.expr:
                    p.value = float(sm[j])


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

    @property
    def summary(self) -> str:
        trend = ""
        if len(self.history) >= 2:
            trend = (f" · mean RMSD {self.history[0]['mean']:.4g} → "
                     f"{self.history[-1]['mean']:.4g}")
        return (f"sequential fit: {len(self.recipes)} spectra, {self.passes} pass"
                f"{'es' if self.passes != 1 else ''}{trend}")


def run_sequential(entries, *, passes: int = 2, start: str = "first",
                   propagate=None, smooth: int = 0, tol=None,
                   progress=None, should_stop=None) -> SeqFitResult:
    """Fit a series by warm-starting each spectrum from its fitted neighbour.

    ``entries`` = list of ``(recipe, ppm, amp, window)`` in **series order** (each
    recipe already carries the model to fit). Passes alternate direction
    (forward, backward, …); ``start`` picks the first direction. ``propagate``
    defaults to *all but amplitude* (positions/widths/quadrupolar carry; each
    amplitude is re-fit fresh). ``smooth`` (window ≥ 3) smooths the parameter
    trajectories between passes. ``progress(pass, k, rmsd)`` fires after each
    spectrum; ``should_stop()`` aborts between spectra.
    """
    if len(entries) < 2:
        raise ValueError("sequential fit needs at least two spectra")
    recipes = [e[0] for e in entries]
    ppms = [np.asarray(e[1], float) for e in entries]
    amps = [np.asarray(e[2], float) for e in entries]
    windows = [e[3] for e in entries]
    n = len(entries)
    if propagate is None:
        propagate = all_but_amplitude(recipes)
    propagate = tuple(propagate)
    base = list(range(n)) if start == "first" else list(range(n - 1, -1, -1))

    history = []
    stopped = False
    for p in range(max(1, passes)):
        order = base if p % 2 == 0 else base[::-1]
        prev = None
        rmsds = [float("nan")] * n
        for k in order:
            if prev is not None:
                seed_from(recipes[k], recipes[prev], propagate)
            fitmod.fit(recipes[k], ppms[k], amps[k],
                       window_ppm=windows[k], tol=tol)
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
                        propagated=propagate)
