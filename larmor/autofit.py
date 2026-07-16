"""Auto Fit and Errors Analysis -- the two dmfit Decomposition tools that
turn a fit from "a number" into "a defended number".

Auto Fit (dmfit Decomposition > Auto Fit): the fit landscape of overlapping
quadrupolar sites is riddled with local minima, so a single gradient run from
one starting guess proves nothing. Restart from many randomized starts inside
the bounds, keep the best chi-square.

Errors Analysis (dmfit Decomposition > Errors Analysis): the covariance matrix
assumes a locally quadratic, well-conditioned chi-square. For strongly
correlated parameters (sigma_Cq vs shift_fwhm, amplitudes of overlapping
lines) that assumption breaks. Scan a parameter across a range, re-fitting
everything else at each step, and read the confidence interval off the real
chi-square profile.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

from larmor import fit as fitmod
from larmor.recipe import Recipe


@dataclass
class AutoFitResult:
    recipe: Recipe                  # best recipe found (also modified in place)
    best_rmsd: float
    trials: list[float]             # rmsd of every trial, best first
    n_improved: int                 # how many restarts beat the initial fit
    result: object = None           # the winning FitResult

    @property
    def summary(self) -> str:
        return (f"best RMSD {self.best_rmsd:.5f} over {len(self.trials)} "
                f"starts ({self.n_improved} beat the plain fit)")


def _perturb(recipe: Recipe, rng: np.random.Generator, spread: float) -> None:
    """Randomize every free parameter around its value, inside its bounds."""
    for site in recipe.sites:
        for p in site.params.values():
            if not p.vary or p.expr:
                continue
            lo = p.min if p.min is not None else -np.inf
            hi = p.max if p.max is not None else np.inf
            scale = abs(p.value) * spread if p.value else spread
            if np.isfinite(lo) and np.isfinite(hi):
                scale = min(scale, 0.5 * (hi - lo))
            val = p.value + rng.normal(0.0, scale or spread)
            p.value = float(np.clip(val, lo + 1e-9 if np.isfinite(lo) else val,
                                    hi - 1e-9 if np.isfinite(hi) else val))


def auto_fit(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
             window_ppm: tuple[float, float] | None = None,
             n_starts: int = 12, spread: float = 0.25, seed: int = 0,
             progress=None) -> AutoFitResult:
    """Multi-start fit. Returns the best recipe; the input recipe is updated.

    `progress(i, n, rmsd_best)` is called after every trial when given.
    """
    rng = np.random.default_rng(seed)
    base = json.dumps(recipe.to_dict())

    # trial 0: the fit from the user's own starting point
    best_result = fitmod.fit(recipe, exp_ppm, exp_amp, window_ppm=window_ppm)
    best_rmsd = best_result.rmsd
    best_dict = json.dumps(recipe.to_dict())
    trials = [best_rmsd]
    n_improved = 0
    if progress:
        progress(1, n_starts + 1, best_rmsd)

    for i in range(n_starts):
        trial = Recipe.from_dict(json.loads(base))
        _perturb(trial, rng, spread)
        try:
            res = fitmod.fit(trial, exp_ppm, exp_amp, window_ppm=window_ppm)
        except Exception:
            continue                     # a wild start can be unsimulatable
        trials.append(res.rmsd)
        if res.rmsd < best_rmsd - 1e-12:
            best_rmsd = res.rmsd
            best_dict = json.dumps(trial.to_dict())
            best_result = res
            n_improved += 1
        if progress:
            progress(i + 2, n_starts + 1, best_rmsd)

    # write the winner back into the caller's recipe object
    winner = json.loads(best_dict)
    recipe.sites = Recipe.from_dict(winner).sites
    recipe.fit_rmsd = best_rmsd
    recipe.fit_window_ppm = winner.get("fit_window_ppm")
    note = (f"auto fit: best of {len(trials)} starts, RMSD {best_rmsd:.5f}"
            + (f"; {n_improved} restart(s) beat the initial fit -- the "
               "landscape has local minima" if n_improved else
               "; no restart improved on the initial fit"))
    if note not in recipe.notes:
        recipe.notes.append(note)
    return AutoFitResult(recipe=recipe, best_rmsd=best_rmsd,
                         trials=sorted(trials), n_improved=n_improved,
                         result=best_result)


@dataclass
class ErrorProfile:
    site: int
    param: str
    values: np.ndarray              # scanned values
    chi2: np.ndarray                # chi-square at each (others re-fitted)
    best_value: float
    chi2_min: float
    ci68: tuple[float | None, float | None]
    ci95: tuple[float | None, float | None]
    notes: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        lo, hi = self.ci68
        if lo is None and hi is None:
            return (f"s{self.site}.{self.param} = {self.best_value:.4g} "
                    "(1σ interval not bracketed in the scanned range)")
        lo_s = f"{lo:.4g}" if lo is not None else "<scan"
        hi_s = f"{hi:.4g}" if hi is not None else ">scan"
        return (f"s{self.site}.{self.param} = {self.best_value:.4g} "
                f"[1σ: {lo_s} … {hi_s}]")


def _crossings(x: np.ndarray, y: np.ndarray, level: float, best: float):
    """Where the profile crosses `level`, on each side of the minimum."""
    lo = hi = None
    imin = int(np.argmin(y))
    # left branch
    for i in range(imin, 0, -1):
        if y[i - 1] >= level >= y[i]:
            f = (level - y[i]) / (y[i - 1] - y[i] or 1.0)
            lo = float(x[i] + f * (x[i - 1] - x[i]))
            break
    # right branch
    for i in range(imin, len(x) - 1):
        if y[i + 1] >= level >= y[i]:
            f = (level - y[i]) / (y[i + 1] - y[i] or 1.0)
            hi = float(x[i] + f * (x[i + 1] - x[i]))
            break
    return lo, hi


def error_profile(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                  site: int, param: str,
                  window_ppm: tuple[float, float] | None = None,
                  n_points: int = 15, span: float = 3.0,
                  progress=None) -> ErrorProfile:
    """chi-square profile of one parameter (dmfit's Errors Analysis).

    The parameter is fixed at each scanned value while EVERY other free
    parameter is re-fitted, so correlations are absorbed rather than ignored.
    `span` = how many stderr (or 25% of the value if no stderr) to scan each
    way. Confidence intervals come from the delta-chi-square rule for one
    parameter of interest: 1.00 for 1σ, 3.84 for 2σ (95%).
    """
    base = json.dumps(recipe.to_dict())
    p0 = recipe.sites[site].params[param]
    center = p0.value
    step = p0.stderr if p0.stderr else abs(center) * 0.25 or 0.25
    lo_v, hi_v = center - span * step, center + span * step
    if p0.min is not None:
        lo_v = max(lo_v, p0.min)
    if p0.max is not None:
        hi_v = min(hi_v, p0.max)
    values = np.linspace(lo_v, hi_v, n_points)

    chi2 = []
    notes = []
    for k, v in enumerate(values):
        trial = Recipe.from_dict(json.loads(base))
        tp = trial.sites[site].params[param]
        tp.value = float(v)
        tp.vary = False               # fixed at the scan point
        tp.expr = None
        try:
            res = fitmod.fit(trial, exp_ppm, exp_amp, window_ppm=window_ppm)
            chi2.append(float(res.lmfit_result.chisqr))
        except Exception:
            chi2.append(np.nan)
        if progress:
            progress(k + 1, n_points, float(v))
    chi2 = np.array(chi2)
    ok = np.isfinite(chi2)
    if ok.sum() < 3:
        raise RuntimeError("chi-square profile failed: too few valid points")
    values, chi2 = values[ok], chi2[ok]

    chi2_min = float(np.min(chi2))
    best = float(values[int(np.argmin(chi2))])
    ci68 = _crossings(values, chi2, chi2_min + 1.00, best)
    ci95 = _crossings(values, chi2, chi2_min + 3.84, best)
    if ci68[0] is None or ci68[1] is None:
        notes.append("1σ not bracketed — widen `span`; the parameter may be "
                     "poorly determined")
    return ErrorProfile(site=site, param=param, values=values, chi2=chi2,
                        best_value=best, chi2_min=chi2_min,
                        ci68=ci68, ci95=ci95, notes=notes)
