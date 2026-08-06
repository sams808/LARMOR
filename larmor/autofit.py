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
            # only chisqr is read below -- this scan point's own covariance/
            # error bars are never used (the WHOLE profile's shape is the
            # error estimate), so skip the errorbar-rescue retry
            res = fitmod.fit(trial, exp_ppm, exp_amp, window_ppm=window_ppm,
                             compute_errorbars=False)
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


# --------------------------------------------------------------------------
# Monte-Carlo errors (dmfit "Errors ▸ Monte Carlo"; pydmfit errorsMonteCarlo.py)
#
# A parametric bootstrap: take the best fit, add synthetic Gaussian noise at the
# residual level to the *model*, re-fit, and repeat N times. The spread of each
# parameter across the trials is its uncertainty. Unlike the covariance matrix
# this captures non-linearity and parameter correlations; unlike the χ² profile
# it does every parameter at once and yields a full distribution (histogram).
# dmfit/pydmfit report each parameter as mean ± σ with σ = sqrt(var) and a
# percentage σ/mean·100 — reproduced here.

@dataclass
class MCParam:
    site: int
    param: str
    label: str                      # e.g. "s0.Cq_MHz"
    best: float                     # best-fit value
    mean: float                     # mean over the MC trials
    std: float                      # sqrt(var) over the trials (the MC error)
    values: np.ndarray = field(default_factory=lambda: np.empty(0))

    @property
    def pct(self) -> float:
        return abs(self.std / self.mean) * 100.0 if self.mean else float("nan")


@dataclass
class MonteCarloResult:
    trials: int
    n_ok: int
    noise: float                    # σ of the synthetic noise (data units)
    seed: int
    params: list[MCParam] = field(default_factory=list)

    @property
    def summary(self) -> str:
        return (f"Monte-Carlo errors from {self.n_ok}/{self.trials} synthetic "
                f"refits · noise σ = {self.noise:.4g}")

    def report(self) -> str:
        lines = [self.summary, ""]
        w = max((len(p.label) for p in self.params), default=8)
        for p in self.params:
            pc = f"{p.pct:.2f}%" if np.isfinite(p.pct) else "—"
            lines.append(f"{p.label:<{w}}  {p.mean:12.6g} ± {p.std:.4g}   ({pc})")
        return "\n".join(lines)


def monte_carlo_errors(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                       window_ppm: tuple[float, float] | None = None,
                       n_trials: int = 200, seed: int = 0,
                       noise: float | None = None, progress=None,
                       should_stop=None) -> MonteCarloResult:
    """Estimate parameter errors by Monte-Carlo (synthetic-noise refits).

    The recipe is fitted once to fix the best fit and estimate the noise level
    (residual std over the window, unless `noise` is given). Then `n_trials`
    synthetic spectra = best-fit model + Gaussian(0, noise) are each re-fitted
    from the best fit; the std of each free parameter over the trials is its
    error. `progress(k, n_trials)` is called per trial; `should_stop()` truthy
    aborts early (returns what was collected).
    """
    from larmor import engine

    rng = np.random.default_rng(seed)
    exp_ppm = np.asarray(exp_ppm, float)
    exp_amp = np.asarray(exp_amp, float)

    # 1. lock the best fit + its model on the experimental axis. A kernel model
    # (Czjzek family) simulates on its OWN grid regardless of exp_ppm, so the
    # model must be interpolated onto exp_ppm before it can be compared to
    # exp_amp — the same pattern larmor.fit uses for its residual.
    best = Recipe.from_dict(json.loads(json.dumps(recipe.to_dict())))
    fitmod.fit(best, exp_ppm, exp_amp, window_ppm=window_ppm)
    mx, model_raw, _ = engine.simulate(best, exp_ppm=exp_ppm)
    model = np.interp(exp_ppm, mx, np.asarray(model_raw, float))

    # 2. noise level: residual std inside the fit window
    if window_ppm:
        lo, hi = min(window_ppm), max(window_ppm)
        m = (exp_ppm >= lo) & (exp_ppm <= hi)
    else:
        m = np.ones(exp_ppm.shape, bool)
    resid = exp_amp[m] - model[m]
    sigma = float(noise) if noise else float(np.std(resid))
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.std(exp_amp)) * 1e-3 or 1.0

    # 3. free parameters to track
    tracked = [(i, pn, f"s{i}.{pn}")
               for i, s in enumerate(best.sites)
               for pn, p in s.params.items() if p.vary and not p.expr]
    best_vals = {(i, pn): float(best.sites[i].params[pn].value)
                 for i, pn, _ in tracked}
    collected: dict = {(i, pn): [] for i, pn, _ in tracked}

    # 4. MC trials — refit model + synthetic noise, starting from the best fit
    base_best = json.dumps(best.to_dict())
    n_ok = 0
    for k in range(n_trials):
        if should_stop is not None and should_stop():
            break
        synth = model + rng.normal(0.0, sigma, size=model.shape)
        trial = Recipe.from_dict(json.loads(base_best))
        try:
            # the MC estimate IS the spread of .value across trials -- each
            # trial's own covariance/error bars are never read, so skip the
            # errorbar-rescue retry (worth up to 2x on every one of n_trials)
            fitmod.fit(trial, exp_ppm, synth, window_ppm=window_ppm,
                      compute_errorbars=False)
        except Exception:
            if progress:
                progress(k + 1, n_trials)
            continue
        for i, pn, _ in tracked:
            collected[(i, pn)].append(float(trial.sites[i].params[pn].value))
        n_ok += 1
        if progress:
            progress(k + 1, n_trials)

    # 5. per-parameter statistics (mean ± sqrt(var), matching dmfit/pydmfit)
    params = []
    for i, pn, label in tracked:
        vals = np.asarray(collected[(i, pn)], float)
        if vals.size:
            mean = float(np.mean(vals))
            std = float(np.sqrt(np.var(vals)))
        else:
            mean, std = best_vals[(i, pn)], float("nan")
        params.append(MCParam(site=i, param=pn, label=label,
                              best=best_vals[(i, pn)], mean=mean, std=std,
                              values=vals))
    return MonteCarloResult(trials=n_trials, n_ok=n_ok, noise=sigma, seed=seed,
                            params=params)
