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
import time
from dataclasses import dataclass, field

import numpy as np

from larmor import fit as fitmod
from larmor.recipe import Recipe

#: ``parallel="auto"``: the process pool is worth starting only when the
#: sequential work it would replace is longer than its own start-up. That
#: start-up is ~15 s on Windows (a spawned interpreter per worker, each
#: importing mrsimulator and loading a kernel) -- measured on a single-line
#: 23Na Czjzek profile: 16.8 s through the pool against 2.2 s sequentially,
#: with nothing on screen for the first 15 s. The first item always runs in
#: this process and is timed; the rest go to the pool when their estimated
#: sequential time exceeds the threshold (a warm pool costs only dispatch).
POOL_BREAKEVEN_S = 20.0
POOL_BREAKEVEN_WARM_S = 4.0


def _map_items(fn, items: list, *, parallel, max_workers=None, should_stop=None,
               on_result=None, executor=None, heartbeat=None) -> tuple[list, bool]:
    """``parallel_map`` with the ``"auto"`` strategy on top. Returns
    ``(results, ran_parallel)``. ``parallel`` is True / False / "auto"."""
    from larmor import parallel as par

    if parallel != "auto" or not items:
        use = bool(parallel)
        res = par.parallel_map(fn, items, max_workers=max_workers,
                               should_stop=should_stop, on_result=on_result,
                               use_processes=use, executor=executor,
                               heartbeat=heartbeat)
        return res, use and (executor is not None
                             or len(items) >= par.MIN_ITEMS_FOR_PROCESSES)
    results: list = [None] * len(items)
    if should_stop is not None and should_stop():
        return results, False
    t0 = time.perf_counter()
    try:
        first = fn(items[0])
    except Exception:
        first = None
    t_first = time.perf_counter() - t0
    results[0] = first
    if on_result:
        on_result(0, first)
    if heartbeat is not None:
        heartbeat()
    rest = items[1:]
    if not rest:
        return results, False
    threshold = (POOL_BREAKEVEN_WARM_S if (executor is not None or par.pool_is_warm())
                 else POOL_BREAKEVEN_S)
    use = len(rest) * t_first > threshold

    def _shifted(i, r):
        if on_result:
            on_result(i + 1, r)

    results[1:] = par.parallel_map(fn, rest, max_workers=max_workers,
                                   should_stop=should_stop, on_result=_shifted,
                                   use_processes=use, executor=executor,
                                   heartbeat=heartbeat)
    return results, use and (executor is not None
                             or len(rest) >= par.MIN_ITEMS_FOR_PROCESSES)


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
    #: degrees of freedom of the full fit (points in the window minus free
    #: parameters) and the residual variance chi2_min / dof it implies --
    #: the unit the delta-chi-square rule is stated in
    dof: int = 0
    noise_var: float = 0.0
    level68: float = 0.0            # chi2 at the 1-sigma crossing
    level95: float = 0.0            # chi2 at the 2-sigma (95 %) crossing
    #: Σ residual² of the recipe as it was handed in (no refit), and whether
    #: it sits within the 1σ level of the scan's minimum -- False means the
    #: scan found a deeper minimum than the fit on screen
    chi2_fit: float | None = None
    fit_at_minimum: bool = True
    #: did the scan points run through the process pool (``parallel="auto"``
    #: decides from the first point's duration)
    ran_parallel: bool = False

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


def _profile_point_worker(item):
    """One χ² profile scan point: fix `param` at `value`, refit everything
    else free, return chisqr. Module-level (not a closure) so it can be
    pickled and sent to a worker process -- see larmor/parallel.py."""
    base_json, site, param, value, exp_ppm, exp_amp, window_ppm = item
    trial = Recipe.from_dict(json.loads(base_json))
    tp = trial.sites[site].params[param]
    tp.value = float(value)
    tp.vary = False                   # fixed at the scan point
    tp.expr = None
    try:
        # only chisqr is read below -- this scan point's own covariance/
        # error bars are never used (the WHOLE profile's shape is the
        # error estimate), so skip the errorbar-rescue retry
        res = fitmod.fit(trial, exp_ppm, exp_amp, window_ppm=window_ppm,
                         compute_errorbars=False)
        return float(res.lmfit_result.chisqr)
    except Exception:
        return None


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
                  progress=None, should_stop=None, parallel="auto",
                  max_workers: int | None = None, executor=None,
                  heartbeat=None) -> ErrorProfile:
    """chi-square profile of one parameter (dmfit's Errors Analysis).

    The parameter is fixed at each scanned value while EVERY other free
    parameter is re-fitted, so correlations are absorbed rather than ignored.
    `span` = how many stderr (or 25% of the value if no stderr) to scan each
    way. Confidence intervals come from the delta-chi-square rule for one
    parameter of interest: 1.00 for 1σ, 3.84 for 2σ (95%) -- in units of
    the residual variance. The fit minimises an UNWEIGHTED sum of squares,
    so chi2 carries the data's intensity units; the rule holds for
    residuals divided by the noise sigma, which is unknown and is
    estimated the way lmfit scales its covariance: chi2_min / dof, with
    dof = points in the window - free parameters. The levels are therefore
    chi2_min (1 + 1.00 / dof) and chi2_min (1 + 3.84 / dof). Stating them
    as chi2_min + 1.00 / + 3.84 (as this function did until 0.12.1) puts
    them a part in 1e13 above the minimum on real-intensity data, so the
    interpolated crossings collapsed onto the best value and the tool
    reported a zero-width interval for every parameter.

    Every scan point is an independent refit -- ``parallel=True`` runs them
    across a process pool (larmor.parallel) instead of one at a time, which
    is where nearly all of this function's time goes for a nucleus/model
    whose fit itself isn't instant; ``"auto"`` (the default) times the first
    point in this process and uses the pool only when the remaining points
    would take longer than the pool's own start-up (``POOL_BREAKEVEN_S``);
    ``False`` never leaves this process. ``executor``: reuse an
    already-running pool instead of starting one just for this call -- pass
    one when calling this many times back-to-back (see
    ``batchfit.batch_error_analysis``). ``heartbeat()`` is called while
    waiting (see ``parallel.parallel_map``) so a dialog can stay alive.

    ``notes`` says when the scan found a LOWER χ² than the recipe's own
    (``fit_at_minimum`` False): the reported best value and interval then
    belong to that deeper minimum, and the fit on screen should be run
    again before its numbers are quoted.
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

    notes = []
    items = [(base, site, param, float(v), exp_ppm, exp_amp, window_ppm)
             for v in values]

    def _cb(i, _r):
        if progress:
            progress(i + 1, n_points, float(values[i]))

    raw, ran_parallel = _map_items(_profile_point_worker, items, parallel=parallel,
                                   max_workers=max_workers, should_stop=should_stop,
                                   on_result=_cb, executor=executor,
                                   heartbeat=heartbeat)
    chi2 = np.array([np.nan if c is None else c for c in raw])
    ok = np.isfinite(chi2)
    if ok.sum() < 3:
        raise RuntimeError("chi-square profile failed: too few valid points")
    values, chi2 = values[ok], chi2[ok]

    chi2_min = float(np.min(chi2))
    best = float(values[int(np.argmin(chi2))])
    dof = _profile_dof(recipe, exp_ppm, window_ppm)
    noise_var = chi2_min / dof if chi2_min > 0 else 0.0
    level68 = chi2_min + 1.00 * noise_var
    level95 = chi2_min + 3.84 * noise_var
    ci68 = _crossings(values, chi2, level68, best)
    ci95 = _crossings(values, chi2, level95, best)
    step_v = abs(float(values[1] - values[0])) if values.size > 1 else 0.0
    if ci68[0] is None or ci68[1] is None:
        notes.append("1σ not bracketed — widen `span`; the parameter may be "
                     "poorly determined")
    elif ci68 == (best, best) or (ci68[1] - ci68[0]) < step_v:
        # both crossings lie between the minimum and its neighbours: the
        # interval is pure interpolation across one step; a narrower span
        # resolves it (the old 5 %-of-a-step test let 90 % of such cases
        # through unflagged)
        notes.append("1σ interval narrower than the scan step — rerun with "
                     "a smaller `span` to resolve it")
    # the profile refits everything else at every point: when it finds a
    # lower χ² than the fit reached, the fit was not at the minimum, and the
    # best value / interval reported here belong to the DEEPER minimum, not
    # to the numbers on screen. Say so rather than let the two disagree in
    # silence (seen on a single-line 23Na Czjzek fit: σ(Cq) 2.16 in the
    # table, 1.35 from the profile, χ² 16 % lower)
    chi2_fit = _chi2_at(recipe, exp_ppm, exp_amp, window_ppm)
    fit_at_minimum = True
    if chi2_fit is not None and chi2_fit > level68:
        fit_at_minimum = False
        drop = 100.0 * (chi2_fit - chi2_min) / chi2_fit if chi2_fit else 0.0
        notes.append(f"the scan reached a χ² {drop:.1f} % below the fit's own "
                     f"({param} = {best:.4g} there): the fit had not converged "
                     "— Fit again (Auto fit finds it), then rescan")
    return ErrorProfile(site=site, param=param, values=values, chi2=chi2,
                        best_value=best, chi2_min=chi2_min,
                        ci68=ci68, ci95=ci95, notes=notes, dof=dof,
                        noise_var=noise_var, level68=level68,
                        level95=level95, chi2_fit=chi2_fit,
                        fit_at_minimum=fit_at_minimum,
                        ran_parallel=ran_parallel)


def _chi2_at(recipe: Recipe, exp_ppm, exp_amp, window_ppm) -> float | None:
    """Σ residual² of ``recipe`` AS GIVEN (no refit) over the fit window --
    the χ² the profile's levels are compared with. None when the model
    cannot be simulated."""
    from larmor import engine

    try:
        x = np.asarray(exp_ppm, float)
        y = np.asarray(exp_amp, float)
        window_ppm = window_ppm or recipe.fit_window_ppm
        if window_ppm is not None and len(window_ppm) == 2:
            hi, lo = max(window_ppm), min(window_ppm)
            sel = (x >= lo) & (x <= hi)
        else:
            sel = np.ones(x.shape, bool)
        mx, model, _ = engine.simulate(recipe, exp_ppm=x)
        yi = np.interp(x[sel], np.asarray(mx, float), np.asarray(model, float))
        return float(np.sum((yi - y[sel]) ** 2))
    except Exception:                                    # noqa: BLE001
        return None


def _profile_dof(recipe: Recipe, exp_ppm: np.ndarray,
                 window_ppm: tuple[float, float] | None) -> int:
    """Degrees of freedom of the full fit: data points inside the fit
    window minus the parameters it varies (linked ones are not free)."""
    x = np.asarray(exp_ppm, float)
    window_ppm = window_ppm or recipe.fit_window_ppm     # as fit.fit resolves it
    if window_ppm is not None and len(window_ppm) == 2:
        hi, lo = max(window_ppm), min(window_ppm)
        n_pts = int(np.count_nonzero((x >= lo) & (x <= hi)))
    else:
        n_pts = int(x.size)
    n_free = sum(1 for s in recipe.sites for p in s.params.values()
                 if p.vary and not p.expr)
    return max(n_pts - n_free, 1)


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
    #: per-trial window integrals of EVERY site (n_ok x n_sites, draw order),
    #: re-integrated on each refit so family sums / ratios take their error
    #: from the spread of per-trial sums (larmor.families); None when a
    #: trial's integration failed
    site_integrals: np.ndarray | None = None
    #: the (hi, lo) ppm window those integrals were taken over
    window_ppm: tuple | None = None
    #: did the trials run through the process pool (``parallel="auto"``)
    ran_parallel: bool = False
    #: parameters whose refitted best differs from the recipe's value by
    #: more than the Monte-Carlo σ: the fit on screen was not at the
    #: minimum the trials scatter around (labels)
    moved: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        s = (f"Monte-Carlo errors from {self.n_ok}/{self.trials} synthetic "
             f"refits · noise σ = {self.noise:.4g}")
        if self.moved:
            s += (" · the refit moved " + ", ".join(self.moved)
                  + " by more than σ: the fit on screen had not converged — "
                  "Fit again, then rerun")
        return s

    def report(self) -> str:
        lines = [self.summary, ""]
        w = max((len(p.label) for p in self.params), default=8)
        for p in self.params:
            pc = f"{p.pct:.2f}%" if np.isfinite(p.pct) else "—"
            lines.append(f"{p.label:<{w}}  {p.mean:12.6g} ± {p.std:.4g}   ({pc})")
        return "\n".join(lines)


def _mc_trial_worker(item):
    """One Monte-Carlo trial: refit the model against ONE synthetic noisy
    spectrum, return ``(fitted value of every tracked free parameter, the
    per-site window integrals of the refit or None)``. The integrals feed the
    Monte-Carlo basis of the family sums (larmor.families): amplitude is the
    PEAK HEIGHT for most models, so summing per-trial amplitudes would be
    wrong -- every trial is re-integrated instead. Module-level (not a
    closure) so it can be pickled and sent to a worker process -- see
    larmor/parallel.py."""
    base_json, exp_ppm, synth_amp, window_ppm, tracked = item
    trial = Recipe.from_dict(json.loads(base_json))
    try:
        # the MC estimate IS the spread of .value across trials -- each
        # trial's own covariance/error bars are never read, so skip the
        # errorbar-rescue retry (worth up to 2x on every one of n_trials)
        fitmod.fit(trial, exp_ppm, synth_amp, window_ppm=window_ppm,
                  compute_errorbars=False)
    except Exception:
        return None
    values = {(i, pn): float(trial.sites[i].params[pn].value)
              for i, pn, _ in tracked}
    try:
        from larmor.quantify import site_integrals
        integrals = [float(v) for v in site_integrals(trial, window_ppm)[0]]
    except Exception:
        integrals = None            # disables the basis, not the run
    return values, integrals


def monte_carlo_errors(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                       window_ppm: tuple[float, float] | None = None,
                       n_trials: int = 200, seed: int = 0,
                       noise: float | None = None, progress=None,
                       should_stop=None, parallel="auto",
                       max_workers: int | None = None,
                       executor=None, heartbeat=None) -> MonteCarloResult:
    """Estimate parameter errors by Monte-Carlo (synthetic-noise refits).

    The recipe is fitted once to fix the best fit and estimate the noise level
    (residual std over the window, unless `noise` is given). Then `n_trials`
    synthetic spectra = best-fit model + Gaussian(0, noise) are each re-fitted
    from the best fit; the std of each free parameter over the trials is its
    error. `progress(k, n_trials)` is called per completed trial (k counts
    completions, not draw order, when `parallel=True`); `should_stop()`
    truthy aborts early (returns what was collected).

    Every trial is an independent refit -- ``parallel=True`` runs them across
    a process pool (larmor.parallel) instead of one at a time; the synthetic
    noise draws themselves stay a single up-front SEQUENTIAL loop over `rng`
    so the trial set (and therefore the result, given a fixed seed) doesn't
    depend on how work happens to be scheduled across worker processes.
    ``"auto"`` (the default) times the first trial in this process and uses
    the pool only when the rest would outlast its start-up
    (``POOL_BREAKEVEN_S``). ``executor``: reuse an already-running pool
    instead of starting one just for this call. ``heartbeat()`` is called
    while waiting so a dialog can stay alive.

    The best fit the trials scatter around is a REFIT of ``recipe``; when
    that refit moves a parameter by more than the Monte-Carlo σ the result's
    ``moved`` names it (and ``summary`` says so): the numbers on screen were
    not at the minimum.
    """
    from larmor import engine

    rng = np.random.default_rng(seed)
    exp_ppm = np.asarray(exp_ppm, float)
    exp_amp = np.asarray(exp_amp, float)

    # 1. lock the best fit + its model on the experimental axis. A kernel model
    # (Czjzek family) simulates on its OWN grid regardless of exp_ppm, so the
    # model must be interpolated onto exp_ppm before it can be compared to
    # exp_amp — the same pattern larmor.fit uses for its residual.
    given = {(i, pn): float(p.value) for i, s in enumerate(recipe.sites)
             for pn, p in s.params.items() if p.vary and not p.expr}
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

    # 4. MC trials — refit model + synthetic noise, starting from the best fit.
    # The noise draws are generated up front in one sequential pass over
    # `rng` (NOT inside the possibly-parallel worker) so the trial set --
    # and therefore the whole result, for a fixed seed -- is independent of
    # execution order/worker scheduling.
    base_best = json.dumps(best.to_dict())
    synths = [model + rng.normal(0.0, sigma, size=model.shape)
             for _ in range(n_trials)]
    items = [(base_best, exp_ppm, synth, window_ppm, tracked) for synth in synths]

    done = 0

    def _cb(_i, _r):
        nonlocal done
        done += 1
        if progress:
            progress(done, n_trials)

    raw, ran_parallel = _map_items(_mc_trial_worker, items, parallel=parallel,
                                   max_workers=max_workers, should_stop=should_stop,
                                   on_result=_cb, executor=executor,
                                   heartbeat=heartbeat)
    n_ok = 0
    integral_rows: list = []
    for r in raw:
        if r is None:
            continue
        values, ints = r
        for key, val in values.items():
            collected[key].append(val)
        integral_rows.append(ints)
        n_ok += 1
    # per-trial site integrals, stacked in draw order (parallel_map returns
    # item order) when every ok trial has them; the window is whatever the
    # trials integrated -- the fit's own (hi, lo), resolved exactly as the
    # worker's site_integrals() did
    site_ints = None
    if integral_rows and all(x is not None for x in integral_rows):
        try:
            site_ints = np.asarray(integral_rows, float)
        except ValueError:
            site_ints = None
    mc_window = None
    try:
        from larmor.quantify import site_integrals as _site_integrals
        mc_window = _site_integrals(best, window_ppm)[1]
    except Exception:
        site_ints = None

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
    # the refit's best against the values handed in: a move beyond the MC σ
    # means the trials scatter around a minimum the screen does not show
    moved = [label for i, pn, label in tracked
             if (i, pn) in given and np.isfinite(best_vals[(i, pn)])
             and abs(best_vals[(i, pn)] - given[(i, pn)])
             > max(next((p.std for p in params if p.label == label), 0.0),
                   1e-12 * max(1.0, abs(given[(i, pn)])))
             and np.isfinite(next((p.std for p in params if p.label == label),
                                  float("nan")))]
    return MonteCarloResult(trials=n_trials, n_ok=n_ok, noise=sigma, seed=seed,
                            params=params, site_integrals=site_ints,
                            window_ppm=mc_window, ran_parallel=ran_parallel,
                            moved=moved)
