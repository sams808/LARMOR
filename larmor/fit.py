"""lmfit-based refinement of a LARMOR recipe against an experimental spectrum.

Parameter names, bounds and defaults come from the model registry
(larmor.models); the residual reuses cached kernels/lineshapes so a fit runs
in seconds. Fitted values AND their standard errors are written back into the
recipe's Param objects, and constraint expressions (Param.expr) are honored
with full error propagation.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass

import numpy as np
import lmfit

from larmor import models as model_registry
from larmor.engine import (grid_restrictable, make_context, simulate_site,
                           site_width_margin)
from larmor.recipe import Recipe

# user-facing constraint syntax: s<index>.<recipe param name>, e.g.
# "0.5 * s0.amplitude" -- translated to lmfit's internal "0.5 * s0_amp"
_EXPR_REF = re.compile(r"\bs(\d+)\.([A-Za-z_][A-Za-z0-9_]*)")


class ConstraintError(ValueError):
    """A constraint expression referenced an unknown site or parameter."""


def ftol_from_pct(pct) -> float | None:
    """Map a user 'completion threshold' — the % change in the residual stdev
    below which the fit is considered converged — to a scipy least_squares
    ``ftol``. Since stdev ∝ √cost, a relative change ``δ`` in stdev is ~½ the
    relative change in the cost, so ``ftol ≈ 2·pct/100``. ``None``/0 → solver
    default (fit to machine tolerance)."""
    try:
        p = float(pct)
    except (TypeError, ValueError):
        return None
    if p <= 0:
        return None
    return max(1e-15, 2.0 * p / 100.0)


def _tol_kws(tol) -> dict:
    """least_squares stop tolerances from a completion threshold (``ftol`` on the
    cost, matched ``xtol``); empty when no threshold is set."""
    ft = ftol_from_pct(tol)
    return {} if ft is None else {"ftol": ft, "xtol": ft}


def _key(site, pname: str) -> str:
    return model_registry.get(site.model).key_of(pname)


def _lmfit_name(i: int, site, pname: str) -> str:
    return f"s{i}_{_key(site, pname)}"


def translate_expr(expr: str, recipe: Recipe) -> str:
    """Turn 's0.amplitude'-style references into lmfit parameter names."""

    def repl(m: re.Match) -> str:
        i, pname = int(m.group(1)), m.group(2)
        if i >= len(recipe.sites):
            raise ConstraintError(
                f"constraint {expr!r}: site s{i} does not exist "
                f"(recipe has {len(recipe.sites)} sites)")
        if pname not in recipe.sites[i].params:
            valid = ", ".join(recipe.sites[i].params)
            raise ConstraintError(
                f"constraint {expr!r}: s{i} has no parameter {pname!r} "
                f"(valid: {valid})")
        return _lmfit_name(i, recipe.sites[i], pname)

    return _EXPR_REF.sub(repl, expr)


@dataclass
class FitResult:
    recipe: Recipe
    lmfit_result: lmfit.minimizer.MinimizerResult
    x_ppm: np.ndarray
    y_exp: np.ndarray
    y_fit: np.ndarray
    per_site: list[np.ndarray]
    rmsd: float
    frozen_sites: list[str] = None
    #: user-facing names (s0.sigma_Cq_MHz, ...) of parameters that finished
    #: the fit pinned at a min/max bound -- usually a sign that a constraint
    #: or starting model is fighting the data
    at_bounds: list[str] = None

    @property
    def report(self) -> str:
        return lmfit.fit_report(self.lmfit_result)


def _make_params(recipe: Recipe) -> lmfit.Parameters:
    params = lmfit.Parameters()
    # pass 1: every parameter exists as a plain value, so that pass-2
    # expressions can reference any of them regardless of site order
    for i, site in enumerate(recipe.sites):
        try:
            pdefs = {pd.name: pd for pd in model_registry.get(site.model).params}
        except Exception:
            pdefs = {}
        for pname, p in site.params.items():
            # fall back to the MODEL's physical bounds when the recipe omits them,
            # so e.g. a width can never be fitted negative even if an old recipe
            # saved it with min=None (a released width otherwise blows up)
            pd = pdefs.get(pname)
            pmin = p.min if p.min is not None else (pd.min if pd else None)
            pmax = p.max if p.max is not None else (pd.max if pd else None)
            # bounds are meaningless for a FIXED parameter (lmfit never moves
            # it) but lmfit's own Parameter still validates min != max
            # unconditionally -- widen to unbounded rather than crash on a
            # degenerate recipe bound like min=max=0 (the "exclude this
            # component" signature batch fit writes onto a locked amplitude,
            # or simply an old recipe with an accidental min==max)
            if not p.vary:
                pmin, pmax = None, None
            params.add(
                _lmfit_name(i, site, pname),
                value=p.value,
                vary=p.vary,
                min=pmin if pmin is not None else -np.inf,
                max=pmax if pmax is not None else np.inf,
            )
    # pass 2: attach constraint expressions
    for i, site in enumerate(recipe.sites):
        for pname, p in site.params.items():
            if p.expr:
                name = _lmfit_name(i, site, pname)
                try:
                    params[name].expr = translate_expr(p.expr, recipe)
                    # evaluate now so a broken expression fails loudly here,
                    # not deep inside the minimizer
                    params.update_constraints()
                except ConstraintError:
                    raise
                except Exception as exc:
                    raise ConstraintError(
                        f"constraint {p.expr!r} on s{i}.{pname} is invalid: {exc}"
                    ) from exc
    return params


def _apply_params(recipe: Recipe, params: lmfit.Parameters) -> None:
    for i, site in enumerate(recipe.sites):
        for pname in site.params:
            lp = params[_lmfit_name(i, site, pname)]
            site.params[pname].value = float(lp.value)
            site.params[pname].stderr = (
                float(lp.stderr) if lp.stderr is not None else None)


def _model(recipe: Recipe, params: lmfit.Parameters, ctx,
           ) -> tuple[np.ndarray, list[np.ndarray]]:
    _apply_params(recipe, params)
    per_site = [simulate_site(s, ctx) for s in recipe.sites]
    return np.sum(per_site, axis=0), per_site


def fit(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
        window_ppm: tuple[float, float] | None = None,
        kernel=None, iter_cb=None, tol=None, frame_cb=None,
        frame_every: int = 10, compute_errorbars: bool = True) -> FitResult:
    """Refine `recipe` against (exp_ppm, exp_amp). Modifies recipe in place.

    `kernel` is accepted for backward compatibility and ignored; kernels are
    cached process-wide and resolved automatically. `tol` is the completion
    threshold (% change in the residual stdev, see ``ftol_from_pct``).
    `frame_cb(x_ppm, y_model, iteration)`, if given, is called each iteration with
    the current full model curve — for a live fit animation.

    `compute_errorbars`: when the primary least_squares pass can't get a valid
    covariance (an ill-conditioned Jacobian — a parameter at/near a bound, or
    degenerate with another free one), the default (True) reruns the WHOLE
    optimization a second time with a different algorithm to try to rescue
    error bars. That is a good tradeoff for one interactive fit but doubles the
    cost of every caller that never reads the resulting stderr anyway (batch
    fit's initial pass, each Monte-Carlo trial, each chi-square-profile scan
    point) — pass False there.
    """
    zones = [z for z in (recipe.fit_zones or []) if z and len(z) == 2]
    if zones:
        # dmfit-style Zones: residual evaluated on the union of the regions
        sel = np.zeros(exp_ppm.shape, dtype=bool)
        for zhi, zlo in zones:
            sel |= (exp_ppm >= min(zhi, zlo)) & (exp_ppm <= max(zhi, zlo))
        hi = max(max(z) for z in zones)
        lo = min(min(z) for z in zones)
    else:
        window = window_ppm or recipe.fit_window_ppm
        if window is None:
            window = (float(np.max(exp_ppm)), float(np.min(exp_ppm)))
        hi, lo = max(window), min(window)
        sel = (exp_ppm >= lo) & (exp_ppm <= hi)
    xw, yw = exp_ppm[sel], exp_amp[sel]

    # a full-grid context is always needed for the RETURNED model curve (so the
    # displayed fit still spans the whole experiment, not just the fit window
    # -- unchanged from before). Models that tolerate it (engine.
    # grid_restrictable -- pointwise formulas, or ones that simulate on their
    # OWN independent grid like the Czjzek kernel family) also get a SEPARATE
    # context restricted to the window + a generous margin, used ONLY during
    # optimisation -- so every one of the thousands of Jacobian-probe
    # evaluations simulates a fraction of the points instead of the whole
    # spectrum. Models excluded from the allowlist (quad_ct/quad_first/
    # quad_csa/csa_*) keep simulating on the full grid exactly as before.
    ctx_full = make_context(recipe, exp_ppm=exp_ppm)
    if grid_restrictable(recipe):
        margin = site_width_margin(recipe.sites)
        sel_pad = (exp_ppm >= lo - margin) & (exp_ppm <= hi + margin)
        ctx = make_context(recipe, exp_ppm=exp_ppm[sel_pad])
    else:
        ctx = ctx_full

    params = _make_params(recipe)

    # A site whose center lies outside the fit window is unconstrained by the
    # data (its parameters would wander and make the covariance singular), so
    # freeze it. Expression-linked parameters follow their master and stay.
    frozen: list[str] = []
    for i, site in enumerate(recipe.sites):
        center = site.params.get("isotropic_chemical_shift_ppm")
        if center is None:
            continue          # spans the window (e.g. a background spectrum)
        pos = center.value
        if not (lo <= pos <= hi):
            for pname, p in site.params.items():
                if not p.expr:
                    params[_lmfit_name(i, site, pname)].vary = False
            if not site.params["amplitude"].expr:
                params[_lmfit_name(i, site, "amplitude")].value = 0.0
            frozen.append(site.label or f"site-{i}")

    # analytic global amplitude pre-scale so the optimizer starts on-scale
    y0, _ = _model(recipe, params, ctx)
    y0w = np.interp(xw, ctx.x_ppm, y0)
    denom = float(y0w @ y0w)
    if denom > 0:
        scale = float(yw @ y0w) / denom
        for i, site in enumerate(recipe.sites):
            amp = params[_lmfit_name(i, site, "amplitude")]
            if amp.vary:
                amp.value *= scale

    def residual(p):
        y, _ = _model(recipe, p, ctx)
        return np.interp(xw, ctx.x_ppm, y) - yw

    _fs = {"n": 0}

    def _main_cb(p, it, resid, *a, **k):
        # emit the live model curve only every `frame_every` iterations (and the
        # first couple) — computing/redrawing every iteration would slow the fit;
        # the FINAL model is always drawn by the caller when the fit finishes
        _fs["n"] += 1
        if frame_cb is not None and (_fs["n"] <= 2
                                     or _fs["n"] % max(1, frame_every) == 0):
            try:
                # the full-grid context, not the (possibly window-restricted)
                # fit context -- already throttled to every frame_every-th
                # iteration, so the extra width costs nothing noticeable, and
                # the animated curve's span then matches the final result's
                ym, _ = _model(recipe, p, ctx_full)
                frame_cb(np.asarray(ctx_full.x_ppm, float), np.asarray(ym, float),
                         _fs["n"])
            except Exception:
                pass
        return iter_cb(p, it, resid, *a, **k) if iter_cb is not None else None

    with warnings.catch_warnings():
        # lmfit computes a covariance-based stderr for every free parameter
        # INSIDE minimize() itself, unconditionally -- with several
        # correlated/near-degenerate parameters (routine for overlapping
        # quadrupolar sites, and the norm rather than the exception once a
        # χ² profile or Monte-Carlo pass has most parameters FIXED at each
        # point) that covariance is frequently not positive-definite, and
        # lmfit sqrt()s its diagonal regardless, printing a RuntimeWarning
        # per ill-conditioned parameter. Harmless (this fit's own use of
        # that stderr is a NaN either way, and neither the caller nor
        # compute_errorbars ever reads it when it isn't wanted) but at
        # hundreds-to-thousands of fits per error-analysis run, purely
        # noisy -- silence just these two known-benign warnings (the stderr
        # sqrt, and the correlation divide that follows it for the same
        # ill-conditioned covariance), narrowly scoped to this call so an
        # unrelated RuntimeWarning elsewhere still surfaces.
        warnings.filterwarnings("ignore", message="invalid value encountered in sqrt",
                                category=RuntimeWarning)
        warnings.filterwarnings("ignore",
                                message="invalid value encountered in scalar divide",
                                category=RuntimeWarning)
        result = lmfit.minimize(residual, params, method="least_squares",
                                iter_cb=(_main_cb if (frame_cb or iter_cb) else None),
                                **_tol_kws(tol))

    def _at_bounds(res) -> list[str]:
        names = []
        for n, p in res.params.items():
            if not p.vary:
                continue
            span = max(1.0, abs(p.value))
            if (np.isfinite(p.min) and (p.value - p.min) < 1e-3 * span) or \
               (np.isfinite(p.max) and (p.max - p.value) < 1e-3 * span):
                names.append(n)
        return names

    at_bounds_internal = _at_bounds(result)
    if compute_errorbars and not result.errorbars:
        # covariance didn't come out of least_squares; Levenberg-Marquardt from
        # the solution usually recovers it
        retry_params = result.params.copy()
        # parameters pinned at a bound have a one-sided derivative that breaks
        # the covariance -- hold them and report errors conditional on that
        for n in at_bounds_internal:
            retry_params[n].vary = False
        # a site whose amplitude collapsed to ~zero leaves its remaining
        # parameters without any influence on the residual: pin the whole
        # site for the covariance pass
        amp_names = [_lmfit_name(i, s, "amplitude")
                     for i, s in enumerate(recipe.sites)]
        amp_scale = max((abs(retry_params[n].value) for n in amp_names),
                        default=1.0)
        for i, site in enumerate(recipe.sites):
            if abs(retry_params[amp_names[i]].value) <= 1e-6 * amp_scale:
                for pname in site.params:
                    retry_params[_lmfit_name(i, site, pname)].vary = False
        retry = lmfit.minimize(residual, retry_params, method="leastsq",
                               iter_cb=iter_cb)
        if retry.errorbars:
            result = retry
    _apply_params(recipe, result.params)

    # user-facing names, e.g. "s0.sigma_Cq_MHz"
    key_to_name = {}
    for i, site in enumerate(recipe.sites):
        for pname in site.params:
            key_to_name[_lmfit_name(i, site, pname)] = f"s{i}.{pname}"
    at_bounds = [key_to_name.get(n, n) for n in at_bounds_internal]

    # the returned curve always spans the FULL experimental axis, whether or
    # not the optimisation itself ran on a restricted grid
    y_fit, per_site = _model(recipe, result.params, ctx_full)
    y_fit_w = np.interp(xw, ctx_full.x_ppm, y_fit)
    rmsd = float(np.sqrt(np.mean((y_fit_w - yw) ** 2)) / (yw.max() or 1.0))

    recipe.fit_window_ppm = (hi, lo)
    recipe.fit_rmsd = rmsd
    if frozen:
        note = f"sites frozen (center outside fit window {hi}..{lo} ppm): " + ", ".join(frozen)
        if note not in recipe.notes:
            recipe.notes.append(note)
    if at_bounds:
        note = ("parameters finished at a bound (check constraints/starting "
                "model; uncertainties are conditional on them): " + ", ".join(at_bounds))
        if note not in recipe.notes:
            recipe.notes.append(note)
    return FitResult(recipe=recipe, lmfit_result=result, x_ppm=ctx_full.x_ppm,
                     y_exp=exp_amp, y_fit=y_fit, per_site=per_site, rmsd=rmsd,
                     frozen_sites=frozen, at_bounds=at_bounds)
