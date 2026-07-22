"""lmfit-based refinement of a LARMOR recipe against an experimental spectrum.

Parameter names, bounds and defaults come from the model registry
(larmor.models); the residual reuses cached kernels/lineshapes so a fit runs
in seconds. Fitted values AND their standard errors are written back into the
recipe's Param objects, and constraint expressions (Param.expr) are honored
with full error propagation.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import lmfit

from larmor import models as model_registry
from larmor.engine import make_context, simulate_site
from larmor.recipe import Recipe

# user-facing constraint syntax: s<index>.<recipe param name>, e.g.
# "0.5 * s0.amplitude" -- translated to lmfit's internal "0.5 * s0_amp"
_EXPR_REF = re.compile(r"\bs(\d+)\.([A-Za-z_][A-Za-z0-9_]*)")


class ConstraintError(ValueError):
    """A constraint expression referenced an unknown site or parameter."""


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
        for pname, p in site.params.items():
            params.add(
                _lmfit_name(i, site, pname),
                value=p.value,
                vary=p.vary,
                min=p.min if p.min is not None else -np.inf,
                max=p.max if p.max is not None else np.inf,
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
        kernel=None, iter_cb=None) -> FitResult:
    """Refine `recipe` against (exp_ppm, exp_amp). Modifies recipe in place.

    `kernel` is accepted for backward compatibility and ignored; kernels are
    cached process-wide and resolved automatically.
    """
    ctx = make_context(recipe, exp_ppm=exp_ppm)
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

    result = lmfit.minimize(residual, params, method="least_squares",
                            iter_cb=iter_cb)

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
    if not result.errorbars:
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

    y_fit, per_site = _model(recipe, result.params, ctx)
    y_fit_w = np.interp(xw, ctx.x_ppm, y_fit)
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
    return FitResult(recipe=recipe, lmfit_result=result, x_ppm=ctx.x_ppm,
                     y_exp=exp_amp, y_fit=y_fit, per_site=per_site, rmsd=rmsd,
                     frozen_sites=frozen, at_bounds=at_bounds)
