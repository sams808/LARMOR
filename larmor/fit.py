"""lmfit-based refinement of a LARMOR recipe against an experimental spectrum.

The residual reuses the precomputed Czjzek kernel, so a full fit runs in
seconds. Fitted values AND their standard errors are written back into the
recipe's Param objects -- the uncertainty dmfit never reported.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import lmfit

from larmor.engine import Axis, CzjzekKernel, build_kernel, needs_kernel, simulate_site
from larmor.recipe import Recipe

#: recipe param name -> lmfit-safe suffix
_PARAM_KEYS = {
    "isotropic_chemical_shift_ppm": "pos",
    "sigma_Cq_MHz": "sigma",
    "shift_fwhm_ppm": "fwhm",
    "amplitude": "amp",
    "gl": "gl",
}
_KEY_TO_PARAM = {v: k for k, v in _PARAM_KEYS.items()}

# user-facing constraint syntax: s<index>.<recipe param name>, e.g.
# "0.5 * s0.amplitude" -- translated to lmfit's internal "0.5 * s0_amp"
_EXPR_REF = re.compile(r"\bs(\d+)\.([A-Za-z_][A-Za-z0-9_]*)")


class ConstraintError(ValueError):
    """A constraint expression referenced an unknown site or parameter."""


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
        return f"s{i}_{_PARAM_KEYS[pname]}"

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
                f"s{i}_{_PARAM_KEYS[pname]}",
                value=p.value,
                vary=p.vary,
                min=p.min if p.min is not None else -np.inf,
                max=p.max if p.max is not None else np.inf,
            )
    # pass 2: attach constraint expressions
    for i, site in enumerate(recipe.sites):
        for pname, p in site.params.items():
            if p.expr:
                name = f"s{i}_{_PARAM_KEYS[pname]}"
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
            lp = params[f"s{i}_{_PARAM_KEYS[pname]}"]
            site.params[pname].value = float(lp.value)
            site.params[pname].stderr = (
                float(lp.stderr) if lp.stderr is not None else None)


def _model(recipe: Recipe, params: lmfit.Parameters, kernel: CzjzekKernel,
           ) -> tuple[np.ndarray, list[np.ndarray]]:
    _apply_params(recipe, params)
    per_site = [simulate_site(s, kernel) for s in recipe.sites]
    return np.sum(per_site, axis=0), per_site


def fit(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
        window_ppm: tuple[float, float] | None = None,
        kernel: "CzjzekKernel | Axis | None" = None) -> FitResult:
    """Refine `recipe` against (exp_ppm, exp_amp). Modifies recipe in place."""
    if kernel is None:
        if needs_kernel(recipe):
            kernel = build_kernel(recipe.nucleus, recipe.larmor_frequency_MHz,
                                  recipe.spin_rate_Hz)
        else:
            # analytic-only recipe (e.g. spin-1/2): fit straight on the data axis
            order = np.argsort(exp_ppm)
            kernel = Axis(x_ppm=np.asarray(exp_ppm)[order])
    window = window_ppm or recipe.fit_window_ppm
    if window is None:
        window = (float(np.max(exp_ppm)), float(np.min(exp_ppm)))
    hi, lo = max(window), min(window)
    sel = (exp_ppm >= lo) & (exp_ppm <= hi)
    xw, yw = exp_ppm[sel], exp_amp[sel]

    params = _make_params(recipe)

    # A site whose center lies outside the fit window is unconstrained by the
    # data (its parameters would wander and make the covariance singular), so
    # freeze it. Typical case: dmfit's ad-hoc Gauss/Lor sideband lines --
    # LARMOR's kernel simulates sidebands physically, so fitting those lines
    # again would double-count them anyway.
    frozen: list[str] = []
    for i, site in enumerate(recipe.sites):
        pos = site.params["isotropic_chemical_shift_ppm"].value
        if not (lo <= pos <= hi):
            # freeze only the FREE parameters; expression-linked ones follow
            # their master parameter and are already constrained
            for pname, p in site.params.items():
                if not p.expr:
                    params[f"s{i}_{_PARAM_KEYS[pname]}"].vary = False
            if not site.params["amplitude"].expr:
                # inert outside the window: zero its amplitude for this fit
                params[f"s{i}_amp"].value = 0.0
            frozen.append(site.label or f"site-{i}")

    # analytic global amplitude pre-scale so the optimizer starts on-scale
    y0, _ = _model(recipe, params, kernel)
    y0w = np.interp(xw, kernel.x_ppm, y0)
    denom = float(y0w @ y0w)
    if denom > 0:
        scale = float(yw @ y0w) / denom
        for i in range(len(recipe.sites)):
            if params[f"s{i}_amp"].vary:
                params[f"s{i}_amp"].value *= scale

    def residual(p):
        y, _ = _model(recipe, p, kernel)
        return np.interp(xw, kernel.x_ppm, y) - yw

    result = lmfit.minimize(residual, params, method="least_squares")

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
        amp_scale = max((abs(retry_params[f"s{i}_amp"].value)
                         for i in range(len(recipe.sites))), default=1.0)
        for i, site in enumerate(recipe.sites):
            if abs(retry_params[f"s{i}_amp"].value) <= 1e-6 * amp_scale:
                for pname in site.params:
                    retry_params[f"s{i}_{_PARAM_KEYS[pname]}"].vary = False
        retry = lmfit.minimize(residual, retry_params, method="leastsq")
        if retry.errorbars:
            result = retry
    _apply_params(recipe, result.params)

    # user-facing names, e.g. "s0.sigma_Cq_MHz"
    at_bounds = []
    for n in at_bounds_internal:
        i, suffix = n.split("_", 1)
        at_bounds.append(f"{i}.{_KEY_TO_PARAM[suffix]}")

    y_fit, per_site = _model(recipe, result.params, kernel)
    y_fit_w = np.interp(xw, kernel.x_ppm, y_fit)
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
    return FitResult(recipe=recipe, lmfit_result=result, x_ppm=kernel.x_ppm,
                     y_exp=exp_amp, y_fit=y_fit, per_site=per_site, rmsd=rmsd,
                     frozen_sites=frozen, at_bounds=at_bounds)
