"""lmfit-based refinement of a LARMOR recipe against an experimental spectrum.

The residual reuses the precomputed Czjzek kernel, so a full fit runs in
seconds. Fitted values AND their standard errors are written back into the
recipe's Param objects -- the uncertainty dmfit never reported.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import lmfit

from larmor.engine import CzjzekKernel, build_kernel, simulate_site
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

    @property
    def report(self) -> str:
        return lmfit.fit_report(self.lmfit_result)


def _make_params(recipe: Recipe) -> lmfit.Parameters:
    params = lmfit.Parameters()
    for i, site in enumerate(recipe.sites):
        for pname, p in site.params.items():
            params.add(
                f"s{i}_{_PARAM_KEYS[pname]}",
                value=p.value,
                vary=p.vary,
                min=p.min if p.min is not None else -np.inf,
                max=p.max if p.max is not None else np.inf,
            )
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
        kernel: CzjzekKernel | None = None) -> FitResult:
    """Refine `recipe` against (exp_ppm, exp_amp). Modifies recipe in place."""
    if kernel is None:
        kernel = build_kernel(recipe.nucleus, recipe.larmor_frequency_MHz,
                              recipe.spin_rate_Hz)
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
            for pname in site.params:
                params[f"s{i}_{_PARAM_KEYS[pname]}"].vary = False
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
    if not result.errorbars:
        # covariance didn't come out of least_squares; Levenberg-Marquardt from
        # the solution usually recovers it
        retry = lmfit.minimize(residual, result.params, method="leastsq")
        if retry.errorbars:
            result = retry
    _apply_params(recipe, result.params)

    y_fit, per_site = _model(recipe, result.params, kernel)
    y_fit_w = np.interp(xw, kernel.x_ppm, y_fit)
    rmsd = float(np.sqrt(np.mean((y_fit_w - yw) ** 2)) / (yw.max() or 1.0))

    recipe.fit_window_ppm = (hi, lo)
    recipe.fit_rmsd = rmsd
    if frozen:
        note = f"sites frozen (center outside fit window {hi}..{lo} ppm): " + ", ".join(frozen)
        if note not in recipe.notes:
            recipe.notes.append(note)
    return FitResult(recipe=recipe, lmfit_result=result, x_ppm=kernel.x_ppm,
                     y_exp=exp_amp, y_fit=y_fit, per_site=per_site, rmsd=rmsd,
                     frozen_sites=frozen)
