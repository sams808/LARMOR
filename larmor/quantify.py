"""Quantification: per-site integrals and relative populations.

This is dmfit's results table -- the thing that actually gets published for a
glass: each site's integrated intensity as a percentage of the total, with an
uncertainty.

Error model: for a fixed lineshape, a site's integral is proportional to its
amplitude parameter, so the fractional error on the integral is taken from the
fractional error on the amplitude (first-order approximation; shape-parameter
covariance is neglected and this is stated in the output).
"""
from __future__ import annotations

import numpy as np

from larmor.engine import make_context, simulate_site
from larmor.recipe import Recipe


def quantify(recipe: Recipe, window_ppm: tuple[float, float] | None = None,
             ) -> dict:
    """Integrate every site over the window. Returns a JSON-friendly table."""
    ctx = make_context(recipe)
    window = window_ppm or recipe.fit_window_ppm or \
        (float(ctx.x_ppm.max()), float(ctx.x_ppm.min()))
    hi, lo = max(window), min(window)
    sel = (ctx.x_ppm >= lo) & (ctx.x_ppm <= hi)

    rows = []
    for i, site in enumerate(recipe.sites):
        y = simulate_site(site, ctx)
        integral = float(np.trapezoid(y[sel], ctx.x_ppm[sel]))
        amp = site.params["amplitude"]
        rel_err = (amp.stderr / amp.value) if (amp.stderr and amp.value) else None
        pos = site.params["isotropic_chemical_shift_ppm"]
        rows.append({
            "site": f"s{i}",
            "label": site.label or site.model,
            "model": site.model,
            "position_ppm": pos.value,
            "position_err": pos.stderr,
            "integral": integral,
            "integral_err": abs(integral) * rel_err if rel_err is not None else None,
        })

    total = sum(abs(r["integral"]) for r in rows) or 1.0
    for r in rows:
        r["fraction_pct"] = 100.0 * abs(r["integral"]) / total
        r["fraction_err_pct"] = (
            100.0 * r["integral_err"] / total
            if r["integral_err"] is not None else None)

    return {
        "window_ppm": [hi, lo],
        "rows": rows,
        "note": "fraction errors are first-order (amplitude covariance only; "
                "lineshape-parameter covariance neglected)",
    }
