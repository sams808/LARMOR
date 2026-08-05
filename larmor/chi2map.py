"""χ² surface over a pair of fitted parameters.

Scanning two parameters on a grid (holding the rest at their fitted values) and
mapping the residual sum-of-squares shows *which* parameters the data actually
determines: a round basin = well-determined, a long diagonal valley = the pair
trades off (degenerate), exactly what the correlation number quantifies. Qt-free
and testable.
"""
from __future__ import annotations

import numpy as np

from larmor import engine
from larmor.recipe import Recipe


def varying_params(recipe: dict) -> list[tuple]:
    """(site_index, param_name, label) for every non-fixed, non-linked, non-gl
    parameter — the candidates for a χ² map axis."""
    out = []
    for i, s in enumerate(recipe.get("sites", [])):
        for pn, p in s.get("params", {}).items():
            if pn == "gl":
                continue
            if isinstance(p, dict) and (p.get("vary") is False or p.get("expr")):
                continue
            out.append((i, pn, f"s{i} {s.get('label') or s.get('model')}: {pn}"))
    return out


def chi2_surface(recipe: dict, exp_ppm, exp_amp, window,
                 axis_a: tuple, axis_b: tuple, span: float = 0.25, n: int = 15):
    """χ²(a, b) on an n×n grid around the fitted values of the two parameters.

    ``axis_a``/``axis_b`` are ``(site_index, param_name)``. Returns
    ``(A, B, Z, (a0, b0))`` with Z[ib, ia] = Σ residual² at (A[ia], B[ib])."""
    rec = Recipe.from_dict(recipe) if isinstance(recipe, dict) else recipe
    exp_ppm = np.asarray(exp_ppm, float)
    exp_amp = np.asarray(exp_amp, float)
    if window:
        hi, lo = max(window), min(window)
    else:
        hi, lo = float(exp_ppm.max()), float(exp_ppm.min())
    sel = (exp_ppm >= lo) & (exp_ppm <= hi)
    xw, yw = exp_ppm[sel], exp_amp[sel]

    pa = rec.sites[axis_a[0]].params[axis_a[1]]
    pb = rec.sites[axis_b[0]].params[axis_b[1]]
    a0, b0 = float(pa.value), float(pb.value)
    da = abs(a0) * span or span
    db = abs(b0) * span or span
    A = np.linspace(a0 - da, a0 + da, n)
    B = np.linspace(b0 - db, b0 + db, n)
    Z = np.zeros((n, n))
    for ib, bv in enumerate(B):
        pb.value = bv
        for ia, av in enumerate(A):
            pa.value = av
            x, y, _ = engine.simulate(rec, exp_ppm=exp_ppm)
            yi = np.interp(xw, x, y)
            Z[ib, ia] = float(np.sum((yi - yw) ** 2))
    pa.value, pb.value = a0, b0
    return A, B, Z, (a0, b0)
