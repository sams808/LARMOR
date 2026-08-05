"""Residual diagnostics: is a fit's residual *white noise*, or is there
unmodelled structure the RMSD is hiding?

A small RMSD can still come from a systematically wrong model — the residual then
shows long runs of the same sign (a Wald–Wolfowitz *runs test* catches this) and
positive lag-1 autocorrelation. Flagging that tells the human "add/adjust a line"
even when the fit looks numerically fine. Qt-free and testable.
"""
from __future__ import annotations

import math

import numpy as np


def runs_test(residual) -> dict:
    """Wald–Wolfowitz runs test on the sign of the (median-centred) residual.

    Returns z (standard normal; strongly negative → too few runs → structure),
    the two-sided p-value, the run count and a ``structured`` flag.
    """
    r = np.asarray(residual, float)
    r = r[np.isfinite(r)]
    signs = np.sign(r - np.median(r))
    signs = signs[signs != 0]
    n = signs.size
    if n < 12:
        return {"z": 0.0, "p": 1.0, "runs": 0, "n": n, "structured": False}
    n_pos = int(np.sum(signs > 0))
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return {"z": 0.0, "p": 1.0, "runs": 1, "n": n, "structured": True}
    runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
    exp = 1.0 + 2.0 * n_pos * n_neg / n
    var = (2.0 * n_pos * n_neg * (2.0 * n_pos * n_neg - n)) / (n * n * (n - 1))
    if var <= 0:
        return {"z": 0.0, "p": 1.0, "runs": runs, "n": n, "structured": False}
    z = (runs - exp) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return {"z": z, "p": p, "runs": runs, "n": n, "structured": abs(z) > 3.0}


def lag1_autocorr(residual) -> float:
    """Lag-1 autocorrelation of the residual (≈1 = smoothly structured, 0 = white)."""
    r = np.asarray(residual, float)
    r = r[np.isfinite(r)]
    if r.size < 3:
        return 0.0
    r = r - r.mean()
    denom = float(np.sum(r * r))
    if denom <= 0:
        return 0.0
    return float(np.sum(r[1:] * r[:-1]) / denom)


def residual_structure(residual) -> dict:
    """Combine the runs test and lag-1 autocorrelation into a single verdict."""
    rt = runs_test(residual)
    ac = lag1_autocorr(residual)
    structured = bool(rt["structured"] or ac > 0.4)
    msg = ""
    if structured:
        why = []
        if rt["structured"]:
            why.append(f"runs test z={rt['z']:.1f}")
        if ac > 0.4:
            why.append(f"lag-1 autocorr {ac:.2f}")
        msg = ("residual is structured (" + ", ".join(why) +
               ") — the model may be missing a component or a lineshape is wrong")
    return {"structured": structured, "message": msg,
            "runs_z": rt["z"], "lag1": ac}
