"""Physical sanity checks for a fitted model.

A least-squares fit will happily return numbers that are *unphysical* — an
asymmetry η outside [0, 1], a negative line width, a negative-amplitude (negative
population) line, or a site whose centre sits outside the fit window (so its
parameters are unconstrained). This module flags those so the human can judge
whether the fit is physically meaningful — LARMOR's paramount concern. It only
reports; it never silently "corrects" a value (that stays a human decision).

Qt-free and testable.
"""
from __future__ import annotations

import numpy as np

#: parameter names that are an **asymmetry** η and must lie in [0, 1]
_ETA = {"eta", "eta_q", "etaq", "eta_cs", "etacs"}
#: names that are a **width / spread** and must be strictly positive
_WIDTH_EXACT = {"sigma_cq_mhz", "sigma_zeta_ppm", "deta"}


def _is_width(name: str) -> bool:
    low = name.lower()
    return ("fwhm" in low or "width" in low or low in _WIDTH_EXACT)


def check_recipe(recipe, window=None, *, tol: float = 1e-6) -> list[dict]:
    """Return a list of ``{site, label, param, message}`` physical warnings.

    ``window`` (hi, lo) ppm, if given, flags any site centred outside it.
    """
    warns: list[dict] = []
    lo = hi = None
    if window and len(window) == 2:
        lo, hi = min(window), max(window)

    def add(i, label, param, message):
        warns.append({"site": i, "label": label, "param": param, "message": message})

    for i, site in enumerate(recipe.sites):
        label = site.label or site.model
        for pn, p in site.params.items():
            try:
                v = float(p.value)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(v):
                add(i, label, pn, f"{pn} is not finite")
                continue
            low = pn.lower()
            if pn == "amplitude":
                if v < -tol:
                    add(i, label, pn, f"negative amplitude ({v:.3g}) — a "
                                      "negative-population line")
            elif low in _ETA:
                if not (-tol <= v <= 1 + tol):
                    add(i, label, pn, f"η = {v:.3g} outside [0, 1]")
            elif pn == "gl":
                if not (-tol <= v <= 1 + tol):
                    add(i, label, pn, f"Gauss/Lorentz mix = {v:.3g} outside [0, 1]")
            elif _is_width(low):
                if v <= tol:
                    add(i, label, pn, f"{pn} = {v:.3g} ≤ 0 — a non-physical width")
        centre = site.params.get("isotropic_chemical_shift_ppm")
        if centre is not None and lo is not None:
            c = float(centre.value)
            if not (lo <= c <= hi):
                add(i, label, "isotropic_chemical_shift_ppm",
                    f"δiso = {c:.3g} ppm is outside the fit window "
                    f"[{lo:.3g}, {hi:.3g}] — unconstrained")
    return warns


def summarize(warns: list[dict]) -> str:
    """A one-line, human-readable summary of the warnings (empty if none)."""
    if not warns:
        return ""
    parts = [f"{w['label']}: {w['message']}" for w in warns]
    return "⚠ physical check — " + "; ".join(parts)
