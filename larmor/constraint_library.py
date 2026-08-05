"""Reusable constraint sets: capture the links / bounds / fixes on a model and
re-apply them to another one.

The physics you trust — "tie every δ_iso to site 0", "C_Q(field 2) = C_Q(field
1)", "fix the Gauss/Lorentz mix", "bound η to [0, 1]" — is worth naming once and
re-using across a project instead of re-typing per fit. A captured set is a plain
dict (per site index → per param → {expr, min, max, vary}); the desktop layer
stores named sets in QSettings. Qt-free and testable.
"""
from __future__ import annotations

import numpy as np


def capture(recipe) -> dict:
    """Extract the non-default constraints (links, finite bounds, fixes) of a
    recipe, keyed by site index and parameter name."""
    sites = []
    for s in recipe.sites:
        entry = {}
        for pn, p in s.params.items():
            c = {}
            if getattr(p, "expr", None):
                c["expr"] = p.expr
            if p.min is not None and np.isfinite(p.min):
                c["min"] = float(p.min)
            if p.max is not None and np.isfinite(p.max):
                c["max"] = float(p.max)
            if not p.vary:
                c["vary"] = False
            if c:
                entry[pn] = c
        sites.append(entry)
    return {"sites": sites}


def apply(recipe, cset: dict) -> list[str]:
    """Apply a captured set to ``recipe`` in place (matched by site index +
    parameter name; missing sites/params are skipped). Returns the params set."""
    applied = []
    for i, entry in enumerate(cset.get("sites", [])):
        if i >= len(recipe.sites):
            break
        for pn, c in entry.items():
            p = recipe.sites[i].params.get(pn)
            if p is None:
                continue
            if "expr" in c:
                p.expr = c["expr"]
            if "min" in c:
                p.min = c["min"]
            if "max" in c:
                p.max = c["max"]
            if "vary" in c:
                p.vary = bool(c["vary"])
            applied.append(f"s{i}.{pn}")
    return applied


def describe(cset: dict) -> str:
    """A short human summary of what a constraint set does."""
    n_links = n_bounds = n_fixed = 0
    for entry in cset.get("sites", []):
        for c in entry.values():
            n_links += bool(c.get("expr"))
            n_bounds += ("min" in c or "max" in c)
            n_fixed += (c.get("vary") is False)
    bits = []
    if n_links:
        bits.append(f"{n_links} link{'s' * (n_links != 1)}")
    if n_bounds:
        bits.append(f"{n_bounds} bounded")
    if n_fixed:
        bits.append(f"{n_fixed} fixed")
    return ", ".join(bits) or "empty"
