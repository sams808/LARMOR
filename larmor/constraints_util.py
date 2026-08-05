"""Utilities for keeping parameter constraints (Param.expr) consistent when the
set of sites changes.

A constraint is stored as an expression over site references, e.g.
``"s3.isotropic_chemical_shift_ppm + 5.3"`` (site 3's δiso plus 5.3). When a site
is deleted the remaining sites are renumbered, so every reference must be remapped
— otherwise a link like ``E = D + 5.3`` can silently become ``E = E + 5.3`` (a
self-reference), which recurses forever when the fitter evaluates it. This module
does that remapping (and drops links that can no longer be satisfied). Qt-free and
testable.
"""
from __future__ import annotations

import re

_SITE_REF = re.compile(r"\bs(\d+)\.")


def site_refs(expr: str) -> set[int]:
    """The set of site indices a constraint expression references."""
    return {int(m.group(1)) for m in _SITE_REF.finditer(expr or "")}


def references_self(expr: str, site_idx: int, param: str) -> bool:
    """True if a constraint references its OWN (site, parameter) — a direct
    self-reference that would recurse forever at fit time."""
    return bool(re.search(rf"\bs{site_idx}\.{re.escape(param)}\b", expr or ""))


def remap_exprs_after_delete(sites: list, deleted_idx: int) -> list[str]:
    """Fix constraints in ``sites`` (mutated in place) after site ``deleted_idx``
    was removed. Drops any expr that referenced the deleted site or that would
    become a self-reference; shifts ``s<k>`` with ``k>deleted_idx`` down by one.
    Returns the ``"s<i>.<param>"`` labels of the dropped constraints."""
    dropped: list[str] = []
    for new_i, site in enumerate(sites):
        for pname, p in (site.get("params", {}) or {}).items():
            expr = p.get("expr") if isinstance(p, dict) else None
            if not expr:
                continue
            refs = site_refs(expr)
            if deleted_idx in refs:                 # target is gone
                p["expr"] = None
                dropped.append(f"s{new_i}.{pname}")
                continue

            def _shift(m):
                k = int(m.group(1))
                return f"s{k - 1 if k > deleted_idx else k}."

            new_expr = _SITE_REF.sub(_shift, expr)
            if new_i in site_refs(new_expr):        # became a self-reference
                p["expr"] = None
                dropped.append(f"s{new_i}.{pname}")
            else:
                p["expr"] = new_expr
    return dropped


def remap_exprs_after_move(sites: list, old_to_new: dict) -> None:
    """Remap constraints after the sites were reordered. ``sites`` is already in
    its NEW order but each ``s<k>`` still refers to an OLD index; ``old_to_new``
    maps old index → new index. A reorder is a bijection, so nothing is dropped
    (references just follow the sites they point at)."""
    def _map(m):
        return f"s{old_to_new.get(int(m.group(1)), int(m.group(1)))}."

    for site in sites:
        for p in (site.get("params", {}) or {}).values():
            if isinstance(p, dict) and p.get("expr"):
                p["expr"] = _SITE_REF.sub(_map, p["expr"])

