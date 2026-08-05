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
_REF_FULL = re.compile(r"\bs(\d+)\.([A-Za-z_][A-Za-z0-9_]*)")


def param_refs(expr: str) -> list[tuple[int, str]]:
    """The (site index, parameter name) pairs a constraint expression references."""
    return [(int(m.group(1)), m.group(2)) for m in _REF_FULL.finditer(expr or "")]


def _find_cycle(graph: dict):
    """Return a node on a dependency cycle, or None (constraint graph = nodes are
    (site, param) with an expr; edges are the params they reference)."""
    color: dict = {}

    def dfs(node):
        color[node] = 1                          # grey (on the stack)
        for ref in graph.get(node, []):
            if ref not in graph:                 # references a free (leaf) param
                continue
            if color.get(ref) == 1:              # back-edge → cycle
                return ref
            if color.get(ref) is None and (r := dfs(ref)):
                return r
        color[node] = 2                          # black (done)
        return None

    for node in graph:
        if color.get(node) is None and (r := dfs(node)):
            return r
    return None


def sanitize_constraints(sites: list) -> list[str]:
    """Remove constraints that would make the fitter loop forever (recursion):
    direct self-references, references to a non-existent site, and cross-line
    cycles (A→B→A). Mutates ``sites`` in place; returns the dropped labels.

    Repairs a recipe loaded from a file that carries a broken link — the everyday
    cause of 'maximum recursion depth exceeded' during a fit."""
    n = len(sites)
    dropped: list[str] = []

    def _graph():
        g = {}
        for i, s in enumerate(sites):
            for pn, p in (s.get("params", {}) or {}).items():
                if isinstance(p, dict) and p.get("expr"):
                    g[(i, pn)] = param_refs(p["expr"])
        return g

    # 1) self-references and dangling (out-of-range) references
    for (i, pn), refs in _graph().items():
        if any((ri == i and rp == pn) or ri >= n or ri < 0 for ri, rp in refs):
            sites[i]["params"][pn]["expr"] = None
            dropped.append(f"s{i}.{pn}")

    # 2) cross-line cycles — break one link per cycle until none remain
    for _ in range(n * 4 + len(dropped) + 1):    # bounded
        cyc = _find_cycle(_graph())
        if cyc is None:
            break
        i, pn = cyc
        sites[i]["params"][pn]["expr"] = None
        dropped.append(f"s{i}.{pn}")
    return dropped


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

