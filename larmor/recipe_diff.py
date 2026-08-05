"""Compare two fits parameter-by-parameter — a quick sanity-check of the current
model against a previously-saved fit or a literature dmfit result.

Qt-free and testable; the desktop layer renders the rows in a table.
"""
from __future__ import annotations


def _val(site: dict | None, pname: str):
    if not site:
        return None
    p = site.get("params", {}).get(pname)
    if p is None:
        return None
    return float(p.get("value")) if isinstance(p, dict) else float(p)


def recipe_diff(current: dict, reference: dict) -> list[dict]:
    """Site-by-site, parameter-by-parameter differences between two recipes.

    Sites are matched by index. Each row is ``{site, label, param, current,
    reference, delta, delta_pct}`` (``None`` where a value is absent)."""
    cs = current.get("sites", []) or []
    rs = reference.get("sites", []) or []
    rows: list[dict] = []
    for i in range(max(len(cs), len(rs))):
        c = cs[i] if i < len(cs) else None
        r = rs[i] if i < len(rs) else None
        names: list[str] = []
        for site in (c, r):
            for pn in (site.get("params", {}) if site else {}):
                if pn != "gl" and pn not in names:
                    names.append(pn)
        label = ((c or r) or {}).get("label") or ((c or r) or {}).get("model", "")
        for pn in names:
            cv, rv = _val(c, pn), _val(r, pn)
            delta = (cv - rv) if (cv is not None and rv is not None) else None
            dpct = (100.0 * delta / rv) if (delta is not None and rv) else None
            rows.append({"site": i, "label": label, "param": pn,
                         "current": cv, "reference": rv,
                         "delta": delta, "delta_pct": dpct})
    return rows
