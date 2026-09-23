"""Derived parameter status: † held fixed, ‡ finished at a bound, § linked.

A results table must never present a held value as a fitted one. The status
of every parameter is DERIVED from what the recipe already stores -- ``vary``,
``expr``, and the value against its effective bound -- never persisted: it is
retroactive on every saved recipe, self-clearing on edit, and needs no recipe
schema change (``Param(**p)`` raises on unknown keys in released versions).

The at-bound test is the fit's own (``fit._at_bounds``): a FREE parameter
within ``AT_BOUND_REL_TOL`` · max(1, |value|) of a finite bound, where the
bound is the recipe's min/max or, when the recipe omits it, the model's
physical bound from the registry (``fit._make_params``' fallback). fit.py
imports ``bound_side`` / ``effective_bounds`` back from here, so the fit and
every table agree by construction. Verified read-only on the 32 PBi Final2
recipes: the derived at-bound set equals the fit's own note in 32/32.

Qt-free and lmfit-free: imports only ``larmor.models`` (the registry loads in
~0.2 s without mrsimulator or lmfit) and ``larmor.recipe``, so methods.py and
io/export.py can format a saved recipe without loading the engine.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from larmor import models as model_registry
from larmor.recipe import Recipe

__all__ = [
    "AT_BOUND_REL_TOL", "AT_BOUND_NOTE_PREFIX", "KINDS", "MARK", "LATEX",
    "CSV_COLUMNS", "ParamStatus", "effective_bounds", "bound_side",
    "param_status", "recipe_statuses", "site_constraints", "csv_fields",
    "footnote", "summary",
]

#: relative tolerance of the at-bound test -- the fit's own 1e-3 rule
AT_BOUND_REL_TOL = 1e-3
#: the recipe note fit() writes (and prunes before rewriting) when
#: parameters finished at a bound; the text after the prefix is unchanged
AT_BOUND_NOTE_PREFIX = "parameters finished at a bound"
KINDS = ("free", "fixed", "linked", "at_bound")
#: the glyphs on screen and in Markdown ('' for free)
MARK = {"free": "", "fixed": "†", "at_bound": "‡", "linked": "§"}
#: the same glyphs as LaTeX superscripts ('' for free)
LATEX = {"free": "", "fixed": r"$^{\dagger}$", "at_bound": r"$^{\ddagger}$",
         "linked": r"$^{\S}$"}
#: the five status columns appended to every long CSV (batch dialog, CLI)
CSV_COLUMNS = ("vary", "min", "max", "expr", "at_bound")

_ORDER = ("fixed", "at_bound", "linked")
_SIDE_WORD = {"min": "lower", "max": "upper"}


@dataclass(frozen=True)
class ParamStatus:
    """The derived status of one parameter."""

    kind: str                  # free | fixed | linked | at_bound
    side: str = ""             # 'min' | 'max' when kind == 'at_bound'
    bound: float | None = None  # the bound hit, when kind == 'at_bound'
    expr: str = ""             # the expression, when kind == 'linked'

    @property
    def marker(self) -> str:
        return MARK.get(self.kind, "")

    @property
    def latex(self) -> str:
        return LATEX.get(self.kind, "")

    @property
    def csv_flag(self) -> str:
        """'fixed' | 'linked' | 'at_min' | 'at_max' | '' -- words, never
        glyphs, so spreadsheet filters work."""
        if self.kind == "at_bound":
            return f"at_{self.side}"
        return self.kind if self.kind in ("fixed", "linked") else ""

    @property
    def word(self) -> str:
        """'held fixed' | 'at its lower bound 10.37' | 'linked: 0.19 * s0.amplitude' | ''"""
        if self.kind == "fixed":
            return "held fixed"
        if self.kind == "at_bound":
            return f"at its {_SIDE_WORD.get(self.side, self.side)} bound {_g(self.bound)}"
        if self.kind == "linked":
            return f"linked: {self.expr}"
        return ""


# ------------------------------------------------------------------ inputs
def _field(p, name: str, default=None):
    """Read ``name`` from a recipe.Param, a ``{value, stderr, vary, min, max,
    expr}`` dict (the desktop's shape) or anything else (-> default)."""
    if isinstance(p, dict):
        return p.get(name, default)
    return getattr(p, name, default)


def _is_param_like(p) -> bool:
    return isinstance(p, dict) or hasattr(p, "value")


def _finite(x) -> bool:
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _g(v) -> str:
    return "" if v is None else f"{float(v):.4g}"


def _num(v) -> str:
    return "" if v is None else f"{float(v):.8g}"


# -------------------------------------------------------------------- rule
def effective_bounds(model: str | None, pname: str, p
                     ) -> tuple[float | None, float | None]:
    """The bounds the fit enforces: the recipe's min/max, else the model's
    physical bounds from the registry (``ParamDef.min/max``); any registry
    failure (unknown model, odd parameter) leaves the recipe bounds alone --
    exactly ``fit._make_params``' fallback."""
    pmin, pmax = _field(p, "min"), _field(p, "max")
    if pmin is None or pmax is None:
        try:
            pdefs = {pd.name: pd for pd in model_registry.get(model).params}
        except Exception:
            pdefs = {}
        pd = pdefs.get(pname)
        if pmin is None:
            pmin = pd.min if pd else None
        if pmax is None:
            pmax = pd.max if pd else None
    return pmin, pmax


def bound_side(value, lo, hi, tol: float = AT_BOUND_REL_TOL) -> str | None:
    """'min' when ``value - lo < tol * max(1, |value|)``, 'max' likewise for
    ``hi``, else None. None and ±inf bounds never match (lmfit's unbounded
    parameters carry ±inf). The fit's own test, moved here."""
    if not _finite(value):
        return None
    v = float(value)
    span = max(1.0, abs(v))
    if _finite(lo) and (v - float(lo)) < tol * span:
        return "min"
    if _finite(hi) and (float(hi) - v) < tol * span:
        return "max"
    return None


def param_status(model: str | None, pname: str, p) -> ParamStatus:
    """Derive the status of one parameter from a recipe.Param, a
    ``{value, stderr, vary, min, max, expr}`` dict, or a bare number (which
    carries no constraint information and reads as free). Precedence:
    linked (expr) > fixed (vary False) > at a bound (free and within the
    tolerance of an effective bound) > free."""
    if p is None or not _is_param_like(p):
        return ParamStatus("free")
    expr = _field(p, "expr") or ""
    if expr:
        return ParamStatus("linked", expr=str(expr))
    vary = _field(p, "vary", True)
    if vary is not None and not vary:
        return ParamStatus("fixed")
    lo, hi = effective_bounds(model, pname, p)
    side = bound_side(_field(p, "value"), lo, hi)
    if side is None:
        return ParamStatus("free")
    return ParamStatus("at_bound", side=side,
                       bound=float(lo if side == "min" else hi))


# ------------------------------------------------------------- per recipe
def _sites(recipe):
    """Yield ``(model, label, params)`` from a Recipe or a recipe dict."""
    if isinstance(recipe, Recipe):
        for s in recipe.sites:
            yield s.model, s.label or "", s.params
        return
    for s in (recipe or {}).get("sites", []) or []:
        if isinstance(s, dict):
            yield s.get("model"), s.get("label") or "", s.get("params") or {}
        else:
            yield getattr(s, "model", None), getattr(s, "label", "") or "", \
                getattr(s, "params", {}) or {}


def recipe_statuses(recipe) -> dict[tuple[int, str], ParamStatus]:
    """``{(site_index, param_name): ParamStatus}`` for a Recipe or dict."""
    out: dict[tuple[int, str], ParamStatus] = {}
    for i, (model, _label, params) in enumerate(_sites(recipe)):
        for pname, p in params.items():
            out[(i, pname)] = param_status(model, pname, p)
    return out


def site_constraints(recipe, i: int) -> dict:
    """The constraints of site ``i`` as ``{'fixed': [pname], 'linked':
    [(pname, expr)], 'at_bound': [(pname, side, bound)]}`` -- the report's
    one-line-per-sample SI paragraph."""
    out: dict = {"fixed": [], "linked": [], "at_bound": []}
    for k, (model, _label, params) in enumerate(_sites(recipe)):
        if k != i:
            continue
        for pname, p in params.items():
            st = param_status(model, pname, p)
            if st.kind == "fixed":
                out["fixed"].append(pname)
            elif st.kind == "linked":
                out["linked"].append((pname, st.expr))
            elif st.kind == "at_bound":
                out["at_bound"].append((pname, st.side, st.bound))
    return out


# ------------------------------------------------------------- vocabularies
def csv_fields(p, st=None) -> list[str]:
    """The five CSV_COLUMNS cells for one parameter: ``['True'|'False'|'',
    min, max, expr, side]`` with '' for None and numbers as %.8g. ``p`` is a
    recipe.Param, a param dict, or a table row carrying ``vary/min/max/expr/
    at_bound``; ``st`` a ParamStatus, a side string, or None (then the side
    is read from ``p['at_bound']``)."""
    if isinstance(st, ParamStatus):
        side = st.side
    elif st is None:
        side = _field(p, "at_bound", "") or ""
    else:
        side = str(st)
    vary = _field(p, "vary")
    return ["" if vary is None else str(bool(vary)),
            _num(_field(p, "min")), _num(_field(p, "max")),
            str(_field(p, "expr") or ""), str(side)]


def footnote(entries: list[tuple[str, ParamStatus]], style: str = "latex") -> str:
    """One footnote line for a table: only the kinds present, in the order
    fixed / at a bound / linked; each at-bound entry with its bound, each
    linked entry with its expression; '' when every status is free.
    ``style`` 'latex' uses $\\dagger$ $\\ddagger$ $\\S$, 'text' the glyphs."""
    glyph = ({"fixed": r"$\dagger$", "at_bound": r"$\ddagger$", "linked": r"$\S$"}
             if style == "latex" else MARK)
    parts = []
    for kind in _ORDER:
        hits = [(label, st) for label, st in entries if st.kind == kind]
        if not hits:
            continue
        if kind == "fixed":
            parts.append(f"{glyph[kind]} held fixed.")
        elif kind == "at_bound":
            items = "; ".join(
                f"{label} ({_SIDE_WORD.get(st.side, st.side)} bound {_g(st.bound)})"
                for label, st in hits)
            parts.append(f"{glyph[kind]} finished at a bound: {items}.")
        else:
            items = "; ".join(f"{label} = {st.expr}" for label, st in hits)
            parts.append(f"{glyph[kind]} linked: {items}.")
    return " ".join(parts)


def summary(recipe) -> str:
    """'1 fixed · 2 at a bound · 1 linked' (kinds present only) or ''."""
    counts = {"fixed": 0, "at_bound": 0, "linked": 0}
    for st in recipe_statuses(recipe).values():
        if st.kind in counts:
            counts[st.kind] += 1
    words = {"fixed": "fixed", "at_bound": "at a bound", "linked": "linked"}
    return " · ".join(f"{counts[k]} {words[k]}" for k in _ORDER if counts[k])
