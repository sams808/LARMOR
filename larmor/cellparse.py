"""Parse what a user types into a fit-parameter cell.

Components are addressed by LETTERS (A, B, C, …), like dmfit — never by
number. A cell accepts:

    62.6266              a plain value
    A                    link: equal to component A's SAME parameter
    A+20                 link: A's value + 20 (in the parameter's own unit)
    A+20kHz              link: A's value + 20 kHz, converted to this unit
    A-1.5kHz             link: A's value − 1.5 kHz
    0.5B   /  0.5*B      link: half of component B's value (ratios)
    2A+10                link: 2·A + 10
    [0..100]             constrain to min 0, max 100 (value unchanged)
    62.6 [0..100]        set the value AND the bounds
    A+20 [50..80]        link AND bounds

Bounds accept  [lo..hi] , [lo:hi] , [lo,hi] ; a blank side means unbounded
(e.g. [0..] or [..100]). The parser is pure (no Qt) so it is unit-tested.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# a coefficient, a letter run, an optional signed offset, an optional unit
_LINK_RE = re.compile(
    r"^\s*(?P<coef>[+-]?\d*\.?\d+)?\s*\*?\s*"
    r"(?P<letter>[A-Za-z]+)"
    r"\s*(?P<off>[+-]\s*\d*\.?\d+)?\s*"
    r"(?P<unit>ppm|khz|hz|mhz)?\s*$", re.IGNORECASE)
_BOUNDS_RE = re.compile(r"\[\s*(?P<lo>[^\],:]*?)\s*(?:\.\.|:|,)\s*(?P<hi>[^\]]*?)\s*\]")
_NUM_RE = re.compile(r"^[+-]?\d*\.?\d+([eE][+-]?\d+)?$")
_NUM_UNIT_RE = re.compile(
    r"^(?P<num>[+-]?\d*\.?\d+([eE][+-]?\d+)?)\s*(?P<unit>ppm|khz|hz|mhz)$",
    re.IGNORECASE)


def _abs_in_native(value: float, unit: str, param_unit: str,
                   larmor_MHz: float) -> float:
    """An absolute value typed with a unit (e.g. '300Hz', '2ppm'), converted to
    the parameter's own unit — lets a width/position be entered in Hz OR ppm."""
    u, pu = unit.lower(), (param_unit or "").lower()
    if pu == "ppm" and larmor_MHz <= 0:
        raise ValueError("no Larmor frequency set — cannot convert Hz↔ppm")
    hz = {"hz": value, "khz": value * 1e3, "mhz": value * 1e6,
          "ppm": value * larmor_MHz}[u]
    if pu == "ppm":
        return hz / larmor_MHz
    if pu == "mhz":
        return hz / 1e6
    if pu == "hz":
        return hz
    raise ValueError(f"can't enter {unit} for a value in {param_unit or '—'}")


# --------------------------------------------------------------------------
def index_to_letter(i: int) -> str:
    """0→A, 25→Z, 26→AA (bijective base-26)."""
    i = int(i) + 1
    s = ""
    while i > 0:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def letter_to_index(s: str) -> int:
    n = 0
    for ch in s.strip().upper():
        if not ("A" <= ch <= "Z"):
            raise ValueError(f"not a component letter: {s!r}")
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _to_native(offset: float, unit: str | None, param_unit: str,
               larmor_MHz: float) -> float:
    """Convert a typed offset into the parameter's own unit."""
    if not unit or unit.lower() == param_unit.lower():
        return offset
    u = unit.lower()
    hz = {"hz": 1.0, "khz": 1e3, "mhz": 1e6}.get(u)
    if hz is None:
        return offset
    hz *= offset
    pu = param_unit.lower()
    if pu == "ppm":
        if larmor_MHz <= 0:
            raise ValueError("no Larmor frequency set — cannot convert Hz→ppm "
                             "(Process ▸ Experiment parameters…)")
        return hz / larmor_MHz            # Hz / (MHz) = ppm
    if pu == "mhz":
        return hz / 1e6                   # Hz → MHz  (e.g. Cq)
    if pu == "khz":
        return hz / 1e3
    raise ValueError(f"can't add {unit} to a value in {param_unit}")


def _fmt(x: float) -> str:
    return f"{x:g}"


@dataclass
class CellResult:
    """What to apply to the Param after parsing a cell. Fields left as the
    sentinel `_UNSET` are not touched."""

    value: float | None = None
    set_value: bool = False
    expr: str | None = None
    set_expr: bool = False
    min: float | None = None
    set_min: bool = False
    max: float | None = None
    set_max: bool = False
    error: str | None = None


def parse_cell(text: str, *, param_name: str, param_unit: str,
               this_index: int, n_sites: int, larmor_MHz: float) -> CellResult:
    """Parse a cell entry for the parameter `param_name` of site `this_index`.

    Links always reference the SAME parameter of the target component.
    """
    text = (text or "").strip()
    res = CellResult()

    # 1) pull off an optional bounds suffix
    mb = _BOUNDS_RE.search(text)
    if mb:
        lo, hi = mb.group("lo").strip(), mb.group("hi").strip()
        try:
            res.set_min = True
            res.min = None if lo in ("", "-", "-inf", "−inf") else float(lo)
            res.set_max = True
            res.max = None if hi in ("", "inf", "+inf") else float(hi)
        except ValueError:
            return CellResult(error=f"bad bounds: {mb.group(0)}")
        if res.min is not None and res.max is not None and res.min >= res.max:
            return CellResult(error="min must be < max")
        text = (text[:mb.start()] + text[mb.end():]).strip()

    if not text:                          # bounds-only edit
        return res

    # 2) a plain number un-links and sets the value
    if _NUM_RE.match(text):
        res.set_value = True
        res.value = float(text)
        res.set_expr = True               # clear any existing link
        res.expr = None
        return res

    # 2b) a plain number WITH a unit (e.g. 300Hz, 1.5kHz, 2ppm) -> native unit
    mu = _NUM_UNIT_RE.match(text)
    if mu:
        try:
            res.value = _abs_in_native(float(mu.group("num")), mu.group("unit"),
                                       param_unit, larmor_MHz)
        except (ValueError, KeyError) as exc:
            return CellResult(error=str(exc))
        res.set_value = True
        res.set_expr = True
        res.expr = None
        return res

    # 3) a letter link
    m = _LINK_RE.match(text)
    if not m:
        return CellResult(error=f"cannot read {text!r} — type a number, or a "
                                "link like A+20, A+20kHz, or 0.5B")
    try:
        idx = letter_to_index(m.group("letter"))
    except ValueError as exc:
        return CellResult(error=str(exc))
    if idx >= n_sites:
        return CellResult(error=f"component {m.group('letter').upper()} does "
                                f"not exist ({n_sites} lines)")
    if idx == this_index:
        return CellResult(error="a parameter cannot link to itself")

    coef = float(m.group("coef")) if m.group("coef") else 1.0
    off_native = 0.0
    if m.group("off"):
        off = float(m.group("off").replace(" ", ""))
        try:
            off_native = _to_native(off, m.group("unit"), param_unit, larmor_MHz)
        except ValueError as exc:
            return CellResult(error=str(exc))

    ref = f"s{idx}.{param_name}"
    expr = ref if coef == 1.0 else f"{_fmt(coef)}*{ref}"
    if abs(off_native) > 1e-12:
        expr += f" + {_fmt(off_native)}" if off_native >= 0 else f" - {_fmt(-off_native)}"
    res.set_expr = True
    res.expr = expr
    return res


# --------------------------------------------------------------------------
_EXPR_LINK = re.compile(
    r"^\s*(?:(?P<coef>[-+]?\d*\.?\d+)\s*\*\s*)?"
    r"s(?P<idx>\d+)\.(?P<name>[A-Za-z_][A-Za-z0-9_]*)"
    r"\s*(?P<off>[+-]\s*\d*\.?\d+)?\s*$")


def format_link(expr: str, param_name: str) -> str:
    """Turn an internal expr back into the friendly letter form for display.

    `s0.isotropic_chemical_shift_ppm + 20` → 'A+20';
    `0.5*s1.amplitude` → '0.5B'. Unrecognized exprs are shown as-is.
    """
    m = _EXPR_LINK.match(expr or "")
    if not m or m.group("name") != param_name:
        return expr
    letter = index_to_letter(int(m.group("idx")))
    coef = m.group("coef")
    out = (letter if not coef or float(coef) == 1.0
           else f"{_fmt(float(coef))}{letter}")
    if m.group("off"):
        off = float(m.group("off").replace(" ", ""))
        out += f"+{_fmt(off)}" if off >= 0 else f"-{_fmt(-off)}"
    return out
