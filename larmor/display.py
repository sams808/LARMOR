"""Display-only transforms of the 1D spectrum view -- Qt-free.

Three things the desktop draws differently from what it stores, all pure
array arithmetic so the view (larmor.desktop.plot) and the overlay cockpit
(larmor.desktop.mw_overlays) apply them from one place and a test can check
the numbers without a window:

* **Y normalisation of the whole plotting area** (View > Y axis). One factor,
  computed from the ACTIVE spectrum -- its maximum, its trapezoid area over
  the axis, or the area of a ppm region -- multiplies everything drawn for
  that spectrum: experiment, model, components, residual, animation frames
  and the paddles (whose drags map back through the inverse). The recipe,
  the fit and every data export stay in raw units.
* **The per-overlay display transform** of the Datasets dock: a scale
  factor, an x shift and a y offset (a fraction of the active spectrum's
  displayed span) plus the global stack offset, applied to a compared
  spectrum's arrays at draw time only. Under a normalisation mode each
  overlay is normalised by its OWN maximum / area first, so shapes compare
  across spectra of different intensity; the scale then acts on top. The
  scale is therefore a DISPLAY-space number and 'match height' fills it
  through :func:`match_scale_display`, not with the raw peak ratio.
* **The Full view**: the data's x extent with a small margin and a y range
  with room for the residual strip below zero, computed from the arrays
  (never from pyqtgraph's auto-range, which the model items would widen),
  plus the rule that sends a view back to it when it no longer shows the
  data at all or has zoomed out far beyond it.
"""
from __future__ import annotations

import numpy as np

#: the Y-axis display modes, in menu order
Y_MODES = ("raw", "max", "area", "region")

#: the per-overlay display transform and its neutral values -- the keys an
#: overlay dict, a workspace snapshot and a project bundle carry (a file
#: written before they existed reads back with these)
OVERLAY_DEFAULTS = {"scale": 1.0, "shift": 0.0, "yoff": 0.0}

#: plain names of the modes (the desktop adds accelerators)
Y_MODE_LABELS = {
    "raw": "Raw intensity",
    "max": "Normalise to maximum",
    "area": "Normalise to area",
    "region": "Normalise to area of a region…",
}

#: fraction of the data span left free beyond each end of the x axis by the
#: Full view (a curve that ends at the last point is not cut visibly)
FULL_MARGIN_FRAC = 0.02
#: the residual strip sits at -0.10 x max; the Full view leaves this much of
#: the maximum below zero for it
RESID_ROOM_FRAC = 0.12
#: padding of the Full view's y range, a fraction of its height
Y_PAD_FRAC = 0.08
#: a requested view wider than this many data spans goes back to the Full view
SNAP_WIDTH_FACTOR = 3.0
#: the mouse may pan this many data spans beyond the data in x ...
LIMIT_X_SPANS = 1.0
#: ... and this many in y (stacked overlays and a tall first guess need room)
LIMIT_Y_SPANS = 2.0

_trapezoid = getattr(np, "trapezoid", None) or getattr(np, "trapz")


def _finite_pair(x, y):
    """(x, y) as float arrays trimmed to a common length and to finite points."""
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    n = min(x.size, y.size)
    x, y = x[:n], y[:n]
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok]


# ------------------------------------------------------------- Y normalisation
def area(x, y, region=None) -> float:
    """Trapezoid area of y over x, positive for a descending ppm axis too
    (integrated over |dx|), restricted to ``region = (hi, lo)`` in axis
    units when given. 0.0 when fewer than two finite points fall in the
    region."""
    x, y = _finite_pair(x, y)
    if region is not None:
        hi, lo = max(region), min(region)
        keep = (x >= lo) & (x <= hi)
        x, y = x[keep], y[keep]
    if x.size < 2:
        return 0.0
    order = np.argsort(x)
    return float(abs(_trapezoid(y[order], x[order])))


def y_factor(mode: str, x, y, region=None) -> float:
    """The multiplier that takes raw y to display y under ``mode``: 1/max for
    "max", 1/area for "area", 1/area(region) for "region" and 1.0 for "raw".
    Also 1.0 whenever the reference is not positive and finite -- no data,
    an all-zero or inverted trace, a region with no points -- so a
    degenerate spectrum is drawn raw rather than blown up."""
    if mode not in Y_MODES or mode == "raw":
        return 1.0
    if y is None:
        return 1.0
    ya = np.asarray(y, float).ravel()
    if not ya.size:
        return 1.0
    if mode == "max":
        ref = float(np.nanmax(ya)) if np.isfinite(ya).any() else 0.0
    elif mode == "area":
        ref = area(x, ya) if x is not None else 0.0
    else:                                   # region
        if region is None or x is None:
            return 1.0
        ref = area(x, ya, region)
    if not np.isfinite(ref) or ref <= 0.0:
        return 1.0
    return 1.0 / ref


def _fmt_ppm(v: float) -> str:
    return f"{float(v):g}".replace("-", "−")


def format_region(region) -> str:
    """'100…−20 ppm' (hi first, a real minus sign); '' for None."""
    if region is None:
        return ""
    hi, lo = max(region), min(region)
    return f"{_fmt_ppm(hi)}…{_fmt_ppm(lo)} ppm"


def parse_region(text) -> tuple[float, float] | None:
    """The inverse of the persisted 'hi,lo' form; None when it is not two
    finite, distinct numbers."""
    if not text:
        return None
    try:
        a, b = (float(v) for v in str(text).split(","))
    except (TypeError, ValueError):
        return None
    if not (np.isfinite(a) and np.isfinite(b)) or a == b:
        return None
    return (max(a, b), min(a, b))


def y_axis_label(mode: str, region=None) -> str:
    """The left-axis label for a mode: 'intensity' raw, 'intensity
    (normalised to max)', '(normalised to area)', '(normalised to area
    100…−20 ppm)'."""
    if mode == "max":
        return "intensity (normalised to max)"
    if mode == "area":
        return "intensity (normalised to area)"
    if mode == "region" and region is not None:
        return f"intensity (normalised to area {format_region(region)})"
    return "intensity"


# ------------------------------------------------------------- overlays
def overlay_display(ppm, amp, *, scale: float = 1.0, shift: float = 0.0,
                    yoff: float = 0.0, stack: float = 0.0, span: float = 1.0,
                    mode: str = "raw", region=None):
    """(x, y) as DRAWN for one compared spectrum; the stored arrays are never
    modified (new arrays come back even for the defaults).

    ``x = ppm + shift`` and ``y = normalise(amp) * scale + (yoff + stack) *
    span``, where ``normalise`` is :func:`y_factor` of ``mode`` computed from
    this overlay's OWN (shifted) trace -- under "max" every overlay reaches
    ``scale`` (1.0 by default) -- and ``span`` is the active spectrum's
    displayed span, so ``yoff`` and the ``stack`` offset are fractions of
    what is on screen."""
    x = np.asarray(ppm, float) + float(shift)
    y = np.asarray(amp, float)
    f = y_factor(mode, x, y, region)
    return x, y * (f * float(scale)) + (float(yoff) + float(stack)) * float(span)


def _peak(y) -> float:
    """The trace's maximum, 0.0 when it has no finite point."""
    a = np.asarray(y, float).ravel()
    if not a.size or not np.isfinite(a).any():
        return 0.0
    return float(np.nanmax(a))


def match_scale(active_amp, amp) -> float:
    """The scale that brings an overlay's RAW maximum to the active
    spectrum's raw maximum; 1.0 when either maximum is not positive."""
    pa, pb = _peak(active_amp), _peak(amp)
    if pa <= 0.0 or pb <= 0.0:
        return 1.0
    return pa / pb


def match_scale_display(active_amp, ppm, amp, *, active_factor: float = 1.0,
                        mode: str = "raw", region=None,
                        shift: float = 0.0) -> float:
    """The 'match height' factor as :func:`overlay_display` consumes it: the
    × that makes this overlay's DRAWN peak equal the active spectrum's drawn
    peak under ``mode``.

    ``active_factor`` is what the active trace is drawn with (the view's
    y_scale); the overlay's own factor is derived here exactly as
    ``overlay_display`` does -- from its SHIFTED axis, so a lined-up
    reference matches over the same region. The result is therefore
    ``(max_active · f_active) / (max_overlay · f_overlay)``: the raw peak
    ratio under "raw" (where both factors are 1.0), exactly 1.0 under "max"
    (both traces already reach 1.0) and the displayed-peak ratio under
    "area" / "region". 1.0 whenever a reference is not positive and finite.
    """
    x = np.asarray(ppm, float) + float(shift)
    fb = y_factor(mode, x, amp, region)
    fa = float(active_factor)
    if not np.isfinite(fa) or fa <= 0.0:
        fa = 1.0
    pa, pb = _peak(active_amp) * fa, _peak(amp) * fb
    if pa <= 0.0 or pb <= 0.0 or not np.isfinite(pa / pb):
        return 1.0
    return pa / pb


def overlay_badge(scale: float = 1.0, shift: float = 0.0,
                  yoff: float = 0.0) -> str:
    """The Datasets label's suffix for a transformed overlay -- '×2.5 ·
    +1.2 ppm · ↑0.3' -- listing only the non-default parts; '' when the
    overlay is drawn as stored."""
    bits = []
    if abs(float(scale) - 1.0) > 1e-12:
        bits.append(f"×{float(scale):.4g}")
    if abs(float(shift)) > 1e-12:
        bits.append(f"{float(shift):+.4g} ppm")
    if abs(float(yoff)) > 1e-12:
        arrow = "↑" if yoff > 0 else "↓"
        bits.append(f"{arrow}{abs(float(yoff)):.3g}")
    return " · ".join(bits)


# ------------------------------------------------------------- the Full view
def full_extents(x, y, margin_frac: float = FULL_MARGIN_FRAC,
                 resid_frac: float = RESID_ROOM_FRAC,
                 pad_frac: float = Y_PAD_FRAC):
    """``((xlo, xhi), (ylo, yhi))`` of the Full view from the data arrays
    alone: the x extent widened by ``margin_frac`` of its span on each side,
    the y extent with ``resid_frac`` of the maximum kept free below zero for
    the residual strip and ``pad_frac`` of the height as padding. None when
    there is no finite point."""
    x, y = _finite_pair(x, y)
    if not x.size:
        return None
    xlo, xhi = float(x.min()), float(x.max())
    span = xhi - xlo
    m = margin_frac * span if span > 0 else max(abs(xlo) * 0.01, 1.0)
    ylo, yhi = float(y.min()), float(y.max())
    ylo = min(ylo, -resid_frac * yhi)          # room for the residual strip
    pad = pad_frac * ((yhi - ylo) or 1.0)
    return (xlo - m, xhi + m), (ylo - pad, yhi + pad)


def snap_back(requested, data, width_factor: float = SNAP_WIDTH_FACTOR) -> bool:
    """Should a REQUESTED x range ``(a, b)`` be replaced by the Full view?
    True when it does not intersect the data extent ``(lo, hi)`` at all --
    the user panned or zoomed completely outside the spectrum -- or is wider
    than ``width_factor`` data spans. Either pair in either order; a zoom
    inside the data is never touched."""
    try:
        lo, hi = float(min(data)), float(max(data))
        a, b = float(min(requested)), float(max(requested))
    except (TypeError, ValueError):
        return False
    span = hi - lo
    if not np.isfinite([lo, hi, a, b]).all() or span <= 0.0:
        return False
    if b <= lo or a >= hi:
        return True
    return (b - a) > width_factor * span


def view_limits(x_bounds, y_bounds, x_spans: float = LIMIT_X_SPANS,
                y_spans: float = LIMIT_Y_SPANS) -> dict:
    """``ViewBox.setLimits`` keywords for a generous envelope around what is
    drawn: ``x_spans`` data spans beyond each end in x, ``y_spans`` in y --
    the mouse cannot pan into nowhere, yet every zoom inside stays free."""
    xlo, xhi = float(min(x_bounds)), float(max(x_bounds))
    ylo, yhi = float(min(y_bounds)), float(max(y_bounds))
    sx = (xhi - xlo) or max(abs(xlo) * 0.01, 1.0)
    sy = (yhi - ylo) or max(abs(yhi) * 0.01, 1.0)
    return {"xMin": xlo - x_spans * sx, "xMax": xhi + x_spans * sx,
            "yMin": ylo - y_spans * sy, "yMax": yhi + y_spans * sy}
