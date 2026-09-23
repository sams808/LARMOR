"""Drag-to-phase gesture arithmetic (TopSpin-style), Qt-free.

The SpectrumView turns a canvas drag into pixel deltas; everything that turns
pixels into degrees lives here, so the sensitivities, the dominant-axis lock
and the pivot re-expression are unit-tested without a window and retuned in
one place.

Conventions: drag right = +p0, drag up = +p1 (screen-relative, independent of
the inverted ppm axis); p0 wraps into [-180, 180] exactly like the Processing
panel's quick-step buttons; p1 is left to the panel's spin-box clamp.
"""
from __future__ import annotations

from dataclasses import dataclass

#: zero-order degrees per horizontal logical pixel
P0_DEG_PER_PX = 0.25
#: first-order degrees per vertical logical pixel (screen-up positive)
P1_DEG_PER_PX = 1.0
#: Shift multiplier (fine control)
FINE = 0.1
#: travel (px, either axis) before the dominant axis is locked
LOCK_PX = 4.0


def wrap_p0(deg: float) -> float:
    """Wrap a zero-order phase into [-180, 180].

    The ProcessingPanel rule: 180 and -180 stay where they are, 190 -> -170,
    -190 -> 170, 540 -> 180.
    """
    v = float(deg)
    while v > 180.0:
        v -= 360.0
    while v < -180.0:
        v += 360.0
    return v


def compensate_p0(p0: float, p1: float, old_frac: float,
                  new_frac: float) -> float:
    """The zero order that leaves the spectrum unchanged when the p1 pivot
    moves from ``old_frac`` to ``new_frac`` (0..1 fractions along the axis).

    ``processing.op_phase`` applies ``p0 + p1 * (idx - pivot_frac)``, so
    ``op_phase(y, p0, p1, old) == op_phase(y, p0 + p1*(new - old), p1, new)``
    exactly; the result is wrapped like every other p0.
    """
    return wrap_p0(float(p0) + float(p1) * (float(new_frac) - float(old_frac)))


@dataclass
class DragGesture:
    """Cumulative dp0 / dp1 (degrees) of one button-down .. release gesture.

    Travel is withheld until ``LOCK_PX`` in either direction; the gesture then
    locks to the dominant axis (horizontal -> p0, vertical -> p1) for the rest
    of the drag and the withheld pixels are applied, so nothing is lost and a
    wobble on the other axis never leaks. ``fine`` is applied per event, so
    Shift may be pressed or released mid-drag without a jump.
    """
    lock: str | None = None      # None until decided, then "p0" | "p1"
    dp0: float = 0.0             # cumulative degrees since button-down
    dp1: float = 0.0
    travel_x: float = 0.0        # withheld pixels while lock is None
    travel_y: float = 0.0

    def move(self, dx_px: float, dy_up_px: float,
             fine: bool = False) -> tuple[float, float] | None:
        """Feed one pointer step (``dy_up_px`` positive for a drag UP the
        screen). Returns the cumulative ``(dp0, dp1)``, or None while the
        gesture is still below the lock threshold."""
        k = FINE if fine else 1.0
        if self.lock is None:
            self.travel_x += float(dx_px)
            self.travel_y += float(dy_up_px)
            if max(abs(self.travel_x), abs(self.travel_y)) < LOCK_PX:
                return None
            self.lock = ("p0" if abs(self.travel_x) >= abs(self.travel_y)
                         else "p1")
            # apply the withheld travel with the CURRENT modifier
            dx_px, dy_up_px = self.travel_x, self.travel_y
            self.travel_x = self.travel_y = 0.0
        if self.lock == "p0":
            self.dp0 += P0_DEG_PER_PX * k * float(dx_px)
        else:
            self.dp1 += P1_DEG_PER_PX * k * float(dy_up_px)
        return (self.dp0, self.dp1)

    @property
    def changed(self) -> bool:
        """Has the gesture produced a phase change (locked and non-zero)?"""
        return self.lock is not None and bool(self.dp0 or self.dp1)
