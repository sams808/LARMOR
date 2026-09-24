"""Plain physical axes: never let pyqtgraph SI-prefix a unit.

pyqtgraph's AxisItem rescales its tick values to the displayed range and
moves the factor into the unit label: a static 81Br spectrum spanning
+-5787 ppm is drawn as ticks 6 ... -6 under "chemical shift (kppm)" (Sam,
2026-09-24, on a WURST-QCPMG dataset). "kppm" is not a unit, and neither is
"mMHz" (a sub-MHz Czjzek C_Q axis) or "m1/K" (1000/T above 1000 K). Without
a unit string the same machinery still scales any axis whose whole range
fits within +-1 -- normalised integrals, dS/S0, RMSDs, population fractions
are ticked 0 ... 800 under "(x0.001)".

``enableAutoSIPrefix(False)`` alone is not enough: pyqtgraph recomputes the
prefix inside that very call and never again once disabled, so an axis
that already showed +-6000 keeps a stale 'k' for good (this is what a
theme switch does to a plot that re-runs the call). plain_units() also
resets the scale, so it is correct whenever it is called.
"""
from __future__ import annotations

from collections.abc import Iterable, Iterator

import pyqtgraph as pg
from PySide6.QtWidgets import QWidget

#: the axes a plot labels by default (the 2D map puts F1 on "right")
DEFAULT_AXES = ("bottom", "left")
_ALL_AXES = ("bottom", "left", "right", "top")


def _plot_item(plot) -> pg.PlotItem:
    if isinstance(plot, pg.PlotItem):
        return plot
    get = getattr(plot, "getPlotItem", None)          # PlotWidget
    if callable(get):
        return get()
    raise TypeError(f"not a PlotWidget / PlotItem: {type(plot).__name__}")


def plain_axis(ax: pg.AxisItem) -> None:
    """Disable the SI prefix on ONE axis and undo a prefix already applied."""
    ax.enableAutoSIPrefix(False)
    if ax.autoSIPrefixScale != 1.0 or ax.labelUnitPrefix:
        ax.autoSIPrefixScale = 1.0
        ax.labelUnitPrefix = ""
        ax._updateLabel()                             # what pyqtgraph itself calls


def plain_units(*plots, axes: Iterable[str] = DEFAULT_AXES) -> None:
    """Call right after creating a plot: its ``axes`` keep the unit they are
    labelled with -- ticks in ppm / Hz / MHz / s / %, never kppm or (x0.001).

    Accepts PlotWidgets (via getPlotItem) and PlotItems (what
    GraphicsLayoutWidget.addPlot returns)."""
    for plot in plots:
        pi = _plot_item(plot)
        for name in axes:
            plain_axis(pi.getAxis(name))


def all_axes(widget: QWidget) -> Iterator[pg.AxisItem]:
    """Every AxisItem of every PlotItem under ``widget`` (the widget itself
    included when it is a plot), for a test to assert the state of a whole
    dialog. PlotWidget and GraphicsLayoutWidget are both GraphicsViews."""
    views = list(widget.findChildren(pg.GraphicsView))
    if isinstance(widget, pg.GraphicsView):
        views.insert(0, widget)
    seen: set[int] = set()
    for view in views:
        scene = view.scene()
        if scene is None:
            continue
        for item in scene.items():
            if isinstance(item, pg.PlotItem):
                for name in _ALL_AXES:
                    ax = item.getAxis(name)
                    if id(ax) not in seen:
                        seen.add(id(ax))
                        yield ax
