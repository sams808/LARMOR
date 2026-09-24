"""A physical axis never carries an SI prefix (larmor.desktop.axes).

Sam, 2026-09-24: "when opening some data (my WURST-QCPMG for example), the
x axis in the plot area is broken: it says from 6 to -6 ppm where the data
range from about 6000 to -6000 ppm." pyqtgraph's AxisItem had rescaled the
ticks and quietly relabelled the unit "kppm". Every plot of every dialog now
goes through plain_units(); this file pins the two pyqtgraph behaviours the
helper exists for (a prefix frozen by a late enableAutoSIPrefix(False); the
(x0.001) scaling of a unit-less axis whose whole range fits within +-1), the
tick strings on the real 81Br spectrum, the main view, and the axis state of
every dialog that owns a plot.
"""
import os

import numpy as np
import pytest

from conftest import MAGLAB_81BR, require

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
import pyqtgraph as pg  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


#: the reproduction range: MAGLAB_81BR / "30" spans +-5787 ppm
WIDE_PPM = np.linspace(5787.0, -5787.0, 4096)
WIDE_AMP = np.exp(-(WIDE_PPM / 800.0) ** 2)


def _shown(qapp, widget, w=700, h=450):
    """Offscreen widgets keep the (0, 1) view range until they are painted."""
    widget.resize(w, h)
    widget.show()
    qapp.processEvents()
    return widget


def _tick_strings(ax, px=2000):
    """The major tick strings pyqtgraph draws for the axis' current range
    (an inverted NMR axis stores it high -> low)."""
    lo, hi = sorted(ax.range)
    spacing, values = ax.tickValues(lo, hi, px)[0]
    return ax.tickStrings(values, ax.autoSIPrefixScale * ax.scale, spacing)


def _spans_wideline(ax):
    lo, hi = sorted(ax.range)
    return lo < -5000 and hi > 5000


def _tick_span(ax):
    vals = [abs(float(s)) for s in _tick_strings(ax)]
    return max(vals)


def _labelled(widget):
    from larmor.desktop.axes import all_axes
    return sorted((ax.orientation, ax.labelText, ax.labelUnits)
                  for ax in all_axes(widget) if ax.label.isVisible())


def _assert_plain(widget, n_plots: int):
    """Every bottom / left axis and every labelled axis under ``widget`` has
    the prefix off; the walker saw ``n_plots`` plots (4 axes each)."""
    from larmor.desktop.axes import all_axes
    axes = list(all_axes(widget))
    assert len(axes) == 4 * n_plots, (len(axes), n_plots)
    bad = [(ax.orientation, ax.labelText, ax.labelUnits) for ax in axes
           if ax.autoSIPrefix
           and (ax.label.isVisible() or ax.orientation in ("bottom", "left"))]
    assert bad == [], f"axes still SI-prefixed: {bad}"


# --------------------------------------------------------------------------
# the helper itself
# --------------------------------------------------------------------------
def test_plain_units_takes_plotwidgets_and_plotitems(qapp):
    from larmor.desktop.axes import all_axes, plain_units

    pw = pg.PlotWidget()
    glw = pg.GraphicsLayoutWidget()
    p1 = glw.addPlot(row=0, col=0)
    p2 = glw.addPlot(row=1, col=0)
    for ax in (*all_axes(pw), *all_axes(glw)):
        assert ax.autoSIPrefix is True                     # pyqtgraph's default
    plain_units(pw, p1)                                    # bottom + left
    plain_units(p2, axes=("bottom", "left", "right"))
    state = {(id(ax.linkedView()), ax.orientation): ax.autoSIPrefix
             for ax in (*all_axes(pw), *all_axes(glw))}
    off = [k for k, v in state.items() if not v]
    assert len(off) == 2 + 2 + 3
    assert all(o in ("bottom", "left") for _, o in off
               if _ != id(p2.getViewBox()))
    assert state[(id(p2.getViewBox()), "right")] is False
    assert state[(id(p2.getViewBox()), "top")] is True   # not asked for
    with pytest.raises(TypeError):
        plain_units(pw.getViewBox())                       # no axes to speak of
    # all_axes: the plot widget itself, and every plot under a plain QWidget
    assert len(list(all_axes(pw))) == 4 and len(list(all_axes(glw))) == 8


def test_plain_units_resets_a_prefix_pyqtgraph_already_applied(qapp):
    """The reproduction, and the reason the helper is more than a one-liner:
    a ppm axis shown at +-5787 reads 6 ... -6 under "(kppm)"; a bare
    enableAutoSIPrefix(False) at that point recomputes and FREEZES the 'k'
    (pyqtgraph 0.14: updateAutoSIPrefix ignores the flag it was just given
    and is never called again once it is off). plain_units() undoes it."""
    from larmor.desktop.axes import plain_units

    pw = pg.PlotWidget()
    pw.setLabel("bottom", "chemical shift", units="ppm")
    pw.plot(WIDE_PPM, WIDE_AMP)
    _shown(qapp, pw)
    ax = pw.getPlotItem().getAxis("bottom")
    assert _spans_wideline(ax)                               # the view did autorange
    assert ax.labelUnitPrefix == "k" and _tick_span(ax) <= 10          # the bug
    assert "(kppm)" in ax.labelString()
    ax.enableAutoSIPrefix(False)                            # the naive fix ...
    assert ax.autoSIPrefixScale == 0.001 and _tick_span(ax) <= 10      # ... frozen
    plain_units(pw)
    qapp.processEvents()
    assert ax.autoSIPrefix is False and ax.autoSIPrefixScale == 1.0
    assert ax.labelUnitPrefix == "" and "(ppm)" in ax.labelString()
    assert _tick_span(ax) >= 4000
    # and it stays plain through later range changes
    pw.getViewBox().setXRange(-6500.0, 6500.0, padding=0)
    qapp.processEvents()
    assert ax.autoSIPrefixScale == 1.0 and _tick_span(ax) >= 4000


def test_unitless_axis_within_unity_stays_in_its_own_numbers(qapp):
    """Without a unit string pyqtgraph still scales an axis whose whole range
    fits within +-1 -- dS/S0, a normalised integral, an RMSD, a population
    fraction -- to ticks 0 ... 800 under "(x0.001)". Not on a plain axis."""
    from larmor.desktop.axes import plain_units

    x = np.linspace(0.0, 10.0, 50)
    y = 0.8 * np.linspace(0.0, 1.0, 50)
    control = pg.PlotWidget()
    control.setLabel("left", "ΔS/S₀")
    control.plot(x, y)
    _shown(qapp, control)
    cax = control.getPlotItem().getAxis("left")
    assert cax.autoSIPrefixScale == 1000.0 and "(x0.001)" in cax.labelString()
    assert _tick_span(cax) >= 100                            # 800, not 0.8

    plain = pg.PlotWidget()
    plain_units(plain)
    plain.setLabel("left", "ΔS/S₀")
    plain.plot(x, y)
    _shown(qapp, plain)
    pax = plain.getPlotItem().getAxis("left")
    assert pax.autoSIPrefixScale == 1.0 and "(x" not in pax.labelString()
    assert 0.5 <= _tick_span(pax) <= 1.0
    # the label reads exactly what was asked for
    assert "ΔS/S₀" in pax.labelString() and pax.labelString().count("(") == 0


def test_units_that_are_not_si_units_never_get_a_prefix(qapp):
    """1000/T in 1/K above 1000 K (a melt), a sub-MHz Czjzek C_Q axis: the
    prefixes pyqtgraph would invent are "m1/K" and "mMHz"."""
    from larmor.desktop.axes import plain_units

    for units, x in (("1/K", np.linspace(0.65, 0.95, 20)),
                     ("MHz", np.linspace(0.0, 0.6, 20))):
        pw = pg.PlotWidget()
        plain_units(pw)
        pw.setLabel("bottom", "axis", units=units)
        pw.plot(x, x)
        _shown(qapp, pw)
        ax = pw.getPlotItem().getAxis("bottom")
        assert ax.labelUnitPrefix == "" and f"({units})" in ax.labelString()
        assert _tick_span(ax) < 2.0, _tick_strings(ax)


# --------------------------------------------------------------------------
# the main view
# --------------------------------------------------------------------------
def test_spectrum_view_shows_a_wideline_spectrum_in_thousands_of_ppm(qapp):
    """Pinned: the main SpectrumView disables the prefix at construction and
    a +-5787 ppm spectrum reads -6000 ... 6000 under "(ppm)"."""
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    v.set_experiment(WIDE_PPM, WIDE_AMP)
    _shown(qapp, v)
    v.getPlotItem().getViewBox().autoRange()
    qapp.processEvents()
    ax = v.getPlotItem().getAxis("bottom")
    assert ax.autoSIPrefix is False and ax.autoSIPrefixScale == 1.0
    assert "(ppm)" in ax.labelString() and "kppm" not in ax.labelString()
    assert _tick_span(ax) >= 4000, _tick_strings(ax)


@pytest.mark.xfail(strict=True, reason=(
    "larmor/desktop/plot.py SpectrumView.apply_theme() re-calls "
    "enableAutoSIPrefix(False) on the bottom axis; run after a wideline "
    "spectrum is displayed (every theme switch, mw_menus._apply_theme_live) "
    "that call freezes the 'k' prefix: ticks 6 ... -6 under 'chemical shift "
    "(kppm)'. Fix: replace the call with plain_units(self) from "
    "larmor.desktop.axes (resets the scale), then drop this marker."))
def test_spectrum_view_keeps_thousands_of_ppm_after_a_theme_reapply(qapp):
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    v.set_experiment(WIDE_PPM, WIDE_AMP)
    _shown(qapp, v)
    v.getPlotItem().getViewBox().autoRange()
    qapp.processEvents()
    ax = v.getPlotItem().getAxis("bottom")
    assert _tick_span(ax) >= 4000
    v.apply_theme()                                   # what a theme switch does
    qapp.processEvents()
    assert ax.autoSIPrefixScale == 1.0 and "kppm" not in ax.labelString()
    assert _tick_span(ax) >= 4000, _tick_strings(ax)


# --------------------------------------------------------------------------
# the real dataset
# --------------------------------------------------------------------------
def test_81br_static_spectrum_is_ticked_in_thousands_of_ppm(qapp):
    """MAGLAB_81BR / "30": the WCPMG 81Br spectrum Sam opened, 65536 points
    over +-5787 ppm. In a batch-fit cell, in the main view, and on a plain
    ppm axis given the data first and the helper second."""
    from larmor.desktop.axes import plain_units
    from larmor.desktop.batchfit_dialog import BatchFitDialog
    from larmor.desktop.plot import SpectrumView
    from larmor.loader import load_any

    expno = require(MAGLAB_81BR / "30")
    ppm, amp, rec, _meta, _warns = load_any(str(expno))
    ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
    assert ppm.size == 65536 and rec["nucleus"] == "81Br"
    assert ppm.max() > 5500 and ppm.min() < -5500

    dlg = BatchFitDialog(None, [str(expno)], None)
    try:
        assert len(dlg._cells) == 1
        cell = _shown(qapp, dlg._cells[0]["plot"])
        ax = cell.getPlotItem().getAxis("bottom")
        assert ax.autoSIPrefix is False and ax.autoSIPrefixScale == 1.0
        assert _spans_wideline(ax)                          # inverted: high -> low
        assert _tick_span(ax) >= 4000, _tick_strings(ax)
        _assert_plain(dlg, n_plots=1)
    finally:
        dlg.close()

    v = SpectrumView()
    v.set_experiment(ppm, amp)
    _shown(qapp, v)
    v.getPlotItem().getViewBox().autoRange()
    qapp.processEvents()
    ax = v.getPlotItem().getAxis("bottom")
    assert "(ppm)" in ax.labelString() and _tick_span(ax) >= 4000

    late = pg.PlotWidget()
    late.setLabel("bottom", "chemical shift", units="ppm")
    late.plot(ppm, amp)
    _shown(qapp, late)
    ax = late.getPlotItem().getAxis("bottom")
    assert ax.labelUnitPrefix == "k" and _tick_span(ax) <= 10        # the report
    plain_units(late)
    qapp.processEvents()
    assert "(ppm)" in ax.labelString() and _tick_span(ax) >= 4000
