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


# --------------------------------------------------------------------------
# every dialog that owns a plot: constructed offscreen with minimal inputs,
# every bottom / left / labelled axis plain, the labels the walker must see
# --------------------------------------------------------------------------
def _csv_series(tmp_path, positions=(13.0, 15.0, 17.0)):
    """Three 11B spectra as CSVs + the shared model (as test_seqfit_ui)."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 60, 500)
    paths = []
    for k, pos in enumerate(positions):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                    sites=[SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(pos),
                        "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100),
                        "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        p = tmp_path / f"s{k}.csv"
        with open(p, "w", encoding="utf-8") as f:
            f.write("# nucleus = 11B\n# larmor_MHz = 160\n")
            for xi, yi in zip(x, m):
                f.write(f"{xi:.4f} {yi:.4f}\n")
        paths.append(str(p))
    return paths, _recipe_11b()


def _recipe_11b():
    return {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
            "sites": [{"model": "gauss_lor", "label": "A", "params": {
                "isotropic_chemical_shift_ppm": {"value": 12.0, "min": 0, "max": 30},
                "shift_fwhm_ppm": {"value": 5.0, "min": 0.1},
                "amplitude": {"value": 80.0, "min": 0},
                "gl": {"value": 1.0, "vary": False}}}]}


def _map_2d():
    from larmor import twod
    f2 = np.linspace(-50, 50, 200); f1 = np.linspace(-30, 30, 80)
    Z = (np.exp(-((f2[None, :] - 10) / 3) ** 2)
         * np.exp(-((f1[:, None] - 5) / 3) ** 2))
    return twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="1H", larmor_MHz=500.0)


def test_batchfit_grid_cells_are_plain(qapp, tmp_path):
    from larmor.desktop.batchfit_dialog import BatchFitDialog
    paths, model = _csv_series(tmp_path)
    dlg = BatchFitDialog(None, paths, model)
    try:
        _assert_plain(dlg, n_plots=3)
        assert _labelled(dlg) == []                 # the cells carry no axis label
    finally:
        dlg.close()


def test_seqfit_dialog_plots_are_plain(qapp, tmp_path):
    from larmor.desktop.seqfit_dialog import SeqFitDialog
    paths, model = _csv_series(tmp_path)
    dlg = SeqFitDialog(None, paths, model)
    try:
        _assert_plain(dlg, n_plots=3)               # spectrum, RMSD, trajectory
        assert _labelled(dlg) == [("bottom", "spectrum", ""),
                                  ("bottom", "spectrum", ""),
                                  ("left", "RMSD", "")]
    finally:
        dlg.close()


def test_cofit_dialog_result_plots_are_plain(qapp):
    from types import SimpleNamespace

    from larmor import twod
    from larmor.desktop.cofit_dialog import CofitDialog
    from larmor.recipe import Recipe

    x = np.linspace(-50, 120, 400)
    amp = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    base = {"kind": "1d", "label": "MAS", "ppm": x, "amp": amp,
            "nucleus": "27Al", "larmor": 195.5}
    f2 = np.linspace(-40, 110, 60); f1 = np.linspace(-30, 90, 50)
    Z = np.exp(-0.5 * (((f2[None, :] - 55) / 12) ** 2
                       + ((f1[:, None] - 60) / 10) ** 2))
    d2 = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="27Al", larmor_MHz=195.5)
    recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 195.5, "sites": [
        {"model": "czjzek", "label": "AlIV", "params": {
            "isotropic_chemical_shift_ppm": {"value": 60.0, "vary": True,
                                             "min": 0, "max": 120},
            "sigma_Cq_MHz": {"value": 2.0, "vary": True, "min": 0.2, "max": 8},
            "shift_fwhm_ppm": {"value": 12.0, "vary": True, "min": 1, "max": 30},
            "line_fwhm_ppm": {"value": 4.0, "vary": True, "min": 0},
            "amplitude": {"value": 1.0, "vary": True, "min": 0}}}]}
    dlg = CofitDialog(None, recipe, base)
    try:
        dlg.datasets = [base, {"kind": "2d", "label": "MQMAS", "data2d": d2,
                               "nucleus": "27Al", "larmor": 195.5}]
        dlg._result = SimpleNamespace(
            recipes=[Recipe.from_dict(recipe)], rmsd=[0.031, 0.048],
            per_dataset=[{"kind": "1d", "x": x, "y_fit": 1.05 * amp},
                         {"kind": "2d", "f2": f2, "f1": f1, "z_fit": Z * 0.98,
                          "per_site": [Z * 0.98]}])
        dlg._plot_result()
        qapp.processEvents()
        _assert_plain(dlg, n_plots=2)
        assert _labelled(dlg) == [("bottom", "F2 (ppm)", ""), ("bottom", "ppm", ""),
                                  ("left", "F1 (ppm)", "")]
    finally:
        dlg.close()


def test_contour_2d_view_map_projections_and_phasing_traces_are_plain(qapp):
    """The map (F2 on the bottom, F1 on the RIGHT axis), both projections,
    and the per-pick phasing traces the view builds on demand."""
    from larmor.desktop.twod_view import Contour2DView

    v = Contour2DView()
    v.set_data(_map_2d().normalized(), "syn")
    qapp.processEvents()
    _assert_plain(v, n_plots=3)
    assert _labelled(v) == [("bottom", "F2 (ppm)", ""), ("right", "F1", "")]
    v._picks = [10, 20]; v._pick_axis = "f2"; v._pivot = 10.0
    v._enter_phasing()
    qapp.processEvents()
    _assert_plain(v, n_plots=5)
    assert ("bottom", "ppm", "") in _labelled(v)


def test_twod_dialog_plots_are_plain(qapp):
    from larmor.desktop.twod_dialog import TwoDDialog
    dlg = TwoDDialog(None, None)
    try:
        _assert_plain(dlg, n_plots=3)
        assert _labelled(dlg) == [("bottom", "F2 (ppm)", ""),
                                  ("right", "F1 (ppm)", "")]
    finally:
        dlg.close()


def test_tool_dialogs_are_plain(qapp):
    from larmor.desktop.tool_dialogs import ErrorsDialog, RedorDialog
    x = np.linspace(-20, 60, 500)
    cases = ((RedorDialog(None, None),
              [("bottom", "recoupling time / s", ""), ("left", "ΔS/S₀", "")]),
             (ErrorsDialog(None, _recipe_11b(), x, np.ones_like(x), (60.0, -20.0)),
              [("bottom", "parameter value", ""), ("left", "χ²", "")]))
    for dlg, labels in cases:
        try:
            _assert_plain(dlg, n_plots=1)
            assert _labelled(dlg) == labels
        finally:
            dlg.close()


def test_montecarlo_dialog_is_plain(qapp):
    from larmor.desktop.montecarlo_dialog import MonteCarloDialog
    x = np.linspace(-20, 60, 500)
    dlg = MonteCarloDialog(None, _recipe_11b(), x, np.ones_like(x), (60.0, -20.0))
    try:
        _assert_plain(dlg, n_plots=1)
        assert _labelled(dlg) == [("bottom", "parameter value", ""),
                                  ("left", "count", "")]
    finally:
        dlg.close()


def test_series_plot_dialog_subplots_are_plain(qapp):
    from larmor.batchfit import BatchFitResult
    from larmor.desktop.series_plot import SeriesPlotDialog
    from larmor.recipe import Param, Recipe, SiteModel

    recs = [Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                   sample=f"g{k}", sites=[
                       SiteModel(model="gauss_lor", label="A", params={
                           "isotropic_chemical_shift_ppm": Param(pos),
                           "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100),
                           "gl": Param(1.0, vary=False)}),
                       SiteModel(model="gauss_lor", label="B", params={
                           "isotropic_chemical_shift_ppm": Param(2.0),
                           "shift_fwhm_ppm": Param(3.0), "amplitude": Param(50),
                           "gl": Param(1.0, vary=False)})])
            for k, pos in enumerate((15.0, 15.3, 14.7))]
    res = BatchFitResult(recipes=recs, labels=[f"g{k}" for k in range(3)],
                         rmsd=[0.0] * 3, per_dataset=[], shared=(), released=())
    dlg = SeriesPlotDialog(None, res)
    try:
        n = len(dlg._subplots)
        assert n >= 4
        _assert_plain(dlg, n_plots=n)
        labels = _labelled(dlg)
        assert len(labels) == n and all(o == "left" for o, _, _ in labels)
        assert ("left", "population % (integral)", "") in labels
    finally:
        dlg.close()


def test_vt_dialog_axis_stays_in_1_per_k(qapp):
    """1000/T for a melt (T > 1000 K) lies below 1: pyqtgraph would relabel
    the axis 'm1/K' and tick it 700 ... 900."""
    from larmor.desktop.vt_dialog import VtDialog
    dlg = VtDialog(None)
    try:
        _assert_plain(dlg, n_plots=1)
        assert _labelled(dlg) == [("bottom", "1000/T", "1/K"),
                                  ("left", "ln(rate)", "")]
        inv_t = 1000.0 / np.array([1123.0, 1223.0, 1323.0, 1423.0])
        dlg.plot.plot(inv_t, -3.0 * inv_t)
        _shown(qapp, dlg)
        ax = dlg.plot.getPlotItem().getAxis("bottom")
        assert max(ax.range) < 1.0                      # the range pyqtgraph scales
        assert "(1/K)" in ax.labelString()
        assert _tick_span(ax) < 1.0, _tick_strings(ax)
    finally:
        dlg.close()


def test_czjzek_dist_dialog_axis_stays_in_mhz(qapp):
    """Zoomed to the sub-MHz part of a 23Na distribution the C_Q axis would
    read 'mMHz' with ticks 100 ... 600."""
    from larmor.desktop.axes import all_axes
    from larmor.desktop.czjzek_dist_dialog import CzjzekDistDialog
    r = {"nucleus": "23Na", "sites": [{"model": "czjzek", "label": "Na",
                                       "params": {"sigma_Cq_MHz": {"value": 0.05}}}]}
    dlg = CzjzekDistDialog(None, r)
    try:
        _shown(qapp, dlg)
        _assert_plain(dlg, n_plots=1)
        assert _labelled(dlg) == [("bottom", "C_Q", "MHz"), ("left", "P(C_Q)", "")]
        ax = next(a for a in all_axes(dlg) if a.orientation == "bottom")
        ax.linkedView().setXRange(0.0, 0.6, padding=0)
        qapp.processEvents()
        assert "(MHz)" in ax.labelString() and ax.labelUnitPrefix == ""
        assert _tick_span(ax) < 1.0, _tick_strings(ax)
    finally:
        dlg.close()


def test_satrec_dialog_plots_are_plain(qapp):
    from larmor.desktop.satrec_dialog import SatrecDialog
    dlg = SatrecDialog(None, None)
    try:
        _assert_plain(dlg, n_plots=2)                   # spectrum + build-up
        assert _labelled(dlg) == [("bottom", "delay (s)", ""),
                                  ("bottom", "shift", "ppm"),
                                  ("left", "integral (norm.)", "")]
    finally:
        dlg.close()


def test_baseline_dialog_plots_are_plain(qapp):
    from larmor.desktop.baseline_dialog import BaselineDialog
    x = np.linspace(-20, 60, 500)
    dlg = BaselineDialog(None, x, np.exp(-((x - 15) / 4) ** 2))
    try:
        _assert_plain(dlg, n_plots=2)
        assert _labelled(dlg) == [("bottom", "shift", "ppm"),
                                  ("left", "corrected", ""),
                                  ("left", "intensity", "")]
    finally:
        dlg.close()
