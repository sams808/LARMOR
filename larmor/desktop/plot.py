"""Central spectrum view: pyqtgraph plot with draggable site markers and
click-to-add. All rendering is direct numpy -> GPU-backed canvas: no browser,
no serialization, instant interaction even on 32k-point spectra."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont

from larmor import cellparse
from larmor.desktop import theme
from larmor.desktop.axes import plain_units
from larmor.phasedrag import DragGesture

#: fallback categorical palette (used before a theme is applied / in headless code)
SITE_COLORS = theme.LIGHT_SERIES


def site_color(i: int) -> str:
    """Categorical colour for site i, from the active theme's series palette."""
    series = theme.active().series or SITE_COLORS
    return series[i % len(series)]


class ScaledAxis(pg.AxisItem):
    """A bottom axis that can DISPLAY the ppm-space data in kHz or MHz.

    Wideline patterns are read and reported in kHz from the reference, so the
    axis must speak that language -- but every plotted coordinate (curves,
    paddles, zones, markers) stays in ppm, the app's one internal unit. Only
    the tick VALUES are transformed, and the ticks are chosen round in the
    display unit (pyqtgraph's own setScale() keeps ppm-round ticks, which
    lands kHz labels like 43.2 / 86.4)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._factor = 1.0                    # display value = ppm * factor

    def set_factor(self, factor: float):
        self._factor = float(factor) if factor else 1.0
        self.picture = None
        self.update()

    def tickValues(self, minVal, maxVal, size):
        f = self._factor
        if f == 1.0:
            return super().tickValues(minVal, maxVal, size)
        vals = super().tickValues(minVal * f, maxVal * f, size)
        return [(spacing / f, [v / f for v in ticks])
                for spacing, ticks in vals]

    def tickStrings(self, values, scale, spacing):
        f = self._factor
        if f == 1.0:
            return super().tickStrings(values, scale, spacing)
        return super().tickStrings([v * f for v in values], scale, spacing * f)


#: display-axis units: label text and the factor from (ppm, sfo_MHz) -- the
#: ppm -> Hz convention is nu = delta * SFO (convert.ppm_to_Hz)
AXIS_UNITS = ("ppm", "kHz", "MHz")


def axis_factor(unit: str, sfo_MHz: float) -> float:
    """display value = ppm * factor. 1.0 for ppm or when SFO is unknown."""
    if unit == "kHz" and sfo_MHz > 0:
        return sfo_MHz / 1000.0
    if unit == "MHz" and sfo_MHz > 0:
        return sfo_MHz / 1e6
    return 1.0


class AnchoredViewBox(pg.ViewBox):
    """A ViewBox whose wheel-zoom keeps the data point **under the cursor** fixed.

    pyqtgraph's default anchors full-plot zoom at the mouse, but single-axis zoom
    (scrolling over an axis) computes its centre from the axis-item-local position
    mapped through the *viewbox* transform — which is wrong, so the axis zooms
    about its middle. Re-deriving the centre from the scene position fixes both:
    scroll on the x-axis at 100 ppm and it zooms about 100 ppm; likewise for y.

    It is also where a TopSpin drag-to-phase gesture is intercepted: pyqtgraph
    delivers a drag to the ViewBox only when no item (pivot, paddle, ruler,
    anchor, zone) accepted it first, so routing a plain left-button canvas
    drag to ``phase_drag_handler`` leaves every item-owned gesture, the axis
    drags (axis=0/1), the middle/right buttons and Ctrl+left (pan) untouched.
    """

    #: callable(MouseDragEvent) armed by SpectrumView.set_phase_drag_mode;
    #: None = the ViewBox pans/zooms as usual
    phase_drag_handler = None

    def phase_drag_takes(self, ev, axis=None) -> bool:
        """Does this drag belong to the phase gesture (left button on the
        canvas, no Ctrl, a handler armed)?"""
        return (self.phase_drag_handler is not None and axis is None
                and ev.button() == Qt.LeftButton
                and not (ev.modifiers() & Qt.ControlModifier))

    def mouseDragEvent(self, ev, axis=None):
        if self.phase_drag_takes(ev, axis):
            ev.accept()
            self.phase_drag_handler(ev)
            return
        super().mouseDragEvent(ev, axis)

    def wheelEvent(self, ev, axis=None):
        if axis in (0, 1):
            mask = [False, False]
            mask[axis] = self.state["mouseEnabled"][axis]
        else:
            mask = self.state["mouseEnabled"][:]
        if not any(mask):
            ev.ignore()
            return
        s = 1.02 ** (ev.delta() * self.state["wheelScaleFactor"])
        s = [(None if m is False else s) for m in mask]
        center = self.mapSceneToView(ev.scenePos())     # the point under the cursor
        self._resetTarget()
        self.scaleBy(s, center)
        ev.accept()
        self.sigRangeChangedManually.emit(mask)


class SpectrumView(pg.PlotWidget):
    """Experiment + model + components + residual, with dmfit-style paddles."""

    add_requested = Signal(float, float)      # (ppm, amplitude) from a click
    exit_add_mode = Signal()                  # right-click while placing lines
    marker_moved = Signal(int, float)         # legacy: (site index, new ppm)
    paddle_moved = Signal(int, float, float, float)   # index, pos, amp, fwhm
    paddle_released = Signal(int)
    cursor_moved = Signal(float, float)       # live x/y for the status bar
    file_dropped = Signal(str)                # a data file dragged onto the plot
    files_dropped_overlay = Signal(list)      # Shift + drop: overlay, keep the fit
    #: View > Overlays off: the compared spectra stay listed but are not drawn
    _overlays_hidden = False
    calibrate_picked = Signal(float)          # snapped peak ppm to reference
    measure_changed = Signal(float, float)    # two ppm cursors (ruler)
    #: a REAL frequency-domain spectrum was placed on the canvas
    #: (set_experiment) -- the workbench's hook to re-validate its display
    #: state (FID / imaginary channel) against the new data
    experiment_set = Signal()
    #: drag-to-phase: cumulative (dp0, dp1) in degrees since button-down
    phase_dragged = Signal(float, float)
    phase_drag_released = Signal()            # the button came up
    pivot_moved = Signal(float)               # new pivot fraction after a drag

    def __init__(self, parent=None):
        t = theme.active()
        super().__init__(parent, background=t.plot_bg, viewBox=AnchoredViewBox(),
                         axisItems={"bottom": ScaledAxis(orientation="bottom")})
        pi = self.getPlotItem()
        pi.invertX(True)                              # ppm convention
        self._axis_unit = "ppm"                       # display unit only
        self._axis_sfo_MHz = 0.0
        tick_font = QFont()
        tick_font.setPointSize(9)
        for name in ("bottom", "left"):
            pi.getAxis(name).setStyle(tickFont=tick_font, tickLength=-5)
        pi.showAxis("top"); pi.getAxis("top").setStyle(showValues=False)
        pi.showAxis("right"); pi.getAxis("right").setStyle(showValues=False)
        pi.getViewBox().setMouseMode(pg.ViewBox.PanMode)
        pi.setContentsMargins(6, 10, 6, 6)

        leg = self.addLegend(offset=(12, 10), labelTextSize="9pt")
        self._legend = leg

        self._exp = self.plot([], [], pen=pg.mkPen(t.experiment, width=1.4),
                              name="experiment", antialias=True)
        # Every MODEL-side item (model, residual, components, fit animation,
        # paddles, markers) is added with ignoreBounds=True: pyqtgraph's
        # auto-range unions the bounds of every added item, and the model is
        # simulated on the Czjzek KERNEL axis (engine.make_context: at least
        # 150 kHz wide, 1.25x the data span otherwise), so each simulation
        # or fit result used to stretch the x axis past the data -- and the
        # fit-animation label, pinned to the view corner, fed the range back
        # into itself frame after frame (2500, 5000, 7000 ppm on a long fit).
        # Only the experiment (and compared overlays) define the auto range;
        # the model is also drawn only across the data (_model_mask).
        self._model = self._add_model_curve(pen=pg.mkPen(t.model, width=1.8),
                                            name="model")
        self._resid = self._add_model_curve(pen=pg.mkPen(t.resid, width=1.0),
                                            name="residual")
        # faint zero line for the offset residual strip
        self._resid_zero = pg.InfiniteLine(
            angle=0, pen=pg.mkPen(t.resid_zero, width=1, style=Qt.DotLine))
        self._resid_zero.setVisible(False)
        pi.addItem(self._resid_zero, ignoreBounds=True)
        pi.getAxis("left").setStyle(tickTextWidth=48, autoExpandTextSpace=False)
        # peak-mode downsampling + clip-to-view: a zero-filled wideline
        # spectrum is 64k points, and repainting every one of them on every
        # pan made the window crawl. "peak" keeps min AND max per screen
        # column, so noise spikes and sharp horns stay honest.
        for c in (self._exp, self._model, self._resid):
            self._tune_curve(c)
        self._components: list[pg.PlotDataItem] = []
        self._markers: list[pg.InfiniteLine] = []
        self._add_mode: str | None = None
        self.show_components = True
        self.show_residual = True
        # component names: pinned at each maximum (View > Component labels)
        # or shown for the component under the cursor when not pinned
        self.show_labels = False
        self._comp_labels: list[pg.TextItem] = []
        self._comp_data: list[tuple] = []      # (site, x, y, text) visible
        self._hover_label: pg.TextItem | None = None

        self.scene().sigMouseClicked.connect(self._on_click)
        self.scene().sigMouseMoved.connect(self._on_move)
        self._paddles: list = []
        self.setAcceptDrops(True)                # drag a spectrum onto the plot

        # manual baseline: draggable anchors + live PCHIP preview
        self._bl_mode = False
        self._bl_anchors: list[pg.TargetItem] = []
        self._bl_curve = self.plot([], [], pen=pg.mkPen(t.baseline, width=1.4,
                                                        style=Qt.DashLine))
        # dmfit-style fit zones
        self._zones: list[pg.LinearRegionItem] = []
        # calibrate (click a peak) and measure (two-cursor ruler)
        self._cal_mode = False
        self._measure_lines: list[pg.InfiniteLine] = []
        # comparison overlays (other datasets drawn behind the active one)
        self._overlay_items: list[pg.PlotDataItem] = []
        # phase pivot (TopSpin-style): p1 rotates about this reference point
        self._pivot: pg.InfiniteLine | None = None
        # drag-to-phase mode (set_phase_drag_mode) and the gesture in flight
        self._phase_mode = False
        self._pd_gesture: DragGesture | None = None

        # display domain: "freq" (ppm axis, every item) or "time" (the FID on
        # a non-inverted ms axis, every ppm item hidden). The real frequency
        # trace is remembered so the pivot, calibrate snapping and the return
        # from the FID never read the displayed channel or a ms axis.
        self._domain = "freq"
        self._freq_x: np.ndarray | None = None
        self._freq_y: np.ndarray | None = None
        self._freq_ranges = None            # (xRange, yRange) saved on entry
        self._fid_aq: float | None = None   # AQ (ms) of the FID last drawn
        self._trace_label = "experiment"

        # live fit animation: a bright morphing curve + fading ghost trail
        self._anim_ghosts: list[pg.PlotDataItem] = []
        self._anim_main: pg.PlotDataItem | None = None
        self._anim_hist: list = []
        self._anim_label: pg.TextItem | None = None

        self.apply_theme()          # axis pens, labels, grid, legend from theme
        # LARMOR export + "send to Plotting studio" (and drop the freeze-prone
        # native pyqtgraph exporter)
        try:
            from larmor.desktop.plot_menu import attach_plot_menu
            attach_plot_menu(self, title="spectrum")
        except Exception:
            pass

    #: fraction of the data span drawn beyond each end of the experiment
    #: when the model is masked to the data (a small margin, so a curve that
    #: spills just past the last point is not cut visibly)
    MODEL_MARGIN_FRAC = 0.02

    def _add_model_curve(self, pen=None, name=None, **kw) -> pg.PlotDataItem:
        """A curve that never takes part in auto-range (see __init__)."""
        item = pg.PlotDataItem([], [], pen=pen, name=name, antialias=True, **kw)
        self.getPlotItem().addItem(item, ignoreBounds=True)
        return item

    def _model_mask(self, x, exp_x=None):
        """Boolean mask of the model points inside the experiment's x range
        (plus MODEL_MARGIN_FRAC of its span on each side); None when there is
        no experiment to mask against."""
        ref = exp_x if exp_x is not None and len(exp_x) else self._freq_x
        if ref is None or not len(ref) or x is None:
            return None
        ref = np.asarray(ref, float)
        lo, hi = float(np.min(ref)), float(np.max(ref))
        m = self.MODEL_MARGIN_FRAC * (hi - lo)
        xa = np.asarray(x, float)
        return (xa >= lo - m) & (xa <= hi + m)

    @staticmethod
    def _masked(mask, *arrays):
        if mask is None:
            return arrays
        return tuple(None if a is None else np.asarray(a)[mask] for a in arrays)

    # ---------------------------------------------------------------- fit animation
    def start_fit_animation(self):
        """Prepare the animated-fit overlay (call before a fit begins)."""
        if self._domain == "time":          # ppm curves never go on the ms axis
            return
        t = theme.active()
        if self._anim_main is None:
            self._anim_ghosts = []
            for _ in range(3):                      # a short fading trail
                g = self._add_model_curve(pen=pg.mkPen(t.accent, width=1))
                g.setZValue(40)
                self._anim_ghosts.append(g)
            self._anim_main = self._add_model_curve(
                pen=pg.mkPen(t.accent, width=2.2))
            self._anim_main.setZValue(45)
            self._anim_label = pg.TextItem(color=t.accent, anchor=(0, 0))
            self._anim_label.setZValue(46)
            # pinned to the view corner each frame: with bounds it drove the
            # auto-range feedback loop described in __init__
            self.getPlotItem().addItem(self._anim_label, ignoreBounds=True)
        self._anim_hist = []
        self._anim_label.setText("")
        for it in (*self._anim_ghosts, self._anim_main):
            it.setData([], [])
        self._anim_label.setVisible(True)

    def set_fit_frame(self, x, y, iteration: int, rms: float | None = None):
        """Show the current model curve, pushing the previous ones into a fading
        trail — so convergence (or divergence) is visible as it happens."""
        if self._domain == "time":
            return
        if self._anim_main is None:
            self.start_fit_animation()
        from PySide6.QtGui import QColor
        x = np.asarray(x, float); y = np.asarray(y, float)
        x, y = self._masked(self._model_mask(x), x, y)     # data range only
        self._anim_hist.append((x, y))
        self._anim_hist = self._anim_hist[-4:]      # main + up to 3 ghosts
        ghosts = self._anim_hist[:-1]
        base = QColor(theme.active().accent)
        alphas = [45, 80, 120]
        n = len(self._anim_ghosts)
        for gi, g in enumerate(self._anim_ghosts):
            idx = gi - (n - len(ghosts))            # newest ghost is brightest
            if 0 <= idx < len(ghosts):
                gx, gy = ghosts[idx]
                c = QColor(base); c.setAlpha(alphas[min(gi, 2)])
                g.setPen(pg.mkPen(c, width=1))
                g.setData(gx, gy)
            else:
                g.setData([], [])
        self._anim_main.setData(x, y)
        # a small live read-out anchored to the top-left of the (inverted-x) view
        try:
            vb = self.getPlotItem().getViewBox()
            xr, yr = vb.viewRange()
            self._anim_label.setPos(max(xr), max(yr))
        except Exception:
            pass
        txt = f"● fitting — iter {iteration}"
        if rms is not None and np.isfinite(rms):
            txt += f" · rms {rms:.3g}"
        self._anim_label.setText(txt)

    def stop_fit_animation(self):
        """Clear the animated-fit overlay (the final model takes over)."""
        for it in (*self._anim_ghosts, self._anim_main):
            if it is not None:
                it.setData([], [])
        if self._anim_label is not None:
            self._anim_label.setText(""); self._anim_label.setVisible(False)
        self._anim_hist = []

    def apply_theme(self):
        """(Re)apply the active colour theme to the plot's persistent items.
        Dynamic items (site markers, components, paddles) pick up the new series
        colours on the next redraw."""
        t = theme.active()
        self.setBackground(t.plot_bg)
        pi = self.getPlotItem()
        axis_pen = pg.mkPen(t.axis, width=1.2)
        for name in ("bottom", "left"):
            ax = pi.getAxis(name)
            ax.setPen(axis_pen)
            ax.setTextPen(pg.mkPen(t.axis))
        pi.getAxis("top").setPen(pg.mkPen(t.axis_minor))
        pi.getAxis("right").setPen(pg.mkPen(t.axis_minor))
        label_style = {"color": t.axis, "font-size": "10pt"}
        self._apply_axis_label(label_style)
        # never let pyqtgraph SI-prefix a ppm axis ("kppm" is not a unit)
        # plain_units, not enableAutoSIPrefix(False): pyqtgraph recomputes the
        # prefix INSIDE that call and never again once off, so re-applying
        # it from apply_theme() while a +-6000 ppm spectrum was shown froze a
        # "k" prefix -- ticks 6 ... -6 under a "(kppm)" label (Sam, 2026-09-24)
        plain_units(pi, axes=("bottom",))
        self.setLabel("left", "intensity", **label_style)
        self.showGrid(x=True, y=True, alpha=t.grid_alpha)
        self._exp.setPen(pg.mkPen(t.experiment, width=1.4))
        self._model.setPen(pg.mkPen(t.model, width=1.8))
        self._resid.setPen(pg.mkPen(t.resid, width=1.0))
        self._resid_zero.setPen(pg.mkPen(t.resid_zero, width=1, style=Qt.DotLine))
        self._bl_curve.setPen(pg.mkPen(t.baseline, width=1.4, style=Qt.DashLine))
        if self._legend is not None:
            self._legend.setLabelTextColor(t.text)
            try:
                r, g, b = theme._rgb(t.legend_bg)
                self._legend.setBrush(pg.mkBrush(r, g, b, 235))
                self._legend.setPen(pg.mkPen(t.border_soft))
            except Exception:
                pass
        ph = getattr(self, "_placeholder", None)
        if ph is not None and ph.isVisible():
            ph.setStyleSheet(f"color: {t.text_dim}; font-size: 13px; "
                             "background: transparent;")

    @staticmethod
    def _tune_curve(item):
        """Downsampling + view clipping for a data-carrying curve."""
        item.setDownsampling(auto=True, method="peak")
        item.setClipToView(True)

    # ---------- axis display unit ----------
    def set_axis_unit(self, unit: str, sfo_MHz: float):
        """Display the bottom axis in ppm, kHz or MHz. Data, paddles, zones
        and every plotted coordinate stay in ppm; only tick labels change."""
        self._axis_unit = unit if unit in AXIS_UNITS else "ppm"
        self._axis_sfo_MHz = float(sfo_MHz or 0.0)
        ax = self.getPlotItem().getAxis("bottom")
        if isinstance(ax, ScaledAxis):
            ax.set_factor(axis_factor(self._axis_unit, self._axis_sfo_MHz))
        t = theme.active()
        self._apply_axis_label({"color": t.axis, "font-size": "10pt"})

    @property
    def domain(self) -> str:
        """"freq" while the spectrum is shown, "time" while the FID is."""
        return self._domain

    def set_trace_label(self, text: str):
        """Rename the experiment curve's legend entry ('experiment',
        'experiment (imag)', 'FID (real)', ...) -- also what an export says."""
        if text == self._trace_label:
            return
        self._trace_label = text
        try:
            lab = self._legend.getLabel(self._exp)
        except Exception:
            lab = None
        if lab is not None:
            lab.setText(text)

    def _apply_axis_label(self, label_style: dict):
        if self._domain == "time":         # a theme change keeps the ms label
            self.setLabel("bottom", "time", units="ms", **label_style)
            return
        unit = self._axis_unit
        if unit != "ppm" and self._axis_sfo_MHz > 0:
            self.setLabel("bottom", "frequency offset", units=unit,
                          **label_style)
        else:
            self.setLabel("bottom", "chemical shift", units="ppm",
                          **label_style)

    # ---------- drag & drop ----------
    def dragEnterEvent(self, ev):
        if ev.mimeData().hasUrls():
            ev.acceptProposedAction()
        else:
            super().dragEnterEvent(ev)

    def dragMoveEvent(self, ev):
        if ev.mimeData().hasUrls():
            ev.acceptProposedAction()
        else:
            super().dragMoveEvent(ev)

    def dropEvent(self, ev):
        urls = ev.mimeData().urls()
        if urls:
            paths = [u.toLocalFile() for u in urls if u.toLocalFile()]
            if paths and (ev.modifiers() & Qt.ShiftModifier):
                # Shift held: compare on top of the active spectrum, the
                # fit stays; every dropped file becomes an overlay
                self.files_dropped_overlay.emit(paths)
            else:
                self.file_dropped.emit(paths[0] if paths else urls[0].toLocalFile())
            ev.acceptProposedAction()
        else:
            super().dropEvent(ev)

    def _on_move(self, scene_pos):
        vb = self.getPlotItem().getViewBox()
        if self.sceneBoundingRect().contains(scene_pos):
            p = vb.mapSceneToView(scene_pos)
            self.cursor_moved.emit(float(p.x()), float(p.y()))
            if not self.show_labels and self._domain == "freq":
                self._hover_component(float(p.x()))

    # ---------- add mode ----------
    def set_add_mode(self, model_name: str | None):
        self._add_mode = model_name
        self.setCursor(Qt.CrossCursor if model_name else Qt.ArrowCursor)
        # while placing lines, right-click exits the mode instead of opening the
        # viewbox menu — so suppress that menu until the mode ends. NOTE:
        # pyqtgraph's setMenuEnabled(False) DESTROYS the ViewBoxMenu object, and
        # setMenuEnabled(True) builds a brand-new default one from scratch — so
        # re-enabling silently wiped our custom "Export figure…" / "Send to
        # Plotting studio" items (added once at construction) the first time add
        # mode was ever toggled off. Re-attach them every time the menu comes
        # back so they survive every add/exit cycle, not just the first.
        self.getPlotItem().getViewBox().setMenuEnabled(model_name is None)
        if model_name is None:
            try:
                from larmor.desktop.plot_menu import attach_plot_menu
                attach_plot_menu(self, title="spectrum")
            except Exception:
                pass

    def _on_click(self, ev):
        if ev.button() == Qt.RightButton:
            if self._add_mode is not None:      # right-click once to stop adding
                ev.accept()
                self.exit_add_mode.emit()
            return
        if ev.button() != Qt.LeftButton:
            return
        if self._domain == "time":
            return          # no add-site / calibrate / baseline picks on a ms axis
        vb = self.getPlotItem().getViewBox()
        if not self.sceneBoundingRect().contains(ev.scenePos()):
            return
        p = vb.mapSceneToView(ev.scenePos())
        if self._cal_mode:
            self.calibrate_picked.emit(self._snap_peak(float(p.x())))
            ev.accept()
            return
        if self._bl_mode:
            self._add_baseline_anchor(float(p.x()), float(p.y()))
            ev.accept()
            return
        if self._add_mode is not None:
            self.add_requested.emit(float(p.x()), abs(float(p.y())))
            ev.accept()

    # ---------- calibrate (reference a peak) ----------
    def set_calibrate_mode(self, on: bool):
        self._cal_mode = on
        self.setCursor(Qt.CrossCursor if on else Qt.ArrowCursor)

    def _snap_peak(self, x_click: float) -> float:
        """Snap a click to the nearest local maximum of the experiment (the
        real frequency trace, whatever channel is displayed)."""
        x, y = self._freq_x, self._freq_y
        if x is None or not len(x):
            return x_click
        span = 0.01 * (float(np.max(x)) - float(np.min(x)))
        near = np.abs(x - x_click) <= max(span, 1e-9)
        if not near.any():
            return x_click
        idx = np.where(near)[0]
        return float(x[idx[int(np.argmax(y[idx]))]])

    # ---------- measure (two-cursor ruler) ----------
    def set_measure_mode(self, on: bool):
        for ln in self._measure_lines:
            self.removeItem(ln)
        self._measure_lines.clear()
        if not on:
            return
        if self._domain == "time" and self._freq_ranges is not None:
            (x0, x1) = self._freq_ranges[0]     # ppm cursors, even from FID view
        else:
            (x0, x1), _ = self.getPlotItem().getViewBox().viewRange()
        lo, hi = min(x0, x1), max(x0, x1)
        mc = theme.active().measure
        for frac in (0.35, 0.65):
            ln = pg.InfiniteLine(pos=lo + frac * (hi - lo), angle=90,
                                 movable=True,
                                 pen=pg.mkPen(mc, width=1.4),
                                 hoverPen=pg.mkPen(mc, width=2.4),
                                 label="{value:.1f}",
                                 labelOpts={"color": mc, "position": 0.92})
            ln.sigPositionChanged.connect(self._emit_measure)
            self.addItem(ln)
            ln.setVisible(self._domain == "freq")
            self._measure_lines.append(ln)
        self._emit_measure()

    def _emit_measure(self, *_):
        if len(self._measure_lines) == 2:
            self.measure_changed.emit(float(self._measure_lines[0].value()),
                                      float(self._measure_lines[1].value()))

    # ---------- phase pivot (TopSpin-style first-order phasing) ----------
    def show_phase_pivot(self, on: bool):
        """Show/hide the draggable p1 pivot. Defaults to the tallest peak so
        first-order phasing leaves that peak in phase (as TopSpin does)."""
        if not on and self._phase_mode:
            # the dock's visibilityChanged(False) must not remove the pivot
            # mid-phasing: phase_pivot_frac() would fall back to 0.5 and p1
            # would silently re-rotate about the centre at the next move
            return
        if not on:
            if self._pivot is not None:
                self.removeItem(self._pivot)
                self._pivot = None
            return
        if self._pivot is not None:
            return
        # the tallest REAL peak, whatever channel (or the FID) is displayed
        x, y = self._freq_x, self._freq_y
        if x is None or not len(x):
            return
        px = (float(x[int(np.argmax(y))]) if y is not None and len(y)
              else float(np.median(x)))
        pv = theme.active().pivot
        self._pivot = pg.InfiniteLine(
            pos=px, angle=90, movable=True,
            pen=pg.mkPen(pv, width=1.4, style=Qt.DashLine),
            hoverPen=pg.mkPen(pv, width=2.4),
            label="pivot", labelOpts={"color": pv, "position": 0.06})
        self.addItem(self._pivot)
        self._pivot.setVisible(self._domain == "freq")
        self._pivot.sigPositionChangeFinished.connect(
            lambda *_: self.pivot_moved.emit(self.phase_pivot_frac()))

    def phase_pivot_frac(self) -> float:
        """The pivot as a 0..1 fraction along the data (op_phase convention).
        Read on the frequency axis -- with the FID shown, _exp.xData is
        milliseconds and would give a garbage pivot on every live tick."""
        x = self._freq_x
        if self._pivot is None or x is None or not len(x):
            return 0.5
        idx = int(np.argmin(np.abs(np.asarray(x) - float(self._pivot.value()))))
        return idx / max(len(x) - 1, 1)

    # ---------- drag to phase (TopSpin gesture) ----------
    def set_phase_drag_mode(self, on: bool):
        """Arm / release the drag-to-phase gesture: a left-button drag on empty
        canvas is routed to the phasedrag arithmetic (horizontal = p0,
        vertical = p1 about the pivot) instead of panning. The pivot is
        created at once when data exists (the dock's visibility event never
        fires offscreen or with a hidden parent) and survives a dock hide
        while the mode is on (see show_phase_pivot)."""
        self._phase_mode = bool(on)
        self._pd_gesture = None
        if on:
            self.show_phase_pivot(True)
        self._arm_phase_drag()

    def phase_drag_active(self) -> bool:
        return self._phase_mode

    def _arm_phase_drag(self):
        """The ViewBox takes canvas drags only while the mode is on AND the
        frequency axis is shown -- phasing a displayed FID is meaningless, so
        a left drag pans the FID as usual and the cursor says so."""
        live = self._phase_mode and self._domain == "freq"
        vb = self.getPlotItem().getViewBox()
        vb.phase_drag_handler = self._phase_drag if live else None
        if live:
            self.setCursor(Qt.SizeAllCursor)
        elif self._phase_mode or self.cursor().shape() == Qt.SizeAllCursor:
            self.setCursor(Qt.ArrowCursor)

    def _phase_drag(self, ev):
        """One MouseDragEvent of the gesture (called by AnchoredViewBox).
        Scene units are logical pixels in a pg.GraphicsView; scene y grows
        downward, so a drag UP is +p1."""
        if ev.isStart() or self._pd_gesture is None:
            self._pd_gesture = DragGesture()
        d = ev.scenePos() - ev.lastScenePos()
        fine = bool(ev.modifiers() & Qt.ShiftModifier)
        res = self._pd_gesture.move(float(d.x()), -float(d.y()), fine)
        if res is not None:
            self.phase_dragged.emit(*res)
        if ev.isFinish():
            self._pd_gesture = None
            self.phase_drag_released.emit()

    # ---------- manual baseline ----------
    def set_baseline_mode(self, on: bool):
        self._bl_mode = on
        self.setCursor(Qt.PointingHandCursor if on else Qt.ArrowCursor)

    def _add_baseline_anchor(self, x: float, y: float):
        # the anchor's fill follows the theme (a white dot on a dark plot
        # read as a stray data point); the baseline-coloured ring is the mark
        t = pg.TargetItem(pos=(x, y), size=11, movable=True,
                          pen=pg.mkPen(theme.active().baseline, width=1.5),
                          brush=pg.mkBrush(*theme._rgb(theme.active().base), 220))
        t.sigPositionChanged.connect(lambda *_: self._update_baseline_curve())
        self.addItem(t)
        self._bl_anchors.append(t)
        self._update_baseline_curve()

    def baseline_anchors(self) -> list[tuple[float, float]]:
        return sorted(((float(t.pos().x()), float(t.pos().y()))
                       for t in self._bl_anchors), key=lambda a: a[0])

    def clear_baseline(self):
        for t in self._bl_anchors:
            self.removeItem(t)
        self._bl_anchors.clear()
        self._bl_curve.setData([], [])

    def baseline_curve(self, x: np.ndarray) -> np.ndarray | None:
        """Evaluate the anchor spline on x (PCHIP, edge-constant outside)."""
        pts = self.baseline_anchors()
        if len(pts) < 2:
            return None
        ax = np.array([p[0] for p in pts])
        ay = np.array([p[1] for p in pts])
        from scipy.interpolate import PchipInterpolator

        f = PchipInterpolator(ax, ay, extrapolate=False)
        out = f(x)
        out = np.where(np.isnan(out) & (x < ax[0]), ay[0], out)
        out = np.where(np.isnan(out) & (x > ax[-1]), ay[-1], out)
        return np.nan_to_num(out)

    def _update_baseline_curve(self):
        x = self._freq_x
        if x is None or not len(x):
            return
        y = self.baseline_curve(np.asarray(x))
        if y is None:
            self._bl_curve.setData([], [])
        else:
            self._bl_curve.setData(x, y)

    # ---------- fit zones ----------
    def set_zones(self, zones: list, on_change=None):
        """zones: list of [hi_ppm, lo_ppm]; draggable teal regions."""
        for r in self._zones:
            self.removeItem(r)
        self._zones.clear()
        zr, zg, zb = theme._rgb(theme.active().measure)
        for z in zones or []:
            region = pg.LinearRegionItem(values=(min(z), max(z)),
                                         brush=pg.mkBrush(zr, zg, zb, 26),
                                         hoverBrush=pg.mkBrush(zr, zg, zb, 45),
                                         pen=pg.mkPen(theme.active().measure, width=1))
            region.setZValue(-5)
            if on_change:
                region.sigRegionChangeFinished.connect(
                    lambda *_: on_change(self.zone_values()))
            self.addItem(region)
            region.setVisible(self._domain == "freq")
            self._zones.append(region)

    def zone_values(self) -> list:
        vals = []
        for r in self._zones:
            a, b = r.getRegion()
            vals.append([max(a, b), min(a, b)])
        return vals

    # ---------- literature shift-range overlay (assignment guide) ----------
    def set_ref_ranges(self, ranges: list[dict] | None, citation: str = "",
                       positions: list[dict] | None = None):
        """Draw (or clear, with None/[]) labeled, NON-interactive shaded
        δiso spans — the literature assignment guide (View ▸ Literature
        shift ranges). Each entry: {label, lo_ppm, hi_ppm, quad, note};
        the quadrupolar note (P_Q/C_Q — not a shift-axis quantity) goes
        into the label's second line and the hover tooltip."""
        for item in getattr(self, "_ref_items", []):
            self.removeItem(item)
        self._ref_items: list = []
        if not ranges and not positions:
            return
        self._draw_ref_positions(positions or [])
        if not ranges:
            return
        t = theme.active()
        rr, rg, rb = theme._rgb(t.pivot)
        for k, r in enumerate(ranges):
            region = pg.LinearRegionItem(
                values=(float(r["lo_ppm"]), float(r["hi_ppm"])),
                movable=False,
                brush=pg.mkBrush(rr, rg, rb, 18),
                pen=pg.mkPen(t.pivot, width=1, style=Qt.DotLine))
            region.setZValue(-20)                 # behind data, zones, model
            tip = f"{r['label']}: {r['lo_ppm']:g} … {r['hi_ppm']:g} ppm"
            if r.get("quad"):
                tip += f"\n{r['quad']}"
            if r.get("note"):
                tip += f"\n{r['note']}"
            # per-entry citation when the dataset carries one (different
            # nuclei come from different primary papers), else the fallback
            ref = r.get("_citation") or citation
            if ref:
                tip += f"\n[{ref}]"
            region.setToolTip(tip)
            label = pg.TextItem(
                r["label"] + (f"\n{r['quad'].split(';')[0]}" if r.get("quad")
                              else ""),
                color=t.pivot, anchor=(0.5, 0.0))
            label.setZValue(-19)
            # stagger label heights so adjacent ranges don't overwrite
            # each other; positions refresh with the view via ViewBox signal
            self.addItem(region)
            self.addItem(label)
            region.setVisible(self._domain == "freq")
            label.setVisible(self._domain == "freq")
            self._ref_items += [region, label]
            self._place_ref_label(label, r, k)
        vb = self.getPlotItem().getViewBox()
        vb.sigRangeChanged.connect(self._refresh_ref_labels)
        self._ref_meta = list(ranges)

    def _draw_ref_positions(self, positions: list[dict]):
        """Reported single positions (e.g. crystalline fluorides) as dotted
        ticks with a compound label, cleared together with the ranges."""
        t = theme.active()
        for k, r in enumerate(positions):
            line = pg.InfiniteLine(
                pos=float(r["ppm"]), angle=90, movable=False,
                pen=pg.mkPen(t.text_dim, width=1, style=Qt.DotLine))
            line.setZValue(-18)
            tip = f"{r['label']}: {r['ppm']:g} ppm"
            if r.get("note"):
                tip += f"\n{r['note']}"
            if r.get("_citation"):
                tip += f"\n{r['_citation']}"
            line.setToolTip(tip)
            line.label = pg.InfLineLabel(
                line, r["label"], position=0.97 - 0.04 * (k % 4),
                color=t.text_dim, movable=False)
            self.addItem(line)
            line.setVisible(self._domain == "freq")
            self._ref_items.append(line)

    def _place_ref_label(self, label, r: dict, k: int):
        vb = self.getPlotItem().getViewBox()
        (_x0, _x1), (y0, y1) = vb.viewRange()
        frac = 0.97 - 0.07 * (k % 3)              # 3-step stagger
        label.setPos((float(r["lo_ppm"]) + float(r["hi_ppm"])) / 2.0,
                     y0 + frac * (y1 - y0))

    def _refresh_ref_labels(self, *_):
        items = getattr(self, "_ref_items", [])
        meta = getattr(self, "_ref_meta", [])
        labels = [it for it in items if isinstance(it, pg.TextItem)]
        for k, (label, r) in enumerate(zip(labels, meta)):
            self._place_ref_label(label, r, k)

    # ---------- onboarding placeholder (empty canvas) ----------
    def set_placeholder(self, text: str | None):
        """Show faint centred guidance on an empty canvas, or clear it (None).
        Mouse-transparent, so drag-and-drop and clicks pass straight through."""
        if not text:
            if getattr(self, "_placeholder", None) is not None:
                self._placeholder.hide()
            return
        if getattr(self, "_placeholder", None) is None:
            from PySide6.QtWidgets import QLabel
            lbl = QLabel(self)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setWordWrap(True)
            lbl.setAttribute(Qt.WA_TransparentForMouseEvents)
            self._placeholder = lbl
        self._placeholder.setText(text)
        self._placeholder.setStyleSheet(
            f"color: {theme.active().text_dim}; font-size: 13px; "
            "background: transparent;")
        self._position_placeholder()
        self._placeholder.show()
        self._placeholder.raise_()

    def _position_placeholder(self):
        lbl = getattr(self, "_placeholder", None)
        if lbl is None:
            return
        w = min(max(self.width() - 40, 200), 480)
        h = 130
        lbl.setGeometry((self.width() - w) // 2, (self.height() - h) // 2, w, h)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self._position_placeholder()

    # ---------- data ----------
    def set_experiment(self, x: np.ndarray, y: np.ndarray):
        """A REAL frequency-domain spectrum arrived: remember it, leave the
        FID view if it was on, draw it, and tell the workbench
        (experiment_set) so a stale FID / imaginary display is dropped."""
        self._freq_x, self._freq_y = x, y
        if self._domain == "time":
            self._leave_time_domain()
        if x is not None and len(x):
            self.set_placeholder(None)          # real data arrived — hide the hint
            if self._phase_mode:
                # data arriving while drag-to-phase is armed: the dock is
                # already visible, so nothing else would create the pivot
                self.show_phase_pivot(True)     # no-op once it exists
        self._exp.setData(x, y)
        self.set_trace_label("experiment")
        self.experiment_set.emit()

    def set_channel_trace(self, y: np.ndarray, label: str):
        """Draw another channel (imag, |S|) of the SAME spectrum on the ppm
        axis. The remembered real trace is untouched (pivot default and
        calibrate snapping keep using it) and experiment_set is not emitted."""
        if self._domain == "time":
            self._leave_time_domain()
        self._exp.setData(self._freq_x, y)
        self.set_trace_label(label)

    def set_fid(self, t_ms: np.ndarray, y: np.ndarray, label: str):
        """Show a time-domain trace (ms, not inverted) instead of the
        spectrum; every ppm item is hidden until set_experiment /
        set_channel_trace brings the frequency axis back."""
        pi = self.getPlotItem()
        t_ms = np.asarray(t_ms, float)
        aq = float(t_ms[-1]) if t_ms.size else 0.0
        entering = self._domain == "freq"
        if entering:
            self._freq_ranges = pi.getViewBox().viewRange()
            self._domain = "time"
            pi.invertX(False)
            ax = pi.getAxis("bottom")
            if isinstance(ax, ScaledAxis):
                ax.set_factor(1.0)
            th = theme.active()
            self._apply_axis_label({"color": th.axis, "font-size": "10pt"})
            self._set_freq_items_visible(False)
            self._arm_phase_drag()          # no phasing of a FID: left drag pans
        self._exp.setData(t_ms, np.asarray(y, float))
        self.set_trace_label(label)
        if entering or self._fid_aq != aq:      # new AQ (zf / TDeff): refit
            pi.enableAutoRange()
        self._fid_aq = aq

    def _leave_time_domain(self):
        """Frequency axis back: ppm items visible, unit/label/factor and the
        pre-toggle zoom restored."""
        pi = self.getPlotItem()
        self._domain = "freq"
        self._fid_aq = None
        pi.invertX(True)
        self.set_axis_unit(self._axis_unit, self._axis_sfo_MHz)
        self._set_freq_items_visible(True)
        self._arm_phase_drag()              # the gesture is back with the ppm axis
        if self._freq_ranges is not None:
            (x0, x1), (y0, y1) = self._freq_ranges
            pi.setXRange(min(x0, x1), max(x0, x1), padding=0)
            pi.setYRange(min(y0, y1), max(y0, y1), padding=0)
            self._freq_ranges = None

    def _set_freq_items_visible(self, on: bool):
        """Hide / show everything that lives in ppm (the experiment curve is
        the one item that switches meaning)."""
        for it in (self._model, self._resid, self._bl_curve):
            it.setVisible(on)
        groups = (self._components, self._paddles, self._markers,
                  self._overlay_items, self._zones, self._measure_lines,
                  self._bl_anchors, self._comp_labels,
                  getattr(self, "_ref_items", []))
        for group in groups:
            for it in group:
                it.setVisible(on)
        if self._pivot is not None:
            self._pivot.setVisible(on)
        if self._hover_label is not None and not on:
            self._hover_label.setVisible(False)
        if on:
            # the zero line belongs to a drawn residual strip only
            rx = self._resid.xData
            self._resid_zero.setVisible(rx is not None and len(rx) > 0)
        else:
            self._resid_zero.setVisible(False)
            self.stop_fit_animation()

    def set_title(self, text: str):
        self.getPlotItem().setTitle(
            f"<span style='color:{theme.active().text}; font-size:11pt'>{text}</span>"
            if text else None)

    def set_model(self, x, total, per_site, labels, hidden: set[int],
                  exp_x=None, exp_y=None):
        if x is None:
            self._model.setData([], [])
            self._resid.setData([], [])
            for c in self._components:
                self.removeItem(c)
            self._components.clear()
            self._comp_data = []
            self._refresh_comp_labels()
            return
        # the model is simulated on the (wider) kernel axis: draw it only
        # across the experiment (+ a small margin) so nothing model-side ever
        # reaches past the data -- the view range is the data's and the user's
        mask = self._model_mask(x, exp_x)
        if mask is not None and mask.sum() < 2:
            mask = None                      # nothing of the model in range
        x, total = self._masked(mask, x, total)
        per_site = list(self._masked(mask, *per_site)) if per_site else []
        self._model.setData(x, total)

        # residual, offset below zero as a dedicated strip
        if self.show_residual and exp_x is not None and len(exp_x):
            yi = np.interp(exp_x, x, total)
            offset = -0.10 * float(np.max(exp_y)) if len(exp_y) else 0.0
            self._resid.setData(exp_x, (exp_y - yi) + offset)
            self._resid_zero.setPos(offset)
            self._resid_zero.setVisible(True)
        else:
            self._resid_zero.setVisible(False)
            self._resid.setData([], [])

        # components: reuse items, add/remove as needed (created hidden while
        # the FID is displayed -- every simulation allocates new ones)
        while len(self._components) < len(per_site):
            item = self._add_model_curve()
            self._tune_curve(item)
            item.setVisible(self._domain == "freq")
            self._components.append(item)
        while len(self._components) > len(per_site):
            self.removeItem(self._components.pop())
        self._comp_data = []
        xa = np.asarray(x, float)
        order = np.argsort(xa)
        for i, ys in enumerate(per_site):
            item = self._components[i]
            if self.show_components and i not in hidden:
                col = pg.mkColor(site_color(i))
                item.setPen(pg.mkPen(col, width=1.3, style=Qt.DashLine))
                fill = pg.mkColor(col); fill.setAlpha(28)
                item.setData(x, ys, fillLevel=0.0, fillBrush=pg.mkBrush(fill))
                name = labels[i] if i < len(labels) and labels[i] else ""
                text = cellparse.index_to_letter(i) + (f" \u00b7 {name}" if name else "")
                ya = np.asarray(ys, float)
                if ya.shape == xa.shape and ya.size:
                    self._comp_data.append((i, xa[order], ya[order], text))
            else:
                item.setData([], [])
        self._refresh_comp_labels()

    # ---------- component labels ----------
    def set_show_labels(self, on: bool):
        """Pin every visible component's letter and name at its maximum
        (View > Component labels). Off, the name still appears for the
        component under the cursor."""
        self.show_labels = bool(on)
        if self.show_labels and self._hover_label is not None:
            self._hover_label.setVisible(False)
        self._refresh_comp_labels()

    def _refresh_comp_labels(self):
        for it in self._comp_labels:
            self.removeItem(it)
        self._comp_labels.clear()
        if not self.show_labels:
            return
        for i, cx, cy, text in self._comp_data:
            j = int(np.argmax(cy))
            if not np.isfinite(cy[j]) or cy[j] <= 0:
                continue
            lab = pg.TextItem(text, color=site_color(i), anchor=(0.5, 1.0))
            lab.setZValue(30)
            self.addItem(lab)
            lab.setPos(float(cx[j]), float(cy[j]))
            lab.setVisible(self._domain == "freq")
            self._comp_labels.append(lab)

    def _hover_component(self, x: float):
        """Name the tallest component under the cursor (labels not pinned)."""
        best = None
        for i, cx, cy, text in self._comp_data:
            if x < cx[0] or x > cx[-1]:
                continue
            yi = float(np.interp(x, cx, cy))
            top = float(np.max(cy))
            if top <= 0 or yi < 0.05 * top:
                continue
            if best is None or yi > best[1]:
                best = (i, yi, text)
        if best is None:
            if self._hover_label is not None:
                self._hover_label.setVisible(False)
            return
        i, yi, text = best
        if self._hover_label is None:
            self._hover_label = pg.TextItem("", anchor=(0.5, 1.0))
            self._hover_label.setZValue(31)
            self.addItem(self._hover_label)
        self._hover_label.setColor(pg.mkColor(site_color(i)))
        self._hover_label.setText(text)
        self._hover_label.setPos(x, yi)
        self._hover_label.setVisible(True)

    # ---------- comparison overlays ----------
    def set_overlays(self, overlays: list[tuple]):
        """overlays: [(x, y, color, label), ...] drawn behind the active data."""
        for it in self._overlay_items:
            self.removeItem(it)
        self._overlay_items.clear()
        for x, y, color, label in overlays:
            item = self.plot(x, y, pen=pg.mkPen(color, width=1.1),
                             name=label, antialias=True)
            self._tune_curve(item)
            item.setZValue(-10)
            item.setVisible(self._domain == "freq" and not self._overlays_hidden)
            self._overlay_items.append(item)

    def set_overlays_hidden(self, hidden: bool):
        """View > Overlays: hide or show every compared spectrum at once
        without removing it from the Datasets dock."""
        self._overlays_hidden = bool(hidden)
        for it in self._overlay_items:
            it.setVisible(self._domain == "freq" and not self._overlays_hidden)

    # ---------- markers (legacy InfiniteLine API kept for tests) ----------
    def set_markers(self, positions: list[tuple[int, float, bool]]):
        """positions: [(site_index, ppm, draggable), ...] for visible sites."""
        for m in self._markers:
            self.removeItem(m)
        self._markers.clear()
        for idx, ppm, draggable in positions:
            line = pg.InfiniteLine(
                pos=ppm, angle=90, movable=draggable,
                pen=pg.mkPen(site_color(idx), width=1.3, style=Qt.DashLine),
                hoverPen=pg.mkPen(site_color(idx), width=2.5),
            )
            line.site_index = idx
            if draggable:
                line.sigPositionChangeFinished.connect(self._marker_done)
            self.getPlotItem().addItem(line, ignoreBounds=True)
            line.setVisible(self._domain == "freq")
            self._markers.append(line)

    def _marker_done(self, line):
        self.marker_moved.emit(line.site_index, float(line.value()))

    # ---------- dmfit-style paddles ----------
    def set_paddles(self, states: list[tuple[int, float, float, float, bool]]):
        """states: [(site_index, pos_ppm, amp, fwhm_ppm, movable), ...]."""
        from larmor.desktop.paddle import Paddle

        for p in self._paddles:
            self.removeItem(p)
        self._paddles.clear()
        for idx, pos, amp, fwhm, movable in states:
            pad = Paddle(idx, site_color(idx), pos, amp, fwhm, movable)
            pad.moved.connect(self.paddle_moved)
            pad.released.connect(self.paddle_released)
            # a paddle is a model handle: it must not pull the auto range
            # (a linked sideband copy can sit outside the data)
            self.getPlotItem().addItem(pad, ignoreBounds=True)
            pad.setVisible(self._domain == "freq")
            self._paddles.append(pad)

    def show_paddles(self, on: bool):
        for p in self._paddles:
            p.setVisible(on and self._domain == "freq")

    def current_xrange(self) -> tuple[float, float]:
        """(hi, lo) of the displayed frequency window in ppm -- while the FID
        is shown, the window saved on entry (callers hand it to dialogs as a
        ppm range)."""
        if self._domain == "time":
            if self._freq_ranges is not None:
                x0, x1 = self._freq_ranges[0]
            elif self._freq_x is not None and len(self._freq_x):
                x0, x1 = float(np.min(self._freq_x)), float(np.max(self._freq_x))
            else:
                x0 = x1 = 0.0
            return (max(x0, x1), min(x0, x1))
        (x0, x1), _ = self.getPlotItem().getViewBox().viewRange()
        return (max(x0, x1), min(x0, x1))
