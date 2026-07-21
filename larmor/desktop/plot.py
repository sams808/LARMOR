"""Central spectrum view: pyqtgraph plot with draggable site markers and
click-to-add. All rendering is direct numpy -> GPU-backed canvas: no browser,
no serialization, instant interaction even on 32k-point spectra."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont

SITE_COLORS = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
               "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", "#ff7f0e"]


def site_color(i: int) -> str:
    return SITE_COLORS[i % len(SITE_COLORS)]


class SpectrumView(pg.PlotWidget):
    """Experiment + model + components + residual, with dmfit-style paddles."""

    add_requested = Signal(float, float)      # (ppm, amplitude) from a click
    marker_moved = Signal(int, float)         # legacy: (site index, new ppm)
    paddle_moved = Signal(int, float, float, float)   # index, pos, amp, fwhm
    paddle_released = Signal(int)
    cursor_moved = Signal(float, float)       # live x/y for the status bar
    file_dropped = Signal(str)                # a data file dragged onto the plot
    calibrate_picked = Signal(float)          # snapped peak ppm to reference
    measure_changed = Signal(float, float)    # two ppm cursors (ruler)

    def __init__(self, parent=None):
        super().__init__(parent, background="#fcfdfc")
        pi = self.getPlotItem()
        pi.invertX(True)                              # ppm convention
        # professional axis styling: dark ink axes, tick labels, light grid
        axis_pen = pg.mkPen("#37424a", width=1.2)
        tick_font = QFont()
        tick_font.setPointSize(9)
        for name in ("bottom", "left"):
            ax = pi.getAxis(name)
            ax.setPen(axis_pen)
            ax.setTextPen(pg.mkPen("#37424a"))
            ax.setStyle(tickFont=tick_font, tickLength=-5)
        pi.getAxis("top").setPen(pg.mkPen("#c5ccc6"))
        pi.getAxis("right").setPen(pg.mkPen("#c5ccc6"))
        pi.showAxis("top"); pi.getAxis("top").setStyle(showValues=False)
        pi.showAxis("right"); pi.getAxis("right").setStyle(showValues=False)
        label_style = {"color": "#37424a", "font-size": "10pt"}
        self.setLabel("bottom", "chemical shift", units="ppm", **label_style)
        self.setLabel("left", "intensity", **label_style)
        self.showGrid(x=True, y=True, alpha=0.08)
        pi.getViewBox().setMouseMode(pg.ViewBox.PanMode)
        pi.setContentsMargins(6, 10, 6, 6)

        leg = self.addLegend(offset=(12, 10), labelTextSize="9pt",
                             brush=pg.mkBrush(255, 255, 255, 235),
                             pen=pg.mkPen("#d7dcd9"), labelTextColor="#16202a")
        self._legend = leg

        self._exp = self.plot([], [], pen=pg.mkPen("#1a2831", width=1.4),
                              name="experiment", antialias=True)
        self._model = self.plot([], [], pen=pg.mkPen("#c0392b", width=1.8),
                                name="model", antialias=True)
        self._resid = self.plot([], [], pen=pg.mkPen("#9aa5ab", width=1.0),
                                name="residual", antialias=True)
        # faint zero line for the offset residual strip
        self._resid_zero = pg.InfiniteLine(
            angle=0, pen=pg.mkPen("#d0d6d1", width=1, style=Qt.DotLine))
        self._resid_zero.setVisible(False)
        self.addItem(self._resid_zero)
        pi.getAxis("left").setStyle(tickTextWidth=48, autoExpandTextSpace=False)
        self._components: list[pg.PlotDataItem] = []
        self._markers: list[pg.InfiniteLine] = []
        self._add_mode: str | None = None
        self.show_components = True
        self.show_residual = True

        self.scene().sigMouseClicked.connect(self._on_click)
        self.scene().sigMouseMoved.connect(self._on_move)
        self._paddles: list = []
        self.setAcceptDrops(True)                # drag a spectrum onto the plot

        # manual baseline: draggable anchors + live PCHIP preview
        self._bl_mode = False
        self._bl_anchors: list[pg.TargetItem] = []
        self._bl_curve = self.plot([], [], pen=pg.mkPen("#c88a1e", width=1.4,
                                                        style=Qt.DashLine))
        # dmfit-style fit zones
        self._zones: list[pg.LinearRegionItem] = []
        # calibrate (click a peak) and measure (two-cursor ruler)
        self._cal_mode = False
        self._measure_lines: list[pg.InfiniteLine] = []
        # comparison overlays (other datasets drawn behind the active one)
        self._overlay_items: list[pg.PlotDataItem] = []

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
            self.file_dropped.emit(urls[0].toLocalFile())
            ev.acceptProposedAction()
        else:
            super().dropEvent(ev)

    def _on_move(self, scene_pos):
        vb = self.getPlotItem().getViewBox()
        if self.sceneBoundingRect().contains(scene_pos):
            p = vb.mapSceneToView(scene_pos)
            self.cursor_moved.emit(float(p.x()), float(p.y()))

    # ---------- add mode ----------
    def set_add_mode(self, model_name: str | None):
        self._add_mode = model_name
        self.setCursor(Qt.CrossCursor if model_name else Qt.ArrowCursor)

    def _on_click(self, ev):
        if ev.button() != Qt.LeftButton:
            return
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
        """Snap a click to the nearest local maximum of the experiment."""
        x, y = self._exp.xData, self._exp.yData
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
        (x0, x1), _ = self.getPlotItem().getViewBox().viewRange()
        lo, hi = min(x0, x1), max(x0, x1)
        for frac in (0.35, 0.65):
            ln = pg.InfiniteLine(pos=lo + frac * (hi - lo), angle=90,
                                 movable=True,
                                 pen=pg.mkPen("#0e7c86", width=1.4),
                                 hoverPen=pg.mkPen("#0e7c86", width=2.4),
                                 label="{value:.1f}",
                                 labelOpts={"color": "#0e7c86",
                                            "position": 0.92})
            ln.sigPositionChanged.connect(self._emit_measure)
            self.addItem(ln)
            self._measure_lines.append(ln)
        self._emit_measure()

    def _emit_measure(self, *_):
        if len(self._measure_lines) == 2:
            self.measure_changed.emit(float(self._measure_lines[0].value()),
                                      float(self._measure_lines[1].value()))

    # ---------- manual baseline ----------
    def set_baseline_mode(self, on: bool):
        self._bl_mode = on
        self.setCursor(Qt.PointingHandCursor if on else Qt.ArrowCursor)

    def _add_baseline_anchor(self, x: float, y: float):
        t = pg.TargetItem(pos=(x, y), size=11, movable=True,
                          pen=pg.mkPen("#c88a1e", width=1.5),
                          brush=pg.mkBrush(255, 255, 255, 220))
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
        x = self._exp.xData
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
        for z in zones or []:
            region = pg.LinearRegionItem(values=(min(z), max(z)),
                                         brush=pg.mkBrush(14, 124, 134, 26),
                                         hoverBrush=pg.mkBrush(14, 124, 134, 45),
                                         pen=pg.mkPen("#0e7c86", width=1))
            region.setZValue(-5)
            if on_change:
                region.sigRegionChangeFinished.connect(
                    lambda *_: on_change(self.zone_values()))
            self.addItem(region)
            self._zones.append(region)

    def zone_values(self) -> list:
        vals = []
        for r in self._zones:
            a, b = r.getRegion()
            vals.append([max(a, b), min(a, b)])
        return vals

    # ---------- data ----------
    def set_experiment(self, x: np.ndarray, y: np.ndarray):
        self._exp.setData(x, y)

    def set_title(self, text: str):
        self.getPlotItem().setTitle(
            f"<span style='color:#37424a; font-size:11pt'>{text}</span>"
            if text else None)

    def set_model(self, x, total, per_site, labels, hidden: set[int],
                  exp_x=None, exp_y=None):
        if x is None:
            self._model.setData([], [])
            self._resid.setData([], [])
            for c in self._components:
                self.removeItem(c)
            self._components.clear()
            return
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

        # components: reuse items, add/remove as needed
        while len(self._components) < len(per_site):
            item = self.plot([], [])
            self._components.append(item)
        while len(self._components) > len(per_site):
            self.removeItem(self._components.pop())
        for i, ys in enumerate(per_site):
            item = self._components[i]
            if self.show_components and i not in hidden:
                col = pg.mkColor(site_color(i))
                item.setPen(pg.mkPen(col, width=1.3, style=Qt.DashLine))
                fill = pg.mkColor(col); fill.setAlpha(28)
                item.setData(x, ys, fillLevel=0.0, fillBrush=pg.mkBrush(fill))
            else:
                item.setData([], [])

    # ---------- comparison overlays ----------
    def set_overlays(self, overlays: list[tuple]):
        """overlays: [(x, y, color, label), ...] drawn behind the active data."""
        for it in self._overlay_items:
            self.removeItem(it)
        self._overlay_items.clear()
        for x, y, color, label in overlays:
            item = self.plot(x, y, pen=pg.mkPen(color, width=1.1),
                             name=label, antialias=True)
            item.setZValue(-10)
            self._overlay_items.append(item)

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
            self.addItem(line)
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
            self.addItem(pad)
            self._paddles.append(pad)

    def show_paddles(self, on: bool):
        for p in self._paddles:
            p.setVisible(on)

    def current_xrange(self) -> tuple[float, float]:
        (x0, x1), _ = self.getPlotItem().getViewBox().viewRange()
        return (max(x0, x1), min(x0, x1))
