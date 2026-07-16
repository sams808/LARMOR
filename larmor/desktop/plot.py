"""Central spectrum view: pyqtgraph plot with draggable site markers and
click-to-add. All rendering is direct numpy -> GPU-backed canvas: no browser,
no serialization, instant interaction even on 32k-point spectra."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal

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

    def __init__(self, parent=None):
        super().__init__(parent, background="w")
        self.getPlotItem().invertX(True)              # ppm convention
        self.setLabel("bottom", "shift / ppm")
        self.setLabel("left", "intensity")
        self.showGrid(x=True, y=True, alpha=0.15)
        self.getPlotItem().getViewBox().setMouseMode(pg.ViewBox.PanMode)
        leg = self.addLegend(offset=(8, 8), labelTextSize="8pt")
        leg.setBrush(pg.mkBrush(255, 255, 255, 200))

        self._exp = self.plot([], [], pen=pg.mkPen("#16202a", width=1.2),
                              name="experiment")
        self._model = self.plot([], [], pen=pg.mkPen("#d62728", width=1.6),
                                name="model")
        self._resid = self.plot([], [], pen=pg.mkPen("#8a969e", width=0.9))
        self._components: list[pg.PlotDataItem] = []
        self._markers: list[pg.InfiniteLine] = []
        self._add_mode: str | None = None
        self.show_components = True
        self.show_residual = True

        self.scene().sigMouseClicked.connect(self._on_click)
        self.scene().sigMouseMoved.connect(self._on_move)
        self._paddles: list = []

        # manual baseline: draggable anchors + live PCHIP preview
        self._bl_mode = False
        self._bl_anchors: list[pg.TargetItem] = []
        self._bl_curve = self.plot([], [], pen=pg.mkPen("#c88a1e", width=1.4,
                                                        style=Qt.DashLine))
        # dmfit-style fit zones
        self._zones: list[pg.LinearRegionItem] = []

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
        if self._bl_mode:
            self._add_baseline_anchor(float(p.x()), float(p.y()))
            ev.accept()
            return
        if self._add_mode is not None:
            self.add_requested.emit(float(p.x()), abs(float(p.y())))
            ev.accept()

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

        # residual, offset below zero
        if self.show_residual and exp_x is not None and len(exp_x):
            yi = np.interp(exp_x, x, total)
            offset = -0.08 * float(np.max(exp_y)) if len(exp_y) else 0.0
            self._resid.setData(exp_x, (exp_y - yi) + offset)
        else:
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
                item.setPen(pg.mkPen(site_color(i), width=1.0,
                                     style=Qt.DashLine))
                item.setData(x, ys)
            else:
                item.setData([], [])

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
