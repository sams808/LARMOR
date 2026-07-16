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
        if self._add_mode is None or ev.button() != Qt.LeftButton:
            return
        vb = self.getPlotItem().getViewBox()
        if not self.sceneBoundingRect().contains(ev.scenePos()):
            return
        p = vb.mapSceneToView(ev.scenePos())
        self.add_requested.emit(float(p.x()), abs(float(p.y())))
        ev.accept()

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
