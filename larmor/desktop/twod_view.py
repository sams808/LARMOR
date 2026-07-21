"""Embeddable 2D contour view for the main window's central stack.

Shows any 2D dataset (real MQMAS or a pseudo-2D relaxation array) as a
log-contour map with F1/F2 projections and a diagonal, and lets the user pull
a 1D trace out of it for fitting: the F2 skyline/sum projection, or a single
row picked with a draggable horizontal cursor. That 1D trace is emitted to the
fit workbench, so a 2D opens, displays, and THEN the user chooses what to do.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox, QDoubleSpinBox, QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
    QWidget,
)


class Contour2DView(QWidget):
    #: (ppm, amp, label) of a 1D trace the user wants to fit
    slice_to_fit = Signal(object, object, str)

    def __init__(self):
        super().__init__()
        self.data = None
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)

        bar = QHBoxLayout()
        self.title = QLabel("2D dataset")
        self.title.setStyleSheet("font-weight: 600; color: #16202a;")
        bar.addWidget(self.title, 1)
        bar.addWidget(QLabel("levels"))
        self.nlevels = QDoubleSpinBox()
        self.nlevels.setRange(4, 40); self.nlevels.setValue(12)
        self.nlevels.valueChanged.connect(self._redraw)
        bar.addWidget(self.nlevels)
        bar.addWidget(QLabel("floor ×σ"))
        self.floor = QDoubleSpinBox()
        self.floor.setRange(1.0, 40.0); self.floor.setValue(8.0)
        self.floor.valueChanged.connect(self._redraw)
        bar.addWidget(self.floor)
        v.addLayout(bar)

        # contour with top (F2) and left (F1) projections
        self.glw = pg.GraphicsLayoutWidget()
        self.glw.setBackground("#fcfdfc")
        self.p_top = self.glw.addPlot(row=0, col=1)
        self.p_top.setMaximumHeight(80); self.p_top.hideAxis("bottom")
        self.p_main = self.glw.addPlot(row=1, col=1)
        self.p_left = self.glw.addPlot(row=1, col=0)
        self.p_left.setMaximumWidth(80); self.p_left.hideAxis("left")
        self.p_main.setLabel("bottom", "F2 (ppm)")
        self.p_main.setLabel("right", "F1")
        self.p_main.getViewBox().invertX(True)
        self.p_main.getViewBox().invertY(True)
        self.p_top.setXLink(self.p_main)
        self.p_left.setYLink(self.p_main)
        v.addWidget(self.glw, 1)
        self._slice_line = None

        pick = QHBoxLayout()
        pick.addWidget(QLabel("Send a 1D trace to fitting:"))
        self.btnSky = QPushButton("F2 skyline →")
        self.btnSky.clicked.connect(lambda: self._emit_projection("skyline"))
        self.btnSum = QPushButton("F2 sum →")
        self.btnSum.clicked.connect(lambda: self._emit_projection("sum"))
        self.btnRow = QPushButton("row at cursor →")
        self.btnRow.setToolTip("drag the dashed horizontal line to a row, then "
                               "send that F2 slice to the fit workbench")
        self.btnRow.clicked.connect(self._emit_row)
        pick.addWidget(self.btnSky)
        pick.addWidget(self.btnSum)
        pick.addWidget(self.btnRow)
        pick.addStretch(1)
        self.hint = QLabel("")
        self.hint.setStyleSheet("color: #5a6871;")
        pick.addWidget(self.hint)
        v.addLayout(pick)

    # ------------------------------------------------------------------
    def set_data(self, data, title: str = ""):
        self.data = data.normalized() if hasattr(data, "normalized") else data
        self.title.setText(title or "2D dataset")
        f1_kind = "arrayed (relaxation)" if getattr(data, "notes", None) and \
            any("pseudo" in n or "arrayed" in n for n in data.notes) else "F1"
        self.hint.setText(f"F1 = {f1_kind}")
        self._redraw()

    def _redraw(self):
        if self.data is None:
            return
        d = self.data
        for p in (self.p_main, self.p_top, self.p_left):
            p.clear()
        z = d.z
        edge = max(1, min(z.shape) // 20)
        frame = np.concatenate([z[:edge].ravel(), z[-edge:].ravel(),
                                z[:, :edge].ravel(), z[:, -edge:].ravel()])
        floor = max(self.floor.value() * float(frame.std() or 1e-6), 0.02)
        n = int(self.nlevels.value())
        top = float(np.nanmax(z)) or 1.0
        levels = np.logspace(np.log10(floor), np.log10(top), n)
        cmap = pg.colormap.get("viridis")
        colors = cmap.getLookupTable(0.0, 1.0, n)

        f2, f1 = d.f2_ppm, d.f1_ppm
        tr = pg.QtGui.QTransform()
        tr.translate(f2[0], f1[0])
        tr.scale((f2[-1] - f2[0]) / max(z.shape[1] - 1, 1),
                 (f1[-1] - f1[0]) / max(z.shape[0] - 1, 1))
        for lvl, col in zip(levels, colors):
            iso = pg.IsocurveItem(data=z.T, level=lvl,
                                  pen=pg.mkPen(tuple(int(c) for c in col), width=1))
            iso.setTransform(tr)
            self.p_main.addItem(iso)

        lo = max(min(f2), min(f1)); hi = min(max(f2), max(f1))
        self.p_main.plot([lo, hi], [lo, hi],
                         pen=pg.mkPen("#b9c1bc", style=Qt.DashLine))
        self.p_top.plot(f2, d.z.max(axis=0), pen=pg.mkPen("#0e7c86"))
        self.p_left.plot(d.z.max(axis=1), f1, pen=pg.mkPen("#0e7c86"))

        self._slice_line = pg.InfiniteLine(
            pos=float(f1[len(f1) // 2]), angle=0, movable=True,
            pen=pg.mkPen("#d62728", width=1.3, style=Qt.DashLine))
        self.p_main.addItem(self._slice_line)

    # ------------------------------------------------------------------
    def _emit_projection(self, mode: str):
        if self.data is None:
            return
        d = self.data
        proj = (d.z.max(axis=0) if mode == "skyline" else d.z.sum(axis=0))
        self.slice_to_fit.emit(np.asarray(d.f2_ppm), np.asarray(proj),
                               f"F2 {mode} of {self.title.text()}")

    def _emit_row(self):
        if self.data is None or self._slice_line is None:
            return
        d = self.data
        y = float(self._slice_line.value())
        i = int(np.argmin(np.abs(d.f1_ppm - y)))
        self.slice_to_fit.emit(np.asarray(d.f2_ppm), np.asarray(d.z[i]),
                               f"F1={d.f1_ppm[i]:.1f} slice of {self.title.text()}")
