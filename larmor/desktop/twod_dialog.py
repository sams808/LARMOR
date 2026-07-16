"""2D MQMAS viewer and fit dialog: contour display with projections."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
)


class TwoDDialog(QDialog):
    def __init__(self, parent, expno: str | None):
        super().__init__(parent)
        self.setWindowTitle("2D MQMAS")
        self.resize(900, 680)
        self.expno = expno
        self.data = None
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.lbl = QLabel(expno or "no 2D EXPNO selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Open 2D EXPNO…")
        btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(btn)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("levels"))
        self.nlevels = QDoubleSpinBox()
        self.nlevels.setRange(4, 40); self.nlevels.setValue(12)
        self.nlevels.valueChanged.connect(self._redraw)
        opts.addWidget(self.nlevels)
        opts.addWidget(QLabel("floor ×σ"))
        self.floor = QDoubleSpinBox()
        self.floor.setRange(1.0, 30.0); self.floor.setValue(8.0)
        self.floor.valueChanged.connect(self._redraw)
        opts.addWidget(self.floor)
        opts.addWidget(QLabel("shear"))
        self.shear = QDoubleSpinBox()
        self.shear.setRange(-3.0, 3.0); self.shear.setDecimals(3)
        self.shear.setToolTip("manual shear factor (mrsimulator MQMAS data is "
                              "already sheared; 0 = none)")
        opts.addWidget(self.shear)
        btnShear = QPushButton("Apply shear")
        btnShear.clicked.connect(self._apply_shear)
        opts.addWidget(btnShear)
        opts.addStretch(1)
        v.addLayout(opts)

        # contour with left (F1) and top (F2) projections
        self.glw = pg.GraphicsLayoutWidget()
        self.glw.setBackground("w")
        self.p_top = self.glw.addPlot(row=0, col=1)
        self.p_top.setMaximumHeight(90)
        self.p_top.hideAxis("bottom")
        self.p_main = self.glw.addPlot(row=1, col=1)
        self.p_left = self.glw.addPlot(row=1, col=0)
        self.p_left.setMaximumWidth(90)
        self.p_left.hideAxis("left")
        self.p_main.setLabel("bottom", "F2 (ppm)")
        self.p_main.setLabel("right", "F1 (ppm)")
        self.p_main.getViewBox().invertX(True)
        self.p_main.getViewBox().invertY(True)
        self.p_top.setXLink(self.p_main)
        self.p_left.setYLink(self.p_main)
        v.addWidget(self.glw, 1)

        self.res = QLabel("open a processed 2D dataset (2rr)")
        self.res.setStyleSheet("color: #37424a;")
        v.addWidget(self.res)

        if expno:
            self._load(expno)

    def _pick(self):
        p = QFileDialog.getExistingDirectory(self, "Bruker 2D EXPNO")
        if p:
            self._load(p)

    def _load(self, expno):
        from larmor import twod

        try:
            self.data = twod.read_bruker_2d(expno).normalized()
        except Exception as exc:
            self.res.setText(f"failed to read 2D: {exc}")
            return
        self.expno = expno
        self.lbl.setText(expno)
        self.res.setText(
            f"{self.data.nucleus or '?'} · "
            f"{self.data.z.shape[0]}×{self.data.z.shape[1]} points")
        self._redraw()

    def _apply_shear(self):
        from larmor import twod

        if self.data is None or self.shear.value() == 0:
            return
        self.data = twod.shear(self.data, self.shear.value())
        self._redraw()
        self.res.setText(f"sheared by {self.shear.value():g}")

    def _redraw(self):
        if self.data is None:
            return
        d = self.data
        self.p_main.clear()
        self.p_top.clear()
        self.p_left.clear()

        # contour lines via IsocurveItem at log-spaced levels
        z = d.z
        edge = max(1, min(z.shape) // 20)
        frame = np.concatenate([z[:edge].ravel(), z[-edge:].ravel(),
                                z[:, :edge].ravel(), z[:, -edge:].ravel()])
        floor = max(self.floor.value() * float(frame.std()), 0.02)
        n = int(self.nlevels.value())
        levels = np.logspace(np.log10(floor), 0, n)
        cmap = pg.colormap.get("viridis")
        colors = cmap.getLookupTable(0.0, 1.0, n)

        # map data-index space to ppm for the isocurves
        f2, f1 = d.f2_ppm, d.f1_ppm
        tr = pg.QtGui.QTransform()
        tr.translate(f2[0], f1[0])
        tr.scale((f2[-1] - f2[0]) / z.shape[1], (f1[-1] - f1[0]) / z.shape[0])
        for lvl, col in zip(levels, colors):
            iso = pg.IsocurveItem(data=z.T, level=lvl,
                                  pen=pg.mkPen(tuple(int(c) for c in col), width=1))
            iso.setTransform(tr)
            self.p_main.addItem(iso)

        # diagonal + projections
        lo = max(min(f2), min(f1)); hi = min(max(f2), max(f1))
        self.p_main.plot([lo, hi], [lo, hi],
                         pen=pg.mkPen("#888", style=Qt.DashLine))
        self.p_top.plot(f2, d.projection("f2"), pen=pg.mkPen("#0e7c86"))
        self.p_left.plot(d.projection("f1"), f1, pen=pg.mkPen("#0e7c86"))
