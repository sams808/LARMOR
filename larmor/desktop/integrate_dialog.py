"""Integration / measurement dialog: drag regions over a spectrum and read the
integral, percent, centre of mass and FWHM of each — dmfit/TopSpin integration,
ssNake FWHM / Centre of Mass / Integrals."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QDialog, QHBoxLayout, QLabel, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)

from larmor import measure as M

REGION_COLORS = ["#0e7c86", "#d62728", "#2ca02c", "#9467bd", "#e377c2",
                 "#8c564b", "#17becf", "#bcbd22"]


class IntegralsDialog(QDialog):
    def __init__(self, parent, ppm, amp):
        super().__init__(parent)
        self.setWindowTitle("Integrals & measurements")
        self.resize(900, 620)
        self.ppm = np.asarray(ppm, float)
        self.amp = np.asarray(amp, float)
        self.regions: list[pg.LinearRegionItem] = []

        v = QVBoxLayout(self)
        bar = QHBoxLayout()
        self.btnAdd = QPushButton("+ Add region")
        self.btnAdd.clicked.connect(self._add_region)
        self.btnClear = QPushButton("Clear")
        self.btnClear.clicked.connect(self._clear)
        self.btnCsv = QPushButton("Copy CSV")
        self.btnCsv.clicked.connect(self._csv)
        bar.addWidget(self.btnAdd); bar.addWidget(self.btnClear)
        bar.addWidget(self.btnCsv); bar.addStretch(1)
        self.hint = QLabel("drag the shaded edges; percents are over all regions")
        self.hint.setStyleSheet("color: #5a6871;")
        bar.addWidget(self.hint)
        v.addLayout(bar)

        self.plot = pg.PlotWidget(background="#fcfdfc")
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "shift", units="ppm")
        self.plot.plot(self.ppm, self.amp, pen=pg.mkPen("#1a2831", width=1.2))
        v.addWidget(self.plot, 1)

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["region (ppm)", "integral", "%", "centre (ppm)", "FWHM (ppm)"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setMaximumHeight(180)
        v.addWidget(self.table)

        if self.ppm.size:
            self._add_region()

    # ------------------------------------------------------------------
    def _add_region(self):
        if not self.ppm.size:
            return
        i = len(self.regions)
        col = REGION_COLORS[i % len(REGION_COLORS)]
        (x0, x1), _ = self.plot.getPlotItem().getViewBox().viewRange()
        lo, hi = sorted((x0, x1))
        span = hi - lo
        c = pg.mkColor(col); c.setAlpha(40)
        r = pg.LinearRegionItem(values=(lo + 0.4 * span, lo + 0.6 * span),
                                brush=pg.mkBrush(c), pen=pg.mkPen(col, width=1.4))
        r.sigRegionChanged.connect(self._recompute)
        self.plot.addItem(r)
        self.regions.append(r)
        self._recompute()

    def _clear(self):
        for r in self.regions:
            self.plot.removeItem(r)
        self.regions.clear()
        self._recompute()

    def _ranges(self):
        return [tuple(r.getRegion()) for r in self.regions]

    def _rows(self):
        return M.integrate_regions(self.ppm, self.amp, self._ranges())

    def _recompute(self, *_):
        rows = self._rows()
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            hi, lo = row["range"]
            vals = [f"{hi:.2f} … {lo:.2f}", f"{row['integral']:.4g}",
                    f"{row['percent']:.1f}", f"{row['centre']:.2f}",
                    f"{row['fwhm']:.2f}"]
            for j, val in enumerate(vals):
                self.table.setItem(i, j, QTableWidgetItem(val))

    def _csv(self):
        lines = ["hi_ppm,lo_ppm,integral,percent,centre_ppm,fwhm_ppm"]
        for row in self._rows():
            hi, lo = row["range"]
            lines.append(f"{hi:.3f},{lo:.3f},{row['integral']:.6g},"
                         f"{row['percent']:.3f},{row['centre']:.3f},{row['fwhm']:.3f}")
        QApplication.clipboard().setText("\n".join(lines))
        self.hint.setText("copied CSV to clipboard")
