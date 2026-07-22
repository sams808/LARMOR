"""Variable-temperature dialog: enter (T, rate) points and fit Arrhenius or VFT."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout,
)

from larmor import vt


class VtDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("Variable temperature — Arrhenius / VFT")
        self.resize(760, 560)
        v = QVBoxLayout(self)

        bar = QHBoxLayout()
        bar.addWidget(QLabel("model"))
        self.mode = QComboBox(); self.mode.addItems(["Arrhenius", "VFT"])
        bar.addWidget(self.mode)
        b1 = QPushButton("+ row"); b1.clicked.connect(lambda: self._add_row())
        b2 = QPushButton("Paste (T rate per line)"); b2.clicked.connect(self._paste)
        self.btnFit = QPushButton("Fit"); self.btnFit.clicked.connect(self._fit)
        bar.addWidget(b1); bar.addWidget(b2); bar.addWidget(self.btnFit)
        bar.addStretch(1)
        v.addLayout(bar)

        body = QHBoxLayout()
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["T (K)", "rate"])
        self.table.setMaximumWidth(220)
        body.addWidget(self.table)
        self.plot = pg.PlotWidget(background="#fcfdfc")
        self.plot.setLabel("bottom", "1000/T", units="1/K")
        self.plot.setLabel("left", "ln(rate)")
        body.addWidget(self.plot, 1)
        v.addLayout(body)

        self.res = QLabel(""); self.res.setStyleSheet(
            "font-weight: 700; color: #0a5a62;")
        v.addWidget(self.res)
        for _ in range(3):
            self._add_row()

    def _add_row(self, T="", rate=""):
        r = self.table.rowCount(); self.table.insertRow(r)
        self.table.setItem(r, 0, QTableWidgetItem(str(T)))
        self.table.setItem(r, 1, QTableWidgetItem(str(rate)))

    def _paste(self):
        txt = QApplication.clipboard().text()
        self.table.setRowCount(0)
        for line in txt.splitlines():
            parts = line.replace(",", " ").split()
            if len(parts) >= 2:
                self._add_row(parts[0], parts[1])

    def _data(self):
        T, k = [], []
        for r in range(self.table.rowCount()):
            try:
                a = float(self.table.item(r, 0).text())
                b = float(self.table.item(r, 1).text())
            except (ValueError, AttributeError):
                continue
            if b > 0:
                T.append(a); k.append(b)
        return np.array(T), np.array(k)

    def _fit(self):
        T, k = self._data()
        if T.size < 3:
            self.res.setText("need at least 3 (T, rate) points with rate > 0")
            return
        try:
            if self.mode.currentText() == "Arrhenius":
                fit = vt.fit_arrhenius(T, k)
                self.res.setText(f"Ea = {fit['Ea_kJmol']:.2f} kJ/mol   ·   "
                                 f"A = {fit['A']:.3g}   ·   rmsd(ln k) "
                                 f"{fit['rmsd_lnk']:.3g}")
            else:
                fit = vt.fit_vft(T, k)
                self.res.setText(f"A = {fit['A']:.3g}   ·   B = {fit['B_K']:.1f} K"
                                 f"   ·   T0 = {fit['T0_K']:.1f} K")
        except Exception as exc:
            self.res.setText(f"fit failed: {exc}")
            return
        self.plot.clear()
        self.plot.plot(1000.0 / T, np.log(k), pen=None, symbol="o",
                       symbolBrush="#0e7c86", symbolSize=7)
        tt = np.linspace(T.min(), T.max(), 200)
        self.plot.plot(1000.0 / tt, np.log(fit["curve"](tt)),
                       pen=pg.mkPen("#d62728", width=1.6))
