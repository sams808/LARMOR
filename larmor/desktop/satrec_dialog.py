"""Saturation-recovery (T1) analysis dialog: pick an EXPNO, LARMOR processes
every vdlist slice, integrates, fits T1 and shows the build-up curve."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
)


class SatrecDialog(QDialog):
    def __init__(self, parent, expno: str | None):
        super().__init__(parent)
        self.setWindowTitle("Saturation recovery — automatic T1")
        self.resize(820, 560)
        self.expno = expno
        self.result = None

        v = QVBoxLayout(self)
        top = QHBoxLayout()
        self.lbl = QLabel(expno or "no EXPNO selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btnPick = QPushButton("Choose EXPNO…")
        btnPick.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(btnPick)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("EM lb (Hz)"))
        self.lb = QDoubleSpinBox()
        self.lb.setRange(0, 1e5)
        self.lb.setValue(100)
        opts.addWidget(self.lb)
        self.stretched = QCheckBox("stretched exponential (β)")
        opts.addWidget(self.stretched)
        self.magnitude = QCheckBox("magnitude mode")
        self.magnitude.setToolTip("phase-insensitive; use if phasing is unstable")
        opts.addWidget(self.magnitude)
        self.btnRun = QPushButton("Analyze")
        self.btnRun.setDefault(True)
        self.btnRun.clicked.connect(self._run)
        opts.addWidget(self.btnRun)
        opts.addStretch(1)
        v.addLayout(opts)

        self.plot = pg.PlotWidget(background="w")
        self.plot.setLogMode(x=True, y=False)
        self.plot.setLabel("bottom", "recovery delay / s")
        self.plot.setLabel("left", "normalized integral")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        v.addWidget(self.plot, 1)

        bottom = QHBoxLayout()
        self.res_lbl = QLabel("")
        self.res_lbl.setStyleSheet("font-weight: 700; color: #0a5a62; font-size: 14px;")
        bottom.addWidget(self.res_lbl, 1)
        self.btnCsv = QPushButton("Copy CSV")
        self.btnCsv.setEnabled(False)
        self.btnCsv.clicked.connect(self._csv)
        bottom.addWidget(self.btnCsv)
        v.addLayout(bottom)

    def _pick(self):
        path = QFileDialog.getExistingDirectory(
            self, "Bruker EXPNO with ser + vdlist (read-only)")
        if path:
            self.expno = path
            self.lbl.setText(path)

    def _run(self):
        if not self.expno:
            return
        p = Path(self.expno)
        if not (p / "ser").exists() or not (p / "vdlist").exists():
            self.res_lbl.setText("this EXPNO has no ser + vdlist "
                                 "(not a pseudo-2D saturation recovery)")
            return
        self.btnRun.setEnabled(False)
        self.res_lbl.setText("processing slices…")
        QApplication.processEvents()
        try:
            from larmor import satrec

            self.result = satrec.analyze(
                p, lb_hz=self.lb.value(),
                stretched=self.stretched.isChecked(),
                mode="magnitude" if self.magnitude.isChecked() else "phase")
        except Exception as exc:
            self.res_lbl.setText(f"failed: {exc}")
            self.btnRun.setEnabled(True)
            return
        r = self.result
        self.plot.clear()
        pos = r.delays_s > 0
        self.plot.plot(r.delays_s[pos], r.integrals[pos],
                       pen=None, symbol="o", symbolSize=7,
                       symbolBrush=None, symbolPen=pg.mkPen("#0e7c86", width=1.5))
        tt = np.logspace(np.log10(max(r.delays_s[pos].min(), 1e-4)),
                         np.log10(r.delays_s.max() * 1.3), 300)
        self.plot.plot(tt, r.curve(tt), pen=pg.mkPen("#d62728", width=1.6))
        notes = ("   ·   " + " · ".join(r.notes)) if r.notes else ""
        self.res_lbl.setText(r.summary + notes)
        self.btnCsv.setEnabled(True)
        self.btnRun.setEnabled(True)

    def _csv(self):
        if not self.result:
            return
        r = self.result
        lines = ["delay_s,integral_norm"]
        lines += [f"{d},{i}" for d, i in zip(r.delays_s, r.integrals)]
        lines.append(f"# {r.summary}  window {r.window_ppm[0]:.1f}..{r.window_ppm[1]:.1f} ppm")
        QApplication.clipboard().setText("\n".join(lines))
        self.res_lbl.setText(r.summary + "   ·   copied as CSV")
