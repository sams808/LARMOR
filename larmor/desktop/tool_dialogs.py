"""Tool dialogs: REDOR and Errors Analysis (the DFT .magres import lives in
magres_dialog.py)."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from larmor.desktop import theme
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QFileDialog, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QVBoxLayout,
)


class RedorDialog(QDialog):
    def __init__(self, parent, expno: str | None):
        super().__init__(parent)
        self.setWindowTitle("REDOR — dipolar coupling & distance")
        self.resize(760, 520)
        self.expno = expno
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.lbl = QLabel(expno or "no EXPNO selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Choose EXPNO…")
        btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(btn)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("observed"))
        self.iso1 = QLineEdit("13C"); self.iso1.setFixedWidth(60)
        opts.addWidget(self.iso1)
        opts.addWidget(QLabel("dephased by"))
        self.iso2 = QLineEdit("15N"); self.iso2.setFixedWidth(60)
        opts.addWidget(self.iso2)
        opts.addWidget(QLabel("regime"))
        self.regime = QComboBox()
        self.regime.addItems(["auto", "short", "pair"])
        opts.addWidget(self.regime)
        self.btnRun = QPushButton("Analyze")
        self.btnRun.setDefault(True)
        self.btnRun.clicked.connect(self._run)
        opts.addWidget(self.btnRun)
        opts.addStretch(1)
        v.addLayout(opts)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.setLabel("bottom", "recoupling time / s")
        self.plot.setLabel("left", "ΔS/S₀")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        v.addWidget(self.plot, 1)

        self.res = QLabel("")
        self.res.setStyleSheet(f"font-weight: 700; color: {theme.active().accent}; font-size: 14px;")
        self.res.setWordWrap(True)
        v.addWidget(self.res)

    def _pick(self):
        p = QFileDialog.getExistingDirectory(self, "EXPNO with redor.txt")
        if p:
            self.expno = p
            self.lbl.setText(p)

    def _run(self):
        if not self.expno:
            return
        from larmor import redor

        pair = (self.iso1.text().strip(), self.iso2.text().strip())
        try:
            res = redor.analyze_expno(self.expno, pair=pair,
                                      regime=self.regime.currentText())
        except Exception as exc:
            self.res.setText(f"failed: {exc}")
            return
        self.plot.clear()
        self.plot.plot(res.ntr_s, res.ds_s0, pen=None, symbol="o", symbolSize=7,
                       symbolBrush=None, symbolPen=pg.mkPen("#0e7c86", width=1.5))
        tt = np.linspace(res.ntr_s.min(), res.ntr_s.max(), 200)
        self.plot.plot(tt, res.curve(tt), pen=pg.mkPen("#d62728", width=1.6))
        self.res.setText(res.summary + "   ·   " + " · ".join(res.notes))


class ErrorsDialog(QDialog):
    def __init__(self, parent, recipe: dict, ppm, amp, window):
        super().__init__(parent)
        self.setWindowTitle("Errors Analysis — χ² profile")
        self.resize(720, 520)
        self.recipe, self.ppm, self.amp, self.window = recipe, ppm, amp, window
        v = QVBoxLayout(self)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("site"))
        self.site = QComboBox()
        for i, s in enumerate(recipe["sites"]):
            self.site.addItem(f"s{i} — {s.get('label') or s['model']}", i)
        opts.addWidget(self.site)
        opts.addWidget(QLabel("parameter"))
        self.param = QComboBox()
        opts.addWidget(self.param)
        self.btnRun = QPushButton("Scan")
        self.btnRun.setDefault(True)
        self.btnRun.clicked.connect(self._run)
        opts.addWidget(self.btnRun)
        opts.addStretch(1)
        v.addLayout(opts)
        self.site.currentIndexChanged.connect(self._fill_params)
        self._fill_params()

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.setLabel("bottom", "parameter value")
        self.plot.setLabel("left", "χ²")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        v.addWidget(self.plot, 1)
        self.res = QLabel("")
        self.res.setStyleSheet(f"font-weight: 700; color: {theme.active().accent};")
        v.addWidget(self.res)

    def _fill_params(self):
        i = self.site.currentData()
        self.param.clear()
        varying = [n for n, p in self.recipe["sites"][i]["params"].items()
                   if p.get("vary", True) and not p.get("expr")]
        self.param.addItems(varying)

    def _run(self):
        from larmor import autofit
        from larmor.recipe import Recipe

        i = self.site.currentData()
        pname = self.param.currentText()
        if not pname:
            return
        self.res.setText("scanning…")
        QApplication.processEvents()
        try:
            prof = autofit.error_profile(
                Recipe.from_dict(self.recipe), self.ppm, self.amp,
                site=i, param=pname, window_ppm=self.window, parallel=True)
        except Exception as exc:
            self.res.setText(f"failed: {exc}")
            return
        self.plot.clear()
        self.plot.plot(prof.values, prof.chi2, pen=pg.mkPen("#0e7c86", width=1.5),
                       symbol="o", symbolSize=5)
        # 1sigma and 2sigma levels, in units of the residual variance
        # chi2_min / dof (see autofit.error_profile)
        for lvl, col in ((prof.level68, "#d62728"),
                         (prof.level95, "#c88a1e")):
            line = pg.InfiniteLine(pos=lvl, angle=0,
                                   pen=pg.mkPen(col, style=Qt.DashLine))
            self.plot.addItem(line)
        self.res.setText(prof.summary + ("   ·   " + " · ".join(prof.notes)
                                         if prof.notes else ""))
