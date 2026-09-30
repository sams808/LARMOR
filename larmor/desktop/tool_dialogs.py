"""Tool dialogs: REDOR and Errors Analysis (the DFT .magres import lives in
magres_dialog.py)."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from larmor.desktop import theme
from larmor.desktop.axes import plain_units
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QHBoxLayout,
    QLabel, QLineEdit, QProgressBar, QPushButton, QSpinBox, QVBoxLayout,
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
        plain_units(self.plot)
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
    """Errors Analysis (dmfit): the χ² profile of one parameter, every other
    free parameter refitted at each scan point.

    Runs synchronously (the scan is a plain function), but never freezes:
    the profile's ``heartbeat`` pumps the event loop while the process pool
    starts and between points, so the progress bar moves and Stop works.
    Before 0.15.1 the pool's ~15 s start-up showed as a dead "scanning…"
    label in a modal window, which on a single-line fit was the whole run
    (2 s of actual work) -- the tool looked broken.
    """

    def __init__(self, parent, recipe: dict, ppm, amp, window):
        super().__init__(parent)
        self.setWindowTitle("Errors Analysis — χ² profile")
        self.resize(760, 560)
        self.recipe, self.ppm, self.amp, self.window = recipe, ppm, amp, window
        self._stop = False
        self.prof = None
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
        opts.addWidget(QLabel("span"))
        self.span = QDoubleSpinBox()
        self.span.setRange(0.2, 30.0); self.span.setDecimals(1)
        self.span.setSingleStep(0.5); self.span.setValue(3.0)
        self.span.setToolTip("how far to scan each way, in multiples of the "
                             "parameter's error bar (25 % of its value when the "
                             "fit gave none); the 1σ crossings must fall inside "
                             "the scan -- widen when 'not bracketed', narrow when "
                             "'narrower than the scan step'")
        opts.addWidget(self.span)
        opts.addWidget(QLabel("points"))
        self.points = QSpinBox()
        self.points.setRange(5, 61); self.points.setValue(15)
        self.points.setToolTip("scan points (each is a full refit of the other "
                               "parameters)")
        opts.addWidget(self.points)
        self.btnRun = QPushButton("Scan")
        self.btnRun.setDefault(True)
        self.btnRun.clicked.connect(self._run)
        opts.addWidget(self.btnRun)
        self.btnStop = QPushButton("Stop")
        self.btnStop.setEnabled(False)
        self.btnStop.setToolTip("keep the points scanned so far (at least three "
                                "are needed for an interval)")
        self.btnStop.clicked.connect(lambda: setattr(self, "_stop", True))
        opts.addWidget(self.btnStop)
        opts.addStretch(1)
        v.addLayout(opts)
        self.site.currentIndexChanged.connect(self._fill_params)
        self._fill_params()

        self.prog = QProgressBar()
        self.prog.setValue(0)
        self.prog.setFormat("%v / %m points")
        v.addWidget(self.prog)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        plain_units(self.plot)
        self.plot.setLabel("bottom", "parameter value")
        self.plot.setLabel("left", "χ²")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        v.addWidget(self.plot, 1)
        self.res = QLabel("")
        self.res.setWordWrap(True)
        self.res.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.res.setStyleSheet(f"font-weight: 700; color: {theme.active().accent};")
        v.addWidget(self.res)
        self.note = QLabel("")
        self.note.setWordWrap(True)
        self.note.setStyleSheet(f"color: {theme.active().text_dim};")
        v.addWidget(self.note)

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
        n = int(self.points.value())
        self._stop = False
        self.prof = None
        self.btnRun.setEnabled(False)
        self.btnStop.setEnabled(True)
        self.prog.setRange(0, n)
        self.prog.setValue(0)
        self.res.setText("scanning — each point is a refit of the other parameters…")
        self.note.setText("")
        QApplication.processEvents()

        def prog(k, ntot, _value):
            self.prog.setValue(k)
            QApplication.processEvents()

        try:
            prof = autofit.error_profile(
                Recipe.from_dict(self.recipe), self.ppm, self.amp,
                site=i, param=pname, window_ppm=self.window,
                n_points=n, span=float(self.span.value()),
                progress=prog, should_stop=lambda: self._stop,
                heartbeat=QApplication.processEvents, parallel="auto")
        except Exception as exc:
            self.res.setText(("stopped — " if self._stop else "failed: ") + str(exc))
            self.btnRun.setEnabled(True)
            self.btnStop.setEnabled(False)
            return
        self.btnRun.setEnabled(True)
        self.btnStop.setEnabled(False)
        self.prof = prof
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
        # where the fit on screen sits: a vertical marker at its value
        try:
            fitted = float(self.recipe["sites"][i]["params"][pname]["value"])
            self.plot.addItem(pg.InfiniteLine(
                pos=fitted, angle=90,
                pen=pg.mkPen(theme.active().text_dim, style=Qt.DashLine)))
        except (KeyError, TypeError, ValueError):
            pass
        self.res.setText(prof.summary + ("   ·   stopped early" if self._stop else ""))
        notes = list(prof.notes)
        if notes:
            self.note.setText(" · ".join(notes))
            self.note.setStyleSheet(
                f"color: {'#A8570F' if not prof.fit_at_minimum else theme.active().text_dim};")
