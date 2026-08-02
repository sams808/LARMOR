"""Monte-Carlo error estimation (dmfit 'Errors ▸ Monte Carlo').

Adds synthetic noise (at the residual level) to the best-fit model, re-fits N
times, and reports each parameter as mean ± σ with a percentage and a
distribution histogram — the same statistics as dmfit / pydmfit
(errorsMonteCarlo.py). Complements the covariance stderr and the χ² profile.
The physics is in larmor.autofit.monte_carlo_errors; this is only the UI.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QDialogButtonBox, QHBoxLayout, QLabel,
    QProgressBar, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
    QVBoxLayout,
)

from larmor.desktop import theme


class MonteCarloDialog(QDialog):
    def __init__(self, parent, recipe: dict, ppm, amp, window):
        super().__init__(parent)
        self.setWindowTitle("Monte-Carlo errors — synthetic-noise refits")
        self.resize(760, 620)
        self.recipe, self.ppm, self.amp, self.window = recipe, ppm, amp, window
        self._result = None
        self._stop = False
        v = QVBoxLayout(self)

        intro = QLabel(
            "Adds Gaussian noise at the residual level to the best-fit model and "
            "re-fits N times; the spread of each parameter is its error. Captures "
            "correlations and non-linearity that the covariance matrix misses "
            "(dmfit ▸ Errors ▸ Monte Carlo).")
        intro.setWordWrap(True)
        intro.setStyleSheet(f"color:{theme.active().text_dim};")
        v.addWidget(intro)

        row = QHBoxLayout()
        row.addWidget(QLabel("trials"))
        self.n = QSpinBox(); self.n.setRange(10, 5000); self.n.setValue(200)
        self.n.setToolTip("more trials = smoother error estimate, slower")
        row.addWidget(self.n)
        row.addWidget(QLabel("seed"))
        self.seed = QSpinBox(); self.seed.setRange(0, 10_000); self.seed.setValue(0)
        row.addWidget(self.seed)
        self.btnRun = QPushButton("Run"); self.btnRun.setDefault(True)
        self.btnRun.clicked.connect(self._run)
        row.addWidget(self.btnRun)
        self.btnStop = QPushButton("Stop"); self.btnStop.setEnabled(False)
        self.btnStop.clicked.connect(lambda: setattr(self, "_stop", True))
        row.addWidget(self.btnStop)
        row.addStretch(1)
        v.addLayout(row)

        self.prog = QProgressBar(); self.prog.setValue(0)
        v.addWidget(self.prog)

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["parameter", "best", "mean", "± σ (MC)", "σ %"])
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table, 1)

        hrow = QHBoxLayout()
        hrow.addWidget(QLabel("histogram"))
        self.pick = QComboBox()
        self.pick.currentIndexChanged.connect(self._draw_hist)
        hrow.addWidget(self.pick, 1)
        v.addLayout(hrow)
        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.setLabel("bottom", "parameter value")
        self.plot.setLabel("left", "count")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        self.plot.setMinimumHeight(160)
        v.addWidget(self.plot)

        self.status = QLabel("")
        self.status.setStyleSheet(f"font-weight:600; color:{theme.active().accent};")
        v.addWidget(self.status)

        bb = QDialogButtonBox(QDialogButtonBox.Close | QDialogButtonBox.Help)
        self.btnApply = bb.addButton("Use as fit errors",
                                     QDialogButtonBox.ApplyRole)
        self.btnApply.setEnabled(False)
        self.btnApply.setToolTip("write the MC σ into each parameter's error "
                                 "in the fit table")
        self.btnApply.clicked.connect(self._apply_errors)
        self.btnCopy = bb.addButton("Copy report", QDialogButtonBox.ActionRole)
        self.btnCopy.setEnabled(False)
        self.btnCopy.clicked.connect(self._copy)
        bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.helpRequested.connect(self._help)
        v.addWidget(bb)

    # ------------------------------------------------------------------
    def _run(self):
        from larmor import autofit
        from larmor.recipe import Recipe

        self._stop = False
        self.btnRun.setEnabled(False); self.btnStop.setEnabled(True)
        self.status.setText("running…")
        n = self.n.value()
        self.prog.setRange(0, n)

        def prog(k, ntot):
            self.prog.setValue(k)
            if k % 5 == 0 or k == ntot:
                QApplication.processEvents()

        try:
            self._result = autofit.monte_carlo_errors(
                Recipe.from_dict(self.recipe), self.ppm, self.amp,
                window_ppm=self.window, n_trials=n, seed=self.seed.value(),
                progress=prog, should_stop=lambda: self._stop)
        except Exception as exc:
            self.status.setText(f"failed: {exc}")
            self.btnRun.setEnabled(True); self.btnStop.setEnabled(False)
            return

        self.btnRun.setEnabled(True); self.btnStop.setEnabled(False)
        self._fill_table()
        self.btnApply.setEnabled(True); self.btnCopy.setEnabled(True)
        self.status.setText(self._result.summary
                            + ("  ·  stopped early" if self._stop else ""))

    def _fill_table(self):
        r = self._result
        self.table.setRowCount(len(r.params))
        self.pick.blockSignals(True); self.pick.clear()
        for i, p in enumerate(r.params):
            pc = f"{p.pct:.2f}" if np.isfinite(p.pct) else "—"
            for c, val in enumerate((p.label, f"{p.best:.5g}", f"{p.mean:.5g}",
                                     f"{p.std:.4g}", pc)):
                it = QTableWidgetItem(val)
                if c:
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.table.setItem(i, c, it)
            self.pick.addItem(p.label, i)
        self.pick.blockSignals(False)
        self.table.resizeColumnsToContents()
        self._draw_hist()

    def _draw_hist(self):
        self.plot.clear()
        if self._result is None or self.pick.currentIndex() < 0:
            return
        p = self._result.params[self.pick.currentData()]
        vals = p.values
        if vals.size < 2 or p.std <= 0:
            return
        counts, edges = np.histogram(vals, bins=min(20, max(6, vals.size // 8)))
        r, g, b = theme._rgb(theme.active().measure)
        self.plot.plot(edges, counts, stepMode="center",
                       fillLevel=0, brush=(r, g, b, 90),
                       pen=pg.mkPen(theme.active().measure, width=1))
        # Gaussian(mean, σ) overlay, scaled to the histogram
        xs = np.linspace(edges[0], edges[-1], 200)
        g = np.exp(-0.5 * ((xs - p.mean) / p.std) ** 2)
        g = g / g.max() * counts.max() if counts.max() else g
        self.plot.plot(xs, g, pen=pg.mkPen(theme.active().model, width=1.6))
        self.plot.addItem(pg.InfiniteLine(
            pos=p.best, angle=90,
            pen=pg.mkPen(theme.active().text_dim, style=Qt.DashLine)))

    def _apply_errors(self):
        if self._result is None:
            return
        for mp in self._result.params:
            try:
                self.recipe["sites"][mp.site]["params"][mp.param]["stderr"] = \
                    float(mp.std)
            except (KeyError, IndexError, ValueError):
                pass
        main = self.parent()
        if main is not None and hasattr(main, "lines_table"):
            main.lines_table.rebuild(self.recipe, getattr(main, "hidden", set()))
            main.statusBar().showMessage(
                "Monte-Carlo σ written to the fit table")
        self.status.setText(self._result.summary + "  ·  errors applied")

    def _copy(self):
        if self._result is not None:
            QApplication.clipboard().setText(self._result.report())
            self.status.setText("report copied to clipboard")

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "spectra-1d", "1D spectra — processing & fitting")
