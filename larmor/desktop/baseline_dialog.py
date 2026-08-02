"""Interactive iterative baseline corrector (Yon et al. 2020).

A LARMOR re-implementation of the workflow of the *MY Baseline Corrector* GUI:
tune the parameters and see the estimated baseline and the corrected spectrum
update live, then Apply. The physics is in larmor.baseline; this is only the UI.

Reference (shown in the dialog and required when you publish results):
    M. Yon, F. Fayon, D. Massiot, V. Sarou-Kanian, "Iterative baseline
    correction algorithm for dead time truncated one-dimensional solid-state MAS
    NMR spectra", Solid State Nucl. Magn. Reson. 110, 101699 (2020),
    doi:10.1016/j.ssnmr.2020.101699 ; github.com/maximeYon/Baseline_Corrector
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
    QHBoxLayout, QLabel, QSpinBox, QVBoxLayout,
)

from larmor.baseline import iterative_baseline
from larmor.desktop import theme


class BaselineDialog(QDialog):
    """Preview + apply the iterative dead-time baseline correction."""

    def __init__(self, parent, ppm: np.ndarray, amp: np.ndarray):
        super().__init__(parent)
        self.setWindowTitle("Iterative baseline — Yon et al. 2020")
        self.resize(820, 640)
        self.ppm = np.asarray(ppm, float)
        self.amp = np.asarray(amp, float)

        v = QVBoxLayout(self)
        intro = QLabel(
            "For a <b>rolling baseline from receiver dead time</b> in "
            "pulse-acquire MAS, where a polynomial fails. Baseline points are "
            "picked automatically by a histogram filter (noise band, excluding "
            "peaks and negative spikes); a smoothing spline is fit and subtracted "
            "iteratively. Optionally restrict the baseline to broad (dead-time) "
            "components in the time domain.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        # ---- plots: spectrum + baseline (top), corrected (bottom) ----
        t = theme.active()
        glw = pg.GraphicsLayoutWidget()
        glw.setBackground(t.plot_bg)
        self.p_top = glw.addPlot(row=0, col=0)
        self.p_top.showGrid(x=True, y=True, alpha=0.12)
        self.p_top.setLabel("left", "intensity")
        self.p_top.addLegend(offset=(-10, 10))
        self.p_bot = glw.addPlot(row=1, col=0)
        self.p_bot.showGrid(x=True, y=True, alpha=0.12)
        self.p_bot.setLabel("left", "corrected")
        self.p_bot.setLabel("bottom", "shift", units="ppm")
        self.p_bot.setXLink(self.p_top)
        glw.ci.layout.setRowStretchFactor(0, 3)
        glw.ci.layout.setRowStretchFactor(1, 2)
        v.addWidget(glw, 1)

        self.c_raw = self.p_top.plot(pen=pg.mkPen(t.experiment, width=1),
                                     name="spectrum")
        self.c_base = self.p_top.plot(pen=pg.mkPen(t.baseline, width=1.6),
                                      name="baseline")
        self.c_corr = self.p_bot.plot(pen=pg.mkPen(t.measure, width=1),
                                      name="corrected")
        for p in (self.p_top, self.p_bot):
            p.setXRange(self.ppm.max(), self.ppm.min())   # NMR: decreasing ppm

        # ---- controls ----
        form = QFormLayout()
        self.sp_dead = QSpinBox()
        self.sp_dead.setRange(0, 8192)
        self.sp_dead.setSingleStep(2)
        self.sp_dead.setValue(0)
        self.sp_dead.setToolTip(
            "Dead-time restriction: time-domain points to keep (≈ 2·DE/DW).\n"
            "0 = plain iterative histogram baseline (no dead-time restriction).")
        form.addRow("Dead-time points (0 = off):", self.sp_dead)

        self.sp_smooth = QDoubleSpinBox()
        self.sp_smooth.setDecimals(3)
        self.sp_smooth.setRange(0.001, 1000.0)
        self.sp_smooth.setValue(1.0)
        self.sp_smooth.setToolTip("Smoothing-spline stiffness (larger = smoother "
                                  "baseline). 1.0 suits a slow rolling baseline.")
        form.addRow("Smoothness:", self.sp_smooth)

        self.sp_thr = QDoubleSpinBox()
        self.sp_thr.setDecimals(2)
        self.sp_thr.setRange(0.2, 3.0)
        self.sp_thr.setSingleStep(0.1)
        self.sp_thr.setValue(1.0)
        self.sp_thr.setToolTip("Histogram threshold width (dmfit AdvanceFilter/100). "
                               "Smaller keeps fewer points as baseline.")
        form.addRow("Threshold factor:", self.sp_thr)
        v.addLayout(form)

        self.status = QLabel("")
        self.status.setStyleSheet("color:#555")
        v.addWidget(self.status)

        cite = QLabel(
            "Method: M. Yon, F. Fayon, D. Massiot, V. Sarou-Kanian, "
            "<i>Solid State Nucl. Magn. Reson.</i> <b>110</b>, 101699 (2020), "
            "doi:10.1016/j.ssnmr.2020.101699 — cite this if you publish results.")
        cite.setWordWrap(True)
        cite.setStyleSheet("color:#666; font-size:11px")
        v.addWidget(cite)

        bb = QDialogButtonBox(QDialogButtonBox.Apply | QDialogButtonBox.Cancel
                              | QDialogButtonBox.Help)
        bb.button(QDialogButtonBox.Apply).clicked.connect(self.accept)
        bb.rejected.connect(self.reject)
        bb.helpRequested.connect(lambda: self._help())
        row = QHBoxLayout()
        row.addWidget(bb)
        v.addLayout(row)

        # live preview (debounced so dragging a spinbox stays smooth)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(120)
        self._timer.timeout.connect(self._preview)
        for w in (self.sp_dead, self.sp_smooth, self.sp_thr):
            w.valueChanged.connect(lambda *_: self._timer.start())
        self._preview()

    # ------------------------------------------------------------------
    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "processing-reference", "Processing reference")

    def params(self) -> dict:
        return {"dead_time_pts": int(self.sp_dead.value()),
                "smoothness": float(self.sp_smooth.value()),
                "threshold_factor": float(self.sp_thr.value())}

    def _preview(self):
        try:
            r = iterative_baseline(self.amp, **self.params())
        except Exception as exc:  # noqa: BLE001 - keep the dialog alive
            self.status.setText(f"preview failed: {exc}")
            return
        self.c_raw.setData(self.ppm, self.amp)
        self.c_base.setData(self.ppm, r.baseline)
        self.c_corr.setData(self.ppm, r.corrected)
        self.status.setText(
            f"{r.n_iter} iterations, "
            f"{'converged' if r.converged else 'stopped at max_iter'} · "
            f"baseline peak = {np.abs(r.baseline).max():.3g} "
            f"({100 * np.abs(r.baseline).max() / (np.abs(self.amp).max() or 1):.1f}% "
            "of the spectrum)")
