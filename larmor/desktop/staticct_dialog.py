"""Read (C_Q, η, δiso) off the static CT pattern on screen — drag three
markers, no fit.

Two solid lines go on the pattern's two divergent horns; the dashed one goes
on the pattern's informative outer limit (the step-like edge on the far side
— the other limit usually coincides with a horn and is refused with a hint).
The reading updates live; one button seeds a quad_ct site from it.
Physics and inversion live in larmor.staticct (Qt-free, tested)."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QDialog, QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
)

from larmor.desktop import theme


class StaticCtDialog(QDialog):
    #: (Cq_MHz, eta, delta_iso_ppm) — "add as quad_ct line"
    seed_site = Signal(float, float, float)

    def __init__(self, parent, ppm, amp, nucleus: str, spin: float,
                 larmor_MHz: float):
        super().__init__(parent)
        self.setWindowTitle(f"Read static pattern — {nucleus or '?'}")
        self.resize(920, 560)
        self.ppm = np.asarray(ppm, float)
        self.amp = np.asarray(amp, float)
        self.spin = float(spin)
        self.larmor_MHz = float(larmor_MHz)
        self.reading = None

        v = QVBoxLayout(self)
        t = theme.active()
        self.plot = pg.PlotWidget(background=t.plot_bg)
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "chemical shift", units="ppm")
        self.plot.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
        self.plot.plot(self.ppm, self.amp,
                       pen=pg.mkPen(t.experiment, width=1.2))
        v.addWidget(self.plot, 1)

        lo, hi = float(self.ppm.min()), float(self.ppm.max())
        span = hi - lo
        mk = dict(angle=90, movable=True)
        self.horn_a = pg.InfiniteLine(pos=lo + 0.35 * span,
                                      pen=pg.mkPen(t.accent, width=2), **mk)
        self.horn_b = pg.InfiniteLine(pos=lo + 0.60 * span,
                                      pen=pg.mkPen(t.accent, width=2), **mk)
        self.edge = pg.InfiniteLine(pos=lo + 0.90 * span,
                                    pen=pg.mkPen(t.model, width=2,
                                                 style=Qt.DashLine), **mk)
        for line, label in ((self.horn_a, "horn"), (self.horn_b, "horn"),
                            (self.edge, "edge")):
            line.label = pg.InfLineLabel(line, label, position=0.92,
                                         color=t.text)
            self.plot.addItem(line)
            line.sigPositionChanged.connect(self._recompute)

        self.res = QLabel("")
        self.res.setStyleSheet(f"color: {t.accent}; font-weight: 600;")
        self.res.setWordWrap(True)
        v.addWidget(self.res)
        hint = QLabel(
            "Solid lines on the two divergent horns; the dashed line on the "
            "pattern's step-like outer limit (the side away from the taller "
            "horn — the other limit sits ON a horn and carries no extra "
            "information). δ_CG plus this reading is what a single-field "
            "static pattern can honestly give.")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color: {t.text_dim};")
        v.addWidget(hint)

        bottom = QHBoxLayout()
        bottom.addStretch(1)
        self.btnSeed = QPushButton("Add as quad_ct line")
        self.btnSeed.setEnabled(False)
        self.btnSeed.clicked.connect(self._seed)
        btnClose = QPushButton("Close")
        btnClose.clicked.connect(self.reject)
        bottom.addWidget(self.btnSeed)
        bottom.addWidget(btnClose)
        v.addLayout(bottom)
        self._recompute()

    def _recompute(self, *_):
        from larmor.staticct import read_cq_eta

        out = read_cq_eta(float(self.horn_a.value()),
                          float(self.horn_b.value()),
                          float(self.edge.value()),
                          self.spin, self.larmor_MHz)
        self.reading = out if out.ok else None
        self.btnSeed.setEnabled(out.ok)
        if out.ok:
            self.res.setText(
                f"C_Q = {out.Cq_MHz:.2f} MHz   η = {out.eta:.2f}   "
                f"P_Q = {out.Pq_MHz:.2f} MHz   δiso = "
                f"{out.delta_iso_ppm:.1f} ppm   "
                f"(feature mismatch {out.resid_ppm:.1f} ppm)")
        else:
            self.res.setText(out.message or "place the three markers")

    def _seed(self):
        if self.reading is not None:
            r = self.reading
            self.seed_site.emit(r.Cq_MHz, r.eta, r.delta_iso_ppm)
            self.accept()
