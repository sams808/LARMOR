"""Herzfeld–Berger sideband analysis: (ζ, η) from the INTENSITIES of a MAS
sideband manifold — a measurement, not a lineshape fit.

Drag the centreband line onto the isotropic peak, choose how many sidebands
each side to integrate; the dialog integrates each window, fits the
normalised intensities to the exact Herzfeld–Berger pattern for the
spectrometer's ν₀ and ν_rot, and reports the tensor in every convention. One
button seeds a csa_mas site so the reading becomes a starting point for a
full lineshape fit. Physics and inversion live in larmor.herzfeld_berger
(Qt-free, tested)."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QDialog, QDoubleSpinBox, QHBoxLayout, QLabel, QPushButton, QSpinBox,
    QTableWidget, QTableWidgetItem, QVBoxLayout,
)

from larmor.desktop import theme


class HerzfeldBergerDialog(QDialog):
    #: (zeta_ppm, eta, delta_iso_ppm, amplitude) — "add as csa_mas line"
    seed_site = Signal(float, float, float, float)

    def __init__(self, parent, ppm, amp, nucleus: str, larmor_MHz: float,
                 spin_rate_Hz: float, spin: float = 0.5):
        super().__init__(parent)
        self.setWindowTitle(f"Herzfeld–Berger sideband analysis — {nucleus or '?'}")
        self.resize(980, 640)
        self.ppm = np.asarray(ppm, float)
        self.amp = np.asarray(amp, float)
        self.larmor_MHz = float(larmor_MHz)
        self.spin_rate_Hz = float(spin_rate_Hz)
        self.step_ppm = self.spin_rate_Hz / self.larmor_MHz
        self.fit = None
        self.measured: dict = {}

        t = theme.active()
        v = QVBoxLayout(self)
        self.plot = pg.PlotWidget(background=t.plot_bg)
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "chemical shift", units="ppm")
        self.plot.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
        self.plot.plot(self.ppm, self.amp, pen=pg.mkPen(t.experiment, width=1.2))
        v.addWidget(self.plot, 1)

        # centreband: start on the tallest peak (the user drags it if a
        # sideband happens to be taller)
        c0 = float(self.ppm[int(np.argmax(self.amp))]) if self.amp.size else 0.0
        self.centre = pg.InfiniteLine(pos=c0, angle=90, movable=True,
                                      pen=pg.mkPen(t.accent, width=2))
        self.centre.label = pg.InfLineLabel(self.centre, "centreband",
                                            position=0.92, color=t.text)
        self.plot.addItem(self.centre)
        self.centre.sigPositionChanged.connect(self._measure)
        self._regions: list[pg.LinearRegionItem] = []
        self._stems = pg.BarGraphItem(x=[], height=[], width=0.3 * self.step_ppm,
                                      brush=pg.mkBrush(t.model), pen=None)
        self._stems.setOpacity(0.55)
        self.plot.addItem(self._stems)

        row = QHBoxLayout()
        row.addWidget(QLabel("sidebands each side"))
        self.n_side = QSpinBox(); self.n_side.setRange(1, 30); self.n_side.setValue(3)
        self.n_side.valueChanged.connect(self._measure)
        row.addWidget(self.n_side)
        row.addWidget(QLabel("integration half-width (ppm)"))
        self.half = QDoubleSpinBox(); self.half.setDecimals(3)
        self.half.setRange(0.001, max(self.step_ppm, 0.002))
        self.half.setValue(0.35 * self.step_ppm)
        self.half.valueChanged.connect(self._measure)
        row.addWidget(self.half)
        row.addWidget(QLabel(f"ν_rot = {self.spin_rate_Hz / 1000:.2f} kHz = "
                             f"{self.step_ppm:.3g} ppm"))
        row.addStretch(1)
        self.btnFit = QPushButton("Fit ζ, η from the intensities")
        self.btnFit.clicked.connect(self._fit)
        row.addWidget(self.btnFit)
        v.addLayout(row)

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["N", "position (ppm)",
                                              "measured", "Herzfeld–Berger"])
        self.table.setMaximumHeight(170)
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table)

        self.res = QLabel("")
        self.res.setStyleSheet(f"color: {t.accent}; font-weight: 600;")
        self.res.setWordWrap(True)
        v.addWidget(self.res)
        hint = ("Integrated intensities of the centreband and ±N sidebands, "
                "normalised over the orders shown, are matched to the exact "
                "sideband pattern of a shielding tensor (Herzfeld & Berger 1980)"
                ". Sign of ζ follows from the +N/−N asymmetry; η is poorly "
                "determined when fewer than two sidebands carry intensity.")
        if spin > 0.5:
            hint += (" This nucleus is quadrupolar: its sidebands also carry "
                     "the satellite/first-order manifold, so read ζ and η as "
                     "an upper bound on the CSA, not a pure measurement.")
        lab = QLabel(hint); lab.setWordWrap(True)
        lab.setStyleSheet(f"color: {t.text_dim};")
        v.addWidget(lab)

        bottom = QHBoxLayout(); bottom.addStretch(1)
        self.btnSeed = QPushButton("Add as CSA line (csa_mas)")
        self.btnSeed.setEnabled(False)
        self.btnSeed.clicked.connect(self._seed)
        btnClose = QPushButton("Close"); btnClose.clicked.connect(self.reject)
        bottom.addWidget(self.btnSeed); bottom.addWidget(btnClose)
        v.addLayout(bottom)
        self._measure()

    # ------------------------------------------------------------------
    def orders(self) -> list[int]:
        n = int(self.n_side.value())
        return list(range(-n, n + 1))

    def _measure(self, *_):
        from larmor.herzfeld_berger import measure_manifold

        c = float(self.centre.value())
        hw = float(self.half.value())
        raw = measure_manifold(self.ppm, self.amp, c, self.larmor_MHz,
                               self.spin_rate_Hz, int(self.n_side.value()),
                               half_width_ppm=hw)
        tot = sum(raw.values())
        self.measured = {n: (v / tot if tot > 0 else 0.0) for n, v in raw.items()}
        t = theme.active()
        for r in self._regions:
            self.plot.removeItem(r)
        self._regions = []
        for n in self.orders():
            pos = c + n * self.step_ppm
            reg = pg.LinearRegionItem(values=(pos - hw, pos + hw), movable=False,
                                      brush=pg.mkBrush(*theme._rgb(t.pivot), 22),
                                      pen=pg.mkPen(t.pivot, width=1, style=Qt.DotLine))
            reg.setZValue(-5)
            self.plot.addItem(reg)
            self._regions.append(reg)
        self.fit = None
        self.btnSeed.setEnabled(False)
        self._stems.setOpts(x=[], height=[])
        self._fill_table()
        self.res.setText("drag the centreband onto the isotropic peak, then Fit")

    def _fill_table(self):
        c = float(self.centre.value())
        orders = self.orders()
        self.table.setRowCount(len(orders))
        for i, n in enumerate(orders):
            vals = [f"{n:+d}", f"{c + n * self.step_ppm:.2f}",
                    f"{self.measured.get(n, 0.0):.4f}",
                    f"{self.fit.fitted[n]:.4f}" if self.fit and n in self.fit.fitted
                    else ""]
            for j, txt in enumerate(vals):
                item = QTableWidgetItem(txt)
                item.setFlags(Qt.ItemIsEnabled)
                self.table.setItem(i, j, item)

    def _fit(self):
        from larmor import convert as C
        from larmor.herzfeld_berger import fit_sidebands

        out = fit_sidebands(self.measured, self.larmor_MHz, self.spin_rate_Hz)
        if not out.ok:
            self.fit = None
            self.res.setText(out.message)
            self.btnSeed.setEnabled(False)
            return
        self.fit = out
        c = float(self.centre.value())
        d11, d22, d33 = C.csa_principal_from_haeberlen(c, out.zeta_ppm, out.eta)
        span, skew = C.csa_span_skew(d11, d22, d33)
        self.res.setText(
            f"ζ = {out.zeta_ppm:.1f} ppm (shielding; δaniso = {-out.zeta_ppm:.1f})"
            f"   η = {out.eta:.2f}   δiso = {c:.2f} ppm   |   "
            f"δ11 / δ22 / δ33 = {d11:.1f} / {d22:.1f} / {d33:.1f} ppm   "
            f"Ω = {span:.1f} ppm   κ = {skew:.2f}   |   "
            f"RMS misfit of intensities {out.rms:.4f}"
            + (f"   — {out.message}" if out.message else ""))
        # fitted intensities as stems, scaled to the measured centreband height
        hw = float(self.half.value())
        m0 = (self.ppm >= c - hw) & (self.ppm <= c + hw)
        peak0 = float(np.max(self.amp[m0])) if m0.any() else float(np.max(self.amp))
        f0 = out.fitted.get(0, 0.0)
        scale = peak0 / f0 if f0 > 0 else 0.0
        xs = [c + n * self.step_ppm for n in out.orders]
        hs = [out.fitted[n] * scale for n in out.orders]
        self._stems.setOpts(x=xs, height=hs, width=0.3 * self.step_ppm)
        self._fill_table()
        self.btnSeed.setEnabled(True)

    def _seed(self):
        if self.fit is None:
            return
        amp = float(np.max(np.abs(self.amp))) if self.amp.size else 1.0
        self.seed_site.emit(self.fit.zeta_ppm, self.fit.eta,
                            float(self.centre.value()), amp)
        self.accept()
