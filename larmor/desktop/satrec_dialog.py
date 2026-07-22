"""Guided relaxation (T1/T2) analysis — the TopSpin guided-T1 workflow.

Steps, top to bottom:
  1. open a pseudo-2D (ser) — the arrayed relaxation experiment;
  2. LARMOR processes every slice and shows the LAST slice (fully relaxed,
     most signal) as the reference spectrum;
  3. the user drags one or more integration ZONES over the peaks of interest;
  4. each zone becomes a build-up curve — the points are plotted (delay vs
     integral), one series per zone;
  5. the user clicks any aberrant point to exclude it;
  6. Fit gives T1 (or T2) per zone on the surviving points.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QLabel, QPushButton, QSpinBox, QSplitter, QVBoxLayout, QWidget,
)

ZONE_COLORS = ["#0e7c86", "#d62728", "#2ca02c", "#9467bd", "#e377c2"]

#: the model actually fitted, shown above the build-up like TopSpin's report
FORMULAE = {
    "satrec": "I(t) = I₀·(1 − f·e<sup>−(t/T1)<sup>β</sup></sup>)",
    "invrec": "I(t) = I₀·(1 − 2f·e<sup>−(t/T1)<sup>β</sup></sup>)",
    "cpmg": "I(t) = I₀·e<sup>−(t/T2)<sup>β</sup></sup>",
    "t1rho": "I(t) = I₀·e<sup>−(t/T1ρ)<sup>β</sup></sup>",
}


class SatrecDialog(QDialog):
    def __init__(self, parent, expno: str | None):
        super().__init__(parent)
        self.setWindowTitle("Guided relaxation — T1 / T2")
        self.resize(980, 760)
        self.expno = expno
        self.x_ppm = None
        self.slices = None
        self.delays = None
        self.zones: list[pg.LinearRegionItem] = []
        self.keep: dict[int, np.ndarray] = {}   # per-zone kept-point mask
        self.results: dict[int, dict] = {}
        self._scatter: dict[int, pg.ScatterPlotItem] = {}

        v = QVBoxLayout(self)

        # ---- source + options bar ----
        top = QHBoxLayout()
        self.lbl = QLabel(expno or "no relaxation EXPNO selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Choose EXPNO…")
        btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(QLabel("kind"))
        self.kind = QComboBox()
        self.kind.addItems(["auto", "satrec", "invrec", "cpmg", "t1rho"])
        top.addWidget(self.kind)
        top.addWidget(QLabel("EM lb (Hz)"))
        self.lb = QDoubleSpinBox(); self.lb.setRange(0, 1e5); self.lb.setValue(100)
        top.addWidget(self.lb)
        self.magnitude = QCheckBox("magnitude")
        top.addWidget(self.magnitude)
        self.stretched = QCheckBox("stretched β")
        top.addWidget(self.stretched)
        top.addWidget(QLabel("slices"))
        self.firstSlice = QSpinBox(); self.firstSlice.setRange(1, 1)
        self.firstSlice.setToolTip("first slice to use (1-based)")
        top.addWidget(self.firstSlice)
        top.addWidget(QLabel("to"))
        self.lastSlice = QSpinBox(); self.lastSlice.setRange(1, 1)
        self.lastSlice.setToolTip("last slice to use — trim early-stopped or "
                                  "already-relaxed acquisitions")
        top.addWidget(self.lastSlice)
        self.firstSlice.valueChanged.connect(self._apply_slice_range)
        self.lastSlice.valueChanged.connect(self._apply_slice_range)
        self.btnProcess = QPushButton("Process slices")
        self.btnProcess.clicked.connect(self._process)
        top.addWidget(self.btnProcess)
        v.addLayout(top)
        self._all_delays = self._all_slices = None

        # ---- two stacked plots ----
        split = QSplitter(Qt.Vertical)

        spec_box = QWidget(); sb = QVBoxLayout(spec_box); sb.setContentsMargins(0, 0, 0, 0)
        zrow = QHBoxLayout()
        zrow.addWidget(QLabel("Integration zones:"))
        self.btnAddZone = QPushButton("+ Add zone")
        self.btnAddZone.clicked.connect(self._add_zone)
        self.btnClearZones = QPushButton("Clear")
        self.btnClearZones.clicked.connect(self._clear_zones)
        zrow.addWidget(self.btnAddZone); zrow.addWidget(self.btnClearZones)
        zrow.addWidget(QLabel("  (drag the shaded edges; the last, most-relaxed "
                              "slice is shown)"))
        zrow.addStretch(1)
        sb.addLayout(zrow)
        self.spec_plot = pg.PlotWidget(background="#fcfdfc")
        self.spec_plot.getPlotItem().invertX(True)
        self.spec_plot.setLabel("bottom", "shift", units="ppm")
        self.spec_plot.showGrid(x=True, y=True, alpha=0.1)
        sb.addWidget(self.spec_plot)
        split.addWidget(spec_box)

        build_box = QWidget(); bb = QVBoxLayout(build_box); bb.setContentsMargins(0, 0, 0, 0)
        brow = QHBoxLayout()
        brow.addWidget(QLabel("Build-up — click a point to exclude an outlier"))
        self.logx = QCheckBox("log delay")
        self.logx.toggled.connect(self._set_logx)
        brow.addWidget(self.logx)
        self.btnAuto = QPushButton("Fit view")
        self.btnAuto.setToolTip("rescale to the points; then zoom/pan freely")
        self.btnAuto.clicked.connect(lambda: self.build_plot.getViewBox()
                                     .autoRange(padding=0.06))
        brow.addWidget(self.btnAuto)
        brow.addStretch(1)
        self.btnFit = QPushButton("Fit T1 / T2")
        self.btnFit.setDefault(True)
        self.btnFit.clicked.connect(self._fit)
        self.btnCsv = QPushButton("Copy CSV"); self.btnCsv.setEnabled(False)
        self.btnCsv.clicked.connect(self._csv)
        brow.addWidget(self.btnFit); brow.addWidget(self.btnCsv)
        bb.addLayout(brow)
        self.build_plot = pg.PlotWidget(background="#fcfdfc")
        # NB: never pass units= on this axis. pyqtgraph's auto-SI-prefix turns
        # a log axis into nonsense ("delay (e27 s)"); keep the unit in the text.
        self.build_plot.setLabel("bottom", "delay (s)")
        self.build_plot.setLabel("left", "integral (norm.)")
        self.build_plot.showGrid(x=True, y=True, alpha=0.15)
        bb.addWidget(self.build_plot)
        split.addWidget(build_box)
        split.setSizes([380, 320])
        v.addWidget(split, 1)

        self.res = QLabel("")
        self.res.setStyleSheet("font-weight: 700; color: #0a5a62; font-size: 13px;")
        self.res.setWordWrap(True)
        v.addWidget(self.res)

        if expno:
            self._process()

    # ------------------------------------------------------------------
    def _pick(self):
        p = QFileDialog.getExistingDirectory(
            self, "Bruker EXPNO with ser + vdlist/vclist")
        if p:
            self.expno = p
            self.lbl.setText(p)
            self._process()

    def _kind(self) -> str | None:
        k = self.kind.currentText()
        return None if k == "auto" else k

    def _process(self):
        if not self.expno:
            return
        from larmor import satrec, series

        p = Path(self.expno)
        if not (p / "ser").exists():
            self.res.setText("this EXPNO has no 'ser' — not an arrayed "
                             "relaxation experiment")
            return
        self.btnProcess.setEnabled(False)
        self.res.setText("processing slices…")
        QApplication.processEvents()
        try:
            delays, self._src = series.read_delays(p)
            mode = "magnitude" if self.magnitude.isChecked() else "phase"
            self.x_ppm, slices = satrec.process_slices(
                p, lb_hz=self.lb.value(), mode=mode)
            n = min(len(delays), slices.shape[0])
            self._all_delays, self._all_slices = delays[:n], slices[:n]
        except Exception as exc:
            self.res.setText(f"failed: {exc}")
            self.btnProcess.setEnabled(True)
            return
        # let the user restrict the slice range (early-stopped / already-relaxed)
        self.firstSlice.blockSignals(True); self.lastSlice.blockSignals(True)
        self.firstSlice.setRange(1, n); self.lastSlice.setRange(1, n)
        self.firstSlice.setValue(1); self.lastSlice.setValue(n)
        self.firstSlice.blockSignals(False); self.lastSlice.blockSignals(False)
        self._apply_slice_range()
        self.btnProcess.setEnabled(True)

    def _apply_slice_range(self, *_):
        if self._all_slices is None:
            return
        lo = self.firstSlice.value() - 1
        hi = max(lo + 2, self.lastSlice.value())
        self.delays = self._all_delays[lo:hi]
        self.slices = self._all_slices[lo:hi]
        self.spec_plot.clear()
        self.spec_plot.plot(self.x_ppm, self.slices[-1],
                            pen=pg.mkPen("#1a2831", width=1.2))
        self._readd_zone_items()
        pk = self._robust_peak_ppm()
        span = self.x_ppm.max() - self.x_ppm.min()
        self.spec_plot.setXRange(pk + 0.18 * span, pk - 0.18 * span, padding=0)
        if not self.zones:
            self._add_zone()
        self.res.setText(
            f"{self.slices.shape[0]} of {self._all_slices.shape[0]} slices · "
            f"delays from {self._src} · drag the zone over your peak, then Fit")
        if self.results:
            self._fit()                  # refit on the new slice range

    def _robust_peak_ppm(self) -> float:
        """Peak of the relaxed slice after removing its rolling baseline and
        ignoring the noisy outer edges, so neither an edge spike nor a
        baseline offset steals the default zone."""
        from scipy.ndimage import uniform_filter1d

        ref = self.slices[-1].astype(float)
        ref = ref - np.median(ref)               # kill the baseline offset
        ref = np.maximum(ref, 0.0)               # phased spectrum: signal is +
        n = ref.size
        guard = max(3, n // 20)
        smooth = uniform_filter1d(ref, max(3, n // 500))
        smooth[:guard] = 0.0
        smooth[-guard:] = 0.0
        if smooth.max() <= 0:
            return float(self.x_ppm[n // 2])
        return float(self.x_ppm[int(np.argmax(smooth))])

    # ------------------------------------------------------------------
    def _add_zone(self):
        if self.x_ppm is None:
            return
        i = len(self.zones)
        col = ZONE_COLORS[i % len(ZONE_COLORS)]
        # default: a window around the robust peak of the relaxed slice
        pk = self._robust_peak_ppm()
        span = 0.06 * (self.x_ppm.max() - self.x_ppm.min())
        c = pg.mkColor(col); c.setAlpha(40)
        region = pg.LinearRegionItem(values=(pk - span, pk + span),
                                     brush=pg.mkBrush(c),
                                     pen=pg.mkPen(col, width=1.4))
        self.spec_plot.addItem(region)
        self.zones.append(region)

    def _clear_zones(self):
        for z in self.zones:
            self.spec_plot.removeItem(z)
        self.zones.clear()
        self.build_plot.clear()
        self.results.clear()
        self._scatter.clear()

    def _readd_zone_items(self):
        for z in self.zones:
            self.spec_plot.addItem(z)

    def _zone_ranges(self) -> list[tuple[float, float]]:
        out = []
        for z in self.zones:
            a, b = z.getRegion()
            out.append((max(a, b), min(a, b)))
        return out

    # ------------------------------------------------------------------
    def _fit(self):
        if self.slices is None or not self.zones:
            return
        from larmor import series

        zones = self._zone_ranges()
        integrals = series.integrate_zones(self.x_ppm, self.slices, zones)
        self.build_plot.clear()
        self._scatter.clear()
        self.results.clear()
        summaries = []
        for zi, ig in enumerate(integrals):
            col = ZONE_COLORS[zi % len(ZONE_COLORS)]
            if ig[-1] < 0:
                ig = -ig
            keep = self.keep.get(zi)
            if keep is None or len(keep) != len(ig):
                keep = np.ones(len(ig), bool)
                self.keep[zi] = keep
            norm = np.abs(ig).max() or 1.0
            self._plot_points(zi, ig / norm, col)
            try:
                r = series.fit_buildup(self.delays, ig, keep=keep,
                                       kind=self._kind() or "satrec",
                                       stretched=self.stretched.isChecked())
            except Exception as exc:
                summaries.append(f"zone {zi + 1}: {exc}")
                continue
            self.results[zi] = r
            if self.logx.isChecked():
                tt = np.logspace(
                    np.log10(max(self.delays[self.delays > 0].min(), 1e-4)),
                    np.log10(self.delays.max() * 1.3), 400)
            else:
                tt = np.linspace(0.0, self.delays.max() * 1.05, 400)
            self.build_plot.plot(tt, r["curve"](tt) / (r["norm"] / norm),
                                 pen=pg.mkPen(col, width=1.6))
            name = {"satrec": "T1", "invrec": "T1", "cpmg": "T2",
                    "t1rho": "T1ρ"}.get(r["kind"], "τ")
            err = f" ± {r['tau_err']:.3g}" if r["tau_err"] else ""
            b = (f", β={r['beta']:.2f}" if r["beta"] != 1.0 else "")
            summaries.append(f"<span style='color:{col}'>zone {zi + 1}</span>: "
                             f"{name} = {r['tau']:.4g}{err} s{b}")
        head = FORMULAE.get(self._kind() or "satrec", "")
        joined = "   ·   ".join(summaries)
        self.res.setText(f"<span style='color:#8a97a0'>{head}</span>&nbsp;&nbsp;"
                         f"&nbsp;&nbsp;{joined}" if head else joined)
        # rescale to the points, then leave the user free to zoom / pan
        self.build_plot.getViewBox().autoRange(padding=0.06)
        self.btnCsv.setEnabled(bool(self.results))

    def _set_logx(self, on: bool):
        self.build_plot.setLogMode(x=on, y=False)
        if self.results:
            self._fit()                       # redraw with the right t-sampling
        else:
            self.build_plot.getViewBox().autoRange(padding=0.06)

    def _plot_points(self, zi: int, y_norm: np.ndarray, col: str):
        keep = self.keep[zi]
        pos = self.delays > 0
        spots = []
        for k in range(len(self.delays)):
            if not pos[k]:
                continue
            on = keep[k]
            spots.append({"pos": (self.delays[k], y_norm[k]), "data": (zi, k),
                          "brush": pg.mkBrush(col) if on else pg.mkBrush(230, 230, 230),
                          "pen": pg.mkPen(col, width=1.4),
                          "size": 9 if on else 7,
                          "symbol": "o" if on else "x"})
        sc = pg.ScatterPlotItem(pxMode=True)
        sc.addPoints(spots)
        sc.sigClicked.connect(self._point_clicked)
        self.build_plot.addItem(sc)
        self._scatter[zi] = sc

    def _point_clicked(self, scatter, points):
        if not points:
            return
        zi, k = points[0].data()
        self.keep[zi][k] = not self.keep[zi][k]     # toggle exclusion
        self._fit()                                  # refit on the survivors

    def _csv(self):
        lines = ["zone,delay_s,integral_norm,kept"]
        zones = self._zone_ranges()
        from larmor import series

        integrals = series.integrate_zones(self.x_ppm, self.slices, zones)
        for zi, ig in enumerate(integrals):
            if ig[-1] < 0:
                ig = -ig
            norm = np.abs(ig).max() or 1.0
            keep = self.keep.get(zi, np.ones(len(ig), bool))
            for k in range(len(self.delays)):
                lines.append(f"{zi + 1},{self.delays[k]:g},"
                             f"{ig[k] / norm:.5f},{int(keep[k])}")
        for zi, r in self.results.items():
            lines.append(f"# zone {zi + 1}: tau={r['tau']:g} s "
                         f"({'stretched ' if r['beta'] != 1 else ''}{r['kind']})")
        QApplication.clipboard().setText("\n".join(lines))
        self.res.setText(self.res.text() + "   ·   copied CSV")
