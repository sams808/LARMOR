"""Embeddable 2D contour view for the main window's central stack.

Shows any 2D dataset (real MQMAS/DQ-SQ or a pseudo-2D relaxation array) as a
contour map with F1/F2 projections and a diagonal, and brings the everyday
TopSpin / ssNake 2D gestures into LARMOR:

  * contour levels, and a positive / negative / both sign choice;
  * interactive **2D phasing** the TopSpin way -- click 1 or 2 reference peaks,
    the selected rows (F2) or columns (F1) are shown full-width as 1D traces
    with the pivot marked, phase them with p0/p1 and ±90/180° steps, then apply
    to every row/column;
  * a manual **shear** for MQMAS data processed elsewhere;
  * a live cursor readout (F2, F1, intensity);
  * pulling a 1D trace out for fitting (F2 skyline/sum, or the row under a
    draggable cursor).
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox, QDoubleSpinBox, QHBoxLayout, QLabel, QMenu, QPushButton, QSlider,
    QStackedWidget, QToolButton, QVBoxLayout, QWidget,
)


class Contour2DView(QWidget):
    #: (ppm, amp, label) of a 1D trace the user wants to fit
    slice_to_fit = Signal(object, object, str)
    #: (f2_ppm, f1_ppm) where the user clicked to place a 2D site
    add_requested = Signal(float, float)
    #: ask the app to load a 1D spectrum to overlay on a projection ("f2"/"f1")
    load_1d_for_projection = Signal(str)

    def __init__(self):
        super().__init__()
        self.data = None            # currently displayed (committed) data
        self._orig = None           # as loaded, never mutated
        self._committed = None      # phases the user has applied so far
        self._add_model = None      # model name while placing a 2D site
        self._model = None          # fitted overlay: (z, f2_ppm, f1_ppm)
        #: HMQC: an external 1D per projection axis, {"f2":(ppm,amp),"f1":...}
        self._hmqc = {"f2": None, "f1": None}
        self._hmqc_scale = {"f2": 1.0, "f1": 1.0}
        self._hmqc_color = {"f2": "#e8832a", "f1": "#6a4fb0"}
        self._picks: list[int] = []  # reference row/column indices (1 or 2)
        self._pick_axis = "f2"
        self._pivot = None          # pivot ppm on the phased axis
        self._pref = []             # per-pick (coords, raw trace, curve, plot)
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)

        # ---- contour / display bar ----
        bar = QHBoxLayout()
        self.title = QLabel("2D dataset")
        self.title.setStyleSheet("font-weight: 600; color: #16202a;")
        bar.addWidget(self.title, 1)
        bar.addWidget(QLabel("contours"))
        self.sign = QComboBox()
        self.sign.addItems(["positive", "negative", "both"])
        self.sign.currentTextChanged.connect(self._redraw)
        bar.addWidget(self.sign)
        bar.addWidget(QLabel("levels"))
        self.nlevels = QDoubleSpinBox()
        self.nlevels.setRange(3, 40); self.nlevels.setValue(12)
        self.nlevels.valueChanged.connect(self._redraw)
        bar.addWidget(self.nlevels)
        bar.addWidget(QLabel("floor ×σ"))
        self.floor = QDoubleSpinBox()
        self.floor.setRange(0.5, 40.0); self.floor.setSingleStep(0.5)
        self.floor.setValue(8.0)
        self.floor.valueChanged.connect(self._redraw)
        bar.addWidget(self.floor)
        bar.addWidget(QLabel("shear"))
        self.shearv = QDoubleSpinBox()
        self.shearv.setRange(-3.0, 3.0); self.shearv.setDecimals(3)
        self.shearv.setToolTip("manual shear F1' = F1 + k·F2 (0 = none; MQMAS "
                               "from mrsimulator is already sheared)")
        bar.addWidget(self.shearv)
        self.btnShear = QPushButton("Apply shear")
        self.btnShear.clicked.connect(self._apply_shear)
        bar.addWidget(self.btnShear)
        self.btnPhase = QPushButton("Phase 2D")
        self.btnPhase.setCheckable(True)
        self.btnPhase.toggled.connect(self._toggle_phase)
        bar.addWidget(self.btnPhase)
        self.btnCal = QPushButton("Calibrate")
        self.btnCal.setCheckable(True)
        self.btnCal.setToolTip("click a peak, then set its F2 / F1 ppm")
        bar.addWidget(self.btnCal)
        self.btnMeasure = QPushButton("Measure")
        self.btnMeasure.setCheckable(True)
        self.btnMeasure.setToolTip("two draggable markers: ΔF2 / ΔF1 in ppm and Hz")
        self.btnMeasure.toggled.connect(self._toggle_measure)
        bar.addWidget(self.btnMeasure)
        self.btnHmqc = QPushButton("HMQC")
        self.btnHmqc.setCheckable(True)
        self.btnHmqc.setToolTip("overlay 1D spectra on the projections and "
                                "extract un-correlated features")
        self.btnHmqc.toggled.connect(lambda on: self.hmqcbar.setVisible(on))
        bar.addWidget(self.btnHmqc)
        self.btnOps = QToolButton(); self.btnOps.setText("2D ops ▾")
        self.btnOps.setPopupMode(QToolButton.InstantPopup)
        ops = QMenu(self.btnOps)
        ops.addAction("Transpose F1↔F2", lambda: self._op("transpose"))
        ops.addAction("Reverse F2", lambda: self._op("rev_f2"))
        ops.addAction("Reverse F1", lambda: self._op("rev_f1"))
        ops.addAction("Symmetrize (diagonal)", lambda: self._op("symmetrize"))
        ops.addSeparator()
        ops.addAction("Diagonal → fit", lambda: self._op("diagonal"))
        ops.addAction("F2 skyline → fit", lambda: self._emit_projection("skyline"))
        ops.addAction("F2 sum → fit", lambda: self._emit_projection("sum"))
        self.btnOps.setMenu(ops)
        bar.addWidget(self.btnOps)
        v.addLayout(bar)
        self._mtargets: list = []

        # ---- HMQC bar (hidden until HMQC) ----
        self.hmqcbar = QWidget()
        hb = QHBoxLayout(self.hmqcbar); hb.setContentsMargins(0, 0, 0, 0)
        hb.addWidget(QLabel("overlay 1D on projection:"))
        b2 = QPushButton("F2 1D…"); b2.clicked.connect(
            lambda: self.load_1d_for_projection.emit("f2"))
        b1 = QPushButton("F1 1D…"); b1.clicked.connect(
            lambda: self.load_1d_for_projection.emit("f1"))
        hb.addWidget(b2); hb.addWidget(b1)
        hb.addWidget(QLabel("scale"))
        self.hmqcScale = QDoubleSpinBox()
        self.hmqcScale.setRange(0.0, 1e6); self.hmqcScale.setDecimals(3)
        self.hmqcScale.setValue(1.0)
        self.hmqcScale.setToolTip("scale the correlated (projection) part before "
                                  "subtracting it from the 1D")
        self.hmqcScale.valueChanged.connect(self._on_hmqc_scale)
        hb.addWidget(self.hmqcScale)
        self.btnFitScale = QPushButton("fit scale")
        self.btnFitScale.setToolTip("least-squares scale of the projection to the "
                                    "1D over the visible range (robust; not just "
                                    "peak-match)")
        self.btnFitScale.clicked.connect(self._fit_scale)
        hb.addWidget(self.btnFitScale)
        self.btnUncF2 = QPushButton("uncorrelated F2 →")
        self.btnUncF2.setToolTip("1D − scale·(F2 projection) → features that do "
                                 "NOT correlate, sent to the fit workbench")
        self.btnUncF2.clicked.connect(lambda: self._emit_uncorrelated("f2"))
        self.btnUncF1 = QPushButton("uncorrelated F1 →")
        self.btnUncF1.clicked.connect(lambda: self._emit_uncorrelated("f1"))
        hb.addWidget(self.btnUncF2); hb.addWidget(self.btnUncF1)
        hb.addStretch(1)
        self.hmqcbar.setVisible(False)
        v.addWidget(self.hmqcbar)

        # ---- phase bar (hidden until Phase 2D) ----
        self.phasebar = QWidget()
        pb = QHBoxLayout(self.phasebar); pb.setContentsMargins(0, 0, 0, 0)
        pb.addWidget(QLabel("of"))
        self.pdim = QComboBox()
        self.pdim.addItems(["rows (F2)", "columns (F1)"])
        self.pdim.currentIndexChanged.connect(self._clear_picks)
        pb.addWidget(self.pdim)
        pb.addWidget(QLabel("using"))
        self.npoints = QComboBox()
        self.npoints.addItems(["1 point", "2 points"])
        self.npoints.currentIndexChanged.connect(self._clear_picks)
        pb.addWidget(self.npoints)
        self.btnGoPhase = QPushButton("Phase ▶")
        self.btnGoPhase.setToolTip("show the selected rows/columns and phase them")
        self.btnGoPhase.clicked.connect(self._enter_phasing)
        pb.addWidget(self.btnGoPhase)
        pb.addWidget(QLabel("p0"))
        self.p0 = QSlider(Qt.Horizontal); self.p0.setRange(-180, 180)
        self.p0.setMaximumWidth(120)
        self.p0v = QDoubleSpinBox(); self.p0v.setRange(-180, 180)
        self.p0.valueChanged.connect(self.p0v.setValue)
        self.p0v.valueChanged.connect(lambda x: self.p0.setValue(int(x)))
        self.p0v.valueChanged.connect(self._phase_changed)
        pb.addWidget(self.p0); pb.addWidget(self.p0v)
        pb.addWidget(QLabel("p1"))
        self.p1 = QSlider(Qt.Horizontal); self.p1.setRange(-720, 720)
        self.p1.setMaximumWidth(120)
        self.p1v = QDoubleSpinBox(); self.p1v.setRange(-720, 720)
        self.p1.valueChanged.connect(self.p1v.setValue)
        self.p1v.valueChanged.connect(lambda x: self.p1.setValue(int(x)))
        self.p1v.valueChanged.connect(self._phase_changed)
        pb.addWidget(self.p1); pb.addWidget(self.p1v)
        for lbl, d in (("−90°", -90.0), ("+90°", 90.0), ("180°", 180.0)):
            b = QPushButton(lbl); b.setMaximumWidth(52)
            b.clicked.connect(lambda _=False, d=d: self._nudge_p0(d))
            pb.addWidget(b)
        self.btnPhaseApply = QPushButton("Apply")
        self.btnPhaseApply.clicked.connect(self._apply_phase)
        self.btnRepick = QPushButton("Re-pick")
        self.btnRepick.clicked.connect(self._repick)
        self.btnPhaseReset = QPushButton("Reset")
        self.btnPhaseReset.clicked.connect(self._reset_phase)
        pb.addWidget(self.btnPhaseApply); pb.addWidget(self.btnRepick)
        pb.addWidget(self.btnPhaseReset)
        self.phasehint = QLabel("")
        self.phasehint.setStyleSheet("color: #0e7c86;")
        pb.addWidget(self.phasehint); pb.addStretch(1)
        self.phasebar.setVisible(False)
        v.addWidget(self.phasebar)

        # ---- stacked display: contour (0) and phasing traces (1) ----
        self.stack = QStackedWidget()
        # contour with top (F2) and left (F1) projections
        self.glw = pg.GraphicsLayoutWidget()
        self.glw.setBackground("#fcfdfc")
        self.p_top = self.glw.addPlot(row=0, col=1)
        self.p_top.setMaximumHeight(80); self.p_top.hideAxis("bottom")
        self.p_main = self.glw.addPlot(row=1, col=1)
        self.p_left = self.glw.addPlot(row=1, col=0)
        self.p_left.setMaximumWidth(80); self.p_left.hideAxis("left")
        self.p_main.setLabel("bottom", "F2 (ppm)")
        self.p_main.setLabel("right", "F1")
        self.p_main.getViewBox().invertX(True)
        self.p_main.getViewBox().invertY(True)
        # the projections must share the main plot's inversion, else they mirror
        # and pan the opposite way to the contour
        self.p_top.getViewBox().invertX(True)
        self.p_left.getViewBox().invertY(True)
        self.p_top.setXLink(self.p_main)
        self.p_left.setYLink(self.p_main)
        self.p_main.scene().sigMouseClicked.connect(self._on_click)
        self.p_main.scene().sigMouseMoved.connect(self._on_move)
        self.stack.addWidget(self.glw)
        # phasing traces
        self.phase_glw = pg.GraphicsLayoutWidget()
        self.phase_glw.setBackground("#fcfdfc")
        self.stack.addWidget(self.phase_glw)
        v.addWidget(self.stack, 1)
        self._slice_line = None

        # ---- pick-to-fit row ----
        pick = QHBoxLayout()
        pick.addWidget(QLabel("Send a 1D trace to fitting:"))
        self.btnSky = QPushButton("F2 skyline →")
        self.btnSky.clicked.connect(lambda: self._emit_projection("skyline"))
        self.btnSum = QPushButton("F2 sum →")
        self.btnSum.clicked.connect(lambda: self._emit_projection("sum"))
        self.btnRow = QPushButton("row at cursor →")
        self.btnRow.setToolTip("drag the dashed horizontal line to a row, then "
                               "send that F2 slice to the fit workbench")
        self.btnRow.clicked.connect(self._emit_row)
        pick.addWidget(self.btnSky)
        pick.addWidget(self.btnSum)
        pick.addWidget(self.btnRow)
        pick.addStretch(1)
        self.cursor = QLabel("")
        self.cursor.setStyleSheet("color: #5a6871; font-family: Consolas, monospace;")
        pick.addWidget(self.cursor)
        self.hint = QLabel("")
        self.hint.setStyleSheet("color: #5a6871;")
        pick.addWidget(self.hint)
        v.addLayout(pick)

    # ------------------------------------------------------------------
    def set_data(self, data, title: str = ""):
        d = data.normalized() if hasattr(data, "normalized") else data
        self._orig = d
        self._committed = d
        self.data = d
        self.title.setText(title or "2D dataset")
        f1_kind = "arrayed (relaxation)" if getattr(data, "notes", None) and \
            any("pseudo" in n or "arrayed" in n for n in data.notes) else "F1"
        self.hint.setText(f"F1 = {f1_kind}")
        # reset phase/pick state WITHOUT triggering extra redraws (each contour
        # pass over a large 2D is expensive; we redraw exactly once below)
        self.btnPhase.blockSignals(True); self.btnPhase.setChecked(False)
        self.btnPhase.blockSignals(False)
        self.phasebar.setVisible(False)
        self._hmqc = {"f2": None, "f1": None}
        self._hmqc_scale = {"f2": 1.0, "f1": 1.0}
        self._picks = []; self._pivot = None; self._pref = []
        for s in (self.p0v, self.p1v):
            s.blockSignals(True); s.setValue(0); s.blockSignals(False)
        self.stack.setCurrentWidget(self.glw)
        self._redraw()

    # ---------- workspace state (for the workspace manager) ----------
    def get_state(self) -> dict:
        return {
            "orig": self._orig, "committed": self._committed, "data": self.data,
            "hmqc": dict(self._hmqc), "hmqc_scale": dict(self._hmqc_scale),
            "hmqc_color": dict(self._hmqc_color), "model": self._model,
            "sign": self.sign.currentText(),
            "nlevels": self.nlevels.value(), "floor": self.floor.value(),
            "title": self.title.text(),
        }

    def set_state(self, s: dict):
        self._orig = s["orig"]; self._committed = s["committed"]
        self.data = s["data"]
        self._hmqc = dict(s["hmqc"]); self._hmqc_scale = dict(s["hmqc_scale"])
        self._hmqc_color = dict(s["hmqc_color"]); self._model = s["model"]
        self.title.setText(s["title"])
        for w, val in ((self.nlevels, s["nlevels"]), (self.floor, s["floor"])):
            w.blockSignals(True); w.setValue(val); w.blockSignals(False)
        self.sign.blockSignals(True); self.sign.setCurrentText(s["sign"])
        self.sign.blockSignals(False)
        # drop any transient interaction state, then draw
        self.btnPhase.blockSignals(True); self.btnPhase.setChecked(False)
        self.btnPhase.blockSignals(False); self.phasebar.setVisible(False)
        for b in (self.btnCal, self.btnMeasure, self.btnHmqc):
            b.blockSignals(True); b.setChecked(False); b.blockSignals(False)
        self.hmqcbar.setVisible(False)
        self._toggle_measure(False)
        self._picks = []; self._pivot = None; self._pref = []
        self.stack.setCurrentWidget(self.glw)
        self._redraw()

    # ---------- phasing: pick 1-2 peaks on the contour ----------
    def _axis(self) -> str:
        return "f2" if self.pdim.currentIndex() == 0 else "f1"

    def _max_picks(self) -> int:
        return 1 if self.npoints.currentIndex() == 0 else 2

    def _toggle_phase(self, on: bool):
        self.phasebar.setVisible(on)
        self.stack.setCurrentWidget(self.glw)     # always show the map here
        if on:
            self.phasehint.setText("click %d reference peak(s), then Phase ▶"
                                   % self._max_picks())
        self._clear_picks()

    def _clear_picks(self, *_):
        self._picks = []
        self._pick_axis = self._axis()
        self._pivot = None
        for s in (self.p0v, self.p1v):
            s.blockSignals(True); s.setValue(0); s.blockSignals(False)
        if self.btnPhase.isChecked():
            self.stack.setCurrentWidget(self.glw)
            self.phasehint.setText("click %d reference peak(s), then Phase ▶"
                                   % self._max_picks())
        if self.data is not None:
            self._redraw()

    def _repick(self):
        self._picks = []
        self.stack.setCurrentWidget(self.glw)
        self.phasehint.setText("click %d reference peak(s), then Phase ▶"
                               % self._max_picks())
        self._redraw()

    def _on_click(self, ev):
        if self.data is None or ev.button() != Qt.LeftButton:
            return
        if self.btnCal.isChecked():
            self._calibrate_at(self.p_main.getViewBox().mapSceneToView(
                ev.scenePos()))
            ev.accept()
            return
        if self._add_model is not None:
            p = self.p_main.getViewBox().mapSceneToView(ev.scenePos())
            self.add_requested.emit(float(p.x()), float(p.y()))
            ev.accept()
            return
        if not self.btnPhase.isChecked():
            return
        if self.stack.currentWidget() is not self.glw:
            return
        p = self.p_main.getViewBox().mapSceneToView(ev.scenePos())
        d = self.data
        self._pick_axis = self._axis()
        if self._pick_axis == "f2":
            self._pivot = float(p.x())          # pivot on F2
            idx = int(np.argmin(np.abs(d.f1_ppm - p.y())))   # which row
        else:
            self._pivot = float(p.y())          # pivot on F1
            idx = int(np.argmin(np.abs(d.f2_ppm - p.x())))   # which column
        self._picks.append(idx)
        self._picks = self._picks[-self._max_picks():]
        lab = "row" if self._pick_axis == "f2" else "column"
        self.phasehint.setText(
            f"{len(self._picks)}/{self._max_picks()} {lab}(s) picked · "
            f"pivot={self._pivot:.1f} ppm — Phase ▶ when ready")
        ev.accept()
        self._redraw()

    def _enter_phasing(self):
        if self.data is None or not self._picks:
            return
        c = self._committed
        coords = c.f2_ppm if self._pick_axis == "f2" else c.f1_ppm
        self.phase_glw.clear()
        self._pref = []
        exact = c.has_hyper
        for k, idx in enumerate(self._picks):
            raw = c.z[idx] if self._pick_axis == "f2" else c.z[:, idx]
            plot = self.phase_glw.addPlot(row=k, col=0)
            plot.showGrid(x=True, y=True, alpha=0.12)
            plot.getViewBox().invertX(True)
            other = (c.f1_ppm if self._pick_axis == "f2" else c.f2_ppm)[idx]
            lab = "Row" if self._pick_axis == "f2" else "Col"
            tag = "" if exact else "  (Hilbert)"
            plot.setTitle(f"<span style='color:#1f6feb'>{lab} {idx} "
                          f"/ {other:.1f} ppm{tag}</span>")
            curve = plot.plot(coords, np.asarray(raw, float),
                              pen=pg.mkPen("#1f6feb", width=1.2))
            plot.addItem(pg.InfiniteLine(pos=self._pivot, angle=90,
                         pen=pg.mkPen("#d62728", width=1.2)))
            if k == len(self._picks) - 1:
                plot.setLabel("bottom", "ppm")
            self._pref.append((coords, idx, curve))
        self.stack.setCurrentWidget(self.phase_glw)
        self.phasehint.setText(f"pivot={self._pivot:.1f} ppm · "
                               f"ph0={self.p0v.value():.0f} ph1={self.p1v.value():.0f}")

    def _nudge_p0(self, delta: float):
        x = self.p0v.value() + delta
        while x > 180.0:
            x -= 360.0
        while x < -180.0:
            x += 360.0
        self.p0v.setValue(x)

    def _phase_changed(self, *_):
        if not self._pref or self.stack.currentWidget() is not self.phase_glw:
            return
        for coords, idx, curve in self._pref:
            y = self._committed.phase_line(self._pick_axis, idx,
                                           self.p0v.value(), self.p1v.value(),
                                           self._pivot)
            curve.setData(coords, y)
        self.phasehint.setText(f"pivot={self._pivot:.1f} ppm · "
                               f"ph0={self.p0v.value():.0f} ph1={self.p1v.value():.0f}")

    def _apply_phase(self):
        if self._committed is None or not self._pref:
            return
        self._committed = self._committed.phased(
            self._pick_axis, self.p0v.value(), self.p1v.value(), self._pivot)
        self.data = self._committed
        self._picks = []
        for s in (self.p0v, self.p1v):
            s.blockSignals(True); s.setValue(0); s.blockSignals(False)
        self.stack.setCurrentWidget(self.glw)
        self.phasehint.setText("phase applied — pick again or uncheck Phase 2D")
        self._redraw()

    def _reset_phase(self):
        self._committed = self._orig
        self.data = self._orig
        self._picks = []
        for s in (self.p0v, self.p1v):
            s.blockSignals(True); s.setValue(0); s.blockSignals(False)
        self.stack.setCurrentWidget(self.glw)
        self._redraw()

    def _apply_shear(self):
        if self._committed is None or self.shearv.value() == 0:
            return
        from larmor import twod

        self._orig = twod.shear(self._committed, self.shearv.value())
        self._committed = self._orig
        self.data = self._orig
        self._redraw()

    # ---------- 2D operations ----------
    def _op(self, kind: str):
        if self.data is None:
            return
        if kind == "diagonal":
            g, amp = self._committed.diagonal()
            self.slice_to_fit.emit(np.asarray(g), np.asarray(amp),
                                   f"diagonal of {self.title.text()}")
            return
        if kind == "transpose":
            new = self._committed.transposed()
        elif kind == "rev_f2":
            new = self._committed.reversed_axis("f2")
        elif kind == "rev_f1":
            new = self._committed.reversed_axis("f1")
        elif kind == "symmetrize":
            new = self._committed.symmetrized()
        else:
            return
        self._orig = self._committed = self.data = new
        self._redraw()

    # ---------- calibrate + measure ----------
    def _calibrate_at(self, p):
        from PySide6.QtWidgets import QInputDialog

        f2, ok = QInputDialog.getDouble(
            self, "Calibrate F2", f"Set F2 {p.x():.2f} ppm to:",
            float(p.x()), -1e5, 1e5, 3)
        if not ok:
            return
        f1, ok = QInputDialog.getDouble(
            self, "Calibrate F1", f"Set F1 {p.y():.2f} ppm to:",
            float(p.y()), -1e5, 1e5, 3)
        if not ok:
            return
        self._shift_axes(float(f2) - float(p.x()), float(f1) - float(p.y()))
        self.btnCal.setChecked(False)

    def _shift_axes(self, d2: float, d1: float):
        from larmor.twod import Data2D

        def sh(src):
            return Data2D(src.f2_ppm + d2, src.f1_ppm + d1, src.z, src.nucleus,
                          src.larmor_MHz, src.spin_rate_Hz, src.source,
                          list(src.notes), src.ri, src.ir, src.ii)
        self._orig = sh(self._orig)
        self._committed = sh(self._committed)
        self.data = sh(self.data)
        self._redraw()

    def _toggle_measure(self, on: bool):
        for t in self._mtargets:
            self.p_main.removeItem(t)
        self._mtargets = []
        if not on or self.data is None:
            return
        d = self.data
        for fx, fy in ((0.4, 0.55), (0.6, 0.45)):
            x = float(d.f2_ppm.min() + fx * (d.f2_ppm.max() - d.f2_ppm.min()))
            y = float(d.f1_ppm.min() + fy * (d.f1_ppm.max() - d.f1_ppm.min()))
            t = pg.TargetItem(pos=(x, y), size=12, movable=True,
                              pen=pg.mkPen("#0e7c86", width=1.6))
            t.sigPositionChanged.connect(self._emit_measure2d)
            self.p_main.addItem(t)
            self._mtargets.append(t)
        self._emit_measure2d()

    def _emit_measure2d(self, *_):
        if len(self._mtargets) != 2 or self.data is None:
            return
        a, b = (t.pos() for t in self._mtargets)
        d2 = abs(a.x() - b.x()); d1 = abs(a.y() - b.y())
        lar = self.data.larmor_MHz or 0.0
        hz2 = f" ({d2 * lar:.0f} Hz)" if lar else ""
        self.cursor.setText(f"ΔF2 {d2:.2f} ppm{hz2}   ΔF1 {d1:.2f} ppm")

    # ------------------------------------------------------------------
    def _on_move(self, pos):
        if self.btnMeasure.isChecked():
            return                       # keep the Δ readout stable while measuring
        if self.data is None or not self.p_main.sceneBoundingRect().contains(pos):
            return
        p = self.p_main.getViewBox().mapSceneToView(pos)
        d = self.data
        j = int(np.argmin(np.abs(d.f2_ppm - p.x())))
        i = int(np.argmin(np.abs(d.f1_ppm - p.y())))
        self.cursor.setText(f"F2 {p.x():7.1f}  F1 {p.y():7.1f}  z {d.z[i, j]:+.3f}")

    # ---------- fitted-model overlay ----------
    def set_add_mode(self, name):
        self._add_model = name
        self.setCursor(Qt.CrossCursor if name else Qt.ArrowCursor)

    def set_model(self, z, f2_ppm, f1_ppm):
        self._model = (np.asarray(z, float), np.asarray(f2_ppm),
                       np.asarray(f1_ppm)) if z is not None else None
        self._redraw()

    def clear_model(self):
        self._model = None
        self._redraw()

    # ---------- HMQC: overlay 1D on projections + uncorrelated features ----------
    def set_projection_1d(self, axis: str, ppm, amp, color: str | None = None):
        """Store an external 1D for the F2 or F1 projection and peak-match its
        scale to the HMQC sum-projection."""
        ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
        o = np.argsort(ppm)
        self._hmqc[axis] = (ppm[o], amp[o])
        if color:
            self._hmqc_color[axis] = color
        self._hmqc_scale[axis] = self._peak_match(axis)
        self.hmqcScale.blockSignals(True); self.hmqcScale.setValue(1.0)
        self.hmqcScale.blockSignals(False)
        if not self.btnHmqc.isChecked():
            self.btnHmqc.setChecked(True)
        self._redraw()

    def _proj(self, axis: str):
        d = self.data
        coords = d.f2_ppm if axis == "f2" else d.f1_ppm
        return coords, d.projection(axis, "sum")

    def _one_d_on(self, axis: str, coords):
        px, pa = self._hmqc[axis]
        return np.interp(coords, px, pa, left=0.0, right=0.0)

    def _peak_match(self, axis: str) -> float:
        if self.data is None or self._hmqc[axis] is None:
            return 1.0
        coords, proj = self._proj(axis)
        yi = self._one_d_on(axis, coords)
        i = int(np.argmax(proj))
        return float(yi[i] / (proj[i] or 1.0)) if proj[i] else 1.0

    def _eff_scale(self, axis: str) -> float:
        return self._hmqc_scale[axis] * (self.hmqcScale.value() or 1.0)

    def _on_hmqc_scale(self, *_):
        self._redraw()

    def _fit_scale(self):
        """Least-squares scale of each loaded projection to its 1D over the
        visible F2 range — robust when the strongest peak isn't fully correlated."""
        (x0, x1), _ = self.p_main.getViewBox().viewRange()
        lo, hi = min(x0, x1), max(x0, x1)
        for axis in ("f2", "f1"):
            if self._hmqc[axis] is None:
                continue
            coords, proj = self._proj(axis)
            yi = self._one_d_on(axis, coords)
            m = (coords >= lo) & (coords <= hi)
            denom = float(proj[m] @ proj[m])
            if denom > 0:
                self._hmqc_scale[axis] = max(0.0, float(yi[m] @ proj[m]) / denom)
        self.hmqcScale.blockSignals(True); self.hmqcScale.setValue(1.0)
        self.hmqcScale.blockSignals(False)
        self._redraw()

    def _emit_uncorrelated(self, axis: str):
        if self.data is None or self._hmqc[axis] is None:
            return
        coords, proj = self._proj(axis)
        yi = self._one_d_on(axis, coords)
        diff = yi - self._eff_scale(axis) * proj
        self.slice_to_fit.emit(np.asarray(coords), np.asarray(diff),
                               f"un-correlated {axis.upper()} of {self.title.text()}")

    # ------------------------------------------------------------------
    def _contour_levels(self, z):
        edge = max(1, min(z.shape) // 20)
        frame = np.concatenate([z[:edge].ravel(), z[-edge:].ravel(),
                                z[:, :edge].ravel(), z[:, -edge:].ravel()])
        floor = max(self.floor.value() * float(frame.std() or 1e-6), 0.01)
        n = int(self.nlevels.value())
        top = float(np.nanmax(np.abs(z))) or 1.0
        if floor >= top:
            floor = top * 0.1
        return np.logspace(np.log10(floor), np.log10(top), n)

    @staticmethod
    def _decimate(z, f2, f1, cap: int = 512):
        """Stride a large 2D down so marching-squares stays cheap; contours look
        the same but cost a fraction. Returns (zc, f2c, f1c)."""
        s2 = max(1, -(-z.shape[1] // cap))       # ceil division -> under cap
        s1 = max(1, -(-z.shape[0] // cap))
        if s1 == 1 and s2 == 1:
            return z, f2, f1
        return z[::s1, ::s2], f2[::s2], f1[::s1]

    @staticmethod
    def _iso_transform(zc, f2c, f1c):
        tr = pg.QtGui.QTransform()
        tr.translate(f2c[0], f1c[0])
        tr.scale((f2c[-1] - f2c[0]) / max(zc.shape[1] - 1, 1),
                 (f1c[-1] - f1c[0]) / max(zc.shape[0] - 1, 1))
        return tr

    def _redraw(self):
        if self.data is None:
            return
        d = self.data
        for p in (self.p_main, self.p_top, self.p_left):
            p.clear()
        z = d.z
        f2, f1 = d.f2_ppm, d.f1_ppm
        zc, f2c, f1c = self._decimate(z, f2, f1)
        levels = self._contour_levels(zc)
        tr = self._iso_transform(zc, f2c, f1c)

        sign = self.sign.currentText()
        pos_cmap = pg.colormap.get("viridis").getLookupTable(0.0, 1.0, len(levels))
        zt = np.ascontiguousarray(zc.T)
        if sign in ("positive", "both"):
            for lvl, col in zip(levels, pos_cmap):
                iso = pg.IsocurveItem(data=zt, level=lvl,
                                      pen=pg.mkPen(tuple(int(c) for c in col), width=1))
                iso.setTransform(tr); self.p_main.addItem(iso)
        if sign in ("negative", "both"):
            for lvl in levels:
                iso = pg.IsocurveItem(data=zt, level=-lvl,
                                      pen=pg.mkPen("#d62728", width=1))
                iso.setTransform(tr); self.p_main.addItem(iso)

        # fitted-model contours overlaid in orange, dashed
        if self._model is not None:
            mz, mf2, mf1 = self._decimate(*self._model)
            mtop = float(np.nanmax(mz)) or 1.0
            mlev = np.logspace(np.log10(max(0.05 * mtop, 1e-6)),
                               np.log10(mtop), max(4, int(self.nlevels.value()) // 2))
            mtr = self._iso_transform(mz, mf2, mf1)
            mzt = np.ascontiguousarray(mz.T)
            for lvl in mlev:
                iso = pg.IsocurveItem(data=mzt, level=lvl,
                                      pen=pg.mkPen("#e8832a", width=1,
                                                   style=Qt.DashLine))
                iso.setTransform(mtr); self.p_main.addItem(iso)

        lo = max(min(f2), min(f1)); hi = min(max(f2), max(f1))
        self.p_main.plot([lo, hi], [lo, hi],
                         pen=pg.mkPen("#b9c1bc", style=Qt.DashLine))
        # projections — or, in HMQC mode, the scaled sum-projection (correlated)
        # under the external 1D (total), so the gap is the uncorrelated part
        if self._hmqc["f2"] is not None:
            self.p_top.plot(f2, self._eff_scale("f2") * d.projection("f2", "sum"),
                            pen=pg.mkPen("#0e7c86"))
            self.p_top.plot(f2, self._one_d_on("f2", f2),
                            pen=pg.mkPen(self._hmqc_color["f2"], width=1.2))
        else:
            self.p_top.plot(f2, z.max(axis=0), pen=pg.mkPen("#0e7c86"))
        if self._hmqc["f1"] is not None:
            self.p_left.plot(self._eff_scale("f1") * d.projection("f1", "sum"), f1,
                             pen=pg.mkPen("#0e7c86"))
            self.p_left.plot(self._one_d_on("f1", f1), f1,
                             pen=pg.mkPen(self._hmqc_color["f1"], width=1.2))
        else:
            self.p_left.plot(z.max(axis=1), f1, pen=pg.mkPen("#0e7c86"))

        # reference-peak markers while picking for phase
        if self.btnPhase.isChecked() and self._picks:
            for idx in self._picks:
                if self._pick_axis == "f2":
                    self.p_main.addItem(pg.InfiniteLine(
                        pos=float(f1[idx]), angle=0,
                        pen=pg.mkPen("#1f6feb", width=1.2)))
                else:
                    self.p_main.addItem(pg.InfiniteLine(
                        pos=float(f2[idx]), angle=90,
                        pen=pg.mkPen("#1f6feb", width=1.2)))
            if self._pivot is not None:
                ang = 90 if self._pick_axis == "f2" else 0
                self.p_main.addItem(pg.InfiniteLine(
                    pos=self._pivot, angle=ang,
                    pen=pg.mkPen("#d62728", width=1.2, style=Qt.DotLine)))

        self._slice_line = pg.InfiniteLine(
            pos=float(f1[len(f1) // 2]), angle=0, movable=True,
            pen=pg.mkPen("#8a97a0", width=1.1, style=Qt.DashLine))
        self.p_main.addItem(self._slice_line)

    # ------------------------------------------------------------------
    def _emit_projection(self, mode: str):
        if self.data is None:
            return
        d = self.data
        proj = (d.z.max(axis=0) if mode == "skyline" else d.z.sum(axis=0))
        self.slice_to_fit.emit(np.asarray(d.f2_ppm), np.asarray(proj),
                               f"F2 {mode} of {self.title.text()}")

    def _emit_row(self):
        if self.data is None or self._slice_line is None:
            return
        d = self.data
        y = float(self._slice_line.value())
        i = int(np.argmin(np.abs(d.f1_ppm - y)))
        self.slice_to_fit.emit(np.asarray(d.f2_ppm), np.asarray(d.z[i]),
                               f"F1={d.f1_ppm[i]:.1f} slice of {self.title.text()}")
