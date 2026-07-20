"""Open FID… — ssNake-style pre-FT processing of a raw fid/ser.

Load a raw Bruker fid (1D) or ser (2D), build a processing pipeline
interactively (apodization, zero-fill, phase, for 2D the indirect quadrature
mode), watch the spectrum update live, and send the result to the fit
workbench (1D) or the 2D viewer.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QHBoxLayout,
    QLabel, QPushButton, QSpinBox, QVBoxLayout, QWidget,
)


class FidDialog(QDialog):
    #: emitted with (ppm, real_spectrum, meta) when the user accepts a 1D result
    accepted_1d = Signal(object, object, object)
    #: emitted with a twod.Data2D when the user accepts a 2D result
    accepted_2d = Signal(object)

    def __init__(self, parent, path: str | None):
        super().__init__(parent)
        self.setWindowTitle("Open FID — process before Fourier transform")
        self.resize(940, 660)
        self.data = None                # raw NMRData
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.lbl = QLabel(path or "no fid/ser selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Open fid / ser…")
        btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(btn)
        v.addLayout(top)

        # --- windowing / zero-fill ---
        w1 = QHBoxLayout()
        w1.addWidget(QLabel("window"))
        self.wdw = QComboBox()
        self.wdw.addItems(["none", "EM", "GM", "SINE", "QSINE", "TRAF"])
        self.wdw.setCurrentText("EM")
        w1.addWidget(self.wdw)
        w1.addWidget(QLabel("LB (Hz)"))
        self.lb = QDoubleSpinBox(); self.lb.setRange(-1e5, 1e5); self.lb.setValue(50)
        w1.addWidget(self.lb)
        w1.addWidget(QLabel("GB"))
        self.gb = QDoubleSpinBox(); self.gb.setRange(0.001, 1.0)
        self.gb.setDecimals(3); self.gb.setValue(0.1)
        w1.addWidget(self.gb)
        w1.addWidget(QLabel("ZF ×"))
        self.zf = QSpinBox(); self.zf.setRange(1, 16); self.zf.setValue(2)
        w1.addWidget(self.zf)
        w1.addWidget(QLabel("FCOR"))
        self.fcor = QDoubleSpinBox(); self.fcor.setRange(0.0, 2.0)
        self.fcor.setDecimals(2); self.fcor.setValue(0.5)
        w1.addWidget(self.fcor)
        w1.addStretch(1)
        v.addLayout(w1)

        # --- phase / 2D mode ---
        w2 = QHBoxLayout()
        w2.addWidget(QLabel("p0"))
        self.p0 = QDoubleSpinBox(); self.p0.setRange(-180, 180)
        w2.addWidget(self.p0)
        w2.addWidget(QLabel("p1"))
        self.p1 = QDoubleSpinBox(); self.p1.setRange(-720, 720)
        w2.addWidget(self.p1)
        self.btnAuto = QPushButton("Autophase")
        self.btnAuto.clicked.connect(self._autophase)
        w2.addWidget(self.btnAuto)
        self.mode_lbl = QLabel("F1 mode")
        self.mode = QComboBox()
        self.mode.addItems(["States", "States-TPPI", "Echo-Antiecho", "TPPI",
                            "QF"])
        w2.addWidget(self.mode_lbl); w2.addWidget(self.mode)
        w2.addStretch(1)
        v.addLayout(w2)

        self.plot = pg.PlotWidget(background="#fcfdfc")
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "shift", units="ppm")
        self.plot.showGrid(x=True, y=True, alpha=0.1)
        v.addWidget(self.plot, 1)

        bottom = QHBoxLayout()
        self.status = QLabel("")
        self.status.setStyleSheet("color: #37424a;")
        bottom.addWidget(self.status, 1)
        self.btnFt = QPushButton("Transform preview")
        self.btnFt.clicked.connect(self._preview)
        self.btnAccept = QPushButton("Use this spectrum →")
        self.btnAccept.setDefault(True)
        self.btnAccept.clicked.connect(self._accept)
        bottom.addWidget(self.btnFt)
        bottom.addWidget(self.btnAccept)
        v.addLayout(bottom)

        for w in (self.wdw, self.lb, self.gb, self.zf, self.fcor,
                  self.p0, self.p1, self.mode):
            sig = getattr(w, "currentTextChanged", None) or w.valueChanged
            sig.connect(self._preview)

        if path:
            self._load(path)
        self._set_2d_visibility(False)

    # ------------------------------------------------------------------
    def _set_2d_visibility(self, on: bool):
        self.mode_lbl.setVisible(on)
        self.mode.setVisible(on)

    def _pick(self):
        p, _ = QFileDialog.getOpenFileName(
            self, "Open raw Bruker fid or ser", "",
            "Bruker raw (fid ser);;All files (*)")
        if p:
            self._load(p)

    def _load(self, path):
        from larmor.io import bruker

        try:
            self.data = bruker.read(path)
        except Exception as exc:
            self.status.setText(f"failed: {exc}")
            return
        if self.data.domain != "time":
            self.status.setText("that is already a processed spectrum — use "
                                "File ▸ Open instead")
            self.data = None
            return
        self.lbl.setText(path)
        self._set_2d_visibility(self.data.ndim == 2)
        if self.data.ndim == 2:
            self.mode.setCurrentText(self.data.meta.get("fnmode", "States"))
        self.status.setText(self.data.summary + "  ·  " +
                            "; ".join(self.data.warnings))
        self._preview()

    def _ops(self) -> list[dict]:
        ops = []
        if self.fcor.value() != 1.0:
            ops.append({"op": "fcor", "factor": self.fcor.value()})
        w = self.wdw.currentText()
        if w == "EM" and self.lb.value():
            ops.append({"op": "em", "lb_hz": abs(self.lb.value())})
        elif w == "GM":
            ops.append({"op": "gm", "lb_hz": -abs(self.lb.value()),
                        "gb": self.gb.value()})
        elif w in ("SINE", "QSINE"):
            ops.append({"op": "sine", "ssb": 2,
                        "power": 2 if w == "QSINE" else 1})
        elif w == "TRAF":
            ops.append({"op": "traf", "lb_hz": abs(self.lb.value()) or 10.0})
        if self.zf.value() > 1:
            ops.append({"op": "zf", "factor": self.zf.value()})
        return ops

    def _preview(self):
        if self.data is None:
            return
        from larmor import fourier

        try:
            if self.data.ndim == 1:
                ops = self._ops()
                if self.p0.value() or self.p1.value():
                    ops.append({"op": "phase", "p0": self.p0.value(),
                                "p1": self.p1.value()})
                ppm, spec = fourier.ft1d(
                    self.data.data, self.data.axes[0].sw_Hz,
                    self.data.meta["larmor_MHz"], ops=ops)
                self._result = ("1d", ppm, spec)
                self.plot.clear()
                self.plot.plot(ppm, spec.real, pen=pg.mkPen("#1a2831", width=1.2))
            else:
                p = fourier.FT2DParams(
                    f2_ops=self._ops(), f1_ops=[{"op": "sine", "ssb": 2}],
                    mode=self.mode.currentText())
                d2 = fourier.ft2d_from_nmrdata(self.data, p)
                self._result = ("2d", d2)
                self.plot.clear()
                self.plot.plot(d2.f2_ppm, d2.projection("f2"),
                               pen=pg.mkPen("#1a2831", width=1.2))
                self.status.setText("2D transformed — F2 projection shown; "
                                    "'Use this spectrum' opens the 2D viewer")
        except Exception as exc:
            self.status.setText(f"transform failed: {exc}")

    def _autophase(self):
        if self.data is None or self.data.ndim != 1:
            return
        from larmor import fourier, processing as proc

        ppm, spec = fourier.ft1d(self.data.data, self.data.axes[0].sw_Hz,
                                 self.data.meta["larmor_MHz"], ops=self._ops())
        s = proc.Spectrum1D(x_ppm=ppm, y=spec, sfo1_MHz=1.0, sw_Hz=1.0)
        proc.op_autophase(s)
        # recover the p0 the scan found by projecting (approximate display)
        self.status.setText("autophased")
        self.plot.clear()
        self.plot.plot(ppm, s.y.real, pen=pg.mkPen("#1a2831", width=1.2))
        self._result = ("1d", ppm, s.y)

    def _accept(self):
        if not getattr(self, "_result", None):
            return
        if self._result[0] == "1d":
            _, ppm, spec = self._result
            self.accepted_1d.emit(ppm, spec.real, self.data.meta)
        else:
            self.accepted_2d.emit(self._result[1])
        self.accept()
