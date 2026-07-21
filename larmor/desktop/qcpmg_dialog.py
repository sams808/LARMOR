"""QCPMG processing dialog: turn a raw echo train into a spectrum, then send it
to the fit workbench. Coadded envelope (continuous, fittable) or spikelets."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QLabel, QPushButton, QSpinBox, QSplitter, QVBoxLayout, QWidget,
)


class QcpmgDialog(QDialog):
    #: (ppm, amp, meta) of the processed spectrum to fit
    accepted_1d = Signal(object, object, dict)

    def __init__(self, parent, source: str | None):
        super().__init__(parent)
        self.setWindowTitle("QCPMG processing")
        self.resize(940, 680)
        self.source = source
        self.fid = None
        self.meta = {}

        v = QVBoxLayout(self)
        top = QHBoxLayout()
        self.lbl = QLabel(source or "no echo-train FID selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Open echo-train FID…")
        btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1); top.addWidget(btn)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("echo period (pts)"))
        self.period = QSpinBox(); self.period.setRange(4, 100000)
        self.period.setToolTip("auto-detected from the echo-train autocorrelation")
        opts.addWidget(self.period)
        opts.addWidget(QLabel("mode"))
        self.mode = QComboBox()
        self.mode.addItems(["coadded envelope (fit)", "spikelets"])
        self.mode.currentIndexChanged.connect(self._process)
        opts.addWidget(self.mode)
        opts.addWidget(QLabel("LB (Hz)"))
        self.lb = QDoubleSpinBox(); self.lb.setRange(0, 1e5); self.lb.setValue(100)
        opts.addWidget(self.lb)
        self.drop1 = QCheckBox("drop 1st echo"); self.drop1.setChecked(True)
        opts.addWidget(self.drop1)
        self.btnProc = QPushButton("Process")
        self.btnProc.clicked.connect(self._process)
        opts.addWidget(self.btnProc)
        opts.addStretch(1)
        v.addLayout(opts)

        split = QSplitter(Qt.Vertical)
        tw = QWidget(); tv = QVBoxLayout(tw); tv.setContentsMargins(0, 0, 0, 0)
        tv.addWidget(QLabel("echo train (magnitude) — red lines mark the period"))
        self.p_fid = pg.PlotWidget(background="#fcfdfc")
        self.p_fid.setLabel("bottom", "point")
        tv.addWidget(self.p_fid)
        split.addWidget(tw)
        sw = QWidget(); sv = QVBoxLayout(sw); sv.setContentsMargins(0, 0, 0, 0)
        sv.addWidget(QLabel("processed spectrum"))
        self.p_spec = pg.PlotWidget(background="#fcfdfc")
        self.p_spec.getPlotItem().invertX(True)
        self.p_spec.setLabel("bottom", "shift", units="ppm")
        sv.addWidget(self.p_spec)
        split.addWidget(sw)
        split.setSizes([300, 340])
        v.addWidget(split, 1)

        bot = QHBoxLayout()
        self.info = QLabel("")
        self.info.setStyleSheet("color: #0a5a62; font-weight: 600;")
        bot.addWidget(self.info, 1)
        self.btnSend = QPushButton("Send to fit →")
        self.btnSend.setDefault(True)
        self.btnSend.clicked.connect(self._send)
        bot.addWidget(self.btnSend)
        v.addLayout(bot)

        self._ppm = self._spec = None
        if source:
            self._load(source)

    # ------------------------------------------------------------------
    def _pick(self):
        p, _ = QFileDialog.getOpenFileName(self, "Echo-train FID (fid) or EXPNO")
        if not p:
            p = QFileDialog.getExistingDirectory(self, "EXPNO with a QCPMG fid")
        if p:
            self._load(p)

    def _load(self, source):
        from larmor import qcpmg
        from larmor.io import bruker

        try:
            d = bruker.read(source)
            if d.domain != "time" or d.ndim != 1:
                raise ValueError("not a 1D echo-train FID (open the raw fid)")
        except Exception as exc:
            self.info.setText(f"cannot read: {exc}")
            return
        self.source = source
        self.lbl.setText(str(source))
        self.fid = np.asarray(d.data, complex)
        self.meta = dict(d.meta)
        self.period.setValue(max(4, qcpmg.detect_period(self.fid)))
        self._draw_train()
        self._process()

    def _carrier_ppm(self) -> float:
        bf1 = self.meta.get("bf1_MHz") or 0.0
        return (self.meta.get("o1_Hz", 0.0) / bf1) if bf1 else 0.0

    def _draw_train(self):
        self.p_fid.clear()
        mag = np.abs(self.fid)
        self.p_fid.plot(np.arange(mag.size), mag, pen=pg.mkPen("#1a2831"))
        p = self.period.value()
        for k in range(1, max(1, mag.size // p)):
            self.p_fid.addItem(pg.InfiniteLine(pos=k * p, angle=90,
                               pen=pg.mkPen("#d62728", width=0.6,
                                            style=Qt.DotLine)))

    def _process(self):
        if self.fid is None:
            return
        from larmor import qcpmg

        sw = self.meta.get("sw_Hz", 0.0)
        sfo = self.meta.get("larmor_MHz", 0.0)
        carrier = self._carrier_ppm()
        p = self.period.value()
        if self.mode.currentIndex() == 0:
            ppm, spec = qcpmg.coadd_spectrum(
                self.fid, p, sw, sfo, carrier, lb_Hz=self.lb.value(),
                drop_first=1 if self.drop1.isChecked() else 0)
            kind = "coadded envelope"
        else:
            ppm, cspec = qcpmg.spikelet_spectrum(
                self.fid, sw, sfo, carrier, lb_Hz=self.lb.value())
            spec = np.abs(cspec)
            kind = "spikelets"
        self._ppm, self._spec = ppm, spec
        self.p_spec.clear()
        self.p_spec.plot(ppm, spec, pen=pg.mkPen("#1a2831"))
        self._draw_train()
        self.info.setText(
            f"{self.meta.get('nucleus', '?')} · period {p} pts · spikelet "
            f"spacing {qcpmg.spikelet_spacing_ppm(p, sw, sfo):.1f} ppm · {kind}")

    def _send(self):
        if self._ppm is None:
            return
        meta = {
            "expno": self.meta.get("expno", ""),
            "title": (self.meta.get("title", "").splitlines() or ["QCPMG"])[0]
            + " (QCPMG)",
            "nucleus": self.meta.get("nucleus", ""),
            "larmor_MHz": self.meta.get("larmor_MHz", 0.0),
            "masr_Hz": self.meta.get("masr_Hz"),
        }
        self.accepted_1d.emit(np.asarray(self._ppm), np.asarray(self._spec), meta)
        self.accept()
