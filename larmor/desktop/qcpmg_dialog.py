"""QCPMG processing dialog (ssNake-style sum-echo workflow).

Turns a raw echo train into a FITTABLE spectrum:
  * split the train into echoes and read off the echo-top T2' decay;
  * sum the echoes (optionally T2-weighted) and whole-echo process them (swap
    the echo top to t=0) → a clean absorption lineshape you fit with the usual
    quadrupolar models — NOT the spikelet comb, which no smooth model can fit;
  * or view the spikelet spectrum.
Send the result to the fit workbench.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QHBoxLayout,
    QLabel, QPushButton, QSlider, QSpinBox, QSplitter, QVBoxLayout, QWidget,
)


class QcpmgDialog(QDialog):
    #: (ppm, amp, meta) of the processed spectrum to fit
    accepted_1d = Signal(object, object, dict)

    def __init__(self, parent, source: str | None):
        super().__init__(parent)
        self.setWindowTitle("QCPMG processing — sum echo / spikelets / T2")
        self.resize(1080, 720)
        self.source = source
        self.fid = None
        self.meta = {}
        self.T2 = None
        self._ppm = self._spec_raw = None

        v = QVBoxLayout(self)
        top = QHBoxLayout()
        self.lbl = QLabel(source or "no echo-train FID selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Open echo-train FID…"); btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(QLabel("period (pts)"))
        self.period = QSpinBox(); self.period.setRange(4, 1000000)
        self.period.valueChanged.connect(self._recompute)
        top.addWidget(self.period)
        top.addWidget(QLabel("echo top"))
        self.top = QSpinBox(); self.top.setRange(0, 1000000)
        self.top.valueChanged.connect(self._recompute)
        top.addWidget(self.top)
        top.addWidget(btn)
        v.addLayout(top)

        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("mode"))
        self.mode = QComboBox()
        self.mode.addItems(["sum echo (absorption — fit this)", "spikelets"])
        self.mode.currentIndexChanged.connect(self._recompute)
        ctl.addWidget(self.mode)
        self.t2w = QCheckBox("T2 weighting")
        self.t2w.setToolTip("weight each echo by exp(-t/T2) before summing "
                            "(matched filter → best S/N)")
        self.t2w.toggled.connect(self._recompute)
        ctl.addWidget(self.t2w)
        ctl.addWidget(QLabel("GB (Hz)"))
        self.gb = QDoubleSpinBox(); self.gb.setRange(0, 1e5); self.gb.setValue(200)
        self.gb.setToolTip("Gaussian line broadening")
        self.gb.valueChanged.connect(self._recompute)
        ctl.addWidget(self.gb)
        ctl.addWidget(QLabel("p0"))
        self.p0 = QDoubleSpinBox(); self.p0.setRange(-720, 720); self.p0.setDecimals(2)
        self.p0.setWrapping(True); self.p0.valueChanged.connect(self._rephase)
        ctl.addWidget(self.p0)
        ctl.addWidget(QLabel("p1"))
        self.p1 = QDoubleSpinBox(); self.p1.setRange(-36000, 36000); self.p1.setDecimals(2)
        self.p1.valueChanged.connect(self._rephase)
        ctl.addWidget(self.p1)
        ctl.addWidget(QLabel("step"))
        self.pstep = QDoubleSpinBox(); self.pstep.setRange(0.01, 90.0)
        self.pstep.setDecimals(2); self.pstep.setValue(1.0)
        self.pstep.setToolTip("phase increment per click/scroll — lower it for "
                              "fine control (the arrows/scroll move by this much)")
        self.pstep.valueChanged.connect(self._set_pstep)
        ctl.addWidget(self.pstep)
        self.btnAuto = QPushButton("Autophase"); self.btnAuto.clicked.connect(self._autophase)
        ctl.addWidget(self.btnAuto)
        self.btnHelp = QPushButton("Help")
        self.btnHelp.setToolTip("open the QCPMG processing guide")
        self.btnHelp.clicked.connect(self._help)
        ctl.addWidget(self.btnHelp)
        ctl.addStretch(1)
        v.addLayout(ctl)
        self._set_pstep()

        split = QSplitter(Qt.Horizontal)
        # left: echo train + T2 decay
        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(QLabel("echo train (magnitude); red lines = period"))
        self.p_fid = pg.PlotWidget(background="#fcfdfc")
        self.p_fid.setLabel("bottom", "point"); self.p_fid.setMaximumHeight(180)
        lv.addWidget(self.p_fid)
        lv.addWidget(QLabel("echo-top decay → T2 (the evolution time)"))
        self.p_t2 = pg.PlotWidget(background="#fcfdfc")
        self.p_t2.setLabel("bottom", "time", units="s")
        lv.addWidget(self.p_t2)
        self.t2lbl = QLabel(""); self.t2lbl.setStyleSheet(
            "font-weight: 700; color: #6a4fb0;")
        lv.addWidget(self.t2lbl)
        split.addWidget(left)
        # right: spectrum
        right = QWidget(); rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0)
        rv.addWidget(QLabel("processed spectrum"))
        self.p_spec = pg.PlotWidget(background="#fcfdfc")
        self.p_spec.getPlotItem().invertX(True)
        self.p_spec.setLabel("bottom", "shift", units="ppm")
        rv.addWidget(self.p_spec)
        split.addWidget(right)
        split.setSizes([420, 640])
        v.addWidget(split, 1)

        bot = QHBoxLayout()
        self.info = QLabel(""); self.info.setStyleSheet(
            "color: #0a5a62; font-weight: 600;")
        bot.addWidget(self.info, 1)
        self.btnSend = QPushButton("Send to fit →"); self.btnSend.setDefault(True)
        self.btnSend.clicked.connect(self._send)
        bot.addWidget(self.btnSend)
        v.addLayout(bot)

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
        p = max(4, qcpmg.detect_period(self.fid))
        ech = qcpmg.split_echoes(self.fid, p)
        for w, val in ((self.period, p), (self.top, qcpmg.echo_top_point(ech))):
            w.blockSignals(True); w.setValue(val); w.blockSignals(False)
        self._fit_t2()
        self._recompute()

    def _carrier_ppm(self) -> float:
        bf1 = self.meta.get("bf1_MHz") or 0.0
        return (self.meta.get("o1_Hz", 0.0) / bf1) if bf1 else 0.0

    def _fit_t2(self):
        from larmor import qcpmg

        if self.fid is None:
            return
        p = self.period.value()
        ech = qcpmg.split_echoes(self.fid, p)
        top = min(self.top.value(), p - 1)
        tau = p / (self.meta.get("sw_Hz", 1.0) or 1.0)
        decay = qcpmg.echo_decay(ech, top)
        self.T2, curve = qcpmg.fit_t2(tau, decay)
        t = np.arange(decay.size) * tau
        self.p_t2.clear()
        self.p_t2.plot(t, decay / (decay.max() or 1.0), pen=None, symbol="o",
                       symbolSize=5, symbolBrush="#6a4fb0")
        tt = np.linspace(0, t.max() if t.size else 1, 300)
        self.p_t2.plot(tt, curve(tt), pen=pg.mkPen("#6a4fb0", width=1.6))
        lb = 1.0 / (np.pi * self.T2) if self.T2 else 0.0
        self.t2lbl.setText(
            f"T2 = {self.T2 * 1e3:.2f} ms  →  matched apodization "
            f"LB = {lb:.0f} Hz  (= 1/πT2)    "
            f"(τecho = {tau * 1e3:.3f} ms, {ech.shape[0]} echoes)")

    def _draw_train(self):
        self.p_fid.clear()
        mag = np.abs(self.fid)
        self.p_fid.plot(np.arange(mag.size), mag, pen=pg.mkPen("#1a2831"))
        p = self.period.value()
        for k in range(1, max(1, min(30, mag.size // p))):
            self.p_fid.addItem(pg.InfiniteLine(pos=k * p, angle=90,
                               pen=pg.mkPen("#d62728", width=0.6, style=Qt.DotLine)))

    # ------------------------------------------------------------------
    def _recompute(self, *_):
        """Rebuild the (unphased) spectrum from the current split/mode/weight."""
        if self.fid is None:
            return
        from larmor import qcpmg

        sw = self.meta.get("sw_Hz", 0.0); sfo = self.meta.get("larmor_MHz", 0.0)
        carrier = self._carrier_ppm(); p = self.period.value()
        self._fit_t2(); self._draw_train()
        if self.mode.currentIndex() == 0:
            t2w = self.T2 if self.t2w.isChecked() else None
            self._ppm, self._spec_raw = qcpmg.sum_echo_spectrum(
                self.fid, p, sw, sfo, carrier, top=self.top.value(),
                t2_weight_s=t2w, gb_Hz=self.gb.value())
            self._kind = "sum echo (absorption)"
        else:
            self._ppm, self._spec_raw = qcpmg.spikelet_spectrum(
                self.fid, sw, sfo, carrier, lb_Hz=max(self.gb.value(), 1.0))
            self._kind = "spikelets"
        self._rephase()

    def _rephase(self, *_):
        if self._spec_raw is None:
            return
        from larmor import qcpmg

        n = self._spec_raw.size
        ramp = np.arange(n) / n
        ph = np.exp(-1j * (np.deg2rad(self.p0.value())
                           + np.deg2rad(self.p1.value()) * ramp))
        self._spec = np.real(self._spec_raw * ph)
        self.p_spec.clear()
        self.p_spec.plot(self._ppm, self._spec, pen=pg.mkPen("#1a2831"))
        sw = self.meta.get("sw_Hz", 0.0); sfo = self.meta.get("larmor_MHz", 0.0)
        spc = qcpmg.spikelet_spacing_ppm(self.period.value(), sw, sfo)
        self.info.setText(
            f"{self.meta.get('nucleus', '?')} · period {self.period.value()} pts "
            f"· spikelet spacing {spc:.1f} ppm · {self._kind} · "
            f"p0 {self.p0.value()}° p1 {self.p1.value()}°")

    def _set_pstep(self, *_):
        s = self.pstep.value()
        self.p0.setSingleStep(s); self.p1.setSingleStep(s)

    def _autophase(self):
        from larmor import qcpmg

        if self._spec_raw is None:
            return
        p0, p1 = qcpmg.autophase(self._spec_raw)
        for w, val in ((self.p0, p0), (self.p1, p1)):
            w.blockSignals(True); w.setValue(val); w.blockSignals(False)
        self._rephase()

    def _help(self):
        from larmor.desktop.help_dialog import show_help

        show_help(self, "qcpmg", "QCPMG processing guide")

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
