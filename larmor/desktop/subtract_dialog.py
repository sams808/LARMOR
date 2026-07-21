"""Background subtraction: sample − scale·(another spectrum).

Load a background/reference spectrum, scale (and optionally shift) it to match
the sample, preview the difference live, and send the result to the workbench.
Unlike the 'background spectrum' fit component (which fits the amplitude), this
produces a new spectrum you can save and reuse.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QDialog, QDoubleSpinBox, QFileDialog, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QVBoxLayout,
)


class SubtractDialog(QDialog):
    applied = Signal(object, object)          # (ppm, amp) of the difference

    def __init__(self, parent, ppm, amp, meta: dict):
        super().__init__(parent)
        self.setWindowTitle("Subtract a spectrum (background)")
        self.resize(880, 560)
        self.ppm = np.asarray(ppm, float)
        self.amp = np.asarray(amp, float)
        self.meta = meta
        self.bg_ppm = self.bg_amp = None

        v = QVBoxLayout(self)
        top = QHBoxLayout()
        self.lbl = QLabel("no background loaded")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Load background spectrum…")
        btn.clicked.connect(self._load)
        top.addWidget(self.lbl, 1); top.addWidget(btn)
        v.addLayout(top)

        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("scale"))
        self.scale = QDoubleSpinBox(); self.scale.setRange(-1e6, 1e6)
        self.scale.setDecimals(4); self.scale.setValue(1.0)
        self.scale.setSingleStep(0.05)
        self.scale.valueChanged.connect(self._redraw)
        ctl.addWidget(self.scale)
        self.btnAuto = QPushButton("Auto scale (this view)")
        self.btnAuto.setToolTip("least-squares scale over the current x-range")
        self.btnAuto.clicked.connect(self._auto)
        ctl.addWidget(self.btnAuto)
        ctl.addWidget(QLabel("shift (ppm)"))
        self.shift = QDoubleSpinBox(); self.shift.setRange(-1e5, 1e5)
        self.shift.setDecimals(3); self.shift.valueChanged.connect(self._redraw)
        ctl.addWidget(self.shift)
        ctl.addStretch(1)
        v.addLayout(ctl)

        self.plot = pg.PlotWidget(background="#fcfdfc")
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "shift", units="ppm")
        self.plot.addLegend(offset=(10, 10))
        v.addWidget(self.plot, 1)

        bot = QHBoxLayout()
        self.info = QLabel(""); self.info.setStyleSheet("color: #5a6871;")
        bot.addWidget(self.info, 1)
        self.btnApply = QPushButton("Apply → workbench"); self.btnApply.setDefault(True)
        self.btnApply.setEnabled(False); self.btnApply.clicked.connect(self._apply)
        bot.addWidget(self.btnApply)
        v.addLayout(bot)
        self._redraw()

    # ------------------------------------------------------------------
    def _load(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Background spectrum", "",
            "Spectra (*.fxmla *.json 1r *.csv *.txt);;All files (*)")
        if not path:
            return
        try:
            from larmor.desktop.app import _load_any

            bp, ba, *_ = _load_any(path)
        except Exception:
            try:
                from larmor.io import bruker

                d = bruker.read(path)
                if d.ndim != 1 or d.domain != "freq":
                    raise ValueError("not a 1D spectrum")
                bp, ba = np.asarray(d.axes[0].values), np.asarray(d.data, float)
            except Exception as exc:
                QMessageBox.warning(self, "Background", f"cannot read: {exc}")
                return
        self.bg_ppm, self.bg_amp = np.asarray(bp, float), np.asarray(ba, float)
        self.lbl.setText(Path(path).name)
        self.btnApply.setEnabled(True)
        self._redraw()

    def _difference(self):
        from larmor.io import spectra

        return spectra.subtract(self.ppm, self.amp, self.bg_ppm, self.bg_amp,
                                scale=self.scale.value(),
                                shift_ppm=self.shift.value())

    def _auto(self):
        if self.bg_ppm is None:
            return
        from larmor.io import spectra

        (x0, x1), _ = self.plot.getPlotItem().getViewBox().viewRange()
        s = spectra.best_scale(self.ppm, self.amp, self.bg_ppm, self.bg_amp,
                               shift_ppm=self.shift.value(), window=(x0, x1))
        self.scale.setValue(s)

    def _redraw(self):
        self.plot.clear()
        self.plot.plot(self.ppm, self.amp, pen=pg.mkPen("#1a2831", width=1.4),
                       name="sample")
        if self.bg_ppm is not None:
            from larmor.io import spectra

            bg_on = spectra.subtract(self.ppm, np.zeros_like(self.amp),
                                     self.bg_ppm, -self.bg_amp,
                                     scale=self.scale.value(),
                                     shift_ppm=self.shift.value())
            self.plot.plot(self.ppm, bg_on, pen=pg.mkPen("#e8832a", width=1.1),
                           name="scaled background")
            diff = self._difference()
            self.plot.plot(self.ppm, diff, pen=pg.mkPen("#0e7c86", width=1.4),
                           name="difference")
            self.info.setText(f"difference = sample − {self.scale.value():.4g} × "
                              f"background (shift {self.shift.value():+.3g} ppm)")

    def _apply(self):
        if self.bg_ppm is None:
            return
        self.applied.emit(self.ppm, self._difference())
        self.accept()
