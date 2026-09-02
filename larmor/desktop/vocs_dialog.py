"""Stitch a frequency-stepped (VOCS) acquisition into one wideline spectrum.

Add the sub-spectra (one EXPNO per transmitter offset), pick skyline or
average, trim the excitation roll-off, watch the live preview, send the
stitched spectrum to the workbench. The combining itself is larmor.vocs
(Qt-free, tested); this dialog is only the loading and the preview.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QHBoxLayout, QLabel,
    QListWidget, QListWidgetItem, QMessageBox, QPushButton, QVBoxLayout,
)

from larmor.desktop import theme

#: preview palette for the pieces (stitch result draws in the accent colour)
_PIECE_COLORS = ["#7f8fa6", "#8c7ae6", "#4a9d7f", "#b8860b", "#a55b6b",
                 "#5b7fa5", "#7a9a55", "#9a6b8f"]


class VocsDialog(QDialog):
    #: (ppm, amp, note_lines, source_paths)
    applied = Signal(object, object, object, object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Stitch frequency-stepped (VOCS) spectra")
        self.resize(960, 620)
        #: [(path, ppm, amp, meta_summary)]
        self.pieces: list[tuple] = []

        v = QVBoxLayout(self)
        top = QHBoxLayout()
        self.btnAdd = QPushButton("Add sub-spectra…")
        self.btnAdd.clicked.connect(self._add)
        self.btnRemove = QPushButton("Remove selected")
        self.btnRemove.clicked.connect(self._remove)
        top.addWidget(self.btnAdd)
        top.addWidget(self.btnRemove)
        top.addStretch(1)
        top.addWidget(QLabel("combine:"))
        self.method = QComboBox()
        self.method.addItem("skyline (max — the VOCS standard)", "skyline")
        self.method.addItem("average (coverage-weighted mean)", "average")
        self.method.currentIndexChanged.connect(self._recompute)
        top.addWidget(self.method)
        top.addWidget(QLabel("trim per edge:"))
        self.trim = QDoubleSpinBox()
        self.trim.setRange(0.0, 45.0)
        self.trim.setValue(10.0)
        self.trim.setSuffix(" %")
        self.trim.setSingleStep(2.5)
        self.trim.setToolTip(
            "Excitation rolls off toward each piece's edges; trim past the "
            "roll-off, especially for 'average' (a rolled-off edge pulls the "
            "mean down at every junction).")
        self.trim.valueChanged.connect(self._recompute)
        top.addWidget(self.trim)
        v.addLayout(top)

        self.list = QListWidget()
        self.list.setMaximumHeight(110)
        v.addWidget(self.list)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "chemical shift", units="ppm")
        self.plot.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
        v.addWidget(self.plot, 1)

        self.status = QLabel("add two or more sub-spectra to begin")
        self.status.setStyleSheet(f"color: {theme.active().text_dim};")
        v.addWidget(self.status)

        bottom = QHBoxLayout()
        bottom.addStretch(1)
        self.btnApply = QPushButton("→ Workbench")
        self.btnApply.setEnabled(False)
        self.btnApply.clicked.connect(self._apply)
        btnClose = QPushButton("Close")
        btnClose.clicked.connect(self.reject)
        bottom.addWidget(self.btnApply)
        bottom.addWidget(btnClose)
        v.addLayout(bottom)

        self._result = None

    # ------------------------------------------------------------ loading
    def _add(self):
        from larmor.desktop.app import _load_any

        paths, _ = QFileDialog.getOpenFileNames(
            self, "Sub-spectra (one per transmitter offset)", "",
            "Spectra (*.fxmla *.json 1r *.txt *.csv *.dat);;All files (*)")
        for path in paths:
            try:
                ppm, amp, recipe, *_ = _load_any(path)
            except Exception as exc:
                QMessageBox.warning(self, "VOCS stitch",
                                    f"cannot read {path}: {exc}")
                continue
            ppm = np.asarray(ppm, float)
            amp = np.asarray(amp, float)
            nuc = (recipe or {}).get("nucleus", "")
            lar = (recipe or {}).get("larmor_frequency_MHz", 0.0)
            self.pieces.append((path, ppm, amp, f"{nuc} {lar:.2f} MHz"))
        self._refill()

    def _remove(self):
        row = self.list.currentRow()
        if 0 <= row < len(self.pieces):
            self.pieces.pop(row)
            self._refill()

    def _refill(self):
        self.list.clear()
        for i, (path, ppm, _amp, meta) in enumerate(self.pieces):
            col = _PIECE_COLORS[i % len(_PIECE_COLORS)]
            it = QListWidgetItem(
                f"{Path(path).name}   [{ppm.min():.0f} … {ppm.max():.0f} ppm]"
                f"   {meta}")
            it.setForeground(pg.mkColor(col))
            self.list.addItem(it)
        self._recompute()

    # ------------------------------------------------------------ stitching
    def _recompute(self, *_):
        from larmor.vocs import stitch

        self.plot.clear()
        self._result = None
        self.btnApply.setEnabled(False)
        for i, (_path, ppm, amp, _meta) in enumerate(self.pieces):
            self.plot.plot(ppm, amp,
                           pen=pg.mkPen(_PIECE_COLORS[i % len(_PIECE_COLORS)],
                                        width=1.0, style=Qt.DotLine))
        if len(self.pieces) < 2:
            self.status.setText("add two or more sub-spectra to begin")
            return
        try:
            res = stitch([(p[1], p[2]) for p in self.pieces],
                         method=self.method.currentData(),
                         trim_frac=self.trim.value() / 100.0,
                         labels=[Path(p[0]).name for p in self.pieces])
        except ValueError as exc:
            self.status.setText(str(exc))
            return
        self.plot.plot(res.ppm, res.amp,
                       pen=pg.mkPen(theme.active().accent, width=1.8))
        self._result = res
        warn = [n for n in res.notes if n.startswith("WARNING")]
        self.status.setText(warn[0] if warn else res.notes[0])
        self.btnApply.setEnabled(True)

    def _apply(self):
        if self._result is None:
            return
        res = self._result
        self.applied.emit(res.ppm, res.amp, list(res.notes),
                          [p[0] for p in self.pieces])
        self.accept()
