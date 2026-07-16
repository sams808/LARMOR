"""Tool dialogs: REDOR, DFT (.magres) import, Errors Analysis."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QLabel, QLineEdit, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)


class RedorDialog(QDialog):
    def __init__(self, parent, expno: str | None):
        super().__init__(parent)
        self.setWindowTitle("REDOR — dipolar coupling & distance")
        self.resize(760, 520)
        self.expno = expno
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.lbl = QLabel(expno or "no EXPNO selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Choose EXPNO…")
        btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(btn)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("observed"))
        self.iso1 = QLineEdit("13C"); self.iso1.setFixedWidth(60)
        opts.addWidget(self.iso1)
        opts.addWidget(QLabel("dephased by"))
        self.iso2 = QLineEdit("15N"); self.iso2.setFixedWidth(60)
        opts.addWidget(self.iso2)
        opts.addWidget(QLabel("regime"))
        self.regime = QComboBox()
        self.regime.addItems(["auto", "short", "pair"])
        opts.addWidget(self.regime)
        self.btnRun = QPushButton("Analyze")
        self.btnRun.setDefault(True)
        self.btnRun.clicked.connect(self._run)
        opts.addWidget(self.btnRun)
        opts.addStretch(1)
        v.addLayout(opts)

        self.plot = pg.PlotWidget(background="w")
        self.plot.setLabel("bottom", "recoupling time / s")
        self.plot.setLabel("left", "ΔS/S₀")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        v.addWidget(self.plot, 1)

        self.res = QLabel("")
        self.res.setStyleSheet("font-weight: 700; color: #0a5a62; font-size: 14px;")
        self.res.setWordWrap(True)
        v.addWidget(self.res)

    def _pick(self):
        p = QFileDialog.getExistingDirectory(self, "EXPNO with redor.txt")
        if p:
            self.expno = p
            self.lbl.setText(p)

    def _run(self):
        if not self.expno:
            return
        from larmor import redor

        pair = (self.iso1.text().strip(), self.iso2.text().strip())
        try:
            res = redor.analyze_expno(self.expno, pair=pair,
                                      regime=self.regime.currentText())
        except Exception as exc:
            self.res.setText(f"failed: {exc}")
            return
        self.plot.clear()
        self.plot.plot(res.ntr_s, res.ds_s0, pen=None, symbol="o", symbolSize=7,
                       symbolBrush=None, symbolPen=pg.mkPen("#0e7c86", width=1.5))
        tt = np.linspace(res.ntr_s.min(), res.ntr_s.max(), 200)
        self.plot.plot(tt, res.curve(tt), pen=pg.mkPen("#d62728", width=1.6))
        self.res.setText(res.summary + "   ·   " + " · ".join(res.notes))


class MagresDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("Import DFT tensors (.magres)")
        self.resize(720, 480)
        self.result_sites: list[dict] = []
        self._sites = []
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.lbl = QLabel("no file")
        btn = QPushButton("Open .magres…")
        btn.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(btn)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("isotope"))
        self.iso = QComboBox()
        opts.addWidget(self.iso)
        opts.addWidget(QLabel("model"))
        self.model = QComboBox()
        from larmor import models as reg

        self.model.addItems([m["name"] for m in reg.describe_all()
                             if m["name"] in ("quad_ct", "quad_csa", "czjzek")])
        opts.addWidget(self.model)
        opts.addWidget(QLabel("σ_ref (ppm)"))
        self.ref = QDoubleSpinBox(); self.ref.setRange(-1e4, 1e4)
        self.ref.setToolTip("reference shielding to convert to chemical shift")
        opts.addWidget(self.ref)
        opts.addStretch(1)
        v.addLayout(opts)

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(
            ["site", "isotope", "Cq (MHz) / η", "σ_iso (ppm)"])
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table, 1)

        row = QHBoxLayout()
        self.note = QLabel("")
        self.note.setStyleSheet("color: #93a0a8;")
        row.addWidget(self.note, 1)
        add = QPushButton("Add these sites to the fit")
        add.setDefault(True)
        add.clicked.connect(self._accept)
        row.addWidget(add)
        v.addLayout(row)
        self.iso.currentTextChanged.connect(self._refresh)

    def _pick(self):
        p, _ = QFileDialog.getOpenFileName(self, "CASTEP/QE .magres",
                                           "", "magres (*.magres);;All (*)")
        if not p:
            return
        from larmor import dft

        try:
            self._sites = dft.read_magres(p)
            warnings = dft.assign_isotopes(self._sites)
        except Exception as exc:
            self.note.setText(f"failed: {exc}")
            return
        self.lbl.setText(Path(p).name)
        isos = sorted({s.isotope for s in self._sites if s.isotope})
        self.iso.clear()
        self.iso.addItems(isos)
        self.note.setText(" · ".join(warnings) if warnings else
                          f"{len(self._sites)} sites, {len(isos)} isotopes")
        self._refresh()

    def _refresh(self):
        from larmor import dft

        iso = self.iso.currentText()
        sites = dft.sites_for_isotope(self._sites, iso) if iso else []
        self.table.setRowCount(len(sites))
        for r, s in enumerate(sites):
            q = s.quadrupolar()
            sh = s.shielding()
            cells = [s.label, s.isotope,
                     f"{q['Cq_MHz']:.3f} / {q['eta']:.2f}" if q else "—",
                     f"{sh['iso_ppm']:.1f}" if sh else "—"]
            for c, t in enumerate(cells):
                it = QTableWidgetItem(t)
                it.setFlags(Qt.ItemIsEnabled)
                self.table.setItem(r, c, it)

    def _accept(self):
        from larmor import dft

        iso = self.iso.currentText()
        ref = self.ref.value() or None
        sites = dft.sites_for_isotope(self._sites, iso) if iso else []
        self.result_sites = []
        for s in sites:
            sd = s.to_site_dict(model=self.model.currentText(),
                                reference_ppm=ref)
            self.result_sites.append({k: v for k, v in sd.items()
                                      if k != "notes"})
        self.accept()


class ErrorsDialog(QDialog):
    def __init__(self, parent, recipe: dict, ppm, amp, window):
        super().__init__(parent)
        self.setWindowTitle("Errors Analysis — χ² profile")
        self.resize(720, 520)
        self.recipe, self.ppm, self.amp, self.window = recipe, ppm, amp, window
        v = QVBoxLayout(self)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("site"))
        self.site = QComboBox()
        for i, s in enumerate(recipe["sites"]):
            self.site.addItem(f"s{i} — {s.get('label') or s['model']}", i)
        opts.addWidget(self.site)
        opts.addWidget(QLabel("parameter"))
        self.param = QComboBox()
        opts.addWidget(self.param)
        self.btnRun = QPushButton("Scan")
        self.btnRun.setDefault(True)
        self.btnRun.clicked.connect(self._run)
        opts.addWidget(self.btnRun)
        opts.addStretch(1)
        v.addLayout(opts)
        self.site.currentIndexChanged.connect(self._fill_params)
        self._fill_params()

        self.plot = pg.PlotWidget(background="w")
        self.plot.setLabel("bottom", "parameter value")
        self.plot.setLabel("left", "χ²")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        v.addWidget(self.plot, 1)
        self.res = QLabel("")
        self.res.setStyleSheet("font-weight: 700; color: #0a5a62;")
        v.addWidget(self.res)

    def _fill_params(self):
        i = self.site.currentData()
        self.param.clear()
        varying = [n for n, p in self.recipe["sites"][i]["params"].items()
                   if p.get("vary", True) and not p.get("expr")]
        self.param.addItems(varying)

    def _run(self):
        from larmor import autofit
        from larmor.recipe import Recipe

        i = self.site.currentData()
        pname = self.param.currentText()
        if not pname:
            return
        self.res.setText("scanning…")
        QApplication.processEvents()
        try:
            prof = autofit.error_profile(
                Recipe.from_dict(self.recipe), self.ppm, self.amp,
                site=i, param=pname, window_ppm=self.window)
        except Exception as exc:
            self.res.setText(f"failed: {exc}")
            return
        self.plot.clear()
        self.plot.plot(prof.values, prof.chi2, pen=pg.mkPen("#0e7c86", width=1.5),
                       symbol="o", symbolSize=5)
        # 1sigma and 2sigma levels
        for lvl, col in ((prof.chi2_min + 1.0, "#d62728"),
                         (prof.chi2_min + 3.84, "#c88a1e")):
            line = pg.InfiniteLine(pos=lvl, angle=0,
                                   pen=pg.mkPen(col, style=Qt.DashLine))
            self.plot.addItem(line)
        self.res.setText(prof.summary + ("   ·   " + " · ".join(prof.notes)
                                         if prof.notes else ""))
