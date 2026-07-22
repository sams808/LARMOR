"""Co-fit dialog: fit the current model against several datasets at once.

The sites of the current fit are the shared physical model; each dataset (1D
spectrum or 2D MQMAS map) gets its own copy with its own field, and the chosen
parameters (Cq, eta, delta, widths…) are tied across all of them. The classic
use is a multi-field fit or an MQMAS map fit jointly with a 1D MAS spectrum.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox,
    QPlainTextEdit, QPushButton, QScrollArea, QSplitter, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

from larmor.multifit import DEFAULT_SHARE


def _add_contours(pw, z, f2, f1, color, style, n_levels=7, name=""):
    """Draw log-spaced positive contours of z(f1, f2) into a PlotWidget, mapping
    array indices to ppm coordinates (same approach as the main 2D view)."""
    z = np.asarray(z, float)
    # decimate so marching-squares stays cheap on fine grids
    s2 = max(1, -(-z.shape[1] // 512)); s1 = max(1, -(-z.shape[0] // 512))
    zc, f2c, f1c = z[::s1, ::s2], np.asarray(f2)[::s2], np.asarray(f1)[::s1]
    top = float(np.nanmax(zc)) if zc.size else 0.0
    if top <= 0 or f2c.size < 2 or f1c.size < 2:
        return
    levels = np.logspace(np.log10(0.06 * top), np.log10(top), n_levels)
    tr = pg.QtGui.QTransform()
    tr.translate(f2c[0], f1c[0])
    tr.scale((f2c[-1] - f2c[0]) / max(zc.shape[1] - 1, 1),
             (f1c[-1] - f1c[0]) / max(zc.shape[0] - 1, 1))
    zt = np.ascontiguousarray(zc.T)
    for lvl in levels:
        iso = pg.IsocurveItem(data=zt, level=lvl,
                              pen=pg.mkPen(color, width=1, style=style))
        iso.setTransform(tr)
        pw.addItem(iso)


class CofitDialog(QDialog):
    applied = Signal(object)                 # fitted master recipe (dict)

    def __init__(self, parent, base_recipe: dict, base_dataset: dict):
        super().__init__(parent)
        self.setWindowTitle("Co-fit datasets (shared model)")
        self.resize(1140, 620)
        self.base_recipe = base_recipe
        self.datasets: list[dict] = [base_dataset]   # each: kind, label, ...
        self._result = None

        outer = QHBoxLayout(self)
        split = QSplitter(Qt.Horizontal)
        outer.addWidget(split)
        left = QWidget(); v = QVBoxLayout(left); v.setContentsMargins(0, 0, 0, 0)
        n = len(base_recipe.get("sites", []))
        v.addWidget(QLabel(
            f"Shared model: <b>{n} site(s)</b> from the current fit. Add the "
            "other spectra/maps of the same sample; tie the physical "
            "parameters below."))

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["kind", "dataset", "nucleus",
                                              "Larmor (MHz)"])
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table, 1)

        row = QHBoxLayout()
        self.btnAdd = QPushButton("＋ Add dataset…")
        self.btnAdd.clicked.connect(self._add)
        self.btnRemove = QPushButton("Remove selected")
        self.btnRemove.clicked.connect(self._remove)
        row.addWidget(self.btnAdd); row.addWidget(self.btnRemove); row.addStretch(1)
        v.addLayout(row)

        v.addWidget(QLabel("Tie across datasets:"))
        shrow = QHBoxLayout()
        self.share_boxes = {}
        present = {p for s in base_recipe.get("sites", [])
                   for p in s.get("params", {})}
        for name in DEFAULT_SHARE:
            if name in present:
                cb = QCheckBox(name.replace("_", " ").replace(" ppm", "")
                               .replace(" MHz", ""))
                cb.setChecked(True)
                self.share_boxes[name] = cb
                shrow.addWidget(cb)
        shrow.addStretch(1)
        v.addLayout(shrow)

        run = QHBoxLayout()
        self.btnRun = QPushButton("Run co-fit")
        self.btnRun.setDefault(True); self.btnRun.clicked.connect(self._run)
        self.btnApply = QPushButton("Apply shared params to current fit")
        self.btnApply.setEnabled(False); self.btnApply.clicked.connect(self._apply)
        run.addWidget(self.btnRun); run.addWidget(self.btnApply); run.addStretch(1)
        v.addLayout(run)

        self.report = QPlainTextEdit(); self.report.setReadOnly(True)
        self.report.setStyleSheet("font-family: Consolas, monospace; font-size: 11px;")
        v.addWidget(self.report, 1)

        # right: a live plot per dataset (experiment vs shared-model fit)
        self.plot_host = QWidget()
        self._plot_v = QVBoxLayout(self.plot_host)
        self._plot_v.setContentsMargins(4, 4, 4, 4)
        self._plot_hint = QLabel("Run the co-fit to see every dataset overlaid "
                                 "with the shared model.")
        self._plot_hint.setWordWrap(True); self._plot_hint.setAlignment(Qt.AlignCenter)
        self._plot_v.addWidget(self._plot_hint)
        self._plot_v.addStretch(1)
        scroll = QScrollArea(); scroll.setWidgetResizable(True)
        scroll.setWidget(self.plot_host)

        split.addWidget(left)
        split.addWidget(scroll)
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 1)
        split.setSizes([460, 680])

        self._refresh()

    # ------------------------------------------------------------------
    def _refresh(self):
        self.table.setRowCount(len(self.datasets))
        for i, d in enumerate(self.datasets):
            cells = [d["kind"], d["label"], d.get("nucleus", ""),
                     f"{d.get('larmor', 0.0):.3f}"]
            for j, c in enumerate(cells):
                self.table.setItem(i, j, QTableWidgetItem(str(c)))

    def _add(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Add a dataset to co-fit", "",
            "Spectra / maps (*.fxmla *.json 1r 2rr *.txt *.csv);;All files (*)")
        if not path:
            path = QFileDialog.getExistingDirectory(self, "…or a 2D EXPNO/pdata")
        if not path:
            return
        try:
            self.datasets.append(self._load(path))
        except Exception as exc:
            QMessageBox.warning(self, "Add dataset", f"cannot read: {exc}")
            return
        self._refresh()

    def _remove(self):
        r = self.table.currentRow()
        if 0 < r < len(self.datasets):          # never remove the base dataset
            del self.datasets[r]
            self._refresh()

    def _load(self, path: str) -> dict:
        from larmor import twod
        from larmor.io import bruker

        # try a 2D first
        try:
            ref = bruker.resolve(path)
            if ref.ndim == 2 and ref.target == "2rr":
                d = twod.read_bruker_2d(path)
                if d.z.ndim == 2:
                    return {"kind": "2d", "label": Path(path).name, "data2d": d,
                            "nucleus": d.nucleus, "larmor": d.larmor_MHz}
        except Exception:
            pass
        # else a 1D spectrum
        from larmor.desktop.app import _load_any

        ppm, amp, recipe, *_ = _load_any(path)
        return {"kind": "1d", "label": recipe.get("sample") or Path(path).name,
                "ppm": np.asarray(ppm), "amp": np.asarray(amp),
                "nucleus": recipe.get("nucleus", ""),
                "larmor": recipe.get("larmor_frequency_MHz", 0.0)}

    # ------------------------------------------------------------------
    def _clone_recipe(self, ds: dict):
        from larmor.recipe import Recipe

        r = Recipe.from_dict(json.loads(json.dumps(self.base_recipe)))
        if ds.get("nucleus"):
            r.nucleus = ds["nucleus"]
        if ds.get("larmor"):
            r.larmor_frequency_MHz = ds["larmor"]
        return r

    def _run(self):
        from larmor.multifit import fit_cofit

        if len(self.datasets) < 2:
            QMessageBox.information(self, "Co-fit",
                                    "add at least one more dataset")
            return
        if not self.base_recipe.get("sites"):
            QMessageBox.information(self, "Co-fit", "the current fit has no lines")
            return
        entries = []
        for ds in self.datasets:
            r = self._clone_recipe(ds)
            spec = ds["data2d"] if ds["kind"] == "2d" else (ds["ppm"], ds["amp"])
            entries.append((r, spec))
        share = tuple(n for n, cb in self.share_boxes.items() if cb.isChecked())
        self.report.setPlainText("running co-fit… (building any MQMAS kernels)")
        self.btnRun.setEnabled(False)
        from PySide6.QtWidgets import QApplication

        QApplication.processEvents()
        try:
            self._result = fit_cofit(entries, share=share)
        except Exception as exc:
            self.report.setPlainText(f"co-fit failed: {exc}")
            self.btnRun.setEnabled(True)
            return
        self.btnRun.setEnabled(True)
        self.btnApply.setEnabled(True)
        lines = [f"shared: {', '.join(share)}", ""]
        for ds, r, rmsd in zip(self.datasets, self._result.recipes,
                               self._result.rmsd):
            extra = ""
            if ds["kind"] == "2d":
                extra = (f"  ·  F1 ref "
                         f"{getattr(r, 'mqmas_f1_ref_ppm', 0.0):+.1f} ppm")
            lines.append(f"[{ds['kind']}] {ds['label']}: RMSD {rmsd:.4f}{extra}")
        lines.append("")
        for i, site in enumerate(self._result.recipes[0].sites):
            lines.append(f"site {site.label or i} ({site.model}):")
            for pn, p in site.params.items():
                err = f" ± {p.stderr:.4g}" if p.stderr else ""
                lines.append(f"    {pn} = {p.value:.5g}{err}")
        self.report.setPlainText("\n".join(lines))
        self._plot_result()

    # ------------------------------------------------------------------ plots
    def _clear_plots(self):
        while self._plot_v.count():
            it = self._plot_v.takeAt(0)
            w = it.widget()
            if w is not None:
                w.setParent(None)

    def _plot_result(self):
        """One plot per dataset: 1D spectra as experiment-vs-fit overlays, 2D
        maps as experiment contours with the shared-model contours on top."""
        if self._result is None:
            return
        self._clear_plots()
        per = self._result.per_dataset
        for ds, pd, rmsd in zip(self.datasets, per, self._result.rmsd):
            title = f"[{ds['kind']}] {ds['label']} — RMSD {rmsd:.4f}"
            if ds["kind"] == "1d":
                w = self._plot_1d(ds, pd, title)
            else:
                w = self._plot_2d(ds, pd, title)
            self._plot_v.addWidget(w, 1)

    def _plot_1d(self, ds, pd, title):
        pw = pg.PlotWidget(title=title)
        pw.setMinimumHeight(230)
        pw.getViewBox().invertX(True)          # ppm high -> low
        pw.setLabel("bottom", "ppm")
        pw.addLegend(offset=(-10, 10))
        pw.plot(np.asarray(ds["ppm"]), np.asarray(ds["amp"]),
                pen=pg.mkPen("#888", width=1), name="experiment")
        pw.plot(np.asarray(pd["x"]), np.asarray(pd["y_fit"]),
                pen=pg.mkPen("#d62728", width=2), name="fit")
        return pw

    def _plot_2d(self, ds, pd, title):
        pw = pg.PlotWidget(title=title)
        pw.setMinimumHeight(300)
        vb = pw.getViewBox()
        vb.invertX(True); vb.invertY(True); vb.setAspectLocked(False)
        pw.setLabel("bottom", "F2 (ppm)"); pw.setLabel("left", "F1 (ppm)")
        d = ds["data2d"].normalized()
        _add_contours(pw, d.z, d.f2_ppm, d.f1_ppm, "#3b7dd8",
                      Qt.SolidLine, name="experiment")
        _add_contours(pw, np.asarray(pd["z_fit"]), np.asarray(pd["f2"]),
                      np.asarray(pd["f1"]), "#e8832a", Qt.DashLine, name="fit")
        # legend proxies (IsocurveItems don't register in the legend)
        leg = pw.addLegend(offset=(-10, 10))
        leg.addItem(pg.PlotDataItem(pen=pg.mkPen("#3b7dd8", width=2)), "experiment")
        leg.addItem(pg.PlotDataItem(pen=pg.mkPen("#e8832a", width=2,
                                                 style=Qt.DashLine)), "model")
        return pw

    def _apply(self):
        if self._result is None:
            return
        self.applied.emit(self._result.recipes[0].to_dict())
        self.accept()
