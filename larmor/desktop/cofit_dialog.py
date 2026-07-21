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
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox,
    QPlainTextEdit, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout,
    QWidget,
)

from larmor.multifit import DEFAULT_SHARE


class CofitDialog(QDialog):
    applied = Signal(object)                 # fitted master recipe (dict)

    def __init__(self, parent, base_recipe: dict, base_dataset: dict):
        super().__init__(parent)
        self.setWindowTitle("Co-fit datasets (shared model)")
        self.resize(760, 560)
        self.base_recipe = base_recipe
        self.datasets: list[dict] = [base_dataset]   # each: kind, label, ...
        self._result = None

        v = QVBoxLayout(self)
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
            lines.append(f"[{ds['kind']}] {ds['label']}: RMSD {rmsd:.4f}")
        lines.append("")
        for i, site in enumerate(self._result.recipes[0].sites):
            lines.append(f"site {site.label or i} ({site.model}):")
            for pn, p in site.params.items():
                err = f" ± {p.stderr:.4g}" if p.stderr else ""
                lines.append(f"    {pn} = {p.value:.5g}{err}")
        self.report.setPlainText("\n".join(lines))

    def _apply(self):
        if self._result is None:
            return
        self.applied.emit(self._result.recipes[0].to_dict())
        self.accept()
