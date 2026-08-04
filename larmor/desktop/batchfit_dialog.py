"""Batch fit dialog: one shared model fitted to many 1D spectra at once.

Ctrl/Shift-select spectra in the Explorer, hit "Batch fit selected…". This shows
them in a 3×3 grid (tabs for more), fits them with all parameters except
amplitude shared, and can then "release" chosen parameters to drift a little per
spectrum. Front-end for larmor.batchfit; the fit runs in an interruptible worker.
"""
from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFileDialog, QGridLayout, QHBoxLayout, QLabel, QProgressBar, QPushButton,
    QScrollArea, QTabWidget, QVBoxLayout, QWidget,
)

from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS

PER_TAB = 9        # 3×3 grid per tab


class _BatchWorker(QThread):
    done = Signal(object)
    failed = Signal(str)
    progress = Signal(int, float)

    def __init__(self, entries, release, frac):
        super().__init__()
        self.entries, self.release, self.frac = entries, release, frac
        self._stop = False

    def request_stop(self):
        self._stop = True

    def run(self):
        try:
            from larmor.batchfit import batch_fit

            state = {"n": 0}

            def cb(params, it, resid, *a, **k):
                state["n"] += 1
                try:
                    rms = float(np.sqrt(np.mean(np.asarray(resid, float) ** 2)))
                except Exception:
                    rms = float("nan")
                self.progress.emit(state["n"], rms)
                return True if self._stop else None

            res = batch_fit(self.entries, release=self.release,
                            release_frac=self.frac, iter_cb=cb)
            self.done.emit(res)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))


class BatchFitDialog(QDialog):
    def __init__(self, parent, paths, model_recipe: dict | None):
        super().__init__(parent)
        self.setWindowTitle("Batch fit — one shared model, amplitudes per spectrum")
        self.resize(940, 720)
        self._model_sites = ((model_recipe or {}).get("sites") or None)
        self._window = ((model_recipe or {}).get("fit_window_ppm") or None)
        self._result = None
        self._worker = None
        self._cells = []
        self._data = self._load(paths)

        v = QVBoxLayout(self)
        top = QHBoxLayout()
        self.lblModel = QLabel()
        b_model = QPushButton("Model from recipe…")
        b_model.setToolTip("load a saved recipe to use as the shared model")
        b_model.clicked.connect(self._pick_model)
        top.addWidget(self.lblModel, 1)
        top.addWidget(b_model)
        v.addLayout(top)

        self.tabs = QTabWidget()
        v.addWidget(self.tabs, 1)
        self._build_grid()

        # release panel
        rel = QHBoxLayout()
        rel.addWidget(QLabel("Release (let drift per spectrum):"))
        self._relbox = QScrollArea()
        self._relbox.setWidgetResizable(True)
        self._relbox.setMaximumHeight(56)
        holder = QWidget(); self._rellay = QHBoxLayout(holder)
        self._rellay.setContentsMargins(2, 2, 2, 2)
        self._relbox.setWidget(holder)
        rel.addWidget(self._relbox, 1)
        rel.addWidget(QLabel("±"))
        self.frac = QDoubleSpinBox()
        self.frac.setRange(1, 100); self.frac.setValue(10); self.frac.setSuffix(" %")
        rel.addWidget(self.frac)
        v.addLayout(rel)
        self._rel_checks = {}
        self._fill_release_params()

        self.prog = QProgressBar(); v.addWidget(self.prog)
        self.status = QLabel(self._model_status())
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"font-weight:600; color:{theme.active().accent};")
        v.addWidget(self.status)

        bb = QDialogButtonBox(QDialogButtonBox.Close | QDialogButtonBox.Help)
        self.btnFit = bb.addButton("Fit (shared)", QDialogButtonBox.ApplyRole)
        self.btnFit.clicked.connect(lambda: self._run(release=False))
        self.btnRelease = bb.addButton("Release fit", QDialogButtonBox.ApplyRole)
        self.btnRelease.clicked.connect(lambda: self._run(release=True))
        self.btnRelease.setEnabled(False)
        self.btnStop = bb.addButton("Stop", QDialogButtonBox.ResetRole)
        self.btnStop.setEnabled(False)
        self.btnStop.clicked.connect(self._stop)
        self.btnSave = bb.addButton("Save fits & table…", QDialogButtonBox.ActionRole)
        self.btnSave.setEnabled(False)
        self.btnSave.clicked.connect(self._save)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.helpRequested.connect(self._help)
        v.addWidget(bb)
        self._update_fit_enabled()

    # ------------------------------------------------------------------
    def _load(self, paths):
        from larmor.loader import load_any

        data = []
        for p in paths:
            try:
                ppm, amp, rec, meta, warns = load_any(p)
            except Exception:
                continue
            data.append({
                "ppm": np.asarray(ppm, float), "amp": np.asarray(amp, float),
                "nucleus": rec.get("nucleus", ""),
                "larmor": float(rec.get("larmor_frequency_MHz", 0.0) or 0.0),
                "spin": float(rec.get("spin_rate_Hz", 0.0) or 0.0),
                "sample": rec.get("sample") or Path(p).stem, "path": p})
            # if no model yet, adopt a loaded recipe that has sites
            if self._model_sites is None and rec.get("sites"):
                self._model_sites = rec["sites"]
        return data

    def _build_grid(self):
        n = len(self._data)
        for start in range(0, n, PER_TAB):
            page = QWidget(); grid = QGridLayout(page)
            for j in range(PER_TAB):
                k = start + j
                if k >= n:
                    break
                d = self._data[k]
                cell = QWidget(); cv = QVBoxLayout(cell)
                cv.setContentsMargins(2, 2, 2, 2); cv.setSpacing(1)
                title = QLabel(d["sample"])
                title.setStyleSheet("font-size:10px; font-weight:600;")
                plot = pg.PlotWidget(background=theme.active().plot_bg)
                plot.setMenuEnabled(False); plot.hideAxis("left")
                plot.getPlotItem().getViewBox().setMouseEnabled(False, False)
                plot.setMinimumHeight(120)
                exp = plot.plot(d["ppm"], d["amp"],
                                pen=pg.mkPen(theme.active().experiment, width=1))
                model = plot.plot([], [],
                                  pen=pg.mkPen(theme.active().model, width=1.2))
                plot.setXRange(d["ppm"].max(), d["ppm"].min())
                rmsd = QLabel(""); rmsd.setStyleSheet("font-size:9px; color:#888;")
                cv.addWidget(title); cv.addWidget(plot, 1); cv.addWidget(rmsd)
                grid.addWidget(cell, j // 3, j % 3)
                self._cells.append({"exp": exp, "model": model, "rmsd": rmsd})
            self.tabs.addTab(page, f"{start + 1}–{min(start + PER_TAB, n)}")

    def _fill_release_params(self):
        while self._rellay.count():
            it = self._rellay.takeAt(0)
            if it.widget():
                it.widget().deleteLater()
        self._rel_checks = {}
        if not self._model_sites:
            return
        from larmor.batchfit import all_but_amplitude
        from larmor.recipe import Recipe

        names = all_but_amplitude([Recipe.from_dict(
            {"sites": self._model_sites, "nucleus": "", "larmor_frequency_MHz": 0})])
        for pn in names:
            if pn == "gl":
                continue
            c = QCheckBox(PARAM_LABELS.get(pn, pn))
            self._rel_checks[pn] = c
            self._rellay.addWidget(c)
        self._rellay.addStretch(1)

    # ------------------------------------------------------------------
    def _model_status(self) -> str:
        if not self._model_sites:
            return ("no shared model yet — load a recipe (or fit one spectrum "
                    "first, then reopen from a fit)")
        nuclei = {d["nucleus"] for d in self._data if d["nucleus"]}
        warn = ("  ⚠ mixed nuclei" if len(nuclei) > 1 else "")
        return (f"{len(self._data)} spectra · model: "
                f"{len(self._model_sites)} line(s) · sharing all but amplitude"
                + warn)

    def _model_label(self) -> str:
        return ("model: none — load a recipe" if not self._model_sites
                else f"model: {len(self._model_sites)} line(s), sharing all "
                     "parameters except amplitude")

    def _update_fit_enabled(self):
        ok = bool(self._model_sites) and len(self._data) >= 2
        self.btnFit.setEnabled(ok)
        self.lblModel.setText(self._model_label())
        self.status.setText(self._model_status())

    def _pick_model(self):
        from larmor.recipe import Recipe

        path, _ = QFileDialog.getOpenFileName(
            self, "Model recipe", "", "LARMOR recipe (*.json);;All (*)")
        if not path:
            return
        try:
            d = Recipe.load(path).to_dict()
        except Exception as exc:  # noqa: BLE001
            self.status.setText(f"could not load model: {exc}")
            return
        if not d.get("sites"):
            self.status.setText("that recipe has no lines")
            return
        self._model_sites = d["sites"]
        self._window = d.get("fit_window_ppm") or self._window
        self._fill_release_params()
        self._update_fit_enabled()

    def _entries(self):
        from larmor.recipe import Recipe

        out = []
        for d in self._data:
            rec = Recipe.from_dict({
                "nucleus": d["nucleus"], "larmor_frequency_MHz": d["larmor"],
                "spin_rate_Hz": d["spin"], "sample": d["sample"],
                "sites": copy.deepcopy(self._model_sites)})
            out.append((rec, d["ppm"], d["amp"], self._window))
        return out

    def _run(self, release: bool):
        if not self._model_sites or len(self._data) < 2:
            self.status.setText("need a model and at least two spectra")
            return
        rel = tuple(pn for pn, c in self._rel_checks.items() if c.isChecked()) \
            if release else ()
        if release and not rel:
            self.status.setText("tick at least one parameter to release")
            return
        self.btnFit.setEnabled(False); self.btnRelease.setEnabled(False)
        self.btnStop.setEnabled(True); self.btnSave.setEnabled(False)
        self.prog.setRange(0, 0)                       # busy while fitting
        self.status.setText("fitting…")
        self._worker = _BatchWorker(self._entries(), rel, self.frac.value() / 100.0)
        self._worker.progress.connect(
            lambda it, rms: self.status.setText(f"fitting — iter {it} · rms {rms:.3g}"))
        self._worker.done.connect(self._done)
        self._worker.failed.connect(self._failed)
        self._worker.start()

    def _stop(self):
        if self._worker is not None and self._worker.isRunning():
            self._worker.request_stop()
            self.status.setText("stopping — keeping the latest values…")

    def _done(self, result):
        self._result = result
        self.prog.setRange(0, 100); self.prog.setValue(100)
        self.btnStop.setEnabled(False)
        self._update_fit_enabled()
        self.btnRelease.setEnabled(True); self.btnSave.setEnabled(True)
        for k, pd in enumerate(result.per_dataset):
            if k < len(self._cells):
                self._cells[k]["model"].setData(np.asarray(pd["x"], float),
                                                np.asarray(pd["y_fit"], float))
                self._cells[k]["rmsd"].setText(f"RMSD {result.rmsd[k]:.4f}")
        self.status.setText(result.summary)

    def _failed(self, msg):
        self.prog.setRange(0, 100); self.prog.setValue(0)
        self.btnStop.setEnabled(False)
        self._update_fit_enabled()
        self.status.setText(f"batch fit failed: {msg}")

    def _save(self):
        if self._result is None:
            return
        from larmor.recipe import Recipe
        from larmor import batchfit

        folder = QFileDialog.getExistingDirectory(self, "Save fits + table")
        if not folder:
            return
        folder = Path(folder)
        n = 0
        for rec in self._result.recipes:
            slug = "".join(c if c.isalnum() or c in "-_" else "_"
                           for c in (rec.sample or "fit"))[:40] or "fit"
            try:
                Recipe.from_dict(rec.to_dict()).save(folder / f"{slug}.recipe.json")
                n += 1
            except Exception:
                pass
        # a flat CSV of the shared / per-spectrum values
        rows = batchfit.shared_table(self._result)
        import csv
        with open(folder / "batch_table.csv", "w", newline="",
                  encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["scope", "site", "label", "param", "value", "stderr"])
            for r in rows:
                w.writerow([r["scope"], r["site"], r["label"], r["param"],
                            f"{r['value']:.6g}",
                            "" if r["stderr"] is None else f"{r['stderr']:.4g}"])
        self.status.setText(f"saved {n} recipe(s) + batch_table.csv to {folder}")

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "multi-dataset", "Multi-dataset & co-fitting")
