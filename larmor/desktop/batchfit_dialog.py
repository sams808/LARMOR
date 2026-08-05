"""Batch fit dialog: one shared model fitted to many 1D spectra at once.

Ctrl/Shift-select spectra in the Explorer, hit "Batch fit selected…". This shows
them in a 3×3 grid (tabs for more), fits them with all parameters except
amplitude shared, and can then **release** chosen parameters — per-parameter — to
drift a little per spectrum. Front-end for larmor.batchfit; the fit runs in an
interruptible worker (Cancel = revert, Stop = keep the last iteration).

Made to be a proper workbench: NMR-style axes (high→low ppm), mouse zoom,
independent-vs-shared scale, live component curves, an optional per-fit baseline,
a completion threshold, and one-click saving of every fit in LARMOR format.
"""
from __future__ import annotations

import copy
import datetime as _dt
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
    QGridLayout, QHBoxLayout, QInputDialog, QLabel, QMessageBox, QProgressBar,
    QPushButton, QScrollArea, QSpinBox, QTabWidget, QVBoxLayout, QWidget,
)

from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS
from larmor.desktop.plot import site_color

PER_TAB = 9        # 3×3 grid per tab


def _poly_baseline(x, y, order: int) -> np.ndarray:
    """Asymmetric polynomial baseline: fit, then iteratively keep points at/below
    the fit (peaks are positive) — a robust rolling baseline for 1D NMR."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.ones(y.size, bool)
    base = np.zeros_like(y)
    for _ in range(12):
        if m.sum() < order + 2:
            break
        c = np.polyfit(x[m], y[m], order)
        base = np.polyval(c, x)
        resid = y - base
        sd = float(np.std(resid[m])) or 1.0
        m = resid < sd                     # drop points that stick up (signal)
    return base


def estimate_baseline(x, y, kind: str, order: int = 3) -> np.ndarray:
    """Per-spectrum baseline of the chosen kind (returns the baseline to subtract)."""
    if kind == "Polynomial":
        return _poly_baseline(x, y, order)
    if kind == "Iterative (Yon 2020)":
        from larmor.baseline import iterative_baseline
        return iterative_baseline(np.asarray(y, float),
                                  x=np.asarray(x, float)).baseline
    if kind == "Flat (edge median)":
        n = max(3, y.size // 20)
        lvl = float(np.median(np.concatenate([y[:n], y[-n:]])))
        return np.full_like(np.asarray(y, float), lvl)
    return np.zeros_like(np.asarray(y, float))


class _BatchWorker(QThread):
    done = Signal(object, str)          # (result, stop_mode)
    failed = Signal(str)
    progress = Signal(int, float)

    def __init__(self, entries, release, frac, tol=None):
        super().__init__()
        self.entries, self.release, self.frac, self.tol = \
            entries, release, frac, tol
        self._stop = False
        self._mode = ""

    def request_stop(self, mode: str):
        self._mode = mode               # "cancel" (revert) | "stop" (keep last)
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
                            release_frac=self.frac, iter_cb=cb, tol=self.tol)
            self.done.emit(res, self._mode)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))


class BatchFitDialog(QDialog):
    def __init__(self, parent, paths, model_recipe: dict | None):
        super().__init__(parent)
        self.setWindowTitle("Batch fit — one shared model, amplitudes per spectrum")
        self.resize(1000, 780)
        self._model_sites = ((model_recipe or {}).get("sites") or None)
        self._window = ((model_recipe or {}).get("fit_window_ppm") or None)
        self._recipe_tag = ""            # set when a model is loaded from a recipe
        self._result = None
        self._worker = None
        self._cells: list[dict] = []
        self._show_comp = False
        self._shared_scale = False
        self._baseline_kind = "None"
        self._data = self._load(paths)

        v = QVBoxLayout(self)

        # ---- model row ----
        top = QHBoxLayout()
        self.lblModel = QLabel()
        b_model = QPushButton("Model from recipe…")
        b_model.setToolTip("load a saved recipe to use as the shared model")
        b_model.clicked.connect(self._pick_model)
        top.addWidget(self.lblModel, 1)
        top.addWidget(b_model)
        v.addLayout(top)

        # ---- the grid of spectra ----
        self.tabs = QTabWidget()
        v.addWidget(self.tabs, 1)
        self._build_grid()

        # ---- view options ----
        opt = QHBoxLayout()
        self.chkComp = QCheckBox("components")
        self.chkComp.setToolTip("overlay each site's component curve on every fit")
        self.chkComp.toggled.connect(self._toggle_components)
        self.chkShared = QCheckBox("shared scale")
        self.chkShared.setToolTip("put every plot on one common x/y scale "
                                  "(off = each spectrum auto-scales; drag to zoom, "
                                  "right-click-drag or scroll to rescale, "
                                  "right-click ▸ View All to reset)")
        self.chkShared.toggled.connect(self._toggle_shared)
        opt.addWidget(self.chkComp)
        opt.addWidget(self.chkShared)
        opt.addSpacing(16)
        opt.addWidget(QLabel("baseline:"))
        self.lblBaseline = QLabel("none")
        self.lblBaseline.setStyleSheet(f"color:{theme.active().text_dim};")
        b_base = QPushButton("Fit baseline…")
        b_base.setToolTip("estimate and subtract a baseline from every spectrum "
                          "independently before fitting")
        b_base.clicked.connect(self._fit_baseline)
        b_baser = QPushButton("Reset")
        b_baser.setToolTip("restore the raw spectra (remove the baseline)")
        b_baser.clicked.connect(self._reset_baseline)
        opt.addWidget(b_base)
        opt.addWidget(b_baser)
        opt.addWidget(self.lblBaseline, 1)
        v.addLayout(opt)

        # ---- release panel (per-parameter) ----
        rel = QHBoxLayout()
        rel.addWidget(QLabel("Release per spectrum:"))
        self._relbox = QScrollArea()
        self._relbox.setWidgetResizable(True)
        self._relbox.setMaximumHeight(56)
        holder = QWidget(); self._rellay = QHBoxLayout(holder)
        self._rellay.setContentsMargins(2, 2, 2, 2)
        self._relbox.setWidget(holder)
        rel.addWidget(self._relbox, 1)
        rel.addWidget(QLabel("±"))
        self.frac = QDoubleSpinBox()
        self.frac.setRange(0.1, 100); self.frac.setValue(10); self.frac.setSuffix(" %")
        self.frac.setToolTip("how far a released parameter may drift around its "
                             "shared value")
        rel.addWidget(self.frac)
        v.addLayout(rel)
        self._rel_checks: dict = {}
        self._fill_release_params()

        # ---- fit controls ----
        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("Completion threshold:"))
        self.tol = QDoubleSpinBox()
        self.tol.setRange(0.0, 50.0); self.tol.setDecimals(3)
        self.tol.setValue(_saved_tol()); self.tol.setSuffix(" % Δσ")
        self.tol.setToolTip("stop once the residual stdev changes by less than "
                            "this between iterations (0 = fit to full precision)")
        ctl.addWidget(self.tol)
        ctl.addStretch(1)
        v.addLayout(ctl)

        self.prog = QProgressBar(); v.addWidget(self.prog)
        self.status = QLabel(self._model_status())
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"font-weight:600; color:{theme.active().accent};")
        v.addWidget(self.status)

        # ---- buttons ----
        bb = QDialogButtonBox(QDialogButtonBox.Close | QDialogButtonBox.Help)
        self.btnFit = bb.addButton("Fit", QDialogButtonBox.ApplyRole)
        self.btnFit.setToolTip("share all parameters but amplitude; also release "
                               "the ticked parameters (±) per spectrum")
        self.btnFit.clicked.connect(self._run)
        self.btnCancel = bb.addButton("Cancel", QDialogButtonBox.ResetRole)
        self.btnCancel.setToolTip("stop and DISCARD this fit (revert to before it)")
        self.btnCancel.setEnabled(False)
        self.btnCancel.clicked.connect(lambda: self._interrupt("cancel"))
        self.btnStop = bb.addButton("Stop", QDialogButtonBox.ResetRole)
        self.btnStop.setToolTip("stop now and KEEP the latest iteration's values")
        self.btnStop.setEnabled(False)
        self.btnStop.clicked.connect(lambda: self._interrupt("stop"))
        self.btnSave = bb.addButton("Save individual fits…", QDialogButtonBox.ActionRole)
        self.btnSave.setToolTip("write one LARMOR .recipe.json per spectrum")
        self.btnSave.setEnabled(False)
        self.btnSave.clicked.connect(self._save_individual)
        self.btnTable = bb.addButton("Save table…", QDialogButtonBox.ActionRole)
        self.btnTable.setToolTip("write a batch_table.csv of shared / per-spectrum values")
        self.btnTable.setEnabled(False)
        self.btnTable.clicked.connect(self._save_table)
        self.btnSeries = bb.addButton("Series plot…", QDialogButtonBox.ActionRole)
        self.btnSeries.setToolTip("plot how each parameter evolves along the series")
        self.btnSeries.setEnabled(False)
        self.btnSeries.clicked.connect(self._series_plot)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.helpRequested.connect(self._help)
        v.addWidget(bb)
        self._update_fit_enabled()

    # ------------------------------------------------------------------ load
    def _load(self, paths):
        from larmor.loader import load_any

        data = []
        for p in paths:
            try:
                ppm, amp, rec, meta, warns = load_any(p)
            except Exception:
                continue
            ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
            data.append({
                "ppm": ppm, "amp": amp, "amp0": amp.copy(),
                "nucleus": rec.get("nucleus", ""),
                "larmor": float(rec.get("larmor_frequency_MHz", 0.0) or 0.0),
                "spin": float(rec.get("spin_rate_Hz", 0.0) or 0.0),
                "sample": rec.get("sample") or Path(p).stem, "path": p,
                "proc": _proc_number(p)})
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
                title = QLabel(d["sample"])          # sample name, top-left
                title.setStyleSheet("font-size:10px; font-weight:600;")
                title.setToolTip(f"{d['sample']}"
                                 + (f" · proc {d['proc']}" if d["proc"] else "")
                                 + f" · {d['nucleus']}")
                plot = pg.PlotWidget(background=theme.active().plot_bg)
                plot.setMenuEnabled(True)            # right-click ▸ View All etc.
                plot.hideAxis("left")
                plot.getPlotItem().getViewBox().setMouseEnabled(True, True)  # zoom
                plot.setMinimumHeight(120)
                exp = plot.plot(d["ppm"], d["amp"],
                                pen=pg.mkPen(theme.active().experiment, width=1))
                model = plot.plot([], [],
                                  pen=pg.mkPen(theme.active().model, width=1.4))
                plot.setXRange(d["ppm"].max(), d["ppm"].min())   # NMR: high→low
                rmsd = QLabel(""); rmsd.setStyleSheet(
                    f"font-size:9px; color:{theme.active().text_dim};")
                cv.addWidget(title); cv.addWidget(plot, 1); cv.addWidget(rmsd)
                grid.addWidget(cell, j // 3, j % 3)
                self._cells.append({"plot": plot, "exp": exp, "model": model,
                                    "rmsd": rmsd, "comp": [], "title": title})
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
        for pn in names:                             # ALL lineshape params (incl gl)
            c = QCheckBox(PARAM_LABELS.get(pn, pn))
            c.setToolTip(f"let {pn} drift ±{self.frac.value():.0f}% per spectrum")
            self._rel_checks[pn] = c
            self._rellay.addWidget(c)
        self._rellay.addStretch(1)

    # ------------------------------------------------------------------ status
    def _model_status(self) -> str:
        if not self._model_sites:
            return ("no shared model yet — load a recipe (or fit one spectrum "
                    "first, then reopen from a fit)")
        nuclei = {d["nucleus"] for d in self._data if d["nucleus"]}
        warn = ("  ⚠ mixed nuclei" if len(nuclei) > 1 else "")
        return (f"{len(self._data)} spectra · model: {len(self._model_sites)} "
                f"line(s) · sharing all but amplitude" + warn)

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
        self._recipe_tag = Path(path).stem.replace(".recipe", "")
        self._fill_release_params()
        self._update_fit_enabled()

    # ------------------------------------------------------------------ baseline
    def _fit_baseline(self):
        kinds = ["Polynomial", "Iterative (Yon 2020)", "Flat (edge median)"]
        dlg = QDialog(self); dlg.setWindowTitle("Fit baseline")
        lay = QVBoxLayout(dlg)
        row = QHBoxLayout(); row.addWidget(QLabel("Type:"))
        combo = QComboBox(); combo.addItems(kinds); row.addWidget(combo, 1)
        lay.addLayout(row)
        row2 = QHBoxLayout(); row2.addWidget(QLabel("Polynomial order:"))
        order = QSpinBox(); order.setRange(0, 8); order.setValue(3)
        row2.addWidget(order); row2.addStretch(1)
        lay.addLayout(row2)
        combo.currentTextChanged.connect(
            lambda t: order.setEnabled(t == "Polynomial"))
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject)
        lay.addWidget(bb)
        if dlg.exec() != QDialog.Accepted:
            return
        kind, o = combo.currentText(), order.value()
        self._baseline_kind = kind
        for k, d in enumerate(self._data):
            base = estimate_baseline(d["ppm"], d["amp0"], kind, o)
            d["amp"] = d["amp0"] - base
            if k < len(self._cells):
                self._cells[k]["exp"].setData(d["ppm"], d["amp"])
        self.lblBaseline.setText(f"{kind}" + (f" (order {o})"
                                              if kind == "Polynomial" else ""))
        self.status.setText(f"baseline subtracted per spectrum: {kind}")

    def _reset_baseline(self):
        self._baseline_kind = "None"
        for k, d in enumerate(self._data):
            d["amp"] = d["amp0"].copy()
            if k < len(self._cells):
                self._cells[k]["exp"].setData(d["ppm"], d["amp"])
        self.lblBaseline.setText("none")
        self.status.setText("baseline removed — using raw spectra")

    # ------------------------------------------------------------------ view opts
    def _toggle_components(self, on: bool):
        self._show_comp = on
        self._refresh_components()

    def _refresh_components(self):
        for cell in self._cells:
            for it in cell["comp"]:
                cell["plot"].removeItem(it)
            cell["comp"] = []
        if not self._show_comp or self._result is None:
            return
        from larmor import engine
        for k, rec in enumerate(self._result.recipes):
            if k >= len(self._cells):
                break
            cell = self._cells[k]
            try:
                x, _tot, per = engine.simulate(rec, exp_ppm=self._data[k]["ppm"])
            except Exception:
                continue
            for i, ys in enumerate(per):
                it = cell["plot"].plot(
                    x, np.asarray(ys, float),
                    pen=pg.mkPen(site_color(i), width=1, style=Qt.DashLine))
                cell["comp"].append(it)

    def _toggle_shared(self, on: bool):
        self._shared_scale = on
        self._apply_scale()

    def _apply_scale(self):
        if not self._cells:
            return
        if self._shared_scale:
            xmin = min(float(d["ppm"].min()) for d in self._data)
            xmax = max(float(d["ppm"].max()) for d in self._data)
            ymin = min(float(np.min(d["amp"])) for d in self._data)
            ymax = max(float(np.max(d["amp"])) for d in self._data)
            for cell in self._cells:
                vb = cell["plot"].getViewBox()
                vb.setXRange(xmax, xmin, padding=0.02)
                vb.setYRange(ymin, ymax, padding=0.05)
        else:
            for k, cell in enumerate(self._cells):
                d = self._data[k]
                vb = cell["plot"].getViewBox()
                vb.enableAutoRange(axis="y")
                vb.setXRange(float(d["ppm"].max()), float(d["ppm"].min()),
                             padding=0.02)

    # ------------------------------------------------------------------ fit
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

    def _run(self):
        if not self._model_sites or len(self._data) < 2:
            self.status.setText("need a model and at least two spectra")
            return
        rel = tuple(pn for pn, c in self._rel_checks.items() if c.isChecked())
        self._pre = [(c["model"].xData, c["model"].yData) for c in self._cells]
        self.btnFit.setEnabled(False)
        self.btnCancel.setEnabled(True); self.btnStop.setEnabled(True)
        self.btnSave.setEnabled(False); self.btnTable.setEnabled(False)
        self.btnSeries.setEnabled(False)
        self.prog.setRange(0, 0)                       # busy while fitting
        rtxt = (f" · releasing {', '.join(rel)} (±{self.frac.value():.0f}%)"
                if rel else "")
        self.status.setText("fitting…" + rtxt)
        tol = self.tol.value() or None
        _save_tol(self.tol.value())
        self._worker = _BatchWorker(self._entries(), rel,
                                    self.frac.value() / 100.0, tol)
        self._worker.progress.connect(
            lambda it, rms: self.status.setText(
                f"fitting — iter {it} · rms {rms:.3g}" + rtxt))
        self._worker.done.connect(self._done)
        self._worker.failed.connect(self._failed)
        self._worker.start()

    def _interrupt(self, mode: str):
        if self._worker is not None and self._worker.isRunning():
            self._worker.request_stop(mode)
            self.status.setText("cancelling — discarding this fit…" if mode == "cancel"
                                else "stopping — keeping the latest values…")

    def _done(self, result, mode: str = ""):
        self.prog.setRange(0, 100); self.prog.setValue(100 if mode != "cancel" else 0)
        self.btnCancel.setEnabled(False); self.btnStop.setEnabled(False)
        self._update_fit_enabled()
        if mode == "cancel":                           # revert to the pre-fit view
            for k, cell in enumerate(self._cells):
                px, py = self._pre[k] if k < len(getattr(self, "_pre", [])) else (None, None)
                cell["model"].setData(px if px is not None else [],
                                      py if py is not None else [])
            self.status.setText("fit cancelled — reverted")
            return
        self._result = result
        self.btnSave.setEnabled(True); self.btnTable.setEnabled(True)
        self.btnSeries.setEnabled(True)
        for k, pd in enumerate(result.per_dataset):
            if k < len(self._cells):
                self._cells[k]["model"].setData(np.asarray(pd["x"], float),
                                                np.asarray(pd["y_fit"], float))
                self._cells[k]["rmsd"].setText(f"RMSD {result.rmsd[k]:.4f}")
        self._refresh_components()
        note = " (stopped early)" if mode == "stop" else ""
        self.status.setText(result.summary + note)

    def _failed(self, msg):
        self.prog.setRange(0, 100); self.prog.setValue(0)
        self.btnCancel.setEnabled(False); self.btnStop.setEnabled(False)
        self._update_fit_enabled()
        self.status.setText(f"batch fit failed: {msg}")

    # ------------------------------------------------------------------ save
    def _auto_name(self, rec, proc) -> str:
        parts = [rec.sample or "fit", rec.nucleus or ""]
        if self._recipe_tag:
            parts.append(self._recipe_tag)
        parts.append("batch")
        parts.append(_dt.datetime.now().strftime("%Y%m%d_%H%M"))
        slug = "_".join(p for p in parts if p)
        return _slug(slug)

    def _save_individual(self):
        if self._result is None:
            return
        from larmor.recipe import Recipe

        folder = QFileDialog.getExistingDirectory(self, "Save individual fits")
        if not folder:
            return
        folder = Path(folder)
        mode = QMessageBox.question(
            self, "Naming",
            "Name the files automatically?\n\n"
            "Yes — auto (sample_nucleus"
            + ("_recipe" if self._recipe_tag else "") + "_batch_YYYYMMDD_HHMM)\n"
            "No — type a name for each fit.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        recs = self._result.recipes
        n = 0
        for k, rec in enumerate(recs):
            proc = self._data[k]["proc"] if k < len(self._data) else ""
            if mode == QMessageBox.Yes:
                name = self._auto_name(rec, proc)
            else:
                default = self._auto_name(rec, proc)
                label = (f"Name for “{rec.sample}”"
                         + (f" (proc {proc})" if proc else "")
                         + f"   [{k + 1} of {len(recs)}]")
                text, ok = QInputDialog.getText(self, "Fit name", label,
                                                text=default)
                if not ok:
                    break
                name = _slug(text) or default
            try:
                Recipe.from_dict(rec.to_dict()).save(
                    folder / f"{name}.recipe.json")
                n += 1
            except Exception:
                pass
        self.status.setText(f"saved {n} fit(s) to {folder}")

    def _save_table(self):
        if self._result is None:
            return
        import csv
        from larmor import batchfit

        path, _ = QFileDialog.getSaveFileName(
            self, "Save batch table", "batch_table.csv", "CSV (*.csv)")
        if not path:
            return
        rows = batchfit.shared_table(self._result)
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["scope", "site", "label", "param", "value", "stderr"])
            for r in rows:
                w.writerow([r["scope"], r["site"], r["label"], r["param"],
                            f"{r['value']:.6g}",
                            "" if r["stderr"] is None else f"{r['stderr']:.4g}"])
        self.status.setText(f"saved {path}")

    def _series_plot(self):
        if self._result is None:
            return
        from larmor.desktop.series_plot import SeriesPlotDialog
        SeriesPlotDialog(self, self._result).exec()

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "multi-dataset", "Multi-dataset & co-fitting")


# ---------------------------------------------------------------- module helpers
def _slug(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in (s or ""))[:80]


def _proc_number(path: str) -> str:
    parts = Path(path).parts
    if "pdata" in parts:
        i = parts.index("pdata")
        if i + 1 < len(parts):
            return parts[i + 1]
    return ""


def _saved_tol() -> float:
    from PySide6.QtCore import QSettings
    try:
        return float(QSettings("LARMOR", "app").value("fitStdevPct", 0.0) or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _save_tol(v: float):
    from PySide6.QtCore import QSettings
    QSettings("LARMOR", "app").setValue("fitStdevPct", float(v))
