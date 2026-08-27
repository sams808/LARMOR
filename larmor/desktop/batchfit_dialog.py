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
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
    QGridLayout, QHBoxLayout, QInputDialog, QLabel, QMessageBox, QProgressBar,
    QPushButton, QScrollArea, QSpinBox, QTabWidget, QVBoxLayout, QWidget,
)

from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS
from larmor.desktop.plot import site_color

PER_TAB = 9        # 3×3 grid per tab


#: baseline kind (the "Fit baseline…" combo) -> the larmor.processing op that
#: reproduces it, so recording {"op": ..., **kwargs} into a recipe's
#: `processing` and replaying it later (loader.apply_processing, used by the
#: Plotting studio's batch-grid experiment/residual traces) gives back the
#: EXACT same corrected spectrum shown here -- one implementation, not two
#: that could quietly drift apart.
BASELINE_OPS = {
    "Polynomial": ("baseline", ("order",)),
    "Iterative (Yon 2020)": ("iterbaseline", ()),
    "Flat (edge median)": ("flat_baseline", ()),
}


def estimate_baseline(x, y, kind: str, order: int = 3) -> np.ndarray:
    """Per-spectrum baseline of the chosen kind (returns the baseline to
    subtract), computed via the same larmor.processing op that gets recorded
    into the recipe -- see BASELINE_OPS."""
    if kind not in BASELINE_OPS:
        return np.zeros_like(np.asarray(y, float))
    from larmor import processing as proc

    op_name, kw_names = BASELINE_OPS[kind]
    kwargs = {"order": order} if "order" in kw_names else {}
    s = proc.Spectrum1D(x_ppm=np.asarray(x, float),
                        y=np.asarray(y, float).astype(complex),
                        sfo1_MHz=1.0, sw_Hz=0.0, domain="freq")
    corrected = proc.OPS[op_name](s, **kwargs)
    return np.asarray(y, float) - corrected.y.real


def baseline_processing_op(kind: str, order: int = 3) -> list[dict]:
    """The processing-pipeline step(s) equivalent to estimate_baseline's
    correction, for recording onto a batch-fit recipe (see _entries())."""
    if kind not in BASELINE_OPS:
        return []
    op_name, kw_names = BASELINE_OPS[kind]
    step = {"op": op_name}
    if "order" in kw_names:
        step["order"] = order
    return [step]


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
                            release_frac=self.frac, iter_cb=cb, tol=self.tol,
                            should_stop=lambda: self._stop)
            self.done.emit(res, self._mode)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))


class _ErrorWorker(QThread):
    """Run a per-spectrum error analysis (covariance / Monte-Carlo / χ² profile)
    on a finished batch result, off the UI thread and interruptibly."""
    done = Signal(object)                       # the (mutated) result
    failed = Signal(str)
    progress = Signal(int, int, int, int)       # spectrum k, n, sub-step j, tot

    def __init__(self, result, data, method, n_trials, seed, n_points):
        super().__init__()
        self.result, self.data, self.method = result, data, method
        self.n_trials, self.seed, self.n_points = n_trials, seed, n_points
        self._stop = False

    def request_stop(self, mode: str = ""):
        self._stop = True

    def run(self):
        try:
            from larmor.batchfit import batch_error_analysis

            batch_error_analysis(
                self.result, self.data, method=self.method,
                n_trials=self.n_trials, seed=self.seed, n_points=self.n_points,
                progress=lambda k, n, j, tot: self.progress.emit(k, n, j, tot),
                should_stop=lambda: self._stop,
                parallel=True)
            self.done.emit(self.result)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))


class BatchFitDialog(QDialog):
    def __init__(self, parent, paths, model_recipe: dict | None):
        super().__init__(parent)
        self.setWindowTitle("Batch fit — one shared model, amplitudes per spectrum")
        self.resize(1000, 780)
        self._src_paths = [str(p) for p in (paths or [])]
        self._model_sites = ((model_recipe or {}).get("sites") or None)
        self._window = ((model_recipe or {}).get("fit_window_ppm") or None)
        self._recipe_tag = ""            # set when a model is loaded from a recipe
        self._result = None
        self._worker = None
        self._cells: list[dict] = []
        self._show_comp = False
        self._shared_scale = False
        self._baseline_kind = "None"
        self._excluded: dict[int, set[int]] = {}   # cell index -> excluded site indices
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

        # a prominent warning if the selection mixes nuclei (a shared model across
        # different nuclei is meaningless)
        self.warnBanner = QLabel("")
        self.warnBanner.setWordWrap(True)
        self.warnBanner.setStyleSheet(
            "background:#c0392b; color:white; font-weight:600; padding:4px 8px; "
            "border-radius:3px;")
        self.warnBanner.setVisible(False)
        v.addWidget(self.warnBanner)
        self._update_nuclei_warning()

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
        b_tsave = QPushButton("Save setup…")
        b_tsave.setToolTip("save this batch setup (release set, baseline, "
                           "threshold) as a reusable template")
        b_tsave.clicked.connect(self._save_template)
        b_tload = QPushButton("Load setup…")
        b_tload.clicked.connect(self._load_template)
        opt.addWidget(b_tsave)
        opt.addWidget(b_tload)
        v.addLayout(opt)

        # ---- release panel (per-parameter) ----
        rel = QHBoxLayout()
        rel.addWidget(QLabel("Release per spectrum (else held fixed):"))
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

        # ---- error calculation (the methods we have; export with the chosen one) ----
        er = QHBoxLayout()
        er.addWidget(QLabel("<b>Error calculation:</b>"))
        self.errCombo = QComboBox()
        self.errCombo.addItem("Covariance (from the fit)", "covariance")
        self.errCombo.addItem("Monte-Carlo (synthetic-noise refits)", "montecarlo")
        self.errCombo.addItem("χ² profile (error analysis)", "profile")
        self.errCombo.setToolTip(
            "how per-spectrum parameter errors are estimated:\n"
            "• Covariance — the least-squares covariance matrix (instant, from the fit)\n"
            "• Monte-Carlo — refit N synthetic noisy copies; captures correlations "
            "and non-linearity the covariance misses\n"
            "• χ² profile — scan each parameter and refit the rest; a real 1σ "
            "confidence interval")
        self.errCombo.currentIndexChanged.connect(self._on_err_method)
        er.addWidget(self.errCombo)
        self.errNlbl = QLabel("trials")
        er.addWidget(self.errNlbl)
        self.errN = QSpinBox(); self.errN.setRange(5, 5000); self.errN.setValue(200)
        self.errN.setToolTip("Monte-Carlo: number of synthetic refits per spectrum · "
                             "χ² profile: points scanned per parameter")
        er.addWidget(self.errN)
        self.btnErr = QPushButton("Compute errors")
        self.btnErr.setToolTip("estimate errors for every spectrum with the selected "
                               "method (writes them into each fit)")
        self.btnErr.setEnabled(False)
        self.btnErr.clicked.connect(self._compute_errors)
        er.addWidget(self.btnErr)
        self.btnErrCsv = QPushButton("Export CSV…")
        self.btnErrCsv.setToolTip("write a CSV of every fitted parameter with its "
                                  "value and error, using the SELECTED error method "
                                  "(computes it first if needed)")
        self.btnErrCsv.setEnabled(False)
        self.btnErrCsv.clicked.connect(self._export_csv)
        er.addWidget(self.btnErrCsv)
        self.errStatus = QLabel("")
        self.errStatus.setStyleSheet(f"color:{theme.active().text_dim};")
        er.addWidget(self.errStatus, 1)
        v.addLayout(er)
        self._err_worker = None
        self._export_after = False
        self._export_path = None
        self._on_err_method()

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

        autorow = QHBoxLayout()
        self.chkAutoRecipes = QCheckBox(
            "also save individual fits (.recipe.json) next to the CSV")
        self.chkAutoRecipes.setChecked(True)
        self.chkAutoRecipes.setToolTip(
            "so the Plotting studio's batch-grid finds the real saved fits "
            "next to a batch table automatically (highest-fidelity match, "
            "ahead of rebuilding from the CSV's own values) — auto-named, "
            "same as \"Save individual fits…\"'s automatic option")
        autorow.addStretch(1)
        autorow.addWidget(self.chkAutoRecipes)
        v.addLayout(autorow)
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
                "sample": sample_label(p, rec), "path": p,
                "proc": _proc_number(p), "snr": _snr(amp),
                "baseline_ops": []})   # per-spectrum manual baseline (2-point…)
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
                plot.hideAxis("left")
                plot.getPlotItem().invertX(True)     # NMR: ppm runs high → low
                plot.getPlotItem().getViewBox().setMouseEnabled(True, True)  # zoom
                plot.setMinimumHeight(120)
                exp = plot.plot(d["ppm"], d["amp"],
                                pen=pg.mkPen(theme.active().experiment, width=1))
                model = plot.plot([], [],
                                  pen=pg.mkPen(theme.active().model, width=1.4))
                plot.setXRange(d["ppm"].min(), d["ppm"].max())   # inverted → high→low
                plot.scene().sigMouseClicked.connect(
                    lambda ev, kk=k: self._cell_clicked(kk, ev))
                rmsd = QLabel(""); rmsd.setStyleSheet(
                    f"font-size:9px; color:{theme.active().text_dim};")
                cv.addWidget(title); cv.addWidget(plot, 1); cv.addWidget(rmsd)
                grid.addWidget(cell, j // 3, j % 3)
                self._cells.append({"plot": plot, "exp": exp, "model": model,
                                    "rmsd": rmsd, "comp": [], "title": title,
                                    "bl_picking": False, "bl_markers": [],
                                    "bl_line": None})
                self._attach_cell_menu(k)
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
            c.setToolTip(f"off: {pn} is held fixed at the recipe value. "
                         f"on: {pn} is fit per spectrum, ±{self.frac.value():.0f}% "
                         "around the recipe value")
            c.toggled.connect(lambda _=False: self.status.setText(self._model_status()))
            self._rel_checks[pn] = c
            self._rellay.addWidget(c)
        self._rellay.addStretch(1)

    def _update_nuclei_warning(self):
        nuclei = sorted({d["nucleus"] for d in self._data if d["nucleus"]})
        if len(nuclei) > 1:
            self.warnBanner.setText(
                "⚠ Mixed nuclei selected (" + ", ".join(nuclei) + "). A single "
                "shared model cannot describe different nuclei — select spectra of "
                "one nucleus for a meaningful batch fit.")
            self.warnBanner.setVisible(True)
        else:
            self.warnBanner.setVisible(False)

    # ------------------------------------------------------------------ status
    def _model_status(self) -> str:
        if not self._model_sites:
            return ("no model yet — load a recipe (or fit one spectrum first, "
                    "then reopen from a fit)")
        nuclei = {d["nucleus"] for d in self._data if d["nucleus"]}
        warn = ("  ⚠ mixed nuclei" if len(nuclei) > 1 else "")
        rel = [pn for pn, c in getattr(self, "_rel_checks", {}).items()
               if c.isChecked()]
        rtxt = (" · releasing " + ", ".join(rel)) if rel else ""
        return (f"{len(self._data)} spectra · model: {len(self._model_sites)} "
                f"line(s) · everything fixed at the recipe except amplitude"
                + rtxt + warn)

    def _model_label(self) -> str:
        return ("model: none — load a recipe" if not self._model_sites
                else f"model: {len(self._model_sites)} line(s) · held fixed except "
                     "amplitude (tick Release to let a parameter move)")

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
        # a new model may not share the old one's site indices/order --
        # any per-spectrum "Exclude component" picks no longer mean anything
        if self._excluded:
            self._excluded.clear()
            for k in range(len(self._cells)):
                self._update_exclude_title(k)
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
        had_manual = any(d.get("baseline_ops") for d in self._data)
        op = baseline_processing_op(kind, o)
        for k, d in enumerate(self._data):
            base = estimate_baseline(d["ppm"], d["amp0"], kind, o)
            d["amp"] = d["amp0"] - base
            # recorded onto the recipe (see _entries()) so a saved fit
            # reproduces this exact correction from the raw source later --
            # replaces rather than layers on any prior per-spectrum manual
            # pick, since this recomputes fresh from amp0 for every spectrum
            d["baseline_ops"] = list(op)
            if k < len(self._cells):
                self._cells[k]["exp"].setData(d["ppm"], d["amp"])
        self.lblBaseline.setText(f"{kind}" + (f" (order {o})"
                                              if kind == "Polynomial" else ""))
        msg = f"baseline subtracted per spectrum: {kind}"
        if had_manual:
            msg += "  (replaced per-spectrum manual baselines)"
        self.status.setText(msg)

    def _reset_baseline(self):
        self._baseline_kind = "None"
        for k in range(len(self._data)):
            self._end_bg_pick(k)        # abandon any in-progress manual pick
        for k, d in enumerate(self._data):
            d["amp"] = d["amp0"].copy()
            d["baseline_ops"] = []
            if k < len(self._cells):
                self._cells[k]["exp"].setData(d["ppm"], d["amp"])
        self.lblBaseline.setText("none")
        self.status.setText("baseline removed — using raw spectra")

    # ------------------------------------------------------------------ manual
    # per-spectrum 2-point linear baseline (right-click a cell's plot)
    def _attach_cell_menu(self, k: int):
        """(Re)build one cell's right-click menu: Export/Send-to-studio plus
        the per-spectrum baseline actions. Called at construction AND every
        time the menu is re-enabled after a pick, because pyqtgraph's
        setMenuEnabled(True) discards the old ViewBoxMenu and builds a fresh
        DEFAULT one from scratch — silently wiping any custom items unless
        they are re-added (this is also why the baseline options must never
        appear to "vanish": re-attaching keeps them available every time)."""
        from larmor.desktop.plot_menu import attach_plot_menu

        plot = self._cells[k]["plot"]
        attach_plot_menu(plot, title=self._data[k]["sample"], parent=self)
        vb_menu = plot.getPlotItem().getViewBox().menu
        vb_menu.addSeparator()
        act_bg = QAction("Add 2-point linear baseline", vb_menu)
        act_bg.setToolTip(
            "click two points on THIS spectrum (one each side of the peaks); "
            "the straight line through them is subtracted — right-click to "
            "cancel. Available again any time you right-click, including "
            "after applying one, so you can add another or replace it.")
        act_bg.triggered.connect(lambda _=False, kk=k: self._start_bg_pick(kk))
        vb_menu.addAction(act_bg)
        act_bg_clear = QAction("Clear this spectrum's baseline", vb_menu)
        act_bg_clear.triggered.connect(
            lambda _=False, kk=k: self._clear_cell_baseline(kk))
        vb_menu.addAction(act_bg_clear)

        if self._model_sites:
            vb_menu.addSeparator()
            excl_menu = vb_menu.addMenu("Exclude component")
            excl_menu.setToolTip(
                "force this component's amplitude to exactly zero for THIS "
                "spectrum only (e.g. a line that's only real in some samples) "
                "— excluded components are held at zero rather than fit, and "
                "are left out of the exported table/CSV and any plot built "
                "from it, not reported as a fitted zero")
            excluded = self._excluded.get(k, set())
            for i, site in enumerate(self._model_sites):
                label = site.get("label") or f"s{i}"
                act = QAction(label, excl_menu)
                act.setCheckable(True)
                act.setChecked(i in excluded)
                act.toggled.connect(
                    lambda on, kk=k, ii=i: self._toggle_exclude(kk, ii, on))
                excl_menu.addAction(act)

    def _start_bg_pick(self, k: int):
        cell = self._cells[k]
        if cell["bl_picking"]:
            return
        cell["bl_picking"] = True
        cell["bl_line"] = None
        cell["plot"].setCursor(Qt.PointingHandCursor)
        cell["plot"].getPlotItem().getViewBox().setMenuEnabled(False)
        self.status.setText(
            f"click two baseline points on “{self._data[k]['sample']}” — one "
            "each side of the peaks (drag either to adjust; right-click to "
            "cancel, or to Apply once both are placed)")

    def _cell_clicked(self, k: int, ev):
        cell = self._cells[k]
        if not cell["bl_picking"]:
            return
        if ev.button() == Qt.RightButton:
            ev.accept()
            if len(cell["bl_markers"]) >= 2:
                self._confirm_bg_pick(k, ev)   # both points placed: Apply / Cancel
            else:
                self._end_bg_pick(k)
                self.status.setText("2-point background cancelled")
            return
        if ev.button() != Qt.LeftButton:
            return
        if len(cell["bl_markers"]) >= 2:
            return                              # drag an existing point instead
        plot = cell["plot"]
        if not plot.sceneBoundingRect().contains(ev.scenePos()):
            return
        vb = plot.getPlotItem().getViewBox()
        p = vb.mapSceneToView(ev.scenePos())
        marker = pg.TargetItem(
            pos=(float(p.x()), float(p.y())), size=11, movable=True,
            pen=pg.mkPen(theme.active().baseline, width=1.5),
            brush=pg.mkBrush(255, 255, 255, 220))
        marker.sigPositionChanged.connect(lambda *_, kk=k: self._update_bg_preview(kk))
        plot.addItem(marker)
        cell["bl_markers"].append(marker)
        ev.accept()
        if len(cell["bl_markers"]) == 2:
            self._update_bg_preview(k)
            self.status.setText(
                f"drag either point on “{self._data[k]['sample']}” to adjust — "
                "right-click to Apply or Cancel")

    def _bg_points(self, k: int):
        """Current (possibly dragged) positions of the two picked markers."""
        return [(float(m.pos().x()), float(m.pos().y()))
                for m in self._cells[k]["bl_markers"]]

    def _update_bg_preview(self, k: int):
        cell = self._cells[k]
        if len(cell["bl_markers"]) < 2:
            return
        (x1, y1), (x2, y2) = self._bg_points(k)
        if cell["bl_line"] is None:
            cell["bl_line"] = cell["plot"].plot(
                [x1, x2], [y1, y2],
                pen=pg.mkPen(theme.active().baseline, width=1.2, style=Qt.DashLine))
        else:
            cell["bl_line"].setData([x1, x2], [y1, y2])

    def _confirm_bg_pick(self, k: int, ev):
        """Right-click with both points placed: a tiny Apply/Cancel menu, so a
        bad click can be fixed (drag) or the whole pick abandoned, instead of
        committing the moment the second point lands."""
        cell = self._cells[k]
        pos = cell["plot"].mapToGlobal(cell["plot"].mapFromScene(ev.scenePos()))
        choice = self._ask_apply_or_cancel(cell["plot"], pos)
        if choice == "apply":
            self._apply_bg_pick(k)
        elif choice == "cancel":
            self._end_bg_pick(k)
            self.status.setText("2-point background cancelled")

    def _ask_apply_or_cancel(self, plot, global_pos) -> str:
        """Show the tiny confirm menu; returns "apply" or "cancel". Isolated
        from _confirm_bg_pick so tests can stub the UI without exec()'ing a
        real native menu (which blocks under Qt regardless of monkeypatching
        QMenu.exec — it's a wrapped virtual, not a plain Python attribute)."""
        from PySide6.QtWidgets import QMenu

        menu = QMenu(plot)
        act_apply = menu.addAction("Apply this baseline")
        menu.addAction("Cancel")
        chosen = menu.exec(global_pos)
        return "apply" if chosen is act_apply else "cancel"

    def _apply_bg_pick(self, k: int):
        cell = self._cells[k]
        d = self._data[k]
        (x1, y1), (x2, y2) = self._bg_points(k)
        if abs(x2 - x1) < 1e-9:
            self.status.setText(
                "the two points are at the same position — drag one apart, "
                "then right-click to Apply")
            return                               # keep picking; don't discard
        m_ = (y2 - y1) / (x2 - x1)          # line through the two picked points
        base = m_ * (d["ppm"] - x1) + y1
        d["amp"] = d["amp"] - base
        d["baseline_ops"].append(
            {"op": "twopoint_bg", "x1": x1, "y1": y1, "x2": x2, "y2": y2})
        cell["exp"].setData(d["ppm"], d["amp"])
        self._end_bg_pick(k)
        self.status.setText(
            f"2-point background subtracted on “{d['sample']}” "
            "(carried into its saved fit / recipe)")

    def _end_bg_pick(self, k: int):
        cell = self._cells[k]
        cell["bl_picking"] = False
        for m in cell["bl_markers"]:
            cell["plot"].removeItem(m)
        cell["bl_markers"] = []
        if cell["bl_line"] is not None:
            cell["plot"].removeItem(cell["bl_line"])
            cell["bl_line"] = None
        cell["plot"].unsetCursor()
        cell["plot"].getPlotItem().getViewBox().setMenuEnabled(True)
        self._attach_cell_menu(k)   # re-enabling built a fresh blank menu

    def _clear_cell_baseline(self, k: int):
        self._end_bg_pick(k)
        d = self._data[k]
        d["amp"] = d["amp0"].copy()
        d["baseline_ops"] = []
        if k < len(self._cells):
            self._cells[k]["exp"].setData(d["ppm"], d["amp"])
        self.status.setText(
            f"baseline cleared for “{d['sample']}” (back to the raw spectrum)")

    # ------------------------------------------------------------------ exclusion
    def _toggle_exclude(self, k: int, site_idx: int, on: bool):
        s = self._excluded.setdefault(k, set())
        if on:
            s.add(site_idx)
        else:
            s.discard(site_idx)
            if not s:
                self._excluded.pop(k, None)
        self._update_exclude_title(k)
        d = self._data[k]
        label = (self._model_sites[site_idx].get("label") or f"s{site_idx}"
                if self._model_sites else f"s{site_idx}")
        self.status.setText(
            (f"“{label}” excluded for " if on else f"“{label}” restored for ")
            + f"“{d['sample']}” (locked to zero amplitude; left out of the "
              "exported table and any plot built from it)")

    def _update_exclude_title(self, k: int):
        if k >= len(self._cells):
            return
        d = self._data[k]
        excluded = sorted(self._excluded.get(k, ()))
        text = d["sample"]
        if excluded and self._model_sites:
            labels = [self._model_sites[i].get("label") or f"s{i}"
                     for i in excluded if i < len(self._model_sites)]
            text += "  (excluded: " + ", ".join(labels) + ")"
        self._cells[k]["title"].setText(text)

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
                vb.setXRange(xmin, xmax, padding=0.02)   # invertX handles direction
                vb.setYRange(ymin, ymax, padding=0.05)
        else:
            for k, cell in enumerate(self._cells):
                d = self._data[k]
                vb = cell["plot"].getViewBox()
                vb.enableAutoRange(axis="y")
                vb.setXRange(float(d["ppm"].min()), float(d["ppm"].max()),
                             padding=0.02)

    # ------------------------------------------------------------------ fit
    def _entries(self):
        from larmor.recipe import Recipe

        out = []
        for k, d in enumerate(self._data):
            rec = Recipe.from_dict({
                "nucleus": d["nucleus"], "larmor_frequency_MHz": d["larmor"],
                "spin_rate_Hz": d["spin"], "sample": d["sample"],
                # whatever baseline correction is active -- the global "Fit
                # baseline…" tool (see BASELINE_OPS) or a per-spectrum manual
                # 2-point pick -- recorded so the exported fit reproduces the
                # EXACT corrected spectrum against the raw source later (the
                # Plotting studio's batch-grid replays this for "experiment")
                "processing": list(d.get("baseline_ops") or []),
                "source_path": d.get("path", ""),
                "sites": copy.deepcopy(self._model_sites)})
            # "Exclude component" (right-click a cell): lock that site's
            # amplitude at exactly zero for THIS spectrum only -- batchfit.
            # free_amplitudes() recognises and preserves this lock (it's not
            # "an amplitude that fit near zero", it's "not part of this fit")
            for i in self._excluded.get(k, ()):
                if i < len(rec.sites):
                    amp = rec.sites[i].params.get("amplitude")
                    if amp is not None:
                        amp.value, amp.vary = 0.0, False
                        amp.min, amp.max = 0.0, 0.0
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
        self.btnErr.setEnabled(False); self.btnErrCsv.setEnabled(False)
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
        if self._err_worker is not None and self._err_worker.isRunning():
            self._err_worker.request_stop(mode)
            self.status.setText("stopping error analysis — keeping what's computed…")
            return
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
        self.btnErr.setEnabled(True); self.btnErrCsv.setEnabled(True)
        self._refresh_err_status()
        for k, pd in enumerate(result.per_dataset):
            if k < len(self._cells):
                self._cells[k]["model"].setData(np.asarray(pd["x"], float),
                                                np.asarray(pd["y_fit"], float))
                self._cells[k]["rmsd"].setText(f"RMSD {result.rmsd[k]:.4f}")
        self._refresh_components()
        flagged = self._flag_quality(result)
        note = " (stopped early)" if mode == "stop" else ""
        if flagged:
            note += f"  ⚠ {len(flagged)} spectrum/spectra flagged (RMSD outlier " \
                    "or low S/N) — hover the cells"
        self.status.setText(result.summary + note)

    def _flag_quality(self, result) -> list[int]:
        """Flag spectra whose RMSD is an outlier or whose S/N is low, colouring
        their RMSD label red with a reason tooltip (a per-spectrum quality gate)."""
        r = np.asarray(result.rmsd, float)
        med = float(np.median(r))
        mad = float(np.median(np.abs(r - med))) or (0.1 * med + 1e-9)
        hi = med + 3.0 * 1.4826 * mad
        flagged = []
        for k, cell in enumerate(self._cells):
            if k >= len(result.rmsd):
                break
            reasons = []
            if r[k] > hi:
                reasons.append(f"RMSD {r[k]:.3g} is an outlier (>{hi:.3g})")
            snr = self._data[k].get("snr", float("inf"))
            if snr < 20:
                reasons.append(f"low S/N ({snr:.0f})")
            lbl = cell["rmsd"]
            if reasons:
                flagged.append(k)
                lbl.setText(f"⚠ RMSD {r[k]:.4f}")
                lbl.setStyleSheet("font-size:9px; color:#c0392b; font-weight:600;")
                lbl.setToolTip("; ".join(reasons))
            else:
                lbl.setStyleSheet(
                    f"font-size:9px; color:{theme.active().text_dim};")
                lbl.setToolTip(f"S/N ≈ {snr:.0f}" if np.isfinite(snr) else "")
        return flagged

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

    def _save_all_recipes_to(self, folder: Path) -> int:
        """Write one auto-named LARMOR .recipe.json per fitted spectrum into
        `folder` -- the core of "Save individual fits…"'s automatic-naming
        path, also called right after a CSV export (see the "also save
        individual fits…" checkbox) so a batch table normally comes with its
        own recipes sitting right next to it: series_grid.find_recipes_near_csv
        finds them by sample name for free, giving the Plotting studio's
        batch-grid the real saved fit (bounds, vary flags, processing/
        baseline included) instead of falling back to rebuilding a value-only
        recipe from the CSV's own rows."""
        from larmor.recipe import Recipe
        if self._result is None:
            return 0
        n = 0
        for k, rec in enumerate(self._result.recipes):
            proc = self._data[k]["proc"] if k < len(self._data) else ""
            name = self._auto_name(rec, proc)
            try:
                Recipe.from_dict(rec.to_dict()).save(folder / f"{name}.recipe.json")
                n += 1
            except Exception:
                pass
        return n

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
        if mode == QMessageBox.Yes:
            n = self._save_all_recipes_to(folder)
        else:
            recs = self._result.recipes
            n = 0
            for k, rec in enumerate(recs):
                proc = self._data[k]["proc"] if k < len(self._data) else ""
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

        from larmor.desktop.paths import suggest_save_dir
        start = suggest_save_dir(self._src_paths[0] if self._src_paths else None)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save batch table",
            str(Path(start) / "batch_table.csv") if start else "batch_table.csv",
            "CSV (*.csv)")
        if not path:
            return
        rows = batchfit.shared_table(self._result)
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            # model + source_path: so the Plotting studio's batch-grid figure
            # can find each row's spectrum/fit straight from this CSV, even
            # without "Save individual fits…" too
            w.writerow(["scope", "site", "label", "param", "value", "stderr",
                        "model", "source_path"])
            for r in rows:
                w.writerow([r["scope"], r["site"], r["label"], r["param"],
                            f"{r['value']:.6g}",
                            "" if r["stderr"] is None else f"{r['stderr']:.4g}",
                            r.get("model", ""), r.get("source_path", "")])
        msg = f"saved {path}"
        if self.chkAutoRecipes.isChecked():
            n = self._save_all_recipes_to(Path(path).parent)
            msg += f" · {n} individual fit(s) saved alongside it"
        self.status.setText(msg)

    def _save_template(self):
        import json
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QInputDialog
        name, ok = QInputDialog.getText(self, "Save batch setup", "Template name:")
        if not ok or not name.strip():
            return
        tpl = {"release": [pn for pn, c in self._rel_checks.items() if c.isChecked()],
               "frac": self.frac.value(), "tol": self.tol.value(),
               "baseline": self._baseline_kind,
               "components": self.chkComp.isChecked(),
               "shared_scale": self.chkShared.isChecked()}
        s = QSettings("LARMOR", "app")
        lib = json.loads(s.value("batchTemplates", "{}") or "{}")
        lib[name.strip()] = tpl
        s.setValue("batchTemplates", json.dumps(lib))
        self.status.setText(f"saved batch setup “{name.strip()}”")

    def _load_template(self):
        import json
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QInputDialog
        lib = json.loads(QSettings("LARMOR", "app").value(
            "batchTemplates", "{}") or "{}")
        if not lib:
            self.status.setText("no saved batch setups yet"); return
        names = list(lib.keys())
        choice, ok = QInputDialog.getItem(self, "Load batch setup", "Template:",
                                          names, 0, False)
        if not ok:
            return
        tpl = lib[choice]
        for pn, c in self._rel_checks.items():
            c.setChecked(pn in tpl.get("release", []))
        self.frac.setValue(tpl.get("frac", 10))
        self.tol.setValue(tpl.get("tol", 0.0))
        self.chkComp.setChecked(tpl.get("components", False))
        self.chkShared.setChecked(tpl.get("shared_scale", False))
        self.status.setText(f"loaded setup “{choice}” — Fit to apply")

    def _series_plot(self):
        if self._result is None:
            return
        from larmor.desktop.series_plot import SeriesPlotDialog
        SeriesPlotDialog(self, self._result).exec()

    # ------------------------------------------------------------------ errors
    def _on_err_method(self):
        m = self.errCombo.currentData()
        prof = (m == "profile")
        self.errNlbl.setText("points" if prof else "trials")
        self.errNlbl.setVisible(m != "covariance")
        self.errN.setVisible(m != "covariance")
        if prof:
            self.errN.setRange(5, 61)
            if self.errN.value() > 61:
                self.errN.setValue(15)
        elif m == "montecarlo":
            self.errN.setRange(5, 5000)
            if self.errN.value() < 40:
                self.errN.setValue(200)
        self._refresh_err_status()

    def _refresh_err_status(self):
        if self._result is None:
            self.errStatus.setText("")
            return
        sel = self.errCombo.currentData()
        have = sel in getattr(self._result, "error_detail", {})
        name = self.errCombo.currentText().split(" (")[0]
        self.errStatus.setText(
            f"{name} ready to export" if have
            else f"Export will compute {name} first")

    def _compute_errors(self):
        if self._result is None:
            return
        # covariance is NOT a free snapshot: batch_fit's initial pass skips
        # the errorbar-rescue retry for speed, so Param.stderr is commonly
        # None until this actually refits with compute_errorbars=True -- so
        # it goes through the same threaded worker as Monte-Carlo/profile,
        # not a synchronous "just read what's there" shortcut.
        self._start_err_worker(self.errCombo.currentData())

    def _export_csv(self):
        if self._result is None:
            return
        m = self.errCombo.currentData()
        have = m in getattr(self._result, "error_detail", {})
        from larmor.desktop.paths import suggest_save_dir
        start = suggest_save_dir(self._src_paths[0] if self._src_paths else None)
        path, _ = QFileDialog.getSaveFileName(
            self, "Export fit table with errors",
            str(Path(start) / f"batch_fit_{m}.csv") if start else f"batch_fit_{m}.csv",
            "CSV (*.csv)")
        if not path:
            return
        if not have:                     # compute the selected method, then export
            self._export_path = path
            self._export_after = True
            self._start_err_worker(m)
        else:
            self._result.error_method = m
            self._write_err_csv(path)

    def _start_err_worker(self, method: str):
        data = [(d["ppm"], d["amp"], self._window) for d in self._data]
        n = self.errN.value()
        self.btnFit.setEnabled(False)
        self.btnErr.setEnabled(False); self.btnErrCsv.setEnabled(False)
        self.btnSave.setEnabled(False); self.btnTable.setEnabled(False)
        self.btnSeries.setEnabled(False)
        self.btnStop.setEnabled(True); self.btnCancel.setEnabled(False)
        self.prog.setRange(0, len(data)); self.prog.setValue(0)
        name = self.errCombo.currentText().split(" (")[0]
        self.status.setText(f"computing {name} errors…")
        self._err_worker = _ErrorWorker(self._result, data, method, n, 0, n)
        self._err_worker.progress.connect(self._err_progress)
        self._err_worker.done.connect(self._err_done)
        self._err_worker.failed.connect(self._err_failed)
        self._err_worker.start()

    def _err_progress(self, k, n, j, tot):
        self.prog.setValue(min(k, n))
        self.status.setText(f"error analysis — spectrum {min(k + 1, n)}/{n}"
                            + (f" · step {j}/{tot}" if tot > 1 else ""))

    def _post_err_enable(self):
        self.prog.setRange(0, 100)
        self.btnStop.setEnabled(False)
        self._update_fit_enabled()
        self.btnErr.setEnabled(True); self.btnErrCsv.setEnabled(True)
        self.btnSave.setEnabled(True); self.btnTable.setEnabled(True)
        self.btnSeries.setEnabled(True)

    def _err_done(self, result):
        self._post_err_enable()
        self.prog.setValue(100)
        method = getattr(result, "error_method", "covariance")
        # a worker completing "successfully" does NOT mean the numbers are
        # any good -- every refit/scan can fail or come back degenerate
        # (silently, one exception at a time) while the worker itself still
        # finishes and reports done. Count how many parameter-spectrum
        # combinations actually got a usable error before claiming success.
        detail = (getattr(result, "error_detail", {}) or {}).get(method) or []
        total = sum(len(d) for d in detail)
        ok = sum(1 for d in detail for pe in d.values()
                if pe.stderr is not None and np.isfinite(pe.stderr))
        if detail and total and ok == 0:
            self.status.setText(
                f"⚠ {method} produced NO usable errors for any of "
                f"{len(result.recipes)} spectra — every refit/scan failed or "
                "was degenerate (a near-zero-amplitude site, an unidentifiable "
                "released parameter, or too few points for the model). Try "
                "fewer released parameters, a looser release %, or check the "
                "flagged/high-RMSD spectra's fit quality first.")
        elif detail and total and ok < total:
            self.status.setText(
                f"{method} errors: {ok}/{total} parameter-spectrum "
                f"combinations succeeded ({len(result.recipes)} spectra) — "
                "some failed or were degenerate; export still writes what "
                "succeeded (blank for what didn't)")
        else:
            self.status.setText(
                f"{method} errors computed for {len(result.recipes)} spectra")
        self._refresh_err_status()
        if self._export_after:
            self._export_after = False
            self._write_err_csv(self._export_path)

    def _err_failed(self, msg):
        self._post_err_enable()
        self.prog.setValue(0)
        self._export_after = False
        self.status.setText(f"error analysis failed: {msg}")

    def _write_err_csv(self, path):
        if not path or self._result is None:
            return
        import csv
        from larmor import batchfit

        method = self.errCombo.currentData()

        def num(v, fmt=".6g"):
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return ""
            return format(v, fmt)

        rows = batchfit.error_table(self._result, method=method)
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            # model + source_path: so the Plotting studio's batch-grid figure
            # can find each row's spectrum/fit straight from this CSV, even
            # without "Save individual fits…" too
            w.writerow(["scope", "site", "label", "param", "value", "stderr",
                        "error_method", "sigma_pct", "ci68_lo", "ci68_hi",
                        "model", "source_path"])
            for r in rows:
                w.writerow([r["scope"], r["site"], r["label"], r["param"],
                            num(r["value"]), num(r["stderr"]), r["error_method"],
                            num(r["sigma_pct"], ".3g"),
                            num(r["ci68_lo"]), num(r["ci68_hi"]),
                            r.get("model", ""), r.get("source_path", "")])
        msg = f"exported {Path(path).name} · {method} errors"
        if self.chkAutoRecipes.isChecked():
            n = self._save_all_recipes_to(Path(path).parent)
            msg += f" · {n} individual fit(s) saved alongside it"
        self.status.setText(msg)

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "multi-dataset", "Multi-dataset & co-fitting")


# ---------------------------------------------------------------- module helpers
def _slug(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in (s or ""))[:80]


def sample_label(path, rec) -> str:
    """A meaningful sample name for a spectrum: the recipe's sample if it is not
    just the nucleus, else the sample **folder** derived from the path (the first
    ancestor that is not a proc/expno number) — so a title of "31P" becomes the
    real sample directory name."""
    nucleus = (rec.get("nucleus") or "").strip()
    name = (rec.get("sample") or "").strip()
    if name and name.lower() != nucleus.lower():
        return name
    for seg in reversed(Path(path).parts):
        low = seg.lower()
        if low.endswith(".fid"):                 # a Varian dataset folder
            return seg[:-4]
        if seg == "pdata" or seg.isdigit() or low in ("1r", "2rr", "fid", "ser"):
            continue
        return seg
    return name or Path(path).stem


def _snr(amp) -> float:
    """Crude signal-to-noise: peak ÷ RMS of the quiet spectrum edges."""
    a = np.asarray(amp, float)
    n = max(3, a.size // 20)
    noise = float(np.std(np.concatenate([a[:n], a[-n:]]))) or 1.0
    return float(np.max(np.abs(a)) / noise)


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
        return float(QSettings("LARMOR", "app").value("fitStdevPct", 0.1) or 0.0)
    except (TypeError, ValueError):
        return 0.1


def _save_tol(v: float):
    from PySide6.QtCore import QSettings
    QSettings("LARMOR", "app").setValue("fitStdevPct", float(v))
