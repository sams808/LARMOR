"""Batch fit dialog: one shared model fitted to many 1D spectra at once.

Ctrl/Shift-select spectra in the Explorer, hit "Batch fit selected…". This shows
them in a 3×3 grid (tabs for more), fits them with all parameters except
amplitude shared, and can then **release** chosen parameters — per-parameter — to
drift a little per spectrum. Front-end for larmor.batchfit; the fit runs in an
interruptible worker (Cancel = revert, Stop = keep the last iteration).

Made to be a proper workbench: NMR-style axes (high→low ppm), mouse zoom,
independent-vs-shared scale, live component curves, an optional per-fit baseline,
a completion threshold, and one-click saving of every fit in LARMOR format.
A results table under the grid (one row per spectrum, sortable) is linked both
ways to the cells: selecting a row spotlights its spectrum, clicking a spectrum
finds its row.
"""
from __future__ import annotations

import copy
import datetime as _dt
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QAction, QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
    QDoubleSpinBox, QFileDialog, QGridLayout, QHBoxLayout, QInputDialog, QLabel,
    QMessageBox, QProgressBar, QPushButton, QScrollArea, QSpinBox, QSplitter,
    QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from larmor import comparability
from larmor.desktop import theme
from larmor.desktop.comparability_dialog import (
    ComparabilityBar, ComparabilityDialog, ReprocessDialog, _CHECK_CSS_COLOR,
    apply_reprocess, revert_reprocess,
)
from larmor.desktop.panels import PARAM_LABELS
from larmor.desktop.plot import site_color
from larmor.io.scan import disambiguate, sample_label
from larmor.series_table import SeriesTable

PER_TAB = 9        # 3×3 grid per tab

#: a cell's idle sample-name label; _apply_highlight restores it after a spotlight
_TITLE_CSS = "font-size:10px; font-weight:600;"
#: the flagged-RMSD red -- the same literal as the mixed-nuclei banner and the
#: cell label (a warning must read as one on every theme, so not a theme role).
#: Two levels, the fit-health strip's: red = bad (RMSD outlier, mixed nuclei),
#: amber (comparability_dialog._CHECK_CSS_COLOR) = check (a spectrum acquired
#: or processed differently from the series majority)
_FLAG_CSS_COLOR = "#c0392b"


class _NumItem(QTableWidgetItem):
    """QTableWidgetItem sorts its text lexicographically ("9.5" > "10.2");
    this one compares the float stored under ``Qt.UserRole + 1`` so a header
    click orders numbers. A missing or non-finite key sorts after every
    number, so blank cells sink to the bottom of an ascending sort."""

    def __lt__(self, other):
        return _sort_key(self) < _sort_key(other)


def _sort_key(item) -> float:
    v = item.data(Qt.UserRole + 1)
    try:
        f = float(v)
    except (TypeError, ValueError):
        return float("inf")
    return f if np.isfinite(f) else float("inf")


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
        self.resize(1000, 860)
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
        self._hl: int | None = None          # spotlighted spectrum (table row <-> grid cell)
        self._flag_reasons: dict[int, str] = {}   # k -> "RMSD … is an outlier; low S/N"
        self._reprocess: dict | None = None  # {"template", "phase"} once reprocessed from fid
        self._data = self._load(paths)
        # the series' identity (larmor.series_table): names, replicate groups,
        # order and joined composition columns -- edited by "Series table…"
        self._series = SeriesTable.from_spectra(self._data)

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

        # the comparability line: were these spectra acquired and processed
        # alike? (acqus / procs / auditp against the series majority; hidden
        # for CSV/fxmla-only series). Details… opens the table, Reprocess all
        # from fid… rebuilds every member with one common pipeline.
        self.compBar = ComparabilityBar()
        self.compBar.details_requested.connect(self._show_comparability)
        self.compBar.reprocess_requested.connect(self._reprocess_all)
        self.compBar.revert_requested.connect(self._use_topspin_processing)
        self.compBar.set_comparison(self._comparison, self._n_with_fid())
        v.addWidget(self.compBar)

        # ---- the grid of spectra, the results table under it ----
        split = QSplitter(Qt.Vertical)
        self.tabs = QTabWidget()
        split.addWidget(self.tabs)
        self._build_grid()
        self._build_table()                  # self.table: fixed columns until a fit
        split.addWidget(self.table)
        split.setCollapsible(0, False)       # the grid can never be dragged shut
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 1)
        split.setSizes([560, 160])
        v.addWidget(split, 1)

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
        self.chkTable = QCheckBox("table")
        self.chkTable.setChecked(True)
        self.chkTable.setToolTip(
            "show the results table under the grid — one row per spectrum: "
            "click a row to spotlight its spectrum (the others dim), click a "
            "spectrum to find its row; sort any column to bring an outlier to "
            "the top")
        self.chkTable.toggled.connect(self._toggle_table)
        opt.addWidget(self.chkComp)
        opt.addWidget(self.chkShared)
        opt.addWidget(self.chkTable)
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
                                  "(computes it first if needed), plus its vary / "
                                  "min / max / expr / at_bound status columns")
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
        self.btnTable.setToolTip("write a batch_table.csv of shared / per-spectrum "
                                 "values — the same numbers as the table under the "
                                 "grid, in long form, with each row's vary / min / "
                                 "max / expr / at_bound status columns")
        self.btnTable.setEnabled(False)
        self.btnTable.clicked.connect(self._save_table)
        self.btnBundle = bb.addButton("Publication bundle…", QDialogButtonBox.ActionRole)
        self.btnBundle.setToolTip(
            "write everything about this batch to one folder: batch_table.csv "
            "(+ the error table if computed), one .recipe.json and one _curves.csv "
            "per spectrum (ppm, experiment exactly as fitted, model, residual, "
            "components), manifest.csv (source file + SHA-256, EXPNO, NS, D1, SF/SR, "
            "processing, fit window, RMSD) and README.txt with the Methods paragraph "
            "and software versions")
        self.btnBundle.setEnabled(False)
        self.btnBundle.clicked.connect(self._export_bundle)
        self.btnSeries = bb.addButton("Series plot…", QDialogButtonBox.ActionRole)
        self.btnSeries.setToolTip("plot how each parameter evolves along the series")
        self.btnSeries.setEnabled(False)
        self.btnSeries.clicked.connect(self._series_plot)
        self.btnSeriesTable = bb.addButton("Series table…", QDialogButtonBox.ActionRole)
        self.btnSeriesTable.setToolTip(
            "names, replicate groups and order of the series (the grid, the results "
            "table, the CSV rows and the saved recipes follow) and composition "
            "columns joined from a CSV for the Series plot's x axis")
        self.btnSeriesTable.clicked.connect(self._edit_series)
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
        # No auto-default button anywhere in the dialog. Enter/Return inside
        # the table is ignored by QAbstractItemView and lands on the dialog,
        # which clicks whichever button Qt promoted to default on show -- the
        # FIRST auto-default push button in the focus chain ("Model from
        # recipe…" here; Close inside the button box) -- so results could be
        # discarded or a file dialog opened by a stray Enter. Much likelier
        # now that a cell click hands the table focus. Space still presses a
        # focused button.
        for b in self.findChildren(QPushButton):
            b.setAutoDefault(False)
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
                # the TopSpin arrays, kept for "Use TopSpin processing" after a
                # reprocess from fid
                "ppm_src": ppm.copy(), "amp_src": amp.copy(),
                "nucleus": rec.get("nucleus", ""),
                "larmor": float(rec.get("larmor_frequency_MHz", 0.0) or 0.0),
                "spin": float(rec.get("spin_rate_Hz", 0.0) or 0.0),
                "sample": sample_label(p, rec), "path": p,
                # where the name came from (io/scan through the loader): the
                # sample folder and the title's first line, for the Series table
                "folder": str((rec.get("provenance") or {}).get("sample_folder", "")),
                "title": str((rec.get("provenance") or {}).get("title", "")),
                "proc": _proc_number(p), "snr": _snr(amp),
                # acqus / procs / auditp (None for CSV / fxmla), the reprocess
                # chain and its EXPNO once "Reprocess all from fid…" ran
                "params": comparability.read_params(p),
                "proc_ops": [], "expno": "",
                "baseline_ops": []})   # per-spectrum manual baseline (2-point…)
        # two folders of one glass (or two EXPNOs of one folder) must never
        # share a label: labels become table scopes and recipe file names
        for d, lab in zip(data, disambiguate([d["sample"] for d in data],
                                             [d["path"] for d in data])):
            d["group"] = d["sample"]             # the base label: replicates share it
            d["sample"] = lab
            if self._model_sites is None and rec.get("sites"):
                self._model_sites = rec["sites"]
        self._comparison = comparability.compare(
            [d["params"] for d in data], [d["sample"] for d in data])
        return data

    def _recompare(self):
        """The series compared as it stands now: a member rebuilt from its
        fid reads with the common template's shape keys (and no abs), so
        only members still carrying their TopSpin processing stay flagged."""
        params = []
        for d in self._data:
            p = d.get("params")
            if p is not None and d.get("proc_ops") and self._reprocess:
                p = comparability.as_reprocessed(p, self._reprocess["template"])
            params.append(p)
        return comparability.compare(params, [d["sample"] for d in self._data])

    def _n_with_fid(self) -> int:
        return sum(1 for d in self._data
                   if d.get("params") is not None and d["params"].has_fid)

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
                title.setStyleSheet(_TITLE_CSS)
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
                                    "title_css": _TITLE_CSS,
                                    "bl_picking": False, "bl_markers": [],
                                    "bl_line": None})
                self._attach_cell_menu(k)
            self.tabs.addTab(page, f"{start + 1}–{min(start + PER_TAB, n)}")
        self._refresh_cell_marks()

    def _refresh_cell_marks(self):
        """The comparability marker on the cells: a deviating spectrum's
        title ends with ⚠ in amber and its tooltip says what differs from
        the series majority -- the RMSD flag's discoverability path, in the
        'check' colour. Re-run after a reprocess / revert."""
        for k, cell in enumerate(self._cells):
            sig = self._comparison.signature(k)
            ink = (_FLAG_CSS_COLOR if self._comparison.level_of(k) == "bad"
                   else _CHECK_CSS_COLOR)
            cell["title_css"] = _TITLE_CSS + (f" color:{ink};" if sig else "")
            title = cell.get("title")
            if title is not None and k < len(self._data):
                d = self._data[k]
                tip = (f"{d['sample']}"
                       + (f" · proc {d['proc']}" if d.get("proc") else "")
                       + f" · {d.get('nucleus', '')}")
                if sig:
                    tip += "\n" + "\n".join(self._comparison.details(k))
                title.setToolTip(tip)
            self._update_exclude_title(k)
        self._apply_highlight()

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
        """Estimate and subtract a baseline from every spectrum's amp0. amp0
        is the TopSpin spectrum until "Reprocess all from fid…" replaces it
        (apply_reprocess resets amp0), so after a reprocess the baseline is
        estimated on -- and Reset restores -- the reprocessed spectrum, not
        the 1r; "Use TopSpin processing" on the comparability bar does that."""
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
        """Baseline off: back to amp0 (the reprocessed spectrum after a
        reprocess from fid, else the TopSpin spectrum)."""
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

    # ------------------------------------------------------------------ comparability
    def _show_comparability(self):
        ComparabilityDialog(self, self._comparison).exec()

    def _reprocess_all(self):
        dlg = ReprocessDialog(self, self._comparison, self._n_with_fid())
        if dlg.exec() == QDialog.Accepted:
            self._apply_reprocess(dlg.template(), dlg.phase())

    def _clear_result(self):
        """The fit describes spectra that no longer exist (reprocessed or
        reverted): drop it and disable what depends on it, exactly as _run
        does until _done."""
        self._result = None
        self._flag_reasons.clear()
        t = theme.active()
        for cell in self._cells:
            cell["model"].setData([], [])
            for it in cell["comp"]:
                cell["plot"].removeItem(it)
            cell["comp"] = []
            cell["rmsd"].setText("")
            cell["rmsd"].setStyleSheet(f"font-size:9px; color:{t.text_dim};")
            cell["rmsd"].setToolTip("")
        for b in (self.btnSave, self.btnTable, self.btnBundle, self.btnSeries,
                  self.btnErr, self.btnErrCsv):
            b.setEnabled(False)
        self._refresh_err_status()

    def _apply_reprocess(self, template: dict, phase: str = "own"):
        """Rebuild every member with a fid through comparability.reprocess
        (loader.apply_processing: the path a saved recipe replays), record
        the chain per spectrum, drop the stale fit and baselines, and turn
        the comparability line into the applied pipeline. The testable core
        of "Reprocess all from fid…"."""
        n, notes = apply_reprocess(self._data, self._cells, template, phase,
                                   self.status.setText)
        self._reprocess = {"template": dict(template), "phase": phase}
        self._baseline_kind = "None"
        self.lblBaseline.setText("none")
        self._clear_result()
        self._comparison = self._recompare()
        self._refresh_cell_marks()
        self.compBar.set_reprocessed(template, n, len(self._data), phase)
        self._fill_table(None)
        self._apply_scale()
        ph = {"own": "each keeps its own phase", "auto": "autophased",
              "none": "unphased"}.get(phase, phase)
        msg = (f"reprocessed {n} spectra from fid with one common pipeline ({ph}) "
               "— Fit again; Fit baseline… still applies")
        if notes:
            msg += "  ⚠ " + "; ".join(notes)
        self.status.setText(msg)

    def _use_topspin_processing(self):
        """Back to the 1r arrays (ppm_src / amp_src) and the original flags."""
        revert_reprocess(self._data, self._cells)
        self._reprocess = None
        self._comparison = self._recompare()
        self._baseline_kind = "None"
        self.lblBaseline.setText("none")
        self._clear_result()
        self._refresh_cell_marks()
        self.compBar.set_comparison(self._comparison, self._n_with_fid())
        self._fill_table(None)
        self._apply_scale()
        self.status.setText("spectra as TopSpin processed them (1r) — Fit again")

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
            # a plain left click spotlights this spectrum and finds its table
            # row -- never a drag-zoom (pyqtgraph delivers a click only for a
            # button that never became a drag); a right click keeps falling
            # through to the ViewBox menu. Only ev.button() is read here.
            if ev.button() == Qt.LeftButton:
                self._select_spectrum(k)
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
        if self._comparison.signature(k):        # acquired/processed differently
            text += " ⚠"
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
        self._apply_highlight()       # fresh curves inherit the dim/spotlight state

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

    # ------------------------------------------------------------------ table
    # The results table under the grid and the two-way spotlight link between
    # its rows and the cells. Row identity is ALWAYS the spectrum index k,
    # stored under Qt.UserRole on the "#" item of each row (see _row_k): labels
    # may collide (two procs of one sample) and rows move when a header is
    # clicked, so neither the label nor the row position can name a spectrum.
    def _build_table(self):
        t = self.table = QTableWidget(0, 0)
        t.setEditTriggers(QAbstractItemView.NoEditTriggers)
        t.setSelectionBehavior(QAbstractItemView.SelectRows)
        t.setSelectionMode(QAbstractItemView.SingleSelection)
        t.setAlternatingRowColors(True)
        t.verticalHeader().setVisible(False)     # it renumbers by position after a sort
        t.horizontalHeader().setStretchLastSection(False)
        # setSortingEnabled(True) sorts by the header's current indicator every
        # time it is re-enabled after a refill (Qt's default: column 0
        # descending): start in series order; a header click then takes over
        # and persists across refills
        t.horizontalHeader().setSortIndicator(0, Qt.AscendingOrder)
        t.setMinimumHeight(80)
        t.itemSelectionChanged.connect(self._on_table_selection)
        self._fill_table(None)

    def _table_method(self, result) -> str:
        """The error method the table shows: the combo's pick once it has been
        computed, else the covariance stderr already on the recipes -- the same
        choice _refresh_err_status reports, because error_table labels an
        unstored method while silently falling back to Param.stderr."""
        m = self.errCombo.currentData()
        return m if m in (getattr(result, "error_detail", {}) or {}) else "covariance"

    def _fill_table(self, result):
        """(Re)build the table: the fixed columns (#, sample, proc, S/N, RMSD)
        always, plus one column per site x parameter once ``result`` exists.
        Screen and Save table… / Export CSV… come from the same batchfit rows,
        so they cannot disagree. Keeps the current spotlight."""
        from larmor import batchfit

        t = self.table
        t.blockSignals(True)
        t.setSortingEnabled(False)           # filling while sorting is on scrambles rows
        has_proc = any(d.get("proc") for d in self._data)
        # the comparability column exists only when there are Bruker
        # parameters to compare, so CSV series keep their headers
        has_params = self._comparison.level != "none"
        fixed = (["#", "sample"] + (["proc"] if has_proc else [])
                 + (["comparability"] if has_params else []) + ["S/N", "RMSD"])
        if result is None:
            cols, cells = [], {}
        else:
            rows = batchfit.error_table(result, self._table_method(result))
            cols, cells = batchfit.pivot_by_spectrum(rows, len(result.recipes))
        headers = fixed + [
            f"s{i} {label}\n" + ("population %" if pn == "population_pct"
                                 else PARAM_LABELS.get(pn, pn))
            for i, label, pn in cols]
        t.setRowCount(0)
        t.setColumnCount(len(headers))
        t.setHorizontalHeaderLabels(headers)
        t.setRowCount(len(self._data))
        for k, d in enumerate(self._data):
            items = []
            it = _NumItem(str(k + 1))
            it.setData(Qt.UserRole, k)                        # THE identity
            it.setData(Qt.UserRole + 1, float(k))
            it.setToolTip(str(d.get("path", ""))
                          + (f" · proc {d['proc']}" if d.get("proc") else ""))
            items.append(it)
            items.append(QTableWidgetItem(str(d.get("sample", ""))))
            if has_proc:
                items.append(QTableWidgetItem(str(d.get("proc", ""))))
            if has_params:
                items.append(self._comparability_item(k, d))
            snr = d.get("snr")
            snr = float(snr) if _finite(snr) else float("nan")
            it = _NumItem(f"{snr:.0f}" if np.isfinite(snr) else "")
            it.setData(Qt.UserRole + 1, snr)
            items.append(it)
            if result is None or k >= len(result.rmsd):
                items.append(_NumItem("—"))
            else:
                rmsd = float(result.rmsd[k])
                it = _NumItem(f"{rmsd:.4f}")
                it.setData(Qt.UserRole + 1, rmsd)
                if k in self._flag_reasons:                   # same gate as the cell label
                    it.setText(f"⚠ {rmsd:.4f}")
                    it.setForeground(QColor(_FLAG_CSS_COLOR))
                    it.setToolTip(self._flag_reasons[k])
                items.append(it)
            for col in cols:
                i, _label, pn = col
                row = cells.get((k, col))
                if row is None:
                    it = _NumItem("")
                    rec = result.recipes[k] if k < len(result.recipes) else None
                    if (rec is not None and i < len(rec.sites)
                            and batchfit.is_zeroed_out(
                                rec.sites[i].params.get("amplitude"))):
                        it.setToolTip("excluded for this spectrum")
                else:
                    value, stderr = row.get("value"), row.get("stderr")
                    txt = _num(value, ".4g")
                    if _finite(stderr):
                        txt += f" ± {float(stderr):.2g}"
                    side = row.get("at_bound")          # N5: ‡ finished at a bound
                    if side:
                        txt += " ‡"
                    it = _NumItem(txt)
                    if _finite(value):
                        it.setData(Qt.UserRole + 1, float(value))
                    it.setForeground(QColor(site_color(i)))
                    tip = f"{pn} · {row.get('error_method', '')}"
                    if _finite(row.get("sigma_pct")):
                        tip += f" · σ {float(row['sigma_pct']):.2g} %"
                    if side:
                        bound = row.get(side)
                        tip += (f" · ‡ finished at its {'lower' if side == 'min' else 'upper'} "
                                f"bound ({_num(bound, '.4g') if _finite(bound) else '?'}): "
                                "the shared value does not fit this spectrum — widen "
                                "the release % or fit the series sequentially")
                    it.setToolTip(tip)
                items.append(it)
            for c, it in enumerate(items):
                t.setItem(k, c, it)
        t.setSortingEnabled(True)                # re-sorts by the current indicator
        t.resizeColumnsToContents()
        if self._hl is not None:
            r = self._row_of(self._hl)
            if r is not None:
                t.selectRow(r)
        t.blockSignals(False)
        self._apply_highlight()

    def _comparability_item(self, k: int, d: dict) -> QTableWidgetItem:
        """'reprocessed' / the deviation signature in amber (red for bad) /
        '✓' -- a plain item, so sorting groups the three kinds together."""
        if d.get("proc_ops"):
            it = QTableWidgetItem("reprocessed")
            if self._reprocess:
                it.setToolTip(comparability.describe_template(
                    self._reprocess["template"], self._reprocess.get("phase", "own")))
            return it
        sig = self._comparison.signature(k)
        if not sig:
            return QTableWidgetItem("✓")
        it = QTableWidgetItem(sig)
        it.setForeground(QColor(_FLAG_CSS_COLOR if self._comparison.level_of(k) == "bad"
                                else _CHECK_CSS_COLOR))
        it.setToolTip("\n".join(self._comparison.details(k)))
        return it

    def _row_k(self, row: int):
        """Spectrum index of a table row, read from Qt.UserRole on its "#" item
        -- never the row position (rows move when sorted) or the sample label
        (two procs of one sample share it). None for an empty row."""
        it = self.table.item(row, 0)
        return None if it is None else it.data(Qt.UserRole)

    def _row_of(self, k: int):
        """The table row currently showing spectrum k (layout- and sort-agnostic)."""
        for r in range(self.table.rowCount()):
            if self._row_k(r) == k:
                return r
        return None

    def _on_table_selection(self):
        rows = self.table.selectionModel().selectedRows()
        self._highlight_cell(self._row_k(rows[0].row()) if rows else None)

    def _select_spectrum(self, k: int):
        """Cell -> table: select spectrum k's row, scroll it to the centre and
        hand the table focus, so ↑/↓ then step the spotlight along the series."""
        r = self._row_of(k)
        if r is not None:
            self.table.selectRow(r)          # fires _on_table_selection unless already selected
            it = self.table.item(r, 0)
            if it is not None:
                self.table.scrollToItem(it, QAbstractItemView.PositionAtCenter)
            if not self.table.isHidden():
                self.table.setFocus()
        self._highlight_cell(k)              # idempotent; also covers r is None

    def _highlight_cell(self, k):
        """Table -> cells: spotlight spectrum k (None clears) and bring its
        grid page forward."""
        self._hl = k
        self._apply_highlight()
        if k is not None and self.tabs.count():
            self.tabs.setCurrentIndex(k // PER_TAB)
        self._selection_status()

    def _apply_highlight(self):
        """Render self._hl onto the cells: the spotlighted cell gets a 2 px
        accent frame, an accent title and wider experiment/model pens; every
        other cell's curves (components included) fade to 35 % opacity --
        theme-agnostic, and no item is rebuilt. Every key is read with .get()
        and None-checked because tests stub cells with partial dicts."""
        t = theme.active()
        for j, cell in enumerate(self._cells):
            sel = (self._hl == j)
            op = 0.35 if (self._hl is not None and not sel) else 1.0
            exp = cell.get("exp")
            if exp is not None:
                exp.setOpacity(op)
                exp.setPen(pg.mkPen(t.experiment, width=2.0 if sel else 1.0))
            model = cell.get("model")
            if model is not None:
                model.setOpacity(op)
                model.setPen(pg.mkPen(t.model, width=2.2 if sel else 1.4))
            for it in cell.get("comp") or []:
                it.setOpacity(op)
            plot = cell.get("plot")
            if plot is not None:
                plot.getViewBox().setBorder(
                    pg.mkPen(t.accent, width=2) if sel else None)
            title = cell.get("title")
            if title is not None:
                # idle: the cell's own css (amber when it deviates from the
                # series majority); the accent spotlight still wins while selected
                title.setStyleSheet(
                    _TITLE_CSS + f" background:{t.accent}; color:{t.accent_text}; "
                    "padding:0 3px; border-radius:2px;" if sel
                    else cell.get("title_css", _TITLE_CSS))

    def _selection_status(self):
        """Status line for the spotlight -- the discoverability path for the
        keyboard and the flag reasons; clearing restores the fit summary."""
        k = self._hl
        if k is None:
            self.status.setText(self._result.summary if self._result is not None
                                else self._model_status())
            return
        d = self._data[k] if k < len(self._data) else {}
        parts = [f"spectrum {k + 1}/{len(self._data)} “{d.get('sample', '')}”"]
        if self._result is not None and k < len(self._result.rmsd):
            parts.append(f"RMSD {self._result.rmsd[k]:.4f}")
        if k in self._flag_reasons:
            parts.append("⚠ " + self._flag_reasons[k])
        sig = self._comparison.signature(k)
        if sig:                                   # differs from the series majority
            parts.append("⚠ " + sig)
        if _finite(d.get("snr")):
            parts.append(f"S/N {float(d['snr']):.0f}")
        self.status.setText(" · ".join(parts)
                            + " — ↑/↓ step through the series · Esc clears")

    def _toggle_table(self, on: bool):
        self.table.setVisible(on)

    def keyPressEvent(self, ev):
        # the first Esc clears the spotlight (intercepted only while one
        # exists); the second reaches QDialog and closes, as before
        if ev.key() == Qt.Key_Escape and self._hl is not None:
            self.table.clearSelection()      # -> _on_table_selection -> _highlight_cell(None)
            self._highlight_cell(None)       # also when the table held no selection
            ev.accept()
            return
        super().keyPressEvent(ev)

    # ------------------------------------------------------------------ fit
    def _entries(self):
        from larmor.recipe import Recipe

        out = []
        for k, d in enumerate(self._data):
            rec = Recipe.from_dict({
                "nucleus": d["nucleus"], "larmor_frequency_MHz": d["larmor"],
                "spin_rate_Hz": d["spin"], "sample": d["sample"],
                # the reprocess-from-fid chain (if any) then whatever baseline
                # correction is active -- the global "Fit baseline…" tool (see
                # BASELINE_OPS) or a per-spectrum manual 2-point pick --
                # recorded so the exported fit reproduces the EXACT spectrum
                # fitted here from the raw source later (the Plotting studio's
                # batch-grid and a reopened recipe replay this for "experiment")
                "processing": (list(d.get("proc_ops") or [])
                               + list(d.get("baseline_ops") or [])),
                "processing_from_raw": bool(d.get("proc_ops")),
                # a chain that starts in the time domain must point at the
                # EXPNO: loader.apply_processing checks bruker.is_expno(source)
                # before proc.from_bruker_fid, so the 1r path would make every
                # saved recipe of a reprocessed batch fail to replay
                "source_path": (d.get("expno") or d.get("path", "")),
                # where the display name came from (the Series table's folder /
                # title columns) and the replicate group it belongs to
                "provenance": {"sample_folder": d.get("folder", ""),
                               "title": d.get("title", ""),
                               "series_group": d.get("group", d["sample"])},
                # a member NOT reprocessed that differs from the series
                # majority carries the fact into its saved fit
                "notes": ([comparability.caveat_note(self._comparison, k)]
                          if (not d.get("proc_ops")
                              and comparability.caveat_note(self._comparison, k))
                          else []),
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
        self.btnBundle.setEnabled(False); self.btnSeries.setEnabled(False)
        self.btnErr.setEnabled(False); self.btnErrCsv.setEnabled(False)
        self.btnSeriesTable.setEnabled(False)          # _done / _err_done index cells by k
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
        self.btnSeriesTable.setEnabled(True)
        if mode == "cancel":                           # revert to the pre-fit view
            for k, cell in enumerate(self._cells):
                px, py = self._pre[k] if k < len(getattr(self, "_pre", [])) else (None, None)
                cell["model"].setData(px if px is not None else [],
                                      py if py is not None else [])
            self.status.setText("fit cancelled — reverted")
            return
        self._result = result
        self.btnSave.setEnabled(True); self.btnTable.setEnabled(True)
        self.btnBundle.setEnabled(True); self.btnSeries.setEnabled(True)
        self.btnErr.setEnabled(True); self.btnErrCsv.setEnabled(True)
        self._refresh_err_status()
        for k, pd in enumerate(result.per_dataset):
            if k < len(self._cells):
                self._cells[k]["model"].setData(np.asarray(pd["x"], float),
                                                np.asarray(pd["y_fit"], float))
                self._cells[k]["rmsd"].setText(f"RMSD {result.rmsd[k]:.4f}")
        self._refresh_components()
        flagged = self._flag_quality(result)
        self._fill_table(result)
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
        self._flag_reasons.clear()            # the table's RMSD column reads these
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
                self._flag_reasons[k] = "; ".join(reasons)
                lbl.setText(f"⚠ RMSD {r[k]:.4f}")
                lbl.setStyleSheet(
                    f"font-size:9px; color:{_FLAG_CSS_COLOR}; font-weight:600;")
                lbl.setToolTip(self._flag_reasons[k])
            else:
                lbl.setStyleSheet(
                    f"font-size:9px; color:{theme.active().text_dim};")
                lbl.setToolTip(f"S/N ≈ {snr:.0f}" if np.isfinite(snr) else "")
        return flagged

    def _failed(self, msg):
        self.prog.setRange(0, 100); self.prog.setValue(0)
        self.btnCancel.setEnabled(False); self.btnStop.setEnabled(False)
        self._update_fit_enabled()
        self.btnSeriesTable.setEnabled(True)
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
        from larmor import batchfit

        from larmor.desktop.paths import suggest_save_dir
        start = suggest_save_dir(self._src_paths[0] if self._src_paths else None)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save batch table",
            str(Path(start) / "batch_table.csv") if start else "batch_table.csv",
            "CSV (*.csv)")
        if not path:
            return
        # batchfit.SHARED_HEADER: model + source_path so the Plotting studio's
        # batch-grid figure can find each row's spectrum/fit straight from
        # this CSV, even without "Save individual fits…" too
        batchfit.write_shared_csv(self._result, path)
        msg = f"saved {path}"
        msg += self._write_series_sidecar(path)
        if self.chkAutoRecipes.isChecked():
            n = self._save_all_recipes_to(Path(path).parent)
            msg += f" · {n} individual fit(s) saved alongside it"
        self.status.setText(msg)

    def _export_bundle(self):
        """Publication bundle…: everything about this batch into one folder
        (larmor.io.bundle.write_bundle) -- the tables, one recipe copy and
        one _curves.csv per spectrum with the experiment EXACTLY as fitted
        (d["amp"], after any baseline; d["amp0"] as experiment_raw when a
        baseline was applied), manifest.csv and README.txt. The error table
        for the currently selected method is included only when it has
        been computed; nothing is computed here."""
        if self._result is None:
            return
        from larmor.io import bundle
        from larmor.desktop.paths import suggest_save_dir

        start = suggest_save_dir(self._src_paths[0] if self._src_paths else None)
        folder = QFileDialog.getExistingDirectory(
            self, "Publication bundle — choose an output folder", start)
        if not folder:
            return
        folder = Path(folder)
        if (folder / "manifest.csv").exists():
            ans = QMessageBox.question(
                self, "Replace bundle?",
                f"“{folder.name}” already holds a bundle (manifest.csv). "
                "Replace its files?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if ans != QMessageBox.Yes:
                return
        data = [(d["ppm"], d["amp"], self._window) for d in self._data]
        raw = [d["amp0"] if d.get("baseline_ops") else None for d in self._data]
        m = self.errCombo.currentData()
        method = m if m in (getattr(self._result, "error_detail", {}) or {}) else None
        try:
            res = bundle.write_bundle(
                self._result, data, folder, kind="batch",
                source_paths=[d["path"] for d in self._data], raw=raw,
                error_method=method)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Bundle failed", str(exc))
            return
        self.status.setText(res.summary + (
            f" · {len(res.warnings)} warning(s), see README.txt" if res.warnings else ""))

    def _save_template(self):
        import json
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QInputDialog
        name, ok = QInputDialog.getText(self, "Save batch setup", "Template name:")
        if not ok or not name.strip():
            return
        tpl = self._template_dict()
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
        self._apply_template(lib[choice])
        self.status.setText(f"loaded setup “{choice}” — Fit to apply")

    def _template_dict(self) -> dict:
        """The reusable setup -- release set, drift %, threshold, baseline
        kind and the view options: what Save setup… stores, and the core of
        a project's batch session (session_state)."""
        return {"release": [pn for pn, c in self._rel_checks.items() if c.isChecked()],
                "frac": self.frac.value(), "tol": self.tol.value(),
                "baseline": self._baseline_kind,
                "components": self.chkComp.isChecked(),
                "shared_scale": self.chkShared.isChecked()}

    def _apply_template(self, tpl: dict):
        for pn, c in self._rel_checks.items():
            c.setChecked(pn in tpl.get("release", []))
        self.frac.setValue(tpl.get("frac", 10))
        self.tol.setValue(tpl.get("tol", 0.0))
        self.chkComp.setChecked(tpl.get("components", False))
        self.chkShared.setChecked(tpl.get("shared_scale", False))

    def _series_plot(self):
        if self._result is None:
            return
        from larmor.desktop.series_plot import SeriesPlotDialog
        SeriesPlotDialog(self, self._result, series=self._series).exec()

    # ------------------------------------------------------------------ series table
    def _edit_series(self):
        """Series table…: names, replicate groups, order and joined columns of
        the series in one editor (larmor.desktop.series_table_dialog). OK
        renames in place; when rows moved, the grid is rebuilt in the new
        order with the fitted curves, exclusions, flags and spotlight following."""
        if any(c.get("bl_picking") for c in self._cells):
            self.status.setText("finish or cancel the 2-point baseline pick first "
                                "(right-click that cell) — the series table rebuilds "
                                "the grid")
            return
        from larmor.desktop.series_table_dialog import SeriesTableDialog

        busy = ((self._worker is not None and self._worker.isRunning())
                or (self._err_worker is not None and self._err_worker.isRunning()))
        dlg = SeriesTableDialog(self, self._series, locked=busy)
        if dlg.exec() != QDialog.Accepted:
            return
        self._apply_series(dlg.table(), dlg.perm())

    def _apply_series(self, table: SeriesTable, perm: list):
        """Adopt an edited series table. ``perm`` (new position -> current
        index) reorders the spectra first (_reorder); then every display name
        is written where labels live -- d["sample"], the fitted result's
        labels and its recipes' sample -- so the cells, the Results table, the
        CSV scopes, the auto-named recipes and the Series plot read the new
        names without a refit."""
        n = len(self._data)
        perm = list(perm)
        if perm != list(range(n)):
            self._reorder(perm)
        self._series = table
        for k, (d, row) in enumerate(zip(self._data, table.rows)):
            d["sample"] = row.display_name
            d["group"] = row.group
            if self._result is not None and k < len(self._result.labels):
                self._result.labels[k] = row.display_name
                self._result.recipes[k].sample = row.display_name
        self._comparison = self._recompare()          # its per-spectrum labels
        self._refresh_cell_marks()                    # titles, tooltips, exclusions
        for k in range(len(self._cells)):
            self._rebuild_cell_menu(k)                # the export title reads the name
        self._fill_table(self._result)
        self._selection_status()
        if self._hl is None:
            self.status.setText(self.status.text() + " · series table applied")

    def _reorder(self, perm: list):
        """Permute the spectra (``perm[new] = old``): _data / _src_paths and
        every index-keyed structure (_excluded, _flag_reasons, _hl), rebuild
        the grid, then re-pair a finished result by source_path
        (batchfit.align_result) and repaint it through _done -- the sequence
        apply_session_state already takes for a reopened session."""
        from larmor import batchfit

        self._data = [self._data[k] for k in perm]
        loaded = [str(d["path"]) for d in self._data]
        self._src_paths = loaded + [p for p in self._src_paths if p not in set(loaded)]
        where = {old: new for new, old in enumerate(perm)}
        self._excluded = {where[k]: v for k, v in self._excluded.items() if k in where}
        self._flag_reasons = {where[k]: v for k, v in self._flag_reasons.items()
                              if k in where}
        self._hl = where.get(self._hl) if self._hl is not None else None
        self._comparison = self._recompare()
        self.tabs.clear()
        self._cells = []
        self._build_grid()
        self._apply_scale()
        if self._result is not None:
            res, _dropped = batchfit.align_result(self._result,
                                                  [d["path"] for d in self._data])
            if len(res.per_dataset) != len(res.recipes):
                res.per_dataset = batchfit.result_curves(res, [d["ppm"] for d in self._data])
            self._done(res)
        else:
            self._fill_table(None)
        for k in self._excluded:
            self._rebuild_cell_menu(k)
        if self._hl is not None and self.tabs.count():
            self.tabs.setCurrentIndex(self._hl // PER_TAB)
        self._apply_highlight()

    def _write_series_sidecar(self, csv_path) -> str:
        """``<stem>_series.csv`` beside a table export: position, name, group,
        folder, title, source_path, every Series-table column and the RMSD --
        the wide companion of the long CSV, which never gains rows (the
        studio's rebuild iterates its scope rows). Returns the status note."""
        from larmor import series_table

        side = Path(csv_path).with_name(Path(csv_path).stem + "_series.csv")
        extra = {"RMSD": list(self._result.rmsd)} if self._result is not None else None
        try:
            series_table.write_series_csv(self._series, side, extra=extra)
        except OSError as exc:
            return f" · series table not written ({exc})"
        return f" · series table {side.name}"

    # ------------------------------------------------------------------ session
    # A project bundle (larmor.project) keeps a batch-fit session as a row of
    # the Workspaces dock: the setup, the spectra, the model, per-spectrum
    # baselines / exclusions, the error method and the fitted result. The
    # dialog is modal, so "saving the session" means recording its state when
    # it closes (app.py _remember_batch) and rebuilding it from that state
    # when the row is reopened (_open_batch_session -> apply_session_state).
    def session_state(self) -> dict:
        """Everything needed to reopen this dialog as it is. The result goes
        through batchfit.result_to_dict (recipes, RMSDs, error detail -- the
        curves are recomputed on reopen). Pass through project.json_safe
        before json.dumps."""
        from larmor import batchfit

        per = {}
        for k, d in enumerate(self._data):
            per[d["path"]] = {
                "baseline_ops": [dict(o) for o in (d.get("baseline_ops") or [])],
                # the reprocess-from-fid chain and its EXPNO, replayed BEFORE
                # the baseline and the stored result on reopen
                "proc_ops": [dict(o) for o in (d.get("proc_ops") or [])],
                "expno": d.get("expno", ""),
                "excluded": sorted(self._excluded.get(k, ()))}
        return {**self._template_dict(),
                "paths": list(self._src_paths),
                # the Series table: rows in paths order, each with its source_path
                "series": self._series.to_dict(),
                "reprocess": copy.deepcopy(self._reprocess),
                "model_sites": copy.deepcopy(self._model_sites),
                "window": list(self._window) if self._window else None,
                "recipe_tag": self._recipe_tag,
                "baseline_kind": self._baseline_kind,
                "per_spectrum": per,
                "error_method": self.errCombo.currentData(),
                "error_n": self.errN.value(),
                "auto_recipes": self.chkAutoRecipes.isChecked(),
                "table": self.chkTable.isChecked(),
                "result": (batchfit.result_to_dict(self._result)
                           if self._result is not None else None)}

    def apply_session_state(self, state: dict) -> list[str]:
        """Restore a session_state() dict onto a dialog built over
        ``state["paths"]`` and the model. Baseline ops are replayed through
        larmor.processing exactly as BASELINE_OPS promises; exclusions go
        through _toggle_exclude; the error method is set BEFORE the count
        (the combo handler clamps the count's range); the fitted result is
        re-paired with the spectra that actually loaded (batchfit.align_result),
        its curves recomputed and handed to _done, so cells, RMSD labels,
        buttons, the table and the error status repaint through the live-fit
        path. Returns the notes (spectra that could not be reloaded, whose
        results were dropped); they are also appended to the status line."""
        from larmor import batchfit

        notes: list[str] = []
        self._recipe_tag = state.get("recipe_tag") or ""
        self._apply_template(state)
        if state.get("series"):
            # the Series table back BEFORE the result block: names, groups and
            # columns pair by source_path, so _done paints the restored labels;
            # a spectrum that did not reload simply loses its row
            self._apply_series(
                SeriesTable.from_dict(state["series"]).aligned_to(
                    [d["path"] for d in self._data], fresh=self._series),
                list(range(len(self._data))))
        self._baseline_kind = (state.get("baseline_kind")
                               or state.get("baseline") or "None")
        per = state.get("per_spectrum") or {}
        order = None
        for k, d in enumerate(self._data):
            ps = per.get(d["path"]) or {}
            pops = [dict(o) for o in (ps.get("proc_ops") or [])]
            if pops:
                # the reprocess-from-fid chain first: amp0 becomes the
                # reprocessed spectrum, the baseline replay below starts from it
                p = d.get("params")
                if p is not None and p.has_fid:
                    try:
                        ppm, amp, _n = comparability.reprocess(p, pops)
                    except Exception as exc:  # noqa: BLE001
                        notes.append(f"could not reprocess {d['sample']} from its "
                                     f"fid ({exc}) — TopSpin spectrum kept")
                    else:
                        d["ppm"], d["amp"], d["amp0"] = ppm, amp, amp.copy()
                        d["proc_ops"], d["expno"] = pops, p.expno
                        d["snr"] = _snr(amp)
                        if k < len(self._cells):
                            self._cells[k]["exp"].setData(ppm, amp)
                            self._cells[k]["plot"].setXRange(float(ppm.min()),
                                                             float(ppm.max()))
                else:
                    notes.append(f"could not reprocess {d['sample']} from its fid "
                                 "— TopSpin spectrum kept")
            ops = [dict(o) for o in (ps.get("baseline_ops") or [])]
            if ops:
                from larmor import processing as proc

                s = proc.apply(proc.from_processed(d["ppm"], d["amp0"],
                                                   d["larmor"] or 1.0), ops)
                d["amp"] = np.asarray(s.y.real, float)
                d["baseline_ops"] = ops
                if k < len(self._cells):
                    self._cells[k]["exp"].setData(d["ppm"], d["amp"])
                for o in ops:
                    if "order" in o:
                        order = o["order"]
            excluded = [int(i) for i in (ps.get("excluded") or [])]
            for i in excluded:
                self._toggle_exclude(k, i, True)
            if excluded:
                self._rebuild_cell_menu(k)   # checkable entries follow the state
        if state.get("reprocess"):
            self._reprocess = copy.deepcopy(state["reprocess"])
            self._comparison = self._recompare()
            self._refresh_cell_marks()
            self.compBar.set_reprocessed(
                self._reprocess.get("template") or {},
                sum(1 for d in self._data if d.get("proc_ops")), len(self._data),
                self._reprocess.get("phase", "own"))
            self._fill_table(None)
        kind = self._baseline_kind
        if kind in (None, "", "None"):
            self.lblBaseline.setText("none")
        else:
            self.lblBaseline.setText(
                kind + (f" (order {order})"
                        if kind == "Polynomial" and order is not None else ""))
        idx = self.errCombo.findData(state.get("error_method"))
        if idx >= 0:
            self.errCombo.setCurrentIndex(idx)
        if state.get("error_n") is not None:
            self.errN.setValue(int(state["error_n"]))
        self.chkAutoRecipes.setChecked(bool(state.get("auto_recipes", True)))
        self.chkTable.setChecked(bool(state.get("table", True)))
        self._result = None
        res_d = state.get("result")
        if res_d:
            full = batchfit.result_from_dict(res_d)
            res, _dropped = batchfit.align_result(full, [d["path"] for d in self._data])
            if res.recipes and len(res.recipes) == len(self._data):
                res.per_dataset = batchfit.result_curves(
                    res, [d["ppm"] for d in self._data])
                self._done(res)
            elif res.recipes:
                notes.append("the fitted results could not be paired with the "
                             "loaded spectra — run Fit again")
        loaded = {d["path"] for d in self._data}
        missing = [p for p in (state.get("paths") or []) if p not in loaded]
        if missing:
            notes.insert(0, f"{len(missing)} spectrum/spectra could not be "
                            "reloaded — their results were dropped: "
                            + ", ".join(Path(p).name for p in missing))
        if self._result is None:
            self.status.setText(self._model_status())
        if notes:
            self.status.setText(self.status.text() + "  ⚠ " + "; ".join(notes))
        return notes

    def session_title(self) -> str:
        t = f"batch fit · {len(self._data)} spectra"
        return t + (f" · {self._recipe_tag}" if self._recipe_tag else "")

    def _rebuild_cell_menu(self, k: int):
        """Rebuild one cell's right-click menu from scratch so its checkable
        "Exclude component" entries reflect self._excluded (pyqtgraph's
        setMenuEnabled(True) builds a fresh default menu, which
        _attach_cell_menu then populates -- the same sequence _end_bg_pick
        relies on)."""
        if k >= len(self._cells) or self._cells[k]["bl_picking"]:
            return
        vb = self._cells[k]["plot"].getPlotItem().getViewBox()
        vb.setMenuEnabled(False)
        vb.setMenuEnabled(True)
        self._attach_cell_menu(k)

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
        if self._result is not None:          # the table shows what Export CSV… writes
            self._fill_table(self._result)

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
        self.btnBundle.setEnabled(False); self.btnSeries.setEnabled(False)
        self.btnSeriesTable.setEnabled(False)
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
        self.btnBundle.setEnabled(True); self.btnSeries.setEnabled(True)
        self.btnSeriesTable.setEnabled(True)

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
        self._fill_table(result)              # parameter cells switch to "value ± err"
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
        from larmor import batchfit

        method = self.errCombo.currentData()
        # batchfit.ERROR_HEADER: model + source_path so the Plotting studio's
        # batch-grid figure can find each row's spectrum/fit straight from
        # this CSV, even without "Save individual fits…" too
        batchfit.write_error_csv(self._result, path, method)
        msg = f"exported {Path(path).name} · {method} errors"
        msg += self._write_series_sidecar(path)
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


def _finite(v) -> bool:
    """True for a real, finite number (None / NaN / non-numeric -> False)."""
    try:
        return v is not None and bool(np.isfinite(float(v)))
    except (TypeError, ValueError):
        return False


def _num(v, fmt: str) -> str:
    return format(float(v), fmt) if _finite(v) else ""


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
