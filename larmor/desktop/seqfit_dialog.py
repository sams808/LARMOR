"""Sequential (forward–backward) fit dialog.

Fit a series one spectrum at a time with full per-spectrum control, carrying each
fitted spectrum's parameters forward to seed the next — then, optionally, let it
run itself: N passes sweeping forward and backward, with live RMSD evolution and
trajectory smoothing. Front-end for :mod:`larmor.seqfit`; 1D only.
"""
from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QInputDialog, QLabel, QMessageBox, QProgressBar, QPushButton,
    QScrollArea, QSpinBox, QVBoxLayout, QWidget,
)

from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS
from larmor.desktop.plot import site_color
from larmor.desktop.batchfit_dialog import _slug, _proc_number, _saved_tol, _save_tol
import datetime as _dt


class _SeqWorker(QThread):
    done = Signal(object, str)
    failed = Signal(str)
    step = Signal(int, int, float)          # (pass, spectrum index, rmsd)

    def __init__(self, entries, passes, start, propagate, smooth, tol):
        super().__init__()
        self.entries, self.passes, self.start = entries, passes, start
        self.propagate, self.smooth, self.tol = propagate, smooth, tol
        self._stop = False
        self._mode = ""

    def request_stop(self, mode: str):
        self._mode = mode
        self._stop = True

    def run(self):
        try:
            from larmor.seqfit import run_sequential
            res = run_sequential(
                self.entries, passes=self.passes, start=self.start,
                propagate=self.propagate, smooth=self.smooth, tol=self.tol,
                progress=lambda p, k, r: self.step.emit(p, k, r),
                should_stop=lambda: self._stop)
            self.done.emit(res, self._mode)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))


class SeqFitDialog(QDialog):
    def __init__(self, parent, paths, model_recipe: dict | None):
        super().__init__(parent)
        self.setWindowTitle("Sequential fit — forward / backward series sweep")
        self.resize(1140, 780)
        self._model_sites = ((model_recipe or {}).get("sites") or None)
        self._window = ((model_recipe or {}).get("fit_window_ppm") or None)
        self._recipe_tag = ""
        self._result = None
        self._worker = None
        self._cur = 0
        self._data = self._load(paths)
        self._recipes = [self._seed_recipe(d) for d in self._data]

        root = QHBoxLayout(self)

        # ---------------- left: controls ----------------
        left = QWidget(); lv = QVBoxLayout(left); left.setMaximumWidth(470)
        root.addWidget(left)

        nav = QHBoxLayout()
        self.btnPrev = QPushButton("◀ Prev"); self.btnPrev.clicked.connect(self._prev)
        self.btnNext = QPushButton("Next ▶"); self.btnNext.clicked.connect(self._next)
        self.lblNav = QLabel(""); self.lblNav.setStyleSheet("font-weight:600;")
        nav.addWidget(self.btnPrev); nav.addWidget(self.lblNav, 1); nav.addWidget(self.btnNext)
        lv.addLayout(nav)

        self.chkSeedMove = QCheckBox("seed from the spectrum I came from when I move")
        self.chkSeedMove.setChecked(True)
        lv.addWidget(self.chkSeedMove)

        from larmor.desktop.table import LinesTable
        self.table = LinesTable()
        self.table.edited.connect(self._resim_current)
        self.table.compute.connect(self._resim_current)
        self.table.fit.connect(self._fit_current)
        lv.addWidget(self.table, 1)

        frow = QHBoxLayout()
        self.btnFit = QPushButton("Fit current")
        self.btnFit.setStyleSheet("font-weight:600;")
        self.btnFit.clicked.connect(self._fit_current)
        self.btnFitNext = QPushButton("Fit → seed next ▶")
        self.btnFitNext.clicked.connect(self._fit_then_next)
        frow.addWidget(self.btnFit); frow.addWidget(self.btnFitNext)
        lv.addLayout(frow)

        lv.addWidget(QLabel("Carry these parameters to the next spectrum:"))
        self._propbox = QScrollArea(); self._propbox.setWidgetResizable(True)
        self._propbox.setMaximumHeight(52)
        holder = QWidget(); self._proplay = QHBoxLayout(holder)
        self._proplay.setContentsMargins(2, 2, 2, 2)
        self._propbox.setWidget(holder)
        lv.addWidget(self._propbox)
        self._prop_checks: dict = {}
        self._fill_prop_params()

        # auto panel
        auto = QWidget(); av = QVBoxLayout(auto)
        auto.setStyleSheet(f"QWidget{{border:1px solid {theme.active().border};"
                           "border-radius:4px;}}")
        arow = QHBoxLayout()
        arow.addWidget(QLabel("Passes:"))
        self.cbPasses = QComboBox(); self.cbPasses.addItems(["1", "2", "4", "8", "16"])
        self.cbPasses.setCurrentText("2")
        self.cbPasses.setToolTip("each pass sweeps one direction; 2 = forward then "
                                 "backward, 4 = F/B/F/B, …")
        arow.addWidget(self.cbPasses)
        arow.addWidget(QLabel("start:"))
        self.cbStart = QComboBox(); self.cbStart.addItems(["first", "last"])
        arow.addWidget(self.cbStart)
        arow.addWidget(QLabel("smooth:"))
        self.smooth = QSpinBox(); self.smooth.setRange(0, 15); self.smooth.setValue(0)
        self.smooth.setToolTip("moving-average window applied to each parameter's "
                               "trajectory between passes (0 = off)")
        arow.addWidget(self.smooth)
        av.addLayout(arow)
        brow = QHBoxLayout()
        self.btnAuto = QPushButton("Auto ⇄ forward–backward fit")
        self.btnAuto.setStyleSheet("font-weight:600;")
        self.btnAuto.clicked.connect(self._run_auto)
        self.btnCancel = QPushButton("Cancel"); self.btnCancel.setEnabled(False)
        self.btnCancel.clicked.connect(lambda: self._interrupt("cancel"))
        self.btnStop = QPushButton("Stop"); self.btnStop.setEnabled(False)
        self.btnStop.clicked.connect(lambda: self._interrupt("stop"))
        brow.addWidget(self.btnAuto, 1); brow.addWidget(self.btnCancel); brow.addWidget(self.btnStop)
        av.addLayout(brow)
        trow = QHBoxLayout()
        trow.addWidget(QLabel("Completion threshold:"))
        self.tol = QDoubleSpinBox(); self.tol.setRange(0.0, 50.0); self.tol.setDecimals(3)
        self.tol.setValue(_saved_tol()); self.tol.setSuffix(" % Δσ")
        trow.addWidget(self.tol); trow.addStretch(1)
        av.addLayout(trow)
        lv.addWidget(auto)

        self.prog = QProgressBar(); lv.addWidget(self.prog)
        self.status = QLabel("fit each spectrum, or run the auto sweep")
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"font-weight:600; color:{theme.active().accent};")
        lv.addWidget(self.status)

        # ---------------- right: plots ----------------
        right = QWidget(); rv = QVBoxLayout(right); root.addWidget(right, 1)
        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.setMinimumHeight(300)
        self.plot.getPlotItem().getViewBox().setMouseEnabled(True, True)
        self.plot.getPlotItem().invertX(True)          # NMR: ppm high → low
        self.plot.hideAxis("left")
        from larmor.desktop.plot_menu import attach_plot_menu
        attach_plot_menu(self.plot, title="spectrum", parent=self)
        self._curExp = self.plot.plot([], [], pen=pg.mkPen(theme.active().experiment, width=1))
        self._curModel = self.plot.plot([], [], pen=pg.mkPen(theme.active().model, width=1.6))
        self._curComp: list = []
        self.plotTitle = QLabel(""); self.plotTitle.setStyleSheet("font-weight:600;")
        rv.addWidget(self.plotTitle)
        rv.addWidget(self.plot, 3)

        evo = QHBoxLayout()
        self.rmsPlot = pg.PlotWidget(background=theme.active().plot_bg)
        self.rmsPlot.setLabel("bottom", "spectrum"); self.rmsPlot.setLabel("left", "RMSD")
        self.rmsPlot.showGrid(x=True, y=True, alpha=0.2)
        self._rmsCurve = self.rmsPlot.plot([], [], pen=pg.mkPen(theme.active().accent, width=2),
                                           symbol="o", symbolSize=6)
        evo.addWidget(self.rmsPlot, 1)
        trajbox = QWidget(); tv = QVBoxLayout(trajbox); tv.setContentsMargins(0, 0, 0, 0)
        prow = QHBoxLayout(); prow.addWidget(QLabel("trajectory:"))
        self.cbTraj = QComboBox(); self.cbTraj.currentIndexChanged.connect(self._draw_traj)
        prow.addWidget(self.cbTraj, 1)
        tv.addLayout(prow)
        self.trajPlot = pg.PlotWidget(background=theme.active().plot_bg)
        self.trajPlot.setLabel("bottom", "spectrum")
        self.trajPlot.showGrid(x=True, y=True, alpha=0.2)
        self._trajCurve = self.trajPlot.plot([], [], pen=pg.mkPen(theme.active().pivot, width=2),
                                             symbol="o", symbolSize=6)
        tv.addWidget(self.trajPlot)
        evo.addWidget(trajbox, 1)
        rv.addLayout(evo, 2)

        # ---------------- buttons ----------------
        bb = QDialogButtonBox(QDialogButtonBox.Close | QDialogButtonBox.Help)
        self.btnSave = bb.addButton("Save individual fits…", QDialogButtonBox.ActionRole)
        self.btnSave.clicked.connect(self._save_individual)
        self.btnSeries = bb.addButton("Series plot…", QDialogButtonBox.ActionRole)
        self.btnSeries.clicked.connect(self._series_plot)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.helpRequested.connect(self._help)
        rv.addWidget(bb)

        self._live_rmsd = [float("nan")] * len(self._data)
        self._fill_traj_options()
        self._show_current()

    # ------------------------------------------------------------------ data
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
                "sample": rec.get("sample") or Path(p).stem, "path": p,
                "proc": _proc_number(p)})
            if self._model_sites is None and rec.get("sites"):
                self._model_sites = rec["sites"]
        return data

    def _seed_recipe(self, d) -> dict:
        return {"nucleus": d["nucleus"], "larmor_frequency_MHz": d["larmor"],
                "spin_rate_Hz": d["spin"], "sample": d["sample"],
                "fit_window_ppm": self._window,
                "sites": copy.deepcopy(self._model_sites or [])}

    def _fill_prop_params(self):
        while self._proplay.count():
            it = self._proplay.takeAt(0)
            if it.widget():
                it.widget().deleteLater()
        self._prop_checks = {}
        if not self._model_sites:
            return
        from larmor.batchfit import all_but_amplitude
        from larmor.recipe import Recipe
        names = all_but_amplitude([Recipe.from_dict(
            {"sites": self._model_sites, "nucleus": "", "larmor_frequency_MHz": 0})])
        for pn in names:
            c = QCheckBox(PARAM_LABELS.get(pn, pn)); c.setChecked(True)
            self._prop_checks[pn] = c
            self._proplay.addWidget(c)
        self._proplay.addStretch(1)

    def _propagate(self):
        return tuple(pn for pn, c in self._prop_checks.items() if c.isChecked())

    def _fill_traj_options(self):
        self.cbTraj.blockSignals(True)
        self.cbTraj.clear()
        from larmor.recipe import Recipe
        rec = Recipe.from_dict(self._recipes[0])
        for i, s in enumerate(rec.sites):
            for pn in s.params:
                if pn == "gl":
                    continue
                self.cbTraj.addItem(f"s{i} {s.label or s.model}: {PARAM_LABELS.get(pn, pn)}",
                                    (i, pn))
        self.cbTraj.blockSignals(False)

    # ------------------------------------------------------------------ nav
    def _show_current(self):
        d = self._data[self._cur]
        self.lblNav.setText(f"spectrum {self._cur + 1} / {len(self._data)}  ·  "
                            f"{d['sample']}"
                            + (f"  (proc {d['proc']})" if d["proc"] else ""))
        self.plotTitle.setText(d["sample"])
        self.table.rebuild(self._recipes[self._cur], set())
        self.btnPrev.setEnabled(self._cur > 0)
        self.btnNext.setEnabled(self._cur < len(self._data) - 1)
        self._resim_current()

    def _prev(self):
        if self._cur > 0:
            src = self._cur
            self._cur -= 1
            if self.chkSeedMove.isChecked():
                self._seed_from(src)
            self._show_current()

    def _next(self):
        if self._cur < len(self._data) - 1:
            src = self._cur
            self._cur += 1
            if self.chkSeedMove.isChecked():
                self._seed_from(src)
            self._show_current()

    def _seed_from(self, src_idx):
        from larmor.recipe import Recipe
        from larmor import seqfit
        src = Recipe.from_dict(self._recipes[src_idx])
        dst = Recipe.from_dict(self._recipes[self._cur])
        seqfit.seed_from(dst, src, self._propagate())
        self._recipes[self._cur] = dst.to_dict()

    # ------------------------------------------------------------------ sim/fit
    def _resim_current(self):
        from larmor import engine
        from larmor.recipe import Recipe
        d = self._data[self._cur]
        self._curExp.setData(d["ppm"], d["amp"])
        self.plot.getPlotItem().getViewBox().setXRange(
            float(d["ppm"].min()), float(d["ppm"].max()), padding=0.02)
        for it in self._curComp:
            self.plot.removeItem(it)
        self._curComp = []
        try:
            rec = Recipe.from_dict(self._recipes[self._cur])
            x, total, per = engine.simulate(rec, exp_ppm=d["ppm"])
        except Exception:
            self._curModel.setData([], []); return
        self._curModel.setData(x, total)
        for i, ys in enumerate(per):
            it = self.plot.plot(x, np.asarray(ys, float),
                                pen=pg.mkPen(site_color(i), width=1, style=Qt.DashLine))
            self._curComp.append(it)

    def _fit_current(self):
        from larmor import fit as fitmod
        from larmor.recipe import Recipe
        d = self._data[self._cur]
        rec = Recipe.from_dict(self._recipes[self._cur])
        try:
            fitmod.fit(rec, d["ppm"], d["amp"], window_ppm=self._window,
                       tol=self.tol.value() or None)
        except Exception as exc:  # noqa: BLE001
            self.status.setText(f"fit failed: {exc}"); return
        _save_tol(self.tol.value())
        self._recipes[self._cur] = rec.to_dict()
        from larmor.seqfit import _rmsd
        self._live_rmsd[self._cur] = _rmsd(rec, d["ppm"], d["amp"], self._window)
        self.table.rebuild(self._recipes[self._cur], set())
        self._resim_current(); self._draw_rms(); self._draw_traj()
        self.status.setText(f"fitted spectrum {self._cur + 1} · "
                            f"RMSD {self._live_rmsd[self._cur]:.4g}")

    def _fit_then_next(self):
        self._fit_current()
        if self._cur < len(self._data) - 1:
            self._next()

    # ------------------------------------------------------------------ auto
    def _entries(self):
        from larmor.recipe import Recipe
        return [(Recipe.from_dict(self._recipes[k]), d["ppm"], d["amp"], self._window)
                for k, d in enumerate(self._data)]

    def _run_auto(self):
        if len(self._data) < 2 or not self._model_sites:
            self.status.setText("need a model and at least two spectra"); return
        self.btnAuto.setEnabled(False); self.btnFit.setEnabled(False)
        self.btnCancel.setEnabled(True); self.btnStop.setEnabled(True)
        self.prog.setRange(0, 0)
        self._pre_recipes = copy.deepcopy(self._recipes)
        _save_tol(self.tol.value())
        self._worker = _SeqWorker(
            self._entries(), int(self.cbPasses.currentText()),
            self.cbStart.currentText(), self._propagate(),
            self.smooth.value(), self.tol.value() or None)
        self._worker.step.connect(self._on_step)
        self._worker.done.connect(self._auto_done)
        self._worker.failed.connect(self._auto_failed)
        self.status.setText("sweeping…")
        self._worker.start()

    def _on_step(self, p, k, rmsd):
        self._live_rmsd[k] = rmsd
        self._draw_rms()
        finite = [v for v in self._live_rmsd if np.isfinite(v)]
        self.status.setText(f"pass {p + 1} · spectrum {k + 1} · RMSD {rmsd:.4g}"
                            + (f" · mean {np.mean(finite):.4g}" if finite else ""))

    def _interrupt(self, mode):
        if self._worker is not None and self._worker.isRunning():
            self._worker.request_stop(mode)
            self.status.setText("cancelling…" if mode == "cancel" else "stopping…")

    def _auto_done(self, result, mode):
        self.prog.setRange(0, 100); self.prog.setValue(0 if mode == "cancel" else 100)
        self.btnAuto.setEnabled(True); self.btnFit.setEnabled(True)
        self.btnCancel.setEnabled(False); self.btnStop.setEnabled(False)
        if mode == "cancel":
            self._recipes = self._pre_recipes
            self.status.setText("sweep cancelled — reverted"); self._show_current()
            return
        self._result = result
        self._recipes = [r.to_dict() for r in result.recipes]
        self._live_rmsd = list(result.rmsd)
        self._show_current(); self._draw_rms(); self._draw_traj()
        means = " → ".join(f"{h['mean']:.4g}" for h in result.history)
        self.status.setText(result.summary + f"   (pass means: {means})")

    def _auto_failed(self, msg):
        self.prog.setRange(0, 100); self.prog.setValue(0)
        self.btnAuto.setEnabled(True); self.btnFit.setEnabled(True)
        self.btnCancel.setEnabled(False); self.btnStop.setEnabled(False)
        self.status.setText(f"sweep failed: {msg}")

    # ------------------------------------------------------------------ evolution plots
    def _draw_rms(self):
        x = np.arange(1, len(self._data) + 1)
        y = np.array(self._live_rmsd, float)
        self._rmsCurve.setData(x[np.isfinite(y)], y[np.isfinite(y)])

    def _draw_traj(self):
        data = self.cbTraj.currentData()
        if not data:
            return
        i, pn = data
        from larmor.recipe import Recipe
        vals = []
        for d in self._recipes:
            try:
                vals.append(float(Recipe.from_dict(d).sites[i].params[pn].value))
            except (KeyError, IndexError):
                vals.append(np.nan)
        x = np.arange(1, len(vals) + 1)
        self._trajCurve.setData(x, np.array(vals, float))
        self.trajPlot.setLabel("left", self.cbTraj.currentText())

    # ------------------------------------------------------------------ save
    def _auto_name(self, rec, proc):
        parts = [rec.sample or "fit", rec.nucleus or ""]
        if self._recipe_tag:
            parts.append(self._recipe_tag)
        parts.append("seq")
        parts.append(_dt.datetime.now().strftime("%Y%m%d_%H%M"))
        return _slug("_".join(p for p in parts if p))

    def _save_individual(self):
        from larmor.recipe import Recipe
        folder = QFileDialog.getExistingDirectory(self, "Save individual fits")
        if not folder:
            return
        folder = Path(folder)
        mode = QMessageBox.question(
            self, "Naming", "Name the files automatically?\n\n"
            "Yes — auto (sample_nucleus_seq_YYYYMMDD_HHMM)\n"
            "No — type a name for each fit.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        recs = [Recipe.from_dict(d) for d in self._recipes]
        n = 0
        for k, rec in enumerate(recs):
            proc = self._data[k]["proc"] if k < len(self._data) else ""
            if mode == QMessageBox.Yes:
                name = self._auto_name(rec, proc)
            else:
                default = self._auto_name(rec, proc)
                text, ok = QInputDialog.getText(
                    self, "Fit name",
                    f"Name for “{rec.sample}”" + (f" (proc {proc})" if proc else "")
                    + f"   [{k + 1} of {len(recs)}]", text=default)
                if not ok:
                    break
                name = _slug(text) or default
            try:
                rec.save(folder / f"{name}.recipe.json"); n += 1
            except Exception:
                pass
        self.status.setText(f"saved {n} fit(s) to {folder}")

    def _series_plot(self):
        from larmor.recipe import Recipe
        from larmor.seqfit import SeqFitResult
        from larmor.desktop.series_plot import SeriesPlotDialog
        recs = [Recipe.from_dict(d) for d in self._recipes]
        res = self._result or SeqFitResult(
            recipes=recs, labels=[d["sample"] for d in self._data],
            rmsd=list(self._live_rmsd), per_dataset=[], history=[], passes=0,
            propagated=self._propagate())
        SeriesPlotDialog(self, res).exec()

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "multi-dataset", "Multi-dataset & co-fitting")
