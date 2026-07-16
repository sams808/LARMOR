"""LARMOR desktop: main window wiring the backend modules to Qt."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QDockWidget, QFileDialog, QInputDialog, QLabel, QMainWindow,
    QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, QTabWidget,
    QToolBar, QVBoxLayout, QWidget, QHBoxLayout, QPlainTextEdit,
)

from larmor import models as model_registry
from larmor.desktop.panels import ProcessingPanel, SitesPanel
from larmor.desktop.plot import SpectrumView
from larmor.recipe import Recipe

APP_STYLE = """
QMainWindow, QDockWidget { background: #eef0ee; }
QToolBar { background: #ffffff; border-bottom: 1px solid #d7dcd9; spacing: 4px; padding: 3px; }
QToolButton { padding: 4px 8px; border-radius: 5px; }
QToolButton:hover { background: #e2f0f0; }
QToolButton:checked { background: #0e7c86; color: white; }
QFrame#siteCard { background: white; border: 1px solid #d7dcd9; border-radius: 6px; }
QTableWidget { background: white; }
QStatusBar { background: #ffffff; border-top: 1px solid #d7dcd9; }
"""


class FitWorker(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, recipe_dict, ppm, amp, window):
        super().__init__()
        self.recipe_dict, self.ppm, self.amp, self.window = \
            recipe_dict, ppm, amp, window

    def run(self):
        try:
            from larmor import fit as fitmod

            recipe = Recipe.from_dict(self.recipe_dict)
            result = fitmod.fit(recipe, self.ppm, self.amp,
                                window_ppm=self.window)
            self.done.emit(result)
        except Exception as exc:
            self.failed.emit(str(exc))


class SimWorker(QThread):
    """Latest-wins simulation thread: UI stays fluid during kernel builds."""

    done = Signal(object, object, object)   # x, total, per_site
    failed = Signal(str)

    def __init__(self, recipe_dict, exp_ppm):
        super().__init__()
        self.recipe_dict, self.exp_ppm = recipe_dict, exp_ppm

    def run(self):
        try:
            from larmor import engine
            from larmor import fit as fitmod

            recipe = Recipe.from_dict(self.recipe_dict)
            params = fitmod._make_params(recipe)      # resolves links
            fitmod._apply_params(recipe, params)
            x, total, per_site = engine.simulate(recipe, exp_ppm=self.exp_ppm)
            self.done.emit(x, total, per_site)
        except Exception as exc:
            self.failed.emit(str(exc))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LARMOR")
        self.resize(1360, 840)

        self.source_path: str | None = None
        self.recipe: dict | None = None
        self.exp_ppm = np.array([])
        self.exp_amp = np.array([])
        self.hidden: set[int] = set()
        self.undo_stack: list[str] = []
        self.redo_stack: list[str] = []
        self._sim_worker: SimWorker | None = None
        self._sim_pending = False
        self._fit_worker: FitWorker | None = None
        self._last_quant = None

        self.view = SpectrumView()
        self.setCentralWidget(self.view)
        self.view.add_requested.connect(self.add_site_at)
        self.view.marker_moved.connect(self.marker_moved)

        self._build_toolbar()
        self._build_side_dock()
        self._build_results_dock()
        self.statusBar().showMessage("Open a dmfit .fxmla, a LARMOR recipe, "
                                     "or a Bruker EXPNO folder to begin")

        self._sim_timer = QTimer(self)
        self._sim_timer.setSingleShot(True)
        self._sim_timer.setInterval(150)
        self._sim_timer.timeout.connect(self._simulate_now)

        self._restore_session()

    # ------------------------------------------------------------- toolbar
    def _build_toolbar(self):
        tb = QToolBar("main")
        tb.setMovable(False)
        self.addToolBar(tb)

        def act(text, slot, shortcut=None, tip=None):
            a = QAction(text, self)
            if shortcut:
                a.setShortcut(QKeySequence(shortcut))
            if tip:
                a.setToolTip(tip)
            a.triggered.connect(slot)
            tb.addAction(a)
            return a

        act("Open file…", self.open_file, "Ctrl+O",
            "dmfit .fxmla or LARMOR .recipe.json")
        act("Open EXPNO…", self.open_expno, "Ctrl+Shift+O",
            "Bruker experiment folder (read-only)")
        self.actSave = act("Save recipe", self.save_recipe, "Ctrl+S")
        tb.addSeparator()
        self.actFit = act("Fit", self.run_fit, "F5", "run the fit (F5)")
        self.actQuant = act("Quantify", self.run_quantify, "F6",
                            "integrals / fractions table (F6)")
        act("Figure…", self.open_figure_dialog, None,
            "export a publication figure")
        tb.addSeparator()
        self.actUndo = act("↩", self.undo, "Ctrl+Z", "undo")
        self.actRedo = act("↪", self.redo, "Ctrl+Y", "redo")
        tb.addSeparator()

        tb.addWidget(QLabel(" Add site: "))
        self._model_actions = {}
        for m in model_registry.describe_all():
            a = QAction(m["label"], self)
            a.setCheckable(True)
            a.setToolTip(m["description"] + "  (then click on the spectrum)")
            a.triggered.connect(
                lambda checked, name=m["name"]: self._set_add_mode(
                    name if checked else None))
            tb.addAction(a)
            self._model_actions[m["name"]] = a

        tb.addSeparator()
        act("window ⇐ zoom", self.window_from_zoom, None,
            "use the current x-zoom as the fit window")
        self._update_enabled()

    def _set_add_mode(self, name):
        for n, a in self._model_actions.items():
            a.setChecked(n == name)
        self.view.set_add_mode(name)
        if name:
            self.statusBar().showMessage(
                f"click on the spectrum to place a {name} site (Esc to cancel)")

    def keyPressEvent(self, ev):
        if ev.key() == Qt.Key_Escape:
            self._set_add_mode(None)
        super().keyPressEvent(ev)

    # ------------------------------------------------------------- docks
    def _build_side_dock(self):
        dock = QDockWidget("Model", self)
        dock.setFeatures(QDockWidget.DockWidgetMovable)
        tabs = QTabWidget()
        self.sites_panel = SitesPanel()
        self.sites_panel.changed.connect(self.on_params_changed)
        self.sites_panel.structure.connect(self.on_site_structure)
        tabs.addTab(self.sites_panel, "Sites")
        self.proc_panel = ProcessingPanel()
        self.proc_panel.apply_requested.connect(self.apply_processing)
        self.proc_panel.reset_requested.connect(self.reset_processing)
        tabs.addTab(self.proc_panel, "Processing")
        dock.setWidget(tabs)
        dock.setMinimumWidth(380)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def _build_results_dock(self):
        self.results_dock = QDockWidget("Results", self)
        self.results_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                      QDockWidget.DockWidgetClosable)
        w = QWidget()
        v = QVBoxLayout(w)
        head = QHBoxLayout()
        self.results_summary = QLabel("")
        self.results_summary.setStyleSheet("font-weight: 600;")
        head.addWidget(self.results_summary, 1)
        btnCsv = QPushButton("Copy CSV")
        btnCsv.clicked.connect(self.copy_csv)
        head.addWidget(btnCsv)
        v.addLayout(head)
        self.qtable = QTableWidget(0, 4)
        self.qtable.setHorizontalHeaderLabels(
            ["site", "position (ppm)", "integral", "fraction (%)"])
        self.qtable.horizontalHeader().setStretchLastSection(True)
        self.qtable.verticalHeader().setVisible(False)
        v.addWidget(self.qtable)
        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        self.report.setMaximumHeight(160)
        self.report.setStyleSheet("font-family: Consolas, monospace; font-size: 10px;")
        v.addWidget(self.report)
        self.results_dock.setWidget(w)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.results_dock)
        self.results_dock.hide()

    def _update_enabled(self):
        loaded = self.recipe is not None
        for a in (self.actSave, self.actFit, self.actQuant):
            a.setEnabled(loaded)
        self.actUndo.setEnabled(bool(self.undo_stack))
        self.actRedo.setEnabled(bool(self.redo_stack))

    # ------------------------------------------------------------- loading
    def open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open spectrum or recipe", self._last_dir(),
            "NMR fits (*.fxmla *.json);;All files (*)")
        if path:
            self.load_source(path)

    def open_expno(self):
        path = QFileDialog.getExistingDirectory(
            self, "Open Bruker EXPNO folder (read-only)", self._last_dir())
        if path:
            self.load_source(path)

    def _last_dir(self) -> str:
        return QSettings("LARMOR", "app").value("lastDir", "")

    def load_source(self, path: str):
        self.statusBar().showMessage("loading…")
        QApplication.processEvents()
        try:
            ppm, amp, recipe, meta, warnings = _load_any(path)
        except Exception as exc:
            QMessageBox.warning(self, "Load failed", str(exc))
            self.statusBar().showMessage("load failed")
            return
        self.source_path = path
        self.exp_ppm, self.exp_amp = ppm, amp
        self.recipe = recipe
        self.hidden.clear()
        self.undo_stack.clear()
        self.redo_stack.clear()
        QSettings("LARMOR", "app").setValue("lastDir", str(Path(path).parent))
        self.setWindowTitle(f"LARMOR — {Path(path).name}")
        self.view.set_experiment(ppm, amp)
        self.view.getPlotItem().enableAutoRange()
        # open on the region of interest, not the full sideband manifold
        positions = [s["params"]["isotropic_chemical_shift_ppm"]["value"]
                     for s in recipe["sites"]
                     if "isotropic_chemical_shift_ppm" in s["params"]]
        if positions:
            lo, hi = min(positions) - 120, max(positions) + 120
            self.view.setXRange(lo, hi, padding=0)
            sel = (ppm >= lo) & (ppm <= hi)
            if sel.any():
                self.view.setYRange(float(amp[sel].min()) * 1.1,
                                    float(amp[sel].max()) * 1.15, padding=0)
        msg = meta + ("   ⚠ " + " • ".join(warnings) if warnings else "")
        self.statusBar().showMessage(msg)
        self.sites_panel.rebuild(self.recipe, self.hidden)
        self._update_markers()
        self._update_enabled()
        if self.recipe["sites"]:
            self.request_simulation()
        self._persist_session()

    # ------------------------------------------------------------- sites
    def snapshot(self):
        if self.recipe is None:
            return
        self.undo_stack.append(json.dumps(self.recipe))
        if len(self.undo_stack) > 60:
            self.undo_stack.pop(0)
        self.redo_stack.clear()
        self._update_enabled()
        self._persist_session()

    def undo(self):
        if not self.undo_stack:
            return
        self.redo_stack.append(json.dumps(self.recipe))
        self.recipe = json.loads(self.undo_stack.pop())
        self._after_structural_change()

    def redo(self):
        if not self.redo_stack:
            return
        self.undo_stack.append(json.dumps(self.recipe))
        self.recipe = json.loads(self.redo_stack.pop())
        self._after_structural_change()

    def _after_structural_change(self):
        self.sites_panel.rebuild(self.recipe, self.hidden)
        self._update_markers()
        self._update_enabled()
        self.request_simulation()

    def add_site_at(self, ppm: float, amp: float):
        name = next((n for n, a in self._model_actions.items()
                     if a.isChecked()), None)
        if not name or self.recipe is None:
            return
        self.snapshot()
        m = model_registry.get(name)
        params = {}
        for p in m.params:
            params[p.name] = {"value": p.default, "stderr": None,
                              "vary": p.vary, "min": p.min, "max": p.max,
                              "expr": None}
        params["isotropic_chemical_shift_ppm"]["value"] = ppm
        params["amplitude"]["value"] = amp or 1.0
        n = len(self.recipe["sites"])
        self.recipe["sites"].append(
            {"model": name, "label": f"{m.label.split(' ')[0]}-{n}",
             "params": params})
        self._set_add_mode(None)
        self._after_structural_change()

    def on_site_structure(self, idx: int, action: str):
        if self.recipe is None:
            return
        if action == "remove":
            self.snapshot()
            self.recipe["sites"].pop(idx)
            self.hidden.discard(idx)
        elif action == "duplicate":
            self.snapshot()
            copy = json.loads(json.dumps(self.recipe["sites"][idx]))
            copy["label"] = (copy.get("label") or "site") + "-copy"
            for p in copy["params"].values():
                p["stderr"] = None
            self.recipe["sites"].append(copy)
        elif action == "visibility":
            (self.hidden.discard(idx) if idx in self.hidden
             else self.hidden.add(idx))
        self._after_structural_change()

    def on_params_changed(self):
        self._persist_session()
        self._update_markers()
        self.request_simulation()

    def marker_moved(self, idx: int, ppm: float):
        if self.recipe and idx < len(self.recipe["sites"]):
            self.snapshot()
            self.recipe["sites"][idx]["params"][
                "isotropic_chemical_shift_ppm"]["value"] = ppm
            self.sites_panel.rebuild(self.recipe, self.hidden)
            self.request_simulation()

    def _update_markers(self):
        if not self.recipe:
            self.view.set_markers([])
            return
        marks = []
        for i, s in enumerate(self.recipe["sites"]):
            if i in self.hidden:
                continue
            p = s["params"].get("isotropic_chemical_shift_ppm")
            if p:
                marks.append((i, p["value"], not p.get("expr")))
        self.view.set_markers(marks)

    # ------------------------------------------------------------- simulate
    def request_simulation(self):
        self._sim_timer.start()

    def _simulate_now(self):
        if not self.recipe or not self.recipe["sites"]:
            self.view.set_model(None, None, None, None, self.hidden)
            return
        if self._sim_worker and self._sim_worker.isRunning():
            self._sim_pending = True
            return
        self._sim_worker = SimWorker(json.loads(json.dumps(self.recipe)),
                                     self.exp_ppm)
        self._sim_worker.done.connect(self._sim_done)
        self._sim_worker.failed.connect(
            lambda msg: self.statusBar().showMessage("simulate: " + msg))
        self._sim_worker.start()

    def _sim_done(self, x, total, per_site):
        labels = [s.get("label") or s["model"] for s in self.recipe["sites"]]
        self.view.set_model(x, total, per_site, labels, self.hidden,
                            self.exp_ppm, self.exp_amp)
        if self._sim_pending:
            self._sim_pending = False
            self._sim_timer.start()

    # ------------------------------------------------------------- fit
    def window_from_zoom(self):
        pass  # the window is read from the zoom at fit time; kept as a no-op
              # action so users discover the behavior from the tooltip

    def run_fit(self):
        if not self.recipe or not self.recipe["sites"]:
            self.statusBar().showMessage("add at least one site first")
            return
        if self._fit_worker and self._fit_worker.isRunning():
            return
        self.snapshot()
        hi, lo = self.view.current_xrange()
        self.statusBar().showMessage(
            f"fitting in window {hi:.1f} … {lo:.1f} ppm …")
        self.actFit.setEnabled(False)
        self._fit_worker = FitWorker(json.loads(json.dumps(self.recipe)),
                                     self.exp_ppm, self.exp_amp, (hi, lo))
        self._fit_worker.done.connect(self._fit_done)
        self._fit_worker.failed.connect(self._fit_failed)
        self._fit_worker.start()

    def _fit_failed(self, msg: str):
        self.actFit.setEnabled(True)
        QMessageBox.warning(self, "Fit failed", msg)
        self.statusBar().showMessage("fit failed")

    def _fit_done(self, result):
        self.actFit.setEnabled(True)
        self.recipe = result.recipe.to_dict()
        self.sites_panel.rebuild(self.recipe, self.hidden)
        self._update_markers()
        labels = [s.get("label") or s["model"] for s in self.recipe["sites"]]
        self.view.set_model(result.x_ppm, result.y_fit, result.per_site,
                            labels, self.hidden, self.exp_ppm, self.exp_amp)
        bits = [f"RMSD {result.rmsd:.4f}"]
        if result.frozen_sites:
            bits.append("frozen: " + ", ".join(result.frozen_sites))
        if result.at_bounds:
            bits.append("⚠ at bounds: " + ", ".join(result.at_bounds))
        self.results_summary.setText("   ·   ".join(bits))
        self.report.setPlainText(result.report)
        self.results_dock.show()
        self.statusBar().showMessage(
            "fit done" + ("  ⚠ some parameters at bounds — check constraints"
                          if result.at_bounds else ""))
        self.run_quantify(show=False)
        self._persist_session()

    # ------------------------------------------------------------- quantify
    def run_quantify(self, show: bool = True):
        if not self.recipe or not self.recipe["sites"]:
            return
        from larmor.quantify import quantify

        hi, lo = self.view.current_xrange()
        try:
            q = quantify(Recipe.from_dict(self.recipe), window_ppm=(hi, lo))
        except Exception as exc:
            self.statusBar().showMessage("quantify: " + str(exc))
            return
        self._last_quant = q
        rows = q["rows"]
        self.qtable.setRowCount(len(rows))
        for r, row in enumerate(rows):
            pos = f"{row['position_ppm']:.2f}"
            if row["position_err"]:
                pos += f" ± {row['position_err']:.2f}"
            frac = f"{row['fraction_pct']:.1f}"
            if row["fraction_err_pct"] is not None:
                frac += f" ± {row['fraction_err_pct']:.1f}"
            for c, text in enumerate([f"{row['label']}  ({row['model']})",
                                      pos, f"{row['integral']:.4g}", frac]):
                item = QTableWidgetItem(text)
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.qtable.setItem(r, c, item)
        self.qtable.resizeColumnsToContents()
        if show:
            self.results_dock.show()

    def copy_csv(self):
        if not self._last_quant:
            return
        head = ("site,model,position_ppm,position_err,integral,"
                "integral_err,fraction_pct,fraction_err_pct")
        lines = [head]
        for r in self._last_quant["rows"]:
            lines.append(",".join(str(r.get(k, "") if r.get(k) is not None else "")
                                  for k in ("label", "model", "position_ppm",
                                            "position_err", "integral",
                                            "integral_err", "fraction_pct",
                                            "fraction_err_pct")))
        QApplication.clipboard().setText("\n".join(lines))
        self.statusBar().showMessage("quantification table copied as CSV")

    # ------------------------------------------------------------- processing
    def apply_processing(self, ops: list, use_raw: bool):
        if not self.source_path:
            return
        from larmor import processing as proc
        from larmor.io import bruker

        self.statusBar().showMessage("processing…")
        QApplication.processEvents()
        try:
            if use_raw:
                if not bruker.is_expno(Path(self.source_path)):
                    raise ValueError("raw-fid processing needs a Bruker EXPNO")
                s = proc.from_bruker_fid(self.source_path)
            else:
                sfo1 = self.recipe.get("larmor_frequency_MHz", 0.0) if self.recipe else 0.0
                s = proc.from_processed(self.exp_ppm, self.exp_amp, sfo1)
            s = proc.apply(s, ops)
            if s.domain != "freq":
                raise ValueError("pipeline must end in the frequency domain")
        except Exception as exc:
            QMessageBox.warning(self, "Processing failed", str(exc))
            self.statusBar().showMessage("processing failed")
            return
        order = np.argsort(s.x_ppm)
        self.exp_ppm, self.exp_amp = np.asarray(s.x_ppm)[order], s.y.real[order]
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.request_simulation()
        self.statusBar().showMessage("processing applied")

    def reset_processing(self):
        if self.source_path:
            keep = json.loads(json.dumps(self.recipe)) if self.recipe else None
            self.load_source(self.source_path)
            if keep is not None:
                self.recipe = keep
                self._after_structural_change()

    # ------------------------------------------------------------- recipe io
    def save_recipe(self):
        if not self.recipe:
            return
        default = (self.recipe.get("sample") or "fit").strip()
        default = "".join(c if c.isalnum() or c in "-_" else "_"
                          for c in default)[:40] or "fit"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save recipe", str(Path(self._last_dir()) /
                                     f"{default}.recipe.json"),
            "LARMOR recipe (*.json)")
        if not path:
            return
        target = Path(path)
        for parent in [target.parent, *target.parent.parents]:
            if (parent / "acqus").exists() or (parent / "fid").exists() \
                    or (parent / "ser").exists():
                QMessageBox.warning(
                    self, "Refused",
                    f"{parent} is an instrument data folder — pick another "
                    "location. LARMOR never writes next to raw data.")
                return
        Recipe.from_dict(self.recipe).save(target)
        self.statusBar().showMessage(f"recipe saved: {target}")

    # ------------------------------------------------------------- figure
    def open_figure_dialog(self):
        from larmor.desktop.figure_dialog import FigureDialog

        FigureDialog(self, self.source_path, self.recipe).exec()

    # ------------------------------------------------------------- session
    def _persist_session(self):
        if not self.source_path:
            return
        s = QSettings("LARMOR", "app")
        s.setValue("session/source", self.source_path)
        s.setValue("session/recipe", json.dumps(self.recipe or {}))

    def _restore_session(self):
        s = QSettings("LARMOR", "app")
        src = s.value("session/source", "")
        if not src or not Path(src).exists():
            return
        try:
            self.load_source(src)
            saved = s.value("session/recipe", "")
            if saved:
                recipe = json.loads(saved)
                if recipe.get("sites"):
                    self.recipe = recipe
                    self._after_structural_change()
            self.statusBar().showMessage(f"session restored — {src}")
        except Exception:
            pass


def _load_any(path: str):
    """Load any supported source. Returns (ppm, amp, recipe_dict, meta, warnings)."""
    p = Path(path)
    if p.suffix.lower() == ".json":
        recipe = Recipe.load(p)
        if not recipe.source_path or not Path(recipe.source_path).exists():
            raise ValueError(f"recipe's source data not found: {recipe.source_path}")
        ppm, amp, _, meta, warnings = _load_any(recipe.source_path)
        return ppm, amp, recipe.to_dict(), f"recipe {p.name} | {meta}", warnings

    if p.suffix.lower() in (".fxmla", ".fxml"):
        from larmor.io import fxmla

        dm = fxmla.read(p)
        if dm.spectrum is None or dm.is_2d:
            raise ValueError("no 1D experimental data in this fxmla "
                             "(2D fitting arrives in Phase 2)")
        recipe, warnings = fxmla.to_recipe(dm)
        ppm, amp = dm.spectrum.ppm, dm.spectrum.amplitude
        order = np.argsort(ppm)
        return (ppm[order], amp[order], recipe.to_dict(),
                f"dmfit {dm.version} | {dm.comment}", warnings)

    from larmor.io import bruker

    if bruker.is_expno(p):
        exp = bruker.read_expno(p)
        if exp.processed is None:
            raise ValueError("EXPNO has no processed pdata/1 to display")
        ppm, amp = exp.processed_ppm, exp.processed.astype(float)
        order = np.argsort(ppm)
        recipe = Recipe(
            sample=exp.title.splitlines()[0] if exp.title else "",
            source_kind="bruker", source_path=str(p),
            nucleus=exp.nucleus, larmor_frequency_MHz=exp.sfo1_MHz,
            spin_rate_Hz=exp.masr_Hz or 0.0,
        )
        return (ppm[order], amp[order], recipe.to_dict(),
                exp.summary.splitlines()[1], exp.conflicts)

    raise ValueError(f"unrecognized source: {p} (expected .fxmla, "
                     ".recipe.json, or a Bruker EXPNO folder)")


def main() -> int:
    import pyqtgraph as pg

    pg.setConfigOptions(antialias=True)
    app = QApplication(sys.argv)
    app.setApplicationName("LARMOR")
    app.setStyleSheet(APP_STYLE)
    win = MainWindow()
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
