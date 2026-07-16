"""LARMOR desktop: main window, laid out after dmfit.

Structure (mirroring dmfit's decomposition interface):
  - menu bar: File / Process / Decomposition / View / Models / ?
  - thin top toolbar: zoom shortcuts, undo/redo, add-line model buttons
  - narrow LEFT sidebar: quick view buttons (Full, Sites, Y auto, parts, pad)
  - central spectrum canvas with dmfit-style paddles (drag the square top
    handle = position+amplitude, drag the side circles = width)
  - BOTTOM dock: the fit-parameters spreadsheet (one row per line, pin
    checkbox beside every value) with Compute / Fit / chi2 footer,
    tabbed with the Report (quantification + fit report)
  - RIGHT dock: processing panel
  - status bar: live cursor x/y like dmfit
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QDockWidget, QFileDialog, QLabel, QMainWindow, QMessageBox,
    QPlainTextEdit, QPushButton, QTableWidget, QTableWidgetItem, QToolBar,
    QVBoxLayout, QWidget, QHBoxLayout,
)

from larmor import models as model_registry
from larmor.desktop.panels import ProcessingPanel
from larmor.desktop.plot import SpectrumView
from larmor.desktop.table import LinesTable
from larmor.recipe import Recipe

APP_STYLE = """
QMainWindow { background: #eef0ee; }
QMenuBar { background: #ffffff; border-bottom: 1px solid #d7dcd9; }
QToolBar { background: #ffffff; border-bottom: 1px solid #d7dcd9; spacing: 3px; padding: 2px; }
QToolBar#sidebar { border-right: 1px solid #d7dcd9; border-bottom: none; padding: 2px 1px; }
QToolButton { padding: 3px 7px; border-radius: 4px; }
QToolButton:hover { background: #e2f0f0; }
QToolButton:checked { background: #0e7c86; color: white; }
QDockWidget::title { background: #f6f8f6; padding: 3px 8px; border-bottom: 1px solid #d7dcd9; }
QTableWidget { background: white; gridline-color: #e4e8e5; }
QHeaderView::section { background: #f6f8f6; border: none; border-right: 1px solid #d7dcd9;
                       border-bottom: 1px solid #d7dcd9; padding: 2px 6px; }
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
    done = Signal(object, object, object)
    failed = Signal(str)

    def __init__(self, recipe_dict, exp_ppm):
        super().__init__()
        self.recipe_dict, self.exp_ppm = recipe_dict, exp_ppm

    def run(self):
        try:
            from larmor import engine
            from larmor import fit as fitmod

            recipe = Recipe.from_dict(self.recipe_dict)
            params = fitmod._make_params(recipe)
            fitmod._apply_params(recipe, params)
            x, total, per_site = engine.simulate(recipe, exp_ppm=self.exp_ppm)
            self.done.emit(x, total, per_site)
        except Exception as exc:
            self.failed.emit(str(exc))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LARMOR")
        self.resize(1440, 900)

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
        self._paddle_live = False   # true while a paddle is being dragged

        self.view = SpectrumView()
        self.setCentralWidget(self.view)
        self.view.add_requested.connect(self.add_site_at)
        self.view.paddle_moved.connect(self.on_paddle_moved)
        self.view.paddle_released.connect(self.on_paddle_released)
        self.view.cursor_moved.connect(
            lambda x, y: self.pos_label.setText(f"x: {x:.2f} ppm   y: {y:.4g}"))

        self._build_menus()
        self._build_toolbar()
        self._build_sidebar()
        self._build_bottom_docks()
        self._build_right_dock()

        self.pos_label = QLabel("")
        self.statusBar().addPermanentWidget(self.pos_label)
        self.statusBar().showMessage(
            "File > Open… (dmfit .fxmla / LARMOR recipe) or Open EXPNO…")

        self._sim_timer = QTimer(self)
        self._sim_timer.setSingleShot(True)
        self._sim_timer.setInterval(120)
        self._sim_timer.timeout.connect(self._simulate_now)

        self._restore_session()

    # ------------------------------------------------------------- menus
    def _build_menus(self):
        mb = self.menuBar()

        m_file = mb.addMenu("&File")
        self._add(m_file, "&Open…", self.open_file, "Ctrl+O")
        self._add(m_file, "Open &EXPNO…", self.open_expno, "Ctrl+Shift+O")
        m_file.addSeparator()
        self.actSave = self._add(m_file, "&Save recipe", self.save_recipe, "Ctrl+S")
        self._add(m_file, "Figure…", self.open_figure_dialog)
        m_file.addSeparator()
        self._add(m_file, "E&xit", self.close)

        m_proc = mb.addMenu("&Process")
        self._add(m_proc, "Show processing panel",
                  lambda: self.proc_dock.show())
        self._add(m_proc, "Autophase (ACME)",
                  lambda: self.apply_processing([{"op": "autophase"}], False))
        self._add(m_proc, "Baseline (order 3)",
                  lambda: self.apply_processing([{"op": "baseline", "order": 3}], False))
        self._add(m_proc, "Reset to original", self.reset_processing)

        m_dec = mb.addMenu("&Decomposition")
        self._add(m_dec, "&New fit (clear lines)", self.new_fit)
        m_dec.addSeparator()
        self.actFit = self._add(m_dec, "&Fit", self.run_fit, "F5")
        self._add(m_dec, "&Compute", self.request_simulation, "F9")
        m_dec.addSeparator()
        self.actQuant = self._add(m_dec, "&Report (quantify)", self.run_quantify, "F6")

        m_view = mb.addMenu("&View")
        self.actResid = self._add(m_view, "Residual", self._toggle_resid,
                                  checkable=True, checked=True)
        self.actComp = self._add(m_view, "Components", self._toggle_comp,
                                 checkable=True, checked=True)
        self.actPaddles = self._add(m_view, "Show paddles", self._toggle_paddles,
                                    checkable=True, checked=True)
        m_view.addSeparator()
        self._add(m_view, "Zoom to sites", self.zoom_sites)
        self._add(m_view, "Full spectrum", self.zoom_full)

        m_models = mb.addMenu("&Models")
        self._model_actions = {}
        for m in model_registry.describe_all():
            a = QAction(m["label"], self)
            a.setCheckable(True)
            a.setToolTip(m["description"])
            a.triggered.connect(
                lambda checked, name=m["name"]: self._set_add_mode(
                    name if checked else None))
            m_models.addAction(a)
            self._model_actions[m["name"]] = a

        m_help = mb.addMenu("&?")
        self._add(m_help, "About LARMOR", self._about)

    def _add(self, menu, text, slot, shortcut=None, checkable=False, checked=False):
        a = QAction(text, self)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        a.setCheckable(checkable)
        a.setChecked(checked)
        a.triggered.connect(slot)
        menu.addAction(a)
        return a

    def _about(self):
        QMessageBox.information(
            self, "LARMOR",
            "LARMOR — open successor to dmfit\n"
            "mrsimulator + lmfit + pyqtgraph\n"
            "github.com/sams808/LARMOR")

    # ------------------------------------------------------------- toolbar
    def _build_toolbar(self):
        tb = QToolBar("main")
        tb.setMovable(False)
        self.addToolBar(tb)
        self.actUndo = QAction("↩", self)
        self.actUndo.setShortcut(QKeySequence("Ctrl+Z"))
        self.actUndo.setToolTip("undo")
        self.actUndo.triggered.connect(self.undo)
        self.actRedo = QAction("↪", self)
        self.actRedo.setShortcut(QKeySequence("Ctrl+Y"))
        self.actRedo.setToolTip("redo")
        self.actRedo.triggered.connect(self.redo)
        tb.addAction(self.actUndo)
        tb.addAction(self.actRedo)
        tb.addSeparator()
        tb.addWidget(QLabel(" add line: "))
        for name, act in self._model_actions.items():
            tb.addAction(act)
        self._update_enabled()

    def _build_sidebar(self):
        sb = QToolBar("view")
        sb.setObjectName("sidebar")
        sb.setMovable(False)
        sb.setOrientation(Qt.Vertical)
        self.addToolBar(Qt.LeftToolBarArea, sb)
        for text, tip, slot in [
            ("Full", "full spectrum", self.zoom_full),
            ("Sites", "zoom to the fitted region", self.zoom_sites),
            ("Y a", "autoscale Y in the current X window", self.autoscale_y),
            ("pad", "toggle paddles", lambda: self.actPaddles.trigger()),
            ("parts", "toggle components", lambda: self.actComp.trigger()),
        ]:
            a = QAction(text, self)
            a.setToolTip(tip)
            a.triggered.connect(slot)
            sb.addAction(a)

    # ------------------------------------------------------------- docks
    def _build_bottom_docks(self):
        self.lines_dock = QDockWidget("Fit parameters", self)
        self.lines_dock.setFeatures(QDockWidget.DockWidgetMovable)
        self.lines_table = LinesTable()
        self.lines_table.edited.connect(self.on_params_changed)
        self.lines_table.constraint_edited.connect(self.on_structure_changed)
        self.lines_table.structure.connect(self.on_site_structure)
        self.lines_table.compute.connect(self.request_simulation)
        self.lines_table.fit.connect(self.run_fit)
        self.lines_dock.setWidget(self.lines_table)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.lines_dock)

        self.results_dock = QDockWidget("Report", self)
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
            ["line", "position (ppm)", "integral", "fraction (%)"])
        self.qtable.horizontalHeader().setStretchLastSection(True)
        self.qtable.verticalHeader().setVisible(False)
        v.addWidget(self.qtable)
        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        self.report.setStyleSheet("font-family: Consolas, monospace; font-size: 10px;")
        v.addWidget(self.report)
        self.results_dock.setWidget(w)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.results_dock)
        self.tabifyDockWidget(self.lines_dock, self.results_dock)
        self.lines_dock.raise_()

    def _build_right_dock(self):
        self.proc_dock = QDockWidget("Processing", self)
        self.proc_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                   QDockWidget.DockWidgetClosable)
        self.proc_panel = ProcessingPanel()
        self.proc_panel.apply_requested.connect(self.apply_processing)
        self.proc_panel.reset_requested.connect(self.reset_processing)
        self.proc_dock.setWidget(self.proc_panel)
        self.addDockWidget(Qt.RightDockWidgetArea, self.proc_dock)
        self.proc_dock.hide()

    def _update_enabled(self):
        loaded = self.recipe is not None
        for a in (self.actSave, self.actFit, self.actQuant):
            a.setEnabled(loaded)
        self.actUndo.setEnabled(bool(self.undo_stack))
        self.actRedo.setEnabled(bool(self.redo_stack))

    # ------------------------------------------------------------- view ops
    def zoom_full(self):
        self.view.getPlotItem().enableAutoRange()

    def zoom_sites(self):
        if not self.recipe or not self.recipe["sites"]:
            return
        pos = [s["params"]["isotropic_chemical_shift_ppm"]["value"]
               for s in self.recipe["sites"]
               if "isotropic_chemical_shift_ppm" in s["params"]]
        if pos:
            lo, hi = min(pos) - 120, max(pos) + 120
            self.view.setXRange(lo, hi, padding=0)
            self.autoscale_y()

    def autoscale_y(self):
        if not self.exp_ppm.size:
            return
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        sel = (self.exp_ppm >= min(x0, x1)) & (self.exp_ppm <= max(x0, x1))
        if sel.any():
            lo = float(self.exp_amp[sel].min())
            hi = float(self.exp_amp[sel].max())
            pad = 0.15 * (hi - lo or 1.0)
            self.view.setYRange(lo - pad, hi + pad, padding=0)

    def _toggle_resid(self, on):
        self.view.show_residual = on
        self.request_simulation()

    def _toggle_comp(self, on):
        self.view.show_components = on
        self.request_simulation()

    def _toggle_paddles(self, on):
        self.view.show_paddles(on)

    # ------------------------------------------------------------- loading
    def open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open spectrum or recipe", self._last_dir(),
            "NMR fits (*.fxmla *.fxml *.json);;All files (*)")
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
        self.zoom_full()
        if recipe["sites"]:
            self.zoom_sites()
        msg = meta + ("   ⚠ " + " • ".join(warnings) if warnings else "")
        self.statusBar().showMessage(msg)
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        self._update_enabled()
        if self.recipe["sites"]:
            self.request_simulation()
        self._persist_session()

    def new_fit(self):
        if self.recipe is None:
            return
        self.snapshot()
        self.recipe["sites"] = []
        self.on_structure_changed()

    # ------------------------------------------------------------- undo
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
        self.on_structure_changed()

    def redo(self):
        if not self.redo_stack:
            return
        self.undo_stack.append(json.dumps(self.recipe))
        self.recipe = json.loads(self.redo_stack.pop())
        self.on_structure_changed()

    # ------------------------------------------------------------- sites
    def _set_add_mode(self, name):
        for n, a in self._model_actions.items():
            a.setChecked(n == name)
        self.view.set_add_mode(name)
        if name:
            self.statusBar().showMessage(
                f"click on the spectrum to place a {name} line (Esc to cancel)")

    def keyPressEvent(self, ev):
        if ev.key() == Qt.Key_Escape:
            self._set_add_mode(None)
        super().keyPressEvent(ev)

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
        self.on_structure_changed()

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
            copy["label"] = (copy.get("label") or "line") + "-copy"
            for p in copy["params"].values():
                p["stderr"] = None
            self.recipe["sites"].append(copy)
        elif action == "visibility":
            (self.hidden.discard(idx) if idx in self.hidden
             else self.hidden.add(idx))
        self.on_structure_changed()

    def on_structure_changed(self):
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        self._update_enabled()
        self.request_simulation()

    def on_params_changed(self):
        self._persist_session()
        self._update_paddles()
        self.request_simulation()

    # ------------------------------------------------------------- paddles
    def _update_paddles(self):
        if not self.recipe:
            self.view.set_paddles([])
            return
        states = []
        for i, s in enumerate(self.recipe["sites"]):
            if i in self.hidden:
                continue
            p = s["params"]
            if "isotropic_chemical_shift_ppm" not in p:
                continue
            pos = p["isotropic_chemical_shift_ppm"]
            amp = p.get("amplitude", {"value": 1.0})
            fwhm = p.get("shift_fwhm_ppm", {"value": 1.0})
            movable = not (pos.get("expr") or amp.get("expr"))
            states.append((i, pos["value"], amp["value"], fwhm["value"], movable))
        self.view.set_paddles(states)
        self.view.show_paddles(self.actPaddles.isChecked())

    def on_paddle_moved(self, idx, pos, amp, fwhm):
        if not self.recipe or idx >= len(self.recipe["sites"]):
            return
        if not self._paddle_live:
            self.snapshot()
            self._paddle_live = True
        p = self.recipe["sites"][idx]["params"]
        p["isotropic_chemical_shift_ppm"]["value"] = pos
        if not p.get("amplitude", {}).get("expr"):
            p["amplitude"]["value"] = amp
        if "shift_fwhm_ppm" in p and not p["shift_fwhm_ppm"].get("expr"):
            p["shift_fwhm_ppm"]["value"] = fwhm
        self.request_simulation()

    def on_paddle_released(self, idx):
        self._paddle_live = False
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._persist_session()

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
    def run_fit(self):
        if not self.recipe or not self.recipe["sites"]:
            self.statusBar().showMessage("add at least one line first")
            return
        if self._fit_worker and self._fit_worker.isRunning():
            return
        self.snapshot()
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        hi, lo = max(x0, x1), min(x0, x1)
        self.statusBar().showMessage(f"fitting in {hi:.1f} … {lo:.1f} ppm …")
        self.lines_table.btnFit.setEnabled(False)
        self._fit_worker = FitWorker(json.loads(json.dumps(self.recipe)),
                                     self.exp_ppm, self.exp_amp, (hi, lo))
        self._fit_worker.done.connect(self._fit_done)
        self._fit_worker.failed.connect(self._fit_failed)
        self._fit_worker.start()

    def _fit_failed(self, msg: str):
        self.lines_table.btnFit.setEnabled(True)
        QMessageBox.warning(self, "Fit failed", msg)
        self.statusBar().showMessage("fit failed")

    def _fit_done(self, result):
        self.lines_table.btnFit.setEnabled(True)
        self.recipe = result.recipe.to_dict()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        labels = [s.get("label") or s["model"] for s in self.recipe["sites"]]
        self.view.set_model(result.x_ppm, result.y_fit, result.per_site,
                            labels, self.hidden, self.exp_ppm, self.exp_amp)
        self.lines_table.set_chi2(f"RMSD {result.rmsd:.4f}")
        bits = [f"RMSD {result.rmsd:.4f}"]
        if result.frozen_sites:
            bits.append("frozen: " + ", ".join(result.frozen_sites))
        if result.at_bounds:
            bits.append("⚠ at bounds: " + ", ".join(result.at_bounds))
        self.results_summary.setText("   ·   ".join(bits))
        self.report.setPlainText(result.report)
        self.statusBar().showMessage(
            "fit done" + ("  ⚠ parameters at bounds — check constraints"
                          if result.at_bounds else ""))
        self.run_quantify(show=False)
        self._persist_session()

    # ------------------------------------------------------------- quantify
    def run_quantify(self, show: bool = True):
        if not self.recipe or not self.recipe["sites"]:
            return
        from larmor.quantify import quantify

        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        try:
            q = quantify(Recipe.from_dict(self.recipe),
                         window_ppm=(max(x0, x1), min(x0, x1)))
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
            self.results_dock.raise_()

    def copy_csv(self):
        if not self._last_quant:
            return
        head = ("line,model,position_ppm,position_err,integral,"
                "integral_err,fraction_pct,fraction_err_pct")
        lines = [head]
        for r in self._last_quant["rows"]:
            lines.append(",".join(str(r.get(k, "") if r.get(k) is not None else "")
                                  for k in ("label", "model", "position_ppm",
                                            "position_err", "integral",
                                            "integral_err", "fraction_pct",
                                            "fraction_err_pct")))
        QApplication.clipboard().setText("\n".join(lines))
        self.statusBar().showMessage("report table copied as CSV")

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
                self.on_structure_changed()

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
                    self.on_structure_changed()
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
