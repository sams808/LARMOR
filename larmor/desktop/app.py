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
QMainWindow { background: #f0f2f0; }
QMenuBar { background: #fbfcfb; color: #16202a; border-bottom: 1px solid #cfd6d1; }
QMenuBar::item { padding: 4px 10px; background: transparent; }
QMenuBar::item:selected { background: #dcebe9; border-radius: 4px; }
QMenu { background: #ffffff; color: #16202a; border: 1px solid #b9c1bc; }
QMenu::item { padding: 4px 26px 4px 18px; }
QMenu::item:selected { background: #dcebe9; }
QMenu::separator { height: 1px; background: #e0e5e1; margin: 4px 8px; }
QToolBar { background: #fbfcfb; border-bottom: 1px solid #cfd6d1; spacing: 3px; padding: 3px; }
QToolBar#sidebar { border-right: 1px solid #cfd6d1; border-bottom: none; padding: 3px 2px; }
QToolButton { padding: 4px 9px; border-radius: 4px; color: #16202a; border: 1px solid transparent; }
QToolButton:hover { background: #dcebe9; border-color: #b7cfcb; }
QToolButton:checked { background: #0e7c86; color: #ffffff; }
QDockWidget { color: #16202a; }
QDockWidget::title { background: #e7ebe8; padding: 4px 8px; border-top: 1px solid #cfd6d1; }
QTableWidget { background: #ffffff; color: #16202a; gridline-color: #e0e5e1;
               alternate-background-color: #f6f8f6; selection-background-color: #dcebe9;
               selection-color: #16202a; }
QHeaderView::section { background: #eef1ee; color: #37424a; font-weight: 600;
                       border: none; border-right: 1px solid #d3d9d4;
                       border-bottom: 1px solid #c5ccc6; padding: 3px 6px; }
QTableCornerButton::section { background: #eef1ee; border: none; }
QDoubleSpinBox, QSpinBox, QLineEdit { color: #16202a; background: transparent;
                                      selection-background-color: #bcdcd9; }
QPushButton { color: #16202a; background: #ffffff; border: 1px solid #aab4ad;
              border-radius: 4px; padding: 4px 14px; }
QPushButton:hover { background: #dcebe9; }
QPushButton:default { background: #0e7c86; color: #ffffff; border-color: #0e7c86; }
QCheckBox { color: #16202a; }
QLabel { color: #16202a; }
QTabBar::tab { background: #e7ebe8; color: #37424a; padding: 4px 14px;
               border: 1px solid #cfd6d1; border-bottom: none;
               border-top-left-radius: 4px; border-top-right-radius: 4px; }
QTabBar::tab:selected { background: #ffffff; color: #0a5a62; font-weight: 600; }
QStatusBar { background: #fbfcfb; color: #37424a; border-top: 1px solid #cfd6d1; }
QPlainTextEdit { background: #ffffff; color: #16202a; }
QScrollBar:vertical { background: #eef1ee; width: 12px; }
QScrollBar::handle:vertical { background: #c5ccc6; border-radius: 5px; min-height: 30px; }
QScrollBar:horizontal { background: #eef1ee; height: 12px; }
QScrollBar::handle:horizontal { background: #c5ccc6; border-radius: 5px; min-width: 30px; }
"""


def _light_palette():
    """A complete light palette so the app renders identically whatever the
    OS theme (Windows dark mode was bleeding white-on-white text through)."""
    from PySide6.QtGui import QColor, QPalette

    p = QPalette()
    c = QColor
    p.setColor(QPalette.Window, c("#f0f2f0"))
    p.setColor(QPalette.WindowText, c("#16202a"))
    p.setColor(QPalette.Base, c("#ffffff"))
    p.setColor(QPalette.AlternateBase, c("#f6f8f6"))
    p.setColor(QPalette.Text, c("#16202a"))
    p.setColor(QPalette.PlaceholderText, c("#93a0a8"))
    p.setColor(QPalette.Button, c("#fbfcfb"))
    p.setColor(QPalette.ButtonText, c("#16202a"))
    p.setColor(QPalette.ToolTipBase, c("#ffffff"))
    p.setColor(QPalette.ToolTipText, c("#16202a"))
    p.setColor(QPalette.Highlight, c("#0e7c86"))
    p.setColor(QPalette.HighlightedText, c("#ffffff"))
    p.setColor(QPalette.Link, c("#0a5a62"))
    for group in (QPalette.Disabled,):
        p.setColor(group, QPalette.WindowText, c("#9aa5ab"))
        p.setColor(group, QPalette.Text, c("#9aa5ab"))
        p.setColor(group, QPalette.ButtonText, c("#9aa5ab"))
    return p


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


class Fit2DWorker(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, recipe_dict, data2d, method):
        super().__init__()
        self.recipe_dict, self.data2d, self.method = \
            recipe_dict, data2d, method

    def run(self):
        try:
            from larmor import twod

            recipe = Recipe.from_dict(self.recipe_dict)
            result = twod.fit_2d(recipe, self.data2d, method=self.method)
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


class ClickableLabel(QLabel):
    """A QLabel that emits doubleClicked (used for the experiment strip)."""
    doubleClicked = Signal()

    def mouseDoubleClickEvent(self, ev):
        self.doubleClicked.emit()
        super().mouseDoubleClickEvent(ev)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LARMOR")
        self.resize(1440, 900)

        self.source_path: str | None = None
        self.recipe: dict | None = None
        self.exp_ppm = np.array([])
        self.exp_amp = np.array([])
        #: unprocessed workbench spectrum the pipeline is (re)applied from, so
        #: live processing reflects ABSOLUTE settings instead of compounding
        self._proc_base: tuple[np.ndarray, np.ndarray] | None = None
        self.hidden: set[int] = set()
        self.undo_stack: list[str] = []
        self.redo_stack: list[str] = []
        self._sim_worker: SimWorker | None = None
        self._sim_pending = False
        self._fit_worker: FitWorker | None = None
        self._last_quant = None
        self._paddle_live = False   # true while a paddle is being dragged
        self._last_model = None     # (x, total) of the latest simulation
        self._first_sim = False     # autoscale Y once the first model arrives

        # central area holds a 1D spectrum view AND a 2D contour view; the
        # loader switches between them so ANY dataset opens with a basic
        # display and the user then picks what to do with it
        from PySide6.QtWidgets import QStackedWidget

        from larmor.desktop.twod_view import Contour2DView

        self.view = SpectrumView()
        self.view2d = Contour2DView()
        self.central_stack = QStackedWidget()
        self.central_stack.addWidget(self.view)      # index 0: 1D
        self.central_stack.addWidget(self.view2d)    # index 1: 2D
        self.setCentralWidget(self.central_stack)

        self.view.add_requested.connect(self.add_site_at)
        self.view.paddle_moved.connect(self.on_paddle_moved)
        self.view.paddle_released.connect(self.on_paddle_released)
        self.view.file_dropped.connect(self.load_source)
        self.view.cursor_moved.connect(
            lambda x, y: self.pos_label.setText(f"x: {x:.2f} ppm   y: {y:.4g}"))
        self.view.calibrate_picked.connect(self.on_calibrate_picked)
        self.view.measure_changed.connect(self.on_measure_changed)
        self.view2d.slice_to_fit.connect(self._trace_to_workbench)

        self._build_menus()
        self._build_toolbar()
        self._build_sidebar()
        self._build_explorer_dock()
        self._build_datasets_dock()
        self._build_workspaces_dock()
        self._build_bottom_docks()
        self._build_right_dock()

        self.exp_label = ClickableLabel("")
        self.exp_label.setStyleSheet("color: #0a5a62; font-weight: 600;")
        self.exp_label.setToolTip("double-click to edit the experiment "
                                  "parameters (nucleus, Larmor, νrot)")
        self.exp_label.setCursor(Qt.PointingHandCursor)
        self.exp_label.doubleClicked.connect(self.edit_experiment)
        self.statusBar().addPermanentWidget(self.exp_label)
        self.pos_label = QLabel("")
        self.statusBar().addPermanentWidget(self.pos_label)
        self.statusBar().showMessage(
            "File > Open… (dmfit .fxmla / LARMOR recipe) or Open EXPNO…")

        self._sim_timer = QTimer(self)
        self._sim_timer.setSingleShot(True)
        self._sim_timer.setInterval(120)
        self._sim_timer.timeout.connect(self._simulate_now)
        self._data2d = None            # the 2D dataset currently on the map
        self._fit2d_worker = None
        self.view2d.add_requested.connect(self.add_site_2d)
        self.view2d.load_1d_for_projection.connect(self.load_projection_1d)
        # workspace manager (TopSpin-style: open a 2D / extract a trace -> a new
        # workspace you can switch between, close, or save)
        self.workspaces: list[dict] = []
        self.active_ws: int | None = None
        self._ws_mode = "auto"         # "auto" | "reuse" | "new"
        self._in_load_source = False

        # give the spectrum the majority of the height; keep the parameter
        # dock compact so it does not swallow half the window when nearly empty
        self.resizeDocks([self.lines_dock], [260], Qt.Vertical)
        self.lines_dock.setMinimumHeight(150)

        self._restore_session()

    # ------------------------------------------------------------- menus
    def _build_menus(self):
        mb = self.menuBar()

        m_file = mb.addMenu("&File")
        self._add(m_file, "&Open…  (spectrum / recipe / 1r / 2rr)",
                  self.open_file, "Ctrl+O")
        self._add(m_file, "Open &sample…  (list all its spectra)",
                  self.open_sample, "Ctrl+Shift+S")
        self._add(m_file, "Open &EXPNO / folder…", self.open_expno,
                  "Ctrl+Shift+O")
        self.m_recent = m_file.addMenu("Open &recent")
        self._rebuild_recent()
        self._add(m_file, "Open &FID…  (process before FT)", self.open_fid,
                  "Ctrl+F")
        m_file.addSeparator()
        self.actSave = self._add(m_file, "&Save recipe", self.save_recipe, "Ctrl+S")
        self._add(m_file, "Save fit &as…  (txt / csv / json / dmfit)",
                  self.save_fit_as, "Ctrl+Shift+E")
        self._add(m_file, "Save s&pectrum as…  (CSV, reopenable in LARMOR)",
                  self.save_spectrum)
        m_file.addSeparator()
        self._add(m_file, "&Copy plot to clipboard  (with all lines)",
                  self.copy_plot, "Ctrl+Shift+C")
        self._add(m_file, "Save plot &image…  (png / svg)", self.save_plot_image)
        self._add(m_file, "Figure…", self.open_figure_dialog)
        m_file.addSeparator()
        self._add(m_file, "E&xit", self.close)

        m_proc = mb.addMenu("&Process")
        self.actExp = self._add(m_proc, "&Experiment parameters… (νrot, B0, nucleus)",
                                self.edit_experiment)
        m_proc.addSeparator()
        self._add(m_proc, "Show processing panel",
                  lambda: self.proc_dock.show())
        self._add(m_proc, "Processing s&teps…  (remove a step)",
                  self.edit_processing_steps)
        self._add(m_proc, "Autophase (ACME)",
                  lambda: self.apply_processing([{"op": "autophase"}], False))
        self._add(m_proc, "Baseline auto (order 3)",
                  lambda: self.apply_processing([{"op": "baseline", "order": 3}], False))
        self._add(m_proc, "Subtract &averages  (offset from the edges)",
                  lambda: self.apply_processing([{"op": "subtract_avg"}], False))
        self._add(m_proc, "Reset to original", self.reset_processing)
        m_proc.addSeparator()
        self._add(m_proc, "&Calibrate axis…  (click a peak, set its ppm)",
                  self.start_calibrate)
        self.actMeasure = self._add(m_proc, "&Measure Δ (ppm / Hz)",
                                    self.toggle_measure, checkable=True)
        m_proc.addSeparator()
        self._add(m_proc, "Subtract a spectrum (&background)…",
                  self.open_subtract)

        m_dec = mb.addMenu("&Decomposition")
        self._add(m_dec, "&New fit (clear lines)", self.new_fit)
        self._add(m_dec, "Add background &spectrum…  (fit another spectrum)",
                  self.add_background_spectrum)
        self._add(m_dec, "Add a line at every &peak…  (auto peak-pick)",
                  self.autopick_lines)
        self._add(m_dec, "Predict at another &field…  (what at X T?)",
                  self.predict_at_field)
        self._add(m_dec, "Add f&unction line…  (y = f(x; a,b,c,d))",
                  self.add_function_line)
        m_dec.addSeparator()
        self.actFit = self._add(m_dec, "&Fit", self.run_fit, "F5")
        self.actAuto = self._add(m_dec, "&Auto Fit (multi-start)…",
                                 self.run_auto_fit)
        self.actErrors = self._add(m_dec, "&Errors Analysis (χ² profile)…",
                                   self.run_errors_analysis)
        self._add(m_dec, "Co-&fit datasets…  (shared model, 1D + MQMAS)",
                  self.open_cofit)
        self._add(m_dec, "&Compute", self.request_simulation, "F9")
        self._add(m_dec, "Computing &parameters…  (kernel resolution)",
                  self.edit_computing_params)
        m_dec.addSeparator()
        self._add(m_dec, "Add fit &zone", self.add_zone)
        self._add(m_dec, "Clear zones", self.clear_zones)
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
        self._add(m_view, "&Back to 2D map", self.back_to_2d, "Ctrl+2")
        self._add(m_view, "Zoom to sites", self.zoom_sites)
        self._add(m_view, "Full spectrum", self.zoom_full)

        m_models = mb.addMenu("&Models")
        self._model_actions = {}
        for m in model_registry.describe_all():
            if m["name"] == "spectrum":
                continue          # added via Decomposition ▸ Add background…
            a = QAction(m["label"], self)
            a.setCheckable(True)
            a.setToolTip(m["description"])
            a.triggered.connect(
                lambda checked, name=m["name"]: self._set_add_mode(
                    name if checked else None))
            m_models.addAction(a)
            self._model_actions[m["name"]] = a

        m_tools = mb.addMenu("&Tools")
        self._add(m_tools, "&Integrals && measurements…  (integral, %, FWHM, CoM)",
                  self.open_integrals)
        m_tools.addSeparator()
        self._add(m_tools, "Relaxation / series (T1, T2)…", self.open_satrec)
        self._add(m_tools, "Per-site relaxation…  (uses the current fit)",
                  self.open_per_site_relaxation)
        self._add(m_tools, "QCPMG (echo train → spectrum)…", self.open_qcpmg)
        self._add(m_tools, "Variable temperature (Arrhenius / VFT)…", self.open_vt)
        self._add(m_tools, "REDOR (dipolar coupling)…", self.open_redor)
        self._add(m_tools, "Import DFT tensors (.magres)…", self.open_magres)
        m_tools.addSeparator()
        self._add(m_tools, "2D MQMAS viewer/fit…", self.open_twod)
        m_tools.addSeparator()
        self._add(m_tools, "Multi-dataset fit (CLI): larmor multifit a.json b.json",
                  lambda: None).setEnabled(False)

        m_util = mb.addMenu("&Utilities")
        self._add(m_util, "&NMR table…  (Larmor frequencies)", self.open_nmr_table)
        self._add(m_util, "&Conversion tools…  (shift / Cq / dipolar)",
                  self.open_convert)

        m_help = mb.addMenu("&?")
        m_man = m_help.addMenu("User &manuals")
        for name, title in (("spectra-1d", "1D spectra — processing & fitting"),
                            ("mqmas", "MQMAS (2D)"),
                            ("correlation-hmqc", "HMQC & correlation"),
                            ("relaxation", "Relaxation (T1/T2)"),
                            ("qcpmg", "QCPMG")):
            self._add(m_man, title,
                      lambda _=False, n=name, t=title: self._open_manual(n, t))
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

    def _rebuild_recent(self):
        self.m_recent.clear()
        paths = QSettings("LARMOR", "app").value("recent", []) or []
        if isinstance(paths, str):
            paths = [paths]
        for p in paths[:12]:
            act = self.m_recent.addAction(f"{Path(p).name}   —   {Path(p).parent}")
            act.triggered.connect(lambda _=False, pp=p: self.load_source(pp))
        if not paths:
            a = self.m_recent.addAction("(none yet)"); a.setEnabled(False)

    def _add_recent(self, path: str):
        s = QSettings("LARMOR", "app")
        paths = s.value("recent", []) or []
        if isinstance(paths, str):
            paths = [paths]
        paths = [p for p in paths if p != path]
        paths.insert(0, path)
        s.setValue("recent", paths[:12])
        self._rebuild_recent()

    def _active_plot_widget(self):
        return (self.view2d.glw if self.central_stack.currentWidget() is self.view2d
                else self.view)

    def copy_plot(self):
        """Copy the current plot (spectrum + model + all component lines, as
        shown) to the clipboard as an image."""
        pix = self._active_plot_widget().grab()
        QApplication.clipboard().setPixmap(pix)
        self.statusBar().showMessage("plot copied to clipboard")

    def save_plot_image(self):
        w = self._active_plot_widget()
        path, _ = QFileDialog.getSaveFileName(
            self, "Save plot image", str(Path(self._last_dir()) / "plot.png"),
            "PNG image (*.png);;SVG vector (*.svg)")
        if not path:
            return
        if path.lower().endswith(".svg"):
            try:
                from pyqtgraph.exporters import SVGExporter

                scene = (self.view.getPlotItem().scene()
                         if w is self.view else self.view2d.p_main.scene())
                SVGExporter(scene).export(path)
            except Exception as exc:
                QMessageBox.warning(self, "Save image", f"SVG export failed: {exc}")
                return
        else:
            w.grab().save(path)
        self.statusBar().showMessage(f"plot saved to {Path(path).name}")

    def open_integrals(self):
        from larmor.desktop.integrate_dialog import IntegralsDialog

        if not self.exp_ppm.size:
            self.statusBar().showMessage("open a 1D spectrum first")
            return
        IntegralsDialog(self, self.exp_ppm, self.exp_amp).exec()

    def open_nmr_table(self):
        from larmor.desktop.utilities import NmrTableDialog

        h1 = 400.0
        if self.recipe and self.recipe.get("larmor_frequency_MHz") and self.recipe.get("nucleus"):
            try:
                from larmor import nuclei as N

                iso = next(i for i in N.all_isotopes()
                           if i.symbol == self.recipe["nucleus"])
                # back out the magnet's ¹H frequency from this nucleus
                h1 = self.recipe["larmor_frequency_MHz"] * N.GAMMA_1H / abs(iso.gamma_MHz_T)
            except Exception:
                pass
        NmrTableDialog(self, h1).exec()

    def open_convert(self):
        from larmor.desktop.utilities import ConvertDialog

        sfo = (self.recipe.get("larmor_frequency_MHz", 100.0)
               if self.recipe else 100.0) or 100.0
        ConvertDialog(self, sfo).exec()

    def edit_processing_steps(self):
        from larmor.desktop.dialogs import ProcessingStepsDialog

        ops = (self.recipe.get("processing") if self.recipe else None) or []
        if not ops:
            self.statusBar().showMessage("no processing steps applied yet")
            return
        dlg = ProcessingStepsDialog(self, ops)
        if dlg.exec():
            self.apply_processing(dlg.result_ops(),
                                  bool(self.recipe.get("processing_from_raw")))

    def edit_computing_params(self):
        from larmor.desktop.dialogs import ComputingParamsDialog

        if ComputingParamsDialog(self).exec():
            self.statusBar().showMessage(
                "computing parameters updated — kernels rebuild on the next fit")
            self.request_simulation()

    def _open_manual(self, name: str, title: str):
        from larmor.desktop.help_dialog import show_help

        show_help(self, name, title)

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
        tb.setIconSize(tb.iconSize())
        self.addToolBar(tb)
        self.actUndo = QAction("↩  Undo", self)
        self.actUndo.setShortcut(QKeySequence("Ctrl+Z"))
        self.actUndo.setToolTip("undo (Ctrl+Z)")
        self.actUndo.triggered.connect(self.undo)
        self.actRedo = QAction("↪  Redo", self)
        self.actRedo.setShortcut(QKeySequence("Ctrl+Y"))
        self.actRedo.setToolTip("redo (Ctrl+Y)")
        self.actRedo.triggered.connect(self.redo)
        tb.addAction(self.actUndo)
        tb.addAction(self.actRedo)
        tb.addSeparator()
        lab = QLabel("  Add line ")
        lab.setStyleSheet("color:#5a6871; font-weight:600;")
        tb.addWidget(lab)
        for name, act in self._model_actions.items():
            tb.addAction(act)
        self._update_enabled()

    def _build_sidebar(self):
        sb = QToolBar("view")
        sb.setObjectName("sidebar")
        sb.setMovable(False)
        sb.setOrientation(Qt.Vertical)
        sb.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self.addToolBar(Qt.LeftToolBarArea, sb)
        for text, tip, slot in [
            ("Full", "show the full spectrum", self.zoom_full),
            ("Sites", "zoom to the fitted region", self.zoom_sites),
            ("Y-fit", "autoscale Y in the current X window", self.autoscale_y),
            ("Paddles", "toggle the on-spectrum paddles", lambda: self.actPaddles.trigger()),
            ("Parts", "toggle the component curves", lambda: self.actComp.trigger()),
            ("Resid.", "toggle the residual", lambda: self.actResid.trigger()),
            ("2D map", "back to the 2D contour map (Ctrl+2)", self.back_to_2d),
        ]:
            a = QAction(text, self)
            a.setToolTip(tip)
            a.triggered.connect(slot)
            sb.addAction(a)

    def _build_explorer_dock(self):
        from larmor.desktop.explorer import ExplorerPanel

        self.explorer_dock = QDockWidget("Explorer", self)
        self.explorer_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                       QDockWidget.DockWidgetClosable)
        self.explorer = ExplorerPanel()
        self._proj_pick_axis = None      # HMQC: awaiting an Explorer pick
        self.explorer.open_requested.connect(self._explorer_open)
        self.explorer_dock.setWidget(self.explorer)
        self.explorer_dock.setMinimumWidth(230)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.explorer_dock)

    def _build_workspaces_dock(self):
        from larmor.desktop.workspaces import WorkspacePanel

        self.ws_dock = QDockWidget("Workspaces", self)
        self.ws_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                 QDockWidget.DockWidgetClosable)
        self.ws_panel = WorkspacePanel()
        self.ws_panel.switch.connect(self.switch_workspace)
        self.ws_panel.close.connect(self.close_workspace)
        self.ws_panel.save.connect(self.save_workspace)
        self.ws_dock.setWidget(self.ws_panel)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.ws_dock)
        self.tabifyDockWidget(self.explorer_dock, self.ws_dock)
        self.explorer_dock.raise_()

    def _build_datasets_dock(self):
        from larmor.desktop.datasets import DatasetsPanel

        self._overlays: list[dict] = []
        self.datasets_dock = QDockWidget("Datasets", self)
        self.datasets_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                       QDockWidget.DockWidgetClosable)
        self.datasets_panel = DatasetsPanel()
        self.datasets_panel.add_requested.connect(self.add_overlay_dialog)
        self.datasets_panel.make_active.connect(self.overlay_make_active)
        self.datasets_panel.remove.connect(self.overlay_remove)
        self.datasets_panel.visibility_changed.connect(self.overlay_visibility)
        self.datasets_panel.offset_changed.connect(lambda _: self._refresh_overlays())
        self.datasets_dock.setWidget(self.datasets_panel)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.datasets_dock)
        self.tabifyDockWidget(self.explorer_dock, self.datasets_dock)
        self.explorer_dock.raise_()

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
        self.proc_panel.baseline_mode.connect(self._baseline_mode)
        self.proc_panel.baseline_apply.connect(self.apply_manual_baseline)
        self.proc_panel.baseline_clear.connect(self.view.clear_baseline)
        self.proc_dock.setWidget(self.proc_panel)
        self.addDockWidget(Qt.RightDockWidgetArea, self.proc_dock)
        self.proc_dock.hide()

    def _update_enabled(self):
        loaded = self.recipe is not None
        for a in (self.actSave, self.actFit, self.actQuant):
            a.setEnabled(loaded)
        self.actUndo.setEnabled(bool(self.undo_stack))
        self.actRedo.setEnabled(bool(self.redo_stack))

    # ------------------------------------------------------------- experiment
    def _update_exp_label(self):
        if not self.recipe:
            self.exp_label.setText("")
            return
        nu = self.recipe.get("spin_rate_Hz", 0.0) or 0.0
        mas = f"νrot {nu:.0f} Hz" if nu else "static"
        sr = self.recipe.get("sr_hz", 0.0) or 0.0
        sr_txt = f" · SR {sr:.0f} Hz" if sr else ""
        self.exp_label.setText(
            f"{self.recipe.get('nucleus', '?')} · "
            f"{self.recipe.get('larmor_frequency_MHz', 0):.3f} MHz · {mas}{sr_txt}")

    def edit_experiment(self):
        if self.recipe is None:
            return
        from larmor.desktop.dialogs import ExperimentDialog

        self.snapshot()
        old_sr = self.recipe.get("sr_hz", 0.0) or 0.0
        dlg = ExperimentDialog(self, self.recipe)
        if dlg.exec():
            # a changed SR re-references the ppm axis by ΔSR / SFO1
            new_sr = self.recipe.get("sr_hz", 0.0) or 0.0
            larmor = self.recipe.get("larmor_frequency_MHz", 0.0) or 0.0
            if larmor and abs(new_sr - old_sr) > 1e-9 and self.exp_ppm.size:
                d_ppm = (new_sr - old_sr) / larmor
                self.exp_ppm = self.exp_ppm + d_ppm
                if self._proc_base is not None:
                    self._proc_base = (self._proc_base[0] + d_ppm,
                                       self._proc_base[1])
                self.view.set_experiment(self.exp_ppm, self.exp_amp)
            self._update_exp_label()
            self.statusBar().showMessage(
                "experiment updated — re-simulating (a new spin rate builds "
                "a new kernel once)")
            self.request_simulation()
            self._persist_session()
        else:
            self.undo_stack.pop()   # dialog cancelled: drop the snapshot

    # ------------------------------------------------------------- baseline
    def _baseline_mode(self, on: bool):
        self.view.set_baseline_mode(on)
        self._set_add_mode(None)
        if on:
            self.statusBar().showMessage(
                "manual baseline: click to place anchors, drag to shape; "
                "then 'Subtract'")

    def apply_manual_baseline(self):
        base = self.view.baseline_curve(self.exp_ppm)
        if base is None:
            self.statusBar().showMessage("place at least 2 baseline anchors first")
            return
        self.exp_amp = self.exp_amp - base
        self.view.clear_baseline()
        self.proc_panel.btnBlPick.setChecked(False)
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.request_simulation()
        self.statusBar().showMessage(
            "manual baseline subtracted ('Reset to original' undoes)")

    # ------------------------------------------------------------- zones
    def add_zone(self):
        if self.recipe is None:
            return
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        width = abs(x1 - x0) * 0.25
        center = (x0 + x1) / 2.0
        zones = self.recipe.get("fit_zones") or []
        zones.append([center + width / 2.0, center - width / 2.0])
        self.recipe["fit_zones"] = zones
        self._sync_zones()
        self.statusBar().showMessage(
            "fit zone added — drag its edges; the fit uses the union of all "
            "zones instead of the zoom")

    def clear_zones(self):
        if self.recipe is None:
            return
        self.recipe["fit_zones"] = []
        self._sync_zones()

    def _sync_zones(self):
        zones = (self.recipe or {}).get("fit_zones") or []
        self.view.set_zones(zones, on_change=self._zones_changed)
        self._persist_session()

    def _zones_changed(self, values: list):
        if self.recipe is not None:
            self.recipe["fit_zones"] = values
            self._persist_session()

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
        """Fit Y to everything visible in the current X window: experiment,
        model, and the residual offset below zero."""
        if not self.exp_ppm.size:
            return
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        lo_x, hi_x = min(x0, x1), max(x0, x1)
        sel = (self.exp_ppm >= lo_x) & (self.exp_ppm <= hi_x)
        if not sel.any():
            return
        lo = float(self.exp_amp[sel].min())
        hi = float(self.exp_amp[sel].max())
        if self._last_model is not None:
            mx, my = self._last_model
            msel = (mx >= lo_x) & (mx <= hi_x)
            if msel.any():
                hi = max(hi, float(np.max(my[msel])))
                lo = min(lo, float(np.min(my[msel])))
        lo = min(lo, -0.12 * hi)          # room for the offset residual
        pad = 0.08 * (hi - lo or 1.0)
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
            self, "Open spectrum, recipe, or Bruker file", self._last_dir(),
            "All supported (*.fxmla *.fxml *.json 1r 2rr fid ser);;"
            "dmfit / recipe (*.fxmla *.fxml *.json);;"
            "Bruker processed (1r 2rr);;Bruker raw (fid ser);;All files (*)")
        if path:
            self.load_source(path)

    def open_expno(self):
        path = QFileDialog.getExistingDirectory(
            self, "Open Bruker EXPNO or pdata folder (read-only)",
            self._last_dir())
        if path:
            self.load_source(path)

    def open_sample(self):
        path = QFileDialog.getExistingDirectory(
            self, "Choose a sample folder (lists all its spectra)",
            self._last_dir())
        if path:
            self.explorer.load_sample(path)
            self.explorer_dock.show()
            self.explorer_dock.raise_()
            self.statusBar().showMessage(
                "sample scanned — double-click a spectrum in the Explorer "
                "to open it")

    def open_fid(self):
        from larmor.desktop.fid_dialog import FidDialog

        path = None
        if self.source_path and Path(self.source_path).is_dir() and \
                ((Path(self.source_path) / "fid").exists() or
                 (Path(self.source_path) / "ser").exists()):
            path = str(Path(self.source_path) /
                       ("ser" if (Path(self.source_path) / "ser").exists()
                        else "fid"))
        dlg = FidDialog(self, path)
        dlg.accepted_1d.connect(self._fid_to_workbench)
        dlg.accepted_2d.connect(self._fid_to_2d)
        dlg.exec()

    def _fid_to_workbench(self, ppm, amp, meta):
        """A 1D spectrum processed from a raw fid becomes the working data."""
        from larmor.recipe import Recipe

        order = np.argsort(ppm)
        self.exp_ppm, self.exp_amp = np.asarray(ppm)[order], np.asarray(amp)[order]
        self._proc_base = None
        self.source_path = meta.get("expno", "")
        self.recipe = Recipe(
            sample=(meta.get("title", "").splitlines() or [""])[0],
            source_kind="bruker", source_path=meta.get("expno", ""),
            nucleus=meta.get("nucleus", ""),
            larmor_frequency_MHz=meta.get("larmor_MHz", 0.0),
            spin_rate_Hz=meta.get("masr_Hz") or 0.0).to_dict()
        self.hidden.clear(); self.undo_stack.clear(); self.redo_stack.clear()
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.view.set_title(self.recipe.get("sample") or "processed FID")
        self._last_model = None
        self._first_sim = True
        self.zoom_full()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles(); self._update_exp_label(); self._update_enabled()
        self.statusBar().showMessage(
            "spectrum from processed FID loaded — add lines and fit")

    def _fid_to_2d(self, data2d):
        from larmor.desktop.twod_dialog import TwoDDialog

        dlg = TwoDDialog(self, None)
        dlg.data = data2d.normalized()
        dlg.lbl.setText(data2d.source or "processed 2D FID")
        dlg._redraw()
        dlg.exec()

    def _last_dir(self) -> str:
        return QSettings("LARMOR", "app").value("lastDir", "")

    # ------------------------------------------------------- workspaces
    def _doc_title(self) -> str:
        if self.central_stack.currentWidget() is self.view2d:
            return self.view2d.title.text() or "2D"
        if self.recipe:
            return (self.recipe.get("sample")
                    or (Path(self.source_path).name if self.source_path else "spectrum"))
        return "empty"

    def _snapshot_doc(self) -> dict:
        is2d = self.central_stack.currentWidget() is self.view2d
        snap = {"kind": "2d" if is2d else "1d",
                "source_path": self.source_path,
                "recipe": json.loads(json.dumps(self.recipe)) if self.recipe else None,
                "hidden": set(self.hidden)}
        if is2d:
            snap["data2d"] = self._data2d
            snap["view2d"] = self.view2d.get_state()
            snap["fittable"] = getattr(self, "_data2d_fittable", False)
        else:
            snap["exp_ppm"] = np.array(self.exp_ppm, copy=True)
            snap["exp_amp"] = np.array(self.exp_amp, copy=True)
            snap["proc_base"] = self._proc_base
            snap["overlays"] = [dict(o) for o in self._overlays]
        return snap

    def _apply_doc(self, snap: dict):
        self.source_path = snap["source_path"]
        self.recipe = (json.loads(json.dumps(snap["recipe"]))
                       if snap["recipe"] else None)
        self.hidden = set(snap["hidden"])
        self.undo_stack.clear(); self.redo_stack.clear()
        if snap["kind"] == "2d":
            self._data2d = snap["data2d"]
            self._data2d_fittable = snap.get("fittable", False)
            self.view2d.set_state(snap["view2d"])
            self.central_stack.setCurrentWidget(self.view2d)
        else:
            self.exp_ppm = snap["exp_ppm"]; self.exp_amp = snap["exp_amp"]
            self._proc_base = snap["proc_base"]
            self._overlays = list(snap["overlays"])
            self.central_stack.setCurrentWidget(self.view)
            self.view.set_experiment(self.exp_ppm, self.exp_amp)
            self.view.set_title(self.recipe.get("sample", "") if self.recipe else "")
            self.lines_table.rebuild(self.recipe, self.hidden)
            self._update_paddles(); self._refresh_overlays(); self._update_sn()
            if self.recipe and self.recipe.get("sites"):
                self.request_simulation()
            else:
                self.view.set_model(None, None, None, None, self.hidden)
        self._update_exp_label(); self._update_enabled()

    def _sync_active(self):
        if self.active_ws is None or not (0 <= self.active_ws < len(self.workspaces)):
            return
        ws = self.workspaces[self.active_ws]
        ws["snap"] = self._snapshot_doc()
        ws["kind"] = ws["snap"]["kind"]
        ws["has_fit"] = bool(self.recipe and self.recipe.get("sites"))
        ws["title"] = self._doc_title()

    def _register_ws(self, kind: str):
        mode, self._ws_mode = self._ws_mode, "auto"
        entry = {"snap": self._snapshot_doc(), "kind": kind,
                 "title": self._doc_title(),
                 "has_fit": bool(self.recipe and self.recipe.get("sites"))}
        reuse = False
        if mode == "reuse":
            reuse = self.active_ws is not None
        elif mode == "auto" and self.active_ws is not None:
            cur = self.workspaces[self.active_ws]
            reuse = (kind == "1d" and cur["kind"] == "1d" and not cur["has_fit"])
        if reuse:
            self.workspaces[self.active_ws] = entry
        else:
            self.workspaces.append(entry)
            self.active_ws = len(self.workspaces) - 1
        self._refresh_ws_panel()

    def _refresh_ws_panel(self):
        items = []
        for ws in self.workspaces:
            icon = "▦" if ws["kind"] == "2d" else ("⤳" if ws["has_fit"] else "∿")
            items.append((icon, ws["title"]))
        self.ws_panel.rebuild(items, self.active_ws if self.active_ws is not None
                              else -1)

    def switch_workspace(self, i: int):
        if i == self.active_ws or not (0 <= i < len(self.workspaces)):
            return
        self._sync_active()
        self.active_ws = i
        self._apply_doc(self.workspaces[i]["snap"])
        self._refresh_ws_panel()
        self.statusBar().showMessage(f"workspace: {self.workspaces[i]['title']}")

    def close_workspace(self, i: int):
        if not (0 <= i < len(self.workspaces)):
            return
        del self.workspaces[i]
        if not self.workspaces:
            self.active_ws = None
            self._refresh_ws_panel()
            return
        if self.active_ws == i:
            self.active_ws = min(i, len(self.workspaces) - 1)
            self._apply_doc(self.workspaces[self.active_ws]["snap"])
        elif self.active_ws is not None and self.active_ws > i:
            self.active_ws -= 1
        self._refresh_ws_panel()

    def save_workspace(self, i: int):
        if not (0 <= i < len(self.workspaces)):
            return
        self.switch_workspace(i)
        if self.central_stack.currentWidget() is self.view2d:
            self.save_recipe()
        elif self.recipe and self.recipe.get("sites"):
            self.save_recipe()
        else:
            self.save_spectrum()

    def load_source(self, path: str, keep_fit: bool | None = None):
        """Open ANY dataset with a basic display, then let the user choose a
        process. 1D spectra go to the fit workbench; 2D datasets show a contour
        map; raw fid/ser get a quick preview. Nothing is ever rejected."""
        self.statusBar().showMessage("loading…")
        QApplication.processEvents()
        self._sync_active()               # persist the outgoing document
        self._in_load_source = True
        try:
            self._load_source_body(path, keep_fit)
        finally:
            self._in_load_source = False
        if Path(path).exists():
            self._add_recent(path)

    def _load_source_body(self, path: str, keep_fit: bool | None):
        # First try the 1D-fittable path (dmfit / recipe / 1D processed).
        try:
            ppm, amp, recipe, meta, warnings = _load_any(path)
        except ValueError:
            # not a 1D-fittable source: it may be a 2D or a raw fid/ser.
            handled = self._load_nonfittable(path)
            if handled:
                return
            # fall through to report the original error
            try:
                _load_any(path)
            except Exception as exc:
                QMessageBox.warning(self, "Load failed", str(exc))
                self.statusBar().showMessage("load failed")
            return
        except Exception as exc:
            QMessageBox.warning(self, "Load failed", str(exc))
            self.statusBar().showMessage("load failed")
            return

        self.central_stack.setCurrentWidget(self.view)   # 1D workbench

        # dmfit behaviour: if a fit is already on screen and the new source
        # brings no fit of its own, offer to keep the current lines
        existing = self.recipe.get("sites") if self.recipe else None
        incoming = recipe.get("sites")
        if keep_fit is None and existing and not incoming:
            btn = QMessageBox.question(
                self, "Keep fit parameters?",
                "A fit is already open. Keep the current lines and fit them "
                "against the new spectrum?\n\n"
                "Yes = keep the lines · No = start empty",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            keep_fit = btn == QMessageBox.Yes
        if keep_fit and existing and not incoming:
            recipe["sites"] = existing
            # carry the experiment parameters the new source knows
            for k in ("nucleus", "larmor_frequency_MHz", "spin_rate_Hz"):
                recipe[k] = recipe.get(k) or self.recipe.get(k)

        self.source_path = path
        self.exp_ppm, self.exp_amp = ppm, amp
        self._proc_base = None
        self.recipe = recipe
        self.hidden.clear()
        self.undo_stack.clear()
        self.redo_stack.clear()
        QSettings("LARMOR", "app").setValue("lastDir", str(Path(path).parent))
        self.setWindowTitle(f"LARMOR — {Path(path).name}")
        self.view.set_experiment(ppm, amp)
        if not self.recipe["sites"]:
            # drop any model curve left over from the previous spectrum
            self.view.set_model(None, None, None, None, self.hidden)
        sample = recipe.get("sample") or Path(path).name
        self.view.set_title(sample)
        self._last_model = None
        self._first_sim = True
        self.zoom_full()
        if recipe["sites"]:
            self.zoom_sites()
        self._update_sn()
        msg = meta + ("   ⚠ " + " • ".join(warnings) if warnings else "")
        self.statusBar().showMessage(msg)
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        self._update_exp_label()
        self._sync_zones()
        self._update_enabled()
        if self.recipe["sites"]:
            self.request_simulation()
        self._refresh_overlays()
        self._register_ws("1d")
        self._persist_session()

    def _load_nonfittable(self, path: str) -> bool:
        """Display a 2D dataset or a raw fid/ser without rejecting it.
        Returns True if it handled the path."""
        from larmor.io import bruker

        try:
            ref = bruker.resolve(path)
            data = bruker.read(path)
        except Exception:
            return False

        self.source_path = path
        QSettings("LARMOR", "app").setValue("lastDir", str(Path(path).parent))
        self.setWindowTitle(f"LARMOR — {Path(path).name}")

        if data.ndim == 2 and data.domain == "freq":
            self._show_2d(_nmrdata_to_data2d(data), Path(path).name,
                          "2D spectrum")
            kind = "arrayed/relaxation" if data.is_pseudo2d else "MQMAS 2D"
            self.statusBar().showMessage(
                f"{kind} displayed — Tools ▸ 2D MQMAS to fit · Tools ▸ "
                "Relaxation for a series · or send a 1D trace to fitting below")
            return True

        if data.ndim == 2 and data.domain == "time":
            from larmor import fourier

            d2 = fourier.ft2d_from_nmrdata(data, fourier.FT2DParams(
                f2_ops=[{"op": "fcor", "factor": 0.5}, {"op": "em", "lb_hz": 100}]))
            self._show_2d(d2, Path(path).name, "raw 2D (quick preview)")
            self.statusBar().showMessage(
                "raw 2D preview — File ▸ Open FID for full processing · "
                "Tools ▸ Relaxation for a T1/T2 series")
            return True

        if data.ndim == 1 and data.domain == "time":
            # quick magnitude-FT preview so the user sees a spectrum at once
            from larmor import fourier

            ppm, spec = fourier.ft1d(
                data.data, data.axes[0].sw_Hz, data.meta["larmor_MHz"],
                ops=[{"op": "fcor", "factor": 0.5}, {"op": "em", "lb_hz": 100}])
            order = np.argsort(ppm)
            spec = np.abs(spec)              # magnitude preview (phase-free)
            self._display_1d(ppm[order], spec[order],
                             data.meta.get("nucleus", ""),
                             data.meta["larmor_MHz"], data.meta.get("masr_Hz"),
                             Path(path).name + " (FID preview)", str(ref.expno))
            self.statusBar().showMessage(
                "raw FID preview (magnitude) — Process ▸ Open FID to apodize, "
                "phase and transform properly, then fit")
            return True
        return False

    # ------------------------------------------------------------------- S/N
    def _update_sn(self):
        """Signal-to-noise: peak signal over the RMS of a signal-free region.
        The noise region is taken from the quiet outer edges of the spectrum
        (robust to where the peak sits), matching TopSpin's sino spirit."""
        y = self.exp_amp
        if y is None or y.size < 20:
            self.lines_table.set_sn("")
            return
        y = np.asarray(y, float)
        edge = max(5, y.size // 10)
        noise_region = np.concatenate([y[:edge], y[-edge:]])
        noise = float(np.std(noise_region - np.median(noise_region)))
        signal = float(np.max(np.abs(y - np.median(noise_region))))
        sn = signal / noise if noise > 0 else 0.0
        self.lines_table.set_sn(f"S/N {sn:,.0f}" if sn else "")

    # --------------------------------------------------------- overlays
    def add_overlay_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Add a spectrum to compare", self._last_dir(),
            "Spectra (*.fxmla *.json 1r 2rr *.txt *.csv);;All files (*)")
        if not path:
            return
        try:
            recipe, ppm, amp, *_ = _load_any(path)
            label = recipe.get("sample") or Path(path).name
        except Exception:
            try:
                from larmor.io import bruker

                d = bruker.read(path)
                if d.ndim != 1 or d.domain != "freq":
                    raise ValueError("not a 1D spectrum")
                ppm, amp = np.asarray(d.axes[0].values), np.asarray(d.data, float)
                label = Path(path).name
            except Exception as exc:
                QMessageBox.warning(self, "Add overlay", f"cannot read: {exc}")
                return
        self._add_overlay(label, np.asarray(ppm), np.asarray(amp), path)

    def _add_overlay(self, label, ppm, amp, source=""):
        from larmor.desktop.datasets import overlay_color

        self._overlays.append({
            "label": label, "ppm": np.asarray(ppm), "amp": np.asarray(amp),
            "color": overlay_color(len(self._overlays)), "visible": True,
            "source": source})
        self._refresh_overlays()
        self.datasets_dock.raise_()

    def overlay_remove(self, i: int):
        if 0 <= i < len(self._overlays):
            del self._overlays[i]
            self._refresh_overlays()

    def overlay_visibility(self, i: int, on: bool):
        if 0 <= i < len(self._overlays):
            self._overlays[i]["visible"] = on
            self._refresh_overlays()

    def overlay_make_active(self, i: int):
        if not (0 <= i < len(self._overlays)):
            return
        ov = self._overlays[i]
        if not ov.get("source"):
            return
        prev = None
        if self.exp_ppm.size and self.source_path:
            prev = (self.recipe.get("sample") if self.recipe else None,
                    self.exp_ppm.copy(), self.exp_amp.copy(), self.source_path)
        del self._overlays[i]
        self._ws_mode = "reuse"           # swap the active spectrum in place
        self.load_source(ov["source"])
        if prev and prev[3] != self.source_path:
            self._add_overlay(prev[0] or Path(prev[3]).name, prev[1], prev[2],
                              prev[3])
        self._refresh_overlays()

    def _refresh_overlays(self):
        if not hasattr(self, "datasets_panel"):
            return
        span = 1.0
        if self.exp_amp.size:
            span = float(np.nanmax(self.exp_amp) - np.nanmin(self.exp_amp)) or 1.0
        step = self.datasets_panel.offset.value() * span
        drawn = []
        for k, ov in enumerate(self._overlays):
            if ov.get("visible", True):
                drawn.append((ov["ppm"], ov["amp"] + step * (k + 1),
                              ov["color"], ov["label"]))
        self.view.set_overlays(drawn)
        label = ""
        if self.recipe is not None:
            label = self.recipe.get("sample") or (
                Path(self.source_path).name if self.source_path else "")
        self.datasets_panel.rebuild(label, self._overlays)

    #: colour assigned to each HMQC projection axis (matches the overlay + the
    #: Explorer highlight)
    PROJ_COLOR = {"f2": "#e8832a", "f1": "#6a4fb0"}

    def load_projection_1d(self, axis: str):
        """HMQC: arm a pick — the next spectrum clicked in the Explorer (or via
        its Browse…) is overlaid on the F2/F1 projection and highlighted."""
        self._proj_pick_axis = axis
        self.explorer_dock.show(); self.explorer_dock.raise_()
        self.statusBar().showMessage(
            f"click a spectrum in the Explorer for the {axis.upper()} projection "
            "(use Browse… to reach one elsewhere)")

    def _explorer_open(self, path: str):
        """Route an Explorer activation: an HMQC projection pick, else a load."""
        axis = self._proj_pick_axis
        if axis is None:
            self.load_source(path)
            return
        self._proj_pick_axis = None
        try:
            ppm, amp = self._read_1d(path)
        except Exception as exc:
            QMessageBox.warning(self, "Projection 1D", f"cannot read: {exc}")
            return
        self.view2d.set_projection_1d(axis, ppm, amp,
                                      color=self.PROJ_COLOR[axis])
        self.explorer.highlight(path, self.PROJ_COLOR[axis])
        self.statusBar().showMessage(
            f"{axis.upper()} 1D overlaid — adjust scale, then 'uncorrelated "
            f"{axis.upper()} →' for the non-correlated features")

    def _read_1d(self, path: str):
        """A 1D (ppm, amp) from any supported source."""
        try:
            ppm, amp, *_ = _load_any(path)
            return np.asarray(ppm), np.asarray(amp)
        except Exception:
            from larmor.io import bruker

            d = bruker.read(path)
            if d.ndim != 1 or d.domain != "freq":
                raise ValueError("not a 1D spectrum")
            return np.asarray(d.axes[0].values), np.asarray(d.data, float)

    def back_to_2d(self):
        """Switch to the most recent 2D workspace (the map you came from)."""
        if self.central_stack.currentWidget() is self.view2d:
            return
        for i in range(len(self.workspaces) - 1, -1, -1):
            if self.workspaces[i]["kind"] == "2d":
                self.switch_workspace(i)
                return
        self.statusBar().showMessage("no 2D workspace is open")

    def _show_2d(self, data2d, title: str, kind: str):
        from larmor.recipe import Recipe

        if not self._in_load_source:
            self._sync_active()
        self._data2d = data2d
        self.view2d.set_data(data2d, f"{kind} — {title}")
        self.view2d.clear_model()
        self.central_stack.setCurrentWidget(self.view2d)
        # a genuine spectroscopic 2D (not a relaxation array) is fittable, so
        # give it a recipe if there isn't a fit already in progress
        pseudo = any("pseudo" in n or "arrayed" in n for n in data2d.notes)
        self._data2d_fittable = not pseudo
        if not pseudo and (self.recipe is None or not self.recipe.get("sites")):
            self.recipe = Recipe(
                sample=title, source_kind="bruker", source_path=self.source_path,
                nucleus=data2d.nucleus,
                larmor_frequency_MHz=data2d.larmor_MHz).to_dict()
            self.hidden.clear()
            self.lines_table.rebuild(self.recipe, self.hidden)
            self._update_exp_label(); self._update_enabled()
        self._register_ws("2d")

    def _display_1d(self, ppm, amp, nucleus, larmor, masr, title, expno):
        """Put a bare 1D spectrum on the workbench (no fit), ready to fit."""
        from larmor.recipe import Recipe

        if not self._in_load_source:
            self._sync_active()
        self.central_stack.setCurrentWidget(self.view)
        self.exp_ppm, self.exp_amp = np.asarray(ppm), np.asarray(amp)
        self._proc_base = None
        self.recipe = Recipe(sample=title, source_kind="bruker",
                             source_path=expno, nucleus=nucleus,
                             larmor_frequency_MHz=larmor,
                             spin_rate_Hz=masr or 0.0).to_dict()
        self.hidden.clear(); self.undo_stack.clear(); self.redo_stack.clear()
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.view.set_model(None, None, None, None, self.hidden)  # clear old model
        self.view.set_title(title)
        self._last_model = None; self._first_sim = True
        self.zoom_full()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles(); self._update_exp_label(); self._update_enabled()
        self._refresh_overlays(); self._update_sn()
        self._register_ws("1d")

    def _trace_to_workbench(self, ppm, amp, label):
        """A 1D trace pulled out of the 2D view becomes a NEW workspace, so the
        2D map stays open in its own workspace."""
        nuc = self.view2d.data.nucleus if self.view2d.data else ""
        larmor = self.view2d.data.larmor_MHz if self.view2d.data else 0.0
        self._ws_mode = "new"
        self._display_1d(ppm, amp, nuc, larmor, None, label,
                         self.source_path or "")
        self.statusBar().showMessage(f"{label} — new workspace; add lines and Fit")

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
        self.view2d.set_add_mode(name)
        if name:
            self.statusBar().showMessage(
                f"placing {name} lines — click to add as many as you want; "
                f"click {name} again (or press Esc) to stop")

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
        # stay in placement mode: drop as many lines as wanted, click the model
        # again (or Esc) to leave the mode
        self.on_structure_changed()
        self.statusBar().showMessage(
            f"added {name} #{n} — click to add more, or click {name} again "
            "(Esc) to stop")

    def add_function_line(self):
        """Add a user y(x; a,b,c,d) expression line (ssNake Function fit)."""
        from PySide6.QtWidgets import QInputDialog

        if self.recipe is None:
            self.statusBar().showMessage("load a spectrum first")
            return
        expr, ok = QInputDialog.getText(
            self, "Function line",
            "y = f(x; a, b, c, d)   (numpy: exp, sin, sqrt, pi …):",
            text="a * exp(-((x - b) / c)**2) + d")
        if not ok or not expr.strip():
            return
        self.snapshot()
        m = model_registry.get("function")
        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        n = len(self.recipe["sites"])
        self.recipe["sites"].append(
            {"model": "function", "label": f"fn-{n}", "func": expr.strip(),
             "params": params})
        self.on_structure_changed()
        self.statusBar().showMessage(f"added function line: {expr.strip()}")

    def autopick_lines(self):
        """Peak-pick the spectrum and drop a Gauss/Lorentz line at each peak."""
        from PySide6.QtWidgets import QInputDialog

        if self.recipe is None or not self.exp_ppm.size:
            self.statusBar().showMessage("load a spectrum first")
            return
        thr, ok = QInputDialog.getDouble(
            self, "Add lines at peaks", "Threshold (% of max):", 5.0, 0.1, 100.0, 1)
        if not ok:
            return
        from larmor import processing as proc

        span = float(np.ptp(self.exp_ppm)) if self.exp_ppm.size else 0.0
        peaks = proc.pick_peaks(self.exp_ppm, self.exp_amp,
                                threshold_frac=thr / 100.0,
                                min_sep_ppm=span / 100.0)
        if not peaks:
            self.statusBar().showMessage("no peaks above that threshold")
            return
        self.snapshot()
        m = model_registry.get("gauss_lor")
        for pk in peaks:
            params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                               "min": p.min, "max": p.max, "expr": None}
                      for p in m.params}
            params["isotropic_chemical_shift_ppm"]["value"] = pk["ppm"]
            params["amplitude"]["value"] = abs(pk["height"])
            params["shift_fwhm_ppm"]["value"] = max(pk.get("fwhm_ppm") or 2.0, 0.2)
            n = len(self.recipe["sites"])
            self.recipe["sites"].append(
                {"model": "gauss_lor", "label": f"pk-{n}", "params": params})
        self.on_structure_changed()
        self.statusBar().showMessage(f"added {len(peaks)} lines at peaks")

    def predict_at_field(self):
        """Simulate the current model at a different field (teaching/planning:
        quadrupolar 2nd-order width ~ 1/B0, chemical shift constant in ppm)."""
        import json
        from PySide6.QtWidgets import QInputDialog

        from larmor import engine, nuclei as N
        from larmor.recipe import Recipe

        if not (self.recipe and self.recipe.get("sites")):
            self.statusBar().showMessage("build a model first")
            return
        nuc = self.recipe.get("nucleus", "")
        cur = self.recipe.get("larmor_frequency_MHz", 0.0) or 0.0
        try:
            iso = next(i for i in N.all_isotopes() if i.symbol == nuc)
            cur_h1 = cur * N.GAMMA_1H / abs(iso.gamma_MHz_T)
        except Exception:
            QMessageBox.warning(self, "Predict", f"unknown nucleus {nuc!r}")
            return
        h1, ok = QInputDialog.getDouble(
            self, "Predict at another field", "Target ¹H frequency (MHz):",
            round(cur_h1) or 400.0, 10.0, 1700.0, 1)
        if not ok:
            return
        new_larmor = h1 * abs(iso.gamma_MHz_T) / N.GAMMA_1H
        rec = Recipe.from_dict(json.loads(json.dumps(self.recipe)))
        rec.larmor_frequency_MHz = new_larmor
        x, total, _ = engine.simulate(
            rec, exp_ppm=self.exp_ppm if self.exp_ppm.size else None)
        self._ws_mode = "new"
        self._display_1d(x, total, nuc, new_larmor, rec.spin_rate_Hz,
                         f"{self.recipe.get('sample') or 'model'} @ {h1:.0f} MHz ¹H",
                         "")
        self.statusBar().showMessage(
            f"predicted at {new_larmor:.1f} MHz ({nuc}, ¹H {h1:.0f})")

    def add_background_spectrum(self):
        """Add another measured spectrum as a fit component (background /
        impurity / reference), scaled by amplitude and shiftable in ppm."""
        if self.recipe is None:
            self.statusBar().showMessage("load a spectrum to fit first")
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Background spectrum to fit", self._last_dir(),
            "Spectra (*.fxmla *.json 1r 2rr *.txt *.csv);;All files (*)")
        if not path:
            return
        try:
            _r, ppm, amp, *_ = _load_any(path)
        except Exception:
            try:
                from larmor.io import bruker

                d = bruker.read(path)
                if d.ndim != 1 or d.domain != "freq":
                    raise ValueError("not a 1D spectrum")
                ppm, amp = np.asarray(d.axes[0].values), np.asarray(d.data, float)
            except Exception as exc:
                QMessageBox.warning(self, "Background spectrum", f"cannot read: {exc}")
                return
        ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
        amp = amp / (np.max(np.abs(amp)) or 1.0)         # unit peak
        self.snapshot()
        m = model_registry.get("spectrum")
        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        # start the amplitude near the data's peak so it is on-scale
        params["amplitude"]["value"] = float(np.max(np.abs(self.exp_amp))
                                              if self.exp_amp.size else 1.0)
        n = len(self.recipe["sites"])
        self.recipe["sites"].append({
            "model": "spectrum", "label": f"bg-{n}",
            "ref": {"ppm": ppm.tolist(), "amp": amp.tolist()},
            "params": params})
        self.on_structure_changed()
        self.statusBar().showMessage(
            f"added background spectrum '{Path(path).name}' — Fit scales and "
            "shifts it")

    def add_site_2d(self, f2_ppm: float, f1_ppm: float):
        """Place a 2D site from a click on the contour: the isotropic shift
        starts at the clicked F1 (isotropic) position."""
        name = next((n for n, a in self._model_actions.items()
                     if a.isChecked()), None)
        if not name or self.recipe is None:
            return
        self.snapshot()
        m = model_registry.get(name)
        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        if "isotropic_chemical_shift_ppm" in params:
            params["isotropic_chemical_shift_ppm"]["value"] = float(f1_ppm)
        if "amplitude" in params:
            params["amplitude"]["value"] = 1.0
        n = len(self.recipe["sites"])
        self.recipe["sites"].append(
            {"model": name, "label": f"{m.label.split(' ')[0]}-{n}",
             "params": params})
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_enabled()
        self.statusBar().showMessage(
            f"added {name} at F1≈{f1_ppm:.1f} ppm — click to add more, or click "
            f"{name} again (Esc) to stop; Decomposition ▸ Fit to fit the 2D")

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
        # the 2D map has its own (explicit) fit path; skip the live 1D sim there
        if self.central_stack.currentWidget() is self.view2d:
            return
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
        self._last_model = (np.asarray(x), np.asarray(total))
        self.view.set_model(x, total, per_site, labels, self.hidden,
                            self.exp_ppm, self.exp_amp)
        if self._first_sim:
            self._first_sim = False
            self.autoscale_y()
        if self._sim_pending:
            self._sim_pending = False
            self._sim_timer.start()

    # ------------------------------------------------------------- fit
    def run_fit(self):
        if not self.recipe or not self.recipe["sites"]:
            self.statusBar().showMessage("add at least one line first")
            return
        if self.central_stack.currentWidget() is self.view2d:
            self.run_fit_2d()
            return
        if self._fit_worker and self._fit_worker.isRunning():
            return
        self.snapshot()
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        hi, lo = max(x0, x1), min(x0, x1)
        zones = self.recipe.get("fit_zones") or []
        self.statusBar().showMessage(
            f"fitting in {len(zones)} zone(s) …" if zones
            else f"fitting in {hi:.1f} … {lo:.1f} ppm …")
        self.lines_table.btnFit.setEnabled(False)
        self._fit_worker = FitWorker(json.loads(json.dumps(self.recipe)),
                                     self.exp_ppm, self.exp_amp, (hi, lo))
        self._fit_worker.done.connect(self._fit_done)
        self._fit_worker.failed.connect(self._fit_failed)
        self._fit_worker.start()

    # ------------------------------------------------------------- 2D fit
    _MODELS_2D = {"czjzek", "ext_czjzek", "quad_ct", "quad_csa"}

    def run_fit_2d(self, method: str = "3QMAS"):
        if self._data2d is None or not getattr(self, "_data2d_fittable", False):
            self.statusBar().showMessage(
                "this 2D is a relaxation array, not an MQMAS map — not fittable")
            return
        bad = [s["model"] for s in self.recipe["sites"]
               if s["model"] not in self._MODELS_2D]
        if bad:
            QMessageBox.warning(
                self, "2D fit", "these models have no 2D implementation: "
                + ", ".join(sorted(set(bad)))
                + "\n2D supports: " + ", ".join(sorted(self._MODELS_2D)))
            return
        if self._fit2d_worker and self._fit2d_worker.isRunning():
            return
        self.snapshot()
        self.statusBar().showMessage(
            "fitting the 2D (building the MQMAS kernel, first time is slow)…")
        self.lines_table.btnFit.setEnabled(False)
        self._fit2d_worker = Fit2DWorker(
            json.loads(json.dumps(self.recipe)), self._data2d, method)
        self._fit2d_worker.done.connect(self._fit2d_done)
        self._fit2d_worker.failed.connect(self._fit_failed)
        self._fit2d_worker.start()

    def _fit2d_done(self, result):
        self.lines_table.btnFit.setEnabled(True)
        self.recipe = result.recipe.to_dict()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self.view2d.set_model(result.z_fit, result.kernel.f2_ppm,
                              result.kernel.f1_ppm)
        self.lines_table.set_chi2(f"RMSD {result.rmsd:.4f}")
        self.results_summary.setText(f"2D MQMAS fit · RMSD {result.rmsd:.4f}")
        self.report.setPlainText(result.report)
        self.statusBar().showMessage(f"2D fit done · RMSD {result.rmsd:.4f}")
        self._persist_session()

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
        self._last_model = (np.asarray(result.x_ppm), np.asarray(result.y_fit))
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
        # only the 1D workbench is processed live; ignore while a 2D map is up
        if self.central_stack.currentWidget() is not self.view:
            return
        if not use_raw and not self.exp_ppm.size:
            return
        from larmor import processing as proc
        from larmor.io import bruker

        live = getattr(self.proc_panel, "chkLive", None)
        live = bool(live and live.isChecked())
        if not live:
            self.statusBar().showMessage("processing…")
            QApplication.processEvents()
        try:
            if use_raw:
                if not bruker.is_expno(Path(self.source_path)):
                    raise ValueError("raw-fid processing needs a Bruker EXPNO")
                s = proc.from_bruker_fid(self.source_path)
            else:
                # apply the pipeline from the UNPROCESSED baseline every time, so
                # a live slider shows the absolute phase rather than compounding
                if self._proc_base is None:
                    self._proc_base = (self.exp_ppm.copy(), self.exp_amp.copy())
                base_ppm, base_amp = self._proc_base
                sfo1 = self.recipe.get("larmor_frequency_MHz", 0.0) if self.recipe else 0.0
                s = proc.from_processed(base_ppm, base_amp, sfo1)
            s = proc.apply(s, ops)
            if s.domain != "freq":
                raise ValueError("pipeline must end in the frequency domain")
        except Exception as exc:
            if not live:                       # never nag on every keystroke
                QMessageBox.warning(self, "Processing failed", str(exc))
            self.statusBar().showMessage("processing failed")
            return
        order = np.argsort(s.x_ppm)
        self.exp_ppm, self.exp_amp = np.asarray(s.x_ppm)[order], s.y.real[order]
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self._update_sn()
        # remember the pipeline in the recipe: saving the fit then saves the
        # processing that produced the spectrum it was fitted against
        if self.recipe is not None:
            self.recipe["processing"] = list(ops)
            self.recipe["processing_from_raw"] = bool(use_raw)
        self.request_simulation()
        self.statusBar().showMessage(
            f"processing applied ({len(ops)} step(s), stored in the recipe)")

    # ------------------------------------------------------------- calibrate
    def start_calibrate(self):
        if not self.exp_ppm.size:
            self.statusBar().showMessage("load a spectrum first")
            return
        self.view.set_calibrate_mode(True)
        self.statusBar().showMessage(
            "calibrate: click the peak whose shift you want to set")

    def on_calibrate_picked(self, peak_ppm: float):
        from PySide6.QtWidgets import QInputDialog

        self.view.set_calibrate_mode(False)
        target, ok = QInputDialog.getDouble(
            self, "Calibrate axis", f"Set the peak at {peak_ppm:.2f} ppm to:",
            peak_ppm, -100000.0, 100000.0, 3)
        if not ok:
            return
        delta = float(target) - float(peak_ppm)
        if abs(delta) < 1e-12:
            return
        self.exp_ppm = self.exp_ppm + delta
        if self._proc_base is not None:
            self._proc_base = (self._proc_base[0] + delta, self._proc_base[1])
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        larmor = self.recipe.get("larmor_frequency_MHz", 0.0) if self.recipe else 0.0
        if self.recipe is not None:
            self.recipe["calibration_ppm"] = \
                self.recipe.get("calibration_ppm", 0.0) + delta
        self.request_simulation()
        hz = delta * larmor
        self.statusBar().showMessage(
            f"axis shifted by {delta:+.3f} ppm"
            + (f"  (SR {hz:+.1f} Hz)" if larmor else ""))

    # --------------------------------------------------------------- measure
    def toggle_measure(self, on: bool):
        self.view.set_measure_mode(on)
        if not on:
            self.statusBar().showMessage("")

    def on_measure_changed(self, p1: float, p2: float):
        larmor = self.recipe.get("larmor_frequency_MHz", 0.0) if self.recipe else 0.0
        dppm = abs(p1 - p2)
        msg = f"Δ = {dppm:.3f} ppm"
        if larmor:
            msg += f"   {dppm * larmor:.1f} Hz"
        msg += f"   ({p1:.2f} → {p2:.2f} ppm)"
        self.statusBar().showMessage(msg)

    def reset_processing(self):
        if self.source_path:
            keep = json.loads(json.dumps(self.recipe)) if self.recipe else None
            self._ws_mode = "reuse"        # reload in place, same workspace
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

    def _guard_write(self, target: Path) -> bool:
        for parent in [target.parent, *target.parent.parents]:
            if (parent / "acqus").exists() or (parent / "fid").exists() \
                    or (parent / "ser").exists():
                QMessageBox.warning(
                    self, "Refused",
                    f"{parent} is an instrument data folder — pick another "
                    "location. LARMOR never writes next to raw data.")
                return False
        return True

    def save_fit_as(self):
        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("nothing to export — add and fit lines first")
            return
        from larmor.io import export

        default = (self.recipe.get("sample") or "fit").strip()
        default = "".join(c if c.isalnum() or c in "-_" else "_"
                          for c in default)[:40] or "fit"
        filters = ";;".join(
            f"{name} (*.{ext})" for name, (ext, _) in export.FORMATS.items())
        path, chosen = QFileDialog.getSaveFileName(
            self, "Save fit as", str(Path(self._last_dir()) / default), filters)
        if not path:
            return
        # figure out the format from the chosen filter (or the extension)
        fmt = None
        for name, (ext, fn) in export.FORMATS.items():
            if name == chosen or path.lower().endswith("." + ext):
                fmt = (ext, fn)
                break
        if fmt is None:
            fmt = ("json", export.export_json)
        ext, fn = fmt
        target = Path(path)
        if target.suffix.lower() != "." + ext:
            target = target.with_suffix("." + ext)
        if not self._guard_write(target):
            return
        recipe = Recipe.from_dict(self.recipe)
        try:
            if ext in ("txt",):
                fn(recipe, self.exp_ppm, self.exp_amp, target)
            elif ext == "fxmla":
                fn(recipe, self.exp_ppm, self.exp_amp, target)
            else:                        # csv, json take (recipe, path)
                fn(recipe, target)
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
            return
        self.statusBar().showMessage(f"fit exported: {target}")

    def open_figure_dialog(self):
        from larmor.desktop.figure_dialog import FigureDialog

        FigureDialog(self, self.source_path, self.recipe).exec()

    def open_subtract(self):
        from larmor.desktop.subtract_dialog import SubtractDialog

        if not self.exp_ppm.size:
            self.statusBar().showMessage("load a spectrum first")
            return
        meta = {"nucleus": self.recipe.get("nucleus", "") if self.recipe else "",
                "larmor_MHz": (self.recipe.get("larmor_frequency_MHz", 0.0)
                               if self.recipe else 0.0)}
        dlg = SubtractDialog(self, self.exp_ppm, self.exp_amp, meta)
        dlg.applied.connect(self._subtract_applied)
        dlg.exec()

    def _subtract_applied(self, ppm, amp):
        self.snapshot()
        self.exp_ppm, self.exp_amp = np.asarray(ppm), np.asarray(amp)
        self._proc_base = None
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        if self.recipe is not None:
            base = self.recipe.get("sample") or "spectrum"
            if "− background" not in base and "minus" not in base:
                self.recipe["sample"] = base + " − background"
            self.view.set_title(self.recipe["sample"])
        if not (self.recipe and self.recipe.get("sites")):
            self.view.set_model(None, None, None, None, self.hidden)
        else:
            self.request_simulation()
        self._update_sn(); self._refresh_overlays()
        self.statusBar().showMessage(
            "background subtracted — File ▸ Save spectrum as… to keep the result")

    def save_spectrum(self):
        from larmor.io import spectra

        if not self.exp_ppm.size:
            self.statusBar().showMessage("no spectrum to save")
            return
        default = (self.recipe.get("sample") if self.recipe else "") or "spectrum"
        default = "".join(c if c.isalnum() or c in "-_ " else "_"
                          for c in default).strip()[:60] or "spectrum"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save spectrum", str(Path(self._last_dir()) / f"{default}.csv"),
            "LARMOR spectrum (*.csv);;Text (*.txt)")
        if not path:
            return
        meta = {}
        if self.recipe:
            meta = {"nucleus": self.recipe.get("nucleus", ""),
                    "larmor_MHz": self.recipe.get("larmor_frequency_MHz", 0.0),
                    "spin_rate_Hz": self.recipe.get("spin_rate_Hz", 0.0),
                    "sample": self.recipe.get("sample", "")}
        spectra.write_csv(path, self.exp_ppm, self.exp_amp, meta)
        self.statusBar().showMessage(f"spectrum saved to {Path(path).name} "
                                     "(reopen it with File ▸ Open…)")

    def open_cofit(self):
        from larmor.desktop.cofit_dialog import CofitDialog

        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage(
                "set up the shared fit (add lines) on one dataset first")
            return
        # the current dataset becomes co-fit dataset 0
        if self.central_stack.currentWidget() is self.view2d and self._data2d:
            base = {"kind": "2d", "label": Path(self.source_path or "2D").name,
                    "data2d": self._data2d, "nucleus": self._data2d.nucleus,
                    "larmor": self._data2d.larmor_MHz}
        else:
            base = {"kind": "1d",
                    "label": self.recipe.get("sample") or "current",
                    "ppm": self.exp_ppm, "amp": self.exp_amp,
                    "nucleus": self.recipe.get("nucleus", ""),
                    "larmor": self.recipe.get("larmor_frequency_MHz", 0.0)}
        dlg = CofitDialog(self, self.recipe, base)
        dlg.applied.connect(self._cofit_applied)
        dlg.exec()

    def _cofit_applied(self, recipe_dict):
        self.snapshot()
        self.recipe = recipe_dict
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        if self.central_stack.currentWidget() is self.view:
            self.request_simulation()
        self.statusBar().showMessage("co-fit shared parameters applied")

    def open_per_site_relaxation(self):
        """Decompose every relaxation slice on the CURRENT fit's lineshapes → a
        T1/T2 per site (not per integration window). Needs a fitted recipe and a
        relaxation ser."""
        from PySide6.QtCore import Qt

        if not (self.recipe and self.recipe.get("sites")):
            self.statusBar().showMessage(
                "fit the most-relaxed slice first — its lines define the sites")
            return
        expno = None
        if self.source_path:
            try:
                from larmor.io import bruker

                ref = bruker.resolve(self.source_path)
                if (ref.expno / "ser").exists():
                    expno = str(ref.expno)
            except Exception:
                expno = None
        if expno is None:
            expno = QFileDialog.getExistingDirectory(
                self, "Relaxation EXPNO (ser + vdlist)")
        if not expno:
            return
        from larmor.recipe import Recipe
        from larmor import series

        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.statusBar().showMessage("decomposing every slice on the fitted lines…")
        QApplication.processEvents()
        try:
            results = series.analyze_per_site(expno, Recipe.from_dict(self.recipe))
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            QMessageBox.warning(self, "Per-site relaxation", str(exc))
            return
        QApplication.restoreOverrideCursor()
        msg = "\n".join(r.summary for r in results) or "no results"
        QMessageBox.information(self, "Per-site relaxation (τ per site)", msg)
        self.statusBar().showMessage("per-site relaxation done")

    def open_vt(self):
        from larmor.desktop.vt_dialog import VtDialog

        VtDialog(self).exec()

    def open_qcpmg(self):
        from larmor.desktop.qcpmg_dialog import QcpmgDialog

        # prefer the raw fid of whatever EXPNO is loaded
        source = None
        if self.source_path:
            try:
                from larmor.io import bruker

                ref = bruker.resolve(self.source_path)
                if (ref.expno / "fid").exists():
                    source = str(ref.expno / "fid")
            except Exception:
                source = None
        dlg = QcpmgDialog(self, source)
        dlg.accepted_1d.connect(self._fid_to_workbench)
        dlg.exec()

    def open_satrec(self):
        from larmor.desktop.satrec_dialog import SatrecDialog

        # resolve whatever is loaded (a 2rr/ser file, a pdata or EXPNO folder)
        # up to the EXPNO that owns the raw ser
        expno = None
        if self.source_path:
            try:
                from larmor.io import bruker

                ref = bruker.resolve(self.source_path)
                if (ref.expno / "ser").exists():
                    expno = str(ref.expno)
            except Exception:
                expno = None
        SatrecDialog(self, expno).exec()

    def open_redor(self):
        from larmor.desktop.tool_dialogs import RedorDialog

        expno = self.source_path if (self.source_path and
                                     Path(self.source_path).is_dir()) else None
        RedorDialog(self, expno).exec()

    def open_magres(self):
        from larmor.desktop.tool_dialogs import MagresDialog

        dlg = MagresDialog(self)
        if dlg.exec() and dlg.result_sites and self.recipe is not None:
            self.snapshot()
            for sd in dlg.result_sites:
                self.recipe["sites"].append(sd)
            self.on_structure_changed()
            self.statusBar().showMessage(
                f"added {len(dlg.result_sites)} site(s) from DFT tensors")

    def open_twod(self):
        from larmor.desktop.twod_dialog import TwoDDialog

        expno = self.source_path if (self.source_path and
                                     Path(self.source_path).is_dir()) else None
        TwoDDialog(self, expno).exec()

    # ------------------------------------------------------------- auto fit
    def run_auto_fit(self):
        if not self.recipe or not self.recipe["sites"]:
            self.statusBar().showMessage("add at least one line first")
            return
        from PySide6.QtWidgets import QInputDialog

        n, ok = QInputDialog.getInt(self, "Auto Fit",
                                    "number of random restarts:", 12, 2, 100)
        if not ok:
            return
        self.snapshot()
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        from larmor import autofit
        from larmor.recipe import Recipe

        self.statusBar().showMessage(f"auto fit: {n} restarts…")
        QApplication.processEvents()
        try:
            r = Recipe.from_dict(self.recipe)
            res = autofit.auto_fit(r, self.exp_ppm, self.exp_amp,
                                   window_ppm=(max(x0, x1), min(x0, x1)),
                                   n_starts=n)
        except Exception as exc:
            QMessageBox.warning(self, "Auto Fit failed", str(exc))
            return
        self.recipe = r.to_dict()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        self.lines_table.set_chi2(f"RMSD {res.best_rmsd:.4f}")
        self.request_simulation()
        self.run_quantify(show=False)
        self.statusBar().showMessage(res.summary)

    def run_errors_analysis(self):
        if not self.recipe or not self.recipe["sites"]:
            return
        from larmor.desktop.tool_dialogs import ErrorsDialog

        ErrorsDialog(self, self.recipe, self.exp_ppm, self.exp_amp,
                     self.view.current_xrange()).exec()

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
    """Load any supported source (shared with the CLI and figure studio)."""
    from larmor.loader import load_any

    return load_any(path)


def _nmrdata_to_data2d(data):
    """Convert a frequency-domain 2D NMRData into a twod.Data2D for display."""
    from larmor.twod import Data2D

    f1, f2 = data.axes
    h = data.hyper or {}
    d = Data2D(f2_ppm=np.asarray(f2.values), f1_ppm=np.asarray(f1.values),
               z=np.asarray(data.data, float), nucleus=data.nucleus,
               larmor_MHz=data.meta.get("larmor_MHz", 0.0),
               source=data.source,
               ri=h.get("ri"), ir=h.get("ir"), ii=h.get("ii"))
    d.notes = list(data.warnings)
    if data.is_pseudo2d:
        d.notes.append("pseudo-2D (arrayed)")
    return d


def asset_path(name: str) -> str:
    """Locate a bundled asset, whether running from source or a frozen exe."""
    here = Path(__file__).resolve()
    for base in (here.parent.parent.parent / "assets",     # repo/assets
                 here.parent.parent / "assets",            # larmor/assets
                 Path(getattr(sys, "_MEIPASS", "")) / "assets"):  # PyInstaller
        p = base / name
        if p.is_file():
            return str(p)
    return ""


def main() -> int:
    import time

    import pyqtgraph as pg
    from PySide6.QtGui import QFont, QIcon, QPixmap
    from PySide6.QtWidgets import QSplashScreen

    pg.setConfigOptions(antialias=True, background="#fcfdfc", foreground="#37424a")
    app = QApplication(sys.argv)
    app.setApplicationName("LARMOR")
    app.setStyle("Fusion")            # deterministic rendering on any OS theme
    app.setPalette(_light_palette())
    for family in ("Segoe UI", "Inter", "Roboto", "Helvetica Neue", "Arial"):
        f = QFont(family, 9)
        if f.exactMatch() or family == "Arial":
            app.setFont(f)
            break
    app.setStyleSheet(APP_STYLE)

    icon = asset_path("larmor_logo.png")
    if icon:
        app.setWindowIcon(QIcon(icon))

    # splash goes up FIRST, heavy imports/build happen behind it (like PRISM)
    splash = None
    shown_at = 0.0
    splash_png = asset_path("larmor_splash.png")
    if splash_png:
        splash = QSplashScreen(QPixmap(splash_png))
        splash.show()
        shown_at = time.time()
        app.processEvents()

    win = MainWindow()
    if icon:
        win.setWindowIcon(QIcon(icon))
    win.show()
    if splash is not None:
        while time.time() - shown_at < 2.0:      # keep the logo up briefly
            app.processEvents()
            time.sleep(0.02)
        splash.finish(win)
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
