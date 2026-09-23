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

Layout (G5). This module is the facade: ``MainWindow`` holds ``__init__``
(window state and build order), ``keyPressEvent`` and ``closeEvent``; every
other method is defined on one mixin in ``larmor/desktop/mw_<block>.py`` --
menus, chrome, files, session, overlays, editing, sidebands, fitting,
processing, cofit, tools -- and the QThreads live in
``larmor/desktop/workers.py``. A new window method goes in the mixin that
owns the state it touches, never here. Mixins define no ``__init__``, no
``Signal`` and no Qt event override, never import this module, and keep
pairwise-disjoint method names; ``tests/test_app_split.py`` enforces the
rules, freezes the member surface and pins the menu tree. Names other code
imports from here (``_load_any``, ``TUTORIALS``, the workers) are re-exported
through ``__all__``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import QFileSystemWatcher, QSettings, Qt, QTimer, Signal
from PySide6.QtWidgets import (QApplication, QLabel, QMainWindow, QVBoxLayout,
                               QWidget)

from larmor.desktop import theme
from larmor.desktop.plot import SpectrumView
from larmor.desktop.fithealth_strip import FitHealthStrip
from larmor.desktop.workers import (_emit_progress, humanize_error, _fit_tol,
                                    FitWorker, Fit2DWorker, KernelWarmWorker,
                                    SimWorker)
from larmor.desktop.mw_files import _FilesMixin, _load_any
from larmor.desktop.mw_overlays import _OverlaysMixin
from larmor.desktop.mw_editing import _EditingMixin
from larmor.desktop.mw_sidebands import _SidebandsMixin
from larmor.desktop.mw_session import _SessionMixin
from larmor.desktop.mw_processing import _ProcessingMixin
from larmor.desktop.mw_fitting import _FittingMixin
from larmor.desktop.mw_cofit import _CofitMixin
from larmor.desktop.mw_tools import _ToolsMixin
from larmor.desktop.mw_chrome import _ChromeMixin
from larmor.desktop.mw_menus import _MenusMixin, TUTORIALS

__all__ = [
    'MainWindow', 'main', 'asset_path', 'TUTORIALS', 'FitWorker',
    'Fit2DWorker', 'KernelWarmWorker', 'SimWorker', '_emit_progress',
    'humanize_error', '_fit_tol', '_load_any']


class ClickableLabel(QLabel):
    """A QLabel that emits doubleClicked (used for the experiment strip)."""
    doubleClicked = Signal()

    def mouseDoubleClickEvent(self, ev):
        self.doubleClicked.emit()
        super().mouseDoubleClickEvent(ev)


class MainWindow(_MenusMixin, _ChromeMixin, _FilesMixin, _SessionMixin,
                 _OverlaysMixin, _EditingMixin, _SidebandsMixin, _FittingMixin,
                 _ProcessingMixin, _CofitMixin, _ToolsMixin, QMainWindow):
    def __init__(self):
        super().__init__()
        # kill pyqtgraph's crash-prone native export dialog app-wide (LARMOR has
        # its own exporter on every plot's right-click menu)
        from larmor.desktop.plot_menu import disable_native_export_globally
        disable_native_export_globally()
        self.setWindowTitle("LARMOR")
        # Fit the initial window to the screen — a fixed 1440×900 overflows small
        # laptops, pushing corners (and the resize grips) off-screen.
        self._fit_to_screen(1440, 900)

        self.source_path: str | None = None
        # File > Watch the source file: reload on change, debounced (a
        # spectrometer writes 1r in several steps; TopSpin replaces the file,
        # which drops a per-file watch, so the parent folder is watched too)
        self._watcher = QFileSystemWatcher(self)
        self._watcher.fileChanged.connect(self._watched_changed)
        self._watcher.directoryChanged.connect(self._watched_changed)
        self._watch_timer = QTimer(self)
        self._watch_timer.setSingleShot(True)
        self._watch_timer.setInterval(600)
        self._watch_timer.timeout.connect(self._reload_watched)
        self._watch_path: str | None = None
        self.recipe: dict | None = None
        self.exp_ppm = np.array([])
        self.exp_amp = np.array([])
        #: unprocessed workbench spectrum the pipeline is (re)applied from, so
        #: live processing reflects ABSOLUTE settings instead of compounding
        self._proc_base: tuple[np.ndarray, np.ndarray] | None = None
        #: complex frequency-domain Spectrum1D from the last apply_processing
        #: (the imag / |S| display channels); None until a pipeline ran
        self._proc_spec = None
        #: the windowed, zero-filled FID captured just before the pipeline's
        #: LAST ft -- what the transform sees (the FID display); None when the
        #: chain has no ft
        self._proc_fid = None
        self._proc_apply_count = 0     # completed applies (a bail-out detector)
        self.hidden: set[int] = set()
        self.undo_stack: list[str] = []
        self.redo_stack: list[str] = []
        self._sim_worker: SimWorker | None = None
        self._sim_pending = False
        self._fit_worker: FitWorker | None = None
        self._last_quant = None
        self._paddle_live = False   # true while a paddle is being dragged
        # drag-to-phase (F2): a gesture in progress, the panel's p0/p1 when
        # it began, and the pivot fraction the last apply used (so a pivot
        # move while phasing can re-express p0 instead of jumping)
        self._phase_live = False
        self._phase_start = (0.0, 0.0)
        self._phase_pivot_frac_last = 0.5
        self._last_model = None     # (x, total) of the latest simulation
        self._first_sim = False     # autoscale Y once the first model arrives
        # fit health (see _health_from_result / _health_live): the verdict on
        # screen, the verdict of the last fit it derives from, and that fit's
        # lmfit result (Parameter correlations) -- all dropped together by
        # _health_reset whenever the active 1D document changes
        self._last_lmfit = None
        self._health = None         # fithealth.Health currently shown (fit or live)
        self._health_fit = None     # fithealth.Health of the last fit

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
        self._build_cofit_page()                     # index 2: co-fit split
        # the stack sits in a container so the fit-health strip can run
        # full-width directly under the spectrum (outside the plot, so it never
        # reaches a figure export); every page check keeps using
        # central_stack.currentWidget(), nothing reads centralWidget()
        container = QWidget()
        lay = QVBoxLayout(container)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.central_stack, 1)
        self.health_strip = FitHealthStrip()
        lay.addWidget(self.health_strip)
        self.setCentralWidget(container)
        self._build_progress()

        self.view.add_requested.connect(self.add_site_at)
        self.view.exit_add_mode.connect(lambda: self._set_add_mode(None))
        self.view.paddle_moved.connect(self.on_paddle_moved)
        self.view.paddle_released.connect(self.on_paddle_released)
        self.view.phase_dragged.connect(self.on_phase_dragged)
        self.view.phase_drag_released.connect(self.on_phase_drag_released)
        self.view.pivot_moved.connect(self._on_pivot_moved)
        self.view.file_dropped.connect(self.load_source)
        self.view.cursor_moved.connect(
            lambda x, y: self.pos_label.setText(
                f"x: {self._format_x(x)}   y: {y:.4g}"))
        self.view.calibrate_picked.connect(self.on_calibrate_picked)
        self.view.measure_changed.connect(self.on_measure_changed)
        # every 'a real spectrum arrived' site (loads, workspace switch, undo,
        # baseline tools, calibrate, SR, ...) re-validates the FID / channel
        # display through this one hook -- no per-site edits
        self.view.experiment_set.connect(self._on_view_experiment_set)
        self.view2d.slice_to_fit.connect(self._trace_to_workbench)

        self._build_menus()
        self._build_toolbar()
        self._build_sidebar()
        self._build_explorer_dock()
        self._build_datasets_dock()
        self._build_workspaces_dock()
        self._build_bottom_docks()
        self._build_right_dock()
        self._build_panels_menu()

        # fit-health strip: every click is navigation, never a recipe edit
        self.health_strip.pill.clicked.connect(self.show_fit_health)
        self.health_strip.open_report.connect(
            lambda: (self.results_dock.show(), self.results_dock.raise_()))
        self.health_strip.show_residual.connect(self._health_show_residual)
        self.health_strip.focus_param.connect(self._health_focus_param)
        self.health_strip.open_correlations.connect(self.show_correlations)
        self.health_strip.open_errors.connect(self.run_errors_analysis)
        self.health_strip.open_help.connect(
            lambda: self._open_manual("spectra-1d",
                                      "1D spectra — processing & fitting"))
        self.central_stack.currentChanged.connect(self._health_show)
        self.health_strip.setVisible(False)      # shown once a 1D recipe has sites

        self.exp_label = ClickableLabel("")
        self.exp_label.setStyleSheet(
            f"color: {theme.active().accent}; font-weight: 600;")
        self.exp_label.setToolTip("double-click to edit the experiment "
                                  "parameters (nucleus, Larmor, νrot)")
        self.exp_label.setCursor(Qt.PointingHandCursor)
        self.exp_label.doubleClicked.connect(self.edit_experiment)
        self.statusBar().addPermanentWidget(self.exp_label)
        self.pos_label = QLabel("")
        self.statusBar().addPermanentWidget(self.pos_label)
        # a red MAS-rate warning at the bottom-right when the spin rate had to
        # be guessed or the sources disagreed (see _update_mas_label); double-
        # clicking it opens the experiment parameters, same as the exp label
        self.mas_label = ClickableLabel("")
        self.mas_label.setStyleSheet(
            "background:#c0392b; color:white; font-weight:600; "
            "padding:1px 8px; border-radius:3px;")
        self.mas_label.setVisible(False)
        self.mas_label.setCursor(Qt.PointingHandCursor)
        self.mas_label.doubleClicked.connect(self.edit_experiment)
        self.statusBar().addPermanentWidget(self.mas_label)
        self.statusBar().showMessage(
            "File > Open… (dmfit .fxmla / LARMOR recipe) or Open EXPNO…")

        self._sim_timer = QTimer(self)
        self._sim_timer.setSingleShot(True)
        self._sim_timer.setInterval(120)
        self._sim_timer.timeout.connect(self._simulate_now)
        # a busy cue if a simulation runs long (the first Czjzek/Amorphous fit
        # builds a kernel — seconds; without this the app looks frozen)
        self._busy = False
        self._busy_timer = QTimer(self)
        self._busy_timer.setSingleShot(True)
        self._busy_timer.setInterval(450)
        self._busy_timer.timeout.connect(self._sim_busy_on)
        # crash-safe autosave: snapshot the session every few minutes (on top of
        # the per-action save) so a crash never costs more than a couple of minutes
        self._autosave_timer = QTimer(self)
        self._autosave_timer.setInterval(180000)          # 3 minutes
        self._autosave_timer.timeout.connect(self._autosave_tick)
        self._autosave_timer.start()
        # debounced per-action session write (see _persist_session)
        self._session_timer = QTimer(self)
        self._session_timer.setSingleShot(True)
        self._session_timer.setInterval(5000)
        self._session_timer.timeout.connect(self._flush_session)
        self._data2d = None            # the 2D dataset currently on the map
        self._fit2d_worker = None
        self.view2d.add_requested.connect(self.add_site_2d)
        self.view2d.load_1d_for_projection.connect(self.load_projection_1d)
        self.view2d.overlay_1d_request.connect(self.overlay_1d_on_2d)
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

    def keyPressEvent(self, ev):
        if ev.key() == Qt.Key_Escape:
            banner = getattr(self, "ssb_banner", None)
            if (banner is not None and not banner.isHidden()
                    and self.central_stack.currentWidget() is self.view):
                self._dismiss_sideband_offer()   # the transient offer first
                return
            self.proc_panel.btnDrag.setChecked(False)   # leave drag-to-phase
            self._set_add_mode(None)
        super().keyPressEvent(ev)

    def closeEvent(self, ev):
        self._flush_session()
        # Join the worker threads before the window (their Python owner) goes
        # away: a QThread destroyed while running aborts the whole process
        # ("QThread: Destroyed while thread is still running") -- seen as test
        # runs that printed every dot and died before the summary, or hung,
        # right after a fixture window that had started a kernel pre-build.
        for name in ("_warm_worker", "_sim_worker", "_fit_worker", "_fit2d_worker"):
            w = getattr(self, name, None)
            if w is None or not hasattr(w, "isRunning"):
                continue
            try:
                if w.isRunning():
                    if hasattr(w, "request_stop"):
                        w.request_stop("stop")
                    w.wait(30_000)
            except RuntimeError:            # already deleted on the C++ side
                pass
        try:
            from larmor.parallel import shutdown_shared_pool
            shutdown_shared_pool()
        except Exception:
            pass
        super().closeEvent(ev)


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


def _crash_log_path() -> str:
    try:
        return os.path.join(os.path.expanduser("~"), "larmor_crash.log")
    except Exception:
        return "larmor_crash.log"


def _install_faulthandler() -> str:
    """Make a C-level crash (segfault) leave a traceback instead of a silent
    "LARMOR stopped". A console run dumps to stderr; the frozen, windowed
    exe has ``sys.stderr is None`` -- ``faulthandler.enable()`` then raises
    RuntimeError, which killed every frozen build at start-up until the
    0.13.0 smoke test caught it -- so it dumps to the crash log instead.
    faulthandler writes to ONE target, hence the either/or. Returns where."""
    import faulthandler

    try:
        if sys.stderr is not None:
            faulthandler.enable()
            return "stderr"
        faulthandler.enable(open(_crash_log_path(), "w"))
        return _crash_log_path()
    except Exception:
        return "unavailable"


def main() -> int:
    import time

    _install_faulthandler()
    _log = _crash_log_path()                 # named in the start-up error text

    import pyqtgraph as pg
    from PySide6.QtGui import QFont, QIcon, QPixmap
    from PySide6.QtWidgets import QSplashScreen

    pg.setConfigOptions(antialias=True)
    app = QApplication(sys.argv)
    app.setApplicationName("LARMOR")
    app.setStyle("Fusion")            # deterministic rendering on any OS theme
    pt = int(QSettings("LARMOR", "app").value("fontPt", 9) or 9)
    for family in ("Segoe UI", "Inter", "Roboto", "Helvetica Neue", "Arial"):
        f = QFont(family, pt)
        if f.exactMatch() or family == "Arial":
            app.setFont(f)
            break
    # apply the saved colour theme (palette + stylesheet + pyqtgraph config).
    # A hidden "aesthetic" override (View ▸ Theme ▸ More styles…) takes
    # precedence when set, so the app opens back into whichever style —
    # normal or aesthetic — was live when it last closed.
    settings = QSettings("LARMOR", "app")
    override = settings.value("appearanceOverride", "")
    if override and override in theme.AESTHETIC_THEMES:
        theme.apply(app, override)
    else:
        saved = settings.value("theme", theme.DEFAULT)
        theme.apply(app, saved if saved in theme.THEMES else theme.DEFAULT)

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

    try:
        win = MainWindow()
        if icon:
            win.setWindowIcon(QIcon(icon))
        win.show()
        win._maybe_show_welcome()      # first-run canvas hint (returning users: no-op)
        if splash is not None:
            while time.time() - shown_at < 2.0:      # keep the logo up briefly
                app.processEvents()
                time.sleep(0.02)
            splash.finish(win)
        return app.exec()
    except Exception:
        import traceback
        traceback.print_exc()
        print(f"\nLARMOR hit an error during start-up. If it was a hard crash, "
              f"a diagnostic was written to:\n  {_log}\n"
              "Please send that file (or this text) to the developer.\n",
              file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
