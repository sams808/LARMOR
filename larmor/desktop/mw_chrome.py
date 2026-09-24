"""Main-window mixin: window furniture.

The five dock builders and ``_fit_to_screen``, the status-bar progress bar
with its Stop / Cancel buttons, the experiment / MAS status labels, the view
operations (zoom, autoscale, residual / components / labels / paddles
toggles) and copy / save of the current plot.

Owned state: the docks and their panels (``explorer``, ``ws_panel``,
``datasets_panel``, ``lines_table`` / ``lines_stack``, ``qtable`` /
``report`` / ``results_summary``, ``proc_panel``), the ``progress`` bar with
``btnStopFit`` / ``btnCancelFit`` and ``_prog_label`` / ``_prog_prev_sdev``;
``exp_label`` / ``mas_label`` / ``pos_label`` are created in ``__init__``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings, QTimer, Qt
from PySide6.QtWidgets import (QApplication, QDockWidget, QHBoxLayout, QLabel,
                               QMessageBox, QPlainTextEdit, QProgressBar,
                               QPushButton, QSizePolicy, QStackedWidget,
                               QTableWidget, QVBoxLayout, QWidget)

from larmor.desktop import theme
from larmor.desktop.panels import ProcessingPanel
from larmor.desktop.table import LinesTable


class _ChromeMixin:
    """Docks, status bar, progress bar, view operations."""

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
        from larmor.desktop.export_dialog import export_pyqtgraph
        item = (self.view.getPlotItem() if w is self.view else self.view2d.p_main)
        try:
            path = export_pyqtgraph(self, item, "plot")
        except Exception as exc:
            QMessageBox.warning(self, "Save image", f"export failed: {exc}")
            return
        if path:
            self.statusBar().showMessage(f"plot saved to {Path(path).name}")

    def _build_explorer_dock(self):
        from larmor.desktop.explorer import ExplorerPanel

        self.explorer_dock = QDockWidget("Explorer", self)
        self.explorer_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                       QDockWidget.DockWidgetClosable)
        self.explorer = ExplorerPanel()
        self._proj_pick_axis = None      # HMQC: awaiting an Explorer pick
        self._pending_apply_model = None  # a model awaiting an Explorer data pick
        self.explorer.open_requested.connect(self._explorer_open)
        self.explorer.batch_requested.connect(self.run_batch_fit)
        self.explorer.inventory_requested.connect(self.open_session_inventory)
        self.explorer.overlay_requested.connect(self.add_overlay_path)
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
        self.ws_panel.activate.connect(self.open_workspace_entry)
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
        self.datasets_panel.color_changed.connect(self.overlay_set_color)
        self.datasets_panel.offset_changed.connect(lambda _: self._refresh_overlays())
        self.datasets_panel.compare_requested.connect(self.compare_overlays)
        self.datasets_panel.match_changed.connect(lambda _: self._refresh_overlays())
        self.datasets_dock.setWidget(self.datasets_panel)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.datasets_dock)
        # stack the left docks as tabs so they share one footprint (kinder on
        # laptop screens); Explorer is the default front tab. Workspaces joins
        # the group when it is built (next).
        self.tabifyDockWidget(self.explorer_dock, self.datasets_dock)
        self.explorer_dock.raise_()

    # ------------------------------------------------------------- docks
    def _build_bottom_docks(self):
        self.lines_dock = QDockWidget("Fit parameters", self)
        self.lines_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                    QDockWidget.DockWidgetClosable)
        self.lines_table = LinesTable()
        self.lines_table.edited.connect(self.on_params_changed)
        self.lines_table.constraint_edited.connect(self.on_structure_changed)
        self.lines_table.family_edited.connect(self.on_family_changed)
        self.lines_table.structure.connect(self.on_site_structure)
        self.lines_table.compute.connect(self.request_simulation)
        self.lines_table.fit.connect(self.run_fit)
        # a stack so co-fit can swap in its split 1D|2D tables in the same dock
        self.lines_stack = QStackedWidget()
        self.lines_stack.addWidget(self.lines_table)     # index 0: normal fit
        self.lines_stack.addWidget(self.cofit_tables)    # index 1: co-fit (built earlier)
        self.lines_dock.setWidget(self.lines_stack)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.lines_dock)

        self.results_dock = QDockWidget("Report", self)
        self.results_dock.setFeatures(QDockWidget.DockWidgetMovable |
                                      QDockWidget.DockWidgetClosable)
        w = QWidget()
        v = QVBoxLayout(w)
        head = QHBoxLayout()
        self.results_summary = QLabel("")
        self.results_summary.setStyleSheet("font-weight: 600;")
        # wrap + shrinkable so a long fit summary never forces the Report dock —
        # and thus the whole window — wider than the screen (it sits full-width in
        # a bottom dock; a non-wrapping label would push its min to the text width)
        self.results_summary.setWordWrap(True)
        self.results_summary.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        head.addWidget(self.results_summary, 1)
        btnCsv = QPushButton("Copy CSV")
        btnCsv.clicked.connect(self.copy_csv)
        head.addWidget(btnCsv)
        btnTex = QPushButton("Copy LaTeX")
        btnTex.setToolTip("copy a LaTeX results table (sites × parameters + "
                          "population %) to the clipboard")
        btnTex.clicked.connect(self.copy_latex)
        head.addWidget(btnTex)
        btnMeth = QPushButton("Copy methods")
        btnMeth.setToolTip("copy a paper-ready Experimental paragraph — spectrometer and "
                           "field, probe, MAS rate, pulse and flip angle, recycle delay, "
                           "scans, referencing, processing, the fit and the software "
                           "versions (a CSV/dmfit source gets the fit sentence only)")
        btnMeth.clicked.connect(self.copy_methods)
        head.addWidget(btnMeth)
        btnBundle = QPushButton("Publication bundle…")
        btnBundle.setToolTip("write a figure (png/pdf/svg) + LaTeX table + CSV + "
                             "methods sentence + report.md for this fit, in one click")
        btnBundle.clicked.connect(self.export_publication_bundle)
        head.addWidget(btnBundle)
        v.addLayout(head)
        self.qtable = QTableWidget(0, 5)
        self.qtable.setHorizontalHeaderLabels(
            ["line", "position (ppm)", "integral", "fraction (%)",
             "outside window (%)"])
        self.qtable.horizontalHeaderItem(4).setToolTip(
            "share of the line's simulated area outside the integration "
            "window — biases its population low; the fit-health strip flags "
            "> 2 % and a click on that chip widens the window")
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
        self.proc_panel.twopoint_mode.connect(self._twopoint_mode)
        self.proc_panel.twopoint_apply.connect(self.apply_twopoint_bg)
        self.proc_panel.twopoint_clear.connect(self.view.clear_baseline)
        self.proc_panel.phase_drag_mode.connect(self._phase_drag_mode)
        self.proc_panel.view_changed.connect(self._on_proc_view_changed)
        self.proc_dock.setWidget(self.proc_panel)
        self.addDockWidget(Qt.RightDockWidgetArea, self.proc_dock)
        self.proc_dock.hide()
        # show the TopSpin-style phase pivot line whenever the panel is open
        self.proc_dock.visibilityChanged.connect(self.view.show_phase_pivot)

    def _fit_to_screen(self, want_w: int, want_h: int):
        """Size the window to at most the available screen (minus a margin) and
        centre it, so it never opens partly off-screen on a small display."""
        from PySide6.QtWidgets import QApplication
        scr = (self.screen() or QApplication.primaryScreen())
        if scr is None:
            self.resize(want_w, want_h)
            return
        g = scr.availableGeometry()
        w = min(want_w, g.width() - 40)
        h = min(want_h, g.height() - 60)
        self.resize(w, h)
        self.move(g.left() + (g.width() - w) // 2,
                  g.top() + (g.height() - h) // 2)

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
        origin = self._mas_origin_sentence()
        tip = ((f"{mas} — {origin} · " if origin else "")
               + "double-click to edit the experiment parameters")
        acq = self.recipe.get("acquisition") or {}
        if acq:
            from larmor import acquisition

            tip += "\n" + "\n".join(acquisition.summary_lines(
                acq, self.recipe.get("software"), self.recipe.get("spin_rate_Hz"),
                self.recipe.get("mas_uncertain")))
        self.exp_label.setToolTip(tip)
        self._update_mas_label()
        self._apply_axis_unit()          # SFO may have changed with the dataset

    _MAS_SOURCE_NAMES = {"acqus": "acqus MASR", "title": "title",
                         "booking": "booking sidecar"}

    def _mas_origin_sentence(self) -> str:
        """Where the recipe's νrot came from, from provenance['mas_rate']
        (larmor.masrate); '' when the spectrum carries no block (CSV…)."""
        block = ((self.recipe or {}).get("provenance") or {}).get("mas_rate")
        if not block:
            return ""
        from larmor import masrate

        names = self._MAS_SOURCE_NAMES
        present = [k for k in ("acqus", "title", "booking")
                   if block.get(f"{k}_Hz") is not None]
        src = str(block.get("source") or "")
        if src == "confirmed" and block.get("confirmed"):
            key = (block.get("session") or "", block.get("rotor") or "",
                   block.get("nucleus") or "")
            return (f"confirmed on {str(block['confirmed'])[:10]} for "
                    f"{masrate.describe_key(key)}")
        if src == "all":
            return "all sources agree"
        if "+" in src:
            pair = src.split("+")
            text = f"{names[pair[0]]} and {names[pair[1]]} agree"
            odd = [k for k in present if k not in pair]
            return text + (f"; {names[odd[0]]} outvoted" if odd else "")
        if src == "highest":
            return ("highest of disagreeing sources"
                    + ("" if block.get("uncertain")
                       else ", confirmed in Experiment parameters"))
        if src == "fallback":
            return "assumed (no source found)"
        if src == "static":
            return "static from acqus / title / pulse program"
        if src in names:
            return f"from the {names[src]} only"
        return ""

    def _update_mas_label(self):
        """Red bottom-right MAS indicator when the spin rate was guessed or the
        sources disagreed. Cleared once the user confirms it (Experiment dialog)."""
        rate = (self.recipe or {}).get("spin_rate_Hz", 0.0) or 0.0
        if self.recipe and self.recipe.get("mas_uncertain"):
            self.mas_label.setText("⚠ assumed static — check!" if rate == 0
                                   else f"⚠ MAS {rate:.0f} Hz — check!")
            block = (self.recipe.get("provenance") or {}).get("mas_rate") or {}
            present = [k for k in ("acqus", "title", "booking")
                       if block.get(f"{k}_Hz") is not None]
            if present and rate:
                from larmor import masrate

                listed = " · ".join(
                    f"{k} {masrate.format_hz(float(block[f'{k}_Hz']))} Hz"
                    for k in present)
                src = str(block.get("source") or "")
                if src == "highest":
                    why = f"{listed} disagree; LARMOR took the highest. "
                elif len(present) == 1:
                    why = (f"{listed} is the only source"
                           + (" and the title's unit had to be repaired"
                              if block.get("title_repaired") else "") + ". ")
                else:
                    why = f"{listed}; a static hint flagged the rate. "
                if block.get("title_note"):
                    why += str(block["title_note"]) + ". "
                self.mas_label.setToolTip(
                    why + "Double-click here (or the experiment strip) to "
                    "compare them, measure from the spinning sidebands, and "
                    "confirm once for this session.")
            else:
                self.mas_label.setToolTip(
                    ("The dataset looks static (MASR recorded as 0, or a "
                     "wideline pulse program) but no second source confirms "
                     "it, so LARMOR assumed static (0 Hz). "
                     if rate == 0 else
                     "The MAS rate was missing or the acqus/title sources "
                     "disagreed, so LARMOR guessed (highest found, or 35714 Hz). ")
                    + "Double-click here (or the experiment strip) to open "
                      "Experiment parameters and clear this warning.")
            self.mas_label.setVisible(True)
        else:
            self.mas_label.setVisible(False)

    # ------------------------------------------------------------- view ops
    def zoom_full(self):
        self.view.getPlotItem().enableAutoRange()

    def zoom_sites(self):
        if self.view.domain == "time":       # a ppm window on the ms axis
            return
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
        if self.view.domain == "time":
            # the FID: its ppm-window selection would pick experiment points
            # whose ppm happens to fall in the ms range
            self.view.getPlotItem().getViewBox().enableAutoRange(y=True)
            return
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

    def _toggle_labels(self, on):
        QSettings("LARMOR", "app").setValue("compLabels", bool(on))
        self.view.set_show_labels(on)

    def _toggle_paddles(self, on):
        self.view.show_paddles(on)

    # ------------------------------------------------------------- progress
    def _build_progress(self):
        """A live fit-progress bar in the status bar: iteration count and the
        current residual RMS, so a long fit shows what it is doing."""
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setFixedWidth(340)
        self.progress.setTextVisible(True)
        self._style_progress()
        self.progress.setVisible(False)
        self.statusBar().addPermanentWidget(self.progress)
        self._prog_label = "fitting"
        # interrupt a running fit: keep the latest iteration, or revert entirely
        self._active_fit_worker = None
        self.btnStopFit = QPushButton("⏹ Stop (keep)")
        self.btnStopFit.setToolTip("stop the optimiser now and keep the latest "
                                   "iteration's parameter values")
        self.btnStopFit.clicked.connect(lambda: self._interrupt_fit("stop"))
        self.btnCancelFit = QPushButton("✖ Cancel (revert)")
        self.btnCancelFit.setToolTip("abort the fit and restore the parameters "
                                     "to what they were before it started")
        self.btnCancelFit.clicked.connect(lambda: self._interrupt_fit("cancel"))
        for b in (self.btnStopFit, self.btnCancelFit):
            b.setVisible(False)
            self.statusBar().addPermanentWidget(b)

    def _show_fit_buttons(self, on: bool):
        self.btnStopFit.setVisible(on)
        self.btnCancelFit.setVisible(on)

    def _interrupt_fit(self, mode: str):
        w = self._active_fit_worker
        if w is not None and w.isRunning():
            w.request_stop(mode)
            self.btnStopFit.setEnabled(False)
            self.btnCancelFit.setEnabled(False)
            self.statusBar().showMessage(
                "stopping the fit — keeping the latest values…" if mode == "stop"
                else "cancelling the fit — reverting…")

    def _style_progress(self):
        t = theme.active()
        self.progress.setStyleSheet(f"""
            QProgressBar {{
                border: 1px solid {t.border}; border-radius: 7px;
                background: {t.base}; height: 16px;
                color: {t.text}; font-size: 11px; text-align: center;
            }}
            QProgressBar::chunk {{
                border-radius: 6px;
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 #16a3a3, stop:0.55 #2fbf8f, stop:1 #7fd66a);
            }}
        """)

    def _progress_start(self, label: str):
        self._prog_label = label
        self._prog_prev_sdev = None          # for the live Δσ% (dmfit varsdev%)
        # a fit's iteration count has no known total, so run the bar as a busy
        # (indeterminate) marquee — it never sticks near the end. The iteration
        # count and the residual stdev in the text give the real feedback.
        self.progress.setRange(0, 0)
        self.progress.setFormat(f"{label} — starting…")
        self.progress.setVisible(True)
        QApplication.processEvents()

    def _progress_tick(self, it: int, rms: float):
        # dmfit-style live read-out: the residual stdev and its % change per
        # iteration (varsdev%). Convergence stops the fit once |Δσ%| drops below
        # the completion threshold (Fit ▸ Fit settings).
        prev = getattr(self, "_prog_prev_sdev", None)
        if prev and prev > 0 and rms == rms:
            dpct = 100.0 * (rms - prev) / prev
            self.progress.setFormat(
                f"{self._prog_label} — iter {it} · sdev {rms:.5g} · "
                f"Δσ {dpct:+.3f}%")
        else:
            self.progress.setFormat(
                f"{self._prog_label} — iter {it} · sdev {rms:.5g}")
        self._prog_prev_sdev = rms

    def _progress_end(self, ok: bool = True):
        self.progress.setRange(0, 100)
        self.progress.setValue(100)
        self.progress.setFormat(("done" if ok else "stopped") + "  %p%")
        self._show_fit_buttons(False)
        self.btnStopFit.setEnabled(True)
        self.btnCancelFit.setEnabled(True)
        QTimer.singleShot(1200, lambda: self.progress.setVisible(False))
