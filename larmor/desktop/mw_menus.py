"""Main-window mixin: the menu bar, toolbar and sidebar.

``_build_menus`` / ``_add`` and everything that exists only to populate or
persist menu state: recent files and apply-recipe rebuilders, the toolbar
and sidebar, the Czjzek-display / panels / axis-unit / theme / text-size
submenus and their setters, ``_update_enabled``, the Help slots (manuals,
tutorials, command palette, About, More). ``TUTORIALS`` (the Help >
Tutorials entries) lives here and is re-exported by ``app.py``.

Owned state: every ``act*`` action, ``m_recent`` / ``m_apply_recipe``,
``_model_actions``, the toolbar / sidebar widgets.
"""
from __future__ import annotations

import os
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (QFileDialog, QHBoxLayout, QLabel, QMessageBox,
                               QPushButton, QToolBar, QVBoxLayout)

from larmor import models as model_registry
from larmor.desktop import theme
from larmor.recipe import Recipe

#: (file stem under docs/tutorials, menu title) -- the Help ▸ Tutorials entries.
#: tests/test_tutorials.py holds the same list and checks the files ship.
TUTORIALS = (
    ("01-first-fit-27Al-czjzek", "1 · A first fit — ²⁷Al Czjzek"),
    ("02-constraints", "2 · Constraints — fix, bound, link"),
    ("03-figures", "3 · Plotting studio figures"),
    ("04-batch-fitting", "4 · Batch fitting a composition series"),
    ("05-mqmas", "5 · MQMAS — 2D processing & fitting"),
    ("06-error-analysis", "6 · Error analysis — covariance, Monte-Carlo, χ² profile"),
    ("07-static-81Br-wcpmg", "7 · Static wideline — ⁸¹Br WURST-CPMG"),
)


class _MenusMixin:
    """Menu bar, toolbar, sidebar and the menu-state setters."""

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
        self._add(m_file, "Open &Varian / Agilent…  (.fid folder)",
                  self.open_varian)
        self._add(m_file, "O&verlay a spectrum…  (compare on top of the active one, "
                          "keeps the fit; Shift + drop does the same)",
                  self.add_overlay_dialog, "Ctrl+Shift+A")
        self.actWatch = self._add(
            m_file, "&Watch the source file  (auto-reload when it changes, "
                    "keep the fit)", self._toggle_watch, checkable=True)
        self.actWatch.setToolTip(
            "re-open the current spectrum whenever its file is rewritten -- "
            "e.g. while acquiring on the spectrometer -- keeping the fit "
            "model so it follows the growing signal")
        m_file.addSeparator()
        self._add(m_file, "Open pro&ject…  (spectra, 2D maps, figures, batch fits)",
                  self.open_project)
        self._add(m_file, "Save projec&t…  (whole session: spectra, 2D maps, "
                          "figures, batch fits)",
                  self.save_project, "Ctrl+Alt+S")
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
        # display projections of the pipeline result (F6): the FID the
        # transform sees, and the real / imaginary / |S| channel
        from PySide6.QtGui import QActionGroup

        self.actTimeDomain = self._add(
            m_proc, "&FID ⇄ spectrum  (time ↔ frequency, re-apodize)",
            self._toggle_time_domain, "Ctrl+T", checkable=True)
        self.actTimeDomain.setToolTip(
            "show the windowed FID the transform sees and re-apply the window "
            "functions live; press again to return to the spectrum")
        m_chan = m_proc.addMenu("Display &channel")
        chan_group = QActionGroup(self)
        chan_group.setExclusive(True)
        self.actChannel = {}
        for name, label in (
                ("real", "&Real"),
                ("imag", "&Imaginary  (inspect while phasing)"),
                ("magnitude", "&Magnitude  |S|  (display only — the pipeline "
                              "op is the panel checkbox)")):
            a = self._add(m_chan, label,
                          lambda _=False, n=name: self._set_channel(n),
                          checkable=True, checked=(name == "real"))
            chan_group.addAction(a)
            self.actChannel[name] = a
        m_chan.addSeparator()
        self._add(m_chan, "C&ycle channel  (real → imag → |S|)",
                  self._cycle_channel, "Ctrl+I")
        self._add(m_proc, "Processing s&teps…  (remove a step)",
                  self.edit_processing_steps)
        self._add(m_proc, "Autophase (ACME)",
                  lambda: self.apply_processing([{"op": "autophase"}], False))
        # non-checkable: the panel's Drag-to-phase button is the single
        # source of truth for the mode, this entry toggles it
        self._add(m_proc, "&Drag to phase  (← → p0, ↑ ↓ p1 about the pivot)",
                  self.start_phase_drag, "Ctrl+P")
        self._add(m_proc, "2-point background…  (pick two flat points)",
                  self.start_twopoint_bg)
        self._add(m_proc, "Baseline auto (order 3)",
                  lambda: self.apply_processing([{"op": "baseline", "order": 3}], False))
        self._add(m_proc, "Baseline iterative (dead-time; Yon 2020)…",
                  self.apply_iterbaseline)
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
        self._add(m_proc, "&WURST excitation profile…  (divide out the sweep)",
                  self.open_wurst_correct)

        m_dec = mb.addMenu("&Decomposition")
        # --- build the model ---
        self._add(m_dec, "&New fit (clear lines)", self.new_fit)
        m_models = m_dec.addMenu("Add &line  (pick a model)")
        self._model_actions = {}
        for m in model_registry.describe_all():
            if m["name"] == "spectrum":
                continue          # added via "Add background spectrum" below
            a = QAction(m["label"], self)
            a.setCheckable(True)
            a.setToolTip(m["description"])
            a.triggered.connect(
                lambda checked, name=m["name"]: self._set_add_mode(
                    name if checked else None))
            m_models.addAction(a)
            self._model_actions[m["name"]] = a
        self.m_apply = m_dec.addMenu("&Apply recipe  (a recent fit → this data)")
        self._rebuild_apply_recipe()
        self._add(m_dec, "Add a line at every &peak…  (auto peak-pick)",
                  self.autopick_lines)
        self._add(m_dec, "&Label lines from literature ranges  (Al[4], BO3, …)",
                  self.label_from_literature)
        self._add(m_dec, "Add spinning &sidebands…  (of a fitted line)",
                  self.add_sidebands)
        self.actDetectSsb = self._add(
            m_dec, "&Detect spinning sidebands  (find the ±νrot repeat)",
            self.detect_sidebands, "Ctrl+Shift+D")
        self.actDetectSsb.setToolTip(
            "look for the repeat of the spectrum at ±νrot (around the "
            "recorded rate, then 1–80 kHz) and offer a linked sideband "
            "manifold or a shifted copy of the spectrum in one click")
        self._add(m_dec, "Add f&unction line…  (y = f(x; a,b,c,d))",
                  self.add_function_line)
        self._add(m_dec, "Add background &spectrum…  (fit another spectrum)",
                  self.add_background_spectrum)
        self._add(m_dec, "Add a &copy of this spectrum…  (shifted — satellite "
                         "/ sideband manifold)",
                  self.add_current_spectrum_line)
        self._add(m_dec, "Add fit &zone", self.add_zone)
        self._add(m_dec, "Clear zones", self.clear_zones)
        m_dec.addSeparator()
        # --- run ---
        self._add(m_dec, "&Simulate  (recompute the model)",
                  self.request_simulation, "F9")
        self.actFit = self._add(m_dec, "&Fit", self.run_fit, "F5")
        self.actAuto = self._add(m_dec, "&Auto Fit (multi-start)…",
                                 self.run_auto_fit)
        m_dec.addSeparator()
        # --- analyze ---
        self.actQuant = self._add(m_dec, "&Report (quantify)", self.run_quantify, "F6")
        # F5 Fit / F6 Report / F7 Health form one cluster
        self.actHealth = self._add(m_dec, "Fit &health details…  (why this verdict)",
                                   self.show_fit_health, "F7")
        self.actErrors = self._add(m_dec, "&Errors Analysis (χ² profile)…",
                                   self.run_errors_analysis)
        # kept as attributes: the fit-health strip's menu lists these SAME
        # QAction objects (no duplicated labels or slots)
        self.actMC = self._add(m_dec, "Monte-&Carlo errors…  (synthetic-noise refits)",
                               self.run_monte_carlo)
        self.actCorr = self._add(m_dec, "Parameter correlations…  (from the last fit)",
                                 self.show_correlations)
        self._add(m_dec, "Compare with a saved fit…  (parameter diff)",
                  self.compare_with_saved_fit)
        self._add(m_dec, "Czjzek distribution P(C_Q)…  (what σ stands for)",
                  self.show_czjzek_dist)
        self.actChi2 = self._add(m_dec, "χ² map (parameter pair)…  (is the pair determined?)",
                                 self.show_chi2_map)
        m_dec.addSeparator()
        # --- advanced / configuration (rarely touched) ---
        m_adv = m_dec.addMenu("Ad&vanced")
        self._add(m_adv, "Co-&fit datasets…  (shared model, 1D + MQMAS)",
                  self.open_cofit)
        self._add(m_adv, "Computing &parameters…  (kernel resolution)",
                  self.edit_computing_params)
        self._add(m_adv, "Fit completion &threshold…  (Δσ % to stop)",
                  self.edit_fit_tol)
        m_adv.addSeparator()
        self._add(m_adv, "Restrict around current values…  (glass protocol, "
                  "Edén 2023)", self.restrict_glass_protocol)
        self._add(m_adv, "Save &constraints as…  (reusable link/bound set)",
                  self.save_constraint_set)
        self._add(m_adv, "Apply saved co&nstraints…", self.apply_constraint_set)
        self._add(m_adv, "MQMAS F1 &reference…  (isotropic-axis align)",
                  self.edit_mqmas_f1_ref)
        self._add(m_adv, "Predict at another &field…  (what at X T?)",
                  self.predict_at_field)

        m_view = mb.addMenu("&View")
        self.m_view = m_view                 # panels submenu filled once docks exist
        self.actResid = self._add(m_view, "Residual", self._toggle_resid,
                                  checkable=True, checked=True)
        self.actComp = self._add(m_view, "Components", self._toggle_comp,
                                 checkable=True, checked=True)
        self.actOverlaysVisible = self._add(
            m_view, "O&verlays  (the compared spectra)", self._toggle_overlays_visible,
            "Ctrl+Shift+V", checkable=True, checked=True)
        self.actOverlaysVisible.setToolTip(
            "show or hide every compared spectrum at once; add one with "
            "File > Overlay a spectrum, Shift + drop on the plot, or the "
            "Explorer's right-click")
        self._add(m_view, "Clear overlays", self.clear_overlays)
        self.actLabels = self._add(
            m_view, "Component &labels  (pin names on the plot)",
            self._toggle_labels, checkable=True,
            checked=bool(QSettings("LARMOR", "app").value(
                "compLabels", False, type=bool)))
        self.actLabels.setToolTip(
            "write each component's letter and name at its maximum; when "
            "off, hovering a component still shows its name")
        self.view.set_show_labels(self.actLabels.isChecked())
        self.actPaddles = self._add(m_view, "Show paddles", self._toggle_paddles,
                                    checkable=True, checked=True)
        self.actAnimateFit = QAction("Animate fits", self)
        self.actAnimateFit.setCheckable(True)
        self.actAnimateFit.setToolTip("draw the model curve as it converges during "
                                      "a 1D fit (a fading trail shows the last few "
                                      "iterations) — watch convergence or divergence")
        self.actAnimateFit.setChecked(bool(QSettings("LARMOR", "app").value(
            "animateFit", True, type=bool)))
        self.actAnimateFit.toggled.connect(
            lambda on: QSettings("LARMOR", "app").setValue("animateFit", bool(on)))
        m_view.addAction(self.actAnimateFit)
        self.actScrollNudge = QAction("Scroll &nudges fit values", self)
        self.actScrollNudge.setCheckable(True)
        self.actScrollNudge.setToolTip("when on, scrolling over a parameter cell "
                                       "nudges its value (off by default so a "
                                       "stray scroll never changes a fit)")
        self.actScrollNudge.setChecked(bool(QSettings("LARMOR", "app").value(
            "scrollNudge", False, type=bool)))
        # setChecked() above fires no signal (connect comes next), so push the
        # saved state into the table module DIRECTLY — otherwise a remembered
        # "on" showed checked here but didn't actually nudge until re-toggled
        from larmor.desktop import table as _table
        _table.set_scroll_nudge(self.actScrollNudge.isChecked())
        self.actScrollNudge.toggled.connect(self._toggle_scroll_nudge)
        m_view.addAction(self.actScrollNudge)
        self.actRefRanges = QAction("&Literature shift ranges  (Edén 2023)",
                                    self)
        self.actRefRanges.setCheckable(True)
        self.actRefRanges.setToolTip(
            "shade the typical literature δiso range of each species for the "
            "current nucleus (labels carry the P_Q/C_Q ranges) — an "
            "assignment guide, sourced from Edén 2023")
        self.actRefRanges.setChecked(bool(QSettings("LARMOR", "app").value(
            "refRanges", False, type=bool)))
        self.actRefRanges.toggled.connect(self._toggle_ref_ranges)
        m_view.addAction(self.actRefRanges)
        self.actSsbOffer = QAction("&Offer spinning-sideband detection on load",
                                   self)
        self.actSsbOffer.setCheckable(True)
        self.actSsbOffer.setToolTip(
            "when a loaded 1D spectrum repeats at ±νrot, a banner over the "
            "plot offers the linked manifold or a shifted copy in one click; "
            "turn off for spectra whose periodic structure is not sidebands "
            "(Decomposition ▸ Detect spinning sidebands still works on demand)")
        self.actSsbOffer.setChecked(bool(QSettings("LARMOR", "app").value(
            "ssbAutoOffer", True, type=bool)))
        self.actSsbOffer.toggled.connect(self._toggle_ssb_offer)
        m_view.addAction(self.actSsbOffer)
        self._build_czjzek_display_menu(m_view)
        m_view.addSeparator()
        self._build_theme_menu(m_view)
        self._build_textsize_menu(m_view)
        m_view.addSeparator()
        self._build_axis_unit_menu(m_view)
        m_view.addSeparator()
        self._add(m_view, "&Back to 2D map", self.back_to_2d, "Ctrl+2")
        self._add(m_view, "Zoom to sites", self.zoom_sites)
        self._add(m_view, "Full spectrum", self.zoom_full)

        m_tools = mb.addMenu("&Tools")
        m_tools.addSection("Analysis")
        self._add(m_tools, "&Integrals && measurements…  (integral, %, FWHM, CoM)",
                  self.open_integrals)
        self._add(m_tools, "&Batch fit report…  (publication table + plots)",
                  self.run_batch_report)
        self._add(m_tools, "E&xperimental section…  (paragraph + Table S1 from acqus / "
                           "procs / title)", self.open_acquisition_table)
        self._add(m_tools, "&Session inventory…  (a month folder as a sample × nucleus "
                           "grid; production EXPNO picks)", self.open_session_inventory)
        self._add(m_tools, "Batch &fit spectra…  (one shared model, 1D)",
                  lambda: self.explorer._batch_clicked())
        self._add(m_tools, "Se&quential fit…  (forward–backward series sweep, 1D)",
                  self.run_seq_fit)
        self._add(m_tools, "Relaxation / series (T1, T2)…", self.open_satrec)
        self._add(m_tools, "Per-site relaxation…  (uses the current fit)",
                  self.open_per_site_relaxation)
        m_tools.addSection("Advanced experiments")
        self._add(m_tools, "QCPMG (echo train → spectrum)…", self.open_qcpmg)
        self._add(m_tools, "Stitch frequency-stepped (&VOCS) spectra…",
                  self.open_vocs)
        self._add(m_tools, "&Herzfeld–Berger sideband analysis (ζ, η)…  (CSA from "
                           "sideband intensities)", self.open_herzfeld_berger)
        self._add(m_tools, "Referencing a&udit…  (SR of a session against its "
                           "¹H adamantane reference)", self.open_referencing_audit)
        self._add(m_tools, "&Read static pattern (C_Q, η)…  (three markers, "
                           "no fit)", self.open_staticct)
        self._add(m_tools, "QCPMG: infinite-field δiso (2 fields)…",
                  self.open_qcpmg_fields)
        self._add(m_tools, "QCPMG: batch infinite-field δiso…",
                  self.open_qcpmg_batch_fields)
        self._add(m_tools, "Variable temperature (Arrhenius / VFT)…", self.open_vt)
        self._add(m_tools, "REDOR (dipolar coupling)…", self.open_redor)
        m_tools.addSection("Import & 2D")
        self._add(m_tools, "Import DFT tensors (.magres)…", self.open_magres)
        self._add(m_tools, "2D MQMAS viewer/fit…", self.open_twod)
        self._add(m_tools, "Multi-dataset fit (CLI): larmor multifit a.json b.json",
                  lambda: None).setEnabled(False)
        m_tools.addSection("Reference")
        self._add(m_tools, "&NMR table…  (Larmor frequencies)", self.open_nmr_table)
        self._add(m_tools, "&Conversion tools…  (shift / Cq / dipolar)",
                  self.open_convert)

        m_plot = mb.addMenu("&Plotting")
        self._add(m_plot, "Plotting &studio…  (build any figure)",
                  self.open_plotting_studio)
        self._add(m_plot, "Plot &current spectrum…", self.plot_current_spectrum)
        self._add(m_plot, "New &2D contour plot…",
                  lambda: self.open_plotting_studio({"kind": "2d", "path": ""}))

        m_help = mb.addMenu("&Help")
        self.actPalette = self._add(m_help, "&Command palette…  (find any menu entry)",
                                    self.open_command_palette, "Ctrl+Shift+P")
        m_help.addSeparator()
        m_man = m_help.addMenu("User &manuals")
        for name, title in (
                ("getting-started", "Getting started"),
                ("spectra-1d", "1D spectra — processing & fitting"),
                ("lineshapes", "Lineshapes — models & physics"),
                ("glass-fitting", "Fitting glasses for publication (Edén 2023)"),
                ("shift-ranges", "Literature shift ranges — data & sources"),
                ("processing-reference", "Processing reference"),
                ("2d-processing", "2D processing"),
                ("mqmas", "MQMAS (2D)"),
                ("correlation-hmqc", "HMQC & correlation"),
                ("relaxation", "Relaxation (T1/T2)"),
                ("qcpmg", "QCPMG"),
                ("multi-dataset", "Multi-dataset & co-fitting"),
                ("dft-tensors", "DFT tensors — import & shift calibration")):
            self._add(m_man, title,
                      lambda _=False, n=name, t=title: self._open_manual(n, t))
        m_tut = m_help.addMenu("&Tutorials")
        for name, title in TUTORIALS:
            self._add(m_tut, title,
                      lambda _=False, n=name, t=title: self._open_tutorial(n, t))
        self._add(m_help, "About LARMOR", self._about)
        self._add(m_help, "More…", self._show_more)

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

    # -------- Apply recipe: re-use a recent fit's model on the open data ------
    def _add_recent_recipe(self, path: str):
        s = QSettings("LARMOR", "app")
        paths = s.value("recentRecipes", []) or []
        if isinstance(paths, str):
            paths = [paths]
        paths = [p for p in paths if p != path]
        paths.insert(0, path)
        s.setValue("recentRecipes", paths[:10])
        self._rebuild_apply_recipe()

    def _rebuild_apply_recipe(self):
        if not hasattr(self, "m_apply"):
            return
        self.m_apply.clear()
        # always let the user reach a recipe anywhere on disk, not just recents
        browse = self.m_apply.addAction("Browse for recipe…")
        browse.triggered.connect(self.apply_recipe_browse)
        self.m_apply.addSeparator()
        paths = QSettings("LARMOR", "app").value("recentRecipes", []) or []
        if isinstance(paths, str):
            paths = [paths]
        paths = [p for p in paths if Path(p).exists()][:10]
        if not paths:
            a = self.m_apply.addAction("(no recent recipes)")
            a.setEnabled(False)
            return
        for p in paths:
            a = self.m_apply.addAction(Path(p).name)
            a.setToolTip(p)
            a.triggered.connect(lambda _=False, path=p: self.apply_recipe(path))

    def apply_recipe_browse(self):
        """Pick a recipe/dmfit file from anywhere and apply its lines to the data."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Apply recipe (its lines → current data)", self._last_dir(),
            "Fits (*.json *.fxmla *.fxml);;LARMOR recipe (*.json);;All (*)")
        if path:
            self.apply_recipe(path)

    def _recipe_model_from(self, path: str) -> dict | None:
        """The model (sites + nucleus/field) from a LARMOR .json or a dmfit
        .fxml/.fxmla — data not required. Returns None (and warns) on failure."""
        try:
            low = path.lower()
            if low.endswith((".fxml", ".fxmla")):
                from larmor.io import fxmla
                recipe, _warns = fxmla.to_recipe(fxmla.read(path))
                return recipe.to_dict()
            return Recipe.load(path).to_dict()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Apply recipe", f"could not load: {exc}")
            return None

    def _open_manual(self, name: str, title: str):
        from larmor.desktop.help_dialog import show_help

        show_help(self, name, title)

    def open_command_palette(self):
        """Help ▸ Command palette (Ctrl+Shift+P): fuzzy-search every menu and
        toolbar command and run the chosen one exactly as from its menu."""
        from larmor.desktop.palette_dialog import show_palette

        show_palette(self)

    def _open_tutorial(self, name: str, title: str):
        from larmor.desktop.help_dialog import show_help

        show_help(self, name, title, kind="tutorial")

    def _about(self):
        from PySide6.QtWidgets import QDialog, QTextBrowser

        from larmor import __version__

        html = f"""
        <h2 style="margin-bottom:2px;">LARMOR <span style="font-size:14px;
        color:#8a9099; font-weight:normal;">v{__version__}</span></h2>
        <p style="color:#4a4f58; margin-top:0;"><i>An open desktop successor to
        dmfit — solid-state NMR lineshape fitting for disordered solids.</i></p>

        <p>Developed by <b>Sam Soudani</b>, in the <b>McCloy</b> group at
        Washington State University, with support from the
        U.S.&nbsp;Department&nbsp;of&nbsp;Energy&nbsp;(DOE).</p>
        <p><a href="https://github.com/sams808/LARMOR">github.com/sams808/LARMOR</a>
        &nbsp;·&nbsp; MIT&nbsp;License</p>

        <h3 style="margin-bottom:2px;">Built on</h3>
        <ul style="margin-top:2px;">
          <li><b>mrsimulator</b> — solid-state NMR spectrum simulation
              (D.&nbsp;J.&nbsp;Srivastava, P.&nbsp;J.&nbsp;Grandinetti, et&nbsp;al.),
              <a href="https://github.com/deepanshs/mrsimulator">deepanshs/mrsimulator</a></li>
          <li><b>lmfit</b> — non-linear least-squares fitting (M.&nbsp;Newville
              et&nbsp;al., doi:10.5281/zenodo.11813)</li>
          <li><b>NumPy</b> — Harris et&nbsp;al., <i>Nature</i> <b>585</b>, 357 (2020)</li>
          <li><b>SciPy</b> — Virtanen et&nbsp;al., <i>Nature Methods</i> <b>17</b>, 261 (2020)</li>
          <li><b>nmrglue</b> — Helmus &amp; Jaroniec, <i>J.&nbsp;Biomol.&nbsp;NMR</i>
              <b>55</b>, 355 (2013)</li>
          <li><b>csdmpy / CSDM</b> — Srivastava et&nbsp;al., <i>PLOS&nbsp;ONE</i>
              <b>15</b>, e0225953 (2020)</li>
          <li><b>Matplotlib</b> — Hunter, <i>Comput.&nbsp;Sci.&nbsp;Eng.</i>
              <b>9</b>, 90 (2007)</li>
          <li><b>PySide6 / Qt&nbsp;for&nbsp;Python</b> (The Qt Company) and
              <b>pyqtgraph</b> (L.&nbsp;Campagnola) — the interface</li>
        </ul>

        <h3 style="margin-bottom:2px;">Inspired by</h3>
        <ul style="margin-top:2px;">
          <li><b>dmfit</b> — D.&nbsp;Massiot et&nbsp;al., <i>Magn.&nbsp;Reson.&nbsp;Chem.</i>
              <b>40</b>, 70 (2002)</li>
          <li><b>ssNake</b> — S.&nbsp;G.&nbsp;J.&nbsp;van&nbsp;Meerten et&nbsp;al.,
              <i>J.&nbsp;Magn.&nbsp;Reson.</i> <b>301</b>, 56 (2019)</li>
          <li><b>Bruker TopSpin</b> — processing &amp; display conventions</li>
        </ul>
        """
        dlg = QDialog(self)
        dlg.setWindowTitle("About LARMOR")
        dlg.resize(560, 560)
        v = QVBoxLayout(dlg)
        v.setContentsMargins(6, 6, 6, 6)
        tb = QTextBrowser()
        tb.setOpenExternalLinks(True)
        tb.setHtml(html)
        v.addWidget(tb)
        row = QHBoxLayout()
        row.addStretch(1)
        ok = QPushButton("Close")
        ok.clicked.connect(dlg.accept)
        row.addWidget(ok)
        v.addLayout(row)
        dlg.exec()

    def _show_more(self):
        """Open the 'More…' card."""
        try:
            from larmor.xfact import show_fact

            self._more_dlg = show_fact(parent=self)
        except Exception as exc:
            self.statusBar().showMessage(f"…never mind ({exc})")

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
        lab.setStyleSheet("font-weight:600;")     # colour from the theme palette
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
        # a checkable mirror of Process ▸ FID ⇄ spectrum (kept in sync by
        # _mirror_view_actions with setChecked -- no triggered emission)
        self.sbFid = QAction("FID", self)
        self.sbFid.setCheckable(True)
        self.sbFid.setToolTip("show the windowed FID / back to the spectrum "
                              "(Ctrl+T)")
        self.sbFid.triggered.connect(self._toggle_time_domain)
        sb.addAction(self.sbFid)
        sb.addSeparator()
        # a short-labelled sidebar mirror of View ▸ Scroll nudges fit values
        self.sbScroll = QAction("Scroll", self)
        self.sbScroll.setCheckable(True)
        self.sbScroll.setChecked(self.actScrollNudge.isChecked())
        self.sbScroll.setToolTip(self.actScrollNudge.toolTip())
        self.sbScroll.toggled.connect(self.actScrollNudge.setChecked)
        sb.addAction(self.sbScroll)

    def _toggle_scroll_nudge(self, on: bool):
        from larmor.desktop import table as _table

        _table.set_scroll_nudge(on)
        QSettings("LARMOR", "app").setValue("scrollNudge", bool(on))
        if getattr(self, "sbScroll", None) is not None:
            self.sbScroll.setChecked(on)            # keep the sidebar in sync
        if getattr(self, "lines_table", None) and self.recipe:
            self.lines_table.rebuild(self.recipe, self.hidden)   # refresh the hint
        self.statusBar().showMessage(
            "scroll over a fit value to nudge it: "
            + ("ON" if on else "off (default)"))

    def _toggle_ref_ranges(self, on: bool):
        QSettings("LARMOR", "app").setValue("refRanges", bool(on))
        self._update_ref_ranges()
        self.statusBar().showMessage(
            "literature shift ranges: " + ("ON — shaded spans are typical "
            "ranges from Edén 2023, hover a span for the P_Q/C_Q note"
            if on else "off"))

    def _update_ref_ranges(self):
        """(Re)draw the literature-range overlay for the CURRENT nucleus —
        called on toggle and whenever the active 1D document changes."""
        from larmor import refranges

        if not hasattr(self, "view"):
            return
        on = getattr(self, "actRefRanges", None) is not None and \
            self.actRefRanges.isChecked()
        nucleus = (self.recipe or {}).get("nucleus", "") if self.recipe else ""
        ranges = refranges.ranges_for(nucleus) if on else []
        positions = refranges.positions_for(nucleus) if on else []
        for r in ranges + positions:           # per-entry primary citation
            r["_citation"] = refranges.citation_for(r)
        self.view.set_ref_ranges(ranges, refranges.CITATION, positions)
        if on and nucleus and not ranges:
            self.statusBar().showMessage(
                f"no literature ranges compiled for {nucleus} yet — see "
                "Help ▸ Literature shift ranges for what is, and "
                "larmor/refranges.py to extend")

    def _build_czjzek_display_menu(self, parent):
        """View ▸ Czjzek width display — pick which of the four literature
        conventions the fit table shows (and accepts as typed input) for a
        Czjzek site's width. Storage/fitting/exports are always σ; this is
        display-deep only (see larmor.desktop.table.CZJZEK_DISPLAYS)."""
        from PySide6.QtGui import QActionGroup

        from larmor.desktop import table as _table

        m = parent.addMenu("Czjzek &width display")
        m.setToolTip("how the fit table quotes a Czjzek distribution's width "
                     "— saved fits and CSVs always store σ")
        # restore the saved choice into the module BEFORE building the actions
        # (setChecked below fires no signal — same init pattern as scroll-nudge)
        saved = str(QSettings("LARMOR", "app").value("czjzekDisplay", "sigma")
                    or "sigma")
        _table.set_czjzek_display(saved)
        group = QActionGroup(self)
        group.setExclusive(True)
        current = _table.czjzek_display_mode()
        for key, (label, k, desc) in _table.CZJZEK_DISPLAYS.items():
            act = QAction(f"{label}  —  {desc}", self)
            act.setCheckable(True)
            act.setChecked(key == current)
            act.triggered.connect(
                lambda _=False, kk=key: self._set_czjzek_display(kk))
            group.addAction(act)
            m.addAction(act)
        self._czjzek_display_group = group

    def _set_czjzek_display(self, mode: str):
        from larmor.desktop import table as _table

        _table.set_czjzek_display(mode)
        QSettings("LARMOR", "app").setValue("czjzekDisplay", mode)
        if getattr(self, "lines_table", None) and self.recipe:
            self.lines_table.rebuild(self.recipe, self.hidden)
        label = _table.CZJZEK_DISPLAYS[_table.czjzek_display_mode()][0]
        self.statusBar().showMessage(
            f"Czjzek width column now shows {label} — stored fits/CSVs "
            "always keep σ")

    def _build_panels_menu(self):
        """View ▸ Panels: show / hide (and reopen) every dock. Each dock also has
        its own close ✕ now, so a panel can be dismissed and brought back."""
        m = self.m_view.addMenu("&Panels")
        for dock in (self.explorer_dock, self.datasets_dock, self.ws_dock,
                     self.lines_dock, self.results_dock, self.proc_dock):
            m.addAction(dock.toggleViewAction())
        # the strip under the spectrum is not a dock; its switch is remembered
        # (small laptops), the actScrollNudge QSettings pattern
        self.actHealthStrip = QAction("Fit &health strip", self)
        self.actHealthStrip.setCheckable(True)
        self.actHealthStrip.setToolTip("the verdict line under the spectrum: "
                                       "residual, physical values, degenerate "
                                       "pairs, bounds, error bars")
        self.actHealthStrip.setChecked(bool(QSettings("LARMOR", "app").value(
            "fitHealthStrip", True, type=bool)))
        self.actHealthStrip.toggled.connect(self._toggle_health_strip)
        m.addAction(self.actHealthStrip)

    def _toggle_health_strip(self, on: bool):
        QSettings("LARMOR", "app").setValue("fitHealthStrip", bool(on))
        self._health_show()

    def _update_enabled(self):
        loaded = self.recipe is not None
        for a in (self.actSave, self.actFit, self.actQuant):
            a.setEnabled(loaded)
        self.actUndo.setEnabled(bool(self.undo_stack))
        self.actRedo.setEnabled(bool(self.redo_stack))
        if getattr(self, "actHealth", None) is not None:
            self.actHealth.setEnabled(self._health is not None)

    # ------------------------------------------------------------- axis unit
    def _build_axis_unit_menu(self, m_view):
        """View > Axis unit: ppm (default), kHz, MHz. Wideline patterns are
        read and reported in kHz from the reference; display only — every
        internal coordinate stays ppm."""
        from PySide6.QtGui import QActionGroup

        from larmor.desktop.plot import AXIS_UNITS

        self._axis_unit = str(QSettings("LARMOR", "app").value(
            "axisUnit", "ppm") or "ppm")
        if self._axis_unit not in AXIS_UNITS:
            self._axis_unit = "ppm"
        m = m_view.addMenu("Axis &unit")
        group = QActionGroup(self)
        group.setExclusive(True)
        self._axis_unit_actions = {}
        for unit, label in (("ppm", "δ (&ppm)"),
                            ("kHz", "&kHz  (offset from 0 ppm)"),
                            ("MHz", "&MHz  (offset from 0 ppm)")):
            a = QAction(label, self)
            a.setCheckable(True)
            a.setChecked(unit == self._axis_unit)
            a.triggered.connect(lambda _=False, u=unit: self._set_axis_unit(u))
            group.addAction(a)
            m.addAction(a)
            self._axis_unit_actions[unit] = a

    def _set_axis_unit(self, unit: str):
        self._axis_unit = unit
        if not os.environ.get("LARMOR_NO_SESSION"):
            QSettings("LARMOR", "app").setValue("axisUnit", unit)
        self._apply_axis_unit()

    def _apply_axis_unit(self):
        sfo = float((self.recipe or {}).get("larmor_frequency_MHz", 0.0) or 0.0)
        self.view.set_axis_unit(getattr(self, "_axis_unit", "ppm"), sfo)

    def _format_x(self, x_ppm: float) -> str:
        """A cursor/readout position in the display unit (always with ppm)."""
        if self.view.domain == "time":       # the FID curve's data is in ms
            return f"{x_ppm:.3f} ms"
        sfo = float((self.recipe or {}).get("larmor_frequency_MHz", 0.0) or 0.0)
        unit = getattr(self, "_axis_unit", "ppm")
        if unit == "kHz" and sfo > 0:
            return f"{x_ppm * sfo / 1e3:.3f} kHz ({x_ppm:.2f} ppm)"
        if unit == "MHz" and sfo > 0:
            return f"{x_ppm * sfo / 1e6:.5f} MHz ({x_ppm:.2f} ppm)"
        return f"{x_ppm:.2f} ppm"

    # ------------------------------------------------------------- theme
    def _build_theme_menu(self, parent):
        from PySide6.QtGui import QActionGroup

        m = parent.addMenu("&Theme")
        self._theme_group = QActionGroup(self)
        self._theme_group.setExclusive(True)
        current = theme.active().name
        for name in theme.names():
            act = m.addAction(name)
            act.setCheckable(True)
            act.setChecked(name == current)
            act.triggered.connect(lambda _=False, n=name: self._set_theme(n))
            self._theme_group.addAction(act)

        # hidden, just-for-fun styles — tucked into a submenu (not a sibling
        # of the 10 normal presets) so they never clutter the everyday,
        # professional-looking Theme list, but applied live just like any
        # other theme once picked.
        m.addSeparator()
        more = m.addMenu("More styles…")
        more.setToolTip("just-for-fun styles, applied immediately")
        self._aesthetic_group = QActionGroup(self)
        self._aesthetic_group.setExclusive(True)
        current_override = QSettings("LARMOR", "app").value("appearanceOverride", "")
        act_normal = more.addAction("Normal")
        act_normal.setCheckable(True)
        act_normal.setChecked(not current_override)
        act_normal.triggered.connect(lambda _=False: self._set_aesthetic_override(""))
        self._aesthetic_group.addAction(act_normal)
        for name in theme.aesthetic_names():
            act = more.addAction(name)
            act.setCheckable(True)
            act.setChecked(name == current_override)
            act.triggered.connect(
                lambda _=False, n=name: self._set_aesthetic_override(n))
            self._aesthetic_group.addAction(act)

    def _set_aesthetic_override(self, name: str):
        """Switch to a hidden "aesthetic" style (or back to "Normal", which
        restores whichever normal theme was last active) — applied live, the
        same as picking any theme from the main list. Kept as a *separate*
        settings key (``appearanceOverride``) from the main list's ``theme``
        key so "Normal" always knows what to restore, and so picking a
        normal theme later cleanly drops any aesthetic override."""
        settings = QSettings("LARMOR", "app")
        settings.setValue("appearanceOverride", name)
        effective = name or settings.value("theme", theme.DEFAULT)
        self._apply_theme_live(effective)
        # reflect the switch in the main list's checkmarks too: none checked
        # while a genuine aesthetic is active, the matching entry when "Normal"
        for act in self._theme_group.actions():
            act.setChecked((not name) and act.text() == effective)
        self.statusBar().showMessage(f"style: {name or effective}")

    def _build_textsize_menu(self, parent):
        from PySide6.QtGui import QActionGroup

        m = parent.addMenu("Text &size")
        cur = int(QSettings("LARMOR", "app").value("fontPt", 9) or 9)
        self._textsize_group = QActionGroup(self)
        self._textsize_group.setExclusive(True)
        for label, pt in (("Small", 8), ("Normal", 9), ("Large", 11),
                          ("Larger", 13)):
            act = m.addAction(f"{label}  ({pt} pt)")
            act.setCheckable(True)
            act.setChecked(pt == cur)
            act.triggered.connect(lambda _=False, p=pt: self._set_text_size(p))
            self._textsize_group.addAction(act)

    def _set_text_size(self, pt: int):
        """Set the application font size (persists; applies to most widgets at
        once, fully on next launch)."""
        from PySide6.QtWidgets import QApplication

        app = QApplication.instance()
        f = app.font()
        f.setPointSize(int(pt))
        app.setFont(f)
        QSettings("LARMOR", "app").setValue("fontPt", int(pt))
        self.statusBar().showMessage(f"text size: {pt} pt")

    def _apply_theme_live(self, name: str):
        """Apply a theme (normal or hidden aesthetic) to the running app and
        every already-built widget — shared by both the main Theme list and
        the hidden "More styles…" submenu, which differ only in what they
        persist to QSettings."""
        from PySide6.QtWidgets import QApplication

        theme.apply(QApplication.instance(), name)
        # persistent hand-coloured label + the status-bar progress bar
        self.exp_label.setStyleSheet(
            f"color: {theme.active().accent}; font-weight: 600;")
        if hasattr(self, "progress"):
            self._style_progress()
        # the verdict pill and chips take their colours from the theme's
        # signal roles; the request_simulation() below then re-runs
        # _health_live, whose signature check leaves a fit verdict untouched
        self.health_strip.apply_theme()
        # re-theme both plot canvases
        self.view.apply_theme()
        if hasattr(self.view2d, "apply_theme"):
            self.view2d.apply_theme()
        # rebuild the dynamic, series-coloured surfaces
        if self.recipe:
            try:
                self.lines_table.rebuild(self.recipe, self.hidden)
            except Exception:
                pass
            self.request_simulation()
            self._update_paddles()

    def _set_theme(self, name: str):
        """Switch to a normal theme live, remember it, and drop any hidden
        aesthetic override (picking from the visible list is an explicit
        "no, use THIS one" that should always win)."""
        settings = QSettings("LARMOR", "app")
        settings.setValue("theme", name)
        if settings.value("appearanceOverride", ""):
            settings.setValue("appearanceOverride", "")
            if hasattr(self, "_aesthetic_group"):
                for act in self._aesthetic_group.actions():
                    act.setChecked(act.text() == "Normal")
        self._apply_theme_live(name)
        self.statusBar().showMessage(f"theme: {name}")

    def _maybe_show_welcome(self):
        """First-run only: guide a brand-new user on the empty canvas. Anyone who
        has opened data before (a non-empty 'recent' list) never sees it — so it
        does not bother returning users."""
        recent = QSettings("LARMOR", "app").value("recent", []) or []
        if recent or self.exp_ppm is not None:
            return
        self.view.set_placeholder(
            "Open a spectrum to begin\n\n"
            "File ▸ Open  (Ctrl+O)   ·   or drag a Bruker folder / file onto "
            "the plot\n\n"
            "New to LARMOR?   ? ▸ User manuals ▸ Getting started")
