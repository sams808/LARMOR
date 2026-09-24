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
from PySide6.QtWidgets import (QFileDialog, QHBoxLayout, QMessageBox,
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
        """Nine menus a student can scan: File · Edit · Process · Fit · Series
        · Tools · View · Plotting · Help. Labels are short; the explainer of
        every row lives in its tooltip / status tip (menus show tooltips, the
        command palette lists them). Families sit in submenus, groups are
        separated. Every QAction keeps its slot, attribute name and shortcut
        (tests/test_app_split.py pins the tree and the shortcut set)."""
        from PySide6.QtGui import QActionGroup

        mb = self.menuBar()

        # ------------------------------------------------------------ File
        m_file = self._menu(mb, "&File")
        self._add(m_file, "&Open…", self.open_file, "Ctrl+O",
                  tip="open a spectrum, a saved fit, a 1r / 2rr or a project")
        self._add(m_file, "Open &sample…", self.open_sample, "Ctrl+Shift+S",
                  tip="list every spectrum of a sample folder in the Explorer")
        self._add(m_file, "Open &EXPNO / folder…", self.open_expno, "Ctrl+Shift+O",
                  tip="open one Bruker experiment folder")
        self._add(m_file, "Open &FID…", self.open_fid, "Ctrl+F",
                  tip="process a raw fid / ser before the Fourier transform")
        self._add(m_file, "Open &Varian / Agilent…", self.open_varian,
                  tip="a Varian / Agilent .fid folder")
        self.m_recent = self._menu(m_file, "Open &recent")
        self._rebuild_recent()
        m_file.addSeparator()
        self._add(m_file, "O&verlay a spectrum…", self.add_overlay_dialog,
                  "Ctrl+Shift+A",
                  tip="compare a spectrum on top of the active one; the fit is "
                      "kept (Shift + drop on the plot does the same)")
        self.actWatch = self._add(
            m_file, "&Watch the source file", self._toggle_watch, checkable=True,
            tip="re-open the current spectrum whenever its file is rewritten -- "
                "e.g. while acquiring on the spectrometer -- keeping the fit "
                "model so it follows the growing signal")
        m_file.addSeparator()
        self._add(m_file, "Open pro&ject…", self.open_project,
                  tip="a whole session: spectra, 2D maps, figures, batch fits")
        self._add(m_file, "Save projec&t…", self.save_project, "Ctrl+Alt+S",
                  tip="save the whole session: spectra, 2D maps, figures, batch fits")
        m_file.addSeparator()
        self.actSave = self._add(m_file, "&Save fit", self.save_recipe, "Ctrl+S",
                                 tip="write the fit recipe (.recipe.json) next to "
                                     "the data -- the reproducible analysis")
        self._add(m_file, "Save fit &as…", self.save_fit_as, "Ctrl+Shift+E",
                  tip="txt / csv / json / dmfit .fxmla")
        m_export = self._menu(m_file, "E&xport",
                              tip="figures, tables and spectra out of LARMOR")
        self._add(m_export, "Save s&pectrum as CSV…", self.save_spectrum,
                  tip="the processed spectrum as a CSV with a metadata header, "
                      "reopenable in LARMOR")
        self._add(m_export, "&Figure…", self.open_figure_dialog,
                  tip="a publication figure of the current fit")
        self._add(m_export, "Save plot &image…", self.save_plot_image,
                  tip="png / svg of the plot as shown")
        self._add(m_export, "&Copy plot to clipboard", self.copy_plot, "Ctrl+Shift+C",
                  tip="the plot with all its lines, as an image")
        m_export.addSeparator()
        self._add(m_export, "Copy report &table as CSV", self.copy_csv,
                  tip="the last Report's table (positions, integrals, fractions)")
        self._add(m_export, "Copy &LaTeX table", self.copy_latex,
                  tip="the fit as a LaTeX table with uncertainties")
        self._add(m_export, "&Publication bundle…", self.export_publication_bundle,
                  tip="figure + LaTeX table + CSV + methods sentence + report.md "
                      "into one folder")
        m_file.addSeparator()
        self._add(m_file, "&Quit", self.close)

        # ------------------------------------------------------------ Edit
        m_edit = self._menu(mb, "&Edit")
        self.actUndo = self._add(m_edit, "↩  Undo", self.undo, "Ctrl+Z",
                                 tip="undo the last model or processing change")
        self.actRedo = self._add(m_edit, "↪  Redo", self.redo, "Ctrl+Y",
                                 tip="redo the change just undone")
        m_edit.addSeparator()
        self._add(m_edit, "&New fit", self.new_fit,
                  tip="clear every line and start the model again")
        m_edit.addSeparator()
        m_models = self._menu(m_edit, "Add &line",
                              tip="pick a lineshape model, then click on the "
                                  "spectrum to place the line")
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
            self._model_actions[m["name"]] = a
        for section, names in self._MODEL_GROUPS:
            m_models.addSection(section)
            for name in names:
                if name in self._model_actions:
                    m_models.addAction(self._model_actions[name])
        grouped = {n for _, names in self._MODEL_GROUPS for n in names}
        rest = [n for n in self._model_actions if n not in grouped]
        if rest:                  # a model the groups do not know yet
            m_models.addSection("More")
            for name in rest:
                m_models.addAction(self._model_actions[name])
        self._add(m_edit, "Add a line at every &peak…", self.autopick_lines,
                  tip="auto peak-pick, one line per maximum")
        self._add(m_edit, "Add f&unction line…", self.add_function_line,
                  tip="y = f(x; a, b, c, d) -- a free function as a line")
        self._add(m_edit, "Add background &spectrum…", self.add_background_spectrum,
                  tip="fit another spectrum as a scaled component")
        self.m_apply = self._menu(m_edit, "&Apply recipe",
                                  tip="re-use a recent fit's model on this data")
        self._rebuild_apply_recipe()
        m_edit.addSeparator()
        m_ssb = self._menu(m_edit, "Spinning &sidebands",
                           tip="find and model the ±νrot repeats of a spectrum")
        self.actDetectSsb = self._add(
            m_ssb, "&Detect spinning sidebands", self.detect_sidebands, "Ctrl+Shift+D",
            tip="look for the repeat of the spectrum at ±νrot (around the "
                "recorded rate, then 1–80 kHz) and offer a linked sideband "
                "manifold or a shifted copy of the spectrum in one click")
        self._add(m_ssb, "Add spinning &sidebands…", self.add_sidebands,
                  tip="a linked ±νrot manifold of a fitted line")
        self._add(m_ssb, "Add a &copy of this spectrum…",
                  self.add_current_spectrum_line,
                  tip="a shifted copy of the spectrum itself -- a satellite or "
                      "sideband manifold")
        m_ssb.addSeparator()
        self.actSsbOffer = QAction("&Offer spinning-sideband detection on load",
                                   self)
        self.actSsbOffer.setCheckable(True)
        self.actSsbOffer.setToolTip(
            "when a loaded 1D spectrum repeats at ±νrot, a banner over the "
            "plot offers the linked manifold or a shifted copy in one click; "
            "turn off for spectra whose periodic structure is not sidebands "
            "(Detect spinning sidebands still works on demand)")
        self.actSsbOffer.setChecked(bool(QSettings("LARMOR", "app").value(
            "ssbAutoOffer", True, type=bool)))
        self.actSsbOffer.toggled.connect(self._toggle_ssb_offer)
        m_ssb.addAction(self.actSsbOffer)
        m_con = self._menu(m_edit, "&Constraints",
                           tip="what the fit may vary: literature ranges, "
                               "glass protocol, saved link / bound sets")
        self._add(m_con, "&Label lines from literature ranges",
                  self.label_from_literature,
                  tip="name each line by the species whose δiso range it falls "
                      "in (Al[4], BO3, …)")
        self._add(m_con, "&Restrict around current values…",
                  self.restrict_glass_protocol,
                  tip="bound every parameter around its current value -- the "
                      "glass protocol (Edén 2023)")
        self._add(m_con, "Save &constraints as…", self.save_constraint_set,
                  tip="a reusable link / bound set")
        self._add(m_con, "Apply saved co&nstraints…", self.apply_constraint_set,
                  tip="a link / bound set saved earlier")
        m_edit.addSeparator()
        self._add(m_edit, "Add fit &zone", self.add_zone,
                  tip="a ppm window the fit is restricted to")
        self._add(m_edit, "Clear zones", self.clear_zones)

        # --------------------------------------------------------- Process
        m_proc = self._menu(mb, "&Process")
        self.actExp = self._add(m_proc, "&Experiment parameters…", self.edit_experiment,
                                tip="νrot, B0, nucleus, 90° pulse -- and where "
                                    "each value came from")
        self._add(m_proc, "Processing s&teps…", self.edit_processing_steps,
                  tip="the recorded processing chain; remove a step")
        self._add(m_proc, "Show processing &panel",
                  lambda: self.proc_dock.show(),
                  tip="source, display, phase, baseline and reference controls")
        m_proc.addSeparator()
        m_phase = self._menu(m_proc, "&Phase")
        self._add(m_phase, "&Autophase (ACME)",
                  lambda: self.apply_processing([{"op": "autophase"}], False),
                  tip="automatic p0 / p1 by entropy minimisation")
        # non-checkable: the panel's Drag-to-phase button is the single
        # source of truth for the mode, this entry toggles it
        self._add(m_phase, "&Drag to phase", self.start_phase_drag, "Ctrl+P",
                  tip="TopSpin-style: drag on the spectrum -- ← → p0, ↑ ↓ p1 "
                      "about the pivot line")
        m_base = self._menu(m_proc, "&Baseline")
        self._add(m_base, "&Polynomial (order 3)",
                  lambda: self.apply_processing([{"op": "baseline", "order": 3}], False),
                  tip="automatic polynomial baseline; the panel sets the order")
        self._add(m_base, "&Iterative (dead-time; Yon 2020)…", self.apply_iterbaseline,
                  tip="the dead-time baseline roll removed iteratively")
        self._add(m_base, "&2-point background…", self.start_twopoint_bg,
                  tip="pick two flat points; the straight line through them "
                      "is subtracted")
        self._add(m_base, "Subtract &averages",
                  lambda: self.apply_processing([{"op": "subtract_avg"}], False),
                  tip="remove the offset read from the spectrum edges")
        m_ref = self._menu(m_proc, "&Reference")
        self._add(m_ref, "&Calibrate axis…", self.start_calibrate,
                  tip="click a peak and type its ppm; the fit follows")
        self.actMeasure = self._add(m_ref, "&Measure Δ (ppm / Hz)",
                                    self.toggle_measure, checkable=True,
                                    tip="drag between two points to read their "
                                        "distance")
        self._add(m_ref, "Referencing a&udit…", self.open_referencing_audit,
                  tip="the SR of a whole session against its ¹H adamantane "
                      "reference")
        m_alg = self._menu(m_proc, "Re&gion / algebra",
                           tip="measure regions, combine spectra")
        self._add(m_alg, "&Integrals && measurements…", self.open_integrals,
                  tip="integral, %, FWHM and centre of mass of dragged regions")
        self._add(m_alg, "Subtract a spectrum (&background)…", self.open_subtract,
                  tip="remove a measured background spectrum")
        self._add(m_alg, "&WURST excitation profile…", self.open_wurst_correct,
                  tip="divide out the sweep's amplitude envelope")
        self._add(m_alg, "Stitch frequency-stepped (&VOCS) spectra…", self.open_vocs,
                  tip="one pattern from several offset acquisitions")
        m_proc.addSeparator()
        # display projections of the pipeline result: the FID the transform
        # sees, and the real / imaginary / |S| channel
        self.actTimeDomain = self._add(
            m_proc, "&FID ⇄ spectrum", self._toggle_time_domain, "Ctrl+T",
            checkable=True,
            tip="show the windowed FID the transform sees and re-apply the window "
                "functions live; press again to return to the spectrum")
        m_chan = self._menu(m_proc, "Display &channel",
                            tip="which channel of the complex spectrum the plot "
                                "shows -- display only, the fit uses the real part")
        chan_group = QActionGroup(self)
        chan_group.setExclusive(True)
        self.actChannel = {}
        for name, label, tip in (
                ("real", "&Real", "the real channel -- what the fit uses"),
                ("imag", "&Imaginary", "inspect the dispersion while phasing"),
                ("magnitude", "&Magnitude |S|",
                 "display only -- the destructive pipeline op is the panel's "
                 "checkbox")):
            a = self._add(m_chan, label,
                          lambda _=False, n=name: self._set_channel(n),
                          checkable=True, checked=(name == "real"), tip=tip)
            chan_group.addAction(a)
            self.actChannel[name] = a
        m_chan.addSeparator()
        self._add(m_chan, "C&ycle channel", self._cycle_channel, "Ctrl+I",
                  tip="real → imaginary → |S| → real")
        m_proc.addSeparator()
        self._add(m_proc, "Reset to &original", self.reset_processing,
                  tip="drop every processing step and reload the source")

        # ------------------------------------------------------------- Fit
        m_fit = self._menu(mb, "F&it")
        self._add(m_fit, "&Simulate", self.request_simulation, "F9",
                  tip="recompute the model without fitting")
        self.actFit = self._add(m_fit, "&Fit", self.run_fit, "F5",
                                tip="least-squares fit of the free parameters")
        self.actAuto = self._add(m_fit, "&Auto fit…", self.run_auto_fit,
                                 tip="multi-start refits from randomised starting "
                                     "values; keeps the best")
        m_fit.addSeparator()
        self.actQuant = self._add(m_fit, "&Report", self.run_quantify, "F6",
                                  tip="quantify: integrals, fractions and "
                                      "uncertainties of every line")
        # F5 Fit / F6 Report / F7 Health form one cluster
        self.actHealth = self._add(m_fit, "Fit &health details…", self.show_fit_health,
                                   "F7", tip="why the fit-health strip says what "
                                             "it says")
        m_err = self._menu(m_fit, "&Errors",
                           tip="uncertainties beyond the covariance matrix")
        self.actErrors = self._add(m_err, "χ² &profile…", self.run_errors_analysis,
                                   tip="Errors Analysis: refit while one parameter "
                                       "is stepped -- its confidence interval")
        # kept as attributes: the fit-health strip's menu lists these SAME
        # QAction objects (no duplicated labels or slots)
        self.actMC = self._add(m_err, "Monte-&Carlo errors…", self.run_monte_carlo,
                               tip="synthetic-noise refits -- a parametric bootstrap")
        self.actCorr = self._add(m_err, "Parameter correlations…",
                                 self.show_correlations,
                                 tip="the correlation matrix of the last fit")
        self.actChi2 = self._add(m_err, "χ² map (parameter pair)…", self.show_chi2_map,
                                 tip="χ² over a plane of two parameters -- is the "
                                     "pair determined?")
        self._add(m_fit, "Compare with a saved fit…", self.compare_with_saved_fit,
                  tip="parameter-by-parameter diff against a saved recipe")
        m_fit.addSeparator()
        self._add(m_fit, "Co-&fit datasets…", self.open_cofit,
                  tip="one shared model over several datasets, 1D + MQMAS")
        self._add(m_fit, "&Predict at another field…", self.predict_at_field,
                  tip="re-simulate the current model at another B0")
        m_mq = self._menu(m_fit, "&MQMAS")
        self._add(m_mq, "2D MQMAS &viewer / fit…", self.open_twod,
                  tip="the 2D contour viewer and its fit")
        self._add(m_mq, "MQMAS F1 &reference…", self.edit_mqmas_f1_ref,
                  tip="pin the isotropic-axis alignment (dmfit style)")
        m_set = self._menu(m_fit, "Fit se&ttings")
        self._add(m_set, "Computing &parameters…", self.edit_computing_params,
                  tip="kernel resolution and simulation grid")
        self._add(m_set, "Fit completion &threshold…", self.edit_fit_tol,
                  tip="the Δσ % below which a fit stops")
        self.actAnimateFit = QAction("&Animate fits", self)
        self.actAnimateFit.setCheckable(True)
        self.actAnimateFit.setToolTip("draw the model curve as it converges during "
                                      "a 1D fit (a fading trail shows the last few "
                                      "iterations) — watch convergence or divergence")
        self.actAnimateFit.setChecked(bool(QSettings("LARMOR", "app").value(
            "animateFit", True, type=bool)))
        self.actAnimateFit.toggled.connect(
            lambda on: QSettings("LARMOR", "app").setValue("animateFit", bool(on)))
        m_set.addAction(self.actAnimateFit)

        # ---------------------------------------------------------- Series
        m_ser = self._menu(mb, "&Series")
        self._add(m_ser, "Batch &fit spectra…",
                  lambda: self.explorer._batch_clicked(),
                  tip="one shared model over the spectra selected in the "
                      "Explorer (1D)")
        self._add(m_ser, "Se&quential fit…", self.run_seq_fit,
                  tip="a forward–backward sweep along a series, each fit "
                      "starting from its neighbour (1D)")
        self._add(m_ser, "&Batch fit report…", self.run_batch_report,
                  tip="publication table and plots from a batch fit")
        m_ser.addSeparator()
        self._add(m_ser, "&Session inventory…", self.open_session_inventory,
                  tip="a month folder as a sample × nucleus grid, with the "
                      "production EXPNO picked per block")
        self._add(m_ser, "E&xperimental section…", self.open_acquisition_table,
                  tip="a methods paragraph and Table S1 from acqus / procs / title")
        self._add(m_ser, "&Compare acquisition parameters…", self.compare_overlays,
                  tip="acqus / procs of the active spectrum and every overlay "
                      "side by side")

        # ----------------------------------------------------------- Tools
        m_tools = self._menu(mb, "&Tools")
        self._add(m_tools, "&NMR table…", self.open_nmr_table,
                  tip="Larmor frequencies, spins and abundances")
        self._add(m_tools, "&Conversion tools…", self.open_convert,
                  tip="shift, Cq and dipolar conversions")
        m_tools.addSeparator()
        self._add(m_tools, "&Herzfeld–Berger sideband analysis…",
                  self.open_herzfeld_berger,
                  tip="ζ and η of the CSA from spinning-sideband intensities")
        self._add(m_tools, "&Read static pattern (C_Q, η)…", self.open_staticct,
                  tip="three markers on a static lineshape, no fit")
        self._add(m_tools, "Czjzek &distribution P(C_Q)…", self.show_czjzek_dist,
                  tip="what a Czjzek σ stands for, drawn")
        m_tools.addSeparator()
        self._add(m_tools, "Import &DFT tensors (.magres)…", self.open_magres,
                  tip="shielding and EFG tensors from a DFT run, calibrated to "
                      "shifts")
        m_rel = self._menu(m_tools, "Re&laxation")
        self._add(m_rel, "Relaxation / series (T1, T2)…", self.open_satrec,
                  tip="a guided saturation-recovery / echo-train workflow")
        self._add(m_rel, "Per-site relaxation…", self.open_per_site_relaxation,
                  tip="T1 of each fitted line, using the current model")
        self._add(m_rel, "Variable temperature (Arrhenius / VFT)…", self.open_vt,
                  tip="activation energy from a τ(T) series")
        m_qc = self._menu(m_tools, "&QCPMG")
        self._add(m_qc, "QCPMG (echo train → spectrum)…", self.open_qcpmg,
                  tip="process an echo train into a spikelet or envelope spectrum")
        self._add(m_qc, "Infinite-field δiso (2 fields)…", self.open_qcpmg_fields,
                  tip="extrapolate δiso from two fields")
        self._add(m_qc, "Batch infinite-field δiso…", self.open_qcpmg_batch_fields,
                  tip="the two-field extrapolation for a whole series")
        self._add(m_tools, "R&EDOR (dipolar coupling)…", self.open_redor,
                  tip="dipolar coupling from a REDOR dephasing curve")

        # ------------------------------------------------------------ View
        m_view = self._menu(mb, "&View")
        self.m_view = m_view
        self.actResid = self._add(m_view, "Residual", self._toggle_resid,
                                  checkable=True, checked=True,
                                  tip="the data − model trace under the spectrum")
        self.actComp = self._add(m_view, "Components", self._toggle_comp,
                                 checkable=True, checked=True,
                                 tip="each line's own curve")
        self.actLabels = self._add(
            m_view, "Component &labels", self._toggle_labels, checkable=True,
            checked=bool(QSettings("LARMOR", "app").value(
                "compLabels", False, type=bool)),
            tip="write each component's letter and name at its maximum; when "
                "off, hovering a component still shows its name")
        self.view.set_show_labels(self.actLabels.isChecked())
        self.actPaddles = self._add(m_view, "&Paddles", self._toggle_paddles,
                                    checkable=True, checked=True,
                                    tip="the on-spectrum handles that move a line")
        self.actOverlaysVisible = self._add(
            m_view, "O&verlays", self._toggle_overlays_visible,
            "Ctrl+Shift+V", checkable=True, checked=True,
            tip="show or hide every compared spectrum at once; add one with "
                "File ▸ Overlay a spectrum, Shift + drop on the plot, or the "
                "Explorer's right-click")
        self._add(m_view, "Clear overlays", self.clear_overlays,
                  tip="drop every compared spectrum")
        self.actRefRanges = QAction("&Literature shift ranges", self)
        self.actRefRanges.setCheckable(True)
        self.actRefRanges.setToolTip(
            "shade the typical literature δiso range of each species for the "
            "current nucleus (labels carry the P_Q/C_Q ranges) — an "
            "assignment guide, sourced from Edén 2023")
        self.actRefRanges.setChecked(bool(QSettings("LARMOR", "app").value(
            "refRanges", False, type=bool)))
        self.actRefRanges.toggled.connect(self._toggle_ref_ranges)
        m_view.addAction(self.actRefRanges)
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
        m_view.addSeparator()
        m_zoom = self._menu(m_view, "&Zoom")
        self._add(m_zoom, "&Full spectrum", self.zoom_full)
        self._add(m_zoom, "Zoom to &sites", self.zoom_sites,
                  tip="the fitted region")
        self._add(m_zoom, "&Back to 2D map", self.back_to_2d, "Ctrl+2",
                  tip="return to the parent 2D map of an extracted trace")
        m_view.addSeparator()
        self._build_axis_unit_menu(m_view)
        self._build_czjzek_display_menu(m_view)
        m_view.addSeparator()
        # filled by _build_panels_menu once the docks exist
        self.m_panels = self._menu(m_view, "&Panels",
                                   tip="show or hide each dock")
        self._build_theme_menu(m_view)
        self._build_textsize_menu(m_view)

        # -------------------------------------------------------- Plotting
        m_plot = self._menu(mb, "&Plotting")
        self._add(m_plot, "Plotting &studio…", self.open_plotting_studio,
                  tip="build any figure from open spectra, fits and series")
        self._add(m_plot, "Plot &current spectrum…", self.plot_current_spectrum,
                  tip="the spectrum on screen, in the studio")
        self._add(m_plot, "New &2D contour plot…",
                  lambda: self.open_plotting_studio({"kind": "2d", "path": ""}))

        # ------------------------------------------------------------ Help
        m_help = self._menu(mb, "&Help")
        self.actPalette = self._add(m_help, "&Command palette…",
                                    self.open_command_palette, "Ctrl+Shift+P",
                                    tip="find and run any menu entry by typing "
                                        "part of its name")
        m_help.addSeparator()
        m_man = self._menu(m_help, "User &manuals")
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
        m_tut = self._menu(m_help, "&Tutorials")
        for name, title in TUTORIALS:
            self._add(m_tut, title,
                      lambda _=False, n=name, t=title: self._open_tutorial(n, t))
        m_help.addSeparator()
        self._add(m_help, "About LARMOR", self._about)
        self._add(m_help, "More…", self._show_more)

    #: Edit ▸ Add line sections, by model name (registry order inside each)
    _MODEL_GROUPS = (
        ("Simple", ("gauss_lor", "gl_norm", "voigt", "jmultiplet")),
        ("Quadrupolar", ("quad_ct", "quad_first", "quad_csa")),
        ("Disordered", ("czjzek", "ext_czjzek", "czjzek_d", "czjzek_corr",
                        "amorphous")),
        ("CSA", ("csa_mas", "csa_czjzek")),
        ("Other", ("sidebands", "exchange2", "function")),
    )

    def _menu(self, parent, title, tip=None):
        """A (sub)menu that shows its rows' tooltips on hover."""
        m = parent.addMenu(title)
        m.setToolTipsVisible(True)
        if tip:
            m.menuAction().setToolTip(tip)
            m.menuAction().setStatusTip(tip)
        return m

    def _add(self, menu, text, slot, shortcut=None, checkable=False, checked=False,
             tip=None):
        a = QAction(text, self)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        a.setCheckable(checkable)
        a.setChecked(checked)
        if tip:
            a.setToolTip(tip)
            a.setStatusTip(tip)
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
        if not hasattr(self, "actUndo"):        # Edit ▸ Undo / Redo own them
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
        # the "＋ Add line" split button, quick buttons and placing label
        # (larmor/desktop/addline_toolbar.py) -- the same model QActions
        from larmor.desktop.addline_toolbar import build_addline_toolbar
        self._addline_toolbar = build_addline_toolbar(self, tb)
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
        its own close ✕ now, so a panel can be dismissed and brought back. The
        submenu itself is created by _build_menus (so it sits with Theme and
        Text size); this fills it once the docks exist."""
        m = self.m_panels
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
        if getattr(self, "_addline_toolbar", None) is not None:
            self._addline_toolbar.apply_theme()
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
