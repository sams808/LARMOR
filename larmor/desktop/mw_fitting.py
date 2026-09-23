"""Main-window mixin: simulate, fit, judge, report.

The simulation debounce and ``SimWorker``, the 1D and 2D fit runs and their
completions, the fit-health verdict (``_health_*``, rendered by the strip
under the plot), quantification, the CSV / LaTeX / Methods clipboard actions,
the publication bundle and auto fit.

Owned state: ``_sim_worker`` / ``_sim_pending`` / ``_busy`` (the timers are
created in ``__init__``), ``_fit_worker`` / ``_fit2d_worker`` /
``_active_fit_worker``, ``_anim_last_ms`` / ``_anim_last_rms``,
``_last_quant``, ``_last_model``, ``_first_sim``, ``_last_lmfit``,
``_health`` / ``_health_fit``.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QFileDialog, QMessageBox,
                               QTableWidgetItem)

from larmor import fithealth
from larmor.desktop.workers import (Fit2DWorker, FitWorker, SimWorker,
                                    humanize_error)
from larmor.recipe import Recipe


class _FittingMixin:
    """Simulation, fitting, fit health, quantification, exports."""

    # ------------------------------------------------------------- simulate
    def _kernel_progress_tick(self, done: int, total: int):
        """A cold Czjzek kernel build in a worker (~23 s wideline) reports
        its chunks here -- previously indistinguishable from a hang."""
        self.statusBar().showMessage(
            f"building lineshape kernel — {done}/{total}"
            + ("  (done)" if done >= total else " …"), 4000)

    def request_simulation(self):
        if getattr(self, "_cofit", None) is not None and self._cofit_active():
            self._cofit_timer.start()             # live preview on both panels
            return
        self._sim_timer.start()

    def _simulate_now(self):
        # the 2D map has its own (explicit) fit path; skip the live 1D sim there
        if self.central_stack.currentWidget() in (self.view2d, self.cofit_page):
            return
        if not self.recipe or not self.recipe["sites"]:
            self.view.set_model(None, None, None, None, self.hidden)
            self._health_reset()          # New fit / last line deleted: no verdict
            return
        if self._sim_worker and self._sim_worker.isRunning():
            self._sim_pending = True
            return
        self._sim_worker = SimWorker(json.loads(json.dumps(self.recipe)),
                                     self.exp_ppm)
        self._sim_worker.done.connect(self._sim_done)
        self._sim_worker.failed.connect(self._on_sim_failed)
        self._sim_worker.kernel_progress.connect(self._kernel_progress_tick)
        self._busy_timer.start()          # arm the busy cue for a slow sim
        self._sim_worker.start()

    def _sim_busy_on(self):
        """Fired only if a sim outlives the busy-timer: show a wait cursor and
        say why (usually a first-time kernel build)."""
        from PySide6.QtGui import QCursor
        from PySide6.QtWidgets import QApplication

        self._busy = True
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        sites = (self.recipe.get("sites") if self.recipe else None) or []
        needs_kernel = any(s.get("model") in ("czjzek", "czjzek_d", "czjzek_corr",
                                              "ext_czjzek", "amorphous")
                           for s in sites)
        self.statusBar().showMessage(
            "building the lineshape kernel (first Czjzek/Amorphous fit is slow, "
            "then instant)…" if needs_kernel else "computing…")

    def _sim_busy_off(self):
        self._busy_timer.stop()
        if self._busy:
            from PySide6.QtWidgets import QApplication
            QApplication.restoreOverrideCursor()
            self._busy = False

    def _on_sim_failed(self, msg):
        self._sim_busy_off()
        self.statusBar().showMessage("simulate: " + msg)

    def _sim_done(self, x, total, per_site):
        self._sim_busy_off()
        labels = [s.get("label") or s["model"] for s in self.recipe["sites"]]
        self._last_model = (np.asarray(x), np.asarray(total))
        self.view.set_model(x, total, per_site, labels, self.hidden,
                            self.exp_ppm, self.exp_amp)
        # the single debounced hook for table edits, paddles, undo/redo,
        # processing, workspace switches and theme re-sims; during a paddle
        # drag only the signature compare runs (the release does a full pass)
        self._health_live(full=not self._paddle_live)
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
        if self.view.domain == "time":
            # the fit window is the viewbox X range (below): in FID view that
            # would be a millisecond window. Covers the menu, F5 and the
            # table's Fit button in one place.
            self.statusBar().showMessage(
                "return to the spectrum first (FID ⇄ spectrum, Ctrl+T) — the "
                "fit window is read from the frequency axis")
            return
        if self._fit_worker and self._fit_worker.isRunning():
            return
        self._sanitize_constraints_before_fit()
        self.snapshot()
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        hi, lo = max(x0, x1), min(x0, x1)
        zones = self.recipe.get("fit_zones") or []
        self.statusBar().showMessage(
            f"fitting in {len(zones)} zone(s) …" if zones
            else f"fitting in {hi:.1f} … {lo:.1f} ppm …")
        self.lines_table.btnFit.setEnabled(False)
        self._progress_start("1D fit")
        animate = bool(getattr(self, "actAnimateFit", None)
                       and self.actAnimateFit.isChecked())
        self._fit_worker = FitWorker(json.loads(json.dumps(self.recipe)),
                                     self.exp_ppm, self.exp_amp, (hi, lo),
                                     animate=animate)
        self._fit_worker.progress.connect(self._progress_tick)
        self._fit_worker.kernel_progress.connect(self._kernel_progress_tick)
        self._fit_worker.done.connect(self._fit_done)
        self._fit_worker.failed.connect(self._fit_failed)
        if animate:
            self._anim_last_ms = 0.0
            self._anim_last_rms = float("nan")
            self._fit_worker.progress.connect(
                lambda it, rms: setattr(self, "_anim_last_rms", rms))
            self._fit_worker.frame.connect(self._fit_frame)
            self.view.start_fit_animation()
        self._active_fit_worker = self._fit_worker
        self._show_fit_buttons(True)
        self._fit_worker.start()

    def _fit_frame(self, x, y, iteration: int):
        """Draw one animation frame. Frames are already emitted only every ~10
        iterations (larmor.fit frame_every), so no per-iteration redraw cost."""
        self.view.set_fit_frame(x, y, iteration,
                                getattr(self, "_anim_last_rms", float("nan")))

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
        self._sanitize_constraints_before_fit()
        self.snapshot()
        self.statusBar().showMessage(
            "fitting the 2D (building the MQMAS kernel, first time is slow)…")
        self.lines_table.btnFit.setEnabled(False)
        self._progress_start("2D MQMAS fit")
        self._fit2d_worker = Fit2DWorker(
            json.loads(json.dumps(self.recipe)), self._data2d, method)
        self._fit2d_worker.progress.connect(self._progress_tick)
        self._fit2d_worker.done.connect(self._fit2d_done)
        self._fit2d_worker.failed.connect(self._fit_failed)
        self._active_fit_worker = self._fit2d_worker
        self._show_fit_buttons(True)
        self._fit2d_worker.start()

    def _fit2d_done(self, result, stop_mode: str = ""):
        self._active_fit_worker = None
        self.lines_table.btnFit.setEnabled(True)
        if stop_mode == "cancel":
            self._progress_end(False)
            self.statusBar().showMessage("2D fit cancelled — parameters unchanged")
            return
        self._progress_end(True)
        self.recipe = result.recipe.to_dict()
        self.lines_table.rebuild(self.recipe, self.hidden)
        from larmor import twod

        self.view2d.set_model(result.z_fit, result.kernel.f2_ppm,
                              twod.mqmas_f1_axis(result.kernel, result.recipe),
                              per_site=result.per_site)
        self.lines_table.set_chi2(f"RMSD {result.rmsd:.4f}")
        f1ref = getattr(result.recipe, "mqmas_f1_ref_ppm", 0.0)
        held = not getattr(result.recipe, "mqmas_f1_ref_vary", True)
        refmsg = f" · F1 ref {f1ref:+.1f} ppm{' (fixed)' if held else ''}"
        self.results_summary.setText(f"2D MQMAS fit · RMSD {result.rmsd:.4f}{refmsg}")
        self.report.setPlainText(
            f"MQMAS F1 isotropic-axis reference offset: {f1ref:+.2f} ppm"
            f"{' (held fixed)' if held else ' (auto-fitted)'}\n\n" + result.report)
        self.statusBar().showMessage(f"2D fit done · RMSD {result.rmsd:.4f}{refmsg}")
        self._persist_session()

    def _sanitize_constraints_before_fit(self):
        """Repair invalid links (self-references / cycles) so the fit can't loop
        forever — the common cause of a recipe that errors on load-then-fit."""
        if not self.recipe or not self.recipe.get("sites"):
            return
        from larmor.constraints_util import sanitize_constraints
        dropped = sanitize_constraints(self.recipe["sites"])
        if dropped:
            self.lines_table.rebuild(self.recipe, self.hidden)
            QMessageBox.information(
                self, "Fixed invalid line links",
                "This model had constraint(s) that reference a line itself or "
                "form a loop, which would make the fit run forever. LARMOR "
                "removed them so you can fit:\n\n  " + ", ".join(dropped)
                + "\n\nRe-add the links you meant (each must reference ANOTHER "
                "line).")

    def _fit_failed(self, msg: str):
        self._active_fit_worker = None
        self.lines_table.btnFit.setEnabled(True)
        self.view.stop_fit_animation()
        self._progress_end(False)
        QMessageBox.warning(self, "Fit failed", humanize_error(msg))
        self.statusBar().showMessage("fit failed")

    @staticmethod
    def _residual_noise_ratio(result):
        """Delegate to fithealth.residual_noise_ratio (the body moved to the
        Qt-free core with the fit-health strip); kept so callers holding a
        duck-typed result still work. To remove when app.py is split."""
        return fithealth.residual_noise_ratio(getattr(result, "y_exp", None),
                                              getattr(result, "y_fit", None))

    # ------------------------------------------------------------- fit health
    def _health_from_result(self, result):
        """Judge a finished fit (plain or Auto Fit): one fithealth.Health from
        the FitResult, the quantify rows of THIS fit and the window's data;
        remembered as the reference the live passes compare against."""
        rl = getattr(result, "lmfit_result", None)
        self._last_lmfit = rl                              # for correlations
        rows = (self._last_quant or {}).get("rows")
        h = fithealth.assess(
            result.recipe, result.y_exp, result.y_fit, lmfit_result=rl,
            at_bounds=result.at_bounds or [], frozen=result.frozen_sites or [],
            window=getattr(result.recipe, "fit_window_ppm", None),
            rmsd=result.rmsd, quant_rows=rows, fitted=True, ppm=self.exp_ppm,
            x_fit=getattr(result, "x_ppm", None))
        h.recipe_sig = fithealth.recipe_signature(self.recipe)
        h.data_sig = fithealth.data_signature(self.exp_ppm, self.exp_amp)
        self._health_fit = h
        self._health_apply(h)
        return h

    def _health_live(self, full: bool = True):
        """Re-judge the live model after the debounced simulation landed.
        ``full=False`` (a paddle drag) only compares the cheap signatures and
        marks the pill; the full pass interpolates the model onto the
        experimental axis and reruns the residual + physical checks, carrying
        the fit's covariance flags forward as stale."""
        if (self._last_model is None or not len(self.exp_ppm)
                or not (self.recipe and self.recipe.get("sites"))):
            return
        base = self._health_fit
        if not full:
            if base is None or (self._health is not None and self._health.stale):
                return
            if (fithealth.recipe_signature(self.recipe) != base.recipe_sig
                    or fithealth.data_signature(self.exp_ppm, self.exp_amp)
                    != base.data_sig):
                self.health_strip.set_stale_hint(True)
            return
        x, tot = self._last_model
        xp = np.asarray(self.exp_ppm)
        if x.shape == xp.shape and np.array_equal(x, xp):
            y = tot
        else:
            if x.size > 1 and x[0] > x[-1]:
                x, tot = x[::-1], tot[::-1]
            y = np.interp(xp, x, tot)     # kernel models simulate on their own axis
        h = fithealth.reassess_live(base, self.recipe, self.exp_amp, y, ppm=xp,
                                    window=self.recipe.get("fit_window_ppm"))
        if h is not self._health and h != self._health:
            self._health_apply(h)

    def _health_apply(self, h):
        self._health = h
        self.health_strip.set_health(h)
        if h.fitted and not h.stale:
            # the Report header documents the LAST FIT; live edits do not
            # rewrite it (the strip carries the live state)
            self.results_summary.setText(h.summary())
            self.results_summary.setToolTip(h.summary_tooltip())
        self._update_enabled()
        self._health_show()

    def _health_show(self, *_):
        on = (self.central_stack.currentWidget() is self.view
              and bool(self.recipe and self.recipe.get("sites"))
              and getattr(self, "actHealthStrip", None) is not None
              and self.actHealthStrip.isChecked())
        self.health_strip.setVisible(on)

    def _health_reset(self):
        self._health = None
        self._health_fit = None
        self._last_lmfit = None
        self.health_strip.set_health(None)
        self._update_enabled()
        self._health_show()

    def show_fit_health(self, *_):
        """Decomposition ▸ Fit health details (F7) and the pill click: every
        flag plus the analysis tools, anchored under the pill."""
        if self._health is None:
            self.statusBar().showMessage("run a fit first (F5)")
            return
        self.health_strip.show_details(
            [self.actCorr, self.actErrors, self.actMC, self.actChi2])

    def _health_focus_param(self, i: int, pname: str):
        self.lines_dock.show()
        self.lines_dock.raise_()
        self.lines_table.select_param(i, pname)

    def _health_show_residual(self):
        # QAction.setChecked emits toggled, not triggered, so the slot would
        # not run: trigger() (the sidebar's own route), guarded to never
        # toggle the residual OFF
        if not self.actResid.isChecked():
            self.actResid.trigger()
        self.view.show_residual = True

    def _fit_done(self, result, stop_mode: str = ""):
        self._active_fit_worker = None
        self.lines_table.btnFit.setEnabled(True)
        self.view.stop_fit_animation()          # the final model replaces the trail
        if stop_mode == "cancel":
            # the worker fitted a COPY, so self.recipe is still the pre-fit state
            self._progress_end(False)
            self.statusBar().showMessage("fit cancelled — parameters unchanged")
            return
        self._progress_end(True)
        self.recipe = result.recipe.to_dict()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        labels = [s.get("label") or s["model"] for s in self.recipe["sites"]]
        self._last_model = (np.asarray(result.x_ppm), np.asarray(result.y_fit))
        self.view.set_model(result.x_ppm, result.y_fit, result.per_site,
                            labels, self.hidden, self.exp_ppm, self.exp_amp)
        # populations first, so the verdict's population rule reads THIS
        # fit's rows (a failed quantify must not leave the previous fit's)
        self._last_quant = None
        self.run_quantify(show=False)
        # one verdict -- residual within the noise / structured, physical
        # values, degenerate pairs, bounds, covariance, populations -- shown
        # by the strip under the spectrum, the Report header and the status
        # bar; the three dialogs stay the detail views (fithealth.assess)
        h = self._health_from_result(result)
        self.lines_table.set_chi2(h.chi_text())
        self.report.setPlainText(result.report)
        self.statusBar().showMessage(
            ("fit stopped — kept the latest iteration values"
             if stop_mode == "stop" else "fit done") + h.status_suffix())
        self._remember_site_defaults()   # per-nucleus smart defaults for next time
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

    def copy_latex(self):
        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("no fit to tabulate")
            return
        from larmor import methods
        tex = methods.latex_table(self.recipe, self._last_quant,
                                  caption=(self.recipe.get("sample") or ""))
        QApplication.clipboard().setText(tex)
        self.statusBar().showMessage("LaTeX table copied to clipboard")

    def copy_methods(self):
        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("no fit to describe")
            return
        from larmor import methods
        QApplication.clipboard().setText(methods.methods_sentence(self.recipe))
        self.statusBar().showMessage("methods sentence copied to clipboard")

    def export_publication_bundle(self):
        """One click: figure (png/pdf/svg) + LaTeX table + CSV + methods sentence
        + report.md for the current fit, into a chosen folder."""
        if not self.recipe or not self.recipe.get("sites") or not len(self.exp_ppm):
            self.statusBar().showMessage("fit a spectrum first")
            return
        folder = QFileDialog.getExistingDirectory(
            self, "Publication bundle — choose an output folder", self._last_dir())
        if not folder:
            return
        folder = Path(folder)
        from larmor import methods, engine, figures
        from larmor.recipe import Recipe
        from larmor.desktop.plot import site_color

        if self._last_quant is None:
            self.run_quantify(show=False)
        (folder / "methods.txt").write_text(
            methods.methods_sentence(self.recipe), encoding="utf-8")
        (folder / "table.tex").write_text(
            methods.latex_table(self.recipe, self._last_quant,
                                caption=(self.recipe.get("sample") or "")),
            encoding="utf-8")
        # a model-overlay figure from the current fit
        rec = Recipe.from_dict(self.recipe)
        x, total, per = engine.simulate(rec, exp_ppm=self.exp_ppm)
        traces = [{"data": {"x": list(map(float, self.exp_ppm)),
                            "y": list(map(float, self.exp_amp))},
                   "label": "experiment", "color": "#333333"},
                  {"data": {"x": list(map(float, x)), "y": list(map(float, total))},
                   "label": "model", "color": "#d1495b"}]
        for i, ys in enumerate(per):
            s = rec.sites[i]
            traces.append({"data": {"x": list(map(float, x)),
                                    "y": list(map(float, ys))},
                           "label": s.label or s.model, "linestyle": "--",
                           "color": site_color(i)})
        spec = {"kind": "1d", "traces": traces, "style": "article",
                "title": self.recipe.get("sample") or "",
                "xlabel": figures.nucleus_xlabel(self.recipe.get("nucleus", ""))}
        try:
            figures.export(spec, folder / "figure", formats=("png", "pdf", "svg"),
                           dpi=300)
        except Exception as exc:  # noqa: BLE001
            self.statusBar().showMessage(f"figure export failed: {exc}")
        # report.md tying it together
        rows = (self._last_quant or {}).get("rows", [])
        md = [f"# {self.recipe.get('sample') or 'Fit'} — {self.recipe.get('nucleus','')}",
              "", methods.methods_sentence(self.recipe), "",
              "| site | position (ppm) | fraction (%) |", "|---|---|---|"]
        for r in rows:
            md.append(f"| {r.get('label')} | {r.get('position_ppm')} | "
                      f"{r.get('fraction_pct', 0):.1f} |")
        md += ["", "![figure](figure.png)", ""]
        (folder / "report.md").write_text("\n".join(md), encoding="utf-8")
        self.statusBar().showMessage(f"publication bundle written to {folder}")

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
        # the winning FitResult goes through the same path as a plain fit, so
        # an Auto Fit gets the same verdict, chi text and correlations
        self._last_quant = None
        self.run_quantify(show=False)
        if res.result is not None:
            h = self._health_from_result(res.result)
            self.lines_table.set_chi2(h.chi_text())
        else:                              # defensive: auto_fit always sets it
            self._health_reset()
            self.lines_table.set_chi2(f"RMSD {res.best_rmsd:.4f}")
        self.request_simulation()
        self.statusBar().showMessage(res.summary)
