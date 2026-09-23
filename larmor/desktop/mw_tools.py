"""Main-window mixin: the Tools / Analysis menu slots.

Each method opens a dialog with the current state or applies one small
recipe edit handed back by a dialog (constraint sets, glass protocol, static
CT and Herzfeld-Berger seeds go through ``_EditingMixin``; the referencing
correction, VOCS stitching, the DFT tensor import (``_magres_add_sites``:
sites + note + ``provenance['dft_import']``) and the batch / series / error
tools are here).
The dialogs themselves live in their own ``larmor.desktop.*_dialog`` modules.

Owned state: the kept-alive non-modal dialogs (``_qcpmg_dlg``,
``_qcpmg_batch_dlg``, ``_ref_audit``, ``_inventory_dlg``, ``_acq_dlg``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

from larmor.desktop.workers import _fit_tol
from larmor.recipe import Recipe


class _ToolsMixin:
    """Dialog launchers and their small recipe edits."""

    def show_chi2_map(self):
        """χ² surface over a chosen parameter pair (basin vs degenerate valley)."""
        if not self.recipe or not self.recipe.get("sites") or not len(self.exp_ppm):
            self.statusBar().showMessage("fit a spectrum first")
            return
        from larmor.chi2map import varying_params
        if len(varying_params(self.recipe)) < 2:
            self.statusBar().showMessage("need at least two free parameters")
            return
        from larmor.desktop.chi2map_dialog import Chi2MapDialog
        (x0, x1), _ = self.view.getPlotItem().getViewBox().viewRange()
        Chi2MapDialog(self, self.recipe, self.exp_ppm, self.exp_amp,
                      (max(x0, x1), min(x0, x1))).exec()

    def compare_with_saved_fit(self):
        """Load a reference fit (recipe/dmfit) and show a parameter diff table."""
        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("fit a spectrum first")
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Reference fit to compare against", self._last_dir(),
            "Fits (*.json *.fxmla *.fxml);;All (*)")
        if not path:
            return
        ref = self._recipe_model_from(path)
        if ref is None or not ref.get("sites"):
            self.statusBar().showMessage("that file has no fitted lines")
            return
        from larmor.desktop.diff_dialog import RecipeDiffDialog
        RecipeDiffDialog(self, self.recipe, ref, Path(path).name).exec()

    def apply_recipe(self, path: str):
        """Load a saved recipe and drop ITS lines onto the currently open data —
        a starting model reused from a similar sample (the fit then refines)."""
        if self.exp_ppm is None or not self.recipe:
            self.statusBar().showMessage("open a spectrum first")
            return
        d = self._recipe_model_from(path)
        if d is None:
            return
        if not d.get("sites"):
            self.statusBar().showMessage("that recipe has no lines")
            return
        cur = self.recipe.get("nucleus", "")
        if d.get("nucleus") and cur and d["nucleus"] != cur:
            if QMessageBox.question(
                    self, "Different nucleus",
                    f"That recipe is {d['nucleus']} but the open spectrum is "
                    f"{cur}. Apply its lines anyway?") != QMessageBox.Yes:
                return
        self.snapshot()
        self.recipe["sites"] = d["sites"]
        self.on_structure_changed()
        self._add_recent_recipe(path)
        self.statusBar().showMessage(
            f"applied {len(d['sites'])} line(s) from {Path(path).name} — now Fit")

    def open_integrals(self):
        from larmor.desktop.integrate_dialog import IntegralsDialog

        if not self.exp_ppm.size:
            self.statusBar().showMessage("open a 1D spectrum first")
            return
        IntegralsDialog(self, self.exp_ppm, self.exp_amp,
                        sfo_MHz=float((self.recipe or {}).get(
                            "larmor_frequency_MHz", 0.0) or 0.0)).exec()

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

    def edit_computing_params(self):
        from larmor.desktop.dialogs import ComputingParamsDialog

        if ComputingParamsDialog(self).exec():
            self.statusBar().showMessage(
                "computing parameters updated — kernels rebuild on the next fit")
            self.request_simulation()

    def edit_mqmas_f1_ref(self):
        """View / set / fix the MQMAS isotropic-axis (F1) reference offset — the
        alignment between mrsimulator's kernel and the experiment's F1 axis.
        Auto-fitted by default; fix it (dmfit-style) to hold your value."""
        from PySide6.QtWidgets import QInputDialog

        if not self.recipe:
            self.statusBar().showMessage("open a 2D MQMAS map first")
            return
        cur = float(self.recipe.get("mqmas_f1_ref_ppm", 0.0))
        val, ok = QInputDialog.getDouble(
            self, "MQMAS F1 reference",
            "Isotropic-axis (F1) reference offset [ppm].\n"
            "This aligns the model's F1 axis to the experiment's convention.\n\n"
            "Cancel = keep auto-fitting it; OK = hold it fixed at this value.",
            cur, -80.0, 80.0, 2)
        if ok:
            self.recipe["mqmas_f1_ref_ppm"] = float(val)
            self.recipe["mqmas_f1_ref_vary"] = False
            self.statusBar().showMessage(
                f"MQMAS F1 reference fixed at {val:+.2f} ppm — Fit to apply "
                "(uncheck by re-running with auto)")
        else:
            self.recipe["mqmas_f1_ref_vary"] = True
            self.statusBar().showMessage("MQMAS F1 reference set to auto-fit")

    def edit_fit_tol(self):
        """Set the global fit completion threshold — the % change in the residual
        stdev below which every fit (1D, 2D, co-fit, batch) stops. 0 = full
        precision (the solver's own tolerance)."""
        from PySide6.QtWidgets import QInputDialog

        cur = _fit_tol()
        val, ok = QInputDialog.getDouble(
            self, "Fit completion threshold",
            "Stop a fit once the residual stdev (sdev) changes by less than (%):\n"
            "default 0.1 % ≈ dmfit's 1.0e-3;  0 = fit to full precision.", cur,
            0.0, 50.0, 3)
        if ok:
            QSettings("LARMOR", "app").setValue("fitStdevPct", float(val))
            self.statusBar().showMessage(
                "fit completion threshold: "
                + ("full precision" if val <= 0 else f"Δσ < {val:g}%"))

    def save_constraint_set(self):
        """Capture the current model's links/bounds/fixes as a named, reusable
        constraint set (stored in QSettings, applyable to any model)."""
        import json
        from PySide6.QtWidgets import QInputDialog
        from larmor import constraint_library as clib

        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("no model to capture constraints from")
            return
        cset = clib.capture(Recipe.from_dict(self.recipe))
        name, ok = QInputDialog.getText(self, "Save constraints",
                                        f"Name for this set ({clib.describe(cset)}):")
        if not ok or not name.strip():
            return
        s = QSettings("LARMOR", "app")
        lib = json.loads(s.value("constraintLibrary", "{}") or "{}")
        lib[name.strip()] = cset
        s.setValue("constraintLibrary", json.dumps(lib))
        self.statusBar().showMessage(f"saved constraint set “{name.strip()}” "
                                     f"({clib.describe(cset)})")

    def apply_constraint_set(self):
        """Apply a saved constraint set to the current model."""
        import json
        from PySide6.QtWidgets import QInputDialog
        from larmor import constraint_library as clib

        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("open/build a model first")
            return
        lib = json.loads(QSettings("LARMOR", "app").value(
            "constraintLibrary", "{}") or "{}")
        if not lib:
            self.statusBar().showMessage("no saved constraint sets yet")
            return
        items = [f"{k}  ({clib.describe(v)})" for k, v in lib.items()]
        choice, ok = QInputDialog.getItem(self, "Apply constraints",
                                          "Constraint set:", items, 0, False)
        if not ok:
            return
        name = list(lib.keys())[items.index(choice)]
        self.snapshot()
        rec = Recipe.from_dict(self.recipe)
        applied = clib.apply(rec, lib[name])
        self.recipe = rec.to_dict()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self.on_structure_changed()
        self.statusBar().showMessage(f"applied “{name}” to {len(applied)} parameter(s)")

    def restrict_glass_protocol(self):
        """Decomposition ▸ Advanced ▸ Restrict around current values — apply
        Edén 2023 §8.3's restricted-range recommendation in one step:
        δiso ± a window around the current value for every site, a physical
        FWHM floor for analytic (spin-½ style) peaks, amplitudes left free.
        Restricted, never fixed — vary flags are untouched. Undoable."""
        from PySide6.QtWidgets import QInputDialog
        from larmor.constraints_util import restrict_glass_protocol as _restrict

        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("open/build a model first")
            return
        win, ok = QInputDialog.getDouble(
            self, "Restrict around current values (Edén 2023 §8.3)",
            "δiso window: ± ppm around each site's CURRENT shift\n"
            "(peak FWHM of Gauss/Lorentz-family sites also gets a ≥4 ppm\n"
            "amorphous floor; amplitudes stay free; nothing is fixed):",
            3.0, 0.1, 50.0, 1)
        if not ok:
            return
        self.snapshot()
        notes = _restrict(self.recipe["sites"], shift_window_ppm=float(win))
        if not notes:
            self.statusBar().showMessage(
                "nothing to restrict (all shifts/widths linked or pinned)")
            return
        self.lines_table.rebuild(self.recipe, self.hidden)
        self.statusBar().showMessage(
            f"restricted {len(notes)} parameter(s) — δiso ±{win:g} ppm around "
            "current values (Edén 2023 §8.3; Edit ▸ Undo to revert)")

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

    def show_correlations(self):
        from larmor.desktop.correlation_dialog import CorrelationDialog

        lm = getattr(self, "_last_lmfit", None)
        if lm is None:
            self.statusBar().showMessage("run a fit first to see correlations")
            return
        CorrelationDialog(self, lm).exec()

    def show_czjzek_dist(self):
        from larmor.desktop.czjzek_dist_dialog import CzjzekDistDialog

        if not (self.recipe and any(
                s.get("model") in ("czjzek", "czjzek_d", "czjzek_corr",
                                   "ext_czjzek")
                for s in self.recipe.get("sites", []))):
            self.statusBar().showMessage("no Czjzek sites in the current fit")
            return
        CzjzekDistDialog(self, self.recipe).exec()

    def open_qcpmg_batch_fields(self):
        """Many samples x several fields in one grid: drop the processed
        spectra in, supervise each window, extrapolate them all."""
        from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

        nuc = self.recipe.get("nucleus", "") if self.recipe else ""
        dlg = QcpmgBatchFieldsDialog(self, nuc)
        self._qcpmg_batch_dlg = dlg          # non-modal: keep it alive
        dlg.show(); dlg.raise_(); dlg.activateWindow()

    def open_qcpmg_fields(self):
        from larmor.desktop.qcpmg_fields_dialog import shared_fields_dialog

        cur = None
        if self.recipe and self.exp_ppm is not None and self.exp_ppm.size:
            cur = (self.recipe.get("larmor_frequency_MHz", 0.0),
                   np.asarray(self.exp_ppm), np.asarray(self.exp_amp))
        nuc = self.recipe.get("nucleus", "") if self.recipe else ""
        # one persistent, NON-modal instance: fields sent from QCPMG
        # processing sessions accumulate here until Compute
        dlg = shared_fields_dialog(self, nuc, cur)
        dlg.show(); dlg.raise_(); dlg.activateWindow()

    def open_staticct(self):
        """Tools > Read static pattern: C_Q, η, δiso from the positions of
        the two horns and the informative edge of a static CT powder
        pattern — a measurement, not a fit."""
        from larmor.desktop.staticct_dialog import StaticCtDialog
        from larmor.nuclei import all_isotopes

        if not self.exp_ppm.size:
            self.statusBar().showMessage("load a spectrum first")
            return
        nucleus = (self.recipe or {}).get("nucleus", "")
        lar = float((self.recipe or {}).get("larmor_frequency_MHz", 0.0) or 0.0)
        iso = next((i for i in all_isotopes() if i.symbol == nucleus), None)
        if iso is None or iso.spin < 1.0 or lar <= 0:
            QMessageBox.warning(
                self, "Read static pattern",
                "Needs a half-integer quadrupolar nucleus and a Larmor "
                "frequency — set them in Experiment parameters.")
            return
        dlg = StaticCtDialog(self, self.exp_ppm, self.exp_amp, nucleus,
                             iso.spin, lar)
        dlg.seed_site.connect(self._staticct_seed)
        dlg.exec()

    def open_herzfeld_berger(self):
        """Tools > Herzfeld-Berger sideband analysis: (zeta, eta) from the
        integrated intensities of a MAS sideband manifold -- a measurement
        that works where a full csa_mas lineshape fit is not yet trusted."""
        from larmor.desktop.herzfeld_berger_dialog import HerzfeldBergerDialog
        from larmor.nuclei import all_isotopes

        if not self.exp_ppm.size:
            self.statusBar().showMessage("load a spectrum first")
            return
        nucleus = (self.recipe or {}).get("nucleus", "")
        lar = float((self.recipe or {}).get("larmor_frequency_MHz", 0.0) or 0.0)
        nur = float((self.recipe or {}).get("spin_rate_Hz", 0.0) or 0.0)
        if lar <= 0 or nur <= 0:
            QMessageBox.warning(
                self, "Herzfeld–Berger sideband analysis",
                "Needs a MAS spin rate and a Larmor frequency — a static "
                "spectrum has no sidebands to analyse. Set them in the "
                "experiment parameters (double-click the header).")
            return
        iso = next((i for i in all_isotopes() if i.symbol == nucleus), None)
        dlg = HerzfeldBergerDialog(self, self.exp_ppm, self.exp_amp, nucleus,
                                   lar, nur, iso.spin if iso else 0.5)
        dlg.seed_site.connect(self._hb_seed)
        dlg.exec()

    def open_referencing_audit(self):
        """Tools > Referencing audit: compare every SR of a session with the
        value its 1H adamantane reference gives by indirect (Xi) referencing;
        export a TopSpin sr list, log old and new values, or re-reference the
        open spectrum. Non-modal, so the workbench stays usable."""
        from larmor.desktop.referencing_dialog import ReferencingAuditDialog
        from larmor.referencing import session_root

        start = ""
        if self.source_path and Path(self.source_path).exists():
            start = str(session_root(self.source_path))
        else:
            start = str(QSettings("LARMOR", "app").value("lastDir", "") or "")
        dlg = getattr(self, "_ref_audit", None)
        if dlg is None:
            dlg = ReferencingAuditDialog(self, start, self.source_path)
            dlg.apply_requested.connect(self._apply_sr_correction)
            self._ref_audit = dlg
        else:
            dlg.current_source = self.source_path
            if start and not dlg.folder.text():
                dlg.folder.setText(start)
        dlg.show()
        dlg.raise_()

    def open_acquisition_table(self, paths=None):
        """Tools > Experimental section: Table S1 and the Experimental
        paragraph of the Explorer's selected spectra (else the open EXPNO),
        with every parameter that varies across the set highlighted and
        printed as a range. Needs no fit. Non-modal, kept alive; the open
        fit's confirmed MAS rate and referencing audit count as evidence."""
        from larmor.desktop.acquisition_dialog import AcquisitionTableDialog
        from larmor.referencing import session_root

        paths = [p for p in (paths or self.explorer.selected_spectra()) if p]
        if not paths and self.source_path and Path(self.source_path).exists():
            paths = [self.source_path]
        hint = ""
        if self.source_path and Path(self.source_path).exists():
            hint = str(session_root(self.source_path))
        else:
            hint = str(QSettings("LARMOR", "app").value("lastDir", "") or "")
        spin = {}
        if self.recipe and self.source_path and self.recipe.get("source_kind") == "bruker":
            spin[self.source_path] = (self.recipe.get("spin_rate_Hz"),
                                      self.recipe.get("mas_uncertain"))
        dlg = getattr(self, "_acq_dlg", None)
        if dlg is None:
            dlg = AcquisitionTableDialog(self, paths, hint, spin, self.recipe)
            self._acq_dlg = dlg
        else:
            dlg.spin_rates.update(spin)
            dlg.recipe = self.recipe
            if hint and not dlg.session_hint:
                dlg.session_hint = hint
            dlg.add_paths(paths)
        dlg.show()
        dlg.raise_()

    def open_session_inventory(self, folder=None):
        """Tools > Session inventory: one month folder as a sample × nucleus
        grid with the production EXPNO pre-picked per block (the highest
        EXPNO with a pdata/1/1r, demoted for a short NS or a setup / failed
        title), title-vs-folder flags, and one-action hand-off of the picks
        to Batch fit or Sequential fit. Non-modal and kept alive; also reached
        from the Explorer's folder context menu, which passes the folder
        (a QAction passes a bool, which is ignored)."""
        from larmor.desktop.inventory_dialog import SessionInventoryDialog
        from larmor.inventory import inventory_root

        folder = folder if isinstance(folder, str) and folder else None
        if folder:
            start = str(inventory_root(folder))
        elif self.source_path and Path(self.source_path).exists():
            start = str(inventory_root(self.source_path))
        else:
            start = str(QSettings("LARMOR", "app").value("lastDir", "") or "")
        dlg = getattr(self, "_inventory_dlg", None)
        if dlg is None:
            dlg = SessionInventoryDialog(self, start)
            dlg.batch_requested.connect(self.run_batch_fit)
            dlg.seq_requested.connect(self.run_seq_fit)
            dlg.open_requested.connect(self.load_source)
            self._inventory_dlg = dlg
        elif start and (folder or not dlg.folder.text()):
            dlg.folder.setText(start)
        dlg.show()
        dlg.raise_()

    def _apply_sr_correction(self, expno_path: str, new_sr: float, old_sr: float,
                             note: str = ""):
        """Re-reference the OPEN spectrum to a corrected SR (from the audit):
        the axis moves by -(new - old) / SF, the recipe records the new SR and
        the provenance keeps the old one; the data folder is untouched."""
        from larmor.referencing import axis_shift_ppm

        if self.recipe is None or not self.exp_ppm.size or not self.source_path:
            self.statusBar().showMessage("open the spectrum to re-reference first")
            return
        src, target = Path(self.source_path), Path(expno_path)
        if src != target and target not in src.parents:
            self.statusBar().showMessage(
                f"the audit row is {target.name}, the open spectrum is {src.name} "
                "— open that spectrum first")
            return
        larmor = float(self.recipe.get("larmor_frequency_MHz", 0.0) or 0.0)
        if larmor <= 0:
            return
        self.snapshot(with_axis=True)
        d_ppm = axis_shift_ppm(old_sr, new_sr, larmor)
        self.exp_ppm = self.exp_ppm + d_ppm
        if self._proc_base is not None:
            self._proc_base = (self._proc_base[0] + d_ppm, self._proc_base[1])
        self.recipe["sr_hz"] = float(new_sr)
        prov = self.recipe.setdefault("provenance", {}) or {}
        prov["referencing"] = {"old_sr_hz": float(old_sr), "new_sr_hz": float(new_sr),
                               "axis_shift_ppm": float(d_ppm), "note": note}
        self.recipe["provenance"] = prov
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self._update_exp_label()
        self.request_simulation()
        self.statusBar().showMessage(
            f"re-referenced: SR {old_sr:.2f} → {new_sr:.2f} Hz, axis moved "
            f"{d_ppm:+.3f} ppm — {note}" if note else
            f"re-referenced: SR {old_sr:.2f} → {new_sr:.2f} Hz, axis moved {d_ppm:+.3f} ppm",
            12000)

    def open_vocs(self):
        """Tools > Stitch frequency-stepped (VOCS): sub-spectra acquired at
        stepped transmitter offsets combine into one wideline pattern."""
        from larmor.desktop.vocs_dialog import VocsDialog

        dlg = VocsDialog(self)
        dlg.applied.connect(self._vocs_to_workbench)
        dlg.exec()

    def _vocs_to_workbench(self, ppm, amp, notes, sources):
        first = sources[0] if sources else ""
        name = Path(first).parent.name or "VOCS"
        self._display_1d(np.asarray(ppm, float), np.asarray(amp, float),
                         (self.recipe or {}).get("nucleus", ""),
                         (self.recipe or {}).get("larmor_frequency_MHz", 0.0),
                         0.0, f"{name} (VOCS stitch)", first)
        # a stitched pattern is a DERIVED spectrum: record how it was made
        self.recipe["notes"] = list(notes)
        self.recipe["provenance"] = {
            "vocs_sources": [str(s_) for s_ in sources],
            "vocs_note": notes[0] if notes else ""}
        # stitched wideline patterns are static experiments in practice, but
        # the pieces' metadata was not inspected here -- leave the rate as
        # the experiment dialog's problem and flag it for confirmation
        self.recipe["mas_uncertain"] = True
        self._update_exp_label()
        self.statusBar().showMessage(
            f"stitched {len(sources)} sub-spectra — set nucleus/field/rate "
            "in Experiment parameters, then fit")

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
        # NON-modal, and kept alive on self: QCPMG processing is a long
        # session, and exec() froze the whole application behind it -- including
        # the infinite-field window it can open
        dlg = QcpmgDialog(self, source)
        dlg.accepted_1d.connect(self._fid_to_workbench)
        self._qcpmg_dlg = dlg
        dlg.show(); dlg.raise_(); dlg.activateWindow()

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
        """Tools > Import DFT tensors (.magres): spin-aware seeding from
        computed tensors, one site per crystallographic position, shieldings
        converted through a fitted calibration line (or σ_ref). The dialog
        opens without a spectrum too -- its Calibration tab is useful alone;
        Add stays disabled until a spectrum is open."""
        from larmor.desktop.magres_dialog import MagresDialog

        rec = self.recipe or {}
        exp_max = (float(np.max(np.abs(self.exp_amp)))
                   if self.exp_amp is not None and self.exp_amp.size else 0.0)
        dlg = MagresDialog(self, nucleus=rec.get("nucleus", "") or "",
                           recipe=self.recipe,
                           spin_rate_Hz=float(rec.get("spin_rate_Hz", 0.0) or 0.0),
                           exp_max=exp_max)
        if dlg.exec() and dlg.result is not None:
            self._magres_add_sites(dlg.result)

    def _magres_add_sites(self, imp):
        """Apply a MagresImport: one undo snapshot, the master amplitude
        seeded from the spectrum, the sites appended as-is (model / label /
        params only), one recipe note, ``provenance['dft_import']`` and a
        status line naming the file, the conversion and the lock ratios."""
        if self.recipe is None:
            self.statusBar().showMessage("load a spectrum first")
            return
        self.snapshot()
        sites = self.recipe.setdefault("sites", [])
        exp_max = (float(np.max(np.abs(self.exp_amp)))
                   if self.exp_amp is not None and self.exp_amp.size else 0.0)
        if imp.sites and exp_max > 0:
            master = imp.sites[0]["params"]
            amp = master.get("amplitude")
            if amp is not None and not amp.get("expr"):
                height = exp_max / max(int(imp.n_groups), 1)
                if imp.area_amplitude:
                    fwhm = float(master.get("shift_fwhm_ppm", {}).get("value", 1.0)
                                 or 1.0)
                    height *= fwhm * 1.064          # area of a unit-height Gaussian
                amp["value"] = float(height)
        for sd in imp.sites:
            sites.append(sd)
        self.recipe.setdefault("notes", []).append(imp.note)
        prov = dict(self.recipe.get("provenance") or {})
        prov["dft_import"] = imp.provenance
        self.recipe["provenance"] = prov
        self.on_structure_changed()
        msg = (f"added {len(imp.sites)} site(s) from "
               f"{Path(imp.provenance.get('file', '')).name} "
               f"(δ via {imp.conversion})")
        if imp.lock_amplitude and len(imp.ratios) > 1:
            msg += "; populations locked " + ":".join(str(r) for r in imp.ratios)
        if imp.share_width and len(imp.sites) > 1:
            msg += "; one shared linewidth"
        self.statusBar().showMessage(msg, 12000)

    def open_twod(self):
        from larmor.desktop.twod_dialog import TwoDDialog

        expno = self.source_path if (self.source_path and
                                     Path(self.source_path).is_dir()) else None
        TwoDDialog(self, expno).exec()

    def run_errors_analysis(self):
        if not self.recipe or not self.recipe["sites"]:
            return
        from larmor.desktop.tool_dialogs import ErrorsDialog

        ErrorsDialog(self, self.recipe, self.exp_ppm, self.exp_amp,
                     self.view.current_xrange()).exec()

    def run_monte_carlo(self):
        if not self.recipe or not self.recipe["sites"]:
            self.statusBar().showMessage("fit a model first")
            return
        if self.exp_ppm is None:
            self.statusBar().showMessage("open a spectrum first")
            return
        from larmor.desktop.montecarlo_dialog import MonteCarloDialog

        MonteCarloDialog(self, self.recipe, self.exp_ppm, self.exp_amp,
                         self.view.current_xrange()).exec()

    def run_batch_report(self):
        from larmor.desktop.batch_dialog import BatchReportDialog

        start = str(QSettings("LARMOR", "app").value("lastDir", "") or "")
        BatchReportDialog(self, start).exec()

    def run_batch_fit(self, paths):
        """Batch-fit several spectra (from the Explorer) with one shared model."""
        paths = [p for p in (paths or []) if p]
        if len(paths) < 2:
            self.statusBar().showMessage(
                "Ctrl/Shift-select at least two spectra in the Explorer first")
            return
        from larmor.desktop.batchfit_dialog import BatchFitDialog

        model = self.recipe if (self.recipe and self.recipe.get("sites")) else None
        dlg = BatchFitDialog(self, paths, model)
        dlg.exec()
        self._remember_batch(dlg)

    def run_seq_fit(self, paths=None):
        """Sequential (forward-backward) fit of the Explorer-selected series."""
        paths = [p for p in (paths or self.explorer.selected_spectra()) if p]
        if len(paths) < 2:
            self.statusBar().showMessage(
                "Ctrl/Shift-select at least two spectra in the Explorer first")
            return
        from larmor.desktop.seqfit_dialog import SeqFitDialog
        model = self.recipe if (self.recipe and self.recipe.get("sites")) else None
        SeqFitDialog(self, paths, model).exec()

    def open_plotting_studio(self, spec=None, *, ws_index=None):
        """Open the Plotting studio, optionally seeded with a figure spec.

        The studio is modal, so "saving figures" means record-on-close: the
        figure is kept in the Workspaces dock when the studio closes with a
        spec that differs from the one it opened with (a previewed-and-closed
        studio leaves nothing behind); a figure reopened from the dock
        (``ws_index``) is always written back to its row."""
        from larmor.desktop.plotting_studio import PlottingStudio
        # QAction.triggered passes a bool -- only a dict is a seed
        dlg = PlottingStudio(self, spec if isinstance(spec, dict) else None)
        seed = dlg._spec()
        dlg.exec()
        cur = dlg._spec()
        if ws_index is not None or cur != seed:
            self._remember_figure(cur, ws_index)

    def plot_current_spectrum(self):
        """Seed the Plotting studio with the current spectrum (+ model if fitted)."""
        traces = []
        if self.exp_ppm is not None and len(self.exp_ppm):
            traces.append({"data": {"x": list(map(float, self.exp_ppm)),
                                    "y": list(map(float, self.exp_amp))},
                           "label": "experiment"})
        spec = {"kind": "1d", "traces": traces}
        self.open_plotting_studio(spec)     # untouched seed -> nothing kept
