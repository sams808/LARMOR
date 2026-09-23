"""Main-window mixin: editing the model.

Undo / redo (``snapshot``, ``_capture_state``, ``_restore_state``), the
add-line mode, the per-nucleus seeding table ``_NUCLEUS_START`` and its
helpers, every ``add_*`` entry point (sites, function lines, autopick,
literature labels, background and current-spectrum components, 2D sites),
the seeds handed back by the static-CT and Herzfeld-Berger dialogs, site
structure edits (move / delete with constraint remapping) and the paddles.

Owned state: ``undo_stack`` / ``redo_stack``, ``_paddle_live``, the add mode
held on the view.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QCheckBox, QFileDialog, QLabel, QMessageBox

from larmor import cellparse, models as model_registry
from larmor.desktop.mw_files import _load_any


class _EditingMixin:
    """Undo / redo, adding and restructuring sites, paddles."""

    def new_fit(self):
        if self.recipe is None:
            return
        self.snapshot()
        self.recipe["sites"] = []
        self.on_structure_changed()

    # ------------------------------------------------------------- undo
    @staticmethod
    def _recipe_copy(rec):
        """Deep copy of a recipe dict that SHARES each site's "ref" trace.

        A "spectrum" component carries its whole reference spectrum in
        site["ref"] (16k+ points). json.dumps of that on EVERY user action --
        twice, once for undo and once for the session -- measured 29 ms per
        click and ~40 MB across a 60-deep undo stack. The trace is immutable
        by contract (created once in _append_spectrum_site, read by
        engine._render_spectrum, never edited in place -- replace it, never
        mutate it), so undo states can all point at the same list."""
        if not rec:
            return rec
        head = {k: v for k, v in rec.items() if k != "sites"}
        out = json.loads(json.dumps(head))
        out["sites"] = []
        for site in rec.get("sites", []):
            sc = json.loads(json.dumps(
                {k: v for k, v in site.items() if k != "ref"}))
            if "ref" in site:
                sc["ref"] = site["ref"]              # shared, immutable
            out["sites"].append(sc)
        return out

    def _capture_state(self, with_axis=False) -> dict:
        snap = {"recipe": self._recipe_copy(self.recipe)}
        if with_axis and self.exp_ppm is not None and len(self.exp_ppm):
            snap["ppm"] = np.array(self.exp_ppm, float)
            snap["amp"] = np.array(self.exp_amp, float)
            snap["base"] = self._proc_base
        return snap

    def _restore_state(self, snap):
        if isinstance(snap, str):                 # legacy recipe-only snapshot
            self.recipe = json.loads(snap)
            return
        rec = snap["recipe"]
        # copy again on the way OUT, so later edits to the live recipe can
        # never reach the states still sitting in the undo/redo stacks
        self.recipe = (json.loads(rec) if isinstance(rec, str)
                       else self._recipe_copy(rec))
        if "ppm" in snap:                         # calibrate/SR changed the axis
            self.exp_ppm = snap["ppm"]
            self.exp_amp = snap["amp"]
            self._proc_base = snap["base"]
            self.view.set_experiment(self.exp_ppm, self.exp_amp)
            self._update_exp_label()
            # the p0/p1 controls follow the restored recipe (a drag-to-phase
            # undo would otherwise leave the undone numbers in the panel and
            # the next slider nudge would silently re-apply them); silent,
            # and a no-op for calibrate / SR / 2-point snapshots
            self.proc_panel.sync_phase_from(
                self.recipe.get("processing") if self.recipe else None)

    def snapshot(self, with_axis=False):
        """Push an undo state. ``with_axis=True`` also captures the experiment
        axis so an axis-changing op (calibrate / SR) is fully reversible."""
        if self.recipe is None:
            return
        self.undo_stack.append(self._capture_state(with_axis))
        if len(self.undo_stack) > 60:
            self.undo_stack.pop(0)
        self.redo_stack.clear()
        self._update_enabled()
        self._persist_session()

    def undo(self):
        if not self.undo_stack:
            return
        tgt = self.undo_stack.pop()
        self.redo_stack.append(
            self._capture_state(with_axis=isinstance(tgt, dict) and "ppm" in tgt))
        self._restore_state(tgt)
        self.on_structure_changed()

    def redo(self):
        if not self.redo_stack:
            return
        tgt = self.redo_stack.pop()
        self.undo_stack.append(
            self._capture_state(with_axis=isinstance(tgt, dict) and "ppm" in tgt))
        self._restore_state(tgt)
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

    #: sensible starting quadrupolar parameters by nucleus (glasses/materials);
    #: only a starting point — always refined by the fit
    _NUCLEUS_START = {
        "27Al": {"sigma_Cq_MHz": 1.5, "Cq_MHz": 3.0, "eta": 0.6, "shift_fwhm_ppm": 12.0},
        "11B":  {"sigma_Cq_MHz": 1.0, "Cq_MHz": 2.6, "eta": 0.2, "shift_fwhm_ppm": 3.0},
        "23Na": {"sigma_Cq_MHz": 1.0, "Cq_MHz": 2.0, "eta": 0.6, "shift_fwhm_ppm": 8.0},
        "17O":  {"sigma_Cq_MHz": 0.9, "Cq_MHz": 3.5, "eta": 0.3, "shift_fwhm_ppm": 15.0},
        "35Cl": {"sigma_Cq_MHz": 1.6, "Cq_MHz": 3.0, "eta": 0.7, "shift_fwhm_ppm": 30.0},
        "71Ga": {"sigma_Cq_MHz": 2.0, "Cq_MHz": 5.0, "eta": 0.5, "shift_fwhm_ppm": 15.0},
        "93Nb": {"sigma_Cq_MHz": 2.0, "Cq_MHz": 4.0, "eta": 0.5, "shift_fwhm_ppm": 20.0},
    }

    def _seed_nucleus_defaults(self, name: str, params: dict):
        """Pre-fill starting values for a new site from the nucleus:
        first what the user last fitted for this nucleus+model, else the built-in
        quadrupolar starting table."""
        nucleus = (self.recipe or {}).get("nucleus", "")
        # 1) what worked last time for this nucleus + model (remembered on fit)
        remembered = self._remembered_site_defaults(nucleus, name)
        for k, val in (remembered or {}).items():
            if k in params:
                params[k]["value"] = val
        # 2) built-in quadrupolar starting points (only if not remembered)
        if name in ("czjzek", "czjzek_d", "czjzek_corr", "ext_czjzek",
                    "quad_ct", "quad_first", "quad_csa", "csa_czjzek"):
            for k, val in (self._NUCLEUS_START.get(nucleus) or {}).items():
                if k in params and not (remembered and k in remembered):
                    params[k]["value"] = val

    def _seed_from_data(self, name: str, params: dict, centre_ppm: float):
        """Measure the starting width from the SPECTRUM. A per-nucleus table
        cannot cover every isotope: for 81Br the defaults gave a 19 kHz
        needle on a 400 kHz pattern, which the fit has no gradient to escape
        from. Skipped when the user's own remembered values already apply."""
        if self.exp_ppm is None or not len(self.exp_ppm):
            return
        nucleus = (self.recipe or {}).get("nucleus", "")
        if self._remembered_site_defaults(nucleus, name):
            return                     # what worked last time wins
        try:
            from larmor import estimate
            vals = estimate.start_values(
                name, self.exp_ppm, self.exp_amp, nucleus,
                float((self.recipe or {}).get("larmor_frequency_MHz", 0.0) or 0.0),
                centre_ppm=centre_ppm,
                spin_rate_Hz=float(
                    (self.recipe or {}).get("spin_rate_Hz", 0.0) or 0.0))
        except Exception:                                  # noqa: BLE001
            return
        for k, v in vals.items():
            if k in params:
                params[k]["value"] = float(v)

    @staticmethod
    def _remembered_site_defaults(nucleus: str, model: str) -> dict:
        import json
        if not nucleus:
            return {}
        lib = json.loads(QSettings("LARMOR", "app").value("siteDefaults", "{}")
                         or "{}")
        return lib.get(f"{nucleus}/{model}", {})

    def _remember_site_defaults(self):
        """After a fit, remember each site's shape parameters (not position or
        amplitude) per nucleus+model, so the next new site starts from them."""
        import json
        if os.environ.get("LARMOR_NO_SESSION"):     # don't pollute settings in tests
            return
        nucleus = (self.recipe or {}).get("nucleus", "")
        if not nucleus or not self.recipe.get("sites"):
            return
        lib = json.loads(QSettings("LARMOR", "app").value("siteDefaults", "{}")
                         or "{}")
        for s in self.recipe["sites"]:
            keep = {}
            for pn, p in s.get("params", {}).items():
                if pn in ("amplitude", "isotropic_chemical_shift_ppm", "gl"):
                    continue
                v = p.get("value") if isinstance(p, dict) else p
                if v is not None:
                    keep[pn] = float(v)
            if keep:
                lib[f"{nucleus}/{s.get('model')}"] = keep
        QSettings("LARMOR", "app").setValue("siteDefaults", json.dumps(lib))

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
        self._seed_nucleus_defaults(name, params)
        self._seed_from_data(name, params, ppm)
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

    def label_from_literature(self):
        """Decomposition > Label lines from literature ranges: name every
        auto-labelled line after the literature species whose band holds its
        position (View > Literature shift ranges data); user-typed labels are
        kept and reported."""
        from larmor import refranges

        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("add or fit lines first")
            return
        nucleus = self.recipe.get("nucleus", "")
        if not (refranges.ranges_for(nucleus) or refranges.positions_for(nucleus)):
            self.statusBar().showMessage(
                f"no literature ranges compiled for {nucleus or '?'} — see "
                "Help ▸ Literature shift ranges")
            return
        changed, kept, unmatched = [], [], []
        self.snapshot()
        for i, s in enumerate(self.recipe["sites"]):
            p = s.get("params", {}).get("isotropic_chemical_shift_ppm")
            if not p:
                continue
            hit = refranges.assign(nucleus, float(p["value"]))
            letter = cellparse.index_to_letter(i)
            if hit is None:
                unmatched.append(letter)
            elif refranges.is_auto_label(s.get("label")):
                s["label"] = hit["label"]
                changed.append(f"{letter}={hit['label']}")
            else:
                kept.append(f"{letter} ({s['label']} → {hit['label']}?)")
        self.on_structure_changed()
        msg = (f"labelled {len(changed)}: " + ", ".join(changed)) if changed \
            else "no auto-labelled line falls in a literature band"
        if kept:
            msg += "  ·  kept your own: " + ", ".join(kept)
        if unmatched:
            msg += "  ·  outside every band: " + ", ".join(unmatched)
        self.statusBar().showMessage(msg, 12000)

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
            # load_any returns (ppm, amp, recipe, meta, warnings) -- unpacking it
            # as (recipe, ppm, amp, ...) put the recipe DICT into `amp` and threw
            # a TypeError out of the slot, so no background could ever be added
            ppm, amp, *_ = _load_any(path)
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
        try:
            ppm = np.asarray(ppm, float).ravel()
            amp = np.asarray(amp, float).ravel()
        except (TypeError, ValueError) as exc:
            QMessageBox.warning(self, "Background spectrum", f"cannot read: {exc}")
            return
        if ppm.size < 2 or amp.size != ppm.size:
            QMessageBox.warning(self, "Background spectrum",
                                "that file does not hold a 1D spectrum")
            return
        self.snapshot()
        # start the amplitude near the data's peak so it is on-scale
        self._append_spectrum_site(
            ppm, amp, f"bg-{len(self.recipe['sites'])}",
            amplitude=(float(np.max(np.abs(self.exp_amp)))
                       if self.exp_amp.size else 1.0))
        self.statusBar().showMessage(
            f"added background spectrum '{Path(path).name}' — Fit scales and "
            "shifts it")

    def _append_spectrum_site(self, ppm, amp, label: str, amplitude: float,
                              shift_ppm: float = 0.0, shift_vary: bool = True):
        """Append a "spectrum" component carrying (ppm, amp) as its reference
        trace, normalised to unit peak. Shared by BOTH entry points -- a file
        on disk and the spectrum currently on screen -- so the two can never
        drift apart in how the site is built."""
        ppm = np.asarray(ppm, float)
        amp = np.asarray(amp, float)
        amp = amp / (np.max(np.abs(amp)) or 1.0)         # unit peak
        m = model_registry.get("spectrum")
        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        params["amplitude"]["value"] = float(amplitude)
        params["shift_ppm"]["value"] = float(shift_ppm)
        params["shift_ppm"]["vary"] = bool(shift_vary)
        self.recipe["sites"].append({
            "model": "spectrum", "label": label,
            "ref": {"ppm": ppm.tolist(), "amp": amp.tolist()},
            "params": params})
        self.on_structure_changed()

    def add_current_spectrum_line(self):
        """Add the spectrum CURRENTLY ON SCREEN as a fit component: a rigidly
        shifted copy of the data itself.

        That is how a satellite-transition or spinning-sideband manifold is
        removed without parameterising a lineshape for it -- the manifold
        repeats the whole pattern at ±νrot, so a shifted copy of the measured
        spectrum models it directly.

        Unlike "Add background spectrum…" this takes the PROCESSED trace as
        displayed (baseline, phase, ppm referencing, QCPMG output). Re-opening
        the file from disk would give the unprocessed one, which is not what
        has to cancel. The copy is a snapshot: re-processing the workbench
        afterwards does not update it.
        """
        from PySide6.QtWidgets import (QComboBox, QDialog, QDialogButtonBox,
                                       QDoubleSpinBox, QFormLayout)

        if self.recipe is None or self.exp_amp.size < 2:
            self.statusBar().showMessage("open a spectrum first")
            return
        nu_r = float(self.recipe.get("spin_rate_Hz", 0) or 0)
        lar = float(self.recipe.get("larmor_frequency_MHz", 0) or 0)
        nur_ppm = (nu_r / lar) if (nu_r > 0 and lar > 0) else 0.0

        dlg = QDialog(self)
        dlg.setWindowTitle("Add a copy of this spectrum")
        form = QFormLayout(dlg)
        combo = QComboBox()
        if nur_ppm:
            combo.addItem(f"+νrot   ({nur_ppm:+.4g} ppm)", nur_ppm)
            combo.addItem(f"−νrot   ({-nur_ppm:+.4g} ppm)", -nur_ppm)
        combo.addItem("custom", None)
        form.addRow("shift the copy by", combo)
        span = float(np.ptp(self.exp_ppm)) if self.exp_ppm.size else 1e4
        spin = QDoubleSpinBox()
        spin.setRange(-10.0 * max(span, 1.0), 10.0 * max(span, 1.0))
        spin.setDecimals(4)
        spin.setSuffix(" ppm")
        spin.setValue(nur_ppm)
        form.addRow("shift", spin)
        combo.currentIndexChanged.connect(
            lambda i: (spin.setValue(float(combo.itemData(i)))
                       if combo.itemData(i) is not None else None))

        hold = QCheckBox("hold the shift fixed while fitting")
        hold.setChecked(True)
        hold.setToolTip(
            "A copy of the SAME data is degenerate with the rest of the model "
            "when its amplitude AND its shift are both free — the copy alone "
            "reproduces the whole spectrum at amplitude 1, shift 0. Holding "
            "the shift at the sideband spacing keeps it a sideband component.")
        form.addRow(hold)
        form.addRow(QLabel(
            f"νrot = {nur_ppm:.4g} ppm  ({nu_r / 1000:.1f} kHz at {lar:.2f} MHz)"
            if nur_ppm else
            "no MAS rate on this dataset — the shift below is used as typed"))
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        form.addRow(bb)
        if dlg.exec() != QDialog.Accepted:
            return

        shift = float(spin.value())
        fixed = hold.isChecked()
        self.snapshot()
        # a satellite/sideband manifold is a FRACTION of the centre band, so
        # start well below the data's peak rather than at it (as the on-disk
        # background path does, where the reference is a different spectrum)
        self._append_spectrum_site(
            self.exp_ppm, self.exp_amp, f"copy-{len(self.recipe['sites'])}",
            amplitude=0.1 * float(np.max(np.abs(self.exp_amp))),
            shift_ppm=shift, shift_vary=not fixed)
        self.statusBar().showMessage(
            f"added a copy of this spectrum shifted by {shift:+.4g} ppm — Fit "
            + ("scales it (shift held)" if fixed else "scales and re-shifts it"))

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
            # remap the hidden set and every constraint (expr) that referenced a
            # site by index, so a link like "E = D+5.3" never silently becomes a
            # self-reference "E = E+5.3" (which recurses forever at fit time)
            self.hidden = {i - 1 if i > idx else i for i in self.hidden}
            dropped = self._remap_exprs_after_delete(idx)
            if dropped:
                self.statusBar().showMessage(
                    "removed a line — dropped now-invalid constraint(s): "
                    + ", ".join(dropped))
        elif action == "duplicate":
            self.snapshot()
            copy = json.loads(json.dumps(self.recipe["sites"][idx]))
            copy["label"] = (copy.get("label") or "line") + "-copy"
            for p in copy["params"].values():
                p["stderr"] = None
            self.recipe["sites"].append(copy)
        elif action in ("move_up", "move_down"):
            self._move_site(idx, -1 if action == "move_up" else +1)
            return
        elif action == "visibility":
            (self.hidden.discard(idx) if idx in self.hidden
             else self.hidden.add(idx))
        self.on_structure_changed()

    def _move_site(self, idx: int, delta: int):
        """Reorder a line in the table (user comfort), remapping every constraint
        reference and the hidden set so links keep pointing at the right lines."""
        from larmor.constraints_util import remap_exprs_after_move
        sites = self.recipe["sites"]
        j = idx + delta
        if not (0 <= idx < len(sites) and 0 <= j < len(sites)):
            return
        self.snapshot()
        sites[idx], sites[j] = sites[j], sites[idx]
        old_to_new = {i: i for i in range(len(sites))}
        old_to_new[idx], old_to_new[j] = j, idx
        remap_exprs_after_move(sites, old_to_new)
        self.hidden = {old_to_new[i] for i in self.hidden}
        self.on_structure_changed()

    def _remap_exprs_after_delete(self, deleted_idx: int) -> list:
        """After removing site ``deleted_idx``, fix every remaining constraint's
        ``s<k>.param`` references: drop any that pointed at the deleted site or
        would become a self-reference, and shift ``k>deleted_idx`` down by one.
        Returns the labels of dropped constraints."""
        from larmor.constraints_util import remap_exprs_after_delete
        return remap_exprs_after_delete(self.recipe.get("sites", []), deleted_idx)

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
        self._health_live(full=True)          # the one full pass after a drag
        self._persist_session()

    def _staticct_seed(self, cq_MHz: float, eta: float, diso_ppm: float):
        """Turn a static-pattern reading into a quad_ct starting site."""
        if self.recipe is None:
            return
        self.snapshot()
        m = model_registry.get("quad_ct")
        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        params["isotropic_chemical_shift_ppm"]["value"] = float(diso_ppm)
        params["Cq_MHz"]["value"] = float(cq_MHz)
        params["eta"]["value"] = float(eta)
        if self.exp_amp.size:
            params["amplitude"]["value"] = float(np.max(np.abs(self.exp_amp)))
        n = len(self.recipe["sites"])
        self.recipe["sites"].append({"model": "quad_ct",
                                     "label": f"read-{n}", "params": params})
        self.on_structure_changed()
        self.statusBar().showMessage(
            f"added quad_ct from the reading: C_Q {cq_MHz:.2f} MHz, "
            f"η {eta:.2f}, δiso {diso_ppm:.1f} ppm — refine with Fit if "
            "the lineshape supports it")

    def _hb_seed(self, zeta_ppm: float, eta: float, diso_ppm: float,
                 amplitude: float):
        """Turn a Herzfeld-Berger reading into a csa_mas starting site."""
        if self.recipe is None:
            return
        self.snapshot()
        m = model_registry.get("csa_mas")
        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        params["isotropic_chemical_shift_ppm"]["value"] = float(diso_ppm)
        params["zeta_ppm"]["value"] = float(zeta_ppm)
        params["eta"]["value"] = float(eta)
        params["amplitude"]["value"] = float(amplitude)
        n = len(self.recipe["sites"])
        self.recipe["sites"].append({"model": "csa_mas", "label": f"HB-{n}",
                                     "params": params})
        self.on_structure_changed()
        self.statusBar().showMessage(
            f"added csa_mas from the sideband reading: ζ {zeta_ppm:.1f} ppm, "
            f"η {eta:.2f}, δiso {diso_ppm:.2f} ppm — refine with Fit")
