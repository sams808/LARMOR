"""Main-window mixin: spinning sidebands.

The manual ``add_sidebands`` dialog and the auto-detect path (F3): the
detector run after a load or a processing change, the offer banner over the
plot, the user's choice (manifold / copies / model line / dismiss), the
rate update and the remembered dismissal. Physics lives in
``larmor.sidebands``; the banner widget in ``larmor.desktop.sideband_offer``.

Owned state: ``ssb_banner`` (the offer banner), ``_ssb_detection`` (the last
detection) and ``_ssb_dismissed``; the ``ssbAutoOffer`` QSettings flag behind
the menu toggle.
"""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QLabel, QMessageBox

from larmor import models as model_registry


class _SidebandsMixin:
    """Manual and auto-detected spinning-sideband handling."""

    # -------- spinning sidebands: clone a line at pos ± n·νrot ---------------
    def add_sidebands(self):
        """Add lines at ± n·νrot of a fitted line (same parameters, shifted),
        for models that do not already generate their own sideband manifold."""
        import copy

        from PySide6.QtWidgets import (
            QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout,
            QSpinBox)

        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("add or fit a line first")
            return
        nu_r = float(self.recipe.get("spin_rate_Hz", 0) or 0)
        lar = float(self.recipe.get("larmor_frequency_MHz", 0) or 0)
        if nu_r <= 0 or lar <= 0:
            QMessageBox.warning(
                self, "Spinning sidebands",
                "A MAS spin rate and Larmor frequency are needed to place "
                "sidebands — set them in the experiment parameters (double-click "
                "the header).")
            return
        nur_ppm = nu_r / lar
        skip = {"sidebands", "csa_mas", "csa_czjzek", "quad_first"}
        eligible = [(i, s) for i, s in enumerate(self.recipe["sites"])
                    if s.get("model") not in skip]
        if not eligible:
            QMessageBox.information(
                self, "Spinning sidebands",
                "No eligible lines — the present models already include their own "
                "sidebands (use those, or add a Gauss/Lor, Czjzek, … line).")
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Add spinning sidebands")
        form = QFormLayout(dlg)
        combo = QComboBox()
        for i, s in eligible:
            combo.addItem(f"s{i} — {s.get('label') or s['model']}", i)
        # the fit table's line menu preselects its line (add_sidebands_for_line)
        pre = getattr(self, "_ssb_preselect", None)
        self._ssb_preselect = None
        if pre is not None and combo.findData(pre) >= 0:
            combo.setCurrentIndex(combo.findData(pre))
        form.addRow("line", combo)
        fwd = QSpinBox(); fwd.setRange(0, 20); fwd.setValue(1)
        form.addRow("sidebands forward (+νrot)", fwd)
        bwd = QSpinBox(); bwd.setRange(0, 20); bwd.setValue(1)
        form.addRow("sidebands backward (−νrot)", bwd)
        form.addRow(QLabel(f"spacing νrot = {nur_ppm:.3g} ppm  "
                           f"({nu_r / 1000:.1f} kHz)"))
        link = QCheckBox("link to the parent line: position = parent ± k·νrot, "
                         "every shape parameter equal to the parent's; only "
                         "the amplitude is free")
        link.setChecked(True)
        link.setToolTip("a sideband IS the parent line displaced by k·νrot: "
                        "fitting its position and width independently only "
                        "adds parameters the data cannot tell apart. Untick "
                        "for free copies (the previous behaviour).")
        form.addRow(link)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject)
        form.addRow(bb)
        if dlg.exec() != QDialog.Accepted:
            return

        base_i = combo.currentData()
        base = self.recipe["sites"][base_i]
        bname = base.get("label") or base["model"]
        self.snapshot()
        added = []
        for sign, cnt in ((1, fwd.value()), (-1, bwd.value())):
            for k in range(1, cnt + 1):
                s = copy.deepcopy(base)
                s["label"] = f"{bname}{sign * k:+d}sb"
                pos = s["params"]["isotropic_chemical_shift_ppm"]
                pos["value"] = float(pos["value"] + sign * k * nur_ppm)
                pos["stderr"] = None
                if link.isChecked():
                    for pname, p in s["params"].items():
                        p["stderr"] = None
                        if pname == "amplitude":
                            p["expr"] = None
                            p["vary"] = True
                            p["value"] = float(p["value"]) * 0.3 ** k
                        elif pname == "isotropic_chemical_shift_ppm":
                            p["expr"] = (f"s{base_i}.isotropic_chemical_shift_ppm"
                                         f" {'+' if sign > 0 else '-'} "
                                         f"{k * nur_ppm:.6g}")
                        else:
                            p["expr"] = f"s{base_i}.{pname}"
                added.append(s)
        if not added:
            return
        self.recipe["sites"].extend(added)
        self.on_structure_changed()
        self.statusBar().showMessage(
            f"added {len(added)} sideband line(s) of {bname} at ±νrot"
            + (" — position and shape linked to the parent, amplitudes free"
               if link.isChecked() else ""))

    def add_sidebands_for_line(self, idx: int):
        """The fit table's right-click ▸ Add spinning sidebands…: the same
        dialog as Decomposition ▸ Add spinning sidebands…, line ``idx``
        preselected."""
        self._ssb_preselect = int(idx)
        self.add_sidebands()

    # -------- spinning sidebands: detect the ±νrot repeat and offer it (F3)
    #: models that render their own manifold -- never offer sidebands on top
    #: of them (the set add_sidebands skips)
    _SSB_SELF_MODELS = frozenset({"sidebands", "csa_mas", "csa_czjzek",
                                  "quad_first"})

    def _ssb_key(self):
        """Identity of the trace on screen for the dismiss memory: (source,
        size, sum, rounded rate). Re-processing or a rate change makes a new
        key; a workspace switch back to the same trace does not."""
        if self.recipe is None or self.exp_amp is None or not self.exp_amp.size:
            return None
        return (self.source_path, int(self.exp_amp.size),
                float(np.sum(self.exp_amp)),
                round(float(self.recipe.get("spin_rate_Hz", 0.0) or 0.0)))

    def _model_has_manifold(self) -> bool:
        """The model already carries a manifold: a self-sidebanding model, a
        line whose position is an expression (linked sidebands) or a shifted
        copy of the spectrum -- i.e. the results of both offer actions and of
        Add spinning sidebands..., without any label heuristic."""
        for s in (self.recipe or {}).get("sites", []) or []:
            if s.get("model") in self._SSB_SELF_MODELS:
                return True
            p = s.get("params", {}) or {}
            if (p.get("isotropic_chemical_shift_ppm") or {}).get("expr"):
                return True
            shift = (p.get("shift_ppm") or {}).get("value", 0.0) or 0.0
            if s.get("model") == "spectrum" and float(shift) != 0.0:
                return True
        return False

    def _run_sideband_detection(self, *, scan: bool = False):
        """sidebands.detect on the trace on screen with the recipe's Larmor
        frequency and rate; None when there is nothing to run on."""
        from larmor import sidebands

        if self.recipe is None or self.exp_amp is None or self.exp_amp.size < 64:
            return None
        lar = float(self.recipe.get("larmor_frequency_MHz", 0.0) or 0.0)
        nu = float(self.recipe.get("spin_rate_Hz", 0.0) or 0.0)
        return sidebands.detect(self.exp_ppm, self.exp_amp, lar, nu, scan=scan)

    def _maybe_offer_sidebands(self):
        """Evaluate the offer for the active 1D document -- hooked into
        _update_sn, so load, workspace switch, re-processing, FID processing,
        subtraction and the experiment / calibrate dialogs all pass here.
        Every failed guard hides a previous offer. The detector runs around
        the recorded rate; it scans only when that rate is flagged uncertain
        (the red MAS pill is an open question the data can answer); a
        confirmed-static recipe never runs it."""
        act = getattr(self, "actSsbOffer", None)
        if (act is None or not act.isChecked() or not hasattr(self, "view")
                or self.central_stack.currentWidget() is not self.view
                or self.recipe is None or self.exp_amp is None
                or self.exp_amp.size < 64):
            self._dismiss_sideband_offer(remember=False)
            return
        rate = float(self.recipe.get("spin_rate_Hz", 0.0) or 0.0)
        uncertain = bool(self.recipe.get("mas_uncertain"))
        if rate <= 0 and not uncertain:          # confirmed static: nothing runs
            self._dismiss_sideband_offer(remember=False)
            return
        if self._model_has_manifold():
            self._dismiss_sideband_offer(remember=False)
            return
        key = self._ssb_key()
        if key is not None and key == getattr(self, "_ssb_dismissed", None):
            self._dismiss_sideband_offer(remember=False)
            return
        det = self._run_sideband_detection(scan=uncertain)
        if det is not None and det.ok:
            self._show_sideband_offer(det)
        else:
            self._dismiss_sideband_offer(remember=False)

    def detect_sidebands(self):
        """Decomposition > Detect spinning sidebands: the offer on demand --
        the ±2 % window around the recorded rate, then a 1–80 kHz scan -- or
        the reason nothing was found, in the status bar."""
        if self.recipe is None or self.exp_amp is None or self.exp_amp.size < 2:
            self.statusBar().showMessage("load a spectrum first")
            return
        if self.central_stack.currentWidget() is not self.view:
            self.statusBar().showMessage(
                "sideband detection works on the 1D workbench")
            return
        self._ssb_dismissed = None               # the user asked: forget ✕
        det = self._run_sideband_detection(scan=True)
        if det is not None and det.ok:
            self._show_sideband_offer(det)
            self.statusBar().showMessage(
                "spinning sidebands found — Return adds the linked manifold, "
                "Esc dismisses")
        else:
            self._dismiss_sideband_offer(remember=False)
            why = det.message if det is not None else "too few points"
            self.statusBar().showMessage(f"no ±νrot repeat found — {why}")

    def _show_sideband_offer(self, det):
        from larmor.desktop.sideband_offer import SidebandBanner

        banner = getattr(self, "ssb_banner", None)
        if banner is None:
            banner = SidebandBanner(self.view)
            banner.manifold.connect(
                lambda: self._apply_sideband_choice("linked"))
            banner.model_line.connect(
                lambda: self._apply_sideband_choice("model"))
            banner.copy.connect(lambda: self._apply_sideband_choice("copy"))
            banner.dismissed.connect(self._dismiss_sideband_offer)
            self.ssb_banner = banner
        self._ssb_detection = det
        rate = float(self.recipe.get("spin_rate_Hz", 0.0) or 0.0)
        banner.set_detection(det, rate)
        banner.show_guides(det)
        banner.show()
        banner.raise_()

    def _dismiss_sideband_offer(self, remember: bool = True):
        """Hide the banner and its guides. ``remember`` (the user's ✕ / Esc)
        keeps the trace's key so the same trace is not re-offered until the
        data or the rate changes; automatic dismissals pass False."""
        if remember and getattr(self, "_ssb_detection", None) is not None:
            self._ssb_dismissed = self._ssb_key()
        self._ssb_detection = None
        banner = getattr(self, "ssb_banner", None)
        if banner is not None:
            banner.dismiss()

    def _toggle_ssb_offer(self, on: bool):
        QSettings("LARMOR", "app").setValue("ssbAutoOffer", bool(on))
        self.statusBar().showMessage(
            "spinning-sideband offer: "
            + ("ON — a banner appears when a loaded spectrum repeats at ±νrot"
               if on else
               "off (Edit ▸ Spinning sidebands ▸ Detect spinning sidebands still "
               "works on demand)"))
        self._maybe_offer_sidebands()

    @staticmethod
    def _fresh_params(model_name: str) -> dict:
        """{name: {value, stderr, vary, min, max, expr}} from a model's
        ParamDefs -- the add_site_at idiom."""
        return {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                         "min": p.min, "max": p.max, "expr": None}
                for p in model_registry.get(model_name).params}

    def _ssb_set_rate_if_due(self, det, force: bool = False) -> bool:
        """Write the measured νrot into the recipe when the recorded rate is
        missing or flagged uncertain (the red pill asked the question), or
        when ``force`` (the `sidebands` model line has no spacing parameter
        of its own). A certain rate is otherwise never touched -- the ⚠ in
        the banner is the user's cue to double-click the experiment strip."""
        recorded = float(self.recipe.get("spin_rate_Hz", 0.0) or 0.0)
        if not (force or recorded <= 0 or self.recipe.get("mas_uncertain")):
            return False
        self.recipe["spin_rate_Hz"] = float(round(det.nu_rot_Hz))
        self.recipe["mas_uncertain"] = False
        self._update_exp_label()
        return True

    def _apply_sideband_choice(self, choice: str, det=None):
        """All three banner actions ('linked' / 'copy' / 'model'): ONE
        snapshot carrying the axis -- so Ctrl+Z reverts the manifold and a
        written rate together and refreshes the experiment strip -- then one
        recipe mutation and one structure refresh."""
        from larmor import sidebands

        det = det if det is not None else getattr(self, "_ssb_detection", None)
        if det is None or self.recipe is None:
            return
        self.snapshot(with_axis=True)
        recorded = float(self.recipe.get("spin_rate_Hz", 0.0) or 0.0)
        spacing = float(det.spacing_ppm)
        force = False
        if choice == "copy":
            n = self._add_sideband_copies(det)
            msg = (f"added {n} shifted cop{'y' if n == 1 else 'ies'} of this "
                   f"spectrum at ±{spacing:.4g} ppm (shift held) — Fit scales "
                   "them")
        elif choice == "model":
            n_ssb, ratio = self._add_sideband_model_line(det)
            msg = f"added a `sidebands` line (n_ssb {n_ssb}, r {ratio:.2f})"
            force = (recorded > 0 and abs(det.nu_rot_Hz - recorded)
                     > sidebands.RATE_MISMATCH * recorded)
        else:
            n, parent = self._add_sideband_manifold(det)
            msg = (f"added {n} linked sideband line(s) of {parent} at ±k·νrot "
                   f"({spacing:.4g} ppm) — Fit scales them")
        if self._ssb_set_rate_if_due(det, force=force):
            msg += (f" · νrot set to "
                    f"{sidebands.format_hz(self.recipe['spin_rate_Hz'])} Hz "
                    "(a new spin rate builds a new kernel once)")
        # the manifold itself now stops the offer (_model_has_manifold), so
        # the dismiss key is left to the user's own ✕ / Esc
        self._dismiss_sideband_offer(remember=False)
        self.on_structure_changed()
        self.statusBar().showMessage(msg)

    def _sideband_parent_index(self, det) -> int:
        """An existing eligible line sitting on the centreband (within half
        its FWHM), else a new gauss_lor centreband seeded from the picked
        peak (position, FWHM, height). Returns its index."""
        tol = max(0.5 * float(det.centre_fwhm_ppm), 1e-9)
        best = None
        for i, s in enumerate(self.recipe["sites"]):
            if s.get("model") in self._SSB_SELF_MODELS or s.get("model") == "spectrum":
                continue
            p = s.get("params", {}) or {}
            pos = p.get("isotropic_chemical_shift_ppm")
            if not pos or pos.get("expr") or "amplitude" not in p:
                continue
            d = abs(float(pos["value"]) - float(det.centre_ppm))
            if d <= tol and (best is None or d < best[0]):
                best = (d, i)
        if best is not None:
            return best[1]
        params = self._fresh_params("gauss_lor")
        params["isotropic_chemical_shift_ppm"]["value"] = float(det.centre_ppm)
        params["shift_fwhm_ppm"]["value"] = float(max(det.centre_fwhm_ppm, 0.1))
        params["amplitude"]["value"] = float(det.centre_height)
        n = len(self.recipe["sites"])
        label = model_registry.get("gauss_lor").label.split(" ")[0]
        self.recipe["sites"].append(
            {"model": "gauss_lor", "label": f"{label}-{n}", "params": params})
        return n

    def _add_sideband_manifold(self, det):
        """Linked lines, one copy of the parent per MATCHED order: position =
        parent ± k·νrot as a constraint (the LinkPositionDialog form, so the
        table shows 'A+124.6' / 'A-124.6' and the paddle is not draggable),
        every other parameter tied to the parent's, only the amplitude free
        and seeded from the measured tooth -- the linking Add spinning
        sidebands... applies, without its dialog. Sideband intensities follow
        the tensor (Herzfeld-Berger), not a ratio, so amplitudes are never
        tied. Returns (count, parent label)."""
        import copy

        p_idx = self._sideband_parent_index(det)
        parent = self.recipe["sites"][p_idx]
        base_name = parent.get("label") or parent["model"]
        p_pos = float(parent["params"]["isotropic_chemical_shift_ppm"]["value"])
        p_amp = float(parent["params"].get("amplitude", {}).get("value", 1.0) or 1.0)
        added = []
        for o in sorted(det.matched(), key=lambda o: (abs(o.k), -o.k)):
            off = float(o.k * det.spacing_ppm)
            s = copy.deepcopy(parent)
            s["label"] = f"{base_name}{o.k:+d}sb"
            for name, prm in s["params"].items():
                prm["stderr"] = None
                if name == "amplitude":
                    prm["expr"] = None
                    prm["vary"] = True
                    frac = (o.height / det.centre_height
                            if det.centre_height > 0 and o.height > 0
                            else 0.3 ** abs(o.k))
                    prm["value"] = p_amp * float(frac)
                elif name == "isotropic_chemical_shift_ppm":
                    prm["value"] = p_pos + off
                    prm["expr"] = (f"s{p_idx}.isotropic_chemical_shift_ppm "
                                   f"{'+' if off >= 0 else '-'} {abs(off):.6g}")
                else:
                    prm["expr"] = f"s{p_idx}.{name}"
            added.append(s)
        self.recipe["sites"].extend(added)
        return len(added), base_name

    def _add_sideband_copies(self, det) -> int:
        """One `spectrum` component per detected side at ±νrot, shift held,
        amplitude = the measured ±1 tooth height (the reference is unit-peak,
        so that value reproduces the tooth). One copy per SIDE, not per
        order: a copy at +νrot already carries the higher orders of the
        pattern. Returns the count."""
        n = 0
        for sign in det.sides():
            tooth = next(o for o in det.orders if o.k == sign and o.matched)
            self._append_spectrum_site(
                self.exp_ppm, self.exp_amp, f"copy{sign:+d}sb",
                amplitude=float(max(tooth.height, 0.0)),
                shift_ppm=float(sign * det.spacing_ppm), shift_vary=False)
            n += 1
        return n

    def _add_sideband_model_line(self, det):
        """One empirical `sidebands` site seeded from the detection: centre
        position / FWHM / height, r = the measured height ratio, n_ssb = the
        highest matched order. Returns (n_ssb, ratio)."""
        from larmor import sidebands

        params = self._fresh_params("sidebands")
        params["isotropic_chemical_shift_ppm"]["value"] = float(det.centre_ppm)
        params["shift_fwhm_ppm"]["value"] = float(max(det.centre_fwhm_ppm, 0.05))
        params["amplitude"]["value"] = float(det.centre_height)
        ratio = float(sidebands.seed_ratio(det))
        n_ssb = int(max(det.n_orders, 1))
        params["ssb_ratio"]["value"] = ratio
        params["n_ssb"]["value"] = float(n_ssb)
        params["gl"]["value"] = 1.0
        n = len(self.recipe["sites"])
        self.recipe["sites"].append(
            {"model": "sidebands", "label": f"ssb-{n}", "params": params})
        return n_ssb, ratio
