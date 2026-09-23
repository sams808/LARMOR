"""Main-window mixin: the workbench side of ``larmor.processing``.

``apply_processing`` and the display refresh (FID / spectrum toggle, real /
imag / magnitude channels), drag-to-phase, calibrate / measure / reset,
the manual, two-point and iterative baseline tools, zones, the experiment
parameters dialog, the processing-steps editor, WURST correction and
spectrum subtraction.

Owned state: ``_proc_base``, ``_proc_spec``, ``_proc_fid``,
``_proc_apply_count``, ``_phase_live`` / ``_phase_start`` /
``_phase_pivot_frac_last``.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QApplication, QMessageBox


class _ProcessingMixin:
    """Live processing, phasing, calibration, baseline tools."""

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

    def edit_experiment(self):
        if self.recipe is None:
            return
        from larmor.desktop.dialogs import ExperimentDialog

        self.snapshot(with_axis=True)         # an SR change re-references the axis
        old_sr = self.recipe.get("sr_hz", 0.0) or 0.0
        dlg = ExperimentDialog(self, self.recipe)
        if dlg.exec():
            self.recipe["mas_uncertain"] = False     # user confirmed the params
            # a changed SR re-references the ppm axis: delta = (nu - SF) / SF,
            # so a LARGER SR moves every peak to LOWER ppm (TopSpin's
            # direction; this used to add the shift instead)
            new_sr = self.recipe.get("sr_hz", 0.0) or 0.0
            larmor = self.recipe.get("larmor_frequency_MHz", 0.0) or 0.0
            if larmor and abs(new_sr - old_sr) > 1e-9 and self.exp_ppm.size:
                from larmor.referencing import axis_shift_ppm
                d_ppm = axis_shift_ppm(old_sr, new_sr, larmor)
                self.exp_ppm = self.exp_ppm + d_ppm
                if self._proc_base is not None:
                    self._proc_base = (self._proc_base[0] + d_ppm,
                                       self._proc_base[1])
                self.view.set_experiment(self.exp_ppm, self.exp_amp)
                # move the fitted sites with the re-referenced axis (see calibrate)
                self._shift_recipe_positions(self.recipe, d_ppm)
                if self.recipe.get("sites"):
                    self.lines_table.rebuild(self.recipe, self.hidden)
                    self.request_simulation()
            self._update_exp_label()
            self._maybe_offer_sidebands()     # νrot / SR may have changed
            self.statusBar().showMessage(
                "experiment updated — re-simulating (a new spin rate builds "
                "a new kernel once)")
            self.request_simulation()
            self._persist_session()
        else:
            self.undo_stack.pop()   # dialog cancelled: drop the snapshot

    # ------------------------------------------------------------- baseline
    def _baseline_mode(self, on: bool):
        if on:
            self.proc_panel.btnTpPick.setChecked(False)   # not the 2-point picker
            self.proc_panel.btnDrag.setChecked(False)     # nor drag-to-phase
        self.view.set_baseline_mode(on)
        self._set_add_mode(None)
        if on:
            self.statusBar().showMessage(
                "manual baseline (dmfit-style): click to add anchor points, "
                "drag to shape — it is subtracted automatically when you turn "
                "'Pick anchors' back off")
        elif len(self.view.baseline_anchors()) >= 2:
            # dmfit behaviour: exiting anchor mode auto-applies the correction
            self.apply_manual_baseline()

    def start_twopoint_bg(self):
        """Process-menu entry point: reveal the processing panel and arm its
        two-point picker so the user can click two baseline points. The panel's
        Pick / Subtract / Clear buttons then drive it (same flow as the panel)."""
        self.proc_dock.show()
        self.proc_dock.raise_()
        self.proc_panel.btnTpPick.setChecked(True)     # -> _twopoint_mode(True)

    def _twopoint_mode(self, on: bool):
        """Pick two baseline points; subtract the straight line through them."""
        self.proc_panel.btnBlPick.setChecked(False)   # not the anchor baseline
        if on:
            self.proc_panel.btnDrag.setChecked(False)  # nor drag-to-phase
        self.view.set_baseline_mode(on)
        self._set_add_mode(None)
        if on:
            self.statusBar().showMessage(
                "2-point background: click TWO baseline points (one each side of "
                "the peaks); the straight line through them is subtracted when you "
                "turn 'Pick 2 points' back off")
        elif len(self.view.baseline_anchors()) >= 2:
            self.apply_twopoint_bg()

    def apply_twopoint_bg(self):
        """Subtract the straight line through the first and last picked points — a
        flat/tilted 2-point background, extrapolated across the whole spectrum."""
        pts = self.view.baseline_anchors()
        if len(pts) < 2 or self.exp_ppm is None or not len(self.exp_ppm):
            self.statusBar().showMessage("click two baseline points first")
            return
        (x1, y1), (x2, y2) = pts[0], pts[-1]
        if abs(x2 - x1) < 1e-9:
            self.statusBar().showMessage("pick two points at different positions")
            return
        m = (y2 - y1) / (x2 - x1)                      # line through the two points
        base = m * (self.exp_ppm - x1) + y1
        self.snapshot(with_axis=True)                 # undoable (changes the data)
        self.exp_amp = self.exp_amp - base
        if self._proc_base is not None:               # keep live-processing in sync
            bx, by = self._proc_base
            self._proc_base = (bx, by - (m * (bx - x1) + y1))
        self.view.clear_baseline()
        self.proc_panel.btnTpPick.setChecked(False)
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.request_simulation()
        self.statusBar().showMessage(
            "2-point (linear) background subtracted — 'Reset to original' / Undo")

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

    def apply_iterbaseline(self):
        """Iterative dead-time baseline correction (Yon et al. 2020) — interactive
        dialog with a live preview, then apply as a recorded processing step."""
        if self.exp_ppm is None or self.exp_amp is None:
            self.statusBar().showMessage("open a spectrum first")
            return
        from larmor.desktop.baseline_dialog import BaselineDialog

        dlg = BaselineDialog(self, self.exp_ppm, self.exp_amp)
        if dlg.exec() != dlg.Accepted:
            return
        self.apply_processing(
            [{"op": "iterbaseline", **dlg.params()}], False)
        self.statusBar().showMessage(
            "iterative baseline (Yon et al. 2020) applied — "
            "'Reset to original' undoes")

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

    # ------------------------------------------------------------- processing
    def apply_processing(self, ops: list, use_raw: bool):
        if not self.source_path:
            return
        # only the 1D workbench is processed live; ignore while a 2D map is up
        if self.central_stack.currentWidget() is not self.view:
            return
        if not use_raw and not self.exp_ppm.size:
            return
        import dataclasses

        from larmor import processing as proc
        from larmor.io import bruker

        live = getattr(self.proc_panel, "chkLive", None)
        live = bool(live and live.isChecked())
        if not live:
            self.statusBar().showMessage("processing…")
            QApplication.processEvents()
        try:
            if use_raw:
                src = Path(self.source_path)
                if not bruker.is_expno(src):
                    # a dropped `fid` FILE leaves source_path at the file (the
                    # preview); resolve it to its EXPNO instead of refusing
                    try:
                        src = bruker.resolve(src).expno
                    except Exception:
                        raise ValueError(
                            "raw-fid processing needs a Bruker EXPNO")
                s = proc.from_bruker_fid(str(src))
            else:
                # apply the pipeline from the UNPROCESSED baseline every time, so
                # a live slider shows the absolute phase rather than compounding
                if self._proc_base is None:
                    self._proc_base = (self.exp_ppm.copy(), self.exp_amp.copy())
                base_ppm, base_amp = self._proc_base
                sfo1 = self.recipe.get("larmor_frequency_MHz", 0.0) if self.recipe else 0.0
                s = proc.from_processed(base_ppm, base_amp, sfo1)
            # TopSpin-style p1 pivot: rotate first-order phase about the pivot
            # line (default: the tallest peak), not the fixed spectrum centre
            piv = self.view.phase_pivot_frac()
            ops = [dict(o, pivot_frac=piv) if o.get("op") == "phase" else o
                   for o in ops]
            # split at the LAST ft: the state just before it is the windowed,
            # zero-filled FID the transform sees (the FID display). Copied --
            # ops mutate the Spectrum1D in place and op_ft reassigns y / x_ppm
            # / domain on the same object, so a reference would become the
            # spectrum.
            k = max((i for i, o in enumerate(ops) if o.get("op") == "ft"),
                    default=None)
            if k is None:
                s = proc.apply(s, ops)
                fid = None
            else:
                s = proc.apply(s, ops[:k])
                fid = dataclasses.replace(s, y=np.array(s.y, copy=True))
                s = proc.apply(s, ops[k:])
            if s.domain != "freq":
                raise ValueError("pipeline must end in the frequency domain")
        except Exception as exc:
            self._proc_spec = self._proc_fid = None   # never show a stale result
            if not live:                       # never nag on every keystroke
                QMessageBox.warning(self, "Processing failed", str(exc))
            self.statusBar().showMessage("processing failed")
            if self.proc_panel.view_state() != ("freq", "real"):
                self._display_fallback()       # the canvas may show an old FID
            return
        order = np.argsort(s.x_ppm)
        self.exp_ppm, self.exp_amp = np.asarray(s.x_ppm)[order], s.y.real[order]
        self._proc_spec = proc.Spectrum1D(
            x_ppm=np.asarray(s.x_ppm, float)[order],
            y=np.asarray(s.y, complex)[order],
            sfo1_MHz=s.sfo1_MHz, sw_Hz=s.sw_Hz)
        self._proc_fid = fid
        self._proc_apply_count += 1
        self._update_sn()
        # remember the pipeline in the recipe: saving the fit then saves the
        # processing that produced the spectrum it was fitted against
        if self.recipe is not None:
            self.recipe["processing"] = list(ops)
            self.recipe["processing_from_raw"] = bool(use_raw)
        self.request_simulation()
        self.statusBar().showMessage(
            f"processing applied ({len(ops)} step(s), stored in the recipe)")
        # draw whichever projection the panel selects (last, so a FID /
        # channel hint replaces the generic status line)
        self._refresh_display()

    # ------------------------------------------------------- drag to phase (F2)
    def start_phase_drag(self):
        """Process ▸ Drag to phase (Ctrl+P): reveal the Processing panel and
        toggle its Drag-to-phase button -- the single source of truth for the
        mode; its toggled signal runs _phase_drag_mode."""
        self.proc_dock.show()
        self.proc_dock.raise_()
        self.proc_panel.btnDrag.setChecked(not self.proc_panel.btnDrag.isChecked())

    def _phase_drag_mode(self, on: bool):
        """Arm / release the TopSpin gesture on the plot (horizontal drag =
        p0, vertical = p1 about the pivot line). Entering shows the dock (the
        pivot and the p0/p1 numbers are the visible feedback), releases the
        two baseline pickers and add mode, and ticks 'Hilbert first' in pdata
        mode: the workbench spectrum is real-only (exp_amp = y.real and
        from_processed zero-fills the imaginary part), so without the
        reconstructed imaginary channel a phase step would merely scale the
        spectrum by cos(p0). Hilbert is never unticked on exit."""
        pp = self.proc_panel
        if on:
            self.proc_dock.show()
            self.proc_dock.raise_()
            pp.btnBlPick.setChecked(False)
            pp.btnTpPick.setChecked(False)
            self._set_add_mode(None)
            hilb = ""
            if not pp.rb_raw.isChecked() and not pp.chkHilbert.isChecked():
                pp.arm_hilbert()            # silent: the first move applies it
                hilb = (" · Hilbert first switched on (the workbench spectrum "
                        "has no imaginary part)")
            self.view.set_phase_drag_mode(True)
            self._phase_pivot_frac_last = self.view.phase_pivot_frac()
            self._phase_live = False
            self.statusBar().showMessage(
                "drag to phase: ← → p0 · ↑ ↓ p1 about the pivot · Shift = fine "
                "· Ctrl+drag = pan · Esc to stop" + hilb)
            return
        self.view.set_phase_drag_mode(False)
        # the pivot survives while the panel is open, goes with a closed panel
        self.view.show_phase_pivot(not self.proc_dock.isHidden())
        self._phase_live = False
        p0, p1 = pp.phase_values()
        self.statusBar().showMessage(
            f"phase kept: p0 {p0:+.2f}°  p1 {p1:+.1f}° — stored in the recipe "
            "as a phase step; Undo reverts each drag")

    def on_phase_dragged(self, dp0: float, dp1: float):
        """Cumulative gesture deltas from the view. One undo snapshot on the
        FIRST move of a gesture (the paddle idiom), then the panel's p0/p1
        follow the pointer and -- with 'live' ticked -- re-apply through the
        same op-list builder the sliders use, so Hilbert / magnitude / SR and
        the recipe record are identical to typing the numbers."""
        if (self.recipe is None or not self.exp_ppm.size
                or self.central_stack.currentWidget() is not self.view
                or self.view.domain != "freq"):
            return
        pp = self.proc_panel
        if not self._phase_live:
            self.snapshot(with_axis=True)
            self._phase_start = pp.phase_values()
            self._phase_live = True
        live = pp.chkLive.isChecked()
        pp.set_phase(self._phase_start[0] + dp0, self._phase_start[1] + dp1,
                     apply=live)
        if live:
            self._phase_pivot_frac_last = self.view.phase_pivot_frac()
        p0, p1 = pp.phase_values()
        # after set_phase, so it replaces apply_processing's generic line
        self.statusBar().showMessage(
            f"phase: p0 {p0:+.2f}°   p1 {p1:+.1f}°   "
            "(Shift = fine · Ctrl+drag = pan · Esc to stop)")

    def on_phase_drag_released(self):
        pp = self.proc_panel
        if self._phase_live and not pp.chkLive.isChecked():
            pp.set_phase(*pp.phase_values())     # the one deferred apply
            self._phase_pivot_frac_last = self.view.phase_pivot_frac()
        self._phase_live = False
        self._persist_session()

    def _on_pivot_moved(self, new_frac: float):
        """The pivot line was dragged. While phasing, p0 is re-expressed so
        the spectrum on screen does not change (op_phase applies
        p0 + p1*(idx - pivot_frac)) and the recipe records the new pivot with
        the re-expressed p0; no undo entry -- dragging the pivot back is the
        undo. With the mode off the pre-existing behaviour stays (the pivot
        governs the next apply)."""
        if not self.view.phase_drag_active():
            self._phase_pivot_frac_last = new_frac
            return
        from larmor.phasedrag import compensate_p0

        pp = self.proc_panel
        p0, p1 = pp.phase_values()
        old, self._phase_pivot_frac_last = self._phase_pivot_frac_last, new_frac
        if p1 == 0.0:
            return                              # the pivot is irrelevant
        pp.set_phase(compensate_p0(p0, p1, old, new_frac), p1,
                     apply=pp.chkLive.isChecked())

    # ------------------------------------------- display projections (F6)
    def _proc_spec_valid(self) -> bool:
        """Is the complex pipeline result the spectrum on the workbench? Exact
        (same code path produced both; ~50 us on 64k points), so a bypassing
        edit of the exp arrays (2-point / manual baseline, WURST, subtract,
        calibrate, SR, undo, workspace switch, any load) can never leave a
        stale imaginary channel or FID on screen."""
        sp = self._proc_spec
        return (sp is not None and sp.y.size == self.exp_amp.size
                and np.array_equal(sp.x_ppm, self.exp_ppm)
                and np.array_equal(sp.y.real, self.exp_amp))

    def _mirror_view_actions(self):
        """Menu + sidebar follow the panel (setChecked: no triggered)."""
        domain, channel = self.proc_panel.view_state()
        self.actTimeDomain.setChecked(domain == "time")
        self.sbFid.setChecked(domain == "time")
        self.actChannel[channel].setChecked(True)

    def _display_fallback(self, hint: str | None = None):
        """Controls and canvas back to the real spectrum."""
        self.proc_panel.reset_view()
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self._mirror_view_actions()
        if hint:
            self.statusBar().showMessage(hint)

    def _on_view_experiment_set(self):
        """A real spectrum was placed on the canvas. With the FID or another
        channel selected: drop the selection when the complex result no
        longer matches the data (new spectrum, bypassing edit); re-draw the
        channel when the same spectrum was merely re-placed (workspace
        switch back, reload of an identical file)."""
        pp = getattr(self, "proc_panel", None)
        if pp is None:
            return
        domain, channel = pp.view_state()
        if (domain, channel) == ("freq", "real"):
            return
        if domain == "time" or not self._proc_spec_valid():
            pp.reset_view()
            if not self._proc_spec_valid():
                self._proc_spec = self._proc_fid = None
            self._mirror_view_actions()
            return
        self._refresh_display()

    def _refresh_display(self):
        """Draw the projection the panel selects of the last pipeline result:
        (freq, real) the real spectrum as always; (freq, imag | magnitude)
        that channel of the complex result on the same ppm axis; (time, ch)
        the windowed FID captured before the last ft. What the result cannot
        honour falls back to the real spectrum with a hint."""
        from larmor import processing as proc

        pp = self.proc_panel
        domain, channel = pp.view_state()
        if (domain, channel) == ("freq", "real"):
            self.view.set_experiment(self.exp_ppm, self.exp_amp)
            self._mirror_view_actions()
            return
        if not self._proc_spec_valid():
            self._display_fallback(
                "display reset to the real spectrum — Apply processing "
                "rebuilds the complex channel")
            return
        if domain == "time":
            fid = self._proc_fid
            if fid is None:
                self._display_fallback(
                    "the time domain needs a transform in the pipeline — pick "
                    "raw fid, or tick re-apodize (needs the Larmor frequency)")
                return
            t_ms = proc.time_axis_s(fid) * 1e3
            label = {"real": "FID (real)", "imag": "FID (imag)",
                     "magnitude": "|FID|"}[channel]
            self.view.set_fid(t_ms, proc.channel_view(fid.y, channel), label)
            self.statusBar().showMessage(
                f"time domain — {fid.y.size} points, AQ {t_ms[-1]:.2f} ms · "
                "LB / WDW / ZF re-apply live · FID ⇄ spectrum (Ctrl+T) returns")
        else:
            self.view.set_channel_trace(
                proc.channel_view(self._proc_spec.y, channel),
                "experiment (imag)" if channel == "imag" else "|experiment|")
            msg = (("showing the imaginary channel" if channel == "imag"
                    else "showing |S| (display only)")
                   + " — the fit always uses the real spectrum")
            if not np.any(self._proc_spec.y.imag):
                msg += " · this spectrum has no imaginary channel: tick Hilbert first"
            self.statusBar().showMessage(msg)
        self._mirror_view_actions()

    def _on_proc_view_changed(self, domain: str, channel: str):
        """The panel (or its menu / sidebar mirrors) asked for a display
        projection. A non-default one needs a valid complex result: when
        there is none (fresh data, a bypassing edit) the panel is first synced
        to the RECORDED chain so the forced Apply replays that chain -- not
        the panel's leftovers -- then the time domain arms re-apodization for
        a pdata source (no ft in the chain) and the imaginary channel arms
        Hilbert (a real-only spectrum has none)."""
        pp = self.proc_panel
        want_default = (domain, channel) == ("freq", "real")
        if (self.central_stack.currentWidget() is not self.view
                or not self.exp_ppm.size or self.recipe is None):
            pp.reset_view()
            self._mirror_view_actions()
            if not want_default:
                self.statusBar().showMessage("load a 1D spectrum first")
            return
        if want_default:
            self._refresh_display()
            return
        need_apply = False
        if not self._proc_spec_valid():
            ok = pp.sync_from_ops(self.recipe.get("processing") or [],
                                  bool(self.recipe.get("processing_from_raw")))
            if not ok:
                self._display_fallback(
                    "the recorded pipeline has steps the panel cannot drive "
                    "(lp / shift_fid / whole-echo) — Process ▸ Processing steps")
                return
            need_apply = True
        if domain == "time" and not pp.chain_has_ft():
            pp.arm_reapodize()
            need_apply = True
            self.statusBar().showMessage(
                "re-apodizing a processed spectrum: Hilbert-reconstructed "
                "imaginary channel, the window compounds with the one already "
                "applied (raw-fid mode / Open FID are exact)")
        elif (channel != "real" and not pp.chain_has_ft()
                and not pp.chkHilbert.isChecked()):
            pp.arm_hilbert()
            need_apply = True
        if need_apply or (domain == "time" and self._proc_fid is None):
            n0 = self._proc_apply_count
            pp.btnApply.click()       # -> _emit([]) -> apply_processing -> _refresh_display
            if (self._proc_apply_count == n0
                    and pp.view_state() != ("freq", "real")):
                # the pipeline did not run (no source file, a 2D map up)
                self._display_fallback(
                    "processing could not run — the display stays on the "
                    "real spectrum")
        else:
            self._refresh_display()

    def _toggle_time_domain(self, on: bool):
        """Process ▸ FID ⇄ spectrum (Ctrl+T) and the sidebar FID button:
        reveal the Processing dock and drive its button, which emits the
        view_changed the workbench listens to (same flow as the panel)."""
        pp = self.proc_panel
        if on:
            self.proc_dock.show()
            self.proc_dock.raise_()
        if pp.btnDomain.isChecked() != bool(on):
            pp.btnDomain.setChecked(bool(on))        # -> view_changed
        self._mirror_view_actions()          # the trigger never stays out of sync

    def _set_channel(self, name: str):
        pp = self.proc_panel
        rb = {"real": pp.rb_real, "imag": pp.rb_imag, "magnitude": pp.rb_mag}[name]
        if not rb.isChecked():
            rb.setChecked(True)                      # -> view_changed
        self._mirror_view_actions()

    def _cycle_channel(self):
        from larmor.processing import CHANNELS

        cur = self.proc_panel.view_state()[1]
        self._set_channel(CHANNELS[(CHANNELS.index(cur) + 1) % len(CHANNELS)])

    # ------------------------------------------------------------- calibrate
    def start_calibrate(self):
        if not self.exp_ppm.size:
            self.statusBar().showMessage("load a spectrum first")
            return
        self.view.set_calibrate_mode(True)
        self.statusBar().showMessage(
            "calibrate: click the peak whose shift you want to set")

    @staticmethod
    def _shift_recipe_positions(recipe, delta: float):
        """Shift every site's absolute ppm position by delta. Used when the axis
        is re-referenced (calibrate / SR) so the fitted model tracks the peaks
        instead of being left behind. Widths/η/σ/amplitude are unaffected."""
        for s in recipe.get("sites", []):
            params = s.get("params", {})
            for pname in ("isotropic_chemical_shift_ppm", "shift_ppm"):
                p = params.get(pname)
                if isinstance(p, dict) and "value" in p:
                    p["value"] = float(p["value"]) + delta

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
        self.snapshot(with_axis=True)         # calibration is now undoable
        self.exp_ppm = self.exp_ppm + delta
        if self._proc_base is not None:
            self._proc_base = (self._proc_base[0] + delta, self._proc_base[1])
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        larmor = self.recipe.get("larmor_frequency_MHz", 0.0) if self.recipe else 0.0
        if self.recipe is not None:
            # re-referencing relabels the axis; move the fitted sites the same way
            # so the model stays on the peaks (was: only exp_ppm shifted -> desync)
            self.recipe["calibration_ppm"] = \
                self.recipe.get("calibration_ppm", 0.0) + delta
            self._shift_recipe_positions(self.recipe, delta)
            if self.recipe.get("sites"):
                self.lines_table.rebuild(self.recipe, self.hidden)
        self.request_simulation()
        hz = delta * larmor
        self.statusBar().showMessage(
            f"axis shifted by {delta:+.3f} ppm"
            + (f"  (SR {hz:+.1f} Hz)" if larmor else ""))
        self._maybe_offer_sidebands()         # the axis moved under the guides

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
            dhz = dppm * larmor
            msg += (f"   {dhz / 1000.0:.3f} kHz" if dhz >= 1000.0
                    else f"   {dhz:.1f} Hz")
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

    def open_wurst_correct(self):
        """Process > WURST excitation profile: divide the current spectrum by
        the computed WURST-N sweep weighting so intensities across a swept
        wideline pattern are comparable. The amplitude half of swept-pulse
        physics; autophase's p2 is the phase half. Recorded in the recipe's
        provenance."""
        from PySide6.QtWidgets import (QDialog, QDialogButtonBox,
                                       QDoubleSpinBox, QFormLayout)

        from larmor.processing import wurst_profile

        if not self.exp_ppm.size:
            self.statusBar().showMessage("load a spectrum first")
            return
        sfo = float((self.recipe or {}).get("larmor_frequency_MHz", 0.0) or 0.0)
        if sfo <= 0:
            QMessageBox.warning(self, "WURST profile",
                                "the Larmor frequency is needed — set it in "
                                "Experiment parameters first")
            return
        dlg = QDialog(self)
        dlg.setWindowTitle("WURST excitation profile")
        form = QFormLayout(dlg)
        span_khz = float(np.ptp(self.exp_ppm)) * sfo / 1e3
        centre = QDoubleSpinBox(); centre.setRange(-1e6, 1e6)
        centre.setDecimals(1); centre.setSuffix(" ppm")
        centre.setValue(float(0.5 * (self.exp_ppm.max() + self.exp_ppm.min())))
        form.addRow("sweep centre (carrier)", centre)
        sweep = QDoubleSpinBox(); sweep.setRange(1.0, 100000.0)
        sweep.setDecimals(1); sweep.setSuffix(" kHz")
        sweep.setValue(round(span_khz, 1))
        sweep.setToolTip("the WURST pulse's total sweep width (e.g. 2000 kHz "
                         "for a 2 MHz WCPMG sweep)")
        form.addRow("sweep width", sweep)
        order = QDoubleSpinBox(); order.setRange(2.0, 200.0)
        order.setDecimals(0); order.setValue(80.0)
        order.setToolTip("the N in WURST-N (envelope 1 − |cos(πt/τ)|^N)")
        form.addRow("WURST order N", order)
        floor = QDoubleSpinBox(); floor.setRange(1.0, 100.0)
        floor.setDecimals(0); floor.setSuffix(" %"); floor.setValue(10.0)
        floor.setToolTip("never divide by less than this fraction of full "
                         "excitation — dividing further amplifies edge noise")
        form.addRow("correction floor", floor)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject)
        form.addRow(bb)
        if dlg.exec() != QDialog.Accepted:
            return
        self.snapshot(with_axis=True)       # undoable (changes the data)
        w = wurst_profile(self.exp_ppm, sfo, float(centre.value()),
                          float(sweep.value()), n=float(order.value()),
                          floor=float(floor.value()) / 100.0)
        self.exp_amp = self.exp_amp / w
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        if self.recipe is not None:
            prov = dict(self.recipe.get("provenance") or {})
            prov["wurst_correct"] = {
                "centre_ppm": float(centre.value()),
                "sweep_kHz": float(sweep.value()),
                "n": float(order.value()),
                "floor": float(floor.value()) / 100.0}
            self.recipe["provenance"] = prov
        if self.recipe and self.recipe.get("sites"):
            self.request_simulation()
        self._update_sn(); self._refresh_overlays()
        self.statusBar().showMessage(
            f"WURST profile divided out ({sweep.value():.0f} kHz sweep, "
            f"N={order.value():.0f}) — File ▸ Save spectrum as… to keep it")

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
