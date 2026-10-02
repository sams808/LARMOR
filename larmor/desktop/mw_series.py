"""Main-window mixin: the series mode of the workbench.

Series ▸ Sequential fit… used to open a dialog of its own (a cramped lines
table, Prev / Next, a synchronous Fit current that froze the window). It now
turns the selected spectra into ordinary workspaces -- one per member, in the
given order, each tagged ``ws["series"]`` -- so every tool of the main window
(add / remove lines, bounds and links, paddles, the processing panel, Auto
fit, undo, the fit-health strip) applies to every member, and a thin series
bar above the plot (``larmor.desktop.series_bar.SeriesBar``) walks the
series, carries a model from one neighbour to the next
(``larmor.seriesmode.carry_into``), chains Fit → next through the window's
own FitWorker and runs the forward / backward auto sweep in a ``SeqWorker``
(``larmor.seqfit.run_sequential``, the function the CLI runs too). The
Series table, the Series plot, the Acquisition table, Save all fits and the
publication bundle live on the bar.

Owned state: ``_series`` (the active ``larmor.seriesmode.SeriesSpec`` or
None), ``series_bar`` (built once by ``_build_series_bar``, inserted above
the central stack), ``_seq_worker`` (the ``SeqWorker`` of a running auto
sweep), ``_series_chain`` (the FitWorker a Fit → next is waiting on),
``_series_sweep_passes`` (for the progress text). The member tags are a
projection of ``_series`` (``_series_retag``) and travel with the project
bundle; ``_series_restore_from_tags`` rebuilds the series after open_project.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QFileDialog, QInputDialog, QMessageBox

from larmor import fithealth
from larmor import seriesmode as sm
from larmor.desktop.windowtray import show_tool_window
from larmor.desktop.workers import humanize_error
from larmor.recipe import Recipe


class _SeriesMixin:
    """The series bar, the walk, the carry rules, Fit → next and the sweep."""

    # ------------------------------------------------------------- build
    def _build_series_bar(self):
        from larmor.desktop.series_bar import SeriesBar

        self._series = None
        self._seq_worker = None
        self._series_chain = None
        self._series_sweep_passes = 0
        self.series_bar = SeriesBar()
        # above the spectrum, inside the central page (not a dock): the
        # container's layout holds [series bar, central stack, health strip]
        self.centralWidget().layout().insertWidget(0, self.series_bar)
        bar = self.series_bar
        bar.member_clicked.connect(self.series_go)
        bar.prev_requested.connect(self.series_prev)
        bar.next_requested.connect(self.series_next)
        bar.fit_next_requested.connect(self.series_fit_then_next)
        bar.sweep_requested.connect(self.series_auto_sweep)
        bar.stop_requested.connect(lambda: self._interrupt_fit("stop"))
        bar.seed_on_move_toggled.connect(self.series_set_seed_on_move)
        bar.carry_changed.connect(self.series_set_carry)
        bar.carry_menu_opening.connect(self._series_fill_carry_menu)
        bar.copy_model_requested.connect(self.series_copy_model)
        bar.table_requested.connect(self.series_table)
        bar.plot_requested.connect(self.series_plot)
        bar.acquisition_requested.connect(self.series_acquisition_table)
        bar.save_all_requested.connect(self.series_save_all)
        bar.bundle_requested.connect(self.series_bundle)
        bar.details_requested.connect(self.series_show_comparability)
        bar.end_requested.connect(self.end_series)
        self.central_stack.currentChanged.connect(self._series_page_changed)
        bar.setVisible(False)

    # ------------------------------------------------------------- lookups
    def _series_rows(self) -> list:
        """``[(workspace index, tag)]`` of the active series' members in
        series order. A member whose workspace was closed drops out of the
        series here (the dock is the truth about what is open)."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return []
        by_key: dict = {}
        for i, ws in enumerate(self.workspaces):
            tag = ws.get("series")
            if (ws.get("kind") == "1d" and isinstance(tag, dict)
                    and tag.get("id") == spec.id):
                by_key.setdefault(tag.get("key"), i)
        gone = [m.key for m in spec.members if m.key not in by_key]
        for key in gone:
            spec.remove(key)
        if gone:
            self._series_retag()
        return [(by_key[m.key], self.workspaces[by_key[m.key]]["series"])
                for m in spec.members]

    def _series_retag(self):
        """Write the member tags from the spec (index, name, group, options,
        fit state, the Series table on the first member)."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        by_key = {ws["series"].get("key"): ws for ws in self.workspaces
                  if isinstance(ws.get("series"), dict)
                  and ws["series"].get("id") == spec.id}
        for k, m in enumerate(spec.members):
            ws = by_key.get(m.key)
            if ws is not None:
                ws["series"] = spec.tag_for(k)

    def _series_current_index(self):
        """The series position of the active workspace, None when the active
        document is not a member."""
        spec = getattr(self, "_series", None)
        if spec is None or self.active_ws is None:
            return None
        if not (0 <= self.active_ws < len(self.workspaces)):
            return None
        tag = self.workspaces[self.active_ws].get("series")
        if not isinstance(tag, dict) or tag.get("id") != spec.id:
            return None
        return spec.index_of(tag.get("key"))

    def _series_recipe_at(self, ws_i: int):
        if ws_i == self.active_ws:
            return self.recipe
        return (self.workspaces[ws_i].get("snap") or {}).get("recipe")

    def _series_member_data(self) -> list:
        """``[{key, name, recipe, ppm, amp, source_path, ws}]`` in series
        order from the workspace snapshots (the active one synced first)."""
        spec = self._series
        self._sync_active()
        out = []
        for k, (ws_i, _tag) in enumerate(self._series_rows()):
            snap = self.workspaces[ws_i].get("snap") or {}
            m = spec.members[k]
            out.append({"key": m.key, "name": m.name,
                        "recipe": snap.get("recipe") or {},
                        "ppm": snap.get("exp_ppm"), "amp": snap.get("exp_amp"),
                        "source_path": str(snap.get("source_path")
                                           or m.source_path or ""),
                        "ws": ws_i})
        return out

    def _series_open_1d(self) -> list:
        """Indices of the open 1D workspaces holding a spectrum, dock order."""
        self._sync_active()
        out = []
        for i, ws in enumerate(self.workspaces):
            snap = ws.get("snap") or {}
            ppm = snap.get("exp_ppm")
            if ws.get("kind") == "1d" and ppm is not None and len(ppm):
                out.append(i)
        return out

    def _series_busy(self) -> bool:
        w = getattr(self, "_seq_worker", None)
        if w is not None and w.isRunning():
            self.statusBar().showMessage(
                "the auto sweep is running — Stop (keep) or Cancel (revert) first")
            return True
        fw = getattr(self, "_fit_worker", None)
        if fw is not None and fw.isRunning():
            self.statusBar().showMessage(
                "a fit is still running — wait for it, or Stop / Cancel it")
            return True
        return False

    # ------------------------------------------------------------- the bar
    def _series_refresh_bar(self):
        """Redraw the bar from the spec and the workspaces: one button per
        member with its status dot (``seriesmode.member_status`` on the live
        recipe of the active member, the snapshot of the others; 'edited'
        when the recipe's signature left the one its last fit stamped),
        the current member, the comparability chip, the page visibility.
        Called from ``_refresh_ws_panel`` and the fit-health hooks; never
        re-simulates anything."""
        bar = getattr(self, "series_bar", None)
        spec = getattr(self, "_series", None)
        if bar is None:
            return
        if spec is None:
            bar.setVisible(False)
            return
        rows = self._series_rows()
        if not rows:
            self._series = None
            bar.setVisible(False)
            self._update_enabled()
            return
        members = []
        renamed = False
        for k, (ws_i, _tag) in enumerate(rows):
            m = spec.members[k]
            # a Rename… in the Workspaces dock reaches the member (the Series
            # table sets the same custom title the other way round)
            custom = self.workspaces[ws_i].get("custom_title")
            if custom and str(custom) != m.name:
                m.name = str(custom)
                renamed = True
            rec = self._series_recipe_at(ws_i)
            stale = bool(m.fit_sig) and sm.signature_key(rec) != m.fit_sig
            st = sm.member_status(rec, health_stale=stale, failed=m.failed)
            r = sm.rmsd_of(rec)
            flag = spec.comparability_flag(k)
            tip = m.name + (f" · RMSD {r:.4g}" if r is not None else " · not fitted")
            if m.source_path:
                tip += f" · {m.source_path}"
            if flag:
                tip += f"\n⚠ {flag} — unlike the series majority"
            members.append({"name": m.name, "status": st, "tip": tip, "flag": bool(flag)})
        if renamed:
            spec.refresh_comparison()
            self._series_retag()
        bar.set_members(members)
        cur = self._series_current_index()
        bar.set_current(-1 if cur is None else cur)
        bar.set_comparison(spec.comparison)
        bar.setVisible(self.central_stack.currentWidget() is self.view)

    def _series_page_changed(self, *_):
        self._series_refresh_bar()

    def _series_on_health(self, h):
        """Hook from ``_health_apply`` / ``_health_reset``: a fresh (not
        stale) fit verdict on a member stamps the recipe signature its dot
        refers to; any verdict change redraws the dots."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        k = self._series_current_index()
        if k is not None and h is not None and getattr(h, "fitted", False) \
                and not getattr(h, "stale", False):
            sig = sm.signature_key(self.recipe)
            if sig:
                m = spec.members[k]
                m.fit_sig, m.failed = sig, False
                self._series_retag()
        self._series_refresh_bar()

    def _series_on_fit_failed(self):
        """Hook from ``_fit_failed``: the member's dot turns red and a
        Fit → next chain stops."""
        self._series_chain = None
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        k = self._series_current_index()
        if k is not None:
            spec.members[k].failed = True
            self._series_retag()
        self._series_refresh_bar()

    # ------------------------------------------------------------- start / end
    def run_seq_fit(self, paths=None):
        """Series ▸ Sequential fit… (also the Explorer selection and the
        Session inventory's picks): the spectra become a series in the main
        window. With fewer than two paths but two or more 1D spectra open,
        the open spectra become the series instead."""
        if isinstance(paths, (list, tuple)):
            paths = [p for p in paths if p]
        else:
            paths = [p for p in self.explorer.selected_spectra() if p]
        if len(paths) >= 2:
            model = self.recipe if (self.recipe and self.recipe.get("sites")) else None
            self.start_series(paths, model)
            return
        if len(self._series_open_1d()) >= 2:
            self.start_series_from_workspaces()
            return
        self.statusBar().showMessage(
            "a series needs two or more spectra — Ctrl/Shift-select them in the "
            "Explorer, pick them in the Session inventory, or open them first")

    def start_series(self, paths, model=None):
        """Open ``paths`` as one workspace each (in this order; a path already
        open as a 1D workspace is adopted rather than reloaded), tag them as
        a series, apply ``model`` (the recipe the caller had on screen) to the
        FIRST member only, switch to it and show the bar."""
        paths = [str(p) for p in (paths or []) if p]
        if not paths:
            self.statusBar().showMessage("no spectra to make a series of")
            return
        if self._series_busy():
            return
        if self._series is not None:
            self.end_series(quiet=True)
        model_amp = sm.amp_max(self.exp_amp) if model else 0.0
        indices: list = []
        failed = 0
        for p in paths:
            i = self._series_find_open(p)
            if i is None:
                n0 = len(self.workspaces)
                self._ws_mode = "new"
                try:
                    self.load_source(p, keep_fit=False)
                finally:
                    self._ws_mode = "auto"
                i = self.active_ws
                if (len(self.workspaces) == n0 or i is None
                        or self.workspaces[i].get("kind") != "1d"):
                    failed += 1
                    continue
            if i not in indices:
                indices.append(i)
        if not indices:
            self.statusBar().showMessage("none of the spectra could be opened as a 1D "
                                         "spectrum — no series started")
            return
        self._series_begin(indices, model=model, model_amp=model_amp)
        if failed:
            self.statusBar().showMessage(
                self.statusBar().currentMessage()
                + f" · {failed} path(s) could not be opened and were left out")

    def start_series_from_workspaces(self, indices=None):
        """Series ▸ Series from open spectra: the open 1D workspaces (dock
        order, or the given indices) become a series."""
        if self._series_busy():
            return
        if isinstance(indices, (list, tuple)):
            rows = [int(i) for i in indices if isinstance(i, int)]
            rows = [i for i in rows if i in set(self._series_open_1d())]
        else:
            rows = self._series_open_1d()
        if len(rows) < 2:
            self.statusBar().showMessage(
                "open at least two 1D spectra first (File ▸ Open…, or Ctrl/Shift-"
                "select them in the Explorer and Series ▸ Sequential fit…)")
            return
        self._series_begin(rows)

    def _series_find_open(self, path: str):
        """An open 1D workspace showing ``path`` that is not in a series."""
        for i, ws in enumerate(self.workspaces):
            if ws.get("kind") != "1d" or isinstance(ws.get("series"), dict):
                continue
            src = (self.source_path if i == self.active_ws
                   else (ws.get("snap") or {}).get("source_path"))
            if src and str(src) == str(path):
                return i
        return None

    def _series_begin(self, indices: list, model=None, model_amp: float = 0.0):
        """The common tail of start_series / start_series_from_workspaces:
        the spec from the workspaces (names from scan.sample_label, a
        Rename… of the dock kept), the tags and titles, the caller's model
        into the first member, the switch to it, the bar."""
        if self._series is not None:
            self.end_series(quiet=True)
        self._sync_active()
        wss = [self.workspaces[i] for i in indices]
        snaps = [ws.get("snap") or {} for ws in wss]
        paths = [str(s.get("source_path") or "") for s in snaps]
        recs = [s.get("recipe") or {} for s in snaps]
        names = [ws.get("custom_title") for ws in wss]
        spec = sm.SeriesSpec.build(paths, recs, names=names)
        for k, rec in enumerate(recs):
            # a member opened from a saved fit starts 'fitted' and is tracked
            # for edits from its loaded values
            if rec.get("sites") and sm.rmsd_of(rec) is not None:
                spec.members[k].fit_sig = sm.signature_key(rec)
        self._series = spec
        for k, ws in enumerate(wss):
            ws["series"] = spec.tag_for(k)
            ws["custom_title"] = spec.members[k].name
            ws["title"] = spec.members[k].name
        note = ""
        if isinstance(model, dict) and model.get("sites"):
            snap = snaps[0]
            new, note = sm.carry_into(snap.get("recipe") or {}, model, None,
                                      sm.amp_max(snap.get("exp_amp")), model_amp,
                                      src_name="the model on screen")
            snap["recipe"] = new
            wss[0]["has_fit"] = bool(new.get("sites"))
            if indices[0] == self.active_ws:
                self.recipe = new
                self.on_structure_changed()
        self.series_go(0, seed=False)
        self._series_retag()
        self._refresh_ws_panel()
        self._update_enabled()
        self._series_fill_carry_menu()
        self._series_refresh_bar()
        n = spec.n
        self.statusBar().showMessage(
            f"series of {n} spectra — click a name or ◀ ▶ to walk it (the next "
            "spectrum takes this one's model), Fit → next fits and moves on, "
            "Auto sweep fits them all" + (f" · {note}" if note else ""), 15000)

    def end_series(self, quiet=False):
        """Series ▸ End series and the bar's ✕: the bar goes, the tags go,
        the spectra stay open."""
        w = getattr(self, "_seq_worker", None)
        if w is not None and w.isRunning():
            self.statusBar().showMessage("the auto sweep is running — Stop it first")
            return
        spec = getattr(self, "_series", None)
        self._series = None
        self._series_chain = None
        if spec is not None:
            for ws in self.workspaces:
                tag = ws.get("series")
                if isinstance(tag, dict) and tag.get("id") == spec.id:
                    ws.pop("series", None)
        bar = getattr(self, "series_bar", None)
        if bar is not None:
            bar.setVisible(False)
        self._refresh_ws_panel()
        self._update_enabled()
        if not quiet:
            self.statusBar().showMessage(
                "series ended — the spectra stay open in the Workspaces dock")

    # ------------------------------------------------------------- walking
    def series_go(self, k, seed=None):
        """Switch to member ``k``; with seed-on-move (or ``seed=True``) the
        spectrum landed on takes the model of the member just left
        (``seriesmode.carry_into``: a copy with scaled amplitudes into an
        empty one, carried values into one with its own lines), behind one
        undo snapshot."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        if self._series_busy():
            return
        rows = self._series_rows()
        try:
            k = int(k)
        except (TypeError, ValueError):
            return
        if not (0 <= k < len(rows)):
            return
        src_k = self._series_current_index()
        src_i = self.active_ws if src_k is not None else None
        dst_i, _tag = rows[k]
        if dst_i != self.active_ws:
            self.switch_workspace(dst_i)
        note = ""
        want = spec.options.seed_on_move if seed is None else bool(seed)
        if want and src_i is not None and src_i != dst_i and src_k is not None:
            note = self._series_seed_from_ws(src_i, spec.members[src_k].name)
        self._series_refresh_bar()
        self.statusBar().showMessage(
            f"series {k + 1}/{len(rows)}: {spec.members[k].name}"
            + (f" — {note}" if note else ""))

    def series_prev(self):
        k = self._series_current_index()
        if k is None:
            self.series_go(0, seed=False)
        elif k > 0:
            self.series_go(k - 1)

    def series_next(self):
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        k = self._series_current_index()
        if k is None:
            self.series_go(0, seed=False)
        elif k < spec.n - 1:
            self.series_go(k + 1)
        else:
            self.statusBar().showMessage("last spectrum of the series")

    def _series_seed_from_ws(self, src_i: int, src_name: str = "") -> str:
        """Carry workspace ``src_i``'s model into the active member."""
        spec = self._series
        snap = self.workspaces[src_i].get("snap") or {}
        src_rec = snap.get("recipe") or {}
        if not src_rec.get("sites") or self.recipe is None:
            return ""
        carry = spec.carry_for([self.recipe, src_rec])
        new, note = sm.carry_into(self.recipe, src_rec, carry, sm.amp_max(self.exp_amp),
                                  sm.amp_max(snap.get("exp_amp")), src_name=src_name)
        if new == self.recipe:
            return note
        self.snapshot()
        self.recipe = new
        self.on_structure_changed()
        return note

    # ------------------------------------------------------------- Fit → next
    def series_fit_then_next(self):
        """Fit the current member through the window's own Fit (FitWorker,
        progress bar, Stop / Cancel, animation), then move to the next member
        and seed it. The one-shot slot is connected AFTER run_fit's own, so
        ``_fit_done`` has landed the result before the move."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        if self._series_current_index() is None:
            self.statusBar().showMessage(
                "switch to a spectrum of the series first (click it in the bar)")
            return
        if self._series_busy():
            return
        before = getattr(self, "_fit_worker", None)
        self.run_fit()
        w = getattr(self, "_fit_worker", None)
        if w is None or w is before:
            return                        # run_fit refused and said why
        self._series_chain = w
        w.done.connect(self._series_chain_done)

    def _series_chain_done(self, result, stop_mode: str = ""):
        w = self._series_chain
        self._series_chain = None
        if w is None or w is not getattr(self, "_fit_worker", None):
            return
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        if stop_mode:
            self.statusBar().showMessage(
                ("fit stopped" if stop_mode == "stop" else "fit cancelled")
                + " — the chain stops here")
            return
        k = self._series_current_index()
        if k is None:
            return
        done = spec.members[k].name
        r = sm.rmsd_of(self.recipe)
        head = f"fitted {done}" + (f" (RMSD {r:.4g})" if r is not None else "")
        if k >= spec.n - 1:
            self.statusBar().showMessage(
                head + " — last spectrum of the series; Table… / Plot… / Save ▾ next")
            return
        self.series_go(k + 1)
        self.statusBar().showMessage(head + " → " + self.statusBar().currentMessage())

    # ------------------------------------------------------------- auto sweep
    def series_auto_sweep(self, passes=None, start=None, smooth=None):
        """Auto sweep ▾ Run: ``seqfit.run_sequential`` over every member in a
        SeqWorker; progress on the bar and the status-bar progress bar, whose
        Stop (keep) / Cancel (revert) buttons drive the worker."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        if self._series_busy():
            return
        data = self._series_member_data()
        if len(data) < 2:
            self.statusBar().showMessage("a sweep needs at least two spectra")
            return
        opts = spec.options
        passes = int(passes) if passes else int(opts.passes)
        start = str(start) if start in sm.START_CHOICES else str(opts.start)
        smooth = int(smooth) if smooth is not None and str(smooth).isdigit() else int(opts.smooth)
        try:
            entries = sm.entries_for_sweep(data, default_window=self.view.current_xrange())
        except ValueError as exc:
            self.statusBar().showMessage(str(exc), 15000)
            return
        carry = spec.carry_for([d["recipe"] for d in data])
        opts.passes, opts.start, opts.smooth = passes, start, smooth
        self._series_sweep_passes = passes
        from larmor.desktop.workers import SeqWorker, _fit_tol

        w = SeqWorker(entries, passes, start, carry, smooth, _fit_tol() or None)
        w.step.connect(self._series_sweep_step)
        w.done.connect(self._series_sweep_done)
        w.failed.connect(self._series_sweep_failed)
        self._seq_worker = w
        self._active_fit_worker = w
        self._show_fit_buttons(True)
        self._progress_start("series sweep")
        self.lines_table.btnFit.setEnabled(False)
        self.series_bar.set_sweep_running(True)
        self.series_bar.set_sweep_progress("starting…")
        self._series_retag()
        self.statusBar().showMessage(
            f"auto sweep: {passes} pass{'es' if passes != 1 else ''} from the {start} "
            f"spectrum over {len(data)} spectra …")
        w.start()

    def _series_sweep_step(self, p: int, k: int, rmsd: float):
        spec = getattr(self, "_series", None)
        n = spec.n if spec is not None else 0
        text = f"pass {p + 1}/{self._series_sweep_passes} · {k + 1}/{n} · RMSD {rmsd:.4g}"
        self.series_bar.set_sweep_progress(text)
        self.progress.setFormat("series sweep — " + text)

    def _series_sweep_cleanup(self, ok: bool):
        self._active_fit_worker = None
        self._progress_end(ok)
        self.lines_table.btnFit.setEnabled(True)
        self.series_bar.set_sweep_running(False)

    def _series_sweep_done(self, result, stop_mode: str = ""):
        """Apply the sweep: every member the sweep fitted takes its fitted
        recipe (``seriesmode.apply_sweep_result`` by series order) and a fit
        verdict built from the sweep's own curves -- the on-screen member
        live (one undo snapshot), the others in their workspace snapshots.
        Stop keeps what was fitted; Cancel changes nothing."""
        self._series_sweep_cleanup(stop_mode != "cancel")
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        if stop_mode == "cancel":
            self.statusBar().showMessage("sweep cancelled — every spectrum keeps its model")
            return
        rows = self._series_rows()
        data = self._series_member_data()
        if len(data) != len(getattr(result, "recipes", []) or []):
            self.statusBar().showMessage(
                "the series changed while the sweep ran — its result was not applied")
            return
        try:
            applied = sm.apply_sweep_result(data, result)
        except ValueError as exc:
            self.statusBar().showMessage(f"sweep result not applied: {exc}")
            return
        n_ok = 0
        for k, ((ws_i, _tag), upd) in enumerate(zip(rows, applied)):
            if upd is None:
                continue
            rec = upd["recipe"]
            m = spec.members[k]
            m.fit_sig, m.failed = sm.signature_key(rec), False
            n_ok += 1
            ws = self.workspaces[ws_i]
            snap = ws.get("snap") or {}
            h = self._series_health_for(rec, snap, upd)
            if ws_i == self.active_ws:
                self.snapshot()
                self.recipe = rec
                self.lines_table.rebuild(self.recipe, self.hidden)
                self._update_paddles()
                self._last_lmfit = None
                self._last_mc = self._last_mc_sig = None
                self._last_quant = None
                self.run_quantify(show=False)
                if h is not None:
                    self._health_fit = h
                    self._health_apply(h)
                    self.lines_table.set_chi2(h.chi_text())
                else:
                    self._health_reset()
                self.request_simulation()
            else:
                snap["recipe"] = rec
                snap["health"] = snap["health_fit"] = h
                snap["lmfit"] = None
                ws["has_fit"] = bool(rec.get("sites"))
        self._series_retag()
        self._refresh_ws_panel()
        self._series_refresh_bar()
        self._persist_session()
        means = [h.get("mean") for h in (getattr(result, "history", None) or [])]
        means = [m for m in means if m is not None and np.isfinite(m)]
        msg = result.summary
        if means:
            msg += " · pass means " + " → ".join(f"{m:.4g}" for m in means)
        if stop_mode == "stop":
            msg += f" · stopped early, {n_ok} of {len(rows)} fitted"
        self.statusBar().showMessage(msg, 20000)

    @staticmethod
    def _series_health_for(rec: dict, snap: dict, upd: dict):
        """A fit-health verdict for a swept member from the sweep's own model
        curve (no re-simulation): ``fithealth.assess`` on the data axis."""
        try:
            ppm = np.asarray(snap.get("exp_ppm"), float)
            amp = np.asarray(snap.get("exp_amp"), float)
            y = upd.get("y_fit")
            x = upd.get("x")
            if y is None or x is None:
                return None
            h = fithealth.assess(rec, amp, np.asarray(y, float),
                                 window=rec.get("fit_window_ppm"), rmsd=upd.get("rmsd"),
                                 fitted=True, ppm=ppm, x_fit=np.asarray(x, float))
            h.recipe_sig = fithealth.recipe_signature(rec)
            h.data_sig = fithealth.data_signature(ppm, amp)
            return h
        except Exception:                                  # noqa: BLE001
            return None

    def _series_sweep_failed(self, msg: str):
        self._series_sweep_cleanup(False)
        QMessageBox.warning(self, "Series sweep failed", humanize_error(msg))
        self.statusBar().showMessage("series sweep failed")

    # ------------------------------------------------------------- carry menu
    def _series_fill_carry_menu(self):
        """Before the Carry menu shows (and at series start): the checklist
        from every parameter present on the members, the ticks from the
        options, the sweep form from the options."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        recs = [d["recipe"] for d in self._series_member_data()]
        cands = sm.carry_candidates(recs)
        self.series_bar.set_options(spec.options.seed_on_move, cands,
                                    spec.carry_for(recs), spec.options.passes,
                                    spec.options.start, spec.options.smooth)

    def series_set_seed_on_move(self, on: bool):
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        spec.options.seed_on_move = bool(on)
        self._series_retag()
        self.statusBar().showMessage(
            "seed on move: " + ("on — the spectrum you move to takes this one's model"
                                if on else "off — moving leaves every spectrum as it is"))

    def series_set_carry(self, names):
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        spec.set_carry(tuple(names or ()), self.series_bar.carry_candidates())
        self._series_retag()
        from larmor.desktop.panels import PARAM_LABELS
        labels = [PARAM_LABELS.get(n, n) for n in (names or ())]
        self.statusBar().showMessage(
            "carried between spectra: " + (", ".join(labels) if labels else "nothing")
            + " (amplitudes are always re-fitted)")

    def series_copy_model(self, replace: bool = False):
        """Carry ▾ Copy this model to every spectrum (``replace=False``:
        empty members take a copy with scaled amplitudes, members with their
        own lines keep them and receive the carried values) / Replace every
        spectrum's lines with this model (``replace=True``)."""
        spec = getattr(self, "_series", None)
        if spec is None or self._series_busy():
            return
        k = self._series_current_index()
        if k is None:
            self.statusBar().showMessage("switch to the spectrum whose model to copy first")
            return
        self._sync_active()
        rows = self._series_rows()
        src_i = rows[k][0]
        src_snap = self.workspaces[src_i].get("snap") or {}
        src_rec = src_snap.get("recipe") or {}
        if not src_rec.get("sites"):
            self.statusBar().showMessage("this spectrum has no lines to copy")
            return
        src_amp = sm.amp_max(src_snap.get("exp_amp"))
        src_name = spec.members[k].name
        n_copied = n_seeded = 0
        for j, (ws_i, _tag) in enumerate(rows):
            if ws_i == src_i:
                continue
            ws = self.workspaces[ws_i]
            snap = ws.get("snap") or {}
            dst = snap.get("recipe") or {}
            had_lines = bool(dst.get("sites"))
            if replace and had_lines:
                dst = {**dst, "sites": []}
            carry = None if (replace or not had_lines) else spec.carry_for([dst, src_rec])
            new, _note = sm.carry_into(dst, src_rec, carry, sm.amp_max(snap.get("exp_amp")),
                                       src_amp, src_name=src_name)
            if had_lines and not replace:
                n_seeded += 1
            else:
                n_copied += 1
                spec.members[j].fit_sig, spec.members[j].failed = "", False
            snap["recipe"] = new
            ws["has_fit"] = bool(new.get("sites"))
        self._series_retag()
        self._refresh_ws_panel()
        self._series_refresh_bar()
        msg = f"model of {src_name} copied into {n_copied} spectr{'um' if n_copied == 1 else 'a'}"
        if n_seeded:
            msg += (f"; {n_seeded} with their own lines kept them and received the "
                    "carried values")
        self.statusBar().showMessage(msg, 12000)

    # ------------------------------------------------------------- outputs
    def series_table(self):
        """Table…: the Series table dialog; OK renames members, regroups and
        reorders the series (tags, bar and dock titles follow)."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        from PySide6.QtWidgets import QDialog

        from larmor.desktop.series_table_dialog import SeriesTableDialog

        w = getattr(self, "_seq_worker", None)
        busy = w is not None and w.isRunning()
        dlg = SeriesTableDialog(self, spec.series_table(), locked=busy)
        if dlg.exec() != QDialog.Accepted:
            return
        self._series_apply_table(dlg.table(), dlg.perm())

    def _series_apply_table(self, table, perm):
        """The dialog-free half of series_table (tests call it): ``perm``
        (new position -> current index) reorders the members; names go to
        the member, the workspace title and the recipe's sample."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        rows_before = self._series_rows()
        moved = spec.apply_table(table, perm)
        order = [rows_before[k][0] for k in perm]
        for k, ws_i in enumerate(order):
            m = spec.members[k]
            ws = self.workspaces[ws_i]
            ws["custom_title"] = m.name
            ws["title"] = m.name
            snap = ws.get("snap") or {}
            if isinstance(snap.get("recipe"), dict):
                snap["recipe"]["sample"] = m.name
            if ws_i == self.active_ws and isinstance(self.recipe, dict):
                self.recipe["sample"] = m.name
                self.view.set_title(m.name)
        self._series_retag()
        self._refresh_ws_panel()
        self._series_refresh_bar()
        self.statusBar().showMessage(f"series table applied — {spec.n} spectra"
                                     + (" in a new order" if moved else ""))

    def series_plot(self):
        """Plot…: the Series plot over the members that carry lines."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        from larmor.desktop.series_plot import SeriesPlotDialog
        from larmor.seqfit import SeqFitResult

        data = self._series_member_data()
        keep = [k for k, d in enumerate(data) if (d["recipe"] or {}).get("sites")]
        if not keep:
            self.statusBar().showMessage(
                "no lines on any spectrum yet — add a model, or Carry ▾ Copy this "
                "model to every spectrum")
            return
        recs = [Recipe.from_dict(data[k]["recipe"]) for k in keep]
        rmsd = [sm.rmsd_of(data[k]["recipe"]) for k in keep]
        res = SeqFitResult(
            recipes=recs, labels=[spec.members[k].name for k in keep],
            rmsd=[float("nan") if r is None else r for r in rmsd], per_dataset=[],
            history=[], passes=0, propagated=spec.carry_for([d["recipe"] for d in data]))
        table = spec.series_table()
        if len(keep) != len(data):
            table = table.permuted(keep)
        show_tool_window(SeriesPlotDialog(self, res, series=table))
        if len(keep) != len(data):
            self.statusBar().showMessage(
                f"series plot over {len(keep)} of {len(data)} spectra — the ones "
                "without lines are left out")

    def series_acquisition_table(self):
        """Acquisition…: the Experimental-section window over the members,
        with their confirmed spin rates (no fit needed)."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        from larmor.desktop.acquisition_dialog import AcquisitionTableDialog

        data = self._series_member_data()
        paths = [d["source_path"] for d in data if d["source_path"]]
        spin = {}
        for d in data:
            rec = d["recipe"] or {}
            acq = rec.get("acquisition") or {}
            spin[d["source_path"]] = (float(rec.get("spin_rate_Hz") or 0.0),
                                      bool(rec.get("mas_uncertain") or acq.get("mas_uncertain")))
        show_tool_window(AcquisitionTableDialog(self, paths, spin_rates=spin))

    def series_save_all(self):
        """Save ▾ Save all fits…: one .recipe.json per member into a folder,
        auto-named sample_nucleus_seq_YYYYMMDD_HHMM or typed one by one."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        from larmor.desktop.paths import suggest_save_dir

        data = self._series_member_data()
        if not data:
            return
        start = suggest_save_dir(data[0]["source_path"], self._last_dir())
        folder = QFileDialog.getExistingDirectory(self, "Save all fits of the series", start)
        if not folder:
            return
        mode = QMessageBox.question(
            self, "Naming", "Name the files automatically?\n\n"
            "Yes — auto (sample_nucleus_seq_YYYYMMDD_HHMM)\n"
            "No — type a name for each fit.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        folder = Path(folder)
        n = unfitted = 0
        for k, d in enumerate(data):
            if not d["recipe"]:
                continue
            rec = Recipe.from_dict(d["recipe"])
            rec.sample = spec.members[k].name
            if not rec.source_path and d["source_path"]:
                rec.source_path = d["source_path"]
            name = sm.auto_name(rec, spec.members[k].name)
            if mode != QMessageBox.Yes:
                text, ok = QInputDialog.getText(
                    self, "Fit name",
                    f"Name for “{spec.members[k].name}”   [{k + 1} of {len(data)}]",
                    text=name)
                if not ok:
                    break
                name = sm.slug(text) or name
            if sm.rmsd_of(d["recipe"]) is None:
                unfitted += 1
            try:
                rec.save(folder / f"{name}.recipe.json")
                n += 1
            except Exception as exc:                       # noqa: BLE001
                self.statusBar().showMessage(f"{name}: not saved — {exc}")
        self.statusBar().showMessage(
            f"saved {n} fit(s) to {folder}"
            + (f" — {unfitted} of them not fitted yet" if unfitted else ""), 12000)

    def series_bundle(self):
        """Save ▾ Publication bundle…: ``larmor.io.bundle.write_bundle``
        (kind 'seq') over the members' live recipes and the arrays as
        fitted; a member never fitted gets a manifest row 'not fitted'."""
        spec = getattr(self, "_series", None)
        if spec is None:
            return
        from larmor import batchfit
        from larmor.desktop.paths import suggest_save_dir
        from larmor.io import bundle

        data = self._series_member_data()
        if not data:
            return
        recs = []
        for k, d in enumerate(data):
            r = Recipe.from_dict(d["recipe"]) if d["recipe"] else Recipe()
            r.sample = spec.members[k].name
            recs.append(r)
        rmsd = [sm.rmsd_of(d["recipe"]) for d in data]
        res = batchfit.BatchFitResult(
            recipes=recs, labels=spec.names(),
            rmsd=[float("nan") if r is None else r for r in rmsd], per_dataset=[],
            shared=(), released=batchfit.all_but_amplitude(recs))
        start = suggest_save_dir(data[0]["source_path"], self._last_dir())
        folder = QFileDialog.getExistingDirectory(
            self, "Publication bundle — choose an output folder", start)
        if not folder:
            return
        folder = Path(folder)
        if (folder / "manifest.csv").exists():
            ans = QMessageBox.question(
                self, "Replace bundle?",
                f"“{folder.name}” already holds a bundle (manifest.csv). Replace its files?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if ans != QMessageBox.Yes:
                return
        arrays = [(d["ppm"], d["amp"], (d["recipe"] or {}).get("fit_window_ppm"))
                  for d in data]
        try:
            out = bundle.write_bundle(res, arrays, folder, kind="seq",
                                      source_paths=[d["source_path"] for d in data])
        except Exception as exc:                           # noqa: BLE001
            QMessageBox.warning(self, "Bundle failed", str(exc))
            return
        self.statusBar().showMessage(out.summary + (
            f" · {len(out.warnings)} warning(s), see README.txt" if out.warnings else ""),
            15000)

    def series_show_comparability(self):
        """The chip: every acqus / procs parameter of every member."""
        spec = getattr(self, "_series", None)
        if spec is None or spec.comparison is None:
            return
        from larmor.desktop.comparability_dialog import ComparabilityDialog

        show_tool_window(ComparabilityDialog(self, spec.comparison))

    # ------------------------------------------------------------- persistence
    def _series_restore_from_tags(self):
        """After open_project: the 1D entries whose tags share one series id
        (the largest group) become the active series again, with their
        names, order, options, fit state and Series table; the acqus /
        procs parameters are re-read for the comparability chip."""
        if getattr(self, "_series", None) is not None and self._series_rows():
            return
        groups: dict = {}
        for i, ws in enumerate(self.workspaces):
            tag = ws.get("series")
            if ws.get("kind") == "1d" and isinstance(tag, dict) and tag.get("id"):
                groups.setdefault(str(tag["id"]), []).append(i)
        if not groups:
            return
        _sid, idx = max(groups.items(), key=lambda kv: len(kv[1]))
        tags = [self.workspaces[i]["series"] for i in idx]
        params = []
        for i in idx:
            snap = self.workspaces[i].get("snap") or {}
            p = snap.get("source_path") or self.workspaces[i]["series"].get("source_path")
            try:
                from larmor import comparability
                params.append(comparability.read_params(p) if p else None)
            except Exception:                              # noqa: BLE001
                params.append(None)
        try:
            spec = sm.SeriesSpec.from_tags(tags, params)
        except ValueError:
            return
        # the sources may have moved with the project: the snapshots know
        for i in idx:
            tag = self.workspaces[i]["series"]
            m = spec.member(tag.get("key"))
            src = (self.workspaces[i].get("snap") or {}).get("source_path")
            if m is not None and src:
                m.source_path = str(src)
        self._series = spec
        self._series_retag()
        self._refresh_ws_panel()
        self._update_enabled()
        self._series_fill_carry_menu()
        self._series_refresh_bar()
