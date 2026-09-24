"""Main-window mixin: workspaces, the project bundle and the session file.

Workspace rows (register / switch / close / save), the ``.larproj.json``
project bundle (``save_project`` / ``open_project`` over ``larmor.project``),
the figure / batch session rows, and the crash-recovery session file
(``_persist_session`` debounced, ``_flush_session``, ``_autosave_tick``,
``_restore_session``). ``closeEvent`` stays in ``app.py`` and calls
``_flush_session``.

Owned state: ``workspaces``, ``active_ws``, ``_ws_mode``, the session and
autosave timers.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QFileDialog, QMessageBox

from larmor import display, project
from larmor.project import DOC_KINDS, SESSION_KINDS


class _SessionMixin:
    """Workspaces, project bundle, session file."""

    # ------------------------------------------------------- workspaces
    def _doc_title(self) -> str:
        if self.central_stack.currentWidget() is self.view2d:
            return self.view2d.title.text() or "2D"
        if self.recipe:
            return (self.recipe.get("sample")
                    or (Path(self.source_path).name if self.source_path else "spectrum"))
        return "empty"

    # ---------- project (the whole session: spectra, 2D maps, figures, batch fits) ----------
    def save_project(self):
        """Save every workspace -- 1D spectra (processing + fit + overlays by
        reference), 2D maps by reference with their recorded processing, kept
        figures and batch-fit sessions -- as one reopenable project file. The
        schema and the per-kind entry builders live in larmor.project."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save project", self._suggest_dir(),
            "LARMOR project (*.larproj.json)")
        if not path:
            return
        self._sync_active()
        bundle, dropped = project.build_bundle(
            self.workspaces, self.active_ws, str(Path(path).parent))
        Path(path).write_text(json.dumps(bundle), encoding="utf-8")
        self._add_recent(path)
        msg = f"project saved — {project.summary(bundle)}"
        if dropped:
            msg += (f"  ({dropped} workspace(s) without a reloadable source "
                    "not included)")
        self.statusBar().showMessage(msg)

    def open_project(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open project", self._last_dir(),
            "LARMOR project (*.larproj.json *.json)")
        if not path:
            return
        self._open_project_path(path)

    def _open_project_path(self, path: str):
        """Open a project bundle (also reached by File ▸ Open…, drag-drop and
        Open recent through load_source): migrate it (larmor.project), ask
        Replace / Add / Cancel when a session with a fit or a non-1D
        workspace is already open, restore every entry by kind, then report
        everything that was migrated, relocated or missing in ONE box -- and
        only when there is something to say."""
        try:
            data, notes = project.load_bundle(path)
        except Exception as exc:
            QMessageBox.warning(self, "Open project", f"cannot read: {exc}")
            return
        wss = data.get("workspaces", [])
        if not wss:
            self.statusBar().showMessage("no workspaces in this project")
            return
        offset = 0
        if any(ws["has_fit"] or ws["kind"] != "1d" for ws in self.workspaces):
            # today's silent `self.workspaces = []` destroyed unsaved work
            mode = self._confirm_open_mode()
            if mode == "cancel":
                self.statusBar().showMessage("open project cancelled")
                return
            if mode == "add":
                self._sync_active()
                offset = len(self.workspaces)
            else:
                self.workspaces = []
                self.active_ws = None
        else:
            self.workspaces = []            # bare unfitted spectra: replace
            self.active_ws = None
        project_dir = str(Path(path).parent)
        for w in wss:
            kind = w.get("kind", "1d")
            title = w.get("title") or kind
            if kind == "1d":
                self._restore_1d_entry(w, notes)
            elif kind == "2d":
                self._restore_2d_entry(w, project_dir, notes)
            elif kind == "figure":
                spec = w.get("spec")
                if isinstance(spec, dict):
                    self._add_session_entry("figure", title, {"spec": spec})
                else:
                    notes.append(f"figure '{title}': no spec — skipped")
            elif kind == "batch":
                st = w.get("state")
                if not isinstance(st, dict):
                    notes.append(f"batch session '{title}': no state — skipped")
                    continue
                st, missing = project.relocate_batch_state(st, project_dir)
                if missing:
                    notes.append(
                        f"batch session '{title}': {len(missing)} spectrum "
                        "file(s) not found — " + ", ".join(missing))
                self._add_session_entry("batch", title, {"state": st})
            else:
                notes.append(f"workspace '{title}' of unknown kind '{kind}' "
                             "— skipped")
        act = data.get("active")
        if (isinstance(act, int) and 0 <= act + offset < len(self.workspaces)
                and self.workspaces[act + offset]["kind"] in DOC_KINDS):
            self.switch_workspace(act + offset)
        self._refresh_ws_panel()
        self._add_recent(path)
        self.statusBar().showMessage(f"project opened — {project.summary(data)}")
        self._report_project_notes(notes)

    def _confirm_open_mode(self) -> str:
        """Open project into a non-empty session: 'replace' | 'add' | 'cancel'.
        Its own method so tests patch it (a modal box hangs offscreen)."""
        box = QMessageBox(self)
        box.setWindowTitle("Open project")
        box.setIcon(QMessageBox.Question)
        box.setText("A session is already open.")
        box.setInformativeText("Replace it with the project, or add the "
                               "project's workspaces to it?")
        b_rep = box.addButton("Replace current session", QMessageBox.AcceptRole)
        b_add = box.addButton("Add to current session", QMessageBox.ActionRole)
        b_cancel = box.addButton(QMessageBox.Cancel)
        box.setEscapeButton(b_cancel)
        box.exec()
        clicked = box.clickedButton()
        if clicked is b_rep:
            return "replace"
        if clicked is b_add:
            return "add"
        return "cancel"

    def _report_project_notes(self, notes: list):
        """Everything worth knowing about an open -- schema migrated, sources
        relocated or missing, overlays / projections dropped, 'run Fit' -- in
        one box, only when there are notes. Its own method so tests patch it."""
        if notes:
            QMessageBox.information(self, "Project opened", "\n".join(notes))

    def _restore_1d_entry(self, w: dict, notes: list):
        """A 1D entry: embedded spectrum, recipe, hidden sites and overlays by
        reference -- the v1 body, with the missing-overlay report collected
        into the notes instead of a status message."""
        self.exp_ppm = np.asarray(w.get("exp_ppm", []), float)
        self.exp_amp = np.asarray(w.get("exp_amp", []), float)
        self.recipe = w.get("recipe")
        self.hidden = set(w.get("hidden", []))
        self.source_path = w.get("source_path")
        self._proc_base = None
        self._overlays = []
        self.central_stack.setCurrentWidget(self.view)
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.view.set_title(w.get("title") or (
            self.recipe.get("sample") if self.recipe else ""))
        self.lines_table.rebuild(self.recipe, self.hidden)
        if self.recipe and self.recipe.get("sites"):
            self.request_simulation()
        missing_overlays = []
        for ov in w.get("overlays", []):
            src = ov.get("source", "")
            try:
                label, ppm, amp, info = self._read_overlay_source(src)
            except Exception:
                missing_overlays.append(ov.get("label") or src)
                continue
            self._add_overlay(ov.get("label") or label, ppm, amp, src, info)
            new = self._overlays[-1]
            new["visible"] = bool(ov.get("visible", True))
            if ov.get("color"):
                new["color"] = ov["color"]
            # the display transform (scale / shift / yoff) -- neutral values
            # for a bundle written before it existed
            for key, default in display.OVERLAY_DEFAULTS.items():
                new[key] = project.num_or(ov.get(key), default)
        if w.get("overlays"):
            self._refresh_overlays()          # reflect restored visible/color
        if missing_overlays:
            notes.append("overlays could not be relocated: "
                         + ", ".join(missing_overlays))
        self._ws_mode = "new"
        self._register_ws("1d")
        self._update_paddles(); self._update_exp_label(); self._update_enabled()

    def _restore_2d_entry(self, w: dict, project_dir: str, notes: list):
        """A 2D entry, BY REFERENCE: relocate the source (as saved, then
        project-relative), reopen it through the production path
        (_load_nonfittable -> _show_2d -> _register_ws), then lay the saved
        recipe, contour settings, recorded processing (twod.replay_ops on the
        fresh load) and projection overlays on top. The fitted model overlay
        is not re-simulated (the MQMAS kernel build takes seconds); the
        parameters come back and the notes say to run Fit."""
        title = w.get("title") or "2D"
        src = project.relocate(w.get("source_path"), project_dir,
                               w.get("source_rel"))
        if src is None:
            notes.append(f"2D map '{title}': source not found "
                         f"({w.get('source_path')}) — skipped")
            return
        self.recipe = None       # _show_2d seeds a blank recipe only when none
        self._ws_mode = "new"
        self._in_load_source = True
        try:
            ok = self._load_nonfittable(src)
        finally:
            self._in_load_source = False
        if not ok:
            self._ws_mode = "auto"
            notes.append(f"2D map '{title}': source could not be read ({src}) "
                         "— skipped")
            return
        self.recipe = w.get("recipe") or self.recipe
        self.hidden = set(w.get("hidden", []))
        self._data2d_fittable = bool(
            w.get("fittable", getattr(self, "_data2d_fittable", False)))
        view = w.get("view") or {}
        try:
            self.view2d.apply_persisted_view(view)
        except Exception as exc:
            notes.append(f"2D map '{title}': processing could not be replayed "
                         f"({exc}) — showing the raw map")
        for axis, pj in (view.get("projections") or {}).items():
            s1 = project.relocate(pj.get("source", ""), project_dir)
            try:
                if s1 is None:
                    raise FileNotFoundError(pj.get("source", ""))
                ppm, amp = self._read_1d(s1)
                self.view2d.set_projection_1d(
                    axis, ppm, amp,
                    color=pj.get("color") or self.PROJ_COLOR.get(axis),
                    source=s1, scale=pj.get("scale"))
            except Exception:
                notes.append(f"2D map '{title}': {axis.upper()} projection "
                             f"overlay not found ({pj.get('source')})")
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._sync_active()
        self._update_exp_label(); self._update_enabled()
        if self.recipe and self.recipe.get("sites"):
            notes.append(f"2D map '{title}': fit parameters restored — run Fit "
                         "to redraw the model overlay")

    def _snapshot_doc(self) -> dict:
        is2d = self.central_stack.currentWidget() is self.view2d
        snap = {"kind": "2d" if is2d else "1d",
                "source_path": self.source_path,
                "recipe": json.loads(json.dumps(self.recipe)) if self.recipe else None,
                "hidden": set(self.hidden)}
        if is2d:
            snap["data2d"] = self._data2d
            snap["view2d"] = self.view2d.get_state()
            snap["fittable"] = getattr(self, "_data2d_fittable", False)
        else:
            snap["exp_ppm"] = np.array(self.exp_ppm, copy=True)
            snap["exp_amp"] = np.array(self.exp_amp, copy=True)
            snap["proc_base"] = self._proc_base
            snap["overlays"] = [dict(o) for o in self._overlays]
            # in-memory only (like the arrays above; never serialised): the
            # verdict and covariance of this workspace's last fit, so that
            # switching back restores them instead of reading 'no fit yet'
            snap["health"] = self._health
            snap["health_fit"] = self._health_fit
            snap["lmfit"] = self._last_lmfit
        return snap

    def _apply_doc(self, snap: dict):
        self.source_path = snap["source_path"]
        self.recipe = (json.loads(json.dumps(snap["recipe"]))
                       if snap["recipe"] else None)
        self.hidden = set(snap["hidden"])
        self.undo_stack.clear(); self.redo_stack.clear()
        if snap["kind"] == "2d":
            self._data2d = snap["data2d"]
            self._data2d_fittable = snap.get("fittable", False)
            self.view2d.set_state(snap["view2d"])
            self.central_stack.setCurrentWidget(self.view2d)
            self._health_reset()          # a 2D workspace carries no 1D verdict
        else:
            self.exp_ppm = snap["exp_ppm"]; self.exp_amp = snap["exp_amp"]
            self._proc_base = snap["proc_base"]
            self._overlays = list(snap["overlays"])
            self.central_stack.setCurrentWidget(self.view)
            self.view.set_experiment(self.exp_ppm, self.exp_amp)
            self.view.set_title(self.recipe.get("sample", "") if self.recipe else "")
            self.lines_table.rebuild(self.recipe, self.hidden)
            self._update_paddles(); self._refresh_overlays(); self._update_sn()
            # _update_sn dropped the verdict; put this workspace's own back.
            # The re-simulation below finds equal signatures and keeps it.
            self._health = snap.get("health")
            self._health_fit = snap.get("health_fit")
            self._last_lmfit = snap.get("lmfit")
            self.health_strip.set_health(self._health)
            self._health_show()
            if self.recipe and self.recipe.get("sites"):
                self.request_simulation()
            else:
                self.view.set_model(None, None, None, None, self.hidden)
        self._update_exp_label(); self._update_enabled()

    def _sync_active(self):
        if self.active_ws is None or not (0 <= self.active_ws < len(self.workspaces)):
            return
        ws = self.workspaces[self.active_ws]
        ws["snap"] = self._snapshot_doc()
        ws["kind"] = ws["snap"]["kind"]
        ws["has_fit"] = bool(self.recipe and self.recipe.get("sites"))
        ws["title"] = self._doc_title()

    def _register_ws(self, kind: str):
        mode, self._ws_mode = self._ws_mode, "auto"
        entry = {"snap": self._snapshot_doc(), "kind": kind,
                 "title": self._doc_title(),
                 "has_fit": bool(self.recipe and self.recipe.get("sites"))}
        reuse = False
        if mode == "reuse":
            reuse = self.active_ws is not None
        elif mode == "auto" and self.active_ws is not None:
            cur = self.workspaces[self.active_ws]
            reuse = (kind == "1d" and cur["kind"] == "1d" and not cur["has_fit"])
        if reuse:
            self.workspaces[self.active_ws] = entry
        else:
            self.workspaces.append(entry)
            self.active_ws = len(self.workspaces) - 1
        self._refresh_ws_panel()

    def _add_session_entry(self, kind: str, title: str, payload: dict) -> int:
        """Append a figure / batch-session row to the Workspaces dock (the
        same list as the documents, so row == index holds). Never touches
        active_ws: a session row is never the displayed document."""
        self.workspaces.append({"snap": {"kind": kind, **payload}, "kind": kind,
                                "title": title, "has_fit": False})
        self._refresh_ws_panel()
        return len(self.workspaces) - 1

    def _refresh_ws_panel(self):
        items = []
        for ws in self.workspaces:
            icon = ({"2d": "▦", "figure": "◫", "batch": "☷"}.get(ws["kind"])
                    or ("⤳" if ws["has_fit"] else "∿"))
            items.append((icon, ws["title"]))
        self.ws_panel.rebuild(items, self.active_ws if self.active_ws is not None
                              else -1)

    def switch_workspace(self, i: int):
        if i == self.active_ws or not (0 <= i < len(self.workspaces)):
            return
        if self.workspaces[i]["kind"] in SESSION_KINDS:
            # a figure / batch row only selects: opening its modal dialog is
            # an explicit gesture (Enter, double-click, Save) -- see
            # open_workspace_entry. active_ws never points at a session row.
            self.statusBar().showMessage(
                f"{self.workspaces[i]['title']} — Enter / double-click (or "
                "Save) opens it")
            return
        self._sync_active()
        self.active_ws = i
        self._apply_doc(self.workspaces[i]["snap"])
        self._refresh_ws_panel()
        self.statusBar().showMessage(f"workspace: {self.workspaces[i]['title']}")

    def close_workspace(self, i: int):
        if not (0 <= i < len(self.workspaces)):
            return
        del self.workspaces[i]
        if not self.workspaces:
            self.active_ws = None
            self._refresh_ws_panel()
            return
        if self.active_ws == i:
            # the nearest DOCUMENT takes over; session rows cannot be displayed
            docs = [j for j, ws in enumerate(self.workspaces)
                    if ws["kind"] in DOC_KINDS]
            if not docs:
                self.active_ws = None
            else:
                self.active_ws = min(docs, key=lambda j: abs(j - i))
                self._apply_doc(self.workspaces[self.active_ws]["snap"])
        elif self.active_ws is not None and self.active_ws > i:
            self.active_ws -= 1
        self._refresh_ws_panel()

    def open_workspace_entry(self, i: int):
        """Enter / double-click on a Workspaces row: a document switches, a
        figure / batch session reopens in its dialog."""
        if not (0 <= i < len(self.workspaces)):
            return
        if self.workspaces[i]["kind"] in SESSION_KINDS:
            self._open_session_entry(i)
        else:
            self.switch_workspace(i)

    def save_workspace(self, i: int):
        if not (0 <= i < len(self.workspaces)):
            return
        if self.workspaces[i]["kind"] in SESSION_KINDS:
            self._open_session_entry(i)      # the Save button reopens a session
            return
        self.switch_workspace(i)
        if self.central_stack.currentWidget() is self.view2d:
            self.save_recipe()
        elif self.recipe and self.recipe.get("sites"):
            self.save_recipe()
        else:
            self.save_spectrum()

    # ------------------------------------------- figures & batch sessions
    # Figures and batch-fit sessions are rows of the Workspaces dock (kind
    # "figure" / "batch") so File ▸ Save project… stores them. Their dialogs
    # are modal, so a row is recorded when the dialog closes and updated in
    # place when it is reopened from the dock and closed again.
    def _open_session_entry(self, i: int):
        """Reopen a figure (Plotting studio, or the Figure studio for a kind
        the studio cannot edit) or a batch session in its dialog."""
        ws = self.workspaces[i]
        if ws["kind"] == "figure":
            spec = ws["snap"].get("spec") or {}
            from larmor.desktop.plotting_studio import EDITABLE_KINDS
            if spec.get("kind", "1d") in EDITABLE_KINDS:
                self.open_plotting_studio(spec, ws_index=i)
            else:
                from larmor.desktop.figure_dialog import FigureDialog
                dlg = FigureDialog(self, self.source_path, self.recipe)
                dlg.spec.setPlainText(json.dumps(spec, indent=2))
                dlg.exec()
                try:
                    new = json.loads(dlg.spec.toPlainText())
                except ValueError:
                    new = None
                if isinstance(new, dict):
                    self._remember_figure(new, i)
        elif ws["kind"] == "batch":
            self._open_batch_session(ws["snap"].get("state") or {}, i)
        self._refresh_ws_panel()

    def _remember_figure(self, spec, ws_index=None):
        """Keep a drawable figure spec as a Workspaces row -- updating row
        ``ws_index`` in place when given, else appending. An empty spec (no
        traces / EXPNO / panels / categories) is never kept."""
        from larmor.desktop.plotting_studio import spec_is_empty

        if not isinstance(spec, dict) or spec_is_empty(spec):
            return
        spec = project.json_safe(spec)          # tuples -> lists, like the file
        title = spec.get("title") or f"{spec.get('kind', '1d')} figure"
        if (ws_index is not None and 0 <= ws_index < len(self.workspaces)
                and self.workspaces[ws_index]["kind"] == "figure"):
            ws = self.workspaces[ws_index]
            ws["snap"]["spec"] = spec
            ws["title"] = title
        else:
            self._add_session_entry("figure", title, {"spec": spec})
        self._refresh_ws_panel()
        self.statusBar().showMessage(
            "figure kept in the Workspaces dock — File ▸ Save project… stores it")

    def _remember_batch(self, dlg, ws_index=None):
        """Keep a batch-fit dialog's session as a Workspaces row once it has a
        result (a configured-but-unfitted dialog adds nothing). A row reopened
        from the dock is updated in place; if its dialog closed WITHOUT a
        result (e.g. every spectrum file was missing) the stored session --
        results included -- is left as it was rather than overwritten."""
        if dlg._result is None:
            if ws_index is not None:
                self.statusBar().showMessage(
                    "batch session left unchanged — no fit result to record")
            return
        state = dlg.session_state()
        title = dlg.session_title()
        if (ws_index is not None and 0 <= ws_index < len(self.workspaces)
                and self.workspaces[ws_index]["kind"] == "batch"):
            ws = self.workspaces[ws_index]
            ws["snap"]["state"] = state
            ws["title"] = title
        else:
            self._add_session_entry("batch", title, {"state": state})
        self._refresh_ws_panel()
        self.statusBar().showMessage(
            "batch session kept in the Workspaces dock — File ▸ Save project… "
            "stores it")

    def _open_batch_session(self, state: dict, ws_index: int):
        """Rebuild a batch-fit dialog from a stored session (spectra, model,
        settings, results re-aligned to the spectra that load) and record it
        back into its row on close."""
        from larmor.desktop.batchfit_dialog import BatchFitDialog

        model = ({"sites": state.get("model_sites"),
                  "fit_window_ppm": state.get("window")}
                 if state.get("model_sites") else None)
        dlg = BatchFitDialog(self, list(state.get("paths") or []), model)
        notes = dlg.apply_session_state(state)
        if notes:
            self.statusBar().showMessage("; ".join(notes))
        dlg.exec()
        self._remember_batch(dlg, ws_index)

    # ------------------------------------------------------------- session
    @staticmethod
    def _session_file() -> Path:
        """The crash-recovery session lives in a FILE, not the registry.

        It used to be two QSettings values -- on Windows that is
        HKEY_CURRENT_USER, and with a spectrum component in the recipe the
        value was a 0.66 MB JSON string written into the registry on every
        user action. A registry is no place for megabytes of float lists."""
        base = os.environ.get("LOCALAPPDATA") or str(Path.home())
        d = Path(base) / "LARMOR"
        d.mkdir(parents=True, exist_ok=True)
        return d / "session.json"

    def _persist_session(self):
        """Schedule a debounced session write (5 s after the LAST action) --
        serialising the whole recipe on every click measured 15 ms each.
        Also newly gated on LARMOR_NO_SESSION: the write path was unguarded,
        so every TEST that called snapshot() wrote the developer's real
        session keys."""
        if not self.source_path or os.environ.get("LARMOR_NO_SESSION"):
            return
        self._session_timer.start()

    def _flush_session(self):
        """Write the session file now (close, autosave tick, debounce fire)."""
        if not self.source_path or os.environ.get("LARMOR_NO_SESSION"):
            return
        try:
            target = self._session_file()
            tmp = target.with_suffix(".json.tmp")
            tmp.write_text(json.dumps({"source": self.source_path,
                                       "recipe": self.recipe or {}}),
                           encoding="utf-8")
            tmp.replace(target)                     # atomic on the same volume
        except OSError:
            pass

    def _autosave_tick(self):
        """Periodic crash-safe autosave (see the autosave QTimer)."""
        if self.source_path and not (self._active_fit_worker
                                     and self._active_fit_worker.isRunning()):
            self._flush_session()
            self.statusBar().showMessage("session autosaved", 1500)

    def _restore_session(self):
        # tests (and clean-room launches) opt out so they never inherit a
        # developer's real saved recipe
        if os.environ.get("LARMOR_NO_SESSION"):
            return
        src, recipe = "", {}
        try:
            f = self._session_file()
            if f.exists():
                d = json.loads(f.read_text(encoding="utf-8"))
                src, recipe = d.get("source", ""), d.get("recipe", {}) or {}
        except (OSError, ValueError):
            pass
        if not src:
            # migrate a pre-0.12 registry session once, then remove the keys
            s = QSettings("LARMOR", "app")
            src = s.value("session/source", "")
            saved = s.value("session/recipe", "")
            try:
                recipe = json.loads(saved) if saved else {}
            except ValueError:
                recipe = {}
            if src:
                s.remove("session/source")
                s.remove("session/recipe")
        if not src or not Path(src).exists():
            return
        try:
            self.load_source(src)
            if recipe.get("sites"):
                self.recipe = recipe
                self.on_structure_changed()
            self.statusBar().showMessage(f"session restored — {src}")
        except Exception:
            pass
