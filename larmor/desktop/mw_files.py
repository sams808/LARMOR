"""Main-window mixin: every way data enters or leaves the window.

Open dialogs (file / EXPNO / Varian / sample / FID), ``load_source`` and its
body, the explorer and 2D navigation (``_explorer_open``, ``_show_2d``,
``back_to_2d``), the 1D display path, File > Watch, the kernel pre-build
trigger and the save side (recipe, fit-as, spectrum, figure). The module
functions ``_load_any`` and ``_nmrdata_to_data2d`` live here; ``app.py``
re-exports ``_load_any`` for the dialogs that import it lazily.

Owned state: ``source_path``, ``_watch_path`` (the watcher and its timer are
created in ``__init__``), ``_warmed_key`` / ``_warm_worker``, ``_data2d`` /
``_data2d_fittable``, ``_pending_apply_model``, ``_in_load_source``.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

from larmor.desktop.workers import KernelWarmWorker, humanize_error
from larmor.recipe import Recipe


def _load_any(path: str, replay: bool = True):
    """Load any supported source (shared with the CLI and figure studio).
    ``replay=False`` returns a recipe's source data WITHOUT its recorded
    processing (the unprocessed base the live pipeline re-applies from)."""
    from larmor.loader import load_any

    return load_any(path, replay=replay)


def keep_fit_prompt(existing_nucleus, incoming_nucleus) -> tuple[str, bool]:
    """Text and default answer of the "Keep fit parameters?" question.

    Keeping the open lines for a new spectrum is the dmfit habit for the SAME
    nucleus (a composition series). Across nuclei it puts, say, 27Al Czjzek
    sites on a 31P spectrum -- the prompt used to default to Yes and swap
    the nucleus silently -- so the default flips to No and the text says why.
    """
    a = str(existing_nucleus or "").strip()
    b = str(incoming_nucleus or "").strip()
    if not a or not b or a == b:
        return ("A fit is already open. Keep the current lines and fit them "
                "against the new spectrum?\n\n"
                "Yes = keep the lines · No = start empty", True)
    return (f"A fit is already open, but it is a {a} model and this spectrum "
            f"is {b}. Keeping the lines would fit {a} sites to {b} data.\n\n"
            "No = start empty (recommended) · Yes = keep the lines anyway", False)


def _nmrdata_to_data2d(data):
    """Convert a frequency-domain 2D NMRData into a twod.Data2D for display."""
    from larmor.twod import Data2D

    f1, f2 = data.axes
    h = data.hyper or {}
    d = Data2D(f2_ppm=np.asarray(f2.values), f1_ppm=np.asarray(f1.values),
               z=np.asarray(data.data, float), nucleus=data.nucleus,
               larmor_MHz=data.meta.get("larmor_MHz", 0.0),
               source=data.source,
               ri=h.get("ri"), ir=h.get("ir"), ii=h.get("ii"))
    d.notes = list(data.warnings)
    if data.is_pseudo2d:
        d.notes.append("pseudo-2D (arrayed)")
    return d


class _FilesMixin:
    """Opening, loading, navigating and saving data (see module docstring)."""

    # ------------------------------------------------------------- watch
    def _toggle_watch(self, on: bool):
        if on and not self.source_path:
            self.actWatch.setChecked(False)
            self.statusBar().showMessage("open a spectrum first, then watch it")
            return
        self._retarget_watch()
        self.statusBar().showMessage(
            f"watching {Path(self.source_path).name}: the spectrum reloads "
            "when the file changes, the fit is kept" if on
            else "stopped watching the source file")

    def _watch_targets(self) -> list[str]:
        """The file(s)/folder(s) whose change means new data: the source
        itself and its folder (a rewritten file loses its own watch)."""
        p = Path(self.source_path) if self.source_path else None
        if p is None or not p.exists():
            return []
        if p.is_dir():                     # an EXPNO: the processed data lives below
            cands = [p / "pdata" / "1" / "1r", p / "pdata" / "1" / "2rr",
                     p / "fid", p / "ser"]
            files = [str(c) for c in cands if c.exists()]
            return files + [str(c.parent) for c in cands if c.exists()] + [str(p)]
        return [str(p), str(p.parent)]

    def _retarget_watch(self):
        if self._watcher.files():
            self._watcher.removePaths(self._watcher.files())
        if self._watcher.directories():
            self._watcher.removePaths(self._watcher.directories())
        self._watch_path = None
        if getattr(self, "actWatch", None) is None or not self.actWatch.isChecked():
            return
        targets = self._watch_targets()
        if targets:
            self._watcher.addPaths(targets)
            self._watch_path = self.source_path

    def _watched_changed(self, *_):
        if self._watch_path:
            self._watch_timer.start()

    def _reload_watched(self):
        path = self._watch_path
        if not path or not Path(path).exists():
            return
        try:
            self.load_source(path, keep_fit=True)
            self.statusBar().showMessage(
                f"reloaded {Path(path).name} (watching for changes)")
        except Exception as exc:      # the writer may be mid-way: retry later
            self.statusBar().showMessage(f"reload failed, will retry: {exc}")
            self._watch_timer.start()
        finally:
            self._retarget_watch()    # a replaced file needs a fresh watch

    # ------------------------------------------------------------- loading
    def open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open spectrum, recipe, or Bruker file", self._last_dir(),
            "All supported (*.fxmla *.fxml *.json 1r 2rr fid ser);;"
            "dmfit / recipe (*.fxmla *.fxml *.json);;"
            "Bruker processed (1r 2rr);;Bruker raw (fid ser);;All files (*)")
        if path:
            self.load_source(path)

    def open_expno(self):
        path = QFileDialog.getExistingDirectory(
            self, "Open Bruker EXPNO or pdata folder (read-only)",
            self._last_dir())
        if path:
            self.load_source(path)

    def open_varian(self):
        path = QFileDialog.getExistingDirectory(
            self, "Open a Varian / Agilent .fid folder (procpar + fid)",
            self._last_dir())
        if path:
            self.load_source(path)

    def open_sample(self):
        path = QFileDialog.getExistingDirectory(
            self, "Choose a sample folder (lists all its spectra)",
            self._last_dir())
        if path:
            self.explorer.load_sample(path)
            self.explorer_dock.show()
            self.explorer_dock.raise_()
            self.statusBar().showMessage(
                "sample scanned — double-click a spectrum in the Explorer "
                "to open it")

    def open_fid(self):
        from larmor.desktop.fid_dialog import FidDialog

        path = None
        if self.source_path and Path(self.source_path).is_dir() and \
                ((Path(self.source_path) / "fid").exists() or
                 (Path(self.source_path) / "ser").exists()):
            path = str(Path(self.source_path) /
                       ("ser" if (Path(self.source_path) / "ser").exists()
                        else "fid"))
        dlg = FidDialog(self, path)
        dlg.accepted_1d.connect(self._fid_to_workbench)
        dlg.accepted_2d.connect(self._fid_to_2d)
        dlg.exec()

    def _fid_to_workbench(self, ppm, amp, meta):
        """A 1D spectrum processed from a raw fid becomes the working data."""
        from larmor.recipe import Recipe

        order = np.argsort(ppm)
        self.exp_ppm, self.exp_amp = np.asarray(ppm)[order], np.asarray(amp)[order]
        self._proc_base = None
        self.source_path = meta.get("expno", "")
        # the FULL processing record (every qcpmg_* key) rides along so a
        # saved fit of a QCPMG spectrum still says how it was made
        provenance = {k: v for k, v in meta.items() if k.startswith("qcpmg_")}
        spin_rate = meta.get("spin_rate_Hz") or meta.get("masr_Hz") or 0.0
        mas_uncertain = bool(meta.get("mas_uncertain", False))
        if "mas_sources" in meta:
            # a spectrum processed from an EXPNO's fid: the same three-source
            # resolution and per-session confirmation as File > Open (QCPMG
            # meta carries no sources block and keeps the values above)
            from larmor import masrate

            mas = masrate.resolve_for_load(meta, meta.get("expno", ""))
            spin_rate, mas_uncertain = mas["spin_rate_Hz"], mas["mas_uncertain"]
            provenance["mas_rate"] = mas["provenance"]
            if mas["note"]:
                self.statusBar().showMessage("⚠ " + mas["note"])
        self.recipe = Recipe(
            sample=(meta.get("title", "").splitlines() or [""])[0],
            source_kind="bruker", source_path=meta.get("expno", ""),
            nucleus=meta.get("nucleus", ""),
            larmor_frequency_MHz=meta.get("larmor_MHz", 0.0),
            spin_rate_Hz=spin_rate,
            mas_uncertain=mas_uncertain,
            # the Open FID dialog records the chain that produced its result
            # (window, zf, ft, phase | autophase): the fit is reproducible and
            # FID ⇄ spectrum replays it from the instrument fid. QCPMG's meta
            # has no such record and keeps an empty (pdata) chain, as before.
            processing=list(meta.get("processing") or []),
            processing_from_raw=bool(meta.get("processing_from_raw", False)),
            provenance=provenance).to_dict()
        self.hidden.clear(); self.undo_stack.clear(); self.redo_stack.clear()
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.view.set_title(self.recipe.get("sample") or "processed FID")
        self._last_model = None
        self._first_sim = True
        self.zoom_full()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles(); self._update_exp_label(); self._update_enabled()
        self._update_sn()          # S/N readout + the sideband offer (F3)
        self.statusBar().showMessage(
            "spectrum from processed FID loaded — add lines and fit")

    def _fid_to_2d(self, data2d):
        from larmor.desktop.twod_dialog import TwoDDialog

        dlg = TwoDDialog(self, None)
        dlg.data = data2d.normalized()
        dlg.lbl.setText(data2d.source or "processed 2D FID")
        dlg._redraw()
        dlg.exec()

    def _last_dir(self) -> str:
        return QSettings("LARMOR", "app").value("lastDir", "")

    def _suggest_dir(self) -> str:
        """Save dialogs follow the CURRENT dataset (its sample folder),
        not whichever folder the last file dialog happened to visit."""
        from larmor.desktop.paths import suggest_save_dir
        return suggest_save_dir(getattr(self, "source_path", None),
                                self._last_dir())

    def _resolve_open_path(self, path: str):
        """If `path` is a bare EXPNO with several processed datasets (pdata/N),
        ask the user which one. Returns the chosen pdata/N path, the original
        path, or None if the user cancelled."""
        from larmor.io import bruker

        p = Path(path)
        try:
            if not bruker.is_expno(p):
                return path
        except Exception:
            return path
        pdata = p / "pdata"
        if not pdata.is_dir():
            return path
        procs = sorted((d for d in pdata.iterdir()
                        if d.is_dir() and d.name.isdigit()
                        and ((d / "1r").exists() or (d / "2rr").exists())),
                       key=lambda d: int(d.name))
        if len(procs) <= 1:
            return path

        def _label(d: Path) -> str:
            kind = "2D" if (d / "2rr").exists() else "1D"
            title = ""
            t = d / "title"
            if t.exists():
                try:
                    title = t.read_text(errors="replace").splitlines()[0].strip()
                except Exception:
                    pass
            return f"proc {d.name}  ({kind})" + (f"  ·  {title}" if title else "")

        from PySide6.QtWidgets import QInputDialog
        items = [_label(d) for d in procs]
        choice, ok = QInputDialog.getItem(
            self, "Choose processed data",
            f"'{p.name}' has {len(procs)} processed datasets — which one?",
            items, 0, False)
        if not ok:
            return None
        return str(procs[items.index(choice)])

    def _offer_fxml_no_data(self, path: str):
        """A dmfit fit with a model but no embedded spectrum: let the user apply
        the model to the current data, pick data in the Explorer, or cancel."""
        model = self._recipe_model_from(path)
        if model is None or not model.get("sites"):
            self.statusBar().showMessage("that fit has no lines to apply")
            return
        box = QMessageBox(self)
        box.setWindowTitle("dmfit fit — no embedded data")
        box.setIcon(QMessageBox.Question)
        box.setText(f"“{Path(path).name}” contains a model ({len(model['sites'])} "
                    "line(s)) but no experimental data.")
        box.setInformativeText("Apply its lines to the spectrum on screen, or pick "
                               "experimental data in the Explorer to fit it on?")
        has_current = self.recipe is not None and self.exp_ppm is not None \
            and len(self.exp_ppm) > 0
        b_cur = box.addButton("Apply to current data", QMessageBox.AcceptRole)
        b_pick = box.addButton("Pick data in Explorer…", QMessageBox.ActionRole)
        box.addButton(QMessageBox.Cancel)
        b_cur.setEnabled(bool(has_current))
        box.exec()
        clicked = box.clickedButton()
        if clicked is b_cur and has_current:
            self.snapshot()
            self.recipe["sites"] = model["sites"]
            self.on_structure_changed()
            self.statusBar().showMessage(
                f"applied {len(model['sites'])} line(s) from {Path(path).name} "
                "— now Fit")
        elif clicked is b_pick:
            self._pending_apply_model = model
            self.explorer_dock.show(); self.explorer_dock.raise_()
            self.statusBar().showMessage(
                "pick the experimental data in the Explorer to fit this model on")

    def load_source(self, path: str, keep_fit: bool | None = None):
        """Open ANY dataset with a basic display, then let the user choose a
        process. 1D spectra go to the fit workbench; 2D datasets show a contour
        map; raw fid/ser get a quick preview. Nothing is ever rejected."""
        if str(path).lower().endswith(".larproj.json"):
            # a project bundle: File ▸ Open…, drag-drop and Open recent all
            # land here
            self._open_project_path(str(path))
            return
        path = self._resolve_open_path(path)
        if path is None:                      # user cancelled the proc chooser
            self.statusBar().showMessage("open cancelled")
            return
        self.statusBar().showMessage("loading…")
        QApplication.processEvents()
        self._sync_active()               # persist the outgoing document
        self._in_load_source = True
        try:
            self._load_source_body(path, keep_fit)
        finally:
            self._in_load_source = False
        if Path(path).exists():
            self._add_recent(path)

    def _load_source_body(self, path: str, keep_fit: bool | None):
        # First try the 1D-fittable path (dmfit / recipe / 1D processed).
        try:
            ppm, amp, recipe, meta, warnings = _load_any(path)
        except ValueError as exc:
            # a dmfit fit that carries a MODEL but no embedded data → let the user
            # attach it to a spectrum instead of just failing
            if (Path(path).suffix.lower() in (".fxml", ".fxmla")
                    and "experimental data" in str(exc).lower()):
                self._offer_fxml_no_data(path)
                return
            # not a 1D-fittable source: it may be a 2D or a raw fid/ser.
            handled = self._load_nonfittable(path)
            if handled:
                return
            # fall through to report the original error
            try:
                _load_any(path)
            except Exception as exc:
                QMessageBox.warning(self, "Load failed", humanize_error(exc))
                self.statusBar().showMessage("load failed")
            return
        except Exception as exc:
            QMessageBox.warning(self, "Load failed", humanize_error(exc))
            self.statusBar().showMessage("load failed")
            return

        self.central_stack.setCurrentWidget(self.view)   # 1D workbench

        # dmfit behaviour: if a fit is already on screen and the new source
        # brings no fit of its own, offer to keep the current lines
        existing = self.recipe.get("sites") if self.recipe else None
        incoming = recipe.get("sites")
        if keep_fit is None and existing and not incoming:
            text, default_yes = keep_fit_prompt(self.recipe.get("nucleus"),
                                                recipe.get("nucleus"))
            btn = QMessageBox.question(
                self, "Keep fit parameters?", text,
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes if default_yes else QMessageBox.No)
            keep_fit = btn == QMessageBox.Yes
        if keep_fit and existing and not incoming:
            recipe["sites"] = existing
            # carry the experiment parameters the new source knows
            for k in ("nucleus", "larmor_frequency_MHz", "spin_rate_Hz"):
                recipe[k] = recipe.get(k) or self.recipe.get(k)

        self.source_path = path
        self._retarget_watch()
        self.exp_ppm, self.exp_amp = ppm, amp
        self._proc_base = None
        if (Path(path).suffix.lower() == ".json" and recipe.get("processing")
                and not recipe.get("processing_from_raw")):
            # a reopened recipe arrives ALREADY replayed; the live pipeline
            # must re-apply from the unprocessed source, else the first panel
            # touch (or FID ⇄ spectrum) compounds the recorded chain on top
            # of its own result (p0 40 became 80)
            try:
                b_ppm, b_amp, *_ = _load_any(path, replay=False)
                self._proc_base = (np.asarray(b_ppm, float),
                                   np.asarray(b_amp, float))
            except Exception:
                pass
        self.recipe = recipe
        self.hidden.clear()
        self.undo_stack.clear()
        self.redo_stack.clear()
        QSettings("LARMOR", "app").setValue("lastDir", str(Path(path).parent))
        self.setWindowTitle(f"LARMOR — {Path(path).name}")
        self.view.set_experiment(ppm, amp)
        if not self.recipe["sites"]:
            # drop any model curve left over from the previous spectrum
            self.view.set_model(None, None, None, None, self.hidden)
        sample = recipe.get("sample") or Path(path).name
        self.view.set_title(sample)
        self._last_model = None
        self._first_sim = True
        self.zoom_full()
        if recipe["sites"]:
            self.zoom_sites()
        self._update_sn()
        msg = meta + ("   ⚠ " + " • ".join(warnings) if warnings else "")
        self.statusBar().showMessage(msg)
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles()
        self._update_exp_label()
        self._sync_zones()
        self._update_enabled()
        if self.recipe["sites"]:
            self.request_simulation()
        self._refresh_overlays()
        self._register_ws("1d")
        self._persist_session()
        self._warm_kernel()             # background pre-build (B6)

    def _load_nonfittable(self, path: str) -> bool:
        """Display a 2D dataset or a raw fid/ser without rejecting it.
        Returns True if it handled the path."""
        from larmor.io import bruker

        try:
            ref = bruker.resolve(path)
            data = bruker.read(path)
        except Exception:
            return False

        self.source_path = path
        QSettings("LARMOR", "app").setValue("lastDir", str(Path(path).parent))
        self.setWindowTitle(f"LARMOR — {Path(path).name}")

        if data.ndim == 2 and data.domain == "freq":
            self._show_2d(_nmrdata_to_data2d(data), Path(path).name,
                          "2D spectrum")
            kind = "arrayed/relaxation" if data.is_pseudo2d else "MQMAS 2D"
            self.statusBar().showMessage(
                f"{kind} displayed — Tools ▸ 2D MQMAS to fit · Tools ▸ "
                "Relaxation for a series · or send a 1D trace to fitting below")
            return True

        if data.ndim == 2 and data.domain == "time":
            from larmor import fourier

            d2 = fourier.ft2d_from_nmrdata(data, fourier.FT2DParams(
                f2_ops=[{"op": "fcor", "factor": 0.5}, {"op": "em", "lb_hz": 100}]))
            self._show_2d(d2, Path(path).name, "raw 2D (quick preview)")
            self.statusBar().showMessage(
                "raw 2D preview — File ▸ Open FID for full processing · "
                "Tools ▸ Relaxation for a T1/T2 series")
            return True

        if data.ndim == 1 and data.domain == "time":
            # quick magnitude-FT preview so the user sees a spectrum at once
            from larmor import fourier

            ppm, spec = fourier.ft1d(
                data.data, data.axes[0].sw_Hz, data.meta["larmor_MHz"],
                ops=[{"op": "fcor", "factor": 0.5}, {"op": "em", "lb_hz": 100}])
            order = np.argsort(ppm)
            spec = np.abs(spec)              # magnitude preview (phase-free)
            self._display_1d(ppm[order], spec[order],
                             data.meta.get("nucleus", ""),
                             data.meta["larmor_MHz"],
                             data.meta.get("spin_rate_Hz")
                             or data.meta.get("masr_Hz"),
                             Path(path).name + " (FID preview)", str(ref.expno))
            # the exact chain the preview applied, recorded so the first
            # FID ⇄ spectrum syncs the panel to raw / EM 100 / magnitude and
            # shows the TRUE fid; unticking 'magnitude' and phasing turns the
            # preview into a real workbench
            self.recipe["processing"] = [
                {"op": "fcor", "factor": 0.5}, {"op": "em", "lb_hz": 100},
                {"op": "ft", "offset_ppm": 0.0}, {"op": "magnitude"}]
            self.recipe["processing_from_raw"] = True
            self.statusBar().showMessage(
                "raw FID preview (magnitude, EM 100 Hz) — FID ⇄ spectrum "
                "(Ctrl+T) to re-apodize, untick magnitude and phase in the "
                "Processing panel, or File ▸ Open FID for the full dialog")
            return True
        return False

    # ------------------------------------------------------------------- S/N
    def _update_sn(self):
        """Signal-to-noise: peak signal over the RMS of a signal-free region.
        The noise region is taken from the quiet outer edges of the spectrum
        (robust to where the peak sits), matching TopSpin's sino spirit."""
        # this method runs at every active-1D-document change (see the comment
        # below), which is exactly when the previous fit's verdict and
        # covariance stop describing what is on screen: drop them here, once,
        # rather than at every call site (Parameter correlations used to show
        # another spectrum's covariance after a load or a workspace switch)
        self._health_reset()
        y = self.exp_amp
        if y is None or y.size < 20:
            self.lines_table.set_sn("")
            self._maybe_offer_sidebands()   # a tiny document hides an old offer
            return
        y = np.asarray(y, float)
        edge = max(5, y.size // 10)
        noise_region = np.concatenate([y[:edge], y[-edge:]])
        noise = float(np.std(noise_region - np.median(noise_region)))
        signal = float(np.max(np.abs(y - np.median(noise_region))))
        sn = signal / noise if noise > 0 else 0.0
        self.lines_table.set_sn(f"S/N {sn:,.0f}" if sn else "")
        # _update_sn runs at every active-1D-document change (load, workspace
        # switch, background subtraction, make-active) — exactly when the
        # literature-range overlay must follow the nucleus, so refresh it
        # here rather than duplicating the call at all five call sites.
        # The spinning-sideband offer follows the same rule: one hook here
        # covers load, workspace switch, re-processing, FID processing,
        # WURST division and background subtraction.
        self._update_ref_ranges()
        self._maybe_offer_sidebands()

    def _explorer_open(self, path: str):
        """Route an Explorer activation: a pending model to attach, an HMQC
        projection pick, else a normal load."""
        if getattr(self, "_pending_apply_model", None) is not None:
            model = self._pending_apply_model
            self._pending_apply_model = None
            self.load_source(path)
            if self.recipe is not None and model.get("sites"):
                self.snapshot()
                self.recipe["sites"] = model["sites"]
                self.on_structure_changed()
                self.statusBar().showMessage(
                    f"applied {len(model['sites'])} line(s) to "
                    f"{Path(path).name} — now Fit")
            return
        axis = self._proj_pick_axis
        if axis is None:
            self.load_source(path)
            return
        self._proj_pick_axis = None
        try:
            ppm, amp = self._read_1d(path)
        except Exception as exc:
            QMessageBox.warning(self, "Projection 1D", f"cannot read: {exc}")
            return
        self.view2d.set_projection_1d(axis, ppm, amp,
                                      color=self.PROJ_COLOR[axis], source=path)
        self.explorer.highlight(path, self.PROJ_COLOR[axis])
        self.statusBar().showMessage(
            f"{axis.upper()} 1D overlaid — adjust scale, then 'uncorrelated "
            f"{axis.upper()} →' for the non-correlated features")

    def _read_1d(self, path: str):
        """A 1D (ppm, amp) from any supported source."""
        try:
            ppm, amp, *_ = _load_any(path)
            return np.asarray(ppm), np.asarray(amp)
        except Exception:
            from larmor.io import bruker

            d = bruker.read(path)
            if d.ndim != 1 or d.domain != "freq":
                raise ValueError("not a 1D spectrum")
            return np.asarray(d.axes[0].values), np.asarray(d.data, float)

    def back_to_2d(self):
        """Switch to the most recent 2D workspace (the map you came from)."""
        if self.central_stack.currentWidget() is self.view2d:
            return
        for i in range(len(self.workspaces) - 1, -1, -1):
            if self.workspaces[i]["kind"] == "2d":
                self.switch_workspace(i)
                return
        self.statusBar().showMessage("no 2D workspace is open")

    def _show_2d(self, data2d, title: str, kind: str):
        from larmor.recipe import Recipe

        if not self._in_load_source:
            self._sync_active()
        self._data2d = data2d
        # let the 2D view seed its export dialogs from the dataset it shows
        self.view2d.source_path = getattr(self, "source_path", "") or ""
        self.view2d.set_data(data2d, f"{kind} — {title}")
        self.view2d.clear_model()
        self.central_stack.setCurrentWidget(self.view2d)
        # _show_2d does not pass through _update_sn; the outgoing 1D
        # workspace's verdict was saved by _sync_active() above
        self._health_reset()
        # a genuine spectroscopic 2D (not a relaxation array) is fittable, so
        # give it a recipe if there isn't a fit already in progress
        pseudo = any("pseudo" in n or "arrayed" in n for n in data2d.notes)
        self._data2d_fittable = not pseudo
        if not pseudo and (self.recipe is None or not self.recipe.get("sites")):
            self.recipe = Recipe(
                sample=title, source_kind="bruker", source_path=self.source_path,
                nucleus=data2d.nucleus,
                larmor_frequency_MHz=data2d.larmor_MHz).to_dict()
            self.hidden.clear()
            self.lines_table.rebuild(self.recipe, self.hidden)
            self._update_exp_label(); self._update_enabled()
        self._register_ws("2d")

    def _warm_kernel(self):
        """Start (at most one) background kernel pre-build for the current
        dataset -- see KernelWarmWorker."""
        if os.environ.get("LARMOR_NO_KERNEL_WARM"):
            return                      # tests / benchmarking opt-out
        rec = self.recipe or {}
        nucleus = rec.get("nucleus", "")
        lar = float(rec.get("larmor_frequency_MHz", 0.0) or 0.0)
        if not nucleus or lar <= 0 or not len(self.exp_ppm):
            return
        key = (nucleus, round(lar, 3),
               round(float(rec.get("spin_rate_Hz", 0.0) or 0.0)),
               len(self.exp_ppm))
        if getattr(self, "_warmed_key", None) == key:
            return
        w = getattr(self, "_warm_worker", None)
        if w is not None and w.isRunning():
            return                      # one at a time; next load re-checks
        self._warmed_key = key
        self._warm_worker = KernelWarmWorker(
            nucleus, lar, rec.get("spin_rate_Hz", 0.0) or 0.0, self.exp_ppm)
        self._warm_worker.start()

    def _display_1d(self, ppm, amp, nucleus, larmor, masr, title, expno):
        """Put a bare 1D spectrum on the workbench (no fit), ready to fit."""
        from larmor.recipe import Recipe

        if not self._in_load_source:
            self._sync_active()
        self.central_stack.setCurrentWidget(self.view)
        self.exp_ppm, self.exp_amp = np.asarray(ppm), np.asarray(amp)
        self._proc_base = None
        self.recipe = Recipe(sample=title, source_kind="bruker",
                             source_path=expno, nucleus=nucleus,
                             larmor_frequency_MHz=larmor,
                             spin_rate_Hz=masr or 0.0).to_dict()
        self.hidden.clear(); self.undo_stack.clear(); self.redo_stack.clear()
        self.view.set_experiment(self.exp_ppm, self.exp_amp)
        self.view.set_model(None, None, None, None, self.hidden)  # clear old model
        self.view.set_title(title)
        self._last_model = None; self._first_sim = True
        self.zoom_full()
        self.lines_table.rebuild(self.recipe, self.hidden)
        self._update_paddles(); self._update_exp_label(); self._update_enabled()
        self._refresh_overlays(); self._update_sn()
        self._register_ws("1d")
        self._warm_kernel()             # background pre-build (B6)

    def _trace_to_workbench(self, ppm, amp, label):
        """A 1D trace pulled out of the 2D view becomes a NEW workspace, so the
        2D map stays open in its own workspace."""
        nuc = self.view2d.data.nucleus if self.view2d.data else ""
        larmor = self.view2d.data.larmor_MHz if self.view2d.data else 0.0
        self._ws_mode = "new"
        self._display_1d(ppm, amp, nuc, larmor, None, label,
                         self.source_path or "")
        self.statusBar().showMessage(f"{label} — new workspace; add lines and Fit")

    # ------------------------------------------------------------- recipe io
    def save_recipe(self):
        if not self.recipe:
            return
        default = (self.recipe.get("sample") or "fit").strip()
        default = "".join(c if c.isalnum() or c in "-_" else "_"
                          for c in default)[:40] or "fit"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save recipe", str(Path(self._suggest_dir()) /
                                     f"{default}.recipe.json"),
            "LARMOR recipe (*.json)")
        if not path:
            return
        target = Path(path)
        if not self._guard_write(target):
            return
        Recipe.from_dict(self.recipe).save(target)
        self._add_recent_recipe(str(target))
        self.statusBar().showMessage(f"recipe saved — {target.name}")
        self.statusBar().showMessage(f"recipe saved: {target}")

    def _guard_write(self, target: Path) -> bool:
        """Saving is allowed anywhere — including next to raw data, which is
        where dmfit keeps fits and where the Explorer lists them. The one
        thing a save may never do is REPLACE an acquired file itself."""
        from larmor.desktop.paths import is_instrument_file
        if is_instrument_file(target):
            QMessageBox.warning(
                self, "Refused",
                f"{target.name} is an instrument file (the measurement "
                "itself) — pick another file name. Saving new files next to "
                "the data is fine.")
            return False
        return True

    def save_fit_as(self):
        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage("nothing to export — add and fit lines first")
            return
        from larmor.io import export

        default = (self.recipe.get("sample") or "fit").strip()
        default = "".join(c if c.isalnum() or c in "-_" else "_"
                          for c in default)[:40] or "fit"
        filters = ";;".join(
            f"{name} (*.{ext})" for name, (ext, _) in export.FORMATS.items())
        path, chosen = QFileDialog.getSaveFileName(
            self, "Save fit as", str(Path(self._suggest_dir()) / default), filters)
        if not path:
            return
        # figure out the format from the chosen filter (or the extension)
        fmt = None
        for name, (ext, fn) in export.FORMATS.items():
            if name == chosen or path.lower().endswith("." + ext):
                fmt = (ext, fn)
                break
        if fmt is None:
            fmt = ("json", export.export_json)
        ext, fn = fmt
        target = Path(path)
        if target.suffix.lower() != "." + ext:
            target = target.with_suffix("." + ext)
        if not self._guard_write(target):
            return
        recipe = Recipe.from_dict(self.recipe)
        try:
            if ext in ("txt",):
                fn(recipe, self.exp_ppm, self.exp_amp, target)
            elif ext == "fxmla":
                fn(recipe, self.exp_ppm, self.exp_amp, target)
            else:                        # csv, json take (recipe, path)
                fn(recipe, target)
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
            return
        self.statusBar().showMessage(f"fit exported: {target}")

    def open_figure_dialog(self):
        """The Figure studio (JSON editor). A spec that changed from what the
        dialog opened with is kept in the Workspaces dock (_remember_figure)."""
        from larmor.desktop.figure_dialog import FigureDialog

        dlg = FigureDialog(self, self.source_path, self.recipe)
        seed = dlg.spec.toPlainText()
        dlg.exec()
        text = dlg.spec.toPlainText()
        if text == seed:
            return
        try:
            spec = json.loads(text)
        except ValueError:
            return                      # not valid JSON: nothing to keep
        self._remember_figure(spec)

    def save_spectrum(self):
        from larmor.io import spectra

        if not self.exp_ppm.size:
            self.statusBar().showMessage("no spectrum to save")
            return
        default = (self.recipe.get("sample") if self.recipe else "") or "spectrum"
        default = "".join(c if c.isalnum() or c in "-_ " else "_"
                          for c in default).strip()[:60] or "spectrum"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save spectrum", str(Path(self._suggest_dir()) / f"{default}.csv"),
            "LARMOR spectrum (*.csv);;Text (*.txt)")
        if not path:
            return
        meta = {}
        if self.recipe:
            meta = {"nucleus": self.recipe.get("nucleus", ""),
                    "larmor_MHz": self.recipe.get("larmor_frequency_MHz", 0.0),
                    "spin_rate_Hz": self.recipe.get("spin_rate_Hz", 0.0),
                    "sample": self.recipe.get("sample", "")}
        spectra.write_csv(path, self.exp_ppm, self.exp_amp, meta)
        self.statusBar().showMessage(f"spectrum saved to {Path(path).name} "
                                     "(reopen it with File ▸ Open…)")
