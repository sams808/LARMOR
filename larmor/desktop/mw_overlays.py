"""Main-window mixin: the Datasets-dock overlay cockpit and the 1D-on-2D
projection overlays (``PROJ_COLOR`` is the per-axis colour table).

Owned state: the ``_overlays`` list (created by the datasets dock builder)
and ``_proj_pick_axis`` (which 2D axis a picked 1D file is meant for).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QFileDialog, QMessageBox

from larmor.desktop.mw_files import _load_any


class _OverlaysMixin:
    """Reference / comparison overlays on the 1D and 2D views."""

    # --------------------------------------------------------- overlays
    def add_overlay_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Add a spectrum to compare", self._last_dir(),
            "Spectra (*.fxmla *.json 1r 2rr *.txt *.csv);;All files (*)")
        if not path:
            return
        try:
            label, ppm, amp, info = self._read_overlay_source(path)
        except Exception as exc:
            QMessageBox.warning(self, "Add overlay", f"cannot read: {exc}")
            return
        self._add_overlay(label, ppm, amp, path, info)

    def add_overlay_path(self, path: str) -> bool:
        """Overlay one more spectrum on the active one -- Shift + drop, the
        Explorer's right-click, File > Overlay a spectrum... The active
        spectrum and its fit stay exactly as they are; the new curve is
        drawn behind them in its own colour and listed in the Datasets
        dock (colour, scale, shift, offset, match height, remove, make
        active)."""
        try:
            label, ppm, amp, info = self._read_overlay_source(path)
        except Exception as exc:                          # noqa: BLE001
            self.statusBar().showMessage(
                f"cannot overlay {Path(path).name}: {exc}")
            return False
        self._add_overlay(label, ppm, amp, path, info)
        act = getattr(self, "actOverlaysVisible", None)
        if act is not None and not act.isChecked():
            act.setChecked(True)                          # an added overlay is meant to be seen
            self.view.set_overlays_hidden(False)          # setChecked emits no triggered()
        self.statusBar().showMessage(
            f"overlaid {label} -- {len(self._overlays)} compared spectrum(s); "
            "Datasets dock: colour, scale, shift, offset, match height, remove, "
            "make active")
        return True

    def add_overlay_paths(self, paths):
        """Shift + drop of several files: every one becomes an overlay."""
        n = sum(1 for p in paths if self.add_overlay_path(p))
        if len(paths) > 1:
            self.statusBar().showMessage(
                f"overlaid {n} of {len(paths)} dropped spectra -- "
                f"{len(self._overlays)} compared spectrum(s)")
        return n

    def clear_overlays(self):
        n = len(self._overlays)
        self._overlays.clear()
        self._refresh_overlays()
        self.statusBar().showMessage(
            f"removed {n} overlay(s)" if n else "no overlays to remove")

    def _toggle_overlays_visible(self, on: bool):
        """View > Overlays (Ctrl+Shift+V): hide or show the compared spectra
        without losing them."""
        self.view.set_overlays_hidden(not on)
        if self._overlays:
            self.statusBar().showMessage(
                f"{len(self._overlays)} overlay(s) " + ("shown" if on else "hidden"))

    def compare_overlays(self):
        """Datasets ▸ Compare acquisition…: acqus / procs / auditp of the
        active spectrum and every overlay with a source, in the shared
        comparability table (larmor.comparability)."""
        from larmor import comparability
        from larmor.desktop.comparability_dialog import ComparabilityDialog

        paths, labels = self._overlay_sources()
        cmp = comparability.compare(
            [comparability.read_params(p) if p else None for p in paths], labels)
        if cmp.level == "none":
            self.statusBar().showMessage(
                "no Bruker acquisition files among these spectra — nothing to compare")
            return
        ComparabilityDialog(self, cmp).exec()

    def _overlay_sources(self):
        """(paths, labels) of the active spectrum then every overlay, in the
        Datasets dock's order (index 0 = active; '' when it has no source)."""
        active = ""
        if getattr(self, "recipe", None):
            active = self.recipe.get("sample") or ""
        src = getattr(self, "source_path", "") or ""
        if not active and src:
            active = Path(src).name
        paths, labels = [src], [active or "active"]
        for ov in self._overlays:
            paths.append(ov.get("source") or "")
            labels.append(str(ov.get("label", "")))
        return paths, labels

    def _overlay_comparison(self):
        """The comparison of active + overlays, recomputed only when the
        tuple of sources changes (a colour or offset change must not re-read a
        dozen JCAMP files); None when it cannot be built."""
        from larmor import comparability

        paths, labels = self._overlay_sources()
        srcs = tuple(paths)
        if srcs != getattr(self, "_cmp_srcs", None):
            self._cmp_srcs = srcs
            try:
                self._cmp_ov = comparability.compare(
                    [comparability.read_params(p) if p else None for p in paths], labels)
            except Exception:  # noqa: BLE001 -- a read failure never breaks overlays
                self._cmp_ov = None
        return getattr(self, "_cmp_ov", None)

    @staticmethod
    def _read_overlay_source(path: str):
        """(label, ppm, amp) for any of the overlay-file-picker's supported
        sources -- shared by add_overlay_dialog and open_project (restoring
        saved overlays), so both go through the exact same, exactly-once-
        fixed loading logic."""
        try:
            # load_any's real return order is (ppm, amp, recipe, meta, warnings)
            # -- this used to be unpacked as (recipe, ppm, amp, ...) here,
            # silently mis-assigning `recipe` to the ppm array; recipe.get(...)
            # then always raised AttributeError, so EVERY overlay source
            # _load_any actually supports (.recipe.json, .fxmla,
            # .csv/.txt/.dat) fell through to the Bruker-only fallback below
            # and failed there too -- only a raw Bruker 1r/2rr path ever
            # worked, "by accident", via that fallback. Found while wiring
            # overlay restoration into save_project/open_project.
            ppm, amp, recipe, *_ = _load_any(path)
            label = recipe.get("sample") or Path(path).name
            info = {"nucleus": recipe.get("nucleus", ""),
                    "larmor_MHz": float(recipe.get("larmor_frequency_MHz", 0.0)
                                        or 0.0),
                    "npts": int(np.asarray(ppm).size), "title": ""}
            return label, np.asarray(ppm), np.asarray(amp), info
        except Exception:
            from larmor.io import bruker

            d = bruker.read(path)
            if d.ndim != 1 or d.domain != "freq":
                raise ValueError("not a 1D spectrum") from None
            meta = d.meta or {}
            # "1r" says nothing -- name the overlay by sample folder + EXPNO
            parts = Path(path).parts
            sample = parts[-5] if len(parts) >= 5 else ""      # .../S/EXP/pdata/N/1r
            expno = str(meta.get("expno", "") or "")
            label = " · ".join(x for x in (sample, expno) if x) or Path(path).name
            info = {"nucleus": meta.get("nucleus", ""),
                    "larmor_MHz": float(meta.get("larmor_MHz", 0.0) or 0.0),
                    "npts": int(np.asarray(d.data).size),
                    "title": (meta.get("title", "") or "").splitlines()[0]
                    if meta.get("title") else ""}
            return (label, np.asarray(d.axes[0].values),
                    np.asarray(d.data, float), info)

    def _add_overlay(self, label, ppm, amp, source="", info=None):
        """Append a compared spectrum. Besides its arrays the dict carries
        the display transform -- ``scale`` (x1), ``shift`` (ppm) and
        ``yoff`` (a fraction of the active span) -- applied at draw time
        only; with 'match height' ticked the scale starts at the matching
        factor."""
        from larmor import display
        from larmor.desktop.datasets import overlay_color

        ov = {"label": label, "ppm": np.asarray(ppm), "amp": np.asarray(amp),
              "color": overlay_color(len(self._overlays)), "visible": True,
              "source": source, **display.OVERLAY_DEFAULTS, **(info or {})}
        match = getattr(self.datasets_panel, "match", None)
        if match is not None and match.isChecked():
            ov["scale"] = display.match_scale(self.exp_amp, ov["amp"])
        self._overlays.append(ov)
        self._refresh_overlays()
        self.datasets_dock.raise_()

    def overlay_set_color(self, i: int, color: str):
        if 0 <= i < len(self._overlays) and color:
            self._overlays[i]["color"] = color
            self._refresh_overlays()

    def overlay_remove(self, i: int):
        if 0 <= i < len(self._overlays):
            del self._overlays[i]
            self._refresh_overlays()

    def overlay_visibility(self, i: int, on: bool):
        if 0 <= i < len(self._overlays):
            self._overlays[i]["visible"] = on
            self._refresh_overlays()

    def overlay_set_scale(self, i: int, value: float):
        """Datasets row: the display scale (x). A hand-typed factor means
        'match height' no longer describes the scales, so the box unticks."""
        if 0 <= i < len(self._overlays):
            self._overlays[i]["scale"] = float(value)
            self._uncheck_match()
            self._refresh_overlays()

    def overlay_set_shift(self, i: int, ppm: float):
        """Datasets row: shift the overlay along the axis (ppm), e.g. to
        line up a reference peak. Display only."""
        if 0 <= i < len(self._overlays):
            self._overlays[i]["shift"] = float(ppm)
            self._refresh_overlays()

    def overlay_set_yoff(self, i: int, frac: float):
        """Datasets row: raise the overlay by a fraction of the active
        spectrum's span, on top of the global stack offset."""
        if 0 <= i < len(self._overlays):
            self._overlays[i]["yoff"] = float(frac)
            self._refresh_overlays()

    def overlay_reset(self, i: int):
        """Right-click > Reset: drawn as stored again (x1, no shift, no
        offset)."""
        from larmor import display

        if 0 <= i < len(self._overlays):
            self._overlays[i].update(display.OVERLAY_DEFAULTS)
            self._uncheck_match()
            self._refresh_overlays()

    def overlay_match_height(self, on: bool):
        """Datasets > match height: ticked, every overlay's scale becomes
        the factor that brings its maximum to the active spectrum's --
        written into its row, editable afterwards; unticked, back to x1.
        Display only: the stored arrays and every export are untouched."""
        from larmor import display

        for ov in self._overlays:
            ov["scale"] = (display.match_scale(self.exp_amp, ov["amp"])
                           if on else 1.0)
        self._refresh_overlays()

    def _uncheck_match(self):
        m = getattr(getattr(self, "datasets_panel", None), "match", None)
        if m is not None and m.isChecked():
            m.blockSignals(True)
            m.setChecked(False)
            m.blockSignals(False)

    def overlay_make_active(self, i: int):
        if not (0 <= i < len(self._overlays)):
            return
        ov = self._overlays[i]
        if not ov.get("source"):
            return
        prev = None
        if self.exp_ppm.size and self.source_path:
            prev = (self.recipe.get("sample") if self.recipe else None,
                    self.exp_ppm.copy(), self.exp_amp.copy(), self.source_path)
        del self._overlays[i]
        self._ws_mode = "reuse"           # swap the active spectrum in place
        self.load_source(ov["source"])
        if prev and prev[3] != self.source_path:
            self._add_overlay(prev[0] or Path(prev[3]).name, prev[1], prev[2],
                              prev[3])
        self._refresh_overlays()

    def _refresh_overlays(self):
        if not hasattr(self, "datasets_panel"):
            return
        from larmor import display

        # every overlay is drawn through larmor.display.overlay_display --
        # normalised by its OWN trace under View > Y axis, then its scale,
        # shift and offsets -- the stored arrays untouched; the offsets are
        # fractions of the active spectrum's DISPLAYED span
        span = 1.0
        if self.exp_amp.size:
            span = float(np.nanmax(self.exp_amp) - np.nanmin(self.exp_amp)) or 1.0
        span *= self.view.y_scale()
        step = self.datasets_panel.offset.value()
        mode, region = self.view.y_mode()
        drawn = []
        for k, ov in enumerate(self._overlays):
            if ov.get("visible", True):
                x, y = display.overlay_display(
                    ov["ppm"], ov["amp"], scale=ov.get("scale", 1.0),
                    shift=ov.get("shift", 0.0), yoff=ov.get("yoff", 0.0),
                    stack=step * (k + 1), span=span, mode=mode, region=region)
                drawn.append((x, y, ov["color"], ov["label"]))
        self.view.set_overlays(drawn)
        label, detail = "", ""
        if self.recipe is not None:
            label = self.recipe.get("sample") or (
                Path(self.source_path).name if self.source_path else "")
            nuc = self.recipe.get("nucleus", "")
            mhz = float(self.recipe.get("larmor_frequency_MHz", 0.0) or 0.0)
            bits = [b for b in (nuc, f"{mhz:.1f} MHz" if mhz else "",
                                f"{self.exp_ppm.size} pts"
                                if self.exp_ppm.size else "") if b]
            detail = " · ".join(bits)
        # the comparability marker: an overlay acquired / processed unlike
        # the compared set's majority gets '⚠ LB 100 Hz' on its detail line
        cmp = self._overlay_comparison()
        for i, ov in enumerate(self._overlays):
            ov["comparability"] = cmp.signature(i + 1) if cmp is not None else ""
        if cmp is not None and cmp.signature(0):
            detail = (detail + " · " if detail else "") + "⚠ " + cmp.signature(0)
        self.datasets_panel.rebuild(label, self._overlays, detail)

    #: colour assigned to each HMQC projection axis (matches the overlay + the
    #: Explorer highlight)
    PROJ_COLOR = {"f2": "#e8832a", "f1": "#6a4fb0"}

    def overlay_1d_on_2d(self, axis: str, source: str):
        """Superpose a 1D spectrum on the F2 (direct) or F1 (indirect) projection
        of the 2D map — the easy path: the current 1D or a file, no Explorer."""
        if source == "current":
            if self.exp_ppm is None or not np.asarray(self.exp_ppm).size:
                self.statusBar().showMessage(
                    "no current 1D spectrum — open one first, or use "
                    "Overlay 1D ▸ From file…")
                return
            ppm, amp = np.asarray(self.exp_ppm), np.asarray(self.exp_amp)
            src = self._current_1d_source()
        else:                                    # from a file
            path, _ = QFileDialog.getOpenFileName(
                self, f"1D spectrum to overlay on {axis.upper()}", "",
                "Spectra (*.fxmla *.json 1r *.txt *.csv);;All files (*)")
            if not path:
                path = QFileDialog.getExistingDirectory(
                    self, "…or a 1D EXPNO/pdata folder")
            if not path:
                return
            try:
                ppm, amp = self._read_1d(path)
            except Exception as exc:
                QMessageBox.warning(self, "Overlay 1D", f"cannot read: {exc}")
                return
            src = path
        # source= lets a project bundle save the overlay by reference (an
        # empty source is dropped at save, like a sourceless 1D overlay)
        self.view2d.set_projection_1d(axis, ppm, amp,
                                      color=self.PROJ_COLOR[axis], source=src)
        self.statusBar().showMessage(
            f"1D superposed on the {axis.upper()} projection — tune 'scale' / "
            "'fit scale' in the overlay bar")

    def _current_1d_source(self) -> str:
        """The file behind the 1D arrays held in exp_ppm / exp_amp. While a
        2D map is displayed, source_path is the MAP's path and the 1D arrays
        belong to the last 1D document, so that document's snapshot is looked
        up by content; '' when nothing on disk matches."""
        if self.central_stack.currentWidget() is self.view:
            return self.source_path or ""
        x, y = np.asarray(self.exp_ppm), np.asarray(self.exp_amp)
        for ws in reversed(self.workspaces):
            snap = ws.get("snap") or {}
            if ws.get("kind") != "1d" or not snap.get("source_path"):
                continue
            if (np.array_equal(np.asarray(snap.get("exp_ppm")), x)
                    and np.array_equal(np.asarray(snap.get("exp_amp")), y)):
                return str(snap["source_path"])
        return ""

    def load_projection_1d(self, axis: str):
        """HMQC: arm a pick — the next spectrum clicked in the Explorer (or via
        its Browse…) is overlaid on the F2/F1 projection and highlighted."""
        self._proj_pick_axis = axis
        self.explorer_dock.show(); self.explorer_dock.raise_()
        self.statusBar().showMessage(
            f"click a spectrum in the Explorer for the {axis.upper()} projection "
            "(use Browse… to reach one elsewhere)")
