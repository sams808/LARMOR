"""Main-window mixin: the co-fit page (two datasets, tied parameters).

The central-stack page built by ``_build_cofit_page`` and its whole state
machine: the tie table, per-dataset tables and views, dataset loading,
simulation, the co-fit run and applying the result back to the main recipe.

Owned state: every ``cofit_*`` / ``_cofit_*`` attribute (page, views,
tables, recipes, ties, timers).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import (QApplication, QCheckBox, QFileDialog,
                               QHBoxLayout, QLabel, QMessageBox, QPushButton,
                               QSplitter, QVBoxLayout, QWidget)

from larmor.desktop.mw_files import _load_any
from larmor.desktop.plot import SpectrumView
from larmor.desktop.table import LinesTable
from larmor.desktop.workers import _fit_tol, humanize_error
from larmor.recipe import Recipe


class _CofitMixin:
    """The co-fit split page and its state machine."""

    # ------------------------------------------------------------- co-fit page
    def _build_cofit_page(self):
        """The co-fit central page: 1D (left) and 2D (right) side by side, each
        with its OWN parameter table under it. Parameters are per-dataset so they
        can be decorrelated; a tie bar chooses which ones are shared in the fit."""
        from larmor.desktop.twod_view import Contour2DView

        page = QWidget(); v = QVBoxLayout(page); v.setContentsMargins(2, 2, 2, 2)
        v.setSpacing(3)
        # -- control bar --
        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("<b>Co-fit</b>"))
        self.btnCofitAdd = QPushButton("＋ Add / replace dataset…")
        self.btnCofitAdd.clicked.connect(self._cofit_add_dataset)
        self.btnCofitPrev = QPushButton("Preview")
        self.btnCofitPrev.clicked.connect(lambda: self._cofit_simulate())
        self.btnCofitRun = QPushButton("Run co-fit")
        self.btnCofitRun.clicked.connect(self.run_cofit_fit)
        self.btnCofitClose = QPushButton("Close co-fit")
        self.btnCofitClose.clicked.connect(self.close_cofit)
        for b in (self.btnCofitAdd, self.btnCofitPrev, self.btnCofitRun,
                  self.btnCofitClose):
            ctl.addWidget(b)
        ctl.addStretch(1)
        self.cofit_rmsd = QLabel("")
        self.cofit_rmsd.setStyleSheet("font-weight:600;")   # theme palette colour
        ctl.addWidget(self.cofit_rmsd)
        v.addLayout(ctl)
        # -- tie bar (rebuilt from the model's actual parameters) --
        self.cofit_tie_bar = QWidget()
        self._cofit_tie_lay = QHBoxLayout(self.cofit_tie_bar)
        self._cofit_tie_lay.setContentsMargins(2, 0, 2, 0)
        self._cofit_tie = {}                         # param -> QCheckBox
        v.addWidget(self.cofit_tie_bar)

        # -- two plot columns (1D | 2D); the per-dataset parameter tables live in
        #    the bottom "Fit parameters" dock, split the same way (below) --
        self.cofit_view1d = SpectrumView()
        self.cofit_view2d = Contour2DView()
        self.cofit_table1d = LinesTable()
        self.cofit_table2d = LinesTable()
        for tbl, which in ((self.cofit_table1d, 1), (self.cofit_table2d, 2)):
            tbl.edited.connect(lambda w=which: self._cofit_on_edit(w))
            tbl.constraint_edited.connect(lambda w=which: self._cofit_on_edit(w))
            tbl.structure.connect(lambda r, a, w=which: self._cofit_struct(w, r, a))
            # a multi-selection Remove / Delete: highest row first, so the
            # lower indices stay valid, through the same per-row path
            tbl.remove_lines.connect(
                lambda rows, w=which: [self._cofit_struct(w, r, "remove")
                                       for r in sorted(rows, reverse=True)])
            # sidebands act on the workbench recipe, not on a co-fit copy
            tbl.sidebands_requested.connect(lambda _i: self.statusBar().showMessage(
                "co-fit tables keep both recipes identical — add sidebands on "
                "the workbench fit, then open the co-fit again"))
            tbl.compute.connect(lambda: self._cofit_simulate())   # footer Compute
            tbl.fit.connect(self.run_cofit_fit)                   # footer Fit

        for vw in (self.cofit_view1d, self.cofit_view2d):
            vw.setMinimumWidth(300)
        split = QSplitter(Qt.Horizontal)            # 1D | 2D plots, side by side
        split.addWidget(self.cofit_view1d)
        split.addWidget(self.cofit_view2d)
        split.setChildrenCollapsible(False)
        split.setStretchFactor(0, 1); split.setStretchFactor(1, 1)
        self.cofit_split = split
        v.addWidget(split, 1)

        # the parameter tables, split 1D | 2D to mirror the plots; swapped into
        # the bottom dock while co-fitting (see _build_bottom_docks / _set_cofit_dock)
        def _tblcol(table, header):
            w = QWidget(); cv = QVBoxLayout(w); cv.setContentsMargins(0, 0, 0, 0)
            cv.setSpacing(2)
            lab = QLabel(header); lab.setStyleSheet("font-weight:600;")
            cv.addWidget(lab); cv.addWidget(table, 1)
            w.setMinimumWidth(280)
            return w
        self.cofit_tables = QSplitter(Qt.Horizontal)
        self.cofit_tables.addWidget(_tblcol(self.cofit_table1d, "1D parameters"))
        self.cofit_tables.addWidget(_tblcol(self.cofit_table2d, "2D parameters"))
        self.cofit_tables.setChildrenCollapsible(False)
        self.cofit_tables.setStretchFactor(0, 1)
        self.cofit_tables.setStretchFactor(1, 1)
        self.cofit_page = page
        self.central_stack.addWidget(page)           # index 2
        self._cofit = None
        self._cofit_last = None                      # stash for quick re-entry
        self.central_stack.currentChanged.connect(self._on_central_changed)
        self._cofit_timer = QTimer(self); self._cofit_timer.setSingleShot(True)
        self._cofit_timer.setInterval(200)
        self._cofit_timer.timeout.connect(self._cofit_simulate_now)

    # short, model-aware labels for the tie bar
    _COFIT_LABEL = {
        "isotropic_chemical_shift_ppm": "δiso", "sigma_Cq_MHz": "σCq",
        "shift_fwhm_ppm": "dCS", "line_fwhm_ppm": "line", "Cq_MHz": "Cq",
        "eta": "η", "eta_q": "ηq", "eps": "eps", "zeta_ppm": "ζ",
        "eta_cs": "ηcs", "gl": "G/L", "gauss_fwhm_ppm": "G",
        "lorentz_fwhm_ppm": "L", "czjzek_d": "d",
        "shift_slope_ppm_per_MHz": "dδ/dC_Q", "split_ppm": "Δδ",
        "pop_a": "p(A)", "k_ex_hz": "k_ex",
    }

    def _cofit_tieable(self) -> list:
        """Parameters that actually influence the current model's lineshape (the
        union of the sites' params, minus amplitude) — so the tie bar never
        offers parameters the selected lineshape doesn't have."""
        st = self._cofit or {}
        seen, order = set(), []
        for site in (st.get("r1") or {}).get("sites", []):
            for p in site.get("params", {}):
                if p != "amplitude" and p not in seen:
                    seen.add(p); order.append(p)
        return order

    def _cofit_tie_rebuild(self):
        lay = self._cofit_tie_lay
        while lay.count():
            w = lay.takeAt(0).widget()
            if w is not None:
                w.setParent(None)
        self._cofit_tie = {}
        lay.addWidget(QLabel("Tie 1D↔2D:"))
        tie = (self._cofit or {}).get("tie", set())
        for p in self._cofit_tieable():
            cb = QCheckBox(self._COFIT_LABEL.get(p, p))
            cb.setChecked(p in tie)
            cb.toggled.connect(lambda on, name=p: self._cofit_tie_toggled(name, on))
            self._cofit_tie[p] = cb; lay.addWidget(cb)
        lay.addStretch(1)

    def _cofit_tie_toggled(self, param: str, on: bool):
        st = self._cofit
        if not st:
            return
        tie = st.setdefault("tie", set())
        if on:
            tie.add(param)
            self._cofit_copy_param(st["r1"], st["r2"], param)   # sync 2D <- 1D
            self._cofit_rebuild_tables()
        else:
            tie.discard(param)
        self._cofit_simulate()

    @staticmethod
    def _cofit_copy_param(src, dst, param):
        for i, s in enumerate(src.get("sites", [])):
            if i < len(dst.get("sites", [])) and param in s.get("params", {}) \
                    and param in dst["sites"][i].get("params", {}):
                dst["sites"][i]["params"][param]["value"] = \
                    s["params"][param]["value"]

    def _cofit_home(self):
        """A stable identity for the workspace a co-fit was launched from, so a
        stash is only resumed when the user returns to the *same* dataset — not
        after they open something else (which would otherwise restore the stale
        co-fit over the new data)."""
        if self.active_ws is not None and 0 <= self.active_ws < len(self.workspaces):
            return id(self.workspaces[self.active_ws])
        return ("path", self.source_path)

    def open_cofit(self):
        if not self.recipe or not self.recipe.get("sites"):
            self.statusBar().showMessage(
                "set up the fit (add lines) on one dataset first")
            return
        self._sync_active()
        home = self._cofit_home()
        prev = self._cofit_last
        if (prev and prev.get("home") == home
                and (prev.get("d1") or prev.get("d2"))):
            st = prev                                 # resume the paused co-fit
        else:
            st = {"d1": None, "d2": None, "home": home}
            base = json.loads(json.dumps(self.recipe))
            st["r1"] = json.loads(json.dumps(base))   # independent per-dataset
            st["r2"] = json.loads(json.dumps(base))   # copies of the model
            st["tie"] = set(self._default_tie(base))
            if self.central_stack.currentWidget() is self.view2d and self._data2d:
                st["d2"] = (self._data2d, Path(self.source_path or "2D").name)
            elif self.exp_ppm is not None and np.asarray(self.exp_ppm).size:
                st["d1"] = (np.asarray(self.exp_ppm), np.asarray(self.exp_amp),
                            self.recipe.get("sample") or "current")
        self._cofit = st
        self.central_stack.setCurrentWidget(self.cofit_page)
        self._set_cofit_dock(True)
        self._cofit_tie_rebuild()
        self._cofit_rebuild_tables()
        self._cofit_split_even()
        self._cofit_refresh_panels()
        if st.get("d1") and st.get("d2"):
            self.statusBar().showMessage(
                "co-fit: both datasets loaded — untie the parameters that differ, "
                "then Preview / Run co-fit")
            self._cofit_simulate()
            return
        need = "a 2D MQMAS map" if st["d2"] is None else "a 1D MAS spectrum"
        self.statusBar().showMessage(
            f"co-fit: model from the current fit — add {need} (＋ Add dataset)")
        self._cofit_add_dataset()

    @staticmethod
    def _default_tie(recipe) -> list:
        from larmor.multifit import DEFAULT_SHARE
        present = {p for s in recipe.get("sites", []) for p in s.get("params", {})}
        return [p for p in DEFAULT_SHARE if p in present]

    def _on_central_changed(self, *_):
        """Leaving the co-fit page (loading data, switching workspace, …) exits
        co-fit mode. The datasets are stashed so Fit ▸ Co-fit datasets resumes
        without re-picking the files."""
        if getattr(self, "_cofit", None) is None:
            return
        if self.central_stack.currentWidget() is self.cofit_page:
            return
        self._set_cofit_dock(False)
        self._cofit_last = self._cofit
        self._cofit = None
        self.statusBar().showMessage(
            "co-fit paused (another dataset is showing) — Fit ▸ "
            "Co-fit datasets returns to it with both datasets still loaded")

    def _cofit_split_even(self):
        """Give the 1D and 2D columns half the page each. Qt ignores setSizes on
        a splitter not yet laid out, so re-apply on the next event-loop pass."""
        def apply():
            w = self.cofit_split.width() or self.cofit_page.width() or 1000
            self.cofit_split.setSizes([w // 2, w - w // 2])
        apply()
        QTimer.singleShot(0, apply)

    def _set_cofit_dock(self, on: bool):
        """While co-fitting, the bottom "Fit parameters" dock shows the split
        1D | 2D parameter tables (mirroring the plot split); otherwise it shows
        the normal single Fit-Parameters table."""
        self.lines_stack.setCurrentWidget(
            self.cofit_tables if on else self.lines_table)
        self.lines_dock.setWindowTitle(
            "Co-fit parameters" if on else "Fit parameters")
        if on:
            self.lines_dock.show()
            self.lines_dock.raise_()

            def even():
                w = self.cofit_tables.width() or self.lines_dock.width() or 1000
                self.cofit_tables.setSizes([w // 2, w - w // 2])
            even()
            QTimer.singleShot(0, even)

    def _cofit_rebuild_tables(self):
        st = self._cofit or {}
        if st.get("r1"):
            self.cofit_table1d.rebuild(st["r1"], self.hidden)
        if st.get("r2"):
            self.cofit_table2d.rebuild(st["r2"], self.hidden)

    def _cofit_on_edit(self, which: int):
        """A value/constraint changed in one table. Push every TIED parameter to
        the other recipe (last edit wins for tied params) and re-simulate."""
        st = self._cofit
        if not st:
            return
        src, dst = (st["r1"], st["r2"]) if which == 1 else (st["r2"], st["r1"])
        for p in st.get("tie", set()):
            self._cofit_copy_param(src, dst, p)
        other = self.cofit_table2d if which == 1 else self.cofit_table1d
        other.rebuild(dst, self.hidden)
        self._cofit_simulate()

    def _cofit_struct(self, which: int, row: int, action: str):
        """Structure edits (remove/duplicate/visibility) keep BOTH recipes
        identical — co-fit needs the same sites/models in each."""
        st = self._cofit
        if not st:
            return
        if action == "visibility":
            (self.hidden.discard(row) if row in self.hidden
             else self.hidden.add(row))
        else:
            for r in (st["r1"], st["r2"]):
                sites = r.get("sites", [])
                if action == "remove" and row < len(sites):
                    sites.pop(row)
                elif action == "duplicate" and row < len(sites):
                    copy = json.loads(json.dumps(sites[row]))
                    copy["label"] = (copy.get("label") or "line") + "-copy"
                    sites.append(copy)
            self.hidden.clear()
        self._cofit_rebuild_tables()
        self._cofit_simulate()

    def _cofit_add_dataset(self):
        if self._cofit is None:
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Add a dataset to co-fit", self._last_dir(),
            "Spectra / maps (*.fxmla *.json 1r 2rr *.txt *.csv);;All files (*)")
        if not path:
            path = QFileDialog.getExistingDirectory(self, "…or a 2D EXPNO/pdata")
        if not path:
            return
        try:
            kind, payload = self._cofit_load(path)
        except Exception as exc:
            QMessageBox.warning(self, "Add dataset", f"cannot read: {exc}")
            return
        self._cofit[kind] = payload
        self._cofit_refresh_panels()
        self._cofit_simulate()

    def _cofit_load(self, path: str):
        from larmor import twod
        from larmor.io import bruker

        try:
            ref = bruker.resolve(path)
            if ref.ndim == 2 and ref.target == "2rr":
                d = twod.read_bruker_2d(path)
                if d.z.ndim == 2:
                    return "d2", (d, Path(path).name)
        except Exception:
            pass
        ppm, amp, *_ = _load_any(path)
        return "d1", (np.asarray(ppm), np.asarray(amp), Path(path).name)

    def _cofit_refresh_panels(self):
        st = self._cofit or {}
        if st.get("d1"):
            ppm, amp, label = st["d1"]
            self.cofit_view1d.set_experiment(ppm, amp)
            self.cofit_view1d.set_title(f"1D — {label}")
        if st.get("d2"):
            d, label = st["d2"]
            self.cofit_view2d.set_data(d, f"2D — {label}")

    def _cofit_active(self) -> bool:
        return (self._cofit is not None
                and self.central_stack.currentWidget() is self.cofit_page)

    def _cofit_simulate(self):
        if self._cofit_active():
            self._cofit_timer.start()

    def _cofit_simulate_now(self):
        """Live overlay of EACH dataset's own model on its panel (no fit)."""
        st = self._cofit
        if not st:
            return
        if st.get("d1") and st.get("r1", {}).get("sites"):
            from larmor import engine
            ppm, amp, _ = st["d1"]
            r1 = st["r1"]
            labels = [s.get("label") or s["model"] for s in r1["sites"]]
            rec = Recipe.from_dict(json.loads(json.dumps(r1)))
            try:
                x, total, per = engine.simulate(rec, exp_ppm=ppm)
                self.cofit_view1d.set_model(x, total, per, labels, self.hidden,
                                            ppm, amp)
            except Exception as exc:
                self.statusBar().showMessage(f"co-fit 1D sim: {exc}")
        if st.get("d2") and st.get("r2", {}).get("sites"):
            self._cofit_sim_2d(st["d2"][0], st["r2"])

    def _cofit_sim_2d(self, d2, recipe_dict):
        from larmor import twod
        rec = Recipe.from_dict(json.loads(json.dumps(recipe_dict)))
        d = d2.normalized()
        kernel = twod._kernel_for(rec, d)
        total, per = twod.simulate_2d(rec, kernel)
        if getattr(rec, "mqmas_f1_ref_vary", True):     # auto-align β for preview
            from scipy.interpolate import RegularGridInterpolator
            itp = RegularGridInterpolator((d.f1_ppm, d.f2_ppm), d.z,
                                          bounds_error=False, fill_value=0.0)
            G1, G2 = np.meshgrid(kernel.f1_ppm, kernel.f2_ppm, indexing="ij")
            mf = total.ravel(); mn = np.sqrt((mf * mf).sum()) or 1.0
            best = (-1.0, 0.0)
            for b in np.linspace(-40, 40, 81):
                ev = itp(np.stack([(G1 + b).ravel(), G2.ravel()], -1)).ravel()
                den = np.sqrt((ev * ev).sum()) * mn
                cc = float((ev * mf).sum()) / den if den > 0 else 0.0
                if cc > best[0]:
                    best = (cc, float(b))
            rec.mqmas_f1_ref_ppm = best[1]
        per_disp = [p if i not in self.hidden else np.zeros_like(p)
                    for i, p in enumerate(per)]
        self.cofit_view2d.set_model(np.sum(per_disp, axis=0), kernel.f2_ppm,
                                    twod.mqmas_f1_axis(kernel, rec),
                                    per_site=per_disp)

    def run_cofit_fit(self):
        from larmor.multifit import fit_cofit

        st = self._cofit
        if not st or not st.get("d1") or not st.get("d2"):
            self.statusBar().showMessage("co-fit needs both a 1D and a 2D dataset")
            return
        self.snapshot()
        share = tuple(sorted(st.get("tie", set())))
        entries = [(Recipe.from_dict(json.loads(json.dumps(st["r1"]))),
                    (st["d1"][0], st["d1"][1])),
                   (Recipe.from_dict(json.loads(json.dumps(st["r2"]))),
                    st["d2"][0])]
        self.btnCofitRun.setEnabled(False)
        self._progress_start("co-fit (building MQMAS kernel…)")

        def _cb(params, it, resid, *a, **k):
            self._prog_label = "co-fit"
            try:
                rms = float(np.sqrt(np.mean(np.asarray(resid, float) ** 2)))
            except Exception:
                rms = float("nan")
            self._progress_tick(it if isinstance(it, int) else 0, rms)
            QApplication.processEvents()
        try:
            result = fit_cofit(entries, share=share, iter_cb=_cb, tol=_fit_tol())
        except Exception as exc:
            self._progress_end(False)
            self.btnCofitRun.setEnabled(True)
            QMessageBox.warning(self, "Co-fit", humanize_error(f"co-fit failed: {exc}"))
            return
        self._progress_end(True)
        self.btnCofitRun.setEnabled(True)
        st["r1"] = result.recipes[0].to_dict()
        st["r2"] = result.recipes[1].to_dict()
        self._cofit_rebuild_tables()
        self._cofit_simulate_now()
        self.cofit_rmsd.setText(
            "RMSD  " + " · ".join(f"{k} {r:.4f}" for k, r in
                                  zip(("1D", "2D"), result.rmsd)))
        self.statusBar().showMessage(
            "co-fit done — tied parameters shared, untied ones fit independently")

    def _cofit_apply_to_main(self):
        """Adopt the 1D recipe as the main fit when closing (so the workbench
        keeps the shared model)."""
        st = self._cofit or self._cofit_last
        if st and st.get("r1", {}).get("sites"):
            self.recipe = json.loads(json.dumps(st["r1"]))

    def close_cofit(self):
        self._cofit_apply_to_main()
        self._set_cofit_dock(False)
        self._cofit = None
        self._cofit_last = None                # deliberate close: do not resume
        back = self.view2d if (self._data2d is not None
                               and getattr(self, "_data2d_fittable", False)) \
            else self.view
        self.central_stack.setCurrentWidget(back)
        if self.recipe and self.recipe.get("sites"):
            self.lines_table.rebuild(self.recipe, self.hidden)
        self.statusBar().showMessage("co-fit closed")
