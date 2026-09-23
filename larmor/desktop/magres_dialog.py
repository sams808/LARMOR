"""Tools > Import DFT tensors (.magres): computed shielding / EFG tensors
become starting sites -- spin-aware, one site per crystallographic position,
shieldings converted through a FITTED calibration line.

Two tabs. **Sites**: open the magres, pick the isotope (pre-selected to the
open spectrum's nucleus, its spin shown) and a model the tensor can seed for
that spin (a spin-1/2 nucleus never sees a quadrupolar model), choose the
shielding -> shift conversion (a calibration line δ = a·σ + b, or the single
σ_ref with slope −1), and read every grouped site with its predicted shift
and its uncertainty before adding. **Calibration**: build that line from
reference compounds -- their magres files and the shifts measured in LARMOR
(a saved fit, the current fit or a typed value) -- with residuals, the
(a, b) covariance, the degrees-of-freedom warning and a check that every
reference was computed with the same [calculation] settings.

Every number comes from larmor.dft / larmor.shiftcal (Qt-free, tested); this
module renders and collects choices. File inputs all have a programmatic
entry (``_load``, ``cal.add_reference_magres``, ``cal.delta_from_recipe``,
``cal.load``) so tests never open a file dialog; refusals go to labels, never
to message boxes; the remembered folders go through paths.remember_dir.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QHBoxLayout,
    QLabel, QPushButton, QRadioButton, QTableWidget, QTableWidgetItem,
    QTabWidget, QVBoxLayout, QWidget,
)

from larmor.desktop import theme

#: settings keys for the remembered folders (paths.remember_dir, gated)
MAGRES_DIR_KEY = "magresDir"
SHIFTCAL_DIR_KEY = "shiftcalDir"

_AMBER = "#c88a1e"
_RED = "#d62728"

#: the sign statement every user must be able to read off the screen
SIGN_NOTE = ("ζ is the shielding anisotropy as csa_mas stores it (δaniso = −ζ); "
             "with slope a the seeded anisotropy is −a·ζ_σ.   δ = a·σ + b.")
ZETA_TIP = ("ζ is the shielding anisotropy handed to mrsimulator; δaniso = −ζ "
            "(same sign as Tools ▸ Herzfeld–Berger)")
CQ_TIP_HALF = ("spin-1/2: no quadrupolar interaction; EFG records, if any, "
               "are ignored")
NO_CONVERSION = ("set σ_ref or load/fit a calibration line — computed "
                 "shieldings are not chemical shifts")
TWO_POINT = ("2 points: zero residual degrees of freedom — the line is exact "
             "by construction; the ± reflects the reference-shift errors only.")


def _spin_txt(spin: float | None) -> str:
    if spin is None:
        return "?"
    return str(Fraction(spin).limit_denominator(2))


def _ro(text: str, tip: str = "") -> QTableWidgetItem:
    it = QTableWidgetItem(text)
    it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
    if tip:
        it.setToolTip(tip)
    return it


def _fmt_pm(value: float | None, err: float | None, nd: int = 1) -> str:
    if value is None:
        return "—"
    s = f"{value:.{nd}f}"
    return s if err is None else f"{s} ± {err:.{nd}f}"


@dataclass
class MagresImport:
    """What the dialog hands back: recipe-ready site dicts and the record."""
    sites: list[dict]
    note: str
    provenance: dict
    model: str
    n_groups: int
    ratios: list[int]
    conversion: str            # short text for the status bar
    lock_amplitude: bool
    share_width: bool
    area_amplitude: bool       # amplitude is an area (gl_norm)


# ======================================================================
class _CalibrationTab(QWidget):
    """Reference compounds -> the fitted line δ = a·σ + b."""

    #: emitted with the ShiftCalibration when 'Use this line' is pressed
    use_line = Signal(object)

    def __init__(self, parent, isotope: str = "", recipe: dict | None = None):
        super().__init__(parent)
        self.isotope = isotope or ""
        self.recipe = recipe
        self.points: list = []
        self.calibration = None
        self._filling = False
        self._recipe_sites: list[tuple[str, float, float | None]] = []
        self._recipe_name = ""
        t = theme.active()
        v = QVBoxLayout(self)

        intro = QLabel(
            "A GIPAW shielding σ sits on a scale set by the calculation; the "
            "measured shift δ is relative to a reference. Fit δ = a·σ + b over "
            "reference compounds measured under the SAME referencing as the "
            "sample (a ≠ −1 in practice). Two references make the line exact "
            "by construction — three or more give residuals that mean "
            "something.")
        intro.setWordWrap(True)
        intro.setStyleSheet(f"color:{t.text_dim};")
        v.addWidget(intro)

        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(
            ["reference", "magres", "σ_iso (ppm)", "δ_exp (ppm)", "± (ppm)",
             "source", "residual (ppm)"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setToolTip("δ_exp and ± are editable: double-click to type "
                              "a value (a row without ± makes the whole fit "
                              "unweighted)")
        self.table.itemChanged.connect(self._on_item_changed)
        v.addWidget(self.table, 1)

        row = QHBoxLayout()
        b_add = QPushButton("Add reference (.magres)…")
        b_add.clicked.connect(self._pick_magres)
        b_rec = QPushButton("δ from a saved fit…")
        b_rec.setToolTip("a .recipe.json: the largest-amplitude site's δ ± "
                         "stderr is proposed for the selected row")
        b_rec.clicked.connect(self._pick_recipe)
        self.b_cur = QPushButton("δ from the current fit")
        self.b_cur.setEnabled(bool(recipe and recipe.get("sites")))
        self.b_cur.clicked.connect(lambda: self.delta_from_current(self.recipe))
        self.site_pick = QComboBox()
        self.site_pick.setToolTip("which site of that fit is the reference line")
        self.site_pick.setVisible(False)
        self.site_pick.currentIndexChanged.connect(self._site_picked)
        b_del = QPushButton("Remove")
        b_del.clicked.connect(self.remove_selected)
        for b in (b_add, b_rec, self.b_cur, self.site_pick, b_del):
            row.addWidget(b)
        row.addStretch(1)
        v.addLayout(row)

        self.plot = pg.PlotWidget(background=t.plot_bg)
        self.plot.setLabel("bottom", "σ_iso (computed)", units="ppm")
        self.plot.setLabel("left", "δ_exp (measured)", units="ppm")
        for ax in ("bottom", "left"):
            self.plot.getPlotItem().getAxis(ax).enableAutoSIPrefix(False)
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        self.plot.setMinimumHeight(150)
        v.addWidget(self.plot, 1)

        self.res = QLabel("")
        self.res.setWordWrap(True)
        self.res.setTextInteractionFlags(Qt.TextSelectableByMouse)
        v.addWidget(self.res)
        self.mismatch_lbl = QLabel("")
        self.mismatch_lbl.setWordWrap(True)
        self.mismatch_lbl.setStyleSheet(f"color:{_RED};")
        self.mismatch_lbl.setVisible(False)
        v.addWidget(self.mismatch_lbl)
        self.chk_mismatch = QCheckBox(
            "Fit despite differing [calculation] settings (recorded as a flag)")
        self.chk_mismatch.setVisible(False)
        v.addWidget(self.chk_mismatch)

        bottom = QHBoxLayout()
        self.btnFit = QPushButton("Fit line")
        self.btnFit.clicked.connect(self.fit)
        b_save = QPushButton("Save calibration…")
        b_save.clicked.connect(self._pick_save)
        b_load = QPushButton("Load…")
        b_load.clicked.connect(self._pick_load)
        self.btnUse = QPushButton("Use this line")
        self.btnUse.setEnabled(False)
        self.btnUse.clicked.connect(self._use)
        bottom.addWidget(self.btnFit); bottom.addWidget(b_save)
        bottom.addWidget(b_load); bottom.addStretch(1)
        bottom.addWidget(self.btnUse)
        v.addLayout(bottom)

    # ------------------------------------------------------------ points
    def add_point(self, pt) -> None:
        self.points.append(pt)
        self._invalidate()
        self._rebuild_table()

    def add_reference_magres(self, path: str, isotope: str | None = None) -> list:
        """One row per equivalent group of ``isotope`` in the file (δ_exp
        left for the user to fill). Returns the new points."""
        from larmor import dft
        from larmor.shiftcal import CalibrationPoint

        mf = dft.read_magres_file(path)
        dft.assign_isotopes(mf.sites)
        isos = sorted({s.isotope for s in mf.sites if s.isotope})
        iso = isotope or self.isotope
        if not iso or iso not in isos:
            iso = self.isotope if self.isotope in isos else (isos[0] if isos else "")
        if not iso:
            raise ValueError(f"{mf.name}: no isotope could be assigned")
        self.isotope = self.isotope or iso
        groups = dft.group_equivalent(dft.sites_for_isotope(mf.sites, iso))
        base = mf.calc.prefix or Path(mf.path).name.split(".")[0]
        new = []
        for g in groups:
            sh = g.shielding()
            if not sh:
                continue
            name = base if len(groups) == 1 else f"{base} {g.label}"
            new.append(CalibrationPoint(
                name=name, sigma_iso=float(sh["iso_ppm"]), magres=mf.path,
                n_atoms=int(g.multiplicity), calc=mf.calc))
        for pt in new:
            self.points.append(pt)
        self._invalidate()
        self._rebuild_table()
        if new:
            self.table.selectRow(self.table.rowCount() - 1)
        return new

    def _target_row(self) -> int:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if rows:
            return rows[0]
        return len(self.points) - 1

    def _apply_delta(self, row: int, delta: float, err: float | None,
                     source: str) -> None:
        if not (0 <= row < len(self.points)):
            self.res.setText("add a reference row first")
            self.res.setStyleSheet(f"color:{_AMBER};")
            return
        p = self.points[row]
        p.delta_exp = float(delta)
        p.delta_err = None if err is None else float(err)
        p.source = source
        self._invalidate()
        self._rebuild_table()
        self.table.selectRow(row)

    @staticmethod
    def _sites_of(recipe: dict) -> list[tuple[str, float, float | None, float]]:
        out = []
        for i, s in enumerate(recipe.get("sites", [])):
            p = s.get("params", {}).get("isotropic_chemical_shift_ppm")
            if not p:
                continue
            a = s.get("params", {}).get("amplitude", {})
            out.append((s.get("label") or f"s{i}", float(p.get("value", 0.0)),
                        p.get("stderr"), float(a.get("value", 0.0) or 0.0)))
        return out

    def _use_recipe(self, recipe: dict, name: str, site_index: int | None):
        sites = self._sites_of(recipe)
        if not sites:
            raise ValueError(f"{name}: no site with an isotropic shift")
        if site_index is None:
            site_index = int(np.argmax([s[3] for s in sites]))
        label, delta, err, _amp = sites[site_index]
        self._recipe_sites = [(s[0], s[1], s[2]) for s in sites]
        self._recipe_name = name
        self.site_pick.blockSignals(True)
        self.site_pick.clear()
        for lab, d, e in self._recipe_sites:
            self.site_pick.addItem(f"{lab}: {_fmt_pm(d, e, 2)} ppm")
        self.site_pick.setCurrentIndex(site_index)
        self.site_pick.blockSignals(False)
        self.site_pick.setVisible(len(sites) > 1)
        nuc = recipe.get("nucleus", "")
        sr = recipe.get("sr_hz")
        src = f"fit: {name} {label}" + (f" ({nuc}" if nuc else "") + \
            (f", SR {sr:.1f} Hz)" if (nuc and sr) else (")" if nuc else ""))
        self._apply_delta(self._target_row(), delta,
                          None if err is None else float(err), src)
        return delta, err

    def delta_from_recipe(self, path: str, site_index: int | None = None):
        """Propose the largest-amplitude site's δ ± stderr of a saved fit for
        the selected row (a combo then picks another site)."""
        from larmor.recipe import Recipe

        rec = Recipe.load(path).to_dict()
        return self._use_recipe(rec, Path(path).name, site_index)

    def delta_from_current(self, recipe: dict | None,
                           site_index: int | None = None):
        if not recipe or not recipe.get("sites"):
            self.res.setText("no fit is open")
            self.res.setStyleSheet(f"color:{_AMBER};")
            return None
        return self._use_recipe(recipe, "current fit", site_index)

    def _site_picked(self, idx: int):
        if 0 <= idx < len(self._recipe_sites):
            lab, d, e = self._recipe_sites[idx]
            self._apply_delta(self._target_row(), d, e,
                              f"fit: {self._recipe_name} {lab}")

    def remove_selected(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()},
                      reverse=True)
        for r in rows:
            if 0 <= r < len(self.points):
                self.points.pop(r)
        self._invalidate()
        self._rebuild_table()

    def _on_item_changed(self, item: QTableWidgetItem):
        if self._filling or item.column() not in (3, 4):
            return
        r = item.row()
        if not (0 <= r < len(self.points)):
            return
        txt = item.text().strip().replace("−", "-").replace(",", ".")
        p = self.points[r]
        try:
            val = float(txt) if txt else None
        except ValueError:
            self._rebuild_table()
            return
        if item.column() == 3:
            p.delta_exp = val
            if not p.source or p.source.startswith("typed"):
                p.source = "typed"
        else:
            p.delta_err = val if (val is None or val > 0) else None
        self._invalidate()
        self._rebuild_table()

    # ------------------------------------------------------------ table
    def _invalidate(self):
        self.calibration = None
        self.btnUse.setEnabled(False)

    def _rebuild_table(self):
        self._filling = True
        try:
            self.table.setRowCount(len(self.points))
            res = (self.calibration.residuals if self.calibration else [])
            fitted = [p for p in self.points if p.ready] if self.calibration else []
            for r, p in enumerate(self.points):
                self.table.setItem(r, 0, _ro(p.name))
                self.table.setItem(r, 1, _ro(Path(p.magres).name if p.magres
                                             else "—", p.magres))
                self.table.setItem(r, 2, _ro(f"{p.sigma_iso:.4f}",
                                             f"mean over {p.n_atoms} atom(s)"))
                d = QTableWidgetItem("" if p.delta_exp is None
                                     else f"{p.delta_exp:.3f}")
                e = QTableWidgetItem("" if p.delta_err is None
                                     else f"{p.delta_err:.3f}")
                for it in (d, e):
                    it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable
                                | Qt.ItemIsEditable)
                self.table.setItem(r, 3, d)
                self.table.setItem(r, 4, e)
                self.table.setItem(r, 5, _ro(p.source or ("—" if p.delta_exp is None
                                                          else "typed")))
                rtxt = ""
                if p in fitted:
                    rtxt = f"{res[fitted.index(p)]:+.3f}"
                self.table.setItem(r, 6, _ro(rtxt))
        finally:
            self._filling = False

    # ------------------------------------------------------------ fit
    def fit(self):
        from larmor import dft, shiftcal

        ready = [p for p in self.points if p.ready]
        iso = self.isotope or ""
        mism = shiftcal.mismatches(ready, dft.element_of(iso))
        self.mismatch_lbl.setVisible(bool(mism))
        self.chk_mismatch.setVisible(bool(mism))
        if mism:
            self.mismatch_lbl.setText(
                "references computed with different [calculation] settings: "
                + "; ".join(mism))
            if not self.chk_mismatch.isChecked():
                self.res.setText("the line was not fitted: the references "
                                 "above do not share one shielding scale")
                self.res.setStyleSheet(f"color:{_RED};")
                self._invalidate()
                self._plot()
                return None
        try:
            cal = shiftcal.fit_calibration(ready, iso, allow_mismatch=bool(mism))
        except ValueError as exc:
            self.res.setText(str(exc))
            self.res.setStyleSheet(f"color:{_RED};")
            self._invalidate()
            self._plot()
            return None
        self.calibration = cal
        self.btnUse.setEnabled(True)
        self._rebuild_table()
        txt = cal.describe()
        txt += "\nresiduals: " + ", ".join(
            f"{p.name} {r:+.2f}" for p, r in zip(cal.points, cal.residuals))
        colour = theme.active().accent
        if cal.dof == 0:
            txt += ("\n" + TWO_POINT + " Add references (cryolite, sulphohalite "
                    "and chiolite standards are measured in the 2025-12 / "
                    "2026-03 sessions).")
            colour = _AMBER
        for f in cal.flags:
            if f != shiftcal.FLAG_ZERO_DOF:
                txt += "\n" + f
                if f == shiftcal.FLAG_UNWEIGHTED or "no covariance" in f:
                    colour = _AMBER
        self.res.setText(txt)
        self.res.setStyleSheet(f"color:{colour}; font-weight:600;")
        self._plot()
        return cal

    def _plot(self):
        self.plot.clear()
        t = theme.active()
        pts = [p for p in self.points if p.ready]
        if not pts:
            return
        x = np.array([p.sigma_iso for p in pts])
        y = np.array([p.delta_exp for p in pts])
        e = np.array([p.delta_err if p.delta_err is not None else 0.0 for p in pts])
        if np.any(e > 0):
            self.plot.addItem(pg.ErrorBarItem(x=x, y=y, height=2 * e, beam=0.0,
                                              pen=pg.mkPen(t.text_dim)))
        self.plot.plot(x, y, pen=None, symbol="o", symbolSize=8,
                       symbolBrush=pg.mkBrush(t.accent), symbolPen=None)
        cal = self.calibration
        if cal is not None:
            span = max(float(x.max() - x.min()), 1.0)
            xx = np.linspace(x.min() - 0.1 * span, x.max() + 0.1 * span, 100)
            yy = cal.slope * xx + cal.intercept
            if cal.dof >= 1 and cal.cov is not None:
                ee = np.array([cal.delta_err(v) or 0.0 for v in xx])
                band = pg.FillBetweenItem(
                    pg.PlotDataItem(xx, yy - ee), pg.PlotDataItem(xx, yy + ee),
                    brush=pg.mkBrush(*theme._rgb(t.model), 40))
                self.plot.addItem(band)
            self.plot.plot(xx, yy, pen=pg.mkPen(t.model, width=1.6))
            for xi, yi, r in zip(x, y, cal.residuals):
                self.plot.plot([xi, xi], [yi - r, yi],
                               pen=pg.mkPen(_RED, width=1.2))

    # ------------------------------------------------------------ io
    def save(self, path: str) -> None:
        if self.calibration is None:
            self.fit()
        if self.calibration is None:
            return
        self.calibration.save(path)

    def load(self, path: str):
        from larmor.shiftcal import ShiftCalibration

        cal = ShiftCalibration.load(path)
        self.points = list(cal.points)
        if cal.isotope:
            self.isotope = cal.isotope
        self.calibration = cal
        self.btnUse.setEnabled(True)
        self._rebuild_table()
        self.res.setText(f"{Path(path).name}: {cal.describe()}"
                         + ("\n" + TWO_POINT if cal.dof == 0 else ""))
        self.res.setStyleSheet(
            f"color:{_AMBER if cal.dof == 0 else theme.active().accent}; "
            "font-weight:600;")
        self._plot()
        return cal

    def _use(self):
        if self.calibration is None:
            self.fit()
        if self.calibration is not None:
            self.use_line.emit(self.calibration)

    # ------------------------------------------------------------ pickers
    def _pick_magres(self):
        from larmor.desktop.paths import remember_dir, remembered_dir

        p, _ = QFileDialog.getOpenFileName(
            self, "reference compound: CASTEP/QE .magres",
            remembered_dir(MAGRES_DIR_KEY), "magres (*.magres);;All (*)")
        if not p:
            return
        remember_dir(MAGRES_DIR_KEY, p)
        try:
            self.add_reference_magres(p)
        except Exception as exc:  # noqa: BLE001
            self.res.setText(f"failed: {exc}")
            self.res.setStyleSheet(f"color:{_RED};")

    def _pick_recipe(self):
        p, _ = QFileDialog.getOpenFileName(
            self, "measured reference shift: LARMOR fit", "",
            "LARMOR recipe (*.json);;All (*)")
        if not p:
            return
        try:
            self.delta_from_recipe(p)
        except Exception as exc:  # noqa: BLE001
            self.res.setText(f"failed: {exc}")
            self.res.setStyleSheet(f"color:{_RED};")

    def _pick_save(self):
        from larmor.desktop.paths import remember_dir, remembered_dir

        if self.calibration is None and self.fit() is None:
            return
        name = f"{self.isotope or 'shift'}.shiftcal.json"
        p, _ = QFileDialog.getSaveFileName(
            self, "Save calibration line",
            str(Path(remembered_dir(SHIFTCAL_DIR_KEY) or ".") / name),
            "LARMOR shift calibration (*.shiftcal.json);;All (*)")
        if not p:
            return
        remember_dir(SHIFTCAL_DIR_KEY, p)
        self.save(p)

    def _pick_load(self):
        from larmor.desktop.paths import remember_dir, remembered_dir

        p, _ = QFileDialog.getOpenFileName(
            self, "Load calibration line", remembered_dir(SHIFTCAL_DIR_KEY),
            "LARMOR shift calibration (*.shiftcal.json *.json);;All (*)")
        if not p:
            return
        remember_dir(SHIFTCAL_DIR_KEY, p)
        try:
            self.load(p)
        except Exception as exc:  # noqa: BLE001
            self.res.setText(f"failed: {exc}")
            self.res.setStyleSheet(f"color:{_RED};")


# ======================================================================
class MagresDialog(QDialog):
    def __init__(self, parent, nucleus: str = "", recipe: dict | None = None,
                 spin_rate_Hz: float = 0.0, exp_max: float = 0.0):
        super().__init__(parent)
        self.setWindowTitle("Import DFT tensors (.magres)")
        self.resize(840, 580)
        self.nucleus = nucleus or ""
        self.recipe = recipe
        self.spin_rate_Hz = float(spin_rate_Hz or 0.0)
        self.exp_max = float(exp_max or 0.0)
        self.mf = None                       # dft.MagresFile
        self.calibration = None              # the fitted line in use
        self.result: MagresImport | None = None
        self._groups: list = []
        self._mismatch: list[str] = []
        self._file_notes: list[str] = []
        t = theme.active()

        v = QVBoxLayout(self)
        self.tabs = QTabWidget()
        v.addWidget(self.tabs, 1)
        self.sites_tab = QWidget()
        self.tabs.addTab(self.sites_tab, "Sites")
        self.cal = _CalibrationTab(self, isotope=self.nucleus, recipe=recipe)
        self.tabs.addTab(self.cal, "Calibration")
        self.cal.use_line.connect(self._use_line)
        self._build_sites_tab(t)

        bottom = QHBoxLayout()
        self.btnHelp = QPushButton("Help")
        self.btnHelp.clicked.connect(self._help)
        bottom.addWidget(self.btnHelp)
        bottom.addStretch(1)
        b_cancel = QPushButton("Cancel")
        b_cancel.clicked.connect(self.reject)
        bottom.addWidget(b_cancel)
        self.btnAdd = QPushButton("Add sites to the fit")
        self.btnAdd.setDefault(True)
        self.btnAdd.clicked.connect(self._accept)
        bottom.addWidget(self.btnAdd)
        v.addLayout(bottom)
        self._update_add_enabled()

    # ------------------------------------------------------------ build
    def _build_sites_tab(self, t):
        v = QVBoxLayout(self.sites_tab)

        top = QHBoxLayout()
        self.lbl = QLabel("no file")
        self.lbl.setStyleSheet("font-weight: 600;")
        btn = QPushButton("Open .magres…")
        btn.clicked.connect(self._pick)
        self.calc_lbl = QLabel("")
        self.calc_lbl.setStyleSheet(f"color: {t.text_dim};")
        self.calc_lbl.setToolTip("the [calculation] header: code, functional, "
                                 "cutoffs — what a calibration line must match")
        top.addWidget(self.lbl, 1)
        top.addWidget(self.calc_lbl)
        top.addWidget(btn)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("isotope"))
        self.iso = QComboBox()
        self.iso.setMinimumWidth(130)
        opts.addWidget(self.iso)
        opts.addWidget(QLabel("model"))
        self.model = QComboBox()
        self.model.setToolTip("only models a computed tensor can seed for this "
                              "spin; a spin-1/2 nucleus never sees a "
                              "quadrupolar model")
        opts.addWidget(self.model, 1)
        v.addLayout(opts)

        conv = QVBoxLayout()
        hdr = QLabel("Shift conversion")
        hdr.setStyleSheet("font-weight: 600;")
        conv.addWidget(hdr)
        r1 = QHBoxLayout()
        self.rb_cal = QRadioButton("Calibration line δ = a·σ + b")
        self.rb_cal.setChecked(True)
        self.cal_lbl = QLabel("none loaded")
        self.cal_lbl.setStyleSheet(f"color: {t.text_dim};")
        self.cal_lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.btn_cal_browse = QPushButton("Load…")
        self.btn_cal_browse.clicked.connect(self._browse_calibration)
        b_build = QPushButton("Build one → Calibration tab")
        b_build.clicked.connect(lambda: self.tabs.setCurrentWidget(self.cal))
        r1.addWidget(self.rb_cal); r1.addWidget(self.cal_lbl, 1)
        r1.addWidget(self.btn_cal_browse); r1.addWidget(b_build)
        conv.addLayout(r1)
        r2 = QHBoxLayout()
        self.rb_ref = QRadioButton("σ_ref only (δ = σ_ref − σ, slope −1)")
        self.ref = QDoubleSpinBox()
        self.ref.setRange(-1e4, 1e4)
        self.ref.setDecimals(2)
        self.ref.setSpecialValueText("—")       # the minimum reads as blank
        self.ref.setValue(self.ref.minimum())
        self.ref.setToolTip("absolute shielding of the reference compound "
                            "(ppm): type it — a blank value never converts")
        self.ref.setKeyboardTracking(False)
        r2.addWidget(self.rb_ref); r2.addWidget(self.ref); r2.addStretch(1)
        conv.addLayout(r2)
        v.addLayout(conv)

        grp = QHBoxLayout()
        self.chk_group = QCheckBox("Group symmetry-equivalent atoms")
        self.chk_group.setChecked(True)
        grp.addWidget(self.chk_group)
        grp.addWidget(QLabel("tolerance (ppm)"))
        self.tol = QDoubleSpinBox()
        self.tol.setRange(0.001, 50.0)
        self.tol.setDecimals(3)
        self.tol.setValue(0.05)
        self.tol.setKeyboardTracking(False)
        grp.addWidget(self.tol)
        self.chk_lock = QCheckBox("Lock relative populations to multiplicity")
        grp.addWidget(self.chk_lock)
        self.chk_width = QCheckBox("One shared linewidth")
        self.chk_width.setChecked(True)
        grp.addWidget(self.chk_width)
        grp.addStretch(1)
        v.addLayout(grp)

        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(
            ["site", "atoms", "σ_iso (ppm)", "δ_pred ± (ppm)",
             "ζ (ppm) / η_CS", "C_Q (MHz) / η_Q", "note"])
        self.table.horizontalHeaderItem(4).setToolTip(ZETA_TIP)
        self.table.horizontalHeaderItem(3).setToolTip(
            "δ = a·σ + b with ± from cov(a, b); for a 2-point line the ± "
            "reflects the reference-shift errors only")
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table, 1)

        self.note = QLabel(SIGN_NOTE)
        self.note.setWordWrap(True)
        self.note.setStyleSheet(f"color: {t.text_dim};")
        v.addWidget(self.note)
        self.warn = QLabel("")
        self.warn.setWordWrap(True)
        v.addWidget(self.warn)
        self.chk_override = QCheckBox("Apply anyway (recorded in provenance)")
        self.chk_override.setVisible(False)
        v.addWidget(self.chk_override)

        self.iso.currentIndexChanged.connect(self._fill_models)
        self.model.currentIndexChanged.connect(self._on_model_changed)
        self.rb_cal.toggled.connect(self._refresh)
        self.rb_ref.toggled.connect(self._refresh)
        self.ref.valueChanged.connect(self._ref_typed)
        self.chk_group.toggled.connect(self._refresh)
        self.tol.valueChanged.connect(self._refresh)
        self.chk_lock.toggled.connect(self._update_add_enabled)
        self.chk_width.toggled.connect(self._update_add_enabled)
        self.chk_override.toggled.connect(self._update_add_enabled)

    # ------------------------------------------------------------ file
    def _pick(self):
        from larmor.desktop.paths import remember_dir, remembered_dir

        p, _ = QFileDialog.getOpenFileName(
            self, "CASTEP/QE .magres", remembered_dir(MAGRES_DIR_KEY),
            "magres (*.magres);;All (*)")
        if not p:
            return
        remember_dir(MAGRES_DIR_KEY, p)
        self._load(p)

    def _load(self, path: str) -> None:
        from larmor import dft

        try:
            mf = dft.read_magres_file(path)
            warnings = dft.assign_isotopes(mf.sites)
        except Exception as exc:  # noqa: BLE001
            self.warn.setText(f"failed: {exc}")
            self.warn.setStyleSheet(f"color: {_RED};")
            return
        self.mf = mf
        self.lbl.setText(mf.name)
        self.lbl.setToolTip(mf.path)
        self.calc_lbl.setText(mf.calc.describe())
        isos = sorted({s.isotope for s in mf.sites if s.isotope})
        self._file_notes = list(warnings)
        if self.nucleus and self.nucleus not in isos:
            self._file_notes.append(
                f"the open spectrum is {self.nucleus}; this file has "
                + (", ".join(isos) or "no assignable isotope"))
        self.iso.blockSignals(True)
        self.iso.clear()
        for iso in isos:
            self.iso.addItem(f"{iso} · I = {_spin_txt(dft.spin_of(iso))}", iso)
        if self.nucleus in isos:
            self.iso.setCurrentIndex(isos.index(self.nucleus))
        self.iso.blockSignals(False)
        if not self.cal.isotope and self.iso.currentData():
            self.cal.isotope = self.iso.currentData()
        self._fill_models()

    def _spin(self) -> float | None:
        from larmor import dft

        iso = self.iso.currentData()
        return dft.spin_of(iso) if iso else None

    def _fill_models(self, *_):
        from larmor import dft

        spin = self._spin()
        names = dft.seedable_models(spin)
        default = dft.default_model(
            spin, self.spin_rate_Hz if self.recipe is not None else None)
        self.model.blockSignals(True)
        self.model.clear()
        for n in names:
            self.model.addItem(dft.model_choice_label(n), n)
        if default in names:
            self.model.setCurrentIndex(names.index(default))
        self.model.blockSignals(False)
        self._on_model_changed()

    def current_model(self) -> str:
        return self.model.currentData() or ""

    def _on_model_changed(self, *_):
        from larmor import dft
        from larmor.estimate import _WIDTH_KEY

        model = self.current_model()
        area = model in dft._AREA_AMPLITUDE
        self.chk_lock.blockSignals(True)
        if area:
            self.chk_lock.setText("Lock relative populations to multiplicity "
                                  "(amplitude links)")
            self.chk_lock.setChecked(True)
        else:
            self.chk_lock.setText("Lock relative populations to multiplicity "
                                  "(heights — approximate)")
            self.chk_lock.setChecked(False)
        self.chk_lock.setToolTip(
            "amplitude_k = (m_k / m_0) × amplitude_0 — exact populations for "
            "an area model (gl_norm); for a height model the ratio equals a "
            "population ratio only for equal-shape lines, so read populations "
            "from the Report")
        self.chk_lock.blockSignals(False)
        wkey, is_cq = _WIDTH_KEY.get(model, (None, True))
        self.chk_width.setVisible(bool(wkey) and not is_cq)
        if wkey and not is_cq:
            self.chk_width.setText(f"One shared linewidth ({wkey})")
        self._refresh()

    # ------------------------------------------------------------ conversion
    def _ref_set(self) -> bool:
        return self.ref.value() > self.ref.minimum()

    def _ref_typed(self, *_):
        if self._ref_set():
            self.rb_ref.setChecked(True)
        self._refresh()

    def _conversion(self):
        from larmor import shiftcal

        if self.rb_cal.isChecked():
            return self.calibration
        if self.rb_ref.isChecked() and self._ref_set():
            return shiftcal.from_sigma_ref(float(self.ref.value()),
                                           self.iso.currentData() or "")
        return None

    def _use_line(self, cal):
        self.calibration = cal
        self.rb_cal.setChecked(True)
        src = Path(cal.file).name if getattr(cal, "file", "") else "fitted here"
        self.cal_lbl.setText(f"{src} · {cal.describe()}")
        self.cal_lbl.setToolTip("\n".join(cal.flags) if cal.flags else "")
        self.tabs.setCurrentWidget(self.sites_tab)
        self._refresh()

    def _browse_calibration(self):
        from larmor.desktop.paths import remember_dir, remembered_dir

        p, _ = QFileDialog.getOpenFileName(
            self, "Load calibration line", remembered_dir(SHIFTCAL_DIR_KEY),
            "LARMOR shift calibration (*.shiftcal.json *.json);;All (*)")
        if not p:
            return
        remember_dir(SHIFTCAL_DIR_KEY, p)
        self.load_calibration(p)

    def load_calibration(self, path: str):
        try:
            cal = self.cal.load(path)
        except Exception as exc:  # noqa: BLE001
            self.warn.setText(f"failed: {exc}")
            self.warn.setStyleSheet(f"color: {_RED};")
            return None
        self._use_line(cal)
        return cal

    # ------------------------------------------------------------ table
    def _refresh(self, *_):
        from larmor import dft

        self.table.setRowCount(0)
        self._groups = []
        if self.mf is not None:
            iso = self.iso.currentData()
            sites = dft.sites_for_isotope(self.mf.sites, iso) if iso else []
            if self.chk_group.isChecked():
                groups = dft.group_equivalent(sites, tol_ppm=float(self.tol.value()))
            else:
                groups = list(sites)
            self._groups = groups
            cal = self._conversion()
            spin = self._spin()
            half = spin is not None and spin <= 0.5
            self.table.setRowCount(len(groups))
            self.table.horizontalHeaderItem(5).setToolTip(
                CQ_TIP_HALF if half else "C_Q from the EFG with the isotope's "
                                         "quadrupole moment (mrsimulator)")
            for r, g in enumerate(groups):
                sh, q = g.shielding(), g.quadrupolar()
                mult = f"×{g.multiplicity}" if g.multiplicity > 1 else "1"
                sig = f"{sh['iso_ppm']:.2f}" if sh else "—"
                if sh and cal is not None:
                    pred = _fmt_pm(cal.delta(sh["iso_ppm"]),
                                   cal.delta_err(sh["iso_ppm"]), 1)
                else:
                    pred = "—"
                zeta = (f"{sh['zeta_ppm']:.1f} / {sh['eta']:.2f}" if sh else "—")
                cq = "—" if half else (f"{q['Cq_MHz']:.3f} / {q['eta']:.2f}"
                                       if q else "—")
                note = "; ".join(n for n in g.notes)
                cells = [g.label, mult, sig, pred, zeta, cq, note]
                for c, txt in enumerate(cells):
                    tip = ", ".join(g.members) if (c == 1 and g.members) else ""
                    self.table.setItem(r, c, _ro(txt, tip))
        n = len(self._groups)
        self.btnAdd.setText(f"Add {n} site{'s' if n != 1 else ''} to the fit"
                            if n else "Add sites to the fit")
        self._update_add_enabled()

    def _update_add_enabled(self, *_):
        from larmor import dft, shiftcal

        blocking: list[str] = []
        soft: list[str] = list(self._file_notes)
        cal = self._conversion()
        if self.mf is None:
            blocking.append("open a .magres file")
        elif not self._groups:
            blocking.append("the file has no site of this isotope")
        if cal is None:
            blocking.append(NO_CONVERSION)
        if self.recipe is None:
            blocking.append("load a spectrum first — the sites are added to "
                            "the open fit")
        self._mismatch = []
        if cal is not None and self.mf is not None:
            self._mismatch = shiftcal.compatible(
                cal.calc, self.mf.calc, dft.element_of(self.iso.currentData() or ""))
        self.chk_override.setVisible(bool(self._mismatch))
        if self._mismatch:
            msg = ("the calibration line was computed with different "
                   "[calculation] settings than this file: "
                   + ", ".join(self._mismatch))
            if self.chk_override.isChecked():
                soft.append(msg + " (applied anyway — recorded)")
            else:
                blocking.append(msg)
        if cal is not None and cal.kind == "fit" and cal.dof == 0:
            soft.append(TWO_POINT)
        if blocking:
            self.warn.setText(" · ".join(blocking))
            self.warn.setStyleSheet(f"color: {_RED};")
        elif soft:
            self.warn.setText(" · ".join(soft))
            self.warn.setStyleSheet(f"color: {_AMBER};")
        else:
            self.warn.setText("")
        self.btnAdd.setEnabled(not blocking)
        self.btnAdd.setToolTip("; ".join(blocking) if blocking else
                               "append the sites below to the open fit (undoable)")

    # ------------------------------------------------------------ accept
    def _accept(self):
        from larmor import dft

        self._update_add_enabled()
        if not self.btnAdd.isEnabled() or self.mf is None:
            return
        iso = self.iso.currentData() or ""
        model = self.current_model()
        cal = self._conversion()
        lock = self.chk_lock.isChecked()
        width = self.chk_width.isChecked() and not self.chk_width.isHidden()
        first = len((self.recipe or {}).get("sites", []))
        try:
            dicts, notes = dft.sites_to_recipe_dicts(
                self._groups, model, calibration=cal, first_index=first,
                lock_amplitude=lock, share_width=width)
        except ValueError as exc:
            self.warn.setText(str(exc))
            self.warn.setStyleSheet(f"color: {_RED};")
            return
        override = list(self._mismatch) if self.chk_override.isChecked() else []
        prov = dft.import_provenance(self.mf, iso, model, self._groups, cal,
                                     calc_override=override,
                                     lock_amplitude=lock, share_width=width)
        if cal.kind == "sigma_ref":
            conv = f"σ_ref = {cal.intercept:.1f} ppm"
        else:
            src = Path(cal.file).name if cal.file else "calibration line"
            sa = cal.slope_err
            conv = f"{src}, a = {cal.slope:.4f}" + \
                (f" ± {sa:.4f}" if sa is not None else "")
        note = (f"DFT import: {len(dicts)} site(s) from {self.mf.name} "
                f"({self.mf.calc.describe()}) as {model}, δ via {conv}; "
                + "; ".join(notes))
        self.result = MagresImport(
            sites=dicts, note=note, provenance=prov, model=model,
            n_groups=len(self._groups),
            ratios=[int(g.multiplicity) for g in self._groups],
            conversion=conv, lock_amplitude=lock, share_width=width,
            area_amplitude=model in dft._AREA_AMPLITUDE)
        self.accept()

    def _help(self):
        from larmor.desktop.help_dialog import show_help

        show_help(self, "dft-tensors", "DFT tensors — import & shift calibration")
