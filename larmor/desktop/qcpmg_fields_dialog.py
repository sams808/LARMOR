"""Two-field (or multi-field) infinite-field extrapolation of the isotropic
chemical shift from the central-transition centre of gravity — for QCPMG data
of half-integer quadrupolar nuclei measured at more than one magnetic field
(Sandland et al. 2004; Baasner et al. 2014). See larmor.qcpmg_fields."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from larmor.desktop import theme
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QHBoxLayout, QLabel,
    QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from larmor.qcpmg_fields import (
    FieldPoint, centre_of_gravity, infinite_field_diso,
)


def _spin_of(nucleus: str) -> float:
    try:
        from mrsimulator.spin_system.isotope import ISOTOPE_DATA
        d = ISOTOPE_DATA.get(nucleus or "")
        if d:
            return (d["spin_multiplicity"] - 1) / 2.0
    except Exception:
        pass
    return 1.5


class QcpmgFieldsDialog(QDialog):
    def __init__(self, parent, nucleus: str = "", current=None):
        super().__init__(parent)
        self.setWindowTitle("QCPMG — infinite-field δiso (2+ fields)")
        self.resize(760, 680)
        self._nucleus = nucleus or "35Cl"
        self._current = current            # (larmor_MHz, ppm, amp) of the open spectrum
        v = QVBoxLayout(self)

        intro = QLabel(
            "Extrapolate the central-transition centre of gravity to infinite "
            "field to remove the second-order quadrupolar shift and obtain the "
            "true isotropic shift δiso and C_Q (Sandland 2004 Eq. 1 / Baasner "
            "2014 Fig. 6). Enter δcg at each field, or grab it from the open "
            "spectrum. In the large-C_Q limit CT-selective and non-selective "
            "fields can be combined (see the Lineshapes/QCPMG manual).")
        intro.setWordWrap(True)
        intro.setStyleSheet(f"color:{theme.active().text_dim};")
        v.addWidget(intro)

        top = QHBoxLayout()
        top.addWidget(QLabel(f"nucleus <b>{self._nucleus}</b>  ·  spin I ="))
        self.spin = QDoubleSpinBox(); self.spin.setDecimals(1)
        self.spin.setRange(1.5, 4.5); self.spin.setSingleStep(1.0)
        self.spin.setValue(_spin_of(self._nucleus))
        top.addWidget(self.spin)
        top.addSpacing(16)
        top.addWidget(QLabel("η (assumed)"))
        self.eta = QDoubleSpinBox(); self.eta.setRange(0.0, 1.0)
        self.eta.setSingleStep(0.05); self.eta.setValue(0.7)
        self.eta.setToolTip("η is not determined by two centres of gravity; "
                            "0.7 is the conventional choice (Stebbins & Du 2002)")
        top.addWidget(self.eta)
        top.addStretch(1)
        v.addLayout(top)

        # per-field table: Larmor (MHz), δcg (ppm), ±err, FWHM (ppm), CT-selective
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["Larmor ν₀ (MHz)", "δcg (ppm)", "± err (ppm)",
             "FWHM (ppm)", "CT-selective"])
        self.table.horizontalHeader().setStretchLastSection(True)
        # a QTableWidget's own minimum is ~1.5 rows, so at the default dialog
        # size the plot's fixed minimum squeezed the table nearly out of view
        self.table.setMinimumHeight(150)             # header + ~4 field rows
        v.addWidget(self.table, 1)

        row = QHBoxLayout()
        b_add = QPushButton("＋ Add field")
        b_add.clicked.connect(lambda: self._add_row())
        b_del = QPushButton("Remove selected")
        b_del.clicked.connect(self._del_row)
        b_ds = QPushButton("Add from datasets…")
        b_ds.setToolTip("pick the processed spectra (1r) measured at each "
                        "field: δcg ± σ and FWHM are read off automatically, "
                        "and selecting the row shows the spectrum with a "
                        "draggable band to supervise the values")
        b_ds.clicked.connect(self._pick_datasets)
        self.b_cur = QPushButton("δcg from open spectrum (visible range)")
        self.b_cur.setToolTip("centre of gravity of the currently open spectrum "
                              "over the visible x-range — zoom to the CT band first")
        self.b_cur.clicked.connect(self._from_current)
        self.b_cur.setEnabled(self._current is not None)
        row.addWidget(b_add); row.addWidget(b_del); row.addWidget(b_ds)
        row.addWidget(self.b_cur)
        row.addStretch(1)
        v.addLayout(row)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.setLabel("bottom", "1 / ν₀²", units="MHz⁻²")
        self.plot.setLabel("left", "δcg", units="ppm")
        self.plot.getPlotItem().getAxis("left").enableAutoSIPrefix(False)
        self.plot.showGrid(x=True, y=True, alpha=0.15)
        self.plot.setMinimumHeight(140)
        v.addWidget(self.plot, 1)        # table and plot share extra height

        self.result = QLabel("add at least two fields, then Compute")
        self.result.setStyleSheet(f"font-weight:600; color:{theme.active().accent};")
        self.result.setWordWrap(True)
        v.addWidget(self.result)

        self.wresult = QLabel(
            "Two-field linewidth split (Sandland Eq. 2): also fill FWHM (ppm) "
            "at both fields, then 'Split W_q / W_csd'.")
        self.wresult.setStyleSheet(f"color:{theme.active().text_dim};")
        self.wresult.setWordWrap(True)
        self.wresult.setTextFormat(Qt.RichText)
        v.addWidget(self.wresult)

        bb = QDialogButtonBox()
        b_comp = bb.addButton("Compute δiso", QDialogButtonBox.ApplyRole)
        b_comp.clicked.connect(self._compute)
        b_w = bb.addButton("Split W_q / W_csd", QDialogButtonBox.ApplyRole)
        b_w.clicked.connect(self._compute_widths)
        bb.addButton(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.addButton(QDialogButtonBox.Help).clicked.connect(self._help)
        v.addWidget(bb)

        # dataset supervision: rows added from a dataset keep their spectrum,
        # and selecting one shows it with a draggable band
        self._ds: dict[int, dict] = {}
        self._ds_seq = 0
        self._region = None
        self._cg_line = None
        self.table.itemSelectionChanged.connect(self._show_selected_dataset)

        # seed two rows; prefill the first from the open spectrum's field
        self._add_row(self._current[0] if self._current else 0.0)
        self._add_row()

    # ------------------------------------------------- datasets (supervised)
    def _pick_datasets(self):
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Choose the processed spectrum (1r) measured at each field",
            "", "Bruker processed (1r);;All files (*)")
        errors = []
        for p in paths:
            try:
                from larmor.io import bruker
                d = bruker.read(p)
                if d.ndim != 1 or d.domain != "freq":
                    raise ValueError("not a processed 1D spectrum")
                ppm = np.asarray(d.axes[0].values, float)
                amp = np.asarray(d.data, float)
                self.add_dataset_spectrum(
                    float(d.meta.get("larmor_MHz", 0.0) or 0.0), ppm, amp)
            except Exception as exc:                          # noqa: BLE001
                errors.append(f"{p}: {exc}")
        if errors:
            QMessageBox.warning(self, "Some datasets were skipped",
                                "\n".join(errors))

    def add_dataset_spectrum(self, larmor_MHz: float, ppm, amp,
                             window=None, magnitude: bool = False) -> int:
        """Add one field's spectrum: δcg ± σ and FWHM are computed over the
        given window (or an automatic one) and written into a new row; the
        spectrum stays attached so selecting the row shows it for
        supervision."""
        from larmor import qcpmg
        ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
        if window is not None:
            lo, hi = float(min(window)), float(max(window))
        else:
            hi, lo = qcpmg.cg_window(ppm, amp)
        if not (np.isfinite(hi) and np.isfinite(lo)) or hi <= lo:
            span = float(ppm.max() - ppm.min())
            mid = float(ppm.min()) + span / 2.0
            lo, hi = mid - span / 6.0, mid + span / 6.0
        self._ds_seq += 1
        ds_id = self._ds_seq
        self._ds[ds_id] = {"ppm": ppm, "amp": amp, "window": (lo, hi),
                           "magnitude": bool(magnitude)}
        self._add_row(larmor_MHz)
        r = self.table.rowCount() - 1
        self.table.item(r, 0).setData(Qt.UserRole, ds_id)
        if magnitude:
            it = self.table.item(r, 1)
            if it is not None:
                it.setToolTip("δcg measured on a MAGNITUDE (mc) spectrum")
        self._apply_ds_values(r, ds_id)
        self._warn_mixed_modes()
        self.table.selectRow(r)
        return r

    def _warn_mixed_modes(self):
        """Sandland's extrapolation compares centres of gravity ACROSS
        fields; one measured in magnitude and one in absorption are not the
        same observable, so say so rather than fitting them together
        silently."""
        modes = {bool(d.get("magnitude")) for d in self._ds.values()}
        if len(modes) > 1:
            self.wresult.setText(
                "<span style='color:#c0392b'>⚠ these fields mix magnitude "
                "(mc) and absorption δcg values — they are different "
                "observables; reprocess them the same way before "
                "extrapolating.</span>")

    def _apply_ds_values(self, r: int, ds_id: int):
        from larmor import qcpmg
        d = self._ds[ds_id]
        lo, hi = d["window"]
        cg, sigma = qcpmg.centre_of_gravity(d["ppm"], d["amp"], (hi, lo),
                                            jitter_frac=0.10)
        fw_ppm = qcpmg.fwhm_hz(d["ppm"], d["amp"], 1.0, (hi, lo))
        if np.isfinite(cg):
            self.table.setItem(r, 1, QTableWidgetItem(f"{cg:.2f}"))
            self.table.setItem(r, 2, QTableWidgetItem(f"{max(sigma, 0.1):.1f}"))
        self.table.setItem(r, 3, QTableWidgetItem(f"{fw_ppm:.2f}"))
        if self._cg_line is not None and np.isfinite(cg):
            self._cg_line.setValue(cg)

    def _row_ds_id(self, r: int):
        it = self.table.item(r, 0)
        return it.data(Qt.UserRole) if it is not None else None

    def _show_selected_dataset(self):
        r = self.table.currentRow()
        ds_id = self._row_ds_id(r) if r >= 0 else None
        if ds_id is None or ds_id not in self._ds:
            return
        d = self._ds[ds_id]
        self.plot.clear()
        self._region = self._cg_line = None
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "shift", units="ppm")
        self.plot.setLabel("left", "intensity", units="")
        self.plot.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
        self.plot.setTitle("dataset — drag the band edges; δcg / FWHM in the "
                           "row follow", color=theme.active().text_dim,
                           size="9pt")
        self.plot.plot(d["ppm"], d["amp"],
                       pen=pg.mkPen(theme.active().experiment, width=1.2))
        lo, hi = d["window"]
        self._region = pg.LinearRegionItem(values=(lo, hi), movable=True)
        self.plot.addItem(self._region)
        self._region.sigRegionChangeFinished.connect(
            lambda *_: self._region_moved(r, ds_id))
        self._cg_line = pg.InfiniteLine(pos=0.0, angle=90, movable=False,
                                        pen=pg.mkPen(theme.active().pivot,
                                                     style=Qt.DashLine))
        self.plot.addItem(self._cg_line)
        self._apply_ds_values(r, ds_id)

    def _region_moved(self, r: int, ds_id: int):
        if self._region is None or ds_id not in self._ds:
            return
        a, b = self._region.getRegion()
        self._ds[ds_id]["window"] = (min(a, b), max(a, b))
        if r < self.table.rowCount() and self._row_ds_id(r) == ds_id:
            self._apply_ds_values(r, ds_id)

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "qcpmg", "QCPMG")

    def _add_row(self, larmor=0.0):
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, 0, QTableWidgetItem(f"{larmor:g}" if larmor else ""))
        self.table.setItem(r, 1, QTableWidgetItem(""))
        self.table.setItem(r, 2, QTableWidgetItem("5"))
        self.table.setItem(r, 3, QTableWidgetItem(""))          # FWHM (ppm)
        chk = QCheckBox(); chk.setChecked(True)
        w = QWidget(); lay = QHBoxLayout(w); lay.setContentsMargins(0, 0, 0, 0)
        lay.setAlignment(Qt.AlignCenter); lay.addWidget(chk)
        self.table.setCellWidget(r, 4, w)

    def _del_row(self):
        r = self.table.currentRow()
        if r >= 0:
            self.table.removeRow(r)

    def _from_current(self):
        if self._current is None:
            return
        larmor, ppm, amp = self._current
        (x0, x1) = self.parent().view.getPlotItem().getViewBox().viewRange()[0] \
            if hasattr(self.parent(), "view") else (None, None)
        cg = centre_of_gravity(ppm, amp, x0, x1)
        # write into the first empty δcg cell (or a new row)
        for r in range(self.table.rowCount()):
            if not (self.table.item(r, 1) and self.table.item(r, 1).text().strip()):
                self.table.setItem(r, 0, QTableWidgetItem(f"{larmor:g}"))
                self.table.setItem(r, 1, QTableWidgetItem(f"{cg:.2f}"))
                return
        self._add_row(larmor)
        self.table.setItem(self.table.rowCount() - 1, 1,
                           QTableWidgetItem(f"{cg:.2f}"))

    def _points(self) -> list[FieldPoint]:
        pts = []
        for r in range(self.table.rowCount()):
            try:
                nu = float(self.table.item(r, 0).text())
                dcg = float(self.table.item(r, 1).text())
            except (AttributeError, ValueError):
                continue
            try:
                err = float(self.table.item(r, 2).text())
            except (AttributeError, ValueError):
                err = 5.0
            w = self.table.cellWidget(r, 4)
            sel = w.findChild(QCheckBox).isChecked() if w else True
            pts.append(FieldPoint(nu, dcg, err, sel))
        return pts

    def _fields_fwhm(self):
        """(larmor, fwhm_ppm) for rows that have both filled — for Eq. 2."""
        out = []
        for r in range(self.table.rowCount()):
            try:
                nu = float(self.table.item(r, 0).text())
                fw = float(self.table.item(r, 3).text())
                out.append((nu, fw))
            except (AttributeError, ValueError):
                continue
        return out

    def _compute_widths(self):
        from larmor.qcpmg_fields import two_field_widths
        fw = self._fields_fwhm()
        if len(fw) < 2:
            self.wresult.setText("enter the FWHM (ppm) at two fields")
            return
        (n1, f1), (n2, f2) = fw[0], fw[1]
        ws = two_field_widths(n1, f1, n2, f2)
        if not ws.ok:
            self.wresult.setText("⚠ " + ws.note)
            return
        self.wresult.setText(
            f"quadrupolar width W_q = {ws.wq_lo_ppm:.1f} ppm (at {min(n1,n2):.0f} "
            f"MHz) / {ws.wq_hi_ppm:.1f} ppm (at {max(n1,n2):.0f} MHz)  ·  "
            f"chemical-shift-distribution width W_csd = <b>{ws.wcsd_ppm:.1f} "
            f"ppm</b> (field-independent)")

    def _compute(self):
        pts = self._points()
        if len(pts) < 2:
            self.result.setText("need δcg at ≥ 2 fields (fill Larmor + δcg)")
            return
        try:
            res = infinite_field_diso(pts, spin=self.spin.value(),
                                      eta=self.eta.value())
        except Exception as exc:
            self.result.setText(f"cannot extrapolate: {exc}")
            return
        self.plot.clear()
        # the plot may be in dataset-supervision mode (inverted ppm axis)
        self._region = self._cg_line = None
        self.plot.getPlotItem().invertX(False)
        self.plot.setTitle(None)
        self.plot.setLabel("bottom", "1 / ν₀²", units="MHz⁻²")
        self.plot.setLabel("left", "δcg", units="ppm")
        x = np.array([1.0 / p.larmor_MHz ** 2 for p in pts])
        y = np.array([p.dcg_ppm for p in pts])
        self.plot.plot(x, y, pen=None, symbol="o", symbolBrush="#1f6feb",
                       symbolSize=9)
        xs = np.linspace(0.0, float(x.max()) * 1.05, 50)
        self.plot.plot(xs, res.line(xs),
                       pen=pg.mkPen("#c0392b", width=1.6, style=Qt.DashLine))
        self.plot.plot([0.0], [res.delta_iso_ppm], pen=None, symbol="star",
                       symbolBrush="#c0392b", symbolSize=14)
        self.result.setText(
            f"δiso = <b>{res.delta_iso_ppm:.1f} ± {res.delta_iso_err_ppm:.1f} "
            f"ppm</b>  (intercept, 1/ν₀²→0)   ·   "
            f"C_Q = {res.cq_MHz:.2f} ± {res.cq_err_MHz:.2f} MHz   ·   "
            f"P_Q = {res.pq_MHz:.2f} MHz   (η = {res.eta:g} assumed)")


#: the one shared instance — fields sent from several QCPMG processing
#: sessions accumulate here until the user computes the extrapolation
_shared: QcpmgFieldsDialog | None = None


def shared_fields_dialog(parent=None, nucleus: str = "",
                         current=None) -> QcpmgFieldsDialog:
    """Get (or create) the persistent infinite-field dialog. Non-modal by
    design: process one field's dataset, send it here, process the next,
    send it too — then Compute."""
    global _shared
    if _shared is not None:
        try:
            _shared.isVisible()               # raises once the C++ side died
        except RuntimeError:
            _shared = None
    if _shared is None:
        _shared = QcpmgFieldsDialog(parent, nucleus, current)
    elif current is not None:
        _shared._current = current
        _shared.b_cur.setEnabled(True)
    return _shared
