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
        self.resize(720, 560)
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
        v.addWidget(self.table, 1)

        row = QHBoxLayout()
        b_add = QPushButton("＋ Add field")
        b_add.clicked.connect(lambda: self._add_row())
        b_del = QPushButton("Remove selected")
        b_del.clicked.connect(self._del_row)
        self.b_cur = QPushButton("δcg from open spectrum (visible range)")
        self.b_cur.setToolTip("centre of gravity of the currently open spectrum "
                              "over the visible x-range — zoom to the CT band first")
        self.b_cur.clicked.connect(self._from_current)
        self.b_cur.setEnabled(self._current is not None)
        row.addWidget(b_add); row.addWidget(b_del); row.addWidget(self.b_cur)
        row.addStretch(1)
        v.addLayout(row)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.setLabel("bottom", "1 / ν₀²", units="MHz⁻²")
        self.plot.setLabel("left", "δcg", units="ppm")
        self.plot.showGrid(x=True, y=True, alpha=0.15)
        self.plot.setMinimumHeight(180)
        v.addWidget(self.plot)

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

        # seed two rows; prefill the first from the open spectrum's field
        self._add_row(self._current[0] if self._current else 0.0)
        self._add_row()

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
