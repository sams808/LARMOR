"""Referencing audit window (Tools > Referencing audit): scan one session
folder, find its 1H adamantane reference spectra, check them, and compare
every other acquisition's stored SR with the value indirect referencing
would give. Outputs a TopSpin-ready ``sr`` list plus a CSV, appends the old
and new values to the permanent log, and can re-reference the spectrum open
in the workbench. Logic lives in larmor.referencing (Qt-free, tested)."""
from __future__ import annotations

from pathlib import Path

import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QDoubleSpinBox, QFileDialog, QHBoxLayout,
    QHeaderView, QLabel, QLineEdit, QMessageBox, QPushButton, QSplitter,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from larmor import referencing as R
from larmor.desktop import theme

_VERDICT_COLOR = {"ok": "#2C6A4E", "off": "#9E3324", "unreferenced": "#9E3324",
                  "no reference": "#A8570F", "unprocessed": "#6A7C80",
                  "1H reference": "#0B6A71", "1H": "#6A7C80"}


class ReferencingAuditDialog(QDialog):
    #: (expno_path, new_sr_hz, old_sr_hz, note) -- re-reference the open spectrum
    apply_requested = Signal(str, float, float, str)

    def __init__(self, parent=None, start_dir: str | None = None,
                 current_source: str | None = None):
        super().__init__(parent)
        self.setWindowTitle("Referencing audit — SR against the session's ¹H reference")
        self.setModal(False)
        self.resize(1180, 720)
        self.current_source = current_source
        self.acqs: list = []
        self.rows: list = []
        self.checks: dict = {}
        t = theme.active()
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        top.addWidget(QLabel("session folder (one month, e.g. …/DATA/2026-05)"))
        self.folder = QLineEdit(start_dir or "")
        top.addWidget(self.folder, 1)
        b = QPushButton("Browse…"); b.clicked.connect(self._browse); top.addWidget(b)
        self.btnScan = QPushButton("Scan"); self.btnScan.clicked.connect(self.scan)
        top.addWidget(self.btnScan)
        v.addLayout(top)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("adamantane ¹H at"))
        self.ada = QDoubleSpinBox(); self.ada.setDecimals(2); self.ada.setRange(-5.0, 10.0)
        self.ada.setValue(R.ADAMANTANE_1H_PPM); self.ada.setSuffix(" ppm")
        opts.addWidget(self.ada)
        opts.addWidget(QLabel("tolerance"))
        self.tol = QDoubleSpinBox(); self.tol.setDecimals(3); self.tol.setRange(0.001, 5.0)
        self.tol.setValue(R.DEFAULT_TOL_PPM); self.tol.setSuffix(" ppm")
        opts.addWidget(self.tol)
        self.chkPeak = QCheckBox("derive each ¹H reference from its adamantane "
                                 "peak (instead of the stored SF)")
        opts.addWidget(self.chkPeak)
        opts.addStretch(1)
        v.addLayout(opts)
        for w in (self.ada, self.tol):
            w.valueChanged.connect(self._recompute)
        self.chkPeak.toggled.connect(self._recompute)

        self.status = QLabel("pick the session folder and press Scan")
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color: {t.text_dim};")
        v.addWidget(self.status)

        split = QSplitter(Qt.Horizontal)
        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(QLabel("¹H reference spectra found in the session"))
        self.refs = QTableWidget(0, 6)
        self.refs.setHorizontalHeaderLabels(["sample/EXPNO", "date", "SF (MHz)",
                                             "SR (Hz)", "peak (ppm)", "check"])
        self.refs.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.refs.horizontalHeader().setStretchLastSection(True)
        self.refs.setSelectionBehavior(QTableWidget.SelectRows)
        self.refs.itemSelectionChanged.connect(self._show_ref)
        lv.addWidget(self.refs, 1)
        self.plot = pg.PlotWidget(background=t.plot_bg)
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "¹H shift", units="ppm")
        self.plot.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
        lv.addWidget(self.plot, 1)
        split.addWidget(left)

        right = QWidget(); rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0)
        rv.addWidget(QLabel("every acquisition against its nearest ¹H reference"))
        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(["sample/EXPNO", "nucleus", "date",
                                              "stored SR (Hz)", "expected SR (Hz)",
                                              "Δ (ppm)", "verdict", "note"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.itemSelectionChanged.connect(self._row_selected)
        rv.addWidget(self.table, 1)
        split.addWidget(right)
        split.setSizes([480, 700])
        v.addWidget(split, 1)

        bottom = QHBoxLayout()
        self.chkAll = QCheckBox("list every EXPNO, not only the flagged ones")
        bottom.addWidget(self.chkAll)
        bottom.addStretch(1)
        self.btnCopy = QPushButton("Copy TopSpin sr list")
        self.btnCopy.clicked.connect(self._copy)
        self.btnExport = QPushButton("Export list (CSV + TopSpin text)…")
        self.btnExport.clicked.connect(self.export)
        self.btnApply = QPushButton("Re-reference the open spectrum")
        self.btnApply.setToolTip("adds the SR correction to the spectrum open in the "
                                 "workbench (recipe provenance keeps the old value); "
                                 "the instrument folder is not touched")
        self.btnApply.clicked.connect(self._apply)
        btnClose = QPushButton("Close"); btnClose.clicked.connect(self.close)
        for w in (self.btnCopy, self.btnExport, self.btnApply, btnClose):
            bottom.addWidget(w)
        v.addLayout(bottom)
        self._enable(False)

    # ------------------------------------------------------------------
    def _enable(self, on: bool):
        for w in (self.btnCopy, self.btnExport):
            w.setEnabled(on)
        self.btnApply.setEnabled(on and self._current_row() is not None)

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Session folder (one month)",
                                             self.folder.text() or str(Path.home()))
        if d:
            self.folder.setText(d)

    def scan(self):
        folder = self.folder.text().strip()
        if not folder or not Path(folder).is_dir():
            self.status.setText("that folder does not exist")
            return
        root = R.session_root(folder)
        self.folder.setText(str(root))
        self.status.setText(f"scanning {root} …")
        QApplication.processEvents()
        self.acqs = R.scan_session(root)
        refs = R.find_references(self.acqs)
        self.checks = {r.path: R.check_reference(r, float(self.ada.value())) for r in refs}
        self._recompute()
        prev = R.previous_audits(root)
        warn = (f"  ·  ⚠ this session was audited before ({prev[-1]['time'][:16]}, "
                f"{len(prev)} record(s) in the log) — the audit is meant to run once"
                if prev else "")
        n1h = len(refs)
        self.status.setText(
            f"{len(self.acqs)} acquisitions in {len({a.sample for a in self.acqs})} "
            f"sample folders · {n1h} referenced ¹H spectrum(s)" + warn
            if self.acqs else "no EXPNO found — pick the month folder that holds the sample folders")
        self._enable(bool(self.acqs))

    def _recompute(self, *_):
        if not self.acqs:
            return
        ada = float(self.ada.value())
        if self.chkPeak.isChecked():
            for path, c in self.checks.items():
                if c.peak_ppm is not None:
                    pass
        ref_sf = {}
        if self.chkPeak.isChecked():
            for path, c in self.checks.items():
                if c.peak_ppm is not None:
                    ref_sf[path] = R.sf_from_peak(c.acq.sf_MHz, c.peak_ppm, ada)
        self.rows = R.audit(self.acqs, tol_ppm=float(self.tol.value()), ada_ppm=ada,
                            ref_sf=ref_sf or None)
        self._fill_refs()
        self._fill_rows()

    def _fill_refs(self):
        refs = R.find_references(self.acqs)
        self.refs.setRowCount(len(refs))
        for i, r in enumerate(refs):
            c = self.checks.get(r.path)
            vals = [r.label, r.date_iso, f"{r.sf_MHz:.6f}" if r.sf_MHz else "",
                    f"{r.sr_hz:.2f}" if r.sr_hz is not None else "",
                    f"{c.peak_ppm:.2f}" if c and c.peak_ppm is not None else "",
                    c.status if c else ""]
            for j, txt in enumerate(vals):
                it = QTableWidgetItem(txt)
                it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                it.setData(Qt.UserRole, r.path)
                if j == 5 and c is not None:
                    it.setForeground(QColor("#2C6A4E" if c.ok else "#9E3324"))
                self.refs.setItem(i, j, it)
        if refs and self.refs.currentRow() < 0:
            self.refs.selectRow(0)

    def _fill_rows(self):
        shown = [r for r in self.rows if self.chkAll.isChecked() or r.flagged
                 or r.verdict in ("no reference", "unprocessed")]
        self.table.setRowCount(len(shown))
        for i, r in enumerate(shown):
            vals = [r.acq.label, r.acq.nucleus, r.acq.date_iso,
                    f"{r.acq.sr_hz:.2f}" if r.acq.sr_hz is not None else "—",
                    f"{r.expected_sr_hz:.2f}" if r.expected_sr_hz is not None else "—",
                    f"{r.delta_ppm:+.3f}" if r.delta_ppm is not None else "—",
                    r.verdict, "; ".join([r.note] + r.notes).strip("; ")]
            for j, txt in enumerate(vals):
                it = QTableWidgetItem(txt)
                it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                it.setData(Qt.UserRole, r.acq.path)
                if j == 6:
                    it.setForeground(QColor(_VERDICT_COLOR.get(r.verdict, "#6A7C80")))
                self.table.setItem(i, j, it)
        counts = R.summary(self.rows)
        self.setWindowTitle("Referencing audit — " + "  ·  ".join(
            f"{k}: {n}" for k, n in sorted(counts.items())))
        self.btnApply.setEnabled(self._current_row() is not None)

    def _show_ref(self):
        sel = self.refs.selectedItems()
        if not sel:
            return
        c = self.checks.get(sel[0].data(Qt.UserRole))
        self.plot.clear()
        if c is None or c.ppm is None:
            return
        t = theme.active()
        self.plot.plot(c.ppm, c.amp, pen=pg.mkPen(t.experiment, width=1.2))
        ada = float(self.ada.value())
        self.plot.addItem(pg.InfiniteLine(pos=ada, angle=90,
                                          pen=pg.mkPen(t.accent, width=2)))
        if c.peak_ppm is not None:
            self.plot.addItem(pg.InfiniteLine(pos=c.peak_ppm, angle=90,
                                              pen=pg.mkPen(t.model, width=1,
                                                           style=Qt.DashLine)))
        span = 8.0
        self.plot.setXRange(ada + span, ada - span)

    def _current_row(self):
        """The audit row of the spectrum open in the workbench, if listed."""
        if not self.current_source or not self.rows:
            return None
        src = Path(self.current_source)
        for r in self.rows:
            p = Path(r.acq.path)
            if src == p or p in src.parents:
                return r if r.expected_sr_hz is not None and not r.acq.is_1h else None
        return None

    def _row_selected(self):
        pass

    # ------------------------------------------------------------------
    def _text(self) -> str:
        return R.topspin_list(self.rows, only_flagged=not self.chkAll.isChecked(),
                              session=self.folder.text())

    def _copy(self):
        QApplication.clipboard().setText(self._text())
        self.status.setText("TopSpin sr list copied to the clipboard")

    def export(self):
        if not self.rows:
            return
        root = Path(self.folder.text())
        default = str(Path.home() / f"referencing_audit_{root.name}.csv")
        path, _ = QFileDialog.getSaveFileName(self, "Export the audit", default,
                                              "CSV (*.csv)")
        if not path:
            return
        p = Path(path)
        try:
            R.to_csv(self.rows, p)
            txt = p.with_name(p.stem + "_topspin.txt")
            txt.write_text(self._text(), encoding="utf-8")
            log = R.append_log(self.rows, root, ada_ppm=float(self.ada.value()),
                               tol_ppm=float(self.tol.value()), action="export")
        except OSError as exc:
            QMessageBox.warning(self, "Export", f"could not write: {exc}")
            return
        self.status.setText(f"wrote {p.name} and {txt.name}; old and new values "
                            f"appended to the log {log}")

    def _apply(self):
        r = self._current_row()
        if r is None:
            QMessageBox.information(self, "Re-reference",
                                    "Open one of the listed spectra in the workbench first.")
            return
        note = (f"SR {r.acq.sr_hz:.2f} → {r.expected_sr_hz:.2f} Hz from the ¹H reference "
                f"{r.ref.label} (adamantane at {self.ada.value():g} ppm, Ξ indirect)")
        try:
            R.append_log([r], Path(self.folder.text()), ada_ppm=float(self.ada.value()),
                         tol_ppm=float(self.tol.value()), action="apply")
        except OSError:
            pass
        self.apply_requested.emit(r.acq.path, float(r.expected_sr_hz),
                                  float(r.acq.sr_hz or 0.0), note)
        self.status.setText("correction applied to the open spectrum (undoable)")
