"""Batch infinite-field delta_iso: many samples, several fields, one table.

Set how many samples and how many fields, drop the processed spectra (the
.csv files stage 5 of the QCPMG dialog writes) onto the grid, and every cell
is measured the same way -- centre of gravity and FWHM over a window you can
still supervise cell by cell. Then one click extrapolates every sample to
infinite field and writes the report and the figures.

The measurement of a cell and the extrapolation itself are the same core
functions the single-sample dialog uses (larmor.qcpmg, larmor.qcpmg_fields),
so the two cannot disagree.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog, QHBoxLayout,
    QHeaderView, QLabel, QMessageBox, QPlainTextEdit, QSpinBox, QSplitter,
    QTableWidget, QTableWidgetItem, QVBoxLayout,
)

from larmor.desktop import theme

_log = logging.getLogger(__name__)


class _DropTable(QTableWidget):
    """A grid that accepts spectrum files dropped onto individual cells."""

    def __init__(self, on_drop):
        super().__init__()
        self._on_drop = on_drop
        self.setAcceptDrops(True)
        self.viewport().setAcceptDrops(True)
        self.setDragDropMode(QTableWidget.DropOnly)

    def dragEnterEvent(self, ev):
        if ev.mimeData().hasUrls():
            ev.acceptProposedAction()

    def dragMoveEvent(self, ev):
        if ev.mimeData().hasUrls():
            ev.acceptProposedAction()

    def dropEvent(self, ev):
        if not ev.mimeData().hasUrls():
            return
        pos = ev.position().toPoint() if hasattr(ev, "position") else ev.pos()
        item = self.indexAt(pos)
        if not item.isValid():
            return
        paths = [u.toLocalFile() for u in ev.mimeData().urls() if u.isLocalFile()]
        if paths:
            # several files at once fill the row from the dropped column on,
            # which is how a sample's fields usually arrive
            self._on_drop(item.row(), item.column(), paths)
            ev.acceptProposedAction()


class QcpmgBatchFieldsDialog(QDialog):
    def __init__(self, parent=None, nucleus: str = ""):
        super().__init__(parent)
        self.setWindowTitle("QCPMG — batch infinite-field δiso")
        self.resize(1080, 800)
        self.setMinimumSize(680, 520)
        self.setSizeGripEnabled(True)
        self.setWindowFlags(self.windowFlags() | Qt.Window
                            | Qt.WindowMinimizeButtonHint
                            | Qt.WindowMaximizeButtonHint)
        self.setModal(False)
        self._nucleus = nucleus or ""
        #: (row, col) -> {"path", "ppm", "amp", "larmor", "window", "cg",
        #:                "sigma", "fwhm"}
        self.cells: dict[tuple[int, int], dict] = {}
        self._results: dict = {}
        self._region = None
        self._sel = None

        v = QVBoxLayout(self)
        intro = QLabel(
            "Set the grid, then <b>drop the processed spectra onto the "
            "cells</b> (the .csv files the QCPMG dialog writes with “Save as "
            "dataset…”), or double-click a cell to browse. Each field is one "
            "column; each sample one row. Select a cell to check and drag its "
            "δcg window.")
        intro.setWordWrap(True)
        intro.setStyleSheet(f"color:{theme.active().text_dim};")
        v.addWidget(intro)

        top = QHBoxLayout()
        self.nSamples = QSpinBox(); self.nSamples.setRange(1, 200)
        self.nSamples.setValue(3)
        self.nSamples.valueChanged.connect(self._resize_grid)
        self.nFields = QSpinBox(); self.nFields.setRange(2, 12)
        self.nFields.setValue(2)
        self.nFields.valueChanged.connect(self._resize_grid)
        self.spin = QDoubleSpinBox(); self.spin.setRange(1.5, 4.5)
        self.spin.setSingleStep(1.0); self.spin.setDecimals(1)
        self.spin.setValue(_spin_of(self._nucleus))
        self.eta = QDoubleSpinBox(); self.eta.setRange(0.0, 1.0)
        self.eta.setSingleStep(0.05); self.eta.setValue(0.7)
        self.eta.setToolTip("η is not determined by centres of gravity; 0.7 "
                            "is the conventional choice")
        self.lblNuc = QLabel(f"nucleus <b>{self._nucleus or '—'}</b>")
        for w in ("samples", self.nSamples, "  fields", self.nFields,
                  "   ", self.lblNuc, "  spin I", self.spin, "  η", self.eta):
            top.addWidget(QLabel(w) if isinstance(w, str) else w)
        top.addStretch(1)
        v.addLayout(top)

        split = QSplitter(Qt.Vertical)
        self.table = _DropTable(self._drop_files)
        self.table.itemSelectionChanged.connect(self._show_cell)
        self.table.itemDoubleClicked.connect(self._browse_cell)
        self.table.itemChanged.connect(self._cell_edited)
        split.addWidget(self.table)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "shift", units="ppm")
        self.plot.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
        self.plot.setMinimumHeight(140)
        split.addWidget(self.plot)

        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        f = self.report.font(); f.setFamily("Consolas"); self.report.setFont(f)
        self.report.setMinimumHeight(90)
        split.addWidget(self.report)
        split.setSizes([320, 240, 200])
        v.addWidget(split, 1)

        self.msg = QLabel("")
        self.msg.setWordWrap(True)
        self.msg.setStyleSheet(f"color:{theme.active().accent}; font-weight:600;")
        v.addWidget(self.msg)

        bb = QDialogButtonBox()
        b_comp = bb.addButton("Compute all", QDialogButtonBox.ApplyRole)
        b_comp.clicked.connect(self._compute)
        self.btnReport = bb.addButton("Export report…", QDialogButtonBox.ActionRole)
        self.btnReport.clicked.connect(self._export_report)
        self.btnFig = bb.addButton("Export figures…", QDialogButtonBox.ActionRole)
        self.btnFig.setToolTip("the merged figure with every sample, plus one "
                               "figure per sample, as .png + .svg + .pdf")
        self.btnFig.clicked.connect(self._export_figures)
        bb.addButton(QDialogButtonBox.Close).clicked.connect(self.close)
        bb.addButton(QDialogButtonBox.Help).clicked.connect(self._help)
        for b in (self.btnReport, self.btnFig):
            b.setEnabled(False)
        v.addWidget(bb)

        self._resize_grid()

    # ------------------------------------------------------------- grid
    def _resize_grid(self):
        t = self.table
        t.blockSignals(True)
        rows, fields = self.nSamples.value(), self.nFields.value()
        t.setRowCount(rows)
        t.setColumnCount(1 + fields)
        t.setHorizontalHeaderLabels(
            ["sample"] + [f"field {i + 1}" for i in range(fields)])
        t.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        for r in range(rows):
            if t.item(r, 0) is None:
                t.setItem(r, 0, QTableWidgetItem(f"sample {r + 1}"))
            for c in range(1, 1 + fields):
                if t.item(r, c) is None:
                    it = QTableWidgetItem("")
                    it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                    t.setItem(r, c, it)
        # cells outside the new grid are gone
        self.cells = {k: v for k, v in self.cells.items()
                      if k[0] < rows and k[1] < 1 + fields}
        t.blockSignals(False)
        self._refresh_headers()

    def _refresh_headers(self):
        """Name each field column by the Larmor frequency its spectra carry."""
        for c in range(1, self.table.columnCount()):
            fs = [d["larmor"] for (r, cc), d in self.cells.items() if cc == c]
            label = f"field {c}"
            if fs:
                label += f"\n{np.median(fs):.3f} MHz"
                if max(fs) - min(fs) > 0.05:
                    label += " ⚠"
            self.table.setHorizontalHeaderItem(c, QTableWidgetItem(label))

    def _cell_edited(self, item):
        if item.column() == 0:
            self._set_dirty()

    def _set_dirty(self):
        self._results = {}
        for b in (self.btnReport, self.btnFig):
            b.setEnabled(False)

    # ------------------------------------------------------------ loading
    def _browse_cell(self, item):
        if item.column() == 0:
            return
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Choose the processed spectrum for this cell", "",
            "Spectra (*.csv *.txt 1r);;All files (*)")
        if paths:
            self._drop_files(item.row(), item.column(), paths)

    def _drop_files(self, row: int, col: int, paths: list[str]):
        if col == 0:                       # dropped on the name column
            col = 1
        errors = []
        for k, p in enumerate(paths):
            c = col + k
            if c >= self.table.columnCount():
                errors.append(f"{Path(p).name}: no field column left for it")
                continue
            try:
                self._load_cell(row, c, p)
            except Exception as exc:                          # noqa: BLE001
                errors.append(f"{Path(p).name}: {exc}")
        self._refresh_headers()
        self._set_dirty()
        if errors:
            QMessageBox.warning(self, "Some files were not loaded",
                                "\n".join(errors))

    def _load_cell(self, row: int, col: int, path: str):
        from larmor import qcpmg
        from larmor.loader import load_any

        ppm, amp, _recipe, meta, _warn = load_any(path)
        ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
        larmor = 0.0
        for key in ("larmor_MHz", "larmor_frequency_MHz"):
            if isinstance(meta, dict) and meta.get(key):
                larmor = float(meta[key]); break
        if not larmor:                     # a bare csv: read the header LARMOR
            from larmor.io import spectra
            try:
                _, _, m2 = spectra.read_csv(path)
                larmor = float(m2.get("larmor_MHz", 0.0) or 0.0)
                self._nucleus = self._nucleus or str(m2.get("nucleus", ""))
            except Exception:                                 # noqa: BLE001
                pass
        if not larmor:
            raise ValueError("no Larmor frequency in the file — save it from "
                             "the QCPMG dialog, which records one")
        hi, lo = qcpmg.cg_window(ppm, amp)
        if not (np.isfinite(hi) and np.isfinite(lo)) or hi <= lo:
            span = float(ppm.max() - ppm.min())
            mid = float(ppm.min()) + span / 2.0
            lo, hi = mid - span / 6.0, mid + span / 6.0
        self.cells[(row, col)] = {"path": path, "ppm": ppm, "amp": amp,
                                  "larmor": larmor, "window": (lo, hi)}
        self._measure_cell(row, col)
        if self._nucleus:
            self.lblNuc.setText(f"nucleus <b>{self._nucleus}</b>")
            self.spin.setValue(_spin_of(self._nucleus))

    def _measure_cell(self, row: int, col: int):
        from larmor import qcpmg
        d = self.cells.get((row, col))
        if d is None:
            return
        lo, hi = d["window"]
        cg, sigma = qcpmg.centre_of_gravity(d["ppm"], d["amp"], (hi, lo))
        fw_ppm = qcpmg.fwhm_hz(d["ppm"], d["amp"], 1.0, (hi, lo))
        d["cg"], d["sigma"], d["fwhm"] = cg, max(sigma, 0.1), fw_ppm
        it = self.table.item(row, col)
        if it is not None:
            txt = Path(d["path"]).name
            if np.isfinite(cg):
                txt += f"\nδcg {cg:.1f} ± {d['sigma']:.1f}   FWHM {fw_ppm:.0f} ppm"
            else:
                txt += "\n⚠ no usable signal in the window"
            self.table.blockSignals(True)
            it.setText(txt)
            it.setToolTip(d["path"])
            self.table.blockSignals(False)

    # -------------------------------------------------------- supervision
    def _show_cell(self):
        items = self.table.selectedItems()
        if not items:
            return
        r, c = items[0].row(), items[0].column()
        d = self.cells.get((r, c))
        self.plot.clear()
        self._region = None
        self._sel = (r, c)
        if d is None:
            return
        self.plot.plot(d["ppm"], d["amp"],
                       pen=pg.mkPen(theme.active().experiment, width=1.2))
        lo, hi = d["window"]
        self._region = pg.LinearRegionItem(values=(lo, hi), movable=True)
        self.plot.addItem(self._region)
        self._region.sigRegionChangeFinished.connect(self._region_moved)
        if d.get("cg") is not None and np.isfinite(d["cg"]):
            line = pg.InfiniteLine(pos=d["cg"], angle=90, movable=False,
                                   pen=pg.mkPen(theme.active().pivot,
                                                style=Qt.DashLine))
            self.plot.addItem(line)
        self.plot.setTitle(f"{Path(d['path']).name} — {d['larmor']:.3f} MHz",
                           color=theme.active().text_dim, size="9pt")

    def _region_moved(self):
        if self._region is None or self._sel is None:
            return
        d = self.cells.get(self._sel)
        if d is None:
            return
        a, b = self._region.getRegion()
        d["window"] = (min(a, b), max(a, b))
        self._measure_cell(*self._sel)
        self._set_dirty()
        self._show_cell()

    # ------------------------------------------------------------ compute
    def _rows(self):
        from larmor.qcpmg_fields import FieldPoint
        out = []
        for (r, c), d in sorted(self.cells.items()):
            if d.get("cg") is None or not np.isfinite(d["cg"]):
                continue
            name_item = self.table.item(r, 0)
            name = (name_item.text() if name_item else "") or f"sample {r + 1}"
            out.append((name, FieldPoint(d["larmor"], float(d["cg"]),
                                         float(d["sigma"]), None, name)))
        return out

    def _compute(self):
        from larmor.qcpmg_fields import (InfiniteFieldResult, fit_samples,
                                         report_text, two_field_widths)
        rows = self._rows()
        if len(rows) < 2:
            self.msg.setText("load at least two fields for one sample first")
            return
        self._results = fit_samples(rows, spin=self.spin.value(),
                                    eta=self.eta.value())
        widths = {}
        for name in self._results:
            fw = [(d["larmor"], d["fwhm"]) for (r, c), d in sorted(self.cells.items())
                  if (self.table.item(r, 0).text() if self.table.item(r, 0)
                      else f"sample {r + 1}") == name and d.get("fwhm")]
            if len(fw) >= 2:
                widths[name] = two_field_widths(fw[0][0], fw[0][1],
                                                fw[1][0], fw[1][1])
        self._widths = widths
        self.report.setPlainText(report_text(
            self._results, self.spin.value(), self.eta.value(),
            self._nucleus, widths))
        n_ok = sum(1 for r in self._results.values()
                   if isinstance(r, InfiniteFieldResult))
        self.msg.setText(f"{n_ok} of {len(self._results)} samples extrapolated")
        for b in (self.btnReport, self.btnFig):
            b.setEnabled(n_ok > 0)

    # ------------------------------------------------------------- export
    def _figure_spec(self, only: str | None = None) -> dict:
        from larmor.qcpmg_fields import InfiniteFieldResult
        samples = []
        for name, res in self._results.items():
            if not isinstance(res, InfiniteFieldResult):
                continue
            if only is not None and name != only:
                continue
            samples.append({"label": name,
                            "points": [[p.larmor_MHz, p.dcg_ppm, p.dcg_err_ppm]
                                       for p in res.points]})
        return {"kind": "infinite_field", "style": "article",
                "nucleus": self._nucleus, "spin": self.spin.value(),
                "eta": self.eta.value(), "samples": samples}

    def _export_report(self):
        from larmor.desktop.paths import (FIGURE_DIR_KEY, remember_dir,
                                          remembered_dir)
        start = remembered_dir(FIGURE_DIR_KEY)
        seed = str(Path(start) / "infinite_field_report.txt") if start \
            else "infinite_field_report.txt"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export report", seed, "Text (*.txt);;All files (*)")
        if not path:
            return
        remember_dir(FIGURE_DIR_KEY, path)
        Path(path).write_text(self.report.toPlainText(), encoding="utf-8")
        self.msg.setText(f"report written — {Path(path).name}")

    def _export_figures(self):
        from larmor import figures
        from larmor.desktop.paths import (FIGURE_DIR_KEY, remember_dir,
                                          remembered_dir)
        from larmor.qcpmg_fields import InfiniteFieldResult
        start = remembered_dir(FIGURE_DIR_KEY)
        seed = str(Path(start) / "infinite_field") if start else "infinite_field"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export figures (base name — extensions added)", seed,
            "Figure base name (*)")
        if not path:
            return
        remember_dir(FIGURE_DIR_KEY, path)
        base = Path(path).with_suffix("")
        written = []
        try:
            written += figures.export(self._figure_spec(), base)
            for name, res in self._results.items():
                if not isinstance(res, InfiniteFieldResult):
                    continue
                safe = "".join(ch if ch.isalnum() or ch in "-_" else "_"
                               for ch in name)[:40] or "sample"
                written += figures.export(self._figure_spec(only=name),
                                          f"{base}_{safe}")
        except Exception as exc:                              # noqa: BLE001
            _log.exception("infinite-field figure export failed")
            self.msg.setText(f"figure export failed: {exc}")
            return
        self.msg.setText(f"{len(written)} files written — {base.name}.* "
                         f"(merged) and one set per sample")

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "qcpmg", "QCPMG")


def _spin_of(nucleus: str) -> float:
    try:
        from mrsimulator.spin_system.isotope import ISOTOPE_DATA
        d = ISOTOPE_DATA.get(nucleus or "")
        if d:
            return (d["spin_multiplicity"] - 1) / 2.0
    except Exception:                                         # noqa: BLE001
        pass
    return 1.5
