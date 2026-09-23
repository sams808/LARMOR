"""Session inventory window (Tools > Session inventory…): one month folder
read into a sample x nucleus grid with the production EXPNO pre-picked per
block, the roles / reasons / title-vs-folder flags of every EXPNO in a
detail table, a manual pick override, and one-action hand-off of the picks
to Batch fit or Sequential fit. Logic lives in larmor.inventory (Qt-free,
tested); the window renders it and emits signals the main window's existing
slots accept."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMenu, QMessageBox, QPushButton,
    QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout,
)

from larmor import inventory as I
from larmor import referencing as R
from larmor.desktop import theme
from larmor.desktop.referencing_dialog import _VERDICT_COLOR

_ROLE_COLOR = {"production": "#2C6A4E", "failed": "#9E3324", "short": "#6A7C80",
               "setup": "#6A7C80", "2D": "#6A7C80", "arrayed": "#6A7C80",
               "unprocessed": "#6A7C80", "reference": "#0B6A71"}
_BG_OK = QColor(44, 106, 78, 45)
_BG_WARN = QColor(168, 87, 15, 55)
_BG_BAD = QColor(158, 51, 36, 60)
_AMBER = "#A8570F"
DETAIL_COLS = ["pick", "EXPNO", "nucleus", "role", "kind", "NS", "D1 (s)", "1r",
               "procs", "date", "SR", "title", "flags"]
_ROLE_KEY = Qt.UserRole + 1          # (folder, key, nucleus, expno) on detail items


def _pickable(r: I.Row) -> bool:
    return r.info.ndim == 1 and r.info.has_1r and r.nucleus != "1H"


class SessionInventoryDialog(QDialog):
    #: an openable data path -> MainWindow.load_source
    open_requested = Signal(str)
    #: the picks' openable paths, in sample order -> MainWindow.run_batch_fit
    batch_requested = Signal(list)
    #: -> MainWindow.run_seq_fit
    seq_requested = Signal(list)

    def __init__(self, parent=None, start_dir: str | None = None):
        super().__init__(parent)
        self.setWindowTitle("Session inventory — a month folder as a sample × nucleus grid")
        self.setModal(False)
        self.resize(1180, 760)
        self.inv: I.Inventory | None = None
        self._filling = False
        self._sample = None                  # the sample shown in the detail table
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
        self.chkSr = QCheckBox("SR audit (needs a referenced ¹H in the session)")
        self.chkSr.setChecked(True)
        self.chkSr.setToolTip("second pass through the referencing audit: the SR column "
                              "and a red cell where a pick is unreferenced or off")
        opts.addWidget(self.chkSr)
        opts.addSpacing(16)
        opts.addWidget(QLabel("a pick whose NS is below"))
        self.nsFrac = QDoubleSpinBox(); self.nsFrac.setDecimals(2)
        self.nsFrac.setRange(0.01, 0.9); self.nsFrac.setSingleStep(0.05)
        self.nsFrac.setValue(I.DEFAULT_NS_FRAC)
        self.nsFrac.setToolTip(
            "the production spectrum is the highest EXPNO of the block with a "
            "pdata/1/1r; when its NS is below this fraction of the block's maximum it "
            "is a short shot (the last quick check) and the pick falls to the next "
            "EXPNO. Changing the value re-picks the whole grid (manual picks reset).")
        opts.addWidget(self.nsFrac)
        opts.addWidget(QLabel("× the block's maximum is a short shot"))
        opts.addStretch(1)
        opts.addWidget(QLabel("hand-off nucleus"))
        self.cmbNucleus = QComboBox()
        self.cmbNucleus.setToolTip("the picks of this nucleus go to Batch fit / Sequential "
                                   "fit; clicking a grid column header selects it")
        self.cmbNucleus.currentTextChanged.connect(lambda _t: self._enable(self.inv is not None))
        opts.addWidget(self.cmbNucleus)
        v.addLayout(opts)
        self.nsFrac.valueChanged.connect(self._repick)

        self.status = QLabel("pick the session folder and press Scan")
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color: {t.text_dim};")
        v.addWidget(self.status)

        split = QSplitter(Qt.Vertical)
        self.grid = QTableWidget(0, 0)
        self.grid.setSelectionBehavior(QTableWidget.SelectRows)
        self.grid.setSelectionMode(QTableWidget.SingleSelection)
        self.grid.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.grid.horizontalHeader().setStretchLastSection(True)
        self.grid.horizontalHeader().sectionClicked.connect(self._column_clicked)
        self.grid.itemSelectionChanged.connect(self._grid_selected)
        self.grid.cellDoubleClicked.connect(self._grid_double)
        self.grid.setContextMenuPolicy(Qt.CustomContextMenu)
        self.grid.customContextMenuRequested.connect(self._context_menu)
        split.addWidget(self.grid)

        self.detail = QTableWidget(0, len(DETAIL_COLS))
        self.detail.setHorizontalHeaderLabels(DETAIL_COLS)
        self.detail.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.detail.horizontalHeader().setStretchLastSection(True)
        self.detail.setSelectionBehavior(QTableWidget.SelectRows)
        self.detail.itemChanged.connect(self._pick_toggled)
        self.detail.itemDoubleClicked.connect(self._detail_double)
        split.addWidget(self.detail)
        split.setSizes([380, 320])
        v.addWidget(split, 1)

        bottom = QHBoxLayout()
        self.btnBatch = QPushButton("Batch fit picks…")
        self.btnBatch.setToolTip("send the production pick of every sample, in this "
                                 "order, to Batch fit — panels are named by sample")
        self.btnBatch.clicked.connect(self._batch)
        self.btnSeq = QPushButton("Sequential fit picks…")
        self.btnSeq.setToolTip("the same picks to the forward–backward sequential fit")
        self.btnSeq.clicked.connect(self._seq)
        self.btnOpen = QPushButton("Open pick")
        self.btnOpen.setToolTip("open the selected EXPNO (detail table) or the selected "
                                "sample's pick in the workbench")
        self.btnOpen.clicked.connect(self._open)
        self.btnCopy = QPushButton("Copy picks")
        self.btnCopy.setToolTip("sample, nucleus, EXPNO and path of every pick, tab-separated")
        self.btnCopy.clicked.connect(self._copy)
        self.btnExport = QPushButton("Export CSV…")
        self.btnExport.setToolTip("inventory_<month>.csv (every EXPNO) and "
                                  "inventory_<month>_picks.txt")
        self.btnExport.clicked.connect(self.export)
        for w in (self.btnBatch, self.btnSeq, self.btnOpen, self.btnCopy, self.btnExport):
            bottom.addWidget(w)
        bottom.addStretch(1)
        btnClose = QPushButton("Close"); btnClose.clicked.connect(self.close)
        bottom.addWidget(btnClose)
        v.addLayout(bottom)
        for b in self.findChildren(QPushButton):
            b.setAutoDefault(False)
        self._enable(False)

    # ------------------------------------------------------------------ state
    def _enable(self, on: bool):
        for w in (self.btnCopy, self.btnExport, self.btnOpen):
            w.setEnabled(on)
        n = len(self.current_picks()) if on else 0
        self.btnBatch.setEnabled(n >= 2)
        self.btnSeq.setEnabled(n >= 2)

    def current_nucleus(self) -> str:
        return self.cmbNucleus.currentText()

    def current_picks(self) -> list[str]:
        """The openable paths of the hand-off nucleus' picks, in sample order."""
        if self.inv is None:
            return []
        return [r.openable for r in self.inv.picks(self.current_nucleus() or None)
                if r.openable]

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Session folder (one month)",
                                             self.folder.text() or str(Path.home()))
        if d:
            self.folder.setText(d)

    # ------------------------------------------------------------------ scan
    def scan(self):
        folder = self.folder.text().strip()
        if not folder or not Path(folder).is_dir():
            self.status.setText("that folder does not exist")
            return
        root = I.inventory_root(folder)
        self.folder.setText(str(root))
        self._sample = None

        def prog(k, n, name):
            self.status.setText(f"scanning {root.name} · {k + 1}/{n} · {name}")
            QApplication.processEvents()

        self.status.setText(f"scanning {root} …")
        QApplication.processEvents()
        self.inv = I.build(root, ns_frac=float(self.nsFrac.value()), sr_audit=False,
                           progress=prog)
        self._fill_all()
        if self.chkSr.isChecked() and self.inv.rows:
            self.status.setText(f"{self._summary()}  ·  referencing audit …")
            QApplication.processEvents()

            def prog_sr(k, n, name):
                self.status.setText(f"{self._summary()}  ·  referencing audit {k + 1}/{n} {name}")
                QApplication.processEvents()

            acqs = R.scan_session(root, progress=prog_sr)
            if acqs:
                I.join_sr(self.inv, R.audit(acqs))
            self._fill_all()
        self.status.setText(self._summary())
        self.setWindowTitle(f"Session inventory — {root.name}")

    def _summary(self) -> str:
        inv = self.inv
        if inv is None or not inv.rows:
            return "no EXPNO found — pick the month folder that holds the sample folders"
        return (f"{len(inv.rows)} EXPNOs in {inv.n_folders()} sample folders · "
                f"{' '.join(inv.nuclei())} · {len(inv.picks())} picks · "
                f"{inv.n_flags()} flagged")

    def _repick(self, value: float):
        """A new NS fraction re-runs the roles and flags (no rescan)."""
        if self.inv is None:
            return
        self.inv.ns_frac = float(value)
        I.assign_roles(self.inv.rows, float(value))
        I.flag_titles(self.inv)
        self._fill_all()
        self.status.setText(self._summary() + "  ·  re-picked with the new NS fraction")

    def _fill_all(self):
        inv = self.inv
        keep = self.current_nucleus()
        self.cmbNucleus.blockSignals(True)
        self.cmbNucleus.clear()
        nucs = inv.nuclei() if inv else []
        self.cmbNucleus.addItems(nucs)
        if keep in nucs:
            self.cmbNucleus.setCurrentText(keep)
        elif nucs:
            self.cmbNucleus.setCurrentText(next((n for n in nucs if n != "1H"), nucs[0]))
        self.cmbNucleus.blockSignals(False)
        self._fill_grid()
        self._fill_detail(self._sample)
        self._enable(inv is not None and bool(inv.rows))

    # ------------------------------------------------------------------ grid
    @staticmethod
    def _cell_text(r: I.Row | None, n_block: int) -> str:
        if r is None:
            return "—"
        if r.role == "reference":
            return f"{r.expno} · ¹H reference"
        txt = f"{r.expno} · NS {r.ns}"
        if r.d1_s:
            txt += f" · D1 {r.d1_s:g} s"
        if n_block > 1:
            txt += f"  (+{n_block - 1})"
        return txt

    def _fill_grid(self):
        inv = self.inv
        t = theme.active()
        self.grid.blockSignals(True)
        if inv is None:
            self.grid.setRowCount(0); self.grid.setColumnCount(0)
            self.grid.blockSignals(False)
            return
        samples, nucs, labels = inv.samples(), inv.nuclei(), inv.labels()
        blocks = inv.blocks()
        self.grid.setRowCount(len(samples))
        self.grid.setColumnCount(len(nucs))
        self.grid.setHorizontalHeaderLabels(nucs)
        for i, s in enumerate(samples):
            rows_of = [r for r in inv.rows if r.sample_id == s]
            first = rows_of[0] if rows_of else None
            n2d = sum(1 for r in rows_of if r.info.ndim == 2)
            hdr = QTableWidgetItem(labels[s])
            tip = [f"folder: {s[0]}"]
            if first is not None:
                if first.name.rotor:
                    tip.append(f"rotor: {first.name.rotor}")
                if first.date_iso:
                    tip.append(f"date: {first.date_iso}")
                if first.title_line:
                    tip.append(f"title: {first.title_line}")
                tip.append(f"sample name from: {first.name.source}")
            if n2d:
                tip.append(f"+{n2d} 2D/arrayed")
            hdr.setToolTip("\n".join(tip))
            self.grid.setVerticalHeaderItem(i, hdr)
            for j, n in enumerate(nucs):
                r = inv.cell(s, n)
                n_block = len(blocks.get((s[0], s[1], n), []))
                it = QTableWidgetItem(self._cell_text(r, n_block))
                it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                it.setData(Qt.UserRole, (s[0], s[1], n))
                if r is None:
                    it.setForeground(QColor(t.text_dim))
                    it.setToolTip("no 1D production spectrum" if n_block == 0 else
                                  "no pick in this block (every 1D EXPNO is 2D, "
                                  "unprocessed or a reference)")
                elif r.role == "reference":
                    it.setForeground(QColor(_VERDICT_COLOR["1H reference"]))
                    it.setToolTip("the ¹H reference spectrum — never picked")
                else:
                    if r.sr_verdict in ("unreferenced", "off"):
                        it.setBackground(_BG_BAD)
                    elif r.reasons or r.flags:
                        it.setBackground(_BG_WARN)
                    else:
                        it.setBackground(_BG_OK)
                    tip = list(r.reasons) + list(r.flags)
                    if r.sr_verdict:
                        tip.append(f"SR: {r.sr_verdict}" + (f" — {r.sr_note}" if r.sr_note else ""))
                    it.setToolTip("\n".join(tip) or "clean pick: the highest EXPNO of the "
                                  "block with a pdata/1/1r")
                self.grid.setItem(i, j, it)
        self.grid.blockSignals(False)

    def _column_clicked(self, col: int):
        if self.inv is None:
            return
        nucs = self.inv.nuclei()
        if 0 <= col < len(nucs):
            self.cmbNucleus.setCurrentText(nucs[col])

    def _grid_selected(self):
        it = self.grid.currentItem()
        if it is None or self.inv is None:
            return
        folder, key, _n = it.data(Qt.UserRole)
        self._sample = (folder, key)
        self._fill_detail(self._sample)

    def _grid_double(self, row: int, col: int):
        it = self.grid.item(row, col)
        if it is None or self.inv is None:
            return
        folder, key, n = it.data(Qt.UserRole)
        r = self.inv.cell((folder, key), n)
        if r is not None and r.openable:
            self.open_requested.emit(r.openable)

    # ------------------------------------------------------------------ detail
    def _fill_detail(self, sample=None):
        inv = self.inv
        self._filling = True
        try:
            if inv is None:
                self.detail.setRowCount(0)
                return
            rows = [r for r in inv.rows if sample is None or r.sample_id == sample]
            rows.sort(key=lambda r: (r.folder.lower(), r.expno, r.info.expno))
            self.detail.setRowCount(len(rows))
            for i, r in enumerate(rows):
                info = r.info
                vals = ["", str(r.expno), r.nucleus, r.role, info.kind, str(r.ns),
                        f"{r.d1_s:g}" if r.d1_s else "",
                        "✓" if info.has_1r else ("2rr" if info.has_2rr else "—"),
                        str(info.n_procs), r.date_iso, r.sr_verdict, r.title_line,
                        "; ".join(r.flags)]
                tip = "\n".join(list(r.reasons) + list(r.flags)
                                + ([f"SR: {r.sr_note}"] if r.sr_note else []))
                for j, txt in enumerate(vals):
                    it = QTableWidgetItem(txt)
                    flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable
                    if j == 0 and _pickable(r):
                        flags |= Qt.ItemIsUserCheckable
                        it.setCheckState(Qt.Checked if r.pick else Qt.Unchecked)
                    it.setFlags(flags)
                    it.setData(Qt.UserRole, r.openable or r.path)
                    it.setData(_ROLE_KEY, (r.folder, r.sample, r.nucleus, r.expno))
                    if tip:
                        it.setToolTip(tip)
                    if j == 3 and r.role in _ROLE_COLOR:
                        it.setForeground(QColor(_ROLE_COLOR[r.role]))
                    elif j == 10 and r.sr_verdict:
                        it.setForeground(QColor(_VERDICT_COLOR.get(r.sr_verdict, "#6A7C80")))
                    elif j == 12 and r.flags:
                        it.setForeground(QColor(_AMBER))
                    if r.pick:
                        f = it.font(); f.setBold(True); it.setFont(f)
                    self.detail.setItem(i, j, it)
        finally:
            self._filling = False

    def _pick_toggled(self, item: QTableWidgetItem):
        if self._filling or item.column() != 0 or self.inv is None:
            return
        folder, key, nuc, expno = item.data(_ROLE_KEY)
        block = self.inv.candidates((folder, key), nuc)
        cur = next((r for r in block if r.pick), None)
        if item.checkState() == Qt.Checked:
            self.inv.set_pick((folder, key), nuc, expno)
            msg = f"{key} · {nuc}: EXPNO {expno} chosen by hand"
        elif cur is not None and cur.expno == expno:
            self.inv.set_pick((folder, key), nuc, None)
            msg = f"{key} · {nuc}: no pick — the sample is left out of the hand-off"
        else:
            return
        self._sample = (folder, key)
        self._fill_grid()
        self._fill_detail(self._sample)
        self._enable(True)
        self.status.setText(f"{self._summary()}  ·  {msg}")

    def _detail_double(self, item: QTableWidgetItem):
        path = item.data(Qt.UserRole)
        if path:
            self.open_requested.emit(str(path))

    # ------------------------------------------------------------------ actions
    def _batch(self):
        paths = self.current_picks()
        if len(paths) < 2:
            self.status.setText(f"fewer than two {self.current_nucleus()} picks — nothing "
                                "to batch fit")
            return
        self.batch_requested.emit(paths)
        self.status.setText(f"{len(paths)} {self.current_nucleus()} picks sent to Batch fit")

    def _seq(self):
        paths = self.current_picks()
        if len(paths) < 2:
            self.status.setText(f"fewer than two {self.current_nucleus()} picks — nothing "
                                "to fit sequentially")
            return
        self.seq_requested.emit(paths)
        self.status.setText(f"{len(paths)} {self.current_nucleus()} picks sent to "
                            "Sequential fit")

    def _open(self):
        it = self.detail.currentItem()
        if it is not None and it.data(Qt.UserRole):
            self.open_requested.emit(str(it.data(Qt.UserRole)))
            return
        g = self.grid.currentItem()
        if g is not None and self.inv is not None:
            folder, key, n = g.data(Qt.UserRole)
            r = self.inv.cell((folder, key), n)
            if r is not None and r.openable:
                self.open_requested.emit(r.openable)
                return
        self.status.setText("select an EXPNO in the detail table (or a grid cell) first")

    def _copy(self):
        if self.inv is None:
            return
        QApplication.clipboard().setText(I.picks_text(self.inv, self.current_nucleus() or None))
        self.status.setText(f"{len(self.current_picks())} {self.current_nucleus()} picks "
                            "copied to the clipboard (sample, nucleus, EXPNO, path)")

    def export(self):
        if self.inv is None or not self.inv.rows:
            return
        root = self.inv.root
        default = str(Path.home() / f"inventory_{root.name}.csv")
        path, _ = QFileDialog.getSaveFileName(self, "Export the inventory", default,
                                              "CSV (*.csv)")
        if not path:
            return
        p = Path(path)
        try:
            I.to_csv(self.inv, p)
            txt = p.with_name(p.stem + "_picks.txt")
            txt.write_text(I.picks_text(self.inv, self.current_nucleus() or None),
                           encoding="utf-8")
        except OSError as exc:
            QMessageBox.warning(self, "Export", f"could not write: {exc}")
            return
        self.status.setText(f"wrote {p.name} and {txt.name}")

    def _context_menu(self, pos):
        if self.inv is None:
            return
        it = self.grid.itemAt(pos)
        if it is not None:
            self.grid.setCurrentItem(it)
        m = QMenu(self)
        m.addAction("Open pick", self._open)
        m.addAction("Batch fit picks…", self._batch)
        m.addAction("Sequential fit picks…", self._seq)
        m.addSeparator()
        m.addAction("Copy picks", self._copy)
        m.addAction("Export CSV…", self.export)
        m.exec(self.grid.viewport().mapToGlobal(pos))
