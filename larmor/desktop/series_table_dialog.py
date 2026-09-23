"""Series table dialog: the one editor Batch fit and Sequential fit open from
"Series table…".

Rows are the spectra of the series, in the order the fit uses (the grid, the
results table, the CSV rows, the saved-recipe names and -- in the sequential
fit -- the sweep). Columns: ``#`` | ``name`` (editable; names stay unique) |
``group`` (editable; replicates share a group and are what "average
replicates" collapses) | ``folder`` | ``title`` (read-only, from the
loader) | one editable column per numeric metadata column (typed with
**Add column…** or joined from a CSV with **Join CSV…** -- an EPMA
composition table or another batch's long ``batch_table*.csv``, pivoted).

Front-end only: every rule lives in :mod:`larmor.series_table`. Reordering
uses the up / down arrows and **Sort by…** (natural name order or any numeric
column), as the Plotting studio does -- row drag-and-drop in QTableWidget
moves cells, not rows.
"""
from __future__ import annotations

import copy
import csv
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDialog, QDialogButtonBox, QFileDialog,
    QHBoxLayout, QInputDialog, QLabel, QListWidget, QListWidgetItem, QPushButton,
    QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout,
)

from larmor import series_table as stab
from larmor.desktop.comparability_dialog import _CHECK_CSS_COLOR
from larmor.desktop.paths import suggest_save_dir
from larmor.series_table import Match, SeriesColumn, SeriesTable

#: the fixed columns before the metadata ones
_FIXED = ("#", "name", "group", "folder", "title")
_COL_NAME, _COL_GROUP = 1, 2


def _ro(text: str, tip: str = "") -> QTableWidgetItem:
    it = QTableWidgetItem(str(text))
    it.setFlags(it.flags() & ~Qt.ItemIsEditable)
    if tip:
        it.setToolTip(tip)
    return it


class SeriesTableDialog(QDialog):
    """Edit a copy of ``table``; ``table()`` / ``perm()`` hand back the result
    in the visual order (``perm[new_pos] = original index``). ``locked``
    greys the move / sort buttons (a fit or error worker is running -- the
    dialogs rebuild their grid on reorder)."""

    def __init__(self, parent, table: SeriesTable, locked: bool = False):
        super().__init__(parent)
        self.setWindowTitle("Series table — names, groups, order, columns")
        self.resize(940, 560)
        self._orig = copy.deepcopy(table)
        self._table = copy.deepcopy(table)
        self._order = list(range(len(self._table.rows)))
        self._colmap: list = []
        self._join_note = ""

        v = QVBoxLayout(self)
        intro = QLabel(
            "Rows are the series in fit order — move them with the arrows or "
            "<b>Sort by…</b>. Double-click a <b>name</b> or <b>group</b> to edit "
            "(names stay unique; replicates share a group). <b>Join CSV…</b> adds "
            "composition columns from an EPMA table or another batch's CSV.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        self.grid = QTableWidget()
        self.grid.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.grid.setSelectionMode(QAbstractItemView.SingleSelection)
        self.grid.verticalHeader().setVisible(False)
        self.grid.itemChanged.connect(self._on_item_changed)
        v.addWidget(self.grid, 1)

        bar = QHBoxLayout()
        self.btnUp = QToolButton(); self.btnUp.setArrowType(Qt.UpArrow)
        self.btnUp.setToolTip("move the selected spectrum up the series")
        self.btnUp.clicked.connect(lambda: self._move(-1))
        self.btnDown = QToolButton(); self.btnDown.setArrowType(Qt.DownArrow)
        self.btnDown.setToolTip("move the selected spectrum down the series")
        self.btnDown.clicked.connect(lambda: self._move(+1))
        self.btnSort = QPushButton("Sort by…")
        self.btnSort.setToolTip("order the series by name (natural: 0Ca, 1Ca, … 10Ca) "
                                "or by any numeric column")
        self.btnSort.clicked.connect(self._ask_sort)
        self.btnReset = QPushButton("Reset names")
        self.btnReset.setToolTip("back to the loader's names and groups (sample folder "
                                 "without its date / rotor / operator tokens)")
        self.btnReset.clicked.connect(self._reset_names)
        self.btnAddCol = QPushButton("Add column…")
        self.btnAddCol.setToolTip("a numeric column to type in (e.g. CaO (mol%))")
        self.btnAddCol.clicked.connect(lambda: self._add_column())
        self.btnRemoveCol = QPushButton("Remove column…")
        self.btnRemoveCol.clicked.connect(lambda: self._remove_column())
        self.btnJoin = QPushButton("Join CSV…")
        self.btnJoin.setToolTip("join columns from a wide CSV (one row per sample, e.g. "
                                "EPMA compositions) or from another batch's long "
                                "batch_table*.csv (pivoted per sample); the row mapping "
                                "is shown and confirmed first")
        self.btnJoin.clicked.connect(lambda: self._join_csv())
        self.btnSave = QPushButton("Save table…")
        self.btnSave.setToolTip("write this table as <name>.series.json")
        self.btnSave.clicked.connect(lambda: self._save_table())
        self.btnLoad = QPushButton("Load table…")
        self.btnLoad.setToolTip("names, groups, order and columns from a saved "
                                ".series.json (rows pair by source path, else by name)")
        self.btnLoad.clicked.connect(lambda: self._load_table())
        for b in (self.btnUp, self.btnDown, self.btnSort, self.btnReset, self.btnAddCol,
                  self.btnRemoveCol, self.btnJoin, self.btnSave, self.btnLoad):
            bar.addWidget(b)
        bar.addStretch(1)
        v.addLayout(bar)
        if locked:
            for b in (self.btnUp, self.btnDown, self.btnSort):
                b.setEnabled(False)
                b.setToolTip("reordering waits for the running fit / error analysis")

        self.status = QLabel("")
        self.status.setWordWrap(True)
        v.addWidget(self.status)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
        for b in self.findChildren(QPushButton):
            b.setAutoDefault(False)
        self._fill()

    # ------------------------------------------------------------------ result
    def table(self) -> SeriesTable:
        """The edited table in the visual order."""
        return self._table.permuted(self._order)

    def perm(self) -> list:
        """new position -> original index (identity when nothing moved)."""
        return list(self._order)

    # ------------------------------------------------------------------ grid
    def _fill(self):
        g = self.grid
        g.blockSignals(True)
        headers = list(_FIXED)
        self._colmap = []
        for c in self._table.columns:
            headers.append(c.label); self._colmap.append((c, c.key))
            if c.err_key:
                headers.append(f"{c.label} ±"); self._colmap.append((c, c.err_key))
        g.clear()
        g.setColumnCount(len(headers))
        g.setHorizontalHeaderLabels(headers)
        g.setRowCount(len(self._order))
        for j, (c, key) in enumerate(self._colmap):
            tip = f"key: {key}"
            if c.err_key:
                tip += f" · ± {c.err_key}"
            if c.source:
                tip += f" · from {c.source}"
            it = g.horizontalHeaderItem(len(_FIXED) + j)
            if it is not None:
                it.setToolTip(tip)
        for pos, k in enumerate(self._order):
            r = self._table.rows[k]
            num = _ro(str(pos + 1), r.source_path)
            num.setData(Qt.UserRole, k)
            g.setItem(pos, 0, num)
            name = QTableWidgetItem(r.display_name)
            name.setToolTip("the label of the panel, the table row, the CSV scope and "
                            "the saved recipe — kept unique")
            g.setItem(pos, _COL_NAME, name)
            grp = QTableWidgetItem(r.group)
            grp.setToolTip("replicates of one composition share a group; 'average "
                           "replicates' in the Series plot collapses a group")
            g.setItem(pos, _COL_GROUP, grp)
            g.setItem(pos, 3, _ro(r.folder, r.source_path))
            g.setItem(pos, 4, _ro(r.title, r.title))
            for j, (c, key) in enumerate(self._colmap):
                val = r.values.get(key)
                it = QTableWidgetItem("" if val is None else f"{val:.6g}")
                it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                g.setItem(pos, len(_FIXED) + j, it)
        g.resizeColumnsToContents()
        g.blockSignals(False)
        self._update_status()

    def _on_item_changed(self, item: QTableWidgetItem):
        pos, c = item.row(), item.column()
        if pos < 0 or pos >= len(self._order):
            return
        k = self._order[pos]
        row = self._table.rows[k]
        self.grid.blockSignals(True)
        try:
            if c == _COL_NAME:
                applied = self._table.rename(k, item.text())
                item.setText(applied)
                self.grid.item(pos, _COL_GROUP).setText(row.group)
            elif c == _COL_GROUP:
                self._table.set_group(k, item.text())
                item.setText(row.group)
            elif c >= len(_FIXED) and c - len(_FIXED) < len(self._colmap):
                _col, key = self._colmap[c - len(_FIXED)]
                val = stab.as_float(item.text())
                row.values[key] = val
                item.setText("" if val is None else f"{val:.6g}")
        finally:
            self.grid.blockSignals(False)
        self._update_status()

    def _update_status(self):
        n = len(self._table.rows)
        groups = self._table.groups()
        reps = [g for g, idx in groups.items() if len(idx) > 1]
        parts = [f"{n} spectra",
                 "order: " + ("as loaded" if self._order == list(range(n)) else "custom")]
        if reps:
            parts.append(f"{len(reps)} replicate group(s): "
                         + ", ".join(f"{g} ×{len(groups[g])}" for g in reps))
        if self._table.columns:
            parts.append(f"{len(self._table.columns)} column(s)")
        if self._join_note:
            parts.append(self._join_note)
        self.status.setText(" · ".join(parts))

    # ------------------------------------------------------------------ order
    def _move(self, delta: int):
        pos = self.grid.currentRow()
        if pos < 0:
            sel = self.grid.selectionModel().selectedRows()
            pos = sel[0].row() if sel else -1
        new = pos + delta
        if pos < 0 or not 0 <= new < len(self._order):
            return
        self._order[pos], self._order[new] = self._order[new], self._order[pos]
        self._fill()
        self.grid.selectRow(new)
        self.grid.setCurrentCell(new, _COL_NAME)

    def _sort_by(self, key=None, descending: bool = False):
        perm = self.table().sort_perm(key, descending)
        self._order = [self._order[p] for p in perm]
        self._fill()

    def _ask_sort(self):
        cols = self._table.columns
        choices = ["name (natural)"] + [c.label for c in cols]
        choice, ok = QInputDialog.getItem(self, "Sort by", "order the series by:",
                                          choices, 0, False)
        if not ok:
            return
        direction, ok = QInputDialog.getItem(self, "Sort by", "direction:",
                                             ["ascending", "descending"], 0, False)
        if not ok:
            return
        key = None if choice == choices[0] else cols[choices.index(choice) - 1].key
        self._sort_by(key, direction == "descending")

    def _reset_names(self):
        for r, o in zip(self._table.rows, self._orig.rows):
            r.display_name, r.group = o.display_name, o.group
        self._fill()

    # ------------------------------------------------------------------ columns
    def _add_column(self, label: str | None = None):
        if label is None:
            label, ok = QInputDialog.getText(self, "Add column",
                                             "column label (e.g. CaO (mol%)):")
            if not ok:
                return
        label = (label or "").strip()
        if not label:
            return
        self._table.add_column(SeriesColumn(key=label, label=label))
        self._fill()

    def _remove_column(self, key: str | None = None):
        cols = self._table.columns
        if not cols:
            self.status.setText("no columns to remove")
            return
        if key is None:
            choice, ok = QInputDialog.getItem(self, "Remove column", "column:",
                                              [c.label for c in cols], 0, False)
            if not ok:
                return
            key = cols[[c.label for c in cols].index(choice)].key
        self._table.remove_column(key)
        self._fill()

    def _start_dir(self) -> str:
        return suggest_save_dir(self._table.rows[0].source_path if self._table.rows else None)

    def _join_csv(self, path: str | None = None):
        if path is None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Join a CSV — compositions (one row per sample) or a batch table",
                self._start_dir(), "CSV (*.csv *.tsv *.txt);;All (*)")
        if not path:
            return
        try:
            headers, rows, sample = stab.read_wide_csv(path)
            if stab.is_long_batch_csv(headers):
                headers, rows, sample = stab.read_batch_csv_wide(path)
        except (OSError, UnicodeDecodeError, csv.Error) as exc:
            self.status.setText(f"could not read {Path(path).name}: {exc}")
            return
        if not rows or not headers:
            self.status.setText(f"{Path(path).name} has no data rows")
            return
        matches = stab.propose_mapping(self._table, rows, sample)
        dlg = JoinMappingDialog(self, self._table, headers, rows, sample, matches,
                                Path(path).name, order=self._order)
        if dlg.exec() != QDialog.Accepted:
            self.status.setText("join cancelled — nothing changed")
            return
        cols = dlg.chosen_columns()
        if not cols:
            self.status.setText("no column ticked — nothing joined")
            return
        stab.apply_join(self._table, rows, sample, dlg.matches(), cols,
                        tag=dlg.tag(), source=Path(path).name)
        n_match = sum(1 for m in dlg.matches() if m.csv_row is not None)
        self._join_note = (f"matched {n_match}/{len(self._table.rows)} from "
                           f"{Path(path).name}" + (f" ({dlg.tag()})" if dlg.tag() else ""))
        self._fill()

    # ------------------------------------------------------------------ files
    def _save_table(self, path: str | None = None):
        if path is None:
            start = self._start_dir()
            path, _ = QFileDialog.getSaveFileName(
                self, "Save series table",
                str(Path(start) / "series.series.json") if start else "series.series.json",
                "Series table (*.series.json *.json)")
        if not path:
            return
        try:
            self.table().save(path)
        except OSError as exc:
            self.status.setText(f"could not save: {exc}")
            return
        self.status.setText(f"saved {Path(path).name}")

    def _load_table(self, path: str | None = None):
        if path is None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Load series table", self._start_dir(),
                "Series table (*.series.json *.json);;All (*)")
        if not path:
            return
        try:
            loaded = SeriesTable.load(path)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.status.setText(f"could not load {Path(path).name}: {exc}")
            return
        pairs = self._table.merge_from(loaded)
        n = len(self._table.rows)
        # follow the saved order for the rows that paired; the rest keep theirs
        self._order = sorted(range(n),
                             key=lambda k: (pairs[k] if pairs[k] is not None else n + k))
        matched = sum(1 for p in pairs if p is not None)
        self._join_note = f"loaded {matched}/{n} rows from {Path(path).name}"
        self._fill()


class JoinMappingDialog(QDialog):
    """Confirm a CSV join: one combo per series row (the proposed CSV row
    preselected; unmatched rows amber with ``(none)``; an ambiguous row blank
    with its candidates in the tooltip), a checklist of the numeric columns
    (oxide columns ticked by default, else all; the ± partner named beside
    each) and the analysed / nominal / none tag."""

    def __init__(self, parent, table: SeriesTable, headers, csv_rows, sample_header,
                 matches, source_name: str, order=None):
        super().__init__(parent)
        self.setWindowTitle(f"Join {source_name}")
        self.resize(760, 620)
        self._matches = [Match(m.row, m.csv_row, m.how, list(m.candidates)) for m in matches]
        self._csv_names = [(r.get(sample_header) or "").strip() for r in csv_rows]
        order = list(order) if order is not None else list(range(len(table.rows)))
        self._pos_of = {k: pos for pos, k in enumerate(order)}
        n_match = sum(1 for m in self._matches if m.csv_row is not None)

        v = QVBoxLayout(self)
        head = QLabel(
            f"<b>{source_name}</b>: {len(csv_rows)} rows, sample column "
            f"<b>{sample_header}</b> · proposed {n_match}/{len(table.rows)} matches — "
            "exact name or group first, then a normalised key (P5-Bi8-12 = P5Bi8-12). "
            "Pick <i>(none)</i> to leave a spectrum without values; nothing is joined "
            "on Cancel.")
        head.setWordWrap(True)
        v.addWidget(head)

        self.map = QTableWidget(len(table.rows), 3)
        self.map.setHorizontalHeaderLabels(["spectrum", "group", "CSV row"])
        self.map.verticalHeader().setVisible(False)
        self._combos: dict = {}
        for k in order:
            pos = self._pos_of[k]
            r = table.rows[k]
            self.map.setItem(pos, 0, _ro(r.display_name))
            self.map.setItem(pos, 1, _ro(r.group))
            cb = QComboBox()
            cb.addItem("(none)", None)
            for j, nm in enumerate(self._csv_names):
                cb.addItem(nm, j)
            m = self._matches[k]
            cb.setCurrentIndex(m.csv_row + 1 if m.csv_row is not None else 0)
            if m.how == "ambiguous":
                cb.setToolTip("ambiguous — several CSV rows match loosely: "
                              + ", ".join(self._csv_names[j] for j in m.candidates))
            else:
                cb.setToolTip({"exact": "exact match", "normalised":
                               "matched on a normalised key — check it",
                               "none": "no match — pick a row or leave (none)"}.get(m.how, ""))
            cb.currentIndexChanged.connect(lambda _i, kk=k: self._on_pick(kk))
            self.map.setCellWidget(pos, 2, cb)
            self._combos[k] = cb
            self._paint(k)
        self.map.resizeColumnsToContents()
        v.addWidget(self.map, 3)

        v.addWidget(QLabel("Columns to join (a detected ± partner is joined with its value):"))
        self.cols = QListWidget()
        nums = stab.numeric_headers(headers, csv_rows)
        pairs = stab.pair_error_columns(headers)
        err_set = set(pairs.values())
        candidates = [h for h in nums if h != sample_header and h not in err_set]
        oxides = [h for h in candidates if stab.canonical_oxide(h)]
        default_on = set(oxides) if oxides else set(candidates)
        for h in candidates:
            it = QListWidgetItem(h + (f"     (± {pairs[h]})" if h in pairs else ""))
            it.setData(Qt.UserRole, h)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if h in default_on else Qt.Unchecked)
            self.cols.addItem(it)
        v.addWidget(self.cols, 2)

        trow = QHBoxLayout()
        trow.addWidget(QLabel("tag as:"))
        self.tagSel = QComboBox()
        self.tagSel.addItem("analysed", "analysed")
        self.tagSel.addItem("nominal", "nominal")
        self.tagSel.addItem("none", "")
        self.tagSel.setToolTip("prefixes the column labels (analysed P2O5_mol) so a "
                               "measured and a nominal composition never look alike")
        text = " ".join([source_name] + [str(h) for h in headers]).lower()
        self.tagSel.setCurrentIndex(0 if ("analy" in text or "epma" in text) else 2)
        trow.addWidget(self.tagSel)
        trow.addStretch(1)
        v.addLayout(trow)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)

    def _on_pick(self, k: int):
        cb = self._combos[k]
        m = self._matches[k]
        m.csv_row = cb.currentData()
        m.how = "manual" if m.csv_row is not None else "none"
        self._paint(k)

    def _paint(self, k: int):
        pos = self._pos_of[k]
        it = self.map.item(pos, 0)
        if it is None:
            return
        if self._matches[k].csv_row is None:
            it.setBackground(QColor(_CHECK_CSS_COLOR))
        else:
            it.setBackground(QColor(0, 0, 0, 0))

    def matches(self) -> list:
        return list(self._matches)

    def chosen_columns(self) -> list:
        return [self.cols.item(i).data(Qt.UserRole) for i in range(self.cols.count())
                if self.cols.item(i).checkState() == Qt.Checked]

    def tag(self) -> str:
        return self.tagSel.currentData() or ""
