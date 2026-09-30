"""Fit parameter table: one NON-modal table of every line's parameters from
several fits -- the open workspaces and / or saved fit files (LARMOR
``.recipe.json``, dmfit ``.fxml`` / ``.fxmla``).

The table itself is ``larmor.fittable`` (Qt-free); this dialog only renders
it: the wide layout (one row per fit and line, ``value ± error`` with the
† ‡ § status glyphs) or the long one (one row per parameter with value,
error, status word, bounds and link), sortable by clicking a header (rows
move whole, so a fit's cells never come apart), Copy (tab-separated text a
spreadsheet pastes into cells), Export CSV…, Add fits…, Add open workspaces
and Remove selected. Rows keep their order of arrival until a header is
clicked. Opened through ``windowtray.show_tool_window``, so the workbench
stays usable while it is up.
"""
from __future__ import annotations

import re
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QDialog, QFileDialog,
    QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem,
    QVBoxLayout,
)

from larmor import fittable
from larmor.desktop import theme
from larmor.desktop.paths import is_instrument_file, suggest_save_dir

_NUM = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")
_ROLE_ENTRY = Qt.UserRole          # column 0 carries the entry index of its row
_FILTER = "Fits (*.json *.fxmla *.fxml);;LARMOR recipe (*.json);;dmfit (*.fxml *.fxmla);;All (*)"


class _Cell(QTableWidgetItem):
    """A read-only cell that sorts numerically when it starts with a number
    ('63.804 ± 0.58 †' sorts by 63.804), textually otherwise; empty cells
    sort after everything."""

    def __init__(self, text: str):
        super().__init__(text)
        self.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        m = _NUM.match(text)
        self.key = float(m.group(1)) if m else None

    def __lt__(self, other):
        a, b = self.key, getattr(other, "key", None)
        if a is not None and b is not None:
            return a < b
        if a is not None:
            return True                     # numbers before text and blanks
        if b is not None:
            return False
        ta, tb = self.text(), other.text()
        if not ta:
            return False                    # blanks last
        if not tb:
            return True
        return ta < tb


class FitTableDialog(QDialog):
    """``entries``: the ``fittable.FitEntry`` rows to start with;
    ``workspaces``: a callable returning the open workspaces' entries (the
    main window's ``workspace_fit_entries``), None when there is no window
    to ask."""

    def __init__(self, parent=None, entries=None, workspaces=None):
        super().__init__(parent)
        self._workspaces = workspaces
        self.entries: list = list(entries or [])
        self.resize(1080, 520)
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.chkLong = QCheckBox("long layout")
        self.chkLong.setToolTip(
            "one row per parameter: value, error, status word, bounds and "
            "link — the spreadsheet-friendly form. Unticked: one row per fit "
            "and line, value ± error with the glyphs († fixed, ‡ at a bound, "
            "§ linked; a held value prints without an error bar)")
        self.chkLong.toggled.connect(self._layout_changed)
        top.addWidget(self.chkLong)
        top.addStretch(1)
        self.btnAddFiles = QPushButton("Add fits…")
        self.btnAddFiles.setToolTip("add saved fits: LARMOR .recipe.json or dmfit "
                                    ".fxml / .fxmla files (several at once)")
        self.btnAddFiles.clicked.connect(self._add_files)
        self.btnAddWs = QPushButton("Add open workspaces")
        self.btnAddWs.setToolTip("add every open workspace that has a fit "
                                 "(those already listed are skipped)")
        self.btnAddWs.clicked.connect(self._add_workspaces)
        self.btnAddWs.setEnabled(callable(workspaces))
        self.btnRemove = QPushButton("Remove selected")
        self.btnRemove.setToolTip("drop the fits behind the selected rows from "
                                  "the table (nothing on disk changes)")
        self.btnRemove.clicked.connect(self._remove_selected)
        for b in (self.btnAddFiles, self.btnAddWs, self.btnRemove):
            top.addWidget(b)
        v.addLayout(top)

        self.table = QTableWidget(0, 0)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        # no sort indicator until a header is clicked: enabling sorting would
        # otherwise sort by column 0 at once and lose the order of arrival
        self.table.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setToolTip("click a header to sort by that "
                                                 "column (again to reverse)")
        v.addWidget(self.table, 1)

        bottom = QHBoxLayout()
        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color: {theme.active().text_dim};")
        bottom.addWidget(self.status, 1)
        self.btnCopy = QPushButton("Copy")
        self.btnCopy.setToolTip("copy the table as shown (order and layout) as "
                                "tab-separated text — pastes into spreadsheet cells")
        self.btnCopy.clicked.connect(self.copy_tsv)
        self.btnCsv = QPushButton("Export CSV…")
        self.btnCsv.setToolTip("write the table as shown to a UTF-8 CSV file")
        self.btnCsv.clicked.connect(self._export_csv)
        btnClose = QPushButton("Close")
        btnClose.clicked.connect(self.close)
        for b in (self.btnCopy, self.btnCsv, btnClose):
            bottom.addWidget(b)
        v.addLayout(bottom)
        self.rebuild()

    # ----------------------------------------------------------- entries
    @staticmethod
    def _key(e) -> tuple:
        return (e.name, e.source)

    def add_entries(self, entries) -> int:
        """Append entries (order of arrival kept). Returns how many."""
        entries = list(entries or [])
        self.entries.extend(entries)
        self.rebuild()
        return len(entries)

    def add_new_entries(self, entries) -> int:
        """Append only the entries not already listed (same name and source)
        -- 'Add open workspaces' pressed twice adds nothing twice."""
        have = {self._key(e) for e in self.entries}
        fresh = [e for e in (entries or []) if self._key(e) not in have]
        return self.add_entries(fresh)

    def remove_entries(self, indices) -> int:
        drop = {int(i) for i in indices if 0 <= int(i) < len(self.entries)}
        if not drop:
            return 0
        self.entries = [e for i, e in enumerate(self.entries) if i not in drop]
        self.rebuild()
        return len(drop)

    def selected_entry_indices(self) -> set:
        """The entries behind the selected rows (any cell of a row selects
        its fit), read from column 0's data so a sorted view maps back."""
        out = set()
        sel = self.table.selectionModel()
        if sel is None:
            return out
        for ix in sel.selectedRows() + sel.selectedIndexes():
            it = self.table.item(ix.row(), 0)
            if it is not None and it.data(_ROLE_ENTRY) is not None:
                out.add(int(it.data(_ROLE_ENTRY)))
        return out

    # ------------------------------------------------------------- table
    def is_long(self) -> bool:
        return self.chkLong.isChecked()

    def _layout_changed(self, *_):
        # the columns change meaning between layouts: forget the sort
        self.table.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)
        self.rebuild()

    def _built(self):
        """(headers, text rows, owner entry index per row) for the current
        layout; the per-entry row counts come from building each entry
        alone, which is how the union table lays them out too."""
        if self.is_long():
            headers, rows = fittable.build_long(self.entries)
            counts = [len(fittable.build_long([e])[1]) for e in self.entries]
        else:
            headers, rows = fittable.build_wide(self.entries)
            counts = [len(e.sites) for e in self.entries]
        owners = [i for i, n in enumerate(counts) for _ in range(n)]
        text_rows = [[fittable.cell_text(c) for c in r] for r in rows]
        return list(headers), text_rows, owners

    def rebuild(self):
        hh = self.table.horizontalHeader()
        col, order = hh.sortIndicatorSection(), hh.sortIndicatorOrder()
        headers, rows, owners = self._built()
        t = self.table
        t.setSortingEnabled(False)
        t.clear()
        t.setColumnCount(len(headers))
        t.setRowCount(len(rows))
        t.setHorizontalHeaderLabels(headers)
        for r, (row, owner) in enumerate(zip(rows, owners)):
            for c, txt in enumerate(row):
                it = _Cell(txt)
                if c == 0:
                    it.setData(_ROLE_ENTRY, owner)
                    it.setToolTip(self.entries[owner].source or "")
                t.setItem(r, c, it)
        t.resizeColumnsToContents()
        for c in range(t.columnCount()):
            t.setColumnWidth(c, min(t.columnWidth(c), 260))
        t.setSortingEnabled(True)
        if 0 <= col < t.columnCount() and rows:
            t.sortItems(col, order)
        else:
            hh.setSortIndicator(-1, Qt.AscendingOrder)
        n = len(self.entries)
        self.setWindowTitle("Fit parameter table — "
                            + (f"{n} fit{'s' if n != 1 else ''}" if n else "no fits yet"))
        if not n:
            self.status.setText("no fits yet — Add fits… (saved LARMOR recipes or "
                                "dmfit .fxml / .fxmla files) or Add open workspaces")
        elif self.is_long():
            self.status.setText(f"{len(rows)} parameter rows from {n} fit"
                                f"{'s' if n != 1 else ''} · status: fixed / linked / "
                                "at_min / at_max · click a header to sort")
        else:
            self.status.setText(f"{len(rows)} line rows from {n} fit"
                                f"{'s' if n != 1 else ''} · † fixed · ‡ at a bound · "
                                "§ linked · a held value has no error bar · click a "
                                "header to sort")
        self.btnRemove.setEnabled(n > 0)
        self.btnCopy.setEnabled(n > 0)
        self.btnCsv.setEnabled(n > 0)

    def current_table(self) -> tuple[list[str], list[list[str]]]:
        """(headers, rows) exactly as displayed: the current layout, in the
        current (possibly sorted) row order."""
        t = self.table
        headers = [t.horizontalHeaderItem(c).text() if t.horizontalHeaderItem(c)
                   else "" for c in range(t.columnCount())]
        rows = [[(t.item(r, c).text() if t.item(r, c) is not None else "")
                 for c in range(t.columnCount())] for r in range(t.rowCount())]
        return headers, rows

    # ----------------------------------------------------------- outputs
    def copy_tsv(self) -> str:
        text = fittable.to_tsv(*self.current_table())
        QApplication.clipboard().setText(text)
        self.status.setText(f"copied {self.table.rowCount()} row(s) as tab-separated "
                            "text — paste into a spreadsheet")
        return text

    def export_csv(self, path) -> bool:
        p = Path(str(path))
        if p.suffix.lower() != ".csv":
            p = p.with_suffix(".csv")
        if is_instrument_file(p):
            self.status.setText("refused: that name would replace an acquired file")
            return False
        try:
            fittable.to_csv(*self.current_table(), p)
        except OSError as exc:
            self.status.setText(f"could not write: {exc}")
            return False
        self.status.setText(f"wrote {p.name} to {p.parent}")
        return True

    def _start_dir(self) -> str:
        fallback = str(QSettings("LARMOR", "app").value("lastDir", "") or "")
        for e in reversed(self.entries):
            d = suggest_save_dir(e.source, "")
            if d:
                return d
        return fallback or str(Path.home())

    def _export_csv(self):
        if not self.entries:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export the fit parameter table",
            str(Path(self._start_dir()) / "fit_parameters.csv"), "CSV (*.csv)")
        if path:
            self.export_csv(path)

    # ------------------------------------------------------------ inputs
    def _add_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Add fits", self._start_dir(),
                                                _FILTER)
        if files:
            self.add_files(files)

    def add_files(self, paths) -> int:
        """Load fit files into the table; unreadable ones are named in the
        status line, never in a box. Returns how many were added."""
        entries, bad = [], []
        for p in paths:
            try:
                entries.append(fittable.load_fit_file(p))
            except Exception as exc:                          # noqa: BLE001
                bad.append(f"{Path(str(p)).name}: {exc}")
        n = self.add_entries(entries)
        if bad:
            self.status.setText(f"added {n} fit(s); could not read " + "; ".join(bad))
        return n

    def _add_workspaces(self):
        if not callable(self._workspaces):
            return
        entries = list(self._workspaces() or [])
        n = self.add_new_entries(entries)
        if not entries:
            self.status.setText("no open workspace has a fit")
        elif not n:
            self.status.setText("every fitted workspace is already in the table")

    def _remove_selected(self):
        picked = self.selected_entry_indices()
        if not picked:
            self.status.setText("select a row first — its fit is removed")
            return
        self.remove_entries(picked)
