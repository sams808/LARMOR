"""Workspace manager panel: a list of open documents (1D fits, 2D maps) you can
switch between, close, or save — like TopSpin windows / ssNake workspaces —
plus the session rows (figures, batch-fit sessions) a project keeps alongside
them; those open in their own dialog on Enter / double-click (``activate``).

The panel is a thin view: the main window owns the workspace snapshots and does
the heavy lifting. Snapshots are lightweight (data arrays + recipe), and the
display widgets are shared and re-populated on switch, so many open workspaces
cost little. Row == workspace index, always.

Several rows can be selected (Ctrl / Shift, an Extended selection). A plain
click still switches the document; a click or key that EXTENDS the selection
moves Qt's current row without switching, so building a selection never
bounces the active document around. The right-click menu acts on the clicked
row (Switch to / Open, Rename…, Save) and on the whole selection (Close,
Close others, Send to Plotting studio, Fit parameter table…, Overlay on the
active spectrum); the menu's signals carry row indices and the main window
does the work.
"""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QHBoxLayout, QListWidget, QListWidgetItem, QMenu,
    QPushButton, QVBoxLayout, QWidget,
)

#: the per-row facts the main window hands rebuild() (kind, has_fit)
_ROLE_INFO = Qt.UserRole
_DOC_KINDS = ("1d", "2d")
_EXTEND = Qt.ControlModifier | Qt.ShiftModifier


class _WorkspaceList(QListWidget):
    """The rows. ``extending`` is True while a Ctrl / Shift click or key is
    being handled, so the panel can tell a selection gesture from a switch
    (Qt moves the current row in both cases, and ``currentRowChanged``
    carries no modifiers)."""

    def __init__(self):
        super().__init__()
        self.extending = False
        #: the row a currentRowChanged of the current press already switched
        #: to (None when none): itemClicked on release must not switch again
        self.switched = None

    def _handle(self, ev, call):
        prev = self.extending
        self.extending = bool(ev.modifiers() & _EXTEND)
        try:
            call(ev)
        finally:
            self.extending = prev

    def mousePressEvent(self, ev):
        self.switched = None
        self._handle(ev, super().mousePressEvent)

    def mouseReleaseEvent(self, ev):
        self._handle(ev, super().mouseReleaseEvent)

    def keyPressEvent(self, ev):
        self._handle(ev, super().keyPressEvent)


class WorkspacePanel(QWidget):
    switch = Signal(int)
    close = Signal(int)
    save = Signal(int)
    #: Enter / double-click on a row -- an explicit gesture, unlike
    #: currentRowChanged which also fires on arrow-key navigation, so a
    #: figure / batch row never pops its modal dialog while merely browsing
    activate = Signal(int)
    #: right-click ▸ Rename… on a row
    rename = Signal(int)
    #: right-click ▸ Close (the selection) / Close others: the rows to close
    close_many = Signal(list)
    #: right-click ▸ Send to Plotting studio: the selected 1D document rows
    send_to_studio = Signal(list)
    #: right-click ▸ Fit parameter table…: the selected document rows with a fit
    fit_table = Signal(list)
    #: right-click ▸ Overlay on the active spectrum: selected 1D rows, active excluded
    overlay = Signal(list)

    def __init__(self):
        super().__init__()
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4); v.setSpacing(4)
        self.list = _WorkspaceList()
        self.list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.list.setToolTip("click: switch to a document · Ctrl / Shift-click: "
                             "select several · right-click: actions on the "
                             "selection (Plotting studio, Fit parameter table, "
                             "overlay, close)")
        self.list.currentRowChanged.connect(self._row_changed)
        self.list.itemClicked.connect(self._clicked)
        self.list.itemActivated.connect(
            lambda it: self.activate.emit(self.list.row(it)))
        self.list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.list.customContextMenuRequested.connect(self._context_menu)
        v.addWidget(self.list, 1)
        row = QHBoxLayout()
        self.btnClose = QPushButton("Close")
        self.btnClose.clicked.connect(self._close)
        self.btnSave = QPushButton("Save")
        self.btnSave.setToolTip("save this workspace's fit / spectrum")
        self.btnSave.clicked.connect(self._save)
        row.addWidget(self.btnClose); row.addWidget(self.btnSave); row.addStretch(1)
        v.addLayout(row)
        self._guard = False
        self._active = -1

    def rebuild(self, items: list[tuple[str, str]], active: int, info=None):
        """items: list of (icon, title); active: current row; info: an
        optional parallel list of ``{"kind", "has_fit"}`` dicts the
        right-click menu reads (rows without one are treated as documents
        without a fit)."""
        self._guard = True
        self.list.clear()
        for k, (icon, title) in enumerate(items):
            it = QListWidgetItem(f"{icon}  {title}")
            if info is not None and k < len(info) and isinstance(info[k], dict):
                it.setData(_ROLE_INFO, dict(info[k]))
            self.list.addItem(it)
        self._active = active if 0 <= active < self.list.count() else -1
        if self._active >= 0:
            self.list.setCurrentRow(self._active)
        self._guard = False

    # ------------------------------------------------------------ queries
    def row_info(self, row: int) -> dict:
        it = self.list.item(row) if 0 <= row < self.list.count() else None
        d = it.data(_ROLE_INFO) if it is not None else None
        return dict(d) if isinstance(d, dict) else {}

    def selected_rows(self) -> list[int]:
        """The selected rows, sorted."""
        return sorted(self.list.row(it) for it in self.list.selectedItems())

    def _is_doc(self, row: int) -> bool:
        return self.row_info(row).get("kind", "1d") in _DOC_KINDS

    # ------------------------------------------------------------ switching
    def _row_changed(self, row: int):
        if not self._guard and row >= 0 and not self.list.extending:
            self.list.switched = row
            self.switch.emit(row)

    def _clicked(self, item):
        """A plain click on a row that is already Qt's current row (it was
        Ctrl-clicked into the selection a moment ago) must still switch:
        currentRowChanged stays silent then. A row the press already
        switched to is not switched to twice."""
        if self.list.extending or self._guard:
            return
        row = self.list.row(item)
        if row >= 0 and row != self._active and self.list.switched != row:
            self.switch.emit(row)

    def _close(self):
        r = self.list.currentRow()
        if r >= 0:
            self.close.emit(r)

    def _save(self):
        r = self.list.currentRow()
        if r >= 0:
            self.save.emit(r)

    # --------------------------------------------------------- right-click
    def _context_menu(self, pos):
        menu = self.menu_at(pos)
        if menu is not None:
            menu.exec(self.list.viewport().mapToGlobal(pos))

    def menu_at(self, pos) -> QMenu | None:
        """The menu a right-click at viewport ``pos`` gets (built, not
        shown): a click outside the selection first makes the clicked row
        the whole selection; None off any row."""
        item = self.list.itemAt(pos)
        if item is None:
            return None
        if not item.isSelected():            # a right-click outside the selection
            self.list.clearSelection()       # acts on the clicked row alone
            item.setSelected(True)
        return self.build_menu(self.list.row(item), self.selected_rows())

    def build_menu(self, row: int, rows: list[int]) -> QMenu:
        """The right-click menu for a click on ``row`` with ``rows``
        selected (built, not shown -- tests trigger its actions)."""
        rows = sorted({r for r in rows if 0 <= r < self.list.count()} | {row})
        m = QMenu(self)
        m.setToolTipsVisible(True)
        is_doc = self._is_doc(row)
        a = m.addAction("Switch to" if is_doc else "Open")
        a.setToolTip("show this document on the workbench" if is_doc
                     else "reopen this figure / batch session in its dialog")
        a.setEnabled(not (is_doc and row == self._active))
        a.triggered.connect(lambda *_: self.activate.emit(row))
        a = m.addAction("Rename…")
        a.setToolTip("a name of your own for this row (kept when you switch "
                     "documents; empty restores the document's own title)")
        a.triggered.connect(lambda *_: self.rename.emit(row))
        a = m.addAction("Save")
        a.setToolTip("save this workspace's fit / spectrum")
        a.triggered.connect(lambda *_: self.save.emit(row))
        n = len(rows)
        a = m.addAction("Close" if n <= 1 else f"Close {n} selected")
        a.triggered.connect(lambda *_, rs=list(rows): self.close_many.emit(rs))
        others = [r for r in range(self.list.count()) if r not in rows]
        a = m.addAction("Close others")
        a.setToolTip("close every workspace that is not selected")
        a.setEnabled(bool(others))
        a.triggered.connect(lambda *_, rs=others: self.close_many.emit(rs))
        m.addSeparator()

        docs_1d = [r for r in rows if self.row_info(r).get("kind", "1d") == "1d"]
        fitted = [r for r in rows if self._is_doc(r) and self.row_info(r).get("has_fit")]
        k = len(docs_1d)
        a = m.addAction("Send to Plotting studio" if k <= 1 else
                        f"Send to Plotting studio  ({k} spectra)")
        a.setToolTip("open the Plotting studio with every selected 1D "
                     "workspace: its experiment, the fitted total (dashed) "
                     "and each component, in one colour per spectrum")
        a.setEnabled(k >= 1)
        a.triggered.connect(lambda *_, rs=docs_1d: self.send_to_studio.emit(rs))
        f = len(fitted)
        a = m.addAction("Fit parameter table…" if f <= 1 else
                        f"Fit parameter table…  ({f} fits)")
        a.setToolTip("one table of every line's parameters from the selected "
                     "fits (value ± error, † fixed, ‡ at a bound, § linked); "
                     "add saved fit files there too, copy or export it")
        a.setEnabled(f >= 1)
        a.triggered.connect(lambda *_, rs=fitted: self.fit_table.emit(rs))
        ov = [r for r in docs_1d if r != self._active]
        a = m.addAction("Overlay on the active spectrum" if len(ov) <= 1 else
                        f"Overlay on the active spectrum  ({len(ov)})")
        a.setToolTip("draw the selected spectra behind the active one, as "
                     "compared spectra in the Datasets dock (the active fit "
                     "is untouched)")
        a.setEnabled(bool(ov) and self.row_info(self._active).get("kind", "1d") == "1d"
                     and self._active >= 0)
        a.triggered.connect(lambda *_, rs=ov: self.overlay.emit(rs))
        return m
