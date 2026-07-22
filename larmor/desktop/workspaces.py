"""Workspace manager panel: a list of open documents (1D fits, 2D maps) you can
switch between, close, or save — like TopSpin windows / ssNake workspaces.

The panel is a thin view: the main window owns the workspace snapshots and does
the heavy lifting. Snapshots are lightweight (data arrays + recipe), and the
display widgets are shared and re-populated on switch, so many open workspaces
cost little.
"""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QHBoxLayout, QListWidget, QListWidgetItem, QPushButton, QVBoxLayout, QWidget,
)


class WorkspacePanel(QWidget):
    switch = Signal(int)
    close = Signal(int)
    save = Signal(int)

    def __init__(self):
        super().__init__()
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4); v.setSpacing(4)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._row_changed)
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

    def rebuild(self, items: list[tuple[str, str]], active: int):
        """items: list of (icon, title); active: current row."""
        self._guard = True
        self.list.clear()
        for icon, title in items:
            it = QListWidgetItem(f"{icon}  {title}")
            self.list.addItem(it)
        if 0 <= active < self.list.count():
            self.list.setCurrentRow(active)
        self._guard = False

    def _row_changed(self, row: int):
        if not self._guard and row >= 0:
            self.switch.emit(row)

    def _close(self):
        r = self.list.currentRow()
        if r >= 0:
            self.close.emit(r)

    def _save(self):
        r = self.list.currentRow()
        if r >= 0:
            self.save.emit(r)
