"""Left explorer panel: browse folders / samples and open spectra.

Two ways in:
  * "Open sample…" scans one sample folder and lists every spectrum in it,
    each auto-identified (nucleus, 1D/2D, experiment kind);
  * the tree browses the filesystem, flagging sample folders and EXPNOs, and
    expands lazily.
Double-click (or Enter) on a spectrum opens it in the workbench.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QHBoxLayout, QLineEdit, QPushButton, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget,
)

_ROLE_PATH = Qt.UserRole
_ROLE_OPEN = Qt.UserRole + 1        # the openable data path (None for folders)

_NUC_COLOR = {
    "1H": "#4b5760", "19F": "#0e7c86", "27Al": "#1f77b4", "23Na": "#2ca02c",
    "13C": "#8c564b", "31P": "#d62728", "29Si": "#9467bd", "11B": "#e377c2",
    "17O": "#17becf", "7Li": "#bcbd22", "35Cl": "#7f7f7f",
}


class ExplorerPanel(QWidget):
    open_requested = Signal(str)        # openable data path

    def __init__(self):
        super().__init__()
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        row = QHBoxLayout()
        self.btnSample = QPushButton("Open sample…")
        self.btnSample.setToolTip("scan a sample folder and list every "
                                  "spectrum in it")
        self.btnBrowse = QPushButton("Browse…")
        self.btnBrowse.setToolTip("browse a data folder tree")
        row.addWidget(self.btnSample)
        row.addWidget(self.btnBrowse)
        v.addLayout(row)

        self.filter = QLineEdit()
        self.filter.setPlaceholderText("filter (nucleus, kind, expno)…")
        self.filter.textChanged.connect(self._apply_filter)
        v.addWidget(self.filter)

        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setIndentation(14)
        self.tree.itemActivated.connect(self._activated)
        self.tree.itemExpanded.connect(self._expanded)
        v.addWidget(self.tree, 1)

        self.btnSample.clicked.connect(self._open_sample)
        self.btnBrowse.clicked.connect(self._browse)
        self._hl: dict = {}          # color -> highlighted item

    # ------------------------------------------------------------------
    def _iter_items(self):
        stack = [self.tree.topLevelItem(i)
                 for i in range(self.tree.topLevelItemCount())]
        while stack:
            it = stack.pop()
            yield it
            for i in range(it.childCount()):
                stack.append(it.child(i))

    def highlight(self, path: str, color: str):
        """Tint the tree row whose openable data path is ``path`` (used to mark
        the spectrum picked for an HMQC projection). One row per colour."""
        old = self._hl.pop(color, None)
        if old is not None:
            old.setBackground(0, QBrush())
        target = next((it for it in self._iter_items()
                       if it.data(0, _ROLE_OPEN) == path), None)
        if target is not None:
            c = QColor(color); c.setAlpha(60)
            target.setBackground(0, QBrush(c))
            self._hl[color] = target

    # ------------------------------------------------------------------
    def _open_sample(self):
        from PySide6.QtWidgets import QFileDialog

        folder = QFileDialog.getExistingDirectory(self, "Choose a sample folder")
        if folder:
            self.load_sample(folder)

    def _browse(self):
        from PySide6.QtWidgets import QFileDialog

        folder = QFileDialog.getExistingDirectory(self, "Choose a data folder")
        if folder:
            self.load_tree(folder)

    def load_sample(self, folder: str):
        from larmor.io import scan

        self.tree.clear()
        root = QTreeWidgetItem([Path(folder).name])
        root.setData(0, _ROLE_PATH, folder)
        f = root.font(0); f.setBold(True); root.setFont(0, f)
        self.tree.addTopLevelItem(root)
        for info in scan.scan_sample(folder):
            self.tree.addTopLevelItem(self._experiment_item(info))
        root.setExpanded(True)
        self.tree.expandToDepth(0)

    def load_tree(self, folder: str):
        from larmor.io import scan

        self.tree.clear()
        top = self._folder_item(Path(folder).name, folder,
                                is_sample=scan.is_sample_folder(folder))
        self.tree.addTopLevelItem(top)
        top.setExpanded(True)          # triggers lazy population

    # ------------------------------------------------------------------
    def _experiment_item(self, info) -> QTreeWidgetItem:
        it = QTreeWidgetItem([info.label])
        it.setData(0, _ROLE_PATH, info.path)
        it.setData(0, _ROLE_OPEN, info.openable)
        it.setForeground(0, QBrush(QColor(_NUC_COLOR.get(info.nucleus, "#16202a"))))
        tip = (f"{info.nucleus} · {'2D' if info.ndim == 2 else '1D'} · "
               f"{info.kind}\npulse: {info.pulse_program}")
        if info.title:
            tip += f"\n{info.title}"
        avail = [k for k, ok in (("1r", info.has_1r), ("2rr", info.has_2rr),
                                 ("fid", info.has_fid), ("ser", info.has_ser))
                 if ok]
        tip += "\navailable: " + ", ".join(avail)
        it.setToolTip(0, tip)
        return it

    def _folder_item(self, name: str, path: str, is_sample=False,
                     is_expno=False) -> QTreeWidgetItem:
        prefix = "📁 " if not is_sample else "🧪 "
        it = QTreeWidgetItem([prefix + name])
        it.setData(0, _ROLE_PATH, path)
        if not is_expno:
            it.addChild(QTreeWidgetItem(["…"]))   # lazy placeholder
        return it

    def _expanded(self, item: QTreeWidgetItem):
        # lazy-load a folder the first time it opens
        if item.childCount() == 1 and item.child(0).text(0) == "…":
            item.takeChildren()
            self._populate(item)

    def _populate(self, item: QTreeWidgetItem):
        from larmor.io import scan

        path = item.data(0, _ROLE_PATH)
        for entry in scan.list_dir(path):
            if entry.is_expno and entry.info is not None:
                item.addChild(self._experiment_item(entry.info))
            elif entry.is_expno:
                child = QTreeWidgetItem([entry.name])
                child.setData(0, _ROLE_PATH, entry.path)
                item.addChild(child)
            else:
                item.addChild(self._folder_item(entry.name, entry.path,
                                                is_sample=entry.is_sample))

    def _activated(self, item: QTreeWidgetItem, _col: int):
        openable = item.data(0, _ROLE_OPEN)
        if openable:
            self.open_requested.emit(openable)

    def _apply_filter(self, text: str):
        text = text.strip().lower()

        def match(it: QTreeWidgetItem) -> bool:
            vis = text in it.text(0).lower() if text else True
            child_vis = False
            for i in range(it.childCount()):
                child_vis = match(it.child(i)) or child_vis
            it.setHidden(bool(text) and not (vis or child_vis))
            return vis or child_vis

        for i in range(self.tree.topLevelItemCount()):
            match(self.tree.topLevelItem(i))
