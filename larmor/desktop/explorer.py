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
    QCheckBox, QHBoxLayout, QLineEdit, QPushButton, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget,
)

_ROLE_PATH = Qt.UserRole
_ROLE_OPEN = Qt.UserRole + 1        # the openable data path (None for folders)
_ROLE_KIND = Qt.UserRole + 2        # "exp" | "ph_proc" | "ph_fit" (placeholders)

#: fit files browsable under a proc, with a human origin
_FIT_EXT = {".json": "LARMOR recipe", ".fxml": "dmfit fit", ".fxmla": "dmfit fit"}


def _list_fits(folder: Path) -> list[Path]:
    """LARMOR .recipe.json and dmfit .fxml/.fxmla files saved in a folder."""
    out: list[Path] = []
    try:
        for f in sorted(folder.iterdir()):
            n = f.name.lower()
            if f.is_file() and (n.endswith(".recipe.json") or n.endswith(".fxml")
                                or n.endswith(".fxmla")):
                out.append(f)
    except OSError:
        pass
    return out

_NUC_COLOR = {
    "1H": "#4b5760", "19F": "#0e7c86", "27Al": "#1f77b4", "23Na": "#2ca02c",
    "13C": "#8c564b", "31P": "#d62728", "29Si": "#9467bd", "11B": "#e377c2",
    "17O": "#17becf", "7Li": "#bcbd22", "35Cl": "#7f7f7f",
}


class ExplorerPanel(QWidget):
    open_requested = Signal(str)        # openable data path
    batch_requested = Signal(list)      # openable paths for a batch fit

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

        self._show_fits = False      # show the pdata proc + fits layers
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setIndentation(14)
        self.tree.setSelectionMode(QTreeWidget.ExtendedSelection)   # Ctrl/Shift
        self.tree.itemActivated.connect(self._activated)
        self.tree.itemExpanded.connect(self._expanded)
        v.addWidget(self.tree, 1)

        self.btnBatch = QPushButton("Batch fit selected…")
        self.btnBatch.setToolTip("Ctrl/Shift-click several spectra above, then "
                                 "batch-fit them with one shared model "
                                 "(amplitudes free per spectrum)")
        self.btnBatch.clicked.connect(self._batch_clicked)
        v.addWidget(self.btnBatch)

        # a toggle at the bottom: reveal each experiment's pdata proc folders and,
        # under a proc, the fits saved in it (pick the proc you fit on)
        self.chkFits = QCheckBox("Show pdata procs && their fits")
        self.chkFits.setToolTip(
            "under each experiment, list its pdata proc folders; expand a proc "
            "to see the fits in it (LARMOR .recipe.json / dmfit .fxml) — and "
            "double-click a proc to open that processing for fitting")
        self.chkFits.toggled.connect(self._toggle_fits)
        v.addWidget(self.chkFits)

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
        it.setData(0, _ROLE_KIND, "exp")
        if self._show_fits:
            self._add_proc_placeholder(it)
        return it

    def _add_proc_placeholder(self, exp_item: QTreeWidgetItem):
        """Give an experiment a lazy child so it can be expanded to its procs."""
        ph = QTreeWidgetItem(["…"])
        ph.setData(0, _ROLE_KIND, "ph_proc")
        exp_item.addChild(ph)

    def _folder_item(self, name: str, path: str, is_sample=False,
                     is_expno=False) -> QTreeWidgetItem:
        prefix = "📁 " if not is_sample else "🧪 "
        it = QTreeWidgetItem([prefix + name])
        it.setData(0, _ROLE_PATH, path)
        if not is_expno:
            it.addChild(QTreeWidgetItem(["…"]))   # lazy placeholder
        return it

    def _expanded(self, item: QTreeWidgetItem):
        # lazy-load the placeholder child the first time a node opens
        if item.childCount() != 1:
            return
        ch = item.child(0)
        kind = ch.data(0, _ROLE_KIND)
        if kind == "ph_proc":
            item.takeChildren()
            self._populate_procs(item)
        elif kind == "ph_fit":
            item.takeChildren()
            self._populate_fits(item)
        elif ch.text(0) == "…":
            item.takeChildren()
            self._populate(item)

    def _populate_procs(self, exp_item: QTreeWidgetItem):
        """List an experiment's pdata proc folders (each one openable)."""
        expno = Path(exp_item.data(0, _ROLE_PATH) or "")
        pdata = expno / "pdata"
        procs = sorted((d for d in pdata.iterdir() if d.is_dir()),
                       key=lambda d: (not d.name.isdigit(), int(d.name)
                                      if d.name.isdigit() else d.name)) \
            if pdata.is_dir() else []
        if not procs:
            none = QTreeWidgetItem(["(no pdata procs)"])
            none.setDisabled(True)
            exp_item.addChild(none)
            return
        for d in procs:
            openable = None
            for f in ("2rr", "1r"):
                if (d / f).exists():
                    openable = str(d / f)
                    break
            n_fits = len(_list_fits(d))
            label = f"proc {d.name}" + (f"   ({n_fits} fit{'s' * (n_fits != 1)})"
                                        if n_fits else "")
            it = QTreeWidgetItem(["⚙ " + label])
            it.setData(0, _ROLE_PATH, str(d))
            it.setData(0, _ROLE_OPEN, openable)      # double-click opens this proc
            it.setToolTip(0, f"{d}\ndouble-click to fit on this processing")
            ph = QTreeWidgetItem(["…"])
            ph.setData(0, _ROLE_KIND, "ph_fit")
            it.addChild(ph)
            exp_item.addChild(it)

    def _populate_fits(self, proc_item: QTreeWidgetItem):
        """List the fit files saved in a proc folder (name + origin)."""
        d = Path(proc_item.data(0, _ROLE_PATH) or "")
        fits = _list_fits(d)
        if not fits:
            none = QTreeWidgetItem(["(no fits here)"])
            none.setDisabled(True)
            proc_item.addChild(none)
            return
        for f in fits:
            origin = _FIT_EXT.get("".join(f.suffixes[-1:]).lower(), "fit")
            it = QTreeWidgetItem([f"📄 {f.name}"])
            it.setData(0, _ROLE_PATH, str(f))
            it.setData(0, _ROLE_OPEN, str(f))        # double-click opens the fit
            it.setToolTip(0, f"{origin}\n{f}")
            proc_item.addChild(it)

    def _toggle_fits(self, on: bool):
        self._show_fits = on
        # add/remove the proc layer on the experiments already in the tree
        for it in list(self._iter_items()):
            if it.data(0, _ROLE_KIND) != "exp":
                continue
            it.takeChildren()
            it.setExpanded(False)
            if on:
                self._add_proc_placeholder(it)

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

    def _batch_clicked(self):
        # gather the openable spectra from the selected rows (skip fit files and
        # plain folders); de-duplicate, keep tree order
        paths, seen = [], set()
        for it in self.tree.selectedItems():
            op = it.data(0, _ROLE_OPEN)
            if not op or op in seen:
                continue
            low = op.lower()
            if low.endswith((".recipe.json", ".fxml", ".fxmla")):
                continue                       # a fit, not a spectrum
            paths.append(op)
            seen.add(op)
        self.batch_requested.emit(paths)

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
