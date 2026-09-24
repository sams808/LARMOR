"""Left explorer panel: browse folders / samples, pick a proc, open a spectrum.

  * "Open sample…" scans one sample folder and lists every spectrum in it;
  * "Browse…" walks the filesystem, flagging sample folders and EXPNOs;
  * two toggles reveal the pdata **proc** layer and the **fits** saved in each
    proc (both on by default);
  * folders can be **pinned** so they come back next session;
  * double-clicking a spectrum opens it — an experiment with several procs asks
    which proc to fit on;
  * right-click ▸ **Rename…** on a sample folder or an EXPNO gives it a display
    name kept by LARMOR (``larmor.aliases``) or renames the folder on disk.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QHBoxLayout, QInputDialog, QLabel,
    QLineEdit, QMenu, QPushButton, QRadioButton, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget,
)

from larmor import aliases

_ROLE_PATH = Qt.UserRole
_ROLE_OPEN = Qt.UserRole + 1        # the openable data path (None for folders)
_ROLE_KIND = Qt.UserRole + 2        # exp|proc|fit|folder + ph_proc|ph_fit placeholders
_ROLE_BASE = Qt.UserRole + 3        # the row's text without any alias (relabelling)

_FIT_EXT = {".json": "LARMOR recipe", ".fxml": "dmfit fit", ".fxmla": "dmfit fit"}

_NUC_COLOR = {
    "1H": "#4b5760", "19F": "#0e7c86", "27Al": "#1f77b4", "23Na": "#2ca02c",
    "13C": "#8c564b", "31P": "#d62728", "29Si": "#9467bd", "11B": "#e377c2",
    "17O": "#17becf", "7Li": "#bcbd22", "35Cl": "#7f7f7f",
}


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


def _procs_of(expno) -> list[Path]:
    """pdata proc folders of an experiment, sorted (numeric first)."""
    pdata = Path(expno) / "pdata"
    if not pdata.is_dir():
        return []
    try:
        procs = [d for d in pdata.iterdir() if d.is_dir()]
    except OSError:
        return []
    return sorted(procs, key=lambda d: (not d.name.isdigit(),
                                        int(d.name) if d.name.isdigit() else d.name))


def _proc_openable(procdir: Path):
    for f in ("2rr", "1r"):
        if (procdir / f).exists():
            return str(procdir / f)
    return None


class ExplorerPanel(QWidget):
    open_requested = Signal(str)        # openable data path
    batch_requested = Signal(list)      # openable paths for a batch fit
    inventory_requested = Signal(str)   # a month or sample folder -> Session inventory
    overlay_requested = Signal(str)     # an EXPNO -> overlay on the active spectrum
    renamed = Signal(str, str)          # (old path, new path) -- equal for an alias

    def __init__(self):
        super().__init__()
        from PySide6.QtCore import QSettings
        #: the data paths open in LARMOR (the window sets this) -- a folder
        #: holding one of them is never renamed on disk
        self.open_paths = lambda: []
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        row = QHBoxLayout()
        self.btnSample = QPushButton("Open sample…")
        self.btnSample.setToolTip("scan a sample folder and list every spectrum")
        self.btnBrowse = QPushButton("Browse…")
        self.btnBrowse.setToolTip("browse a data folder tree")
        row.addWidget(self.btnSample)
        row.addWidget(self.btnBrowse)
        v.addLayout(row)

        self.filter = QLineEdit()
        self.filter.setPlaceholderText("filter (nucleus, kind, expno)…")
        self.filter.textChanged.connect(self._apply_filter)
        v.addWidget(self.filter)

        self._show_procs = True
        self._show_fits = True
        self._pinned = list(QSettings("LARMOR", "app").value("pinnedFolders", []) or [])
        if isinstance(self._pinned, str):
            self._pinned = [self._pinned]
        # optional display names for pins (path -> name), JSON in QSettings so
        # the round-trip is platform-independent
        try:
            import json
            raw = QSettings("LARMOR", "app").value("pinnedNames", "") or ""
            self._pin_names = {k: v for k, v in json.loads(raw).items()
                               if isinstance(v, str)} if raw else {}
        except Exception:
            self._pin_names = {}

        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setIndentation(14)
        self.tree.setSelectionMode(QTreeWidget.ExtendedSelection)   # Ctrl/Shift
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._context_menu)
        self.tree.itemActivated.connect(self._activated)
        self.tree.itemExpanded.connect(self._expanded)
        v.addWidget(self.tree, 1)

        self.btnBatch = QPushButton("Batch fit selected…")
        self.btnBatch.setToolTip("Ctrl/Shift-click several spectra above, then "
                                 "batch-fit them with one shared model")
        self.btnBatch.clicked.connect(self._batch_clicked)
        v.addWidget(self.btnBatch)

        togs = QHBoxLayout()
        self.chkProcs = QCheckBox("procs")
        self.chkProcs.setChecked(True)
        self.chkProcs.setToolTip("show each experiment's pdata proc folders "
                                 "(when it has more than one)")
        self.chkProcs.toggled.connect(self._toggle_procs)
        self.chkFits = QCheckBox("fits")
        self.chkFits.setChecked(True)
        self.chkFits.setToolTip("show the fits saved in a proc "
                                "(LARMOR .recipe.json / dmfit .fxml)")
        self.chkFits.toggled.connect(self._toggle_fits)
        togs.addWidget(self.chkProcs)
        togs.addWidget(self.chkFits)
        togs.addStretch(1)
        v.addLayout(togs)

        self.btnSample.clicked.connect(self._open_sample)
        self.btnBrowse.clicked.connect(self._browse)
        self._hl: dict = {}
        self._readd_pins()

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
        old = self._hl.pop(color, None)
        if old is not None:
            old.setBackground(0, QBrush())
        target = next((it for it in self._iter_items()
                       if it.data(0, _ROLE_OPEN) == path), None)
        if target is not None:
            c = QColor(color); c.setAlpha(60)
            target.setBackground(0, QBrush(c))
            self._hl[color] = target

    # ---------------- pins ----------------
    def _pin_label(self, path: str) -> str:
        """The pin's own name, else the folder's LARMOR display name (an
        alias set with Rename…), else the folder name."""
        return self._pin_names.get(path) or aliases.display_name(path)

    def _save_pin_names(self):
        import json

        from PySide6.QtCore import QSettings
        QSettings("LARMOR", "app").setValue("pinnedNames",
                                            json.dumps(self._pin_names))

    def _readd_pins(self):
        for p in self._pinned:
            if Path(p).exists():
                it = self._folder_item(self._pin_label(p), p,
                                       is_sample=self._is_sample(p), pinned=True)
                self.tree.addTopLevelItem(it)

    @staticmethod
    def _is_sample(folder) -> bool:
        try:
            from larmor.io import scan
            return scan.is_sample_folder(folder)
        except Exception:
            return False

    def _pin(self, path: str):
        from PySide6.QtCore import QSettings
        if path not in self._pinned:
            self._pinned.insert(0, path)
            QSettings("LARMOR", "app").setValue("pinnedFolders", self._pinned)
            it = self._folder_item(self._pin_label(path), path,
                                   is_sample=self._is_sample(path), pinned=True)
            self.tree.insertTopLevelItem(0, it)

    def _unpin(self, path: str):
        from PySide6.QtCore import QSettings
        self._pinned = [p for p in self._pinned if p != path]
        QSettings("LARMOR", "app").setValue("pinnedFolders", self._pinned)
        if self._pin_names.pop(path, None) is not None:
            self._save_pin_names()
        for i in range(self.tree.topLevelItemCount()):
            it = self.tree.topLevelItem(i)
            if it.data(0, _ROLE_PATH) == path and it.data(0, _ROLE_KIND) == "pin":
                self.tree.takeTopLevelItem(i)
                break

    def set_pin_name(self, path: str, name: str):
        """Give a pinned folder a display name (empty = back to the folder's
        own name). Persisted, so it survives restarts."""
        name = (name or "").strip()
        if name and name != Path(path).name:
            self._pin_names[path] = name
        else:
            self._pin_names.pop(path, None)
        self._save_pin_names()
        for i in range(self.tree.topLevelItemCount()):
            it = self.tree.topLevelItem(i)
            if it.data(0, _ROLE_PATH) == path and it.data(0, _ROLE_KIND) == "pin":
                it.setText(0, "📌 " + self._pin_label(path))

    def _rename_pin(self, path: str):
        name, ok = QInputDialog.getText(
            self, "Rename pinned folder",
            f"Display name for\n{path}\n(leave empty to use the folder name):",
            text=self._pin_label(path))
        if ok:
            self.set_pin_name(path, name)

    def _context_menu(self, pos):
        item = self.tree.itemAt(pos)
        if item is None:
            return
        path = item.data(0, _ROLE_PATH)
        if not path or not Path(path).is_dir():
            return
        m = QMenu(self)
        if (Path(path) / "acqus").exists():          # an EXPNO: it has a story
            m.addAction("Dataset info…", lambda: self._dataset_info(path))
            m.addAction("Overlay on the current spectrum  (compare, keeps the fit)",
                        lambda: self.overlay_requested.emit(path))
            m.addSeparator()
        else:                                        # a month or sample folder
            m.addAction("Session inventory…  (production picks for every sample)",
                        lambda: self.inventory_requested.emit(path))
            m.addSeparator()
        if aliases.is_expno(path) or self._is_sample(path):
            m.addAction("Rename…", lambda: self.rename_item(path))
            m.addSeparator()
        if path in self._pinned:
            m.addAction("Rename pin…", lambda: self._rename_pin(path))
            m.addAction("Unpin folder", lambda: self._unpin(path))
        else:
            m.addAction("📌 Pin folder", lambda: self._pin(path))
        m.exec(self.tree.viewport().mapToGlobal(pos))

    @staticmethod
    def dataset_info_text(path: str) -> str:
        """The dataset's own description — the full TITLE text the way
        TopSpin's title page shows it, plus the key acquisition facts."""
        from larmor.io import bruker
        parts = []
        try:
            exp = bruker.read_expno(path, verify=False)
            parts.append(exp.summary)
        except Exception as exc:                              # noqa: BLE001
            parts.append(f"(acquisition summary unavailable: {exc})")
        # the summary shows only the title's first line; give the whole text
        for procdir in sorted((Path(path) / "pdata").glob("*")):
            tf = procdir / "title"
            if tf.exists():
                try:
                    full = tf.read_text(errors="replace").strip()
                except OSError:
                    continue
                if full and "\n" in full:
                    parts.append(f"\ntitle ({procdir.name}):\n{full}")
                break
        return "\n".join(parts)

    def _dataset_info(self, path: str):
        from PySide6.QtWidgets import (QDialog, QDialogButtonBox,
                                       QPlainTextEdit, QVBoxLayout)
        dlg = QDialog(self)
        dlg.setWindowTitle(f"{Path(path).name} — dataset info")
        dlg.resize(560, 380)
        v = QVBoxLayout(dlg)
        txt = QPlainTextEdit(self.dataset_info_text(path))
        txt.setReadOnly(True)
        f = txt.font(); f.setFamily("Consolas"); txt.setFont(f)
        v.addWidget(txt)
        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(dlg.reject)
        bb.accepted.connect(dlg.accept)
        bb.button(QDialogButtonBox.Close).clicked.connect(dlg.close)
        v.addWidget(bb)
        from larmor.desktop.windowtray import show_tool_window
        show_tool_window(dlg)                    # read-only: non-modal

    # ---------------- loading ----------------
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
        self._readd_pins()
        root = QTreeWidgetItem([aliases.display_name(folder)])
        root.setData(0, _ROLE_PATH, folder)
        root.setData(0, _ROLE_KIND, "folder")
        root.setData(0, _ROLE_BASE, "")
        if aliases.alias_for(folder):
            root.setToolTip(0, folder)
        f = root.font(0); f.setBold(True); root.setFont(0, f)
        self.tree.addTopLevelItem(root)
        for info in scan.scan_sample(folder):
            self.tree.addTopLevelItem(self._experiment_item(info))
        root.setExpanded(True)

    def load_tree(self, folder: str):
        from larmor.io import scan
        self.tree.clear()
        self._readd_pins()
        top = self._folder_item(Path(folder).name, folder,
                                is_sample=scan.is_sample_folder(folder))
        self.tree.addTopLevelItem(top)
        top.setExpanded(True)

    # ---------------- items ----------------
    def _experiment_item(self, info) -> QTreeWidgetItem:
        it = QTreeWidgetItem([self._exp_text(info.label, info.path)])
        it.setData(0, _ROLE_PATH, info.path)
        it.setData(0, _ROLE_OPEN, info.openable)
        it.setData(0, _ROLE_KIND, "exp")
        it.setData(0, _ROLE_BASE, info.label)
        it.setForeground(0, QBrush(QColor(_NUC_COLOR.get(info.nucleus, "#16202a"))))
        tip = (f"{info.nucleus} · {'2D' if info.ndim == 2 else '1D'} · {info.kind}\n"
               f"pulse: {info.pulse_program}")
        if info.title:
            tip += f"\n{info.title}"
        it.setToolTip(0, tip)
        self._reset_exp_children(it)
        return it

    def _reset_exp_children(self, it: QTreeWidgetItem):
        """Give an experiment the right expandable layer: procs when it has more
        than one, else — for a single proc that holds fits — the fits directly
        (no redundant proc level; double-click still opens that proc)."""
        it.takeChildren(); it.setExpanded(False)
        procs = _procs_of(it.data(0, _ROLE_PATH) or "")
        if self._show_procs and len(procs) > 1:
            self._add_ph(it, "ph_proc")
        elif self._show_fits and len(procs) == 1 and _list_fits(procs[0]):
            self._add_ph(it, "ph_expfit")          # single proc → its fits here

    def _add_ph(self, parent: QTreeWidgetItem, kind: str):
        ph = QTreeWidgetItem(["…"])
        ph.setData(0, _ROLE_KIND, kind)
        parent.addChild(ph)

    def _folder_item(self, name: str, path: str, is_sample=False,
                     is_expno=False, pinned=False) -> QTreeWidgetItem:
        prefix = "📌 " if pinned else ("🧪 " if is_sample else "📁 ")
        if not pinned:
            name = aliases.display_name(path)     # the LARMOR display name
        it = QTreeWidgetItem([prefix + name])
        it.setData(0, _ROLE_PATH, path)
        it.setData(0, _ROLE_KIND, "pin" if pinned else "folder")
        it.setData(0, _ROLE_BASE, prefix)
        if pinned or name != Path(path).name:   # an aliased row says where it points
            it.setToolTip(0, path)
        if not is_expno:
            it.addChild(QTreeWidgetItem(["…"]))    # lazy placeholder
        return it

    def _expanded(self, item: QTreeWidgetItem):
        if item.childCount() != 1:
            return
        ch = item.child(0)
        kind = ch.data(0, _ROLE_KIND)
        if kind == "ph_proc":
            item.takeChildren(); self._populate_procs(item)
        elif kind == "ph_fit":
            item.takeChildren(); self._populate_fits(item)
        elif kind == "ph_expfit":                  # single-proc experiment's fits
            item.takeChildren()
            procs = _procs_of(item.data(0, _ROLE_PATH) or "")
            if procs:
                self._add_fit_items(item, procs[0])
        elif ch.text(0) == "…":
            item.takeChildren(); self._populate(item)

    def _populate_procs(self, exp_item: QTreeWidgetItem):
        for d in _procs_of(exp_item.data(0, _ROLE_PATH) or ""):
            n_fits = len(_list_fits(d))
            label = f"proc {d.name}" + (f"   ({n_fits} fit{'s' * (n_fits != 1)})"
                                        if n_fits else "")
            it = QTreeWidgetItem(["⚙ " + label])
            it.setData(0, _ROLE_PATH, str(d))
            it.setData(0, _ROLE_OPEN, _proc_openable(d))
            it.setData(0, _ROLE_KIND, "proc")
            it.setToolTip(0, f"{d}\ndouble-click to fit on this processing")
            if self._show_fits and n_fits:         # only expandable if it has fits
                self._add_ph(it, "ph_fit")
            exp_item.addChild(it)

    def _populate_fits(self, proc_item: QTreeWidgetItem):
        self._add_fit_items(proc_item, Path(proc_item.data(0, _ROLE_PATH) or ""))

    def _add_fit_items(self, parent: QTreeWidgetItem, folder):
        for f in _list_fits(Path(folder)):
            origin = _FIT_EXT.get("".join(f.suffixes[-1:]).lower(), "fit")
            it = QTreeWidgetItem([f"📄 {f.name}"])
            it.setData(0, _ROLE_PATH, str(f))
            it.setData(0, _ROLE_OPEN, str(f))
            it.setData(0, _ROLE_KIND, "fit")
            it.setToolTip(0, f"{origin}\n{f}")
            parent.addChild(it)

    def _toggle_procs(self, on: bool):
        self._show_procs = on
        for it in list(self._iter_items()):
            if it.data(0, _ROLE_KIND) == "exp":
                self._reset_exp_children(it)

    def _toggle_fits(self, on: bool):
        self._show_fits = on
        for it in list(self._iter_items()):
            kind = it.data(0, _ROLE_KIND)
            if kind == "exp":
                # single-proc experiments carry their fits directly — refresh them;
                # multi-proc ones keep their proc layer (handled below)
                if len(_procs_of(it.data(0, _ROLE_PATH) or "")) <= 1:
                    self._reset_exp_children(it)
            elif kind == "proc":
                it.takeChildren(); it.setExpanded(False)
                if on and _list_fits(Path(it.data(0, _ROLE_PATH) or "")):
                    self._add_ph(it, "ph_fit")

    def _populate(self, item: QTreeWidgetItem):
        from larmor.io import scan
        for entry in scan.list_dir(item.data(0, _ROLE_PATH)):
            if entry.is_expno and entry.info is not None:
                item.addChild(self._experiment_item(entry.info))
            elif entry.is_expno:
                child = QTreeWidgetItem([aliases.display_name(entry.path)])
                child.setData(0, _ROLE_PATH, entry.path)
                child.setData(0, _ROLE_BASE, entry.name)
                item.addChild(child)
            else:
                item.addChild(self._folder_item(entry.name, entry.path,
                                                is_sample=entry.is_sample))

    # ---------------- open ----------------
    def _activated(self, item: QTreeWidgetItem, _col: int):
        if item.data(0, _ROLE_KIND) == "exp":
            procs = _procs_of(item.data(0, _ROLE_PATH) or "")
            if len(procs) > 1:                     # ask which proc to fit on
                d = self._ask_proc(procs)
                if d is not None:
                    op = _proc_openable(d)
                    if op:
                        self.open_requested.emit(op)
                return
            if len(procs) == 1:
                op = _proc_openable(procs[0]) or item.data(0, _ROLE_OPEN)
                if op:
                    self.open_requested.emit(op)
                return
        openable = item.data(0, _ROLE_OPEN)
        if openable:
            self.open_requested.emit(openable)

    def _ask_proc(self, procs: list[Path]):
        labels = []
        for d in procs:
            n = len(_list_fits(d))
            labels.append(f"proc {d.name}" + (f"  ({n} fit{'s' * (n != 1)})"
                                              if n else ""))
        choice, ok = QInputDialog.getItem(
            self, "Choose processing",
            "This experiment has several pdata procs — open which one?",
            labels, 0, False)
        if not ok:
            return None
        return procs[labels.index(choice)]

    def selected_spectra(self) -> list[str]:
        """The openable spectra in the selected rows (fits and folders skipped),
        de-duplicated, in tree order."""
        paths, seen = [], set()
        for it in self.tree.selectedItems():
            op = it.data(0, _ROLE_OPEN)
            if not op or op in seen:
                continue
            if op.lower().endswith((".recipe.json", ".fxml", ".fxmla")):
                continue                           # a fit, not a spectrum
            paths.append(op); seen.add(op)
        return paths

    def _batch_clicked(self):
        self.batch_requested.emit(self.selected_spectra())

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

    # ---------------- rename (alias or on disk) ----------------
    @staticmethod
    def _exp_text(base: str, path: str) -> str:
        """An experiment row: ``<expno> · <nucleus> …`` or, aliased,
        ``<alias> (<expno>) · <nucleus> …`` -- the number stays visible."""
        alias = aliases.alias_for(path)
        if not alias:
            return base
        expno = Path(path).name
        if base.startswith(expno):
            return f"{alias} ({expno})" + base[len(expno):]
        return f"{alias} ({base})"

    def _relabel(self, path: str):
        """Rewrite the text of every row that shows ``path`` (an alias was
        set or cleared)."""
        for it in self._iter_items():
            if it.data(0, _ROLE_PATH) != path:
                continue
            kind = it.data(0, _ROLE_KIND)
            base = it.data(0, _ROLE_BASE) or ""
            if kind == "exp":
                it.setText(0, self._exp_text(base, path))
            elif kind in ("folder", "pin"):
                name = (self._pin_label(path) if kind == "pin"
                        else aliases.display_name(path))
                it.setText(0, base + name)
                it.setToolTip(0, path if name != Path(path).name else "")
            elif kind is None and base:                  # an EXPNO without info
                it.setText(0, aliases.display_name(path))

    def _retarget(self, old: str, new: str):
        """After a rename on disk: every stored path under ``old`` now lives
        under ``new`` (rows, openables, pins), and the renamed row is
        relabelled."""
        def moved(p):
            if not p:
                return p
            if p == old:
                return new
            for sep in ("\\", "/"):
                if p.startswith(old + sep):
                    return new + p[len(old):]
            return p

        for it in self._iter_items():
            for role in (_ROLE_PATH, _ROLE_OPEN):
                v = it.data(0, role)
                if isinstance(v, str) and moved(v) != v:
                    it.setData(0, role, moved(v))
            if it.data(0, _ROLE_PATH) == new:
                kind = it.data(0, _ROLE_KIND)
                base = it.data(0, _ROLE_BASE) or ""
                if kind == "exp":
                    it.setData(0, _ROLE_BASE, base.replace(Path(old).name, Path(new).name, 1))
                elif kind is None and base:
                    it.setData(0, _ROLE_BASE, Path(new).name)
        if any(moved(p) != p for p in self._pinned):
            from PySide6.QtCore import QSettings
            self._pinned = [moved(p) for p in self._pinned]
            QSettings("LARMOR", "app").setValue("pinnedFolders", self._pinned)
            self._pin_names = {moved(k): v for k, v in self._pin_names.items()}
            self._save_pin_names()
        self._relabel(new)

    def rename_item(self, path: str):
        """Right-click ▸ Rename… on a sample folder or an EXPNO: a display
        name kept by LARMOR (an alias), or the folder renamed on disk after
        an explicit confirmation that states both paths."""
        from PySide6.QtWidgets import QMessageBox

        dlg = RenameDialog(self, path)
        if dlg.exec() != RenameDialog.Accepted:
            return
        name = dlg.name()
        if not dlg.on_disk():
            aliases.set_alias(path, name)
            self._relabel(path)
            self.renamed.emit(path, path)
            return
        try:
            target = aliases.check_rename(path, name, open_paths=self.open_paths())
        except aliases.RenameError as exc:
            QMessageBox.warning(self, "Rename on disk", str(exc))
            return
        kind = "EXPNO" if aliases.is_expno(path) else "sample folder"
        ans = QMessageBox.question(
            self, "Rename on disk",
            f"Rename this {kind} on disk?\n\n{path}\n→ {target}\n\n"
            "The move is logged in LARMOR's rename log; nothing inside the "
            "folder is changed.",
            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel)
        if ans != QMessageBox.Yes:
            return
        try:
            target = aliases.rename_folder(path, name, open_paths=self.open_paths())
        except (aliases.RenameError, OSError) as exc:
            QMessageBox.warning(self, "Rename on disk", str(exc))
            return
        self._retarget(path, str(target))
        self.renamed.emit(path, str(target))


class RenameDialog(QDialog):
    """One small dialog, two choices: a display name in LARMOR only, or the
    folder renamed on disk (an EXPNO must stay a number)."""

    def __init__(self, parent, path: str):
        super().__init__(parent)
        self.path = path
        self._is_expno = aliases.is_expno(path)
        self.setWindowTitle(f"Rename {Path(path).name}")
        v = QVBoxLayout(self)
        v.addWidget(QLabel(f"<b>{Path(path).name}</b><br>"
                           f"<span style='color:#7a8089'>{Path(path).parent}</span>"))
        self.edit = QLineEdit(aliases.display_name(path))
        self.edit.selectAll()
        self.edit.setPlaceholderText("new name")
        v.addWidget(self.edit)
        self.rbAlias = QRadioButton("Display name in LARMOR only")
        self.rbAlias.setToolTip("kept by LARMOR (aliases.json), used by the "
                                "Explorer, the window title, the Datasets dock, "
                                "the recipe's sample and batch labels; the "
                                "folder on disk is untouched. Leave the folder's "
                                "own name to remove an alias.")
        self.rbDisk = QRadioButton("Rename the folder on disk")
        self.rbDisk.setToolTip("os.rename of the folder after a confirmation "
                               "that states both paths; refused when the target "
                               "exists or a file inside is open in LARMOR; "
                               "logged in rename_log.jsonl")
        self.rbAlias.setChecked(True)
        v.addWidget(self.rbAlias)
        v.addWidget(self.rbDisk)
        self.hint = QLabel("")
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color:#7a8089;")
        v.addWidget(self.hint)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
        self.rbDisk.toggled.connect(self._update_hint)
        self.edit.textChanged.connect(self._update_hint)
        self._update_hint()

    def _update_hint(self, *_):
        if self.rbDisk.isChecked():
            if self._is_expno:
                self.hint.setText("An EXPNO folder must stay a number so TopSpin "
                                  "still reads it; for a descriptive name use "
                                  "the display name instead.")
            else:
                self.hint.setText(f"The folder becomes "
                                  f"{Path(self.path).with_name(self.name() or '…')}.")
        else:
            self.hint.setText("Shown everywhere in LARMOR; the folder keeps its "
                              "name on disk.")

    def name(self) -> str:
        return self.edit.text().strip()

    def on_disk(self) -> bool:
        return self.rbDisk.isChecked()
