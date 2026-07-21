"""dmfit-style fit-parameters table: one ROW per line, one column per
parameter, a pin checkbox beside every value (checked = FIXED, like dmfit),
colored row headers matching the curves, and Compute / Fit / χ² in the
footer. This is the primary model editor, faithful to dmfit's 'Fit
Parameters' window."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QAction, QColor
from PySide6.QtWidgets import (
    QCheckBox, QHBoxLayout, QHeaderView, QInputDialog, QLabel, QLineEdit,
    QMenu, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from larmor import cellparse
from larmor.desktop.plot import site_color


def _param_unit(model: str, key: str) -> str:
    """The declared unit of a parameter (from the model registry)."""
    try:
        from larmor import models

        for pd in models.get(model).params:
            if pd.name == key:
                return pd.unit or ""
    except Exception:
        pass
    return ""

#: column order and short headers, dmfit-style
PARAM_COLUMNS = [
    ("amplitude", "Amplitude"),
    ("isotropic_chemical_shift_ppm", "Position\n(ppm)"),
    ("shift_fwhm_ppm", "Width\n(ppm)"),
    ("gl", "xG/(1-x)L"),
    ("sigma_Cq_MHz", "σ(Cq)\n(MHz)"),
    ("Cq_MHz", "Cq\n(MHz)"),
    ("eta", "η"),
    ("zeta_ppm", "ζ CSA\n(ppm)"),
]


class _Cell(QWidget):
    """A smart value field + pin checkbox (checked = FIXED, dmfit convention).

    The field is a text box: type a number, a link to another component by its
    LETTER (A+20, A+20kHz, 0.5B), or bounds ([0..100]); see larmor.cellparse.
    """

    edited = Signal()
    pinned = Signal()
    error = Signal(str)

    def __init__(self, p: dict, param_name: str, param_unit: str,
                 this_index: int, n_sites: int, larmor_MHz: float):
        super().__init__()
        self.p = p
        self.ctx = dict(param_name=param_name, param_unit=param_unit,
                        this_index=this_index, n_sites=n_sites,
                        larmor_MHz=larmor_MHz)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(2, 0, 2, 0)
        lay.setSpacing(2)
        self.edit = QLineEdit()
        self.edit.setFrame(False)
        self.edit.setText(self._display_text())
        self.pin = QCheckBox()
        self.pin.setToolTip("pin: hold this parameter fixed during the fit")
        self.pin.setChecked(not p.get("vary", True) and not p.get("expr"))
        self.pin.setEnabled(not p.get("expr"))
        self._style()
        lay.addWidget(self.edit, 1)
        lay.addWidget(self.pin)
        self.edit.editingFinished.connect(self._on_edit)
        self.pin.toggled.connect(self._on_pin)

    # ---- display ----
    def _display_text(self) -> str:
        if self.p.get("expr"):
            return cellparse.format_link(self.p["expr"], self.ctx["param_name"])
        v = self.p["value"]
        return f"{v:.5g}" if abs(v) < 1e5 or v == 0 else f"{v:.6g}"

    def _style(self):
        p = self.p
        tips = []
        bounded = p.get("min") is not None or p.get("max") is not None
        css = ""
        if bounded:
            lo = "−∞" if p.get("min") is None else f"{p['min']:g}"
            hi = "+∞" if p.get("max") is None else f"{p['max']:g}"
            tips.append(f"constrained to [{lo}, {hi}]  ·  edit inline as [{lo}..{hi}]")
            css += "QLineEdit { border-left: 3px solid #0e7c86; }"
        if p.get("expr"):
            tips.append("linked: " + p["expr"] + "  ·  type a number to unlink")
            css += "QLineEdit { background: #e2f0f0; }"
        else:
            tips.append("type a value, a link (A+20, A+20kHz, 0.5B), "
                        "or bounds [0..100]")
        if p.get("stderr"):
            tips.append(f"± {p['stderr']:.3g}")
        self.edit.setStyleSheet(css)
        self.edit.setToolTip("  ·  ".join(tips))

    # ---- edits ----
    def _on_edit(self):
        res = cellparse.parse_cell(self.edit.text(), **self.ctx)
        if res.error:
            self.error.emit(res.error)
            self.edit.setText(self._display_text())     # revert
            return
        if res.set_value:
            self.p["value"] = res.value
        if res.set_expr:
            self.p["expr"] = res.expr
            if res.expr:
                self.p["vary"] = True
        if res.set_min:
            self.p["min"] = res.min
        if res.set_max:
            self.p["max"] = res.max
        # clamp the stored value into any new bounds
        if self.p.get("min") is not None and self.p["value"] < self.p["min"]:
            self.p["value"] = self.p["min"]
        if self.p.get("max") is not None and self.p["value"] > self.p["max"]:
            self.p["value"] = self.p["max"]
        self.pin.setEnabled(not self.p.get("expr"))
        self.edit.setText(self._display_text())
        self._style()
        self.edited.emit()

    def _on_pin(self, checked):
        self.p["vary"] = not checked
        self.pinned.emit()

    def wheelEvent(self, ev):  # scroll on the cell nudges the value
        if self.p.get("expr"):
            return
        step = (abs(self.p["value"]) or 1.0) * (
            0.1 if ev.modifiers() & Qt.ShiftModifier else 0.02)
        self.p["value"] += step if ev.angleDelta().y() > 0 else -step
        lo, hi = self.p.get("min"), self.p.get("max")
        if lo is not None:
            self.p["value"] = max(self.p["value"], lo)
        if hi is not None:
            self.p["value"] = min(self.p["value"], hi)
        self.edit.setText(self._display_text())
        self.edited.emit()


class LinesTable(QWidget):
    """The bottom-dock spreadsheet + Compute/Fit footer."""

    edited = Signal()            # any value changed -> resimulate
    structure = Signal(int, str) # (row, "remove"|"duplicate"|"visibility")
    compute = Signal()
    fit = Signal()
    constraint_edited = Signal()

    def __init__(self):
        super().__init__()
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 2, 4, 2)
        v.setSpacing(3)
        self.table = QTableWidget(0, 0)
        self.table.verticalHeader().setDefaultSectionSize(28)
        self.table.verticalHeader().setFixedWidth(24)
        self.table.setAlternatingRowColors(True)
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._context_menu)
        v.addWidget(self.table, 1)

        foot = QHBoxLayout()
        self.btnCompute = QPushButton("Compute")
        self.btnCompute.setToolTip("re-simulate the model with current values")
        self.btnCompute.clicked.connect(self.compute)
        self.btnFit = QPushButton("Fit")
        self.btnFit.setStyleSheet("font-weight: 600;")
        self.btnFit.clicked.connect(self.fit)
        self.chi2 = QLabel("")
        self.chi2.setStyleSheet("font-weight: 600; color: #0a5a62;")
        foot.addWidget(self.btnCompute)
        foot.addWidget(self.btnFit)
        foot.addSpacing(16)
        foot.addWidget(self.chi2)
        foot.addStretch(1)
        v.addLayout(foot)
        # a full-width, always-visible hint (was previously cut off at the edge)
        self.hint = QLabel(
            "Type in a cell:  a value  ·  a bound  [0..100]  ·  a link to "
            "another line by its letter:  A  ·  A+20  ·  A+20kHz (→ppm)  ·  "
            "0.5B  ·  A+20 [50..80].   pin ☑ = fixed   ·   scroll = nudge   ·   "
            "right-click for menus")
        self.hint.setStyleSheet("color: #5a6871; font-size: 10px; padding: 2px 4px;")
        self.hint.setWordWrap(True)
        v.addWidget(self.hint)
        self._recipe: dict | None = None

    # ------------------------------------------------------------------
    def rebuild(self, recipe: dict | None, hidden: set[int]):
        self._recipe = recipe
        t = self.table
        t.blockSignals(True)
        t.clear()
        sites = (recipe or {}).get("sites", [])
        used = [(k, h) for k, h in PARAM_COLUMNS
                if any(k in s["params"] for s in sites)]
        self._used_keys = [k for k, _ in used]
        t.setColumnCount(2 + len(used))
        t.setHorizontalHeaderLabels(["line", "model"] + [h for _, h in used])
        t.setRowCount(len(sites))
        larmor = (recipe or {}).get("larmor_frequency_MHz", 0.0)
        n = len(sites)
        # vertical header shows the component LETTER (A, B, C…), dmfit-style
        t.setVerticalHeaderLabels([cellparse.index_to_letter(i)
                                   for i in range(n)])
        for i, site in enumerate(sites):
            letter = cellparse.index_to_letter(i)
            head = QTableWidgetItem(f"■ {letter} · {site.get('label') or ''}".rstrip(" ·"))
            head.setForeground(QColor(site_color(i)))
            if i in hidden:
                head.setForeground(QColor("#b9c1bc"))
            head.setToolTip("double-click to rename; right-click for actions. "
                            f"Reference this line in a link as “{letter}”.")
            t.setItem(i, 0, head)
            model_item = QTableWidgetItem(site["model"])
            model_item.setFlags(Qt.ItemIsEnabled)
            model_item.setForeground(QColor("#5a6871"))
            t.setItem(i, 1, model_item)
            for c, key in enumerate(self._used_keys, start=2):
                if key in site["params"]:
                    cell = _Cell(site["params"][key], key,
                                 _param_unit(site["model"], key), i, n, larmor)
                    cell.edited.connect(self.edited)
                    cell.pinned.connect(self.edited)
                    cell.error.connect(self._cell_error)
                    t.setCellWidget(i, c, cell)
                else:
                    blank = QTableWidgetItem("")
                    blank.setFlags(Qt.NoItemFlags)
                    blank.setBackground(QColor("#f3f5f3"))
                    t.setItem(i, c, blank)
        hh = t.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.Interactive)
        hh.setDefaultSectionSize(130)
        hh.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setStretchLastSection(True)
        t.blockSignals(False)
        t.itemChanged.connect(self._renamed)

    def _cell_error(self, msg: str):
        # surface a bad link/bounds entry without crashing the edit
        w = self.window()
        if w is not None and hasattr(w, "statusBar"):
            w.statusBar().showMessage("⚠ " + msg, 6000)

    def _renamed(self, item):
        if self._recipe and item.column() == 0:
            i = item.row()
            if i < len(self._recipe["sites"]):
                # strip the "■ A · " prefix, keep the user's own name
                txt = item.text()
                if "·" in txt:
                    txt = txt.split("·", 1)[1]
                self._recipe["sites"][i]["label"] = txt.strip()

    def set_chi2(self, text: str):
        self.chi2.setText(text)

    # ------------------------------------------------------------------
    def _context_menu(self, pos):
        if not self._recipe:
            return
        row = self.table.rowAt(pos.y())
        col = self.table.columnAt(pos.x())
        if row < 0:
            return
        site = self._recipe["sites"][row]
        menu = QMenu(self)
        if col >= 2 and col - 2 < len(self._used_keys):
            key = self._used_keys[col - 2]
            if key in site["params"]:
                p = site["params"][key]
                others = len(self._recipe["sites"]) > 1
                # dmfit-style presets, no expression writing needed
                if key == "isotropic_chemical_shift_ppm" and others:
                    a = QAction("Position: offset from another line… (ppm/Hz)", menu)
                    a.triggered.connect(lambda: self._preset_position(row, p))
                    menu.addAction(a)
                if key == "amplitude" and others:
                    a = QAction("Amplitude: ratio of another line…", menu)
                    a.triggered.connect(lambda: self._preset_ratio(
                        row, p, "amplitude", "Amplitude ratio", 0.5))
                    menu.addAction(a)
                if key == "shift_fwhm_ppm" and others:
                    a = QAction("Width: same as another line…", menu)
                    a.triggered.connect(lambda: self._preset_ratio(
                        row, p, "shift_fwhm_ppm", "Shared width", 1.0))
                    menu.addAction(a)
                a_link = QAction("Custom link expression…", menu)
                a_link.triggered.connect(lambda: self._edit_link(p))
                menu.addAction(a_link)
                if p.get("expr"):
                    a_un = QAction(f"Unlink  (now: {p['expr']})", menu)
                    a_un.triggered.connect(lambda: self._unlink(p))
                    menu.addAction(a_un)
                bounded = (p.get("min") is not None
                           or p.get("max") is not None)
                label = ("Constrain min / max…  ✓" if bounded
                         else "Constrain min / max…")
                a_bounds = QAction(label, menu)
                a_bounds.triggered.connect(
                    lambda _=False, key=key: self._edit_bounds(site, key))
                menu.addAction(a_bounds)
                if bounded:
                    a_free = QAction("Remove constraints", menu)
                    a_free.triggered.connect(
                        lambda _=False, pp=p: self._clear_bounds(pp))
                    menu.addAction(a_free)
                menu.addSeparator()
        a_vis = QAction("Show / hide on plot", menu)
        a_vis.triggered.connect(lambda: self.structure.emit(row, "visibility"))
        a_dup = QAction("Duplicate line", menu)
        a_dup.triggered.connect(lambda: self.structure.emit(row, "duplicate"))
        a_del = QAction("Remove line", menu)
        a_del.triggered.connect(lambda: self.structure.emit(row, "remove"))
        for a in (a_vis, a_dup, a_del):
            menu.addAction(a)
        menu.exec(self.table.viewport().mapToGlobal(pos))

    def _preset_position(self, row: int, p: dict):
        from larmor.desktop.dialogs import LinkPositionDialog

        dlg = LinkPositionDialog(self, self._recipe, row)
        if dlg.exec() and dlg.expr:
            p["expr"] = dlg.expr
            p["vary"] = True
            self.constraint_edited.emit()

    def _preset_ratio(self, row: int, p: dict, param: str, title: str,
                      default: float):
        from larmor.desktop.dialogs import RatioDialog

        dlg = RatioDialog(self, self._recipe, row, param, title, default)
        if dlg.exec() and dlg.expr:
            p["expr"] = dlg.expr
            p["vary"] = True
            self.constraint_edited.emit()

    def _unlink(self, p: dict):
        p["expr"] = None
        self.constraint_edited.emit()

    def _edit_link(self, p: dict):
        text, ok = QInputDialog.getText(
            self, "Link parameter",
            "expression (empty to unlink), e.g. 0.29 * s0.amplitude:",
            text=p.get("expr") or "")
        if ok:
            p["expr"] = text.strip() or None
            if p["expr"]:
                p["vary"] = True
            self.constraint_edited.emit()

    def _edit_bounds(self, site: dict, key: str):
        from larmor.desktop.dialogs import BoundsDialog

        p = site["params"][key]
        from larmor.desktop.panels import PARAM_LABELS as _PL

        dlg = BoundsDialog(self, _PL.get(key, key), p)
        if dlg.exec():
            p["min"], p["max"] = dlg.result_min, dlg.result_max
            if p["min"] is not None and p["value"] < p["min"]:
                p["value"] = p["min"]
            if p["max"] is not None and p["value"] > p["max"]:
                p["value"] = p["max"]
            self.constraint_edited.emit()

    def _clear_bounds(self, p: dict):
        p["min"] = p["max"] = None
        self.constraint_edited.emit()
