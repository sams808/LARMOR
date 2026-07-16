"""dmfit-style fit-parameters table: one ROW per line, one column per
parameter, a pin checkbox beside every value (checked = FIXED, like dmfit),
colored row headers matching the curves, and Compute / Fit / χ² in the
footer. This is the primary model editor, faithful to dmfit's 'Fit
Parameters' window."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QAction, QColor
from PySide6.QtWidgets import (
    QCheckBox, QDoubleSpinBox, QHBoxLayout, QHeaderView, QInputDialog,
    QLabel, QMenu, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from larmor.desktop.plot import site_color

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
    """value spinbox + pin checkbox (checked = FIXED, dmfit convention)."""

    edited = Signal()
    pinned = Signal()

    def __init__(self, p: dict):
        super().__init__()
        self.p = p
        lay = QHBoxLayout(self)
        lay.setContentsMargins(2, 0, 2, 0)
        lay.setSpacing(2)
        self.spin = QDoubleSpinBox()
        self.spin.setDecimals(4)
        self.spin.setRange(-1e12, 1e12)
        self.spin.setKeyboardTracking(False)
        self.spin.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
        self.spin.setButtonSymbols(QDoubleSpinBox.NoButtons)
        self.spin.setFrame(False)
        self.spin.setValue(p["value"])
        self.spin.setEnabled(not p.get("expr"))
        self.pin = QCheckBox()
        self.pin.setToolTip("pin: hold this parameter fixed during the fit")
        self.pin.setChecked(not p.get("vary", True) and not p.get("expr"))
        self.pin.setEnabled(not p.get("expr"))
        if p.get("expr"):
            self.spin.setToolTip("linked: " + p["expr"])
            self.setStyleSheet("background: #e2f0f0;")
        err = p.get("stderr")
        if err:
            self.spin.setToolTip((p.get("expr") or "") + f"  ± {err:.3g}".strip())
        lay.addWidget(self.spin, 1)
        lay.addWidget(self.pin)
        self.spin.valueChanged.connect(self._on_value)
        self.pin.toggled.connect(self._on_pin)

    def _on_value(self, v):
        self.p["value"] = float(v)
        self.edited.emit()

    def _on_pin(self, checked):
        self.p["vary"] = not checked
        self.pinned.emit()

    def wheelEvent(self, ev):  # scroll on the cell nudges the value
        if not self.spin.isEnabled():
            return
        step = (abs(self.p["value"]) or 1.0) * (0.1 if ev.modifiers() & Qt.ShiftModifier else 0.02)
        self.p["value"] += step if ev.angleDelta().y() > 0 else -step
        self.spin.blockSignals(True)
        self.spin.setValue(self.p["value"])
        self.spin.blockSignals(False)
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
        self.hint = QLabel("pin ☑ = fixed · right-click a cell for link / bounds · scroll a cell to nudge")
        self.hint.setStyleSheet("color: #93a0a8; font-size: 10px;")
        foot.addWidget(self.hint)
        v.addLayout(foot)
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
        for i, site in enumerate(sites):
            head = QTableWidgetItem(f"■ {site.get('label') or f's{i}'}")
            head.setForeground(QColor(site_color(i)))
            if i in hidden:
                head.setForeground(QColor("#b9c1bc"))
            head.setToolTip("double-click to rename; right-click for actions")
            t.setItem(i, 0, head)
            model_item = QTableWidgetItem(site["model"])
            model_item.setFlags(Qt.ItemIsEnabled)
            model_item.setForeground(QColor("#5a6871"))
            t.setItem(i, 1, model_item)
            for c, key in enumerate(self._used_keys, start=2):
                if key in site["params"]:
                    cell = _Cell(site["params"][key])
                    cell.edited.connect(self.edited)
                    cell.pinned.connect(self.edited)
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

    def _renamed(self, item):
        if self._recipe and item.column() == 0:
            i = item.row()
            if i < len(self._recipe["sites"]):
                self._recipe["sites"][i]["label"] = item.text().lstrip("■ ").strip()

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
                a_bounds = QAction("Set min / max…", menu)
                a_bounds.triggered.connect(lambda: self._edit_bounds(p))
                menu.addAction(a_bounds)
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

    def _edit_bounds(self, p: dict):
        lo, ok1 = QInputDialog.getText(self, "Bounds", "min (empty = none):",
                                       text="" if p.get("min") is None else str(p["min"]))
        if not ok1:
            return
        hi, ok2 = QInputDialog.getText(self, "Bounds", "max (empty = none):",
                                       text="" if p.get("max") is None else str(p["max"]))
        if not ok2:
            return
        p["min"] = float(lo) if lo.strip() else None
        p["max"] = float(hi) if hi.strip() else None
        self.constraint_edited.emit()
