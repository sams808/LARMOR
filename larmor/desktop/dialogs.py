"""Dialogs: experiment parameters and dmfit-style parameter links."""
from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
    QLabel, QLineEdit, QVBoxLayout,
)


class ExperimentDialog(QDialog):
    """Edit nucleus / Larmor frequency / MAS rate (nu_rot) for the recipe."""

    def __init__(self, parent, recipe: dict):
        super().__init__(parent)
        self.setWindowTitle("Experiment parameters")
        self.recipe = recipe
        form = QFormLayout(self)

        self.nucleus = QLineEdit(recipe.get("nucleus", ""))
        self.nucleus.setPlaceholderText("e.g. 27Al, 23Na, 29Si")
        form.addRow("nucleus", self.nucleus)

        self.larmor = QDoubleSpinBox()
        self.larmor.setDecimals(4)
        self.larmor.setRange(0.1, 2000.0)
        self.larmor.setSuffix(" MHz")
        self.larmor.setValue(recipe.get("larmor_frequency_MHz", 100.0) or 100.0)
        form.addRow("Larmor frequency", self.larmor)

        self.mas = QDoubleSpinBox()
        self.mas.setDecimals(1)
        self.mas.setRange(0.0, 300000.0)
        self.mas.setSuffix(" Hz")
        self.mas.setToolTip("MAS spinning rate nu_rot; 0 = static")
        self.mas.setValue(recipe.get("spin_rate_Hz", 0.0) or 0.0)
        form.addRow("MAS rate (νrot)", self.mas)

        note = QLabel("Changing these re-simulates every line (kernels are "
                      "cached per field / spin rate).")
        note.setWordWrap(True)
        note.setStyleSheet("color: #93a0a8;")
        form.addRow(note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def _accept(self):
        nuc = self.nucleus.text().strip()
        if nuc:
            try:
                from mrsimulator.spin_system.isotope import Isotope

                Isotope(symbol=nuc)   # validates
            except Exception:
                self.nucleus.setStyleSheet("border: 1px solid #b0442e;")
                self.nucleus.setToolTip(f"unknown isotope: {nuc!r}")
                return
        self.recipe["nucleus"] = nuc
        self.recipe["larmor_frequency_MHz"] = float(self.larmor.value())
        self.recipe["spin_rate_Hz"] = float(self.mas.value())
        self.accept()


def _other_lines(recipe: dict, exclude: int) -> list[tuple[int, str]]:
    return [(j, s.get("label") or f"s{j}")
            for j, s in enumerate(recipe.get("sites", [])) if j != exclude]


class LinkPositionDialog(QDialog):
    """'This line sits at a fixed offset (ppm or Hz) from another line.'"""

    def __init__(self, parent, recipe: dict, row: int):
        super().__init__(parent)
        self.setWindowTitle("Position relative to another line")
        self.recipe, self.row = recipe, row
        self.expr: str | None = None
        v = QVBoxLayout(self)
        form = QFormLayout()

        self.ref = QComboBox()
        for j, label in _other_lines(recipe, row):
            self.ref.addItem(f"s{j} — {label}", j)
        form.addRow("reference line", self.ref)

        self.offset = QDoubleSpinBox()
        self.offset.setDecimals(4)
        self.offset.setRange(-1e7, 1e7)
        form.addRow("offset", self.offset)

        self.unit = QComboBox()
        self.unit.addItems(["ppm", "Hz"])
        form.addRow("unit", self.unit)
        v.addLayout(form)

        note = QLabel("Hz offsets are converted with the recipe's Larmor "
                      "frequency. The position follows the reference during "
                      "the fit, with error propagation.")
        note.setWordWrap(True)
        note.setStyleSheet("color: #93a0a8;")
        v.addWidget(note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        v.addWidget(buttons)

    def _accept(self):
        if self.ref.currentIndex() < 0:
            self.reject()
            return
        j = self.ref.currentData()
        off = float(self.offset.value())
        if self.unit.currentText() == "Hz":
            larmor = self.recipe.get("larmor_frequency_MHz", 0.0) or 1.0
            off = off / larmor          # Hz -> ppm
        self.expr = f"s{j}.isotropic_chemical_shift_ppm + {off:.6g}"
        self.accept()


class RatioDialog(QDialog):
    """'This line's <param> is a fixed multiple of another line's.'"""

    def __init__(self, parent, recipe: dict, row: int, param: str,
                 title: str, default_ratio: float = 1.0):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.expr: str | None = None
        self.param = param
        v = QVBoxLayout(self)
        form = QFormLayout()
        self.ref = QComboBox()
        for j, label in _other_lines(recipe, row):
            self.ref.addItem(f"s{j} — {label}", j)
        form.addRow("reference line", self.ref)
        self.ratio = QDoubleSpinBox()
        self.ratio.setDecimals(6)
        self.ratio.setRange(1e-6, 1e6)
        self.ratio.setValue(default_ratio)
        form.addRow("ratio", self.ratio)
        v.addLayout(form)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        v.addWidget(buttons)

    def _accept(self):
        if self.ref.currentIndex() < 0:
            self.reject()
            return
        j = self.ref.currentData()
        r = float(self.ratio.value())
        self.expr = (f"s{j}.{self.param}" if r == 1.0
                     else f"{r:.6g} * s{j}.{self.param}")
        self.accept()
