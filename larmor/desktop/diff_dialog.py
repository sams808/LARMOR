"""Side-by-side parameter diff of the current fit vs a reference fit."""
from __future__ import annotations

from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QHeaderView, QLabel, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)

from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS
from larmor.recipe_diff import recipe_diff


class RecipeDiffDialog(QDialog):
    def __init__(self, parent, current: dict, reference: dict, ref_name: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Compare fits")
        self.resize(680, 480)
        v = QVBoxLayout(self)
        v.addWidget(QLabel(f"Current fit vs <b>{ref_name or 'reference'}</b> "
                           "— Δ = current − reference (sites matched by order)."))
        rows = recipe_diff(current, reference)
        t = QTableWidget(len(rows), 5)
        t.setHorizontalHeaderLabels(["line · parameter", "current", "reference",
                                     "Δ", "Δ %"])
        for r, d in enumerate(rows):
            name = f"{d['label']} · {PARAM_LABELS.get(d['param'], d['param'])}"
            cells = [name,
                     "" if d["current"] is None else f"{d['current']:.4g}",
                     "" if d["reference"] is None else f"{d['reference']:.4g}",
                     "" if d["delta"] is None else f"{d['delta']:+.4g}",
                     "" if d["delta_pct"] is None else f"{d['delta_pct']:+.1f}%"]
            for c, txt in enumerate(cells):
                it = QTableWidgetItem(txt)
                if c == 3 and d["delta"] is not None and abs(d["delta"]) > 0:
                    # tint big relative disagreements
                    big = d["delta_pct"] is not None and abs(d["delta_pct"]) > 10
                    it.setForeground(QColor("#c0392b" if big
                                            else theme.active().text))
                t.setItem(r, c, it)
        t.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        t.verticalHeader().setVisible(False)
        v.addWidget(t, 1)
        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        v.addWidget(bb)
