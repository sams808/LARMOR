"""Multi-dataset cockpit: overlay / stack / compare several spectra and pick
which one is the active fit target (ssNake multiplot / TopSpin multiple display).

The panel is a thin view over a list of dataset dicts owned by the main window;
it emits intent signals and never touches the data itself.
"""
from __future__ import annotations

import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QDoubleSpinBox, QHBoxLayout, QLabel, QPushButton, QScrollArea,
    QVBoxLayout, QWidget,
)

#: overlay palette, distinct from the site colors
OVERLAY_COLORS = ["#e8832a", "#1f77b4", "#2ca02c", "#9467bd", "#8c564b",
                  "#17becf", "#bcbd22", "#d62728", "#7f7f7f"]


def overlay_color(i: int) -> str:
    return OVERLAY_COLORS[i % len(OVERLAY_COLORS)]


class DatasetsPanel(QScrollArea):
    add_requested = Signal()               # load another spectrum to compare
    make_active = Signal(int)              # promote overlay i to the fit target
    remove = Signal(int)                   # drop overlay i
    visibility_changed = Signal(int, bool)
    offset_changed = Signal(float)         # global vertical stack offset

    def __init__(self):
        super().__init__()
        self.setWidgetResizable(True)
        self._host = QWidget()
        self.setWidget(self._host)
        self._v = QVBoxLayout(self._host)
        self._v.setAlignment(Qt.AlignTop)

        head = QHBoxLayout()
        self.btnAdd = QPushButton("＋ Add spectrum to compare…")
        self.btnAdd.clicked.connect(self.add_requested)
        head.addWidget(self.btnAdd)
        self._v.addLayout(head)

        off = QHBoxLayout()
        off.addWidget(QLabel("stack offset"))
        self.offset = QDoubleSpinBox()
        self.offset.setRange(0.0, 5.0); self.offset.setSingleStep(0.1)
        self.offset.setDecimals(2)
        self.offset.setToolTip("shift each overlay up by this fraction for a "
                               "stacked look (0 = overlaid)")
        self.offset.valueChanged.connect(self.offset_changed)
        off.addWidget(self.offset); off.addStretch(1)
        self._v.addLayout(off)

        self._rows = QVBoxLayout()
        self._v.addLayout(self._rows)
        self._v.addStretch(1)

    def rebuild(self, active_label: str, overlays: list[dict]):
        while self._rows.count():
            item = self._rows.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        active = QLabel(f"● active: <b>{active_label or '(none)'}</b>")
        active.setStyleSheet("color: #16202a;")
        active.setWordWrap(True)
        self._rows.addWidget(active)

        for i, ov in enumerate(overlays):
            row = QWidget()
            h = QHBoxLayout(row); h.setContentsMargins(0, 0, 0, 0)
            chk = QCheckBox()
            chk.setChecked(ov.get("visible", True))
            chk.toggled.connect(lambda on, i=i: self.visibility_changed.emit(i, on))
            h.addWidget(chk)
            swatch = QLabel("■")
            swatch.setStyleSheet(f"color: {ov['color']}; font-size: 14px;")
            h.addWidget(swatch)
            lab = QLabel(ov["label"])
            lab.setToolTip(ov.get("source", ""))
            h.addWidget(lab, 1)
            act = QPushButton("active")
            act.setToolTip("make this the spectrum being fitted")
            act.setEnabled(bool(ov.get("source")))
            act.clicked.connect(lambda _=False, i=i: self.make_active.emit(i))
            h.addWidget(act)
            rm = QPushButton("✕"); rm.setMaximumWidth(28)
            rm.clicked.connect(lambda _=False, i=i: self.remove.emit(i))
            h.addWidget(rm)
            self._rows.addWidget(row)
