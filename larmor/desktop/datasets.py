"""Multi-dataset cockpit: overlay / stack / compare several spectra and pick
which one is the active fit target (ssNake multiplot / TopSpin multiple display).

The panel is a thin view over a list of dataset dicts owned by the main window;
it emits intent signals and never touches the data itself.
"""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal

from larmor.desktop import theme
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
    color_changed = Signal(int, str)       # overlay i gets a new hex color

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

    @staticmethod
    def _detail_of(ov: dict) -> str:
        bits = [ov.get("nucleus", "")]
        if ov.get("larmor_MHz"):
            bits.append(f"{float(ov['larmor_MHz']):.1f} MHz")
        if ov.get("npts"):
            bits.append(f"{int(ov['npts'])} pts")
        if ov.get("title"):
            bits.append(str(ov["title"]))
        return " · ".join(b for b in bits if b)

    def _pick_color(self, i: int, current: str):
        from PySide6.QtGui import QColor
        from PySide6.QtWidgets import QColorDialog
        c = QColorDialog.getColor(QColor(current), self, "Overlay color")
        if c.isValid():
            self.color_changed.emit(i, c.name())

    def rebuild(self, active_label: str, overlays: list[dict],
                active_detail: str = ""):
        while self._rows.count():
            item = self._rows.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        t = theme.active()
        active = QLabel(f"● active: <b>{active_label or '(none)'}</b>"
                        + (f"<br><span style='color:{t.text_dim}; "
                           f"font-size:10px;'>{active_detail}</span>"
                           if active_detail else ""))
        active.setTextFormat(Qt.RichText)
        active.setStyleSheet(f"color: {t.text};")
        active.setWordWrap(True)
        active.setToolTip("the spectrum currently being fitted")
        self._rows.addWidget(active)

        if not overlays:
            hint = QLabel("no comparison spectra yet — “＋ Add spectrum to "
                          "compare…” overlays any 1D file on the active one")
            hint.setStyleSheet(f"color: {t.text_dim}; font-size: 10px;")
            hint.setWordWrap(True)
            self._rows.addWidget(hint)
            return

        for i, ov in enumerate(overlays):
            row = QWidget()
            h = QHBoxLayout(row); h.setContentsMargins(0, 2, 0, 2)
            chk = QCheckBox()
            chk.setChecked(ov.get("visible", True))
            chk.setToolTip("show / hide this overlay in the plot")
            chk.toggled.connect(lambda on, i=i: self.visibility_changed.emit(i, on))
            h.addWidget(chk)
            swatch = QPushButton()
            swatch.setFixedSize(18, 18)
            swatch.setToolTip("click to change this overlay's color")
            swatch.setStyleSheet(
                f"background: {ov['color']}; border: 1px solid {t.border}; "
                "border-radius: 3px;")
            swatch.clicked.connect(
                lambda _=False, i=i, c=ov["color"]: self._pick_color(i, c))
            h.addWidget(swatch)
            detail = self._detail_of(ov)
            lab = QLabel(f"{ov['label']}"
                         + (f"<br><span style='color:{t.text_dim}; "
                            f"font-size:10px;'>{detail}</span>" if detail
                            else ""))
            lab.setTextFormat(Qt.RichText)
            tip = ov.get("source", "")
            if ov.get("title"):
                tip = (tip + "\n" if tip else "") + str(ov["title"])
            lab.setToolTip(tip)
            h.addWidget(lab, 1)
            act = QPushButton("active")
            act.setToolTip("make this the spectrum being fitted")
            act.setEnabled(bool(ov.get("source")))
            act.clicked.connect(lambda _=False, i=i: self.make_active.emit(i))
            h.addWidget(act)
            rm = QPushButton("✕"); rm.setMaximumWidth(28)
            rm.setToolTip("remove this overlay (the file is untouched)")
            rm.clicked.connect(lambda _=False, i=i: self.remove.emit(i))
            h.addWidget(rm)
            self._rows.addWidget(row)
