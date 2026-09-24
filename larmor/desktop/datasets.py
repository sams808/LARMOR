"""Multi-dataset cockpit: overlay / stack / compare several spectra and pick
which one is the active fit target (ssNake multiplot / TopSpin multiple display).

The panel is a thin view over a list of dataset dicts owned by the main window;
it emits intent signals and never touches the data itself. Each overlay row
carries that spectrum's display controls -- colour, a scale factor (x), an
x shift in ppm and a y offset (a fraction of the active spectrum's span) --
which the main window applies at draw time (larmor.display.overlay_display);
the stored arrays are never modified. A refresh that changes only values
(a colour, a factor, a visibility) updates the rows in place, so the spin
box being edited keeps its focus.
"""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractSpinBox, QCheckBox, QDoubleSpinBox, QHBoxLayout, QLabel, QMenu,
    QPushButton, QScrollArea, QVBoxLayout, QWidget,
)

from larmor import display
from larmor.desktop import theme

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
    compare_requested = Signal()           # acqus/procs of active + overlays
    match_changed = Signal(bool)           # fill the scales with the matching factor
    scale_changed = Signal(int, float)     # overlay i: display scale factor (x)
    shift_changed = Signal(int, float)     # overlay i: x shift, ppm
    yoff_changed = Signal(int, float)      # overlay i: y offset, fraction of the active span
    reset_requested = Signal(int)          # overlay i: scale 1, shift 0, offset 0

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
        self.btnCompare = QPushButton("Compare acquisition…")
        self.btnCompare.setToolTip(
            "compare acqus / procs of the active spectrum and every overlay — "
            "window, LB, TDeff, SI, phase, D1, NS, pulse… (Bruker sources only)")
        self.btnCompare.setEnabled(False)       # until an overlay with a source exists
        self.btnCompare.clicked.connect(self.compare_requested)
        head.addWidget(self.btnCompare)
        self._v.addLayout(head)

        off = QHBoxLayout()
        off.addWidget(QLabel("stack offset"))
        self.offset = QDoubleSpinBox()
        self.offset.setRange(0.0, 5.0); self.offset.setSingleStep(0.1)
        self.offset.setDecimals(2)
        self.offset.setToolTip("shift each overlay up by this fraction of the "
                               "active spectrum's span for a stacked look "
                               "(0 = overlaid); each row's own ↑ offset adds to it")
        self.offset.valueChanged.connect(self.offset_changed)
        off.addWidget(self.offset)
        self.match = QCheckBox("match height")
        self.match.setToolTip("fill each overlay's × scale with the factor that "
                              "brings its maximum to the active spectrum's — one "
                              "click, editable afterwards; untick for ×1 "
                              "(display only; nothing is written)")
        self.match.toggled.connect(self.match_changed)
        off.addWidget(self.match); off.addStretch(1)
        self._v.addLayout(off)

        self._rows = QVBoxLayout()
        self._v.addLayout(self._rows)
        self._v.addStretch(1)
        self._row_widgets: list[dict] = []     # one dict of widgets per overlay row
        self._active_label: QLabel | None = None
        self._colors: list[str] = []

    # ------------------------------------------------------------- texts
    @staticmethod
    def _detail_of(ov: dict) -> str:
        bits = [ov.get("nucleus", "")]
        if ov.get("larmor_MHz"):
            bits.append(f"{float(ov['larmor_MHz']):.1f} MHz")
        if ov.get("npts"):
            bits.append(f"{int(ov['npts'])} pts")
        if ov.get("title"):
            bits.append(str(ov["title"]))
        if ov.get("comparability"):          # differs from the compared set's majority
            bits.append("⚠ " + str(ov["comparability"]))
        return " · ".join(b for b in bits if b)

    @staticmethod
    def _badge_of(ov: dict) -> str:
        """'×2.5 · +1.2 ppm · ↑0.3' when the overlay is drawn transformed."""
        return display.overlay_badge(ov.get("scale", 1.0), ov.get("shift", 0.0),
                                     ov.get("yoff", 0.0))

    def _label_html(self, ov: dict, t) -> str:
        html = f"{ov['label']}"
        badge = self._badge_of(ov)
        if badge:
            html += (f" <span style='color:{t.accent}; font-size:10px;'>"
                     f"{badge}</span>")
        detail = self._detail_of(ov)
        if detail:
            html += (f"<br><span style='color:{t.text_dim}; font-size:10px;'>"
                     f"{detail}</span>")
        return html

    @staticmethod
    def _tip_of(ov: dict) -> str:
        tip = ov.get("source", "")
        if ov.get("title"):
            tip = (tip + "\n" if tip else "") + str(ov["title"])
        if ov.get("comparability"):
            tip = (tip + "\n" if tip else "") + "⚠ " + str(ov["comparability"])
        return tip

    @staticmethod
    def _active_html(active_label: str, active_detail: str, t) -> str:
        return (f"● active: <b>{active_label or '(none)'}</b>"
                + (f"<br><span style='color:{t.text_dim}; font-size:10px;'>"
                   f"{active_detail}</span>" if active_detail else ""))

    def _pick_color(self, i: int, current: str):
        from PySide6.QtGui import QColor
        from PySide6.QtWidgets import QColorDialog
        c = QColorDialog.getColor(QColor(current), self, "Overlay color")
        if c.isValid():
            self.color_changed.emit(i, c.name())

    # ------------------------------------------------------------- rows
    def rebuild(self, active_label: str, overlays: list[dict],
                active_detail: str = ""):
        """Show the active spectrum and one row per overlay. When the same
        overlays are listed (same sources and labels, same order) the
        existing rows are updated in place -- a value change never steals
        the focus of the spin box being edited."""
        self.btnCompare.setEnabled(any(ov.get("source") for ov in overlays))
        keys = [(ov.get("source", ""), ov.get("label", "")) for ov in overlays]
        if (overlays and self._active_label is not None
                and keys == [w["key"] for w in self._row_widgets]):
            self._update_in_place(active_label, overlays, active_detail)
            return
        self._clear_rows()
        t = theme.active()
        active = QLabel(self._active_html(active_label, active_detail, t))
        active.setTextFormat(Qt.RichText)
        active.setStyleSheet(f"color: {t.text};")
        active.setWordWrap(True)
        active.setToolTip("the spectrum currently being fitted")
        self._rows.addWidget(active)
        self._active_label = active

        if not overlays:
            hint = QLabel("no comparison spectra yet — “＋ Add spectrum to "
                          "compare…” overlays any 1D file on the active one")
            hint.setStyleSheet(f"color: {t.text_dim}; font-size: 10px;")
            hint.setWordWrap(True)
            self._rows.addWidget(hint)
            return

        for i, ov in enumerate(overlays):
            self._colors.append(ov["color"])
            self._rows.addWidget(self._build_row(i, ov, t))

    def _clear_rows(self):
        while self._rows.count():
            item = self._rows.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self._row_widgets = []
        self._active_label = None
        self._colors = []

    @staticmethod
    def _swatch_style(color: str, t) -> str:
        return (f"background: {color}; border: 1px solid {t.border}; "
                "border-radius: 3px;")

    @staticmethod
    def _set_quiet(spin: QDoubleSpinBox, value: float):
        spin.blockSignals(True)
        spin.setValue(float(value))
        spin.blockSignals(False)

    def _update_in_place(self, active_label: str, overlays: list[dict],
                         active_detail: str):
        t = theme.active()
        self._active_label.setText(self._active_html(active_label, active_detail, t))
        for i, (w, ov) in enumerate(zip(self._row_widgets, overlays)):
            self._colors[i] = ov["color"]
            w["chk"].blockSignals(True)
            w["chk"].setChecked(ov.get("visible", True))
            w["chk"].blockSignals(False)
            w["swatch"].setStyleSheet(self._swatch_style(ov["color"], t))
            w["lab"].setText(self._label_html(ov, t))
            w["lab"].setToolTip(self._tip_of(ov))
            self._set_quiet(w["scale"], ov.get("scale", 1.0))
            self._set_quiet(w["shift"], ov.get("shift", 0.0))
            self._set_quiet(w["yoff"], ov.get("yoff", 0.0))
            w["act"].setEnabled(bool(ov.get("source")))
            w["source"] = ov.get("source", "")

    def _build_row(self, i: int, ov: dict, t) -> QWidget:
        row = QWidget()
        v = QVBoxLayout(row)
        v.setContentsMargins(0, 2, 0, 2)
        v.setSpacing(1)
        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        chk = QCheckBox()
        chk.setChecked(ov.get("visible", True))
        chk.setToolTip("show / hide this overlay in the plot")
        chk.toggled.connect(lambda on, i=i: self.visibility_changed.emit(i, on))
        top.addWidget(chk)
        swatch = QPushButton()
        swatch.setFixedSize(18, 18)
        swatch.setToolTip("click to change this overlay's color")
        swatch.setStyleSheet(self._swatch_style(ov["color"], t))
        swatch.clicked.connect(
            lambda _=False, i=i: self._pick_color(i, self._colors[i]))
        top.addWidget(swatch)
        lab = QLabel(self._label_html(ov, t))
        lab.setTextFormat(Qt.RichText)
        lab.setToolTip(self._tip_of(ov))
        top.addWidget(lab, 1)
        act = QPushButton("active")
        act.setToolTip("make this the spectrum being fitted")
        act.setEnabled(bool(ov.get("source")))
        act.clicked.connect(lambda _=False, i=i: self.make_active.emit(i))
        top.addWidget(act)
        rm = QPushButton("✕"); rm.setMaximumWidth(28)
        rm.setToolTip("remove this overlay (the file is untouched)")
        rm.clicked.connect(lambda _=False, i=i: self.remove.emit(i))
        top.addWidget(rm)
        v.addLayout(top)

        # the display transform of this overlay: scale x, shift ppm, y offset
        # (values are applied on Enter / focus-out / a step, not per keystroke)
        ctl = QHBoxLayout()
        ctl.setContentsMargins(24, 0, 0, 0)
        scale = QDoubleSpinBox()
        scale.setPrefix("× ")
        scale.setRange(0.01, 1000.0)
        scale.setDecimals(3)
        scale.setStepType(QAbstractSpinBox.AdaptiveDecimalStepType)   # log-friendly steps
        scale.setKeyboardTracking(False)
        scale.setValue(float(ov.get("scale", 1.0)))
        scale.setToolTip("multiply this overlay's intensity for display "
                         "(match height fills it in); the stored data are untouched")
        scale.valueChanged.connect(lambda v, i=i: self.scale_changed.emit(i, float(v)))
        ctl.addWidget(scale)
        shift = QDoubleSpinBox()
        shift.setSuffix(" ppm")
        shift.setRange(-1e6, 1e6)
        shift.setDecimals(2)
        shift.setSingleStep(0.1)
        shift.setKeyboardTracking(False)
        shift.setValue(float(ov.get("shift", 0.0)))
        shift.setToolTip("shift this overlay along the axis, e.g. to line up a "
                         "reference peak")
        shift.valueChanged.connect(lambda v, i=i: self.shift_changed.emit(i, float(v)))
        ctl.addWidget(shift)
        yoff = QDoubleSpinBox()
        yoff.setPrefix("↑ ")
        yoff.setRange(-5.0, 5.0)
        yoff.setDecimals(2)
        yoff.setSingleStep(0.05)
        yoff.setKeyboardTracking(False)
        yoff.setValue(float(ov.get("yoff", 0.0)))
        yoff.setToolTip("raise (or lower) this overlay by a fraction of the active "
                        "spectrum's span, on top of the global stack offset")
        yoff.valueChanged.connect(lambda v, i=i: self.yoff_changed.emit(i, float(v)))
        ctl.addWidget(yoff)
        ctl.addStretch(1)
        v.addLayout(ctl)

        row.setContextMenuPolicy(Qt.CustomContextMenu)
        row.customContextMenuRequested.connect(
            lambda pos, i=i, r=row: self.row_menu(i).exec(r.mapToGlobal(pos)))
        self._row_widgets.append({
            "key": (ov.get("source", ""), ov.get("label", "")), "row": row,
            "chk": chk, "swatch": swatch, "lab": lab, "scale": scale,
            "shift": shift, "yoff": yoff, "act": act,
            "source": ov.get("source", "")})
        return row

    def row_menu(self, i: int) -> QMenu:
        """The right-click menu of overlay row i: Reset scale / shift /
        offset, Make active (with a source), Remove. Built fresh each time;
        the row's context-menu handler exec()s it."""
        m = QMenu(self)
        reset = m.addAction("Reset scale / shift / offset")
        reset.setToolTip("draw this overlay as stored again: ×1, no shift, no offset")
        reset.triggered.connect(lambda _=False, i=i: self.reset_requested.emit(i))
        act = m.addAction("Make active")
        act.setToolTip("make this the spectrum being fitted")
        has_source = bool(i < len(self._row_widgets) and self._row_widgets[i]["source"])
        act.setEnabled(has_source)
        act.triggered.connect(lambda _=False, i=i: self.make_active.emit(i))
        m.addSeparator()
        rm = m.addAction("Remove")
        rm.setToolTip("remove this overlay (the file is untouched)")
        rm.triggered.connect(lambda _=False, i=i: self.remove.emit(i))
        m.setToolTipsVisible(True)
        return m
