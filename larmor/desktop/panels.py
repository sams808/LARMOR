"""Side panels: sites (dmfit-style parameter cards) and processing."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import (
    QCheckBox, QDoubleSpinBox, QFrame, QGridLayout, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QRadioButton, QScrollArea, QSizePolicy, QSlider,
    QSpinBox, QToolButton, QVBoxLayout, QWidget,
)

from larmor.desktop.plot import site_color

PARAM_LABELS = {
    "isotropic_chemical_shift_ppm": "δiso (ppm)",
    "sigma_Cq_MHz": "σ(Cq) (MHz)",
    "Cq_MHz": "Cq (MHz)",
    "eta": "η",
    "zeta_ppm": "ζ CSA (ppm)",
    "shift_fwhm_ppm": "FWHM (ppm)",
    "amplitude": "amplitude",
    "gl": "g/l fraction",
}


class ParamSpin(QDoubleSpinBox):
    """Spin box tuned for spectroscopy values: wide range, adaptive step."""

    def __init__(self):
        super().__init__()
        self.setDecimals(5)
        self.setRange(-1e12, 1e12)
        self.setKeyboardTracking(False)
        self.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)


class SiteCard(QFrame):
    changed = Signal()          # any parameter/flag edit (snapshot + resim)
    structure = Signal(str)     # "remove" | "duplicate" | "visibility"

    def __init__(self, index: int, site: dict, hidden: bool):
        super().__init__()
        self.index, self.site = index, site
        self.setFrameShape(QFrame.StyledPanel)
        self.setObjectName("siteCard")
        v = QVBoxLayout(self)
        v.setContentsMargins(6, 4, 6, 6)
        v.setSpacing(2)

        head = QHBoxLayout()
        sw = QToolButton()
        sw.setText("■")
        sw.setStyleSheet(f"color: {site_color(index)}; font-size: 14px; border: none;")
        sw.setToolTip("show/hide this site on the plot")
        sw.clicked.connect(lambda: self.structure.emit("visibility"))
        self.name = QLineEdit(site.get("label") or f"s{index}")
        self.name.setFrame(False)
        self.name.setStyleSheet("font-weight: 600;")
        self.name.editingFinished.connect(self._rename)
        tag = QLabel(f"s{index} · {site['model']}")
        tag.setStyleSheet("color: #5a6871; font-size: 10px;")
        bGear = QToolButton(); bGear.setText("⚙")
        bGear.setToolTip("constraints: link expression / min / max")
        bGear.setCheckable(True)
        bGear.toggled.connect(self._toggle_constraints)
        bDup = QToolButton(); bDup.setText("⧉")
        bDup.setToolTip("duplicate site")
        bDup.clicked.connect(lambda: self.structure.emit("duplicate"))
        bDel = QToolButton(); bDel.setText("✕")
        bDel.setToolTip("remove site")
        bDel.clicked.connect(lambda: self.structure.emit("remove"))
        for w in (sw, self.name, tag, bGear, bDup, bDel):
            head.addWidget(w)
        head.setStretch(1, 1)
        v.addLayout(head)

        self.grid = QGridLayout()
        self.grid.setHorizontalSpacing(6)
        self.grid.setVerticalSpacing(2)
        v.addLayout(self.grid)
        self._constraint_rows: list[QWidget] = []
        self._build_rows()
        if hidden:
            self.setStyleSheet("#siteCard { background: #f3f5f3; } * { color: #93a0a8; }")

    def _rename(self):
        self.site["label"] = self.name.text()
        self.changed.emit()

    def _build_rows(self):
        row = 0
        for pname, p in self.site["params"].items():
            lab = QLabel(PARAM_LABELS.get(pname, pname))
            lab.setToolTip(pname + ("  — linked: " + p["expr"] if p.get("expr") else ""))
            if p.get("expr"):
                lab.setText(lab.text() + " ⚭")
                lab.setStyleSheet("color: #0e7c86; font-weight: 600;")
            spin = ParamSpin()
            spin.setValue(p["value"])
            spin.setEnabled(not p.get("expr"))
            spin.valueChanged.connect(
                lambda val, pp=p: (pp.__setitem__("value", float(val)),
                                   self.changed.emit()))
            vary = QCheckBox()
            vary.setChecked(bool(p.get("vary", True)) and not p.get("expr"))
            vary.setEnabled(not p.get("expr"))
            vary.setToolTip("checked = fitted; unchecked = fixed"
                            if not p.get("expr") else "linked — follows its expression")
            vary.toggled.connect(
                lambda on, pp=p: (pp.__setitem__("vary", bool(on)),
                                  self.changed.emit()))
            err = QLabel("± %.3g" % p["stderr"] if p.get("stderr") else "")
            err.setStyleSheet("color: #0a5a62; font-size: 10px;")
            self.grid.addWidget(lab, row, 0)
            self.grid.addWidget(spin, row, 1)
            self.grid.addWidget(vary, row, 2)
            self.grid.addWidget(err, row, 3)
            row += 1

            # constraints row (hidden until ⚙)
            cw = QWidget()
            ch = QHBoxLayout(cw)
            ch.setContentsMargins(0, 0, 0, 2)
            expr = QLineEdit(p.get("expr") or "")
            expr.setPlaceholderText(f"link: 0.5 * s0.{pname}")
            expr.setStyleSheet("font-family: Consolas, monospace; font-size: 10px;")
            expr.editingFinished.connect(
                lambda pp=p, w=expr: self._set_expr(pp, w.text()))
            lo = QLineEdit("" if p.get("min") is None else str(p["min"]))
            lo.setPlaceholderText("min"); lo.setFixedWidth(56)
            lo.setValidator(QDoubleValidator())
            lo.editingFinished.connect(
                lambda pp=p, w=lo: (pp.__setitem__(
                    "min", float(w.text()) if w.text() else None),
                    self.changed.emit()))
            hi = QLineEdit("" if p.get("max") is None else str(p["max"]))
            hi.setPlaceholderText("max"); hi.setFixedWidth(56)
            hi.setValidator(QDoubleValidator())
            hi.editingFinished.connect(
                lambda pp=p, w=hi: (pp.__setitem__(
                    "max", float(w.text()) if w.text() else None),
                    self.changed.emit()))
            ch.addWidget(expr); ch.addWidget(lo); ch.addWidget(hi)
            cw.setVisible(False)
            self.grid.addWidget(cw, row, 0, 1, 4)
            self._constraint_rows.append(cw)
            row += 1

    def _set_expr(self, p: dict, text: str):
        p["expr"] = text.strip() or None
        if p["expr"]:
            p["vary"] = True
        self.changed.emit()

    def _toggle_constraints(self, on: bool):
        for w in self._constraint_rows:
            w.setVisible(on)


class SitesPanel(QScrollArea):
    changed = Signal()
    structure = Signal(int, str)     # site index, action

    def __init__(self):
        super().__init__()
        self.setWidgetResizable(True)
        self._inner = QWidget()
        self._layout = QVBoxLayout(self._inner)
        self._layout.setAlignment(Qt.AlignTop)
        self._layout.setSpacing(6)
        self.setWidget(self._inner)
        self._hint = QLabel(
            "Pick a model in the toolbar, then click on the spectrum to place "
            "a site.\nDrag the dashed marker to move it. Checkbox = fitted; "
            "⚙ = link / bounds.")
        self._hint.setWordWrap(True)
        self._hint.setStyleSheet("color: #93a0a8;")
        self._layout.addWidget(self._hint)

    def rebuild(self, recipe: dict | None, hidden: set[int]):
        while self._layout.count():
            item = self._layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        if not recipe or not recipe.get("sites"):
            self._layout.addWidget(self._hint)
            self._hint = QLabel(self._hint.text())
            self._hint.setWordWrap(True)
            self._hint.setStyleSheet("color: #93a0a8;")
            return
        for i, site in enumerate(recipe["sites"]):
            card = SiteCard(i, site, i in hidden)
            card.changed.connect(self.changed)
            card.structure.connect(
                lambda action, idx=i: self.structure.emit(idx, action))
            self._layout.addWidget(card)


class ProcessingPanel(QWidget):
    apply_requested = Signal(list, bool)   # (ops, use_raw)
    reset_requested = Signal()

    def __init__(self):
        super().__init__()
        v = QVBoxLayout(self)
        v.setAlignment(Qt.AlignTop)

        self.rb_pdata = QRadioButton("TopSpin-processed (pdata)")
        self.rb_raw = QRadioButton("raw fid (EM → ZF → FT)")
        self.rb_pdata.setChecked(True)
        v.addWidget(self.rb_pdata)
        v.addWidget(self.rb_raw)

        raw = QHBoxLayout()
        raw.addWidget(QLabel("EM lb (Hz)"))
        self.lb = QDoubleSpinBox(); self.lb.setRange(0, 1e5); self.lb.setValue(50)
        raw.addWidget(self.lb)
        raw.addWidget(QLabel("ZF ×"))
        self.zf = QSpinBox(); self.zf.setRange(1, 8); self.zf.setValue(2)
        raw.addWidget(self.zf)
        raw.addWidget(QLabel("offset (ppm)"))
        self.off = QDoubleSpinBox(); self.off.setRange(-1e5, 1e5)
        raw.addWidget(self.off)
        v.addLayout(raw)

        v.addWidget(QLabel("<b>Phase</b>"))
        self.btnAuto = QPushButton("Autophase (ACME)")
        v.addWidget(self.btnAuto)
        ph0 = QHBoxLayout()
        ph0.addWidget(QLabel("p0"))
        self.p0 = QSlider(Qt.Horizontal); self.p0.setRange(-180, 180)
        self.p0v = QDoubleSpinBox(); self.p0v.setRange(-180, 180)
        self.p0.valueChanged.connect(self.p0v.setValue)
        self.p0v.valueChanged.connect(lambda v_: self.p0.setValue(int(v_)))
        ph0.addWidget(self.p0); ph0.addWidget(self.p0v)
        v.addLayout(ph0)
        ph1 = QHBoxLayout()
        ph1.addWidget(QLabel("p1"))
        self.p1 = QSlider(Qt.Horizontal); self.p1.setRange(-720, 720)
        self.p1v = QDoubleSpinBox(); self.p1v.setRange(-720, 720)
        self.p1.valueChanged.connect(self.p1v.setValue)
        self.p1v.valueChanged.connect(lambda v_: self.p1.setValue(int(v_)))
        ph1.addWidget(self.p1); ph1.addWidget(self.p1v)
        v.addLayout(ph1)

        bl = QHBoxLayout()
        bl.addWidget(QLabel("<b>Baseline</b> order"))
        self.blOrder = QSpinBox(); self.blOrder.setRange(0, 9); self.blOrder.setValue(3)
        bl.addWidget(self.blOrder)
        self.btnBaseline = QPushButton("Correct baseline")
        bl.addWidget(self.btnBaseline)
        v.addLayout(bl)

        actions = QHBoxLayout()
        self.btnApply = QPushButton("Apply processing")
        self.btnApply.setDefault(True)
        self.btnReset = QPushButton("Reset to original")
        actions.addWidget(self.btnApply); actions.addWidget(self.btnReset)
        v.addLayout(actions)

        note = QLabel("Processing never writes to instrument files — the "
                      "pipeline is applied in memory and the fit uses the result.")
        note.setWordWrap(True)
        note.setStyleSheet("color: #93a0a8;")
        v.addWidget(note)

        self.btnApply.clicked.connect(lambda: self._emit([]))
        self.btnAuto.clicked.connect(lambda: self._emit([{"op": "autophase"}]))
        self.btnBaseline.clicked.connect(
            lambda: self._emit([{"op": "baseline", "order": self.blOrder.value()}]))
        self.btnReset.clicked.connect(self.reset_requested)

    def _emit(self, extra: list[dict]):
        raw = self.rb_raw.isChecked()
        ops: list[dict] = []
        if raw:
            if self.lb.value() > 0:
                ops.append({"op": "em", "lb_hz": self.lb.value()})
            if self.zf.value() > 1:
                ops.append({"op": "zf", "factor": self.zf.value()})
            ops.append({"op": "ft", "offset_ppm": self.off.value()})
        if self.p0v.value() or self.p1v.value():
            ops.append({"op": "phase", "p0": self.p0v.value(), "p1": self.p1v.value()})
        ops.extend(extra)
        self.apply_requested.emit(ops, raw)
