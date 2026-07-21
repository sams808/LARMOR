"""Side panels: sites (dmfit-style parameter cards) and processing."""
from __future__ import annotations

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFrame, QGridLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QRadioButton, QScrollArea, QSizePolicy,
    QSlider, QSpinBox, QToolButton, QVBoxLayout, QWidget,
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
    baseline_mode = Signal(bool)           # pick-anchors toggle
    baseline_apply = Signal()
    baseline_clear = Signal()

    def __init__(self):
        super().__init__()
        v = QVBoxLayout(self)
        v.setAlignment(Qt.AlignTop)

        self.rb_pdata = QRadioButton("TopSpin-processed (pdata)")
        self.rb_raw = QRadioButton("raw fid (EM → ZF → FT)")
        self.rb_pdata.setChecked(True)
        v.addWidget(self.rb_pdata)
        v.addWidget(self.rb_raw)

        # TopSpin-style window function block
        wdw = QHBoxLayout()
        wdw.addWidget(QLabel("WDW"))
        self.wdw = QComboBox()
        self.wdw.addItems(["none", "EM", "GM", "SINE", "QSINE", "TRAF"])
        self.wdw.setCurrentText("EM")
        wdw.addWidget(self.wdw)
        wdw.addWidget(QLabel("LB"))
        self.lb = QDoubleSpinBox(); self.lb.setRange(-1e5, 1e5); self.lb.setValue(50)
        self.lb.setToolTip("Hz; negative for GM (Lorentz-to-Gauss)")
        wdw.addWidget(self.lb)
        wdw.addWidget(QLabel("GB"))
        self.gb = QDoubleSpinBox(); self.gb.setRange(0.001, 1.0); self.gb.setDecimals(3)
        self.gb.setValue(0.1); self.gb.setToolTip("GM: Gaussian max position (0..1)")
        wdw.addWidget(self.gb)
        wdw.addWidget(QLabel("SSB"))
        self.ssb = QDoubleSpinBox(); self.ssb.setRange(0, 64); self.ssb.setValue(2)
        self.ssb.setToolTip("SINE/QSINE: 2 = cosine bell, 0 = pure sine")
        wdw.addWidget(self.ssb)
        v.addLayout(wdw)

        raw = QHBoxLayout()
        raw.addWidget(QLabel("TDeff"))
        self.tdeff = QSpinBox(); self.tdeff.setRange(0, 10_000_000)
        self.tdeff.setToolTip("use only the first TDeff fid points (0 = all)")
        raw.addWidget(self.tdeff)
        raw.addWidget(QLabel("ZF ×"))
        self.zf = QSpinBox(); self.zf.setRange(1, 16); self.zf.setValue(2)
        raw.addWidget(self.zf)
        raw.addWidget(QLabel("FCOR"))
        self.fcor = QDoubleSpinBox(); self.fcor.setRange(0.0, 2.0)
        self.fcor.setDecimals(2); self.fcor.setValue(0.5)
        raw.addWidget(self.fcor)
        raw.addWidget(QLabel("offset (ppm)"))
        self.off = QDoubleSpinBox(); self.off.setRange(-1e5, 1e5)
        raw.addWidget(self.off)
        v.addLayout(raw)

        srrow = QHBoxLayout()
        srrow.addWidget(QLabel("<b>SR</b> (Hz)"))
        self.sr = QDoubleSpinBox(); self.sr.setRange(-1e6, 1e6); self.sr.setDecimals(2)
        self.sr.setToolTip("spectral reference: shifts the ppm axis by SR/SFO1")
        srrow.addWidget(self.sr)
        self.chkMag = QCheckBox("magnitude")
        self.chkMag.setToolTip("phase-insensitive |S| display")
        srrow.addWidget(self.chkMag)
        self.chkHilbert = QCheckBox("Hilbert first")
        self.chkHilbert.setToolTip("rebuild the imaginary part of pdata (1r) "
                                   "so phase correction works on it")
        srrow.addWidget(self.chkHilbert)
        v.addLayout(srrow)

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

        # TopSpin-style quick zero-order phase steps
        quick = QHBoxLayout()
        quick.addWidget(QLabel("p0 step"))
        for lbl, d in (("−90°", -90.0), ("+90°", 90.0), ("180°", 180.0)):
            b = QPushButton(lbl)
            b.setToolTip("add to the zero-order phase (wraps to ±180°)")
            b.clicked.connect(lambda _=False, d=d: self._nudge_p0(d))
            quick.addWidget(b)
        quick.addStretch(1)
        v.addLayout(quick)

        bl = QHBoxLayout()
        bl.addWidget(QLabel("<b>Baseline auto</b> order"))
        self.blOrder = QSpinBox(); self.blOrder.setRange(0, 9); self.blOrder.setValue(3)
        bl.addWidget(self.blOrder)
        self.btnBaseline = QPushButton("Correct")
        bl.addWidget(self.btnBaseline)
        v.addLayout(bl)

        v.addWidget(QLabel("<b>Baseline manual</b> (dmfit-style anchors)"))
        blm = QHBoxLayout()
        self.btnBlPick = QPushButton("Pick anchors")
        self.btnBlPick.setCheckable(True)
        self.btnBlPick.setToolTip("click on the spectrum to place anchor "
                                  "points; drag them to shape the baseline")
        self.btnBlApply = QPushButton("Subtract")
        self.btnBlClear = QPushButton("Clear")
        blm.addWidget(self.btnBlPick)
        blm.addWidget(self.btnBlApply)
        blm.addWidget(self.btnBlClear)
        v.addLayout(blm)

        actions = QHBoxLayout()
        self.chkLive = QCheckBox("live")
        self.chkLive.setChecked(True)
        self.chkLive.setToolTip("re-apply the pipeline to the spectrum as you "
                                "change any control")
        actions.addWidget(self.chkLive)
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
        self.btnBlPick.toggled.connect(self.baseline_mode)
        self.btnBlApply.clicked.connect(self.baseline_apply)
        self.btnBlClear.clicked.connect(self.baseline_clear)

        # ---- live preview: coalesce rapid edits into one re-apply ----
        self._live_timer = QTimer(self)
        self._live_timer.setSingleShot(True)
        self._live_timer.setInterval(120)
        self._live_timer.timeout.connect(lambda: self._emit([]))
        for w in (self.lb, self.gb, self.ssb, self.fcor, self.off, self.sr,
                  self.p0v, self.p1v):
            w.valueChanged.connect(self._schedule_live)
        for w in (self.tdeff, self.zf):
            w.valueChanged.connect(self._schedule_live)
        self.wdw.currentTextChanged.connect(self._schedule_live)
        for w in (self.chkMag, self.chkHilbert, self.rb_raw, self.rb_pdata):
            w.toggled.connect(self._schedule_live)

    def _schedule_live(self, *_):
        if self.chkLive.isChecked():
            self._live_timer.start()

    def _nudge_p0(self, delta: float):
        v_ = self.p0v.value() + delta
        while v_ > 180.0:
            v_ -= 360.0
        while v_ < -180.0:
            v_ += 360.0
        self.p0v.setValue(v_)          # fires valueChanged -> live re-apply

    def _emit(self, extra: list[dict]):
        raw = self.rb_raw.isChecked()
        ops: list[dict] = []
        if raw:
            if self.tdeff.value() > 0:
                ops.append({"op": "tdeff", "points": self.tdeff.value()})
            if self.fcor.value() != 1.0:
                ops.append({"op": "fcor", "factor": self.fcor.value()})
            w = self.wdw.currentText()
            if w == "EM" and self.lb.value():
                ops.append({"op": "em", "lb_hz": abs(self.lb.value())})
            elif w == "GM":
                ops.append({"op": "gm", "lb_hz": -abs(self.lb.value()),
                            "gb": self.gb.value()})
            elif w == "SINE":
                ops.append({"op": "sine", "ssb": self.ssb.value(), "power": 1})
            elif w == "QSINE":
                ops.append({"op": "sine", "ssb": self.ssb.value(), "power": 2})
            elif w == "TRAF":
                ops.append({"op": "traf", "lb_hz": abs(self.lb.value()) or 10.0})
            if self.zf.value() > 1:
                ops.append({"op": "zf", "factor": self.zf.value()})
            ops.append({"op": "ft", "offset_ppm": self.off.value()})
        if not raw and self.chkHilbert.isChecked():
            ops.append({"op": "hilbert"})
        if self.p0v.value() or self.p1v.value():
            ops.append({"op": "phase", "p0": self.p0v.value(), "p1": self.p1v.value()})
        if self.chkMag.isChecked():
            ops.append({"op": "magnitude"})
        if self.sr.value():
            ops.append({"op": "sr", "sr_hz": self.sr.value()})
        ops.extend(extra)
        self.apply_requested.emit(ops, raw)
