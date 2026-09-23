"""Side panels: sites (dmfit-style parameter cards) and processing."""
from __future__ import annotations

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import (
    QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFrame, QGridLayout,
    QHBoxLayout, QLabel, QLineEdit, QPushButton, QRadioButton, QScrollArea,
    QSizePolicy, QSlider, QSpinBox, QToolButton, QVBoxLayout, QWidget,
)

from larmor.desktop import theme
from larmor.desktop.plot import site_color
from larmor.phasedrag import wrap_p0
from larmor.processing import CHANNELS, OPS, TIME_DOMAIN_OPS

PARAM_LABELS = {
    "isotropic_chemical_shift_ppm": "δiso (ppm)",
    "sigma_Cq_MHz": "σ(Cq) (MHz)",
    "Cq_MHz": "Cq (MHz)",
    "Cq_fwhm_MHz": "ΔCq FWHM (MHz)",
    "eta": "η",
    "eta_q": "ηQ",
    "eta_cs": "η CSA",
    "eta_fwhm": "Δη FWHM",
    "eps": "ε (perturbation)",
    "zeta_ppm": "ζ CSA (ppm)",
    "sigma_zeta_ppm": "σ(ζ) (ppm)",
    "shift_fwhm_ppm": "FWHM (ppm)",
    "line_fwhm_ppm": "line broadening (ppm)",
    "gauss_fwhm_ppm": "Gauss FWHM (ppm)",
    "lorentz_fwhm_ppm": "Lorentz FWHM (ppm)",
    "shift_ppm": "shift (ppm)",
    "j_hz": "J (Hz)",
    "n_j": "n couplings",
    "ssb_ratio": "sideband ratio",
    "n_ssb": "n sidebands",
    "amplitude": "amplitude",
    "gl": "g/l fraction",
    "czjzek_d": "d (Czjzek dimension)",
    "shift_slope_ppm_per_MHz": "dδiso/dC_Q (ppm/MHz)",
    "split_ppm": "Δδ A−B (ppm)",
    "pop_a": "p(A)",
    "k_ex_hz": "k_ex (s⁻¹)",
    "a": "a (function)",
    "b": "b (function)",
    "c": "c (function)",
    "d": "d (function)",
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
        tag.setStyleSheet(f"color: {theme.active().text_dim}; font-size: 10px;")
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
            self.setStyleSheet(f"#siteCard {{ background: {theme.active().alt_base}; }} * {{ color: {theme.active().disabled_text}; }}")

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
                lab.setStyleSheet(f"color: {theme.active().accent}; font-weight: 600;")
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
            err.setStyleSheet(f"color: {theme.active().accent}; font-size: 10px;")
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
        self._hint.setStyleSheet(f"color: {theme.active().text_dim};")
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
            self._hint.setStyleSheet(f"color: {theme.active().text_dim};")
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
    twopoint_mode = Signal(bool)           # 2-point background: pick-two toggle
    twopoint_apply = Signal()
    twopoint_clear = Signal()
    phase_drag_mode = Signal(bool)         # TopSpin drag-to-phase toggle
    #: display projection requested: (domain "time"|"freq", channel
    #: "real"|"imag"|"magnitude") -- a view of the pipeline result, never a
    #: pipeline step (see view_state / reset_view)
    view_changed = Signal(str, str)

    #: recorded time-domain steps the panel has no control for
    _UNSYNCABLE = frozenset({"lp", "shift_fid", "swap_echo", "echo_apodize"})
    _WINDOW_OPS = ("em", "gm", "sine", "traf")

    def __init__(self):
        super().__init__()
        #: frequency-domain steps of a recorded chain the widgets cannot express
        #: (baseline, autophase, subtract_avg, ...), re-appended by _emit so a
        #: forced re-apply keeps them (see sync_from_ops)
        self._carried_ops: list[dict] = []
        self._hilbert_before_reapod = False
        # All controls live inside a scroll area so this panel can be made
        # narrow without forcing the main window wider than the screen (a wide
        # row scrolls instead of pushing the whole window past the monitor).
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        outer.addWidget(scroll)
        content = QWidget()
        scroll.setWidget(content)
        self.setMinimumWidth(280)                # usable floor, well under a screen
        v = QVBoxLayout(content)
        v.setAlignment(Qt.AlignTop)

        self.rb_pdata = QRadioButton("TopSpin-processed (pdata)")
        self.rb_raw = QRadioButton("raw fid (EM → ZF → FT)")
        self.rb_pdata.setChecked(True)
        v.addWidget(self.rb_pdata)
        v.addWidget(self.rb_raw)

        # raw-FID window functions are advanced: collapsed by default so pdata
        # users are not faced with them; auto-expands when 'raw fid' is picked
        self.adv_toggle = QToolButton()
        self.adv_toggle.setText("▸ Raw-FID window functions (advanced)")
        self.adv_toggle.setCheckable(True)
        self.adv_toggle.setStyleSheet(
            "QToolButton { border: none; font-weight: 600; }")
        v.addWidget(self.adv_toggle)
        self._adv = QWidget()
        self._adv.setVisible(False)
        adv = QVBoxLayout(self._adv)
        adv.setContentsMargins(8, 0, 0, 0)
        v.addWidget(self._adv)

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
        adv.addLayout(wdw)

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
        adv.addLayout(raw)

        # re-apodize a spectrum that did NOT come from a raw fid (TopSpin 1r,
        # CSV, the FID / QCPMG / VOCS dialogs' output): Hilbert -> IFT -> the
        # window block above -> FT, still applied from the unprocessed base
        self.chkReapod = QCheckBox(
            "re-apodize this spectrum  (Hilbert → IFT → window → FT)")
        self.chkReapod.setToolTip(
            "for a TopSpin-processed 1r, a CSV, or a spectrum sent from the "
            "FID / QCPMG / VOCS dialogs: rebuild the imaginary channel, go back "
            "to the FID, apply the window functions above, transform again — "
            "a reconstruction that compounds with whatever window was already "
            "applied; the raw-fid mode restarts from the instrument file and "
            "is exact. Hilbert first is mandatory (the IFT of a real-only "
            "spectrum is two-sided) and stays locked while this is on.")
        adv.addWidget(self.chkReapod)

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

        # Display: which projection of the pipeline result the canvas shows.
        # The FID button flips to the windowed, zero-filled FID the transform
        # sees (window / LB / GB / ZF re-apply live there); the radios pick the
        # real, imaginary or magnitude channel. Display only: the fit, S/N,
        # save and overlays always use the real frequency-domain spectrum.
        v.addWidget(QLabel("<b>Display</b>"))
        disp = QHBoxLayout()
        self.btnDomain = QPushButton("FID ⇄ spectrum")
        self.btnDomain.setCheckable(True)
        self.btnDomain.setToolTip(
            "show the windowed FID the transform sees — WDW / LB / GB / ZF "
            "re-apply live; press again for the spectrum (Ctrl+T)")
        disp.addWidget(self.btnDomain)
        # REQUIRED: rb_pdata / rb_raw are plain radios on the same content
        # widget and Qt auto-groups sibling radios -- without an explicit
        # group, checking 'imag' would un-check the source radio and flip the
        # pipeline to raw. Buttons in a QButtonGroup leave the sibling group.
        self._chan_group = QButtonGroup(self)
        self.rb_real = QRadioButton("real")
        self.rb_imag = QRadioButton("imag")
        self.rb_mag = QRadioButton("|S|")
        self.rb_real.setToolTip("the real channel — what the fit uses")
        self.rb_imag.setToolTip("the imaginary channel: inspect it while "
                                "phasing (dispersion should vanish under the "
                                "pivot when p0 / p1 are right); Ctrl+I cycles")
        self.rb_mag.setToolTip("|S| for display only — the destructive "
                               "'magnitude' pipeline op is the checkbox above")
        for i, rb in enumerate((self.rb_real, self.rb_imag, self.rb_mag)):
            self._chan_group.addButton(rb, i)
            disp.addWidget(rb)
        self.rb_real.setChecked(True)
        disp.addStretch(1)
        v.addLayout(disp)

        v.addWidget(QLabel("<b>Phase</b>"))
        self.btnAuto = QPushButton("Autophase (ACME)")
        # TopSpin drag-to-phase: a checkable mode the workbench arms on the
        # plot (see MainWindow._phase_drag_mode); this button is the single
        # source of truth for the mode, the Process-menu entry toggles it
        self.btnDrag = QPushButton("Drag to phase")
        self.btnDrag.setCheckable(True)
        self.btnDrag.setToolTip(
            "TopSpin-style: drag on the spectrum — left/right = p0 (0.25°/px), "
            "up/down = p1 (1°/px) about the pivot line; Shift = fine ×0.1, "
            "Ctrl+drag = pan; start the drag on empty canvas (items keep their "
            "own drags); Esc or click again to stop "
            "(Process ▸ Drag to phase, Ctrl+P)")
        phb = QHBoxLayout()
        phb.addWidget(self.btnAuto)
        phb.addWidget(self.btnDrag)
        v.addLayout(phb)
        ph0 = QHBoxLayout()
        ph0.addWidget(QLabel("p0"))
        self.p0 = QSlider(Qt.Horizontal); self.p0.setRange(-180, 180)
        self.p0v = QDoubleSpinBox(); self.p0v.setRange(-180, 180)
        self.p0.valueChanged.connect(
            lambda v_: self._slider_to_spin(self.p0v, v_))
        self.p0v.valueChanged.connect(lambda v_: self.p0.setValue(int(v_)))
        ph0.addWidget(self.p0); ph0.addWidget(self.p0v)
        v.addLayout(ph0)
        ph1 = QHBoxLayout()
        ph1.addWidget(QLabel("p1"))
        self.p1 = QSlider(Qt.Horizontal); self.p1.setRange(-720, 720)
        self.p1v = QDoubleSpinBox(); self.p1v.setRange(-720, 720)
        self.p1.valueChanged.connect(
            lambda v_: self._slider_to_spin(self.p1v, v_))
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

        # 2-point background: pick two baseline points, subtract the straight line
        # through them (removes a flat/tilted background before the auto baseline)
        tp = QHBoxLayout()
        tp.addWidget(QLabel("<b>2-point background</b>"))
        self.btnTpPick = QPushButton("Pick 2 points"); self.btnTpPick.setCheckable(True)
        self.btnTpPick.setToolTip("click two baseline points on the spectrum "
                                  "(one each side of the peaks); the straight line "
                                  "through them is subtracted — turn off to apply")
        self.btnTpApply = QPushButton("Subtract")
        self.btnTpClear = QPushButton("Clear")
        tp.addWidget(self.btnTpPick); tp.addWidget(self.btnTpApply)
        tp.addWidget(self.btnTpClear)
        v.addLayout(tp)

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
        note.setStyleSheet(f"color: {theme.active().text_dim};")
        v.addWidget(note)

        self.btnApply.clicked.connect(lambda: self._emit([]))
        self.btnAuto.clicked.connect(lambda: self._emit([{"op": "autophase"}]))
        self.btnBaseline.clicked.connect(
            lambda: self._emit([{"op": "baseline", "order": self.blOrder.value()}]))
        self.btnReset.clicked.connect(self._clear_carried)   # before the reload
        self.btnReset.clicked.connect(self.reset_requested)
        self.btnBlPick.toggled.connect(self.baseline_mode)
        self.btnBlApply.clicked.connect(self.baseline_apply)
        self.btnBlClear.clicked.connect(self.baseline_clear)
        self.btnTpPick.toggled.connect(self.twopoint_mode)
        self.btnTpApply.clicked.connect(self.twopoint_apply)
        self.btnTpClear.clicked.connect(self.twopoint_clear)
        self.btnDrag.toggled.connect(self.phase_drag_mode)

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
        for w in (self.chkMag, self.chkHilbert, self.rb_raw, self.rb_pdata,
                  self.chkReapod):
            w.toggled.connect(self._schedule_live)
        self.chkReapod.toggled.connect(self._toggle_reapod)
        # re-apodizing is meaningless in raw-fid mode (the window applies to
        # the instrument fid directly)
        self.rb_raw.toggled.connect(lambda on: self.chkReapod.setEnabled(not on))
        self.adv_toggle.toggled.connect(self._toggle_adv)
        # picking raw-FID mode reveals the window-function controls it needs
        self.rb_raw.toggled.connect(
            lambda on: self.adv_toggle.setChecked(True) if on else None)
        # display projection: both controls funnel into one signal
        self.btnDomain.toggled.connect(self._emit_view)
        self._chan_group.idToggled.connect(
            lambda _id, on: self._emit_view() if on else None)

    def _toggle_adv(self, on: bool):
        self._adv.setVisible(on)
        self.adv_toggle.setText(("▾ " if on else "▸ ")
                                + "Raw-FID window functions (advanced)")

    def _schedule_live(self, *_):
        if self.chkLive.isChecked():
            self._live_timer.start()

    @staticmethod
    def _slider_to_spin(spin: QDoubleSpinBox, v: int):
        """An integer slider tick reaches its spin box only when it carries
        new information. The spin box writes int(value) to the slider, whose
        valueChanged echoed the truncated value straight back, so any
        fractional phase written to the spin box (typed, dragged,
        re-expressed by a pivot move) was clobbered whenever the integer
        part changed (35.96 became 35.0)."""
        if int(spin.value()) != int(v):
            spin.setValue(float(v))

    def _nudge_p0(self, delta: float):
        # fires valueChanged -> live re-apply; wraps like the drag gesture
        self.p0v.setValue(wrap_p0(self.p0v.value() + delta))

    # ------------------------------------------------------- drag to phase
    def phase_values(self) -> tuple[float, float]:
        """(p0, p1) as the controls show them (degrees)."""
        return (self.p0v.value(), self.p1v.value())

    def set_phase(self, p0: float, p1: float, apply: bool = True):
        """Write both phase controls (p0 wrapped to +-180, p1 clamped by the
        spin range; the sliders follow through their valueChanged links) and,
        with ``apply``, re-apply the chain at once through the same op-list
        builder the sliders use. The 120 ms debounce the setValue calls just
        armed is stopped, so a drag applies exactly once per event and an
        Undo resync (``apply=False``) applies nothing."""
        self.p0v.setValue(wrap_p0(p0))
        self.p1v.setValue(float(p1))
        self._live_timer.stop()
        if apply:
            self._emit([])

    def sync_phase_from(self, ops):
        """Mirror the LAST phase step of a recorded chain (both 0.0 when it
        has none) into the controls, silently -- the Undo/Redo resync."""
        p0 = p1 = 0.0
        for step in ops or []:
            if step.get("op") == "phase":
                p0 = float(step.get("p0", 0.0))
                p1 = float(step.get("p1", 0.0))
        self.set_phase(p0, p1, apply=False)

    # ---------------------------------------------------------- display state
    def view_state(self) -> tuple[str, str]:
        """(domain, channel) the canvas should show: ("freq"|"time",
        "real"|"imag"|"magnitude")."""
        cid = self._chan_group.checkedId()
        channel = CHANNELS[cid] if 0 <= cid < len(CHANNELS) else "real"
        return ("time" if self.btnDomain.isChecked() else "freq", channel)

    def reset_view(self):
        """Back to (freq, real) SILENTLY -- no view_changed. The workbench
        calls this when new data arrives or a bypassing edit invalidates the
        complex pipeline result."""
        ws = (self.btnDomain, self.rb_real, self.rb_imag, self.rb_mag,
              self._chan_group)
        for w in ws:
            w.blockSignals(True)
        try:
            self.btnDomain.setChecked(False)
            self.rb_real.setChecked(True)
        finally:
            for w in ws:
                w.blockSignals(False)

    def _emit_view(self, *_):
        self.view_changed.emit(*self.view_state())

    def chain_has_ft(self) -> bool:
        """Does the chain _emit builds contain a transform (so a windowed FID
        exists to show)? Raw-fid mode always does; pdata only re-apodizing."""
        return self.rb_raw.isChecked() or self.chkReapod.isChecked()

    def _toggle_reapod(self, on: bool):
        # Hilbert is mandatory for a one-sided FID: force and lock it while
        # re-apodizing, restore the user's own choice afterwards
        if on:
            self._hilbert_before_reapod = self.chkHilbert.isChecked()
            self.chkHilbert.setChecked(True)
            self.chkHilbert.setEnabled(False)
        else:
            self.chkHilbert.setEnabled(True)
            self.chkHilbert.setChecked(self._hilbert_before_reapod)

    def _clear_carried(self, *_):
        self._carried_ops = []

    def arm_reapodize(self):
        """Tick 're-apodize this spectrum' for the caller's own synchronous
        Apply: Hilbert forced + locked, the window block expanded, EM selected
        when no window was (so LB is live at once), and the live timer stopped
        so the caller's Apply is the only one."""
        self.chkReapod.blockSignals(True)
        try:
            self.chkReapod.setChecked(True)
        finally:
            self.chkReapod.blockSignals(False)
        self._toggle_reapod(True)              # signals were blocked: run it
        if self.wdw.currentText() == "none":
            self.wdw.blockSignals(True)
            self.wdw.setCurrentText("EM")
            self.wdw.blockSignals(False)
        self.adv_toggle.setChecked(True)
        self._live_timer.stop()

    def arm_hilbert(self):
        """Tick 'Hilbert first' without a live tick (the caller applies): a
        real-only spectrum has an identically zero imaginary channel."""
        self.chkHilbert.blockSignals(True)
        try:
            self.chkHilbert.setChecked(True)
        finally:
            self.chkHilbert.blockSignals(False)
        self._live_timer.stop()

    def sync_from_ops(self, ops: list[dict], use_raw: bool) -> bool:
        """Set every control from a RECORDED chain so that the next _emit([])
        reproduces it -- the inverse of _emit. Frequency-domain steps the
        widgets cannot express (baseline, autophase, subtract_avg, ...) are
        carried and re-appended in recorded order. Returns False WITHOUT
        touching any widget when the chain holds a time-domain step the panel
        has no control for (lp, shift_fid, whole-echo, zf to an absolute SI,
        two windows, ...). Silent: no live tick, no signal."""
        ops = [dict(o) for o in (ops or [])]
        names = [o.get("op") for o in ops]
        if any(n not in OPS for n in names) or (self._UNSYNCABLE & set(names)):
            return False
        if names.count("ft") > 1 or names.count("ift") > 1:
            return False
        if "ift" in names and use_raw:
            return False
        if sum(n in self._WINDOW_OPS for n in names) > 1:
            return False
        for o in ops:
            if o.get("op") == "zf" and (
                    o.get("si") or not 1 <= int(o.get("factor", 2)) <= 16):
                return False
            if o.get("op") == "sine" and int(o.get("power", 1)) not in (1, 2):
                return False
        # time-domain steps must sit between the ift (if any) and the ft
        td = [i for i, n in enumerate(names)
              if n in TIME_DOMAIN_OPS and n != "ft"]
        if td or "ft" in names:
            if "ft" not in names or not (use_raw or "ift" in names):
                return False
            k_ft = names.index("ft")
            k_ift = names.index("ift") if "ift" in names else -1
            if any(not (k_ift < i < k_ft) for i in td):
                return False

        widgets = (self.wdw, self.lb, self.gb, self.ssb, self.tdeff, self.zf,
                   self.fcor, self.off, self.sr, self.p0, self.p0v, self.p1,
                   self.p1v, self.chkMag, self.chkHilbert, self.chkReapod,
                   self.rb_raw, self.rb_pdata)
        for w in widgets:
            w.blockSignals(True)
        try:
            self._carried_ops = []
            # neutral: EM with LB 0 emits no window step (an empty chain
            # round-trips to an empty chain) and leaves LB live
            self.wdw.setCurrentText("EM"); self.lb.setValue(0.0)
            self.gb.setValue(0.1); self.ssb.setValue(2.0)
            self.tdeff.setValue(0); self.zf.setValue(1); self.fcor.setValue(1.0)
            self.off.setValue(0.0); self.sr.setValue(0.0)
            self.p0.setValue(0); self.p0v.setValue(0.0)
            self.p1.setValue(0); self.p1v.setValue(0.0)
            self.chkMag.setChecked(False)
            self.chkReapod.setChecked(False)
            self.chkHilbert.setEnabled(True); self.chkHilbert.setChecked(False)
            self.rb_raw.setChecked(bool(use_raw))
            self.rb_pdata.setChecked(not use_raw)
            self.chkReapod.setEnabled(not use_raw)
            for o in ops:
                name = o.pop("op")
                if name == "tdeff":
                    self.tdeff.setValue(int(o.get("points", 0)))
                elif name == "fcor":
                    self.fcor.setValue(float(o.get("factor", 0.5)))
                elif name == "em":
                    self.wdw.setCurrentText("EM")
                    self.lb.setValue(abs(float(o.get("lb_hz", 0.0))))
                elif name == "gm":
                    self.wdw.setCurrentText("GM")
                    self.lb.setValue(abs(float(o.get("lb_hz", -10.0))))
                    self.gb.setValue(float(o.get("gb", 0.1)))
                elif name == "sine":
                    self.wdw.setCurrentText(
                        "QSINE" if int(o.get("power", 1)) == 2 else "SINE")
                    self.ssb.setValue(float(o.get("ssb", 2.0)))
                elif name == "traf":
                    self.wdw.setCurrentText("TRAF")
                    self.lb.setValue(abs(float(o.get("lb_hz", 10.0))))
                elif name == "zf":
                    self.zf.setValue(int(o.get("factor", 2)))
                elif name == "ft":
                    self.off.setValue(float(o.get("offset_ppm", 0.0)))
                elif name == "hilbert":
                    self.chkHilbert.setChecked(True)
                elif name == "ift":
                    self.chkReapod.setChecked(True)
                    self._toggle_reapod(True)
                elif name == "phase":            # pivot_frac: the live pivot governs
                    p0 = float(o.get("p0", 0.0)); p1 = float(o.get("p1", 0.0))
                    self.p0v.setValue(p0); self.p0.setValue(int(round(p0)))
                    self.p1v.setValue(p1); self.p1.setValue(int(round(p1)))
                elif name == "magnitude":
                    self.chkMag.setChecked(True)
                elif name == "sr":
                    self.sr.setValue(float(o.get("sr_hz", 0.0)))
                else:
                    self._carried_ops.append({"op": name, **o})
            if use_raw or self.chkReapod.isChecked():
                self.adv_toggle.setChecked(True)
        finally:
            for w in widgets:
                w.blockSignals(False)
        self._live_timer.stop()
        return True

    # ----------------------------------------------------------- the chain
    def _time_domain_ops(self, with_fcor: bool) -> list[dict]:
        """TDeff, (FCOR,) window, ZF -- the panel's time-domain block, in the
        order the raw-fid chain has always used."""
        ops: list[dict] = []
        if self.tdeff.value() > 0:
            ops.append({"op": "tdeff", "points": self.tdeff.value()})
        if with_fcor and self.fcor.value() != 1.0:
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
        return ops

    def _emit(self, extra: list[dict]):
        raw = self.rb_raw.isChecked()
        reapod = not raw and self.chkReapod.isChecked()
        ops: list[dict] = []
        if raw:
            ops += self._time_domain_ops(with_fcor=True)
            ops.append({"op": "ft", "offset_ppm": self.off.value()})
        elif reapod:
            # Hilbert first is mandatory (the ift of a real-only spectrum is
            # two-sided; a one-sided window on it gives a dispersive line); no
            # fcor (halving the first point of an ift'd spectrum injects a DC
            # step instead of removing one) and no offset (op_ft restores the
            # held axis, centre included)
            ops += [{"op": "hilbert"}, {"op": "ift"}]
            ops += self._time_domain_ops(with_fcor=False)
            ops.append({"op": "ft"})
        if not raw and not reapod and self.chkHilbert.isChecked():
            ops.append({"op": "hilbert"})
        if self.p0v.value() or self.p1v.value():
            ops.append({"op": "phase", "p0": self.p0v.value(), "p1": self.p1v.value()})
        if self.chkMag.isChecked():
            ops.append({"op": "magnitude"})
        if self.sr.value():
            ops.append({"op": "sr", "sr_hz": self.sr.value()})
        ops.extend(self._carried_ops)
        ops.extend(extra)
        self.apply_requested.emit(ops, raw)
