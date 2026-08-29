"""QCPMG processing — a guided, stage-by-stage workflow.

Turns a raw echo train into a FITTABLE spectrum, following the ssNake sum-echo
protocol but with the bookkeeping done for you. Six numbered stages, each
showing the plot you expect and the one number you read off it:

  1 Train & split   — the echo period, READ from the pulse program (CNST7)
  2 Echo & top      — all echoes / first-vs-last (the split check); drag the
                      top marker (ssNake's "Pos 147")
  3 Decay & T2      — echo-top intensity vs echo number, fitted; click a point
                      to exclude it
  4 Apodization     — one click applies the matched filter LB = 1/(pi*T2)
  5 Spectrum        — sum echo AND spikelets on one axes (both, not either/or)
  6 Measure         — delta_CG and FWHM in a draggable window, with the
                      window-sensitivity as an honest error bar

The headline readout and "Send to fit →" stay visible from every stage.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from larmor.desktop import theme
from larmor.desktop.plot_menu import attach_plot_menu
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFileDialog, QFrame, QHBoxLayout, QLabel, QPushButton, QScrollArea,
    QSpinBox, QSplitter, QTabWidget, QVBoxLayout, QWidget,
)


_log = logging.getLogger(__name__)


def _plot(title: str = "", ppm_axis: bool = False, height: int | None = None,
          parent=None):
    pw = pg.PlotWidget(background=theme.active().plot_bg)
    # keep each plot's minimum small: tabs stack up to two of them, and their
    # combined minimum is what decides how far the DIALOG can be shrunk
    pw.setMinimumHeight(90)
    if ppm_axis:
        pw.getPlotItem().invertX(True)
        pw.setLabel("bottom", "shift", units="ppm")
        pw.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
    if height:
        pw.setMaximumHeight(height)
    if title:
        pw.setTitle(title, color=theme.active().text_dim, size="9pt")
    try:                       # right-click: Export figure / Send to studio
        attach_plot_menu(pw, title=title or "QCPMG", parent=parent)
    except Exception:          # noqa: BLE001 - a menu is never worth a crash
        pass
    return pw


def _scrolled(w: QWidget) -> QScrollArea:
    """Wrap a stage page so its layout minimum can never block shrinking the
    dialog: below the content's minimum the page scrolls instead."""
    sa = QScrollArea()
    sa.setWidget(w)
    sa.setWidgetResizable(True)
    sa.setFrameShape(QFrame.NoFrame)
    return sa


def _row(*widgets) -> QWidget:
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    for x in widgets:
        h.addWidget(QLabel(x) if isinstance(x, str) else x)
    h.addStretch(1)
    return w


class QcpmgDialog(QDialog):
    #: (ppm, amp, meta) of the processed spectrum to fit
    accepted_1d = Signal(object, object, dict)

    def __init__(self, parent, source: str | None):
        super().__init__(parent)
        self.setWindowTitle("QCPMG processing — guided sum-echo workflow")
        # freely resizable: each stage sits in a scroll area, so no layout
        # minimum can ever block shrinking; the size the user leaves the
        # dialog at is remembered and restored on the next open
        self.setMinimumSize(480, 360)
        self.setSizeGripEnabled(True)
        restored = False
        if not os.environ.get("LARMOR_NO_SESSION"):
            try:
                from PySide6.QtCore import QSettings
                geo = QSettings("LARMOR", "app").value("qcpmgDialogGeometry")
                if geo is not None:
                    restored = bool(self.restoreGeometry(geo))
            except Exception:
                restored = False
        if not restored:
            w, h = 1080, 760
            scr = self.screen()
            if scr is not None:
                avail = scr.availableGeometry()
                w = min(w, avail.width() - 80)
                h = min(h, avail.height() - 100)
            self.resize(max(640, w), max(480, h))
        self.source = source
        self.fid = None
        self.meta: dict = {}
        self.t2 = None                  # qcpmg.T2Fit
        # exclusions are stored as ABSOLUTE train-echo indices, so changing
        # drop_first / nEch cannot silently re-point them at a different echo
        self.excluded: set[int] = set()
        self.keep: np.ndarray | None = None      # derived mask, for plotting
        self._phased = False            # has the spectrum been phased at all?
        self._window_user = False       # did the user place the CG window?
        self._cg = self._sigma = self._fwhm = None
        self._ppm = self._spec_raw = self._spec = None
        self._spk_ppm = self._spk = None
        self._period_src = "none"
        self._carrier = 0.0
        self._referenced = False
        self._loading = False

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(150)
        self._timer.timeout.connect(self._recompute)

        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.lbl = QLabel(source or "no echo-train FID selected")
        self.lbl.setStyleSheet("font-weight: 600;")
        self.lbl.setWordWrap(True)
        b = QPushButton("Open echo-train FID…")
        b.setToolTip("choose the raw 'fid' file (not pdata) of a QCPMG EXPNO")
        b.clicked.connect(self._pick)
        top.addWidget(self.lbl, 1)
        top.addWidget(b)
        v.addLayout(top)

        self.tabs = QTabWidget()
        self.tabs.tabBar().setUsesScrollButtons(True)   # 6 tabs, narrow window
        self.tabs.addTab(_scrolled(self._build_stage1()), "1 · Train && split")
        self.tabs.addTab(_scrolled(self._build_stage2()), "2 · Echo && top")
        self.tabs.addTab(_scrolled(self._build_stage3()), "3 · Decay && T₂")
        self.tabs.addTab(_scrolled(self._build_stage4()), "4 · Apodization")
        self.tabs.addTab(_scrolled(self._build_stage5()), "5 · Spectrum")
        self.tabs.addTab(_scrolled(self._build_stage6()), "6 · Measure")
        v.addWidget(self.tabs, 1)

        # typing '293' digit by digit must not act on the intermediate '29':
        # _on_period_pts's top.setMaximum(28) would clamp a top of 147 down and
        # the debounced recompute would silently run from the wrong point
        for sb in (self.findChildren(QSpinBox)
                   + self.findChildren(QDoubleSpinBox)):
            sb.setKeyboardTracking(False)

        # persistent headline + actions: visible from every stage
        self.res = QLabel("")
        self.res.setWordWrap(True)
        self.res.setTextFormat(Qt.RichText)
        self.res.setStyleSheet(f"color: {theme.active().accent}; font-weight: 600;")
        v.addWidget(self.res)

        bb = QDialogButtonBox(QDialogButtonBox.Help)
        self.btnCsv = bb.addButton("Copy CSV", QDialogButtonBox.ActionRole)
        self.btnCsv.setToolTip("decay points + every processing value, for the lab book")
        self.btnCsv.clicked.connect(self._copy_csv)
        self.btnSend = bb.addButton("Send to fit →", QDialogButtonBox.AcceptRole)
        self.btnSend.setDefault(True)
        self.btnSend.clicked.connect(self._send)
        bb.helpRequested.connect(self._help)
        v.addWidget(bb)

        for w in (self.btnCsv, self.btnSend):
            w.setEnabled(False)
        if source:
            self._load(source)

    # ---------------------------------------------------------------- stages
    def _build_stage1(self) -> QWidget:
        w = QWidget(); lv = QVBoxLayout(w)
        self.period = QSpinBox(); self.period.setRange(4, 10_000_000)
        self.period.setToolTip("points per echo")
        self.periodHz = QDoubleSpinBox(); self.periodHz.setRange(0.01, 1e7)
        self.periodHz.setDecimals(2); self.periodHz.setSuffix(" Hz")
        self.periodHz.setToolTip("spikelet spacing = 1/τecho — the number the "
                                 "pulse program sets (CNST7)")
        self.nEch = QSpinBox(); self.nEch.setRange(1, 100000)
        self.nEch.setToolTip("echoes to use — later ones are mostly noise")
        self.dropFirst = QSpinBox(); self.dropFirst.setRange(0, 100)
        self.dropFirst.setToolTip("discard leading echoes (a partial first echo)")
        for s in (self.period, self.nEch, self.dropFirst):
            s.valueChanged.connect(self._on_period_pts if s is self.period
                                   else self._queue)
        self.periodHz.valueChanged.connect(self._on_period_hz)
        lv.addWidget(_row("period", self.period, "pts   =", self.periodHz,
                          "    echoes", self.nEch, "   drop first", self.dropFirst))
        self.p_train = _plot("echo train (|FID|) — dotted lines mark the period", parent=self)
        self.p_train.setLabel("bottom", "point")
        lv.addWidget(self.p_train, 1)
        self.lbl1 = QLabel(""); self.lbl1.setWordWrap(True)
        lv.addWidget(self.lbl1)
        return w

    def _build_stage2(self) -> QWidget:
        w = QWidget(); lv = QVBoxLayout(w)
        self.top = QSpinBox(); self.top.setRange(0, 10_000_000)
        self.top.setToolTip("point index of the echo top INSIDE each echo "
                            "(ssNake's 'Pos N') — drag the line instead")
        self.top.valueChanged.connect(self._on_top_spin)
        bAuto = QPushButton("Auto (mean of all echoes)")
        bAuto.clicked.connect(self._auto_top)
        lv.addWidget(_row("echo top", self.top, bAuto))
        split = QSplitter(Qt.Vertical)
        self.p_echoes = _plot("all echoes — drag the vertical line onto the top", parent=self)
        self.p_echoes.setLabel("bottom", "point in echo")
        self.p_firstlast = _plot("first vs last echo — features must line up, "
                                 "or the period is wrong", parent=self)
        self.p_firstlast.setLabel("bottom", "point in echo")
        split.addWidget(self.p_echoes); split.addWidget(self.p_firstlast)
        split.setSizes([420, 260])
        lv.addWidget(split, 1)
        self.topLine = pg.InfiniteLine(
            angle=90, movable=True,
            pen=pg.mkPen(theme.active().pivot, width=1.6),
            hoverPen=pg.mkPen(theme.active().accent, width=2.4),
            label="top {value:0.0f}",
            labelOpts={"position": 0.92, "color": theme.active().pivot})
        self.topLine.sigPositionChangeFinished.connect(self._on_top_line)
        self.p_echoes.addItem(self.topLine)
        self.centreLine = pg.InfiniteLine(
            angle=90, movable=False,
            pen=pg.mkPen(theme.active().text_dim, width=1, style=Qt.DotLine))
        self.p_echoes.addItem(self.centreLine)
        self.lbl2 = QLabel(""); self.lbl2.setWordWrap(True)
        lv.addWidget(self.lbl2)
        return w

    def _build_stage3(self) -> QWidget:
        w = QWidget(); lv = QVBoxLayout(w)
        self.decayMode = QComboBox()
        self.decayMode.addItems(["signed real (ssNake)", "magnitude"])
        self.decayMode.setToolTip("ssNake samples the SIGNED real part; "
                                  "magnitude has a rectified noise floor that "
                                  "biases the tail upward")
        self.decayMode.currentIndexChanged.connect(self._queue)
        self.useOffset = QCheckBox("fit a constant offset")
        self.useOffset.setChecked(True)
        self.useOffset.setToolTip("ssNake's model is C + B·exp(-t/T2); "
                                  "dropping C shifted T2 by 6–67% "
                                  "(median 26%) on a real 12-sample set")
        self.useOffset.toggled.connect(self._queue)
        lv.addWidget(_row("decay from", self.decayMode, self.useOffset))
        self.p_decay = _plot("echo-top intensity vs echo number — "
                             "click a point to exclude it, then it refits", parent=self)
        self.p_decay.setLabel("bottom", "time", units="s")
        lv.addWidget(self.p_decay, 1)
        self.lbl3 = QLabel(""); self.lbl3.setWordWrap(True)
        self.lbl3.setTextFormat(Qt.RichText)
        lv.addWidget(self.lbl3)
        return w

    def _build_stage4(self) -> QWidget:
        w = QWidget(); lv = QVBoxLayout(w)
        self.t2w = QCheckBox("T₂-weight the echo sum (matched filter)")
        self.t2w.setChecked(True)
        self.t2w.toggled.connect(self._queue)
        self.btnMatched = QPushButton("Use matched LB")
        self.btnMatched.setToolTip("set the Lorentzian broadening to 1/(π·T₂)")
        self.btnMatched.clicked.connect(self._use_matched)
        self.lb = QDoubleSpinBox(); self.lb.setRange(0, 1e6)
        self.lb.setDecimals(1); self.lb.setSuffix(" Hz")
        self.lb.setToolTip("Lorentzian broadening of the summed echo")
        self.lb.valueChanged.connect(self._queue)
        self.gb = QDoubleSpinBox(); self.gb.setRange(0, 1e6)
        self.gb.setDecimals(1); self.gb.setSuffix(" Hz"); self.gb.setValue(0.0)
        self.gb.setToolTip("Gaussian broadening of the summed echo")
        self.gb.valueChanged.connect(self._queue)
        lv.addWidget(_row(self.t2w, self.btnMatched, "  LB", self.lb,
                          "  GB", self.gb))
        split = QSplitter(Qt.Vertical)
        self.p_apod_echoes = _plot("apodized echoes (ssNake 'Apodised echoes')", parent=self)
        self.p_apod_echoes.setLabel("bottom", "point in echo")
        self.p_apod_decay = _plot("the weighting along the echo dimension "
                                  "(ssNake 'Apodised D1')", parent=self)
        self.p_apod_decay.setLabel("bottom", "time", units="s")
        split.addWidget(self.p_apod_echoes); split.addWidget(self.p_apod_decay)
        split.setSizes([360, 300])
        lv.addWidget(split, 1)
        self.lbl4 = QLabel(""); self.lbl4.setWordWrap(True)
        lv.addWidget(self.lbl4)
        return w

    def _build_stage5(self) -> QWidget:
        w = QWidget(); lv = QVBoxLayout(w)
        self.showSum = QCheckBox("sum echo (absorption — fit this)")
        self.showSum.setChecked(True); self.showSum.toggled.connect(self._redraw_spec)
        self.showSpk = QCheckBox("spikelets")
        self.showSpk.toggled.connect(self._redraw_spec)
        self.norm = QComboBox(); self.norm.addItems(["unit max", "area", "raw"])
        self.norm.currentIndexChanged.connect(self._redraw_spec)
        self.zf = QComboBox(); self.zf.addItems(["2", "4", "8", "16", "32"])
        self.zf.setCurrentText("16"); self.zf.currentIndexChanged.connect(self._queue)
        lv.addWidget(_row(self.showSum, self.showSpk, "  scale", self.norm,
                          "  zero-fill ×", self.zf))
        self.p0 = QDoubleSpinBox(); self.p0.setRange(-720, 720)
        self.p0.setDecimals(2); self.p0.setWrapping(True)
        self.p0.valueChanged.connect(self._rephase)
        self.p1 = QDoubleSpinBox(); self.p1.setRange(-36000, 36000)
        self.p1.setDecimals(2); self.p1.valueChanged.connect(self._rephase)
        self.p1.setToolTip("rarely needed — a correct whole-echo transform is "
                           "already near-pure absorption, so p0 alone usually does it")
        self.pstep = QDoubleSpinBox(); self.pstep.setRange(0.01, 90.0)
        self.pstep.setDecimals(2); self.pstep.setValue(1.0)
        self.pstep.setToolTip("phase increment per click/scroll")
        self.pstep.valueChanged.connect(self._set_pstep)
        bAuto = QPushButton("Autophase"); bAuto.clicked.connect(self._autophase)
        lv.addWidget(_row("p0", self.p0, "p1", self.p1, "step", self.pstep, bAuto))
        self.p_spec = _plot("QCPMG spectrum", ppm_axis=True, parent=self)
        lv.addWidget(self.p_spec, 1)
        self.lbl5 = QLabel(""); self.lbl5.setWordWrap(True)
        self.lbl5.setTextFormat(Qt.RichText)
        lv.addWidget(self.lbl5)
        self._set_pstep()
        return w

    def _build_stage6(self) -> QWidget:
        w = QWidget(); lv = QVBoxLayout(w)
        bAuto = QPushButton("Auto window (first minima)")
        bAuto.clicked.connect(self._auto_window)
        self.jitter = QDoubleSpinBox(); self.jitter.setRange(1.0, 50.0)
        self.jitter.setValue(10.0); self.jitter.setSuffix(" %")
        self.jitter.setToolTip("how far each window edge is jittered to get the "
                               "sensitivity of δCG — replaces integrating by "
                               "hand three times")
        self.jitter.valueChanged.connect(self._measure)
        lv.addWidget(_row(bAuto, "  edge jitter", self.jitter))
        self.p_measure = _plot("central band", ppm_axis=True, parent=self)
        lv.addWidget(self.p_measure, 1)
        self.region = pg.LinearRegionItem(
            brush=pg.mkBrush(*theme._rgb(theme.active().measure), 40),
            pen=pg.mkPen(theme.active().measure, width=1))
        self.region.setZValue(-10)
        self.region.sigRegionChangeFinished.connect(self._window_dragged)
        self.p_measure.addItem(self.region)
        self.cgLine = pg.InfiniteLine(
            angle=90, movable=False,
            pen=pg.mkPen(theme.active().pivot, width=1.6),
            label="δCG {value:0.1f}",
            labelOpts={"position": 0.9, "color": theme.active().pivot})
        self.p_measure.addItem(self.cgLine)
        self.lbl6 = QLabel(""); self.lbl6.setWordWrap(True)
        self.lbl6.setTextFormat(Qt.RichText)
        lv.addWidget(self.lbl6)
        return w

    # ----------------------------------------------------------------- load
    def _pick(self):
        p, _ = QFileDialog.getOpenFileName(
            self, "QCPMG echo-train FID  (pick the raw 'fid' file)", "",
            "Bruker FID (fid ser);;All files (*)")
        if not p:
            p = QFileDialog.getExistingDirectory(
                self, "…or the EXPNO folder containing it")
        if p:
            self._load(p)

    def _load(self, source):
        from larmor import qcpmg
        from larmor.io import bruker

        try:
            src = Path(source)
            # an EXPNO folder resolves to the PROCESSED data, which is not an
            # echo train -- always prefer the raw fid inside it
            if src.is_dir() and (src / "fid").exists():
                src = src / "fid"
            d = bruker.read(str(src))
            if d.domain != "time" or d.ndim != 1:
                raise ValueError(
                    "that is not a 1D echo-train FID — open the raw 'fid' "
                    "file of the QCPMG EXPNO, not its processed data")
        except Exception as exc:                              # noqa: BLE001
            self.res.setText(f"cannot read: {exc}")
            return

        self._loading = True
        try:
            self.source = str(src)
            self.lbl.setText(str(src))
            self.fid = np.asarray(d.data, complex)
            self.meta = dict(d.meta)
            # a fresh dataset starts from a clean slate FIRST, before anything
            # that can raise: a failure below must not leave the old sample's
            # phase/exclusions/window attached to the new sample's fid
            self.excluded.clear(); self.keep = None
            self.t2 = None
            self._phased = False; self._window_user = False
            self._cg = self._sigma = self._fwhm = None
            self._ppm = self._spec = self._spec_raw = None
            self._spk_ppm = self._spk = None
            for w in (self.lb, self.gb, self.p0, self.p1):
                w.blockSignals(True); w.setValue(0.0); w.blockSignals(False)
            self.dropFirst.blockSignals(True)
            self.dropFirst.setValue(0)
            self.dropFirst.blockSignals(False)

            self._carrier, self._referenced = qcpmg.carrier_ppm(self.meta)
            sw = float(self.meta.get("sw_Hz", 0.0) or 0.0)

            per_f, self._period_src = qcpmg.echo_period_from_meta(
                self.meta, n_points=self.fid.size)
            per = int(round(per_f)) if per_f >= 4 else 0
            if per < 4:
                per = qcpmg.detect_period(self.fid)
                self._period_src = "autocorrelation" if per else "none"
            per = max(4, per or 4)
            self.period.setValue(per)
            self.periodHz.setValue(sw / per if per and sw else 0.0)
            ech = qcpmg.split_echoes(self.fid, per)
            self.nEch.setMaximum(max(1, ech.shape[0]))
            # ALL echoes by default: the fitted constant absorbs the noise
            # floor, and truncating the tail biases T2 low (measured -15% when
            # cut at the 3-sigma echo). n_usable is shown as advice instead.
            self.nEch.setValue(ech.shape[0])
            self.dropFirst.setMaximum(max(0, ech.shape[0] - 2))
            self._n_usable = qcpmg.n_usable_echoes(ech)
            self.top.setMaximum(per - 1)
            self.top.setValue(qcpmg.echo_top_point(ech))
        except Exception as exc:                              # noqa: BLE001
            # never leave the PREVIOUS sample's spectrum sendable under the
            # NEW sample's name -- the spectra were nulled above, so all that
            # is left is to make the actions unreachable
            for w in (self.btnCsv, self.btnSend):
                w.setEnabled(False)
            self.res.setText(f"cannot process: {exc}")
            return
        finally:
            self._loading = False
        for w in (self.btnCsv, self.btnSend):
            w.setEnabled(True)
        self._recompute()
        # phase once automatically: on an UNPHASED spectrum the tallest feature
        # is a noise sliver, and stage 6 would report a delta_CG from it
        if self._spec_raw is not None:
            self._autophase()

    # ------------------------------------------------------------ plumbing
    def done(self, r: int):
        """Every way the dialog closes (X, Esc, Send) funnels through here:
        remember the size/position the user chose for the next open."""
        if not os.environ.get("LARMOR_NO_SESSION"):
            try:
                from PySide6.QtCore import QSettings
                QSettings("LARMOR", "app").setValue("qcpmgDialogGeometry",
                                                    self.saveGeometry())
            except Exception:
                pass
        super().done(r)

    def _queue(self, *_):
        if not self._loading and self.fid is not None:
            self._timer.start()

    def _flush(self):
        """Run any recompute still sitting in the 150 ms debounce NOW.

        Send/Copy inside the debounce window otherwise ships the STALE
        spectrum while the provenance records the NEW widget values -- the
        exported meta would lie about how the data was processed."""
        if self._timer.isActive():
            self._timer.stop()
            self._recompute()

    def _apply_period(self, per: int, sw: float):
        """A period change re-splits the whole train: every echo-indexed
        setting from the OLD split (exclusions, echo count) points at a
        different stretch of data under the new one, so they are reset --
        keeping excluded={3} across a period change would silently mask a
        different piece of the FID."""
        # write the QUANTISED spacing back, so the field never shows a
        # spacing the integer period cannot actually realise
        self.periodHz.blockSignals(True)
        self.periodHz.setValue(sw / per if per and sw else 0.0)
        self.periodHz.blockSignals(False)
        self.top.setMaximum(max(0, per - 1))
        self.excluded.clear()
        self.keep = None
        if self.fid is not None and per >= 4:
            n = max(1, self.fid.size // per)
            for w in (self.nEch, self.dropFirst):
                w.blockSignals(True)
            self.nEch.setMaximum(n)
            self.nEch.setValue(n)
            self.dropFirst.setMaximum(max(0, n - 2))
            self.dropFirst.setValue(0)
            for w in (self.nEch, self.dropFirst):
                w.blockSignals(False)
        self._period_src = "manual"
        self._queue()

    def _on_period_pts(self, val):
        if self._loading:
            return
        self._apply_period(int(val), float(self.meta.get("sw_Hz", 0.0) or 0.0))

    def _on_period_hz(self, val):
        if self._loading:
            return
        sw = float(self.meta.get("sw_Hz", 0.0) or 0.0)
        if val and sw:
            per = max(4, int(round(sw / val)))
            if per == self.period.value():
                # same integer period: still snap the field to the spacing
                # that period actually realises, then stop
                self.periodHz.blockSignals(True)
                self.periodHz.setValue(sw / per)
                self.periodHz.blockSignals(False)
                return
            self.period.blockSignals(True)
            self.period.setValue(per)
            self.period.blockSignals(False)
            self._apply_period(per, sw)

    def _on_top_spin(self, val):
        if self._loading:
            return
        self.topLine.blockSignals(True)
        self.topLine.setValue(float(val))
        self.topLine.blockSignals(False)
        self._queue()

    def _on_top_line(self, *_):
        val = int(round(float(self.topLine.value())))
        val = int(np.clip(val, 0, max(0, self.period.value() - 1)))
        self.topLine.blockSignals(True)
        self.topLine.setValue(float(val))
        self.topLine.blockSignals(False)
        if val != self.top.value():
            self.top.blockSignals(True)
            self.top.setValue(val)
            self.top.blockSignals(False)
            self._queue()

    def _auto_top(self):
        from larmor import qcpmg
        if self.fid is None:
            return
        self.top.setValue(qcpmg.echo_top_point(self._echoes()))

    def _use_matched(self):
        if self.t2 is not None and self.t2.ok:
            self.lb.setValue(round(self.t2.lb_Hz, 1))

    def _set_pstep(self, *_):
        s = self.pstep.value()
        self.p0.setSingleStep(s); self.p1.setSingleStep(s)

    def _echoes(self):
        from larmor import qcpmg
        return qcpmg.split_echoes(self.fid, self.period.value(),
                                  n_echoes=self.nEch.value(),
                                  drop_first=self.dropFirst.value())

    # ----------------------------------------------------------- recompute
    def _recompute(self, *_):
        from larmor import qcpmg
        if self.fid is None:
            return
        try:
            t = theme.active()
            sw = float(self.meta.get("sw_Hz", 0.0) or 0.0)
            sfo = float(self.meta.get("larmor_MHz", 0.0) or 0.0)
            per = self.period.value()
            ech = self._echoes()
            top = int(np.clip(self.top.value(), 0, per - 1))
            tau = per / sw if sw else 0.0

            # -- stage 1
            self.p_train.clear()
            mag = np.abs(self.fid)
            step = max(1, mag.size // 6000)
            self.p_train.plot(np.arange(0, mag.size, step), mag[::step],
                              pen=pg.mkPen(t.experiment, width=1))
            for k in range(1, min(ech.shape[0], 400) + 1):
                self.p_train.addItem(pg.InfiniteLine(
                    pos=k * per, angle=90,
                    pen=pg.mkPen(t.text_dim, width=0.6, style=Qt.DotLine)))
            full = qcpmg.split_echoes(self.fid, per)
            if self.nEch.maximum() != full.shape[0]:
                self.nEch.blockSignals(True)
                self.nEch.setMaximum(max(1, full.shape[0]))
                self.nEch.blockSignals(False)
            self._n_usable = qcpmg.n_usable_echoes(ech, top)
            align = qcpmg.split_alignment(ech)
            src = {"CNST7": "read from the pulse program (CNST7)",
                   "CNST8": "read from the pulse program (CNST8)",
                   "MASR": "assumed rotor-synchronised (MASR)",
                   "autocorrelation": "GUESSED from the autocorrelation — check the markers",
                   "manual": "set by hand",
                   "none": "unknown"}.get(self._period_src, self._period_src)
            warn = align < 0.35
            self.lbl1.setText(
                f"{ech.shape[0]} echoes × {per} pts · τecho = {tau * 1e6:,.1f} µs · "
                f"spikelet spacing {self.periodHz.value():,.1f} Hz "
                f"({self.periodHz.value() / sfo if sfo else 0:.2f} ppm) · "
                f"period {src} · alignment {align:.2f} · "
                f"signal above 3σ through echo "
                f"{getattr(self, '_n_usable', ech.shape[0])} "
                f"(all are kept — the fitted constant absorbs the noise floor)"
                + ("   ⚠ echoes do not line up — check the period"
                   if warn else ""))
            self.lbl1.setStyleSheet(f"color: {t.model if warn else t.text_dim};")

            # -- stage 2
            self.p_echoes.clear()
            self.p_echoes.addItem(self.topLine); self.p_echoes.addItem(self.centreLine)
            n_show = min(ech.shape[0], 60)
            for i in range(n_show):
                self.p_echoes.plot(np.abs(ech[i]),
                                   pen=pg.mkPen(t.experiment, width=1))
            self.centreLine.setValue(qcpmg.echo_centre(per))
            self.topLine.blockSignals(True)
            self.topLine.setValue(float(top))
            self.topLine.blockSignals(False)
            first, last = qcpmg.first_last_echo(ech)
            self.p_firstlast.clear()
            self.p_firstlast.plot(first, pen=pg.mkPen(t.experiment, width=1.4),
                                  name="first")
            self.p_firstlast.plot(last, pen=pg.mkPen(t.model, width=1.4),
                                  name="last")
            off = top - qcpmg.echo_centre(per)
            ok_top = abs(off) <= max(2, per // 50)
            self.lbl2.setText(
                f"echo top = {top} · block centre = {qcpmg.echo_centre(per)} · "
                f"offset {off:+d} pt — "
                + ("whole-echo condition met" if ok_top else
                   "the top is far from the centre; the transform will not be "
                   "pure absorption"))
            self.lbl2.setStyleSheet(f"color: {t.text_dim if ok_top else t.model};")

            # -- stage 3
            mode = "real" if self.decayMode.currentIndex() == 0 else "magnitude"
            decay = qcpmg.echo_decay(ech, top, mode)
            drop = self.dropFirst.value()
            self.keep = np.array([(drop + i) not in self.excluded
                                  for i in range(decay.size)], bool)
            tt = np.arange(decay.size) * tau
            # pass the TIMES explicitly: after an exclusion the kept points are
            # no longer 0, tau, 2*tau, ... and re-timing them shifts T2 badly
            self.t2 = qcpmg.fit_t2(tau, decay[self.keep], t_s=tt[self.keep],
                                   offset=self.useOffset.isChecked(), period=per)
            self._draw_decay(tau, decay)

            # -- stage 4
            t2s = self.t2.T2_s if (self.t2.ok and self.t2w.isChecked()) else None
            self.btnMatched.setEnabled(bool(self.t2.ok))
            self.btnMatched.setText(
                f"Use matched LB = {self.t2.lb_Hz:,.1f} Hz" if self.t2.ok
                else "Use matched LB  (T₂ fit failed)")
            apod = qcpmg.apodize_echoes(ech, tau, t2s)
            self.p_apod_echoes.clear()
            for i in range(min(apod.shape[0], 60)):
                self.p_apod_echoes.plot(np.abs(apod[i]),
                                        pen=pg.mkPen(t.experiment, width=1))
            w = qcpmg.echo_weights(ech.shape[0], tau, t2s)
            tt = np.arange(ech.shape[0]) * tau
            self.p_apod_decay.clear()
            self.p_apod_decay.plot(tt, np.abs(qcpmg.echo_decay(ech, top, "magnitude")
                                              / (np.abs(decay).max() or 1.0)),
                                   pen=pg.mkPen(t.text_dim, width=1))
            self.p_apod_decay.plot(tt, w, pen=pg.mkPen(t.model, width=1.6))
            eff = float(w.sum())
            self.lbl4.setText(
                f"matched filter {'ON' if t2s else 'off'} · "
                f"LB {self.lb.value():,.1f} Hz · GB {self.gb.value():,.1f} Hz · "
                f"effective echoes ≈ {eff:.1f} of {ech.shape[0]}")
            self.lbl4.setStyleSheet(f"color: {t.text_dim};")

            # -- stage 5
            zf = int(self.zf.currentText())
            self._ppm, self._spec_raw = qcpmg.sum_echo_spectrum(
                self.fid, per, sw, sfo, self._carrier, top=top,
                n_echoes=self.nEch.value(), drop_first=self.dropFirst.value(),
                t2_weight_s=t2s, lb_Hz=self.lb.value(), gb_Hz=self.gb.value(),
                zf=zf)
            self._spk_ppm, spk = qcpmg.spikelet_spectrum(
                self.fid, sw, sfo, self._carrier,
                lb_Hz=max(self.lb.value(), 1.0), zf=2)
            self._spk = np.abs(spk)
            for w in (self.btnCsv, self.btnSend):
                w.setEnabled(True)
            self._rephase()
        except Exception as exc:                              # noqa: BLE001
            _log.exception("QCPMG recompute failed")
            # do not leave a stale spectrum reachable behind a new, bad setting
            self._ppm = self._spec = self._spec_raw = None
            self._cg = self._sigma = self._fwhm = None
            for w in (self.btnCsv, self.btnSend):
                w.setEnabled(False)
            self.lbl1.setText(f"cannot process: {exc}")
            self.lbl1.setStyleSheet(f"color: {theme.active().model};")
            self._update_headline()

    def _draw_decay(self, tau, decay):
        t = theme.active()
        self.p_decay.clear()
        n = decay.size
        tt = np.arange(n) * tau
        keep = self.keep if (self.keep is not None and self.keep.size == n) \
            else np.ones(n, bool)
        scat = pg.ScatterPlotItem(
            x=tt, y=decay, size=9, pen=None,
            brush=[pg.mkBrush(t.pivot) if k else pg.mkBrush(t.text_dim)
                   for k in keep],
            symbol=["o" if k else "x" for k in keep],
            data=list(range(n)))
        scat.sigClicked.connect(self._decay_clicked)
        self.p_decay.addItem(scat)
        if self.t2 is not None and self.t2.ok:
            grid = np.linspace(0, tt.max() if n else 1.0, 400)
            self.p_decay.plot(grid, self.t2.model(grid),
                              pen=pg.mkPen(t.model, width=1.6))
            f = self.t2
            err = (f"± {f.T2_err_s * 1e3:.3f}" if np.isfinite(f.T2_err_s) else "")
            rel = (100.0 * f.T2_err_s / f.T2_s
                   if f.T2_s and np.isfinite(f.T2_err_s) else float("nan"))
            n_dec = f.T2_s / tau if tau else float("nan")
            # a decay sampled by only a couple of echoes is not a measurement,
            # however confident the fitted number looks
            poor = f.r2 < 0.7 or (np.isfinite(rel) and rel > 25.0)
            warn = ""
            if poor:
                warn = (f"<br><span style='color:{t.model}'>⚠ poorly determined "
                        f"(R² {f.r2:.2f}, ±{rel:.0f}%): the signal decays in only "
                        f"{n_dec:.1f} echoes, so few points carry real signal. "
                        f"Treat T₂ — and the matched LB derived from it — as "
                        f"indicative; a shorter echo period would measure it "
                        f"properly.</span>")
            self.lbl3.setText(
                f"<b>T₂ = {f.T2_s * 1e3:.3f} {err} ms</b> "
                f"(R² = {f.r2:.4f}, {int(keep.sum())} of {n} echoes, "
                f"decays in {n_dec:.1f} echoes) &nbsp;·&nbsp; "
                f"matched apodization <b>LB = {f.lb_Hz:,.1f} Hz</b> = 1/πT₂"
                f"<br><span style='color:{t.text_dim}'>on ssNake's D1 pseudo-axis "
                f"(1 dwell per echo): T₂ = {f.T2_ssnake_s:.4e} s, "
                f"LB = {f.lb_ssnake_Hz:,.0f} Hz — same filter, axis differs by "
                f"the echo length ({self.period.value()})</span>" + warn)
        else:
            self.lbl3.setText(
                f"<span style='color:{t.model}'>T₂ fit did not converge — "
                f"no matched filter offered. Try the signed-real decay, a "
                f"different echo top, or fewer echoes.</span>")

    def _decay_clicked(self, _scatter, points):
        if not points:
            return
        i = points[0].data()
        if not isinstance(i, int):
            return
        abs_i = self.dropFirst.value() + i          # index in the whole train
        self.excluded.symmetric_difference_update({abs_i})
        self._recompute()

    def _window_dragged(self, *_):
        """The user placed the CG window by hand -- from now on a phase nudge
        or a recompute must not silently re-seed it."""
        self._window_user = True
        self._measure()

    # -------------------------------------------------------------- phase
    def _rephase(self, *_):
        from larmor import qcpmg
        if self._spec_raw is None:
            return
        if self.p0.value() or self.p1.value():
            self._phased = True
        self._spec = qcpmg.phase_spectrum(
            self._spec_raw, self.p0.value(), self.p1.value()).real
        self._redraw_spec()
        if self._window_user:
            self._measure()          # keep the window the user placed
        else:
            self._auto_window()

    def _redraw_spec(self, *_):
        from larmor import qcpmg
        if self._spec is None:
            return
        t = theme.active()
        mode = ["max", "area", "raw"][self.norm.currentIndex()]
        a, b = qcpmg.overlay_pair(self._spec, self._spk, mode)
        self.p_spec.clear()
        shown = []
        if self.showSpk.isChecked():
            self.p_spec.plot(self._spk_ppm, b, pen=pg.mkPen(t.model, width=1.0))
            shown.append("spikelets")
        if self.showSum.isChecked():
            self.p_spec.plot(self._ppm, a, pen=pg.mkPen(t.experiment, width=1.6))
            shown.insert(0, "sum echo")
        ref = ("" if self._referenced else
               f"<span style='color:{t.model}'> · carrier from O1/BF1 — no "
               f"procs, the ppm axis may be offset</span>")
        self.lbl5.setText(
            f"{' ⊕ '.join(shown) or 'nothing shown'} · scale {mode} · "
            f"zero-fill ×{self.zf.currentText()} · p0 {self.p0.value():.1f}° "
            f"p1 {self.p1.value():.1f}° · carrier {self._carrier:.2f} ppm{ref}")
        self._update_headline()

    def _autophase(self):
        from larmor import qcpmg
        if self._spec_raw is None:
            return
        p0, p1 = qcpmg.autophase(self._spec_raw)
        for w, val in ((self.p0, p0), (self.p1, p1)):
            w.blockSignals(True); w.setValue(val); w.blockSignals(False)
        self._phased = True
        self._rephase()

    # ------------------------------------------------------------ measure
    def _auto_window(self, *_):
        from larmor import qcpmg
        if self._spec is None:
            return
        hi, lo = qcpmg.cg_window(self._ppm, self._spec)
        if not (np.isfinite(hi) and np.isfinite(lo)) or hi <= lo:
            # a NaN region is un-grabbable, so 'drag it onto the band' would
            # be impossible advice -- seed the middle third instead
            span = float(self._ppm.max() - self._ppm.min())
            mid = float(self._ppm.min()) + span / 2.0
            lo, hi = mid - span / 6.0, mid + span / 6.0
        self.region.blockSignals(True)
        self.region.setRegion((lo, hi))
        self.region.blockSignals(False)
        self._measure()

    def _measure(self, *_):
        from larmor import qcpmg
        if self._spec is None:
            return
        t = theme.active()
        self.p_measure.clear()
        self.p_measure.addItem(self.region); self.p_measure.addItem(self.cgLine)
        self.p_measure.plot(self._ppm, self._spec,
                            pen=pg.mkPen(t.experiment, width=1.4))
        lo, hi = self.region.getRegion()
        if not self._phased:
            self._cg = self._sigma = self._fwhm = None
            self.cgLine.setVisible(False)      # no stale number on the plot
            self.lbl6.setText(
                f"<span style='color:{t.model}'>phase the spectrum first "
                f"(stage 5 → Autophase) — δCG and FWHM measured on an unphased "
                f"spectrum are meaningless.</span>")
            self._update_headline()
            return
        cg, sigma = qcpmg.centre_of_gravity(self._ppm, self._spec, (hi, lo),
                                            jitter_frac=self.jitter.value() / 100.0)
        fw = qcpmg.fwhm_hz(self._ppm, self._spec,
                           float(self.meta.get("larmor_MHz", 0.0) or 0.0), (hi, lo))
        if not np.isfinite(cg):
            self._cg = self._sigma = self._fwhm = None
            self.cgLine.setVisible(False)      # no stale number on the plot
            self.lbl6.setText(
                f"<span style='color:{t.model}'>the window contains no usable "
                f"signal — drag it onto the central band.</span>")
            self._update_headline()
            return
        self.cgLine.setValue(cg)
        self.cgLine.setVisible(True)
        self._cg, self._sigma, self._fwhm = cg, sigma, fw
        sloppy = sigma > 8.0
        sfo = float(self.meta.get("larmor_MHz", 0.0) or 1.0)
        self.lbl6.setText(
            f"<b>δCG = {cg:.1f} ± {sigma:.1f} ppm</b> &nbsp;·&nbsp; "
            f"<b>FWHM = {fw:,.0f} Hz</b> ({fw / sfo:.1f} ppm) &nbsp;·&nbsp; "
            f"window {max(lo, hi):.1f} … {min(lo, hi):.1f} ppm"
            + (f"<br><span style='color:{t.model}'>the window is sensitive — "
               f"drag its edges onto the first minima either side of the "
               f"central band</span>" if sloppy else ""))
        self._update_headline()

    def _update_headline(self):
        bits = []
        if self.t2 is not None and self.t2.ok:
            bits.append(f"T₂ {self.t2.T2_s * 1e3:.2f} ms · LB {self.t2.lb_Hz:,.0f} Hz")
        if getattr(self, "_cg", None) is not None:
            bits.append(f"δCG {self._cg:.1f} ± {self._sigma:.1f} ppm")
            bits.append(f"FWHM {self._fwhm:,.0f} Hz")
        self.res.setText("   ·   ".join(bits))

    # --------------------------------------------------------------- output
    def _copy_csv(self):
        from PySide6.QtWidgets import QApplication
        from larmor import qcpmg
        self._flush()
        # _spec_raw is None whenever the current settings could not be
        # processed (failed load or failed recompute) -- there is nothing
        # whose provenance could honestly be copied, and re-splitting the fid
        # below would raise the very error that put us in this state
        if self.fid is None or self._spec_raw is None:
            return
        per = self.period.value()
        sw = float(self.meta.get("sw_Hz", 0.0) or 0.0)
        tau = per / sw if sw else 0.0
        ech = self._echoes()
        top = int(np.clip(self.top.value(), 0, per - 1))
        mode = "real" if self.decayMode.currentIndex() == 0 else "magnitude"
        decay = qcpmg.echo_decay(ech, top, mode)
        keep = self.keep if (self.keep is not None and
                             self.keep.size == decay.size) else np.ones(decay.size, bool)
        lines = ["echo,t_s,intensity,kept"]
        for i, y in enumerate(decay):
            lines.append(f"{i},{i * tau:.9g},{float(y):.9g},{int(keep[i])}")
        f = self.t2
        lines += [
            f"# source={self.source}",
            f"# period_pts={per}  source={self._period_src}  tau_echo_s={tau:.9g}",
            f"# echo_top={top}  echoes={ech.shape[0]}  drop_first={self.dropFirst.value()}",
            f"# decay_mode={mode}  offset_fitted={self.useOffset.isChecked()}",
        ]
        if f is not None and f.ok:
            lines += [
                f"# T2_physical_s={f.T2_s:.9g}  T2_physical_err_s={f.T2_err_s:.9g}",
                f"# LB_physical_Hz={f.lb_Hz:.9g}  R2={f.r2:.6g}",
                f"# T2_ssnake_D1_s={f.T2_ssnake_s:.9g}  LB_ssnake_D1_Hz={f.lb_ssnake_Hz:.9g}",
            ]
        lines += [
            f"# LB_applied_Hz={self.lb.value():.9g}  GB_applied_Hz={self.gb.value():.9g}",
            f"# t2_weighted={self.t2w.isChecked()}  zerofill={self.zf.currentText()}",
            f"# p0_deg={self.p0.value():.9g}  p1_deg={self.p1.value():.9g}",
            f"# carrier_ppm={self._carrier:.9g}  referenced={self._referenced}",
        ]
        if getattr(self, "_cg", None) is not None:
            lines += [f"# delta_CG_ppm={self._cg:.9g}  sigma_ppm={self._sigma:.9g}",
                      f"# FWHM_Hz={self._fwhm:.9g}"]
        QApplication.clipboard().setText("\n".join(lines))
        self._update_headline()
        self.res.setText(self.res.text() + "   ·   copied CSV")

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "qcpmg", "QCPMG processing guide")

    def _send(self):
        self._flush()
        if self._ppm is None or self._spec is None:
            return
        meta = {
            "expno": self.meta.get("expno", ""),
            "title": (self.meta.get("title", "").splitlines() or ["QCPMG"])[0]
            + " (QCPMG)",
            "nucleus": self.meta.get("nucleus", ""),
            "larmor_MHz": self.meta.get("larmor_MHz", 0.0),
            # a QCPMG sum-echo spectrum has no MAS sideband manifold to model;
            # forwarding a stale acqus MASR would silently add sidebands
            "spin_rate_Hz": 0.0,
            "mas_uncertain": False,
            "qcpmg_period_pts": self.period.value(),
            "qcpmg_echo_top": self.top.value(),
            "qcpmg_n_echoes": self.nEch.value(),
            "qcpmg_lb_Hz": self.lb.value(),
            "qcpmg_gb_Hz": self.gb.value(),
            "qcpmg_p0_deg": self.p0.value(),
            "qcpmg_p1_deg": self.p1.value(),
            "qcpmg_carrier_ppm": self._carrier,
            "qcpmg_referenced": self._referenced,
        }
        if self.t2 is not None and self.t2.ok:
            meta.update({"qcpmg_T2_s": self.t2.T2_s,
                         "qcpmg_T2_err_s": self.t2.T2_err_s,
                         "qcpmg_matched_lb_Hz": self.t2.lb_Hz})
        if getattr(self, "_cg", None) is not None:
            meta.update({"qcpmg_delta_cg_ppm": self._cg,
                         "qcpmg_delta_cg_sigma_ppm": self._sigma,
                         "qcpmg_fwhm_Hz": self._fwhm})
        self.accepted_1d.emit(np.asarray(self._ppm), np.asarray(self._spec), meta)
        self.accept()
