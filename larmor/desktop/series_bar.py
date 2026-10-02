"""The series bar: one row above the spectrum that drives the series mode.

``[Series] [ s0 ● | s1 ● | s2 ● … ] [◀] [▶] [Carry ▾] [Fit → next ▶]
[Auto sweep ▾] [Table…] [Plot…] [Acquisition…] [Save ▾] [⚠ chip] [✕ End]``

Every spectrum of a series is an ordinary workspace of the main window
(``larmor.desktop.mw_series``); this widget only shows the members with a
status dot (grey not fitted · green fitted · amber edited since its fit ·
red the last fit failed), walks them, and emits what the user asked for.
It holds no model state: the window's series mixin feeds ``set_members`` /
``set_options`` / ``set_comparison`` and listens to the signals. Hidden
until a series is active; constructible on its own, so tests drive it
without a MainWindow.
"""
from __future__ import annotations

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtGui import QAction, QColor, QFontMetrics, QIcon, QPainter, QPixmap
from PySide6.QtWidgets import (QButtonGroup, QComboBox, QFrame, QHBoxLayout,
                               QLabel, QMenu, QPushButton, QScrollArea,
                               QSizePolicy, QToolButton, QWidget, QWidgetAction)

from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS
from larmor.seriesmode import PASS_CHOICES, SMOOTH_CHOICES, START_CHOICES

__all__ = ["SeriesBar", "status_color", "STATUS_WORDS"]

#: what the status dot means, for tooltips
STATUS_WORDS = {"unfitted": "not fitted yet", "fitted": "fitted",
                "edited": "edited since its fit (F5 refits)",
                "failed": "the last fit failed"}
#: a member button's text is elided beyond this width (px); the tooltip has it all
NAME_MAX_PX = 150


def status_color(status: str) -> str:
    """The dot colour of a member status from the active theme's signal
    roles: the fit-health strip's green-teal (measure) for fitted, amber
    (baseline) for edited, red (model) for failed, the dim ink for not yet."""
    t = theme.active()
    return {"fitted": t.measure, "edited": t.baseline,
            "failed": t.model}.get(status, t.text_dim)


def _dot_icon(color: str, size: int = 10) -> QIcon:
    pm = QPixmap(size, size)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing, True)
    p.setPen(Qt.NoPen)
    p.setBrush(QColor(color))
    p.drawEllipse(1, 1, size - 2, size - 2)
    p.end()
    return QIcon(pm)


class SeriesBar(QFrame):
    member_clicked = Signal(int)
    prev_requested = Signal()
    next_requested = Signal()
    fit_next_requested = Signal()
    sweep_requested = Signal(int, str, int)      # (passes, start, smooth)
    stop_requested = Signal()
    seed_on_move_toggled = Signal(bool)
    carry_changed = Signal(object)               # tuple of parameter names ticked
    carry_menu_opening = Signal()                # refresh the checklist before it shows
    copy_model_requested = Signal(bool)          # True = replace every member's lines
    table_requested = Signal()
    plot_requested = Signal()
    acquisition_requested = Signal()
    save_all_requested = Signal()
    bundle_requested = Signal()
    details_requested = Signal()
    end_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("seriesBar")
        self.setAttribute(Qt.WA_StyledBackground, True)
        self._members: list = []
        self._buttons: list = []
        self._current = -1
        self._carry_actions: dict = {}
        self._candidates: tuple = ()
        self._running = False

        h = QHBoxLayout(self)
        h.setContentsMargins(8, 3, 8, 3)
        h.setSpacing(6)

        self.lblTitle = QLabel("Series")
        self.lblTitle.setStyleSheet("font-weight:600;")
        self.lblTitle.setToolTip(
            "every spectrum of the series is a workspace of its own: add or "
            "remove lines, set bounds and links, process, fit and undo on each "
            "one with the usual tools; this bar walks the series and runs the "
            "sweeps. Dots: grey not fitted · green fitted · amber edited since "
            "its fit · red the last fit failed")
        h.addWidget(self.lblTitle)

        # the member strip: checkable buttons in a horizontally scrolling area
        self.strip = QScrollArea()
        self.strip.setWidgetResizable(True)
        self.strip.setFrameShape(QFrame.NoFrame)
        self.strip.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.strip.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.strip.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._stripBody = QWidget()
        self._stripLay = QHBoxLayout(self._stripBody)
        self._stripLay.setContentsMargins(0, 0, 0, 0)
        self._stripLay.setSpacing(2)
        self._stripLay.addStretch(1)
        self.strip.setWidget(self._stripBody)
        self._group = QButtonGroup(self)
        self._group.setExclusive(True)
        self._group.buttonClicked.connect(self._member_button_clicked)
        h.addWidget(self.strip, 1)

        self.btnPrev = self._tool("◀", "previous spectrum of the series (seeds it "
                                       "from this one when 'seed on move' is on)")
        self.btnPrev.clicked.connect(self.prev_requested)
        self.btnNext = self._tool("▶", "next spectrum of the series (seeds it from "
                                       "this one when 'seed on move' is on)")
        self.btnNext.clicked.connect(self.next_requested)
        h.addWidget(self.btnPrev)
        h.addWidget(self.btnNext)

        # Carry ▾
        self.btnCarry = self._tool("Carry ▾", "what moves from one spectrum to the "
                                              "next, and the copy-to-all actions")
        self.btnCarry.setPopupMode(QToolButton.InstantPopup)
        self.mCarry = QMenu(self.btnCarry)
        self.mCarry.setToolTipsVisible(True)
        self.mCarry.aboutToShow.connect(self.carry_menu_opening)
        self.actSeedOnMove = QAction("Seed the next spectrum from this one when I move", self)
        self.actSeedOnMove.setCheckable(True)
        self.actSeedOnMove.setChecked(True)
        self.actSeedOnMove.setToolTip(
            "on ◀ ▶ (and after Fit → next) the spectrum you land on takes this "
            "one's model: a spectrum without lines gets a copy (amplitudes scaled "
            "to its height), one with its own lines keeps them and only receives "
            "the ticked parameters on matching lines")
        self.actSeedOnMove.toggled.connect(self.seed_on_move_toggled)
        self.actCarryHead = QAction("Carried parameters:", self)
        self.actCarryHead.setEnabled(False)
        self.actCopyModel = QAction("Copy this model to every spectrum", self)
        self.actCopyModel.setToolTip(
            "spectra without lines get a copy of this model (amplitudes scaled "
            "to each height); spectra with their own lines keep them and receive "
            "the carried values on matching lines")
        self.actCopyModel.triggered.connect(lambda: self.copy_model_requested.emit(False))
        self.actReplaceModel = QAction("Replace every spectrum's lines with this model", self)
        self.actReplaceModel.setToolTip(
            "every other spectrum drops its own lines and takes a copy of this "
            "model (amplitudes scaled to each height) — their fits start over")
        self.actReplaceModel.triggered.connect(lambda: self.copy_model_requested.emit(True))
        self._rebuild_carry_menu()
        self.btnCarry.setMenu(self.mCarry)
        h.addWidget(self.btnCarry)

        self.btnFitNext = self._tool("Fit → next ▶",
                                     "fit this spectrum (the usual Fit: progress bar, "
                                     "Stop / Cancel, animation), then move to the "
                                     "next one and seed it from the result")
        self.btnFitNext.clicked.connect(self.fit_next_requested)
        h.addWidget(self.btnFitNext)

        # Auto sweep ▾ with its small form, progress text and Stop
        self.btnSweep = self._tool("Auto sweep ▾",
                                   "fit every spectrum in turn, each one starting "
                                   "from its fitted neighbour, forward and back")
        self.btnSweep.setPopupMode(QToolButton.InstantPopup)
        self.mSweep = QMenu(self.btnSweep)
        form = QWidget()
        fl = QHBoxLayout(form)
        fl.setContentsMargins(10, 6, 10, 6)
        fl.setSpacing(6)
        fl.addWidget(QLabel("passes"))
        self.cbPasses = QComboBox()
        self.cbPasses.addItems([str(p) for p in PASS_CHOICES])
        self.cbPasses.setCurrentText("2")
        self.cbPasses.setToolTip("each pass sweeps one direction: 2 = forward then "
                                 "back, 4 = F/B/F/B …")
        fl.addWidget(self.cbPasses)
        fl.addWidget(QLabel("start"))
        self.cbStart = QComboBox()
        self.cbStart.addItems(list(START_CHOICES))
        self.cbStart.setToolTip("which end of the series the first pass starts from")
        fl.addWidget(self.cbStart)
        fl.addWidget(QLabel("smooth"))
        self.cbSmooth = QComboBox()
        self.cbSmooth.addItems([str(s) for s in SMOOTH_CHOICES])
        self.cbSmooth.setToolTip("moving-average window applied to each parameter's "
                                 "trajectory between passes (0 = off) so the series "
                                 "does not jitter")
        fl.addWidget(self.cbSmooth)
        self.btnRun = QPushButton("Run")
        self.btnRun.setDefault(True)
        self.btnRun.setToolTip("start the sweep; Stop keeps what was fitted, Cancel "
                               "(beside the progress bar) reverts everything")
        self.btnRun.clicked.connect(self._run_clicked)
        fl.addWidget(self.btnRun)
        wa = QWidgetAction(self.mSweep)
        wa.setDefaultWidget(form)
        self.mSweep.addAction(wa)
        self.btnSweep.setMenu(self.mSweep)
        h.addWidget(self.btnSweep)
        self.lblProgress = QLabel("")
        self.lblProgress.setStyleSheet("font-weight:600;")
        self.lblProgress.setVisible(False)
        h.addWidget(self.lblProgress)
        self.btnStop = self._tool("⏹ Stop", "stop the sweep after the current spectrum "
                                            "and keep what was fitted")
        self.btnStop.clicked.connect(self.stop_requested)
        self.btnStop.setVisible(False)
        h.addWidget(self.btnStop)

        # outputs
        self.btnTable = self._tool("Table…", "the Series table: names, replicate groups, "
                                             "the order of the series and composition "
                                             "columns joined from a CSV (OK renames and "
                                             "reorders the members)")
        self.btnTable.clicked.connect(self.table_requested)
        self.btnPlot = self._tool("Plot…", "the Series plot: every parameter and "
                                           "population across the series, with export")
        self.btnPlot.clicked.connect(self.plot_requested)
        self.btnAcq = self._tool("Acquisition…", "Table S1 and the Experimental "
                                                 "paragraph for these spectra, from "
                                                 "acqus / procs / title — no fit needed")
        self.btnAcq.clicked.connect(self.acquisition_requested)
        self.btnSave = self._tool("Save ▾", "save every member's fit, or write the "
                                            "publication bundle of the series")
        self.btnSave.setPopupMode(QToolButton.InstantPopup)
        self.mSave = QMenu(self.btnSave)
        self.mSave.setToolTipsVisible(True)
        a = self.mSave.addAction("Save all fits…")
        a.setToolTip("one .recipe.json per spectrum into a folder — automatic names "
                     "(sample_nucleus_seq_YYYYMMDD_HHMM) or one typed per fit")
        a.triggered.connect(self.save_all_requested)
        a = self.mSave.addAction("Publication bundle…")
        a.setToolTip("everything about this series in one folder: seq_table.csv, one "
                     ".recipe.json and one _curves.csv per spectrum, manifest.csv "
                     "(source file + SHA-256, EXPNO, NS, D1, SF/SR, processing, fit "
                     "window, RMSD) and README.txt with the Methods paragraph")
        a.triggered.connect(self.bundle_requested)
        self.btnSave.setMenu(self.mSave)
        for b in (self.btnTable, self.btnPlot, self.btnAcq, self.btnSave):
            h.addWidget(b)

        # comparability chip (acquired / processed alike?)
        self.chip = self._tool("", "")
        self.chip.clicked.connect(self.details_requested)
        self.chip.setVisible(False)
        h.addWidget(self.chip)

        self.btnEnd = self._tool("✕ End series", "remove this bar and the series tags; "
                                                 "the spectra stay open")
        self.btnEnd.clicked.connect(self.end_requested)
        h.addWidget(self.btnEnd)

        self.apply_theme()
        self.hide()

    # ------------------------------------------------------------ building
    @staticmethod
    def _tool(text: str, tip: str) -> QToolButton:
        b = QToolButton()
        b.setText(text)
        b.setToolTip(tip)
        b.setAutoRaise(True)
        b.setToolButtonStyle(Qt.ToolButtonTextOnly)
        b.setFocusPolicy(Qt.NoFocus)
        return b

    def apply_theme(self) -> None:
        t = theme.active()
        self.setStyleSheet(
            f"#seriesBar {{ background: {t.surface}; border-bottom: 1px solid {t.border}; }}")
        self.strip.setStyleSheet("QScrollArea { background: transparent; }")
        self._stripBody.setStyleSheet("background: transparent;")
        self._paint_members()
        if self.chip.isVisible():
            self.set_comparison(getattr(self, "_cmp", None))

    def _rebuild_carry_menu(self) -> None:
        m = self.mCarry
        m.clear()
        m.addAction(self.actSeedOnMove)
        m.addSeparator()
        m.addAction(self.actCarryHead)
        for name in self._candidates:
            m.addAction(self._carry_actions[name])
        if not self._candidates:
            a = m.addAction("(no lines on the series yet)")
            a.setEnabled(False)
        m.addSeparator()
        m.addAction(self.actCopyModel)
        m.addAction(self.actReplaceModel)

    # ------------------------------------------------------------ members
    def set_members(self, members: list) -> None:
        """``members``: ``[{"name", "status", "tip", "flag"}, …]`` in series
        order. Buttons are updated in place when the count is unchanged."""
        members = [dict(m) for m in members]
        if len(members) != len(self._buttons):
            for b in self._buttons:
                self._group.removeButton(b)
                self._stripLay.removeWidget(b)
                b.deleteLater()
            self._buttons = []
            for k in range(len(members)):
                b = QToolButton()
                b.setCheckable(True)
                b.setAutoRaise(True)
                b.setFocusPolicy(Qt.NoFocus)
                b.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
                b.setIconSize(QSize(10, 10))
                self._group.addButton(b, k)
                self._stripLay.insertWidget(k, b)
                self._buttons.append(b)
        self._members = members
        self._paint_members()
        if not (0 <= self._current < len(self._buttons)):
            self._current = -1

    def _paint_members(self) -> None:
        fm = QFontMetrics(self.font())
        for k, (b, m) in enumerate(zip(self._buttons, self._members)):
            name = str(m.get("name") or f"spectrum {k + 1}")
            text = fm.elidedText(name, Qt.ElideMiddle, NAME_MAX_PX)
            if m.get("flag"):
                text += " ⚠"
            b.setText(text)
            st = str(m.get("status") or "unfitted")
            b.setIcon(_dot_icon(status_color(st)))
            tip = str(m.get("tip") or name)
            b.setToolTip(f"{tip}\n{STATUS_WORDS.get(st, st)} — click to switch to it")
            f = b.font()
            f.setBold(k == self._current)
            b.setFont(f)

    def set_current(self, k: int) -> None:
        self._current = int(k) if k is not None else -1
        self._group.blockSignals(True)
        try:
            if 0 <= self._current < len(self._buttons):
                self._buttons[self._current].setChecked(True)
                self.strip.ensureWidgetVisible(self._buttons[self._current], 24, 0)
            else:
                checked = self._group.checkedButton()
                if checked is not None:
                    self._group.setExclusive(False)
                    checked.setChecked(False)
                    self._group.setExclusive(True)
        finally:
            self._group.blockSignals(False)
        self._paint_members()
        n = len(self._buttons)
        self.btnPrev.setEnabled(not self._running and n > 0 and self._current != 0)
        self.btnNext.setEnabled(not self._running and n > 0 and self._current != n - 1)

    def _member_button_clicked(self, btn) -> None:
        k = self._group.id(btn)
        if k >= 0:
            self.member_clicked.emit(int(k))

    def member_buttons(self) -> list:
        return list(self._buttons)

    def current(self) -> int:
        return self._current

    def statuses(self) -> list:
        return [str(m.get("status") or "unfitted") for m in self._members]

    # ------------------------------------------------------------ options
    def set_options(self, seed_on_move: bool, candidates, carried, passes: int,
                    start: str, smooth: int) -> None:
        """Mirror the series options: the seed-on-move toggle, the Carry
        checklist (``candidates`` offered, ``carried`` ticked, labels from
        PARAM_LABELS) and the sweep form."""
        self.actSeedOnMove.blockSignals(True)
        self.actSeedOnMove.setChecked(bool(seed_on_move))
        self.actSeedOnMove.blockSignals(False)
        carried = set(carried or ())
        self._candidates = tuple(candidates or ())
        for name in self._candidates:
            a = self._carry_actions.get(name)
            if a is None:
                a = QAction(PARAM_LABELS.get(name, name), self)
                a.setCheckable(True)
                a.setToolTip(f"carry {name} from one spectrum to the next")
                a.toggled.connect(self._carry_toggled)
                self._carry_actions[name] = a
            a.blockSignals(True)
            a.setChecked(name in carried)
            a.blockSignals(False)
        self._rebuild_carry_menu()
        for cb, val in ((self.cbPasses, passes), (self.cbStart, start), (self.cbSmooth, smooth)):
            i = cb.findText(str(val))
            if i >= 0:
                cb.setCurrentIndex(i)

    def carry(self) -> tuple:
        return tuple(n for n in self._candidates if self._carry_actions[n].isChecked())

    def carry_candidates(self) -> tuple:
        return tuple(self._candidates)

    def _carry_toggled(self, *_):
        self.carry_changed.emit(self.carry())

    def sweep_settings(self) -> tuple:
        return (int(self.cbPasses.currentText()), self.cbStart.currentText(),
                int(self.cbSmooth.currentText()))

    def _run_clicked(self):
        self.mSweep.close()
        p, s, m = self.sweep_settings()
        self.sweep_requested.emit(p, s, m)

    # ------------------------------------------------------------ comparability
    def set_comparison(self, cmp) -> None:
        """The comparability verdict as a chip: hidden when there is nothing
        to compare, '✓ alike' dim, '⚠ check' amber, '✖ mixed' red; the full
        sentence in the tooltip, the Details dialog on click."""
        self._cmp = cmp
        level = getattr(cmp, "level", "none") if cmp is not None else "none"
        if level == "none":
            self.chip.setVisible(False)
            return
        t = theme.active()
        summary = ""
        try:
            summary = cmp.summary()
        except Exception:                                 # noqa: BLE001
            pass
        if level == "ok":
            self.chip.setText("✓ alike")
            self.chip.setStyleSheet(f"color:{t.text_dim};")
        elif level == "bad":
            self.chip.setText("✖ mixed")
            self.chip.setStyleSheet(
                f"background:{t.model}; color:{theme.best_text_on(t.model)}; "
                "font-weight:600; padding:1px 6px; border-radius:3px;")
        else:
            self.chip.setText("⚠ check")
            self.chip.setStyleSheet(
                f"background:{t.baseline}; color:{theme.best_text_on(t.baseline)}; "
                "font-weight:600; padding:1px 6px; border-radius:3px;")
        self.chip.setToolTip((summary + "\n\n" if summary else "")
                             + "click for every acqus / procs parameter of every "
                             "spectrum (the deviants in amber)")
        self.chip.setVisible(True)

    # ------------------------------------------------------------ sweep state
    def set_sweep_running(self, on: bool) -> None:
        self._running = bool(on)
        for w in (self.btnCarry, self.btnFitNext, self.btnSweep, self.btnTable,
                  self.btnSave, self.btnEnd):
            w.setEnabled(not on)
        for b in self._buttons:
            b.setEnabled(not on)
        self.lblProgress.setVisible(on)
        self.btnStop.setVisible(on)
        if not on:
            self.lblProgress.setText("")
        self.set_current(self._current)

    def set_sweep_progress(self, text: str) -> None:
        self.lblProgress.setText(str(text))

    def is_running(self) -> bool:
        return self._running
