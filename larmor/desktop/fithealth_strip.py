"""The fit-health strip: one row under the 1D workbench that renders a
``fithealth.Health`` — a coloured verdict pill, a neutral RMSD chip and one
chip per flag — with a details menu (F7) that also lists the analysis tools.

The strip never touches the recipe: every click is navigation (a signal the
main window connects to the residual toggle, a table cell, the correlation /
errors dialogs, the Report dock or the manual). Colours come from the active
theme's contrast-checked signal roles; a horizontally *Ignored* size policy
means the strip can never widen the window, however many flags it shows.
"""
from __future__ import annotations

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtWidgets import QHBoxLayout, QMenu, QSizePolicy, QToolButton, QWidget

from larmor import fithealth
from larmor.desktop import theme

#: chip wording is cut here (the full detail stays in the tooltip)
CHIP_CHARS = 48


def _elide(text: str, limit: int = CHIP_CHARS) -> str:
    return text if len(text) <= limit else text[:limit - 1] + "…"


class FitHealthStrip(QWidget):
    open_report = Signal()
    show_residual = Signal()
    focus_param = Signal(int, str)
    open_correlations = Signal()
    open_errors = Signal()
    open_help = Signal()
    #: the quantitativity chips: widen the integration window to the tails,
    #: Tools ▸ Relaxation on the sibling T1 EXPNO, the 90° pulse field of
    #: Experiment parameters
    widen_window = Signal()
    open_relaxation = Signal()
    enter_flip = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._health = None
        self._stale_hint = False
        self.setObjectName("fitHealthStrip")
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self._lay = QHBoxLayout(self)
        self._lay.setContentsMargins(6, 2, 6, 2)
        self._lay.setSpacing(6)
        self.pill = QToolButton()
        self.pill.setCursor(Qt.PointingHandCursor)
        self.pill.setFocusPolicy(Qt.StrongFocus)
        self.pill.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self._lay.addWidget(self.pill)
        self.rmsd_chip = QToolButton()
        self.rmsd_chip.setCursor(Qt.PointingHandCursor)
        self.rmsd_chip.setFocusPolicy(Qt.StrongFocus)
        self.rmsd_chip.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self.rmsd_chip.clicked.connect(lambda *_: self.open_report.emit())
        self._lay.addWidget(self.rmsd_chip)
        self._lay.addStretch(1)
        self.chips: list[QToolButton] = []
        self.set_health(None)

    # ------------------------------------------------------------ sizing
    def minimumSizeHint(self) -> QSize:        # Qt override: never wider than 120 px
        base = super().minimumSizeHint()
        return QSize(min(base.width(), 120), base.height())

    # ------------------------------------------------------------ state
    def health(self):
        return self._health

    def set_health(self, h) -> None:
        """Render ``h`` (None → the dim 'no fit yet' pill and no chips)."""
        self._health = h
        self._stale_hint = False
        self._render()

    def set_stale_hint(self, on: bool) -> None:
        """During a paddle drag only the pill wording changes ('edited since
        fit'); the chips wait for the full pass on release."""
        if self._health is None or self._stale_hint == bool(on):
            return
        self._stale_hint = bool(on)
        self.pill.setText(self._pill_text())

    def _pill_text(self) -> str:
        if self._health is None:
            return fithealth.NO_FIT_TEXT
        text = self._health.pill_text()
        if self._stale_hint and text != fithealth.NO_FIT_TEXT \
                and not text.endswith(fithealth.STALE_SUFFIX):
            text = text.replace(" Fit: ", " Model: ", 1)
            if text.startswith("✓ "):
                text = text[2:]
            text += fithealth.STALE_SUFFIX
        return text

    # ------------------------------------------------------------ colours
    def colours(self, level: str) -> tuple[str, str]:
        """(background, foreground) of the verdict pill for a level."""
        t = theme.active()
        if level == "bad":
            bg = t.model
        elif level == "check":
            bg = t.baseline
        elif level == "ok":
            bg = t.measure
        else:
            return t.surface, t.text_dim
        return bg, theme.best_text_on(bg)

    def apply_theme(self) -> None:
        """Re-derive every colour from the active theme, in place."""
        hint = self._stale_hint
        self._render()
        if hint:
            self.set_stale_hint(True)

    def _style_chip(self, chip: QToolButton, level: str, stale: bool = False) -> None:
        t = theme.active()
        ink = {"bad": t.model, "check": t.baseline, "info": t.text_dim}.get(level, t.text)
        border = ink if level in ("bad", "check") else t.border_soft
        if stale:
            ink, border = t.disabled_text, t.border_soft
        chip.setStyleSheet(
            f"QToolButton {{ background: {t.base}; color: {ink}; padding: 1px 6px; "
            f"border: 1px solid {border}; border-radius: 3px; }} "
            f"QToolButton:hover, QToolButton:focus {{ background: {t.hover}; "
            f"border-color: {t.accent}; }}")

    # ------------------------------------------------------------ render
    def _render(self) -> None:
        t = theme.active()
        self.setStyleSheet(
            f"#fitHealthStrip {{ background: {t.surface}; "
            f"border-top: 1px solid {t.border_soft}; }}")
        for c in self.chips:
            self._lay.removeWidget(c)
            c.deleteLater()
        self.chips = []
        h = self._health
        if h is None:
            level, tip = "none", "fit the lines (F5) to judge the model"
            self.rmsd_chip.setVisible(False)
        else:
            level, tip = h.level, h.tooltip()
            chi = h.chi_text()
            self.rmsd_chip.setVisible(bool(h.fitted and chi))
            self.rmsd_chip.setText(chi)
            self.rmsd_chip.setToolTip(
                (fithealth.STALE_PREFIX if h.stale else "")
                + "RMSD and reduced χ² of the last fit — click for the fit report")
            self._style_chip(self.rmsd_chip, "neutral", stale=h.stale)
            for f in h.flags:
                chip = self._make_chip(f)
                self.chips.append(chip)
                self._lay.insertWidget(self._lay.count() - 1, chip)
        bg, fg = self.colours(level)
        border = t.border_soft if level == "none" else bg
        self.pill.setText(self._pill_text())
        self.pill.setToolTip(tip)
        self.pill.setStyleSheet(
            f"QToolButton {{ background: {bg}; color: {fg}; font-weight: 600; "
            f"padding: 1px 8px; border: 1px solid {border}; border-radius: 3px; }} "
            f"QToolButton:hover, QToolButton:focus {{ border-color: {t.accent}; }}")

    def _make_chip(self, f) -> QToolButton:
        chip = QToolButton()
        chip.setToolButtonStyle(Qt.ToolButtonTextOnly)
        chip.setText(_elide(f.text))
        chip.setToolTip((fithealth.STALE_PREFIX if f.stale else "") + f.detail)
        chip.setCursor(Qt.PointingHandCursor if f.target else Qt.ArrowCursor)
        chip.setFocusPolicy(Qt.StrongFocus)
        chip.setProperty("flag_kind", f.kind)
        chip.clicked.connect(lambda *_, fl=f: self._activate(fl))
        self._style_chip(chip, f.level, stale=f.stale)
        return chip

    def _activate(self, f) -> None:
        if f.target == "residual":
            self.show_residual.emit()
        elif f.target == "param":
            if f.params:
                i, pname = f.params[0]
                self.focus_param.emit(int(i), str(pname))
        elif f.target == "correlations":
            self.open_correlations.emit()
        elif f.target == "errors":
            self.open_errors.emit()
        elif f.target == "report":
            self.open_report.emit()
        elif f.target == "widen":
            self.widen_window.emit()
        elif f.target == "relaxation":
            self.open_relaxation.emit()
        elif f.target == "flip":
            self.enter_flip.emit()

    # ------------------------------------------------------------ details
    def details_menu(self, extra_actions=()) -> QMenu:
        """Every flag as a line (triggering it is the chip click), then the
        SAME QAction objects the Decomposition menu holds, then the report and
        the manual."""
        h = self._health
        m = QMenu(self)
        kind = "Fit health" if (h is not None and h.fitted and not h.stale) else "Model health"
        m.addSection(f"{kind} — {self._pill_text()}")
        if h is not None and h.flags:
            for f in h.flags:
                a = m.addAction(f"{fithealth.GLYPH[f.level]} {f.text}")
                a.setToolTip(f.detail)
                a.setEnabled(bool(f.target))
                a.triggered.connect(lambda *_, fl=f: self._activate(fl))
        else:
            a = m.addAction("· no flags" if h is not None else "· not fitted yet")
            a.setEnabled(False)
        m.addSeparator()
        for act in extra_actions:
            m.addAction(act)
        if extra_actions:
            m.addSeparator()
        m.addAction("Show fit report").triggered.connect(lambda *_: self.open_report.emit())
        m.addAction("What the flags mean…").triggered.connect(lambda *_: self.open_help.emit())
        return m

    def show_details(self, extra_actions=()) -> None:
        m = self.details_menu(extra_actions)
        m.exec(self.pill.mapToGlobal(self.pill.rect().bottomLeft()))
        m.deleteLater()
