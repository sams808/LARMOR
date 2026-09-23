"""The spinning-sideband OFFER: a small non-modal banner over the spectrum
canvas, with the dotted comb guides it draws on the same view.

``larmor.sidebands.detect`` finds the +-nu_rot repeat; this widget shows the
claim where it can be checked -- the guides sit on the very peaks the
detector matched -- and carries the two one-click actions (a linked manifold
of lines, or a shifted copy of the spectrum) plus, in the manifold button's
drop-down, the empirical ``sidebands`` model line. It is a child of the
SpectrumView, positioned top-centre under the title by an event filter (the
``set_placeholder`` idiom), so plot.py needs no edit and the banner hides
with its view when a 2D document takes the stack. It never steals focus;
Return (with the banner focused) adds the manifold, Escape dismisses.

Constructible with a bare SpectrumView, so it can be tested without a
MainWindow. Physics and the decision live in larmor.sidebands; the recipe
edits live in MainWindow (larmor/desktop/app.py).
"""
from __future__ import annotations

import pyqtgraph as pg
from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtWidgets import (QFrame, QHBoxLayout, QLabel, QMenu, QPushButton,
                               QToolButton, QVBoxLayout)

from larmor import sidebands
from larmor.desktop import theme

#: y of the banner inside the canvas: under the plot title, clear of the
#: top-left legend
BANNER_TOP = 36
#: the sentence wraps beyond this width (px)
LABEL_MAX_WIDTH = 600


class SidebandBanner(QFrame):
    """'Spinning sidebands: nu_rot ... [Add linked manifold v] [Add shifted
    copy] [x]' -- hidden until :meth:`set_detection`."""

    manifold = Signal()        # [Add linked manifold] / Return
    model_line = Signal()      # drop-down: 'as one `sidebands` model line'
    copy = Signal()            # [Add shifted copy]
    dismissed = Signal()       # [x] / Escape

    def __init__(self, view):
        super().__init__(view)
        self.view = view
        self.det = None
        self.guides: list = []
        self.setObjectName("ssbBanner")
        self.setFocusPolicy(Qt.StrongFocus)
        self.setAttribute(Qt.WA_StyledBackground, True)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 6, 8, 6)
        outer.setSpacing(4)
        self.label = QLabel("")
        self.label.setWordWrap(True)
        self.label.setMaximumWidth(LABEL_MAX_WIDTH)
        self.label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        outer.addWidget(self.label)

        row = QHBoxLayout()
        row.setSpacing(6)
        row.addStretch(1)
        self.btnManifold = QToolButton(self)
        self.btnManifold.setText("Add linked manifold")
        self.btnManifold.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self.btnManifold.setPopupMode(QToolButton.MenuButtonPopup)
        self.btnManifold.setToolTip(
            "one Gauss/Lorentz line per detected order, position = centre ± "
            "k·νrot and shape tied to the centreband, amplitudes free "
            "(sideband intensities follow the tensor, not a ratio). Return "
            "does the same.")
        menu = QMenu(self.btnManifold)
        self.actModelLine = menu.addAction(
            "as one `sidebands` model line  (geometric ratio)")
        self.actModelLine.setToolTip(
            "the empirical dmfit-style manifold: one line whose n_ssb and "
            "ratio r are seeded from the detected teeth")
        self.btnManifold.setMenu(menu)
        self.btnManifold.clicked.connect(lambda *_: self.manifold.emit())
        self.actModelLine.triggered.connect(lambda *_: self.model_line.emit())
        self.btnCopy = QPushButton("Add shifted copy", self)
        self.btnCopy.setToolTip(
            "one copy of the spectrum on screen per detected side, shifted by "
            "±νrot (shift held) and scaled to the first sideband -- the copy "
            "carries the whole pattern, higher orders included")
        self.btnCopy.clicked.connect(lambda *_: self.copy.emit())
        self.btnClose = QToolButton(self)
        self.btnClose.setText("✕")
        self.btnClose.setAutoRaise(True)
        self.btnClose.setToolTip(
            "Not sidebands — dismiss (Esc). Turn the offer off under View ▸ "
            "Offer spinning-sideband detection on load")
        self.btnClose.clicked.connect(lambda *_: self.dismissed.emit())
        row.addWidget(self.btnManifold)
        row.addWidget(self.btnCopy)
        row.addWidget(self.btnClose)
        outer.addLayout(row)

        self._apply_theme()
        view.installEventFilter(self)
        self.hide()

    # ------------------------------------------------------------ look
    def _apply_theme(self):
        t = theme.active()
        self.setStyleSheet(
            f"QFrame#ssbBanner {{ background: {t.base}; color: {t.text}; "
            f"border: 1px solid {t.accent}; border-radius: 6px; }} "
            f"QFrame#ssbBanner QLabel {{ color: {t.text}; border: none; "
            f"background: transparent; }}")

    def eventFilter(self, obj, ev):
        if obj is self.view and ev.type() == QEvent.Resize:
            self._reposition()
        return False                          # never swallow the view's events

    def _reposition(self):
        self.adjustSize()
        x = max(0, (self.view.width() - self.width()) // 2)
        self.move(x, BANNER_TOP)

    def showEvent(self, ev):
        super().showEvent(ev)
        self._reposition()

    # ------------------------------------------------------------ content
    def set_detection(self, det, recipe_rate_Hz: float = 0.0):
        """Show ``det`` in words (sidebands.describe) and remember it."""
        self.det = det
        self._apply_theme()
        self.label.setText(sidebands.describe(det, recipe_rate_Hz))
        self._reposition()

    def show_guides(self, det):
        """A vertical guide at the centreband (solid) and at every MATCHED
        order (dotted), labelled 'centre', '+1', '−1', ..., in the theme's
        pivot colour behind data and paddles."""
        self.clear_guides()
        pv = theme.active().pivot
        marks = [(det.centre_ppm, "centre", Qt.SolidLine)]
        marks += [(o.ppm, f"{o.k:+d}".replace("-", "−"), Qt.DotLine)
                  for o in det.matched()]
        for pos, text, style in marks:
            line = pg.InfiniteLine(
                pos=float(pos), angle=90, movable=False,
                pen=pg.mkPen(pv, width=1, style=style),
                label=text, labelOpts={"color": pv, "position": 0.95,
                                       "movable": False})
            line.setZValue(-15)
            self.view.addItem(line)
            self.guides.append(line)

    def clear_guides(self):
        for g in self.guides:
            try:
                self.view.removeItem(g)
            except Exception:
                pass
        self.guides = []

    def dismiss(self):
        """Hide, drop the guides and forget the detection (no signal)."""
        self.hide()
        self.clear_guides()
        self.det = None

    # ------------------------------------------------------------ keys
    def keyPressEvent(self, ev):
        if ev.key() in (Qt.Key_Return, Qt.Key_Enter):
            ev.accept()
            self.manifold.emit()
            return
        if ev.key() == Qt.Key_Escape:
            ev.accept()
            self.dismissed.emit()
            return
        super().keyPressEvent(ev)
