"""Manuals and tutorials in non-modal help windows, with a text-size bar.

``show_help(parent, name, title, kind)`` is the one entry point -- the Help
menu, the command palette (which triggers the same menu actions) and every
dialog's Help button go through it. The page opens as a non-modal tool
window (``windowtray.TOOL_WINDOW_FLAGS``: minimise / maximise / close
buttons) parented to the caller's window, so the workbench -- or the tool
dialog the button sits in -- stays usable while the manual is read. Open
pages are kept in the owner window's ``_help_windows`` dict keyed by
``(kind, name)``: asking for a page that is already open raises the existing
window; closing a page frees it (``WA_DeleteOnClose``).

Text size: A- / A+ / Reset in the window's toolbar (Ctrl+-, Ctrl+= or Ctrl++,
Ctrl+0; Ctrl + mouse wheel) re-render the page through
``mdrender.render_help_html(md, scale)`` so the typeset equations grow with
the text. The factor is one per application (module state, persisted as
QSettings ``helpTextScale`` unless LARMOR_NO_SESSION), applied to every open
help window at once and to every window opened afterwards. The window
geometry is remembered the same way (``helpWindowGeometry``).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from PySide6.QtCore import QEvent, QSettings, QSize, Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (QDialog, QHBoxLayout, QLabel, QPushButton,
                               QSizePolicy, QTextBrowser, QToolBar,
                               QVBoxLayout, QWidget)

from larmor.desktop.windowtray import TOOL_WINDOW_FLAGS, is_alive

__all__ = ["help_path", "tutorial_path", "show_help", "HelpWindow",
           "ZOOM_STEPS", "text_scale", "set_text_scale"]

#: the text-size factors A- / A+ step through (1.0 = the application font)
ZOOM_STEPS = (0.7, 0.8, 0.9, 1.0, 1.15, 1.3, 1.5, 1.75, 2.0, 2.5)
_SCALE_KEY = "helpTextScale"
_GEOMETRY_KEY = "helpWindowGeometry"

_scale: float | None = None          # the session's factor, read lazily
_open: list = []                     # every HelpWindow alive


def help_path(name: str) -> Path | None:
    """Locate larmor/help/<name>.md, whether from source or a frozen build."""
    here = Path(__file__).resolve()
    for base in (here.parent.parent / "help",                 # larmor/help
                 Path(getattr(sys, "_MEIPASS", "")) / "larmor" / "help"):
        p = base / f"{name}.md"
        if p.is_file():
            return p
    return None


def tutorial_path(name: str) -> Path | None:
    """Locate docs/tutorials/<name>.md, whether from source or a frozen build
    (packaging/larmor.spec bundles the folder as docs/tutorials)."""
    here = Path(__file__).resolve()
    for base in (here.parents[2] / "docs" / "tutorials",
                 Path(getattr(sys, "_MEIPASS", "")) / "docs" / "tutorials"):
        p = base / f"{name}.md"
        if p.is_file():
            return p
    return None


# --------------------------------------------------------------- settings
def _settings() -> QSettings:
    return QSettings("LARMOR", "app")


def _persist() -> bool:
    return not os.environ.get("LARMOR_NO_SESSION")


def text_scale() -> float:
    """The current help text-size factor (from QSettings on first use)."""
    global _scale
    if _scale is None:
        _scale = 1.0
        if _persist():
            try:
                v = float(_settings().value(_SCALE_KEY, 1.0) or 1.0)
                if ZOOM_STEPS[0] <= v <= ZOOM_STEPS[-1]:
                    _scale = v
            except (TypeError, ValueError):
                pass
    return _scale


def set_text_scale(scale: float) -> float:
    """Set the factor for the session, persist it, re-render every open help
    window; returns the clamped value."""
    global _scale
    _scale = float(min(max(float(scale), ZOOM_STEPS[0]), ZOOM_STEPS[-1]))
    if _persist():
        try:
            _settings().setValue(_SCALE_KEY, _scale)
        except Exception:
            pass
    for w in list(_open):
        if is_alive(w):
            w.apply_scale(_scale)
    return _scale


def _step(scale: float, direction: int) -> float:
    """The next factor above (+1) or below (-1) ``scale``."""
    if direction > 0:
        for s in ZOOM_STEPS:
            if s > scale + 1e-6:
                return s
        return ZOOM_STEPS[-1]
    for s in reversed(ZOOM_STEPS):
        if s < scale - 1e-6:
            return s
    return ZOOM_STEPS[0]


# ----------------------------------------------------------------- window
class HelpWindow(QDialog):
    """One manual / tutorial page: a QTextBrowser under a text-size bar."""

    def __init__(self, parent, name: str, title: str, kind: str, markdown: str):
        super().__init__(parent)
        self.page = (kind, name)
        self.markdown = markdown
        self.scale = text_scale()
        self.html = ""
        self._registry: dict | None = None
        self.setWindowTitle(title)
        self.setWindowFlags(TOOL_WINDOW_FLAGS)
        self.setModal(False)
        self.setAttribute(Qt.WA_DeleteOnClose, True)
        self.setSizeGripEnabled(True)

        v = QVBoxLayout(self)
        bar = QToolBar()
        bar.setMovable(False)
        bar.setIconSize(QSize(14, 14))
        bar.setToolButtonStyle(Qt.ToolButtonTextOnly)
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        bar.addWidget(spacer)
        bar.addWidget(QLabel("Text size "))
        self.actSmaller = QAction("A−", self)
        self.actSmaller.setShortcuts([QKeySequence("Ctrl+-"),
                                      QKeySequence(QKeySequence.ZoomOut)])
        self.actSmaller.setToolTip("smaller text (Ctrl+−, or Ctrl + mouse wheel)")
        self.actSmaller.triggered.connect(self.zoom_out)
        self.actLarger = QAction("A+", self)
        self.actLarger.setShortcuts([QKeySequence("Ctrl+="), QKeySequence("Ctrl++"),
                                     QKeySequence(QKeySequence.ZoomIn)])
        self.actLarger.setToolTip("larger text (Ctrl+=, or Ctrl + mouse wheel)")
        self.actLarger.triggered.connect(self.zoom_in)
        self.actReset = QAction("Reset", self)
        self.actReset.setShortcut(QKeySequence("Ctrl+0"))
        self.actReset.setToolTip("default text size (Ctrl+0)")
        self.actReset.triggered.connect(self.zoom_reset)
        bar.addAction(self.actSmaller)
        bar.addAction(self.actLarger)
        bar.addAction(self.actReset)
        self.lblScale = QLabel("100 %")
        self.lblScale.setToolTip("the text size of every help window; kept "
                                 "for the next launch")
        bar.addWidget(self.lblScale)
        v.addWidget(bar)

        self.tb = QTextBrowser()
        self.tb.setOpenExternalLinks(True)
        self.tb.setStyleSheet("QTextBrowser { background: #ffffff; padding: 6px 10px; }")
        self.tb.viewport().installEventFilter(self)     # Ctrl + wheel
        v.addWidget(self.tb, 1)
        row = QHBoxLayout()
        row.addStretch(1)
        btn = QPushButton("Close")
        btn.setAutoDefault(False)
        btn.clicked.connect(self.accept)
        row.addWidget(btn)
        v.addLayout(row)

        self._render()
        restored = False
        if _persist():
            try:
                geo = _settings().value(_GEOMETRY_KEY)
                if geo is not None:
                    restored = bool(self.restoreGeometry(geo))
            except Exception:
                restored = False
        if not restored:
            self.resize(880, 720)
        _open.append(self)
        self.destroyed.connect(lambda *_: self._forget())

    # ------------------------------------------------------------ zooming
    def zoom_in(self):
        set_text_scale(_step(self.scale, +1))

    def zoom_out(self):
        set_text_scale(_step(self.scale, -1))

    def zoom_reset(self):
        set_text_scale(1.0)

    def apply_scale(self, scale: float):
        """Re-render at ``scale`` (called for every open window by
        ``set_text_scale``)."""
        if abs(scale - self.scale) < 1e-9 and self.html:
            return
        self.scale = float(scale)
        self._render()

    def eventFilter(self, obj, ev):
        if ev.type() == QEvent.Wheel and (ev.modifiers() & Qt.ControlModifier):
            dy = ev.angleDelta().y()
            if dy > 0:
                self.zoom_in()
            elif dy < 0:
                self.zoom_out()
            return True                 # not the browser's own font zoom
        return super().eventFilter(obj, ev)

    def _render(self):
        vb = self.tb.verticalScrollBar()
        frac = vb.value() / vb.maximum() if vb.maximum() > 0 else 0.0
        try:
            from PySide6.QtWidgets import QApplication

            from larmor.desktop.mdrender import help_css, render_help_html

            html = render_help_html(self.markdown, self.scale)
            doc = self.tb.document()
            f = QApplication.font()
            if f.pointSizeF() > 0:
                f.setPointSizeF(f.pointSizeF() * self.scale)
            doc.setDefaultFont(f)
            doc.setDefaultStyleSheet(help_css(self.scale))
            self.tb.setHtml(html)                  # typeset the equations
            self.html = html
        except Exception:
            try:
                self.tb.setMarkdown(self.markdown)
            except Exception:
                self.tb.setPlainText(self.markdown)
        self.lblScale.setText(f"{round(self.scale * 100):d} %")
        self.actSmaller.setEnabled(self.scale > ZOOM_STEPS[0] + 1e-6)
        self.actLarger.setEnabled(self.scale < ZOOM_STEPS[-1] - 1e-6)
        if frac > 0:
            def back():
                if is_alive(self):
                    sb = self.tb.verticalScrollBar()
                    sb.setValue(int(round(frac * sb.maximum())))
            back()
            QTimer.singleShot(0, back)

    # ----------------------------------------------------------- plumbing
    def _forget(self):
        _open[:] = [w for w in _open if w is not self]
        reg = self._registry
        if reg is not None:
            for k, w in list(reg.items()):
                if w is self:
                    reg.pop(k, None)

    def done(self, r: int):
        """Every way the window closes (X, Esc, Close) funnels through here:
        remember the geometry for the next help window, leave the registry."""
        if _persist():
            try:
                _settings().setValue(_GEOMETRY_KEY, self.saveGeometry())
            except Exception:
                pass
        self._forget()
        super().done(r)


def show_help(parent, name: str, title: str = "Help",
              kind: str = "manual") -> HelpWindow:
    """Open (or raise) the manual / tutorial ``name`` in a non-modal window
    owned by ``parent``'s window; ``kind`` is "manual" (larmor/help) or
    "tutorial" (docs/tutorials). Returns the window."""
    p = tutorial_path(name) if kind == "tutorial" else help_path(name)
    missing = "Tutorial not found." if kind == "tutorial" else "Manual not found."
    text = (p.read_text(encoding="utf-8") if p
            else f"# {name}\n\n{missing}")
    owner = None
    if parent is not None:
        try:
            owner = parent.window()
        except RuntimeError:
            owner = None
    reg = None
    if owner is not None:
        reg = getattr(owner, "_help_windows", None)
        if reg is None:
            reg = {}
            owner._help_windows = reg
    key = (kind, name)
    existing = reg.get(key) if reg is not None else None
    if existing is not None and is_alive(existing):
        if existing.windowState() & Qt.WindowMinimized:
            existing.setWindowState(existing.windowState() & ~Qt.WindowMinimized)
        existing.show()
        existing.raise_()
        existing.activateWindow()
        return existing
    w = HelpWindow(owner, name, title, kind, text)
    if reg is not None:
        reg[key] = w
        w._registry = reg
    w.show()
    w.raise_()
    w.activateWindow()
    return w
