"""Non-modal tool windows and the "open windows" bar of the left sidebar.

Window furniture from the 0.14.2 field-use feedback, all about secondary
windows getting in the way of -- or getting lost behind -- the workbench:

* ``show_tool_window(dlg)`` shows a read-only viewer or a calculator dialog
  (the NMR table, a manual, a correlation matrix, the Czjzek distribution ...)
  as a NON-modal top-level window with minimise / maximise / close buttons,
  kept alive on its owner window and freed when closed, so the rest of the
  application stays usable while it is open. Dialogs whose caller reads a
  value back (``if dlg.exec(): ...``) stay modal and never come through here.

* ``OpenWindowsBar`` -- a compact strip at the bottom of the left sidebar
  listing EVERY open tool window of the session, oldest first, as a small
  button with the elided title (full title in the tooltip) under a header
  line "3 windows open"; a second QCPMG dialog opened by mistake is visible
  at a glance. Clicking a button raises the window (or restores it when it
  is minimised); a right-click offers Restore / Minimise / Close / Close all;
  the strip hides itself when no tool window is open.

  Tracking is an application-level event filter (installed once by
  ``install_tray`` from ``_build_sidebar``), so no dialog registers by hand:
  QEvent.Show adds a top-level, non-modal Window / Dialog other than the main
  window; Close, a programmatic Hide and ``destroyed`` drop it; tooltips,
  menus and popups (their window types) and non-window widgets (the co-fit
  split panels) are never listed; a modal dialog blocks the application while
  it is up, so it cannot be a forgotten window and is not listed either.

* Minimising. On Windows an owned top-level window that is minimised shrinks
  to a small floating title bar at the bottom-left of the SCREEN. The same
  filter watches WindowStateChange: a listed window that becomes minimised is
  hidden instead and its bar entry is drawn dimmed; the button (or the
  Restore item) brings it back with ``showNormal()`` / ``showMaximized()``.
  A minimised MODAL dialog would leave the application blocked behind that
  floating bar, so it is restored at once instead.

``tray_minimize(window)`` / ``tray_restore(window)`` drive the same code
without the window manager, for scripts and tests (offscreen Qt does not
deliver a WindowStateChange the way a real window manager does).
"""
from __future__ import annotations

import shiboken6
from PySide6.QtCore import QEvent, QObject, QSize, Qt, QTimer
from PySide6.QtGui import QFontMetrics
from PySide6.QtWidgets import (QApplication, QFrame, QLabel, QMenu,
                               QSizePolicy, QStyle, QToolButton, QVBoxLayout,
                               QWidget)

__all__ = ["TOOL_WINDOW_FLAGS", "is_alive", "show_tool_window",
           "OpenWindowsBar", "install_tray", "current_tray", "open_windows",
           "tray_minimize", "tray_restore"]

#: a real window (not a Qt.Dialog): title bar, system menu, minimise /
#: maximise / close buttons -- so it can be put aside while the main window
#: is used, and so the bar has a minimise to catch
TOOL_WINDOW_FLAGS = (Qt.Window | Qt.WindowTitleHint | Qt.WindowSystemMenuHint
                     | Qt.WindowMinMaxButtonsHint | Qt.WindowCloseButtonHint)

#: window types the bar lists; Tool (floating docks), Popup (menus, combo
#: drop-downs), ToolTip and SplashScreen windows are left alone
_LISTED_TYPES = (Qt.Window, Qt.Dialog)

#: parentless tool windows have no owner to be kept alive on
_orphans: list = []


def is_alive(obj) -> bool:
    try:
        return obj is not None and shiboken6.isValid(obj)
    except Exception:
        return False


# ----------------------------------------------------------- tool windows
def show_tool_window(dlg, owner=None):
    """Show ``dlg`` as a non-modal tool window and return it.

    ``owner`` (default: the window of the dialog's parent) keeps the dialog
    alive in its ``_tool_windows`` list until the dialog is closed; closing
    frees it (``WA_DeleteOnClose``). The dialog's own Close / accept / reject
    plumbing keeps working -- ``done()`` closes a non-modal QDialog too."""
    if owner is None:
        p = dlg.parentWidget()
        owner = p.window() if p is not None else None
    dlg.setWindowFlags(TOOL_WINDOW_FLAGS)
    if hasattr(dlg, "setModal"):
        dlg.setModal(False)
    dlg.setAttribute(Qt.WA_DeleteOnClose, True)
    keep = _orphans
    if owner is not None:
        keep = getattr(owner, "_tool_windows", None)
        if keep is None:
            keep = []
            owner._tool_windows = keep
    if not any(k is dlg for k in keep):
        keep.append(dlg)

        def forget(*_a, _keep=keep, _dlg=dlg):
            _keep[:] = [k for k in _keep if k is not _dlg]

        dlg.destroyed.connect(forget)
        if hasattr(dlg, "finished"):
            dlg.finished.connect(forget)
    dlg.show()
    dlg.raise_()
    dlg.activateWindow()
    return dlg


# --------------------------------------------------------------------- bar
class _Entry:
    __slots__ = ("window", "button", "title", "minimized", "maximized", "slot")

    def __init__(self, window, title: str):
        self.window = window
        self.title = title
        self.button = None
        self.minimized = False
        self.maximized = False
        self.slot = None


class _BarStrip(QWidget):
    """The strip at the bottom of the sidebar: a stretch, then a panel (a
    thin rule, the "N windows open" header, one button per window) that is
    hidden while no tool window is open. The strip itself stays as the
    sidebar's spacer."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("windowTray")
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addStretch(1)
        self.panel = QWidget()
        self.panel.setObjectName("windowTrayPanel")
        lay = QVBoxLayout(self.panel)
        lay.setContentsMargins(0, 4, 0, 0)
        lay.setSpacing(2)
        rule = QFrame()
        rule.setFrameShape(QFrame.HLine)
        rule.setFrameShadow(QFrame.Sunken)
        lay.addWidget(rule)
        self.header = QLabel("")
        self.header.setObjectName("windowTrayHeader")
        self.header.setAlignment(Qt.AlignHCenter)
        self.header.setWordWrap(True)
        lay.addWidget(self.header)
        self._lay = lay
        self._buttons: list[QToolButton] = []
        outer.addWidget(self.panel)
        self.panel.hide()

    def add_button(self, btn: QToolButton):
        self._buttons.append(btn)
        self._lay.addWidget(btn)

    def remove_button(self, btn: QToolButton):
        self._buttons = [b for b in self._buttons if b is not btn]
        # a tracked window's destroyed() can fire AFTER the main window (and
        # this strip's layout) is gone -- at application exit or when a test
        # closes the window first; nothing is left to update then
        if not shiboken6.isValid(self._lay) or not shiboken6.isValid(btn):
            return
        self._lay.removeWidget(btn)
        btn.hide()
        btn.setParent(None)
        btn.deleteLater()

    def buttons(self) -> list:
        return list(self._buttons)

    def set_count(self, n: int):
        if n <= 0:
            self.panel.hide()
            return
        self.header.setText(f"{n} window{'s' if n != 1 else ''} open")
        self.header.setToolTip(
            f"{n} tool window{'s' if n != 1 else ''} open in this session -- "
            "click one to bring it to the front; right-click for Restore / "
            "Minimise / Close / Close all")
        try:
            from larmor.desktop import theme
            self.header.setStyleSheet(
                f"color: {theme.active().text_dim}; font-size: 10px;")
        except Exception:
            pass
        self.panel.show()


class OpenWindowsBar(QObject):
    """Application-level event filter + the sidebar strip (see module doc)."""

    def __init__(self, main_window, toolbar=None, parent=None):
        super().__init__(parent if parent is not None else main_window)
        self._main = main_window
        self._toolbar = toolbar
        self._entries: list[_Entry] = []
        self.widget = _BarStrip()
        self.strip_action = None

    # ------------------------------------------------------------ queries
    def entries(self) -> list:
        """The open tool windows, oldest first (dead wrappers dropped)."""
        dead = [e for e in self._entries if not is_alive(e.window)]
        for e in dead:
            self._untrack(e, dead=True)
        return [e.window for e in self._entries]

    def count(self) -> int:
        return len(self.entries())

    def is_minimized(self, window) -> bool:
        e = self._entry(window)
        return bool(e is not None and e.minimized)

    def button_for(self, window) -> QToolButton | None:
        e = self._entry(window)
        return e.button if e is not None else None

    def is_shown(self) -> bool:
        """Whether the strip currently shows anything (False = no window)."""
        return not self.widget.panel.isHidden()

    def _entry(self, window) -> _Entry | None:
        for e in self._entries:
            if e.window is window:
                return e
        return None

    def _is_candidate(self, obj) -> bool:
        """A top-level Window / Dialog other than the main window."""
        if not isinstance(obj, QWidget) or not is_alive(obj):
            return False
        try:
            if not obj.isWindow() or obj is self._main:
                return False
            return obj.windowType() in _LISTED_TYPES
        except RuntimeError:            # wrapper outliving its C++ object
            return False

    # ------------------------------------------------------- event filter
    def eventFilter(self, obj, ev):
        t = ev.type()
        if t == QEvent.Show:
            if self._is_candidate(obj):
                self._shown(obj)
        elif t == QEvent.WindowStateChange:
            if self._is_candidate(obj):
                self._state_changed(obj)
        elif t == QEvent.Close:
            e = self._entry(obj)
            if e is not None:
                self._closed(e)
        elif t == QEvent.Hide:
            e = self._entry(obj)
            if e is not None and not e.minimized and not ev.spontaneous():
                self._untrack(e)          # hidden by code: not open any more
        elif t == QEvent.WindowTitleChange:
            e = self._entry(obj)
            if e is not None:
                e.title = obj.windowTitle() or e.title
                self._refresh(e)
        return False

    def _shown(self, w):
        e = self._entry(w)
        if e is None:
            if w.isModal():
                return                      # blocks the app: cannot be forgotten
            self._track(w)
            return
        if e.minimized:
            # app code re-shows a trayed window (a kept-alive dialog's
            # dlg.show()) or restore() runs: leave the minimised state, and
            # clear the bit BEFORE the platform show so it does not come
            # back as a floating title bar (QShowEvent precedes show_sys)
            e.minimized = False
            self._clear_minimized(w)
        self._refresh(e)

    def _state_changed(self, w):
        state = w.windowState()
        e = self._entry(w)
        if not (state & Qt.WindowMinimized):
            if e is not None and e.minimized:   # restored by other means
                e.minimized = False
                self._refresh(e)
            return
        if e is not None and e.minimized:
            return                          # already handled
        if w.isModal():
            # a minimised modal would leave the application blocked behind
            # a floating title bar: bring it straight back
            QTimer.singleShot(0, lambda: self._clear_minimized(w, show=True))
            return
        self.minimize(w)

    def _closed(self, e: _Entry):
        w = e.window
        self._untrack(e)
        self._clear_minimized(w)
        # a closeEvent may be vetoed (ignore()): the window is then still
        # visible a moment later and belongs back in the bar
        QTimer.singleShot(0, lambda: self._recheck(w))

    def _recheck(self, w):
        if (self._entry(w) is None and self._is_candidate(w)
                and not w.isModal() and w.isVisible()):
            self._track(w)

    @staticmethod
    def _clear_minimized(w, show: bool = False):
        if not is_alive(w):
            return
        state = w.windowState()
        if state & Qt.WindowMinimized:
            w.setWindowState(state & ~Qt.WindowMinimized)
        if show:
            if state & Qt.WindowMaximized:
                w.showMaximized()
            else:
                w.showNormal()
            w.raise_()
            w.activateWindow()

    # ------------------------------------------------------------ actions
    def activate(self, w):
        """The button click: restore a minimised window, raise a visible one."""
        e = self._entry(w)
        if e is not None and e.minimized:
            self.restore(w)
            return
        if is_alive(w):
            w.show()
            w.raise_()
            w.activateWindow()

    def minimize(self, w):
        """Hide ``w`` and draw its entry dimmed (a window not yet listed --
        shown before the bar existed -- is listed first). Idempotent."""
        if not is_alive(w):
            return
        e = self._entry(w)
        if e is None:
            e = self._track(w)
        if e.minimized:
            return
        e.minimized = True
        e.maximized = bool(w.windowState() & Qt.WindowMaximized)
        w.hide()
        self._refresh(e)

    def restore(self, w):
        """Bring a minimised window back (normal or maximised as it was)."""
        e = self._entry(w)
        if e is None or not is_alive(w):
            return
        e.minimized = False
        w.setWindowState(w.windowState() & ~Qt.WindowMinimized)
        if e.maximized:
            w.showMaximized()
        else:
            w.showNormal()
        w.raise_()
        w.activateWindow()
        self._refresh(e)

    def close_window(self, w):
        """Close a listed window (its entry goes with it)."""
        e = self._entry(w)
        if e is not None:
            self._untrack(e)
            self._clear_minimized(w)
        if is_alive(w):
            w.close()

    def close_all(self):
        for w in list(self.entries()):
            self.close_window(w)

    def context_menu(self, w) -> QMenu:
        """The Restore / Minimise / Close / Close all menu of an entry
        (built, not shown)."""
        e = self._entry(w)
        m = QMenu(self.widget)
        m.addAction("Restore").triggered.connect(lambda: self.activate(w))
        a_min = m.addAction("Minimise")
        a_min.triggered.connect(lambda: self.minimize(w))
        a_min.setEnabled(not (e is not None and e.minimized))
        m.addAction("Close").triggered.connect(lambda: self.close_window(w))
        m.addSeparator()
        m.addAction("Close all").triggered.connect(self.close_all)
        return m

    # ----------------------------------------------------------- plumbing
    def _track(self, w) -> _Entry:
        e = _Entry(w, w.windowTitle() or w.objectName() or "window")
        self._entries.append(e)
        e.button = self._make_button(e)
        self.widget.add_button(e.button)

        def gone(*_a, _e=e):
            self._untrack(_e, dead=True)

        e.slot = gone
        w.destroyed.connect(gone)
        self._refresh(e)
        return e

    def _untrack(self, e: _Entry, dead: bool = False):
        if e not in self._entries:
            return
        self._entries.remove(e)
        if e.button is not None:
            self.widget.remove_button(e.button)
            e.button = None
        if not dead and e.slot is not None and is_alive(e.window):
            try:
                e.window.destroyed.disconnect(e.slot)
            except (RuntimeError, TypeError):
                pass
        e.slot = None
        self.widget.set_count(len(self._entries))

    def _label_width(self) -> int:
        for wdg in (self.widget, self._toolbar):
            if wdg is not None and is_alive(wdg) and wdg.width() > 24:
                return max(40, wdg.width() - 10)
        return 64

    def _refresh(self, e: _Entry):
        """Text, dimming and tooltip of one entry + the header count."""
        btn = e.button
        if btn is None or not is_alive(btn):
            return
        btn.setText(QFontMetrics(btn.font()).elidedText(
            e.title, Qt.ElideRight, self._label_width()))
        # dynamic properties (for tests and stylesheets); "windowTitle" and
        # "minimized" are QWidget's own Q_PROPERTYs, hence the prefix
        btn.setProperty("trayTitle", e.title)
        btn.setProperty("trayMinimized", bool(e.minimized))
        if e.minimized:
            try:
                from larmor.desktop import theme
                dim = theme.active().text_dim
            except Exception:
                dim = None
            btn.setStyleSheet(
                ("QToolButton#trayButton { font-style: italic; "
                 + (f"color: {dim}; " if dim else "") + "}"))
            btn.setToolTip(f"{e.title}\n(minimised) click to restore; "
                           "right-click for Restore / Minimise / Close")
        else:
            btn.setStyleSheet("")
            btn.setToolTip(f"{e.title}\nclick to bring to the front; "
                           "right-click for Restore / Minimise / Close")
        self.widget.set_count(len(self._entries))

    def _make_button(self, e: _Entry) -> QToolButton:
        btn = QToolButton()
        btn.setObjectName("trayButton")
        btn.setAutoRaise(True)
        btn.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        btn.setIcon(btn.style().standardIcon(QStyle.SP_TitleBarNormalButton))
        btn.setIconSize(QSize(12, 12))
        btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        btn.clicked.connect(lambda *_: self.activate(e.window))
        btn.setContextMenuPolicy(Qt.CustomContextMenu)
        btn.customContextMenuRequested.connect(
            lambda pos, b=btn, w=e.window: self.context_menu(w).exec(
                b.mapToGlobal(pos)))
        return btn


# -------------------------------------------------------------- singleton
_bar: OpenWindowsBar | None = None


def install_tray(main_window, toolbar) -> OpenWindowsBar:
    """Create the bar for ``main_window``, put its strip at the bottom of the
    (vertical) sidebar ``toolbar`` and install the filter on the application.
    A previous bar's filter is removed first, so exactly one is active
    (tests build several main windows in one process)."""
    global _bar
    app = QApplication.instance()
    old = _bar
    if old is not None and is_alive(old) and app is not None:
        try:
            app.removeEventFilter(old)
        except RuntimeError:
            pass
    bar = OpenWindowsBar(main_window, toolbar)
    bar.strip_action = toolbar.addWidget(bar.widget)   # a text-less widget action
    if app is not None:
        app.installEventFilter(bar)
    _bar = bar
    return bar


def current_tray() -> OpenWindowsBar | None:
    return _bar if is_alive(_bar) else None


def open_windows() -> list:
    """The tool windows the active bar lists, oldest first."""
    bar = current_tray()
    return bar.entries() if bar is not None else []


def tray_minimize(window) -> bool:
    """Hide ``window`` into the active bar (dimmed entry, no floating bar)."""
    bar = current_tray()
    if bar is None:
        return False
    bar.minimize(window)
    return bar.is_minimized(window)


def tray_restore(window) -> bool:
    """Restore ``window`` from the active bar; False when it was not minimised."""
    bar = current_tray()
    if bar is None or not bar.is_minimized(window):
        return False
    bar.restore(window)
    return True
