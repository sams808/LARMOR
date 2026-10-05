"""Shared plumbing of the desktop sweeps (``gui_menus`` / ``gui_tools``):
one offscreen QApplication under the default theme, a guard that puts the
user's saved settings back exactly as they were, a trap for exceptions
raised inside Qt slots (PySide6 only prints them) and for Qt's own
critical messages, a timer that dismisses every modal dialog the sweep
provokes, and an event pump.

The sweep runs under ``QT_QPA_PLATFORM=offscreen`` (set by the orchestrator
before Qt loads): no window shows, and file / colour / input dialogs are
plain widgets, so ``ModalDismisser`` can reject them -- the native Windows
dialogs could not be closed from code.
"""
from __future__ import annotations

import os
import sys
import time
import traceback

__all__ = ["ensure_app", "SettingsGuard", "ExceptionTrap", "ModalDismisser", "pump",
           "close_tool_windows"]

_SETTINGS = ("LARMOR", "app")


def ensure_app():
    """The one QApplication of the sweep (Fusion style, the default theme
    applied the way ``larmor.desktop.app.main`` does)."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("LARMOR_NO_SESSION", "1")
    os.environ.setdefault("LARMOR_NO_KERNEL_WARM", "1")
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication([sys.argv[0] if sys.argv else "larmor"])
        app.setStyle("Fusion")
    if not getattr(app, "_distcheck_themed", False):
        from larmor.desktop import theme
        theme.apply(app, theme.DEFAULT)
        app._distcheck_themed = True
    return app


class SettingsGuard:
    """Snapshot ``QSettings("LARMOR", "app")`` on entry and restore every key
    on exit, so a sweep that changes the theme, the recent files or any
    option leaves the user's registry exactly as it found it."""

    def __enter__(self):
        from PySide6.QtCore import QSettings
        s = QSettings(*_SETTINGS)
        self._saved = {k: s.value(k) for k in s.allKeys()}
        return self

    def __exit__(self, *exc):
        from PySide6.QtCore import QSettings
        s = QSettings(*_SETTINGS)
        s.clear()
        for k, v in self._saved.items():
            s.setValue(k, v)
        s.sync()
        return False


class ExceptionTrap:
    """Collect exceptions raised in Qt slots (``sys.excepthook``) and Qt's
    critical / fatal messages while active. ``errors`` is a list of strings."""

    _QT_LEVELS = None

    def __init__(self, warnings_too: bool = False):
        self.errors: list = []
        self.warnings: list = []
        self._warnings_too = warnings_too

    def __enter__(self):
        from PySide6.QtCore import QtMsgType, qInstallMessageHandler
        self._old_hook = sys.excepthook
        self._old_handler = qInstallMessageHandler(self._qt)
        self._crit = {QtMsgType.QtCriticalMsg, QtMsgType.QtFatalMsg}
        self._warn = QtMsgType.QtWarningMsg
        sys.excepthook = self._hook
        return self

    def __exit__(self, *exc):
        from PySide6.QtCore import qInstallMessageHandler
        sys.excepthook = self._old_hook
        qInstallMessageHandler(self._old_handler)
        return False

    def _hook(self, exc_type, exc, tb):
        self.errors.append("".join(traceback.format_exception(exc_type, exc, tb)))

    def _qt(self, mode, context, message):
        text = str(message)
        if mode in self._crit:
            self.errors.append("Qt critical: " + text)
        elif mode == self._warn:
            self.warnings.append(text)
            if self._warnings_too:
                self.errors.append("Qt warning: " + text)

    def take(self) -> list:
        """The errors so far, and reset."""
        out, self.errors = self.errors, []
        return out


class ModalDismisser:
    """While active, every modal widget (message box, file / input / colour
    dialog, any QDialog.exec) is rejected within ~100 ms, and its title is
    recorded, so a swept action that asks a question returns at once with
    the answer "cancel"."""

    def __init__(self, interval_ms: int = 100):
        from PySide6.QtCore import QTimer
        self.seen: list = []
        self._timer = QTimer()
        self._timer.setInterval(interval_ms)
        self._timer.timeout.connect(self._tick)

    def __enter__(self):
        self._timer.start()
        return self

    def __exit__(self, *exc):
        self._timer.stop()
        return False

    def _tick(self):
        from PySide6.QtWidgets import QApplication, QDialog, QMessageBox
        w = QApplication.activeModalWidget()
        if w is None:
            return
        self.seen.append(w.windowTitle() or type(w).__name__)
        try:
            if isinstance(w, QMessageBox):
                w.reject()
            elif isinstance(w, QDialog):
                w.reject()
            else:
                w.close()
        except Exception:                                 # noqa: BLE001
            pass


def pump(ms: int = 200) -> None:
    """Process events for ``ms`` milliseconds (without QTest)."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance()
    end = time.perf_counter() + ms / 1000.0
    while time.perf_counter() < end:
        app.processEvents()
        time.sleep(0.005)
    app.processEvents()


def close_tool_windows(keep=()) -> int:
    """Close every top-level widget except the ones in ``keep`` (the main
    window); returns how many were closed."""
    from PySide6.QtWidgets import QApplication
    n = 0
    keep = set(id(w) for w in keep)
    for w in list(QApplication.topLevelWidgets()):
        if id(w) in keep or not w.isVisible():
            continue
        try:
            w.close()
            n += 1
        except Exception:                                 # noqa: BLE001
            pass
    pump(50)
    return n
