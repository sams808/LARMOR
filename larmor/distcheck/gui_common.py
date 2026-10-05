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
    option leaves the user's registry exactly as it found it.

    The snapshot also goes to DISK first (``backup_path``: a ``reg export``
    of ``HKCU\\Software\\LARMOR`` on Windows, a JSON dump elsewhere) under
    ``%LOCALAPPDATA%\\LARMOR\\settings_backup`` (the last 20 kept), because a
    guard that lives only in memory is lost with its process: a test that
    hung in a modal dialog and had to be killed took the user's saved
    recent files, pinned folders and libraries with it (2026-10-05)."""

    backup_path = None

    def __enter__(self):
        from larmor.desktop.prefs import settings
        self.backup_path = backup_settings()
        s = settings()
        self._saved = {k: s.value(k) for k in s.allKeys()}
        return self

    def __exit__(self, *exc):
        from larmor.desktop.prefs import settings
        s = settings()
        s.clear()
        for k, v in self._saved.items():
            s.setValue(k, v)
        s.sync()
        return False


def backup_settings():
    """Write a restorable copy of the user's LARMOR settings and return its
    path (None when nothing could be written). Windows: ``reg export`` of
    ``HKCU\\Software\\LARMOR`` (``reg import <file>`` puts it back after
    ``reg delete HKCU\\Software\\LARMOR /f``); elsewhere a JSON of the keys."""
    import datetime
    import json
    import shutil
    import subprocess
    from pathlib import Path

    from larmor.desktop.prefs import SETTINGS_FILE_ENV, settings

    base = Path(os.environ.get("LOCALAPPDATA") or Path.home()) / "LARMOR" / "settings_backup"
    try:
        base.mkdir(parents=True, exist_ok=True)
    except OSError:
        return None
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + f"_{os.getpid()}"
    try:
        ini = os.environ.get(SETTINGS_FILE_ENV)
        if ini:
            # an .ini store (the test suite, the distribution check): a copy
            src = Path(ini)
            if not src.exists():
                return None
            target = base / f"LARMOR_settings_{stamp}.ini"
            shutil.copyfile(src, target)
        elif sys.platform == "win32":
            target = base / f"LARMOR_settings_{stamp}.reg"
            r = subprocess.run(["reg", "export", r"HKCU\Software\LARMOR", str(target), "/y"],
                               capture_output=True, text=True, timeout=30)
            if r.returncode != 0:            # no key yet: nothing to back up
                target = None
        else:
            s = settings()
            target = base / f"LARMOR_settings_{stamp}.json"
            target.write_text(json.dumps({k: s.value(k) for k in s.allKeys()},
                                         default=str, indent=1), encoding="utf-8")
        old = sorted(base.glob("LARMOR_settings_*"), key=lambda p: p.stat().st_mtime)
        for p in old[:-20]:
            try:
                p.unlink()
            except OSError:
                pass
        return target
    except Exception:                                     # noqa: BLE001
        return None


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


def running_threads() -> list:
    """Every Python-owned QThread still running (FitWorker, SimWorker, a
    dialog's worker...), the main thread excluded."""
    import gc

    import shiboken6
    from PySide6.QtCore import QThread

    main = shiboken6.getCppPointer(QThread.currentThread())[0]
    out = []
    for obj in gc.get_objects():
        if isinstance(obj, QThread):
            try:
                if (shiboken6.isValid(obj) and shiboken6.getCppPointer(obj)[0] != main
                        and obj.isRunning()):
                    out.append(obj)
            except RuntimeError:                          # deleted meanwhile
                pass
    return out


def join_threads(timeout: float = 30.0) -> list:
    """Wait, pumping events, for every running QThread; the class names of
    those still running after ``timeout`` s each."""
    left = []
    for t in running_threads():
        end = time.perf_counter() + timeout
        try:
            while t.isRunning() and time.perf_counter() < end:
                pump(20)
            if t.isRunning():
                left.append(type(t).__name__)
        except RuntimeError:
            pass
    return left


def destroy_all_windows():
    """Wait for every running worker thread, then destroy the leftover
    windows (below). Returns ``(destroyed, workers_waited_for,
    workers_still_running)``: nothing is destroyed while a worker runs --
    Qt aborts the interpreter when a running QThread's owner is deleted
    (which is how the first version of this clean-up crashed after the
    error-tools stage)."""
    waited = [type(t).__name__ for t in running_threads()]
    still = join_threads(30.0) if waited else []
    if still:
        return 0, waited, still
    return _destroy_windows(), waited, []


def _destroy_windows() -> int:
    """Between GUI stages: perform the deletions the stage queued, then
    destroy the parentless windows that were never shown; returns how many
    were destroyed.

    ``deleteLater`` is honoured only when control returns to an event loop,
    and ``processEvents()`` is not that -- so every window a stage closed
    stayed alive, and every theme switch of a later stage restyled all of
    them: the menu sweep took 6 minutes on its own and more than 40 after
    the other stages in one process. ``sendPostedEvents(DeferredDelete)``
    performs the pending deletions.

    What is NOT destroyed: a widget with a parent (it goes with its owner),
    Qt's popups (menus, combo-box drop-downs, tool tips: their owners keep
    pointers to them -- destroying every top-level widget crashed the
    interpreter), and a window that has been SHOWN (it owns a native window
    handle). Deleting the main window a dialog stage had shown made Qt warn
    "shared QObject was deleted directly" and abort the interpreter; the
    same deletion outside the check -- a shown main window closed and
    deleted mid-session, after a fit, themed or not, and LARMOR's own exit
    through app.main() -- went through cleanly every time, so the cause is
    in the harness, and the few shown windows are simply left closed."""
    import gc

    from PySide6.QtCore import QCoreApplication, QEvent, Qt
    from PySide6.QtWidgets import QApplication, QMenu

    app = QApplication.instance()
    if app is None:
        return 0
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    skip_types = {Qt.WindowType.Popup, Qt.WindowType.ToolTip, Qt.WindowType.Desktop,
                  Qt.WindowType.Drawer, Qt.WindowType.Sheet}
    n = 0
    for w in list(QApplication.topLevelWidgets()):
        try:
            if w.parent() is not None or isinstance(w, QMenu):
                continue
            if w.windowType() in skip_types or w.windowHandle() is not None:
                continue
            w.close()
            w.deleteLater()
            n += 1
        except RuntimeError:                              # already gone on the C++ side
            pass
    for _ in range(3):
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        app.processEvents()
    gc.collect()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    return n


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
