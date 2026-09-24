"""The "open windows" bar of the left sidebar (larmor/desktop/windowtray.py).

Sam's field-use feedback: a minimised popup window used to shrink to a
floating title bar at the bottom-left of the screen, and there was no way to
see which tool windows were already open in a session. The bar is an
application-level event filter feeding a strip at the very bottom of the
sidebar; these tests drive it offscreen the way a window manager would
(``setWindowState(Qt.WindowMinimized)`` sends the same QWindowStateChangeEvent)
and through the public helpers ``tray_minimize`` / ``tray_restore``.
"""
import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEvent, Qt  # noqa: E402
from PySide6.QtWidgets import (QApplication, QDialog, QLabel, QMenu,  # noqa: E402
                               QToolBar, QVBoxLayout, QWidget)


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    w.show()
    qapp.processEvents()
    yield w
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    if bar is not None:
        bar.close_all()
    w.close()


def _tool(win, title):
    from larmor.desktop.windowtray import show_tool_window

    d = QDialog(win)
    d.setWindowTitle(title)
    QVBoxLayout(d).addWidget(QLabel(title))
    show_tool_window(d)
    QApplication.processEvents()
    return d


def _flush(qapp):
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    qapp.processEvents()


def test_bar_sits_at_the_bottom_of_the_sidebar_and_starts_hidden(win):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    assert bar is not None and bar is win._window_bar
    sidebar = next(tb for tb in win.findChildren(QToolBar) if tb.objectName() == "sidebar")
    assert bar.widget.parent() is sidebar
    assert sidebar.actions()[-1] is bar.strip_action         # after Scroll
    assert sidebar.widgetForAction(bar.strip_action) is bar.widget
    assert not bar.strip_action.text()                       # not a palette command
    assert bar.count() == 0 and not bar.is_shown()


def test_bar_lists_open_tool_windows_oldest_first_with_a_count(win, qapp):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    d1 = _tool(win, "QCPMG processing — 2702")
    d2 = _tool(win, "QCPMG processing — 2702")      # a duplicate stays visible
    assert bar.entries() == [d1, d2] == windowtray.open_windows()
    assert bar.count() == 2 and bar.is_shown()
    assert bar.widget.header.text() == "2 windows open"
    btns = bar.widget.buttons()
    assert len(btns) == 2
    assert all(b.property("trayTitle") == "QCPMG processing — 2702" for b in btns)
    assert all("QCPMG processing — 2702" in b.toolTip() for b in btns)
    assert not d1.isModal() and not d2.isModal()
    # the main window keeps taking menu actions while they are open
    before = win.actResid.isChecked()
    win.actResid.trigger()
    assert win.actResid.isChecked() != before
    # a title change follows into the entry
    d2.setWindowTitle("Referencing audit")
    qapp.processEvents()
    assert bar.button_for(d2).property("trayTitle") == "Referencing audit"
    d1.close()
    _flush(qapp)
    assert bar.entries() == [d2] and bar.widget.header.text() == "1 window open"


def test_minimising_hides_the_window_and_dims_its_entry(win, qapp):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    d1 = _tool(win, "Session inventory")
    d2 = _tool(win, "Referencing audit")
    d1.setWindowState(Qt.WindowMinimized)                    # what the WM does
    qapp.processEvents()
    assert d1.isHidden()                                     # no floating title bar
    assert bar.is_minimized(d1) and not bar.is_minimized(d2)
    assert bar.count() == 2                                  # still listed
    b1 = bar.button_for(d1)
    assert b1.property("trayMinimized") is True and "italic" in b1.styleSheet()
    assert "minimised" in b1.toolTip()
    b1.click()                                               # restore
    qapp.processEvents()
    assert d1.isVisible() and not (d1.windowState() & Qt.WindowMinimized)
    assert not bar.is_minimized(d1)
    assert bar.button_for(d1).property("trayMinimized") is False
    assert bar.button_for(d1).styleSheet() == ""


def test_public_helpers_tray_minimize_and_tray_restore(win, qapp):
    from larmor.desktop import windowtray

    d = _tool(win, "Herzfeld–Berger")
    assert windowtray.tray_minimize(d) is True
    assert d.isHidden() and windowtray.current_tray().is_minimized(d)
    assert windowtray.tray_restore(d) is True
    assert d.isVisible() and not windowtray.current_tray().is_minimized(d)
    assert windowtray.tray_restore(d) is False               # not minimised now


def test_closing_a_minimised_window_drops_its_entry(win, qapp):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    d = _tool(win, "Experimental section")
    windowtray.tray_minimize(d)
    assert bar.count() == 1
    d.close()
    _flush(qapp)
    assert bar.count() == 0 and not bar.is_shown()


def test_context_menu_restore_minimise_close_and_close_all(win, qapp):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    d1 = _tool(win, "Comparability")
    d2 = _tool(win, "Series evolution")
    m = bar.context_menu(d1)
    labels = [a.text() for a in m.actions() if not a.isSeparator()]
    assert labels == ["Restore", "Minimise", "Close", "Close all"]
    next(a for a in m.actions() if a.text() == "Minimise").trigger()
    qapp.processEvents()
    assert d1.isHidden() and bar.is_minimized(d1)
    m = bar.context_menu(d1)
    assert not next(a for a in m.actions() if a.text() == "Minimise").isEnabled()
    next(a for a in m.actions() if a.text() == "Restore").trigger()
    qapp.processEvents()
    assert d1.isVisible() and not bar.is_minimized(d1)
    next(a for a in bar.context_menu(d1).actions() if a.text() == "Close").trigger()
    _flush(qapp)
    assert bar.entries() == [d2]
    next(a for a in bar.context_menu(d2).actions() if a.text() == "Close all").trigger()
    _flush(qapp)
    assert bar.count() == 0 and not bar.is_shown()


def test_hidden_by_code_leaves_the_bar_and_reshowing_relists_last(win, qapp):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    d1 = _tool(win, "first")
    d2 = _tool(win, "second")
    d1.hide()                                                # a kept-alive dialog closing
    qapp.processEvents()
    assert bar.entries() == [d2]
    d1.show()
    qapp.processEvents()
    assert bar.entries() == [d2, d1]                         # a new opening: last


def test_modal_dialogs_are_never_listed_and_a_minimised_one_is_restored(win, qapp):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    md = QDialog(win)
    md.setWindowTitle("Save fit as")
    md.setModal(True)
    md.show()
    qapp.processEvents()
    assert md not in bar.entries() and bar.count() == 0
    md.setWindowState(Qt.WindowMinimized)
    qapp.processEvents()
    qapp.processEvents()                                     # the 0-ms restore
    assert not (md.windowState() & Qt.WindowMinimized) and md.isVisible()
    assert bar.count() == 0
    md.close()


def test_filter_ignores_popups_tool_windows_and_child_widgets(win):
    from larmor.desktop import windowtray

    bar = windowtray.current_tray()
    assert not bar._is_candidate(QMenu(win))                 # Qt.Popup
    assert not bar._is_candidate(QWidget(win, Qt.Tool))      # floating dock style
    assert not bar._is_candidate(QLabel("x", win))           # not a window
    assert not bar._is_candidate(win)                        # the main window
    d = QDialog(win)
    assert bar._is_candidate(d)


def test_one_bar_is_active_across_main_windows(win, qapp, monkeypatch):
    """A second MainWindow (tests build several per process) takes over the
    application filter; the first one's bar stops collecting."""
    from larmor.desktop import windowtray
    from larmor.desktop.app import MainWindow

    first = windowtray.current_tray()
    other = MainWindow()
    try:
        assert windowtray.current_tray() is other._window_bar is not first
        d = _tool(other, "on the second window")
        assert d in other._window_bar.entries()
        assert d not in first.entries()
    finally:
        other._window_bar.close_all()
        other.close()
