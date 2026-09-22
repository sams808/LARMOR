"""Offscreen tests of the command palette's Qt half (F8).

collect_commands / CommandPalette / show_palette are exercised on a small
synthetic QMainWindow (menus, a main toolbar, a sidebar toolbar), then one
MainWindow test pins the Help entry, the coverage of the live menu tree, the
run-through of a checkable, and -- as a permanent guard the codebase lacked --
the uniqueness of every shortcut across every menu and toolbar action.
"""
import gc
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
import shiboken6  # noqa: E402
from PySide6.QtCore import QEvent, Qt  # noqa: E402
from PySide6.QtGui import QAction, QColor, QKeyEvent, QKeySequence  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import (  # noqa: E402
    QApplication, QDialog, QLabel, QMainWindow, QMenu, QToolBar,
)

from larmor.desktop import theme  # noqa: E402
from larmor.desktop.palette_dialog import (  # noqa: E402
    PLACEHOLDER, CommandPalette, collect_commands, show_palette,
)


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _menu_window():
    """'&File' with a shortcut action, separator, section header, a hinted
    submenu holding one real entry and a disabled placeholder, a hidden
    action, a disabled hint and a checked checkable with a custom tooltip."""
    win = QMainWindow()
    acts = {}
    m_file = win.menuBar().addMenu("&File")
    acts["open"] = m_file.addAction("Open…")
    acts["open"].setShortcut(QKeySequence("Ctrl+O"))
    m_file.addSeparator()
    m_file.addSection("Sec")
    m_recent = m_file.addMenu("Open &recent  (last 12)")
    acts["recent"] = m_recent.addAction("x.fxmla   —   C:/d")
    none = m_recent.addAction("(none yet)")
    none.setEnabled(False)
    hidden = m_file.addAction("Hidden one")
    hidden.setVisible(False)
    acts["cli"] = m_file.addAction("CLI hint")
    acts["cli"].setEnabled(False)
    acts["resid"] = m_file.addAction("Residual")
    acts["resid"].setCheckable(True)
    acts["resid"].setChecked(True)
    acts["resid"].setToolTip("show the residual")
    return win, acts


def _full_window():
    """The menu window plus a main toolbar (a menu duplicate, a toolbar-only
    Undo, a QLabel widget action) and a sidebar toolbar."""
    win, acts = _menu_window()
    tb = QToolBar("main")
    win.addToolBar(tb)
    tb.addAction(acts["open"])
    acts["undo"] = QAction("Undo", win)
    acts["undo"].setShortcut(QKeySequence("Ctrl+Z"))
    tb.addAction(acts["undo"])
    tb.addWidget(QLabel("  Add line "))
    sb = QToolBar("view")
    sb.setObjectName("sidebar")
    win.addToolBar(Qt.LeftToolBarArea, sb)
    sb.addAction("Full")
    return win, acts


def _row(dlg):
    return dlg.tree.indexOfTopLevelItem(dlg.tree.currentItem())


def _key(dlg, key):
    dlg.keyPressEvent(QKeyEvent(QEvent.KeyPress, key, Qt.NoModifier))


# ------------------------------------------------------------ collecting
def test_collect_commands_walks_submenus_and_skips_noise(qapp):
    win, acts = _menu_window()
    cmds = collect_commands(win)
    assert [(c.path, c.label) for c in cmds] == [
        (("File",), "Open…"),
        (("File", "Open recent"), "x.fxmla   —   C:/d"),   # hint + mnemonic stripped
        (("File",), "CLI hint"),
        (("File",), "Residual"),
    ]
    open_ = cmds[0]
    assert open_.action is acts["open"]
    assert open_.shortcut == "Ctrl+O"
    assert open_.checkable is False and open_.enabled is True
    assert open_.display == "File ▸ Open…"
    assert open_.glyph == ""
    assert open_.tooltip == ""                    # Qt's default tooltip (== text)
    assert "Ctrl+O" in open_.search_text
    assert cmds[2].enabled is False
    resid = cmds[3]
    assert resid.checkable and resid.checked
    assert resid.glyph == "☑ "
    assert resid.tooltip == "show the residual"
    labels = [c.label for c in cmds]
    for absent in ("Sec", "Hidden one", "(none yet)"):
        assert absent not in labels


def test_collect_dedupes_toolbar_and_skips_sidebar_and_widget_actions(qapp):
    win, acts = _full_window()
    cmds = collect_commands(win)
    opens = [c for c in cmds if c.label == "Open…"]
    assert len(opens) == 1 and opens[0].path == ("File",)   # menu path wins
    undo = [c for c in cmds if c.label == "Undo"]
    assert len(undo) == 1
    assert undo[0].path == ("Toolbar",) and undo[0].shortcut == "Ctrl+Z"
    assert all(c.label for c in cmds)                       # no QLabel row
    assert not any(c.label == "Full" for c in cmds)         # sidebar skipped
    assert len(cmds) == 5


def test_collect_leaves_wrapper_ownership_alone(qapp):
    """Regression: a walk through QAction.menu() re-parented the window's
    QMenu wrappers under temporary action wrappers; once those died, shiboken
    invalidated the menus and every C++-created action under them (rows
    reported 'already deleted', the toolbar duplicate was no longer
    recognised, and the app's own submenu handles were dead)."""
    win, acts = _menu_window()
    recent_menu = next(m for m in win.menuBar().findChildren(QMenu)
                       if m.title().startswith("Open &recent"))
    cmds = collect_commands(win)
    cmds2 = collect_commands(win)                  # a second open, same window
    gc.collect()
    assert all(shiboken6.isValid(c.action) for c in cmds + cmds2)
    assert shiboken6.isValid(recent_menu)
    assert all(shiboken6.isValid(a) for a in acts.values())
    assert [c.label for c in cmds] == [c.label for c in cmds2]
    # the C++-created recent entry can still be driven
    fired = []
    acts["recent"].triggered.connect(lambda *_: fired.append(1))
    next(c for c in cmds if c.label.startswith("x.fxmla")).action.trigger()
    assert fired == [1]
    del cmds, cmds2
    gc.collect()
    assert shiboken6.isValid(recent_menu) and shiboken6.isValid(acts["recent"])


# ------------------------------------------------------------- the dialog
def test_dialog_filters_and_chooses_on_enter(qapp):
    win, acts = _full_window()
    cmds = collect_commands(win)
    dlg = CommandPalette(None, cmds)
    assert dlg.tree.topLevelItemCount() == len(cmds) == 5
    assert _row(dlg) == 0
    assert dlg.tree.topLevelItem(0).text(0) == "File ▸ Open…"
    assert dlg.tree.topLevelItem(4).text(0) == "Toolbar ▸ Undo"
    dlg.edit.setText("undo")
    assert dlg.tree.topLevelItemCount() == 1
    item = dlg.tree.topLevelItem(0)
    assert "Undo" in item.text(0) and item.text(1) == "Ctrl+Z"
    dlg.edit.returnPressed.emit()
    assert dlg.chosen is acts["undo"]
    assert dlg.result() == QDialog.Accepted


def test_checkable_rows_carry_state_glyph_and_tooltip(qapp):
    win, acts = _full_window()
    dlg = CommandPalette(None, collect_commands(win))
    dlg.edit.setText("resid")
    item = dlg.tree.topLevelItem(0)
    assert item.text(0) == "☑ File ▸ Residual"
    assert item.toolTip(0) == "show the residual"
    acts["resid"].setChecked(False)
    dlg2 = CommandPalette(None, collect_commands(win))   # re-collected per open
    dlg2.edit.setText("resid")
    assert dlg2.tree.topLevelItem(0).text(0) == "☐ File ▸ Residual"


def test_disabled_row_is_grey_and_inert(qapp):
    win, acts = _full_window()
    dlg = CommandPalette(None, collect_commands(win))
    dlg.edit.setText("cli")
    assert dlg.tree.topLevelItemCount() == 1
    item = dlg.tree.topLevelItem(0)
    assert item.foreground(0).color() == QColor(theme.active().disabled_text)
    assert item.foreground(1).color() == QColor(theme.active().disabled_text)
    dlg.edit.returnPressed.emit()
    assert dlg.chosen is None
    assert dlg.result() != QDialog.Accepted
    assert not dlg.note.isHidden()
    assert "not available" in dlg.note.text() and "CLI hint" in dlg.note.text()
    dlg.edit.setText("undo")                       # typing again clears the note
    assert dlg.note.isHidden()


def test_no_match_placeholder_is_inert(qapp):
    win, _ = _full_window()
    dlg = CommandPalette(None, collect_commands(win))
    dlg.edit.setText("zzqqxx")
    assert dlg.tree.topLevelItemCount() == 1
    item = dlg.tree.topLevelItem(0)
    assert item.text(0) == PLACEHOLDER
    assert not (item.flags() & Qt.ItemIsEnabled)
    assert not (item.flags() & Qt.ItemIsSelectable)
    assert dlg.current_command() is None
    dlg.edit.returnPressed.emit()
    assert dlg.chosen is None
    assert dlg.result() != QDialog.Accepted
    dlg.edit.setText("")                           # back to the full map
    assert dlg.tree.topLevelItemCount() == 5 and _row(dlg) == 0


def test_keyboard_navigation(qapp):
    win, _ = _full_window()
    dlg = CommandPalette(None, collect_commands(win))
    n = dlg.tree.topLevelItemCount()
    assert n == 5 and _row(dlg) == 0
    _key(dlg, Qt.Key_Down)
    _key(dlg, Qt.Key_Down)
    assert _row(dlg) == 2
    for _ in range(n):
        _key(dlg, Qt.Key_Down)
    assert _row(dlg) == n - 1                      # clamped at the last row
    _key(dlg, Qt.Key_Up)
    assert _row(dlg) == n - 2
    _key(dlg, Qt.Key_PageUp)
    assert _row(dlg) == 0
    _key(dlg, Qt.Key_PageDown)
    assert _row(dlg) == n - 1
    _key(dlg, Qt.Key_PageUp)
    # the real path: keys typed into the focused line edit reach the dialog
    dlg.show()
    dlg.edit.setFocus()
    QTest.keyClick(dlg.edit, Qt.Key_Down)
    assert _row(dlg) == 1
    QTest.keyClick(dlg.edit, Qt.Key_Up)
    assert _row(dlg) == 0
    dlg.close()


# ---------------------------------------------------------- the launcher
def test_show_palette_triggers_after_close_and_toggles_checkable(qapp, monkeypatch):
    win, acts = _full_window()
    resid = acts["resid"]
    fired = []
    resid.triggered.connect(lambda *_: fired.append(1))

    def fake_exec(self):
        self.chosen = resid
        return QDialog.Accepted

    monkeypatch.setattr(CommandPalette, "exec", fake_exec)
    assert resid.isChecked() is True
    assert show_palette(win) is resid
    assert fired == [1]
    assert resid.isChecked() is False               # trigger() toggled it
    assert win.statusBar().currentMessage() == "File ▸ Residual"

    monkeypatch.setattr(CommandPalette, "exec", lambda self: QDialog.Rejected)
    assert show_palette(win) is None
    assert fired == [1]


# ------------------------------------------------------- the real window
def test_mainwindow_palette_entry_coverage_and_shortcut_uniqueness(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    try:
        help_menu = next(m for m in win.menuBar().findChildren(QMenu)
                         if m.title() == "&Help")
        first = next(a for a in help_menu.actions() if not a.isSeparator())
        assert first is win.actPalette
        assert first.text().startswith("&Command palette")
        assert first.shortcut().toString() == "Ctrl+Shift+P"
        labels = [a.text() for a in help_menu.actions() if a.text()]
        assert labels.index("More…") == labels.index("About LARMOR") + 1

        cmds = collect_commands(win)
        assert len(cmds) >= 110, len(cmds)
        assert all(c.path and c.label for c in cmds)
        by = {(c.path, c.label): c for c in cmds}
        fit = by[(("Decomposition",), "Fit")]
        assert fit.shortcut == "F5" and fit.enabled is False    # no data yet
        undo = by[(("Toolbar",), "↩  Undo")]
        assert undo.shortcut == "Ctrl+Z"
        assert by[(("View", "Panels"), "Explorer")].checkable
        assert (("Help", "User manuals"), "QCPMG") in by
        assert (("Decomposition", "Apply recipe"), "Browse for recipe…") in by
        assert by[(("View",), "Residual")].checked is True
        assert by[(("Help",), "Command palette…  (find any menu entry)"
                   )].shortcut == "Ctrl+Shift+P"
        gauss = [c for c in cmds if c.label == "Gauss/Lorentz"]
        assert len(gauss) == 1                                  # toolbar dedupe
        assert gauss[0].path == ("Decomposition", "Add line")
        assert not any(c.label.startswith("(") for c in cmds)
        shortcuts = [c.shortcut for c in cmds if c.shortcut]
        dupes = sorted({s for s in shortcuts if shortcuts.count(s) > 1})
        assert not dupes, f"shortcut shadowed across menus/toolbars: {dupes}"

        # the walk must leave the window's own wrappers alive: the recent-
        # files and apply-recipe submenus are rebuilt later via .clear()
        gc.collect()
        assert all(shiboken6.isValid(c.action) for c in cmds)
        for handle in (win.m_recent, win.m_apply, win.m_view, win.actFit,
                       win.actUndo, win.actPalette):
            assert shiboken6.isValid(handle)
        win._rebuild_recent()                                  # must not raise
        win._rebuild_apply_recipe()

        # run a checkable through the palette: toggles, status line written
        def fake_exec(self):
            self.chosen = win.actResid
            return QDialog.Accepted

        monkeypatch.setattr(CommandPalette, "exec", fake_exec)
        assert win.actResid.isChecked() is True
        win.open_command_palette()
        assert win.actResid.isChecked() is False
        assert win.statusBar().currentMessage() == "View ▸ Residual"
    finally:
        win.close()
