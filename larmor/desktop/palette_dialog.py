"""Command palette: fuzzy search over every menu and toolbar action.

Help ▸ Command palette… (Ctrl+Shift+P) opens a modal dialog with a search
line above the full list of commands, each written as
"Menu ▸ Submenu ▸ Label  (hint)" with its shortcut in a second column.

The list is not a registry: collect_commands() walks the live QMenuBar and
the main QToolBar at open time, so every action the window holds -- the ones
built by MainWindow._add() and the ~50 built outside it (models, themes,
docks, recent files, recipes) -- is covered by construction, with its current
enabled/checked state. Matching is delegated to the Qt-free larmor.palette.

The chosen QAction is triggered AFTER the dialog's exec() returns, the same
way the sidebar drives menu actions (actResid.trigger()): checkables toggle,
exclusive groups switch, dock toggles show/hide, and any modal the action
opens never nests inside the palette's event loop.
"""
from __future__ import annotations

from dataclasses import dataclass

import shiboken6
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QColor, QKeySequence
from PySide6.QtWidgets import (
    QDialog, QHeaderView, QLabel, QLineEdit, QMainWindow, QMenu, QToolBar,
    QTreeWidget, QTreeWidgetItem, QVBoxLayout,
)

from larmor import palette
from larmor.desktop import theme

__all__ = ["Command", "collect_commands", "CommandPalette", "show_palette"]

PLACEHOLDER = "(no matching command)"


def _ptr(obj) -> int:
    """C++ address of a wrapped Qt object: a stable identity for actions and
    menus that does not depend on which Python wrapper is in hand."""
    return shiboken6.getCppPointer(obj)[0]


@dataclass
class Command:
    """One leaf QAction of the menu/toolbar tree, as the palette shows it."""
    action: QAction
    path: tuple[str, ...]      # menu titles, mnemonic- and hint-stripped
    label: str                 # leaf text, mnemonic-stripped, hint kept
    shortcut: str              # QKeySequence.NativeText, '' when none
    checkable: bool
    checked: bool
    enabled: bool
    tooltip: str               # the action's own tooltip, '' when default

    @property
    def display(self) -> str:
        return " ▸ ".join(self.path + (self.label,))

    @property
    def glyph(self) -> str:
        if self.checkable:
            return "☑ " if self.checked else "☐ "
        return ""

    @property
    def search_text(self) -> str:
        """What the matcher scores: path, label, hint and shortcut."""
        return (self.display + "  " + self.shortcut).rstrip()


def _default_tooltip(text: str) -> str:
    """QAction.toolTip() falls back to the text with mnemonics and ellipses
    removed; a tooltip equal to that carries no extra information."""
    return palette.strip_mnemonic(text).replace("...", "").replace("…", "").strip()


def _command(action: QAction, path: tuple[str, ...], seen: set[int]) -> Command | None:
    if action.isSeparator() or not action.isVisible():   # incl. addSection headers
        return None
    label = palette.strip_mnemonic(action.text())
    if not label:                                     # tb.addWidget(...) actions
        return None
    enabled = action.isEnabled()
    if not enabled and label.startswith("("):         # '(none yet)' placeholders
        return None
    if _ptr(action) in seen:                          # model actions: menu + toolbar
        return None
    seen.add(_ptr(action))
    tip = action.toolTip()
    if tip.strip() in (action.text(), label, _default_tooltip(action.text())):
        tip = ""
    return Command(
        action=action, path=path, label=label,
        shortcut=action.shortcut().toString(QKeySequence.NativeText),
        checkable=action.isCheckable(), checked=action.isChecked(),
        enabled=enabled, tooltip=tip)


def collect_commands(window: QMainWindow) -> list[Command]:
    """Every leaf action of the menu bar (recursively) and of the main
    toolbars, in menu order; the toolbar entries come last under 'Toolbar'
    so an action present in both keeps its menu path. Separators, section
    headers, hidden actions, widget actions and disabled '(…)' placeholders
    are skipped; the live enabled/checked state is read at call time.

    Submenus are resolved through a menuAction() -> QMenu map built from
    findChildren(QMenu), never through QAction.menu(): in PySide6 that call
    re-parents the window's own QMenu wrapper under the (temporary) action
    wrapper that asked, and when the temporary dies shiboken invalidates the
    menu wrapper and every C++-created action under it -- MainWindow.m_recent
    would be dead after one palette open. findChildren(), menuAction() and
    actions() leave wrapper ownership alone (measured on PySide6 6.11.1)."""
    seen: set[int] = set()
    out: list[Command] = []
    submenus = {_ptr(m.menuAction()): m for m in window.findChildren(QMenu)}

    def walk(actions, path: tuple[str, ...]) -> None:
        for a in actions:
            menu = submenus.get(_ptr(a))
            if menu is not None:
                if a.isVisible():
                    walk(menu.actions(), path + (palette.strip_hint(menu.title()),))
                continue
            cmd = _command(a, path, seen)
            if cmd is not None:
                out.append(cmd)

    walk(window.menuBar().actions(), ())
    for tb in window.findChildren(QToolBar):
        # direct children only (a toolbar inside a panel is not a command
        # surface); the sidebar holds short-labelled mirrors of View entries
        if tb.parent() is window and tb.objectName() != "sidebar":
            walk(tb.actions(), ("Toolbar",))
    return out


class CommandPalette(QDialog):
    """Search line + two-column list (command, shortcut) + a hidden note.

    Up/Down move one row, PgUp/PgDn ten, Enter or double-click chooses,
    Esc closes. Disabled commands are listed greyed and inert so the palette
    teaches that a command exists before data is open."""

    def __init__(self, parent, commands: list[Command]):
        super().__init__(parent)
        self.setWindowTitle("Command palette")
        self.resize(640, 420)
        self.commands = list(commands)
        self.chosen: QAction | None = None

        v = QVBoxLayout(self)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(
            "type part of a command — e.g. fit, czj, theme dark, F5")
        v.addWidget(self.edit)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderHidden(True)
        self.tree.setRootIsDecorated(False)
        self.tree.setUniformRowHeights(True)
        self.tree.setAllColumnsShowFocus(True)
        self.tree.setSelectionMode(QTreeWidget.SingleSelection)
        hdr = self.tree.header()
        hdr.setStretchLastSection(False)
        hdr.setSectionResizeMode(0, QHeaderView.Stretch)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        v.addWidget(self.tree, 1)

        self.note = QLabel()
        self.note.setWordWrap(True)
        self.note.setStyleSheet(f"color: {theme.active().text_dim};")
        self.note.hide()
        v.addWidget(self.note)

        self.edit.textChanged.connect(self.refill)
        self.edit.returnPressed.connect(self._choose)
        self.tree.itemActivated.connect(lambda *_: self._choose())
        self.refill("")
        self.edit.setFocus()

    # ------------------------------------------------------------ the list
    def refill(self, text: str) -> None:
        """Re-rank the commands for ``text`` and rebuild the rows."""
        self.note.hide()
        self.tree.clear()
        idx = palette.rank(text, [c.search_text for c in self.commands])
        grey = QColor(theme.active().disabled_text)
        for i in idx:
            c = self.commands[i]
            item = QTreeWidgetItem([c.glyph + c.display, c.shortcut])
            item.setData(0, Qt.UserRole, i)
            if c.tooltip:
                item.setToolTip(0, c.tooltip)
            if not c.enabled:
                item.setForeground(0, grey)
                item.setForeground(1, grey)
            self.tree.addTopLevelItem(item)
        if idx:
            self.tree.setCurrentItem(self.tree.topLevelItem(0))
        else:
            item = QTreeWidgetItem([PLACEHOLDER, ""])
            item.setData(0, Qt.UserRole, -1)
            item.setFlags(item.flags() & ~Qt.ItemIsSelectable & ~Qt.ItemIsEnabled)
            self.tree.addTopLevelItem(item)

    def current_command(self) -> Command | None:
        item = self.tree.currentItem()
        if item is None:
            return None
        i = item.data(0, Qt.UserRole)
        if i is None or int(i) < 0:
            return None
        return self.commands[int(i)]

    def _move(self, step: int) -> None:
        n = self.tree.topLevelItemCount()
        if not n:
            return
        cur = self.tree.currentItem()
        row = self.tree.indexOfTopLevelItem(cur) if cur is not None else 0
        row = max(0, min(n - 1, row + step))
        self.tree.setCurrentItem(self.tree.topLevelItem(row))

    def keyPressEvent(self, ev):
        # the line edit ignores these keys, so they arrive here while typing
        key = ev.key()
        if key == Qt.Key_Down:
            self._move(1)
        elif key == Qt.Key_Up:
            self._move(-1)
        elif key == Qt.Key_PageDown:
            self._move(10)
        elif key == Qt.Key_PageUp:
            self._move(-10)
        else:
            super().keyPressEvent(ev)          # Esc rejects as usual
            return
        ev.accept()

    # ---------------------------------------------------------- the choice
    def _choose(self) -> None:
        cmd = self.current_command()
        if cmd is None:
            return
        if not cmd.enabled:
            self.note.setText(
                f"'{cmd.label}' is not available right now — its menu entry "
                "is disabled (usually: open a spectrum first)")
            self.note.show()
            return
        self.chosen = cmd.action
        self.accept()


def show_palette(window) -> QAction | None:
    """Open the palette over ``window``; run the chosen action after the
    dialog has closed and return it (None when nothing was chosen).

    The status line 'Menu ▸ Label' is written BEFORE trigger() so an action
    that posts its own message (add-mode, scroll-nudge, Czjzek display)
    replaces it, while silent actions keep the confirmation."""
    cmds = collect_commands(window)
    dlg = CommandPalette(window, cmds)
    try:
        if dlg.exec() != QDialog.Accepted or dlg.chosen is None:
            return None
        act = dlg.chosen
    finally:
        dlg.deleteLater()
    display = palette.strip_mnemonic(act.text())
    for c in cmds:
        if c.action is act:
            display = c.display
            break
    window.statusBar().showMessage(display)
    act.trigger()
    return act
