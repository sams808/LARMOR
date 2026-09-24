"""The "＋ Add line" tool on the main toolbar (larmor/desktop/addline_toolbar.py).

A split button showing the CURRENT model, a grouped menu of every model with
descriptive tooltips, quick buttons for the most used models and a "placing …"
label while a placing mode is on -- all arranging the SAME checkable QActions
``MainWindow._model_actions`` holds, so _set_add_mode, the command palette and
the existing tests keep working.
"""
import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtCore import QSettings  # noqa: E402
from PySide6.QtGui import QAction  # noqa: E402
from PySide6.QtWidgets import QApplication, QToolBar  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _main_toolbar(win) -> QToolBar:
    return next(tb for tb in win.findChildren(QToolBar)
                if tb.parent() is win and tb.objectName() != "sidebar")


def test_groups_cover_every_model_and_tooltips_describe_them():
    from larmor import models
    from larmor.desktop import addline_toolbar as at

    grouped = [n for _, names in at.GROUPS for n in names]
    assert len(grouped) == len(set(grouped))
    reg = set(models.REGISTRY) - {"spectrum"}            # spectrum is a dialog
    assert set(grouped) == reg, set(grouped) ^ reg
    assert set(at.QUICK) <= reg and len(at.QUICK) == 5
    for n in grouped:
        assert n in models.REGISTRY, n
    tip = at.model_tooltip("czjzek")
    assert tip.startswith("Czjzek (quad. distribution) — Czjzek distribution")
    assert "Key parameters: position, σ(Cq), width (dCS), amplitude" in tip
    assert "(lb held unless freed)" in tip
    assert "Gauss/Lorentz mix held unless freed" in at.model_tooltip("gauss_lor")
    assert at.model_tooltip("no_such_model") == "no_such_model"


def test_toolbar_holds_the_split_button_quick_buttons_and_the_same_actions(qapp, win):
    from larmor.desktop import addline_toolbar as at

    tbw = win._addline_toolbar
    assert isinstance(tbw, at.AddLineToolbar)
    tb = _main_toolbar(win)
    acts = [a for a in tb.actions() if not a.isSeparator()]
    # Undo, Redo, the split button, five quick buttons, the placing label
    assert [a.text() for a in acts[:2]] == ["↩  Undo", "↪  Redo"]
    quick = [a for a in acts if a in win._model_actions.values()]
    assert [next(n for n, a in win._model_actions.items() if a is q) for q in quick] \
        == list(at.QUICK)
    assert tb.widgetForAction(acts[2]) is tbw.button
    # every model action sits in the split menu exactly once, under its group
    menu_acts = [a for a in tbw.menu.actions() if not a.isSeparator()]
    model_acts = [a for a in menu_acts if a in win._model_actions.values()]
    assert len(model_acts) == len(win._model_actions) == len(set(map(id, model_acts)))
    sections = [a.text() for a in tbw.menu.actions() if a.isSeparator() and a.text()]
    assert sections == [t for t, _ in at.GROUPS]
    assert tbw.action_background.text() == "Spectrum (background)…"
    assert all(a.toolTip() for a in model_acts)
    # win._model_actions is still the name -> checkable QAction dict
    assert all(isinstance(a, QAction) and a.isCheckable()
               for a in win._model_actions.values())
    assert set(win._model_actions) == {n for _, names in at.GROUPS for n in names}


def test_split_button_shows_the_current_model_and_the_placing_state(qapp, win):
    tbw = win._addline_toolbar
    QSettings("LARMOR", "app").remove("addLine/model")
    assert tbw.placing_model() is None and not tbw.placing.isVisibleTo(win)
    assert tbw.button.text().startswith("＋ Add line: ")
    # a model action checked anywhere (menu, quick button, Decomposition menu)
    # becomes the current model and turns the placing label on
    win._model_actions["czjzek"].setChecked(True)
    assert tbw.placing_model() == "czjzek" and tbw.current_model() == "czjzek"
    assert tbw.button.text() == "＋ Add line: Czjzek (quad. distribution)"
    assert tbw.button.isChecked()
    assert tbw.placing.isVisibleTo(win)
    assert tbw.placing.text() == ("placing Czjzek (quad. distribution) — click on "
                                  "the spectrum · Esc to stop")
    assert "Esc to stop" in tbw.button.toolTip()
    # Esc / _set_add_mode(None) leaves the model on the button, label off
    win._set_add_mode(None)
    assert tbw.placing_model() is None and not tbw.placing.isVisibleTo(win)
    assert tbw.current_model() == "czjzek" and not tbw.button.isChecked()
    assert tbw.button.text() == "＋ Add line: Czjzek (quad. distribution)"
    # the main part of the split button toggles that model's placing mode
    tbw.button.click()
    assert win._model_actions["czjzek"].isChecked() and tbw.placing_model() == "czjzek"
    assert win.view._add_mode == "czjzek"
    tbw.button.click()
    assert tbw.placing_model() is None
    assert not any(a.isChecked() for a in win._model_actions.values())
    # the choice is remembered (globally, and per nucleus once one is known)
    assert QSettings("LARMOR", "app").value("addLine/model") == "czjzek"
    from larmor.desktop.addline_toolbar import AddLineToolbar

    fresh = AddLineToolbar(win, QToolBar())
    assert fresh.current_model() == "czjzek"
    QSettings("LARMOR", "app").remove("addLine/model")


def test_palette_lists_every_model_once_through_the_edit_menu(qapp, win):
    """The split button and Edit > Add line share ONE QAction per model, so
    the palette lists each model exactly once, under Edit > Add line."""
    from larmor.desktop.palette_dialog import collect_commands

    cmds = collect_commands(win)
    for name, act in win._model_actions.items():
        hits = [c for c in cmds if c.action is act]
        assert len(hits) == 1, (name, len(hits))
        assert tuple(hits[0].path[:2]) == ("Edit", "Add line"), hits[0].path


def test_theme_switch_restyles_the_button(qapp, win):
    from larmor.desktop import theme

    tbw = win._addline_toolbar
    prev = theme.active().name
    try:
        other = next(n for n in theme.names() if n != prev)
        win._apply_theme_live(other)
        assert theme.active().accent in tbw.button.styleSheet()
    finally:
        win._apply_theme_live(prev)
