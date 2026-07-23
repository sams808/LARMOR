"""The Lineshapes reference manual must ship, cover every registered model, and
be reachable from the ? menu."""
import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_lineshapes_manual_covers_every_model():
    from larmor import models
    from larmor.desktop.help_dialog import help_path

    p = help_path("lineshapes")
    assert p is not None, "larmor/help/lineshapes.md is missing"
    raw = p.read_text(encoding="utf-8")
    flat = " ".join(raw.split())                      # collapse line wraps
    for m in models.describe_all():
        assert f"`{m['name']}`" in raw, f"model {m['name']} not in lineshapes.md"
    for token in ("Czjzek", "d'Espinose", "Central Limit Theorem",
                  "quadrupolar product", "Haeberlen", "mrsimulator",
                  "Sam Soudani"):
        assert token in flat, f"missing: {token}"


def test_lineshapes_item_in_help_menu():
    from PySide6.QtWidgets import QApplication, QMenu

    os.environ["LARMOR_NO_SESSION"] = "1"
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    help_menu = next(m for m in win.menuBar().findChildren(QMenu)
                     if m.title() == "&?")
    labels = [a.text() for a in help_menu.actions() if a.text()]
    assert any("Lineshapes" in x for x in labels)
    win.close()
