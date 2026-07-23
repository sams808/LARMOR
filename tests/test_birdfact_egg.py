"""The birdfact easter egg (vendored from github.com/sams808/XFact): it must
be birds-only, produce a card offline (no network), and open from the ? menu."""
import os

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("requests")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_only_birds_pack_is_vendored():
    from larmor.xfact import list_packs
    assert list_packs() == ["birds"]


def test_offline_card_needs_no_network():
    from larmor.xfact.packs.birds import offline
    card = offline.load_offline_card()
    assert card["category"] == "bird" and card["tier"] == "offline"
    assert card["image_bytes"] and len(card["image_bytes"]) > 1000   # a real jpg
    assert card["name"] and card["headline"] and card["grid"]


def test_more_menu_item_opens_birdfact(monkeypatch):
    from PySide6.QtWidgets import QApplication, QMenu

    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    # force the fully-offline path (no live iNaturalist / cataas calls)
    import larmor.xfact.core.http as http
    monkeypatch.setattr(http, "get_json", lambda *a, **k: None)
    monkeypatch.setattr(http, "get_bytes", lambda *a, **k: None)

    app = QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow
    win = MainWindow()
    help_menu = next(m for m in win.menuBar().findChildren(QMenu)
                     if m.title() == "&?")
    labels = [a.text() for a in help_menu.actions() if a.text()]
    assert "More…" in labels
    # positioned directly below About LARMOR
    assert labels.index("More…") == labels.index("About LARMOR") + 1

    win._birdfact()                      # trigger the egg
    dlg = win._birdfact_dlg
    dlg._worker.wait(4000)               # let the offline fallback resolve
    app.processEvents()
    assert dlg.name_lab.text()           # a bird name got shown
    dlg.close(); win.close()
