"""The ? ▸ More… card (vendored xfact): it must be a single pack, produce a
card offline (no network), and open from the ? menu."""
import os

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("requests")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_single_pack_is_vendored():
    from larmor.xfact import list_packs
    assert list_packs() == ["birds"]


def test_offline_card_needs_no_network():
    from larmor.xfact.packs.birds import offline
    card = offline.load_offline_card()
    assert card["category"] == "bird" and card["tier"] == "offline"
    assert card["image_bytes"] and len(card["image_bytes"]) > 1000   # a real jpg
    assert card["name"] and card["headline"] and card["grid"]


def test_offline_pool_is_varied():
    """The offline fallback must have real variety — it used to be only three
    samples (so a network-blocked user saw the same one or two over and over)."""
    from larmor.xfact.packs.birds import offline
    names = {c["name"] for c in offline.OFFLINE_CARDS}
    assert len(names) >= 10
    for c in offline.OFFLINE_CARDS:                 # every sample's jpg is present
        assert (offline.ASSETS_DIR / c["image_file"]).exists()


def test_more_menu_item_opens_card(monkeypatch):
    from PySide6.QtWidgets import QApplication, QMenu

    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    # force the fully-offline path (no live network calls)
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

    win._show_more()
    dlg = win._more_dlg
    dlg._worker.wait(4000)               # let the offline fallback resolve
    app.processEvents()
    assert dlg.name_lab.text()           # a card got shown
    dlg.close(); win.close()
