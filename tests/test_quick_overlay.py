"""Quick overlays in the main window: Shift + drop, File > Overlay a
spectrum..., View > Overlays / Clear overlays, the Explorer's right-click and
the Datasets dock's 'match height' -- all on top of an active spectrum whose
fit stays untouched."""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QMimeData, QPointF, QUrl, Qt  # noqa: E402
from PySide6.QtGui import QDropEvent  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from larmor.io import spectra  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _active(win):
    ppm = np.linspace(80.0, -20.0, 600)
    amp = 100.0 * np.exp(-((ppm - 30.0) / 6.0) ** 2)
    win._display_1d(ppm, amp, "27Al", 130.3, 20000.0, "active", "active")
    return ppm, amp


def _csv(tmp_path, name, scale):
    ppm = np.linspace(80.0, -20.0, 500)
    amp = scale * np.exp(-((ppm - 10.0) / 4.0) ** 2)
    return spectra.write_csv(tmp_path / name, ppm, amp, {"nucleus": "27Al"}), amp


def _drop(win, paths, modifiers):
    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(p)) for p in paths])
    ev = QDropEvent(QPointF(20.0, 20.0), Qt.CopyAction, mime, Qt.LeftButton, modifiers)
    win.view.dropEvent(ev)


def test_shift_drop_overlays_and_plain_drop_still_loads(win, tmp_path, monkeypatch):
    _active(win)
    p1, _ = _csv(tmp_path, "ref1.csv", 5.0)
    p2, _ = _csv(tmp_path, "ref2.csv", 7.0)
    loaded = []
    monkeypatch.setattr(win, "load_source", lambda path, keep_fit=None: loaded.append(path))
    win.view.file_dropped.disconnect()
    win.view.file_dropped.connect(win.load_source)
    _drop(win, [p1, p2], Qt.ShiftModifier)
    assert loaded == []
    assert [ov["label"] for ov in win._overlays] == ["ref1", "ref2"]
    assert len(win.view._overlay_items) == 2
    assert "2 of 2 dropped" in win.statusBar().currentMessage()
    # the active spectrum is exactly what it was
    assert win.recipe["nucleus"] == "27Al" and win.exp_ppm.size == 600
    _drop(win, [p1], Qt.NoModifier)
    from pathlib import Path
    assert [Path(p) for p in loaded] == [Path(p1)]     # no Shift: the ordinary load path
    assert len(win._overlays) == 2


def test_menu_toggle_clear_and_match_height(win, tmp_path):
    ppm, amp = _active(win)
    path, ref = _csv(tmp_path, "ref.csv", 5.0)
    assert win.add_overlay_path(str(path)) is True
    assert "overlaid ref" in win.statusBar().currentMessage()
    item = win.view._overlay_items[0]
    ref_max = float(ref.max())                  # the sampled maximum (< 5.0); CSV keeps ~7 digits
    assert item.isVisible() and item.yData.max() == pytest.approx(ref_max, rel=1e-5)
    # match height: drawn at the active maximum, stored array untouched
    win.datasets_panel.match.setChecked(True)
    item = win.view._overlay_items[0]
    assert item.yData.max() == pytest.approx(float(amp.max()), rel=1e-9)   # the active peak
    assert win._overlays[0]["amp"].max() == pytest.approx(ref_max, rel=1e-5)
    win.datasets_panel.match.setChecked(False)
    assert win.view._overlay_items[0].yData.max() == pytest.approx(ref_max, rel=1e-5)
    # View > Overlays hides without removing; adding one turns it back on
    win.actOverlaysVisible.trigger()                     # checked -> unchecked, like a click
    assert not win.actOverlaysVisible.isChecked()
    assert not win.view._overlay_items[0].isVisible() and len(win._overlays) == 1
    win.add_overlay_path(str(path))
    assert win.actOverlaysVisible.isChecked()
    assert all(it.isVisible() for it in win.view._overlay_items)
    win.clear_overlays()
    assert win._overlays == [] and win.view._overlay_items == []
    assert "removed 2" in win.statusBar().currentMessage()
    # an unreadable path is a status line, never a dialog
    assert win.add_overlay_path(str(tmp_path / "nope.csv")) is False
    assert "cannot overlay nope.csv" in win.statusBar().currentMessage()


def test_explorer_request_and_menu_entries(win, tmp_path):
    _active(win)
    path, _ = _csv(tmp_path, "ref.csv", 3.0)
    win.explorer.overlay_requested.emit(str(path))
    assert [ov["label"] for ov in win._overlays] == ["ref"]
    texts = [a.text() for m in win.menuBar().findChildren(type(win.menuBar().actions()[0].menu()))
             for a in m.actions()]
    assert any(t.startswith("O&verlay a spectrum") for t in texts)
    assert any(t.startswith("O&verlays") for t in texts) and "Clear overlays" in texts
    shortcuts = {a.shortcut().toString() for m in win.menuBar().findChildren(
        type(win.menuBar().actions()[0].menu())) for a in m.actions() if a.shortcut().toString()}
    assert {"Ctrl+Shift+A", "Ctrl+Shift+V"} <= shortcuts
