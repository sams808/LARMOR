"""Offscreen tests of the Session inventory window (larmor.desktop.
inventory_dialog): scan, grid, detail, the manual pick override, the batch /
sequential hand-off payloads, open / copy / export; the Explorer's folder
context menu; and the MainWindow route (menu row, kept-alive dialog, the
signals reaching run_batch_fit / load_source)."""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtWidgets import (  # noqa: E402
    QApplication, QFileDialog, QMenu, QTreeWidgetItem,
)

from test_inventory import A, D, _month  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _labels(d):
    return [d.grid.verticalHeaderItem(i).text() for i in range(d.grid.rowCount())]


def _expnos(d):
    return [d.detail.item(i, 1).text() for i in range(d.detail.rowCount())]


def test_dialog_scans_fills_grid_and_hands_picks_to_batch(qapp, tmp_path):
    from larmor.desktop.inventory_dialog import SessionInventoryDialog

    root = _month(tmp_path)
    d = SessionInventoryDialog(None, str(root))
    try:
        d.chkSr.setChecked(False)
        assert not d.btnBatch.isEnabled() and not d.btnExport.isEnabled()
        d.scan()
        assert d.grid.rowCount() == 5
        assert [d.grid.horizontalHeaderItem(j).text()
                for j in range(d.grid.columnCount())] == ["27Al", "31P", "1H"]
        labels = _labels(d)
        assert labels == ["Base0Ca", "P1-Bi1-12", "P5-Bi1-12",
                          "P5-Bi8-12 (04272026)", "P5-Bi8-12 (05082026)"]
        ia = labels.index("P5-Bi8-12 (04272026)")
        cell = d.grid.item(ia, 0).text()
        assert cell.startswith("2704") and "NS 512" in cell and "(+3)" in cell
        assert "2799" in d.grid.item(ia, 0).toolTip()          # the demotion, verbatim
        assert d.grid.item(ia, 2).text() == "—"                 # no 1H: dim placeholder
        assert d.grid.item(labels.index("P5-Bi1-12"), 2).text().startswith("1 ·")
        assert "21 EXPNOs in 5 sample folders" in d.status.text()
        assert "6 picks" in d.status.text()
        assert d.windowTitle().endswith("2026-05")
        assert d.cmbNucleus.currentText() == "27Al"             # first non-1H nucleus
        assert len(d.current_picks()) == 4                      # A, B, D, E
        d.cmbNucleus.setCurrentText("31P")
        assert d.btnBatch.isEnabled() and d.btnSeq.isEnabled()
        got, got2 = [], []
        d.batch_requested.connect(lambda paths: got.append(list(paths)))
        d.seq_requested.connect(lambda paths: got2.append(list(paths)))
        d._batch()
        exp = [str(root / A / "3114" / "pdata" / "1" / "1r"),
               str(root / D / "3102" / "pdata" / "1" / "1r")]
        assert got == [exp]
        d._seq()
        assert got2 == [exp]
        assert "sent to Sequential fit" in d.status.text()
        d._column_clicked(0)                                    # header click -> nucleus
        assert d.cmbNucleus.currentText() == "27Al"
        # the NS fraction re-picks the grid without a rescan
        d.nsFrac.setValue(0.01)
        assert d.grid.item(ia, 0).text().startswith("2799")
        d.nsFrac.setValue(0.25)
        assert d.grid.item(ia, 0).text().startswith("2704")
        # an empty folder: the status says what to pick instead
        empty = tmp_path / "empty"
        empty.mkdir()
        d.folder.setText(str(empty))
        d.scan()
        assert "no EXPNO found" in d.status.text() and not d.btnBatch.isEnabled()
    finally:
        d.close()


def test_detail_pick_override_is_exclusive_open_copy_export(qapp, tmp_path, monkeypatch):
    from larmor.desktop.inventory_dialog import SessionInventoryDialog

    root = _month(tmp_path)
    d = SessionInventoryDialog(None, str(root))
    try:
        d.chkSr.setChecked(False)
        d.scan()
        assert d.detail.rowCount() == 21                        # nothing selected: all rows
        ia = _labels(d).index("P5-Bi8-12 (04272026)")
        d.grid.setCurrentCell(ia, 0)
        d._grid_selected()
        ex = _expnos(d)
        assert ex == ["2701", "2702", "2704", "2799", "3102", "3103", "3104", "3112", "3113", "3114"]
        assert d.detail.item(ex.index("2704"), 1).font().bold()
        assert d.detail.item(ex.index("2704"), 0).checkState() == Qt.Checked
        assert d.detail.item(ex.index("2799"), 3).text() == "short"
        assert d.detail.item(ex.index("2701"), 3).text() == "arrayed"
        assert not (d.detail.item(ex.index("2701"), 0).flags() & Qt.ItemIsUserCheckable)
        assert d.detail.item(ex.index("3114"), 6).text() == "300"      # D1 (s)
        # tick 2702: exclusive per block, the grid cell follows
        d.detail.item(ex.index("2702"), 0).setCheckState(Qt.Checked)
        ex = _expnos(d)
        assert d.detail.item(ex.index("2702"), 0).checkState() == Qt.Checked
        assert d.detail.item(ex.index("2704"), 0).checkState() == Qt.Unchecked
        assert d.detail.item(ex.index("2702"), 3).text() == "production"
        assert d.detail.item(ex.index("2704"), 3).text() == "candidate"
        assert d.grid.item(ia, 0).text().startswith("2702")
        assert "chosen by hand" in d.status.text()
        # open the selected detail row
        opened = []
        d.open_requested.connect(opened.append)
        d.detail.setCurrentCell(ex.index("2702"), 1)
        d._open()
        assert opened == [str(root / A / "2702" / "pdata" / "1" / "1r")]
        d._detail_double(d.detail.item(ex.index("3114"), 5))
        assert opened[-1] == str(root / A / "3114" / "pdata" / "1" / "1r")
        d._grid_double(ia, 1)
        assert opened[-1] == str(root / A / "3114" / "pdata" / "1" / "1r")
        # copy: tab-separated picks of the hand-off nucleus
        d.cmbNucleus.setCurrentText("31P")
        d._copy()
        clip = QApplication.clipboard().text()
        assert clip.startswith("# sample\tnucleus\tEXPNO\tpath") and "\t31P\t3114\t" in clip
        assert "\t27Al\t" not in clip
        # export: CSV + picks list next to it
        out = tmp_path / "inventory_2026-05.csv"
        monkeypatch.setattr(QFileDialog, "getSaveFileName",
                            staticmethod(lambda *a, **k: (str(out), "CSV (*.csv)")))
        d.export()
        assert out.exists() and (tmp_path / "inventory_2026-05_picks.txt").exists()
        assert "wrote inventory_2026-05.csv and inventory_2026-05_picks.txt" in d.status.text()
        assert ",production," in out.read_text(encoding="utf-8")
        # untick the pick: the sample leaves the hand-off
        d.cmbNucleus.setCurrentText("27Al")
        n_before = len(d.current_picks())
        ex = _expnos(d)
        d.detail.item(ex.index("2702"), 0).setCheckState(Qt.Unchecked)
        assert len(d.current_picks()) == n_before - 1
        assert d.grid.item(ia, 0).text() == "—" and "left out" in d.status.text()
    finally:
        d.close()


def test_explorer_context_menu_emits_inventory_request_for_folders(qapp, tmp_path, monkeypatch):
    from larmor.desktop import explorer as ex

    root = _month(tmp_path)
    panel = ex.ExplorerPanel()
    seen = {}

    class _Menu(QMenu):                      # no popup loop: record and return
        def exec(self, *a, **k):
            seen["menu"] = self
            return None

    monkeypatch.setattr(ex, "QMenu", _Menu)
    got = []
    panel.inventory_requested.connect(got.append)
    folder = QTreeWidgetItem([root.name])
    folder.setData(0, ex._ROLE_PATH, str(root))
    folder.setData(0, ex._ROLE_KIND, "folder")
    panel.tree.addTopLevelItem(folder)
    monkeypatch.setattr(panel.tree, "itemAt", lambda pos: folder)
    panel._context_menu(panel.tree.rect().center())
    texts = [a.text() for a in seen["menu"].actions() if not a.isSeparator()]
    assert any(t.startswith("Session inventory") for t in texts)
    next(a for a in seen["menu"].actions() if a.text().startswith("Session inventory")).trigger()
    assert got == [str(root)]
    # an EXPNO item offers Dataset info, not the inventory
    expno = QTreeWidgetItem(["2702"])
    expno.setData(0, ex._ROLE_PATH, str(root / A / "2702"))
    expno.setData(0, ex._ROLE_KIND, "exp")
    panel.tree.addTopLevelItem(expno)
    monkeypatch.setattr(panel.tree, "itemAt", lambda pos: expno)
    panel._context_menu(panel.tree.rect().center())
    texts = [a.text() for a in seen["menu"].actions() if not a.isSeparator()]
    assert not any(t.startswith("Session inventory") for t in texts)
    assert any(t.startswith("Dataset info") for t in texts)
    panel.close()


def test_app_menu_opens_inventory_and_routes_to_batch_fit(qapp, tmp_path, monkeypatch):
    os.environ["LARMOR_NO_SESSION"] = "1"
    from larmor.desktop.app import MainWindow

    root = _month(tmp_path)
    win = MainWindow()
    try:
        labels = [a.text() for m in win.menuBar().findChildren(QMenu) for a in m.actions()]
        assert any("Session inventory" in t for t in labels)
        received, loaded = [], []
        monkeypatch.setattr(win, "run_batch_fit", lambda paths: received.append(list(paths)))
        monkeypatch.setattr(win, "load_source", lambda path, *a, **k: loaded.append(path))
        win.open_session_inventory()
        dlg = win._inventory_dlg
        assert dlg is not None and dlg.isVisible()
        win.open_session_inventory(False)                # QAction.triggered's bool
        assert win._inventory_dlg is dlg                  # kept alive, one instance
        win.open_session_inventory(str(root / A / "2702"))
        assert dlg.folder.text() == str(root)             # an EXPNO -> its month
        dlg.batch_requested.emit(["a", "b", "c"])
        assert received == [["a", "b", "c"]]
        dlg.open_requested.emit("x")
        assert loaded == ["x"]
        # the Explorer's folder context menu lands in the same window
        win.explorer.inventory_requested.emit(str(root / A))
        assert win._inventory_dlg is dlg and dlg.folder.text() == str(root)
        # the real hand-off: scan, pick a nucleus, send
        dlg.chkSr.setChecked(False)
        dlg.scan()
        dlg.cmbNucleus.setCurrentText("31P")
        dlg._batch()
        assert received[-1] == [str(root / A / "3114" / "pdata" / "1" / "1r"),
                                str(root / D / "3102" / "pdata" / "1" / "1r")]
    finally:
        win.close()
