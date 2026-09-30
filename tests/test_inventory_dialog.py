"""Offscreen tests of the Session inventory window (larmor.desktop.
inventory_dialog): scan, grid, detail, the manual pick override, the batch /
sequential hand-off payloads, open / copy / export; the Fix… menu (the
remedies of larmor.inventory.fixes_for as actions, the SR override written
and undone, the TopSpin command, the 1H reference, the rename through the
Explorer's RenameDialog, the fid hand-off, pick anyway, dataset info); the
Explorer's folder context menu; and the MainWindow route (menu row,
kept-alive dialog, the signals reaching run_batch_fit / load_source /
open_fid_path / _on_inventory_renamed). Every store is a tmp file."""
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtGui import QColor  # noqa: E402
from PySide6.QtWidgets import (  # noqa: E402
    QApplication, QDialog, QFileDialog, QMenu, QMessageBox, QTreeWidgetItem,
)

from test_inventory import A, B, C, CRY, D, E, _month, _sr_month  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXPNO_3616 = ROOT / "examples" / "pCABS2-4" / "3616"


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def stores(tmp_path, monkeypatch):
    """Own SR-override store, referencing log, alias file and rename log."""
    monkeypatch.setenv("LARMOR_SR_OVERRIDES", str(tmp_path / "store" / "sr_overrides.json"))
    monkeypatch.setenv("LARMOR_REF_LOG", str(tmp_path / "store" / "referencing_log.jsonl"))
    monkeypatch.setenv("LARMOR_ALIASES", str(tmp_path / "store" / "aliases.json"))
    monkeypatch.setenv("LARMOR_RENAME_LOG", str(tmp_path / "store" / "rename_log.jsonl"))
    return tmp_path


def _texts(menu):
    return [a.text() for a in menu.actions() if not a.isSeparator()]


def _action(menu, prefix):
    return next(a for a in menu.actions() if a.text().startswith(prefix))


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


def test_fix_menu_offers_and_applies_the_remedies(qapp, stores, monkeypatch):
    from larmor import referencing as R
    from larmor.desktop import explorer as ex
    from larmor.desktop.inventory_dialog import _BG_BAD, WRENCH, SessionInventoryDialog

    root = _sr_month(stores)
    d = SessionInventoryDialog(None, str(root))
    opened, fids = [], []
    d.open_requested.connect(opened.append)
    d.fid_requested.connect(fids.append)
    try:
        d.scan()                                            # with the SR audit
        assert "1 pick unreferenced/off" in d.status.text()
        assert "right-click a row (or press Fix…)" in d.status.text()
        assert d.btnFix.isEnabled()
        ib = _labels(d).index("P1-Bi1-12")
        nucs = [d.grid.horizontalHeaderItem(j).text() for j in range(d.grid.columnCount())]
        assert nucs == ["23Na", "27Al", "31P", "1H"]
        jal = nucs.index("27Al")
        assert d.grid.item(ib, jal).background().color() == _BG_BAD
        assert "Apply the audited SR in LARMOR" in d.grid.item(ib, jal).toolTip()
        d.grid.setCurrentCell(ib, jal)
        d._grid_selected()
        ex_ = _expnos(d)
        i2702 = ex_.index("2702")
        assert d.detail.item(i2702, 10).text() == "unreferenced"
        assert d.detail.item(i2702, 12).text().startswith(WRENCH)
        assert "Apply the audited SR in LARMOR" in d.detail.item(i2702, 12).toolTip()
        assert not d.detail.item(ex_.index("2701"), 12).text()      # nothing to fix
        d.detail.setCurrentCell(i2702, 1)
        row = d.current_row()
        assert row is not None and row.expno == 2702
        m = d.fix_menu(row)
        assert isinstance(m, QMenu)
        assert _texts(m) == [
            "Apply the audited SR in LARMOR (0.00 → -135.31 Hz)",
            "Copy the TopSpin command:  sr -135.31",
            "Open the ¹H reference (04272026_RS40339_P5-Bi1-12/1)",
            "Apply the audited SR to every flagged pick… (1)",
            "Dataset info…", "Open pick", "Batch fit picks…", "Sequential fit picks…",
            "Copy picks", "Export CSV…"]
        assert "reversible" in _action(m, "Apply the audited SR in LARMOR").toolTip()
        # the TopSpin command lands on the clipboard; the reference opens
        _action(m, "Copy the TopSpin").trigger()
        assert QApplication.clipboard().text() == "sr -135.31" and "TopSpin" in d.status.text()
        _action(m, "Open the ¹H reference").trigger()
        assert opened == [str(root / C / "1" / "pdata" / "1" / "1r")]
        # apply: the store is written, the log appended, the cell repainted green
        _action(m, "Apply the audited SR in LARMOR").trigger()
        ov = R.override_for(root / B / "2702")
        assert ov and ov["new_sr_hz"] == pytest.approx(-135.31, abs=0.02) and ov["old_sr_hz"] == 0.0
        assert ov["nucleus"] == "27Al" and "RS40339" in ov["note"]
        assert [r["action"] for r in R.previous_audits(root)] == ["override"]
        ex_ = _expnos(d)
        cell = d.detail.item(ex_.index("2702"), 10)
        assert cell.text() == "corrected (LARMOR)"
        assert cell.foreground().color() == QColor("#2C6A4E")
        assert "LARMOR applies SR -135.31 Hz" in cell.toolTip()
        assert d.grid.item(ib, jal).background().color() != _BG_BAD
        assert "corrected (LARMOR)" in d.grid.item(ib, jal).toolTip()
        assert "0 picks unreferenced/off" in d.status.text()
        assert "1 SR corrected in LARMOR" in d.status.text() and "Undo" in d.status.text()
        d.detail.setCurrentCell(ex_.index("2702"), 1)
        m2 = d.fix_menu(d.current_row())
        assert _texts(m2)[0] == "Undo the LARMOR SR correction (0.00 → -135.31 Hz)"
        assert not any(t.startswith("Apply the audited SR to every") for t in _texts(m2))
        _action(m2, "Undo the LARMOR").trigger()
        assert R.override_for(root / B / "2702") is None
        assert d.detail.item(_expnos(d).index("2702"), 10).text() == "unreferenced"
        assert [r["action"] for r in R.previous_audits(root)] == ["override", "override-cleared"]
        # the bulk entry: a confirmation that lists the picks, Cancel writes nothing
        asked = []
        monkeypatch.setattr(QMessageBox, "question",
                            staticmethod(lambda *a, **k: (asked.append(a[2]), QMessageBox.Cancel)[1]))
        _action(d.fix_menu(None), "Apply the audited SR to every").trigger()
        assert asked and "EXPNO 2702" in asked[0] and "0.00 → -135.31" in asked[0]
        assert R.override_for(root / B / "2702") is None
        monkeypatch.setattr(QMessageBox, "question",
                            staticmethod(lambda *a, **k: QMessageBox.Yes))
        _action(d.fix_menu(None), "Apply the audited SR to every").trigger()
        assert R.override_for(root / B / "2702")["new_sr_hz"] == pytest.approx(-135.31, abs=0.02)
        assert "1 SR correction(s) recorded" in d.status.text()
        R.clear_override(root / B / "2702")
        d._refresh_sr(str(root / B / "2702"))
        # pick anyway: the short 2703 becomes the block's production spectrum
        d.detail.setCurrentCell(_expnos(d).index("2703"), 1)
        m3 = d.fix_menu(d.current_row())
        assert _texts(m3)[0] == "Pick EXPNO 2703 anyway"
        _action(m3, "Pick EXPNO 2703").trigger()
        ex_ = _expnos(d)
        assert d.detail.item(ex_.index("2703"), 3).text() == "production"
        assert d.detail.item(ex_.index("2702"), 3).text() == "candidate"
        assert d.grid.item(ib, jal).text().startswith("2703") and "chosen by hand" in d.status.text()
        # the fid hand-off: the unprocessed 31P itself and from the pick it outranks
        d.detail.setCurrentCell(ex_.index("3102"), 1)
        m4 = d.fix_menu(d.current_row())
        assert _texts(m4)[0] == "Process the fid of EXPNO 3102 in LARMOR…"
        _action(m4, "Process the fid").trigger()
        assert fids == [str(root / B / "3102" / "fid")] and "Open FID" in d.status.text()
        d.detail.setCurrentCell(ex_.index("3101"), 1)
        assert _texts(d.fix_menu(d.current_row()))[0] == \
            "Process the fid of EXPNO 3102 (NS 140) in LARMOR…"
        # dataset info: the shared read-only window
        seen = []
        monkeypatch.setattr(ex, "show_dataset_info", lambda parent, path: seen.append(path))
        _action(d.fix_menu(d.current_row()), "Dataset info").trigger()
        assert seen == [str(root / B / "3101")]
        # a right-click on the detail table and on a grid cell builds the same menu
        pops = []

        class _Menu(QMenu):
            def exec(self, *a, **k):
                pops.append(_texts(self))
                return None

        import larmor.desktop.inventory_dialog as mod
        monkeypatch.setattr(mod, "QMenu", _Menu)
        r2702 = d.detail.item(_expnos(d).index("2702"), 1)
        monkeypatch.setattr(d.detail, "itemAt", lambda pos: r2702)
        d._detail_context_menu(d.detail.rect().center())
        assert pops[-1][0].startswith("Apply the audited SR in LARMOR")
        monkeypatch.setattr(d.grid, "itemAt", lambda pos: d.grid.item(ib, jal))
        d._context_menu(d.grid.rect().center())
        assert pops[-1][0] == "Dataset info…"       # the cell's pick is the clean 2703 now
        assert "Open pick" in pops[-1] and not any(t.startswith("Apply") for t in pops[-1])
        d.detail.setCurrentCell(_expnos(d).index("2702"), 1)
        d._fix_button()
        assert len(pops) == 3 and pops[-1][0].startswith("Apply the audited SR in LARMOR")
    finally:
        d.close()


def test_fix_menu_rename_alias_and_on_disk(qapp, stores, monkeypatch):
    from larmor import aliases
    from larmor.desktop import explorer as ex
    from larmor.desktop.inventory_dialog import SessionInventoryDialog

    root = _sr_month(stores)
    d = SessionInventoryDialog(None, str(root))
    got = []
    d.renamed.connect(lambda a, b: got.append((a, b)))
    try:
        d.chkSr.setChecked(False)
        d.scan()
        # (a) the rotor flag on E/33 -> a display name for the sample folder
        d.grid.setCurrentCell(_labels(d).index("Base0Ca"), 0)
        d._grid_selected()
        d.detail.setCurrentCell(_expnos(d).index("33"), 1)
        m = d.fix_menu(d.current_row())
        assert _texts(m)[0] == f"Rename the sample folder {E}…"

        def accept_alias(self):
            self.edit.setText("Base zero")
            self.rbAlias.setChecked(True)
            return QDialog.Accepted

        monkeypatch.setattr(ex.RenameDialog, "exec", accept_alias)
        _action(m, "Rename the sample folder").trigger()
        assert aliases.alias_for(root / E) == "Base zero" and (root / E).is_dir()
        assert "Base zero" in _labels(d) and "Base0Ca" not in _labels(d)
        assert got == [(str(root / E), str(root / E))]
        assert "display name 'Base zero' set" in d.status.text()
        assert d._sample == (E, "Base zero") and "33" in _expnos(d)   # the detail followed
        # (b) the nucleus-token flag on Cryolite/16 -> a display name for the EXPNO
        d.grid.setCurrentCell(_labels(d).index("Cryolite"), 0)
        d._grid_selected()
        d.detail.setCurrentCell(_expnos(d).index("16"), 1)
        m = d.fix_menu(d.current_row())
        assert _texts(m)[0] == "Rename EXPNO 16…"

        def accept_expno_alias(self):
            self.edit.setText("was a 1H title")
            self.rbAlias.setChecked(True)
            return QDialog.Accepted

        monkeypatch.setattr(ex.RenameDialog, "exec", accept_expno_alias)
        _action(m, "Rename EXPNO 16").trigger()
        assert aliases.alias_for(root / CRY / "16") == "was a 1H title"
        assert d.detail.item(0, 1).text() == "16 (was a 1H title)"
        assert got[-1] == (str(root / CRY / "16"), str(root / CRY / "16"))
        # (c) on disk: the confirmation states both paths, the session is rescanned
        d.grid.setCurrentCell(_labels(d).index("Base zero"), 0)
        d._grid_selected()
        d.detail.setCurrentCell(_expnos(d).index("33"), 1)
        new_name = "01192026_SR31648_Base0Ca_SS_ALP"          # the rotor the title says

        def accept_disk(self):
            self.edit.setText(new_name)
            self.rbDisk.setChecked(True)
            return QDialog.Accepted

        monkeypatch.setattr(ex.RenameDialog, "exec", accept_disk)
        asked = []
        monkeypatch.setattr(QMessageBox, "question",
                            staticmethod(lambda *a, **k: (asked.append(a[2]), QMessageBox.Yes)[1]))
        _action(d.fix_menu(d.current_row()), "Rename the sample folder").trigger()
        assert asked and str(root / E) in asked[0] and str(root / new_name) in asked[0]
        assert (root / new_name / "33" / "acqus").exists() and not (root / E).exists()
        assert got[-1] == (str(root / E), str(root / new_name))
        assert aliases.alias_for(root / new_name) == "Base zero"     # the alias followed
        assert "Base zero" in _labels(d) and "renamed on disk" in d.status.text()
        assert aliases.rename_log_path().exists()
        r33 = next(r for r in d.inv.rows if r.expno == 33)
        assert r33.folder == new_name and not any("rotor" in f for f in r33.flags)
        # a folder holding an open document is refused before any confirmation
        asked.clear()
        warned = []
        monkeypatch.setattr(QMessageBox, "warning",
                            staticmethod(lambda *a, **k: warned.append(a[2])))
        d.open_paths = lambda: [str(root / new_name / "33" / "pdata" / "1" / "1r")]
        d.grid.setCurrentCell(_labels(d).index("Base zero"), 0)
        d._grid_selected()
        d.detail.setCurrentCell(_expnos(d).index("33"), 1)
        r33.flags.append("title says rotor SR1 but the folder is SR2")   # a flag to act on

        def accept_disk_v2(self):
            self.edit.setText(new_name + "_v2")
            self.rbDisk.setChecked(True)
            return QDialog.Accepted

        monkeypatch.setattr(ex.RenameDialog, "exec", accept_disk_v2)
        _action(d.fix_menu(d.current_row()), "Rename the sample folder").trigger()
        assert warned and "open in LARMOR" in warned[0] and not asked
        assert (root / new_name).exists() and not (root / (new_name + "_v2")).exists()
        assert got[-1] == (str(root / E), str(root / new_name))        # nothing new emitted
    finally:
        d.close()


def test_show_dataset_info_is_a_non_modal_window(qapp, tmp_path):
    from PySide6.QtWidgets import QPlainTextEdit

    from larmor.desktop import explorer as ex

    root = _month(tmp_path)
    w = ex.show_dataset_info(None, str(root / A / "3114"))
    try:
        assert w.windowTitle() == "3114 — dataset info" and not w.isModal()
        txt = w.findChild(QPlainTextEdit)
        assert txt.isReadOnly() and "Rotor RS2427418" in txt.toPlainText()
    finally:
        w.close()


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
        # the Fix… menu's signals: a fid to open_fid_path, a rename to the
        # Explorer and the document relabelling; the open paths are shared
        fids = []
        monkeypatch.setattr(win, "open_fid_path", lambda p: fids.append(p))
        dlg.fid_requested.emit("some/fid")
        assert fids == ["some/fid"]
        assert dlg.open_paths == win._open_source_paths
        relabelled = []
        monkeypatch.setattr(win.explorer, "apply_rename", lambda a, b: relabelled.append((a, b)))
        dlg.renamed.emit(str(root / A), str(root / A))
        assert relabelled == [(str(root / A), str(root / A))]
        assert "display name set" in win.statusBar().currentMessage()
    finally:
        win.close()


def test_open_fid_path_opens_the_fid_dialog_wired_to_the_workbench(qapp, monkeypatch):
    if not (EXPNO_3616 / "fid").exists():
        pytest.skip("examples/pCABS2-4/3616 not shipped here")
    os.environ["LARMOR_NO_SESSION"] = "1"
    from larmor.desktop import fid_dialog
    from larmor.desktop.app import MainWindow

    shown = []
    monkeypatch.setattr(fid_dialog.FidDialog, "exec", lambda self: shown.append(self) or 0)
    win = MainWindow()
    try:
        win.open_fid_path(str(EXPNO_3616 / "fid"))
        assert len(shown) == 1
        d = shown[0]
        assert d.lbl.text() == str(EXPNO_3616 / "fid") and d.data is not None
        assert d.data.ndim == 1 and d.data.domain == "time"
        d._accept()                                   # accepted_1d -> _fid_to_workbench
        assert win.source_path == str(EXPNO_3616)
        assert win.recipe["source_kind"] == "bruker" and win.recipe["nucleus"] == "27Al"
        assert win.exp_ppm.size > 0 and win.recipe["processing_from_raw"]
    finally:
        win.close()
