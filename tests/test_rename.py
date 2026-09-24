"""Explorer ▸ Rename…: a display name kept by LARMOR (an alias) or the folder
renamed on disk.

The alias store (larmor.aliases, Qt-free) is read by scan.sample_name /
sample_label, so the recipe's sample, batch labels and the window title all
follow it; the Explorer tree shows it and offers Rename… on sample folders
and EXPNOs. A rename on disk goes through check_rename (target exists, a
file inside is open, an EXPNO stays a number, a plain name) and is appended
to rename_log.jsonl. Every disk operation here runs on a tmp_path copy of a
fake EXPNO -- never on instrument data.
"""
import json
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from larmor import aliases  # noqa: E402
from larmor.io import scan  # noqa: E402

FOLDER = "01192026_SR31649_Base0Ca_SS_ALP"


def _acqus_text() -> str:
    """A COMPLETE JCAMP acqus. nmrglue's read_jcamp (behind
    comparability.read_params, which the Datasets dock refresh calls) loops
    forever on a truncated ``(0..63)`` array, so every array carries its 64
    values and the file ends with ``##END=``."""
    d = ["5", "0.5", "1"] + ["0"] * 61
    pulses = ["0", "2.5"] + ["0"] * 62
    lines = ["##TITLE= Parameter file, TopSpin 4.1", "##JCAMPDX= 5.0",
             "##DATATYPE= Parameter Values", "##ORIGIN= Bruker BioSpin GmbH",
             "##$NUC1= <27Al>", "##$SFO1= 130.3", "##$BF1= 130.3", "##$O1= 0",
             "##$SW_h= 100000", "##$TD= 4096", "##$NS= 16", "##$RG= 64",
             "##$PULPROG= <zg>", "##$PARMODE= 0", "##$DATE= 1780000000",
             "##$D= (0..63)", " ".join(d), "##$P= (0..63)", " ".join(pulses),
             "##END="]
    return "\n".join(lines) + "\n"


def _fake_expno(root: Path, folder: str = FOLDER, expno: str = "10") -> Path:
    e = root / folder / expno
    (e / "pdata" / "1").mkdir(parents=True)
    (e / "acqus").write_text(_acqus_text())
    (e / "pdata" / "1" / "1r").write_bytes(b"\0" * 8)
    (e / "pdata" / "1" / "title").write_text("27Al MAS\n")
    return e


@pytest.fixture()
def store(tmp_path, monkeypatch):
    """Own alias file and rename log per test -- the developer's
    %LOCALAPPDATA% store is never read or written."""
    monkeypatch.setenv("LARMOR_ALIASES", str(tmp_path / "store" / "aliases.json"))
    monkeypatch.setenv("LARMOR_RENAME_LOG", str(tmp_path / "store" / "rename_log.jsonl"))
    return tmp_path


# ---------------------------------------------------------------- Qt-free core
def test_alias_round_trip_through_scan(store):
    expno = _fake_expno(store)
    sample = expno.parent
    r1 = str(expno / "pdata" / "1" / "1r")
    assert scan.sample_name(expno, "27Al MAS").key == "Base0Ca"
    assert aliases.alias_for(sample) == "" and aliases.display_name(sample) == FOLDER

    aliases.set_alias(sample, "Glass A")
    assert aliases.aliases_path().exists()
    assert aliases.alias_for(sample) == "Glass A"
    # the key is the normalised path: case and separators do not matter
    assert aliases.alias_for(str(sample).replace("\\", "/").upper()) == "Glass A"
    name = scan.sample_name(expno, "27Al MAS")
    assert (name.key, name.source, name.folder) == ("Glass A", "alias", FOLDER)
    # the alias beats the recipe's own sample and the path rule
    assert scan.sample_label(r1, {"sample": "Base0Ca", "nucleus": "27Al"}) == "Glass A"
    assert scan.sample_label(r1, {}) == "Glass A"
    assert aliases.window_label(r1) == "Glass A · 10"

    # an EXPNO alias is more specific than the sample's
    aliases.set_alias(expno, "run 10 (echo)")
    assert scan.sample_label(r1, {"sample": "Base0Ca"}) == "run 10 (echo)"
    assert aliases.window_label(r1) == "run 10 (echo)"
    assert aliases.window_label(store / "elsewhere.csv") == "elsewhere.csv"

    # empty, or the folder's own name, removes the alias
    aliases.set_alias(expno, "")
    aliases.set_alias(sample, FOLDER)
    assert aliases.load() == {}
    assert scan.sample_name(expno, "27Al MAS").key == "Base0Ca"
    assert scan.sample_label(r1, {"sample": "typed"}) == "typed"


def test_rename_folder_on_disk_moves_alias_and_logs(store):
    expno = _fake_expno(store)
    sample = expno.parent
    aliases.set_alias(expno, "first run")
    new = aliases.rename_folder(expno, "11")
    assert new == sample / "11" and new.is_dir() and not expno.exists()
    assert (new / "pdata" / "1" / "1r").read_bytes() == b"\0" * 8   # layout intact
    assert aliases.alias_for(new) == "first run" and aliases.alias_for(expno) == ""
    lines = aliases.rename_log_path().read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    rec = json.loads(lines[0])
    assert rec["action"] == "rename" and rec["kind"] == "expno"
    assert Path(rec["old"]) == expno and Path(rec["new"]) == new
    assert rec["alias_moved"] is True and rec["larmor"]

    # the sample folder: the aliases keyed under it follow, the log grows
    aliases.set_alias(sample, "Glass A")
    new_sample = aliases.rename_folder(sample, "Base0Ca_repacked")
    assert new_sample.is_dir() and (new_sample / "11" / "acqus").exists()
    assert aliases.alias_for(new_sample) == "Glass A"
    assert aliases.alias_for(new_sample / "11") == "first run"
    assert len(aliases.rename_log_path().read_text(encoding="utf-8").splitlines()) == 2
    # renamed TO its alias: the alias is redundant and dropped
    aliases.rename_folder(new_sample, "Glass A")
    assert aliases.alias_for(store / "Glass A") == ""
    assert aliases.alias_for(store / "Glass A" / "11") == "first run"


def test_rename_refusals(store):
    expno = _fake_expno(store)
    sample = expno.parent
    other = _fake_expno(store, expno="12")
    with pytest.raises(aliases.RenameError, match="must stay a number"):
        aliases.check_rename(expno, "echo")
    with pytest.raises(aliases.RenameError, match="already exists"):
        aliases.check_rename(expno, "12")
    with pytest.raises(aliases.RenameError, match="empty"):
        aliases.check_rename(expno, "   ")
    with pytest.raises(aliases.RenameError, match="unchanged"):
        aliases.check_rename(sample, FOLDER)
    with pytest.raises(aliases.RenameError, match="plain folder name"):
        aliases.check_rename(sample, "a/b")
    with pytest.raises(aliases.RenameError, match="plain folder name"):
        aliases.check_rename(sample, "..")
    with pytest.raises(aliases.RenameError, match="is not a folder"):
        aliases.check_rename(store / "nowhere", "x")
    # a file inside is open in LARMOR (the window's own or an overlay's)
    open_1r = str(other / "pdata" / "1" / "1r")
    with pytest.raises(aliases.RenameError, match="open in LARMOR"):
        aliases.check_rename(other, "13", open_paths=[open_1r])
    with pytest.raises(aliases.RenameError, match="open in LARMOR"):
        aliases.check_rename(sample, "Renamed", open_paths=["", open_1r])
    # a sibling open elsewhere does not block
    assert aliases.check_rename(expno, "13", open_paths=[open_1r]) == sample / "13"
    assert expno.exists() and not (sample / "13").exists()    # check_rename moves nothing
    assert not aliases.rename_log_path().exists()


# --------------------------------------------------------------- the Explorer
@pytest.fixture(scope="module")
def qapp():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


def test_explorer_tree_shows_alias_and_offers_rename(qapp, store, monkeypatch):
    from PySide6.QtWidgets import QMenu

    from larmor.desktop import explorer as ex

    expno = _fake_expno(store)
    sample = expno.parent
    aliases.set_alias(sample, "Glass A")
    aliases.set_alias(expno, "first run")
    panel = ex.ExplorerPanel()
    try:
        panel.load_tree(str(store))
        root = panel.tree.topLevelItem(panel.tree.topLevelItemCount() - 1)
        panel._expanded(root)                                   # populate the month
        smp = next(root.child(i) for i in range(root.childCount())
                   if root.child(i).data(0, ex._ROLE_PATH) == str(sample))
        assert smp.text(0) == "🧪 Glass A" and smp.toolTip(0) == str(sample)
        panel._expanded(smp)
        exp = smp.child(0)
        assert exp.data(0, ex._ROLE_KIND) == "exp"
        assert exp.text(0).startswith("first run (10) · 27Al 1D")
        assert exp.data(0, ex._ROLE_BASE).startswith("10 · 27Al 1D")

        # clearing the alias relabels in place
        aliases.set_alias(expno, "")
        panel._relabel(str(expno))
        assert exp.text(0).startswith("10 · 27Al 1D")

        # Rename… is offered on the sample folder and the EXPNO, not on the month
        seen = {}

        class _Menu(QMenu):
            def exec(self, *a, **k):
                seen["menu"] = self

        monkeypatch.setattr(ex, "QMenu", _Menu)
        for item, expect in ((smp, True), (exp, True), (root, False)):
            monkeypatch.setattr(panel.tree, "itemAt", lambda pos, it=item: it)
            panel._context_menu(panel.tree.rect().center())
            texts = [a.text() for a in seen["menu"].actions() if not a.isSeparator()]
            assert ("Rename…" in texts) is expect, (item.text(0), texts)
    finally:
        panel.close()


def test_explorer_rename_dialog_alias_and_disk(qapp, store, monkeypatch):
    from PySide6.QtWidgets import QDialog, QMessageBox

    from larmor.desktop import explorer as ex

    expno = _fake_expno(store)
    sample = expno.parent
    panel = ex.ExplorerPanel()
    got = []
    panel.renamed.connect(lambda a, b: got.append((a, b)))
    try:
        panel.load_sample(str(sample))
        exp = next(it for it in panel._iter_items() if it.data(0, ex._ROLE_KIND) == "exp")
        assert exp.data(0, ex._ROLE_PATH) == str(expno)

        # (a) display name only
        def accept_alias(self):
            self.edit.setText("first run")
            self.rbAlias.setChecked(True)
            return QDialog.Accepted

        monkeypatch.setattr(ex.RenameDialog, "exec", accept_alias)
        panel.rename_item(str(expno))
        assert aliases.alias_for(expno) == "first run" and expno.exists()
        assert exp.text(0).startswith("first run (10)")
        assert got == [(str(expno), str(expno))]

        # (b) on disk, refused while its 1r is open in LARMOR
        warned = []
        monkeypatch.setattr(QMessageBox, "warning",
                            staticmethod(lambda *a, **k: warned.append(a[2])))
        asked = []
        monkeypatch.setattr(QMessageBox, "question",
                            staticmethod(lambda *a, **k: (asked.append(a[2]), QMessageBox.Yes)[1]))

        def accept_disk(self):
            self.edit.setText("11")
            self.rbDisk.setChecked(True)
            return QDialog.Accepted

        monkeypatch.setattr(ex.RenameDialog, "exec", accept_disk)
        panel.open_paths = lambda: [str(expno / "pdata" / "1" / "1r")]
        panel.rename_item(str(expno))
        assert warned and "open in LARMOR" in warned[0] and expno.exists()
        assert not asked                        # refused BEFORE the confirmation

        # (b) on disk, confirmed: moved, alias follows, tree retargeted, signal
        panel.open_paths = lambda: []
        panel.rename_item(str(expno))
        assert asked and str(expno) in asked[0] and str(sample / "11") in asked[0]
        assert (sample / "11" / "acqus").exists() and not expno.exists()
        assert exp.data(0, ex._ROLE_PATH) == str(sample / "11")
        assert exp.data(0, ex._ROLE_OPEN) == str(sample / "11" / "pdata" / "1" / "1r")
        assert exp.text(0).startswith("first run (11) · 27Al 1D")
        assert aliases.alias_for(sample / "11") == "first run"
        assert got[-1] == (str(expno), str(sample / "11"))
        rec = json.loads(aliases.rename_log_path().read_text(encoding="utf-8").splitlines()[-1])
        assert Path(rec["new"]) == sample / "11"

        # cancelling the confirmation moves nothing
        monkeypatch.setattr(QMessageBox, "question",
                            staticmethod(lambda *a, **k: QMessageBox.Cancel))

        def accept_disk_12(self):
            self.edit.setText("12")
            self.rbDisk.setChecked(True)
            return QDialog.Accepted

        monkeypatch.setattr(ex.RenameDialog, "exec", accept_disk_12)
        panel.rename_item(str(sample / "11"))
        assert (sample / "11").exists() and not (sample / "12").exists()
    finally:
        panel.close()


def test_main_window_relabels_open_documents(qapp, store, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    expno = _fake_expno(store)
    sample = expno.parent
    r1 = str(expno / "pdata" / "1" / "1r")
    win = MainWindow()
    try:
        assert win.explorer.open_paths == win._open_source_paths
        win.source_path = r1
        win.recipe = {"sample": "Base0Ca", "nucleus": "27Al", "sites": [],
                      "larmor_frequency_MHz": 130.3}
        win.setWindowTitle("LARMOR — 1r")
        assert r1 in win._open_source_paths()

        # an alias set in the Explorer relabels the title bar and the recipe
        aliases.set_alias(sample, "Glass A")
        win.explorer.renamed.emit(str(sample), str(sample))
        assert win.windowTitle() == "LARMOR — Glass A · 10"
        assert win.recipe["sample"] == "Glass A"
        assert "display name set" in win.statusBar().currentMessage()

        # a rename on disk retargets the open path (the Explorer refuses it
        # while the file is open; the relabel hook itself must still cope)
        new = sample.with_name("01192026_SR31649_GlassA_SS_ALP")
        os.rename(sample, new)                # a bare move: no alias follows
        win.explorer.renamed.emit(str(sample), str(new))
        assert win.source_path == str(new / "10" / "pdata" / "1" / "1r")
        assert win.windowTitle() == "LARMOR — 1r"
        assert win.recipe["sample"] == "GlassA"
        assert "renamed on disk" in win.statusBar().currentMessage()
    finally:
        win.close()
