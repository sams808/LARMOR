"""A public build must open and work with whatever the registry holds: every
saved preference LARMOR reads is set to garbage here, and the window must
still build, show its first-run hint, add a line, rebuild its recent menus
and open its tools. The user's real settings are snapshotted and restored
(QSettings is the live registry; see larmor.distcheck.gui_common)."""
import json
import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtCore import QSettings  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from larmor.distcheck.gui_common import SettingsGuard  # noqa: E402

GARBAGE = {
    "fontPt": "abc", "theme": "NoSuchTheme", "appearanceOverride": "nope",
    "recent": "not-a-list", "recentRecipes": 42, "siteDefaults": "{not json",
    "constraintLibrary": "[", "batchTemplates": "x", "pinnedFolders": "Z:\\nowhere",
    "pinnedNames": "::::", "czjzekDisplay": "??", "yAxisMode": "bogus",
    "yAxisRegion": "bogus", "lastDir": "Z:\\nowhere\\at\\all",
    "plottingStudio/dataRoots": ";;;", "qcpmgDialogGeometry": b"junk",
    "fitStdevPct": "abc", "session/source": "Z:\\nope\\1r", "session/recipe": "{bad",
    "ssbOfferOnLoad": "maybe", "welcomeShown": "x",
}


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def garbage_settings():
    # the guard snapshots the registry to memory AND to a .reg file first;
    # only the keys under test are overwritten (never clear(): a process
    # killed mid-test must leave as little damage as possible)
    with SettingsGuard() as guard:
        assert guard.backup_path is None or guard.backup_path.exists()
        s = QSettings("LARMOR", "app")
        for k, v in GARBAGE.items():
            s.setValue(k, v)
        s.sync()
        yield s


def test_prefs_helpers_survive_garbage(qapp, garbage_settings):
    from larmor.desktop import prefs
    from larmor.desktop.app import saved_font_pt

    assert saved_font_pt() == 9
    assert prefs.json_setting("siteDefaults", {}) == {}
    assert prefs.json_setting("constraintLibrary", {}) == {}
    assert prefs.json_setting("batchTemplates", {}) == {}
    assert prefs.list_setting("recent") == ["not-a-list"]        # a lone string is a list
    assert prefs.list_setting("recentRecipes") == []
    assert prefs.str_setting("yAxisMode", "raw", allowed=("raw", "max")) == "raw"
    # and sane values still read back as themselves
    garbage_settings.setValue("siteDefaults", json.dumps({"27Al/czjzek": {"a": 1}}))
    garbage_settings.setValue("recent", ["a", "b"])
    garbage_settings.setValue("fontPt", 11)
    assert prefs.json_setting("siteDefaults", {}) == {"27Al/czjzek": {"a": 1}}
    assert prefs.list_setting("recent") == ["a", "b"]
    assert saved_font_pt() == 11
    garbage_settings.setValue("fontPt", 400)
    assert saved_font_pt() == 9                                   # out of range: default


def test_window_builds_and_works_with_garbage_in_every_saved_setting(
        qapp, garbage_settings, tmp_path, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_KERNEL_WARM", "1")
    from larmor.desktop.app import MainWindow
    from larmor.distcheck.gui_common import ModalDismisser

    # the Plotting studio opens modally (exec): the dismisser rejects every
    # modal within ~100 ms, so nothing here can hang the test run
    dismiss = ModalDismisser()
    dismiss.__enter__()
    win = MainWindow()
    try:
        win._add_recent = lambda p: None
        win._maybe_show_welcome()                 # the first-run hint reads 'recent'
        win._rebuild_recent()
        win._rebuild_apply_recipe()
        # a spectrum, a line (siteDefaults is parsed on Add line), a simulation
        import numpy as np
        x = np.linspace(-20, 60, 300)
        y = 100 * np.exp(-0.5 * ((x - 15) / 4) ** 2)
        p = tmp_path / "s.csv"
        with open(p, "w", encoding="utf-8") as f:
            f.write("# nucleus = 27Al\n# larmor_MHz = 130.3\n")
            for xi, yi in zip(x, y):
                f.write(f"{xi:.4f} {yi:.4f}\n")
        win.load_source(str(p), keep_fit=False)
        win._model_actions["czjzek"].setChecked(True)
        win.add_site_at(15.0, 100.0)
        assert len(win.recipe["sites"]) == 1
        win._model_actions["gauss_lor"].setChecked(True)
        win.add_site_at(5.0, 10.0)
        assert len(win.recipe["sites"]) == 2
        # menus and tools that parse other saved values
        win.open_plotting_studio()                # plottingStudio/dataRoots
        win.explorer._pinned                      # pinnedFolders was read at build time
        from larmor.desktop.workers import _fit_tol
        assert _fit_tol() == 0.1                  # fitStdevPct garbage -> the default
        for w in list(QApplication.topLevelWidgets()):
            if w is not win and w.isVisible():
                w.close()
        qapp.processEvents()
        assert any("Plotting studio" in t for t in dismiss.seen), dismiss.seen
    finally:
        dismiss.__exit__(None, None, None)
        win.close()
        qapp.processEvents()


def test_corrupt_session_file_is_ignored(qapp, tmp_path, monkeypatch):
    """The crash-recovery session lives in %LOCALAPPDATA%/LARMOR/session.json;
    garbage there (or a source that no longer exists) must not stop start-up."""
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    monkeypatch.setenv("LARMOR_NO_KERNEL_WARM", "1")
    from larmor.desktop.app import MainWindow

    d = tmp_path / "LARMOR"
    d.mkdir()
    for payload in ("{not json", json.dumps({"source": "Z:\\nope\\1r", "recipe": "bad"}),
                    json.dumps({"source": 12, "recipe": {"sites": "x"}})):
        (d / "session.json").write_text(payload, encoding="utf-8")
        with SettingsGuard():
            win = MainWindow()
            try:
                assert win.exp_ppm is None or len(win.exp_ppm) == 0
            finally:
                win.close()
                qapp.processEvents()
