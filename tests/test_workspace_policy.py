"""An opened spectrum stays open, and a saved fit can be put on the open data.

Sam: "Why I cannot open more than 4-5 spectra, I double clicked them but
they are not staying open in the workspace". The Workspaces dock reused the
active 1D workspace whenever it carried no fit, so browsing a sample folder
replaced the previous spectrum each time and only fitted ones accumulated.
"""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    monkeypatch.setenv("LARMOR_NO_KERNEL_WARM", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _csv(path, centre):
    ppm = np.linspace(-60.0, 100.0, 256)
    amp = np.exp(-4 * np.log(2) * ((ppm - centre) / 4.0) ** 2)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# LARMOR spectrum\n# nucleus=11B\n# larmor_MHz=160.46\n"
                "# spin_rate_Hz=20000.0\nppm,intensity\n")
        for x, y in zip(ppm, amp):
            f.write(f"{x},{y}\n")
    return path


def test_every_opened_spectrum_is_its_own_workspace(win, qapp, tmp_path):
    paths = [_csv(tmp_path / f"s{k}.csv", 5.0 * k) for k in range(6)]
    for p in paths:
        win.load_source(str(p), keep_fit=False)
        qapp.processEvents()
    assert len(win.workspaces) == 6
    assert [ws["kind"] for ws in win.workspaces] == ["1d"] * 6
    assert win.active_ws == 5
    assert not any(ws["has_fit"] for ws in win.workspaces)
    # each keeps its own data
    win.switch_workspace(2)
    qapp.processEvents()
    assert win.source_path == str(paths[2])
    assert float(win.exp_ppm[int(np.argmax(win.exp_amp))]) == pytest.approx(10.0, abs=0.5)
    # Reset to original reloads IN PLACE (an explicit reuse), never a 7th --
    # and the processing record goes with the processing: the recipe used to
    # keep listing the steps while the data was back to raw, so a saved fit
    # replayed a chain the user had just dropped (distribution check)
    win.apply_processing([{"op": "phase", "p0": 15.0, "p1": 0.0}], False)
    qapp.processEvents()
    assert [o["op"] for o in win.recipe["processing"]][-1] == "phase"
    assert win.proc_panel.phase_values() == (15.0, 0.0)
    win.reset_processing()
    qapp.processEvents()
    assert len(win.workspaces) == 6 and win.active_ws == 2
    assert win.recipe["processing"] == [] and not win.recipe.get("processing_from_raw")
    assert win.proc_panel.phase_values() == (0.0, 0.0)


def test_file_menu_opens_a_fit_on_the_current_spectrum(win, qapp, tmp_path, monkeypatch):
    from PySide6.QtGui import QAction

    acts = {a.text(): a for a in win.findChildren(QAction)}
    act = acts.get("Open a f&it on this spectrum…")
    assert act is not None and act.toolTip()
    called = []
    monkeypatch.setattr(win, "apply_recipe_browse", lambda *a: called.append(1))
    # the action was wired at construction: reconnect to the patched slot
    act.triggered.disconnect()
    act.triggered.connect(win.apply_recipe_browse)
    act.trigger()
    assert called == [1]


def test_explorer_fit_row_offers_apply_to_the_open_spectrum(qapp, tmp_path):
    from larmor.desktop.explorer import ExplorerPanel

    panel = ExplorerPanel()
    try:
        got = []
        panel.apply_fit_requested.connect(lambda p: got.append(p))
        menu = panel.fit_menu(str(tmp_path / "a.recipe.json"), [str(tmp_path / "a.recipe.json")])
        texts = [a.text() for a in menu.actions() if not a.isSeparator()]
        assert texts[:2] == ["Open", "Apply to the open spectrum"]
        next(a for a in menu.actions() if a.text() == "Apply to the open spectrum").trigger()
        assert got == [str(tmp_path / "a.recipe.json")]
    finally:
        panel.close()


def test_main_window_wires_the_explorer_apply_signal(win, qapp):
    """The Explorer's 'Apply to the open spectrum' puts a saved fit's lines on
    the spectrum on screen: the shipped 27Al recipe onto the shipped 27Al 1r."""
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    win.load_source(str(repo / "examples" / "pCABS2-4" / "3616" / "pdata" / "1" / "1r"),
                    keep_fit=False)
    qapp.processEvents()
    assert win.recipe["sites"] == []
    win.explorer.apply_fit_requested.emit(str(repo / "examples" / "pCABS2-4_27Al.recipe.json"))
    qapp.processEvents()
    assert [s["model"] for s in win.recipe["sites"]] == ["czjzek"] * 3
    assert len(win.workspaces) == 1                      # applied, not opened anew
