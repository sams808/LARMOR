"""The Workspaces dock's multi-selection and right-click menu, the inline
fit traces the studio is seeded with, and -- in ONE MainWindow -- the slots
behind the menu: Send to Plotting studio, Rename…, Close, Overlay on the
active spectrum and the Fit parameter table."""
import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QPoint, Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REC_AL = ROOT / "examples" / "pCABS2-4_27Al.recipe.json"

ITEMS = [("⤳", "glassA"), ("∿", "glassB"), ("▦", "map"), ("◫", "fig")]
INFO = [{"kind": "1d", "has_fit": True}, {"kind": "1d", "has_fit": False},
        {"kind": "2d", "has_fit": True}, {"kind": "figure", "has_fit": False}]


@pytest.fixture(scope="module")
def qapp():
    yield QApplication.instance() or QApplication([])


def _panel(qapp, active=0):
    from larmor.desktop.workspaces import WorkspacePanel

    p = WorkspacePanel()
    p.resize(260, 220)
    p.show()
    p.rebuild(ITEMS, active, INFO)
    qapp.processEvents()
    return p


def _click(p, row, mods=Qt.NoModifier):
    it = p.list.item(row)
    QTest.mouseClick(p.list.viewport(), Qt.LeftButton, mods,
                     p.list.visualItemRect(it).center())


def _action(menu, prefix):
    for a in menu.actions():
        if not a.isSeparator() and a.text().startswith(prefix):
            return a
    raise AssertionError(f"no action starting with {prefix!r}: "
                         f"{[a.text() for a in menu.actions()]}")


# ------------------------------------------------------------ bare panel
def test_plain_click_switches_but_an_extending_click_does_not(qapp):
    p = _panel(qapp)
    got = []
    p.switch.connect(got.append)
    assert p.selected_rows() == [0]
    assert p.row_info(2) == {"kind": "2d", "has_fit": True} and p.row_info(9) == {}

    _click(p, 2, Qt.ControlModifier)          # extend: Qt moves current, no switch
    assert got == []
    assert p.selected_rows() == [0, 2] and p.list.currentRow() == 2
    _click(p, 1, Qt.ShiftModifier)            # range: no switch either
    assert got == []
    assert 1 in p.selected_rows()

    # a plain click on the row that is ALREADY Qt's current row still
    # switches (currentRowChanged stays silent then)
    assert p.list.currentRow() == 1
    _click(p, 1)
    assert got == [1]
    _click(p, 3)                              # a session row: the window decides
    assert got == [1, 3]
    p.rebuild(ITEMS, 1, INFO)                 # the window made row 1 active
    _click(p, 1)                              # clicking the active row: nothing
    assert got == [1, 3]

    # keyboard: Shift+Down extends without switching, a plain Down switches
    p.list.setFocus()
    QTest.keyClick(p.list, Qt.Key_Down, Qt.ShiftModifier)
    assert got == [1, 3] and p.list.currentRow() == 2
    QTest.keyClick(p.list, Qt.Key_Down)
    assert got == [1, 3, 3]


def test_menu_actions_emit_the_right_indices(qapp):
    p = _panel(qapp)
    caught = {}
    for name in ("activate", "rename", "save", "close_many", "send_to_studio",
                 "fit_table", "overlay"):
        getattr(p, name).connect(lambda v, n=name: caught.setdefault(n, []).append(v))

    m = p.build_menu(1, [0, 1, 2, 3])
    texts = [a.text() for a in m.actions() if not a.isSeparator()]
    assert texts[:5] == ["Switch to", "Rename…", "Save", "Close 4 selected",
                         "Close others"]
    assert _action(m, "Switch to").isEnabled()
    assert not _action(m, "Close others").isEnabled()
    assert _action(m, "Send to").text() == "Send to Plotting studio  (2 spectra)"
    assert _action(m, "Fit parameter").text() == "Fit parameter table…  (2 fits)"
    assert _action(m, "Overlay").text() == "Overlay on the active spectrum"
    _action(m, "Send to").trigger()
    assert caught["send_to_studio"] == [[0, 1]]          # the 1D rows only
    _action(m, "Fit parameter").trigger()
    assert caught["fit_table"] == [[0, 2]]               # the fitted docs only
    _action(m, "Overlay").trigger()
    assert caught["overlay"] == [[1]]                    # the active row excluded
    _action(m, "Close 4").trigger()
    assert caught["close_many"] == [[0, 1, 2, 3]]
    _action(m, "Rename").trigger()
    _action(m, "Save").trigger()
    _action(m, "Switch to").trigger()
    assert (caught["rename"], caught["save"], caught["activate"]) == ([1], [1], [1])

    # a lone figure row: Open, nothing to plot / table / overlay
    m2 = p.build_menu(3, [3])
    assert [a.text() for a in m2.actions() if not a.isSeparator()][0] == "Open"
    for prefix in ("Send to", "Fit parameter", "Overlay"):
        assert not _action(m2, prefix).isEnabled()
    assert _action(m2, "Close").text() == "Close"
    _action(m2, "Close others").trigger()
    assert caught["close_many"][-1] == [0, 1, 2]

    # the active row: cannot switch to itself nor overlay itself; one fit
    m3 = p.build_menu(0, [0])
    assert not _action(m3, "Switch to").isEnabled()
    assert not _action(m3, "Overlay").isEnabled()
    assert _action(m3, "Fit parameter").isEnabled()
    assert _action(m3, "Send to").text() == "Send to Plotting studio"
    for a in m3.actions():
        if not a.isSeparator():
            assert a.toolTip() or a.text() in ("Close", "Save")


def test_right_click_outside_the_selection_acts_on_that_row(qapp):
    p = _panel(qapp)
    _click(p, 2, Qt.ControlModifier)
    assert p.selected_rows() == [0, 2]
    pos = p.list.visualItemRect(p.list.item(1)).center()
    m = p.menu_at(pos)                        # row 1 is not selected: it alone
    assert p.selected_rows() == [1]
    assert _action(m, "Close").text() == "Close"
    assert _action(m, "Close others").isEnabled()
    _click(p, 2, Qt.ControlModifier)          # a click inside the selection
    m2 = p.menu_at(pos)                       # keeps it
    assert p.selected_rows() == [1, 2]
    assert _action(m2, "Close").text() == "Close 2 selected"
    assert p.menu_at(QPoint(-5, -5)) is None  # off any row: no menu


# ---------------------------------------------------------- fit traces
def _recipe(sample, pos, model="gauss_lor"):
    return {"nucleus": "27Al", "larmor_frequency_MHz": 130.3, "sample": sample,
            "sites": [{"model": model, "label": "A", "params": {
                "isotropic_chemical_shift_ppm": {"value": pos},
                "shift_fwhm_ppm": {"value": 8},
                "gl": {"value": 0.5},
                "amplitude": {"value": 1.0}}}]}


def test_fit_traces_inline_data_styles_and_failure_note():
    from larmor import figures

    x = np.linspace(-40, 120, 300)
    y = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    traces, note = figures.fit_traces("g", x, y, _recipe("g", 60.0), color="#123456")
    assert note == ""
    assert [t["label"] for t in traces] == ["g · experiment", "g · fit", "g · A"]
    assert all(t["color"] == "#123456" for t in traces)
    assert traces[1]["linestyle"] == "--" and "linestyle" not in traces[0]
    assert traces[2]["linewidth"] == 0.8 and traces[2]["alpha"] == 0.7
    assert len(traces[1]["data"]["x"]) == 300
    peak = int(np.argmax(traces[1]["data"]["y"]))
    assert traces[1]["data"]["x"][peak] == pytest.approx(60.0, abs=1.0)
    # no sites: the experiment alone; a model that cannot simulate: a note
    only, note = figures.fit_traces("g", x, y, None)
    assert len(only) == 1 and note == ""
    bad, note = figures.fit_traces("g", x, y, _recipe("g", 60.0, model="no_such_model"))
    assert len(bad) == 1 and note.startswith("g: ")
    # the spec the studio receives is drawable as it is
    import matplotlib.pyplot as plt
    fig = figures.render_1d({"kind": "1d", "traces": traces})
    assert len(fig.axes[0].lines) == 3
    plt.close(fig)


# ---------------------------------------------------------- main window
def _mkws(w, sample, pos):
    """Register a fitted 1D workspace (one gauss_lor line) on the window."""
    x = np.linspace(-40, 120, 300)
    y = np.exp(-0.5 * ((x - pos) / 8) ** 2)
    w.exp_ppm, w.exp_amp = x, y
    w.recipe = _recipe(sample, pos)
    w.source_path = "src_" + sample
    w.hidden = set()
    w.central_stack.setCurrentWidget(w.view)
    w.view.set_experiment(x, y)
    w.lines_table.rebuild(w.recipe, w.hidden)
    w._ws_mode = "new"
    w._register_ws("1d")


def test_main_window_slots_behind_the_menu(qapp, monkeypatch):
    """The one MainWindow of this file: studio seeding, rename persistence,
    fit entries / table, overlay and close-many."""
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    try:
        _mkws(w, "glassA", 60.0)
        _mkws(w, "glassB", 20.0)
        assert [ws["title"] for ws in w.workspaces] == ["glassA", "glassB"]
        assert w.active_ws == 1
        assert w.ws_panel.row_info(0) == {"kind": "1d", "has_fit": True}

        # ---- Send to Plotting studio: the spec, captured instead of shown
        captured = []
        monkeypatch.setattr(w, "open_plotting_studio",
                            lambda spec=None, **k: captured.append(spec))
        w.send_workspaces_to_studio([0, 1])
        spec = captured[-1]
        assert spec["kind"] == "1d" and spec["title"] == "glassA + 1 more"
        assert [t["label"] for t in spec["traces"]] == [
            "glassA · experiment", "glassA · fit", "glassA · A",
            "glassB · experiment", "glassB · fit", "glassB · A"]
        colours = [t["color"] for t in spec["traces"]]
        assert colours[0] == colours[1] == colours[2]
        assert colours[3] == colours[4] == colours[5] and colours[0] != colours[3]
        assert spec["traces"][1]["linestyle"] == "--"
        assert spec["traces"][2]["alpha"] == 0.7
        assert spec["traces"][3]["data"]["y"] == pytest.approx(
            list(w.workspaces[1]["snap"]["exp_amp"]))
        assert "2 spectrum(s)" in w.statusBar().currentMessage()
        # one workspace: its own title; session rows are never sent
        w.send_workspaces_to_studio([0])
        assert captured[-1]["title"] == "glassA" and len(captured[-1]["traces"]) == 3
        w._add_session_entry("figure", "fig", {"spec": spec})
        w.send_workspaces_to_studio([2])
        assert len(captured) == 2
        assert "select at least one" in w.statusBar().currentMessage()
        # a model that cannot be simulated: the experiment goes, the bar says why
        w.workspaces[0]["snap"]["recipe"]["sites"][0]["model"] = "no_such_model"
        w.send_workspaces_to_studio([0])
        assert len(captured[-1]["traces"]) == 1
        assert "could not simulate" in w.statusBar().currentMessage()
        w.workspaces[0]["snap"]["recipe"]["sites"][0]["model"] = "gauss_lor"

        # ---- Rename… survives switches (custom_title), '' restores
        w.set_workspace_title(0, "my glass")
        assert w.workspaces[0]["title"] == "my glass"
        w.switch_workspace(0)
        w.switch_workspace(1)
        assert w.workspaces[0]["title"] == "my glass"
        assert w.ws_panel.list.item(0).text().endswith("my glass")
        assert [e.name for e in w.workspace_fit_entries([0])] == ["my glass"]
        w.set_workspace_title(0, "")
        assert w.workspaces[0]["title"] == "glassA"
        assert "custom_title" not in w.workspaces[0]

        # ---- the fit entries and the table (tool window captured, not shown)
        entries = w.workspace_fit_entries()
        assert [e.name for e in entries] == ["glassA", "glassB"]
        assert entries[0].source == "src_glassA" and entries[0].nucleus == "27Al"
        shown = []
        monkeypatch.setattr("larmor.desktop.mw_tools.show_tool_window",
                            lambda d, *a, **k: shown.append(d) or d)
        w.open_fit_table(False)                   # a QAction's bool: every fit
        assert shown[-1].table.rowCount() == 2
        assert shown[-1].windowTitle() == "Fit parameter table — 2 fits"
        w.open_fit_table([1])
        assert shown[-1].table.rowCount() == 1
        w.open_fit_table_files([str(REC_AL), str(ROOT / "nope.fxml")])
        assert shown[-1].table.rowCount() == 3
        assert "could not read" in w.statusBar().currentMessage()
        shown[-1]._add_workspaces()               # the dialog asks the window
        assert shown[-1].table.rowCount() == 5

        # ---- Overlay on the active spectrum (active = glassB): glassA only
        w.overlay_workspaces([0, 1])
        assert [o["label"] for o in w._overlays] == ["glassA"]
        assert w._overlays[0]["nucleus"] == "27Al" and w._overlays[0]["npts"] == 300

        # ---- Close (selection) and Close others
        w.close_workspaces([2, 0])                # the figure row and glassA
        assert [ws["title"] for ws in w.workspaces] == ["glassB"]
        assert w.active_ws == 0 and w.ws_panel.list.count() == 1
    finally:
        w.close()
