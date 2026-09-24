"""The Fit-parameters table's line menu and multi-line removal (wip/WB).

Right-click on a row: Add spinning sidebands… (the Decomposition dialog with
the line preselected), Duplicate, Rename…, Hide/Show, Fix all / Free all,
Move, Remove -- and "Remove N selected lines" when several rows are selected.
Delete on a multi-selection removes them all in ONE undo step with every
constraint remapped (a link to a removed line is dropped, one between the
survivors renumbered -- never left to become a self-reference).
"""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication, QInputDialog  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _recipe_dict(n_sites=4, model="gauss_lor"):
    from larmor.recipe import Param, Recipe, SiteModel

    sites = [SiteModel(model=model, label=f"L{i}", params={
        "isotropic_chemical_shift_ppm": Param(10.0 * i),
        "shift_fwhm_ppm": Param(4.0), "amplitude": Param(1.0),
        "gl": Param(0.5, vary=False)}) for i in range(n_sites)]
    return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3,
                  spin_rate_Hz=13030.0, sites=sites).to_dict()


def _texts(menu):
    return [a.text() for a in menu.actions() if not a.isSeparator()]


def _select_rows(t, rows):
    from PySide6.QtCore import QItemSelectionModel as M

    t.table.clearSelection()
    for r in rows:
        t.table.selectionModel().select(t.table.model().index(r, 0),
                                        M.SelectionFlag.Select | M.SelectionFlag.Rows)
    t.table.setCurrentCell(rows[-1], 0, M.SelectionFlag.NoUpdate)


# ------------------------------------------------------------- the table
def test_line_menu_lists_every_line_action_and_keeps_the_cell_menus(qapp):
    from larmor.desktop.table import LinesTable

    t = LinesTable()
    rec = _recipe_dict(3)
    t.rebuild(rec, {2})
    menu = t._build_menu(0, None)
    texts = _texts(menu)
    assert menu.actions()[0].menu().title() == "Family"       # unchanged, first
    assert texts[1].startswith("Add spinning sidebands…")   # the section is a separator
    for want in ("Duplicate line", "Rename line…", "Hide on plot  (Ctrl+H toggles)",
                 "Fix all parameters of the line", "Free all parameters of the line",
                 "Move line up", "Move line down", "Remove line"):
        assert want in texts, (want, texts)
    assert "Change model" not in " ".join(texts)              # no model-change path
    # a hidden line offers Show; the first line cannot move up
    assert "Show on plot  (Ctrl+H toggles)" in _texts(t._build_menu(2, None))
    assert not next(a for a in menu.actions() if a.text() == "Move line up").isEnabled()
    # the cell menus are still there on top of the line actions
    cell_menu = t._build_menu(1, "shift_fwhm_ppm")
    ct = _texts(cell_menu)
    assert ct[0] == "Width: same as another line…" and "Constrain min / max…" in ct
    assert "Remove line" in ct and ct.index("Constrain min / max…") < ct.index("Remove line")
    # a self-sidebanding model greys the sidebands entry
    t.rebuild(_recipe_dict(1, model="csa_mas") | {"sites": [
        {"model": "csa_mas", "label": "c", "params": rec["sites"][0]["params"]}]}, set())
    a = next(a for a in t._build_menu(0, None).actions()
             if a.text().startswith("Add spinning sidebands"))
    assert not a.isEnabled()
    t.close()


def test_line_menu_actions_emit_and_edit(qapp, monkeypatch):
    from larmor.desktop.table import LinesTable

    t = LinesTable()
    rec = _recipe_dict(3)
    t.rebuild(rec, set())
    got, edits, ssb = [], [], []
    t.structure.connect(lambda r, a: got.append((r, a)))
    t.edited.connect(lambda: edits.append(1))
    t.sidebands_requested.connect(ssb.append)
    menu = t._build_menu(1, None)
    by = {a.text(): a for a in menu.actions() if not a.isSeparator()}
    next(a for k, a in by.items() if k.startswith("Add spinning sidebands")).trigger()
    assert ssb == [1]
    by["Duplicate line"].trigger()
    by["Remove line"].trigger()
    assert got == [(1, "duplicate"), (1, "remove")]
    # Fix all pins every non-linked parameter and refreshes the pin boxes in place
    rec["sites"][1]["params"]["amplitude"]["expr"] = "0.5 * s0.amplitude"
    rec["sites"][1]["params"]["amplitude"]["vary"] = True
    t.rebuild(rec, set())
    t._set_all_vary(1, False)
    ps = rec["sites"][1]["params"]
    assert not ps["isotropic_chemical_shift_ppm"]["vary"] and not ps["shift_fwhm_ppm"]["vary"]
    assert ps["amplitude"]["vary"] and ps["amplitude"]["expr"]      # linked: untouched
    col = 2 + t._used_keys.index("shift_fwhm_ppm")
    cell = t.table.cellWidget(1, col)
    assert cell.pin.isChecked() and cell.is_held()
    assert "italic" in cell.edit.styleSheet() and "held" in cell.edit.toolTip()
    assert len(edits) == 1
    menu = t._build_menu(1, None)
    by = {a.text(): a for a in menu.actions() if not a.isSeparator()}
    assert not by["Fix all parameters of the line"].isEnabled()
    by["Free all parameters of the line"].trigger()
    assert all(p["vary"] for p in ps.values())
    assert not cell.pin.isChecked() and "italic" not in cell.edit.styleSheet()
    assert len(edits) == 2
    # Rename… writes the label and the letter cell, and asks for a redraw
    monkeypatch.setattr(QInputDialog, "getText",
                        staticmethod(lambda *a, **k: ("AlO4", True)))
    t._rename(1)
    assert rec["sites"][1]["label"] == "AlO4"
    assert t.table.item(1, 0).text() == "■ B · AlO4"
    assert len(edits) == 3
    monkeypatch.setattr(QInputDialog, "getText",
                        staticmethod(lambda *a, **k: ("nope", False)))
    t._rename(1)
    assert rec["sites"][1]["label"] == "AlO4" and len(edits) == 3
    t.close()


def test_multi_selection_drives_remove_from_menu_and_delete_key(qapp):
    from larmor.desktop.table import LinesTable

    t = LinesTable()
    t.rebuild(_recipe_dict(4), set())
    many, one = [], []
    t.remove_lines.connect(lambda rows: many.append(list(rows)))
    t.structure.connect(lambda r, a: one.append((r, a)))
    _select_rows(t, [0, 2])
    assert t.selected_rows() == [0, 2]
    assert t.selected_rows(3) == [0, 2, 3]                 # the clicked row joins
    menu = t._build_menu(2, None)
    a_del = next(a for a in menu.actions() if a.text().startswith("Remove "))
    assert a_del.text() == "Remove 2 selected lines  (A, C)"
    a_del.trigger()
    assert many == [[0, 2]] and one == []
    QTest.keyClick(t.table, Qt.Key_Delete)
    assert many == [[0, 2], [0, 2]] and one == []
    # a single selection keeps the structure path every embedder handles
    t.table.clearSelection()
    t.table.setCurrentCell(3, 0)
    assert t.selected_rows() == [3]
    assert _texts(t._build_menu(3, None)).count("Remove line") == 1
    QTest.keyClick(t.table, Qt.Key_Delete)
    assert one == [(3, "remove")] and len(many) == 2
    t.close()


# ------------------------------------------------------------ the window
@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _window_with_lines(win, n=5):
    x = np.linspace(-100.0, 100.0, 801)
    y = np.exp(-0.5 * ((x - 10.0) / 3.0) ** 2)
    win._display_1d(x, y, "27Al", 130.3, 13030.0, "t", "x")
    win.recipe["sites"] = _recipe_dict(n)["sites"]
    win.recipe["spin_rate_Hz"] = 13030.0
    win.on_structure_changed()
    return win.recipe["sites"]


def test_remove_sites_takes_every_selected_line_in_one_undo_step(qapp, win):
    sites = _window_with_lines(win, 5)
    # links: B ← A (dropped with A), D ← C (renumbered), E ← B (dropped with B)
    sites[1]["params"]["shift_fwhm_ppm"]["expr"] = "s0.shift_fwhm_ppm"
    sites[3]["params"]["isotropic_chemical_shift_ppm"]["expr"] = \
        "s2.isotropic_chemical_shift_ppm + 5.3"
    sites[4]["params"]["amplitude"]["expr"] = "0.5 * s1.amplitude"
    win.hidden = {1, 4}
    n_undo = len(win.undo_stack)
    win.remove_sites([1, 0, 1, 99, -1])                  # duplicates and junk ignored
    labels = [s["label"] for s in win.recipe["sites"]]
    assert labels == ["L2", "L3", "L4"]
    assert len(win.undo_stack) == n_undo + 1               # one snapshot
    # C↔D survives, renumbered; every link to a removed line is gone
    assert win.recipe["sites"][1]["params"]["isotropic_chemical_shift_ppm"]["expr"] == \
        "s0.isotropic_chemical_shift_ppm + 5.3"
    assert win.recipe["sites"][2]["params"]["amplitude"]["expr"] is None
    assert win.hidden == {2}                               # E followed its line
    from larmor.constraints_util import sanitize_constraints
    assert sanitize_constraints(win.recipe["sites"]) == []
    assert "removed 2 lines" in win.statusBar().currentMessage()
    assert "dropped now-invalid constraint(s)" in win.statusBar().currentMessage()
    assert win.lines_table.table.rowCount() == 3
    win.undo()
    assert [s["label"] for s in win.recipe["sites"]] == ["L0", "L1", "L2", "L3", "L4"]
    assert win.recipe["sites"][4]["params"]["amplitude"]["expr"] == "0.5 * s1.amplitude"
    # the table's multi-selection reaches it through the wired signal
    win.lines_table.remove_lines.emit([3, 4])
    assert [s["label"] for s in win.recipe["sites"]] == ["L0", "L1", "L2"]
    # nothing to do on an empty / all-invalid request
    n_undo = len(win.undo_stack)
    win.remove_sites([])
    win.remove_sites([7])
    assert len(win.undo_stack) == n_undo


def test_line_menu_sidebands_preselects_the_line(qapp, win, monkeypatch):
    from PySide6.QtWidgets import QDialog

    _window_with_lines(win, 3)
    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Accepted)
    win.lines_table.sidebands_requested.emit(2)              # the menu's signal
    sites = win.recipe["sites"]
    assert len(sites) == 5
    plus = next(s for s in sites if s["label"] == "L2+1sb")
    assert plus["params"]["isotropic_chemical_shift_ppm"]["expr"].startswith(
        "s2.isotropic_chemical_shift_ppm + ")
    assert plus["params"]["shift_fwhm_ppm"]["expr"] == "s2.shift_fwhm_ppm"
    assert win._ssb_preselect is None                       # consumed
    # the plain menu entry still defaults to the first eligible line
    win.add_sidebands()
    assert any(s["label"] == "L0+1sb" for s in win.recipe["sites"])


def test_cofit_tables_honour_the_multi_selection_removal(qapp, win):
    """The co-fit tables (two recipes kept identical) take a multi-selection
    Remove through their per-row path, highest row first; the sidebands
    entry only explains where sidebands are added."""
    import json

    _window_with_lines(win, 4)
    r = json.loads(json.dumps(win.recipe))
    win._cofit = {"d1": None, "d2": None, "home": None, "tie": set(),
                  "r1": json.loads(json.dumps(r)), "r2": json.loads(json.dumps(r))}
    win._cofit_rebuild_tables()
    win.cofit_table1d.remove_lines.emit([0, 2])
    for key in ("r1", "r2"):
        assert [s["label"] for s in win._cofit[key]["sites"]] == ["L1", "L3"]
    assert win.cofit_table1d.table.rowCount() == 2
    assert win.cofit_table2d.table.rowCount() == 2
    win.cofit_table2d.sidebands_requested.emit(0)
    assert "co-fit tables keep both recipes identical" in win.statusBar().currentMessage()
    assert len(win._cofit["r1"]["sites"]) == 2            # nothing added
    win._cofit = None
