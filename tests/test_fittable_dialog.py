"""The Fit parameter table dialog renders larmor.fittable: wide / long
layouts, Copy (TSV with headers), Export CSV (round trip), Add fits… / Add
open workspaces / Remove selected, and header sorting that moves whole rows."""
import csv
import os
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from larmor import fittable as F  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REC_AL = ROOT / "examples" / "pCABS2-4_27Al.recipe.json"
REC_B = ROOT / "examples" / "pCABS2-4_11B.recipe.json"
FXML_AL = ROOT / "examples" / "pCABS2-4" / "3616" / "pdata" / "1" / "1r.fxml"


@pytest.fixture(scope="module")
def qapp():
    yield QApplication.instance() or QApplication([])


def _dlg(qapp, entries=None, workspaces=None):
    from larmor.desktop.fittable_dialog import FitTableDialog

    return FitTableDialog(None, entries, workspaces=workspaces)


def _ws_entry(name="ws1"):
    return F.entry_from_recipe(name, {
        "nucleus": "27Al", "larmor_frequency_MHz": 130.3, "fit_rmsd": 0.01,
        "sites": [{"model": "gauss_lor", "label": "A", "params": {
            "isotropic_chemical_shift_ppm": {"value": 60.0, "stderr": 0.2},
            "shift_fwhm_ppm": {"value": 8.0}, "gl": {"value": 0.5, "vary": False},
            "amplitude": {"value": 1.0}}}]}, source="s1")


def test_builds_the_wide_table_in_arrival_order(qapp):
    d = _dlg(qapp, [F.load_fit_file(REC_AL), F.load_fit_file(REC_B)])
    headers, rows = d.current_table()
    exp_h, exp_rows = F.build_wide(d.entries)
    assert headers == exp_h and rows == exp_rows          # untouched, unsorted
    assert d.table.rowCount() == 6
    assert d.table.horizontalHeader().sortIndicatorSection() == -1
    assert d.windowTitle() == "Fit parameter table — 2 fits"
    assert "† fixed" in d.status.text()
    assert d.table.item(3, 0).toolTip() == str(REC_B)     # the row's fit file
    assert not d.btnAddWs.isEnabled()                    # no window to ask


def test_empty_dialog_hints_and_disables(qapp):
    d = _dlg(qapp)
    assert d.windowTitle() == "Fit parameter table — no fits yet"
    assert "Add fits" in d.status.text()
    assert not (d.btnCopy.isEnabled() or d.btnCsv.isEnabled() or d.btnRemove.isEnabled())
    d.add_entries([_ws_entry()])
    assert d.windowTitle() == "Fit parameter table — 1 fit"
    assert d.btnCopy.isEnabled() and d.btnCsv.isEnabled()


def test_copy_has_headers_and_csv_round_trips_both_layouts(qapp, tmp_path):
    d = _dlg(qapp, [F.load_fit_file(REC_AL)])
    text = d.copy_tsv()
    lines = text.split("\n")
    headers, rows = d.current_table()
    assert lines[0].split("\t") == headers and len(lines) == 4
    assert lines[1].split("\t") == rows[0]
    assert QApplication.clipboard().text() == text
    assert "copied 3 row(s)" in d.status.text()
    out = tmp_path / "wide.csv"
    assert d.export_csv(out)
    with open(out, encoding="utf-8", newline="") as f:
        back = list(csv.reader(f))
    assert back[0] == headers and back[1:] == rows
    assert "wide.csv" in d.status.text()
    # the long layout: one row per parameter, status words, numbers as %.8g
    d.chkLong.setChecked(True)
    h2, r2 = d.current_table()
    assert h2 == list(F.LONG_HEADERS) and len(r2) == 12
    assert r2[0][:4] == ["pCABS2-4_27Al", "AlO$_4$", "czjzek", "isotropic_chemical_shift_ppm"]
    assert r2[0][4] == "62.815324" and r2[0][6] == ""
    assert d.export_csv(tmp_path / "long")               # .csv appended
    with open(tmp_path / "long.csv", encoding="utf-8", newline="") as f:
        back2 = list(csv.reader(f))
    assert back2[0] == h2 and back2[1:] == r2
    assert d.copy_tsv().split("\n")[0].split("\t") == h2
    d.chkLong.setChecked(False)
    assert d.current_table()[0] == headers


def test_add_workspaces_files_and_remove_selected(qapp, tmp_path):
    asked = []

    def provider():
        asked.append(1)
        return [_ws_entry()]

    d = _dlg(qapp, workspaces=provider)
    assert d.btnAddWs.isEnabled()
    d._add_workspaces()
    assert len(d.entries) == 1 and asked == [1]
    d._add_workspaces()                                   # already listed: skipped
    assert len(d.entries) == 1 and "already" in d.status.text()
    n = d.add_files([str(REC_B), str(tmp_path / "missing.recipe.json")])
    assert n == 1 and "could not read" in d.status.text()
    assert [e.name for e in d.entries] == ["ws1", "pCABS2-4_11B"]
    assert d.table.rowCount() == 4
    # Remove selected maps a row back to its fit
    d.table.selectRow(2)
    assert d.selected_entry_indices() == {1}
    d._remove_selected()
    assert [e.name for e in d.entries] == ["ws1"] and d.table.rowCount() == 1
    d.table.clearSelection()
    d._remove_selected()                                  # nothing selected: a hint
    assert len(d.entries) == 1 and "select a row" in d.status.text()


def test_sorting_keeps_rows_intact_and_outputs_follow_the_view(qapp):
    entries = [F.load_fit_file(REC_AL), F.load_fit_file(REC_B), F.load_fit_file(FXML_AL)]
    d = _dlg(qapp, entries)
    headers, before = d.current_table()
    pos = headers.index("Position (ppm)")
    d.table.sortItems(pos, Qt.DescendingOrder)
    _, after = d.current_table()
    assert sorted(map(tuple, before)) == sorted(map(tuple, after))   # whole rows
    vals = [float(r[pos].split()[0]) for r in after]
    assert vals == sorted(vals, reverse=True)
    assert after != before
    # numbers sort before blanks (the dmfit fit has no RMSD); Copy follows the view
    rc = headers.index("RMSD")
    d.table.sortItems(rc, Qt.AscendingOrder)
    _, r3 = d.current_table()
    assert [r[rc] for r in r3][:3] == ["0.0037236"] * 3
    assert [r[rc] for r in r3][-3:] == [""] * 3
    assert all(r[0] == "1r" for r in r3[-3:])
    assert d.copy_tsv().split("\n")[1].split("\t") == r3[0]
    # a rebuild keeps the sort; a layout change forgets it
    d.add_entries([_ws_entry()])
    _, r4 = d.current_table()
    assert r4[0][rc] == "0.0037236" and d.table.rowCount() == 10
    d.chkLong.setChecked(True)
    assert d.table.horizontalHeader().sortIndicatorSection() == -1
    assert d.current_table()[1] == [[F.cell_text(c) for c in r]
                                    for r in F.build_long(d.entries)[1]]
    # Remove selected after a sort still removes the right fit
    d.chkLong.setChecked(False)
    d.table.sortItems(rc, Qt.DescendingOrder)
    _, r5 = d.current_table()
    row_1r = next(i for i, r in enumerate(r5) if r[0] == "1r")
    d.table.selectRow(row_1r)
    d._remove_selected()
    assert [e.name for e in d.entries] == ["pCABS2-4_27Al", "pCABS2-4_11B", "ws1"]
