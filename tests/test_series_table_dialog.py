"""SeriesTableDialog / JoinMappingDialog offscreen: rename with uniqueness,
move + perm, sort, typed columns, the CSV join through the confirm dialog,
save / load of a .series.json."""
import csv
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtGui import QColor  # noqa: E402
from PySide6.QtWidgets import QApplication, QDialog, QFileDialog, QInputDialog  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _dialog(names=None, locked=False):
    from larmor.desktop.series_table_dialog import SeriesTableDialog
    from larmor.series_table import SeriesTable
    t = SeriesTable.from_paths(["C:/d/a.csv", "C:/d/b.csv", "C:/d/c.csv"],
                               names=names or ["P5-Bi8-12", "P5-Bi0", "Base0Ca"])
    return SeriesTableDialog(None, t, locked=locked)


def test_series_table_dialog_renames_moves_sorts_and_reports_perm(qapp, monkeypatch):
    dlg = _dialog()
    assert dlg.grid.rowCount() == 3 and dlg.grid.columnCount() == 5
    assert [dlg.grid.horizontalHeaderItem(i).text() for i in range(5)] == \
        ["#", "name", "group", "folder", "title"]
    assert dlg.perm() == [0, 1, 2] and "order: as loaded" in dlg.status.text()
    dlg.grid.item(0, 1).setText("P5Bi8-12")                 # edit the name cell
    assert dlg.table().labels()[0] == "P5Bi8-12"
    assert dlg.table().rows[0].group == "P5Bi8-12"          # a singleton's group follows
    dlg.grid.item(1, 1).setText("Base0Ca")                  # collides -> (2)
    assert dlg.grid.item(1, 1).text() == "Base0Ca (2)"
    assert dlg.table().labels() == ["P5Bi8-12", "Base0Ca (2)", "Base0Ca"]
    dlg.grid.item(1, 2).setText("Base0Ca")                  # make it a replicate
    assert dlg.table().rows[1].group == "Base0Ca"
    assert "1 replicate group(s): Base0Ca ×2" in dlg.status.text()
    dlg.grid.selectRow(2); dlg.grid.setCurrentCell(2, 1)
    dlg._move(-1)
    assert dlg.perm() == [0, 2, 1]
    assert [dlg.grid.item(r, 0).text() for r in range(3)] == ["1", "2", "3"]
    assert [dlg.grid.item(r, 1).text() for r in range(3)] == ["P5Bi8-12", "Base0Ca", "Base0Ca (2)"]
    assert "order: custom" in dlg.status.text()
    dlg._move(-1)                                           # row 1 up
    assert dlg.perm() == [2, 0, 1]
    dlg._sort_by(None)
    assert dlg.table().labels() == ["Base0Ca", "Base0Ca (2)", "P5Bi8-12"]
    assert dlg.perm() == [2, 1, 0]
    answers = iter([("name (natural)", True), ("descending", True)])
    monkeypatch.setattr(QInputDialog, "getItem", staticmethod(lambda *a, **k: next(answers)))
    dlg._ask_sort()
    assert dlg.table().labels() == ["P5Bi8-12", "Base0Ca (2)", "Base0Ca"]
    assert dlg.perm() == [0, 1, 2]
    dlg._reset_names()
    assert dlg.table().labels() == ["P5-Bi8-12", "P5-Bi0", "Base0Ca"]
    assert [r.group for r in dlg.table().rows] == ["P5-Bi8-12", "P5-Bi0", "Base0Ca"]
    locked = _dialog(locked=True)
    assert not locked.btnUp.isEnabled() and not locked.btnDown.isEnabled()
    assert not locked.btnSort.isEnabled() and locked.btnJoin.isEnabled()
    assert "waits" in locked.btnSort.toolTip()


def test_series_table_dialog_typed_columns_edit_sort_and_remove(qapp):
    dlg = _dialog()
    dlg._add_column("CaO (mol%)")
    assert dlg.grid.columnCount() == 6
    assert dlg.grid.horizontalHeaderItem(5).text() == "CaO (mol%)"
    dlg.grid.item(0, 5).setText("2.5")
    dlg.grid.item(1, 5).setText("abc")                      # not a number -> empty
    dlg.grid.item(2, 5).setText("1")
    t = dlg.table()
    assert [r.values["CaO (mol%)"] for r in t.rows] == [2.5, None, 1.0]
    assert dlg.grid.item(1, 5).text() == ""
    assert "1 column(s)" in dlg.status.text()
    dlg._sort_by("CaO (mol%)")
    assert dlg.table().labels() == ["Base0Ca", "P5-Bi8-12", "P5-Bi0"]   # missing last
    dlg._sort_by("CaO (mol%)", descending=True)
    assert dlg.table().labels() == ["P5-Bi8-12", "Base0Ca", "P5-Bi0"]
    dlg._remove_column("CaO (mol%)")
    assert dlg.grid.columnCount() == 5 and dlg.table().columns == []
    dlg._remove_column()
    assert "no columns" in dlg.status.text()


_HEADER = ["sample", "series", "n_matrix", "structure", "SiO2_mol", "SiO2_mol_sd",
           "P2O5_mol", "P2O5_mol_sd", "Bi2O3_mol", "Bi2O3_mol_sd", "method", "Vm", "Vm_u"]
_ROWS = [["P5Bi8-12", "P5", "15", "homogeneous", "48.1", "0.8", "5.2", "0.1", "7.9", "0.3",
          "plain means", "27.3", "0.09"],
         ["P5Bi0", "P5", "10", "homogeneous", "52.0", "0.9", "5.0", "0.1", "0.0", "0.0",
          "plain means", "26.7", "0.06"]]


def _batch_csv(tmp_path, labels):
    from larmor import batchfit
    from larmor.recipe import Recipe, SiteModel, Param
    recs = []
    for k, lab in enumerate(labels):
        a = Param(100.0 + k); a.stderr = 0.5
        recs.append(Recipe(nucleus="31P", larmor_frequency_MHz=202.0, sample=lab, sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(10.0), "shift_fwhm_ppm": Param(6.0),
                "amplitude": a, "gl": Param(1.0, vary=False)})]))
    res = batchfit.BatchFitResult(recipes=recs, labels=list(labels), rmsd=[0.0] * len(recs),
                                  per_dataset=[], shared=(), released=())
    p = tmp_path / "batch_table.csv"
    batchfit.write_shared_csv(res, p)
    return str(p)


def test_join_csv_dialog_prefills_mapping_and_adds_tagged_columns(qapp, tmp_path, monkeypatch):
    from larmor.desktop import series_table_dialog as mod
    from larmor.desktop.comparability_dialog import _CHECK_CSS_COLOR

    epma = tmp_path / "EPMA_analyzed.csv"
    with open(epma, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(_HEADER)
        w.writerows(_ROWS)
    captured = {}

    def fake_exec(self):
        captured["dlg"] = self
        return QDialog.Accepted
    monkeypatch.setattr(mod.JoinMappingDialog, "exec", fake_exec)
    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(epma), "")))
    dlg = _dialog()
    dlg._join_csv()
    jd = captured["dlg"]
    assert jd.tag() == "analysed"                            # 'EPMA' in the file name
    assert jd.chosen_columns() == ["SiO2_mol", "P2O5_mol", "Bi2O3_mol"]   # oxides ticked, Vm not
    names = [jd.cols.item(i).data(Qt.UserRole) for i in range(jd.cols.count())]
    assert "Vm" in names and "n_matrix" in names and "SiO2_mol_sd" not in names
    assert "(± Vm_u)" in jd.cols.item(names.index("Vm")).text()
    assert jd._combos[0].currentText() == "P5Bi8-12" and jd._combos[1].currentText() == "P5Bi0"
    assert jd._combos[2].currentText() == "(none)"
    assert jd.map.item(2, 0).background().color().name() == QColor(_CHECK_CSS_COLOR).name()
    assert jd.map.item(0, 0).background().color().alpha() == 0
    assert [m.how for m in jd.matches()] == ["normalised", "normalised", "none"]
    t = dlg.table()
    assert [c.label for c in t.columns] == ["analysed SiO2_mol", "analysed P2O5_mol",
                                            "analysed Bi2O3_mol"]
    assert t.columns[1].err_key == "P2O5_mol_sd" and t.columns[1].source == "EPMA_analyzed.csv"
    assert t.rows[0].values["P2O5_mol"] == 5.2 and t.rows[0].values["P2O5_mol_sd"] == 0.1
    assert t.rows[2].values["P2O5_mol"] is None
    assert "matched 2/3 from EPMA_analyzed.csv (analysed)" in dlg.status.text()
    assert dlg.grid.columnCount() == 5 + 6
    assert dlg.grid.horizontalHeaderItem(6).text() == "analysed SiO2_mol ±"
    assert dlg.grid.item(0, 7).text() == "5.2"

    # a long LARMOR batch CSV is recognised and offered pivoted, every column ticked
    bpath = _batch_csv(tmp_path, ["P5-Bi0", "Base0Ca", "other"])
    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (bpath, "")))
    dlg2 = _dialog()
    dlg2._join_csv()
    jd2 = captured["dlg"]
    heads = [jd2.cols.item(i).data(Qt.UserRole) for i in range(jd2.cols.count())]
    assert "A amplitude" in heads and "A population_pct" in heads and "scope" not in heads
    assert jd2.tag() == "" and set(jd2.chosen_columns()) == set(heads)
    assert [m.how for m in jd2.matches()] == ["none", "exact", "exact"]
    t2 = dlg2.table()
    assert t2.column("A amplitude").err_key == "A amplitude ±"
    assert t2.rows[1].values["A amplitude"] == 100.0 and t2.rows[2].values["A amplitude"] == 101.0
    assert t2.rows[0].values["A amplitude"] is None
    assert "matched 2/3 from batch_table.csv" in dlg2.status.text()

    # Cancel joins nothing
    monkeypatch.setattr(mod.JoinMappingDialog, "exec", lambda self: QDialog.Rejected)
    dlg3 = _dialog()
    dlg3._join_csv()
    assert dlg3.table().columns == [] and "cancelled" in dlg3.status.text()
    # a manual pick in the confirm dialog is honoured
    monkeypatch.setattr(mod.JoinMappingDialog, "exec", fake_exec)
    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(epma), "")))
    orig_init = mod.JoinMappingDialog.__init__

    def picking_init(self, *a, **k):
        orig_init(self, *a, **k)
        self._combos[2].setCurrentIndex(2)                  # Base0Ca -> P5Bi0 by hand
    monkeypatch.setattr(mod.JoinMappingDialog, "__init__", picking_init)
    dlg4 = _dialog()
    dlg4._join_csv()
    assert dlg4.table().rows[2].values["P2O5_mol"] == 5.0
    assert captured["dlg"].matches()[2].how == "manual"
    assert "matched 3/3" in dlg4.status.text()


def test_save_and_load_table_restore_names_groups_columns_and_order(qapp, tmp_path):
    dlg = _dialog()
    dlg._add_column("CaO")
    for r, v in zip(range(3), ("0", "2", "4")):
        dlg.grid.item(r, 5).setText(v)
    dlg.grid.item(0, 1).setText("renamed")
    dlg.grid.item(1, 2).setText("G")
    dlg.grid.selectRow(1); dlg.grid.setCurrentCell(1, 1)
    dlg._move(-1)
    assert dlg.perm() == [1, 0, 2]
    p = tmp_path / "t.series.json"
    dlg._save_table(str(p))
    assert p.exists() and "saved t.series.json" in dlg.status.text()

    fresh = _dialog()
    fresh._load_table(str(p))
    assert fresh.perm() == [1, 0, 2]                         # the saved order is followed
    assert fresh.table().labels() == dlg.table().labels() == ["P5-Bi0", "renamed", "Base0Ca"]
    assert [r.group for r in fresh.table().rows] == ["G", "renamed", "Base0Ca"]
    assert [r.values["CaO"] for r in fresh.table().rows] == [2.0, 0.0, 4.0]
    assert "loaded 3/3 rows from t.series.json" in fresh.status.text()
    # over a different series only the rows that pair (by name) take values
    from larmor.desktop.series_table_dialog import SeriesTableDialog
    from larmor.series_table import SeriesTable
    other = SeriesTableDialog(None, SeriesTable.from_paths(["x", "y"], names=["Base0Ca", "new"]))
    other._load_table(str(p))
    assert other.table().labels() == ["Base0Ca", "new"]
    assert [r.values["CaO"] for r in other.table().rows] == [4.0, None]
    assert "loaded 1/2 rows" in other.status.text()
    other._load_table(str(tmp_path / "missing.series.json"))
    assert "could not load" in other.status.text()
