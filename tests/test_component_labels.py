"""Small workbench comforts: component names on the plot (pinned or on
hover), keyboard line removal and arrow-key nudging in the parameter table,
and sidebands that are LINKED to their parent line instead of free copies."""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _two_components():
    x = np.linspace(-50.0, 50.0, 401)
    a = np.exp(-0.5 * ((x - 20.0) / 3.0) ** 2)
    b = 0.5 * np.exp(-0.5 * ((x + 10.0) / 5.0) ** 2)
    return x, a, b


def test_component_labels_pin_at_each_maximum(qapp):
    import pyqtgraph as pg
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    x, a, b = _two_components()
    v.set_model(x, a + b, [a, b], ["AlIV", ""], set(), x, a + b)
    assert v._comp_labels == []                      # off by default
    v.set_show_labels(True)
    labels = {lab.toPlainText(): lab.pos() for lab in v._comp_labels}
    assert set(labels) == {"A \u00b7 AlIV", "B"}       # letter + name; letter alone
    assert labels["A \u00b7 AlIV"].x() == pytest.approx(20.0, abs=0.3)
    assert labels["B"].x() == pytest.approx(-10.0, abs=0.3)
    assert all(isinstance(lab, pg.TextItem) for lab in v._comp_labels)

    # a hidden component gets no label; clearing the model clears the labels
    v.set_model(x, a + b, [a, b], ["AlIV", "AlVI"], {1}, x, a + b)
    assert [lab.toPlainText() for lab in v._comp_labels] == ["A \u00b7 AlIV"]
    v.set_model(None, None, [], [], set())
    assert v._comp_labels == []
    v.set_show_labels(False)
    v.close()


def test_hover_names_the_tallest_component_under_the_cursor(qapp):
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    x, a, b = _two_components()
    v.set_model(x, a + b, [a, b], ["AlIV", "AlVI"], set(), x, a + b)
    v._hover_component(20.0)
    assert v._hover_label is not None and v._hover_label.isVisible()
    assert v._hover_label.toPlainText() == "A \u00b7 AlIV"
    v._hover_component(-10.0)
    assert v._hover_label.toPlainText() == "B \u00b7 AlVI"
    v._hover_component(45.0)                         # nothing there
    assert not v._hover_label.isVisible()
    # pinned labels take over: the hover label is hidden
    v._hover_component(20.0)
    v.set_show_labels(True)
    assert not v._hover_label.isVisible()
    v.close()


def _recipe_dict(n_sites=2):
    from larmor.recipe import Param, Recipe, SiteModel

    sites = [SiteModel(model="gauss_lor", label=f"L{i}", params={
        "isotropic_chemical_shift_ppm": Param(10.0 * i),
        "shift_fwhm_ppm": Param(4.0), "amplitude": Param(1.0),
        "gl": Param(0.5)}) for i in range(n_sites)]
    return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3,
                  spin_rate_Hz=13030.0, sites=sites).to_dict()


def test_delete_key_removes_the_selected_line(qapp):
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest
    from larmor.desktop.table import LinesTable

    t = LinesTable()
    got = []
    t.structure.connect(lambda row, action: got.append((row, action)))
    t.rebuild(_recipe_dict(3), set())
    t.table.setCurrentCell(1, 0)
    QTest.keyClick(t.table, Qt.Key_Delete)
    assert got == [(1, "remove")]
    QTest.keyClick(t.table, Qt.Key_D, Qt.ControlModifier)
    assert got[-1] == (1, "duplicate")
    QTest.keyClick(t.table, Qt.Key_H, Qt.ControlModifier)
    assert got[-1] == (1, "visibility")
    # no selection -> no action
    got.clear()
    t.table.setCurrentCell(-1, -1)
    QTest.keyClick(t.table, Qt.Key_Delete)
    assert got == []
    t.close()


def test_arrow_keys_nudge_a_parameter_cell(qapp):
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest
    from larmor.desktop.table import LinesTable

    t = LinesTable()
    edits = []
    t.edited.connect(lambda: edits.append(1))
    rec = _recipe_dict(2)
    rec["sites"][1]["params"]["shift_fwhm_ppm"]["expr"] = "s0.shift_fwhm_ppm"
    t.rebuild(rec, set())
    col = 2 + t._used_keys.index("shift_fwhm_ppm")
    cell = t.table.cellWidget(0, col)
    p = rec["sites"][0]["params"]["shift_fwhm_ppm"]
    QTest.keyClick(cell.edit, Qt.Key_Up)
    assert p["value"] == pytest.approx(4.0 * 1.02)
    QTest.keyClick(cell.edit, Qt.Key_Down, Qt.ShiftModifier)
    assert p["value"] == pytest.approx(4.08 * 0.9)
    QTest.keyClick(cell.edit, Qt.Key_PageUp)
    assert p["value"] == pytest.approx(4.08 * 0.9 * 1.10)
    assert float(cell.edit.text()) == pytest.approx(p["value"], rel=1e-4)
    assert len(edits) == 3
    # bounds clamp; a linked parameter is not nudged
    p["max"] = p["value"]
    QTest.keyClick(cell.edit, Qt.Key_Up)
    assert p["value"] == pytest.approx(p["max"])
    linked = t.table.cellWidget(1, col)
    before = rec["sites"][1]["params"]["shift_fwhm_ppm"]["value"]
    QTest.keyClick(linked.edit, Qt.Key_Up)
    assert rec["sites"][1]["params"]["shift_fwhm_ppm"]["value"] == before
    t.close()


def test_sidebands_are_linked_to_their_parent(qapp, monkeypatch):
    """Decomposition > Add spinning sidebands: position = parent +- k*nu_rot as
    a constraint, every shape parameter tied to the parent, amplitude free."""
    from PySide6.QtWidgets import QDialog
    from larmor.desktop.app import MainWindow
    from larmor.recipe import Recipe

    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Accepted)
    win = MainWindow()
    try:
        x = np.linspace(-100.0, 100.0, 801)
        y = np.exp(-0.5 * ((x - 10.0) / 3.0) ** 2)
        win._display_1d(x, y, "27Al", 130.3, 13030.0, "t", "x")
        win.recipe["sites"] = _recipe_dict(1)["sites"]
        win.recipe["spin_rate_Hz"] = 13030.0
        win.on_structure_changed()
        win.add_sidebands()                          # defaults: +1 / -1, linked
        sites = win.recipe["sites"]
        assert len(sites) == 3
        nur = 13030.0 / 130.3
        plus = next(s for s in sites if s["label"].endswith("+1sb"))
        minus = next(s for s in sites if s["label"].endswith("-1sb"))
        p = plus["params"]
        assert p["isotropic_chemical_shift_ppm"]["expr"] == \
            f"s0.isotropic_chemical_shift_ppm + ({nur:.6g})"
        assert p["isotropic_chemical_shift_ppm"]["value"] == pytest.approx(nur)
        assert minus["params"]["isotropic_chemical_shift_ppm"]["expr"] == \
            f"s0.isotropic_chemical_shift_ppm + ({-nur:.6g})"
        assert p["shift_fwhm_ppm"]["expr"] == "s0.shift_fwhm_ppm"
        assert p["gl"]["expr"] == "s0.gl"
        assert p["amplitude"]["expr"] is None and p["amplitude"]["vary"]
        assert p["amplitude"]["value"] == pytest.approx(0.3)
        # the constraints are valid for the fit translator
        from larmor.fit import translate_expr
        rec = Recipe.from_dict(win.recipe)
        for s in sites[1:]:
            for prm in s["params"].values():
                if prm["expr"]:
                    translate_expr(prm["expr"], rec)
    finally:
        win.close()
