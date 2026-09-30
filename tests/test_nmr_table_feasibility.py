"""The NMR table shows how practical each element is for NMR of glasses
(Youngman 2018, Fig. 1) as the cell fill, with the spin still on the border."""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from larmor import nuclei as N  # noqa: E402


def test_feasibility_partitions_the_whole_periodic_table():
    grid = {s for row in N.PERIODIC_ROWS for s in row if s}
    cats = list(N.FEASIBILITY)
    assert cats == ["favorable", "challenging", "very_difficult", "impractical",
                    "impossible"]
    seen = []
    for c in cats:
        seen += list(N.FEASIBILITY[c])
    assert len(seen) == len(set(seen)), "an element sits in two categories"
    assert set(seen) == grid, sorted(set(seen) ^ grid)
    assert set(N.FEASIBILITY_LABEL) == set(N.FEASIBILITY_NOTE) == set(cats)
    assert "Youngman" in N.FEASIBILITY_SOURCE


@pytest.mark.parametrize("element,category", [
    ("Al", "favorable"), ("Si", "favorable"), ("B", "favorable"),
    ("Na", "favorable"), ("F", "favorable"), ("Pb", "favorable"),
    ("V", "favorable"), ("Sc", "favorable"),
    ("O", "challenging"), ("Sn", "challenging"), ("Cs", "challenging"),
    ("Ga", "challenging"), ("Tl", "challenging"),
    ("Mg", "very_difficult"), ("Ca", "very_difficult"), ("Ge", "very_difficult"),
    ("N", "very_difficult"), ("La", "very_difficult"), ("Nb", "very_difficult"),
    ("Ti", "impractical"), ("Fe", "impractical"), ("Zn", "impractical"),
    ("Cl", "impractical"), ("Br", "impractical"), ("Bi", "impractical"),
    ("Ar", "impossible"), ("Ce", "impossible"), ("Th", "impossible"),
    ("Fr", "impossible"),
])
def test_feasibility_spot_checks_against_the_figure(element, category):
    assert N.feasibility(element) == category


def test_unknown_element_has_no_category():
    assert N.feasibility("Xx") is None


pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def test_table_cells_are_filled_by_category_and_say_why(qapp):
    from larmor.desktop.utilities import FEASIBILITY_FILL, NmrTableDialog

    dlg = NmrTableDialog(None, 400.0)
    try:
        assert set(FEASIBILITY_FILL) == set(N.FEASIBILITY)
        al = dlg._buttons["Al"]
        assert FEASIBILITY_FILL["favorable"] in al.styleSheet()
        assert "Favorable" in al.toolTip() and "Youngman" in al.toolTip()
        assert dlg.cell_fill("O") == FEASIBILITY_FILL["challenging"]
        assert dlg.cell_fill("Ar") == FEASIBILITY_FILL["impossible"]
        # the spin border is still there
        assert "border: 1.5px solid" in al.styleSheet()
        assert dlg.legend.isVisibleTo(dlg)
        # fill off: plain cells again, tooltips keep the category text
        dlg.shade.setCurrentIndex(1)
        assert dlg.cell_fill("Al") is None
        assert "background" not in dlg._buttons["Al"].styleSheet()
        assert not dlg.legend.isVisibleTo(dlg)
        assert "Favorable" in dlg._buttons["Al"].toolTip()
    finally:
        dlg.close()
