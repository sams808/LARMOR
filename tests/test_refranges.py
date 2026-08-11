"""Literature shift-range overlay (View > Literature shift ranges):
the Qt-free data module and the SpectrumView/app wiring. Data values are
SOURCED (Eden 2023) -- these tests pin the structure and sanity, not the
science."""
import os

import pytest

from larmor import refranges


def test_ranges_structure_and_sanity():
    assert set(refranges.REF_RANGES) == {"27Al", "11B", "29Si", "31P",
                                         "17O", "23Na", "25Mg", "19F"}
    for nuc, entries in refranges.REF_RANGES.items():
        assert entries, nuc
        for r in entries:
            assert r["lo_ppm"] < r["hi_ppm"], (nuc, r["label"])
            assert r["label"]
            assert "quad" in r and "note" in r
            # every entry names its primary source, and it resolves
            assert r["ref"] in refranges.REFS, (nuc, r["label"])
            assert refranges.citation_for(r) == refranges.REFS[r["ref"]]
    assert "Ed\u00e9n 2023" in refranges.CITATION


def test_new_nuclei_values_are_sourced_sanely():
    """Spot-check the added nuclei against their primary sources."""
    o = {r["label"]: r for r in refranges.ranges_for("17O")}
    assert o["Si\u2013O\u2013Si (BO)"]["ref"] == "dirken1997"
    assert o["NBO (Si\u2013O\u2013M)"]["ref"] == "du2003"
    assert o["NBO (Si\u2013O\u2013M)"]["hi_ppm"] == 75.0        # K-NBO at 71 fits inside

    na = refranges.ranges_for("23Na")
    assert len(na) == 1 and na[0]["lo_ppm"] == -20.0 and na[0]["hi_ppm"] == 10.0

    mg = {r["label"]: r for r in refranges.ranges_for("25Mg")}
    assert mg["Mg[6]"]["ref"] == "shimoda2007"
    assert mg["Mg[4]/Mg[5]"]["lo_ppm"] == 30.0

    # 19F: spin-1/2 (no quad note), all NEGATIVE vs CFCl3, Baasner ladder
    f = refranges.ranges_for("19F")
    assert len(f) == 5
    assert all(r["hi_ppm"] < 0 for r in f)
    assert all(r["quad"] == "" for r in f)
    centers = sorted((r["lo_ppm"] + r["hi_ppm"]) / 2 for r in f)
    assert centers == [-225.0, -188.0, -168.0, -146.0, -113.0]


def test_ranges_for_normalizes_and_defaults_empty():
    assert refranges.ranges_for("27Al")
    assert refranges.ranges_for(" 27Al ")          # stray whitespace tolerated
    assert refranges.ranges_for("7Li") == []       # not compiled -> empty
    assert refranges.ranges_for(None) == []
    assert refranges.ranges_for("") == []


def test_quadrupolar_notes_present_where_they_matter():
    """The 'what about Cq' half of the feature: quadrupolar nuclei carry
    their P_Q/C_Q ranges as label/tooltip text (a width is not a
    shift-axis quantity, so it can't be a span)."""
    al = {r["label"]: r for r in refranges.ranges_for("27Al")}
    assert all("P_Q" in r["quad"] for r in al.values())
    b = {r["label"]: r for r in refranges.ranges_for("11B")}
    assert "2.4" in b["B[3] (BO3)"]["quad"]        # BO3 C_Q 2.4-2.8 MHz
    assert "0.2" in b["B[4] (BO4)"]["quad"]        # BO4 C_Q 0.2-0.8 MHz


def test_spectrum_view_draws_and_clears_ref_ranges():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    import pyqtgraph as pg
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    v.set_ref_ranges(refranges.ranges_for("27Al"), refranges.CITATION)
    items = v._ref_items
    regions = [i for i in items if isinstance(i, pg.LinearRegionItem)]
    labels = [i for i in items if isinstance(i, pg.TextItem)]
    assert len(regions) == 3 and len(labels) == 3      # Al[4]/Al[5]/Al[6]
    assert all(not r.movable for r in regions)          # guide, not a control
    assert "P_Q" in regions[0].toolTip()
    assert "Ed\u00e9n 2023" in regions[0].toolTip()

    v.set_ref_ranges(None)                              # clears completely
    assert v._ref_items == []
    v.close()
