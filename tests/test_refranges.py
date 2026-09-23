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


def test_19f_fluoride_ladder_positions():
    """The crystalline MF..MF4 ladder (Bureau 1997): single positions vs
    CFCl3, sourced, ordered by cation as the superposition model reads."""
    pos = refranges.positions_for("19F")
    by = {}
    for r in pos:
        assert r["ref"] in refranges.REFS and r["note"], r
        assert -240 < r["ppm"] < 120, r
        by.setdefault(r["label"], []).append(r["ppm"])
    assert by["NaF"] == [-224.0] and by["LiF"] == [-204.0]
    assert by["CaF2"] == [-108.0] and by["BaF2"] == [-14.0]
    assert by["KF"] == [-133.0] and by["CsF"] == [-11.0]
    assert sorted(by["LaF3"]) == [-23.0, 25.0]      # two-site fluoride
    # the alkali ladder runs Cs -> Rb -> K -> Li -> Na towards shielding
    assert by["CsF"] > by["RbF"] > by["KF"] > by["LiF"] > by["NaF"]
    assert refranges.positions_for("27Al") == []
    assert refranges.positions_for(None) == []


def test_assign_reads_a_position_into_a_species():
    """A line's position maps to the literature band that holds it (the
    narrowest of overlapping bands), or to a reported compound within
    POSITION_TOLERANCE_PPM; auto-generated labels are recognised."""
    al = refranges.assign("27Al", 65.0)
    assert al["label"] == "Al[4]" and al["kind"] == "range"
    assert refranges.assign("27Al", 38.0)["label"] == "Al[5]"
    assert refranges.assign("27Al", 200.0) is None
    b = refranges.assign("11B", 0.5)
    assert b["label"].startswith("B[4]")
    f = refranges.assign("19F", -224.0)          # NaF sits inside the F-Na(n) band
    assert f["kind"] == "range" and f["label"] == "F\u2013Na(n)"
    f2 = refranges.assign("19F", -12.0)          # only the ladder reaches here
    assert f2["kind"] == "position" and f2["label"] == "CsF"
    assert refranges.assign("19F", 60.0) is None  # nothing within tolerance
    assert refranges.assign("7Li", 0.0) is None

    for auto in ("", None, "Czjzek-3", "pk-0", "read-1", "HB-2", "line-copy",
                 "Al-copy", "A+1sb", "pk-0-1sb"):
        assert refranges.is_auto_label(auto), auto
    for own in ("AlIV", "Al[4]", "Q3 site", "B4 ring", "s0 main"):
        assert not refranges.is_auto_label(own), own


def test_label_lines_from_literature_in_the_app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ["LARMOR_NO_SESSION"] = "1"
    pytest.importorskip("PySide6")
    import numpy as np
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow
    from larmor.recipe import Param, Recipe, SiteModel

    win = MainWindow()
    try:
        x = np.linspace(-50.0, 120.0, 601)
        win._display_1d(x, np.exp(-((x - 65.0) / 5.0) ** 2), "27Al", 130.3,
                        14000.0, "t", "x")
        rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3,
                     spin_rate_Hz=14000.0, sites=[
            SiteModel(model="gauss_lor", label="pk-0", params={
                "isotropic_chemical_shift_ppm": Param(65.0),
                "shift_fwhm_ppm": Param(5.0), "amplitude": Param(1.0),
                "gl": Param(0.5)}),
            SiteModel(model="gauss_lor", label="my octahedral", params={
                "isotropic_chemical_shift_ppm": Param(5.0),
                "shift_fwhm_ppm": Param(5.0), "amplitude": Param(1.0),
                "gl": Param(0.5)}),
            SiteModel(model="gauss_lor", label="pk-2", params={
                "isotropic_chemical_shift_ppm": Param(110.0),
                "shift_fwhm_ppm": Param(5.0), "amplitude": Param(1.0),
                "gl": Param(0.5)})]).to_dict()
        rec["sites"][1]["family"] = "mine"          # a user tag: never overwritten
        win.recipe["sites"] = rec["sites"]
        win.on_structure_changed()
        win.label_from_literature()
        labels = [s["label"] for s in win.recipe["sites"]]
        assert labels == ["Al[4]", "my octahedral", "pk-2"]
        msg = win.statusBar().currentMessage()
        assert "A=Al[4]" in msg and "kept your own" in msg and "outside" in msg
        # N3: the Al[4] band also fills the EMPTY family tag; the pre-tagged
        # line keeps its own; the line outside every band stays untagged
        fams = [s.get("family", "") for s in win.recipe["sites"]]
        assert fams == ["Al(IV)", "mine", ""]
        assert "tagged 1 family" in msg
        assert [f["family"] for f in win._last_quant["families"]] == ["Al(IV)", "mine"]
        win.undo()
        assert win.recipe["sites"][0]["label"] == "pk-0"
        assert "family" not in win.recipe["sites"][0]
        # a line dropped inside the Al[4] band is pre-tagged and the status says so
        win._model_actions["gauss_lor"].setChecked(True)
        win.add_site_at(62.0, 1.0)
        added = win.recipe["sites"][-1]
        assert added["family"] == "Al(IV)"
        assert "family Al(IV) (literature band" in win.statusBar().currentMessage()
        win.add_site_at(110.0, 1.0)                  # outside every band: no tag
        assert "family" not in win.recipe["sites"][-1]
        assert "family" not in win.statusBar().currentMessage()
        win._model_actions["gauss_lor"].setChecked(False)
    finally:
        win.close()


def test_ranges_carry_families_only_where_the_band_names_one():
    """N3: the 27Al and 11B bands name a structural species, so they seed a
    site's family tag; 17O (overlapping BO/NBO bands), 29Si and 31P (a band
    cannot name a Qn) carry none."""
    assert refranges.family_for("27Al", 65.0) == "Al(IV)"
    assert refranges.family_for("27Al", 37.0) == "Al(V)"
    assert refranges.family_for("27Al", 5.0) == "Al(VI)"
    assert refranges.family_for("11B", 15.0) == "BO3"
    assert refranges.family_for("11B", 0.5) == "BO4"
    assert refranges.family_for("29Si", -90.0) == ""
    assert refranges.family_for("17O", 50.0) == ""
    assert refranges.family_for("31P", 0.0) == ""
    assert refranges.family_for("7Li", 0.0) == "" and refranges.family_for(None, 1.0) == ""
    assert refranges.family_for("27Al", 200.0) == ""          # outside every band
    for nuc, entries in refranges.REF_RANGES.items():
        for r in entries:
            assert ("family" in r) == (nuc in ("27Al", "11B")), (nuc, r["label"])
    # every seeded family is a preset of its nucleus
    from larmor import families
    for nuc in ("27Al", "11B"):
        for r in refranges.REF_RANGES[nuc]:
            assert r["family"] in families.presets_for(nuc)


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

    # 19F: five glass bands plus the crystalline ladder as ticks
    f_pos = refranges.positions_for("19F")
    v.set_ref_ranges(refranges.ranges_for("19F"), refranges.CITATION, f_pos)
    ticks = [i for i in v._ref_items if isinstance(i, pg.InfiniteLine)]
    assert len(ticks) == len(f_pos) and all(not t.movable for t in ticks)
    assert any("NaF" in t.toolTip() for t in ticks)
    v.set_ref_ranges(None)
    assert v._ref_items == []
    v.close()
