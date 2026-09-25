"""Per-overlay display controls in the Datasets dock.

Sam (2026-09): "When adding a spectrum to compare in Datasets, we must be
able to modify colour, scaling and shifting of the added spectrum easily."

Every overlay row has, besides the colour swatch, a scale factor (x), an x
shift (ppm) and a y offset (a fraction of the active span); 'match height'
fills the scales with the matching factor and stays a one-click helper; a
right-click menu offers Reset / Make active / Remove. The transform is
applied at draw time (larmor.display.overlay_display) -- the stored arrays
are never modified -- and travels with the overlay through the workspace
snapshot and the project bundle (neutral values for older files)."""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QDoubleSpinBox, QFileDialog  # noqa: E402

from larmor import display, project  # noqa: E402
from larmor.io import spectra  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _active(win):
    ppm = np.linspace(80.0, -20.0, 600)
    amp = 100.0 * np.exp(-((ppm - 30.0) / 6.0) ** 2)
    win._display_1d(ppm, amp, "27Al", 130.3, 20000.0, "active", "active")
    return ppm, amp


def _csv(tmp_path, name, scale):
    ppm = np.linspace(80.0, -20.0, 500)
    amp = scale * np.exp(-((ppm - 10.0) / 4.0) ** 2)
    return spectra.write_csv(tmp_path / name, ppm, amp, {"nucleus": "27Al"}), amp


def _row(win, i=0):
    return win.datasets_panel._row_widgets[i]


def test_row_controls_change_the_drawn_item_not_the_arrays(win, tmp_path):
    ppm, amp = _active(win)
    path, ref = _csv(tmp_path, "ref.csv", 5.0)
    assert win.add_overlay_path(str(path))
    assert "scale, shift, offset" in win.statusBar().currentMessage()
    ov = win._overlays[0]
    assert (ov["scale"], ov["shift"], ov["yoff"]) == (1.0, 0.0, 0.0)
    ppm0, amp0 = ov["ppm"].copy(), ov["amp"].copy()
    ref_max = float(amp0.max())
    span = float(amp.max() - amp.min())
    row = _row(win)
    for key in ("scale", "shift", "yoff"):
        assert isinstance(row[key], QDoubleSpinBox)
    assert "reference peak" in row["shift"].toolTip()

    row["scale"].setValue(2.5)                          # the spin box drives the dict
    assert ov["scale"] == 2.5
    item = win.view._overlay_items[0]
    assert item.yData.max() == pytest.approx(2.5 * ref_max, rel=1e-9)
    row["shift"].setValue(1.2)
    item = win.view._overlay_items[0]
    assert np.allclose(item.xData, ppm0 + 1.2)
    row["yoff"].setValue(0.3)
    item = win.view._overlay_items[0]
    assert np.allclose(item.yData, 2.5 * amp0 + 0.3 * span)
    # the global stack offset still adds on top
    win.datasets_panel.offset.setValue(0.4)
    item = win.view._overlay_items[0]
    assert np.allclose(item.yData, 2.5 * amp0 + 0.7 * span)
    # stored arrays untouched; the label carries the badge; rows updated in place
    assert np.array_equal(ov["ppm"], ppm0) and np.array_equal(ov["amp"], amp0)
    assert "\u00d72.5 \u00b7 +1.2 ppm \u00b7 \u21910.3" in _row(win)["lab"].text()
    assert _row(win)["scale"] is row["scale"]
    assert row["scale"].value() == 2.5 and row["shift"].value() == 1.2

    # right-click > Reset: as stored again
    menu = win.datasets_panel.row_menu(0)
    texts = [a.text() for a in menu.actions() if not a.isSeparator()]
    assert texts == ["Reset scale / shift / offset", "Make active", "Remove"]
    assert [a for a in menu.actions() if a.text() == "Make active"][0].isEnabled()
    [a for a in menu.actions() if a.text().startswith("Reset")][0].trigger()
    ov = win._overlays[0]
    assert (ov["scale"], ov["shift"], ov["yoff"]) == (1.0, 0.0, 0.0)
    item = win.view._overlay_items[0]
    assert np.allclose(item.yData, amp0 + 0.4 * span) and np.allclose(item.xData, ppm0)
    assert display.overlay_badge(1.0, 0.0, 0.0) == "" and "\u00d7" not in _row(win)["lab"].text()
    assert _row(win)["scale"].value() == 1.0
    # Remove through the menu
    [a for a in menu.actions() if a.text() == "Remove"][0].trigger()
    assert win._overlays == [] and win.view._overlay_items == []


def test_match_height_fills_the_scale_and_stays_a_helper(win, tmp_path):
    ppm, amp = _active(win)
    path, ref = _csv(tmp_path, "ref.csv", 5.0)
    win.add_overlay_path(str(path))
    ref_max = float(win._overlays[0]["amp"].max())
    win.datasets_panel.match.setChecked(True)
    factor = float(amp.max()) / ref_max
    assert win._overlays[0]["scale"] == pytest.approx(factor)
    assert _row(win)["scale"].value() == pytest.approx(factor, rel=1e-3)
    assert win.view._overlay_items[0].yData.max() == pytest.approx(float(amp.max()), rel=1e-9)
    # an overlay added while ticked starts matched
    path2, _ = _csv(tmp_path, "ref2.csv", 20.0)
    win.add_overlay_path(str(path2))
    assert win.view._overlay_items[1].yData.max() == pytest.approx(float(amp.max()), rel=1e-9)
    # a hand-typed scale unticks the box (the scales are no longer 'matched')
    _row(win, 1)["scale"].setValue(3.0)
    assert not win.datasets_panel.match.isChecked()
    assert win._overlays[0]["scale"] == pytest.approx(factor)      # the other keeps its factor
    assert win._overlays[1]["scale"] == 3.0
    # unticking by hand puts every scale back to x1
    win.datasets_panel.match.setChecked(True)
    win.datasets_panel.match.setChecked(False)
    assert [ov["scale"] for ov in win._overlays] == [1.0, 1.0]
    assert win.view._overlay_items[0].yData.max() == pytest.approx(ref_max, rel=1e-9)


def test_match_height_matches_the_drawn_peaks_in_every_y_mode(win, tmp_path):
    """Group A (2026-09-24 review): the stored × was the ratio of the two RAW
    maxima while each overlay is drawn normalised by its OWN trace, so under
    View > Y axis > Normalise to maximum the overlay came out
    max(active)/max(overlay) times the active spectrum -- 20x here. The
    factor is now a display-space one: the drawn peaks agree in every mode,
    and raw mode still stores the raw ratio."""
    ppm, amp = _active(win)
    path, _ref = _csv(tmp_path, "weak.csv", 5.0)
    assert win.add_overlay_path(str(path))
    ov = win._overlays[0]
    raw_ratio = float(amp.max()) / float(ov["amp"].max())
    assert raw_ratio == pytest.approx(20.0, rel=1e-3)
    stored = ov["amp"].copy()

    def drawn():
        return (float(win.view._exp.yData.max()),
                float(win.view._overlay_items[0].yData.max()))

    for mode, region, expect in (("raw", None, raw_ratio),
                                 ("max", None, 1.0),
                                 ("area", None, None),
                                 ("region", (60.0, 0.0), None)):
        win.view.set_y_mode(mode, region)
        win.datasets_panel.match.setChecked(False)
        win.datasets_panel.match.setChecked(True)
        act_peak, ov_peak = drawn()
        assert ov_peak == pytest.approx(act_peak, rel=1e-9), mode
        if expect is not None:
            assert ov["scale"] == pytest.approx(expect, rel=1e-3), mode
        # the box the user reads shows the factor actually in use
        assert _row(win)["scale"].value() == pytest.approx(ov["scale"], rel=2e-3), mode
        # RAW units everywhere but the canvas
        assert np.array_equal(ov["amp"], stored) and np.array_equal(win.exp_amp, amp)
    win.view.set_y_mode("raw")


def test_transform_round_trips_through_snapshot_and_project(win, tmp_path, monkeypatch):
    win._confirm_open_mode = lambda: "replace"
    win._report_project_notes = lambda notes: None
    win._add_recent = lambda p: None
    ppm, amp = _active(win)
    path, ref = _csv(tmp_path, "compare.csv", 5.0)
    win.add_overlay_path(str(path))
    win.overlay_set_scale(0, 2.5)
    win.overlay_set_shift(0, -1.5)
    win.overlay_set_yoff(0, 0.25)
    win._sync_active()                     # the workspace registered by _display_1d
    ws = win.workspaces[win.active_ws]
    # the workspace snapshot carries the whole overlay dict
    snap_ov = ws["snap"]["overlays"][0]
    assert (snap_ov["scale"], snap_ov["shift"], snap_ov["yoff"]) == (2.5, -1.5, 0.25)
    entry = project.entry_1d(ws)
    assert entry["overlays"][0] == {
        "label": "compare", "color": win._overlays[0]["color"], "visible": True,
        "source": str(path), "scale": 2.5, "shift": -1.5, "yoff": 0.25}

    proj = str(tmp_path / "s.larproj.json")
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (proj, "")))
    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (proj, "")))
    win.save_project()
    win.open_project()
    ov = win._overlays[0]
    assert (ov["scale"], ov["shift"], ov["yoff"]) == (2.5, -1.5, 0.25)
    item = win.view._overlay_items[0]
    assert np.allclose(item.xData, ov["ppm"] - 1.5)
    assert item.yData.max() == pytest.approx(2.5 * float(ov["amp"].max()) + 0.25 * float(amp.max() - amp.min()), rel=1e-6)
    assert "\u00d72.5" in _row(win)["lab"].text()

    # an entry written before the keys existed restores neutral values
    old = {"kind": "1d", "title": "old", "source_path": "srcA",
           "recipe": {"nucleus": "27Al", "larmor_frequency_MHz": 130.3, "sites": []},
           "hidden": [], "exp_ppm": ppm.tolist(), "exp_amp": amp.tolist(),
           "overlays": [{"label": "old", "source": str(path), "color": "#000000",
                         "visible": True, "scale": "not a number"}]}
    win._restore_1d_entry(old, notes=[])
    ov = win._overlays[0]
    assert (ov["scale"], ov["shift"], ov["yoff"]) == (1.0, 0.0, 0.0)
    assert np.allclose(win.view._overlay_items[0].xData, ov["ppm"])


def test_panel_rows_offer_the_controls_and_update_in_place(qapp):
    from larmor.desktop.datasets import DatasetsPanel

    p = DatasetsPanel()
    got = []
    p.scale_changed.connect(lambda i, v: got.append(("scale", i, v)))
    p.shift_changed.connect(lambda i, v: got.append(("shift", i, v)))
    p.yoff_changed.connect(lambda i, v: got.append(("yoff", i, v)))
    p.reset_requested.connect(lambda i: got.append(("reset", i)))
    p.make_active.connect(lambda i: got.append(("active", i)))
    p.remove.connect(lambda i: got.append(("remove", i)))
    ovs = [{"label": "a", "color": "#e8832a", "visible": True, "source": "C:/a/1r",
            "scale": 1.0, "shift": 0.0, "yoff": 0.0},
           {"label": "b", "color": "#1f77b4", "visible": False, "source": "",
            "scale": 2.0, "shift": -3.0, "yoff": 0.5}]
    p.rebuild("active", ovs, "27Al")
    rows = p._row_widgets
    assert len(rows) == 2
    assert rows[1]["scale"].value() == 2.0 and rows[1]["shift"].value() == -3.0
    assert rows[1]["yoff"].value() == 0.5 and not rows[1]["chk"].isChecked()
    assert "\u00d72" in rows[1]["lab"].text() and "\u21910.5" in rows[1]["lab"].text()
    assert rows[0]["scale"].minimum() == 0.01 and rows[0]["scale"].maximum() == 1000.0
    assert rows[0]["scale"].prefix().startswith("\u00d7")
    rows[0]["scale"].setValue(4.0)
    rows[1]["shift"].setValue(2.0)
    rows[1]["yoff"].setValue(-0.1)
    assert got == [("scale", 0, 4.0), ("shift", 1, 2.0), ("yoff", 1, -0.1)]
    # the menu of a source-less overlay cannot make it active
    menu = p.row_menu(1)
    assert not [a for a in menu.actions() if a.text() == "Make active"][0].isEnabled()
    [a for a in menu.actions() if a.text().startswith("Reset")][0].trigger()
    [a for a in menu.actions() if a.text() == "Remove"][0].trigger()
    assert got[-2:] == [("reset", 1), ("remove", 1)]
    # the same overlays again (a colour change): rows updated in place, quietly
    got.clear()
    ovs[0]["color"] = "#123456"
    ovs[0]["scale"] = 4.0
    p.rebuild("active", ovs, "27Al")
    assert p._row_widgets[0]["scale"] is rows[0]["scale"]
    assert "#123456" in p._row_widgets[0]["swatch"].styleSheet()
    assert got == []
    # a different list: fresh rows
    p.rebuild("active", ovs[:1], "27Al")
    assert len(p._row_widgets) == 1 and p._row_widgets[0]["scale"] is not rows[0]["scale"]
    p.rebuild("active", [], "")
    assert p._row_widgets == []
