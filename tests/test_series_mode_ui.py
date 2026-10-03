"""The series mode of the workbench (wip/SEQ): Series ▸ Sequential fit… no
longer opens a dialog -- every spectrum of the series is a workspace of the
main window and a series bar above the plot walks it, carries the model,
chains Fit → next through the window's FitWorker and runs the auto sweep in
a SeqWorker. Synthetic 3-spectrum CSV series throughout; every fit thread is
pumped to completion (bounded) before a test ends."""
import csv
import json
import os
import time

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    monkeypatch.setenv("LARMOR_NO_KERNEL_WARM", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    w._add_recent = lambda p: None                 # keep the real Open recent clean
    w._report_project_notes = lambda notes: None   # never a modal box offscreen
    w._confirm_open_mode = lambda: "replace"
    yield w
    _settle(w, qapp)
    w.close()


# ------------------------------------------------------------------ helpers
def _series(tmp_path, amps=(100.0, 50.0, 80.0)):
    """Three 11B gauss_lor spectra at 13 / 15 / 17 ppm (the old
    test_seqfit_ui helper) with DIFFERENT heights, so a carried model has
    to rescale its amplitudes; plus the one-line model a user would have on
    screen."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 60, 500)
    paths = []
    for k, (pos, amp) in enumerate(zip([13.0, 15.0, 17.0], amps)):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                    sites=[SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(pos),
                        "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                        "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        p = tmp_path / f"s{k}.csv"
        with open(p, "w", encoding="utf-8") as f:
            f.write("# nucleus = 11B\n# larmor_MHz = 160\n")
            for xi, yi in zip(x, m):
                f.write(f"{xi:.4f} {yi:.4f}\n")
        paths.append(str(p))
    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 12.0, "min": 0, "max": 30},
                 "shift_fwhm_ppm": {"value": 5.0, "min": 0.1},
                 "amplitude": {"value": 80.0, "min": 0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    return paths, model


def _pump(qapp, cond, timeout=120.0, what="condition"):
    t0 = time.perf_counter()
    while not cond():
        if time.perf_counter() - t0 > timeout:
            raise AssertionError(f"timed out waiting for {what}")
        qapp.processEvents()
        time.sleep(0.01)
    qapp.processEvents()


def _fit_idle(win):
    fw = win._fit_worker
    return (fw is None or not fw.isRunning()) and win._active_fit_worker is None


def _sweep_idle(win):
    w = win._seq_worker
    return (w is None or not w.isRunning()) and win._active_fit_worker is None


def _settle(win, qapp, timeout=30.0):
    """Let the debounced simulation land and its thread finish."""
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < timeout:
        qapp.processEvents()
        sw = win._sim_worker
        if not (win._sim_timer.isActive() or win._sim_pending
                or (sw is not None and sw.isRunning())):
            break
        time.sleep(0.01)
    qapp.processEvents()


def _pos(rec):
    return rec["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"]


def _amp(rec):
    return rec["sites"][0]["params"]["amplitude"]["value"]


def _dock_titles(win):
    return [win.ws_panel.list.item(i).text() for i in range(win.ws_panel.list.count())]


def _start(win, qapp, tmp_path, with_model=True):
    paths, model = _series(tmp_path)
    win.start_series(paths, model if with_model else None)
    _settle(win, qapp)
    return paths, model


# ------------------------------------------------------------------ the bar alone
def test_series_bar_widget_alone(qapp):
    from larmor.desktop.series_bar import SeriesBar, status_color

    bar = SeriesBar()
    got = []
    bar.member_clicked.connect(got.append)
    bar.carry_changed.connect(lambda names: got.append(("carry", names)))
    bar.sweep_requested.connect(lambda p, s, m: got.append(("sweep", p, s, m)))
    bar.seed_requested.connect(lambda d, s: got.append(("seed", d, s)))
    bar.lock_toggled.connect(lambda k, on: got.append(("lock", k, on)))
    bar.remove_requested.connect(lambda k: got.append(("remove", k)))
    bar.rename_requested.connect(lambda k: got.append(("rename", k)))
    bar.set_members([{"name": "a", "status": "unfitted", "tip": "a"},
                     {"name": "b", "status": "fitted", "tip": "b", "flag": True},
                     {"name": "c", "status": "edited", "tip": "c"}])
    btns = bar.member_buttons()
    assert [b.text() for b in btns] == ["a", "b ⚠", "c"]
    assert bar.statuses() == ["unfitted", "fitted", "edited"]
    assert "fitted — click" in btns[1].toolTip() and "edited since" in btns[2].toolTip()
    bar.set_current(1)
    assert btns[1].isChecked() and bar.current() == 1 and btns[1].font().bold()
    assert bar.btnPrev.isEnabled() and bar.btnNext.isEnabled()
    assert bar.actSeedPrev.isEnabled() and bar.actSeedNext.isEnabled()
    assert not bar.actLock.isChecked()
    bar.set_current(2)
    assert not bar.btnNext.isEnabled() and not bar.actSeedNext.isEnabled()
    btns[0].click()
    assert got[-1] == 0
    # a kept member: the mark before its name, the Carry menu's keep entry
    # ticked and its seed entries off; the member menu carries the same
    bar.set_members([{"name": "a", "status": "unfitted"},
                     {"name": "b", "status": "fitted", "locked": True, "flag": True},
                     {"name": "c", "status": "edited"}])
    assert [b.text() for b in bar.member_buttons()] == ["a", "🔒 b ⚠", "c"]
    assert bar.locked() == [False, True, False]
    assert "kept" in bar.member_buttons()[1].toolTip()
    bar.set_current(1)
    assert bar.actLock.isChecked() and not bar.actSeedPrev.isEnabled()
    bar.actLock.setChecked(False)
    assert got[-1] == ("lock", 1, False)
    bar.set_current(0)
    bar.actSeedNext.trigger()
    assert got[-1] == ("seed", 0, 1)
    assert not bar.actSeedPrev.isEnabled()
    menu = bar.build_member_menu(1)
    acts = {a.text(): a for a in menu.actions() if a.text()}
    assert list(acts)[0] == "Switch to b"
    assert acts["Keep this fit  🔒"].isChecked()
    assert not acts["Seed it from the previous spectrum"].isEnabled()      # kept
    acts["Keep this fit  🔒"].setChecked(False)
    assert got[-1] == ("lock", 1, False)
    menu2 = bar.build_member_menu(2)
    acts2 = {a.text(): a for a in menu2.actions() if a.text()}
    assert not acts2["Seed it from the next spectrum"].isEnabled()         # last member
    acts2["Seed it from the previous spectrum"].trigger()
    assert got[-1] == ("seed", 2, 1)
    acts2["Seed it from the current spectrum"].trigger()
    assert got[-1] == ("seed", 2, 0)
    acts2["Remove from the series"].trigger()
    assert got[-1] == ("remove", 2)
    acts2["Rename…"].trigger()
    assert got[-1] == ("rename", 2)
    acts2["Switch to c"].trigger()
    assert got[-1] == 2
    # options: the Carry checklist with PARAM_LABELS, the sweep form
    bar.set_options(False, ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm"),
                    ("isotropic_chemical_shift_ppm",), 4, "last", 3)
    assert not bar.actSeedOnMove.isChecked()
    assert bar.carry() == ("isotropic_chemical_shift_ppm",)
    labels = [a.text() for a in bar.mCarry.actions() if a.isCheckable()]
    assert "δiso (ppm)" in labels and "FWHM (ppm)" in labels
    assert bar.sweep_settings() == (4, "last", 3)
    bar._carry_actions["shift_fwhm_ppm"].setChecked(True)
    assert got[-1] == ("carry", ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm"))
    bar._run_clicked()
    assert got[-1] == ("sweep", 4, "last", 3)
    # the sweep state greys the walk and shows progress + Stop
    bar.set_sweep_running(True)
    bar.set_sweep_progress("pass 1/2 · 2/3 · RMSD 0.01")
    assert bar.lblProgress.isVisibleTo(bar) and bar.btnStop.isVisibleTo(bar)
    assert not bar.btnFitNext.isEnabled() and not btns[0].isEnabled()
    bar.set_sweep_running(False)
    assert bar.btnFitNext.isEnabled() and btns[0].isEnabled()
    # the member count can change in place
    bar.set_members([{"name": "only", "status": "failed"}])
    assert len(bar.member_buttons()) == 1 and bar.statuses() == ["failed"]
    assert status_color("failed") != status_color("fitted")
    bar.deleteLater()


# ------------------------------------------------------------------ start
def test_start_series_tags_three_workspaces_and_shows_the_bar(win, qapp, tmp_path):
    paths, model = _start(win, qapp, tmp_path)
    spec = win._series
    assert spec is not None and spec.n == 3 and spec.names() == ["s0", "s1", "s2"]
    assert len(win.workspaces) == 3 and all(ws["kind"] == "1d" for ws in win.workspaces)
    tags = [ws["series"] for ws in win.workspaces]
    assert [t["index"] for t in tags] == [0, 1, 2] and len({t["id"] for t in tags}) == 1
    assert [t["name"] for t in tags] == ["s0", "s1", "s2"]
    assert win.active_ws == 0 and win.source_path == paths[0]
    bar = win.series_bar
    assert bar.isVisibleTo(win) and len(bar.member_buttons()) == 3
    assert bar.current() == 0 and bar.member_buttons()[0].isChecked()
    # the model on screen landed on the FIRST member only
    assert len(win.recipe["sites"]) == 1 and _pos(win.recipe) == 12.0
    assert not (win.workspaces[1]["snap"]["recipe"] or {}).get("sites")
    assert not (win.workspaces[2]["snap"]["recipe"] or {}).get("sites")
    assert bar.statuses() == ["unfitted"] * 3
    # the dock mirrors the series order
    assert [t.split("  ", 1)[1] for t in _dock_titles(win)] == ["1· s0", "2· s1", "3· s2"]
    assert win.actEndSeries.isEnabled()
    # the bar belongs to the 1D page: it hides with the co-fit page
    win.central_stack.setCurrentWidget(win.cofit_page)
    assert not bar.isVisibleTo(win)
    win.central_stack.setCurrentWidget(win.view)
    assert bar.isVisibleTo(win)


def test_series_next_copies_the_model_into_the_empty_neighbour_with_scaled_amplitudes(
        win, qapp, tmp_path):
    paths, model = _start(win, qapp, tmp_path)
    a0 = _amp(win.recipe)
    win.series_next()
    _settle(win, qapp)
    assert win.active_ws == 1 and win._series_current_index() == 1
    assert win.series_bar.current() == 1
    rec = win.recipe
    assert len(rec["sites"]) == 1 and _pos(rec) == 12.0 and rec["sites"][0]["label"] == "A"
    # member 1 is half as high as member 0: the copied amplitude is halved
    assert _amp(rec) == pytest.approx(a0 * 0.5, rel=0.02)
    assert "1 line copied from s0, amplitudes scaled ×0.50" in win.statusBar().currentMessage()
    assert win.undo_stack                              # the seed is one undo step
    win.undo()
    _settle(win, qapp)
    assert not win.recipe["sites"]
    win.redo()
    _settle(win, qapp)
    assert len(win.recipe["sites"]) == 1
    # seed on move off: moving leaves the empty member empty
    win.series_set_seed_on_move(False)
    assert win.workspaces[1]["series"]["options"]["seed_on_move"] is False
    win.series_next()
    _settle(win, qapp)
    assert win.active_ws == 2 and not win.recipe["sites"]
    win.series_next()
    assert "last spectrum" in win.statusBar().currentMessage()
    win.series_prev()
    _settle(win, qapp)
    assert win.active_ws == 1


def test_added_line_on_a_member_survives_moving_away_and_back(win, qapp, tmp_path):
    paths, model = _start(win, qapp, tmp_path)
    win.series_next()
    _settle(win, qapp)                       # member 1 took a copy (1 line)
    win._model_actions["gauss_lor"].setChecked(True)
    win.add_site_at(30.0, 5.0)
    _settle(win, qapp)
    assert len(win.recipe["sites"]) == 2
    win.series_prev()
    _settle(win, qapp)
    assert win.active_ws == 0 and len(win.recipe["sites"]) == 1   # its own one line
    win.series_next()
    _settle(win, qapp)
    assert win.active_ws == 1 and len(win.recipe["sites"]) == 2   # still two
    assert "keeps its own 2 lines" in win.statusBar().currentMessage()
    # and removing a line there is respected too
    win.remove_sites([1])
    _settle(win, qapp)
    win.series_prev()
    _settle(win, qapp)
    win.series_next()
    _settle(win, qapp)
    assert len(win.recipe["sites"]) == 1


def test_moving_never_changes_a_member_that_has_lines(win, qapp, tmp_path):
    """Sam's case: sample 1 fitted, sample 2 worked on, back to sample 1 --
    its fit must still be there (the first version re-seeded it from the
    neighbour on every move and ruined it)."""
    from larmor.seriesmode import rmsd_of

    paths, model = _start(win, qapp, tmp_path)
    win.series_fit_then_next()                     # fits s0, lands on s1 with a copy
    _pump(qapp, lambda: _fit_idle(win) and win._series_chain is None and win.active_ws == 1,
          what="Fit → next to land on member 2")
    _settle(win, qapp)
    fitted0 = json.loads(json.dumps(win.workspaces[0]["snap"]["recipe"]))
    assert rmsd_of(fitted0) is not None
    # on s1, work as a user does on sample 2: move the copied line, add one
    win.recipe["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] = 25.0
    win.on_params_changed()
    _settle(win, qapp)
    win._model_actions["gauss_lor"].setChecked(True)
    win.add_site_at(28.0, 5.0)
    _settle(win, qapp)
    assert len(win.recipe["sites"]) == 2
    # back to s0: nothing moved, the dot is still green
    win.series_prev()
    _settle(win, qapp)
    assert win.active_ws == 0
    assert win.recipe["sites"] == fitted0["sites"]
    assert rmsd_of(win.recipe) == rmsd_of(fitted0)
    assert win.series_bar.statuses()[0] == "fitted"
    assert "s0 keeps its own 1 line (fitted)" in win.statusBar().currentMessage()
    # forward again: s1 keeps its two lines where the user put them
    win.series_next()
    _settle(win, qapp)
    assert win.active_ws == 1 and len(win.recipe["sites"]) == 2 and _pos(win.recipe) == 25.0
    assert "s1 keeps its own 2 lines" in win.statusBar().currentMessage()
    assert "Carry ▾" in win.statusBar().currentMessage()        # the unfitted hint
    # an empty member still takes a copy of the model you leave
    win.series_next()
    _settle(win, qapp)
    assert win.active_ws == 2 and len(win.recipe["sites"]) == 2
    assert "2 lines copied from s1" in win.statusBar().currentMessage()


def test_seeding_from_a_neighbour_is_explicit_and_pairs_lines_by_label(win, qapp, tmp_path):
    paths, model = _start(win, qapp, tmp_path)          # s0: one line "A" at 12
    win.series_set_seed_on_move(False)
    win.series_go(1)
    _settle(win, qapp)
    own_a = json.loads(json.dumps(model["sites"][0]))   # the same component, named alike
    own_a["params"]["isotropic_chemical_shift_ppm"] = {"value": 14.0, "min": 0, "max": 14.5}
    own_a["params"]["amplitude"]["value"] = 33.0
    own_c = json.loads(json.dumps(model["sites"][0]))
    own_c["label"] = "C"                                 # a component s0 does not have
    own_c["params"]["isotropic_chemical_shift_ppm"] = {"value": 25.0, "min": 0, "max": 40}
    win.recipe["sites"] = [own_c, own_a]                 # its own order
    win.on_structure_changed()
    _settle(win, qapp)
    win.series_set_seed_on_move(True)
    win.series_go(0)
    _settle(win, qapp)
    p = win.recipe["sites"][0]["params"]
    p["isotropic_chemical_shift_ppm"]["value"] = 16.0
    p["shift_fwhm_ppm"]["value"] = 9.0
    p["amplitude"]["value"] = 999.0
    win.on_params_changed()
    _settle(win, qapp)
    win.series_go(1)                                     # a move: untouched
    _settle(win, qapp)
    assert win.recipe["sites"][1]["params"]["isotropic_chemical_shift_ppm"]["value"] == 14.0
    # Carry ▾ Seed this spectrum from the previous one: by label, clipped,
    # amplitudes and the nameless-elsewhere line untouched, one undo step
    win.series_seed(1, 0)
    _settle(win, qapp)
    rec = win.recipe
    assert [s["label"] for s in rec["sites"]] == ["C", "A"]
    c, a = rec["sites"][0]["params"], rec["sites"][1]["params"]
    assert a["isotropic_chemical_shift_ppm"]["value"] == 14.5             # clipped to ITS bound
    assert a["shift_fwhm_ppm"]["value"] == 9.0
    assert a["amplitude"]["value"] == 33.0
    assert c["isotropic_chemical_shift_ppm"]["value"] == 25.0
    msg = win.statusBar().currentMessage()
    assert "1 matching line seeded from s0" in msg and "(C has no counterpart there)" in msg
    assert win.series_bar.statuses()[1] == "unfitted"
    win.undo()
    _settle(win, qapp)
    assert win.recipe["sites"][1]["params"]["isotropic_chemical_shift_ppm"]["value"] == 14.0
    # the same action onto a member that is not on screen lands in its snapshot
    win.series_go(0)
    _settle(win, qapp)
    win.series_seed(1, 0)
    snap = win.workspaces[1]["snap"]["recipe"]
    assert snap["sites"][1]["params"]["isotropic_chemical_shift_ppm"]["value"] == 14.5
    assert win.active_ws == 0
    # no name in common: nothing moves, the status says why
    snap["sites"][1]["label"] = "Q"
    win.series_seed(1, 0)
    assert "shares a name" in win.statusBar().currentMessage()
    assert snap["sites"][1]["params"]["shift_fwhm_ppm"]["value"] == 9.0
    snap["sites"][1]["label"] = "A"
    # the Carry checklist rules the explicit seed too
    win._series_fill_carry_menu()
    cands = win.series_bar.carry_candidates()
    assert "shift_fwhm_ppm" in cands
    win.series_set_carry(tuple(n for n in cands if n != "shift_fwhm_ppm"))
    assert win._series.options.carry_off == ("shift_fwhm_ppm",)
    win.recipe["sites"][0]["params"]["shift_fwhm_ppm"]["value"] = 3.3
    win.on_params_changed()
    _settle(win, qapp)
    win.series_seed(1, 0)
    assert win.workspaces[1]["snap"]["recipe"]["sites"][1]["params"]["shift_fwhm_ppm"]["value"] == 9.0
    # seeding an EMPTY member explicitly is a whole copy
    win.series_seed(2, 0)
    assert len(win.workspaces[2]["snap"]["recipe"]["sites"]) == 1
    assert "1 line copied from s0" in win.statusBar().currentMessage()


def test_keep_this_fit_protects_a_member_from_seeding_copies_and_the_sweep(
        win, qapp, tmp_path):
    paths, model = _start(win, qapp, tmp_path)
    win.series_copy_model()
    win.series_fit_then_next()
    _pump(qapp, lambda: _fit_idle(win) and win._series_chain is None and win.active_ws == 1,
          what="Fit → next")
    _settle(win, qapp)
    fitted0 = json.loads(json.dumps(win.workspaces[0]["snap"]["recipe"]))
    win.series_set_locked(0, True)
    assert win._series.members[0].locked and win.workspaces[0]["series"]["locked"] is True
    assert win.series_bar.locked() == [True, False, False]
    assert win.series_bar.member_buttons()[0].text().startswith("🔒")
    assert "kept 🔒" in win.statusBar().currentMessage()
    # seeding into it refuses; copy-to-all walks past it
    win.series_seed(0, 1)
    assert "is kept 🔒" in win.statusBar().currentMessage()
    assert win.workspaces[0]["snap"]["recipe"] == fitted0
    win.recipe["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] = 20.0
    win.on_params_changed()
    _settle(win, qapp)
    win.series_copy_model()
    assert "1 kept 🔒 left alone" in win.statusBar().currentMessage()
    assert win.workspaces[0]["snap"]["recipe"] == fitted0
    assert _pos(win.workspaces[2]["snap"]["recipe"]) == 20.0
    # the sweep: s0 seeds s1 and is never refitted; the others are fitted
    win.series_auto_sweep(2, "first", 0)
    assert win._seq_worker is not None and win._seq_worker.fixed == (0,)
    assert "1 kept as seed" in win.statusBar().currentMessage()
    _pump(qapp, lambda: _sweep_idle(win), what="the auto sweep")
    _settle(win, qapp)
    assert win.workspaces[0]["snap"]["recipe"] == fitted0
    assert win.series_bar.statuses() == ["fitted"] * 3
    pos = [_pos(win._series_recipe_at(i)) for i in range(3)]
    assert pos == pytest.approx([13.0, 15.0, 17.0], abs=0.4)
    msg = win.statusBar().currentMessage()
    assert "3 spectra, 1 kept" in msg and "1 kept spectrum untouched" in msg
    # release, and a series kept whole cannot sweep
    win.series_set_locked(0, False)
    assert not win._series.members[0].locked and "released" in win.statusBar().currentMessage()
    for k in range(3):
        win.series_set_locked(k, True)
    w = win._seq_worker
    win.series_auto_sweep(2, "first", 0)
    assert win._seq_worker is w and "every spectrum is kept" in win.statusBar().currentMessage()
    # the lock survives a project round trip
    win._sync_active()
    from larmor import project
    bundle, _dropped = project.build_bundle(win.workspaces, win.active_ws, str(tmp_path))
    assert all(w_["series"]["locked"] is True for w_ in bundle["workspaces"])


def test_member_menu_removes_a_spectrum_from_the_series(win, qapp, tmp_path):
    paths, model = _start(win, qapp, tmp_path)
    menu = win.series_bar.build_member_menu(1)
    acts = {a.text(): a for a in menu.actions() if a.text()}
    assert list(acts)[0] == "Switch to s1"
    acts["Keep this fit  🔒"].setChecked(True)
    assert win._series.members[1].locked
    acts["Remove from the series"].trigger()
    qapp.processEvents()
    assert win._series.n == 2 and win._series.names() == ["s0", "s2"]
    assert len(win.workspaces) == 3 and "series" not in win.workspaces[1]
    assert [t.split("  ", 1)[1] for t in _dock_titles(win)] == ["1· s0", "s1", "2· s2"]
    assert [b.text() for b in win.series_bar.member_buttons()] == ["s0", "s2"]
    assert "left the series" in win.statusBar().currentMessage()
    # the walk follows the shorter series
    win.series_next()
    _settle(win, qapp)
    assert win.active_ws == 2
    # the last members out end the series
    win.series_remove_member(1)
    win.series_remove_member(0)
    assert win._series is None and not win.series_bar.isVisibleTo(win)
    assert len(win.workspaces) == 3 and not any("series" in ws for ws in win.workspaces)
    assert "series ended" in win.statusBar().currentMessage()


def test_series_bar_paints_the_current_member_readable_and_cascades_no_bare_style(
        win, qapp, tmp_path):
    """The screenshot bugs: the strip body's bare 'background: transparent'
    cascaded onto the checked member (accent background gone, white name on
    white) and onto the member tool tips (a black box). Every style sheet
    on the way down to a member button now carries a selector, and under
    the application theme (applied here the way main() does; MainWindow
    itself does not) the checked button renders in the accent colour with
    readable text. Measured offscreen on the Light theme: 0 accent pixels
    and 3913 transparent ones with the bare sheet, 3533 accent pixels of
    4225 with the fix."""
    from PySide6.QtGui import QColor
    from larmor.desktop import theme

    paths, model = _start(win, qapp, tmp_path)
    bar = win.series_bar
    btn = bar.member_buttons()[0]
    w = btn
    while w is not None:
        ss = w.styleSheet().strip()
        assert not ss or "{" in ss.split(";")[0], (w.objectName() or type(w).__name__, ss)
        if w is bar:
            break
        w = w.parentWidget()
    old_ss, old_pal = qapp.styleSheet(), qapp.palette()
    theme.apply(qapp, theme.DEFAULT)
    win.show()
    qapp.processEvents()
    try:
        img = btn.grab().toImage()
        t = theme.active()
        acc, txt = QColor(t.accent), QColor(t.accent_text)

        def near(c, ref, tol):
            return (abs(c.red() - ref.red()) + abs(c.green() - ref.green())
                    + abs(c.blue() - ref.blue())) <= tol

        n_acc = n_txt = n_clear = 0
        for y in range(img.height()):
            for x in range(img.width()):
                c = img.pixelColor(x, y)
                if c.alpha() == 0:
                    n_clear += 1
                elif near(c, acc, 40):
                    n_acc += 1
                elif near(c, txt, 60):
                    n_txt += 1
        area = img.width() * img.height()
        assert n_clear == 0, n_clear                    # nothing shows through
        assert n_acc > 0.3 * area, (n_acc, area)
        assert n_txt > 8, n_txt
    finally:
        win.hide()
        qapp.setStyleSheet(old_ss)
        qapp.setPalette(old_pal)


# ------------------------------------------------------------------ fitting
def test_fit_then_next_runs_the_fit_worker_and_lands_on_the_next_member_fitted(
        win, qapp, tmp_path):
    from larmor.seriesmode import rmsd_of

    paths, model = _start(win, qapp, tmp_path)
    assert win.series_bar.statuses()[0] == "unfitted"
    win.series_fit_then_next()
    assert win._fit_worker is not None and win._series_chain is win._fit_worker
    _pump(qapp, lambda: _fit_idle(win) and win._series_chain is None and win.active_ws == 1,
          what="Fit → next to land on member 2")
    _settle(win, qapp)
    rec0 = win.workspaces[0]["snap"]["recipe"]
    assert rmsd_of(rec0) is not None and _pos(rec0) == pytest.approx(13.0, abs=0.4)
    assert win.series_bar.statuses()[0] == "fitted"
    assert win.workspaces[0]["series"]["fit_sig"]
    # member 2 was seeded from the FITTED member 1 (a copy: it had no lines)
    assert len(win.recipe["sites"]) == 1 and _pos(win.recipe) == pytest.approx(13.0, abs=0.4)
    assert win.series_bar.statuses()[1] == "unfitted"
    assert win.statusBar().currentMessage().startswith("fitted s0 (RMSD")
    # an edit on a fitted member turns its dot amber; undo turns it green again
    win.series_set_seed_on_move(False)
    win.series_go(0)
    _settle(win, qapp)
    assert win.series_bar.statuses()[0] == "fitted"
    win.snapshot()
    win.recipe["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] += 1.0
    win.on_params_changed()
    _settle(win, qapp)
    assert win.series_bar.statuses()[0] == "edited"
    win.undo()
    _settle(win, qapp)
    assert win.series_bar.statuses()[0] == "fitted"
    # the chain from the LAST member fits it and stops there
    win.series_go(2)
    _settle(win, qapp)
    win.series_copy_model()                       # nothing to copy: no lines here
    assert "no lines to copy" in win.statusBar().currentMessage()
    win.series_set_seed_on_move(True)
    win.series_go(1)                              # member 1 keeps its lines (source has none)
    _settle(win, qapp)
    assert len(win.recipe["sites"]) == 1
    win.series_fit_then_next()
    _pump(qapp, lambda: _fit_idle(win) and win._series_chain is None and win.active_ws == 2,
          what="Fit → next to land on member 3")
    _settle(win, qapp)
    win.series_fit_then_next()
    _pump(qapp, lambda: _fit_idle(win) and win._series_chain is None,
          what="the last member's fit")
    _settle(win, qapp)
    assert win.active_ws == 2 and "last spectrum" in win.statusBar().currentMessage()
    assert win.series_bar.statuses() == ["fitted", "fitted", "fitted"]
    pos = [_pos(win._series_recipe_at(i)) for i in range(3)]
    assert pos == pytest.approx([13.0, 15.0, 17.0], abs=0.4)


def test_auto_sweep_fits_every_member_in_a_worker(win, qapp, tmp_path):
    from larmor.seriesmode import rmsd_of

    paths, model = _start(win, qapp, tmp_path)
    # a sweep needs a line on every spectrum: the status says so
    win.series_auto_sweep(2, "first", 0)
    assert "needs at least one line" in win.statusBar().currentMessage()
    assert win._seq_worker is None
    win.series_copy_model()
    assert "copied into 2 spectra" in win.statusBar().currentMessage()
    assert all(len(ws["snap"]["recipe"]["sites"]) == 1 for ws in win.workspaces[1:])
    win.series_auto_sweep(2, "first", 0)
    assert win._seq_worker is not None and win.series_bar.is_running()
    assert win._active_fit_worker is win._seq_worker
    assert win.btnStopFit.isVisibleTo(win) and not win.lines_table.btnFit.isEnabled()
    _pump(qapp, lambda: _sweep_idle(win), what="the auto sweep")
    _settle(win, qapp)
    assert not win.series_bar.is_running() and win.lines_table.btnFit.isEnabled()
    assert win.series_bar.statuses() == ["fitted"] * 3
    recs = [win._series_recipe_at(i) for i in range(3)]
    assert [_pos(r) for r in recs] == pytest.approx([13.0, 15.0, 17.0], abs=0.4)
    assert all(rmsd_of(r) is not None and np.isfinite(rmsd_of(r)) for r in recs)
    msg = win.statusBar().currentMessage()
    assert "sequential fit: 3 spectra, 2 passes" in msg and "pass means" in msg
    # the on-screen member got a fit verdict, the others carry theirs in their snapshot
    assert win._health is not None and win._health.fitted and not win._health.stale
    assert all(win.workspaces[i]["snap"]["health"] is not None for i in (1, 2))
    assert win.workspaces[0]["series"]["options"]["passes"] == 2
    assert win.undo_stack                          # the on-screen member can undo the sweep
    # switching to a swept member shows its fitted verdict
    win.series_set_seed_on_move(False)
    win.series_go(2)
    _settle(win, qapp)
    assert win._health is not None and win._health.fitted
    assert win.series_bar.statuses() == ["fitted"] * 3


# ------------------------------------------------------------------ table / plot / outputs
def test_series_table_reorder_renames_and_reorders_members_and_dock_titles(
        win, qapp, tmp_path, monkeypatch):
    paths, model = _start(win, qapp, tmp_path, with_model=False)
    tbl = win._series.series_table().permuted([2, 1, 0])
    tbl.rename(0, "last")
    win._series_apply_table(tbl, [2, 1, 0])
    spec = win._series
    assert spec.names() == ["last", "s1", "s0"]
    assert spec.paths() == [paths[2], paths[1], paths[0]]
    assert win.workspaces[2]["series"]["index"] == 0 and win.workspaces[2]["series"]["name"] == "last"
    assert win.workspaces[0]["series"]["index"] == 2
    assert win.workspaces[2]["title"] == "last"
    assert win.workspaces[2]["snap"]["recipe"]["sample"] == "last"
    titles = [t.split("  ", 1)[1] for t in _dock_titles(win)]
    assert titles == ["3· s0", "2· s1", "1· last"]
    assert [b.text() for b in win.series_bar.member_buttons()] == ["last", "s1", "s0"]
    assert "new order" in win.statusBar().currentMessage()
    # the walk follows the new order
    win.series_go(0)
    _settle(win, qapp)
    assert win.active_ws == 2 and win.recipe["sample"] == "last"
    win.series_next()
    _settle(win, qapp)
    assert win.active_ws == 1
    # a Rename… in the Workspaces dock reaches the member and the bar
    win.set_workspace_title(1, "middle")
    assert spec.names() == ["last", "middle", "s0"]
    assert [b.text() for b in win.series_bar.member_buttons()] == ["last", "middle", "s0"]
    assert win.workspaces[1]["series"]["name"] == "middle"
    # the dialog route: OK with nothing changed keeps everything
    from PySide6.QtWidgets import QDialog
    from larmor.desktop import series_table_dialog as std
    monkeypatch.setattr(std.SeriesTableDialog, "exec", lambda self: QDialog.Accepted)
    win.series_table()
    assert spec.names() == ["last", "s1", "s0"] and "series table applied" in \
        win.statusBar().currentMessage()
    # the series survives a project round trip in its new order and names
    win._sync_active()
    from larmor import project
    bundle, dropped = project.build_bundle(win.workspaces, win.active_ws, str(tmp_path))
    assert dropped == 0 and [w["series"]["index"] for w in bundle["workspaces"]] == [2, 1, 0]


def test_series_plot_and_acquisition_table_open_as_tool_windows(win, qapp, tmp_path, monkeypatch):
    from larmor.desktop import mw_series
    from larmor.desktop.acquisition_dialog import AcquisitionTableDialog
    from larmor.desktop.series_plot import SeriesPlotDialog

    opened = []
    monkeypatch.setattr(mw_series, "show_tool_window", lambda dlg, owner=None: opened.append(dlg))
    paths, model = _start(win, qapp, tmp_path)
    win.series_plot()                                 # only member 0 has lines
    assert isinstance(opened[-1], SeriesPlotDialog) and opened[-1]._labels == ["s0"]
    assert "left out" in win.statusBar().currentMessage()
    win.series_copy_model()
    win.series_plot()
    assert isinstance(opened[-1], SeriesPlotDialog) and opened[-1]._labels == ["s0", "s1", "s2"]
    win.series_acquisition_table()
    assert isinstance(opened[-1], AcquisitionTableDialog)
    for d in opened:
        d.close()
        d.deleteLater()
    qapp.processEvents()


def test_save_all_fits_writes_one_recipe_per_member(win, qapp, tmp_path, monkeypatch):
    import re
    from larmor.recipe import Recipe

    paths, model = _start(win, qapp, tmp_path)
    win.series_copy_model()
    out = tmp_path / "fits"
    out.mkdir()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *a, **k: str(out)))
    monkeypatch.setattr(QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.Yes))
    win.series_save_all()
    files = sorted(out.glob("*.recipe.json"))
    assert len(files) == 3
    recs = [Recipe.load(f) for f in files]
    assert sorted(r.source_path for r in recs) == sorted(paths)
    assert all(len(r.sites) == 1 for r in recs)
    assert all(re.fullmatch(r"s\d_11B_seq_\d{8}_\d{4}", f.name[:-len(".recipe.json")])
               for f in files)
    assert [r.sample for r in recs] == ["s0", "s1", "s2"]
    assert "saved 3 fit(s)" in win.statusBar().currentMessage()
    assert "3 of them not fitted yet" in win.statusBar().currentMessage()


def test_publication_bundle_writes_manifest_after_a_manual_fit(win, qapp, tmp_path, monkeypatch):
    from larmor.recipe import Recipe

    paths, model = _start(win, qapp, tmp_path)
    win.series_copy_model()
    win.series_fit_then_next()                        # member 0 fitted, now on member 1
    _pump(qapp, lambda: _fit_idle(win) and win._series_chain is None and win.active_ws == 1,
          what="Fit → next")
    _settle(win, qapp)
    folder = tmp_path / "bundle"
    folder.mkdir()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *a, **k: str(folder)))
    win.series_bundle()
    names = {p.name for p in folder.iterdir()}
    assert {"seq_table.csv", "manifest.csv", "README.txt"} <= names
    assert len(list(folder.glob("*.recipe.json"))) == 3
    assert len(list(folder.glob("*_curves.csv"))) == 3
    with open(folder / "manifest.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    for k, r in enumerate(rows):
        assert r["source_path"] == paths[k]
        assert Recipe.load(folder / r["recipe_file"]).source_path == paths[k]
    assert rows[0]["note"] == "" and float(rows[0]["rmsd"]) > 0
    assert [r["note"] for r in rows[1:]] == ["not fitted", "not fitted"]
    assert "manifest.csv" in win.statusBar().currentMessage()


# ------------------------------------------------------------------ end / entry points / project
def test_end_series_hides_the_bar_and_keeps_the_workspaces(win, qapp, tmp_path):
    paths, model = _start(win, qapp, tmp_path)
    win.end_series()
    assert win._series is None and not win.series_bar.isVisibleTo(win)
    assert len(win.workspaces) == 3 and not any("series" in ws for ws in win.workspaces)
    assert not win.actEndSeries.isEnabled()
    assert [t.split("  ", 1)[1] for t in _dock_titles(win)] == ["s0", "s1", "s2"]
    assert "spectra stay open" in win.statusBar().currentMessage()
    # the walk slots are inert without a series
    win.series_next()
    win.series_fit_then_next()
    assert win._fit_worker is None


def test_run_seq_fit_uses_the_open_spectra_when_fewer_than_two_paths_are_given(
        win, qapp, tmp_path):
    paths, model = _series(tmp_path)
    for p in paths:
        win.load_source(p, keep_fit=False)
        qapp.processEvents()
    assert len(win.workspaces) == 3 and win._series is None
    win.run_seq_fit()                          # nothing selected in the Explorer
    _settle(win, qapp)
    assert win._series is not None and win._series.n == 3
    assert win._series.paths() == paths and len(win.workspaces) == 3   # adopted, not reloaded
    assert win.active_ws == 0 and win.series_bar.isVisibleTo(win)
    # the inventory / Explorer hand-off: paths already open are adopted too
    win.end_series()
    win.run_seq_fit(paths[1:])
    _settle(win, qapp)
    assert win._series.n == 2 and win._series.paths() == paths[1:]
    assert len(win.workspaces) == 3
    assert win.workspaces[0].get("series") is None
    assert [t.split("  ", 1)[1] for t in _dock_titles(win)] == ["s0", "1· s1", "2· s2"]
    # Series ▸ Series from open spectra takes every open 1D spectrum
    win.end_series()
    win.start_series_from_workspaces()
    assert win._series.n == 3
    # a closed member drops out of the series; the rest re-index
    win.close_workspace(1)
    qapp.processEvents()
    assert win._series.n == 2 and win._series.names() == ["s0", "s2"]
    assert [t.split("  ", 1)[1] for t in _dock_titles(win)] == ["1· s0", "2· s2"]


def test_project_round_trip_restores_the_series(win, qapp, tmp_path):
    from larmor import project

    paths, model = _start(win, qapp, tmp_path)
    win.series_copy_model()
    tbl = win._series.series_table()
    tbl.rename(1, "mid")
    win._series_apply_table(tbl, [0, 1, 2])
    win.series_set_carry(("isotropic_chemical_shift_ppm",))
    win._series.options.passes = 4
    win._series_retag()
    sid = win._series.id
    win._sync_active()
    bundle, dropped = project.build_bundle(win.workspaces, win.active_ws, str(tmp_path))
    assert dropped == 0 and all(w.get("series", {}).get("id") == sid for w in bundle["workspaces"])
    proj = tmp_path / "series.larproj.json"
    proj.write_text(json.dumps(bundle), encoding="utf-8")

    win._open_project_path(str(proj))          # replace mode (stubbed)
    _settle(win, qapp)
    spec = win._series
    assert spec is not None and spec.id == sid and spec.n == 3
    assert spec.names() == ["s0", "mid", "s2"]
    assert spec.paths() == paths
    assert spec.options.passes == 4 and "shift_fwhm_ppm" in spec.options.carry_off
    assert len(win.workspaces) == 3
    assert [ws["series"]["index"] for ws in win.workspaces] == [0, 1, 2]
    assert all(len(ws["snap"]["recipe"]["sites"]) == 1 for ws in win.workspaces)
    assert win.series_bar.isVisibleTo(win) and len(win.series_bar.member_buttons()) == 3
    assert [b.text() for b in win.series_bar.member_buttons()] == ["s0", "mid", "s2"]
    assert [t.split("  ", 1)[1] for t in _dock_titles(win)] == ["1· s0", "2· mid", "3· s2"]
    assert win.actEndSeries.isEnabled()
    # and the walk works on the restored series
    win.series_go(2)
    _settle(win, qapp)
    assert win.active_ws == 2 and win.series_bar.current() == 2
