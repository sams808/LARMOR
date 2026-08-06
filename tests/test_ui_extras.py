"""Regressions for the small UI utilities added for the workflow tweaks:
the scroll-nudge opt-in flag and the Explorer's pdata-proc / fits browsing.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication, QTreeWidgetItem  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def test_scroll_nudge_flag_defaults_off_and_toggles(qapp):
    from larmor.desktop import table
    assert table.scroll_nudge_enabled() is False        # off by default
    table.set_scroll_nudge(True)
    assert table.scroll_nudge_enabled() is True
    table.set_scroll_nudge(False)
    assert table.scroll_nudge_enabled() is False


def test_list_fits_finds_recipes_and_dmfit(tmp_path):
    from larmor.desktop import explorer
    (tmp_path / "a.recipe.json").write_text("{}")
    (tmp_path / "b.fxml").write_text("<x/>")
    (tmp_path / "c.fxmla").write_text("<x/>")
    (tmp_path / "procs").write_text("junk")            # not a fit
    (tmp_path / "1r").write_bytes(b"\0")               # not a fit
    names = [f.name for f in explorer._list_fits(tmp_path)]
    assert names == ["a.recipe.json", "b.fxml", "c.fxmla"]


def test_explorer_proc_and_fit_layers(qapp, tmp_path):
    from larmor.desktop import explorer
    expno = tmp_path / "1118"
    for proc, fits in (("1", []), ("15", ["x.recipe.json", "y.fxml"])):
        d = expno / "pdata" / proc
        d.mkdir(parents=True)
        (d / "1r").write_bytes(b"\0")
        for fn in fits:
            (d / fn).write_text("{}")

    panel = explorer.ExplorerPanel()
    exp = QTreeWidgetItem(["1118"])
    exp.setData(0, explorer._ROLE_PATH, str(expno))
    exp.setData(0, explorer._ROLE_KIND, "exp")
    panel._populate_procs(exp)
    labels = [exp.child(i).text(0) for i in range(exp.childCount())]
    assert any("proc 1" in t for t in labels)
    assert any("proc 15" in t and "2 fit" in t for t in labels)

    proc15 = next(exp.child(i) for i in range(exp.childCount())
                  if "15" in exp.child(i).text(0))
    # opening a proc opens its 1r (the processing to fit on)
    assert proc15.data(0, explorer._ROLE_OPEN).endswith("1r")
    proc15.takeChildren()
    panel._populate_fits(proc15)
    fit_names = [proc15.child(i).text(0) for i in range(proc15.childCount())]
    assert any("x.recipe.json" in t for t in fit_names)
    # a fit is openable (double-click opens it)
    assert proc15.child(0).data(0, explorer._ROLE_OPEN).endswith(".recipe.json")


def test_fit_can_be_interrupted_keeping_last_values(qapp):
    """request_stop makes the iter_cb abort lmfit; the worker still returns a
    result and reports the stop mode (so the app can keep the latest values)."""
    from larmor.desktop.app import FitWorker, _emit_progress
    from larmor.engine import make_context, simulate_site
    from larmor.recipe import Recipe, SiteModel, Param

    # the iter_cb returns True (abort) exactly when should_stop() is truthy
    assert _emit_progress(_FakeSig(), lambda: True)(None, 1, np.zeros(3)) is True
    assert _emit_progress(_FakeSig(), lambda: False)(None, 1, np.zeros(3)) is None

    x = np.linspace(-20, 120, 200)
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="gauss_lor", label="p", params={
            "isotropic_chemical_shift_ppm": Param(60, min=0, max=120),
            "shift_fwhm_ppm": Param(8, min=1, max=40),
            "gl": Param(0.5, min=0, max=1, vary=False),
            "amplitude": Param(1e6, min=0)})])
    ctx = make_context(r, exp_ppm=x)
    y = np.sum([simulate_site(s, ctx) for s in r.sites], axis=0)
    fw = FitWorker(r.to_dict(), x, y * 1.3, (120, -20))
    got = {}
    fw.done.connect(lambda res, mode: got.update(res=res, mode=mode))
    fw.request_stop("stop")            # abort on the first iteration
    fw.run()
    assert got.get("mode") == "stop"
    assert got.get("res") is not None  # last-iteration result is returned


class _FakeSig:
    def emit(self, *a):
        pass


def test_emit_progress_dmfit_style_convergence(qapp):
    """The completion threshold stops the fit once the residual stdev stops
    changing by more than the threshold (dmfit 'sdev not changing' criterion)."""
    from larmor.desktop.app import _emit_progress
    cb = _emit_progress(_FakeSig(), lambda: False, converge_frac=1e-3)  # 0.1%
    assert cb(None, 1, np.full(100, 10.0)) is None       # first iteration
    assert cb(None, 2, np.full(100, 9.98)) is None        # Δ 0.2% > 0.1% → keep going
    assert cb(None, 3, np.full(100, 9.9795)) is True      # Δ ~0.005% < 0.1% → stop
    # with no threshold it never converges on its own
    cb2 = _emit_progress(_FakeSig(), lambda: False, converge_frac=None)
    assert cb2(None, 1, np.full(100, 10.0)) is None
    assert cb2(None, 2, np.full(100, 10.0)) is None


def test_batch_fit_dialog_loads_grid_and_fits(qapp, tmp_path):
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import engine
    from larmor.desktop.batchfit_dialog import BatchFitDialog, _BatchWorker

    x = np.linspace(-20, 60, 600)
    paths = []
    for k, (sh, amp) in enumerate(((0.0, 100), (0.3, 70))):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                    sites=[SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(15.0 + sh),
                        "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                        "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        d = m + np.random.default_rng(k).normal(0, 1.5, x.size)
        p = tmp_path / f"s{k}.csv"
        with open(p, "w", encoding="utf-8") as f:
            f.write("# nucleus = 11B\n# larmor_MHz = 160\n")
            for xi, yi in zip(x, d):
                f.write(f"{xi:.4f} {yi:.4f}\n")
        paths.append(str(p))
    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 14.0, "min": 0, "max": 30},
                 "shift_fwhm_ppm": {"value": 5.0, "min": 0.1},
                 "amplitude": {"value": 80.0, "min": 0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    dlg = BatchFitDialog(None, paths, model)
    assert len(dlg._data) == 2 and dlg.tabs.count() == 1
    assert "isotropic_chemical_shift_ppm" in dlg._rel_checks
    # error-calculation menu offers the methods we have, disabled until a fit
    assert [dlg.errCombo.itemData(i) for i in range(dlg.errCombo.count())] == \
        ["covariance", "montecarlo", "profile"]
    assert not dlg.btnErr.isEnabled() and not dlg.btnErrCsv.isEnabled()

    w = _BatchWorker(dlg._entries(), (), 0.1)
    w.done.connect(dlg._done)
    w.run()                                   # synchronous
    assert dlg._result is not None
    assert dlg.btnSave.isEnabled()
    assert dlg._cells[0]["model"].xData is not None
    # error buttons enable after the fit; covariance CSV exports with the method
    assert dlg.btnErr.isEnabled() and dlg.btnErrCsv.isEnabled()
    out = tmp_path / "errs.csv"
    dlg.errCombo.setCurrentIndex(0)           # covariance
    dlg._write_err_csv(str(out))
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].split(",") == ["scope", "site", "label", "param", "value",
                                   "stderr", "error_method", "sigma_pct",
                                   "ci68_lo", "ci68_hi"]
    amp = [ln for ln in lines[1:] if ln.split(",")[3] == "amplitude"]
    assert amp and all(ln.split(",")[6] == "covariance" for ln in amp)


def test_batch_dialog_per_spectrum_twopoint_baseline(qapp, tmp_path):
    """Right-click 'Add 2-point linear baseline' on one cell: two picks
    subtract the line through them for THAT spectrum only, the correction is
    recorded on its recipe (processing + source_path) so it survives the fit
    and export, and is reproducible on reload; other spectra are untouched."""
    import pyqtgraph as pg
    from PySide6.QtCore import Qt
    from larmor import batchfit
    from larmor.desktop.batchfit_dialog import BatchFitDialog
    from larmor.loader import load_any
    from larmor.recipe import Recipe

    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 14.0},
                 "shift_fwhm_ppm": {"value": 5.0},
                 "amplitude": {"value": 80.0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    dlg = BatchFitDialog(None, [], model)

    x = np.linspace(-20, 60, 300)
    # spectrum 0: a real file on disk (source for the replay check) with a
    # tilted background; spectrum 1: clean, must stay untouched
    y0 = 0.4 * x + 5.0 + 50 * np.exp(-0.5 * ((x - 15) / 6) ** 2)
    p0 = tmp_path / "s0.csv"
    p0.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                  "\n".join(f"{xi:.6f} {yi:.6f}" for xi, yi in zip(x, y0)))
    y1 = 50 * np.exp(-0.5 * ((x - 15) / 6) ** 2)
    dlg._data = [
        {"ppm": x.copy(), "amp": y0.copy(), "amp0": y0.copy(), "nucleus": "11B",
         "larmor": 160.0, "spin": 0.0, "sample": "s0", "path": str(p0),
         "proc": "", "snr": 50, "baseline_ops": []},
        {"ppm": x.copy(), "amp": y1.copy(), "amp0": y1.copy(), "nucleus": "11B",
         "larmor": 160.0, "spin": 0.0, "sample": "s1", "path": "s1.csv",
         "proc": "", "snr": 50, "baseline_ops": []},
    ]
    dlg._cells = [{"plot": pg.PlotWidget(), "exp": None, "model": None,
                  "rmsd": None, "comp": [], "title": None,
                  "bl_picking": False, "bl_markers": [], "bl_line": None}
                 for _ in range(2)]
    for c in dlg._cells:
        c["exp"] = c["plot"].plot([], [])

    # right-click cancel BEFORE two points are placed leaves the spectrum untouched
    dlg._start_bg_pick(0)
    assert dlg._cells[0]["bl_picking"]

    class _RightClick:
        def button(self): return Qt.RightButton
        def scenePos(self): return None
        def accept(self): pass

    dlg._cell_clicked(0, _RightClick())
    assert not dlg._cells[0]["bl_picking"]
    assert np.allclose(dlg._data[0]["amp"], y0)

    # place two points (a WRONG second one), verify it does NOT auto-apply —
    # the whole point is that a bad click must be fixable, not committed instantly
    dlg._start_bg_pick(0)
    cell0 = dlg._cells[0]
    for pos in [(-18.0, 0.4 * -18.0 + 5.0), (58.0, 999.0)]:   # 2nd point is wrong
        m = pg.TargetItem(pos=pos, movable=True)
        m.sigPositionChanged.connect(lambda *_: dlg._update_bg_preview(0))
        cell0["plot"].addItem(m)
        cell0["bl_markers"].append(m)
    dlg._update_bg_preview(0)
    assert cell0["bl_picking"]                          # still armed — not applied
    assert np.allclose(dlg._data[0]["amp"], y0)          # data untouched so far
    assert cell0["bl_line"] is not None                  # live preview shown

    # fix the bad point by DRAGGING it (what the user asked for) instead of
    # having to cancel and restart from scratch
    cell0["bl_markers"][1].setPos(58.0, 0.4 * 58.0 + 5.0)
    assert dlg._bg_points(0)[1] == pytest.approx((58.0, 0.4 * 58.0 + 5.0))

    dlg._apply_bg_pick(0)                                # the explicit confirm step
    assert not dlg._cells[0]["bl_picking"]

    d0, d1 = dlg._data
    edge = np.concatenate([d0["amp"][:15], d0["amp"][-15:]])
    assert abs(float(np.mean(edge))) < 0.5           # tilt removed on s0
    assert d0["baseline_ops"] == [
        {"op": "twopoint_bg", "x1": -18.0, "y1": pytest.approx(-2.2),
         "x2": 58.0, "y2": pytest.approx(28.2)}]
    assert np.allclose(d1["amp"], y1)                # s1 untouched
    assert d1["baseline_ops"] == []

    # the correction is carried into the recipe fed to the fit
    entries = dlg._entries()
    rec0, rec1 = entries[0][0], entries[1][0]
    assert rec0.processing == d0["baseline_ops"]
    assert rec0.source_path == str(p0)
    assert rec1.processing == []

    res = batchfit.batch_fit(entries)
    assert res.recipes[0].processing == d0["baseline_ops"]   # survives the fit

    # saved fit reproduces the corrected spectrum on reload (closed-loop)
    out = tmp_path / "s0.recipe.json"
    Recipe.from_dict(res.recipes[0].to_dict()).save(out)
    _, amp_reloaded, _, _, warns = load_any(str(out))
    assert any("processing step" in w for w in warns)
    assert np.allclose(np.sort(amp_reloaded), np.sort(d0["amp"]), atol=1e-6)

    # "Clear this spectrum's baseline" restores the raw spectrum
    dlg._clear_cell_baseline(0)
    assert d0["baseline_ops"] == []
    assert np.allclose(d0["amp"], y0)


def test_batch_baseline_menu_survives_apply_and_allows_a_second_one(qapp):
    """The right-click baseline options must NOT disappear after applying one
    (the bug report this guards) -- pyqtgraph rebuilds a bare default menu
    every time setMenuEnabled(True) runs, so the cell's custom items (Export /
    Send to studio / the two baseline actions) must be re-attached each time
    picking ends. A second correction must then be addable (and compose with
    the first) or cancellable, not just the very first."""
    import pyqtgraph as pg
    from larmor.desktop.batchfit_dialog import BatchFitDialog

    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 14.0},
                 "shift_fwhm_ppm": {"value": 5.0}, "amplitude": {"value": 80.0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    dlg = BatchFitDialog(None, [], model)
    x = np.linspace(-20, 60, 300)
    y = 0.4 * x + 5.0 + 50 * np.exp(-0.5 * ((x - 15) / 6) ** 2)
    dlg._data = [{"ppm": x.copy(), "amp": y.copy(), "amp0": y.copy(),
                 "nucleus": "11B", "larmor": 160.0, "spin": 0.0, "sample": "s0",
                 "path": "s0.csv", "proc": "", "snr": 50, "baseline_ops": []}]
    plot = pg.PlotWidget()
    dlg._cells = [{"plot": plot, "exp": plot.plot([], []), "model": None,
                  "rmsd": None, "comp": [], "title": None,
                  "bl_picking": False, "bl_markers": [], "bl_line": None}]
    dlg._attach_cell_menu(0)
    vb = plot.getPlotItem().getViewBox()
    expected = [a.text() for a in vb.menu.actions()]
    assert "Add 2-point linear baseline" in expected
    assert "Clear this spectrum's baseline" in expected

    def pick(pts):
        dlg._start_bg_pick(0)
        for p in pts:
            m = pg.TargetItem(pos=p, movable=True)
            plot.addItem(m)
            dlg._cells[0]["bl_markers"].append(m)
        dlg._apply_bg_pick(0)

    pick([(-18.0, 0.4 * -18.0 + 5.0), (58.0, 0.4 * 58.0 + 5.0)])
    assert [a.text() for a in vb.menu.actions()] == expected   # still there
    assert len(dlg._data[0]["baseline_ops"]) == 1

    # right-click again works: a second correction composes with the first
    pick([(-19.0, -0.1), (59.0, 0.1)])
    assert [a.text() for a in vb.menu.actions()] == expected   # still there
    assert len(dlg._data[0]["baseline_ops"]) == 2


def test_batch_baseline_right_click_confirm_apply_and_cancel(qapp, monkeypatch):
    """Once both points are down, right-click offers Apply/Cancel — it must NOT
    silently commit (the bug report this guards): Cancel discards the pick and
    leaves the spectrum untouched; Apply commits it."""
    import pyqtgraph as pg
    from PySide6.QtCore import Qt
    from larmor.desktop.batchfit_dialog import BatchFitDialog

    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 14.0},
                 "shift_fwhm_ppm": {"value": 5.0}, "amplitude": {"value": 80.0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    dlg = BatchFitDialog(None, [], model)
    x = np.linspace(-20, 60, 200)
    y = 0.4 * x + 5.0
    dlg._data = [{"ppm": x.copy(), "amp": y.copy(), "amp0": y.copy(),
                 "nucleus": "11B", "larmor": 160.0, "spin": 0.0, "sample": "s0",
                 "path": "s0.csv", "proc": "", "snr": 50, "baseline_ops": []}]
    plot = pg.PlotWidget()
    dlg._cells = [{"plot": plot, "exp": plot.plot([], []), "model": None,
                  "rmsd": None, "comp": [], "title": None,
                  "bl_picking": False, "bl_markers": [], "bl_line": None}]

    class _RightClick:
        def button(self): return Qt.RightButton
        def scenePos(self):
            from PySide6.QtCore import QPointF
            return QPointF(0.0, 0.0)
        def accept(self): pass

    def _place_two_points():
        dlg._start_bg_pick(0)
        for pos in [(-18.0, 0.4 * -18.0 + 5.0), (58.0, 0.4 * 58.0 + 5.0)]:
            m = pg.TargetItem(pos=pos, movable=True)
            plot.addItem(m)
            dlg._cells[0]["bl_markers"].append(m)

    # Cancel: discards the pick, spectrum stays raw
    _place_two_points()
    monkeypatch.setattr(dlg, "_ask_apply_or_cancel", lambda plot, pos: "cancel")
    dlg._cell_clicked(0, _RightClick())
    assert not dlg._cells[0]["bl_picking"]
    assert np.allclose(dlg._data[0]["amp"], y)
    assert dlg._data[0]["baseline_ops"] == []

    # Apply: commits it
    _place_two_points()
    monkeypatch.setattr(dlg, "_ask_apply_or_cancel", lambda plot, pos: "apply")
    dlg._cell_clicked(0, _RightClick())
    assert not dlg._cells[0]["bl_picking"]
    assert dlg._data[0]["baseline_ops"]
    assert not np.allclose(dlg._data[0]["amp"], y)


def test_batch_baseline_coincident_points_keep_picking_open():
    """Dragging both points to the same x must not silently discard the pick —
    the user needs to be able to keep adjusting, not start over."""
    import pyqtgraph as pg
    from larmor.desktop.batchfit_dialog import BatchFitDialog

    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 14.0},
                 "shift_fwhm_ppm": {"value": 5.0}, "amplitude": {"value": 80.0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    dlg = BatchFitDialog(None, [], model)
    x = np.linspace(-20, 60, 100)
    dlg._data = [{"ppm": x, "amp": x * 0.0, "amp0": x * 0.0, "nucleus": "11B",
                 "larmor": 160.0, "spin": 0.0, "sample": "s0", "path": "s0.csv",
                 "proc": "", "snr": 50, "baseline_ops": []}]
    plot = pg.PlotWidget()
    dlg._cells = [{"plot": plot, "exp": plot.plot([], []), "model": None,
                  "rmsd": None, "comp": [], "title": None,
                  "bl_picking": True, "bl_markers": [], "bl_line": None}]
    for pos in [(5.0, 1.0), (5.0, 2.0)]:      # same x -> a degenerate line
        m = pg.TargetItem(pos=pos, movable=True)
        plot.addItem(m)
        dlg._cells[0]["bl_markers"].append(m)

    dlg._apply_bg_pick(0)
    assert dlg._cells[0]["bl_picking"]        # still open — didn't discard
    assert len(dlg._cells[0]["bl_markers"]) == 2
    assert dlg._data[0]["baseline_ops"] == []


def _expno(tmp_path, procs):
    """Build an EXPNO with the given ``{proc: [fit filenames]}`` pdata layout."""
    expno = tmp_path / "10"
    for proc, fits in procs.items():
        d = expno / "pdata" / proc
        d.mkdir(parents=True)
        (d / "1r").write_bytes(b"\0")
        for fn in fits:
            (d / fn).write_text("{}")
    return expno


def test_procs_toggle_adds_layer_only_when_multiple(qapp, tmp_path):
    from larmor.desktop import explorer
    expno = _expno(tmp_path, {"1": [], "15": ["x.recipe.json"]})
    panel = explorer.ExplorerPanel()
    exp = QTreeWidgetItem(["10"])
    exp.setData(0, explorer._ROLE_KIND, "exp")
    exp.setData(0, explorer._ROLE_PATH, str(expno))
    panel.tree.addTopLevelItem(exp)
    panel.chkProcs.setChecked(False)
    assert exp.childCount() == 0
    panel.chkProcs.setChecked(True)
    assert exp.childCount() == 1                        # >1 proc → expandable


def test_single_proc_no_fits_not_expandable(qapp, tmp_path):
    from larmor.desktop import explorer
    expno = _expno(tmp_path, {"1": []})                 # one proc, no fits
    panel = explorer.ExplorerPanel()
    exp = QTreeWidgetItem(["10"])
    exp.setData(0, explorer._ROLE_KIND, "exp")
    exp.setData(0, explorer._ROLE_PATH, str(expno))
    panel.tree.addTopLevelItem(exp)
    panel.chkProcs.setChecked(False); panel.chkProcs.setChecked(True)
    assert exp.childCount() == 0                        # open directly, nothing to show


def test_single_proc_with_fits_shows_fits_on_experiment(qapp, tmp_path):
    # regression: a single-proc experiment that holds a fit must still show it
    # (no redundant proc layer) — the .fxml directly under the experiment
    from larmor.desktop import explorer
    expno = _expno(tmp_path, {"1": ["P1_31P.fxml"]})
    panel = explorer.ExplorerPanel()
    exp = QTreeWidgetItem(["3102"])
    exp.setData(0, explorer._ROLE_KIND, "exp")
    exp.setData(0, explorer._ROLE_PATH, str(expno))
    panel._reset_exp_children(exp)
    assert exp.childCount() == 1
    assert exp.child(0).data(0, explorer._ROLE_KIND) == "ph_expfit"
    exp.takeChildren()
    panel._add_fit_items(exp, explorer._procs_of(str(expno))[0])
    names = [exp.child(i).text(0) for i in range(exp.childCount())]
    assert any("P1_31P.fxml" in t for t in names)
    assert exp.child(0).data(0, explorer._ROLE_OPEN).endswith(".fxml")


def test_proc_without_fits_not_expandable(qapp, tmp_path):
    from larmor.desktop import explorer
    expno = _expno(tmp_path, {"1": [], "2": ["f.recipe.json"]})
    panel = explorer.ExplorerPanel()
    exp = QTreeWidgetItem(["10"])
    exp.setData(0, explorer._ROLE_PATH, str(expno))
    exp.setData(0, explorer._ROLE_KIND, "exp")
    panel._populate_procs(exp)
    kids = [exp.child(i) for i in range(exp.childCount())]
    proc1 = next(c for c in kids if c.text(0).strip().endswith("proc 1"))
    proc2 = next(c for c in kids if "proc 2" in c.text(0))
    assert proc1.childCount() == 0                      # no fits → not expandable
    assert proc2.childCount() == 1                      # has a fit → expandable
