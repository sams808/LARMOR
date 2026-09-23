"""Regressions for the small UI utilities added for the workflow tweaks:
the scroll-nudge opt-in flag and the Explorer's pdata-proc / fits browsing.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
# this module's figure-export test states LARMOR_NO_SESSION as its
# precondition but used to inherit it from whichever sibling module was
# collected first, so running this file alone read the real QSettings
os.environ.setdefault("LARMOR_NO_SESSION", "1")

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
                                   "ci68_lo", "ci68_hi", "model", "source_path"]
    amp = [ln for ln in lines[1:] if ln.split(",")[3] == "amplitude"]
    assert amp and all(ln.split(",")[6] == "covariance" for ln in amp)


def test_batch_dialog_covariance_errors_go_through_the_threaded_worker(qapp, tmp_path):
    """Regression: covariance is NOT a free snapshot (batch_fit's initial pass
    skips the errorbar-rescue retry for speed) -- Compute errors / Export CSV
    with Covariance selected must go through the same _ErrorWorker as Monte-
    Carlo/profile (so it actually refits with compute_errorbars=True), not the
    old synchronous "just read whatever's already there" shortcut, which left
    every exported row blank whenever the fast fit's covariance was singular."""
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor.desktop.batchfit_dialog import BatchFitDialog, _BatchWorker
    from larmor import batchfit, engine

    # a genuinely degenerate 2-site-per-spectrum setup (see test_batchfit.py's
    # _degenerate_batch_entries): amplitudes are perfectly correlated, so the
    # fast fit's covariance is reliably singular
    x = np.linspace(-30, 30, 400)
    paths = []
    for k in range(2):
        truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(0.0),
                "shift_fwhm_ppm": Param(10.0), "amplitude": Param(80.0),
                "gl": Param(1.0, vary=False)})])
        _, y, _ = engine.simulate(truth, exp_ppm=x)
        d = y + np.random.default_rng(k).normal(0, 0.3, x.size)
        p = tmp_path / f"s{k}.csv"
        p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                     "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, d)))
        paths.append(str(p))
    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "sites": [
        {"model": "gauss_lor", "label": "A", "params": {
            "isotropic_chemical_shift_ppm": {"value": 0.0, "vary": False},
            "shift_fwhm_ppm": {"value": 10.0, "vary": False},
            "amplitude": {"value": 40.0, "min": 0},
            "gl": {"value": 1.0, "vary": False}}},
        {"model": "gauss_lor", "label": "B", "params": {
            "isotropic_chemical_shift_ppm": {"value": 0.0, "vary": False},
            "shift_fwhm_ppm": {"value": 10.0, "vary": False},
            "amplitude": {"value": 40.0, "min": 0},
            "gl": {"value": 1.0, "vary": False}}}]}
    dlg = BatchFitDialog(None, paths, model)
    w = _BatchWorker(dlg._entries(), (), 0.1)
    w.done.connect(dlg._done)
    w.run()
    assert dlg._result is not None
    # (the exact "every stderr is None right after the fit" precondition is
    # proven deterministically at the pure-function level in test_batchfit.py's
    # test_error_analysis_covariance_refits_when_the_fast_fit_had_none; this
    # test's job is the DIALOG WIRING -- that Covariance goes through the
    # worker and ends up with real numbers regardless of the raw starting point)

    dlg.errCombo.setCurrentIndex(0)        # covariance
    dlg._compute_errors()                  # must route through the worker...
    assert dlg._err_worker is not None
    dlg._err_worker.run()                  # ...and run synchronously here
    dlg._err_done(dlg._result)
    rows = batchfit.error_table(dlg._result, method="covariance")
    assert any(r["stderr"] is not None for r in rows)
    assert "NO usable errors" not in dlg.status.text()


def test_batch_worker_request_stop_reaches_batch_fit(qapp):
    """_BatchWorker.request_stop must actually propagate to batch_fit's
    should_stop -- not just abort lmfit's iter_cb (which only affects the ONE
    spectrum being fit when Stop is pressed) -- so spectra later in the batch
    are also skipped, not fit to completion regardless."""
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import engine
    from larmor.desktop.batchfit_dialog import _BatchWorker

    x = np.linspace(-20, 60, 400)

    def spec(sh, amp, seed):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(15.0 + sh),
                "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        return m + np.random.default_rng(seed).normal(0, 1.0, x.size)

    def start(sample):
        return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=sample,
                      sites=[SiteModel(model="gauss_lor", label="A", params={
                          "isotropic_chemical_shift_ppm": Param(14.0, min=0, max=30),
                          "shift_fwhm_ppm": Param(5.0, min=0.1),
                          "amplitude": Param(80.0, min=0),
                          "gl": Param(1.0, vary=False)})])

    entries = [(start(f"g{k}"), x, spec(0.0, 100.0, k), (60.0, -20.0))
              for k in range(3)]
    w = _BatchWorker(entries, (), 0.1)

    orig_run = w.run

    def counting_run():
        # request_stop as soon as the batch is under way, before any spectrum
        # could plausibly finish fitting on its own
        w.request_stop("stop")
        orig_run()

    w.run = counting_run
    results = []
    w.done.connect(lambda res, mode: results.append((res, mode)))
    w.run()

    assert results
    res, mode = results[0]
    assert mode == "stop"
    assert len(res.recipes) == 3          # every entry still present, aligned
    # with stop requested before the loop even starts, nothing was fit
    for k in range(3):
        assert res.recipes[k].sites[0].params["amplitude"].value == \
            pytest.approx(80.0)


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


def test_batch_dialog_exclude_component_and_auto_save_recipes(qapp, tmp_path):
    """Right-click 'Exclude component' locks a site's amplitude to zero for
    ONE spectrum only; the exclusion is carried into _entries(), survives
    the fit, is left out of the exported table, and (with the "also save
    individual fits" checkbox, default on) the CSV export auto-saves a
    .recipe.json per spectrum next to it -- series_grid then finds the
    excluded spectrum's own recipe and reconstructs correctly."""
    from PySide6.QtWidgets import QLabel
    from larmor import batchfit, series_grid
    from larmor.desktop.batchfit_dialog import BatchFitDialog

    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0,
             "sites": [
                 {"model": "gauss_lor", "label": "A", "params": {
                     "isotropic_chemical_shift_ppm": {"value": 14.0},
                     "shift_fwhm_ppm": {"value": 5.0},
                     "amplitude": {"value": 80.0, "min": 0},
                     "gl": {"value": 1.0, "vary": False}}},
                 {"model": "gauss_lor", "label": "B", "params": {
                     "isotropic_chemical_shift_ppm": {"value": -5.0},
                     "shift_fwhm_ppm": {"value": 3.0},
                     "amplitude": {"value": 40.0, "min": 0},
                     "gl": {"value": 1.0, "vary": False}}}]}
    dlg = BatchFitDialog(None, [], model)
    assert dlg.chkAutoRecipes.isChecked()          # on by default

    x = np.linspace(-30, 40, 400)
    from larmor import engine
    from larmor.recipe import Recipe

    def make(sample, seed):
        rec = Recipe.from_dict({**model, "sample": sample})
        _, y, _ = engine.simulate(rec, exp_ppm=x)
        data = y + np.random.default_rng(seed).normal(0, 1.0, x.size)
        p = tmp_path / f"{sample}.csv"
        p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                     "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, data)))
        return p, data

    p0, y0 = make("g0", 0)
    p1, y1 = make("g1", 1)
    dlg._data = [
        {"ppm": x.copy(), "amp": y0, "amp0": y0.copy(), "nucleus": "11B",
         "larmor": 160.0, "spin": 0.0, "sample": "g0", "path": str(p0),
         "proc": "", "snr": 50, "baseline_ops": []},
        {"ppm": x.copy(), "amp": y1, "amp0": y1.copy(), "nucleus": "11B",
         "larmor": 160.0, "spin": 0.0, "sample": "g1", "path": str(p1),
         "proc": "", "snr": 50, "baseline_ops": []},
    ]
    dlg._cells = [{"title": QLabel()} for _ in range(2)]

    dlg._toggle_exclude(1, 1, True)          # exclude site B for g1 only
    assert dlg._excluded == {1: {1}}
    assert "(excluded: B)" in dlg._cells[1]["title"].text()
    dlg._toggle_exclude(1, 1, False)         # toggling off clears it
    assert dlg._excluded == {}
    dlg._toggle_exclude(1, 1, True)

    entries = dlg._entries()
    amp1 = entries[1][0].sites[1].params["amplitude"]
    assert amp1.value == 0.0 and amp1.vary is False
    assert entries[0][0].sites[1].params["amplitude"].vary is True  # g0 unaffected

    dlg._result = batchfit.batch_fit(entries)
    assert dlg._result.recipes[1].sites[1].params["amplitude"].value == 0.0

    csv_path = tmp_path / "batch_table.csv"
    from PySide6.QtWidgets import QFileDialog

    orig_get_save = QFileDialog.getSaveFileName
    QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (str(csv_path), ""))
    try:
        dlg._save_table()
    finally:
        QFileDialog.getSaveFileName = orig_get_save

    assert csv_path.exists()
    saved_recipes = list(tmp_path.glob("*.recipe.json"))
    assert len(saved_recipes) == 2                   # auto-saved alongside the CSV

    panels, warnings = series_grid.load_panels(str(csv_path))
    by_sample = {p.sample: p for p in panels}
    # tier-1 resolution matched the REAL saved recipes (not a CSV rebuild) --
    # full fidelity, so g1's recipe file still carries both sites (site B
    # zeroed, not deleted); what must actually change is what's DRAWN
    assert not by_sample["g0"].reconstructed and not by_sample["g1"].reconstructed
    assert by_sample["g0"].n_sites == 2 and by_sample["g1"].n_sites == 2

    from larmor import figures
    import matplotlib.pyplot as plt
    fig = figures.render({"kind": "batch_grid",
                          "panels": [{"recipe": by_sample["g1"].path}]})
    # site B (excluded, zero-amplitude) never drew a component or legend entry
    assert not fig.legends or "B" not in {t.get_text() for t in fig.legends[0].get_texts()}
    plt.close(fig)


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


def _czjzek_recipe_dict(sigma=1.1815):
    return {"nucleus": "27Al", "larmor_frequency_MHz": 156.28,
            "spin_rate_Hz": 35714.0, "sites": [
                {"model": "czjzek", "label": "Al", "params": {
                    "isotropic_chemical_shift_ppm": {"value": 64.5},
                    "sigma_Cq_MHz": {"value": sigma, "min": 0.1, "max": 4.0},
                    "shift_fwhm_ppm": {"value": 9.6},
                    "line_fwhm_ppm": {"value": 0.5},
                    "amplitude": {"value": 100.0}}}]}


def test_czjzek_display_convention_rescales_cell_and_header(qapp):
    """One fitted sigma, four literature conventions (sigma / 2sigma dmfit
    sCZ_CQ / 4sigma dmfit CQ box / sqrt5 sigma P_Q): the fit table's sigma
    column can display any of them, converts typed values AND bounds back to
    the stored sigma, and the column header names the active convention --
    the stored recipe value must never change with the display mode."""
    from larmor.desktop import table
    from larmor.desktop.table import LinesTable

    # the ambient mode is whatever the developer's own app settings hold (a
    # main-window test earlier in the suite applies the saved preference), so
    # pin the starting state instead of asserting the factory default
    prev = table.czjzek_display_mode()
    table.set_czjzek_display("sigma")
    t = LinesTable()
    try:
        rec = _czjzek_recipe_dict(sigma=1.0)
        for mode, k in [("sigma", 1.0), ("cq2", 2.0), ("dmfit", 4.0),
                        ("pq", 2.0 * 5.0 ** 0.5)]:
            table.set_czjzek_display(mode)
            t.rebuild(rec, set())
            col = 2 + t._used_keys.index("sigma_Cq_MHz")
            cell = t.table.cellWidget(0, col)
            shown = float(cell.edit.text())
            # rel=1e-4: the cell renders with "%.5g"
            assert shown == pytest.approx(1.0 * k, rel=1e-4), mode
            # the header names the convention
            head = t.table.horizontalHeaderItem(col).text()
            assert table.CZJZEK_DISPLAYS[mode][0] in head
            # stored value untouched by the display mode
            assert rec["sites"][0]["params"]["sigma_Cq_MHz"]["value"] == 1.0

        # typing a value in dmfit-CQ mode stores sigma = value / 4
        table.set_czjzek_display("dmfit")
        t.rebuild(rec, set())
        col = 2 + t._used_keys.index("sigma_Cq_MHz")
        cell = t.table.cellWidget(0, col)
        cell.edit.setText("4.7")
        cell._on_edit()
        assert rec["sites"][0]["params"]["sigma_Cq_MHz"]["value"] == \
            pytest.approx(4.7 / 4.0)
        # bounds typed in display units land in sigma units too
        cell.edit.setText("[0..8]")
        cell._on_edit()
        p = rec["sites"][0]["params"]["sigma_Cq_MHz"]
        assert p["min"] == pytest.approx(0.0)
        assert p["max"] == pytest.approx(2.0)
        # and the derived side-label shows sigma when sigma is hidden
        assert "σ" in cell.derived.text()
    finally:
        table.set_czjzek_display(prev)                   # never leak the mode


def test_pin_rename_persists_and_updates_label(qapp, tmp_path):
    """Right-click > Rename pin: the display name sticks to the pin, survives a
    restart (fresh panel), resets on empty input, and is dropped on unpin."""
    from PySide6.QtCore import QSettings

    from larmor.desktop import explorer
    s = QSettings("LARMOR", "app")
    old_pins = s.value("pinnedFolders", [])
    old_names = s.value("pinnedNames", "")
    try:
        s.setValue("pinnedFolders", [])
        s.setValue("pinnedNames", "")
        panel = explorer.ExplorerPanel()
        folder = tmp_path / "GlassSeries_2026"
        folder.mkdir()
        panel._pin(str(folder))
        assert "GlassSeries_2026" in panel.tree.topLevelItem(0).text(0)
        panel.set_pin_name(str(folder), "LAW glasses (MagLab)")
        assert panel.tree.topLevelItem(0).text(0) == "📌 LAW glasses (MagLab)"

        panel2 = explorer.ExplorerPanel()          # "restart"
        it2 = panel2.tree.topLevelItem(0)
        assert "LAW glasses (MagLab)" in it2.text(0)
        assert it2.toolTip(0) == str(folder)       # the real path stays visible
        panel2.set_pin_name(str(folder), "  ")     # empty -> folder's own name
        assert panel2._pin_label(str(folder)) == "GlassSeries_2026"
        panel2.set_pin_name(str(folder), "again")
        panel2._unpin(str(folder))                 # unpin forgets the name

        panel3 = explorer.ExplorerPanel()
        assert panel3._pin_names == {}
    finally:
        s.setValue("pinnedFolders", old_pins)
        s.setValue("pinnedNames", old_names)


def test_plot_menu_prunes_the_broken_average_submenu(qapp):
    """pyqtgraph's 'Plot Options > Average' popped up as a big empty white box
    (a bare QListWidget with nothing in it); it is removed on every plot that
    gets the shared menu."""
    import pyqtgraph as pg

    from larmor.desktop.plot_menu import attach_plot_menu
    pw = pg.PlotWidget()
    attach_plot_menu(pw, title="t")
    texts = [a.text() for a in pw.getPlotItem().ctrlMenu.actions()]
    assert texts and "Average" not in texts


def test_ppm_axes_never_grow_an_si_prefix(qapp):
    """A wide quadrupolar spectrum made pyqtgraph relabel the axis 'kppm' --
    kilo-ppm is not a unit a spectroscopist has ever used."""
    from larmor.desktop.qcpmg_dialog import _plot
    pw = _plot("x", ppm_axis=True)
    assert pw.getPlotItem().getAxis("bottom").autoSIPrefix is False


def test_amorphous_distribution_widths_get_table_columns(qapp):
    """dmfit's Amorphous fits a Gaussian Cq distribution (FWHM_CQ) -- the
    engine always fitted it, but the parameter-table column allowlist did not
    include it, so the value was invisible and uneditable from the table."""
    from larmor.desktop.table import LinesTable

    params = {k: {"value": v, "vary": True} for k, v in [
        ("isotropic_chemical_shift_ppm", 17.4), ("Cq_MHz", 2.6),
        ("eta", 0.2), ("Cq_fwhm_MHz", 0.3), ("eta_fwhm", 0.0),
        ("shift_fwhm_ppm", 5.0), ("line_fwhm_ppm", 0.5), ("gl", 0.0),
        ("amplitude", 62.0)]}
    rec = {"nucleus": "11B", "larmor_frequency_MHz": 160.46,
           "sites": [{"model": "amorphous", "label": "BO3", "params": params}]}
    t = LinesTable()
    t.rebuild(rec, set())
    for key in ("Cq_fwhm_MHz", "eta_fwhm", "line_fwhm_ppm"):
        assert key in t._used_keys, key
        col = 2 + t._used_keys.index(key)
        assert t.table.cellWidget(0, col) is not None, key
    # a recipe with no amorphous site does not grow the extra columns
    t.rebuild({"nucleus": "27Al", "larmor_frequency_MHz": 130.32, "sites": [
        {"model": "gauss_lor", "label": "", "params": {
            k: {"value": 1.0, "vary": True} for k in
            ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm",
             "amplitude", "gl")}}]}, set())
    assert "Cq_fwhm_MHz" not in t._used_keys


def test_suggest_save_dir_follows_the_current_dataset(tmp_path):
    """Save dialogs must propose the CURRENT dataset's own folder (where
    dmfit keeps fits and the Explorer lists them) -- not whichever folder the
    last file dialog visited."""
    from larmor.desktop.paths import suggest_save_dir

    sample = tmp_path / "pCABS2-4"
    proc = sample / "3616" / "pdata" / "1"
    proc.mkdir(parents=True)
    (sample / "3616" / "acqus").write_text("##$PULPROG= <zg>")
    (sample / "3616" / "fid").write_bytes(b"\0")
    (proc / "1r").write_bytes(b"\0")

    assert suggest_save_dir(proc / "1r") == str(proc)      # fit next to 1r
    assert suggest_save_dir(sample / "3616") == str(sample / "3616")
    plain = tmp_path / "exports"; plain.mkdir()
    (plain / "s.csv").write_text("x")
    assert suggest_save_dir(plain / "s.csv") == str(plain)
    # no source -> the caller's fallback
    assert suggest_save_dir(None, "FB") == "FB"
    assert suggest_save_dir(tmp_path / "gone" / "x.csv", "FB") == "FB"


def test_saves_allowed_anywhere_but_never_over_an_acquired_file(tmp_path):
    """Saving next to raw data is allowed (that is the dmfit convention);
    only REPLACING an acquired file itself is refused."""
    from larmor.desktop.paths import is_instrument_file

    proc = tmp_path / "3616" / "pdata" / "1"
    proc.mkdir(parents=True)
    (tmp_path / "3616" / "fid").write_bytes(b"\0")
    (proc / "1r").write_bytes(b"\0")

    assert not is_instrument_file(proc / "myfit.recipe.json")   # new file: fine
    assert not is_instrument_file(tmp_path / "3616" / "out.csv")
    assert is_instrument_file(proc / "1r")                      # the data itself
    assert is_instrument_file(tmp_path / "3616" / "fid")
    # the name alone is not enough -- only an EXISTING acquired file is guarded
    assert not is_instrument_file(tmp_path / "elsewhere" / "1r")


def test_app_save_dialogs_seed_from_the_open_dataset(qapp, tmp_path):
    from larmor.desktop.app import MainWindow

    w = MainWindow.__new__(MainWindow)      # no full app construction needed
    sample = tmp_path / "S1"
    (sample / "10" / "pdata" / "1").mkdir(parents=True)
    (sample / "10" / "acqus").write_text("x")
    w.source_path = str(sample / "10" / "pdata" / "1" / "1r")
    (sample / "10" / "pdata" / "1" / "1r").write_bytes(b"\0")
    assert w._suggest_dir() == str(sample / "10" / "pdata" / "1")


def _one_site_recipe(model_name):
    from larmor import models
    m = models.get(model_name)
    params = {p.name: {"value": p.default, "vary": p.vary,
                       "min": p.min, "max": p.max} for p in m.params}
    return {"nucleus": "27Al", "larmor_frequency_MHz": 130.32,
            "spin_rate_Hz": 20000.0,
            "sites": [{"model": model_name, "label": "", "params": params}]}


def test_every_registered_model_shows_every_parameter(qapp):
    """The Amorphous ΔCq FWHM was fitted but invisible because the table's
    column allowlist did not know the key; 8 of 15 models had the same class
    of hole (Voigt widths, J coupling, sideband ratio, ext-Czjzek eps, the
    two quad+CSA etas, ...). Every parameter of every registered model must
    get a column and a live cell -- including models added in the future,
    via the automatic-header fallback."""
    from larmor.models.base import REGISTRY

    from larmor.desktop.table import LinesTable
    t = LinesTable()
    for name, m in REGISTRY.items():
        t.rebuild(_one_site_recipe(name), set())
        for p in m.params:
            assert p.name in t._used_keys, f"{name}.{p.name} has no column"
            col = 2 + t._used_keys.index(p.name)
            assert t.table.cellWidget(0, col) is not None, f"{name}.{p.name}"
            head = t.table.horizontalHeaderItem(col).text()
            assert head.strip(), f"{name}.{p.name} has a blank header"


def test_mixed_model_recipe_gets_one_cell_per_own_parameter(qapp):
    """Sites of DIFFERENT models in one recipe: the table is the union of
    their parameters, each row carrying cells only for its own model's keys
    (blank holes elsewhere, never a wrong-model cell)."""
    from larmor.models.base import REGISTRY

    from larmor.desktop.table import LinesTable
    names = ["gauss_lor", "voigt", "jmultiplet", "czjzek", "quad_csa",
             "amorphous", "spectrum", "czjzek_d", "czjzek_corr", "exchange2"]
    sites = [_one_site_recipe(n)["sites"][0] for n in names]
    rec = {"nucleus": "27Al", "larmor_frequency_MHz": 130.32,
           "spin_rate_Hz": 20000.0, "sites": sites}
    t = LinesTable()
    t.rebuild(rec, set())
    for row, name in enumerate(names):
        keys = {p.name for p in REGISTRY[name].params}
        for c, key in enumerate(t._used_keys, start=2):
            cell = t.table.cellWidget(row, c)
            if key in keys:
                assert cell is not None, f"row {name} missing cell for {key}"
            else:
                assert cell is None, f"row {name} has a spurious {key} cell"


def test_dataset_info_text_reads_like_topspin_title_page(qapp):
    """Right-click > Dataset info on an EXPNO: the acquisition summary plus
    the FULL title text (the summary alone shows only its first line)."""
    from pathlib import Path

    from larmor.desktop import explorer
    expno = Path(__file__).resolve().parents[1] / "examples" / "pCABS2-4" / "3616"
    if not expno.exists():
        pytest.skip("example data not present")
    text = explorer.ExplorerPanel.dataset_info_text(str(expno))
    assert "27Al" in text and "zg" in text
    assert "MASR" in text                      # spinning rate is on the page
    assert "MAS=26 kHz" in text                # a later line of the full title


def test_datasets_panel_rows_show_detail_and_offer_color(qapp):
    from PySide6.QtWidgets import QLabel

    from larmor.desktop.datasets import DatasetsPanel

    p = DatasetsPanel()
    got = {}
    p.color_changed.connect(lambda i, c: got.update(i=i, c=c))
    p.rebuild("pCABS2-4", [
        {"label": "35Cl_2025-12 · 1", "color": "#e8832a", "visible": True,
         "source": "C:/x/1/pdata/1/1r", "nucleus": "35Cl",
         "larmor_MHz": 78.354, "npts": 4096, "title": "conditions OK"},
        {"label": "bare", "color": "#1f77b4", "visible": True, "source": ""},
    ], active_detail="35Cl · 78.4 MHz · 4096 pts")
    joined = " | ".join(w.text() for w in p._host.findChildren(QLabel))
    assert "78.4 MHz" in joined                 # active detail line
    assert "35Cl" in joined and "4096 pts" in joined
    assert "bare" in joined                     # detail-less overlay still fine
    # the swatch is a button that emits color_changed (picker itself is modal,
    # so emit directly to prove the wiring)
    p.color_changed.emit(1, "#123456")
    assert got == {"i": 1, "c": "#123456"}


def test_datasets_panel_empty_state_hints(qapp):
    from PySide6.QtWidgets import QLabel

    from larmor.desktop.datasets import DatasetsPanel
    p = DatasetsPanel()
    p.rebuild("something", [])
    texts = " ".join(w.text() for w in p._host.findChildren(QLabel))
    assert "no comparison spectra" in texts


def test_figure_exports_remember_their_folder(qapp, tmp_path, monkeypatch):
    """Exporting a second figure must start where the first one landed, not
    back in the LARMOR root. Gated on LARMOR_NO_SESSION so the suite itself
    never touches the real setting."""
    from PySide6.QtCore import QSettings

    from larmor.desktop import paths

    # under LARMOR_NO_SESSION (the test default) nothing is read or written
    assert paths.remembered_dir(paths.FIGURE_DIR_KEY, "FB") == "FB"
    paths.remember_dir(paths.FIGURE_DIR_KEY, tmp_path / "x.png")

    s = QSettings("LARMOR", "app")
    old = s.value(paths.FIGURE_DIR_KEY)
    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    try:
        fig_dir = tmp_path / "figs"; fig_dir.mkdir()
        paths.remember_dir(paths.FIGURE_DIR_KEY, fig_dir / "fig1.png")
        assert paths.remembered_dir(paths.FIGURE_DIR_KEY) == str(fig_dir)
        # a remembered folder that no longer exists falls back cleanly
        paths.remember_dir(paths.FIGURE_DIR_KEY, tmp_path / "gone" / "f.png")
        assert paths.remembered_dir(paths.FIGURE_DIR_KEY, "FB") == "FB"
    finally:
        if old is None:
            s.remove(paths.FIGURE_DIR_KEY)
        else:
            s.setValue(paths.FIGURE_DIR_KEY, old)


def test_computing_params_controls_are_all_wired(qapp, monkeypatch):
    """D1: the Computing-parameters dialog had two decorative controls -- a
    'Cq max (MHz)' whose setting nothing read (the 1D ceiling is the
    automatic ladder) and an 'η steps' that never reached the Czjzek
    render. The knob is gone (replaced by an explanatory label + a live
    cache readout) and η steps now genuinely changes the built kernel."""
    from larmor import engine
    from larmor.desktop.dialogs import ComputingParamsDialog
    from larmor.models.base import SimContext
    from larmor.models.quadrupolar import _render_czjzek

    assert "cq_max_MHz" not in engine.KERNEL_SETTINGS

    dlg = ComputingParamsDialog(None)
    assert not hasattr(dlg, "cqmax")
    old = dict(engine.KERNEL_SETTINGS)
    try:
        dlg.neta.setValue(5)
        dlg._accept()
        assert engine.KERNEL_SETTINGS["n_eta"] == 5

        # the wired path: the kernel the render builds carries the new n_eta
        ctx = SimContext(nucleus="27Al", larmor_MHz=130.32,
                         spin_rate_Hz=12500.0,
                         x_ppm=np.linspace(-120, 160, 512))
        engine.clear_kernel_cache()
        v = {"isotropic_chemical_shift_ppm": 60.0, "sigma_Cq_MHz": 1.5,
             "shift_fwhm_ppm": 5.0, "line_fwhm_ppm": 0.0, "amplitude": 1.0}
        _render_czjzek(v, ctx)
        (key, kernel), = engine._KERNEL_CACHE.items()
        assert kernel.eta_grid.size == 5
    finally:
        engine.KERNEL_SETTINGS.update(old)
        engine.clear_kernel_cache()


# ---------------------------------------------------------------------------
# F5: the batch results table under the grid, index-linked both ways to the
# spectrum cells. Synthetic '# nucleus = 11B' CSVs only (no real data), the
# dialog built alone and the fit run synchronously through its own worker.

def _batch_dialog(tmp_path, n, noise=None):
    """A BatchFitDialog over n synthetic 1-site 11B spectra (amplitude and
    shift vary along the series; ``noise[k]`` is spectrum k's noise sigma,
    1.5 for all by default) with a matching 1-site model loaded."""
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import engine
    from larmor.desktop.batchfit_dialog import BatchFitDialog

    x = np.linspace(-20, 60, 600)
    paths = []
    for k in range(n):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                    sites=[SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(15.0 + 0.1 * k),
                        "shift_fwhm_ppm": Param(6.0),
                        "amplitude": Param(100.0 - 5.0 * k),
                        "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        sigma = 1.5 if noise is None else noise[k]
        d = m + np.random.default_rng(k).normal(0, sigma, x.size)
        p = tmp_path / f"batch{k:02d}.csv"
        p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                     "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, d)),
                     encoding="utf-8")
        paths.append(str(p))
    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 15.0, "min": 0, "max": 30},
                 "shift_fwhm_ppm": {"value": 6.0, "min": 0.1},
                 "amplitude": {"value": 80.0, "min": 0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    return BatchFitDialog(None, paths, model)


def _fit_batch(dlg):
    """Run the batch fit synchronously through the dialog's own worker."""
    from larmor.desktop.batchfit_dialog import _BatchWorker

    w = _BatchWorker(dlg._entries(), (), 0.1)
    w.done.connect(dlg._done)
    w.run()
    assert dlg._result is not None
    return dlg._result


class _Click:
    """Stub of pyqtgraph's MouseClickEvent with only button() / scenePos() /
    accept(), like the _RightClick stubs above: the not-picking branch of
    _cell_clicked must never need more than the button."""

    def __init__(self, button):
        self._button = button

    def button(self):
        return self._button

    def scenePos(self):
        from PySide6.QtCore import QPointF
        return QPointF(0.0, 0.0)

    def accept(self):
        pass


def _headers(t):
    return [t.horizontalHeaderItem(c).text() for c in range(t.columnCount())]


def _selected_rows(t):
    return [i.row() for i in t.selectionModel().selectedRows()]


def test_batch_table_exists_before_fit_and_fills_after(qapp, tmp_path):
    """The table is there from construction (fixed columns, RMSD '—') so the
    link works while setting up; _done refills it with one column per site
    x parameter from the SAME rows Export CSV… writes, identity intact."""
    from PySide6.QtCore import Qt

    dlg = _batch_dialog(tmp_path, 2)
    t = dlg.table
    assert not t.isHidden()
    assert t.rowCount() == 2
    assert _headers(t) == ["#", "sample", "S/N", "RMSD"]     # CSVs: no proc column
    for k in range(2):
        assert t.item(k, 0).data(Qt.UserRole) == k
        assert t.item(k, 0).text() == str(k + 1)
        assert t.item(k, 1).text() == dlg._data[k]["sample"]
        assert t.item(k, 3).text() == "—"
    dlg._cell_clicked(1, _Click(Qt.LeftButton))
    assert dlg._hl == 1 and _selected_rows(t) == [1]

    res = _fit_batch(dlg)
    headers = _headers(t)
    assert headers[:4] == ["#", "sample", "S/N", "RMSD"]
    assert any(h.startswith("s0 A") for h in headers[4:])
    assert any(h.endswith("population %") for h in headers[4:])
    amp_col = headers.index("s0 A\namplitude")
    assert t.rowCount() == 2
    for k in range(2):
        r = dlg._row_of(k)
        assert t.item(r, 0).data(Qt.UserRole) == k
        shown = float(t.item(r, amp_col).text().split(" ±")[0])
        assert shown == pytest.approx(
            res.recipes[k].sites[0].params["amplitude"].value, rel=1e-3)
        assert f"{res.rmsd[k]:.4f}" in t.item(r, 3).text()
    # the spotlight survived the refill
    assert dlg._hl == 1 and _selected_rows(t) == [dlg._row_of(1)]


def test_batch_table_row_selection_spotlights_cell(qapp, tmp_path):
    from PySide6.QtCore import Qt
    from larmor.desktop import theme

    dlg = _batch_dialog(tmp_path, 2)
    res = _fit_batch(dlg)
    c0, c1 = dlg._cells
    dlg.table.selectRow(dlg._row_of(1))
    assert dlg._hl == 1
    assert c1["plot"].getViewBox().border.style() != Qt.NoPen
    assert c0["plot"].getViewBox().border.style() == Qt.NoPen
    assert c1["exp"].opts["pen"].widthF() == pytest.approx(2.0)
    assert c0["exp"].opts["pen"].widthF() == pytest.approx(1.0)
    assert c1["model"].opts["pen"].widthF() == pytest.approx(2.2)
    assert c0["model"].opts["pen"].widthF() == pytest.approx(1.4)
    assert c0["exp"].opacity() == pytest.approx(0.35)
    assert c1["exp"].opacity() == pytest.approx(1.0)
    accent = theme.active().accent
    assert accent in c1["title"].styleSheet()
    assert accent not in c0["title"].styleSheet()
    status = dlg.status.text()
    assert dlg._data[1]["sample"] in status and "RMSD" in status
    assert "Esc clears" in status

    dlg.table.clearSelection()
    assert dlg._hl is None
    for c in (c0, c1):
        assert c["plot"].getViewBox().border.style() == Qt.NoPen
        assert c["exp"].opacity() == pytest.approx(1.0)
        assert c["exp"].opts["pen"].widthF() == pytest.approx(1.0)
        assert accent not in c["title"].styleSheet()
    assert dlg.status.text() == res.summary


def test_batch_cell_click_selects_and_scrolls_to_its_row(qapp, tmp_path):
    from PySide6.QtCore import Qt

    dlg = _batch_dialog(tmp_path, 2)
    _fit_batch(dlg)
    t = dlg.table
    dlg._cell_clicked(1, _Click(Qt.LeftButton))
    assert _selected_rows(t) == [1]
    assert t.currentRow() == 1
    assert dlg._hl == 1
    # a right click outside picking mode keeps falling through to the menu
    dlg._cell_clicked(0, _Click(Qt.RightButton))
    assert _selected_rows(t) == [1] and dlg._hl == 1
    dlg._cell_clicked(0, _Click(Qt.LeftButton))
    assert _selected_rows(t) == [0] and dlg._hl == 0
    assert dlg._cells[1]["plot"].getViewBox().border.style() == Qt.NoPen
    assert dlg._cells[0]["plot"].getViewBox().border.style() != Qt.NoPen


def test_batch_highlight_brings_the_cells_tab_forward(qapp, tmp_path):
    from PySide6.QtCore import Qt

    dlg = _batch_dialog(tmp_path, 10)             # two pages; no fit needed
    assert dlg.tabs.count() == 2 and len(dlg._cells) == 10
    dlg._highlight_cell(9)
    assert dlg.tabs.currentIndex() == 1
    assert dlg._cells[9]["plot"].getViewBox().border.style() != Qt.NoPen
    assert dlg._cells[0]["exp"].opacity() == pytest.approx(0.35)
    dlg._highlight_cell(0)
    assert dlg.tabs.currentIndex() == 0
    dlg._highlight_cell(None)
    for c in dlg._cells:
        assert c["plot"].getViewBox().border.style() == Qt.NoPen
        assert c["exp"].opacity() == pytest.approx(1.0)


def test_batch_table_identity_survives_sorting(qapp, tmp_path):
    """Sort by RMSD descending puts the outlier on top; selecting that row
    must spotlight the outlier's SPECTRUM (identity = Qt.UserRole k), not
    whatever spectrum used to sit at row 0."""
    from PySide6.QtCore import Qt
    from larmor.desktop.batchfit_dialog import _NumItem

    dlg = _batch_dialog(tmp_path, 3, noise=(1.0, 6.0, 1.5))   # the middle one is bad
    res = _fit_batch(dlg)
    t = dlg.table
    worst = int(np.argmax(res.rmsd))
    assert worst == 1
    assert 1 in dlg._flag_reasons and 0 not in dlg._flag_reasons
    rmsd_col = _headers(t).index("RMSD")
    for k in range(3):                       # the flag gate is the cell label's
        txt = t.item(dlg._row_of(k), rmsd_col).text()
        assert txt.startswith("⚠ ") == (k in dlg._flag_reasons)
        assert f"{res.rmsd[k]:.4f}" in txt

    t.sortItems(rmsd_col, Qt.DescendingOrder)
    assert t.item(0, 0).text() == str(worst + 1)
    t.selectRow(0)
    assert dlg._hl == worst
    assert _selected_rows(t) == [0]
    # the sort persists across a refill and the spotlight follows its spectrum
    dlg._fill_table(res)
    assert t.item(0, 0).text() == str(worst + 1)
    assert dlg._hl == worst and _selected_rows(t) == [0]

    # numeric, not lexicographic; blanks sink to the bottom
    a = _NumItem("9.5"); a.setData(Qt.UserRole + 1, 9.5)
    b = _NumItem("10.2"); b.setData(Qt.UserRole + 1, 10.2)
    assert a < b and not (b < a)
    blank = _NumItem("")
    assert a < blank and not (blank < a)


def test_batch_escape_clears_selection_before_closing_and_arrows_step(qapp, tmp_path):
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QDialogButtonBox, QPushButton

    dlg = _batch_dialog(tmp_path, 3)
    dlg.show()                                    # offscreen; key stepping needs a laid-out view
    try:
        # Enter can never silently press a button: nothing is auto-default,
        # so Qt promoted no button to default on show
        assert all(not b.autoDefault()
                   for b in dlg.findChild(QDialogButtonBox).buttons())
        assert not any(b.isDefault() for b in dlg.findChildren(QPushButton))

        dlg._select_spectrum(0)
        assert dlg._hl == 0
        assert dlg.focusWidget() is dlg.table     # so the arrows work at once
        QTest.keyClick(dlg.table, Qt.Key_Down)
        assert dlg._hl == 1
        QTest.keyClick(dlg.table, Qt.Key_Down)
        assert dlg._hl == 2
        QTest.keyClick(dlg.table, Qt.Key_Up)
        assert dlg._hl == 1

        rejected = []
        dlg.rejected.connect(lambda: rejected.append(1))
        QTest.keyClick(dlg.table, Qt.Key_Escape)  # ignored by the view, reaches the dialog
        assert dlg._hl is None and rejected == []
        assert _selected_rows(dlg.table) == []
        QTest.keyClick(dlg, Qt.Key_Escape)        # QDialog's default is untouched
        assert rejected == [1]
    finally:
        dlg.close()


def test_batch_highlight_survives_components_toggle_and_error_refill(qapp, tmp_path):
    dlg = _batch_dialog(tmp_path, 2)
    res = _fit_batch(dlg)
    t = dlg.table
    t.selectRow(0)
    assert dlg._hl == 0
    dlg.chkComp.setChecked(True)                  # rebuilds the component curves
    assert dlg._cells[0]["comp"] and dlg._cells[1]["comp"]
    assert all(it.opacity() == pytest.approx(0.35) for it in dlg._cells[1]["comp"])
    assert all(it.opacity() == pytest.approx(1.0) for it in dlg._cells[0]["comp"])

    # covariance errors through the dialog's own worker path (Compute errors),
    # the thread joined here instead of run through an event loop
    dlg.errCombo.setCurrentIndex(0)
    dlg._compute_errors()
    assert dlg._err_worker is not None and dlg._err_worker.wait(120_000)
    dlg._err_done(dlg._result)
    assert "covariance" in res.error_detail
    assert t.rowCount() == 2
    assert _selected_rows(t) == [0] and dlg._hl == 0
    amp_col = _headers(t).index("s0 A\namplitude")
    cells = [t.item(r, amp_col) for r in range(2)]
    assert any(" ± " in c.text() for c in cells)
    assert all("covariance" in c.toolTip() for c in cells)


def test_batch_table_checkbox_hides_and_shows(qapp, tmp_path):
    from PySide6.QtCore import Qt

    dlg = _batch_dialog(tmp_path, 2)
    assert dlg.chkTable.isChecked() and not dlg.table.isHidden()
    dlg.chkTable.setChecked(False)
    assert dlg.table.isHidden()                   # never isVisible(): the dialog is not shown
    dlg.table.selectRow(1)                        # the link still works while hidden
    assert dlg._hl == 1
    dlg._cell_clicked(0, _Click(Qt.LeftButton))   # ...and never focuses a hidden table
    assert dlg._hl == 0 and dlg.focusWidget() is not dlg.table
    dlg.chkTable.setChecked(True)
    assert not dlg.table.isHidden()


# ---------------------------------------------------------------------------
# N1: Publication bundle… in the batch dialog (larmor.io.bundle behind it).

def _read_bundle_curves(path):
    """A bundle _curves.csv as a structured array, its '# key=value' header
    dropped first (genfromtxt would take the first '#' line as the header)."""
    from pathlib import Path
    rows = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines()
            if ln and not ln.startswith("#")]
    return np.genfromtxt(rows, delimiter=",", names=True)


def test_batch_dialog_publication_bundle_button_writes_the_folder(qapp, tmp_path, monkeypatch):
    import csv
    from PySide6.QtWidgets import QFileDialog, QMessageBox
    from larmor.desktop import batchfit_dialog

    dlg = _batch_dialog(tmp_path, 3)
    assert not dlg.btnBundle.isEnabled()                 # like Save table… before a fit
    _fit_batch(dlg)
    assert dlg.btnBundle.isEnabled()

    folder = tmp_path / "bundle"
    folder.mkdir()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *a, **k: str(folder)))
    dlg._export_bundle()
    names = {p.name for p in folder.iterdir()}
    assert {"manifest.csv", "README.txt", "batch_table.csv"} <= names
    assert len(list(folder.glob("*.recipe.json"))) == 3
    assert len(list(folder.glob("*_curves.csv"))) == 3
    assert not list(folder.glob("batch_fit_*.csv"))      # no error method computed
    with open(folder / "manifest.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    for k, d in enumerate(dlg._data):
        tab = _read_bundle_curves(folder / rows[k]["curves_file"])
        assert np.array_equal(tab["experiment"], d["amp"])       # exactly as fitted
        assert np.array_equal(tab["ppm"], d["ppm"])
        assert "experiment_raw" not in tab.dtype.names            # no baseline applied
        assert rows[k]["source_path"] == d["path"]
        assert rows[k]["source_kind"] == "csv" and len(rows[k]["source_sha256"]) == 64
    assert "manifest.csv" in dlg.status.text() and "3 spectra" in dlg.status.text()

    # an existing bundle: the replace question answered No writes nothing
    before = (folder / "manifest.csv").read_bytes()
    stamp = (folder / "manifest.csv").stat().st_mtime_ns
    monkeypatch.setattr(QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.No))
    dlg._export_bundle()
    assert (folder / "manifest.csv").stat().st_mtime_ns == stamp
    assert (folder / "manifest.csv").read_bytes() == before

    # the button follows Save table… through an error worker
    monkeypatch.setattr(batchfit_dialog._ErrorWorker, "start", lambda self: None)
    dlg._start_err_worker("covariance")
    assert not dlg.btnBundle.isEnabled() and not dlg.btnTable.isEnabled()
    dlg._post_err_enable()
    assert dlg.btnBundle.isEnabled() and dlg.btnTable.isEnabled()


def test_batch_bundle_experiment_is_the_baselined_array(qapp, tmp_path, monkeypatch):
    """A per-spectrum baseline: the bundle's 'experiment' is d['amp'] (what
    the fit saw), 'experiment_raw' is d['amp0'], and only for that spectrum;
    the manifest's processing column carries the recorded op."""
    import csv
    from PySide6.QtWidgets import QFileDialog
    from larmor import batchfit
    from larmor.desktop.batchfit_dialog import estimate_baseline

    dlg = _batch_dialog(tmp_path, 2)
    d0 = dlg._data[0]
    d0["amp0"] = d0["amp0"] + 20.0                       # a raised baseline on spectrum 0
    base = estimate_baseline(d0["ppm"], d0["amp0"], "Flat (edge median)")
    d0["amp"] = d0["amp0"] - base
    d0["baseline_ops"] = [{"op": "flat_baseline"}]
    dlg._cells[0]["exp"].setData(d0["ppm"], d0["amp"])
    res = batchfit.batch_fit(dlg._entries())
    dlg._done(res)
    assert res.recipes[0].processing == [{"op": "flat_baseline"}]
    assert res.recipes[1].processing == []

    folder = tmp_path / "b"
    folder.mkdir()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *a, **k: str(folder)))
    dlg._export_bundle()
    with open(folder / "manifest.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    tab0 = _read_bundle_curves(folder / rows[0]["curves_file"])
    assert list(tab0.dtype.names)[:3] == ["ppm", "experiment", "experiment_raw"]
    assert np.array_equal(tab0["experiment"], d0["amp"])
    assert np.array_equal(tab0["experiment_raw"], d0["amp0"])
    assert not np.array_equal(tab0["experiment"], tab0["experiment_raw"])
    tab1 = _read_bundle_curves(folder / rows[1]["curves_file"])
    assert "experiment_raw" not in tab1.dtype.names
    assert np.array_equal(tab1["experiment"], dlg._data[1]["amp"])
    assert rows[0]["processing"].startswith('[{"op": "')
    assert rows[1]["processing"] == "[]"
