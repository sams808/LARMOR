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
    w = _BatchWorker(dlg._entries(), (), 0.1)
    w.done.connect(dlg._done)
    w.run()                                   # synchronous
    assert dlg._result is not None
    assert dlg.btnSave.isEnabled()
    assert dlg._cells[0]["model"].xData is not None


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
