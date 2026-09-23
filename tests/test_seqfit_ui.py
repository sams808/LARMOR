"""Sequential-fit dialog: navigation, per-spectrum table, auto sweep worker."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _series(tmp_path):
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import engine
    x = np.linspace(-20, 60, 500)
    paths = []
    for k, pos in enumerate([13.0, 15.0, 17.0]):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                    sites=[SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(pos),
                        "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100),
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


def test_seqfit_dialog_builds_and_navigates(qapp, tmp_path):
    from larmor.desktop.seqfit_dialog import SeqFitDialog
    paths, model = _series(tmp_path)
    dlg = SeqFitDialog(None, paths, model)
    assert len(dlg._data) == 3
    assert [dlg.cbPasses.itemText(i) for i in range(dlg.cbPasses.count())] == \
        ["1", "2", "4", "8", "16"]
    assert dlg._cur == 0
    dlg._next()
    assert dlg._cur == 1
    dlg._prev()
    assert dlg._cur == 0


def test_seqfit_edit_does_not_rescale_the_plot(qapp, tmp_path):
    # a manual zoom must survive edits/fits — only navigation rescales the x-axis
    from larmor.desktop.seqfit_dialog import SeqFitDialog
    paths, model = _series(tmp_path)
    dlg = SeqFitDialog(None, paths, model)
    vb = dlg.plot.getPlotItem().getViewBox()
    vb.setXRange(5.0, 25.0, padding=0)                 # user zooms in
    lo, hi = vb.viewRange()[0]
    dlg._resim_current()                                # a re-sim (edit/fit)
    assert vb.viewRange()[0] == pytest.approx([lo, hi], abs=0.01)   # unchanged
    dlg._show_current(rescale=False)                    # after a sweep — keep zoom
    assert vb.viewRange()[0] == pytest.approx([lo, hi], abs=0.01)


def test_seqfit_dialog_fit_current_and_auto(qapp, tmp_path):
    from larmor.desktop.seqfit_dialog import SeqFitDialog, _SeqWorker
    paths, model = _series(tmp_path)
    dlg = SeqFitDialog(None, paths, model)
    dlg._fit_current()
    assert np.isfinite(dlg._live_rmsd[0])

    w = _SeqWorker(dlg._entries(), 2, "first", dlg._propagate(), 0, None)
    w.done.connect(dlg._auto_done)
    w.run()
    assert dlg._result is not None and dlg._result.passes == 2
    pos = [__import__("larmor.recipe", fromlist=["Recipe"]).Recipe
           .from_dict(d).sites[0].params["isotropic_chemical_shift_ppm"].value
           for d in dlg._recipes]
    assert pos[0] == pytest.approx(13.0, abs=0.4)
    assert pos[-1] == pytest.approx(17.0, abs=0.4)


def test_seqfit_dialog_publication_bundle_after_manual_fits(qapp, tmp_path, monkeypatch):
    """Publication bundle… works after a manual Fit current alone (unfitted
    members get note 'not fitted'), recipes carry their source path, and a
    second bundle after the auto sweep has every member fitted."""
    import csv
    from PySide6.QtWidgets import QFileDialog
    from larmor.desktop.seqfit_dialog import SeqFitDialog, _SeqWorker
    from larmor.recipe import Recipe

    paths, model = _series(tmp_path)
    dlg = SeqFitDialog(None, paths, model)
    assert dlg.btnBundle.isEnabled()
    assert [dlg._seed_recipe(d)["source_path"] for d in dlg._data] == paths
    assert [d["source_path"] for d in dlg._recipes] == paths

    dlg._fit_current()                                  # spectrum 0 only
    folder = tmp_path / "b1"
    folder.mkdir()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *a, **k: str(folder)))
    dlg._export_bundle()
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
        assert (folder / r["curves_file"]).is_file()
    fitted = Recipe.from_dict(dlg._recipes[0])
    assert fitted.fit_rmsd is not None
    assert rows[0]["note"] == ""
    assert float(rows[0]["rmsd"]) == pytest.approx(fitted.fit_rmsd, rel=1e-6)
    assert [r["note"] for r in rows[1:]] == ["not fitted", "not fitted"]
    assert "manifest.csv" in dlg.status.text() and "3 spectra" in dlg.status.text()

    w = _SeqWorker(dlg._entries(), 2, "first", dlg._propagate(), 0, None)
    w.done.connect(dlg._auto_done)
    w.run()
    folder2 = tmp_path / "b2"
    folder2.mkdir()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *a, **k: str(folder2)))
    dlg._export_bundle()
    with open(folder2 / "manifest.csv", newline="", encoding="utf-8") as f:
        rows2 = list(csv.DictReader(f))
    assert [r["note"] for r in rows2] == ["", "", ""]
    assert all(float(r["rmsd"]) > 0 for r in rows2)
