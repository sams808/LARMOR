"""Project files: save / reopen the whole session through the versioned
.larproj.json bundle (larmor.project) -- 1D spectra + fits (+ overlays by
reference), 2D maps by reference with their recorded processing, kept figures
and batch-fit sessions. Synthetic data throughout except the one real-data
2D round trip, which skips cleanly when the dataset is absent."""
import copy
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")
import pytest

from conftest import BRUKER_2RR_MQMAS, require

pytest.importorskip("PySide6")


def test_project_roundtrip():
    from PySide6.QtWidgets import QApplication, QFileDialog
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    w._confirm_open_mode = lambda: "replace"       # fitted workspaces are open
    w._report_project_notes = lambda notes: None   # never a modal box offscreen
    w._add_recent = lambda p: None                 # keep the real Open recent clean

    def mkws(sample, pos):
        x = np.linspace(-40, 120, 300)
        y = np.exp(-0.5 * ((x - pos) / 8) ** 2)
        w.exp_ppm, w.exp_amp = x, y
        w.recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
                    "sample": sample, "sites": [
                        {"model": "gauss_lor", "label": "A", "params": {
                            "isotropic_chemical_shift_ppm": {"value": pos},
                            "shift_fwhm_ppm": {"value": 8},
                            "gl": {"value": 0.5},
                            "amplitude": {"value": 1.0}}}]}
        w.source_path = "src_" + sample
        w.hidden = set()
        w.view.set_experiment(x, y)
        w.lines_table.rebuild(w.recipe, w.hidden)
        w._ws_mode = "new"
        w._register_ws("1d")

    mkws("glassA", 60)
    mkws("glassB", 30)
    assert len(w.workspaces) == 2

    proj = tempfile.mktemp(suffix=".larproj.json")
    QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.save_project()
    assert os.path.exists(proj)

    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.open_project()
    assert len(w.workspaces) == 2
    samples = [ws["snap"]["recipe"]["sample"] for ws in w.workspaces]
    assert samples == ["glassA", "glassB"]
    pos = [ws["snap"]["recipe"]["sites"][0]["params"]
           ["isotropic_chemical_shift_ppm"]["value"] for ws in w.workspaces]
    assert pos == [60, 30]
    w.close()


def _write_csv_spectrum(path, nucleus, x, y):
    path.write_text("# nucleus = " + nucleus + "\n# larmor_MHz = 130.3\n"
                    + "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, y)))


def test_add_overlay_dialog_reads_a_recipe_json_source(tmp_path):
    """Regression: add_overlay_dialog used to unpack load_any()'s return
    tuple in the wrong order (recipe, ppm, amp) instead of the real
    (ppm, amp, recipe, ...), so recipe.get(...) always raised and EVERY
    overlay source _load_any actually supports (recipe.json, fxmla,
    csv/txt/dat) silently fell through to the Bruker-only fallback and
    failed there too -- caught while wiring overlay round-trip into
    save_project/open_project."""
    from PySide6.QtWidgets import QApplication, QFileDialog
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    x = np.linspace(-40, 120, 300)
    y = np.exp(-0.5 * ((x - 20) / 8) ** 2)
    csv_path = tmp_path / "overlay_source.csv"
    _write_csv_spectrum(csv_path, "27Al", x, y)

    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (str(csv_path), ""))
    n_before = len(w._overlays)
    w.add_overlay_dialog()
    assert len(w._overlays) == n_before + 1
    ov = w._overlays[-1]
    assert ov["source"] == str(csv_path)
    assert ov["label"] == "overlay_source"        # from the csv's own stem
    # abs=5e-4: the csv round-trips through "%.4f" text, not full precision
    assert np.allclose(ov["ppm"], x, atol=5e-4)
    assert np.allclose(ov["amp"], y, atol=5e-4)
    w.close()


def test_project_roundtrip_includes_overlays(tmp_path):
    """save_project captured overlays in the live snapshot but never wrote
    them into the saved file, so every overlay silently vanished on
    reopen -- fixed alongside the add_overlay_dialog bug above."""
    from PySide6.QtWidgets import QApplication, QFileDialog
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    w._confirm_open_mode = lambda: "replace"
    w._report_project_notes = lambda notes: None
    w._add_recent = lambda p: None
    x = np.linspace(-40, 120, 300)
    y = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    w.exp_ppm, w.exp_amp = x, y
    w.recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
               "sample": "glassA", "sites": []}
    w.source_path = str(tmp_path / "main.csv")
    _write_csv_spectrum(Path(w.source_path), "27Al", x, y)
    w.hidden = set()
    w.view.set_experiment(x, y)
    w.lines_table.rebuild(w.recipe, w.hidden)

    ov_path = tmp_path / "compare.csv"
    ov_y = np.exp(-0.5 * ((x - 30) / 8) ** 2)
    _write_csv_spectrum(ov_path, "27Al", x, ov_y)
    w._add_overlay("compare", x, ov_y, str(ov_path))
    w._overlays[-1]["visible"] = False

    w._ws_mode = "new"
    w._register_ws("1d")

    proj = tempfile.mktemp(suffix=".larproj.json")
    QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.save_project()

    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.open_project()

    assert len(w.workspaces) == 1
    overlays = w.workspaces[0]["snap"]["overlays"]
    assert len(overlays) == 1
    assert overlays[0]["label"] == "compare"
    assert overlays[0]["source"] == str(ov_path)
    assert overlays[0]["visible"] is False
    # peak position/height, not a full-array compare: a Gaussian's tail spans
    # many orders of magnitude and the csv's "%.4f" text rounds far-tail
    # values to exactly 0 -- checking the peak is what actually matters here
    # (that the right FILE came back), not bit-for-bit tail precision
    restored_amp = np.asarray(overlays[0]["amp"])
    restored_ppm = np.asarray(overlays[0]["ppm"])
    assert restored_amp.max() == pytest.approx(ov_y.max(), abs=1e-3)
    assert restored_ppm[restored_amp.argmax()] == pytest.approx(30.0, abs=0.5)
    w.close()


# =========================================================================
# G1: the bundle carries the whole session (2D maps by reference, figures,
# batch-fit sessions), versioned v2 with a v1 -> v2 migration.

@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    """A MainWindow whose modal project prompts are stubbed: Replace on the
    open-into-a-session question, the 'Project opened' notes recorded (one
    list per call) instead of shown, Open recent left alone."""
    monkeypatch.setenv("LARMOR_NO_KERNEL_WARM", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    w._note_calls = []
    w._confirm_open_mode = lambda: "replace"
    w._report_project_notes = lambda notes: w._note_calls.append(list(notes))
    w._add_recent = lambda p: None
    yield w
    w.close()


def _notes_text(w) -> str:
    return "\n".join(n for call in w._note_calls for n in call)


def _mkws(w, sample, pos):
    """Register a fitted 1D workspace (one gauss_lor line) on the window."""
    x = np.linspace(-40, 120, 300)
    y = np.exp(-0.5 * ((x - pos) / 8) ** 2)
    w.exp_ppm, w.exp_amp = x, y
    w.recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
                "sample": sample, "sites": [
                    {"model": "gauss_lor", "label": "A", "params": {
                        "isotropic_chemical_shift_ppm": {"value": pos},
                        "shift_fwhm_ppm": {"value": 8},
                        "gl": {"value": 0.5},
                        "amplitude": {"value": 1.0}}}]}
    w.source_path = "src_" + sample
    w.hidden = set()
    w.central_stack.setCurrentWidget(w.view)
    w.view.set_experiment(x, y)
    w.lines_table.rebuild(w.recipe, w.hidden)
    w._ws_mode = "new"
    w._register_ws("1d")


def _model():
    """A two-site 11B model (as tests/test_batchfit.py's _start)."""
    return {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
            "sites": [
                {"model": "gauss_lor", "label": "A", "params": {
                    "isotropic_chemical_shift_ppm": {"value": 14.0, "min": 0, "max": 30},
                    "shift_fwhm_ppm": {"value": 5.0, "min": 0.1},
                    "amplitude": {"value": 80.0, "min": 0},
                    "gl": {"value": 1.0, "vary": False}}},
                {"model": "gauss_lor", "label": "B", "params": {
                    "isotropic_chemical_shift_ppm": {"value": 1.0, "min": -10, "max": 15},
                    "shift_fwhm_ppm": {"value": 3.5, "min": 0.1},
                    "amplitude": {"value": 50.0, "min": 0},
                    "gl": {"value": 1.0, "vary": False}}}]}


def _csv_spectra(tmp_path, n):
    """n synthetic two-site 11B CSV spectra loader.load_any reads."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 60, 500)
    paths = []
    for k in range(n):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(15.0 + 0.2 * k),
                "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0 - 10.0 * k),
                "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="B", params={
                "isotropic_chemical_shift_ppm": Param(2.0),
                "shift_fwhm_ppm": Param(3.0), "amplitude": Param(40.0 + 10.0 * k),
                "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        y = m + np.random.default_rng(k).normal(0, 1.5, x.size)
        p = tmp_path / f"batch{k:02d}.csv"
        p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n"
                     + "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, y)),
                     encoding="utf-8")
        paths.append(str(p))
    return paths


def _v1_file(path, samples):
    """A v1 bundle exactly as the previous release wrote it: 1D workspaces
    with embedded arrays, no 'kind' on the entries (tests/test_batch.py's
    _make_larproj shape)."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 60, 200)
    workspaces = []
    for name, pos in samples:
        r = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=name,
                   sites=[SiteModel(model="gauss_lor", label="A", params={
                       "isotropic_chemical_shift_ppm": Param(pos),
                       "shift_fwhm_ppm": Param(6.0), "amplitude": Param(80.0),
                       "gl": Param(1.0, vary=False)})])
        _, y, _ = engine.simulate(r, exp_ppm=x)
        workspaces.append({"title": name, "source_path": f"src_{name}",
                           "recipe": r.to_dict(), "hidden": [],
                           "exp_ppm": x.tolist(), "exp_amp": y.tolist()})
    Path(path).write_text(json.dumps({"larmor_project_version": 1, "active": 0,
                                      "workspaces": workspaces}), encoding="utf-8")
    return str(path)


def _save(w, monkeypatch, path):
    from PySide6.QtWidgets import QFileDialog
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(path), "")))
    w.save_project()
    assert Path(path).exists()


def _open(w, monkeypatch, path):
    from PySide6.QtWidgets import QFileDialog
    monkeypatch.setattr(QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(path), "")))
    w.open_project()


def _dumps(recipe) -> str:
    """Canonical text of a Recipe, NaN-tolerant (json writes NaN as NaN)."""
    return json.dumps(recipe.to_dict(), sort_keys=True)


# ------------------------------------------------------------- schema (Qt-free)
def test_project_bundle_v1_migrates_to_v2_with_defaults():
    from larmor import project

    d = {"larmor_project_version": 1, "active": 0,
         "workspaces": [{"title": "a", "recipe": {"sites": []},
                         "exp_ppm": [1, 2], "exp_amp": [0, 1]}]}
    out, notes = project.migrate_bundle(d)
    assert out["workspaces"][0]["kind"] == "1d"
    assert out["larmor_project_version"] == 2
    assert any("project migrated from schema v1 to v2" in n for n in notes)
    assert "kind" not in d["workspaces"][0]          # the caller's dict is untouched

    v2 = {"larmor_project_version": 2, "active": None,
          "workspaces": [{"kind": "figure", "title": "f", "spec": {"kind": "1d"}}]}
    ref = copy.deepcopy(v2)
    out2, notes2 = project.migrate_bundle(v2)
    assert notes2 == [] and out2 == ref

    out3, notes3 = project.migrate_bundle({"larmor_project_version": 3, "workspaces": []})
    assert len(notes3) == 1 and "newer LARMOR" in notes3[0] and out3["workspaces"] == []

    out4, _ = project.migrate_bundle({"workspaces": [{"title": "x"}]})   # no version: v1
    assert out4["workspaces"][0]["kind"] == "1d" and out4["larmor_project_version"] == 2


def test_project_migration_chain_runs_every_hop_through_run_migrations(monkeypatch):
    from larmor import project, recipe

    monkeypatch.setattr(project, "PROJECT_BUNDLE_VERSION", 3)
    monkeypatch.setitem(project._MIGRATIONS, 2, lambda d: {**d, "hop2": True})
    out, notes = project.migrate_bundle(
        {"larmor_project_version": 1, "workspaces": [{"title": "a"}]})
    assert out["workspaces"][0]["kind"] == "1d" and out["hop2"] is True
    assert notes == ["project migrated from schema v1 to v2",
                     "project migrated from schema v2 to v3"]
    assert out["larmor_project_version"] == 3

    def bump(d):
        return {**d, "sample": d.get("sample", "") + " [migrated]"}
    assert recipe.run_migrations({"sample": "x"}, 0, {0: bump}, 1, "recipe") == \
        ({"sample": "x [migrated]"}, ["recipe migrated from schema v0 to v1"])
    # a missing hop stops the chain, as the recipe loader always did
    assert recipe.run_migrations({"a": 1}, 0, {}, 2) == ({"a": 1}, [])


def test_relocate_prefers_absolute_then_project_relative(tmp_path):
    from larmor import project

    proj = tmp_path / "proj"
    f = proj / "S1" / "35" / "pdata" / "1" / "2rr"
    f.parent.mkdir(parents=True)
    f.write_bytes(b"\x00" * 8)
    assert project.relocate(str(f), str(proj)) == str(f)
    assert project.relocate("D:/old/S1/35/pdata/1/2rr", str(proj),
                            rel="S1/35/pdata/1/2rr") == str(f)
    assert project.relocate("D:/old/S1/35/pdata/1/2rr", str(proj)) is None
    assert project.relocate("", str(proj), rel=None) is None
    assert project.relpath_or_none(str(f), str(proj)) == "S1/35/pdata/1/2rr"
    assert project.relpath_or_none("", str(proj)) is None
    if sys.platform == "win32":
        assert project.relpath_or_none("D:/x/y", "C:/z") is None   # other drive
    assert project.json_safe({"a": np.float64(1.5), "b": np.arange(2), "c": {1, 2},
                              "d": (1, 2), "e": np.bool_(True), "f": Path("x")}) == \
        {"a": 1.5, "b": [0, 1], "c": [1, 2], "d": [1, 2], "e": True, "f": "x"}


def test_relocate_batch_state_rewrites_moved_paths(tmp_path):
    """A project moved together with its spectra: the batch session's paths
    relocate through the project-relative copies, and every reference to a
    moved path -- per-spectrum settings and the restored recipes' source_path,
    which is what aligns results to spectra -- follows."""
    from larmor import project

    new_dir = tmp_path / "moved"
    (new_dir / "data").mkdir(parents=True)
    a = new_dir / "data" / "a.csv"
    a.write_text("x")
    old_a, old_b = "D:/old/data/a.csv", "D:/old/data/b.csv"
    state = {"paths": [old_a, old_b], "paths_rel": ["data/a.csv", "data/b.csv"],
             "per_spectrum": {old_a: {"excluded": [0]}, old_b: {"excluded": []}},
             "result": {"recipes": [{"source_path": old_a}, {"source_path": old_b}]}}
    st, missing = project.relocate_batch_state(state, new_dir)
    assert st["paths"] == [str(a), old_b] and missing == [old_b]
    assert st["per_spectrum"] == {str(a): {"excluded": [0]}, old_b: {"excluded": []}}
    assert [r["source_path"] for r in st["result"]["recipes"]] == [str(a), old_b]
    assert state["paths"] == [old_a, old_b]          # input untouched


def test_relocate_batch_state_rewrites_series_row_paths(tmp_path):
    """The Series table's rows pair with spectra by source_path
    (SeriesTable.aligned_to), so a moved project rewrites them with the same
    mapping as per_spectrum and the recipes; the input dict is untouched."""
    from larmor import project

    new_dir = tmp_path / "moved"
    (new_dir / "data").mkdir(parents=True)
    a = new_dir / "data" / "a.csv"
    a.write_text("x")
    old_a, old_b = "D:/old/data/a.csv", "D:/old/data/b.csv"
    rows = [{"source_path": old_a, "display_name": "a", "group": "a", "folder": "",
             "title": "", "values": {"CaO": 1.0}},
            {"source_path": old_b, "display_name": "b", "group": "a", "folder": "",
             "title": "", "values": {"CaO": 2.0}}]
    state = {"paths": [old_a, old_b], "paths_rel": ["data/a.csv", "data/b.csv"],
             "per_spectrum": {old_a: {}, old_b: {}},
             "series": {"version": 1, "rows": rows,
                        "columns": [{"key": "CaO", "label": "CaO", "err_key": None,
                                     "tag": "", "source": ""}]}}
    st, missing = project.relocate_batch_state(state, new_dir)
    assert missing == [old_b]
    assert [r["source_path"] for r in st["series"]["rows"]] == [str(a), old_b]
    assert st["series"]["rows"][1]["values"] == {"CaO": 2.0}
    assert st["series"]["columns"] == state["series"]["columns"]
    assert state["series"]["rows"][0]["source_path"] == old_a     # input untouched
    # nothing moved: the series dict passes through as is
    st2, _ = project.relocate_batch_state({"paths": [str(a)], "series": state["series"]},
                                          new_dir)
    assert st2["series"] is state["series"]


def test_build_bundle_remaps_active_and_counts_drops(tmp_path):
    from larmor import project

    x = np.linspace(0, 1, 5)
    ws_2d_gone = {"kind": "2d", "title": "gone", "has_fit": False,
                  "snap": {"kind": "2d", "source_path": str(tmp_path / "nope" / "2rr"),
                           "recipe": None, "hidden": set(), "view2d": {}}}
    ws_1d = {"kind": "1d", "title": "s", "has_fit": False,
             "snap": {"kind": "1d", "source_path": "p", "recipe": {"sites": []},
                      "hidden": set(), "exp_ppm": x, "exp_amp": x,
                      "overlays": [{"label": "o", "source": ""},
                                   {"label": "k", "source": "f.csv", "color": "#000",
                                    "visible": False}]}}
    ws_fig = {"kind": "figure", "title": "f", "has_fit": False,
              "snap": {"kind": "figure",
                       "spec": {"kind": "1d", "traces": [], "figsize": (3, 2)}}}
    bundle, dropped = project.build_bundle([ws_2d_gone, ws_1d, ws_fig], 1, tmp_path)
    assert dropped == 1 and bundle["active"] == 0     # the active index follows the drop
    assert [w["kind"] for w in bundle["workspaces"]] == ["1d", "figure"]
    assert bundle["workspaces"][0]["overlays"] == [
        {"label": "k", "color": "#000", "visible": False, "source": "f.csv"}]
    assert bundle["workspaces"][1]["spec"]["figsize"] == [3, 2]
    assert bundle["larmor_project_version"] == 2
    assert project.summary(bundle) == "1 spectrum, 1 figure"
    json.dumps(bundle)
    bundle2, _ = project.build_bundle([ws_2d_gone, ws_1d], 0, tmp_path)
    assert bundle2["active"] is None                  # the active one was dropped


def test_batch_loader_skips_non_1d_bundle_entries(tmp_path):
    from larmor import batch, engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 60, 200)
    r = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample="glassA",
               sites=[SiteModel(model="gauss_lor", label="A", params={
                   "isotropic_chemical_shift_ppm": Param(15.0),
                   "shift_fwhm_ppm": Param(6.0), "amplitude": Param(80.0),
                   "gl": Param(1.0, vary=False)})])
    _, y, _ = engine.simulate(r, exp_ppm=x)
    one_d = {"kind": "1d", "title": "glassA", "source_path": "src", "recipe": r.to_dict(),
             "hidden": [], "exp_ppm": x.tolist(), "exp_amp": y.tolist(), "overlays": []}
    two_d = {"kind": "2d", "title": "map", "source_path": "C:/nowhere/2rr",
             "source_rel": None, "recipe": {**r.to_dict(), "sample": "map"},
             "hidden": [], "fittable": True, "view": {"ops": []}}
    fig = {"kind": "figure", "title": "f",
           "spec": {"kind": "1d", "traces": [{"data": {"x": [0, 1], "y": [0, 1]}}]}}
    bat = {"kind": "batch", "title": "b", "state": {"paths": [], "result": None}}
    p = tmp_path / "s.larproj.json"
    p.write_text(json.dumps({"larmor_project_version": 2, "active": 0,
                             "workspaces": [one_d, two_d, fig, bat]}), encoding="utf-8")
    entries, warnings = batch.load_entries([str(p)])
    assert not warnings and len(entries) == 1 and entries[0].sample == "glassA"

    plain = {k: v for k, v in one_d.items() if k != "kind"}   # a v1 file: no 'kind'
    v1 = tmp_path / "v1.larproj.json"
    v1.write_text(json.dumps({"larmor_project_version": 1, "active": 0,
                              "workspaces": [plain, {**plain, "title": "glassB"}]}),
                  encoding="utf-8")
    entries, warnings = batch.load_entries([str(v1)])
    assert not warnings and len(entries) == 2

    # a plain recipe file is NOT a bundle: the reader must say so (None) and
    # fall back to load_any, not migrate it into an empty project
    from larmor import project
    rec = tmp_path / "one.recipe.json"
    rec.write_text(json.dumps(r.to_dict()), encoding="utf-8")
    assert batch._larproj_workspaces(str(rec)) is None
    with pytest.raises(ValueError, match="not a LARMOR project bundle"):
        project.load_bundle(rec)


# ------------------------------------------------------------- 2D replay core
def test_twod_replay_ops_matches_the_direct_operations():
    from larmor import twod

    rng = np.random.default_rng(0)
    f2 = np.linspace(-100, 100, 200)
    f1 = np.linspace(-50, 50, 40)
    z = (np.exp(-((f2[None, :] - 10) / 6) ** 2) * np.exp(-((f1[:, None] - 5) / 6) ** 2)
         + 0.01 * rng.standard_normal((40, 200)))
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=z, nucleus="27Al", larmor_MHz=130.3)
    ops = [{"op": "shear", "factor": 0.5},
           {"op": "phase", "axis": "f2", "p0": 30.0, "p1": 20.0, "pivot": 0.0},
           {"op": "transpose"}, {"op": "rev_f1"},
           {"op": "shift", "d2": 1.0, "d1": -2.0}, {"op": "symmetrize"}]
    got = twod.replay_ops(d, ops)
    ref = (twod.shear(d, 0.5).phased("f2", 30.0, 20.0, 0.0).transposed()
           .reversed_axis("f1").shifted(1.0, -2.0).symmetrized())
    assert np.allclose(got.z, ref.z)
    assert np.allclose(got.f2_ppm, ref.f2_ppm) and np.allclose(got.f1_ppm, ref.f1_ppm)

    s = d.shifted(1.0, -2.0)
    assert np.allclose(s.f2_ppm, f2 + 1.0) and np.allclose(s.f1_ppm, f1 - 2.0)
    assert s.z is d.z
    with pytest.raises(ValueError, match="bogus"):
        twod.replay_ops(d, [{"op": "bogus"}])
    assert twod.replay_ops(d, []) is d


def test_contour_view_records_ops_and_persisted_view_replays(qapp):
    from larmor import project, twod
    from larmor.desktop.twod_view import Contour2DView

    f2 = np.linspace(-100, 100, 64)
    f1 = np.linspace(-40, 40, 16)
    raw = twod.Data2D(f2_ppm=f2, f1_ppm=f1,
                      z=np.random.RandomState(2).rand(16, 64) + 0.1,
                      nucleus="27Al", larmor_MHz=130.3)
    v = Contour2DView()
    v.set_data(raw, "syn"); qapp.processEvents()
    v._op("transpose")
    v._shift_axes(1.0, -2.0)
    v.shearv.setValue(0.3); v._apply_shear()
    v._pick_axis = "f2"; v._pivot = None; v._pref = [object()]
    v.p0v.setValue(30); v._apply_phase()
    assert v._ops[-1] == {"op": "phase", "axis": "f2", "p0": 30.0, "p1": 0.0,
                          "pivot": None}
    v._shift_axes(0.5, 0.0)
    v._reset_phase()
    # the phase before the second calibrate is gone, the calibrate stays
    assert v._ops == [{"op": "transpose"}, {"op": "shift", "d2": 1.0, "d1": -2.0},
                      {"op": "shear", "factor": 0.3},
                      {"op": "shift", "d2": 0.5, "d1": 0.0}]
    assert np.allclose(v.data.z, v._orig.z)

    v.disp.setCurrentText("density"); v.cmap.setCurrentText("magma")
    v.actFlipF1.setChecked(True)
    x1 = v.data.f2_ppm
    amp = np.exp(-((x1 - x1.mean()) / 5) ** 2)
    v.set_projection_1d("f2", x1, amp, color="#e8832a", source="C:/x/1r", scale=2.0)
    assert v._hmqc_scale["f2"] == 2.0
    st = v.get_state()
    assert st["ops"] == v._ops and st["disp"] == "density" and st["cmap"] == "magma"
    assert st["flip"] == {"f2": True, "f1": True}
    assert st["hmqc_source"]["f2"] == "C:/x/1r"
    pv = project.view2d_persisted(st)
    assert set(pv) == {"ops", "sign", "nlevels", "floor", "title", "disp", "cmap",
                       "flip", "projections"}
    assert pv["projections"] == {"f2": {"source": "C:/x/1r", "color": "#e8832a",
                                        "scale": 2.0}}
    json.dumps(pv)

    v2 = Contour2DView()
    v2.set_data(raw, "syn")
    v2.apply_persisted_view(pv); qapp.processEvents()
    assert np.allclose(v2.data.z, v.data.z)
    assert np.allclose(v2.data.f2_ppm, v.data.f2_ppm)
    assert np.allclose(v2.data.f1_ppm, v.data.f1_ppm)
    assert v2._model_ops == ["transpose"] and v2._ops == v._ops
    assert v2.disp.currentText() == "density" and v2.cmap.currentText() == "magma"
    assert v2.nlevels.value() == v.nlevels.value()
    assert v2.sign.currentText() == v.sign.currentText()
    assert v2.floor.value() == v.floor.value()
    assert v2._flip == v._flip and v2.actFlipF1.isChecked()

    # a trailing phase op is live on reopen AND still undone by Reset: _orig
    # follows the pruned log exactly as it would have live
    pv2 = {**pv, "ops": pv["ops"] + [{"op": "phase", "axis": "f2", "p0": 45.0,
                                      "p1": 0.0, "pivot": None}]}
    v3 = Contour2DView()
    v3.set_data(raw, "syn")
    v3.apply_persisted_view(pv2)
    assert not np.allclose(v3.data.z, v3._orig.z)
    assert np.allclose(v3._orig.z, v2.data.z)
    v3._reset_phase()
    assert v3._ops == v2._ops and np.allclose(v3.data.z, v2.data.z)

    # a legacy state dict (no ops / disp / cmap / flip / hmqc_source) applies
    legacy = {k: val for k, val in v.get_state().items()
              if k not in ("ops", "disp", "cmap", "flip", "hmqc_source")}
    v2.set_state(legacy)
    assert v2._ops == [] and v2._hmqc_source == {"f2": "", "f1": ""}


# ------------------------------------------------------------- bundle open / save
def test_open_project_reads_a_v1_file_unchanged(win, monkeypatch, tmp_path):
    """A bundle written by the previous release (no 'kind') still opens:
    every entry is a 1D workspace, the migration is the only note, and no
    Replace / Add / Cancel prompt fires for an empty session."""
    w = win
    path = _v1_file(tmp_path / "old.larproj.json", [("glassA", 15.0), ("glassB", 20.0)])
    w._confirm_open_mode = lambda: pytest.fail("no prompt for an empty session")
    _open(w, monkeypatch, path)
    assert len(w.workspaces) == 2 and all(ws["kind"] == "1d" for ws in w.workspaces)
    assert [ws["snap"]["recipe"]["sample"] for ws in w.workspaces] == ["glassA", "glassB"]
    msg = w.statusBar().currentMessage()
    assert "project opened" in msg and "2 spectra" in msg
    assert len(w._note_calls) == 1 and len(w._note_calls[0]) == 1
    assert "migrated from schema v1 to v2" in w._note_calls[0][0]
    assert w.active_ws == 0 and w.central_stack.currentWidget() is w.view


def test_project_roundtrip_2d_workspace_by_reference(win, monkeypatch, tmp_path):
    from larmor import twod

    path = str(require(BRUKER_2RR_MQMAS))
    w = win
    w.load_source(path, keep_fit=False)
    assert w.central_stack.currentWidget() is w.view2d and len(w.workspaces) == 1
    w.view2d._op("rev_f2")
    w.view2d._shift_axes(1.0, -2.0)
    w.view2d.nlevels.setValue(7)
    w.view2d.sign.setCurrentText("positive")
    w.view2d.disp.setCurrentText("density")
    w.recipe["sites"].append({"model": "czjzek", "label": "Al", "params": {
        "isotropic_chemical_shift_ppm": {"value": 60.0}, "sigma_Cq_MHz": {"value": 2.0},
        "shift_fwhm_ppm": {"value": 8.0}, "amplitude": {"value": 1.0}}})
    w.lines_table.rebuild(w.recipe, w.hidden)

    proj = tmp_path / "p.larproj.json"
    _save(w, monkeypatch, proj)
    d = json.loads(proj.read_text(encoding="utf-8"))
    ws0 = d["workspaces"][0]
    assert ws0["kind"] == "2d" and "z" not in ws0 and "exp_ppm" not in ws0
    assert ws0["source_path"] == path
    assert ws0["view"]["ops"] == [{"op": "rev_f2"}, {"op": "shift", "d2": 1.0, "d1": -2.0}]
    assert ws0["view"]["disp"] == "density"
    assert ws0["recipe"]["sites"][0]["model"] == "czjzek"
    assert proj.stat().st_size < 200_000                 # by reference, no arrays
    assert "1 2D map" in w.statusBar().currentMessage()

    _open(w, monkeypatch, proj)
    assert len(w.workspaces) == 1 and w.workspaces[0]["kind"] == "2d"
    assert w.central_stack.currentWidget() is w.view2d
    assert w.recipe["sites"][0]["model"] == "czjzek" and w._data2d_fittable is True
    assert w.view2d.nlevels.value() == 7 and w.view2d.sign.currentText() == "positive"
    assert w.view2d.disp.currentText() == "density"
    ref = twod.read_bruker_2d(path).normalized()
    assert np.allclose(w.view2d.data.z, np.flip(ref.z, 1))
    assert np.allclose(w.view2d.data.f2_ppm, ref.f2_ppm + 1.0)
    assert np.allclose(w.view2d.data.f1_ppm, ref.f1_ppm - 2.0)
    assert w.view2d._model is None                         # not re-simulated
    assert "run Fit" in _notes_text(w)
    assert w.workspaces[0]["snap"]["view2d"]["ops"] == ws0["view"]["ops"]


def test_open_project_reports_a_missing_2d_source_and_keeps_the_rest(win, monkeypatch, tmp_path):
    w = win
    x = np.linspace(-40, 120, 100)
    y = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    one_d = {"kind": "1d", "title": "glassA", "source_path": "src_a",
             "recipe": {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
                        "sample": "glassA", "sites": []},
             "hidden": [], "exp_ppm": x.tolist(), "exp_amp": y.tolist(), "overlays": []}
    gone = str(tmp_path / "gone" / "35" / "pdata" / "1" / "2rr")
    two_d = {"kind": "2d", "title": "map", "source_path": gone,
             "source_rel": "gone/35/pdata/1/2rr", "recipe": None, "hidden": [],
             "fittable": True, "view": {"ops": []}}
    proj = tmp_path / "p.larproj.json"
    proj.write_text(json.dumps({"larmor_project_version": 2, "active": 1,
                                "workspaces": [one_d, two_d]}), encoding="utf-8")
    _open(w, monkeypatch, proj)
    assert len(w.workspaces) == 1 and w.workspaces[0]["kind"] == "1d"
    assert w.active_ws == 0
    text = _notes_text(w)
    assert "source not found" in text and gone in text

    # a project moved with its data: the project-relative path reaches the
    # loader; an unreadable file there is reported as such, not as missing
    moved = tmp_path / "moved"
    f = moved / "35" / "pdata" / "1" / "2rr"
    f.parent.mkdir(parents=True)
    f.write_bytes(b"not a bruker file")
    proj2 = moved / "p.larproj.json"
    proj2.write_text(json.dumps({"larmor_project_version": 2, "active": 0,
                                 "workspaces": [one_d, {**two_d, "source_rel": "35/pdata/1/2rr"}]}),
                     encoding="utf-8")
    seen = []

    def fake_load(p):
        seen.append(p)
        return False
    monkeypatch.setattr(w, "_load_nonfittable", fake_load)
    w._note_calls.clear()
    _open(w, monkeypatch, proj2)
    assert seen == [str(f)]
    text = _notes_text(w)
    assert "could not be read" in text and str(f) in text
    assert len(w.workspaces) == 1 and w.active_ws == 0
    assert w._ws_mode == "auto"                          # nothing left armed


def test_open_project_offers_replace_add_cancel_when_a_session_is_open(win, monkeypatch, tmp_path):
    w = win
    _mkws(w, "live", 40.0)                     # a fitted workspace is open
    live = w.workspaces[0]
    path = _v1_file(tmp_path / "two.larproj.json", [("glassA", 15.0), ("glassB", 20.0)])

    w._confirm_open_mode = lambda: "add"
    _open(w, monkeypatch, path)
    assert len(w.workspaces) == 3 and w.workspaces[0] is live and live["has_fit"]
    assert w.active_ws == 1                    # saved active 0 + offset 1
    assert [ws["snap"]["recipe"]["sample"] for ws in w.workspaces] == \
        ["live", "glassA", "glassB"]

    w._confirm_open_mode = lambda: "replace"
    _open(w, monkeypatch, path)
    assert [ws["snap"]["recipe"]["sample"] for ws in w.workspaces] == ["glassA", "glassB"]

    before = list(w.workspaces)
    w._confirm_open_mode = lambda: "cancel"
    _open(w, monkeypatch, path)
    assert len(w.workspaces) == 2 and all(a is b for a, b in zip(w.workspaces, before))
    assert "cancelled" in w.statusBar().currentMessage()

    w.workspaces = []
    w.active_ws = None
    w._confirm_open_mode = lambda: pytest.fail("prompt must not fire for an empty session")
    _open(w, monkeypatch, path)
    assert len(w.workspaces) == 2


def test_load_source_routes_a_larproj_file_to_open_project(win, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox

    w = win
    path = _v1_file(tmp_path / "via_open.larproj.json", [("glassA", 15.0), ("glassB", 20.0)])
    w.load_source(path)                        # File ▸ Open… / drag-drop / recent
    assert len(w.workspaces) == 2 and "project opened" in w.statusBar().currentMessage()

    seen = []
    monkeypatch.setattr(w, "_open_project_path", lambda p: seen.append(p))
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: pytest.fail("Load failed box")))
    target = str(tmp_path / "x.larproj.json")
    w.load_source(target)
    assert seen == [target]
    assert "open cancelled" not in w.statusBar().currentMessage()


# ------------------------------------------------------------- batch sessions
def test_batchfit_result_dict_roundtrip_and_alignment():
    from larmor import batchfit, engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 60, 400)
    entries, data = [], []
    for k, (sh, amp) in enumerate(((0.0, 100.0), (0.3, 70.0), (-0.3, 120.0))):
        truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(15.0 + sh),
                "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(truth, exp_ppm=x)
        y = m + np.random.default_rng(k).normal(0, 1.5, x.size)
        start = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=f"g{k}",
                       source_path=f"p{k}", sites=[
                           SiteModel(model="gauss_lor", label="A", params={
                               "isotropic_chemical_shift_ppm": Param(14.0, min=0, max=30),
                               "shift_fwhm_ppm": Param(5.0, min=0.1),
                               "amplitude": Param(80.0, min=0),
                               "gl": Param(1.0, vary=False)})])
        entries.append((start, x, y, (-10.0, 40.0)))
        data.append((x, y, (-10.0, 40.0)))
    res = batchfit.batch_fit(entries)
    batchfit.batch_error_analysis(res, data, method="covariance")

    d = batchfit.result_to_dict(res)
    json.dumps(d)
    back = batchfit.result_from_dict(json.loads(json.dumps(d)))
    assert [_dumps(r) for r in back.recipes] == [_dumps(r) for r in res.recipes]
    assert back.labels == res.labels and back.rmsd == pytest.approx(res.rmsd)
    assert back.shared == res.shared and back.released == res.released
    assert back.release_frac == res.release_frac and back.error_method == "covariance"
    pe0 = res.error_detail["covariance"][0][(0, "amplitude")]
    assert back.error_detail["covariance"][0][(0, "amplitude")].stderr == pe0.stderr
    assert back.per_dataset == [] and back.lmfit_result is None
    key = [(r["scope"], r["site"], r["param"], r["value"]) for r in batchfit.error_table(res)]
    assert [(r["scope"], r["site"], r["param"], r["value"])
            for r in batchfit.error_table(back)] == key

    sub, dropped = batchfit.align_result(back, ["p2", "p0"])
    assert sub.labels == [res.labels[2], res.labels[0]]
    assert sub.rmsd == pytest.approx([res.rmsd[2], res.rmsd[0]])
    assert dropped == ["p1"]
    assert sub.error_detail["covariance"][0] is back.error_detail["covariance"][2]
    curves = batchfit.result_curves(sub, [x[:300], x])
    assert len(curves) == 2
    assert curves[0]["x"].size == curves[0]["y_fit"].size > 0
    assert curves[1]["y_fit"].size == curves[1]["x"].size > 0
    none, dropped = batchfit.align_result(back, ["zz"])
    assert none.recipes == [] and dropped == ["p0", "p1", "p2"]


def test_batchfit_dialog_session_state_roundtrip(qapp, tmp_path):
    from larmor import batchfit, project
    from larmor.desktop.batchfit_dialog import BatchFitDialog, estimate_baseline

    paths = _csv_spectra(tmp_path, 2)
    model = _model()
    dlg = BatchFitDialog(None, paths, model)
    pn = "isotropic_chemical_shift_ppm"
    dlg._rel_checks[pn].setChecked(True)
    dlg.frac.setValue(7); dlg.tol.setValue(0.5)
    dlg.errCombo.setCurrentIndex(dlg.errCombo.findData("montecarlo"))
    dlg.errN.setValue(55)
    dlg._toggle_exclude(1, 1, True)              # site B absent from spectrum 1
    d0 = dlg._data[0]
    base = estimate_baseline(d0["ppm"], d0["amp0"], "Flat (edge median)")
    d0["amp"] = d0["amp0"] - base
    d0["baseline_ops"] = [{"op": "flat_baseline"}]
    dlg._cells[0]["exp"].setData(d0["ppm"], d0["amp"])
    res = batchfit.batch_fit(dlg._entries(), release=(pn,), release_frac=0.07)
    dlg._done(res)

    st = dlg.session_state()
    wire = json.loads(json.dumps(project.json_safe(st)))    # as the file carries it
    assert st["paths"] == paths and st["release"] == [pn]
    assert st["per_spectrum"][paths[1]]["excluded"] == [1]
    assert st["per_spectrum"][paths[0]]["baseline_ops"] == [{"op": "flat_baseline"}]
    assert st["error_method"] == "montecarlo" and st["error_n"] == 55
    assert st["result"] is not None and st["model_sites"] == model["sites"]

    dlg2 = BatchFitDialog(None, paths, {"sites": wire["model_sites"],
                                        "fit_window_ppm": wire["window"]})
    notes = dlg2.apply_session_state(wire)
    assert notes == []
    assert dlg2._rel_checks[pn].isChecked()
    assert not dlg2._rel_checks["shift_fwhm_ppm"].isChecked()
    assert dlg2.frac.value() == 7 and dlg2.tol.value() == 0.5
    assert dlg2.errCombo.currentData() == "montecarlo" and dlg2.errN.value() == 55
    assert dlg2._excluded == {1: {1}}
    assert "excluded" in dlg2._cells[1]["title"].text()
    assert np.allclose(dlg2._data[0]["amp"], dlg._data[0]["amp"])
    assert dlg2.lblBaseline.text() == "none"      # no global kind, a per-spectrum op
    assert dlg2._result is not None
    assert [_dumps(r) for r in dlg2._result.recipes] == [_dumps(r) for r in res.recipes]
    assert dlg2._result.rmsd == pytest.approx(res.rmsd)
    assert dlg2._cells[0]["model"].xData is not None and len(dlg2._cells[0]["model"].xData) > 0
    assert dlg2._cells[0]["rmsd"].text().endswith(f"{res.rmsd[0]:.4f}")
    assert dlg2.btnSave.isEnabled() and dlg2.btnTable.isEnabled() and dlg2.btnErr.isEnabled()
    assert dlg2.table.rowCount() == 2 and dlg2.table.columnCount() > 4
    assert dlg2.session_title().startswith("batch fit · 2 spectra")
    # the state round-trips through a second dialog unchanged
    st2 = dlg2.session_state()
    assert {k: v for k, v in json.loads(json.dumps(project.json_safe(st2))).items()
            if k != "result"} == {k: v for k, v in wire.items() if k != "result"}


def test_batchfit_dialog_session_state_keeps_series_names_order_and_columns(qapp, tmp_path):
    """A renamed, reordered series with a typed column survives the session
    round trip: the rows ride in paths order with their source_path, come
    back before the result is re-paired, and a spectrum that cannot be
    reloaded loses only its own row."""
    from larmor import batchfit, project
    from larmor.desktop.batchfit_dialog import BatchFitDialog
    from larmor.series_table import SeriesColumn

    paths = _csv_spectra(tmp_path, 3)
    model = _model()
    dlg = BatchFitDialog(None, paths, model)
    tbl = dlg._series.permuted([1, 0, 2])
    tbl.rename(1, "glassZ")                               # spectrum 0, now second
    tbl.add_column(SeriesColumn("CaO (mol%)", "CaO (mol%)"), {0: 1.0, 1: 0.0, 2: 2.0})
    dlg._apply_series(tbl, [1, 0, 2])
    assert [d["path"] for d in dlg._data] == [paths[1], paths[0], paths[2]]
    assert [d["sample"] for d in dlg._data] == ["batch01", "glassZ", "batch02"]
    res = batchfit.batch_fit(dlg._entries())
    dlg._done(res)
    assert res.labels == ["batch01", "glassZ", "batch02"]

    st = dlg.session_state()
    wire = json.loads(json.dumps(project.json_safe(st)))
    assert wire["paths"] == [paths[1], paths[0], paths[2]]
    assert [r["source_path"] for r in wire["series"]["rows"]] == wire["paths"]
    assert wire["series"]["rows"][1]["display_name"] == "glassZ"
    assert wire["series"]["columns"][0]["key"] == "CaO (mol%)"
    dlg2 = BatchFitDialog(None, wire["paths"], {"sites": wire["model_sites"],
                                                "fit_window_ppm": wire["window"]})
    assert dlg2.apply_session_state(wire) == []
    assert dlg2._series.labels() == ["batch01", "glassZ", "batch02"]
    assert [d["sample"] for d in dlg2._data] == ["batch01", "glassZ", "batch02"]
    assert [r.values["CaO (mol%)"] for r in dlg2._series.rows] == [1.0, 0.0, 2.0]
    assert dlg2._result.labels == ["batch01", "glassZ", "batch02"]
    assert [r.sample for r in dlg2._result.recipes] == dlg2._result.labels
    assert dlg2._cells[1]["title"].text().startswith("glassZ")
    assert dlg2.table.item(dlg2._row_of(1), 1).text() == "glassZ"
    assert dlg2._result.rmsd == pytest.approx(res.rmsd)
    assert dlg2.session_state()["series"] == wire["series"]

    os.remove(paths[1])
    dlg3 = BatchFitDialog(None, wire["paths"], {"sites": wire["model_sites"]})
    notes = dlg3.apply_session_state(wire)
    assert notes and "could not be reloaded" in notes[0]
    assert dlg3._series.labels() == ["glassZ", "batch02"]
    assert [r.values["CaO (mol%)"] for r in dlg3._series.rows] == [0.0, 2.0]
    assert dlg3._result.labels == ["glassZ", "batch02"]
    assert [d["sample"] for d in dlg3._data] == ["glassZ", "batch02"]


def test_batchfit_dialog_session_aligns_results_when_a_spectrum_is_missing(qapp, tmp_path):
    from larmor import batchfit
    from larmor.desktop.batchfit_dialog import BatchFitDialog

    paths = _csv_spectra(tmp_path, 3)
    model = _model()
    dlg = BatchFitDialog(None, paths, model)
    res = batchfit.batch_fit(dlg._entries())
    dlg._done(res)
    st = dlg.session_state()
    os.remove(paths[1])

    dlg3 = BatchFitDialog(None, paths, {"sites": st["model_sites"],
                                        "fit_window_ppm": st["window"]})
    notes = dlg3.apply_session_state(st)
    assert len(dlg3._data) == 2 and dlg3._result is not None
    assert [r.source_path for r in dlg3._result.recipes] == [paths[0], paths[2]]
    assert dlg3._result.rmsd == pytest.approx([res.rmsd[0], res.rmsd[2]])
    assert notes and "could not be reloaded" in notes[0]
    assert "could not be reloaded" in dlg3.status.text()
    assert dlg3._cells[1]["rmsd"].text().endswith(f"{res.rmsd[2]:.4f}")

    os.remove(paths[2])
    dlg4 = BatchFitDialog(None, paths, {"sites": st["model_sites"]})
    dlg4.apply_session_state(st)
    assert len(dlg4._data) == 1 and dlg4._result is not None
    assert len(dlg4._result.recipes) == 1

    os.remove(paths[0])
    dlg5 = BatchFitDialog(None, paths, {"sites": st["model_sites"]})
    notes5 = dlg5.apply_session_state(st)
    assert dlg5._data == [] and dlg5._result is None and notes5


def test_project_roundtrip_figures_and_batch_sessions(win, monkeypatch, tmp_path):
    w = win
    _mkws(w, "glassA", 60.0)
    spec = {"kind": "1d", "traces": [{"data": {"x": [1, 2, 3], "y": [1, 2, 3]},
                                      "label": "t"}], "title": "Fig A"}
    w._remember_figure(spec)
    w._remember_figure({"kind": "2d", "path": ""})          # empty: never kept
    st = {"paths": [str(tmp_path / "a.csv"), str(tmp_path / "b.csv")], "release": [],
          "frac": 10.0, "tol": 0.1, "baseline": "None", "components": False,
          "shared_scale": False, "model_sites": _model()["sites"], "window": None,
          "recipe_tag": "", "baseline_kind": "None", "per_spectrum": {},
          "error_method": "covariance", "error_n": 200, "auto_recipes": True,
          "table": True, "result": None}
    w._add_session_entry("batch", "batch fit · 2 spectra", {"state": st})
    assert [ws["kind"] for ws in w.workspaces] == ["1d", "figure", "batch"]
    assert w.active_ws == 0
    texts = [w.ws_panel.list.item(i).text() for i in range(w.ws_panel.list.count())]
    assert [t[0] for t in texts] == ["⤳", "◫", "☷"]

    proj = tmp_path / "s.larproj.json"
    _save(w, monkeypatch, proj)
    d = json.loads(proj.read_text(encoding="utf-8"))
    assert [ws["kind"] for ws in d["workspaces"]] == ["1d", "figure", "batch"]
    assert d["workspaces"][1]["spec"] == spec and d["active"] == 0
    saved = d["workspaces"][2]["state"]
    assert {k: v for k, v in saved.items() if k != "paths_rel"} == st
    assert saved["paths_rel"] == ["a.csv", "b.csv"]
    msg = w.statusBar().currentMessage()
    assert "1 spectrum" in msg and "1 figure" in msg and "1 batch session" in msg

    _open(w, monkeypatch, proj)
    assert [ws["kind"] for ws in w.workspaces] == ["1d", "figure", "batch"]
    assert w.workspaces[1]["snap"]["spec"] == spec
    assert w.workspaces[2]["snap"]["state"]["paths"] == st["paths"]
    assert w.active_ws == 0 and w.central_stack.currentWidget() is w.view
    text = _notes_text(w)                     # the missing spectra, reported once
    assert "batch session" in text and "2 spectrum file(s) not found" in text
    assert len(w._note_calls) == 1


def test_studio_close_records_only_changed_specs_and_updates_in_place(win, monkeypatch):
    from larmor.desktop import figure_dialog as fd
    from larmor.desktop import plotting_studio as ps

    w = win
    _mkws(w, "glassA", 60.0)
    n0 = len(w.workspaces)

    def exec_add_trace(self):
        self._push_trace({"data": {"x": [0, 1], "y": [0, 1]}, "label": "t"})
        self.title.setText("T1")
        return 0
    monkeypatch.setattr(ps.PlottingStudio, "exec", exec_add_trace)
    w.open_plotting_studio()
    assert len(w.workspaces) == n0 + 1
    assert w.workspaces[-1]["kind"] == "figure" and w.workspaces[-1]["title"] == "T1"

    idx = len(w.workspaces) - 1
    monkeypatch.setattr(ps.PlottingStudio, "exec",
                        lambda self: (self.title.setText("T2"), 0)[1])
    w.open_workspace_entry(idx)                         # reopened: updated in place
    assert len(w.workspaces) == n0 + 1 and w.workspaces[idx]["title"] == "T2"
    assert w.workspaces[idx]["snap"]["spec"]["title"] == "T2"
    assert len(w.workspaces[idx]["snap"]["spec"]["traces"]) == 1

    monkeypatch.setattr(ps.PlottingStudio, "exec", lambda self: 0)
    w.open_plotting_studio({"kind": "2d", "path": ""})   # empty seed, untouched
    assert len(w.workspaces) == n0 + 1
    w.plot_current_spectrum()                             # seeded, untouched
    assert len(w.workspaces) == n0 + 1
    monkeypatch.setattr(ps.PlottingStudio, "exec",
                        lambda self: (self.title.setText("T3"), 0)[1])
    w.plot_current_spectrum()                             # seeded, then changed
    assert len(w.workspaces) == n0 + 2 and w.workspaces[-1]["title"] == "T3"

    fig_spec = {"kind": "1d", "traces": [{"data": {"x": [0, 1], "y": [1, 0]},
                                          "label": "u"}], "title": "F"}
    monkeypatch.setattr(fd.FigureDialog, "exec",
                        lambda self: (self.spec.setPlainText(json.dumps(fig_spec)), 0)[1])
    w.open_figure_dialog()
    assert len(w.workspaces) == n0 + 3 and w.workspaces[-1]["snap"]["spec"] == fig_spec
    monkeypatch.setattr(fd.FigureDialog, "exec", lambda self: 0)
    w.open_figure_dialog()
    assert len(w.workspaces) == n0 + 3
    assert w.active_ws == 0

    # a figure kind the studio cannot edit reopens in the Figure studio
    w._add_session_entry("figure", "inf", {"spec": {"kind": "infinite_field",
                                                     "recipe": "x.json"}})
    kinds = []
    monkeypatch.setattr(fd.FigureDialog, "exec",
                        lambda self: (kinds.append(json.loads(self.spec.toPlainText())["kind"]), 0)[1])
    w.open_workspace_entry(len(w.workspaces) - 1)
    assert kinds == ["infinite_field"]


def test_run_batch_fit_records_a_session_only_when_it_produced_a_result(win, monkeypatch, tmp_path):
    from larmor import batchfit
    from larmor.desktop import batchfit_dialog as bd

    w = win
    paths = _csv_spectra(tmp_path, 2)
    monkeypatch.setattr(bd.BatchFitDialog, "exec", lambda self: 0)
    w.run_batch_fit(paths)
    assert w.workspaces == []                              # no fit, nothing kept

    def exec_fit(self):
        self._done(batchfit.batch_fit(self._entries()))
        return 0
    monkeypatch.setattr(bd.BatchFitDialog, "exec", exec_fit)
    w.recipe = _model()
    w.run_batch_fit(paths)
    assert len(w.workspaces) == 1 and w.workspaces[0]["kind"] == "batch"
    st = w.workspaces[0]["snap"]["state"]
    assert st["result"] is not None and st["paths"] == paths
    assert w.workspaces[0]["title"].startswith("batch fit · 2 spectra")
    assert w.active_ws is None
    rmsd = list(st["result"]["rmsd"])

    monkeypatch.setattr(bd.BatchFitDialog, "exec", lambda self: 0)
    w.open_workspace_entry(0)                              # reopened, results restored
    assert len(w.workspaces) == 1
    assert w.workspaces[0]["snap"]["state"]["result"]["rmsd"] == pytest.approx(rmsd)

    for p in paths:                                        # every spectrum gone
        os.remove(p)
    w.open_workspace_entry(0)
    assert w.workspaces[0]["snap"]["state"]["result"] is not None   # not overwritten
    assert "unchanged" in w.statusBar().currentMessage()


def test_session_entries_never_reach_apply_doc(win, monkeypatch):
    from larmor.desktop import plotting_studio as ps

    w = win
    _mkws(w, "glassA", 60.0)
    spec = {"kind": "1d", "traces": [{"data": {"x": [0, 1], "y": [0, 1]}, "label": "t"}],
            "title": "F"}
    w._remember_figure(spec)
    assert [ws["kind"] for ws in w.workspaces] == ["1d", "figure"] and w.active_ws == 0

    w.switch_workspace(1)                       # a click only selects a session row
    assert w.active_ws == 0 and "Enter" in w.statusBar().currentMessage()

    w.close_workspace(0)                        # the only document goes away
    assert w.active_ws is None and len(w.workspaces) == 1
    assert w.workspaces[0]["kind"] == "figure"

    _mkws(w, "glassB", 30.0)                    # index 1, active 1
    assert w.active_ws == 1
    n = {"count": 0}
    monkeypatch.setattr(ps.PlottingStudio, "exec",
                        lambda self: n.__setitem__("count", n["count"] + 1) or 0)
    w.save_workspace(0)                         # the Save button reopens a session
    assert n["count"] == 1 and w.active_ws == 1
    w.open_workspace_entry(0)
    assert n["count"] == 2 and w.active_ws == 1
    assert w.ws_panel.list.currentRow() == 1

    w.close_workspace(0)                        # the figure row, while active is 1
    assert w.active_ws == 0 and w.workspaces[0]["kind"] == "1d"
    assert w.recipe["sample"] == "glassB"       # the document was not re-applied


def test_batchfit_session_state_replays_reprocess_on_reopen(qapp, tmp_path, monkeypatch):
    """A batch reprocessed from its fids reopens reprocessed: the session
    state carries the template and each spectrum's chain + EXPNO, and
    apply_session_state replays the chain BEFORE the baseline and the
    stored result; a member whose fid is gone keeps its TopSpin spectrum
    with a note."""
    from conftest import fake_spectrum_params
    from larmor import comparability, project
    from larmor.desktop.batchfit_dialog import BatchFitDialog

    paths = _csv_spectra(tmp_path, 2)
    model = _model()
    expnos = []
    for k in range(2):
        e = tmp_path / f"E{k}" / "24"
        e.mkdir(parents=True)
        (e / "acqus").write_text("##TITLE= stub\n##END=\n", encoding="utf-8")
        expnos.append(str(e))
    fakes = {f"batch0{k}": fake_spectrum_params(f"batch0{k}", lb=100 * k, has_fid=True,
                                                expno=expnos[k]) for k in range(2)}
    monkeypatch.setattr(comparability, "read_params",
                        lambda p, procno=None: fakes.get(Path(str(p)).stem))
    x = np.linspace(-20, 60, 500)
    # a 0.7 pedestal so the flat (edge-median) baseline replayed below has
    # something to remove
    fake_amp = (3.0 * np.exp(-((x - 15.0) / 5.0) ** 2)
                + 1.0 * np.exp(-((x - 2.0) / 3.0) ** 2) + 0.7)
    monkeypatch.setattr(comparability, "reprocess",
                        lambda params, ops: (x.copy(), fake_amp.copy(), ["replayed"]))
    dlg = BatchFitDialog(None, paths, model)
    assert "LB 0 / 100 Hz" in dlg.compBar.text()
    template = {"wdw": 0, "lb_hz": 0.0, "gb": 0.0, "ssb": 0.0, "tdeff": 1024,
                "si": 32768, "fcor": 1.0}
    dlg._apply_reprocess(template, "own")
    ops0 = dlg._data[0]["proc_ops"]
    assert ops0 and dlg._data[0]["expno"] == expnos[0]
    # a baseline over the reprocessed spectrum
    d0 = dlg._data[0]
    d0["amp"] = d0["amp0"] - 0.5
    d0["baseline_ops"] = [{"op": "flat_baseline"}]

    st = dlg.session_state()
    wire = json.loads(json.dumps(project.json_safe(st)))
    assert wire["reprocess"] == {"template": template, "phase": "own"}
    assert wire["per_spectrum"][paths[0]]["proc_ops"] == ops0
    assert wire["per_spectrum"][paths[0]]["expno"] == expnos[0]
    assert wire["per_spectrum"][paths[0]]["baseline_ops"] == [{"op": "flat_baseline"}]

    dlg2 = BatchFitDialog(None, paths, {"sites": wire["model_sites"],
                                        "fit_window_ppm": wire["window"]})
    notes = dlg2.apply_session_state(wire)
    assert notes == []
    assert dlg2._reprocess == wire["reprocess"]
    assert np.allclose(dlg2._data[0]["amp0"], fake_amp)          # the chain replayed
    assert dlg2._data[0]["proc_ops"] == ops0 and dlg2._data[0]["expno"] == expnos[0]
    assert dlg2._data[0]["baseline_ops"] == [{"op": "flat_baseline"}]
    assert np.allclose(dlg2._data[0]["amp"], fake_amp - 0.7, atol=1e-6)   # baseline on top
    assert np.allclose(dlg2._cells[0]["exp"].yData, dlg2._data[0]["amp"])
    assert "reprocessed from fid" in dlg2.compBar.text() and "2/2 spectra" in dlg2.compBar.text()
    assert not dlg2.compBar.btnRevert.isHidden()
    assert not dlg2._cells[1]["title"].text().endswith("⚠")       # LB no longer differs
    rec = dlg2._entries()[0][0]
    assert rec.processing_from_raw and rec.source_path == expnos[0]
    assert rec.processing == ops0 + [{"op": "flat_baseline"}]
    st2 = json.loads(json.dumps(project.json_safe(dlg2.session_state())))
    assert st2["reprocess"] == wire["reprocess"]
    assert st2["per_spectrum"] == wire["per_spectrum"]

    # the fid is gone: the TopSpin spectrum is kept and the note says so
    nofid = {k: fake_spectrum_params(k, lb=100 * i, has_fid=False)
             for i, k in enumerate(("batch00", "batch01"))}
    monkeypatch.setattr(comparability, "read_params",
                        lambda p, procno=None: nofid.get(Path(str(p)).stem))
    dlg3 = BatchFitDialog(None, paths, {"sites": wire["model_sites"],
                                        "fit_window_ppm": wire["window"]})
    notes = dlg3.apply_session_state(wire)
    assert len(notes) == 2 and all("could not reprocess" in n for n in notes)
    assert np.allclose(dlg3._data[1]["amp"], dlg3._data[1]["amp_src"])
    assert dlg3._data[1]["proc_ops"] == [] and dlg3._data[1]["expno"] == ""
    assert "could not reprocess" in dlg3.status.text()
    assert dlg3._cells[1]["title"].text().endswith("⚠")           # still flagged
