"""larmor.provenance: the software stamp, the file-read git commit, the
source hash, the reload checks and the read-out notes -- and the loader /
fit plumbing that writes and reads them."""
import datetime as _dt
import json
import re
import subprocess

import numpy as np
import pytest

from conftest import BRUKER_1R, require
from larmor import fit as fitmod
from larmor import engine, provenance as P
from larmor.loader import load_any
from larmor.recipe import Param, Recipe, SiteModel, sha256_of
from test_acquisition import make_expno


def test_software_stamp_and_git_commit(tmp_path, monkeypatch):
    import importlib.metadata as md

    import larmor

    s = P.software_stamp()
    assert set(s) >= {"larmor", "git_commit", "mrsimulator", "lmfit", "numpy", "scipy",
                      "nmrglue", "python", "platform", "fitted"}
    assert s["larmor"] == larmor.__version__
    assert s["mrsimulator"] == md.version("mrsimulator") and s["scipy"] == md.version("scipy")
    _dt.datetime.fromisoformat(s["fitted"])
    s2 = P.software_stamp()
    assert {k: v for k, v in s.items() if k != "fitted"} == \
        {k: v for k, v in s2.items() if k != "fitted"}
    assert s["git_commit"] == "" or re.fullmatch(r"[0-9a-f]{12}|[0-9a-f]{40}", s["git_commit"])

    # the file reader: never a subprocess
    def boom(*a, **k):
        raise AssertionError("git_commit must not spawn a process")
    monkeypatch.setattr(subprocess, "run", boom)
    c = P.git_commit()
    assert c == "" or re.fullmatch(r"[0-9a-f]{40}", c)
    h = "a" * 40
    # HEAD -> loose ref
    repo = tmp_path / "r1"
    (repo / ".git" / "refs" / "heads").mkdir(parents=True)
    (repo / ".git" / "HEAD").write_text("ref: refs/heads/master\n")
    (repo / ".git" / "refs" / "heads" / "master").write_text(h + "\n")
    assert P.git_commit(repo) == h
    # packed-refs only
    repo2 = tmp_path / "r2"
    (repo2 / ".git").mkdir(parents=True)
    (repo2 / ".git" / "HEAD").write_text("ref: refs/heads/wip/N4\n")
    (repo2 / ".git" / "packed-refs").write_text(
        "# pack-refs with: peeled fully-peeled sorted\n" + "b" * 40 + " refs/heads/wip/N4\n")
    assert P.git_commit(repo2) == "b" * 40
    # detached
    repo3 = tmp_path / "r3"
    (repo3 / ".git").mkdir(parents=True)
    (repo3 / ".git" / "HEAD").write_text("c" * 40 + "\n")
    assert P.git_commit(repo3) == "c" * 40
    # a worktree: .git FILE 'gitdir:' + commondir holding the refs
    main = tmp_path / "main" / ".git"
    (main / "refs" / "heads").mkdir(parents=True)
    (main / "refs" / "heads" / "wt").write_text("d" * 40 + "\n")
    wt_git = main / "worktrees" / "wt"
    wt_git.mkdir(parents=True)
    (wt_git / "HEAD").write_text("ref: refs/heads/wt\n")
    (wt_git / "commondir").write_text("../..\n")
    wt = tmp_path / "wt"
    wt.mkdir()
    (wt / ".git").write_text(f"gitdir: {wt_git}\n")
    assert P.git_commit(wt) == "d" * 40
    assert P.git_commit(tmp_path / "nogit") == ""
    assert P.parse_version("0.13.0") == (0, 13, 0) and P.parse_version("0.14.0rc1") == (0, 14, 0)
    assert P.parse_version("") == ()


def test_recipe_roundtrip_keeps_acquisition_and_software_without_version_bump(tmp_path):
    from larmor.recipe import RECIPE_VERSION

    r = Recipe(nucleus="31P", acquisition={"ns": 14, "d1_s": 300.0},
               software={"larmor": "0.13.0"})
    p = tmp_path / "a.recipe.json"
    r.save(p)
    d = json.loads(p.read_text(encoding="utf-8"))
    assert d["larmor_recipe_version"] == 1 == RECIPE_VERSION
    back = Recipe.load(p)
    assert back.acquisition == {"ns": 14, "d1_s": 300.0} and back.software == {"larmor": "0.13.0"}
    old = {k: v for k, v in d.items() if k not in ("acquisition", "software")}
    r2 = Recipe.from_dict(old)
    assert r2.acquisition == {} and r2.software == {} and r2.notes == []


def _gauss_recipe(pos=10.0):
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(pos), "shift_fwhm_ppm": Param(5.0),
            "amplitude": Param(100.0), "gl": Param(1.0, vary=False)})])


def test_fit_stamps_software_block():
    import larmor
    from larmor.fithealth import recipe_signature

    truth = _gauss_recipe(10.0)
    x = np.linspace(-30, 50, 600)
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    rec = _gauss_recipe(8.0)
    assert rec.software == {}
    before = recipe_signature(rec)
    res = fitmod.fit(rec, x, y)
    assert res.recipe.software["larmor"] == larmor.__version__
    assert "fitted" in res.recipe.software and res.recipe.software["mrsimulator"]
    # the stamp is not part of the 'edited since the fit' signature
    after = recipe_signature(res.recipe)
    assert len(after) == len(before)
    assert all(a[0] == b[0] for a, b in zip(after, before))          # same models
    assert "software" in res.recipe.to_dict()


def _csv(tmp_path, name="s.csv"):
    x = np.linspace(-20, 60, 200)
    y = 100 * np.exp(-0.5 * ((x - 15) / 3) ** 2)
    p = tmp_path / name
    p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n"
                 + "".join(f"{a:.4f} {b:.4f}\n" for a, b in zip(x, y)), encoding="utf-8")
    return p


def test_loader_fills_sha_and_block_and_verifies_on_reload(tmp_path):
    csvp = _csv(tmp_path)
    ppm, amp, rec, meta, warns = load_any(csvp)
    assert rec["source_sha256"] == sha256_of(csvp) and re.fullmatch(r"[0-9a-f]{64}", rec["source_sha256"])
    assert rec["acquisition"] == {} and rec["source_kind"] == "csv"
    rp = tmp_path / "s.recipe.json"
    Recipe.from_dict(rec).save(rp)
    _, _, _, _, w = load_any(rp)
    assert not any(P.SOURCE_CHANGED_PREFIX in x for x in w)
    # the data file changes after the fit -> one status-bar line, no exception
    with open(csvp, "a", encoding="utf-8") as f:
        f.write("70.0000 0.0000\n")
    _, _, _, _, w = load_any(rp)
    assert any(x.startswith(P.SOURCE_CHANGED_PREFIX) for x in w) and len(w) == 1
    assert "s.csv" in w[0] and "refit before publishing" in w[0]
    # a pre-0.14 recipe (no hash) stays silent
    d = dict(rec, source_sha256="")
    rp2 = tmp_path / "old.recipe.json"
    Recipe.from_dict(d).save(rp2)
    _, _, _, _, w = load_any(rp2)
    assert not any(P.SOURCE_CHANGED_PREFIX in x for x in w)
    # SR moved since the fit
    d3 = dict(rec, source_sha256=sha256_of(csvp), sr_hz=50.0)
    rp3 = tmp_path / "sr.recipe.json"
    Recipe.from_dict(d3).save(rp3)
    _, _, _, _, w = load_any(rp3)
    assert any(x.startswith(P.REREFERENCED_PREFIX) for x in w)
    assert "SR 50.00 → 0.00 Hz" in " ".join(w) and "ppm axis moved" in " ".join(w)
    # ... unless the recipe explains it (the audit re-referenced it, or an sr op)
    d4 = dict(d3, provenance={"referencing": {"old_sr_hz": 0.0, "new_sr_hz": 50.0}})
    assert P.verify_source(d4, {"sr_hz": 0.0}) == []
    d5 = dict(d3, processing=[{"op": "sr", "sr_hz": 50.0}])
    assert P.verify_source(d5, {"sr_hz": 0.0}) == []
    assert P.verify_source({}, None) == [] and P.verify_source(None, {"sr_hz": 1.0}) == []
    assert P.carry_source({"source_kind": "csv", "source_sha256": "ab" * 32}) == {
        "source_kind": "csv", "source_sha256": "ab" * 32, "acquisition": {}}
    assert P.carry_source({}) == {"source_kind": "", "source_sha256": "", "acquisition": {}}

    # the real 1r: kind, hash of the exact file, block with the title's flip
    e = require(BRUKER_1R)
    _, _, rec, _, _ = load_any(e)
    assert rec["source_kind"] == "bruker" and rec["source_sha256"] == sha256_of(e)
    a = rec["acquisition"]
    assert a["procno"] == 1 and a["data_file"] == "pdata/1/1r" and a["nucleus"] == "27Al"
    assert a["flip_deg"] == 11.0 and a["ns"] == 56
    assert P.data_file_for(rec) == e
    assert P.source_sha256(rec) == rec["source_sha256"]
    # File > Open FID fits hash the fid
    assert P.data_file_for(rec, from_raw=True) == e.parents[2] / "fid"
    assert P.data_file_for({}) is None and P.source_sha256("Z:/nowhere/1r") == ""


def test_loader_json_reopens_the_fitted_procno(tmp_path):
    e = make_expno(tmp_path / "2026-05", "04272026_RS43009_S1_SS_ALP", 2701, procnos=(1, 2))
    _, amp1, rec1, _, _ = load_any(e)
    _, amp2, rec2, _, w2 = load_any(e / "pdata" / "2")
    assert amp1.max() == pytest.approx(1000.0) and amp2.max() == pytest.approx(2000.0)
    assert rec2["acquisition"]["procno"] == 2 and rec2["acquisition"]["data_file"] == "pdata/2/1r"
    assert rec2["source_sha256"] == sha256_of(e / "pdata" / "2" / "1r")
    assert rec2["source_path"] == str(e)          # source_path stays the EXPNO
    rp = tmp_path / "p2.recipe.json"
    Recipe.from_dict(rec2).save(rp)
    _, amp, rec, _, w = load_any(rp)
    assert amp.max() == pytest.approx(2000.0) and w == []
    assert rec["acquisition"]["procno"] == 2
    rp1 = tmp_path / "p1.recipe.json"
    Recipe.from_dict(rec1).save(rp1)
    _, amp, _, _, w = load_any(rp1)
    assert amp.max() == pytest.approx(1000.0) and w == []
    # a recipe of procno 2 whose 1r is later reprocessed: the hash line
    (e / "pdata" / "2" / "1r").write_bytes((np.zeros(1024, dtype="<i4") + 7).tobytes())
    _, _, _, _, w = load_any(rp)
    assert any(x.startswith(P.SOURCE_CHANGED_PREFIX) and "1r" in x for x in w)


def test_version_notes_only_for_applicable_readout_changes():
    cz = {"software": {"larmor": "0.12.1"}, "sites": [{"model": "czjzek"}]}
    notes = P.version_notes(cz, "0.13.0")
    assert len(notes) == 1 and "fitted with LARMOR 0.12.1" in notes[0]
    assert "Czjzek read-outs" in notes[0] and "0.13.0" in notes[0]
    assert P.version_notes({"software": {"larmor": "0.12.1"},
                            "sites": [{"model": "gauss_lor"}]}, "0.13.0") == []
    assert P.version_notes({"software": {"larmor": "0.13.0"},
                            "sites": [{"model": "czjzek"}]}, "0.13.0") == []
    assert P.version_notes({"sites": [{"model": "czjzek"}]}, "0.13.0") == []
    assert P.version_notes({"software": {"larmor": "0.13.0"},
                            "sites": [{"model": "czjzek"}]}, "0.14.0") == []
    assert P.version_notes({"software": {"larmor": "0.12.1"},
                            "sites": [{"model": "ext_czjzek"}]}, "0.14.0")
    # a newer-than-current recipe is not warned about (a downgrade is not a read-out change)
    assert P.version_notes({"software": {"larmor": "0.15.0"},
                            "sites": [{"model": "czjzek"}]}, "0.13.0") == []
    assert P.version_notes(cz) == P.version_notes(cz, __import__("larmor").__version__)
