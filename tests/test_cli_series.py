"""Headless CLI: `larmor batchfit` and `larmor seqfit` over a series."""
import numpy as np

from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine, cli


def _make_series(tmp_path):
    x = np.linspace(-20, 60, 400)
    paths = []
    for k, pos in enumerate([14.0, 15.0, 16.0]):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                    sites=[SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(pos),
                        "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100),
                        "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        p = tmp_path / f"s{k}.csv"
        p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n"
                     + "".join(f"{xi:.4f} {yi:.4f}\n" for xi, yi in zip(x, m)))
        paths.append(str(p))
    model = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                   sites=[SiteModel(model="gauss_lor", label="A", params={
                       "isotropic_chemical_shift_ppm": Param(12.0, min=0, max=30),
                       "shift_fwhm_ppm": Param(5.0, min=0.1),
                       "amplitude": Param(80.0, min=0),
                       "gl": Param(1.0, vary=False)})])
    mp = tmp_path / "model.recipe.json"
    model.save(mp)
    return paths, str(mp)


def test_cli_batchfit_writes_outputs(tmp_path):
    paths, model = _make_series(tmp_path)
    out = tmp_path / "bout"
    rc = cli.main(["batchfit", *paths, "--model", model, "-o", str(out)])
    assert rc == 0
    assert (out / "batch_table.csv").exists()
    assert len(list(out.glob("*_batch.recipe.json"))) == 3

    # N5: the five status columns end the header; shared rows are held by
    # the batch (vary False), amplitudes free, population rows blank
    import csv
    with open(out / "batch_table.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert list(rows[0].keys())[-5:] == ["vary", "min", "max", "expr", "at_bound"]
    shared = [r for r in rows if r["scope"] == "shared"]
    assert {r["param"] for r in shared} == {"isotropic_chemical_shift_ppm",
                                             "shift_fwhm_ppm", "gl"}
    assert all(r["vary"] == "False" and r["at_bound"] == "" for r in shared)
    assert next(r for r in shared if r["param"] == "shift_fwhm_ppm")["min"] == "0.1"
    amps = [r for r in rows if r["param"] == "amplitude"]
    assert len(amps) == 3
    assert all(r["vary"] == "True" and r["min"] == "0" and r["max"] == "" for r in amps)
    pops = [r for r in rows if r["param"] == "population_pct"]
    assert len(pops) == 3
    assert all(r["vary"] == "" and r["expr"] == "" and r["at_bound"] == "" for r in pops)


def test_cli_seqfit_beats_shared_on_marching_series(tmp_path):
    paths, model = _make_series(tmp_path)
    out = tmp_path / "sout"
    rc = cli.main(["seqfit", *paths, "--model", model, "--passes", "2",
                   "-o", str(out)])
    assert rc == 0
    assert (out / "seq_table.csv").exists()
    recs = sorted(out.glob("*_seq.recipe.json"))
    assert len(recs) == 3
    pos = [Recipe.load(str(r)).sites[0].params["isotropic_chemical_shift_ppm"].value
           for r in recs]
    # each spectrum found its own marching position (12/14/16-ish, not one shared)
    assert max(pos) - min(pos) > 1.0

    # N5: seq_table.csv carries the same status columns; every position is
    # free inside its recipe bounds (0..30) -- no ±window, nothing at a bound
    import csv
    with open(out / "seq_table.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert list(rows[0].keys())[-5:] == ["vary", "min", "max", "expr", "at_bound"]
    prow = [r for r in rows if r["param"] == "isotropic_chemical_shift_ppm"]
    assert len(prow) == 3
    assert all(r["vary"] == "True" and r["at_bound"] == "" for r in prow)
    assert all((r["min"], r["max"]) == ("0", "30") for r in prow)


def test_cli_seqfit_needs_model(tmp_path):
    paths, _ = _make_series(tmp_path)
    # strip the CSV headers so no nucleus/model is inferable? still has no sites →
    # no model available without --model
    rc = cli.main(["seqfit", *paths, "-o", str(tmp_path / "x")])
    assert rc == 1


def test_cli_batchfit_and_seqfit_curves_flag_write_bundle_files(tmp_path, capsys):
    """--curves adds the publication bundle (curves, manifest.csv, README.txt)
    next to the recipes/table the CLI already writes, reusing their stems;
    the CLI table now carries the dialog's 8 columns (model, source_path)."""
    import csv
    from larmor import batchfit

    paths, model = _make_series(tmp_path)
    out = tmp_path / "bout"
    rc = cli.main(["batchfit", *paths, "--model", model, "-o", str(out), "--curves"])
    assert rc == 0
    printed = capsys.readouterr().out
    assert "wrote 3 recipe(s) + batch_table.csv to" in printed
    assert f"wrote 3 curve file(s) + manifest.csv + README.txt to {out}" in printed

    with open(out / "batch_table.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == batchfit.SHARED_HEADER
    per = [r for r in rows[1:] if r[0] != "shared"]
    assert per and all(r[7] in paths for r in per)          # source_path filled
    assert all(r[6] == "gauss_lor" for r in per)
    recipes = sorted(p.name for p in out.glob("*.recipe.json"))
    assert recipes == [f"s{k}_batch.recipe.json" for k in range(3)]   # written once
    for k in range(3):
        assert Recipe.load(out / f"s{k}_batch.recipe.json").source_path == paths[k]

    with open(out / "manifest.csv", newline="", encoding="utf-8") as f:
        man = list(csv.DictReader(f))
    assert len(man) == 3
    assert [m["recipe_file"] for m in man] == [f"s{k}_batch.recipe.json" for k in range(3)]
    assert [m["curves_file"] for m in man] == [f"s{k}_batch_curves.csv" for k in range(3)]
    assert [m["source_path"] for m in man] == paths
    assert all((out / m["curves_file"]).is_file() for m in man)
    assert all(m["source_kind"] == "csv" and len(m["source_sha256"]) == 64 for m in man)
    assert all(m["note"] == "" and float(m["rmsd"]) > 0 for m in man)
    assert (out / "README.txt").is_file()

    # without --curves nothing of the bundle appears
    out2 = tmp_path / "plain"
    assert cli.main(["batchfit", *paths, "--model", model, "-o", str(out2)]) == 0
    assert not (out2 / "manifest.csv").exists()
    assert not list(out2.glob("*_curves.csv")) and not (out2 / "README.txt").exists()

    sout = tmp_path / "sout"
    rc = cli.main(["seqfit", *paths, "--model", model, "--passes", "1",
                   "-o", str(sout), "--curves"])
    assert rc == 0
    assert (sout / "seq_table.csv").is_file() and (sout / "manifest.csv").is_file()
    assert (sout / "README.txt").is_file()
    assert sorted(p.name for p in sout.glob("*_seq_curves.csv")) == \
        [f"s{k}_seq_curves.csv" for k in range(3)]
    assert len(list(sout.glob("*.recipe.json"))) == 3
    with open(sout / "manifest.csv", newline="", encoding="utf-8") as f:
        sman = list(csv.DictReader(f))
    assert [m["recipe_file"] for m in sman] == [f"s{k}_seq.recipe.json" for k in range(3)]
    assert "sequential" in (sout / "README.txt").read_text(encoding="utf-8")
