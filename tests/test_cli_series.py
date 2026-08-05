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


def test_cli_seqfit_needs_model(tmp_path):
    paths, _ = _make_series(tmp_path)
    # strip the CSV headers so no nucleus/model is inferable? still has no sites →
    # no model available without --model
    rc = cli.main(["seqfit", *paths, "-o", str(tmp_path / "x")])
    assert rc == 1
