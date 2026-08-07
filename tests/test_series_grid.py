"""Resolving a batch of fits from a CSV, a folder, or explicit paths, for the
Plotting studio's batch-grid figures (larmor.series_grid)."""
from pathlib import Path

import numpy as np
import pytest

from larmor import engine, series_grid
from larmor.recipe import Recipe, SiteModel, Param


def _make_fit(tmp_path, sample, has_data=True, seed=0):
    """A saved .recipe.json, optionally with a real source spectrum on disk
    (source_path set) so has_data reflects a genuine round-trip, not a flag."""
    sites = [SiteModel(model="gauss_lor", label="A", params={
        "isotropic_chemical_shift_ppm": Param(10.0), "shift_fwhm_ppm": Param(4.0),
        "amplitude": Param(100.0), "gl": Param(1.0, vary=False)})]
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=sample,
                sites=sites)
    if has_data:
        x = np.linspace(-30, 30, 200)
        _, y, _ = engine.simulate(rec, exp_ppm=x)
        data = y + np.random.default_rng(seed).normal(0, 0.5, x.size)
        csv_path = tmp_path / f"{sample}_raw.csv"
        csv_path.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                            "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, data)))
        rec.source_path = str(csv_path)
    p = tmp_path / f"{sample}.recipe.json"
    rec.save(p)
    return p


def test_resolve_paths_accepts_list_folder_or_csv(tmp_path):
    p0 = _make_fit(tmp_path, "g0")
    p1 = _make_fit(tmp_path, "g1")

    paths, warn = series_grid.resolve_paths([str(p0), str(p1)])
    assert paths == [str(p0), str(p1)] and not warn

    paths2, warn2 = series_grid.resolve_paths(str(tmp_path))
    assert set(paths2) == {str(p0), str(p1)} and not warn2

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    paths3, warn3 = series_grid.resolve_paths(str(empty_dir))
    assert paths3 == [] and warn3


def test_find_recipes_near_csv_matches_by_scope(tmp_path):
    p0 = _make_fit(tmp_path, "g0")
    p1 = _make_fit(tmp_path, "g1")
    _make_fit(tmp_path, "g2")            # not referenced by the csv below

    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr\n"
        "g0,s0,A,amplitude,100,2\n"
        "g1,s0,A,amplitude,90,3\n"
        "shared,s0,A,shift_fwhm_ppm,4.0,\n")

    matched = series_grid.find_recipes_near_csv(csv_path)
    assert set(matched) == {str(p0), str(p1)}   # g2 correctly excluded, "shared" ignored


def test_load_panels_reports_has_data_and_falls_back_gracefully(tmp_path):
    _make_fit(tmp_path, "g0", has_data=True)
    _make_fit(tmp_path, "g1", has_data=False)     # no source_path at all

    panels, warnings = series_grid.load_panels(str(tmp_path))
    by_sample = {p.sample: p for p in panels}
    assert by_sample["g0"].has_data is True
    assert by_sample["g1"].has_data is False
    assert by_sample["g0"].models == ("gauss_lor",)
    assert by_sample["g0"].n_sites == 1
    assert not warnings


def test_load_panels_via_csv_auto_retrieves_the_saved_fits(tmp_path):
    """The exact "load the csv, auto-retrieve sample+exp+fit" workflow."""
    _make_fit(tmp_path, "g0", has_data=True, seed=1)
    _make_fit(tmp_path, "g1", has_data=True, seed=2)
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr\n"
        "g0,s0,A,amplitude,100,2\ng1,s0,A,amplitude,90,3\n")

    panels, warnings = series_grid.load_panels(str(csv_path))
    assert {p.sample for p in panels} == {"g0", "g1"}
    assert all(p.has_data for p in panels)
    assert not warnings


def test_load_panels_via_csv_with_no_matching_fits_warns_usefully(tmp_path):
    """An older CSV (no source_path/model columns) with no matching fit next
    to it: the scope isn't silently dropped -- it comes back needs_manual so
    the studio can prompt the user to locate its data directly."""
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text("scope,site,label,param,value,stderr\ng0,s0,A,amplitude,100,2\n")
    panels, warnings = series_grid.load_panels(str(csv_path))
    assert len(panels) == 1
    assert panels[0].sample == "g0" and panels[0].needs_manual and not panels[0].has_data
    assert any("Save individual fits" in w for w in warnings)
    assert any("g0" in w for w in warnings)


def test_load_panels_via_csv_reconstructs_even_a_partial_row_set(tmp_path):
    """A newer CSV with a `model` column but only "amplitude" for a site
    (every other param used a soft render-time default, so the original
    recipe never had a Param for it either -- see
    recipe_from_csv_rows' registry-default fill-in) still reconstructs a
    full fit, not just a data-only panel -- even with zero .recipe.json
    around."""
    x = np.linspace(-20, 20, 50)
    raw = tmp_path / "g0_raw.csv"
    raw.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                   "\n".join(f"{xi:.4f} {xi:.4f}" for xi in x))
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr,model,source_path\n"
        f"g0,s0,A,amplitude,100,2,gauss_lor,{raw}\n")
    panels, warnings = series_grid.load_panels(str(csv_path))
    assert len(panels) == 1
    p = panels[0]
    assert p.sample == "g0" and p.has_data and not p.needs_manual
    assert p.reconstructed and p.path.endswith(".recipe.json")
    assert p.models == ("gauss_lor",)
    # resolved via the source_path hint -- nothing left to ask the user about
    assert not any("couldn't locate data" in w for w in warnings)


def _complete_csv(tmp_path, raw):
    """A batch CSV with a full gauss_lor row set (shared gl + per-scope
    position/width/amplitude, the real shape a batch-fit export has) for
    two scopes sharing one "shared" ladder -- everything needed for a full
    reconstruction."""
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr,model,source_path\n"
        "shared,s0,A,gl,1,,gauss_lor,\n"
        f"g0,s0,A,isotropic_chemical_shift_ppm,10.0,0.1,gauss_lor,{raw}\n"
        f"g0,s0,A,shift_fwhm_ppm,5.0,0.2,gauss_lor,{raw}\n"
        f"g0,s0,A,amplitude,100.0,3.0,gauss_lor,{raw}\n")
    return csv_path


def test_recipe_from_csv_rows_rebuilds_a_complete_fit(tmp_path):
    x = np.linspace(-30, 30, 100)
    raw = tmp_path / "g0_raw.csv"
    raw.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                   "\n".join(f"{xi:.4f} {xi:.4f}" for xi in x))
    csv_path = _complete_csv(tmp_path, raw)
    rows = series_grid.csv_rows_by_scope(csv_path)
    rec = series_grid.recipe_from_csv_rows(rows["shared"], rows["g0"], str(raw))
    assert len(rec.sites) == 1
    site = rec.sites[0]
    assert site.model == "gauss_lor"
    assert site.params["isotropic_chemical_shift_ppm"].value == pytest.approx(10.0)
    assert site.params["isotropic_chemical_shift_ppm"].stderr == pytest.approx(0.1)
    assert site.params["gl"].value == pytest.approx(1.0)
    assert rec.nucleus == "11B" and rec.source_path == str(raw)


def test_recipe_from_csv_rows_fills_missing_params_from_the_model_registry():
    """A site the CSV only gave "amplitude" for (e.g. every other param used
    a soft render-time default and so was never a Param the original recipe
    even carried, like czjzek's line_fwhm_ppm) reconstructs successfully --
    the missing params take the model registry's own bootstrap defaults,
    exactly what a freshly-added site of that model starts from, rather than
    being rejected as "incomplete" (a real bug a mixed-lineshape test caught:
    czjzek's line_fwhm_ppm is exactly this case)."""
    rec = series_grid.recipe_from_csv_rows(
        [], [{"site": "s0", "label": "A", "param": "amplitude",
             "value": "100", "stderr": "", "model": "gauss_lor"}])
    assert len(rec.sites) == 1
    p = rec.sites[0].params
    assert p["amplitude"].value == 100.0
    assert p["isotropic_chemical_shift_ppm"].value == 0.0    # registry default
    assert p["shift_fwhm_ppm"].value == 5.0                  # registry default
    assert p["gl"].value == 1.0                              # registry default


def test_recipe_from_csv_rows_rejects_older_csv_without_model():
    with pytest.raises(ValueError, match="model"):
        series_grid.recipe_from_csv_rows(
            [], [{"site": "s0", "label": "A", "param": "amplitude",
                 "value": "100", "stderr": "", "model": ""}])


def test_load_panels_via_csv_reconstructs_full_fits_when_model_column_present(tmp_path):
    """The exact new workflow: a batch CSV alone (no .recipe.json anywhere)
    is enough for a full deconvolution-grid panel -- experiment AND fit."""
    x = np.linspace(-30, 30, 100)
    raw = tmp_path / "g0_raw.csv"
    raw.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                   "\n".join(f"{xi:.4f} {xi:.4f}" for xi in x))
    csv_path = _complete_csv(tmp_path, raw)
    panels, warnings = series_grid.load_panels(str(csv_path))
    assert len(panels) == 1
    p = panels[0]
    assert p.reconstructed and p.has_data and not p.needs_manual
    assert p.path.endswith(".recipe.json") and Path(p.path).exists()
    assert p.n_sites == 1 and p.models == ("gauss_lor",)

    from larmor import figures
    fig = figures.render({"kind": "batch_grid", "panels": [{"recipe": p.path}]})
    import matplotlib.pyplot as plt
    assert len(fig.axes) >= 1
    plt.close(fig)


def test_resolve_manual_clears_needs_manual(tmp_path):
    from larmor import series_grid as sg
    panel = sg.Panel(path="", sample="g3", nucleus="", models=(), has_data=False,
                     needs_manual=True)
    fixed = sg.resolve_manual(panel, tmp_path / "found_it.csv")
    assert fixed.has_data and not fixed.needs_manual
    assert fixed.data_path == str(tmp_path / "found_it.csv")
    assert panel.needs_manual   # original untouched (dataclasses.replace copies)


def test_recipe_from_csv_rows_drops_a_site_absent_from_this_scope():
    """A site that only exists under "shared" (e.g. excluded/zeroed for this
    sample via batch fit's "Exclude component") must be OMITTED from the
    reconstructed recipe, not treated as an incomplete/broken site."""
    shared = [
        {"scope": "shared", "site": "s0", "label": "A", "param": "gl",
         "value": "1", "stderr": "", "model": "gauss_lor"},
        {"scope": "shared", "site": "s1", "label": "B", "param": "gl",
         "value": "1", "stderr": "", "model": "gauss_lor"},
    ]
    # g0 has both sites; g1 (built below) only reports s0 -- s1 was excluded
    scope_rows_g1 = [
        {"scope": "g1", "site": "s0", "label": "A",
         "param": "isotropic_chemical_shift_ppm", "value": "10.0", "stderr": "",
         "model": "gauss_lor"},
        {"scope": "g1", "site": "s0", "label": "A", "param": "shift_fwhm_ppm",
         "value": "5.0", "stderr": "", "model": "gauss_lor"},
        {"scope": "g1", "site": "s0", "label": "A", "param": "amplitude",
         "value": "100.0", "stderr": "", "model": "gauss_lor"},
    ]
    rec = series_grid.recipe_from_csv_rows(shared, scope_rows_g1)
    assert len(rec.sites) == 1 and rec.sites[0].label == "A"


def test_recipe_from_csv_rows_ignores_population_pct_rows():
    """population_pct (batchfit's derived integrated-population column) is
    not a model parameter -- it must not be treated as one, or fed into the
    site-completeness check, when rebuilding a fit from the CSV."""
    shared = [{"scope": "shared", "site": "s0", "label": "A", "param": "gl",
              "value": "1", "stderr": "", "model": "gauss_lor"}]
    scope_rows = [
        {"scope": "g0", "site": "s0", "label": "A",
         "param": "isotropic_chemical_shift_ppm", "value": "10.0", "stderr": "",
         "model": "gauss_lor"},
        {"scope": "g0", "site": "s0", "label": "A", "param": "shift_fwhm_ppm",
         "value": "5.0", "stderr": "", "model": "gauss_lor"},
        {"scope": "g0", "site": "s0", "label": "A", "param": "amplitude",
         "value": "100.0", "stderr": "", "model": "gauss_lor"},
        {"scope": "g0", "site": "s0", "label": "A", "param": "population_pct",
         "value": "100.0", "stderr": "", "model": "gauss_lor"},
    ]
    rec = series_grid.recipe_from_csv_rows(shared, scope_rows)
    assert len(rec.sites) == 1
    assert "population_pct" not in rec.sites[0].params


def test_load_panels_via_csv_with_excluded_site_reconstructs_the_rest(tmp_path):
    """End-to-end: a real batch-fit CSV export (via batchfit's own table
    builders) with one scope excluding a site reconstructs correctly through
    load_panels -- the excluded site simply isn't part of that scope's panel."""
    from larmor import engine, batchfit
    from larmor.recipe import Recipe, SiteModel, Param

    x = np.linspace(-30, 30, 200)

    def start(sample):
        return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=sample,
                      sites=[
                          SiteModel(model="gauss_lor", label="A", params={
                              "isotropic_chemical_shift_ppm": Param(10.0),
                              "shift_fwhm_ppm": Param(4.0),
                              "amplitude": Param(80.0, min=0),
                              "gl": Param(1.0, vary=False)}),
                          SiteModel(model="gauss_lor", label="B", params={
                              "isotropic_chemical_shift_ppm": Param(-8.0),
                              "shift_fwhm_ppm": Param(3.0),
                              "amplitude": Param(40.0, min=0),
                              "gl": Param(1.0, vary=False)})])

    entries = []
    raws = {}
    for k, sample in enumerate(["g0", "g1"]):
        rec = start(sample)
        _, y, _ = engine.simulate(rec, exp_ppm=x)
        data = y + np.random.default_rng(k).normal(0, 0.5, x.size)
        raw = tmp_path / f"{sample}_raw.csv"
        raw.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                       "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, data)))
        raws[sample] = raw
        rec.source_path = str(raw)
        entries.append((rec, x, data, None))
    entries[1][0].sites[1].params["amplitude"] = Param(
        0.0, vary=False, min=0.0, max=0.0)   # exclude site B for g1

    res = batchfit.batch_fit(entries)
    rows = batchfit.shared_table(res)
    csv_path = tmp_path / "batch_table.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        import csv as _csv
        w = _csv.writer(f)
        w.writerow(["scope", "site", "label", "param", "value", "stderr",
                    "model", "source_path"])
        for r in rows:
            w.writerow([r["scope"], r["site"], r["label"], r["param"],
                       r["value"], r["stderr"] or "", r["model"],
                       r.get("source_path", "")])

    panels, warnings = series_grid.load_panels(str(csv_path))
    by_sample = {p.sample: p for p in panels}
    assert by_sample["g0"].n_sites == 2 and by_sample["g0"].reconstructed
    assert by_sample["g1"].n_sites == 1 and by_sample["g1"].reconstructed
    rec_g1 = Recipe.load(by_sample["g1"].path)
    assert rec_g1.sites[0].label == "A"


def test_csv_rows_by_scope_groups_correctly(tmp_path):
    csv_path = tmp_path / "t.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr\n"
        "g0,s0,A,amplitude,100,2\ng0,s1,B,amplitude,50,1\ng1,s0,A,amplitude,90,3\n")
    grouped = series_grid.csv_rows_by_scope(csv_path)
    assert set(grouped) == {"g0", "g1"}
    assert len(grouped["g0"]) == 2 and len(grouped["g1"]) == 1
