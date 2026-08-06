"""Resolving a batch of fits from a CSV, a folder, or explicit paths, for the
Plotting studio's batch-grid figures (larmor.series_grid)."""
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


def test_load_panels_via_csv_uses_source_path_hint_when_no_recipe_matches(tmp_path):
    """A newer CSV's own source_path column resolves a scope to a data-only
    panel (experiment, no fit) even with zero .recipe.json files around."""
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
    assert p.data_path == str(raw) and p.path == ""
    assert p.models == ("gauss_lor",)
    # resolved via the source_path hint -- nothing left to ask the user about
    assert not any("couldn't locate data" in w for w in warnings)


def test_resolve_manual_clears_needs_manual(tmp_path):
    from larmor import series_grid as sg
    panel = sg.Panel(path="", sample="g3", nucleus="", models=(), has_data=False,
                     needs_manual=True)
    fixed = sg.resolve_manual(panel, tmp_path / "found_it.csv")
    assert fixed.has_data and not fixed.needs_manual
    assert fixed.data_path == str(tmp_path / "found_it.csv")
    assert panel.needs_manual   # original untouched (dataclasses.replace copies)


def test_csv_rows_by_scope_groups_correctly(tmp_path):
    csv_path = tmp_path / "t.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr\n"
        "g0,s0,A,amplitude,100,2\ng0,s1,B,amplitude,50,1\ng1,s0,A,amplitude,90,3\n")
    grouped = series_grid.csv_rows_by_scope(csv_path)
    assert set(grouped) == {"g0", "g1"}
    assert len(grouped["g0"]) == 2 and len(grouped["g1"]) == 1
