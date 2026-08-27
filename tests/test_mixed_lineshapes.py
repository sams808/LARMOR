"""Every fit/batch-fit/reconstruction test elsewhere in this suite uses a
recipe where every site shares ONE model (gauss_lor throughout, or czjzek
throughout, or quad_ct throughout). That's the one place a subtle
"worked for my one dataset, silently wrong for another" bug would hide --
this file exercises a recipe that MIXES three genuinely different models
(gauss_lor, czjzek, quad_ct — three disjoint parameter-name sets, one
kernel-based) through the full pipeline: a single fit, a batch fit (shared
model + per-spectrum "Exclude component" on an arbitrary site), a CSV
export/reconstruction round-trip, and a batch-grid figure.
"""
import numpy as np
import pytest

from larmor import batchfit, engine, series_grid
from larmor import fit as fitmod
from larmor.recipe import Param, Recipe, SiteModel


def _mixed_recipe(sample: str = "") -> Recipe:
    """gauss_lor (no kernel, no quadrupole), czjzek (kernel-based), quad_ct
    (analytic quadrupolar, no kernel) -- three disjoint param-name sets:
    gauss_lor has "gl", czjzek has "sigma_Cq_MHz", quad_ct has "Cq_MHz"/"eta".
    """
    return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sample=sample,
                  sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(20.0, min=0, max=40),
            "shift_fwhm_ppm": Param(4.0, min=0.1),
            "amplitude": Param(50.0, min=0), "gl": Param(1.0, vary=False)}),
        SiteModel(model="czjzek", label="B", params={
            "isotropic_chemical_shift_ppm": Param(60.0, min=40, max=80),
            "sigma_Cq_MHz": Param(2.5, min=0.2, max=8.0),
            "shift_fwhm_ppm": Param(6.0, min=0.1, max=30),
            "amplitude": Param(80.0, min=0)}),
        SiteModel(model="quad_ct", label="C", params={
            "isotropic_chemical_shift_ppm": Param(-10.0, min=-30, max=10),
            "Cq_MHz": Param(3.0, min=0.01, max=40),
            "eta": Param(0.3, min=0, max=1),
            "shift_fwhm_ppm": Param(3.0, min=0.05),
            "amplitude": Param(40.0, min=0)}),
    ])


def _mixed_spectrum(recipe: Recipe, x: np.ndarray, seed: int,
                    scale: tuple[float, float, float] = (1.0, 1.0, 1.0)) -> np.ndarray:
    """A synthetic spectrum from `recipe`, with each site's amplitude
    scaled independently (so a batch of these isn't just one spectrum
    repeated) — kernel-based sites (czjzek) simulate on their own native
    grid, so the standard np.interp-back-onto-x step (also used by
    test_batchfit.py's _czjzek_entries) is needed here too."""
    rec = Recipe.from_dict(recipe.to_dict())
    for site, s in zip(rec.sites, scale):
        site.params["amplitude"].value *= s
    xr, y, _ = engine.simulate(rec, exp_ppm=x)
    y = np.interp(x, np.sort(xr), y[np.argsort(xr)])
    return y + np.random.default_rng(seed).normal(0, 0.4, x.size)


def test_fit_recovers_all_three_models_with_no_cross_talk():
    """larmor.fit.fit(): one recipe, three models, no shared/colliding param
    names -- every site's own parameters converge near their true values."""
    x = np.linspace(-60, 100, 1200)
    truth = _mixed_recipe()
    data = _mixed_spectrum(truth, x, seed=0)

    start = _mixed_recipe()
    start.sites[0].params["isotropic_chemical_shift_ppm"].value = 22.0
    start.sites[1].params["sigma_Cq_MHz"].value = 2.0
    start.sites[2].params["Cq_MHz"].value = 3.5

    res = fitmod.fit(start, x, data)
    assert res.rmsd < 0.05
    assert start.sites[0].params["isotropic_chemical_shift_ppm"].value == pytest.approx(20.0, abs=0.5)
    assert start.sites[0].params["gl"].value == 1.0            # untouched (fixed)
    assert start.sites[1].params["sigma_Cq_MHz"].value == pytest.approx(2.5, abs=0.4)
    assert start.sites[2].params["Cq_MHz"].value == pytest.approx(3.0, abs=0.5)
    assert start.sites[2].params["eta"].value == pytest.approx(0.3, abs=0.15)
    # no model's params leaked onto another site
    assert "gl" not in start.sites[1].params and "gl" not in start.sites[2].params
    assert "Cq_MHz" not in start.sites[0].params and "Cq_MHz" not in start.sites[1].params
    assert "sigma_Cq_MHz" not in start.sites[0].params and "sigma_Cq_MHz" not in start.sites[2].params


def _mixed_batch_entries(n=3, exclude_site_on=None):
    """`n` synthetic spectra sharing one mixed-model recipe (shape held,
    amplitudes free per spectrum, per batchfit's design) -- optionally with
    one site's amplitude pre-locked to exactly zero on ONE spectrum (the
    "Exclude component" signature batchfit.free_amplitudes must respect,
    batchfit_dialog.py's right-click action normally sets this)."""
    x = np.linspace(-60, 100, 900)
    truth = _mixed_recipe()
    entries = []
    for k in range(n):
        scale = (1.0 + 0.3 * k, 1.0 - 0.15 * k, 1.0 + 0.1 * k)
        data = _mixed_spectrum(truth, x, seed=k, scale=scale)
        rec = _mixed_recipe(sample=f"g{k}")
        if exclude_site_on is not None and k == exclude_site_on[0]:
            amp = rec.sites[exclude_site_on[1]].params["amplitude"]
            amp.value, amp.vary, amp.min, amp.max = 0.0, False, 0.0, 0.0
        entries.append((rec, x, data, None))
    return entries


def test_batch_fit_shares_shape_and_frees_amplitude_across_mixed_models():
    entries = _mixed_batch_entries()
    res = batchfit.batch_fit(entries)
    # shape params held fixed at the recipe value (nothing released) for
    # EVERY model, not just whichever model a single-model test happens to use
    for rec in res.recipes:
        assert rec.sites[0].params["shift_fwhm_ppm"].value == pytest.approx(4.0)
        assert rec.sites[1].params["sigma_Cq_MHz"].value == pytest.approx(2.5)
        assert rec.sites[2].params["eta"].value == pytest.approx(0.3)
    # amplitudes adapted per spectrum for every model, tracking the injected scale
    ampA = [r.sites[0].params["amplitude"].value for r in res.recipes]
    ampB = [r.sites[1].params["amplitude"].value for r in res.recipes]
    ampC = [r.sites[2].params["amplitude"].value for r in res.recipes]
    assert ampA[2] > ampA[0]     # gauss_lor amplitude grew (scale 1.0 -> 1.2)
    assert ampB[2] < ampB[0]     # czjzek amplitude shrank (scale 1.0 -> 0.7)
    assert ampC[2] > ampC[0]     # quad_ct amplitude grew (scale 1.0 -> 1.2)
    shared = set(batchfit.all_but_amplitude(res.recipes))
    assert {"gl", "sigma_Cq_MHz", "Cq_MHz", "eta"} <= shared   # every model's shape params


def test_batch_fit_exclusion_works_regardless_of_which_models_site():
    """The zero-amplitude "Exclude component" lock (see is_zeroed_out) must
    be respected identically whether the excluded site is gauss_lor, czjzek,
    or quad_ct -- free_amplitudes() must not special-case by model name."""
    for site_idx, model in enumerate(("gauss_lor", "czjzek", "quad_ct")):
        entries = _mixed_batch_entries(exclude_site_on=(1, site_idx))
        res = batchfit.batch_fit(entries)
        excluded_amp = res.recipes[1].sites[site_idx].params["amplitude"]
        assert excluded_amp.value == 0.0 and not excluded_amp.vary, model
        # the other two spectra's SAME site fit normally (not locked)
        assert res.recipes[0].sites[site_idx].params["amplitude"].vary, model
        assert res.recipes[2].sites[site_idx].params["amplitude"].vary, model


def test_csv_round_trip_preserves_each_sites_own_model(tmp_path):
    """batchfit's table export + series_grid's reconstruction must carry each
    site's OWN model name through -- not assume one model for the whole CSV
    (unlike a converter that assumes one model for the whole file, which could
    safely assume gauss_lor only because that happened to be true there)."""
    entries = _mixed_batch_entries()
    # give each spectrum a resolvable source_path (reconstruction needs one to
    # seed nucleus/larmor_frequency_MHz; without it the recipe stays data-less
    # but should still reconstruct correctly parameter-wise)
    for rec, x, data, _win in entries:
        p = tmp_path / f"{rec.sample}_raw.csv"
        p.write_text("# nucleus = 27Al\n# larmor_MHz = 130.3\n" +
                     "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, data)))
        rec.source_path = str(p)
    res = batchfit.batch_fit(entries)
    rows = batchfit.shared_table(res)

    csv_path = tmp_path / "batch_table.csv"
    import csv as _csv
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = _csv.writer(f)
        w.writerow(["scope", "site", "label", "param", "value", "stderr",
                    "model", "source_path"])
        for r in rows:
            w.writerow([r["scope"], r["site"], r["label"], r["param"],
                       r["value"], r["stderr"] or "", r["model"],
                       r.get("source_path", "")])

    rows_by_scope = series_grid.csv_rows_by_scope(csv_path)
    rebuilt = series_grid.recipe_from_csv_rows(
        rows_by_scope["shared"], rows_by_scope["g0"], rows_by_scope["g0"][0]["source_path"])
    assert [s.model for s in rebuilt.sites] == ["gauss_lor", "czjzek", "quad_ct"]
    assert set(rebuilt.sites[0].params) == {
        "isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "amplitude", "gl"}
    # line_fwhm_ppm: czjzek's own extra Lorentzian-broadening param, never
    # set in _mixed_recipe() (soft default 0.0 at render time) -- filled in
    # from the model registry's own default, not rejected as "missing"
    assert set(rebuilt.sites[1].params) == {
        "isotropic_chemical_shift_ppm", "sigma_Cq_MHz", "shift_fwhm_ppm",
        "line_fwhm_ppm", "amplitude"}
    assert rebuilt.sites[1].params["line_fwhm_ppm"].value == 0.0
    assert set(rebuilt.sites[2].params) == {
        "isotropic_chemical_shift_ppm", "Cq_MHz", "eta", "shift_fwhm_ppm", "amplitude"}

    panels, warnings = series_grid.load_panels(str(csv_path))
    assert all(p.reconstructed for p in panels)
    # resolved via CSV reconstruction -- nothing left to ask the user about,
    # even though the informational "no .recipe.json next to this CSV" note
    # still fires (true: none exist; reconstruction is a separate, successful
    # fallback path, same as series_grid's other CSV-only-resolution tests)
    assert not any("couldn't locate data" in w for w in warnings)


def test_render_batch_grid_handles_a_mixed_model_panel(tmp_path):
    """The figure engine (colors, peak labels, quantify()-based population %)
    must not assume every component shares one model."""
    import matplotlib.pyplot as plt
    from larmor import figures

    x = np.linspace(-60, 100, 900)
    rec = _mixed_recipe(sample="g0")
    data = _mixed_spectrum(rec, x, seed=0)
    fitmod.fit(rec, x, data)     # gives it a real fit_window_ppm for quantify()
    p = tmp_path / "g0_raw.csv"
    p.write_text("# nucleus = 27Al\n# larmor_MHz = 130.3\n" +
                 "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, data)))
    rec.source_path = str(p)
    recipe_path = tmp_path / "g0.recipe.json"
    rec.save(recipe_path)

    fig = figures.render({"kind": "batch_grid",
                          "panels": [{"recipe": str(recipe_path)}],
                          "peak_labels": "position+pct"})
    assert fig.legends and {t.get_text() for t in fig.legends[0].get_texts()} == {
        "A", "B", "C"}
    plt.close(fig)
