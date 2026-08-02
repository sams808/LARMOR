"""Regressions for the batch publication-report tool (larmor.batch).

Guards: (a) the model-aware columns and their error propagation (Czjzek derived
quantities), (b) the full pipeline writes a table + report with populations and
errors, (c) the Monte-Carlo error method is wired, (d) mixed nuclei are flagged.
"""
import csv
import json

import numpy as np
import pytest

from larmor import batch, engine
from larmor.recipe import Recipe, SiteModel, Param


def _write_fit(tmp_path, name, pos1=15.0, nucleus="11B"):
    x = np.linspace(-20, 60, 900)
    r = Recipe(nucleus=nucleus, larmor_frequency_MHz=160.0, spin_rate_Hz=20000.0,
               sample=name, source_kind="csv", sites=[
                   SiteModel(model="gauss_lor", label="BO3", params={
                       "isotropic_chemical_shift_ppm": Param(pos1),
                       "shift_fwhm_ppm": Param(6.0, min=0.1),
                       "amplitude": Param(100.0, min=0.0),
                       "gl": Param(1.0, vary=False)}),
                   SiteModel(model="gauss_lor", label="BO4", params={
                       "isotropic_chemical_shift_ppm": Param(0.0),
                       "shift_fwhm_ppm": Param(3.0, min=0.1),
                       "amplitude": Param(50.0, min=0.0),
                       "gl": Param(1.0, vary=False)})])
    _, model, _ = engine.simulate(r, exp_ppm=x)
    data = model + np.random.default_rng(0).normal(0.0, 2.0, x.size)
    csvp = tmp_path / f"{name}.csv"
    with open(csvp, "w", encoding="utf-8") as f:
        f.write(f"# nucleus = {nucleus}\n# larmor_MHz = 160\n")
        for xi, yi in zip(x, data):
            f.write(f"{xi:.5f} {yi:.5f}\n")
    r.source_path = str(csvp)
    r.fit_window_ppm = (-10.0, 40.0)
    p = tmp_path / f"{name}.recipe.json"
    p.write_text(json.dumps(r.to_dict()), encoding="utf-8")
    return str(p)


def test_czjzek_columns_propagate_errors():
    site = {"model": "czjzek", "params": {
        "isotropic_chemical_shift_ppm": {"value": 60.0, "stderr": 0.2},
        "sigma_Cq_MHz": {"value": 1.6, "stderr": 0.05},
        "shift_fwhm_ppm": {"value": 8.0, "stderr": 0.3}}}
    errs = {"isotropic_chemical_shift_ppm": 0.2, "sigma_Cq_MHz": 0.05,
            "shift_fwhm_ppm": 0.3}
    cols = dict((h, (v, e)) for h, v, e in batch._site_columns(site, errs))
    assert cols["C_Q=2σ (MHz)"] == pytest.approx((3.2, 0.1))
    v, e = cols["√⟨P_Q²⟩ (MHz)"]
    assert v == pytest.approx(np.sqrt(5) * 1.6)
    assert e == pytest.approx(np.sqrt(5) * 0.05)


def test_run_batch_writes_table_and_report(tmp_path):
    paths = [_write_fit(tmp_path, "glassA", 15.0),
             _write_fit(tmp_path, "glassB", 14.5)]
    out = tmp_path / "report"
    res = batch.run_batch(paths, out, error_method="covariance",
                          make_plots=True, formats=("csv", "latex", "markdown"))
    assert res.n_fits == 2 and res.n_sites == 4
    assert (out / "table.csv").exists()
    assert (out / "table.tex").exists()
    assert (out / "report.md").exists()
    assert (out / "figures" / "glassA.png").exists()

    # the CSV holds a value AND an error column for δiso and a population
    with open(out / "table.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 4
    assert any(r["δiso (ppm)"] for r in rows)
    assert any(r["pop (%)"] for r in rows)
    assert any(r["δiso (ppm) err"] for r in rows)          # covariance error present

    md = (out / "report.md").read_text(encoding="utf-8")
    assert "| Sample |" in md and "pop (%)" in md
    assert "figures/glassA.png" in md


@pytest.mark.slow
def test_montecarlo_method_runs(tmp_path):
    paths = [_write_fit(tmp_path, "mcA", 15.0)]
    res = batch.run_batch(paths, tmp_path / "mc", error_method="montecarlo",
                          n_mc=20, make_plots=False, formats=("csv",))
    assert res.error_method == "montecarlo"
    assert res.n_sites == 2


def test_homogeneity_flags_mixed_nuclei(tmp_path):
    a = _write_fit(tmp_path, "al", 60.0, nucleus="27Al")
    b = _write_fit(tmp_path, "bo", 15.0, nucleus="11B")
    entries, _ = batch.load_entries([a, b])
    notes = batch.homogeneity(entries)
    assert any("mixed nuclei" in n for n in notes)
