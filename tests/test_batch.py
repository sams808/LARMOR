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


def _write_fit(tmp_path, name, pos1=15.0, nucleus="11B", bo3_width_start=None,
               fix_bo4_pos=False):
    """A synthetic two-line 11B fit (truth: BO3 width 6 ppm, BO4 at 0 ppm).
    ``bo3_width_start=w`` writes the RECIPE's BO3 width as ``Param(w, min=w)``
    after the data are made from the truth, so a fit from it must stop at
    that floor; ``fix_bo4_pos`` holds the BO4 position. Defaults reproduce
    the original fixture exactly."""
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
    if bo3_width_start is not None:
        r.sites[0].params["shift_fwhm_ppm"] = Param(bo3_width_start, min=bo3_width_start)
    if fix_bo4_pos:
        r.sites[1].params["isotropic_chemical_shift_ppm"].vary = False
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
    # 2σ is dmfit's sCZ_CQ = the Czjzek-paper σ_Cz (not the mode of |C_Q|)
    assert cols["sCZ_CQ=2σ (MHz)"] == pytest.approx((3.2, 0.1))
    # √⟨P_Q²⟩ = √5·σ_Cz = 2√5·σ in the stored (mrsimulator) σ
    v, e = cols["√⟨P_Q²⟩ (MHz)"]
    assert v == pytest.approx(2.0 * np.sqrt(5) * 1.6)
    assert e == pytest.approx(2.0 * np.sqrt(5) * 0.05)


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


def test_run_batch_writes_acquisition_outputs(tmp_path):
    """N4: fits carrying an acquisition block get acquisition.csv / .tex and
    an '## Acquisition' section with ranges and the varying columns; CSV
    sources (no block) write neither and the checkbox-off path writes
    neither either."""
    from pathlib import Path

    from test_acquisition import _block_3102_like

    paths = [_write_fit(tmp_path, "glassA", 15.0), _write_fit(tmp_path, "glassB", 14.5)]
    out = tmp_path / "plain"
    batch.run_batch(paths, out, make_plots=False, formats=("csv", "latex", "markdown"))
    assert not (out / "acquisition.csv").exists() and not (out / "acquisition.tex").exists()
    md = (out / "report.md").read_text(encoding="utf-8")
    assert "## Acquisition" not in md and "no acquisition record for glassA" in md

    for p, d1 in zip(paths, (12.5, 36.0)):
        d = json.loads(Path(p).read_text(encoding="utf-8"))
        d["acquisition"] = _block_3102_like(
            expno_path=f"X:/DATA/2026-01/{d['sample']}/24", sample=d["sample"], nucleus="11B",
            d1_s=d1, spin_rate_Hz=35714.0, mas_uncertain=False)
        d["spin_rate_Hz"], d["mas_uncertain"] = 35714.0, False
        Path(p).write_text(json.dumps(d), encoding="utf-8")
    out2 = tmp_path / "acq"
    res = batch.run_batch(paths, out2, make_plots=False, formats=("csv", "latex", "markdown"))
    assert (out2 / "acquisition.csv").exists() and (out2 / "acquisition.tex").exists()
    assert str(out2 / "acquisition.csv") in res.files and str(out2 / "acquisition.tex") in res.files
    md = (out2 / "report.md").read_text(encoding="utf-8")
    assert "## Acquisition" in md and "12.5–36 s" in md and "Table S1" in md
    assert "Parameters that vary across the series: D1 (s)." in md
    assert "**D1 (s)**" in md
    tex = (out2 / "acquisition.tex").read_text(encoding="utf-8")
    assert max(ord(c) for c in tex) < 0x80 and r"\textbf{12.5}" in tex
    csv_txt = (out2 / "acquisition.csv").read_text(encoding="utf-8")
    assert csv_txt.strip().splitlines()[-1].startswith("# varies")
    out3 = tmp_path / "off"
    res3 = batch.run_batch(paths, out3, make_plots=False, acquisition_table=False)
    assert not (out3 / "acquisition.csv").exists()
    assert "## Acquisition" not in (out3 / "report.md").read_text(encoding="utf-8")
    assert not any("acquisition" in f for f in res3.files)


def test_report_marks_fixed_and_at_bound_cells(tmp_path):
    """N5: every table the batch report writes marks held (†) and at-bound
    (‡) values -- flag columns in table.csv (numbers unmarked), LaTeX
    superscripts plus one legend row in table.tex, glyphs plus the legend
    sentence and a ## Constraints section in report.md."""
    paths = [_write_fit(tmp_path, "glassA", 15.0, bo3_width_start=8.0, fix_bo4_pos=True),
             _write_fit(tmp_path, "glassB", 14.5, bo3_width_start=8.0, fix_bo4_pos=True)]
    out = tmp_path / "report"
    res = batch.run_batch(paths, out, error_method="covariance",
                          make_plots=False, formats=("csv", "latex", "markdown"))
    assert res.n_fits == 2 and res.n_sites == 4

    with open(out / "table.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    bo3 = [r for r in rows if r["label"] == "BO3"]
    bo4 = [r for r in rows if r["label"] == "BO4"]
    assert len(bo3) == 2 and len(bo4) == 2
    assert all(r["dCS/FWHM (ppm) flag"] == "at_min" for r in bo3)
    assert all(r["δiso (ppm) flag"] == "fixed" for r in bo4)
    assert all(r["δiso (ppm) flag"] == "" for r in bo3)
    assert all(r["dCS/FWHM (ppm) flag"] == "" for r in bo4)
    assert all(r["pop (%) flag"] == "" for r in rows)
    for r in rows:                      # numbers stay numbers: no glyph in a value
        for c in ("δiso (ppm)", "dCS/FWHM (ppm)", "pop (%)"):
            assert not any(g in r[c] for g in "†‡§")
            float(r[c])
    assert all(float(r["dCS/FWHM (ppm)"]) == pytest.approx(8.0, abs=0.01) for r in bo3)
    assert all(float(r["δiso (ppm)"]) == 0.0 for r in bo4)
    header = list(rows[0].keys())
    assert header.index("δiso (ppm) err") < header.index("δiso (ppm) flag")

    tex = (out / "table.tex").read_text(encoding="utf-8")
    assert r"$^{\ddagger}$" in tex and r"$^{\dagger}$" in tex and r"\S" not in tex
    tl = tex.splitlines()
    foot = [ln for ln in tl if ln.startswith(r"\multicolumn")]
    assert len(foot) == 1
    assert tl.index(foot[0]) < tl.index(r"\end{tabular}")
    assert r"$\dagger$ held fixed" in foot[0] and r"$\ddagger$ finished at a bound" in foot[0]
    assert "linked" not in foot[0]
    assert "‡" not in tex and "†" not in tex        # Unicode glyphs mapped to LaTeX

    md = (out / "report.md").read_text(encoding="utf-8")
    table_lines = [ln for ln in md.splitlines() if ln.startswith("| glass")]
    assert len(table_lines) == 4
    assert all("‡" in ln for ln in table_lines if "| BO3 |" in ln)
    assert all("†" in ln for ln in table_lines if "| BO4 |" in ln)
    para = next(ln for ln in md.splitlines() if ln.startswith("Errors are ± one standard error"))
    assert "† held fixed · ‡ finished at a bound" in para and "§" not in para
    assert "## Constraints" in md
    cons = md.split("## Constraints", 1)[1].split("## ", 1)[0]
    for sample in ("glassA", "glassB"):
        line = next(ln for ln in cons.splitlines() if ln.startswith(f"- **{sample}**"))
        assert "dCS/FWHM (ppm) (BO3) lower bound 8" in line
        assert "δiso (ppm) (BO4)" in line and "gl (BO3)" in line
        assert line.endswith("linked: —")


def test_report_has_no_markers_for_an_unconstrained_batch(tmp_path):
    """Free values print exactly as before: no glyph, no legend, no
    Constraints section beyond what the recipe actually holds (gl is fixed
    in the fixture, so the section lists just that)."""
    paths = [_write_fit(tmp_path, "freeA", 15.0)]
    out = tmp_path / "report"
    batch.run_batch(paths, out, error_method="covariance",
                    make_plots=False, formats=("csv", "latex", "markdown"))
    with open(out / "table.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert all(r[k] == "" for r in rows for k in r if k.endswith(" flag"))
    tex = (out / "table.tex").read_text(encoding="utf-8")
    assert r"\multicolumn" not in tex and "dagger" not in tex
    md = (out / "report.md").read_text(encoding="utf-8")
    assert "‡" not in md and "Marked values" not in md
    cons = md.split("## Constraints", 1)[1]
    assert "fixed: gl (BO3), gl (BO4) · at a bound: — · linked: —" in cons


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


def _make_larproj(tmp_path, samples: list[tuple[str, float]]) -> str:
    """A minimal stand-in for app.py's save_project() output: workspaces
    with an embedded recipe + exp_ppm/exp_amp, no data on disk needed."""
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
    proj = tmp_path / "session.larproj.json"
    proj.write_text(json.dumps({"larmor_project_version": 1, "active": 0,
                                "workspaces": workspaces}), encoding="utf-8")
    return str(proj)


def test_load_entries_expands_a_project_bundle_into_one_entry_per_workspace(tmp_path):
    """A .larproj.json (app.py's save_project, multiple spectra in one file)
    is not itself a single Recipe -- load_any's plain .json branch can't read
    it (previously: silently skipped with "could not load", despite the
    batch-report help text claiming .larproj was a supported input)."""
    proj = _make_larproj(tmp_path, [("glassA", 15.0), ("glassB", 20.0)])
    entries, warnings = batch.load_entries([proj])
    assert not warnings
    assert len(entries) == 2
    assert {e.sample for e in entries} == {"glassA", "glassB"}
    for e in entries:
        assert e.nucleus == "11B" and e.ppm.size == 200 and e.amp.size == 200


def test_load_entries_skips_an_empty_project_with_a_warning(tmp_path):
    proj = tmp_path / "empty.larproj.json"
    proj.write_text(json.dumps({"larmor_project_version": 1, "workspaces": []}),
                    encoding="utf-8")
    entries, warnings = batch.load_entries([str(proj)])
    assert not entries
    assert any("no spectra in this project" in w for w in warnings)
