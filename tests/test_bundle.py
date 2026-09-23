"""Series publication bundle (larmor.io.bundle): per-spectrum curves with the
experiment exactly as fitted, manifest.csv, README.txt, recipe copies and the
tables -- Qt-free, on the synthetic 11B series of test_batchfit."""
import csv
import json
import re
from pathlib import Path

import numpy as np
import pytest

import larmor
from conftest import BRUKER_1R, require
from larmor import batchfit, methods, seqfit, series_grid
from larmor.io import bundle
from larmor.recipe import Param, Recipe, sha256_of
from test_batchfit import _data_for, _entries, _start


def _read_curves(path):
    """A curves CSV as a structured array (the '# key=value' lines dropped
    first: genfromtxt(names=True) would take '# LARMOR curves' as the header;
    pandas reads the file with comment='#')."""
    rows = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines()
            if ln and not ln.startswith("#")]
    return np.genfromtxt(rows, delimiter=",", names=True)


def _manifest(folder):
    with open(Path(folder) / "manifest.csv", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        return list(r.fieldnames), list(r)


def _header_lines(path):
    return [ln for ln in Path(path).read_text(encoding="utf-8").splitlines()
            if ln.startswith("#")]


def test_write_bundle_writes_every_file_and_the_manifest_contract(tmp_path):
    entries = _entries()
    res = batchfit.batch_fit(entries)
    out = bundle.write_bundle(res, _data_for(entries), tmp_path, kind="batch")

    names = {p.name for p in tmp_path.iterdir()}
    assert {"README.txt", "manifest.csv", "batch_table.csv"} <= names
    for s in ("01_g0", "02_g1", "03_g2"):
        assert f"{s}.recipe.json" in names and f"{s}_curves.csv" in names
    assert not list(tmp_path.glob("batch_fit_*.csv"))     # no error method computed
    assert set(out.files) == names

    header, rows = _manifest(tmp_path)
    assert header == list(bundle.MANIFEST_COLUMNS)
    assert [r["name"] for r in rows] == res.labels
    assert [r["stem"] for r in rows] == ["01_g0", "02_g1", "03_g2"]
    assert [r["index"] for r in rows] == ["1", "2", "3"]
    for r in rows:
        assert (tmp_path / r["recipe_file"]).is_file()
        assert (tmp_path / r["curves_file"]).is_file()
        assert r["processing"] == "[]" and r["n_sites"] == "2"
        assert r["error_method"] == "covariance" and r["note"] == ""
        assert r["larmor_version"] == larmor.__version__
        assert r["mrsimulator_version"]
        assert r["nucleus"] == "11B" and float(r["larmor_MHz"]) == 160.0
        assert (float(r["fit_window_hi_ppm"]), float(r["fit_window_lo_ppm"])) == (40.0, -10.0)
    assert out.kind == "batch" and out.warnings == []
    assert "manifest.csv" in out.summary and "README.txt" in out.summary
    assert out.summary.startswith("publication bundle: 3 spectra")


def test_curves_experiment_is_the_fitted_array_and_rmsd_recomputes(tmp_path):
    entries = _entries()
    res = batchfit.batch_fit(entries)
    bundle.write_bundle(res, _data_for(entries), tmp_path)
    _h, rows = _manifest(tmp_path)
    for k, (rec0, ppm, amp, win) in enumerate(entries):
        tab = _read_curves(tmp_path / rows[k]["curves_file"])
        assert list(tab.dtype.names) == ["ppm", "experiment", "model", "residual",
                                         "s0_A", "s1_B"]
        assert np.array_equal(tab["experiment"], amp)          # bit-exact
        assert np.array_equal(tab["ppm"], ppm)
        assert np.allclose(tab["model"], tab["s0_A"] + tab["s1_B"],
                           atol=1e-6 * tab["model"].max())
        assert np.allclose(tab["residual"], tab["experiment"] - tab["model"],
                           atol=1e-6 * tab["model"].max())
        rec = res.recipes[k]
        r, npts = bundle.curve_rmsd(tab["ppm"], tab["experiment"], tab["model"],
                                    rec.fit_window_ppm)
        assert r == pytest.approx(res.rmsd[k], rel=1e-6)
        assert r == pytest.approx(rec.fit_rmsd, rel=1e-6)
        assert float(rows[k]["rmsd"]) == pytest.approx(r, rel=1e-6)
        hi, lo = max(win), min(win)
        assert npts == int(((ppm >= lo) & (ppm <= hi)).sum())
        assert int(rows[k]["n_points_window"]) == npts
        # a header line per convention, so File > Open reopens the file
        hdr = _header_lines(tmp_path / rows[k]["curves_file"])
        assert hdr[0] == "# LARMOR curves"
        assert "# nucleus=11B" in hdr and "# fit_window_ppm=40.0,-10.0" in hdr
        assert f"# sample={rec.sample}" in hdr


def test_baseline_processing_writes_raw_and_experiment_matches_replay(tmp_path):
    from larmor.loader import apply_processing

    entries = _entries()
    raws = []
    # spectrum 0: a flat offset removed by the recipe's own processing step
    rec, ppm, amp, win = entries[0]
    raw0 = amp + 7.0
    rec.processing = [{"op": "flat_baseline"}]
    ppm_p, amp_p, _notes = apply_processing(rec, ppm, raw0)
    entries[0] = (rec, ppm_p, np.asarray(amp_p, float), win)
    raws.append(raw0)
    raws.append(entries[1][2].copy())      # spectrum 1: raw identical to amp
    raws.append(None)                      # spectrum 2: no raw at all
    res = batchfit.batch_fit(entries)
    assert res.recipes[0].processing == [{"op": "flat_baseline"}]

    bundle.write_bundle(res, _data_for(entries), tmp_path, raw=raws)
    _h, rows = _manifest(tmp_path)

    tab0 = _read_curves(tmp_path / rows[0]["curves_file"])
    assert list(tab0.dtype.names)[:4] == ["ppm", "experiment", "experiment_raw", "model"]
    assert np.array_equal(tab0["experiment_raw"], raw0)
    assert np.array_equal(tab0["experiment"], entries[0][2])     # as fitted, not re-derived
    assert not np.array_equal(tab0["experiment"], raw0)
    assert '# processing=[{"op": "flat_baseline"}]' in _header_lines(
        tmp_path / rows[0]["curves_file"])
    assert rows[0]["processing"] == json.dumps(res.recipes[0].processing)
    assert rows[0]["processing"].startswith('[{"op": "')
    for k in (1, 2):
        tab = _read_curves(tmp_path / rows[k]["curves_file"])
        assert "experiment_raw" not in tab.dtype.names
        assert rows[k]["processing"] == "[]"


def test_excluded_site_is_a_zero_column_and_flagged_and_names_are_unique(tmp_path):
    entries = _entries()[:2]
    for e in entries:
        e[0].sample = "glass"                       # identical labels
    entries[1][0].sites[1].params["amplitude"] = Param(0.0, vary=False, min=0.0, max=0.0)
    res = batchfit.batch_fit(entries)
    assert res.labels == ["glass", "glass"]
    bundle.write_bundle(res, _data_for(entries), tmp_path)
    _h, rows = _manifest(tmp_path)
    assert [r["stem"] for r in rows] == ["01_glass", "02_glass"]

    tab = _read_curves(tmp_path / "02_glass_curves.csv")
    assert list(tab.dtype.names) == ["ppm", "experiment", "model", "residual", "s0_A", "s1_B"]
    assert np.all(tab["s1_B"] == 0.0)
    assert "# excluded=s1_B" in _header_lines(tmp_path / "02_glass_curves.csv")
    assert not any(ln.startswith("# excluded") for ln in
                   _header_lines(tmp_path / "01_glass_curves.csv"))
    assert [r["excluded_sites"] for r in rows] == ["", "s1_B"]
    # the table still omits the excluded site (one amplitude row for s1, not two)
    with open(tmp_path / "batch_table.csv", newline="", encoding="utf-8") as f:
        trows = list(csv.DictReader(f))
    s1_amp = [r for r in trows if r["scope"] == "glass" and r["site"] == "s1"
              and r["param"] == "amplitude"]
    assert len(s1_amp) == 1

    assert len(set(bundle.stems(["g"] * 10))) == 10
    assert bundle.stems(["27Al failed, MAS stopped"]) == ["01_27Al_failed__MAS_stopped"]
    # the dialog / CLI slug rule: alphanumerics (str.isalnum, so 'é' stays),
    # '-' and '_' kept, everything else '_', 'fit' when empty
    assert bundle.stems(["a b", "c/d", "", "é.x"]) == ["01_a_b", "02_c_d", "03_fit", "04_é_x"]
    for s in bundle.stems(["a b", "c/d", "", "é.x"]):
        assert all(c.isalnum() or c in "-_" for c in s), s
    assert bundle.stems([f"s{k}" for k in range(120)])[0] == "001_s0"


def test_seq_wrap_recomputes_normalised_rmsd_and_fills_source_paths(tmp_path):
    entries = _entries()
    src = tmp_path / "src"
    src.mkdir()
    paths = []
    for k, (rec, ppm, amp, win) in enumerate(entries):
        p = src / f"g{k}.csv"
        p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n"
                     + "".join(f"{x!r},{y!r}\n" for x, y in zip(ppm, amp)),
                     encoding="utf-8")
        paths.append(str(p))
        assert not rec.source_path
    sres = seqfit.run_sequential(entries, passes=1)
    recs = list(sres.recipes)
    recs[2] = _start("g2")                          # never fitted (fit_rmsd None)
    wrap = batchfit.BatchFitResult(
        recipes=recs, labels=list(sres.labels), rmsd=list(sres.rmsd),
        per_dataset=[], shared=(), released=batchfit.all_but_amplitude(recs))

    out = bundle.write_bundle(wrap, _data_for(entries), tmp_path / "b", kind="seq",
                              source_paths=paths)
    b = tmp_path / "b"
    _h, rows = _manifest(b)
    for k in (0, 1):
        assert float(rows[k]["rmsd"]) == pytest.approx(recs[k].fit_rmsd, rel=1e-6)
        assert float(rows[k]["rmsd"]) != pytest.approx(sres.rmsd[k], rel=1e-3)   # not the raw one
        assert rows[k]["note"] == ""
    assert rows[2]["note"] == "not fitted" and (b / rows[2]["curves_file"]).is_file()
    assert rows[2]["rmsd"] != ""
    for k in range(3):
        rc = Recipe.load(b / rows[k]["recipe_file"])
        assert rc.source_path == paths[k]
        assert rc.source_sha256 == sha256_of(paths[k])
        assert rows[k]["source_path"] == paths[k]
        assert rows[k]["source_sha256"] == sha256_of(paths[k])
        assert rows[k]["source_kind"] == "csv"
    assert recs[0].source_path == ""                # the in-memory result untouched
    assert (b / "seq_table.csv").is_file() and not (b / "batch_table.csv").exists()
    readme = (b / "README.txt").read_text(encoding="utf-8")
    assert "sequential" in readme
    assert out.warnings == []
    # the table written from the copies carries the source paths
    with open(b / "seq_table.csv", newline="", encoding="utf-8") as f:
        srcs = {r["source_path"] for r in csv.DictReader(f) if r["scope"] != "shared"}
    assert srcs == set(paths)


def test_error_csv_only_when_a_method_was_computed(tmp_path):
    entries = _entries()
    res = batchfit.batch_fit(entries)
    data = _data_for(entries)
    bundle.write_bundle(res, data, tmp_path / "a", error_method="covariance")
    assert not list((tmp_path / "a").glob("batch_fit_*.csv"))

    batchfit.batch_error_analysis(res, data, method="covariance")
    out = bundle.write_bundle(res, data, tmp_path / "b", error_method="covariance")
    p = tmp_path / "b" / "batch_fit_covariance.csv"
    assert p.is_file() and "batch_fit_covariance.csv" in out.summary
    with open(p, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == batchfit.ERROR_HEADER
    assert len(rows) - 1 == len(batchfit.error_table(res, method="covariance"))
    _h, mrows = _manifest(tmp_path / "b")
    assert all(r["error_method"] == "covariance" for r in mrows)

    out = bundle.write_bundle(res, data, tmp_path / "c", error_method="montecarlo")
    assert not list((tmp_path / "c").glob("batch_fit_*.csv"))
    assert (tmp_path / "c" / "manifest.csv").is_file()


def test_source_facts_from_a_real_expno_and_missing_sources_stay_blank(tmp_path):
    blank = bundle.source_facts("C:/nowhere/at/all/1r")
    assert set(blank) == set(bundle.SOURCE_KEYS)
    assert all(v == "" for v in blank.values())
    assert bundle.source_facts("") == blank
    assert bundle.data_file("C:/nowhere/at/all/1r") == (None, "")

    p = tmp_path / "s.csv"
    p.write_text("# nucleus = 11B\n1,2\n2,3\n", encoding="utf-8")
    f = bundle.source_facts(str(p))
    assert f["source_kind"] == "csv" and f["data_file"] == str(p)
    assert f["source_sha256"] == sha256_of(p) and f["expno"] == "" and f["ns"] == ""
    # a recipe stands for its own source
    rec = Recipe(sample="s", source_path=str(p))
    rp = tmp_path / "s.recipe.json"
    rec.save(rp)
    assert bundle.data_file(str(rp)) == (p, "csv")

    require(BRUKER_1R)
    from larmor.io import bruker

    expno = BRUKER_1R.parents[2]
    before = bruker.snapshot(expno)
    f = bundle.source_facts(str(expno))               # the EXPNO dir, as load_any stores it
    assert f["source_kind"] == "bruker"
    assert f["expno"] == "2702" and f["procno"] == "1"
    assert f["ns"] == 56 and f["d1_s"] == 3.0
    assert f["pulprog"] == "zg"
    assert f["sf_MHz"] == pytest.approx(156.2816087, abs=1e-6)
    assert f["bf1_MHz"] == pytest.approx(156.281744)
    assert f["title"] == "27Al failed, MAS stopped"
    df = Path(f["data_file"])
    assert df.name == "1r" and df.parent.name == "1" and df.parent.parent.name == "pdata"
    assert f["source_sha256"] == sha256_of(BRUKER_1R)
    assert re.fullmatch(r"[0-9a-f]{64}", f["source_sha256"])
    assert f["date"]
    assert bundle.source_facts(str(BRUKER_1R)) == f
    assert bruker.snapshot(expno) == before           # read-only


def test_readme_has_methods_versions_and_glossary_and_studio_finds_the_recipes(tmp_path):
    entries = _entries()
    res = batchfit.batch_fit(entries)
    bundle.write_bundle(res, _data_for(entries), tmp_path)
    readme = (tmp_path / "README.txt").read_text(encoding="utf-8")
    assert methods.methods_sentence(res.recipes[0].to_dict(), "covariance") in readme
    assert f"LARMOR {larmor.__version__}" in readme
    assert "mrsimulator" in readme and "lmfit" in readme and "numpy" in readme
    assert "manifest.csv" in readme and "_curves.csv" in readme and ".recipe.json" in readme
    assert "RMSD = sqrt(mean((model - experiment)^2))" in readme
    assert "repr()" in readme and "9 significant digits" in readme
    assert "batch fit of 3 spectra" in readme
    assert "WARNINGS" not in readme

    csv_path = tmp_path / "batch_table.csv"
    found = series_grid.find_recipes_near_csv(csv_path)
    assert sorted(found) == sorted(str(tmp_path / f"{s}.recipe.json")
                                   for s in ("01_g0", "02_g1", "03_g2"))
    panels, _notes = series_grid.load_panels(str(csv_path))
    assert len(panels) == 3 and not any(p.reconstructed for p in panels)
