"""larmor.series_table -- Qt-free identity, order, joined columns, replicate
statistics, OLS and exports of a batch / sequential series.

The folder-name STRINGS are real ones from the NMRFAM tree (2026-01, 2026-03b,
2026-05) but nothing here touches the disk except the one require()-gated
check at the end."""
import csv
import json

import numpy as np
import pytest

from larmor import series_table as st
from larmor.series_table import SeriesColumn, SeriesTable

from conftest import LAW_CA_11B, require

_D = "C:/data/NMRFAM/DATA"
TRIPLE = [f"{_D}/2026-03b/03232026_P5-Bi8-12_SS_ALP/3104/pdata/1/1r",
          f"{_D}/2026-05/04272026_P5-Bi8-12_SS_ALP/3114/pdata/1/1r",
          f"{_D}/2026-05/05082026_P5-Bi8-12_SS_ALP/3102"]
TREE_ORDER = [f"{_D}/2026-01/01192026_SR31649_Base0Ca_SS_ALP/24",
              f"{_D}/2026-01/01202026_SR31648_Base3Ca_SS_ALP/24",
              f"{_D}/2026-01/01202026_SR31649_Base1Ca_SS_ALP/24",
              f"{_D}/2026-01/01202026_SR31649_Base2Ca_SS_ALP/24",
              f"{_D}/2026-01/01202026_SR31649_Base4Ca_SS_ALP/24"]


def _tree_table():
    return SeriesTable.from_paths(TREE_ORDER)


# ---------------------------------------------------------------- identity
def test_from_paths_names_come_from_scan_and_replicates_share_a_group():
    """Names are io/scan's (folder key + date tag when a glass was measured
    twice); the group is the base label, so the P5-Bi8-12 triple averages."""
    t = SeriesTable.from_paths(TRIPLE + ["C:/x/batch00.csv"])
    assert t.labels() == ["P5-Bi8-12 (03232026)", "P5-Bi8-12 (04272026)",
                          "P5-Bi8-12 (05082026)", "batch00"]
    assert [r.group for r in t.rows] == ["P5-Bi8-12"] * 3 + ["batch00"]
    assert t.groups() == {"P5-Bi8-12": [0, 1, 2], "batch00": [3]}
    assert t.has_replicates()
    assert [r.folder for r in t.rows] == ["03232026_P5-Bi8-12_SS_ALP",
                                          "04272026_P5-Bi8-12_SS_ALP",
                                          "05082026_P5-Bi8-12_SS_ALP", ""]
    assert t.paths() == TRIPLE + ["C:/x/batch00.csv"]
    # the dialogs' per-spectrum dicts feed the same constructor
    t2 = SeriesTable.from_spectra([
        {"path": "p0", "sample": "Base0Ca", "group": "Base0Ca",
         "folder": "01192026_SR31649_Base0Ca_SS_ALP", "title": "11B with short tip angle"},
        {"path": "p1", "sample": "Base1Ca"}])
    assert t2.labels() == ["Base0Ca", "Base1Ca"]
    assert t2.rows[0].title == "11B with short tip angle"
    assert t2.rows[0].folder == "01192026_SR31649_Base0Ca_SS_ALP"
    assert t2.rows[1].group == "Base1Ca" and t2.rows[1].folder == ""


def test_rename_keeps_names_unique_and_moves_a_singleton_group():
    t = SeriesTable.from_paths(TRIPLE + [f"{_D}/2026-03b/03232026_P5-Bi0_SS_ALP/3104"])
    assert t.labels()[3] == "P5-Bi0"
    # a hand-typed collision gets the (2) suffix; the name applied is returned
    assert t.rename(3, "P5-Bi8-12 (03232026)") == "P5-Bi8-12 (03232026) (2)"
    assert t.rows[3].group == "P5-Bi8-12 (03232026) (2)"   # was its own group
    # a replicate's group (!= its name) survives a rename
    assert t.rename(1, "second run") == "second run"
    assert t.rows[1].group == "P5-Bi8-12"
    t.set_group(3, "P5-Bi8-12")
    assert t.groups() == {"P5-Bi8-12": [0, 1, 2, 3]}
    t.set_group(3, "")                                       # blank -> its own name
    assert t.rows[3].group == t.rows[3].display_name
    assert st.unique_names(["a", "a", "b", "a", ""]) == ["a", "a (2)", "b", "a (3)",
                                                        "spectrum 5"]


def test_sort_perm_is_natural_stable_puts_missing_last_and_move_keeps_rows_aligned():
    t = _tree_table()
    assert t.labels() == ["Base0Ca", "Base3Ca", "Base1Ca", "Base2Ca", "Base4Ca"]
    assert t.sort_perm(None) == [0, 2, 3, 1, 4]
    assert st.natural_key("10Ca") > st.natural_key("9Ca")
    assert sorted(["Base10Ca", "Base9Ca", "Base1Ca"], key=st.natural_key) == \
        ["Base1Ca", "Base9Ca", "Base10Ca"]
    t.add_column(SeriesColumn("CaO", "CaO (mol%)"), {0: 0, 1: 3, 2: None, 3: 2, 4: 4})
    assert t.sort_perm("CaO") == [0, 3, 1, 4, 2]              # NaN last
    assert t.sort_perm("CaO", descending=True) == [4, 1, 3, 0, 2]
    x, xe = t.x_values("CaO")
    assert x[:2].tolist() == [0.0, 3.0] and np.isnan(x[2]) and np.isnan(xe).all()
    t.move(4, 0)
    assert t.labels()[0] == "Base4Ca"
    assert t.rows[0].source_path == TREE_ORDER[4] and t.rows[0].values["CaO"] == 4.0
    u = _tree_table()
    u.add_column(SeriesColumn("CaO", "CaO (mol%)"), {0: 0, 1: 3, 2: None, 3: 2, 4: 4})
    perm = u.sort_perm(None)
    v = u.permuted(perm)
    u.apply_order(perm)
    assert u.labels() == v.labels() == ["Base0Ca", "Base1Ca", "Base2Ca", "Base3Ca", "Base4Ca"]
    assert [r.values["CaO"] for r in v.rows] == [0.0, None, 2.0, 3.0, 4.0]
    u.remove_column("CaO")
    assert u.columns == [] and "CaO" not in u.rows[0].values


# ---------------------------------------------------------------- CSV join
_EPMA_HEADER = ("sample,series,source,Bi_nom,P_nom,n_matrix,structure,SiO2_mol,"
                "SiO2_mol_sd,P2O5_mol,P2O5_mol_sd,Bi2O3_mol,Bi2O3_mol_sd,method,Vm,Vm_u")
_EPMA_ROWS = [
    ("P5Bi8-12", "P5", "-12", "8", "5", "15", "homogeneous", "48.1", "0.8", "5.2",
     "0.1", "7.9", "0.3", "plain means", "27.3", "0.09"),
    ("P5Bi0", "P5", "parent", "0", "5", "10", "homogeneous", "52.0", "0.9", "5.0",
     "0.1", "0.0", "0.0", "plain means", "26.7", "0.06"),
    ("P0Bi0", "P0", "parent", "0", "0", "10", "homogeneous", "50.6", "0.8", "0.0",
     "0.0", "0.0", "0.0", "shade-reg (p=0.000)", "26.7", "0.06"),
]


def _write_epma(path, delim=","):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=delim)
        w.writerow(_EPMA_HEADER.split(","))
        for r in _EPMA_ROWS:
            w.writerow(r)
    return str(path)


def test_read_wide_csv_sniffs_delimiter_bom_and_sample_column(tmp_path):
    outs = []
    for name, delim in (("a.csv", ","), ("b.csv", ";"), ("c.tsv", "\t")):
        headers, rows, sample = st.read_wide_csv(_write_epma(tmp_path / name, delim))
        assert headers[0] == "sample" and sample == "sample"       # BOM stripped
        outs.append(rows)
    assert outs[0] == outs[1] == outs[2]
    assert [r["sample"] for r in outs[0]] == ["P5Bi8-12", "P5Bi0", "P0Bi0"]
    nums = st.numeric_headers(headers, outs[0])
    assert "method" not in nums and "structure" not in nums and "source" not in nums
    assert {"Bi_nom", "SiO2_mol", "P2O5_mol_sd", "Vm", "Vm_u"} <= set(nums)
    pairs = st.pair_error_columns(headers)
    assert pairs == {"SiO2_mol": "SiO2_mol_sd", "P2O5_mol": "P2O5_mol_sd",
                     "Bi2O3_mol": "Bi2O3_mol_sd", "Vm": "Vm_u"}
    # other sample headers and the first-non-numeric fallback
    (tmp_path / "s.csv").write_text("Scope ;x\ng0;1\ng1;2\n", encoding="utf-8")
    assert st.read_wide_csv(tmp_path / "s.csv")[2] == "Scope"
    (tmp_path / "t.csv").write_text("x,glassname,y\n1,a,2\n3,b,4\n", encoding="utf-8")
    assert st.read_wide_csv(tmp_path / "t.csv")[2] == "glassname"
    assert st.pair_error_columns(["A population_pct", "A population_pct ±", "B amplitude"]) \
        == {"A population_pct": "A population_pct ±"}


def test_propose_mapping_exact_then_group_then_normalised_and_flags_ambiguity():
    t = SeriesTable.from_paths(["p0", "p1", "p2", "p3"],
                               names=["P5-Bi8-12", "P5-Bi0", "Base0Ca", "2Ca12F"])
    rows = [{"sample": s} for s in ("P5Bi8-12", "P5Bi0", "P0Bi0", "2Ca-12F", "P5Bi8-13")]
    m = st.propose_mapping(t, rows, "sample")
    assert [x.how for x in m] == ["normalised", "normalised", "none", "normalised"]
    assert [x.csv_row for x in m] == [0, 1, None, 3]
    # a replicate (tagged name) matches through its group, exactly
    rep = SeriesTable.from_paths(TRIPLE)
    rows2 = [{"sample": "P5-Bi8-12"}, {"sample": "P0Bi0"}]
    m2 = st.propose_mapping(rep, rows2, "sample")
    assert [(x.how, x.csv_row) for x in m2] == [("exact", 0)] * 3
    # exact wins over normalised when both exist
    m3 = st.propose_mapping(t, [{"sample": "P5Bi8-12"}, {"sample": "P5-Bi8-12"}], "sample")
    assert (m3[0].how, m3[0].csv_row) == ("exact", 1)
    # two normalised candidates and no exact one -> ambiguous, no default
    m4 = st.propose_mapping(t, [{"sample": "P5Bi8-12"}, {"sample": "p5 bi8 12"}], "sample")
    assert m4[0].how == "ambiguous" and m4[0].csv_row is None and m4[0].candidates == [0, 1]
    # a batch CSV written before the folder-derived names: its scopes are raw
    # sample folders -> read through io/scan's tokens (how = 'folder')
    old = [{"scope": "04272026_P5-Bi8-12_SS_ALP"}, {"scope": "03232026_P5-Bi0_SS_ALP"},
           {"scope": "shared"}]
    m5 = st.propose_mapping(t, old, "scope")
    assert [(x.how, x.csv_row) for x in m5[:2]] == [("folder", 0), ("folder", 1)]
    assert m5[2].how == "none" and m5[3].how == "none"
    m6 = st.propose_mapping(rep, old, "scope")           # the triple, through its group
    assert [(x.how, x.csv_row) for x in m6] == [("folder", 0)] * 3


def test_apply_join_writes_numbers_only_and_tags_labels(tmp_path):
    headers, rows, sample = st.read_wide_csv(_write_epma(tmp_path / "e.csv"))
    t = SeriesTable.from_paths(["p0", "p1", "p2"], names=["P5-Bi8-12", "P5-Bi0", "Base0Ca"])
    m = st.propose_mapping(t, rows, sample)
    cols = st.apply_join(t, rows, sample, m, ["P2O5_mol", "Bi2O3_mol", "method"],
                         tag="analysed", source="e.csv")
    assert [c.key for c in cols] == ["P2O5_mol", "Bi2O3_mol", "method"]
    assert cols[0] == SeriesColumn("P2O5_mol", "analysed P2O5_mol", "P2O5_mol_sd",
                                   "analysed", "e.csv")
    x, xe = t.x_values("P2O5_mol")
    assert x[:2].tolist() == [5.2, 5.0] and np.isnan(x[2])
    assert xe[:2].tolist() == [0.1, 0.1] and np.isnan(xe[2])
    assert t.rows[0].values["method"] is None            # text -> None, never a string
    assert t.rows[2].values["Bi2O3_mol"] is None         # unmatched row stays empty
    # re-joining the same header replaces the column instead of duplicating it
    st.apply_join(t, rows, sample, m, ["P2O5_mol"], tag="nominal")
    assert [c.label for c in t.columns if c.key == "P2O5_mol"] == ["nominal P2O5_mol"]
    assert len(t.columns) == 3


def _batch_result():
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor.batchfit import BatchFitResult
    recs = []
    for k, amp in enumerate((100.0, 80.0, 120.0)):
        a = Param(amp); a.stderr = 1.0 + k
        recs.append(Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=f"g{k}",
                           sites=[SiteModel(model="gauss_lor", label="A", params={
                               "isotropic_chemical_shift_ppm": Param(15.0),
                               "shift_fwhm_ppm": Param(6.0), "amplitude": a,
                               "gl": Param(1.0, vary=False)}),
                                  SiteModel(model="gauss_lor", label="B", params={
                               "isotropic_chemical_shift_ppm": Param(2.0),
                               "shift_fwhm_ppm": Param(3.0), "amplitude": Param(50.0),
                               "gl": Param(1.0, vary=False)})]))
    return BatchFitResult(recipes=recs, labels=[f"g{k}" for k in range(3)],
                          rmsd=[0.0] * 3, per_dataset=[], shared=("shift_fwhm_ppm",),
                          released=())


def test_read_batch_csv_wide_pivots_a_larmor_long_csv_per_scope(tmp_path):
    from larmor import batchfit
    res = _batch_result()
    path = tmp_path / "batch_table.csv"
    batchfit.write_shared_csv(res, path)
    headers, _rows, _s = st.read_wide_csv(path)
    assert st.is_long_batch_csv(headers)
    assert not st.is_long_batch_csv(_EPMA_HEADER.split(","))
    headers, rows, sample = st.read_batch_csv_wide(path)
    assert sample == "scope" and [r["scope"] for r in rows] == ["g0", "g1", "g2"]
    assert "A amplitude" in headers and "A amplitude ±" in headers
    assert [float(r["A amplitude"]) for r in rows] == [100.0, 80.0, 120.0]
    assert [float(r["A amplitude ±"]) for r in rows] == [1.0, 2.0, 3.0]
    assert "A population_pct" in headers and "B population_pct" in headers
    assert st.pair_error_columns(headers)["A amplitude"] == "A amplitude ±"
    nums = st.numeric_headers(headers, rows)
    assert "A amplitude" in nums and "scope" not in nums
    # joins like a composition table: the 31P bonded-P column of another batch
    t = SeriesTable.from_paths(["p0", "p1", "p2"], names=["g2", "g0", "gX"])
    m = st.propose_mapping(t, rows, sample)
    st.apply_join(t, rows, sample, m, ["A amplitude"], source="batch_table.csv")
    x, xe = t.x_values("A amplitude")
    assert x[:2].tolist() == [120.0, 100.0] and np.isnan(x[2])
    assert xe[:2].tolist() == [3.0, 1.0]


# ---------------------------------------------------------------- statistics
def test_replicate_stats_mean_spread_and_n():
    groups = {"a": [0, 1, 2], "b": [3]}
    r = st.replicate_stats(groups, [1, 2, 3, 5], yerr=[.1, .1, .1, .4],
                           x=[10, 10, 10, 20], xerr=[.5, .5, .5, .2])
    assert r["labels"] == ["a", "b"]
    assert r["y"].tolist() == [2.0, 5.0]
    assert r["yerr"].tolist() == pytest.approx([1.0, 0.4])       # ddof=1; singleton keeps its own
    assert r["n"].tolist() == [3, 1]
    assert r["x"].tolist() == [10.0, 20.0]
    assert r["xerr"].tolist() == pytest.approx([0.5, 0.2])       # identical x -> RMS of ±
    r2 = st.replicate_stats({"a": [0, 1, 2]}, [1, 2, 3], x=[9, 10, 11])
    assert r2["xerr"].tolist() == pytest.approx([1.0])           # differing x -> sample std
    assert r2["yerr"][0] == pytest.approx(1.0)
    r3 = st.replicate_stats({"a": [0, 1]}, [np.nan, 4.0], yerr=[np.nan, np.nan])
    assert r3["y"].tolist() == [4.0] and np.isnan(r3["yerr"][0]) and r3["n"].tolist() == [1]
    assert r3["x"] is None and r3["xerr"] is None


def test_ols_matches_polyfit_and_pearson_r():
    x = np.array([10.0, 20.0, 30.0, 40.0, np.nan, 50.0])
    y = -0.092 * x + 65.0
    y[4] = 12.0                                             # pairs with the NaN x: skipped
    fit = st.ols(x, y)
    assert fit["slope"] == pytest.approx(-0.092) and fit["intercept"] == pytest.approx(65.0)
    assert fit["r"] == pytest.approx(-1.0) and fit["n"] == 5
    assert np.isfinite(fit["slope_err"]) and np.isfinite(fit["intercept_err"])
    two = st.ols([1, 2], [3, 5])
    assert two["slope"] == pytest.approx(2.0) and np.isnan(two["slope_err"]) and two["n"] == 2
    one = st.ols([1], [3])
    assert one["n"] == 1 and all(np.isnan(one[k]) for k in
                                 ("slope", "intercept", "slope_err", "intercept_err", "r"))
    none = st.ols([], [])
    assert none["n"] == 0 and np.isnan(none["slope"])
    lab = st.ols_label(fit)
    assert lab.startswith("slope -0.092 ±") and "r = -1.000" in lab and "n = 5" in lab
    assert "± —" in st.ols_label(two)


# ---------------------------------------------------------------- exports
def _joined_triple():
    t = SeriesTable.from_paths(TRIPLE + [f"{_D}/2026-03b/03232026_P5-Bi0_SS_ALP/3104"])
    t.add_column(SeriesColumn("P2O5_mol", "analysed P2O5_mol", "P2O5_mol_sd", "analysed"),
                 {0: 5.2, 1: 5.2, 2: 5.2, 3: 5.0})
    for r in t.rows:
        r.values["P2O5_mol_sd"] = 0.1
    t.add_column(SeriesColumn("Bi2O3 (mol%)", "Bi2O3 (mol%)"), {0: 7.9, 1: 7.9, 2: 7.9, 3: 0.0})
    t.add_column(SeriesColumn("Vm", "Vm", "Vm_u"), {0: 27.3, 1: 27.3, 2: 27.3, 3: 26.7})
    t.add_column(SeriesColumn("n_matrix", "n_matrix"), {0: 15, 1: 15, 2: 15, 3: 10})
    return t


def test_write_series_csv_is_wide_with_every_column_and_extras(tmp_path):
    t = _joined_triple()
    path = tmp_path / "batch_table_series.csv"
    st.write_series_csv(t, path, extra={"RMSD": [0.01, 0.02, 0.03, 0.04]})
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["position", "name", "group", "folder", "title", "source_path",
                       "analysed P2O5_mol", "analysed P2O5_mol ±", "Bi2O3 (mol%)",
                       "Vm", "Vm ±", "n_matrix", "RMSD"]
    assert rows[1][:4] == ["1", "P5-Bi8-12 (03232026)", "P5-Bi8-12",
                           "03232026_P5-Bi8-12_SS_ALP"]
    assert rows[1][5] == TRIPLE[0] and rows[1][6:8] == ["5.2", "0.1"]
    assert rows[4][1] == "P5-Bi0" and rows[4][-1] == "0.04" and rows[4][10] == ""


def test_dust_rows_maps_oxides_to_canonical_names_and_appends_n4():
    assert st.canonical_oxide("na2o_wt%") == "Na2O"
    assert st.canonical_oxide("Bi2O3 (mol%)") == "Bi2O3"
    assert st.canonical_oxide("P2O5_mol") == "P2O5"
    assert st.canonical_oxide("SiO2 mol %") == "SiO2"
    assert st.canonical_oxide("Na₂O") == "Na2O"
    assert st.canonical_oxide("R'") is None and st.canonical_oxide("Vm") is None
    assert st.canonical_oxide("n_matrix") is None and st.canonical_oxide("Bi_nom") is None
    t = _joined_triple()
    n4 = np.array([0.40, 0.42, 0.44, 0.30])
    n4e = np.array([0.01, 0.01, 0.01, 0.02])
    header, rows, skipped = st.dust_rows(t, n4, n4e, average=False)
    assert header == ["Sample", "P2O5", "Bi2O3", "N4_measured", "N4_measured_err"]
    assert skipped == ["Vm", "n_matrix"]
    assert len(rows) == 4 and rows[0][0] == "P5-Bi8-12 (03232026)"
    assert rows[0][1:] == ["5.2", "7.9", "0.4", "0.01"]
    header, rows, _ = st.dust_rows(t, n4, n4e, average=True)
    assert [r[0] for r in rows] == ["P5-Bi8-12", "P5-Bi0"]
    assert float(rows[0][3]) == pytest.approx(0.42)
    assert float(rows[0][4]) == pytest.approx(np.std([0.40, 0.42, 0.44], ddof=1))
    assert rows[1][3:] == ["0.3", "0.02"]
    # no series table: Sample + N4 only, labelled by the result's names
    header, rows, skipped = st.dust_rows(None, [0.5, 0.6], None, names=["a", "b"])
    assert header == ["Sample", "N4_measured", "N4_measured_err"] and skipped == []
    assert rows == [["a", "0.5", ""], ["b", "0.6", ""]]


def test_species_bar_and_correlation_specs_render():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from larmor import figures
    spec = st.species_bar_spec(["a", "b"], ["s0 A", "s1 B"],
                               np.array([[60.0, 40.0], [np.nan, 100.0]]),
                               xlabel="analysed P2O5 (mol%)", colors=["#123456", None])
    assert spec["kind"] == "species_bar" and spec["categories"] == ["a", "b"]
    assert spec["series"][0] == {"label": "s0 A", "values": [60.0, 0.0], "color": "#123456"}
    assert spec["xlabel"] == "analysed P2O5 (mol%)"
    fig = figures.render(spec)
    assert fig.axes[0].get_xlabel() == "analysed P2O5 (mol%)"
    plt.close(fig)
    x = np.array([1.0, 2.0, 3.0, 4.0]); y = 2.0 * x + 1.0
    traces = st.correlation_traces(x, y, [.1] * 4, [.2] * 4, st.ols(x, y), "s0", "#ff0000")
    assert len(traces) == 2
    assert traces[0]["marker"] == "o" and traces[0]["linestyle"] == "none"
    assert traces[0]["data"]["xerr"] == [.1] * 4 and traces[0]["data"]["yerr"] == [.2] * 4
    line = traces[1]
    assert line["linestyle"] == "--" and line["data"]["x"] == [1.0, 4.0]
    assert line["data"]["y"] == pytest.approx([3.0, 9.0])
    assert "r = 1.000" in line["label"] and "n = 4" in line["label"] and line["label"].startswith("s0 ")
    fig = figures.render({"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
                          "traces": traces, "xlabel": "x"})
    assert len(fig.axes[0].lines) >= 2
    plt.close(fig)


# ---------------------------------------------------------------- persistence
def test_series_table_roundtrip_dict_file_aligned_to_and_merge(tmp_path):
    t = _joined_triple()
    t.rename(3, "parent glass")
    d = json.loads(json.dumps(t.to_dict()))
    back = SeriesTable.from_dict(d)
    assert back == t
    p = tmp_path / "series.series.json"
    t.save(p)
    assert SeriesTable.load(p) == t
    # aligned_to: subset + reorder by source_path, a fresh row for a new path
    fresh = SeriesTable.from_paths([TRIPLE[2], TRIPLE[0], "C:/x/new.csv"],
                                   names=["P5-Bi8-12 (05082026)", "P5-Bi8-12 (03232026)", "new"])
    a = t.aligned_to([TRIPLE[2], TRIPLE[0], "C:/x/new.csv"], fresh=fresh)
    assert a.labels() == ["P5-Bi8-12 (05082026)", "P5-Bi8-12 (03232026)", "new"]
    assert a.rows[2].values == {"P2O5_mol": None, "P2O5_mol_sd": None, "Bi2O3 (mol%)": None,
                                "Vm": None, "Vm_u": None, "n_matrix": None}
    assert a.rows[0].values["P2O5_mol"] == 5.2 and a.columns == t.columns
    b = t.aligned_to(["C:/x/other.csv"])                 # no fresh table: bare row
    assert b.labels() == ["other"] and b.rows[0].group == "other"
    # merge_from: a .series.json over the same spectra restores names / columns,
    # pairing by path, else by (normalised) name
    cur = SeriesTable.from_paths(TRIPLE + ["C:/y/P5Bi0.csv"])
    assert cur.labels()[3] == "P5Bi0"
    saved = SeriesTable.from_dict(t.to_dict())
    saved.rows[3].source_path = "D:/elsewhere/P5-Bi0/1"      # moved -> pairs by norm name
    saved.rows[3].display_name = "P5-Bi0"
    saved.rename(0, "first run")
    assert cur.merge_from(saved) == [0, 1, 2, 3]
    assert cur.labels()[0] == "first run" and cur.labels()[3] == "P5-Bi0"
    lone = SeriesTable.from_paths(["C:/z/unknown.csv", TRIPLE[1]])
    assert lone.merge_from(saved) == [None, 1]
    assert lone.labels() == ["unknown", "P5-Bi8-12 (04272026)"]
    assert lone.rows[0].values["P2O5_mol"] is None and lone.rows[1].values["P2O5_mol"] == 5.2
    assert st.as_float(" 1.5 ") == 1.5 and st.as_float("") is None and st.as_float("x") is None
    assert [c.key for c in cur.columns] == ["P2O5_mol", "Bi2O3 (mol%)", "Vm", "n_matrix"]
    assert cur.rows[3].values["P2O5_mol"] == 5.0


def test_from_paths_on_the_real_2026_01_series():
    require(LAW_CA_11B[0])
    t = SeriesTable.from_paths(LAW_CA_11B)
    assert t.labels() == ["Base0Ca", "Base1Ca", "Base2Ca", "Base3Ca", "Base4Ca"]
    assert all(r.folder.endswith("_SS_ALP") for r in t.rows)
    assert [r.group for r in t.rows] == t.labels()
