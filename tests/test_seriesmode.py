"""The Qt-free series mode: carry rules between neighbours, member statuses,
the series record and its tags, and the round trip into
seqfit.run_sequential and back."""
import copy
import json

import numpy as np
import pytest


def _site(pos, fwhm=6.0, amp=100.0, label="A", pos_bounds=(0, 30), expr=None):
    return {"model": "gauss_lor", "label": label, "params": {
        "isotropic_chemical_shift_ppm": {"value": pos, "stderr": 0.1, "vary": True,
                                         "min": pos_bounds[0], "max": pos_bounds[1],
                                         "expr": expr},
        "shift_fwhm_ppm": {"value": fwhm, "stderr": 0.2, "vary": True, "min": 0.1,
                           "max": None, "expr": None},
        "amplitude": {"value": amp, "stderr": None, "vary": True, "min": 0,
                      "max": None, "expr": None},
        "gl": {"value": 1.0, "stderr": None, "vary": False, "min": 0, "max": 1,
               "expr": None}}}


def _recipe(sites, **extra):
    d = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
         "sample": "x", "sites": sites}
    d.update(extra)
    return d


# ------------------------------------------------------------------ carry
def test_carry_into_empty_target_copies_the_model_and_rescales_amplitudes():
    from larmor.seriesmode import carry_into

    src = _recipe([_site(13.0, amp=100.0, label="A"), _site(20.0, amp=40.0, label="B")],
                  fit_window_ppm=[30.0, -10.0], fit_rmsd=0.01)
    src["sites"][1]["family"] = "BO3"
    src["sites"][1]["params"]["shift_fwhm_ppm"]["expr"] = "s0.shift_fwhm_ppm"
    dst = _recipe([], sample="target")
    new, note = carry_into(dst, src, ("isotropic_chemical_shift_ppm",), 41.0, 100.0,
                           src_name="LAW3Cl1Ca")
    assert dst["sites"] == []                              # inputs untouched
    assert [s["label"] for s in new["sites"]] == ["A", "B"]
    assert new["sites"][1]["family"] == "BO3"
    assert new["sites"][1]["params"]["shift_fwhm_ppm"]["expr"] == "s0.shift_fwhm_ppm"
    amps = [s["params"]["amplitude"]["value"] for s in new["sites"]]
    assert amps == pytest.approx([41.0, 16.4])              # ×0.41
    assert all(p["stderr"] is None for s in new["sites"] for p in s["params"].values())
    assert new["fit_window_ppm"] == [30.0, -10.0]
    assert "fit_rmsd" not in new and new["sample"] == "target"
    assert note == "2 lines copied from LAW3Cl1Ca, amplitudes scaled ×0.41"
    # unknown maxima: no scaling, and the note says nothing about it
    same, note2 = carry_into({}, src, None, 0.0, 100.0)
    assert [s["params"]["amplitude"]["value"] for s in same["sites"]] == [100.0, 40.0]
    assert note2 == "2 lines copied"


def test_carry_into_target_with_lines_seeds_only_carried_params_on_matching_indices():
    from larmor.seriesmode import carry_into

    # the source's two lines are both called "A": labels that identify
    # nothing, so the lines pair by INDEX (the rule unlabelled recipes keep)
    src = _recipe([_site(13.0, fwhm=6.0, amp=100.0), _site(20.0, fwhm=3.0, amp=50.0)])
    # the target has its OWN structure: four lines, the second one linked
    dst = _recipe([_site(12.0, fwhm=5.0, amp=80.0, label="own0"),
                   _site(19.0, fwhm=5.0, amp=10.0, label="own1",
                         expr="s0.isotropic_chemical_shift_ppm + 7"),
                   _site(25.0, fwhm=4.0, amp=5.0, label="own2"),
                   _site(28.0, fwhm=4.0, amp=5.0, label="own3")])
    new, note = carry_into(dst, src, ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm"),
                           80.0, 100.0, src_name="s0")
    assert len(new["sites"]) == 4 and [s["label"] for s in new["sites"]] == \
        ["own0", "own1", "own2", "own3"]
    p0 = new["sites"][0]["params"]
    assert p0["isotropic_chemical_shift_ppm"]["value"] == 13.0
    assert p0["shift_fwhm_ppm"]["value"] == 6.0
    assert p0["amplitude"]["value"] == 80.0                 # not carried
    assert p0["isotropic_chemical_shift_ppm"]["stderr"] is None   # a seeded value is unfitted
    p1 = new["sites"][1]["params"]
    assert p1["isotropic_chemical_shift_ppm"]["value"] == 19.0       # linked: untouched
    assert p1["isotropic_chemical_shift_ppm"]["expr"] == "s0.isotropic_chemical_shift_ppm + 7"
    assert p1["shift_fwhm_ppm"]["value"] == 3.0                      # free: seeded
    assert new["sites"][2]["params"]["isotropic_chemical_shift_ppm"]["value"] == 25.0
    assert new["sites"][3]["params"]["isotropic_chemical_shift_ppm"]["value"] == 28.0
    assert note == ("positions and widths of 2 matching lines seeded from s0; "
                    "this spectrum keeps its own 4 lines (own2, own3 have no "
                    "counterpart there)")
    assert dst["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] == 12.0


def test_carry_into_pairs_lines_by_label_and_names_the_ones_without_a_counterpart():
    """A series fitted one spectrum at a time rarely keeps one model: the
    carry pairs lines by their label, wherever they sit in the list, and a
    line without a namesake on the other side is left exactly as it is."""
    from larmor.seriesmode import carry_into

    src = _recipe([_site(13.0, fwhm=6.0, amp=100.0, label="A"),
                   _site(20.0, fwhm=3.0, amp=50.0, label="B")])
    dst = _recipe([_site(19.0, fwhm=5.0, amp=10.0, label="B"),
                   _site(25.0, fwhm=4.0, amp=5.0, label="C"),
                   _site(12.0, fwhm=5.0, amp=80.0, label="A")])
    new, note = carry_into(dst, src, ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm"),
                           1.0, 1.0, src_name="s0")
    assert [s["label"] for s in new["sites"]] == ["B", "C", "A"]        # own order kept
    b, c, a = (s["params"] for s in new["sites"])
    assert (b["isotropic_chemical_shift_ppm"]["value"], b["shift_fwhm_ppm"]["value"]) == (20.0, 3.0)
    assert (a["isotropic_chemical_shift_ppm"]["value"], a["shift_fwhm_ppm"]["value"]) == (13.0, 6.0)
    assert (c["isotropic_chemical_shift_ppm"]["value"], c["shift_fwhm_ppm"]["value"]) == (25.0, 4.0)
    assert [s["params"]["amplitude"]["value"] for s in new["sites"]] == [10.0, 5.0, 80.0]
    assert note == ("positions and widths of 2 matching lines seeded from s0; this "
                    "spectrum keeps its own 3 lines (C has no counterpart there)")
    # two labelled models with no name in common: nothing moves, and the note
    # says how to make them pair
    other = _recipe([_site(1.0, label="X"), _site(2.0, label="Y")])
    same, note2 = carry_into(other, src, None, 1.0, 1.0, src_name="s0")
    assert same == other
    assert note2.startswith("no line here shares a name with s0's (X, Y vs A, B) — nothing "
                            "seeded; this spectrum keeps its own 2 lines")
    assert "label column" in note2


def test_series_member_lock_round_trips_and_apply_sweep_result_skips_kept_members():
    from larmor.seriesmode import SeriesMember, SeriesSpec, apply_sweep_result

    spec = SeriesSpec(id="s", members=[SeriesMember(key="a", name="a"),
                                       SeriesMember(key="b", name="b", locked=True)])
    tags = [spec.tag_for(0), spec.tag_for(1)]
    assert tags[0]["locked"] is False and tags[1]["locked"] is True
    back = SeriesSpec.from_tags(tags)
    assert [m.locked for m in back.members] == [False, True]
    # a sweep result names the members it kept: they get None, so the
    # window leaves their fit, errors and verdict exactly as they were
    rec = _recipe([_site(10.0)], fit_rmsd=0.01)

    class R:
        recipes = [rec, rec]
        history = [{"rmsd": [0.1, 0.1]}]
        rmsd = [0.1, 0.1]
        per_dataset = []
        fixed = (1,)

    members = [{"recipe": rec, "name": "a"}, {"recipe": rec, "name": "b"}]
    out = apply_sweep_result(members, R())
    assert out[0] is not None and out[1] is None


def test_carry_into_clips_to_the_target_bounds_and_matches_seqfit_seed_from():
    from larmor import seqfit
    from larmor.recipe import Recipe
    from larmor.seriesmode import carry_into

    src = _recipe([_site(45.0, fwhm=0.05, amp=1.0)])
    dst = _recipe([_site(10.0, fwhm=5.0, amp=3.0, pos_bounds=(0, 30))])
    new, _ = carry_into(dst, src, None, 1.0, 1.0)
    p = new["sites"][0]["params"]
    assert p["isotropic_chemical_shift_ppm"]["value"] == 30.0     # clipped to max
    assert p["shift_fwhm_ppm"]["value"] == pytest.approx(0.1)     # clipped to min
    assert p["amplitude"]["value"] == 3.0                          # amplitude never carried by default
    # the same numbers larmor.seqfit.seed_from produces (the sweep's rule)
    d = Recipe.from_dict(dst)
    seqfit.seed_from(d, Recipe.from_dict(src), ("isotropic_chemical_shift_ppm",
                                                 "shift_fwhm_ppm", "gl"))
    ref = d.to_dict()["sites"][0]["params"]
    for pn in ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "gl", "amplitude"):
        assert p[pn]["value"] == pytest.approx(ref[pn]["value"])


def test_carry_into_edge_cases():
    from larmor.seriesmode import carry_into

    dst = _recipe([_site(10.0)])
    new, note = carry_into(dst, _recipe([]), None, 1.0, 1.0)
    assert new == dst and note.startswith("nothing to carry")
    # no parameter in common with the carried set: nothing moves, said plainly
    new, note = carry_into(dst, _recipe([_site(20.0)]), ("sigma_Cq_MHz",), 1.0, 1.0, "B")
    assert new["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] == 10.0
    assert note == "no matching parameters to seed from B; this spectrum keeps its own 1 line"


def test_carry_candidates_and_carry_for_follow_the_unticked_names():
    from larmor.seriesmode import SeriesSpec, carry_candidates

    recs = [_recipe([_site(10.0)]), _recipe([{"model": "czjzek", "label": "C", "params": {
        "isotropic_chemical_shift_ppm": {"value": 60}, "sigma_Cq_MHz": {"value": 2},
        "shift_fwhm_ppm": {"value": 8}, "amplitude": {"value": 1}}}])]
    assert carry_candidates(recs) == ("gl", "isotropic_chemical_shift_ppm",
                                      "shift_fwhm_ppm", "sigma_Cq_MHz")
    spec = SeriesSpec(id="s")
    assert spec.carry_for(recs) == carry_candidates(recs)
    spec.set_carry(["isotropic_chemical_shift_ppm", "sigma_Cq_MHz"], carry_candidates(recs))
    assert spec.options.carry_off == ("gl", "shift_fwhm_ppm")
    assert spec.carry_for(recs) == ("isotropic_chemical_shift_ppm", "sigma_Cq_MHz")
    # an unticked name that is not offered right now stays remembered
    spec.set_carry(["isotropic_chemical_shift_ppm"], ("isotropic_chemical_shift_ppm", "gl"))
    assert spec.options.carry_off == ("gl", "shift_fwhm_ppm")
    assert spec.carry_for([recs[0]]) == ("isotropic_chemical_shift_ppm",)


# ------------------------------------------------------------------ statuses
def test_member_status_and_rmsd_of():
    from larmor.recipe import Recipe
    from larmor.seriesmode import member_status, rmsd_of, signature_key

    assert member_status(None) == "unfitted"
    assert member_status(_recipe([])) == "unfitted"
    assert member_status(_recipe([_site(1.0)])) == "unfitted"          # lines, never fitted
    fitted = _recipe([_site(1.0)], fit_rmsd=0.0123)
    assert member_status(fitted) == "fitted"
    assert member_status(fitted, health_stale=True) == "edited"
    assert member_status(fitted, failed=True) == "failed"
    assert member_status(_recipe([_site(1.0)], fit_rmsd=float("nan")) ) == "unfitted"
    assert rmsd_of(fitted) == 0.0123 and rmsd_of(_recipe([])) is None
    assert rmsd_of(Recipe.from_dict(fitted)) == 0.0123
    # the edited test: a value change moves the signature, a label does not
    k0 = signature_key(fitted)
    assert len(k0) == 16 and signature_key(_recipe([])) == "" and signature_key(None) == ""
    other = copy.deepcopy(fitted)
    other["sites"][0]["label"] = "renamed"
    assert signature_key(other) == k0
    other["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] = 1.5
    assert signature_key(other) != k0
    assert signature_key(Recipe.from_dict(fitted)) == k0


# ------------------------------------------------------------------ the record
def _csv(path, pos, amp=100.0):
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel
    x = np.linspace(-20, 60, 500)
    tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                sites=[SiteModel(model="gauss_lor", label="A", params={
                    "isotropic_chemical_shift_ppm": Param(pos),
                    "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                    "gl": Param(1.0, vary=False)})])
    _, m, _ = engine.simulate(tr, exp_ppm=x)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# nucleus = 11B\n# larmor_MHz = 160\n")
        for xi, yi in zip(x, m):
            f.write(f"{xi:.4f} {yi:.4f}\n")
    return str(path), x, m


def test_series_spec_build_tags_and_restore(tmp_path):
    from larmor.loader import load_any
    from larmor.seriesmode import SeriesOptions, SeriesSpec, proc_number

    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    paths = [_csv(tmp_path / "a" / "glass.csv", 13.0)[0],
             _csv(tmp_path / "b" / "glass.csv", 15.0)[0],
             _csv(tmp_path / "b" / "other.csv", 17.0)[0]]
    recs = [load_any(p)[2] for p in paths]
    spec = SeriesSpec.build(paths, recs, names=[None, None, "typed"])
    assert spec.n == 3 and spec.id.startswith("series-")
    names = spec.names()
    assert names[2] == "typed" and names[0] != names[1]           # disambiguated
    assert names[0].startswith("glass") and names[1].startswith("glass")
    assert [m.group for m in spec.members] == ["glass", "glass", "typed"]
    assert spec.paths() == paths and [m.proc for m in spec.members] == ["", "", ""]
    assert spec.comparison is not None and spec.comparison.level == "none"   # CSV: nothing to compare
    assert spec.comparability_flag(0) == ""
    assert proc_number(r"C:\data\s\12\pdata\3\1r") == "3" and proc_number(paths[0]) == ""
    assert spec.keys() == [f"{spec.id}:{k}" for k in range(3)]
    assert spec.index_of(spec.keys()[1]) == 1 and spec.index_of("nope") is None

    spec.options = SeriesOptions(seed_on_move=False, carry_off=("gl",), passes=4,
                                 start="last", smooth=3)
    spec.members[1].fit_sig = "abcd"
    spec.members[2].failed = True
    spec.table = spec.series_table()
    spec.table.add_column(__import__("larmor.series_table", fromlist=["SeriesColumn"])
                          .SeriesColumn(key="cao", label="CaO (mol%)"),
                          {0: 0.0, 1: 5.0, 2: 10.0})
    tags = [spec.tag_for(k) for k in range(3)]
    json.dumps(tags)                                           # JSON-safe
    assert [t["index"] for t in tags] == [0, 1, 2] and all(t["id"] == spec.id for t in tags)
    assert "table" in tags[0] and "table" not in tags[1]
    assert tags[0]["options"] == {"seed_on_move": False, "carry_off": ["gl"], "passes": 4,
                                  "start": "last", "smooth": 3}
    # restore from shuffled tags: order by index, options, fit state, table columns
    back = SeriesSpec.from_tags([tags[2], tags[0], tags[1]])
    assert back.id == spec.id and back.names() == names and back.keys() == spec.keys()
    assert back.options == spec.options
    assert back.members[1].fit_sig == "abcd" and back.members[2].failed is True
    assert [c.key for c in back.series_table().columns] == ["cao"]
    assert back.series_table().x_values("cao")[0].tolist() == [0.0, 5.0, 10.0]
    # a bad option set falls back to the defaults instead of raising
    o = SeriesOptions.from_dict({"passes": 99, "start": "middle", "smooth": 4})
    assert (o.passes, o.start, o.smooth) == (2, "first", 0)


def test_series_spec_apply_table_reorders_and_renames(tmp_path):
    from larmor.loader import load_any
    from larmor.seriesmode import SeriesSpec

    paths = [_csv(tmp_path / f"s{k}.csv", 13.0 + 2 * k)[0] for k in range(3)]
    spec = SeriesSpec.build(paths, [load_any(p)[2] for p in paths])
    assert spec.names() == ["s0", "s1", "s2"]
    tbl = spec.series_table().permuted([2, 1, 0])
    tbl.rename(0, "last")
    tbl.set_group(1, "mid")
    assert spec.apply_table(tbl, [2, 1, 0]) is True
    assert spec.names() == ["last", "s1", "s0"]
    assert spec.paths() == [paths[2], paths[1], paths[0]]
    assert [m.group for m in spec.members] == ["last", "mid", "s0"]
    assert spec.tag_for(0)["name"] == "last" and spec.tag_for(2)["source_path"] == paths[0]
    assert spec.apply_table(spec.series_table(), [0, 1, 2]) is False
    with pytest.raises(ValueError):
        spec.apply_table(tbl, [0, 0, 1])
    assert spec.remove(spec.keys()[1]) and spec.n == 2 and not spec.remove("nope")


def test_auto_name_and_slug():
    import datetime as dt
    from larmor.seriesmode import auto_name, slug

    when = dt.datetime(2026, 10, 2, 14, 5)
    rec = _recipe([], sample="LAW 3Cl/1Ca")
    assert auto_name(rec, when=when) == "LAW_3Cl_1Ca_11B_seq_20261002_1405"
    assert auto_name(rec, name="member one", when=when) == "member_one_11B_seq_20261002_1405"
    assert auto_name({}, when=when) == "fit_seq_20261002_1405"
    assert slug("a b/c") == "a_b_c" and len(slug("x" * 200)) == 80


# ------------------------------------------------------------------ the sweep
def test_entries_round_trip_through_run_sequential(tmp_path):
    from larmor import seqfit
    from larmor.seriesmode import apply_sweep_result, entries_for_sweep

    series = [_csv(tmp_path / f"s{k}.csv", pos) for k, pos in enumerate([13.0, 15.0, 17.0])]
    model = _recipe([_site(12.0, fwhm=5.0, amp=80.0)])
    members = [{"recipe": copy.deepcopy(model), "ppm": x, "amp": y, "name": f"s{k}"}
               for k, (_p, x, y) in enumerate(series)]
    members[1]["recipe"]["fit_window_ppm"] = [50.0, -10.0]
    entries = entries_for_sweep(members, default_window=(60.0, -20.0))
    assert len(entries) == 3 and entries[0][3] == (60.0, -20.0) and entries[1][3] == (50.0, -10.0)
    assert entries[0][0].sites[0].params["isotropic_chemical_shift_ppm"].value == 12.0
    res = seqfit.run_sequential(entries, passes=2, start="first",
                                propagate=("isotropic_chemical_shift_ppm", "shift_fwhm_ppm"))
    applied = apply_sweep_result(members, res)
    assert all(a is not None for a in applied)
    pos = [a["recipe"]["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"]
           for a in applied]
    assert pos == pytest.approx([13.0, 15.0, 17.0], abs=0.4)
    assert all(np.isfinite(a["rmsd"]) and a["rmsd"] < 1.0 for a in applied)
    assert all(a["recipe"]["fit_rmsd"] is not None for a in applied)
    assert all(len(a["y_fit"]) == len(a["x"]) for a in applied)
    # the members' own dicts were never touched by the sweep
    assert members[0]["recipe"]["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] == 12.0

    # a stop after the first spectrum: only that member comes back fitted
    seen = []
    res2 = seqfit.run_sequential(
        [e for e in entries_for_sweep(members)], passes=2,
        progress=lambda p, k, r: seen.append(k), should_stop=lambda: len(seen) >= 1)
    applied2 = apply_sweep_result(members, res2)
    assert applied2[0] is not None and applied2[1] is None and applied2[2] is None

    members[2]["recipe"]["sites"] = []
    with pytest.raises(ValueError, match="none on s2"):
        entries_for_sweep(members)
    with pytest.raises(ValueError):
        apply_sweep_result(members[:2], res)
