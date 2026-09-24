"""The derived parameter status († fixed · ‡ at a bound · § linked) and its
agreement with the fit's own at-bound rule (larmor.paramstatus)."""
import re
from pathlib import Path

import numpy as np
import pytest

from larmor import paramstatus
from larmor import fit as fitmod
from larmor.paramstatus import ParamStatus, bound_side, param_status
from larmor.recipe import Recipe, SiteModel, Param


def test_bound_side_matches_the_fit_rule():
    """1e-3 · max(1, |value|) below/above a FINITE bound; None / ±inf never."""
    assert bound_side(0.0, 0.0, None) == "min"
    assert bound_side(12.6, 0, 12.6) == "max"
    # the P5-Bi0 27Al sigma is NOT at its bound
    assert bound_side(1.2691513, 1.2172446, 1.4004858) is None
    assert bound_side(10.373903484537355, 10.373903484537353, 11.93556637468276) == "min"
    assert bound_side(3.52, 3.51, None) is None          # 0.01 > 1e-3 * 3.52
    assert bound_side(5.0, None, None) is None
    assert bound_side(5.0, -np.inf, np.inf) is None
    assert bound_side(5.0, float("-inf"), float("inf")) is None
    assert bound_side(float("nan"), 0.0, 1.0) is None
    # the tolerance scales with max(1, |value|)
    assert bound_side(2000.0, 1999.0, None) == "min"
    assert bound_side(0.5, 0.4, None) is None            # span floors at 1
    assert bound_side(0.5, 0.4995, None) == "min"
    # min wins when both sides match (degenerate min == max)
    assert bound_side(0.0, 0.0, 0.0) == "min"


def test_param_status_classifies_param_dict_and_number_inputs():
    m = "gauss_lor"
    assert param_status(m, "shift_fwhm_ppm", Param(5.0)).kind == "free"
    assert param_status(m, "shift_fwhm_ppm", Param(5.0, vary=False)).kind == "fixed"
    st = param_status(m, "shift_fwhm_ppm",
                      Param(5.0, vary=False, expr="s0.shift_fwhm_ppm"))
    assert st.kind == "linked" and st.expr == "s0.shift_fwhm_ppm"   # expr wins
    st = param_status(m, "shift_fwhm_ppm", Param(3.5100000000000002, min=3.51, max=4.29))
    assert st == ParamStatus("at_bound", side="min", bound=3.51)
    st = param_status(m, "shift_fwhm_ppm", Param(7.149999999999857, min=5.85, max=7.15))
    assert st.kind == "at_bound" and st.side == "max" and st.bound == 7.15

    # the desktop's {value, stderr, vary, min, max, expr} dicts read the same
    cases = [
        ({"value": 5.0, "stderr": None, "vary": True, "min": None, "max": None,
          "expr": None}, Param(5.0)),
        ({"value": 5.0, "stderr": 0.0, "vary": False, "min": None, "max": None,
          "expr": None}, Param(5.0, vary=False)),
        ({"value": 5.0, "stderr": None, "vary": False, "min": None, "max": None,
          "expr": "s0.shift_fwhm_ppm"}, Param(5.0, vary=False, expr="s0.shift_fwhm_ppm")),
        ({"value": 3.51, "stderr": None, "vary": True, "min": 3.51, "max": 4.29,
          "expr": None}, Param(3.51, min=3.51, max=4.29)),
    ]
    for d, p in cases:
        assert param_status(m, "shift_fwhm_ppm", d) == param_status(m, "shift_fwhm_ppm", p)
    # a bare number carries no constraint information
    assert param_status(m, "shift_fwhm_ppm", 5.0).kind == "free"
    assert param_status(m, "shift_fwhm_ppm", 0.0).kind == "free"
    assert param_status(m, "shift_fwhm_ppm", None).kind == "free"

    # vocabularies
    flags = {param_status(m, "shift_fwhm_ppm", p).csv_flag for p in (
        Param(5.0), Param(5.0, vary=False), Param(5.0, expr="s0.shift_fwhm_ppm"),
        Param(3.51, min=3.51), Param(7.15, min=5.85, max=7.15))}
    assert flags == {"", "fixed", "linked", "at_min", "at_max"}
    assert (ParamStatus("fixed").marker, ParamStatus("fixed").latex) == ("†", r"$^{\dagger}$")
    assert (ParamStatus("at_bound", "min", 1.0).marker,
            ParamStatus("at_bound", "min", 1.0).latex) == ("‡", r"$^{\ddagger}$")
    assert (ParamStatus("linked", expr="x").marker,
            ParamStatus("linked", expr="x").latex) == ("§", r"$^{\S}$")
    assert (ParamStatus("free").marker, ParamStatus("free").latex) == ("", "")
    assert ParamStatus("fixed").word == "held fixed"
    assert ParamStatus("at_bound", "min", 10.373903).word == "at its lower bound 10.37"
    assert ParamStatus("at_bound", "max", 12.6).word == "at its upper bound 12.6"
    assert ParamStatus("linked", expr="0.19 * s0.amplitude").word == \
        "linked: 0.19 * s0.amplitude"
    assert ParamStatus("free").word == ""


def test_param_status_uses_registry_bounds_like_the_fit():
    # registry ParamDef min = 0.0 for the amplitude, recipe min None
    st = param_status("gauss_lor", "amplitude", Param(0.0))
    assert st == ParamStatus("at_bound", side="min", bound=0.0)
    assert param_status("no_such_model", "amplitude", Param(0.0)).kind == "free"
    assert param_status(None, "amplitude", Param(0.0)).kind == "free"
    assert paramstatus.effective_bounds("gauss_lor", "gl", Param(0.5)) == (0.0, 1.0)
    # an explicit recipe bound is honoured over the registry
    assert paramstatus.effective_bounds("gauss_lor", "gl", Param(0.5, min=0.2)) == (0.2, 1.0)
    assert paramstatus.effective_bounds("gauss_lor", "no_such_param", Param(0.5)) == (None, None)

    # fit._make_params' fallback semantics are unchanged (it now calls the
    # same helper): width min None -> the model's floor 0.1
    from test_param_bounds_fallback import _recipe
    r = _recipe(None)
    p = fitmod._make_params(r)
    name = fitmod._lmfit_name(0, r.sites[0], "shift_fwhm_ppm")
    assert p[name].min == 0.1
    assert p[fitmod._lmfit_name(0, r.sites[0], "amplitude")].min == 0.0
    # a FIXED parameter is still widened to unbounded (lmfit's min != max check)
    r2 = _recipe(None)
    r2.sites[0].params["amplitude"] = Param(0.0, vary=False, min=0.0, max=0.0)
    p2 = fitmod._make_params(r2)
    assert p2[fitmod._lmfit_name(0, r2.sites[0], "amplitude")].min == -np.inf


def _capped_recipe(cap):
    return Recipe(nucleus="27Al", larmor_frequency_MHz=130.32, sites=[
        SiteModel(model="gauss_lor", label="g", params={
            "isotropic_chemical_shift_ppm": Param(0.0, vary=False),
            "shift_fwhm_ppm": Param(8.0, vary=False),
            "amplitude": Param(3.0, min=0.0, max=cap),
            "gl": Param(1.0, vary=False)})])


def test_param_status_agrees_with_fit_at_bounds_and_notes_are_pruned():
    """The capped-amplitude synthetic of tests/test_constraints.py: the
    derived at-bound set equals fit.fit's own list; a refit that no longer
    pins clears the note (and the marker) instead of accumulating notes."""
    x = np.linspace(-60.0, 60.0, 2001)
    y = 10.0 * np.exp(-4 * np.log(2) * (x / 8.0) ** 2)
    recipe = _capped_recipe(5.0)
    result = fitmod.fit(recipe, x, y)
    statuses = paramstatus.recipe_statuses(recipe)
    derived = {f"s{i}.{pn}" for (i, pn), st in statuses.items() if st.kind == "at_bound"}
    assert derived == set(result.at_bounds) == {"s0.amplitude"}
    st = statuses[(0, "amplitude")]
    assert st.side == "max" and st.bound == 5.0
    for pn in ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "gl"):
        assert statuses[(0, pn)].kind == "fixed"
    notes = [n for n in recipe.notes if n.startswith(paramstatus.AT_BOUND_NOTE_PREFIX)]
    assert len(notes) == 1 and "s0.amplitude" in notes[0]
    assert paramstatus.summary(recipe) == "3 fixed · 1 at a bound"

    # a second, identical fit must not add a second copy of the note
    fitmod.fit(recipe, x, y)
    assert sum(n.startswith(paramstatus.AT_BOUND_NOTE_PREFIX) for n in recipe.notes) == 1

    # raise the cap and refit the SAME recipe: unpinned, note gone, marker gone
    recipe.sites[0].params["amplitude"].max = 20.0
    result2 = fitmod.fit(recipe, x, y)
    assert result2.at_bounds == []
    assert paramstatus.recipe_statuses(recipe)[(0, "amplitude")].kind == "free"
    assert not any(n.startswith(paramstatus.AT_BOUND_NOTE_PREFIX) for n in recipe.notes)
    assert paramstatus.summary(recipe) == "3 fixed"


def _final2_shaped_dict():
    """A recipe DICT with the shapes of the PBi Final2 P5-Bi0 31P file: a
    fully fixed site s1 (amplitude value = min = max = 0, stderr 0.0 on every
    held value), s0's width and s3's position at their lower bounds, every
    gl fixed with stderr 0.0. No extra keys anywhere."""
    def P(value, stderr=None, vary=True, mn=None, mx=None, expr=None):
        return {"value": value, "stderr": stderr, "vary": vary,
                "min": mn, "max": mx, "expr": expr}

    def site(label, pos, width, amp, **over):
        params = {"isotropic_chemical_shift_ppm": pos, "shift_fwhm_ppm": width,
                  "amplitude": amp, "gl": P(0.6, 0.0, False)}
        params.update(over)
        return {"model": "gauss_lor", "label": label, "params": params,
                "ref": None, "func": None}

    return {"sample": "P5-Bi0", "nucleus": "31P", "larmor_frequency_MHz": 242.9,
            "spin_rate_Hz": 20000.0, "sites": [
                site("Q0", P(3.9, 0.02, mn=2.79, mx=4.29), P(3.51, None, True, 3.51, 4.29),
                     P(1200.0, 30.0, True, 0.0, None)),
                site("Bi", P(-5.0, 0.0, False), P(9.0, 0.0, False),
                     P(0.0, 0.0, False, 0.0, 0.0)),
                site("Q1", P(-8.0, 0.03, mn=-9.5, mx=-6.5), P(7.0, 0.1, True, 5.85, 7.15),
                     P(900.0, 20.0, True, 0.0, None)),
                site("Q1b", P(2.7900000000023404, None, True, 2.79, 4.29),
                     P(4.0, 0.1, True, 3.51, 4.29), P(300.0, 15.0, True, 0.0, None)),
            ], "notes": [], "larmor_recipe_version": 1}


def test_legacy_final2_shaped_dict_footnote_and_summary():
    d = _final2_shaped_dict()
    st = paramstatus.recipe_statuses(d)
    assert set(st) == {(i, pn) for i in range(4) for pn in
                       ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "amplitude", "gl")}
    assert all(st[(1, pn)].kind == "fixed" for pn in
               ("isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "amplitude", "gl"))
    assert st[(3, "isotropic_chemical_shift_ppm")] == ParamStatus("at_bound", "min", 2.79)
    assert st[(0, "shift_fwhm_ppm")] == ParamStatus("at_bound", "min", 3.51)
    assert all(st[(i, "gl")].kind == "fixed" for i in range(4))
    assert st[(0, "isotropic_chemical_shift_ppm")].kind == "free"
    assert st[(2, "shift_fwhm_ppm")].kind == "free"           # 7.0 vs max 7.15
    assert st[(0, "amplitude")].kind == "free"                # 1200 vs min 0

    c3 = paramstatus.site_constraints(d, 3)
    assert c3["at_bound"] == [("isotropic_chemical_shift_ppm", "min", pytest.approx(2.79))]
    assert c3["fixed"] == ["gl"] and c3["linked"] == []
    c1 = paramstatus.site_constraints(d, 1)
    assert c1["fixed"] == ["isotropic_chemical_shift_ppm", "shift_fwhm_ppm", "amplitude", "gl"]
    assert c1["at_bound"] == []
    # the same dict through Recipe.from_dict reads identically
    assert paramstatus.recipe_statuses(Recipe.from_dict(d)) == st

    entries = [("Bi FWHM (ppm)", st[(1, "shift_fwhm_ppm")]),
               ("Q1b δiso (ppm)", st[(3, "isotropic_chemical_shift_ppm")]),
               ("Q0 FWHM (ppm)", st[(0, "shift_fwhm_ppm")]),
               ("X amplitude", ParamStatus("linked", expr="0.19 * s0.amplitude")),
               ("Q0 δiso (ppm)", st[(0, "isotropic_chemical_shift_ppm")])]
    tex = paramstatus.footnote(entries, "latex")
    assert tex == (r"$\dagger$ held fixed. "
                   r"$\ddagger$ finished at a bound: Q1b δiso (ppm) (lower bound 2.79); "
                   r"Q0 FWHM (ppm) (lower bound 3.51). "
                   r"$\S$ linked: X amplitude = 0.19 * s0.amplitude.")
    txt = paramstatus.footnote(entries, "text")
    assert txt.startswith("† held fixed. ‡ finished at a bound: ") and "§ linked: " in txt
    assert "\\" not in txt
    # only the kinds present, in the order fixed / at bound / linked
    only_linked = paramstatus.footnote([entries[3]])
    assert only_linked == r"$\S$ linked: X amplitude = 0.19 * s0.amplitude."
    assert "held fixed" not in only_linked and "at a bound" not in only_linked
    assert paramstatus.footnote([]) == ""
    assert paramstatus.footnote([entries[4]]) == ""            # all free

    assert paramstatus.summary(d) == "7 fixed · 2 at a bound"   # 4 gl + 3 of site Bi
    from test_methods import _czjzek_recipe
    assert paramstatus.summary(_czjzek_recipe()) == ""
    assert paramstatus.summary(_czjzek_recipe().to_dict()) == ""
    assert paramstatus.summary({"sites": []}) == ""

    # csv_fields: the five long-CSV cells from a Param, a dict or a table row
    assert paramstatus.csv_fields(Param(3.51, min=3.51, max=4.29), ParamStatus("at_bound", "min", 3.51)) \
        == ["True", "3.51", "4.29", "", "min"]
    assert paramstatus.csv_fields(Param(1.0, vary=False), ParamStatus("fixed")) \
        == ["False", "", "", "", ""]
    assert paramstatus.csv_fields({"vary": None, "min": None, "max": None, "expr": "",
                                   "at_bound": ""}) == ["", "", "", "", ""]
    assert paramstatus.csv_fields({"vary": True, "min": 14.925, "max": 15.075,
                                   "expr": "", "at_bound": "max"}) \
        == ["True", "14.925", "15.075", "", "max"]
    assert paramstatus.csv_fields(Param(2.0, expr="0.5 * s0.amplitude"), "") \
        == ["True", "", "", "0.5 * s0.amplitude", ""]


FINAL2 = Path("C:/Users/samso/Desktop/WSU_work/Project_pers/PBi/01_Data/NMR/fitting/Final2")


def test_final2_recipes_derived_status_reproduces_the_fit_notes():
    """READ-ONLY on the 32 PBi Final2 recipes: the derived at-bound set of
    every recipe equals the set its own fit wrote into the recipe note --
    the real-data agreement the design was verified on. Skips when the
    local data is absent."""
    from conftest import require
    require(FINAL2)
    files = sorted(FINAL2.glob("*.recipe.json"))
    if len(files) < 32:
        pytest.skip(f"expected 32 Final2 recipes, found {len(files)}")
    disagreements, n_at, n_fixed, n_linked, n_free = [], 0, 0, 0, 0
    for f in files:
        rec = Recipe.load(f)
        noted: set[str] = set()
        for n in rec.notes:
            if n.startswith(paramstatus.AT_BOUND_NOTE_PREFIX):
                noted |= {t.strip() for t in n.split("):", 1)[1].split(",") if t.strip()}
        derived = set()
        for (i, pn), st in paramstatus.recipe_statuses(rec).items():
            if st.kind == "at_bound":
                derived.add(f"s{i}.{pn}"); n_at += 1
            elif st.kind == "fixed":
                n_fixed += 1
            elif st.kind == "linked":
                n_linked += 1
            else:
                n_free += 1
        if derived != noted:
            disagreements.append((f.name, sorted(derived ^ noted)))
    assert not disagreements, disagreements
    assert n_at > 0 and n_fixed > 0
    assert re.fullmatch(r"\d+ fixed( · \d+ at a bound)?",
                        paramstatus.summary(Recipe.load(files[0])))


def test_default_fixed_lb_is_the_silent_default_kind():
    """A Czjzek lb still at its registry default (pinned, as a fresh site is
    created) is 'default': held (no error bar) but no glyph, no footnote, not
    counted -- a parameter fixed at its model default is not a user
    constraint. Change its value or free it and the ordinary kinds apply."""
    from larmor import models

    assert paramstatus.is_default_fixed("czjzek", "line_fwhm_ppm")
    assert not paramstatus.is_default_fixed("quad_ct", "line_fwhm_ppm")
    assert not paramstatus.is_default_fixed("gauss_lor", "gl")
    assert not paramstatus.is_default_fixed("no_such_model", "line_fwhm_ppm")
    fresh = models.get("czjzek").defaults()
    st = param_status("czjzek", "line_fwhm_ppm", fresh["line_fwhm_ppm"])
    assert st.kind == "default" and st.held
    assert (st.marker, st.latex, st.csv_flag) == ("", "", "fixed")
    assert st.word == "held at the model default"
    # the desktop dict shape and stderr 0.0 read the same
    d = {"value": 0.0, "stderr": 0.0, "vary": False, "min": 0.0, "max": None,
         "expr": None}
    assert param_status("czjzek", "line_fwhm_ppm", d).kind == "default"
    # moved off the default while pinned: a user's fixed value
    assert param_status("czjzek", "line_fwhm_ppm",
                        Param(2.0, vary=False)).kind == "fixed"
    # freed: free, or at its bound like any other free parameter
    assert param_status("czjzek", "line_fwhm_ppm", Param(2.0)).kind == "free"
    assert param_status("czjzek", "line_fwhm_ppm", Param(0.0)).kind == "at_bound"
    # gl of a Gauss/Lorentz is fixed by default too but NOT declared
    # default_fixed: it keeps its dagger
    assert param_status("gauss_lor", "gl", Param(1.0, vary=False)).kind == "fixed"
    # the recipe-level views stay silent about it
    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.0, sites=[
        SiteModel(model="czjzek", label="Al4", params=fresh)])
    assert paramstatus.summary(rec) == ""
    assert paramstatus.site_constraints(rec, 0) == {"fixed": [], "linked": [],
                                                    "at_bound": []}
    assert paramstatus.footnote([("Al4 lb", st)]) == ""
    assert paramstatus.ParamStatus("fixed").held and not paramstatus.ParamStatus("free").held
