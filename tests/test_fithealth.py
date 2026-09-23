"""fithealth: one Qt-free verdict assembled from sanity, identifiability and
the residual diagnostics (no Qt, synthetic data only)."""
import copy
import inspect
import json
import re
import types

import numpy as np

from larmor import fithealth
from larmor.recipe import Param, Recipe, SiteModel

# the arrays of tests/test_fit_diagnostics.py::test_residual_noise_ratio
X = np.linspace(-50, 150, 400)
PEAK = np.exp(-0.5 * ((X - 60) / 8) ** 2)
NOISE = np.random.RandomState(0).normal(0, 0.02, 400)
Y = PEAK + NOISE
WINDOW = (40.0, -10.0)          # test_sanity's window: the 15 ppm site sits inside


def _rec(**over):
    p = {"isotropic_chemical_shift_ppm": Param(15.0), "shift_fwhm_ppm": Param(6.0),
         "amplitude": Param(100.0), "gl": Param(0.5, vary=False)}
    p.update(over)
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0,
                  sites=[SiteModel(model="gauss_lor", label="A", params=p)])


def _lm(names=("s0_isotropic_chemical_shift_ppm", "s0_shift_fwhm_ppm",
               "s0_amplitude"), covar=None, redchi=1.1):
    """A fake lmfit result (the construction of tests/test_correlation_dialog.py)."""
    if covar is None:
        covar = np.eye(len(names))
    return types.SimpleNamespace(var_names=list(names), covar=covar, redchi=redchi)


# unit variances, so corr[0, 2] is -0.95 exactly (a degenerate pair)
_DEGEN_COV = np.array([[1.0, 0.0, -0.95], [0.0, 1.0, 0.1], [-0.95, 0.1, 1.0]])
_DEGEN_NAMES = ("s0_sigma_Cq_MHz", "s0_amplitude", "s1_amplitude")


def test_residual_noise_ratio_moved_intact():
    assert fithealth.residual_noise_ratio(Y, PEAK) < 1.5          # within noise
    assert fithealth.residual_noise_ratio(Y, 0.7 * PEAK) > 3      # structure left
    assert fithealth.residual_noise_ratio(Y[:30], PEAK[:30]) is None   # too short
    assert fithealth.residual_noise_ratio(PEAK, PEAK) is None     # zero edge noise


def test_clean_fit_is_ok():
    h = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01,
                         window=WINDOW)
    assert h.level == "ok" and h.flags == []
    assert h.pill_text() == "✓ Fit: no flags"
    assert h.chi_text() == "RMSD 0.0100 · χ²ᵣ 1.10"
    assert h.summary().startswith(
        "RMSD 0.0100 · χ²ᵣ 1.10   ·   residual within noise")
    tip = h.tooltip()
    assert tip.splitlines()[0].startswith("Fit health — last fit")
    assert "residual ≈ noise" in tip
    assert h.status_suffix() == "  ·  Fit: no flags"


def test_unphysical_and_degenerate_are_bad_and_lead():
    lm = _lm(_DEGEN_NAMES, covar=_DEGEN_COV)
    h = fithealth.assess(_rec(gl=Param(1.4, vary=False)), Y, PEAK, lmfit_result=lm,
                         at_bounds=["s0.shift_fwhm_ppm"], rmsd=0.01, window=WINDOW)
    assert h.level == "bad"
    assert [f.kind for f in h.flags] == ["physical", "degenerate", "at_bounds"]
    assert h.pill_text() == "✗ Fit: not physical · degenerate"
    phys, degen, bounds = h.flags
    assert phys.text == "unphysical ×1" and phys.target == "param"
    assert phys.params == [(0, "gl")]
    assert "Gauss/Lorentz mix = 1.4 outside [0, 1]" in phys.detail
    assert degen.text == "degenerate ×1: s0.σCq↔s1.amp (-0.95)"
    assert degen.target == "correlations" and degen.covariance_based
    assert bounds.text == "at bounds: s0.shift_fwhm_ppm"
    assert bounds.params == [(0, "shift_fwhm_ppm")]
    s = h.summary()                                   # legacy wording pinned
    for legacy in ("⚠ 1 physical warning",
                   "⚠ 1 unidentifiable pair (see Correlations)",
                   "⚠ at bounds: s0.shift_fwhm_ppm"):
        assert legacy in s, (legacy, s)


def test_structured_residual_and_noise_ratio_are_check():
    h = fithealth.assess(_rec(), Y, 0.7 * PEAK, lmfit_result=_lm(), rmsd=0.05,
                         window=WINDOW)
    assert h.level == "check"
    assert [f.kind for f in h.flags] == ["structured", "noise"]
    assert h.noise_ratio > 3 and h.runs_z < -3
    struct, noise = h.flags
    assert re.fullmatch(r"residual \d+\.\d× noise", noise.text)
    assert noise.target == "residual" and struct.target == "residual"
    assert "runs test z=" in struct.detail
    assert h.pill_text() == "⚠ Fit: 2 caveats"
    assert "⚠ structured residual" in h.summary()


def test_missing_covariance_and_population_are_check_and_fit_only():
    nocov = types.SimpleNamespace(var_names=["a"], covar=None, redchi=1.0)
    h = fithealth.assess(_rec(), Y, PEAK, lmfit_result=nocov, rmsd=0.01,
                         window=WINDOW)
    (f,) = h.flags
    assert f.kind == "nocov" and f.text == "no error bars (no covariance)"
    assert f.target == "errors" and f.covariance_based
    assert h.summary().endswith("⚠ no covariance (no error bars)")
    h2 = fithealth.assess(_rec(), Y, PEAK, lmfit_result=None, rmsd=0.01,
                          window=WINDOW)
    assert [f.kind for f in h2.flags] == ["nocov"]

    rows = [{"label": "C", "fraction_pct": 4.0, "fraction_err_pct": 5.0},
            {"label": "A", "fraction_pct": 90.0, "fraction_err_pct": 1.0}]
    h3 = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01,
                          window=WINDOW, quant_rows=rows)
    (pop,) = h3.flags
    assert pop.kind == "population" and pop.text == "population ±≥100 %: C"
    assert pop.target == "report"
    assert "⚠ population ±≥100 %: C" in h3.summary()
    none_rows = [{"label": "C", "fraction_pct": 4.0, "fraction_err_pct": None}]
    assert fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01,
                            window=WINDOW, quant_rows=none_rows).flags == []
    # an unfitted model has no covariance to judge: none of these appear
    h4 = fithealth.assess(_rec(), Y, PEAK, lmfit_result=None, at_bounds=["s0.gl"],
                          quant_rows=rows, fitted=False, window=WINDOW)
    assert not ({"nocov", "degenerate", "at_bounds", "population"} & h4.kinds())


def test_frozen_is_info_and_deduplicates_the_window_warning():
    rec = _rec()
    rec.sites.append(SiteModel(model="gauss_lor", label="C", params={
        "isotropic_chemical_shift_ppm": Param(200.0), "shift_fwhm_ppm": Param(6.0),
        "amplitude": Param(0.0), "gl": Param(0.5, vary=False)}))
    h = fithealth.assess(rec, Y, PEAK, lmfit_result=_lm(), frozen=["C"], rmsd=0.01,
                         window=WINDOW)
    (fr,) = h.flags                                   # one fact, reported once
    assert fr.kind == "frozen" and fr.level == "info" and fr.text == "frozen: C"
    assert fr.params == [(1, "isotropic_chemical_shift_ppm")]
    assert "physical" not in h.kinds()
    assert any("outside the fit window" in w["message"] for w in h.warns)
    assert "frozen: C" in h.summary()
    assert h.level == "ok"                            # info does not count


def test_signatures_and_live_reassessment():
    d = _rec().to_dict()
    d["sites"][0]["ref"] = {"ppm": [1.0, 2.0], "amp": [0.0, 1.0]}
    d["fit_zones"] = [[40.0, -10.0]]
    sig = fithealth.recipe_signature(d)
    d2 = copy.deepcopy(d)
    amp = d2["sites"][0]["params"]["amplitude"]
    amp["stderr"], amp["min"], amp["max"] = 3.0, 0.0, 1e9
    d2["sites"][0]["label"] = "renamed"
    d2["sites"][0]["ref"] = {"ppm": [1.0, 2.0], "amp": [0.0, 1.0]}   # equal list
    assert fithealth.recipe_signature(d2) == sig

    def changed(mutate):
        d3 = copy.deepcopy(d)
        mutate(d3)
        return fithealth.recipe_signature(d3) != sig

    assert changed(lambda x: x["sites"][0]["params"]["amplitude"]
                   .__setitem__("value", 100.0 + 1e-9))
    assert changed(lambda x: x["sites"][0]["params"]["gl"].__setitem__("vary", True))
    assert changed(lambda x: x["sites"][0]["params"]["shift_fwhm_ppm"]
                   .__setitem__("expr", "s0.amplitude"))
    assert changed(lambda x: x.__setitem__("fit_zones", [[40.0, 0.0]]))
    assert changed(lambda x: x["sites"].append(copy.deepcopy(x["sites"][0])))
    assert "json" not in inspect.getsource(fithealth)

    assert fithealth.data_signature(X.copy(), Y.copy()) == fithealth.data_signature(X, Y)
    assert fithealth.data_signature(X + 1.0, Y) != fithealth.data_signature(X, Y)
    assert fithealth.data_signature(np.array([]), np.array([])) == ()

    h_fit = fithealth.assess(d, Y, PEAK, lmfit_result=_lm(_DEGEN_NAMES, _DEGEN_COV),
                             rmsd=0.01, window=WINDOW)
    h_fit.recipe_sig, h_fit.data_sig = sig, fithealth.data_signature(X, Y)
    edited = copy.deepcopy(d)
    edited["sites"][0]["params"]["amplitude"]["value"] = 70.0
    h = fithealth.reassess_live(h_fit, edited, Y, 0.7 * PEAK, ppm=X, window=WINDOW)
    assert h.stale and h.fitted
    assert h.pill_text().endswith(" · edited since fit")
    assert h.pill_text().startswith(("⚠ Model:", "✗ Model:", "Model:"))
    degen = [f for f in h.flags if f.kind == "degenerate"]
    assert degen and degen[0].stale                   # carried from the fit
    assert {"structured", "noise"} <= h.kinds()       # recomputed live
    assert h.rmsd == 0.01 and h.recipe_sig == sig
    # the same values again: the fit verdict itself comes back (identity)
    assert fithealth.reassess_live(h_fit, d, Y, PEAK, ppm=X, window=WINDOW) is h_fit
    # before any fit only the physical checks run
    unphys = copy.deepcopy(d)
    unphys["sites"][0]["params"]["gl"]["value"] = 1.4
    h0 = fithealth.reassess_live(None, unphys, Y, 0.5 * PEAK, ppm=X, window=WINDOW)
    assert h0.level == "bad" and h0.pill_text() == "✗ Model: not physical"
    assert not ({"noise", "structured"} & h0.kinds())
    assert h0.tooltip().splitlines()[0] == "Model health — not fitted yet (F5 to fit)"
    hn = fithealth.reassess_live(None, d, Y, 0.5 * PEAK, ppm=X, window=WINDOW)
    assert hn.level == "none" and hn.pill_text() == fithealth.NO_FIT_TEXT
    assert hn.pill_text() == "no fit yet — F5 fits the lines"


def test_assess_accepts_dict_and_does_not_mutate_it():
    d = _rec().to_dict()
    d["fit_window_ppm"] = [40.0, 20.0]                # a JSON list; 15 ppm is below 20
    before = json.dumps(d)
    h = fithealth.assess(d, Y, PEAK, lmfit_result=_lm(), rmsd=0.01)
    assert json.dumps(d) == before
    phys = [f for f in h.flags if f.kind == "physical"]
    assert phys and "is outside the fit window" in phys[0].detail
    assert fithealth.assess(Recipe.from_dict(d), Y, PEAK, lmfit_result=_lm(),
                            rmsd=0.01).level == h.level == "bad"


def test_pill_text_is_bounded_for_any_flag_load():
    sites = [SiteModel(model="quad_ct", label=f"S{i}", params={
        "isotropic_chemical_shift_ppm": Param(15.0), "Cq_MHz": Param(3.0),
        "eta": Param(1.5), "line_fwhm_ppm": Param(-1.0), "amplitude": Param(100.0)})
        for i in range(12)]
    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.0, sites=sites)
    n = 40
    names = [f"s{i}_amplitude" if i % 2 else f"s{i}_isotropic_chemical_shift_ppm"
             for i in range(n)]
    cov = np.full((n, n), 0.99)
    np.fill_diagonal(cov, 1.0)
    lm = types.SimpleNamespace(var_names=names, covar=cov, redchi=2.0)
    bounds = [f"s{i}.amplitude" for i in range(30)]
    h = fithealth.assess(rec, Y, 0.7 * PEAK, lmfit_result=lm, at_bounds=bounds,
                         rmsd=0.1, window=WINDOW)
    assert len(h.warns) == 24 and len(h.pairs) == 780
    assert len(h.pill_text()) <= 40
    assert h.pill_text() == "✗ Fit: not physical · degenerate"
    assert all(len(f.text) <= 60 for f in h.flags), [f.text for f in h.flags]
    assert h.summary_tooltip().count("↔") == 8            # legacy uni[:8]
    assert fithealth.short_name("s0_amplitude") == "s0.amp"
    assert fithealth.parse_bound_name("s3.sigma_Cq_MHz") == (3, "sigma_Cq_MHz")
    assert fithealth.parse_bound_name("amplitude") is None


# ------------------------------------------------- quantitativity chips
def _facts(t1=True, spin=1.5, title="11B with short tip angle", p90=None, p1=0.425,
           d1=14.0, kind="Single pulse", pulprog="zg"):
    from larmor import quantitativity as Q
    from larmor.satrec import T1Region

    acq = Q.Acquisition(expno="/x/24", nucleus="11B", pulprog=pulprog, kind=kind,
                        ns=256, d1_s=d1, aq_s=0.04, p1_us=p1, plw1_w=100.0,
                        probhd="16_Solenoid (PMAS16)", title=title,
                        p90_us_title=p90, flip_deg_title=None, t1_multiple_claimed=None)
    src = Q.T1Source(expno="/x/23", kind="ct1t2", regions=[
        T1Region(1, 27.382, 5.846, 4.643, (4.643,), 1.207e-2, 7),
        T1Region(2, 5.846, -7.618, 3.896, (3.896,), 1.350e-2, 7)],
        vdlist_max_s=256.0, title_notes_s=[4.64, 3.74]) if t1 else None
    return Q.AcqFacts(acquisition=acq, t1=src, sibling_expno="/x/23", spin=spin,
                      folder="/x", t1_status="ok" if t1 else "missing")


def test_tail_flag_is_check_targets_widen_and_is_carried_stale():
    rows = [{"site": "s0", "label": "AlO6", "fraction_pct": 20.0, "fraction_err_pct": 1.0,
             "tail_outside_pct": 14.0},
            {"site": "s1", "label": "AlO4", "fraction_pct": 80.0, "fraction_err_pct": 1.0,
             "tail_outside_pct": 0.3}]
    h = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                         quant_rows=rows)
    (f,) = h.flags
    assert (f.kind, f.level, f.text, f.target) == \
        ("tail", "check", "tail outside window: AlO6 14 %", "widen")
    assert f.params == [(0, "amplitude")] and f.covariance_based
    assert "14 % of AlO6's simulated area" in f.detail and "limit 2 %" in f.detail
    assert "⚠ tail outside window: AlO6 14 %" in h.summary()
    assert h.level == "check" and h.tail_checked
    assert "tails inside the window" not in h.passing()
    inside = [dict(r, tail_outside_pct=0.4) for r in rows]
    h2 = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                          quant_rows=inside)
    assert h2.flags == [] and "tails inside the window (≤ 2 %)" in h2.passing()
    legacy = [{"label": "C", "fraction_pct": 40.0, "fraction_err_pct": 1.0}]
    h3 = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                          quant_rows=legacy)
    assert h3.flags == [] and not h3.tail_checked
    assert "tails inside" not in " ".join(h3.passing())
    h4 = fithealth.assess(_rec(), Y, PEAK, quant_rows=rows, fitted=False, window=WINDOW)
    assert "tail" not in h4.kinds()
    # a live edit carries the fit's tail verdict, dimmed
    d = _rec().to_dict()
    h.recipe_sig = fithealth.recipe_signature(d)
    h.data_sig = fithealth.data_signature(X, Y)
    edited = copy.deepcopy(d)
    edited["sites"][0]["params"]["amplitude"]["value"] = 70.0
    live = fithealth.reassess_live(h, edited, Y, 0.7 * PEAK, ppm=X, window=WINDOW)
    (tail,) = [f for f in live.flags if f.kind == "tail"]
    assert tail.stale and live.stale and live.tail_checked
    assert fithealth.KIND_ORDER.index("tail") > fithealth.KIND_ORDER.index("population")
    assert fithealth.KIND_ORDER.index("excitation") < fithealth.KIND_ORDER.index("frozen")
    assert fithealth.KIND_ORDER.index("tail") < fithealth.KIND_ORDER.index("recovery")


def test_acquisition_flags_are_live_leveled_and_targeted():
    from larmor import quantitativity as Q

    facts = _facts()
    h = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                         acquisition=facts)
    assert [f.kind for f in h.flags] == ["recovery", "excitation"]
    rec, exc = h.flags
    assert (rec.level, rec.text, rec.target) == \
        ("check", "D1 = 3.0 T1 → 95 % (90° assumed)", "relaxation")
    assert "EXPNO 23" in rec.detail and "4.64 s" in rec.detail
    assert (exc.level, exc.text, exc.target) == ("info", "flip angle unknown (I = 3/2)", "flip")
    assert not rec.covariance_based and not exc.covariance_based
    assert h.pill_text() == "⚠ Fit: 1 caveat"            # info does not count
    assert isinstance(h.acquisition, Q.Check) and h.acquisition.facts is facts
    assert "⚠ D1 = 3.0 T1 → 95 % (90° assumed)" in h.summary()
    assert "acquisition not checked" not in h.tooltip()
    # live: shown before the first fit too
    h0 = fithealth.assess(_rec(), None, None, fitted=False, window=WINDOW, acquisition=facts)
    assert {"recovery", "excitation"} <= h0.kinds() and h0.level == "check"
    # a live edit recomputes the site -> region mapping from the cached facts
    d = _rec().to_dict()
    h.recipe_sig, h.data_sig = fithealth.recipe_signature(d), fithealth.data_signature(X, Y)
    moved = copy.deepcopy(d)
    moved["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] = 1.0
    live = fithealth.reassess_live(h, moved, Y, PEAK, ppm=X, window=WINDOW)
    (rec2,) = [f for f in live.flags if f.kind == "recovery"]
    assert rec2.text == "D1 = 3.6 T1 → 97 % (90° assumed)" and not rec2.stale
    assert live.acquisition.facts is facts                  # reused from prev
    assert live.acquisition.recoveries[0].t1_s == 3.896
    # no usable T1: a grey info chip that never colours the pill
    hn = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                          acquisition=_facts(t1=False))
    (f,) = [f for f in hn.flags if f.kind == "recovery"]
    assert f.level == "info" and f.text == "recycle 14 s — T1 unknown"
    assert hn.level == "ok"
    # spin-1/2: no excitation chip, a passing line instead
    hp = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                          acquisition=_facts(spin=0.5))
    assert "excitation" not in hp.kinds()
    assert "spin-½: any flip angle is quantitative" in hp.passing()
    # the title's P1(90) with a long pulse: an amber excitation chip
    hx = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                          acquisition=_facts(title="P1(90)=2.5", p90=2.5, p1=0.5))
    (f,) = [f for f in hx.flags if f.kind == "excitation"]
    assert f.level == "check" and f.text == "flip 18° > 15° limit (I = 3/2)"
    assert f.target == "flip" and "⚠ flip 18° > 15° limit (I = 3/2)" in hx.summary()
    # a long enough delay and a typed 90-degree pulse: both checks pass
    r = _rec()
    r.provenance = {"quantitativity": {"p90_us": 3.5}}
    hok = fithealth.assess(r, Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW,
                           acquisition=_facts(d1=5 * 4.643))
    assert not ({"recovery", "excitation"} & hok.kinds())
    passing = hok.passing()
    assert any(p.startswith("recycle 23.3 s = 5.0 T1 → ") and "(EXPNO 23)" in p
               for p in passing), passing
    assert "flip 10.9° within the linear regime (≤ 15°)" in passing
    # not a Bruker dataset: nothing judged, and the tooltip says so
    hno = fithealth.assess(_rec(), Y, PEAK, lmfit_result=_lm(), rmsd=0.01, window=WINDOW)
    assert not ({"recovery", "excitation"} & hno.kinds()) and hno.acquisition is None
    assert "acquisition not checked (not a Bruker dataset)" in hno.tooltip()
    assert "acquisition not checked" not in fithealth.assess(
        _rec(), None, None, fitted=False, window=WINDOW).tooltip()


def test_with_quantification_rebuilds_only_quant_flags():
    lm = _lm(_DEGEN_NAMES, covar=_DEGEN_COV)
    rows = [{"site": "s0", "label": "A", "fraction_pct": 50.0, "fraction_err_pct": 1.0,
             "tail_outside_pct": 14.0}]
    prev = fithealth.assess(_rec(), Y, PEAK, lmfit_result=lm, rmsd=0.01, window=WINDOW,
                            quant_rows=rows, acquisition=_facts())
    prev.recipe_sig, prev.data_sig = ("sig",), ("data",)
    assert [f.kind for f in prev.flags] == ["degenerate", "tail", "recovery", "excitation"]
    degen = prev.flags[0]
    fresh = [dict(rows[0], tail_outside_pct=0.2)]
    h = fithealth.with_quantification(prev, _rec(), fresh)
    assert [f.kind for f in h.flags] == ["degenerate", "recovery", "excitation"]
    assert h.flags[0] == degen
    assert (h.recipe_sig, h.data_sig, h.stale, h.rmsd, h.fitted) == \
        (("sig",), ("data",), False, 0.01, True)
    assert h.acquisition.facts is prev.acquisition.facts and h.tail_checked
    assert prev.flags[1].kind == "tail"                     # prev untouched
    h2 = fithealth.with_quantification(prev, _rec().to_dict(), None)
    assert [f.kind for f in h2.flags] == ["degenerate", "recovery", "excitation"]
    assert not h2.tail_checked
    # a typed 90-degree pulse re-judges the recovery without a refit
    r = _rec()
    r.provenance = {"quantitativity": {"p90_us": 3.5, "t1_s": 2.0}}
    h3 = fithealth.with_quantification(prev, r, fresh)
    assert [f.kind for f in h3.flags] == ["degenerate"]
    assert any("recycle" in p for p in h3.passing())
    # on an unfitted verdict nothing fit-only appears
    unfit = fithealth.assess(_rec(), None, None, fitted=False, window=WINDOW,
                             acquisition=_facts())
    h4 = fithealth.with_quantification(unfit, _rec(), rows)
    assert "tail" not in h4.kinds() and "recovery" in h4.kinds()
