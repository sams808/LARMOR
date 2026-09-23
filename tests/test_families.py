"""Site families (larmor.families): summed populations and named ratios with
their error on three stated bases -- covariance (lmfit uvars), Monte-Carlo
(per-trial re-integrated sums) and the flagged independent fallback.

The analytic checks pin the arithmetic: a two-line N4 against the exact
first-order formula over the amplitude covariance, and the Monte-Carlo basis
against the plain spread of the per-trial ratio.
"""
import numpy as np
import pytest

from larmor import families
from larmor.quantify import quantify, site_integrals
from larmor.recipe import Param, Recipe, SiteModel


def _gl(pos, width, amp, label, family="", stderr=None):
    return SiteModel(model="gauss_lor", label=label, family=family, params={
        "isotropic_chemical_shift_ppm": Param(pos),
        "shift_fwhm_ppm": Param(width),
        "amplitude": Param(amp, stderr=stderr),
        "gl": Param(1.0, vary=False)})


def _b_recipe(sites):
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                  sites=sites)


def test_two_site_family_error_matches_analytic_covariance():
    """N4 = b/(a+b) from two equal-width lines: the covariance basis must
    reproduce sqrt((b²σa² + a²σb² − 2ab·cov)/(a+b)⁴) exactly (equal widths
    make integral ∝ amplitude, so the k_i cancel); the flagged independent
    basis gives a different number -- with ρ < 0 the RATIO is known worse
    (the split between the lines is what anticorrelation blurs) while the
    SUM of the pair is known better than its quadrature."""
    from uncertainties import correlated_values

    a, b, sa, sb, rho = 100.0, 50.0, 4.0, 3.0, -0.6
    cov = np.array([[sa ** 2, rho * sa * sb], [rho * sa * sb, sb ** 2]])
    rec = _b_recipe([_gl(15.0, 6.0, a, "A", "BO3", stderr=sa),
                     _gl(0.0, 6.0, b, "B", "BO4", stderr=sb)])
    ua, ub = correlated_values([a, b], cov)
    q = quantify(rec, (40.0, -20.0), uvars={0: ua, 1: ub})
    assert q["family_basis"] == "covariance"
    n4 = next(r for r in q["ratios"] if r["name"] == "N4")
    assert n4["defined"] and n4["value"] == pytest.approx(b / (a + b), rel=1e-9)
    exact = np.sqrt((b * b * sa ** 2 + a * a * sb ** 2 - 2 * a * b * cov[0, 1])
                    / (a + b) ** 4)
    assert n4["err"] == pytest.approx(exact, rel=1e-9)
    # a one-line family IS its site: same fraction, letters and integral
    fams = {f["family"]: f for f in q["families"]}
    assert fams["BO3"]["letters"] == "A" and fams["BO4"]["letters"] == "B"
    for fam, row in (("BO3", q["rows"][0]), ("BO4", q["rows"][1])):
        assert fams[fam]["fraction_pct"] == pytest.approx(row["fraction_pct"], rel=1e-12)
        assert fams[fam]["integral"] == pytest.approx(abs(row["integral"]), rel=1e-12)
    # site rows are byte-identical whatever the basis
    q_ind = quantify(rec, (40.0, -20.0))
    assert q_ind["rows"] == q["rows"]
    assert q_ind["family_basis"] == "independent"
    n4_ind = next(r for r in q_ind["ratios"] if r["name"] == "N4")
    indep = np.sqrt((b * b * sa ** 2 + a * a * sb ** 2) / (a + b) ** 4)
    assert n4_ind["err"] == pytest.approx(indep, rel=1e-9)
    assert n4_ind["err"] != pytest.approx(n4["err"], rel=1e-6)
    assert n4_ind["err"] < n4["err"]                    # ρ < 0: ratio known worse
    assert "independent" in q_ind["family_note"]
    # ...while the SUM of the anticorrelated pair beats its quadrature: tag
    # both as one family beside a third, uncorrelated line
    c, sc = 50.0, 3.0
    cov3 = np.zeros((3, 3))
    cov3[:2, :2] = cov
    cov3[2, 2] = sc ** 2
    rec3 = _b_recipe([_gl(15.0, 6.0, a, "A", "BO3", stderr=sa),
                      _gl(12.0, 6.0, b, "B", "BO3", stderr=sb),
                      _gl(0.0, 6.0, c, "C", "BO4", stderr=sc)])
    u = correlated_values([a, b, c], cov3)
    q3 = quantify(rec3, (40.0, -20.0), uvars=dict(enumerate(u)))
    q3_ind = quantify(rec3, (40.0, -20.0))
    bo3 = next(f for f in q3["families"] if f["family"] == "BO3")
    bo3_ind = next(f for f in q3_ind["families"] if f["family"] == "BO3")
    assert bo3["letters"] == "A+B" and bo3["sites"] == [0, 1]
    assert bo3["fraction_err_pct"] < bo3_ind["fraction_err_pct"]


def test_montecarlo_samples_give_the_trial_spread_and_keep_best_values():
    """Perfectly anticorrelated per-trial integrals (I_A + I_B constant):
    the N4 error IS np.std of the per-trial ratio, a family holding both
    lines has zero spread (smaller than any quadrature), the reported value
    stays the best-fit number (not the trial mean), basis 'montecarlo'."""
    rng = np.random.default_rng(4)
    d = rng.normal(0.0, 2.0, 200)
    samples = np.column_stack([60.0 + d, 40.0 - d])
    integrals, amps = [60.0, 40.0], [10.0, 5.0]
    out = families.summarize(integrals, amps, ["BO3", "BO4"], "11B",
                             sigma=[2.0, 2.0], samples=samples)
    assert out["basis"] == "montecarlo" and "200 trials" in out["note"]
    n4 = next(r for r in out["ratios"] if r["name"] == "N4")
    assert n4["value"] == pytest.approx(0.4)             # best fit, not the mean
    assert n4["err"] == pytest.approx(float(np.std(samples[:, 1] / samples.sum(axis=1))))
    bo4 = next(f for f in out["families"] if f["family"] == "BO4")
    assert bo4["fraction_pct"] == pytest.approx(40.0)
    assert bo4["fraction_err_pct"] == pytest.approx(
        float(np.std(100.0 * samples[:, 1] / samples.sum(axis=1))))
    # both lines in ONE family beside a fixed third line: zero spread, while
    # the independent quadrature of the same column stds is not
    samples3 = np.column_stack([60.0 + d, 40.0 - d, np.full(d.size, 50.0)])
    both = families.summarize([60.0, 40.0, 50.0], [1, 1, 1], ["BO3", "BO3", "BO4"],
                              "11B", samples=samples3)
    ind = families.summarize([60.0, 40.0, 50.0], [1, 1, 1], ["BO3", "BO3", "BO4"],
                             "11B", sigma=[float(np.std(d))] * 2 + [0.0])
    f_mc = next(f for f in both["families"] if f["family"] == "BO3")
    f_ind = next(f for f in ind["families"] if f["family"] == "BO3")
    assert f_mc["fraction_err_pct"] == pytest.approx(0.0, abs=1e-9)
    assert f_ind["fraction_err_pct"] > 0.5
    assert ind["basis"] == "independent"
    # a sample block of the wrong width is ignored, not trusted
    bad = families.summarize(integrals, amps, ["BO3", "BO4"], "11B",
                             sigma=[2.0, 2.0], samples=np.ones((5, 3)))
    assert bad["basis"] == "independent"


def test_ratio_needs_two_denominator_families_untagged_listed_and_weighted_ratios():
    # only BO4 tagged: N4 is not a number, the note names the missing family
    one = families.summarize([60.0, 40.0], [1, 1], ["", "BO4"], "11B")
    n4 = one["ratios"][0]
    assert n4["name"] == "N4" and not n4["defined"] and n4["value"] is None
    assert "BO3" in n4["note"] and one["untagged"] == [0]
    assert [f["family"] for f in one["families"]] == ["BO4"]
    assert one["families"][0]["fraction_pct"] == pytest.approx(40.0)
    # a zero-amplitude tagged family counts as present: N4 = 0.0 is defined
    zero = families.summarize([60.0, 0.0], [1.0, 0.0], ["BO3", "BO4"], "11B")
    n4z = next(r for r in zero["ratios"] if r["name"] == "N4")
    assert n4z["defined"] and n4z["value"] == 0.0
    # a third, untagged line is in the total but in no family
    three = families.summarize([60.0, 20.0, 20.0], [1, 1, 1], ["BO3", "BO4", ""],
                               "11B")
    assert three["untagged"] == [2]
    assert sum(f["fraction_pct"] for f in three["families"]) == pytest.approx(80.0)
    assert "untagged: C" in three["note"]
    # weighted ratios: ⟨CN⟩ Al and ⟨n⟩ Si
    al = families.summarize([50.0, 50.0], [1, 1], ["Al(IV)", "Al(VI)"], "27Al")
    cn = next(r for r in al["ratios"] if r["name"] == "⟨CN⟩ Al")
    assert cn["defined"] and cn["value"] == pytest.approx(5.0)
    assert cn["fmt"].format(cn["value"]) == "5.00"
    si = families.summarize([30.0, 30.0], [1, 1], ["Q2", "Q3"], "29Si")
    n = next(r for r in si["ratios"] if r["name"] == "⟨n⟩ Si")
    assert n["defined"] and n["value"] == pytest.approx(2.5)
    # untagged recipe: nothing to report, basis irrelevant, no note
    none = families.summarize([1.0, 2.0], [1, 1], ["", ""], "11B")
    assert none["families"] == [] and none["ratios"] == [] and none["untagged"] == [0, 1]
    # preset order first, then user tags as first seen
    mixed = families.summarize([1, 1, 1, 1], [1, 1, 1, 1],
                               ["mine", "BO4", "BO3", "mine"], "11B")
    assert [f["family"] for f in mixed["families"]] == ["BO3", "BO4", "mine"]
    assert mixed["families"][2]["letters"] == "A+D"


def test_presets_ratios_and_label_guesses_are_consistent():
    for nuc, defs in families.RATIOS.items():
        presets = set(families.PRESETS[nuc])
        for rd in defs:
            assert set(rd.numerator) <= presets and set(rd.denominator) <= presets, rd
    assert families.presets_for("11B") == ["BO3", "BO4"]
    assert families.presets_for(" 27Al ") == ["Al(IV)", "Al(V)", "Al(VI)"]
    assert families.presets_for("7Li") == [] and families.ratios_for(None) == ()
    assert "bonded" in families.presets_for("31P") and "Q3" in families.presets_for("31P")
    assert families.guess_from_label("11B", "B[4] (BO4)") == "BO4"
    assert families.guess_from_label("11B", "B[3]") == "BO3"
    assert families.guess_from_label("27Al", "AlIV") == "Al(IV)"
    assert families.guess_from_label("27Al", "Al(VI) octahedral") == "Al(VI)"
    assert families.guess_from_label("29Si", "Q3 site") == "Q3"
    assert families.guess_from_label("31P", "P2 chain") == ""
    assert families.guess_from_label("11B", "Czjzek-0") == ""
    assert families.guess_from_label("17O", "BO") == ""      # no aliases compiled
    assert families.DERIVED_PARAMS >= families.GROUP_PARAMS
    assert "population_pct" in families.DERIVED_PARAMS
    assert "population_pct" not in families.GROUP_PARAMS


def test_site_integrals_matches_quantify_and_windows_are_compared_with_tolerance():
    rec = _b_recipe([_gl(15.0, 8.0, 1.0, "A", "BO3"), _gl(0.0, 4.0, 1.0, "B", "BO4")])
    ints, win = site_integrals(rec, (40.0, -20.0))
    q = quantify(rec, (40.0, -20.0))
    assert win == (40.0, -20.0)
    assert list(ints) == pytest.approx([r["integral"] for r in q["rows"]], rel=1e-12)
    # equal peak heights, widths 8 vs 4 ppm: AREAS 2:1, never 50/50
    fams = {f["family"]: f["fraction_pct"] for f in q["families"]}
    assert fams["BO3"] == pytest.approx(66.67, abs=0.3)
    assert fams["BO4"] == pytest.approx(33.33, abs=0.3)
    n4 = next(r for r in q["ratios"] if r["name"] == "N4")
    assert n4["value"] == pytest.approx(1 / 3, abs=0.003)
    # the window resolution order: explicit > recipe.fit_window_ppm > full axis
    rec.fit_window_ppm = (30.0, -10.0)
    assert site_integrals(rec)[1] == (30.0, -10.0)
    assert site_integrals(rec, (-5.0, 25.0))[1] == (25.0, -5.0)   # normalised (hi, lo)
