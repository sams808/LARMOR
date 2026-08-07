"""Regressions for the batch fit (larmor.batchfit): one shared model, many 1D
spectra, amplitudes free per spectrum, with an optional 'release' stage.
"""
import numpy as np
import pytest

from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine, batchfit


def _spec(x, pos_shift, amps, seed):
    tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                sites=[
                    SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(15.0 + pos_shift),
                        "shift_fwhm_ppm": Param(6.0),
                        "amplitude": Param(amps[0]),
                        "gl": Param(1.0, vary=False)}),
                    SiteModel(model="gauss_lor", label="B", params={
                        "isotropic_chemical_shift_ppm": Param(2.0),
                        "shift_fwhm_ppm": Param(3.0),
                        "amplitude": Param(amps[1]),
                        "gl": Param(1.0, vary=False)})])
    _, m, _ = engine.simulate(tr, exp_ppm=x)
    return m + np.random.default_rng(seed).normal(0.0, 1.5, x.size)


def _start(sample):
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                  sample=sample, sites=[
                      SiteModel(model="gauss_lor", label="A", params={
                          "isotropic_chemical_shift_ppm": Param(14.0, min=0, max=30),
                          "shift_fwhm_ppm": Param(5.0, min=0.1),
                          "amplitude": Param(80.0, min=0),
                          "gl": Param(1.0, vary=False)}),
                      SiteModel(model="gauss_lor", label="B", params={
                          "isotropic_chemical_shift_ppm": Param(1.0, min=-10, max=15),
                          "shift_fwhm_ppm": Param(3.5, min=0.1),
                          "amplitude": Param(50.0, min=0),
                          "gl": Param(1.0, vary=False)})])


def _entries():
    x = np.linspace(-20, 60, 800)
    shifts, amps = [0.0, 0.3, -0.3], [[100, 60], [80, 90], [120, 40]]
    return [(_start(f"g{k}"), x, _spec(x, sh, am, k), (-10.0, 40.0))
            for k, (sh, am) in enumerate(zip(shifts, amps))]


def test_all_but_amplitude_excludes_amplitude():
    r = _start("s")
    names = batchfit.all_but_amplitude([r])
    assert "amplitude" not in names
    assert "isotropic_chemical_shift_ppm" in names
    assert "shift_fwhm_ppm" in names


def test_unreleased_shape_stays_fixed_at_recipe_value():
    # nothing released → δiso and width are HELD at the recipe values (14.0 / 5.0)
    # for every spectrum; only the amplitudes are fitted (and differ)
    res = batchfit.batch_fit(_entries())
    posA = [r.sites[0].params["isotropic_chemical_shift_ppm"].value
            for r in res.recipes]
    fwhm = [r.sites[0].params["shift_fwhm_ppm"].value for r in res.recipes]
    ampA = [r.sites[0].params["amplitude"].value for r in res.recipes]
    assert all(p == pytest.approx(14.0) for p in posA)   # NOT fitted — held
    assert all(w == pytest.approx(5.0) for w in fwhm)     # NOT fitted — held
    assert ampA[0] > ampA[1] and ampA[2] > ampA[0]        # amplitudes free/distinct


def test_release_lets_selected_param_move_per_spectrum():
    res = batchfit.batch_fit(
        _entries(), release=("isotropic_chemical_shift_ppm",), release_frac=0.2)
    posA = [r.sites[0].params["isotropic_chemical_shift_ppm"].value
            for r in res.recipes]
    fwhm = [r.sites[0].params["shift_fwhm_ppm"].value for r in res.recipes]
    assert res.released == ("isotropic_chemical_shift_ppm",)
    # released → moves per spectrum, tracking the injected ±0.3 ppm shifts
    assert max(posA) - min(posA) > 0.3
    assert posA[1] > posA[0] > posA[2]
    # width was NOT released → still fixed at the recipe value everywhere
    assert all(w == pytest.approx(5.0) for w in fwhm)


def test_needs_at_least_two_spectra():
    with pytest.raises(ValueError):
        batchfit.batch_fit(_entries()[:1])


def test_shared_table_has_shared_and_per_spectrum_rows():
    res = batchfit.batch_fit(_entries())
    rows = batchfit.shared_table(res)
    assert any(r["scope"] == "shared" and r["param"] == "shift_fwhm_ppm"
               for r in rows)
    assert any(r["param"] == "amplitude" and r["scope"] != "shared" for r in rows)


def test_free_amplitudes_overrides_recipe_locks():
    r = _start("s")
    r.sites[0].params["amplitude"].vary = False        # locked in the recipe
    r.sites[0].params["amplitude"].min = 50.0          # and bounded above zero
    batchfit.free_amplitudes([r])
    amp = r.sites[0].params["amplitude"]
    assert amp.vary is True and amp.min == 0.0          # freed, may reach zero
    # a linked amplitude is left alone
    r.sites[1].params["amplitude"].expr = "0.5 * s0.amplitude"
    batchfit.free_amplitudes([r])
    assert r.sites[1].params["amplitude"].expr == "0.5 * s0.amplitude"


def test_batch_fit_fits_a_locked_amplitude_per_spectrum():
    # the reported bug: a recipe with a LOCKED amplitude must still fit that
    # amplitude per spectrum (and let it adapt), not freeze it
    entries = _entries()
    for rec, *_ in entries:
        rec.sites[0].params["amplitude"].vary = False  # lock site A everywhere
    res = batchfit.batch_fit(entries)
    ampA = [rr.sites[0].params["amplitude"].value for rr in res.recipes]
    assert max(ampA) - min(ampA) > 1.0                 # it adapted per spectrum
    assert all(a >= -1e-6 for a in ampA)               # stayed non-negative


def test_completion_threshold_maps_pct_to_ftol():
    from larmor.fit import ftol_from_pct, _tol_kws
    assert ftol_from_pct(0) is None            # 0 / off -> solver default
    assert ftol_from_pct(None) is None
    assert ftol_from_pct(5) == pytest.approx(0.1)     # 2 * 5/100
    assert _tol_kws(0) == {}
    assert _tol_kws(5)["ftol"] == pytest.approx(0.1)


def _data_for(entries):
    return [(ppm, amp, win) for _rec, ppm, amp, win in entries]


def test_error_analysis_covariance_is_default_and_exportable():
    entries = _entries()
    res = batchfit.batch_fit(entries)
    # this fixture's amplitude-only fit is well-conditioned even without the
    # errorbar-rescue retry, so the raw post-fit stderr already works here --
    # export still goes through batch_error_analysis explicitly (covariance
    # is no longer a free snapshot in general, see the test below)
    batchfit.batch_error_analysis(res, _data_for(entries), method="covariance")
    rows = batchfit.error_table(res, method="covariance")
    assert res.error_method == "covariance"
    amp_rows = [r for r in rows if r["param"] == "amplitude"]
    assert amp_rows and all(r["error_method"] == "covariance" for r in amp_rows)
    # every fitted amplitude has an error and a % error
    assert all(r["stderr"] is not None for r in amp_rows)


def _degenerate_batch_entries():
    """Two IDENTICAL-shape overlapping sites per spectrum: their amplitudes
    are perfectly correlated (only the SUM is determined by the data), so the
    covariance is reliably singular on the fast (no-retry) pass -- a
    controlled trigger, not left to chance like a merely-overlapping fit."""
    x = np.linspace(-30, 30, 500)
    entries = []
    for k in range(3):
        truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(0.0),
                "shift_fwhm_ppm": Param(10.0), "amplitude": Param(80.0),
                "gl": Param(1.0, vary=False)})])
        _, y, _ = engine.simulate(truth, exp_ppm=x)
        data = y + np.random.default_rng(k).normal(0, 0.3, x.size)
        rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=f"g{k}",
                    sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(0.0, vary=False),
                "shift_fwhm_ppm": Param(10.0, vary=False),
                "amplitude": Param(40.0, min=0), "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="B", params={
                "isotropic_chemical_shift_ppm": Param(0.0, vary=False),
                "shift_fwhm_ppm": Param(10.0, vary=False),
                "amplitude": Param(40.0, min=0), "gl": Param(1.0, vary=False)}),
        ])
        entries.append((rec, x, data, None))
    return entries


def test_error_analysis_covariance_refits_when_the_fast_fit_had_none():
    """Regression: batch_fit's initial pass uses compute_errorbars=False for
    speed (see its docstring) -- for a genuinely degenerate released fit this
    leaves every parameter's stderr as None right after the fit (verified: the
    exact mechanism test_engine_fit.py's degenerate-covariance test exercises
    at the single-fit level). Requesting covariance errors must not just
    silently report that emptiness -- it must refit with
    compute_errorbars=True and get the TOTAL (identifiable) amplitude's real
    error, not leave every row blank."""
    entries = _degenerate_batch_entries()
    res = batchfit.batch_fit(entries)
    free = [(i, pn) for rec in res.recipes for i, pn in batchfit._free_params(rec)]
    n_none_raw = sum(1 for rec in res.recipes for i, pn in batchfit._free_params(rec)
                     if rec.sites[i].params[pn].stderr is None)
    assert n_none_raw == len(free)          # confirms the fixture DOES trigger it

    batchfit.batch_error_analysis(res, _data_for(entries), method="covariance")
    rows = batchfit.error_table(res, method="covariance")
    assert any(r["stderr"] is not None for r in rows)


def test_error_analysis_montecarlo_writes_errors_and_detail():
    entries = _entries()
    res = batchfit.batch_fit(entries)
    data = _data_for(entries)
    seen = []
    batchfit.batch_error_analysis(res, data, method="montecarlo", n_trials=40,
                                  seed=1, progress=lambda *a: seen.append(a))
    assert res.error_method == "montecarlo"
    assert "montecarlo" in res.error_detail
    assert "covariance" in res.error_detail        # snapshot preserved
    # MC σ written into the recipes and non-negative
    for rec in res.recipes:
        se = rec.sites[0].params["amplitude"].stderr
        assert se is not None and se >= 0
    rows = batchfit.error_table(res, method="montecarlo")
    assert any(r["param"] == "amplitude" and r["sigma_pct"] is not None
               for r in rows)
    assert seen                                    # progress was reported


def test_error_analysis_profile_gives_confidence_intervals():
    entries = _entries()
    res = batchfit.batch_fit(entries)
    data = _data_for(entries)
    batchfit.batch_error_analysis(res, data, method="profile", n_points=9,
                                  span=3.0)
    assert res.error_method == "profile"
    rows = batchfit.error_table(res, method="profile")
    amp = [r for r in rows if r["param"] == "amplitude"]
    # at least one amplitude has a bracketed 1σ interval around its value
    assert any(r["ci68_lo"] is not None and r["ci68_hi"] is not None
               and r["ci68_lo"] <= r["value"] <= r["ci68_hi"] for r in amp)


def test_switching_export_method_keeps_covariance():
    entries = _entries()
    res = batchfit.batch_fit(entries)
    cov = batchfit.error_table(res, method="covariance")
    cov_amp = next(r["stderr"] for r in cov if r["param"] == "amplitude")
    batchfit.batch_error_analysis(res, _data_for(entries), method="montecarlo",
                                  n_trials=40, seed=2)
    # after MC overwrote Param.stderr, exporting covariance still returns the
    # original least-squares errors (the snapshot), not the MC σ
    cov2 = batchfit.error_table(res, method="covariance")
    cov2_amp = next(r["stderr"] for r in cov2 if r["param"] == "amplitude")
    assert cov2_amp == pytest.approx(cov_amp)


def _czjzek_entries():
    """A batch of synthetic 27Al Czjzek spectra: a kernel-based model, unlike the
    Gauss/Lorentz fixture above -- exercises the batch pipeline end-to-end on the
    lineshape family the user actually runs next (27Al)."""
    x = np.linspace(-40, 120, 700)

    def start(sample):
        return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sample=sample,
                      sites=[SiteModel(model="czjzek", label="Al", params={
                          "isotropic_chemical_shift_ppm": Param(60.0, min=40, max=80),
                          "sigma_Cq_MHz": Param(3.0, min=0.05),
                          "shift_fwhm_ppm": Param(8.0, min=0.1),
                          "line_fwhm_ppm": Param(1.0, min=0.0),
                          "amplitude": Param(100.0, min=0)})])

    def spec(amp, seed):
        tr = start("t")
        tr.sites[0].params["amplitude"].value = amp
        kx, m, _ = engine.simulate(tr)
        y = np.interp(x, np.sort(kx), m[np.argsort(kx)])
        return y + np.random.default_rng(seed).normal(0.0, y.max() * 0.01, x.size)

    return [(start(f"g{k}"), x, spec(amp, k), (120.0, -40.0))
            for k, amp in enumerate([100.0, 80.0, 120.0])]


def test_batch_fit_holds_czjzek_shape_and_frees_amplitude():
    res = batchfit.batch_fit(_czjzek_entries())
    sigma = [r.sites[0].params["sigma_Cq_MHz"].value for r in res.recipes]
    assert all(s == pytest.approx(3.0) for s in sigma)     # shape held fixed
    amps = [r.sites[0].params["amplitude"].value for r in res.recipes]
    assert amps[2] > amps[0] > amps[1]                     # amplitudes track truth


def test_batch_error_analysis_montecarlo_works_on_czjzek():
    # the kernel-model grid-mismatch bug this generalisation surfaced
    # (autofit.monte_carlo_errors) must not resurface here
    entries = _czjzek_entries()
    res = batchfit.batch_fit(entries)
    batchfit.batch_error_analysis(
        res, [(a, b, c) for _, a, b, c in entries], method="montecarlo",
        n_trials=15, seed=1)
    rows = batchfit.error_table(res, method="montecarlo")
    amp_rows = [r for r in rows if r["param"] == "amplitude"]
    assert amp_rows and all(r["stderr"] is not None and r["stderr"] >= 0
                            for r in amp_rows)


def test_snapshot_covariance_normalizes_nan_stderr_to_none():
    """An ill-conditioned covariance makes lmfit report stderr as NaN, not
    None. NaN is truthy, so an unguarded `if se` would treat it as a real
    error -- _snapshot_covariance must normalize it to None."""
    res = batchfit.batch_fit(_entries())
    res.recipes[0].sites[0].params["amplitude"].stderr = float("nan")
    detail = batchfit._snapshot_covariance(res)
    pe = detail[0][(0, "amplitude")]
    assert pe.stderr is None and pe.pct is None


def test_batch_fit_accepts_tol_and_still_runs():
    # a loose threshold still returns a result; released position tracks the data
    res = batchfit.batch_fit(_entries(), tol=1.0,
                             release=("isotropic_chemical_shift_ppm",),
                             release_frac=0.2)
    posA = [r.sites[0].params["isotropic_chemical_shift_ppm"].value
            for r in res.recipes]
    assert posA[0] == pytest.approx(15.0, abs=0.8)


def test_batch_fit_is_independent_per_spectrum_not_a_joint_optimisation():
    """Each spectrum is fit on its own (larmor.fit.fit), not as one combined
    multi-spectrum problem -- verified two ways: (1) the iter_cb / progress
    callback keeps firing continuously across MULTIPLE spectra's fits (the
    existing dialog code relies on this to show one running counter, unchanged
    since before this refactor), and (2) each fitted recipe's rmsd is
    consistent with ONLY its own data (spot-checked against a direct
    single-spectrum larmor.fit.fit call on the same data)."""
    from larmor import fit as fitmod

    entries = _entries()
    calls = []
    res = batchfit.batch_fit(entries, iter_cb=lambda *a, **k: calls.append(1))
    assert len(calls) > len(entries)          # more than one call per spectrum
    assert res.lmfit_result is None           # no single combined result anymore

    # spot check spectrum 1 in isolation: same recipe/data/window as inside
    # batch_fit (shape held fixed, only amplitude free) must give the same rmsd
    rec, ppm, amp, window = entries[1]
    solo = Recipe.from_dict(rec.to_dict())
    for s in solo.sites:
        for pn, p in s.params.items():
            p.vary = (pn == "amplitude")
    fitmod.fit(solo, ppm, amp, window_ppm=window)
    assert res.rmsd[1] == pytest.approx(solo.fit_rmsd, rel=1e-6)


def test_is_zeroed_out_requires_the_full_lock_signature():
    assert batchfit.is_zeroed_out(Param(0.0, vary=False, min=0.0, max=0.0))
    assert not batchfit.is_zeroed_out(Param(0.0, vary=True, min=0.0, max=0.0))
    assert not batchfit.is_zeroed_out(Param(0.0, vary=False, min=0.0, max=None))
    assert not batchfit.is_zeroed_out(Param(1e-9, vary=False, min=0.0, max=0.0))
    assert not batchfit.is_zeroed_out(None)


def test_free_amplitudes_preserves_an_explicit_zero_lock():
    r = _start("s")
    r.sites[1].params["amplitude"] = Param(0.0, vary=False, min=0.0, max=0.0)
    batchfit.free_amplitudes([r])
    amp = r.sites[1].params["amplitude"]
    assert amp.vary is False and amp.value == 0.0     # left locked...
    other = r.sites[0].params["amplitude"]
    assert other.vary is True and other.min == 0.0    # ...unlike a normal site


def test_batch_fit_respects_a_per_spectrum_zero_exclusion():
    """"Exclude component": one spectrum's site B is locked to zero BEFORE
    batch_fit runs -- it must stay exactly zero after fitting, while every
    OTHER spectrum's site B fits normally (the exclusion is per-spectrum,
    not global)."""
    entries = _entries()
    entries[1][0].sites[1].params["amplitude"] = Param(
        0.0, vary=False, min=0.0, max=0.0)
    res = batchfit.batch_fit(entries)
    ampB = [r.sites[1].params["amplitude"].value for r in res.recipes]
    assert ampB[1] == 0.0
    assert not res.recipes[1].sites[1].params["amplitude"].vary
    assert ampB[0] > 1.0 and ampB[2] > 1.0            # the other two fit normally


def test_excluded_site_params_stay_fixed_even_when_that_param_is_released():
    """Regression: excluding site B on one spectrum while ALSO releasing
    'isotropic_chemical_shift_ppm'/'shift_fwhm_ppm' globally (the real-world
    combination that triggered this bug — "release per spectrum" ticked
    alongside an excluded component) must NOT leave site B's position/width
    varying on the excluded spectrum. A zero-amplitude site contributes
    nothing to the residual, so those would be dead-gradient free parameters:
    pure optimiser waste at best, degraded convergence for the real
    parameters at worst — exactly the "worse RMSD, slower fit" symptom
    reported against a live 8-line/17-spectrum batch fit."""
    entries = _entries()
    entries[1][0].sites[1].params["amplitude"] = Param(
        0.0, vary=False, min=0.0, max=0.0)
    res = batchfit.batch_fit(
        entries, release=("isotropic_chemical_shift_ppm", "shift_fwhm_ppm"))

    excluded_site = res.recipes[1].sites[1]
    assert excluded_site.params["amplitude"].value == 0.0
    assert not excluded_site.params["isotropic_chemical_shift_ppm"].vary
    assert not excluded_site.params["shift_fwhm_ppm"].vary

    # the same site on the OTHER (non-excluded) spectra is still released
    for k in (0, 2):
        site = res.recipes[k].sites[1]
        assert site.params["isotropic_chemical_shift_ppm"].vary
        assert site.params["shift_fwhm_ppm"].vary
    # and site A (never excluded, anywhere) is released on every spectrum,
    # including the one where site B was excluded
    for r in res.recipes:
        assert r.sites[0].params["isotropic_chemical_shift_ppm"].vary
        assert r.sites[0].params["shift_fwhm_ppm"].vary

    # quality check: fitting fewer (real) free parameters for the excluded
    # spectrum should not produce a worse fit than fitting it with the SAME
    # exclusion but WITHOUT the (harmless-in-theory, harmful-in-practice)
    # dead columns -- i.e. RMSD should match a control run that never had
    # site B in the model to begin with.
    control = _entries()[1]
    control[0].sites.pop(1)                 # site B never existed at all
    control_res = batchfit.batch_fit(
        [control, _entries()[0], _entries()[2]],
        release=("isotropic_chemical_shift_ppm", "shift_fwhm_ppm"))
    assert res.rmsd[1] == pytest.approx(control_res.rmsd[0], abs=1e-6)


def test_shared_and_error_table_omit_excluded_sites_and_add_population_pct():
    entries = _entries()
    entries[1][0].sites[1].params["amplitude"] = Param(
        0.0, vary=False, min=0.0, max=0.0)
    res = batchfit.batch_fit(entries)

    shared_rows = batchfit.shared_table(res)
    g1_sites = {r["site"] for r in shared_rows if r["scope"] == "g1"}
    g0_sites = {r["site"] for r in shared_rows if r["scope"] == "g0"}
    assert "s1" not in g1_sites                       # excluded for g1...
    assert "s1" in g0_sites                            # ...but present for g0

    pop = [r for r in shared_rows if r["param"] == "population_pct"]
    assert pop                                          # rows were added
    g0_pop = {r["site"]: r["value"] for r in pop if r["scope"] == "g0"}
    assert set(g0_pop) == {"s0", "s1"}
    assert sum(g0_pop.values()) == pytest.approx(100.0, abs=0.5)
    g1_pop = {r["site"]: r["value"] for r in pop if r["scope"] == "g1"}
    assert set(g1_pop) == {"s0"}                        # excluded site has no population row
    assert g1_pop["s0"] == pytest.approx(100.0, abs=0.5)  # the only site left = 100%

    batchfit.batch_error_analysis(res, _data_for(entries), method="covariance")
    err_rows = batchfit.error_table(res, method="covariance")
    err_pop = [r for r in err_rows if r["param"] == "population_pct"]
    assert err_pop and all(r["error_method"] == "covariance" for r in err_pop)
    assert not any(r["scope"] == "g1" and r["site"] == "s1" for r in err_pop)
    # error-method consistency: population stderr differs from a naive
    # "whatever stderr happened to be on the recipe" snapshot once a
    # DIFFERENT method is requested
    batchfit.batch_error_analysis(res, _data_for(entries), method="montecarlo",
                                  n_trials=30, seed=3)
    mc_rows = batchfit.error_table(res, method="montecarlo")
    mc_pop = {(r["scope"], r["site"]): r["stderr"]
             for r in mc_rows if r["param"] == "population_pct"}
    assert all(v is not None and v >= 0 for v in mc_pop.values())


def test_batch_fit_should_stop_gives_every_entry_a_result():
    """A Stop mid-batch must not misalign indices: every entry (fit or not)
    gets a recipe/rmsd/label/per_dataset entry, in order, so the dialog's
    per-cell UI (indexed by k) never goes out of sync."""
    entries = _entries()      # 3 spectra
    seen = {"n": 0}

    def stop_after_first():
        # true once the FIRST spectrum's fit has made real progress, so
        # spectra 1 and 2 are reached only after stop is already requested
        return seen["n"] > 3

    def cb(*a, **k):
        seen["n"] += 1

    res = batchfit.batch_fit(entries, iter_cb=cb, should_stop=stop_after_first)
    assert len(res.recipes) == len(entries)
    assert len(res.rmsd) == len(entries)
    assert len(res.labels) == len(entries)
    assert len(res.per_dataset) == len(entries)
    for pd in res.per_dataset:
        assert len(pd["x"]) == len(pd["y_fit"]) > 0
    # a spectrum reached after the stop keeps its PRE-fit (unmodified) recipe
    unfit_amp = _start("g2").sites[0].params["amplitude"].value
    assert res.recipes[2].sites[0].params["amplitude"].value == \
        pytest.approx(unfit_amp)
