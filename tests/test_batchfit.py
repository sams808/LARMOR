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
    # covariance is available with no extra work (it comes off the fit)
    rows = batchfit.error_table(res, method="covariance")
    assert res.error_method == "covariance"
    amp_rows = [r for r in rows if r["param"] == "amplitude"]
    assert amp_rows and all(r["error_method"] == "covariance" for r in amp_rows)
    # every fitted amplitude has an error and a % error
    assert all(r["stderr"] is not None for r in amp_rows)


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


def test_batch_fit_accepts_tol_and_still_runs():
    # a loose threshold still returns a result; released position tracks the data
    res = batchfit.batch_fit(_entries(), tol=1.0,
                             release=("isotropic_chemical_shift_ppm",),
                             release_frac=0.2)
    posA = [r.sites[0].params["isotropic_chemical_shift_ppm"].value
            for r in res.recipes]
    assert posA[0] == pytest.approx(15.0, abs=0.8)
