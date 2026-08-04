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


def test_shared_stage_ties_shape_frees_amplitude():
    res = batchfit.batch_fit(_entries())
    posA = [r.sites[0].params["isotropic_chemical_shift_ppm"].value
            for r in res.recipes]
    ampA = [r.sites[0].params["amplitude"].value for r in res.recipes]
    # positions identical across spectra (shared), ~15
    assert max(posA) - min(posA) < 1e-3
    assert posA[0] == pytest.approx(15.0, abs=0.3)
    # amplitudes free and distinct (100 / 80 / 120)
    assert ampA[0] > ampA[1] and ampA[2] > ampA[0]


def test_release_lets_selected_param_drift_per_spectrum():
    res = batchfit.batch_fit(
        _entries(), release=("isotropic_chemical_shift_ppm",), release_frac=0.1)
    posA = [r.sites[0].params["isotropic_chemical_shift_ppm"].value
            for r in res.recipes]
    assert res.released == ("isotropic_chemical_shift_ppm",)
    # now distinct, tracking the injected ±0.3 ppm shifts
    assert max(posA) - min(posA) > 0.3
    assert posA[1] > posA[0] > posA[2]
    # widths stayed shared
    fwhm = [r.sites[0].params["shift_fwhm_ppm"].value for r in res.recipes]
    assert max(fwhm) - min(fwhm) < 1e-3


def test_needs_at_least_two_spectra():
    with pytest.raises(ValueError):
        batchfit.batch_fit(_entries()[:1])


def test_shared_table_has_shared_and_per_spectrum_rows():
    res = batchfit.batch_fit(_entries())
    rows = batchfit.shared_table(res)
    assert any(r["scope"] == "shared" and r["param"] == "shift_fwhm_ppm"
               for r in rows)
    assert any(r["param"] == "amplitude" and r["scope"] != "shared" for r in rows)
