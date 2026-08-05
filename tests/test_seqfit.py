"""Sequential (forward-backward) series fitting: warm-start each spectrum from
its fitted neighbour, sweep back to smooth, track evolving parameters."""
import numpy as np
import pytest

from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine, seqfit


def _spec(x, pos, amp, seed):
    tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                sites=[SiteModel(model="gauss_lor", label="A", params={
                    "isotropic_chemical_shift_ppm": Param(pos),
                    "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                    "gl": Param(1.0, vary=False)})])
    _, m, _ = engine.simulate(tr, exp_ppm=x)
    return m + np.random.default_rng(seed).normal(0.0, 1.0, x.size)


def _start(sample):
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                  sample=sample, sites=[SiteModel(model="gauss_lor", label="A",
                      params={
                          "isotropic_chemical_shift_ppm": Param(10.0, min=0, max=30),
                          "shift_fwhm_ppm": Param(5.0, min=0.1),
                          "amplitude": Param(80.0, min=0),
                          "gl": Param(1.0, vary=False)})])


def _entries():
    # a smooth series: position marches 12 -> 18 ppm across 6 spectra
    x = np.linspace(-20, 60, 700)
    positions = np.linspace(12.0, 18.0, 6)
    amps = [100, 90, 80, 85, 95, 110]
    return [(_start(f"g{k}"), x, _spec(x, p, a, k), (-10.0, 40.0))
            for k, (p, a) in enumerate(zip(positions, amps))]


def test_seed_from_copies_values_within_bounds():
    a, b = _start("a"), _start("b")
    b.sites[0].params["isotropic_chemical_shift_ppm"].value = 25.0
    seqfit.seed_from(a, b, ("isotropic_chemical_shift_ppm",))
    assert a.sites[0].params["isotropic_chemical_shift_ppm"].value == 25.0
    # amplitude not in the propagate set -> unchanged
    assert a.sites[0].params["amplitude"].value == 80.0


def test_sequential_recovers_marching_position():
    res = seqfit.run_sequential(_entries(), passes=2, smooth=0)
    pos = [r.sites[0].params["isotropic_chemical_shift_ppm"].value
           for r in res.recipes]
    assert pos[0] == pytest.approx(12.0, abs=0.6)
    assert pos[-1] == pytest.approx(18.0, abs=0.6)
    assert all(pos[i] < pos[i + 1] + 0.5 for i in range(len(pos) - 1))  # monotone-ish


def test_more_passes_and_smoothing_do_not_diverge():
    res = seqfit.run_sequential(_entries(), passes=4, smooth=3)
    assert res.passes == 4
    assert len(res.history) >= 1
    assert np.isfinite(res.history[-1]["mean"])
    # final mean RMSD is comparable to or better than the first pass
    assert res.history[-1]["mean"] <= res.history[0]["mean"] * 1.5


def test_progress_and_stop_callbacks():
    seen = []
    seqfit.run_sequential(_entries(), passes=2,
                          progress=lambda p, k, r: seen.append((p, k, r)))
    assert seen and all(len(t) == 3 for t in seen)

    calls = {"n": 0}

    def stop():
        calls["n"] += 1
        return calls["n"] >= 2

    res = seqfit.run_sequential(_entries(), passes=4, should_stop=stop)
    assert res is not None                       # returns partial result


def test_needs_two_spectra():
    with pytest.raises(ValueError):
        seqfit.run_sequential(_entries()[:1])


def test_direction_alternates():
    res = seqfit.run_sequential(_entries(), passes=2)
    assert res.history[0]["direction"] == "→"
    assert res.history[1]["direction"] == "←"
