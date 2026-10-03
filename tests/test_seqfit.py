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


def test_stop_is_responsive_across_many_passes():
    # Stop must short-circuit the WHOLE sweep, not run all the requested passes
    # (a 16-pass auto sweep was effectively unstoppable before)
    res = seqfit.run_sequential(_entries(), passes=16, should_stop=lambda: True)
    assert len(res.history) == 1                  # stopped after the first check
    assert res is not None


def _two(sample, specs):
    """A recipe with the given (label, position) lines."""
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                  sample=sample, sites=[SiteModel(model="gauss_lor", label=lab, params={
                      "isotropic_chemical_shift_ppm": Param(pos, min=0, max=40),
                      "shift_fwhm_ppm": Param(5.0, min=0.1),
                      "amplitude": Param(80.0, min=0),
                      "gl": Param(1.0, vary=False)}) for lab, pos in specs])


def test_seed_from_pairs_lines_by_label_and_an_empty_set_carries_nothing():
    a = _two("a", [("A", 10.0), ("B", 20.0)])
    b = _two("b", [("B", 22.0), ("A", 12.0), ("C", 30.0)])
    pairs = seqfit.seed_from(a, b, ("isotropic_chemical_shift_ppm",))
    assert pairs == [1, 0]
    assert a.sites[0].params["isotropic_chemical_shift_ppm"].value == 12.0
    assert a.sites[1].params["isotropic_chemical_shift_ppm"].value == 22.0
    # params=() is "nothing": every value stays (None would mean everything)
    seqfit.seed_from(a, _two("c", [("A", 33.0), ("B", 34.0)]), ())
    assert a.sites[0].params["isotropic_chemical_shift_ppm"].value == 12.0
    res = seqfit.run_sequential(_entries()[:2], passes=1, propagate=())
    assert res.propagated == ()


def test_fixed_members_seed_their_neighbours_but_are_never_refitted():
    entries = _entries()
    kept = entries[0][0]
    kept.sites[0].params["isotropic_chemical_shift_ppm"].value = 11.0   # deliberately off
    res = seqfit.run_sequential(entries, passes=2, fixed=(0,))
    assert res.fixed == (0,)
    assert res.recipes[0] is kept
    assert kept.sites[0].params["isotropic_chemical_shift_ppm"].value == 11.0   # untouched
    assert kept.sites[0].params["amplitude"].value == 80.0
    pos = [r.sites[0].params["isotropic_chemical_shift_ppm"].value for r in res.recipes[1:]]
    assert pos[-1] == pytest.approx(18.0, abs=0.6)              # the rest were fitted
    assert all(np.isfinite(h["rmsd"][0]) for h in res.history)  # measured as it stands
    assert "6 spectra, 1 kept" in res.summary
    with pytest.raises(ValueError, match="every spectrum is kept"):
        seqfit.run_sequential(entries, passes=1, fixed=range(6))


def test_smooth_trajectories_pairs_by_label_and_skips_absent_members():
    specs = [[("A", 10.0)], [("A", 14.0), ("C", 30.0)], [("A", 11.0), ("C", 36.0)],
             [("A", 15.0), ("C", 30.0)], [("A", 12.0)]]
    recs = [_two(f"g{k}", s) for k, s in enumerate(specs)]
    seqfit.smooth_trajectories(recs, ("isotropic_chemical_shift_ppm",), list(range(5)), 3)
    a = [r.sites[0].params["isotropic_chemical_shift_ppm"].value for r in recs]
    assert a == pytest.approx([11.3333, 11.6667, 13.3333, 12.6667, 13.0], abs=1e-3)
    c = [r.sites[1].params["isotropic_chemical_shift_ppm"].value for r in recs[1:4]]
    assert c == pytest.approx([32.0, 32.0, 32.0], abs=1e-9)      # smoothed over its 3 members
    assert len(recs[0].sites) == 1 and len(recs[4].sites) == 1   # nothing invented
