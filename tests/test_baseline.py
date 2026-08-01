"""Regressions for the iterative baseline corrector (Yon et al. 2020, SSNMR 110).

Guards the two properties that matter: (1) a rolling baseline under peaks is
removed and the peak AREAS survive; (2) a spectrum that already has a good
baseline is left essentially untouched (the corrector must never eat signal).
"""
import numpy as np
import pytest


def _spectrum(n=4096, seed=0):
    rng = np.random.default_rng(seed)
    x = np.arange(n)
    g = lambda c, w, a: a * np.exp(-0.5 * ((x - c) / w) ** 2)
    peaks = g(1500, 40, 100) + g(2200, 25, 60) + g(2600, 60, 45)
    return x, peaks, rng.normal(0, 1.5, n)


def test_removes_rolling_baseline_and_preserves_peak_area():
    from larmor.baseline import iterative_baseline
    x, peaks, noise = _spectrum()
    # a slow rolling baseline (the dead-time class this method targets)
    roll = 30 * np.sin(2 * np.pi * x / 6000) + 25 * np.exp(-((x - 3200) / 1500) ** 2) + 12
    r = iterative_baseline(peaks + roll + noise)
    assert r.converged
    pk = peaks > 5
    area_off = abs(r.corrected[pk].sum() - peaks[pk].sum()) / peaks[pk].sum()
    assert area_off < 0.03                       # peak area preserved to < 3 %
    base_off = np.sqrt(np.mean((r.baseline - roll) ** 2)) / np.ptp(roll)
    assert base_off < 0.05                       # baseline recovered to < 5 %


def test_leaves_a_good_spectrum_alone():
    """On a spectrum whose baseline is already flat, the estimated baseline must
    be tiny relative to the signal -- the corrector must not invent one."""
    from larmor.baseline import iterative_baseline
    x, peaks, noise = _spectrum()
    r = iterative_baseline(peaks + noise)
    assert np.abs(r.baseline).max() < 0.03 * peaks.max()
    pk = peaks > 5
    assert r.corrected[pk].max() == pytest.approx(peaks[pk].max(), rel=0.05)


def test_dead_time_restriction_runs_and_is_broad():
    """The dead-time (FID-truncation) restriction keeps only broad baseline
    structure: the returned baseline must be smooth (small point-to-point step)."""
    from larmor.baseline import iterative_baseline
    x, peaks, noise = _spectrum()
    roll = 40 * np.exp(-((x - 2048) / 1800) ** 2) + 10
    r = iterative_baseline(peaks + roll + noise, dead_time_pts=64)
    # a broad baseline changes slowly point-to-point
    step = np.max(np.abs(np.diff(r.baseline)))
    assert step < 0.1 * np.ptp(r.baseline) + 1.0
    assert r.baseline.shape == peaks.shape
