"""Multi-experiment correlation engine (architecture; not wired into the app)."""
import numpy as np
import pytest

from larmor import correlate as X


def _obs(label, ppm, amp, kind="1d"):
    return X.Observable(label, ppm, amp, kind)


def test_difference_reproduces_hmqc_case():
    """The two-dataset difference (a 1D minus an HMQC projection) isolates the
    un-correlated feature — the existing HMQC behaviour as a special case."""
    x = np.linspace(-50, 50, 400)
    oned = np.exp(-((x - 10) / 4) ** 2) + 0.7 * np.exp(-((x + 20) / 4) ** 2)
    proj = np.exp(-((x - 10) / 4) ** 2)                 # correlated peak only
    grid, diff = X.difference(_obs("1D", x, oned), [_obs("HMQC", x, proj, "projection")])
    at10 = diff[np.argmin(np.abs(grid - 10))]
    atm20 = diff[np.argmin(np.abs(grid + 20))]
    assert abs(at10) < 0.1          # correlated peak removed
    assert atm20 > 0.5              # un-correlated peak kept


def test_intersection_is_the_common_feature():
    """Three experiments; only one peak appears in all → intersection keeps it."""
    x = np.linspace(-50, 50, 400)
    common = np.exp(-((x - 5) / 4) ** 2)
    a = common + np.exp(-((x + 30) / 4) ** 2)
    b = common + np.exp(-((x + 15) / 4) ** 2)
    c = common + np.exp(-((x - 25) / 4) ** 2)
    grid, inter = X.intersection([_obs("A", x, a), _obs("B", x, b), _obs("C", x, c)])
    assert inter[np.argmin(np.abs(grid - 5))] > 0.8     # common peak survives
    for pos in (-30, -15, 25):
        assert inter[np.argmin(np.abs(grid - pos))] < 0.2   # specific peaks gone


def test_scale_to_is_least_squares():
    x = np.linspace(-10, 10, 200)
    src = np.exp(-(x / 2) ** 2)
    assert X.scale_to(2.5 * src, src) == pytest.approx(2.5, rel=1e-6)


def test_align_on_common_grid():
    a = _obs("a", np.linspace(-10, 10, 50), np.ones(50))
    b = _obs("b", np.linspace(0, 30, 80), np.ones(80))
    grid, amps = X.align([a, b])
    assert grid[0] == pytest.approx(-10.0) and grid[-1] == pytest.approx(30.0)
    assert len(amps) == 2 and amps[0].shape == grid.shape
