"""Residual diagnostics: runs test + lag-1 autocorrelation catch a structured
(systematically wrong) fit even when the RMSD is small."""
import numpy as np

from larmor import diagnostics


def test_white_noise_is_not_flagged():
    rng = np.random.default_rng(0)
    r = rng.normal(0, 1, 500)
    assert not diagnostics.residual_structure(r)["structured"]
    assert abs(diagnostics.runs_test(r)["z"]) < 3.0
    assert abs(diagnostics.lag1_autocorr(r)) < 0.2


def test_structured_residual_is_flagged():
    # a smooth sine (long same-sign runs, high autocorrelation) = clear structure
    x = np.linspace(0, 6 * np.pi, 500)
    r = np.sin(x)
    s = diagnostics.residual_structure(r)
    assert s["structured"]
    assert diagnostics.runs_test(r)["runs"] < 20        # very few runs
    assert diagnostics.lag1_autocorr(r) > 0.9
    assert "structured" in s["message"]


def test_runs_test_handles_short_and_degenerate_input():
    assert diagnostics.runs_test([1, -1, 1])["structured"] is False   # too short
    # degenerate inputs must not crash and return a valid verdict dict
    for bad in (np.ones(50), np.zeros(50), np.array([]), np.array([1.0])):
        out = diagnostics.runs_test(bad)
        assert set(out) >= {"z", "p", "runs", "structured"}


def test_small_rmsd_but_structured():
    # a tiny-amplitude but coherent residual — RMSD looks fine, structure caught
    x = np.linspace(0, 8 * np.pi, 400)
    r = 0.01 * np.sin(x)
    assert diagnostics.residual_structure(r)["structured"]
