"""Identifiability: flag parameter pairs the data cannot separate (|r| >= 0.95)."""
import numpy as np

from larmor import identifiability as idn


class _FakeLM:
    def __init__(self, names, corr):
        self.var_names = names
        d = np.ones(len(names))
        self.covar = np.asarray(corr, float) * np.outer(d, d)


def test_no_pairs_when_uncorrelated():
    lm = _FakeLM(["a", "b"], [[1.0, 0.1], [0.1, 1.0]])
    assert idn.unidentifiable_pairs(lm) == []


def test_flags_degenerate_pair():
    lm = _FakeLM(["pos", "width", "amp"],
                 [[1.0, 0.98, 0.2], [0.98, 1.0, 0.1], [0.2, 0.1, 1.0]])
    uni = idn.unidentifiable_pairs(lm)
    assert len(uni) == 1
    a, b, r = uni[0]
    assert {a, b} == {"pos", "width"} and r > 0.95


def test_threshold_is_configurable():
    lm = _FakeLM(["a", "b"], [[1.0, 0.9], [0.9, 1.0]])
    assert idn.unidentifiable_pairs(lm) == []                 # 0.9 < 0.95
    assert idn.unidentifiable_pairs(lm, thresh=0.85)          # now flagged


def test_missing_covariance_returns_empty():
    class NoCov:
        var_names = ["a"]
        covar = None
    assert idn.unidentifiable_pairs(NoCov()) == []
    assert idn.corr_matrix(NoCov())[1] is None
