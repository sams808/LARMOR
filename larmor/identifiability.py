"""Parameter identifiability from a fit's covariance.

Two parameters with |correlation| ≈ 1 trade off perfectly — the data cannot
separate them, so their individual error bars are meaningless even if the fit
"converged". This module turns the lmfit covariance into a correlation matrix and
lists the **unidentifiable pairs** (|r| ≥ 0.95 by default) so the human knows
which numbers not to over-interpret. Qt-free and testable.
"""
from __future__ import annotations

import numpy as np

IDENTIFIABILITY_THRESHOLD = 0.95


def corr_matrix(lmfit_result):
    """(var_names, correlation matrix) from an lmfit result, or (names, None)."""
    names = list(getattr(lmfit_result, "var_names", []) or [])
    cov = getattr(lmfit_result, "covar", None)
    if not names or cov is None:
        return names, None
    cov = np.asarray(cov, float)
    if cov.ndim != 2 or cov.shape[0] != len(names):
        return names, None
    d = np.sqrt(np.clip(np.diag(cov), 1e-300, None))
    return names, cov / np.outer(d, d)


def unidentifiable_pairs(lmfit_result,
                         thresh: float = IDENTIFIABILITY_THRESHOLD) -> list[tuple]:
    """List of ``(name_i, name_j, r)`` with |r| ≥ ``thresh``, strongest first."""
    names, corr = corr_matrix(lmfit_result)
    if corr is None:
        return []
    out = []
    n = len(names)
    for i in range(n):
        for j in range(i + 1, n):
            r = float(corr[i, j])
            if np.isfinite(r) and abs(r) >= thresh:
                out.append((names[i], names[j], r))
    out.sort(key=lambda t: -abs(t[2]))
    return out
