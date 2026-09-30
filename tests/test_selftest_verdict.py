"""The self-test must see an aborted or barely-run solver.

Its 'parameters moved by more than 1e-9' rule was satisfied by the amplitude
pre-scale and the first Jacobian probe alone, so the 0.14.x self-test passed
on this machine while the completion-threshold callback was aborting every
Fit-button fit after five evaluations. solver_verdict reads the lmfit result
instead.
"""
from types import SimpleNamespace

import numpy as np

from larmor import selftest


def test_solver_verdict_flags_aborts_and_first_sweep_stops():
    ok = SimpleNamespace(aborted=False, nfev=106, message="`ftol` termination condition is satisfied.")
    assert selftest.solver_verdict(ok) == ""
    aborted = SimpleNamespace(aborted=True, nfev=5,
                              message="Fit aborted by user callback. Could not estimate error-bars.")
    v = selftest.solver_verdict(aborted)
    assert "aborted" in v and "user callback" in v
    early = SimpleNamespace(aborted=False, nfev=5, message="`xtol` termination condition is satisfied.")
    v2 = selftest.solver_verdict(early)
    assert "5 evaluations" in v2 and "Jacobian" in v2
    assert selftest.solver_verdict(None) == "no lmfit result was kept"
    assert selftest.MIN_EVALUATIONS == 8


def test_core_fits_fail_on_an_aborted_solver_even_when_values_moved(monkeypatch):
    """A fit whose amplitude pre-scale moved the values but whose solver was
    aborted at once -- exactly the 0.14.x signature -- must FAIL the self-test."""
    from larmor import fit as fitmod

    real_fit = fitmod.fit

    def aborted_fit(rec, x, y, **kw):
        res = real_fit(rec, x, y, **kw)
        lm = SimpleNamespace(aborted=True, nfev=5,
                             message="Fit aborted by user callback. Could not estimate error-bars.")
        res.lmfit_result = lm
        return res

    monkeypatch.setattr(fitmod, "fit", aborted_fit)
    lines = []

    class Say:
        def __call__(self, *parts):
            lines.append(" ".join(str(p) for p in parts))

        def fail(self, stage, exc=None):
            lines.append("FAIL " + stage)

    assert selftest.core_fits(Say()) is False
    assert any(l.startswith("FAIL core gauss_lor: the solver was aborted") for l in lines), lines


def test_core_fits_pass_with_the_real_solver():
    lines = []

    class Say:
        def __call__(self, *parts):
            lines.append(" ".join(str(p) for p in parts))

        def fail(self, stage, exc=None):
            lines.append("FAIL " + stage)

    assert selftest.core_fits(Say()) is True
    assert not any(l.startswith("FAIL") for l in lines), lines
    assert all("solver: ran" in l for l in lines if l.startswith("core ")), lines
    assert np.isfinite(1.0)
