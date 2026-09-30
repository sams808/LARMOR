"""The Fit button honours the completion threshold through the solver, never
through the progress callback.

0.14.0-0.14.4 shipped a callback rule -- 'stop once the residual stdev changes
by less than the threshold between two calls' -- but lmfit calls iter_cb on
every residual EVALUATION, and most evaluations are finite-difference
Jacobian probes that change the residual by a relative 1e-8 (analytic
models) to 1e-4 (kernel models). At the default 0.1 % every Fit-button fit
aborted after ~5 evaluations with the positions untouched ('visibly not
fitting'), while Auto fit -- no threshold, no callback -- converged.
Measured on the shipped 27Al example: nfev 5 / RMSD 0.082 at 0.1 %,
nfev 273 / RMSD 0.0536 at 0.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _model_and_data():
    from larmor.engine import make_context, simulate_site
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20, 120, 400)
    truth = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="gauss_lor", label="p", params={
            "isotropic_chemical_shift_ppm": Param(60, min=0, max=120),
            "shift_fwhm_ppm": Param(8, min=1, max=40),
            "gl": Param(0.5, min=0, max=1, vary=False),
            "amplitude": Param(1e6, min=0)})])
    ctx = make_context(truth, exp_ppm=x)
    y = np.sum([simulate_site(s, ctx) for s in truth.sites], axis=0)
    start = Recipe.from_dict(truth.to_dict())
    # a user's start: 2 ppm off, a quarter too wide, 30 % too small (the
    # abort left exactly these numbers on screen, the amplitude aside)
    start.sites[0].params["isotropic_chemical_shift_ppm"].value = 62.0
    start.sites[0].params["shift_fwhm_ppm"].value = 10.0
    start.sites[0].params["amplitude"].value = 7e5
    return start, x, y


def test_ftol_mapping_is_the_only_threshold_mechanism():
    from larmor.fit import _tol_kws, ftol_from_pct

    assert ftol_from_pct(0.1) == pytest.approx(0.002)
    assert ftol_from_pct(0) is None and ftol_from_pct(None) is None
    # ftol only: an xtol matched to it stopped the solver as soon as the
    # amplitude (|x| ~ 1e6) settled, 1 ppm short of the true position
    assert _tol_kws(0.1) == {"ftol": pytest.approx(0.002)}
    assert "xtol" not in _tol_kws(0.1)
    assert _tol_kws(0) == {}


@pytest.mark.parametrize("threshold", [0.1, 0.0])
def test_fit_worker_moves_the_line_at_the_default_threshold(qapp, monkeypatch,
                                                            threshold):
    from larmor.desktop import workers
    from larmor.desktop.app import FitWorker

    monkeypatch.setattr(workers, "_fit_tol", lambda: threshold)
    start, x, y = _model_and_data()
    fw = FitWorker(start.to_dict(), x, y, (120, -20), animate=False)
    done = {}
    fw.done.connect(lambda res, mode: done.update(res=res, mode=mode))
    fw.failed.connect(lambda msg: done.update(failed=msg))
    fw.run()
    assert "failed" not in done, done.get("failed")
    res = done["res"]
    lm = res.lmfit_result
    assert not getattr(lm, "aborted", False), lm.message
    assert "callback" not in str(lm.message).lower()
    assert lm.nfev > 10                       # more than one Jacobian sweep
    site = res.recipe.sites[0].params
    assert site["isotropic_chemical_shift_ppm"].value == pytest.approx(60.0, abs=0.2)
    assert site["shift_fwhm_ppm"].value == pytest.approx(8.0, abs=0.3)
    assert site["amplitude"].value == pytest.approx(1e6, rel=0.05)
    assert res.rmsd < 0.01
