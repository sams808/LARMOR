import os

import numpy as np

# engine tests build (tiny) kernels; never let them write the real
# on-disk kernel cache under %LOCALAPPDATA% (B5) -- the two disk-cache
# tests below opt back in with their own tmp LOCALAPPDATA
os.environ.setdefault("LARMOR_NO_SESSION", "1")
import pytest

from larmor import engine
from larmor.io import fxmla
from larmor import fit as fitmod
from larmor.recipe import Recipe, SiteModel, Param

from conftest import CAALGLASS, require


def _degenerate_recipe():
    """Two IDENTICAL-shape overlapping sites: their amplitudes are perfectly
    correlated (only the SUM is determined by the data), so the covariance
    matrix is singular and the first least_squares pass reliably fails to get
    error bars -- a controlled trigger for the errorbar-rescue retry."""
    x = np.linspace(-30, 30, 500)
    truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(0.0),
            "shift_fwhm_ppm": Param(10.0), "amplitude": Param(80.0),
            "gl": Param(1.0, vary=False)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    data = y + np.random.default_rng(0).normal(0, 0.3, x.size)
    r = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(0.0, vary=False),
            "shift_fwhm_ppm": Param(10.0, vary=False),
            "amplitude": Param(40.0, min=0), "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="B", params={
            "isotropic_chemical_shift_ppm": Param(0.0, vary=False),
            "shift_fwhm_ppm": Param(10.0, vary=False),
            "amplitude": Param(40.0, min=0), "gl": Param(1.0, vary=False)}),
    ])
    return r, x, data


def test_amplitude_uvars_matches_param_stderr_and_goes_stale_on_edit():
    """fit.amplitude_uvars (N3): the ufloat per site carries the covariance
    stderr lmfit wrote on the Param (and its correlations), returns None once
    an amplitude was edited after the fit, and None without a covariance."""
    from types import SimpleNamespace

    x = np.linspace(-20, 40, 600)
    truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(12.0), "shift_fwhm_ppm": Param(8.0),
            "amplitude": Param(100.0), "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="B", params={
            "isotropic_chemical_shift_ppm": Param(6.0), "shift_fwhm_ppm": Param(6.0),
            "amplitude": Param(60.0), "gl": Param(1.0, vary=False)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    data = y + np.random.default_rng(3).normal(0, 0.5, x.size)
    rec = Recipe.from_dict(truth.to_dict())
    for s in rec.sites:
        s.params["amplitude"].value *= 0.8
    res = fitmod.fit(rec, x, data, window_ppm=(40, -20))
    uv = fitmod.amplitude_uvars(res.lmfit_result, rec)
    assert uv is not None and set(uv) == {0, 1}
    for i in range(2):
        p = rec.sites[i].params["amplitude"]
        assert uv[i].nominal_value == pytest.approx(p.value, rel=1e-12)
        assert uv[i].std_dev == pytest.approx(p.stderr, rel=1e-9)
    # overlapping lines: the two amplitudes are correlated, not independent
    from uncertainties import covariance_matrix
    cm = covariance_matrix([uv[0], uv[1]])
    assert cm[0][1] != 0.0
    # a nudge after the fit makes the covariance stale -> None
    rec.sites[0].params["amplitude"].value *= 1.001
    assert fitmod.amplitude_uvars(res.lmfit_result, rec) is None
    rec.sites[0].params["amplitude"].value /= 1.001
    assert fitmod.amplitude_uvars(res.lmfit_result, rec) is not None
    # no covariance at all (lmfit sets uvars None), or no result -> None
    assert fitmod.amplitude_uvars(SimpleNamespace(uvars=None), rec) is None
    assert fitmod.amplitude_uvars(None, rec) is None


def test_compute_errorbars_false_skips_the_rescue_retry():
    """A degenerate covariance normally triggers a second full optimization
    (leastsq from the converged point) to try to recover error bars.
    compute_errorbars=False must skip that retry entirely -- verified by
    counting actual evaluations, not just checking the final stderr."""
    r1, x, data = _degenerate_recipe()
    calls1 = []
    fitmod.fit(r1, x, data, compute_errorbars=True,
              iter_cb=lambda *a, **k: calls1.append(1))

    r2, x, data = _degenerate_recipe()
    calls2 = []
    fitmod.fit(r2, x, data, compute_errorbars=False,
              iter_cb=lambda *a, **k: calls2.append(1))

    # with the retry skipped, meaningfully fewer evaluations happened
    assert len(calls2) < len(calls1)
    # the fitted (identifiable) TOTAL amplitude is unaffected either way
    total1 = sum(s.params["amplitude"].value for s in r1.sites)
    total2 = sum(s.params["amplitude"].value for s in r2.sites)
    assert total1 == pytest.approx(total2, rel=0.05)


def test_gauss_lor_shapes():
    x = np.linspace(-50, 50, 2001)
    g = engine.gauss_lor(x, 0.0, 10.0, 2.0, gl=1.0)
    l = engine.gauss_lor(x, 0.0, 10.0, 2.0, gl=0.0)
    assert g.max() == pytest.approx(2.0, rel=1e-6)
    assert l.max() == pytest.approx(2.0, rel=1e-6)
    # both peak-normalized with the same FWHM...
    half = np.where(g >= 1.0)[0]
    assert x[half[-1]] - x[half[0]] == pytest.approx(10.0, abs=0.2)
    # ...but the Lorentzian has heavier tails
    tail = np.abs(x) > 25
    assert l[tail].sum() > 5 * g[tail].sum()


@pytest.mark.slow
def test_kernel_and_fit_caalglass():
    """End-to-end: import the dmfit fit, refine with lmfit, beat Phase 0's RMSD."""
    dm = fxmla.read(require(CAALGLASS))
    recipe, _ = fxmla.to_recipe(dm)

    kernel = engine.build_kernel(
        recipe.nucleus, recipe.larmor_frequency_MHz, recipe.spin_rate_Hz)
    assert kernel.K.shape == (80 * 11, 2048)
    assert np.all(np.diff(kernel.x_ppm) > 0)

    # sigma reweighting sanity: larger sigma -> broader lineshape
    y_narrow = kernel.weights(1.0) @ kernel.K
    y_broad = kernel.weights(4.0) @ kernel.K
    width = lambda y: np.sum(y > y.max() / 2)
    assert width(y_broad) > width(y_narrow)

    exp_ppm, exp_amp = dm.spectrum.ppm, dm.spectrum.amplitude
    result = fitmod.fit(recipe, exp_ppm, exp_amp, window_ppm=(150.0, -80.0))

    # refined fit should beat the fixed-parameter Phase 0 replay (RMSD 0.027)
    assert result.rmsd < 0.01

    # dmfit's ad-hoc Gauss/Lor sideband lines (236/208 ppm) sit outside the
    # window and must be frozen, not left to wander
    assert len(result.frozen_sites) == 2

    site1 = recipe.sites[0]
    assert site1.params["isotropic_chemical_shift_ppm"].value == pytest.approx(
        66.2, abs=5.0)
    assert site1.params["sigma_Cq_MHz"].value == pytest.approx(2.0, abs=0.7)
    # the whole point: uncertainties exist now
    assert result.lmfit_result.errorbars
    assert site1.params["sigma_Cq_MHz"].stderr is not None
    assert 0 < site1.params["sigma_Cq_MHz"].stderr < 0.5


def test_grid_restrictable_allowlist():
    from larmor import engine

    def rec(*models):
        return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
            SiteModel(model=m, label=str(i), params={"amplitude": Param(1.0)})
            for i, m in enumerate(models)])

    assert engine.grid_restrictable(rec("gauss_lor"))
    assert engine.grid_restrictable(rec("czjzek"))
    assert engine.grid_restrictable(rec("gauss_lor", "voigt", "czjzek"))
    # the LRU single-site models build their own cached simulation directly
    # from ctx.x_ppm's span and convolution-broaden it -- restricting the grid
    # there would truncate real signal (satellite/sideband manifolds), not
    # just its cost, so they must NOT be treated as restrictable
    for unsafe in ("quad_ct", "quad_first", "quad_csa", "csa_mas", "csa_czjzek"):
        assert not engine.grid_restrictable(rec(unsafe)), unsafe
    # ANY unsafe site anywhere in the recipe disqualifies the whole fit (ctx is
    # shared across all sites -- there's no such thing as a partial restriction)
    assert not engine.grid_restrictable(rec("gauss_lor", "quad_ct"))
    # an unregistered/unaudited future model defaults to NOT restrictable
    assert not engine.grid_restrictable(rec("some_future_model"))


def test_fit_windowed_grid_recovers_truth_and_returns_full_span():
    """A narrow fit window against a MUCH wider experimental sweep must still
    recover the true parameters (the optimisation-time grid restriction can't
    change the answer), and the returned model curve must still span the FULL
    experiment (the restriction is an internal speed optimisation, invisible
    to callers/the UI -- not a change in what gets displayed)."""
    x = np.linspace(-150, 150, 3000)
    truth = Recipe(nucleus="31P", larmor_frequency_MHz=162.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(10.0),
            "shift_fwhm_ppm": Param(5.0), "amplitude": Param(100.0),
            "gl": Param(1.0, vary=False)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    data = y + np.random.default_rng(1).normal(0, 0.5, x.size)

    r = Recipe(nucleus="31P", larmor_frequency_MHz=162.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(9.0, min=-20, max=40),
            "shift_fwhm_ppm": Param(6.0, min=0.1),
            "amplitude": Param(80.0, min=0), "gl": Param(1.0, vary=False)})])
    result = fitmod.fit(r, x, data, window_ppm=(40.0, -20.0))

    p = r.sites[0].params
    assert p["isotropic_chemical_shift_ppm"].value == pytest.approx(10.0, abs=0.1)
    assert p["shift_fwhm_ppm"].value == pytest.approx(5.0, abs=0.2)
    assert p["amplitude"].value == pytest.approx(100.0, rel=0.05)
    # the fit window (±60 ppm) is far narrower than the full ±150 ppm sweep --
    # the returned curve must still cover the whole thing
    assert result.x_ppm.min() <= -149.0 and result.x_ppm.max() >= 149.0
    assert len(result.x_ppm) == len(x)
    assert len(result.y_fit) == len(result.x_ppm)


def test_fit_windowed_grid_matches_unrestricted_for_excluded_model():
    """A quad_ct fit (excluded from the grid-restriction allowlist) must give
    the identical fitted values whether or not grid_restrictable would apply
    to other models -- verifies the exclusion actually prevents any change,
    not just that the recipe was skipped by convention."""
    x = np.linspace(-40, 140, 1200)
    truth = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="quad_ct", label="A", params={
            "isotropic_chemical_shift_ppm": Param(60.0),
            "Cq_MHz": Param(3.0, min=0.01, max=20),
            "eta": Param(0.2, vary=False),
            "shift_fwhm_ppm": Param(2.0, min=0.05),
            "amplitude": Param(100.0, min=0)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    data = y + np.random.default_rng(2).normal(0, y.max() * 0.005, x.size)

    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="quad_ct", label="A", params={
            "isotropic_chemical_shift_ppm": Param(58.0, min=40, max=80),
            "Cq_MHz": Param(3.2, min=0.01, max=20),
            "eta": Param(0.2, vary=False),
            "shift_fwhm_ppm": Param(2.2, min=0.05),
            "amplitude": Param(90.0, min=0)})])
    fitmod.fit(r, x, data, window_ppm=(140.0, -40.0))
    p = r.sites[0].params
    assert p["isotropic_chemical_shift_ppm"].value == pytest.approx(60.0, abs=0.5)
    assert p["Cq_MHz"].value == pytest.approx(3.0, abs=0.3)


def test_kernel_axis_ppm_is_bin_exact_against_a_real_build():
    """make_context derives the kernel axis arithmetically instead of paying
    a full kernel build for it. The formula must be BIN-EXACT (the ppm
    conversion runs against the isotope's reference frequency, not gamma*B0
    -- 0.08 % apart for 27Al, which would be nearly half a ppm at the axis
    edge)."""
    from larmor import engine

    kernel = engine.build_kernel("27Al", 130.32, 12500.0, sw_Hz=150000.0,
                                 npts=256, ref_offset_ppm=30.0,
                                 cq_max_MHz=6.0, n_cq=4, n_eta=3)
    axis = engine.kernel_axis_ppm("27Al", 130.32, sw_Hz=150000.0, npts=256,
                                  ref_offset_ppm=30.0)
    assert axis.shape == kernel.x_ppm.shape
    assert np.allclose(axis, kernel.x_ppm, atol=1e-9, rtol=0)


def test_make_context_does_not_build_a_kernel(monkeypatch):
    """The documented 'two kernel builds' defect: make_context built a full
    kernel and kept only its axis. It must now never call build_kernel."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    def boom(*a, **k):
        raise AssertionError("make_context built a kernel for its axis")

    monkeypatch.setattr(engine, "build_kernel", boom)
    site = SiteModel(model="czjzek", label="c", params={
        "isotropic_chemical_shift_ppm": Param(60.0),
        "sigma_Cq_MHz": Param(2.0), "shift_fwhm_ppm": Param(8.0),
        "amplitude": Param(1.0)})
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.32,
               spin_rate_Hz=12500.0, sites=[site])
    ctx = engine.make_context(r, exp_ppm=np.linspace(-300, 300, 2048))
    assert ctx.x_ppm.size >= 2048


def test_kernel_cache_is_bounded_lru(monkeypatch):
    """B1: the kernel cache evicts least-recently-used entries past its byte
    budget, but never below the keep-min (a fit straddling a Cq-ladder step
    alternates two kernels; evicting one would rebuild every iteration)."""
    from collections import OrderedDict

    from larmor import engine

    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    monkeypatch.setattr(engine, "KERNEL_CACHE_BUDGET_MB", 1.0)   # 1 MB budget

    class FakeKernel:
        def __init__(self, mb):
            self.K = np.zeros(int(mb * 1e6 // 4), dtype=np.float32)
            self.x_ppm = np.zeros(2)

    for i in range(8):                       # 8 x 0.5 MB = 4x the budget
        engine._cache_put(("k", i), FakeKernel(0.5))
    assert len(engine._KERNEL_CACHE) == engine._CACHE_KEEP_MIN
    assert list(engine._KERNEL_CACHE) == [("k", i) for i in (4, 5, 6, 7)]

    # a get refreshes recency: touch the oldest survivor, add one more
    engine.build_kernel  # (real one not needed: exercise the LRU directly)
    engine._KERNEL_CACHE.move_to_end(("k", 4))
    engine._cache_put(("k", 8), FakeKernel(0.5))
    assert ("k", 4) in engine._KERNEL_CACHE      # refreshed -> survived
    assert ("k", 5) not in engine._KERNEL_CACHE  # LRU -> evicted

    info = engine.kernel_cache_info()
    assert info["entries"] == 4 and info["mb"] == pytest.approx(2.0, rel=0.01)


def test_kernel_basis_is_float32():
    """B3: the (Cq, eta) basis ships as float32 -- the grid's own quantisation
    dwarfs 1e-7 relative precision, and the example-fit RMSDs are unchanged
    to <0.1 %. Halves both memory and the per-iteration matmul."""
    from larmor import engine

    k = engine.build_kernel("27Al", 130.32, 12500.0, sw_Hz=150000.0,
                            npts=128, ref_offset_ppm=30.0,
                            cq_max_MHz=6.0, n_cq=3, n_eta=3)
    assert k.K.dtype == np.float32
    y = k.weights(2.0) @ k.K                 # the fit-path multiply
    assert y.dtype == np.float64 and np.all(np.isfinite(y))


def test_kernel_disk_cache_round_trip(monkeypatch, tmp_path):
    """B5: a built kernel persists to disk and a cold in-memory cache loads
    it back bin-identical -- this is what spares a fresh worker process (or
    tomorrow's session) the ~23 s wideline rebuild. Gated off under
    LARMOR_NO_SESSION so the suite itself exercises it only via this
    explicit opt-in."""
    from collections import OrderedDict

    from larmor import engine

    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())

    k1 = engine.build_kernel("27Al", 130.32, 12500.0, sw_Hz=150000.0,
                             npts=128, ref_offset_ppm=30.0,
                             cq_max_MHz=6.0, n_cq=3, n_eta=3)
    files = list((tmp_path / "LARMOR" / "kernels").glob("*.npz"))
    assert len(files) == 1

    # cold memory, warm disk: must come back identical without simulating
    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())

    def boom(*a, **k):
        raise AssertionError("simulated despite a disk hit")

    import mrsimulator
    monkeypatch.setattr(mrsimulator, "Simulator", boom)
    k2 = engine.build_kernel("27Al", 130.32, 12500.0, sw_Hz=150000.0,
                             npts=128, ref_offset_ppm=30.0,
                             cq_max_MHz=6.0, n_cq=3, n_eta=3)
    assert np.array_equal(k1.K, k2.K)
    assert np.array_equal(k1.x_ppm, k2.x_ppm)

    # a corrupt file is discarded and rebuilt, not fatal
    files[0].write_bytes(b"garbage")
    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    monkeypatch.undo()  # restore Simulator for the rebuild
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    k3 = engine.build_kernel("27Al", 130.32, 12500.0, sw_Hz=150000.0,
                             npts=128, ref_offset_ppm=30.0,
                             cq_max_MHz=6.0, n_cq=3, n_eta=3)
    assert np.array_equal(k1.K, k3.K)


def test_kernel_disk_cache_off_under_no_session(monkeypatch, tmp_path):
    from collections import OrderedDict

    from larmor import engine

    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    engine.build_kernel("27Al", 130.32, 12500.0, sw_Hz=150000.0,
                        npts=128, ref_offset_ppm=30.0,
                        cq_max_MHz=6.0, n_cq=3, n_eta=3)
    assert not (tmp_path / "LARMOR" / "kernels").exists()


def test_kernel_build_chunked_is_identical_and_cancellable(monkeypatch):
    """B4: with feedback hooks registered the (Cq, eta) systems simulate in
    chunks -- output must be BIT-identical to the monolithic run, progress
    must count up, and a cancel between chunks must abort cleanly."""
    from collections import OrderedDict

    from larmor import engine

    args = dict(sw_Hz=150000.0, npts=128, ref_offset_ppm=30.0,
                cq_max_MHz=6.0, n_cq=5, n_eta=3)
    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    mono = engine.build_kernel("27Al", 130.32, 12500.0, **args)

    seen = []
    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    with engine.kernel_build_feedback(cb=lambda i, n: seen.append((i, n))):
        chunked = engine.build_kernel("27Al", 130.32, 12500.0, **args)
    assert seen and seen[-1][0] == seen[-1][1] > 1
    assert np.array_equal(mono.K, chunked.K)
    assert np.array_equal(mono.x_ppm, chunked.x_ppm)

    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    calls = {"n": 0}

    def cancel_after_two():
        calls["n"] += 1
        return calls["n"] > 2

    with pytest.raises(engine.KernelBuildCancelled):
        with engine.kernel_build_feedback(cancel=cancel_after_two):
            engine.build_kernel("27Al", 130.32, 12500.0, **args)
    # nothing half-built may enter the caches
    assert len(engine._KERNEL_CACHE) == 0


# ---------------------------------------------------------------------------
# WA (2026-09): the czjzek_d / czjzek_corr fits that took tens of minutes and
# could not be stopped. Root causes, each pinned below:
#   * the optimiser's first trust-region step (scipy trf: radius = the largest
#     parameter, the AMPLITUDE ~3e6) asked _broaden_shift for a 3.15e6 ppm
#     Gaussian -- a 2.4-million-point kernel on a 2560-point axis, minutes per
#     site render;
#   * the error-bar rescue (leastsq) probed its Jacobian at a 1e-5 relative
#     step on a kernel-quantised model and dithered to maxfev (~22 000 evals);
#   * Stop was checked only after a whole residual evaluation / every 8th of a
#     kernel build.
# ---------------------------------------------------------------------------

def _small_kernel_settings(monkeypatch):
    """A tiny kernel so the Czjzek renders below build in well under a second
    and never touch the process-wide cache."""
    from collections import OrderedDict

    monkeypatch.setattr(engine, "_KERNEL_CACHE", OrderedDict())
    monkeypatch.setattr(engine, "KERNEL_SETTINGS",
                        {"npts": 256, "n_cq": 6, "n_eta": 3})


def test_broadening_wider_than_the_axis_is_clamped_and_physical_widths_untouched():
    import time

    from scipy.ndimage import gaussian_filter1d

    from larmor.models.analytic import FWHM_TO_SIGMA
    from larmor.models.quadrupolar import _broaden_shift, _clamp_fwhm

    x = np.linspace(-500.0, 500.0, 2560)
    y = np.exp(-(x / 20.0) ** 2)
    t0 = time.perf_counter()
    runaway = _broaden_shift(x, y, 0.0, 3.15e6)       # the measured request
    assert time.perf_counter() - t0 < 2.0             # was > 40 min
    assert np.array_equal(runaway, _broaden_shift(x, y, 0.0, x[-1] - x[0]))
    assert np.isfinite(runaway).all()
    # a physical width is exactly what it always was
    assert _clamp_fwhm(x, 12.0) == 12.0
    ref = gaussian_filter1d(np.interp(x - 30.0, x, y, left=0.0, right=0.0),
                            12.0 * FWHM_TO_SIGMA / (x[1] - x[0]), mode="constant")
    assert np.array_equal(_broaden_shift(x, y, 30.0, 12.0), ref)


def test_reweight_memo_shared_by_linked_copies_is_bit_identical(monkeypatch):
    from larmor.models import quadrupolar as Q
    from larmor.models.base import SimContext

    _small_kernel_settings(monkeypatch)
    ctx = SimContext("27Al", 130.32, 12500.0, np.linspace(-150.0, 200.0, 256))
    calls = {"n": 0}
    orig = engine.CzjzekKernel.weights

    def counted(self, sigma):
        calls["n"] += 1
        return orig(self, sigma)

    monkeypatch.setattr(engine.CzjzekKernel, "weights", counted)
    v = {"isotropic_chemical_shift_ppm": 60.0, "sigma_Cq_MHz": 1.5,
         "shift_fwhm_ppm": 8.0, "line_fwhm_ppm": 0.0, "amplitude": 2.0,
         "shift_slope_ppm_per_MHz": 1.5, "czjzek_d": 4.0}
    kernel = Q._kernel_for(ctx, 15.0)
    for render, key in ((Q._render_czjzek, "czjzek"),
                        (Q._render_czjzek_corr, "czjzek_corr")):
        kernel.__dict__.pop("_reweight_memo", None)
        calls["n"] = 0
        parent = render(dict(v), ctx)
        copies = [render({**v, "isotropic_chemical_shift_ppm": 60.0 + k * 96.0}, ctx)
                  for k in (-2, -1, 1, 2, 3)]
        assert calls["n"] == 1                 # one reweighting for six renders
        assert list(kernel._reweight_memo)[0][0] == key
        # every copy equals its own uncached render, bit for bit
        for k, y in zip((-2, -1, 1, 2, 3), copies):
            kernel.__dict__.pop("_reweight_memo", None)
            fresh = render({**v, "isotropic_chemical_shift_ppm": 60.0 + k * 96.0}, ctx)
            assert np.array_equal(y, fresh)
        assert parent.max() == pytest.approx(2.0, rel=0.05)   # tiny kernel, interpolated
    # czjzek_d shares the same memo, keyed on (sigma, d)
    kernel.__dict__.pop("_reweight_memo", None)
    a = Q._render_czjzek_d(dict(v), ctx)
    b = Q._render_czjzek_d({**v, "isotropic_chemical_shift_ppm": 10.0}, ctx)
    assert list(kernel._reweight_memo) == [("czjzek_d", 1.5, 4.0)]
    kernel.__dict__.pop("_reweight_memo", None)
    assert np.array_equal(b, Q._render_czjzek_d({**v, "isotropic_chemical_shift_ppm": 10.0}, ctx))
    assert a.max() == pytest.approx(2.0, rel=0.05)


def test_kernel_chunks_are_bounded_in_size():
    assert engine._n_chunks(880) == 10            # 27Al default grid: ~0.7 s each
    assert engine._n_chunks(15) == 8              # small builds keep 8 ticks
    assert engine._n_chunks(14080) == 147         # the 400 MHz ladder step
    assert engine._n_chunks(1) == 1 and engine._n_chunks(0) == 1


def test_same_kernel_requested_from_two_threads_is_built_once(monkeypatch):
    import threading

    _small_kernel_settings(monkeypatch)
    built = {"n": 0}
    orig = engine._build_kernel_uncached

    def counted(*a, **k):
        built["n"] += 1
        return orig(*a, **k)

    monkeypatch.setattr(engine, "_build_kernel_uncached", counted)
    args = dict(sw_Hz=150000.0, npts=128, ref_offset_ppm=30.0,
                cq_max_MHz=6.0, n_cq=5, n_eta=3)
    out = {}

    def go(tag):
        out[tag] = engine.build_kernel("27Al", 130.32, 12500.0, **args)

    ts = [threading.Thread(target=go, args=(t,)) for t in ("warm", "sim")]
    for t in ts:
        t.start()
    for t in ts:
        t.join(60)
    assert built["n"] == 1                        # the race built it twice
    assert out["warm"] is out["sim"]


def test_waiting_for_another_build_still_honours_stop(monkeypatch):
    import threading
    import time

    _small_kernel_settings(monkeypatch)
    args = dict(sw_Hz=150000.0, npts=128, ref_offset_ppm=30.0,
                cq_max_MHz=6.0, n_cq=5, n_eta=3)
    key = ("27Al", 130.32, 12500, 150000, 128, 30.0, 6.0, 5, 3)
    lock = engine._build_lock(key)
    lock.acquire()                                # "another thread is building"
    try:
        flag = {"stop": False}
        threading.Timer(0.3, lambda: flag.__setitem__("stop", True)).start()
        t0 = time.perf_counter()
        with engine.kernel_build_feedback(cancel=lambda: flag["stop"]):
            with pytest.raises(engine.KernelBuildCancelled):
                engine.build_kernel("27Al", 130.32, 12500.0, **args)
        assert time.perf_counter() - t0 < 1.5
    finally:
        lock.release()


def _six_slow_lines(sleep_s: float, monkeypatch):
    """Six gauss_lor sites whose render sleeps: one residual evaluation costs
    6 x sleep_s, the unit of Stop latency before this change."""
    import dataclasses
    import time

    from larmor.models import base as mbase

    model = mbase.REGISTRY["gauss_lor"]
    orig = model.render

    def slow(values, ctx):
        time.sleep(sleep_s)
        return orig(values, ctx)

    monkeypatch.setitem(mbase.REGISTRY, "gauss_lor",
                        dataclasses.replace(model, render=slow))
    x = np.linspace(-20.0, 120.0, 300)
    sites = [SiteModel(model="gauss_lor", label=f"p{i}", params={
        "isotropic_chemical_shift_ppm": Param(20.0 + 15.0 * i, min=-50, max=150),
        "shift_fwhm_ppm": Param(6.0, min=1, max=40),
        "gl": Param(0.5, vary=False),
        "amplitude": Param(1e6, min=0)}) for i in range(6)]
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=sites)
    ctx = engine.make_context(r, exp_ppm=x)
    y = np.sum([engine.simulate_site(s, ctx) for s in r.sites], axis=0)
    return r, x, y * 1.3


def test_stop_between_site_renders_returns_within_a_second_and_keeps_the_last_point(monkeypatch):
    import threading
    import time

    r, x, y = _six_slow_lines(0.08, monkeypatch)        # 0.48 s per evaluation
    flag = {"stop": False, "t_req": None}
    seen = []

    def arm():
        time.sleep(0.15)
        flag["t_req"] = time.perf_counter()
        flag["stop"] = True

    def cb(params, it, resid, *a, **k):
        # the parameter set of every COMPLETED evaluation
        seen.append({n: p.value for n, p in params.items() if p.vary})
        if len(seen) == 3:
            # ask for the stop 0.2 s into the NEXT evaluation (mid-render)
            threading.Thread(target=arm).start()
        return flag["stop"] or None

    with engine.kernel_build_feedback(cancel=lambda: flag["stop"]):
        res = fitmod.fit(r, x, y, window_ppm=(120.0, -20.0), iter_cb=cb)
    t_back = time.perf_counter()
    assert flag["t_req"] is not None
    # one render (0.08 s) to notice + the final full-grid render of the kept
    # point (six renders, 0.48 s); before, a whole evaluation had to finish
    assert t_back - flag["t_req"] < 1.0
    assert res.lmfit_result.aborted
    # the recipe holds the last COMPLETE evaluation's values, not the trial
    # point whose evaluation was abandoned
    last = seen[-1]
    for i, s in enumerate(r.sites):
        for pname, p in s.params.items():
            key = fitmod._lmfit_name(i, s, pname)
            if key in last:
                assert p.value == pytest.approx(last[key], rel=1e-12)
    assert np.isfinite(res.y_fit).all() and res.rmsd < 1.0
