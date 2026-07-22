"""2D MQMAS engine: kernel, simulation, shear and fitting."""
import numpy as np
import pytest

from larmor import twod
from larmor.recipe import Param, Recipe, SiteModel

from conftest import CAALGLASS_MQ, require


def test_hypercomplex_phasing_is_exact():
    """With the imaginary quadrant present, a p0 rotation exactly recovers the
    absorption spectrum (Hilbert can only approximate this)."""
    n2, n1 = 256, 48
    f2 = np.linspace(-100, 100, n2); f1 = np.linspace(-50, 50, n1)
    x = (f2 - 10) / 2.0
    absorption = (1 / (1 + x ** 2))[None, :] * np.ones((n1, 1))
    dispersion = (x / (1 + x ** 2))[None, :] * np.ones((n1, 1))
    th = np.deg2rad(50.0)
    rr = absorption * np.cos(th) - dispersion * np.sin(th)
    ri = absorption * np.sin(th) + dispersion * np.cos(th)
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=rr, ri=ri,
                    ir=np.zeros_like(rr), ii=np.zeros_like(rr))
    assert d.has_hyper
    fixed = d.phased("f2", 50.0, 0.0, pivot_ppm=10.0)
    assert np.max(np.abs(fixed.z - absorption)) < 1e-9
    # phase_line preview matches the applied row exactly
    line = d.phase_line("f2", n1 // 2, 50.0, 0.0, 10.0)
    assert np.allclose(line, fixed.z[n1 // 2], atol=1e-9)


def _small_kernel():
    """A deliberately small kernel: these tests check physics and plumbing,
    not resolution."""
    return twod.build_mqmas_kernel(
        "27Al", 195.483, f2_window=(150.0, -100.0), f1_window=(120.0, -60.0),
        n2=64, n1=48, n_cq=10, n_eta=3, cq_max_MHz=12.0)


def _czjzek_recipe(sigma=2.0, pos=60.0):
    return Recipe(nucleus="27Al", larmor_frequency_MHz=195.483,
                  spin_rate_Hz=0.0,
                  sites=[SiteModel(model="czjzek", label="AlO4", params={
                      "isotropic_chemical_shift_ppm": Param(pos, min=0, max=120),
                      "sigma_Cq_MHz": Param(sigma, min=0.2, max=8.0),
                      "shift_fwhm_ppm": Param(6.0, min=1.0, max=30.0),
                      "amplitude": Param(1.0, min=0.0)})])


# ---------------------------------------------------------------- kernel
@pytest.mark.slow
def test_kernel_shape_and_caching():
    k = _small_kernel()
    assert k.K.shape == (10 * 3, 48, 64)
    assert k.shape == (48, 64)
    assert np.all(np.diff(k.f2_ppm) > 0) and np.all(np.diff(k.f1_ppm) > 0)
    assert np.isfinite(k.K).all()
    # cached: the second call must be the very same object
    assert _small_kernel() is k


def test_unknown_method_rejected():
    with pytest.raises(ValueError, match="unknown 2D method"):
        twod.build_mqmas_kernel("27Al", 195.5, (100, -100), (100, -100),
                                method="17QMAS")


@pytest.mark.slow
def test_kernel_weights_normalize():
    k = _small_kernel()
    w = k.weights(2.0)
    assert w.shape == (30,)
    assert w.sum() == pytest.approx(1.0)
    assert (w >= 0).all()
    we = k.ext_weights(5.0, 0.2, 0.3)
    assert we.sum() == pytest.approx(1.0)


# ---------------------------------------------------------------- physics
@pytest.mark.slow
def test_mqmas_peak_lies_off_the_diagonal_by_the_quadrupolar_shift():
    """The whole point of MQMAS: a quadrupolar site is displaced from the
    F1 = F2 diagonal, and MORE so when Cq (here sigma) grows."""
    k = _small_kernel()

    def peak_of(sigma):
        r = _czjzek_recipe(sigma=sigma, pos=60.0)
        z, _ = twod.simulate_2d(r, k)
        i1, i2 = np.unravel_index(int(np.argmax(z)), z.shape)
        return k.f1_ppm[i1], k.f2_ppm[i2]

    f1_small, f2_small = peak_of(0.5)
    f1_big, f2_big = peak_of(4.0)
    # a small-Cq site sits close to its isotropic shift on both axes
    assert abs(f1_small - 60.0) < 25.0
    # a large-Cq site is pushed further from the diagonal
    assert abs(f1_big - f2_big) > abs(f1_small - f2_small)


@pytest.mark.slow
def test_wider_czjzek_broadens_the_2d_footprint():
    k = _small_kernel()
    def area(sigma):
        z, _ = twod.simulate_2d(_czjzek_recipe(sigma=sigma), k)
        return int((z > 0.5 * z.max()).sum())
    assert area(4.0) > area(1.0)


def _f2f1_correlation(z, k):
    """Intensity-weighted Pearson correlation between the F2 and F1 axes."""
    F2, F1 = np.meshgrid(k.f2_ppm, k.f1_ppm)      # z is (n1, n2)
    p = z / z.sum()
    m2 = (p * F2).sum(); m1 = (p * F1).sum()
    v2 = (p * (F2 - m2) ** 2).sum(); v1 = (p * (F1 - m1) ** 2).sum()
    cov = (p * (F2 - m2) * (F1 - m1)).sum()
    return cov / np.sqrt(v2 * v1)


@pytest.mark.slow
def test_dCS_distribution_elongates_the_peak_along_the_diagonal():
    """dmfit's dCS: a distribution of isotropic shifts moves F2 and F1 together,
    so growing it must raise the F2-F1 correlation toward the +1 diagonal."""
    k = _small_kernel()

    def corr(dcs):
        r = _czjzek_recipe(sigma=1.0, pos=50.0)
        r.sites[0].params["shift_fwhm_ppm"].value = dcs
        z, _ = twod.simulate_2d(r, k)
        return _f2f1_correlation(z, k)

    assert corr(25.0) > corr(3.0) + 0.1
    assert corr(25.0) > 0.6                       # strongly diagonal


@pytest.mark.slow
def test_line_broadening_is_round_not_diagonal():
    """The round point broadening (dmfit wid) must NOT elongate along the
    diagonal the way the CS distribution does -- that is the whole reason they
    are separate parameters."""
    k = _small_kernel()

    def corr_with(dcs, line):
        r = _czjzek_recipe(sigma=1.0, pos=50.0)
        r.sites[0].params["shift_fwhm_ppm"].value = dcs
        r.sites[0].params["line_fwhm_ppm"] = Param(line, min=0.0)
        z, _ = twod.simulate_2d(r, k)
        return _f2f1_correlation(z, k)

    # same tiny CS distribution; adding round broadening barely moves the
    # correlation, while the CS distribution drives it strongly positive
    round_corr = corr_with(2.0, 20.0)
    diag_corr = corr_with(20.0, 0.0)
    assert diag_corr > round_corr + 0.2


def test_czjzek_1d_line_broadening_adds_in_quadrature():
    """1D: dCS and the round line width both blur the single MAS dimension, so
    they combine in quadrature -- and line=0 reproduces the legacy width."""
    from larmor.models.quadrupolar import _czjzek_fwhm

    assert _czjzek_fwhm({"shift_fwhm_ppm": 10.0}) == pytest.approx(10.0)
    assert _czjzek_fwhm({"shift_fwhm_ppm": 10.0, "line_fwhm_ppm": 0.0}) == \
        pytest.approx(10.0)
    assert _czjzek_fwhm({"shift_fwhm_ppm": 8.0, "line_fwhm_ppm": 6.0}) == \
        pytest.approx(10.0)


@pytest.mark.slow
def test_amplitude_scales_linearly():
    k = _small_kernel()
    r1 = _czjzek_recipe()
    z1, _ = twod.simulate_2d(r1, k)
    r2 = _czjzek_recipe()
    r2.sites[0].params["amplitude"].value = 3.0
    z2, _ = twod.simulate_2d(r2, k)
    assert z2.max() == pytest.approx(3.0 * z1.max(), rel=1e-6)


@pytest.mark.slow
def test_unsupported_model_in_2d_is_rejected():
    k = _small_kernel()
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=195.5,
               sites=[SiteModel(model="gauss_lor", params={
                   "isotropic_chemical_shift_ppm": Param(0.0),
                   "shift_fwhm_ppm": Param(5.0), "amplitude": Param(1.0),
                   "gl": Param(1.0)})])
    with pytest.raises(ValueError, match="no 2D implementation"):
        twod.simulate_2d(r, k)


# ---------------------------------------------------------------- data ops
def test_data2d_region_and_projection():
    f2 = np.linspace(-50, 50, 40)
    f1 = np.linspace(-20, 20, 30)
    z = np.zeros((30, 40))
    z[15, 20] = 1.0
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=z)
    r = d.region(f2_range=(-10, 10), f1_range=(-5, 5))
    assert r.z.shape[0] < 30 and r.z.shape[1] < 40
    assert d.projection("f2").shape == (40,)
    assert d.projection("f1").shape == (30,)
    assert d.projection("f2").max() == pytest.approx(1.0)
    assert d.normalized().z.max() == pytest.approx(1.0)


def test_shear_moves_intensity_predictably():
    f2 = np.linspace(-50, 50, 101)
    f1 = np.linspace(-50, 50, 101)
    z = np.zeros((101, 101))
    j = 75                       # f2 = +25 ppm
    i = 50                       # f1 = 0 ppm
    z[i, j] = 1.0
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=z)
    sh = twod.shear(d, factor=1.0, ref_ppm=0.0)
    # F1' = F1 + 1.0 * (F2 - 0) -> the point moves up by 25 ppm in F1
    i_new = int(np.argmax(sh.z[:, j]))
    assert f1[i_new] == pytest.approx(25.0, abs=1.5)
    assert any("sheared" in n for n in sh.notes)


# ---------------------------------------------------------------- fit
@pytest.mark.slow
def test_fit_2d_recovers_a_synthetic_czjzek_site():
    k = _small_kernel()
    truth = _czjzek_recipe(sigma=2.5, pos=65.0)
    z_true, _ = twod.simulate_2d(truth, k)
    rng = np.random.default_rng(0)
    data = twod.Data2D(f2_ppm=k.f2_ppm, f1_ppm=k.f1_ppm,
                       z=z_true + rng.normal(0, 0.004, z_true.shape),
                       nucleus="27Al", larmor_MHz=195.483)

    start = _czjzek_recipe(sigma=1.2, pos=52.0)      # deliberately wrong
    res = twod.fit_2d(start, data, kernel=k)
    p = start.sites[0].params
    assert res.rmsd < 0.05
    assert p["sigma_Cq_MHz"].value == pytest.approx(2.5, abs=0.5)
    assert p["isotropic_chemical_shift_ppm"].value == pytest.approx(65.0, abs=6.0)
    assert res.z_fit.shape == k.shape
    assert len(res.per_site) == 1


@pytest.mark.slow
def test_fit_2d_recovers_the_f1_reference_offset():
    """mrsimulator's kernel and an experimental F1 axis differ by a referencing
    offset; fit_2d fits it (isotropic-axis alignment) and holds it when fixed."""
    k = _small_kernel()
    truth = _czjzek_recipe(sigma=2.0, pos=55.0)
    truth.mqmas_f1_ref_ppm = 15.0                    # data carries a +15 ppm F1 shift
    z, _ = twod.simulate_2d(truth, k)
    data = twod.Data2D(f2_ppm=k.f2_ppm, f1_ppm=k.f1_ppm, z=z,
                       nucleus="27Al", larmor_MHz=195.483)

    auto = _czjzek_recipe(sigma=2.0, pos=55.0)        # ref unknown -> fitted
    twod.fit_2d(auto, data, kernel=k)
    assert auto.mqmas_f1_ref_ppm == pytest.approx(15.0, abs=3.0)

    held = _czjzek_recipe(sigma=2.0, pos=55.0)        # user fixes the referencing
    held.mqmas_f1_ref_ppm = 5.0
    held.mqmas_f1_ref_vary = False
    twod.fit_2d(held, data, kernel=k)
    assert held.mqmas_f1_ref_ppm == 5.0              # untouched by the fit


@pytest.mark.slow
def test_fit_2d_reports_uncertainties():
    k = _small_kernel()
    truth = _czjzek_recipe(sigma=2.0, pos=60.0)
    z, _ = twod.simulate_2d(truth, k)
    rng = np.random.default_rng(1)
    data = twod.Data2D(f2_ppm=k.f2_ppm, f1_ppm=k.f1_ppm,
                       z=z + rng.normal(0, 0.01, z.shape),
                       nucleus="27Al", larmor_MHz=195.483)
    r = _czjzek_recipe(sigma=1.6, pos=57.0)
    res = twod.fit_2d(r, data, kernel=k)
    assert res.lmfit_result.errorbars
    assert r.sites[0].params["sigma_Cq_MHz"].stderr is not None
    assert "sigma" in res.report or "s0_sigma" in res.report


# ---------------------------------------------------------------- real file
def test_caalglass_mq_parses_as_2d():
    """The user's real MQMAS dmfit file: it must be recognized as 2D with
    its Czjzek sites (fitting it needs the full-size kernel; this asserts the
    import path the 2D engine consumes)."""
    from larmor.io import fxmla

    dm = fxmla.read(require(CAALGLASS_MQ))
    assert dm.is_2d
    assert dm.fit_mode == "MQMAS"
    dim = dm.dimensions[0]
    assert dim.nucleus == "27Al"
    assert all(ln.model_name == "CzSimple" for ln in dim.lines)
    recipe, warnings = fxmla.to_recipe(dm)
    assert [s.model for s in recipe.sites] == ["czjzek", "czjzek"]
    assert any("MQMAS" in w for w in warnings)


def test_2d_operations():
    """Transpose, reverse, and diagonal extraction (dmfit 2D ops)."""
    f2 = np.linspace(-40, 80, 60); f1 = np.linspace(20, 60, 30)
    Z = (np.exp(-((f2[None, :] - 30) / 6) ** 2)
         * np.exp(-((f1[:, None] - 40) / 6) ** 2))
    d = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z, nucleus="27Al", larmor_MHz=156.0)
    dt = d.transposed()
    assert dt.z.shape == (60, 30) and np.allclose(dt.z, Z.T)
    assert np.allclose(dt.f2_ppm, f1) and np.allclose(dt.f1_ppm, f2)
    dr = d.reversed_axis("f2")
    assert np.allclose(dr.z, np.flip(Z, 1))
    g, amp = d.diagonal()
    assert g.size == 60 and amp.shape == g.shape


def test_2d_symmetrize():
    f2 = np.linspace(-40, 40, 50); f1 = np.linspace(-40, 40, 50)
    Z = (np.exp(-((f2[None, :] - 10) / 6) ** 2)
         * np.exp(-((f1[:, None] + 10) / 6) ** 2))
    ds = twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=Z).symmetrized()
    assert np.max(np.abs(ds.z - ds.z.T)) < 1e-9            # symmetric now
