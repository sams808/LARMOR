"""Relaxation series suite (satrec/invrec/cpmg, per-site) and REDOR."""
import numpy as np
import pytest

from larmor import redor, series
from larmor.recipe import Param, Recipe, SiteModel


# ---------------------------------------------------------------- series fits
def test_fit_series_satrec():
    t = np.array([0, 0.05, 0.2, 0.5, 1.0, 3.0, 10.0, 30.0])
    y = 1 - np.exp(-t / 2.0)
    out, curve = series.fit_series(t, y, kind="satrec")
    assert out.params["tau"].value == pytest.approx(2.0, rel=0.02)


def test_fit_series_invrec_crosses_zero():
    t = np.array([0.01, 0.5, 1.0, 1.386, 2.0, 5.0, 20.0])
    y = 1 - 2 * np.exp(-t / 2.0)             # zero crossing at T1*ln2
    out, curve = series.fit_series(t, y, kind="invrec")
    assert out.params["tau"].value == pytest.approx(2.0, rel=0.05)
    assert curve(2.0 * np.log(2)) == pytest.approx(0.0, abs=0.05)


def test_fit_series_cpmg_decay():
    t = np.linspace(0.001, 0.05, 12)
    y = np.exp(-t / 0.01)
    out, _ = series.fit_series(t, y, kind="cpmg")
    assert out.params["tau"].value == pytest.approx(0.01, rel=0.03)


def test_fit_series_stretched():
    t = np.logspace(-2, 1.5, 15)
    y = 1 - np.exp(-((t / 3.0) ** 0.6))
    out, _ = series.fit_series(t, y, kind="satrec", stretched=True)
    assert out.params["tau"].value == pytest.approx(3.0, rel=0.1)
    assert out.params["beta"].value == pytest.approx(0.6, abs=0.06)


def test_unknown_kind_rejected():
    import lmfit

    p = lmfit.Parameters()
    p.add("i0", value=1)
    p.add("tau", value=1)
    with pytest.raises(ValueError, match="unknown series kind"):
        series.model_curve("sorcery", p, np.array([1.0]), False)


def test_read_delays_vclist(tmp_path):
    (tmp_path / "vclist").write_text("2\n4\n8\n16\n")
    vals, src = series.read_delays(tmp_path)
    assert src == "vclist"
    assert vals == pytest.approx([2, 4, 8, 16])


# ---------------------------------------------------------------- per site
def test_per_site_separates_two_overlapping_t1():
    """THE point of per-site analysis: two overlapping lines with different
    T1. A single integration window would return one meaningless average."""
    from larmor.engine import make_context, simulate_site

    ppm = np.linspace(-40, 40, 1500)
    recipe = Recipe(nucleus="27Al", larmor_frequency_MHz=195.5, spin_rate_Hz=0,
                    sites=[
        SiteModel(model="gauss_lor", label="fast", params={
            "isotropic_chemical_shift_ppm": Param(6.0),
            "shift_fwhm_ppm": Param(9.0), "amplitude": Param(1.0),
            "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="slow", params={
            "isotropic_chemical_shift_ppm": Param(-6.0),
            "shift_fwhm_ppm": Param(9.0), "amplitude": Param(1.0),
            "gl": Param(1.0, vary=False)}),
    ])
    ctx = make_context(recipe, exp_ppm=ppm)
    base = [simulate_site(s, ctx) for s in recipe.sites]

    delays = np.array([0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 15.0, 40.0])
    t1_true = (0.5, 8.0)
    rng = np.random.default_rng(0)
    slices = np.empty((delays.size, ppm.size))
    for k, t in enumerate(delays):
        a0 = 1 - np.exp(-t / t1_true[0])
        a1 = 1 - np.exp(-t / t1_true[1])
        slices[k] = a0 * base[0] + a1 * base[1] + rng.normal(0, 0.002, ppm.size)

    amps = series.fit_slice_amplitudes(recipe, ppm, slices)
    assert amps.shape == (delays.size, 2)
    results = []
    for i in range(2):
        y = amps[:, i] / amps[:, i].max()
        out, _ = series.fit_series(delays, y, kind="satrec")
        results.append(out.params["tau"].value)
    assert results[0] == pytest.approx(t1_true[0], rel=0.15)
    assert results[1] == pytest.approx(t1_true[1], rel=0.15)
    # and the two are cleanly distinguished
    assert results[1] / results[0] == pytest.approx(16.0, rel=0.3)


def test_fit_slice_amplitudes_is_nonnegative():
    ppm = np.linspace(-20, 20, 500)
    recipe = Recipe(nucleus="27Al", larmor_frequency_MHz=195.5, spin_rate_Hz=0,
                    sites=[SiteModel(model="gauss_lor", params={
                        "isotropic_chemical_shift_ppm": Param(0.0),
                        "shift_fwhm_ppm": Param(5.0), "amplitude": Param(1.0),
                        "gl": Param(1.0, vary=False)})])
    slices = np.random.default_rng(2).normal(0, 1, (3, ppm.size))  # pure noise
    amps = series.fit_slice_amplitudes(recipe, ppm, slices)
    assert (amps >= 0).all()


# ---------------------------------------------------------------- REDOR
def test_redor_universal_curve_landmarks():
    """The powder curve must reproduce the universal REDOR curve of an
    isolated pair (Mueller et al. 1995): ΔS/S0 = 0.5 near lambda = D*N*Tr =
    0.77, a first maximum near lambda = 1.67 with the characteristic
    overshoot (~1.04)."""
    lam = np.linspace(0.001, 2.5, 600)
    ds = redor.redor_pair_curve(1.0, lam)          # D = 1 Hz -> lambda = ntr
    i_max = int(np.argmax(ds))
    assert 1.55 < lam[i_max] < 1.8, f"first max at lambda={lam[i_max]:.2f}"
    assert 1.0 < ds[i_max] < 1.1                   # characteristic overshoot
    half = lam[int(np.argmin(np.abs(ds[:i_max] - 0.5)))]
    assert 0.72 < half < 0.82, f"dS/S0 = 0.5 at lambda={half:.2f}"
    assert ds[0] < 1e-3                            # starts at zero


def test_redor_curve_is_the_exact_isolated_pair_solution():
    """Physics from outside the module: the powder curve equals the closed
    form of Mueller et al. (J. Magn. Reson. A 113, 81, 1995), S/S0 =
    (sqrt2 pi/4) J_1/4(sqrt2 lambda) J_-1/4(sqrt2 lambda); a 13C-15N pair at
    the glycine C'-N distance (2.48 A), its dephasing made from that closed
    form and the CODATA / IUPAC constants, comes back at 2.48 A; the M2 is
    the van Vleck value of one spin-1/2 partner, which Bertmer & Eckert's
    short-time law (4/(3 pi^2)) M2 (N Tr)^2 turns back into the data. The
    curve was once twice too fast in lambda (the I-spin line at
    +-D (3cos^2 - 1) instead of +-(D/2)): every fitted D came out halved,
    every distance 26 % long, and the round-trip tests could not see it."""
    from scipy.special import jv

    def mueller(lam):
        x = np.sqrt(2.0) * np.asarray(lam, float)
        with np.errstate(invalid="ignore"):
            s = (np.sqrt(2.0) * np.pi / 4.0) * jv(0.25, x) * jv(-0.25, x)
        return np.where(x > 0, 1.0 - s, 0.0)

    lam = np.linspace(0.0, 3.0, 301)
    assert np.abs(redor.redor_pair_curve(1.0, lam) - mueller(lam)).max() < 2e-3
    hbar, mu0_4pi = 1.054571817e-34, 1.0e-7
    g13c, g15n = 6.728284e7, -2.712618e7          # rad s^-1 T^-1
    d_true = mu0_4pi * abs(g13c * g15n) * hbar / (2 * np.pi * (2.48e-10) ** 3)   # ~201 Hz
    ntr = np.arange(2, 40, 2) / 10000.0           # 10 kHz MAS, lambda up to 0.76
    ds = mueller(d_true * ntr)
    pair = redor.analyze(ntr, ds, pair=("13C", "15N"), regime="pair")
    assert pair.d_hz == pytest.approx(d_true, rel=0.01)
    assert pair.distance_A == pytest.approx(2.48, abs=0.01)
    short = redor.analyze(ntr, ds, pair=("13C", "15N"), regime="short")
    assert short.d_hz == pytest.approx(d_true, rel=0.05)
    assert pair.m2 == pytest.approx(0.75 * (4.0 / 15.0) * (2 * np.pi * d_true) ** 2, rel=0.02)
    low = ds < 0.05
    assert np.allclose(4.0 / (3.0 * np.pi ** 2) * pair.m2 * ntr[low] ** 2, ds[low], rtol=0.03)


def test_redor_parabola_is_the_small_lambda_limit():
    """The parabola is the second-order expansion of the powder average, so
    it must converge to it as lambda -> 0, and its error must GROW toward the
    quoted validity limit (dS/S0 ~ 0.2) -- that is why the limit exists."""
    lam = np.linspace(1e-3, 0.5, 120)
    exact = redor.redor_pair_curve(1.0, lam)
    para = redor.short_time_curve(1.0, lam)
    rel = np.abs(para - exact) / np.maximum(exact, 1e-12)

    # the decisive check: as lambda -> 0 the parabola IS the powder average
    assert rel[0] < 1e-3
    deep = exact < 0.05          # well inside the regime
    assert rel[deep].max() < 0.02
    edge = exact < 0.2           # at the quoted validity limit: ~8%
    assert rel[edge].max() < 0.085
    # the error grows smoothly with lambda -- the expansion degrades, which is
    # exactly why the 0.2 limit is quoted in the first place
    assert rel[-1] > rel[0]


def test_redor_dipolar_constant_and_distance_roundtrip():
    d = redor.dipolar_constant_hz("13C", "15N", 2.5)
    assert 100 < d < 400          # a 2.5 A C-N pair is a few hundred Hz
    r = redor.distance_angstrom("13C", "15N", d)
    assert r == pytest.approx(2.5, rel=1e-6)
    # 1/r^3 scaling
    assert (redor.dipolar_constant_hz("13C", "15N", 5.0)
            == pytest.approx(d / 8.0, rel=1e-6))


def test_redor_analyze_recovers_d_short_time():
    d_true = 200.0
    ntr = np.linspace(0.0002, 0.0016, 8)          # keeps dS/S0 below ~0.2
    ds = redor.short_time_curve(d_true, ntr)
    res = redor.analyze(ntr, ds, pair=("13C", "15N"), regime="short")
    assert res.d_hz == pytest.approx(d_true, rel=0.02)
    assert res.regime == "short-time parabola"
    assert res.distance_A == pytest.approx(
        redor.distance_angstrom("13C", "15N", d_true), rel=0.02)
    assert "model-free" in " ".join(res.notes)


def test_redor_analyze_full_pair_curve():
    d_true = 350.0
    ntr = np.linspace(0.0002, 0.004, 14)
    ds = redor.redor_pair_curve(d_true, ntr)
    res = redor.analyze(ntr, ds, regime="pair")
    assert res.d_hz == pytest.approx(d_true, rel=0.03)
    assert res.regime == "isolated pair"
    assert "isolated" in " ".join(res.notes)


def test_redor_auto_regime_prefers_short_when_available():
    d_true = 150.0
    ntr = np.linspace(0.0002, 0.0015, 8)
    ds = redor.short_time_curve(d_true, ntr)
    res = redor.analyze(ntr, ds, regime="auto")
    assert res.regime == "short-time parabola"
    assert any("auto-selected" in n for n in res.notes)


def test_redor_txt_parser(tmp_path):
    p = tmp_path / "redor.txt"
    p.write_text(
        "Dataset :\n/data/x\n\nSpinning speed : 35714.000000\n\n"
        "Peak 1 (xy -221.324 ppm)\n\n"
        "    number integral(S0) integral(S*) intensity(S0) intensity(S*)\n\n"
        "         1    -1.00e+10    -9.00e+09             0             0\n"
        "         2    -1.00e+10    -5.00e+09             0             0\n")
    n, ds, masr = redor.read_redor_txt(p)
    assert masr == pytest.approx(35714.0)
    assert n == pytest.approx([1, 2])
    assert ds == pytest.approx([0.1, 0.5])


def test_detect_kind_reads_the_pulse_program(tmp_path):
    """The series kind comes from acqus PULPROG, 1D or pseudo-2D. nmrglue's
    read_acqus_file returns ONE dict keyed by file name; a two-value
    unpacking of it once made detect_kind None for every EXPNO, so every
    kind=None caller (`larmor satrec`, the per-site T1 actions) fitted a
    saturation recovery whatever the experiment -- a CPMG T2 included."""
    def expno(name, prog, two_d=False):
        e = tmp_path / name
        e.mkdir()
        txt = f"##TITLE= Parameter file\n##JCAMPDX= 5.0\n##$PULPROG= <{prog}>\n##$TD= 1024\n##END=\n"
        (e / "acqus").write_text(txt)
        if two_d:
            (e / "acqu2s").write_text(txt)
        return e

    assert series.detect_kind(expno("1", "satrec_zg", two_d=True)) == "satrec"
    assert series.detect_kind(expno("2", "cpmg1d")) == "cpmg"
    assert series.detect_kind(expno("3", "t1ir")) == "invrec"
    assert series.detect_kind(expno("4", "zg")) is None            # unsure: ask, never guess
    assert series.detect_kind(tmp_path / "no-such-expno") is None
