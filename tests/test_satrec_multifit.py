"""Tests for the satrec T1 module, the multi-dataset fit engine, and the
advanced processing operations."""
import numpy as np
import pytest

from larmor import processing, satrec
from larmor.recipe import Param, Recipe, SiteModel

from conftest import EXPNO_1901, require


# ---------------------------------------------------------------- processing
def _fid(n=1024, sw=10000.0, freq=1000.0, t2=0.02):
    t = np.arange(n) / sw
    return processing.Spectrum1D(
        x_ppm=None, y=np.exp(2j * np.pi * freq * t) * np.exp(-t / t2),
        sfo1_MHz=100.0, sw_Hz=sw, domain="time")


@pytest.mark.parametrize("ops", [
    [{"op": "gm", "lb_hz": -20, "gb": 0.3}],
    [{"op": "sine", "ssb": 2}],
    [{"op": "sine", "ssb": 3, "power": 2}],
    [{"op": "traf", "lb_hz": 15}],
    [{"op": "fcor", "factor": 0.5}],
])
def test_windows_keep_peak_position(ops):
    s = processing.apply(_fid(), ops + [{"op": "zf"}, {"op": "ft"}])
    peak = s.x_ppm[np.argmax(np.abs(s.y))]
    assert abs(abs(peak) - 10.0) < 0.3      # 1000 Hz / 100 MHz = 10 ppm


def test_tdeff_truncates():
    s = _fid()
    processing.apply(s, [{"op": "tdeff", "points": 256}])
    assert s.y.size == 256


def test_sr_shifts_axis():
    s = processing.apply(_fid(), [{"op": "ft"}])
    x0 = s.x_ppm[np.argmax(np.abs(s.y))]
    processing.op_sr(s, sr_hz=200.0)        # 200 Hz @ 100 MHz = 2 ppm
    x1 = s.x_ppm[np.argmax(np.abs(s.y))]
    assert x1 - x0 == pytest.approx(2.0, abs=1e-9)


def test_hilbert_enables_rephasing():
    # fcor 0.5 first: an uncorrected first fid point leaves a DC offset that
    # has no Kramers-Kronig partner (that is exactly what TopSpin FCOR fixes)
    ref = processing.apply(_fid(), [{"op": "fcor", "factor": 0.5},
                                    {"op": "ft"}])
    a_max = ref.y.real.max()                # true absorption maximum

    s = processing.apply(_fid(), [{"op": "fcor", "factor": 0.5},
                                  {"op": "ft"}, {"op": "phase", "p0": 70.0}])
    s.y = s.y.real + 0j                     # keep only the real part (like 1r)
    processing.apply(s, [{"op": "hilbert"}, {"op": "autophase"}])
    # note: real.max/|y|.max < 1 even for a perfectly phased off-bin
    # Lorentzian (|y| contains dispersion), so compare to the never-dephased
    # absorption instead
    assert s.y.real.max() > 0.93 * a_max


def test_magnitude():
    s = processing.apply(_fid(), [{"op": "ft"}, {"op": "phase", "p0": 90.0},
                                  {"op": "magnitude"}])
    assert (s.y.real >= 0).all()


# ---------------------------------------------------------------- satrec
@pytest.mark.slow
def test_satrec_real_expno_vs_topspin():
    """Full auto-T1 on the real 19F satrec EXPNO; integrals must correlate
    with TopSpin's own t1ints.txt values."""
    path = require(EXPNO_1901)
    result = satrec.analyze(path, lb_hz=100.0)
    assert result.t1_s > 0
    assert result.integrals.max() == pytest.approx(1.0)
    assert np.all(np.diff(result.delays_s) > 0)
    # monotonic-ish build-up: last integral is the largest
    assert result.integrals[-1] == pytest.approx(result.integrals.max(),
                                                 abs=0.05)

    ts_file = path / "pdata" / "1" / "t1ints.txt"
    if ts_file.exists():
        vals = [float(v) for v in ts_file.read_text().split()]
        ts = []
        i = 1
        while i + 9 <= len(vals) + 1 and vals[i] >= 0:
            ts.append(vals[i + 7])
            i += 9
        ts = np.abs(np.array(ts))
        ts = ts / ts.max()
        n = min(len(ts), len(result.integrals))
        r = np.corrcoef(ts[:n], np.abs(result.integrals[:n]))[0, 1]
        assert r > 0.98, f"correlation with TopSpin integrals too low: {r}"


def test_vdlist_units(tmp_path):
    (tmp_path / "vdlist").write_text("0\n250m\n1\n4s\n100u\n")
    d = satrec.read_vdlist(tmp_path)
    assert d == pytest.approx([0.0, 0.25, 1.0, 4.0, 1e-4])


def test_fit_t1_synthetic():
    t = np.array([0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
    y = 1.0 - np.exp(-t / 2.5)
    out, curve = satrec.fit_t1(t, y)
    assert out.params["t1"].value == pytest.approx(2.5, rel=0.02)
    assert curve(2.5) == pytest.approx(1 - np.exp(-1), rel=0.02)


# ---------------------------------------------------------------- multifit
def _recipe_two_fields(larmor_MHz):
    return Recipe(
        nucleus="27Al", larmor_frequency_MHz=larmor_MHz, spin_rate_Hz=20000.0,
        sites=[SiteModel(model="quad_ct", label="site", params={
            "isotropic_chemical_shift_ppm": Param(52.0),
            "Cq_MHz": Param(4.0, min=0.5, max=12.0),
            "eta": Param(0.15, vary=False),
            "shift_fwhm_ppm": Param(3.0, min=0.5),
            "amplitude": Param(1.0, min=0.0)})],
    )


@pytest.mark.slow
def test_multifield_lifts_degeneracy():
    """Simulate ONE 27Al site at two fields, start the fit from wrong values,
    and require the shared Cq/delta_iso to be recovered."""
    from larmor.engine import make_context, simulate_site

    truth = {"pos": 55.0, "cq": 5.2}
    entries = []
    rng = np.random.default_rng(1)
    for larmor in (104.26, 195.48):          # 9.4 T and 17.6 T
        r = _recipe_two_fields(larmor)
        r.sites[0].params["isotropic_chemical_shift_ppm"].value = truth["pos"]
        r.sites[0].params["Cq_MHz"].value = truth["cq"]
        ctx = make_context(r, exp_ppm=np.linspace(-60, 120, 1500))
        y = simulate_site(r.sites[0], ctx) + rng.normal(0, 0.004, ctx.x_ppm.size)
        # reset to a wrong start
        r.sites[0].params["isotropic_chemical_shift_ppm"].value = 48.0
        r.sites[0].params["Cq_MHz"].value = 3.0
        entries.append((r, ctx.x_ppm.copy(), y))

    from larmor.multifit import fit_multi

    result = fit_multi(entries,
                       share=("isotropic_chemical_shift_ppm", "Cq_MHz", "eta"))
    r0, r1 = result.recipes
    p0 = r0.sites[0].params
    assert p0["Cq_MHz"].value == pytest.approx(truth["cq"], abs=0.15)
    assert p0["isotropic_chemical_shift_ppm"].value == pytest.approx(
        truth["pos"], abs=1.0)
    # the shared parameters are identical across datasets (linked, not copied)
    assert r1.sites[0].params["Cq_MHz"].value == pytest.approx(
        p0["Cq_MHz"].value, rel=1e-9)
    assert max(result.rmsd) < 0.05


@pytest.mark.slow
def test_cofit_mixed_1d_and_2d_mqmas():
    """Co-fit a 1D spectrum and a 2D MQMAS map of the SAME czjzek site with a
    shared sigma_Cq; recover it from a wrong start."""
    from larmor import engine, twod
    from larmor.multifit import fit_cofit

    def czjzek(sigma):
        return Recipe(nucleus="27Al", larmor_frequency_MHz=195.5,
                      sites=[SiteModel(model="czjzek", label="A", params={
                          "isotropic_chemical_shift_ppm": Param(0.0),
                          "sigma_Cq_MHz": Param(sigma, min=0.1, max=8.0),
                          "shift_fwhm_ppm": Param(2.0, min=0.1),
                          "amplitude": Param(1.0, min=0.0)})])

    x = np.linspace(-80, 120, 500)
    xc, tot, _ = engine.simulate(czjzek(3.0), exp_ppm=x)
    y1 = np.interp(x, xc, tot)
    k = twod.build_mqmas_kernel("27Al", 195.5, (120, -80), (60, -40))
    zt, _ = twod.simulate_2d(czjzek(3.0), k)
    d2 = twod.Data2D(f2_ppm=k.f2_ppm, f1_ppm=k.f1_ppm, z=zt, nucleus="27Al",
                     larmor_MHz=195.5)

    res = fit_cofit([(czjzek(1.5), (x, y1)), (czjzek(1.5), d2)])
    assert [pd["kind"] for pd in res.per_dataset] == ["1d", "2d"]
    s0 = res.recipes[0].sites[0].params["sigma_Cq_MHz"].value
    assert s0 == pytest.approx(3.0, abs=0.1)
    assert res.recipes[1].sites[0].params["sigma_Cq_MHz"].value == pytest.approx(
        s0, rel=1e-9)                       # tied across the 1D and 2D
    assert max(res.rmsd) < 0.02


def test_multifit_validates_alignment():
    from larmor.multifit import fit_multi

    r1 = _recipe_two_fields(104.26)
    r2 = _recipe_two_fields(195.48)
    r2.sites[0] = SiteModel(model="gauss_lor", params={
        "isotropic_chemical_shift_ppm": Param(0.0),
        "shift_fwhm_ppm": Param(1.0), "amplitude": Param(1.0),
        "gl": Param(1.0)})
    x = np.linspace(-10, 10, 100)
    with pytest.raises(ValueError, match="model differs"):
        fit_multi([(r1, x, np.zeros(100)), (r2, x, np.zeros(100))])