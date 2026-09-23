"""New processing ops: subtract-averages, scale-SW, inverse FT, real/imag/conj."""
import numpy as np
import pytest

from larmor import processing as P


def test_ft_ift_roundtrip():
    fid = np.random.RandomState(0).randn(256) + 1j * np.random.RandomState(1).randn(256)
    s = P.Spectrum1D(x_ppm=None, y=fid.copy(), sfo1_MHz=100.0, sw_Hz=1e4, domain="time")
    s = P.op_ft(s)
    assert s.domain == "freq"
    s = P.op_ift(s)
    assert s.domain == "time" and s.x_ppm is None
    assert np.max(np.abs(s.y - fid)) < 1e-9


# ---------------------------------------------------------------- F6: ift of a
# processed spectrum (time <-> frequency toggle, re-apodization)
def _lorentzian(x, x0, fwhm):
    g = fwhm / 2.0
    return g ** 2 / ((x - x0) ** 2 + g ** 2)


def _fwhm(x, y):
    """Half-height width walked out from the maximum on the grid."""
    y = np.asarray(y, float)
    i = int(np.argmax(y))
    half = y[i] / 2.0
    li = i
    while li > 0 and y[li] > half:
        li -= 1
    ri = i
    while ri < y.size - 1 and y[ri] > half:
        ri += 1
    return abs(float(x[ri] - x[li]))


_X = np.linspace(-50.0, 50.0, 1024)
_DX = float(_X[1] - _X[0])
_Y = _lorentzian(_X, 12.3, 1.0)          # FWHM 1 ppm = 100 Hz at 100 MHz


def _fid(n=1024, sw=10000.0, f=1000.0, t2=0.02):
    t = np.arange(n) / sw
    y = np.exp(2j * np.pi * f * t) * np.exp(-t / t2)
    return P.Spectrum1D(x_ppm=None, y=y, sfo1_MHz=100.0, sw_Hz=sw,
                        domain="time")


def test_ift_derives_sw_and_ft_restores_the_axis_exactly():
    """from_processed leaves sw_Hz = 0; the old op_ift kept it (every window
    then divided by zero) and op_ft rebuilt the axis about 0 ppm."""
    s = P.from_processed(_X, _Y, 100.0)
    assert s.sw_Hz == 0.0
    s = P.apply(s, [{"op": "hilbert"}, {"op": "ift"}])
    assert s.domain == "time" and s.x_ppm is None
    assert np.allclose(s.x_ppm_hold, _X)
    assert s.sw_Hz == pytest.approx(_DX * 1024 * 100.0, rel=1e-9)
    s = P.apply(s, [{"op": "ft"}])
    assert s.domain == "freq"
    assert np.allclose(s.x_ppm, _X, atol=1e-9)
    assert np.allclose(s.y.real, _Y, atol=1e-8 * _Y.max())
    assert s.x_ppm_hold is None


def test_ift_of_a_processed_spectrum_without_larmor_frequency_is_a_clear_error():
    s = P.from_processed(_X, _Y, 0.0)
    with pytest.raises(ValueError, match="Larmor"):
        P.op_ift(s)


def test_ift_zf_ft_keeps_the_zero_frequency_bin_and_the_peak():
    """fftshift puts 0 Hz at index n//2 on both grids: the zero-filled axis is
    anchored on the held zero bin, not on the mid-span (half a bin off)."""
    s = P.from_processed(_X, _Y, 100.0)
    s = P.apply(s, [{"op": "hilbert"}, {"op": "ift"},
                    {"op": "zf", "factor": 2}, {"op": "ft"}])
    assert s.x_ppm.size == 2048
    assert s.x_ppm[1024] == pytest.approx(_X[512], abs=1e-12)
    d = np.diff(s.x_ppm)
    assert np.all(d > 0) and np.allclose(d, _DX / 2.0)
    assert s.x_ppm[int(np.argmax(s.y.real))] == pytest.approx(12.3, abs=_DX)


def test_reapodize_chain_broadens_by_lb_only_with_hilbert():
    """The panel forces Hilbert before ift: a real-only spectrum's ift is
    two-sided (hermitian), so a one-sided EM window damps the mirrored half
    and loses about half the signal (measured: 0.17 of a 0.33 peak)."""
    s = P.from_processed(_X, _Y, 100.0)
    s = P.apply(s, [{"op": "hilbert"}, {"op": "ift"}])
    e = np.abs(s.y) ** 2
    assert e[:256].sum() > 0.9 * e.sum()                  # one-sided fid
    assert abs(s.y[-1]) < 1e-6 * abs(s.y[0])
    s = P.apply(s, [{"op": "em", "lb_hz": 200}, {"op": "ft"}])
    # 1 ppm + 200 Hz / 100 MHz = 3 ppm; the grid step is 0.098 ppm
    assert _fwhm(s.x_ppm, s.y.real) == pytest.approx(3.0, abs=2 * _DX)
    peak_ok = float(s.y.real.max())
    assert peak_ok == pytest.approx(1.0 / 3.0, rel=0.1)   # area conserved

    bad = P.from_processed(_X, _Y, 100.0)
    bad = P.apply(bad, [{"op": "ift"}])
    assert abs(bad.y[-1]) > 0.5 * abs(bad.y[1])          # two-sided
    bad = P.apply(bad, [{"op": "em", "lb_hz": 200}, {"op": "ft"}])
    assert float(bad.y.real.max()) < 0.6 * peak_ok        # half the signal gone


def test_ift_accepts_a_descending_axis():
    s = P.from_processed(_X[::-1], _Y[::-1], 100.0)
    s = P.apply(s, [{"op": "hilbert"}, {"op": "ift"}, {"op": "ft"}])
    order = np.argsort(s.x_ppm)
    assert np.allclose(s.x_ppm[order], _X, atol=1e-9)
    assert np.allclose(s.y.real[order], _Y, atol=1e-8)


def test_channel_view_and_time_domain_ops_constant():
    g = np.exp(-(_X / 5.0) ** 2)
    y = (1 + 2j) * g
    assert np.allclose(P.channel_view(y, "real"), g)
    assert np.allclose(P.channel_view(y, "imag"), 2 * g)
    assert np.allclose(P.channel_view(y, "magnitude"), np.sqrt(5.0) * g)
    with pytest.raises(ValueError):
        P.channel_view(y, "phase")
    assert P.CHANNELS == ("real", "imag", "magnitude")
    assert P.TIME_DOMAIN_OPS <= set(P.OPS)
    assert "ft" in P.TIME_DOMAIN_OPS and "ift" not in P.TIME_DOMAIN_OPS
    assert P.DOMAIN_AGNOSTIC_OPS <= set(P.OPS)
    assert not (P.DOMAIN_AGNOSTIC_OPS & P.TIME_DOMAIN_OPS)
    f = _fid()
    assert P.time_axis_s(f)[-1] == pytest.approx((f.y.size - 1) / f.sw_Hz)


def test_fcor_does_not_mutate_the_caller_array():
    """ft1d wraps a complex fid with np.asarray (no copy); an in-place fcor
    halved the Open-FID dialog's instrument data on every preview."""
    from larmor import fourier

    fid = _fid().y
    first = complex(fid[0])
    fourier.ft1d(fid, 10000.0, 100.0, ops=[{"op": "fcor", "factor": 0.5}])
    fourier.ft1d(fid, 10000.0, 100.0, ops=[{"op": "fcor", "factor": 0.5}])
    assert fid[0] == first
    s = _fid()
    y_in = s.y
    out = P.op_fcor(s, 0.5)
    assert y_in[0] == first and out.y[0] == first * 0.5


def test_chain_start_domain():
    em, ft = {"op": "em", "lb_hz": 50}, {"op": "ft"}
    hil, ift = {"op": "hilbert"}, {"op": "ift"}
    assert P.chain_start_domain([em, ft]) == "time"
    assert P.chain_start_domain([{"op": "phase", "p0": 10}]) == "freq"
    assert P.chain_start_domain([hil, ift, em, ft]) == "freq"
    assert P.chain_start_domain([]) == "freq"
    assert P.chain_start_domain(None) == "freq"
    assert P.chain_start_domain([{"op": "tdeff", "points": 8}, {"op": "zf"},
                                 ft, {"op": "phase"}]) == "time"
    assert P.chain_start_domain([ift, em, ft]) == "freq"
    # the first domain-restricted op decides; agnostic ops are skipped
    assert P.chain_start_domain([{"op": "baseline"}, em]) == "freq"
    assert P.chain_start_domain([{"op": "scale", "factor": 2}, em, ft]) == "time"


def test_subtract_averages():
    x = np.linspace(-50, 50, 1000)
    y = np.exp(-(x / 5) ** 2) + 0.3
    s = P.op_subtract_avg(P.from_processed(x, y + 0j, 100.0))
    assert abs(np.mean(s.y[:100].real)) < 1e-6


def test_scale_sw_stretches_about_centre():
    x = np.linspace(-50, 50, 500)
    s = P.op_scale_sw(P.from_processed(x, x + 0j, 100.0), 2.0)
    assert (s.x_ppm[-1] - s.x_ppm[0]) == pytest.approx(200.0, rel=1e-9)
    assert 0.5 * (s.x_ppm[0] + s.x_ppm[-1]) == pytest.approx(0.0, abs=1e-9)


def test_real_imag_conj():
    x = np.linspace(-10, 10, 100)
    y = np.exp(-(x / 2) ** 2) * (1 + 1j)
    assert np.allclose(P.op_real(P.from_processed(x, y.copy(), 1.0)).y.imag, 0)
    assert np.allclose(P.op_imag(P.from_processed(x, y.copy(), 1.0)).y.real, y.imag)
    assert np.allclose(P.op_conj(P.from_processed(x, y.copy(), 1.0)).y, np.conj(y))


def test_ops_are_registered():
    for name in ("subtract_avg", "scale_sw", "ift", "real", "imag", "conj",
                "twopoint_bg"):
        assert name in P.OPS


def test_twopoint_bg_subtracts_the_line_through_the_two_points():
    x = np.linspace(-20, 20, 400)
    peak = 50.0 * np.exp(-0.5 * (x / 2.0) ** 2)
    tilt = 0.5 * x + 3.0                    # a tilted linear background
    s = P.op_twopoint_bg(P.from_processed(x, peak + tilt, 100.0),
                         x1=-18.0, y1=0.5 * -18.0 + 3.0,
                         x2=18.0, y2=0.5 * 18.0 + 3.0)
    edge = np.concatenate([s.y.real[:20], s.y.real[-20:]])
    assert abs(float(np.mean(edge))) < 0.1                  # background gone
    assert s.y.real.max() == pytest.approx(50.0, abs=0.5)   # peak preserved


def test_twopoint_bg_rejects_coincident_points():
    x = np.linspace(-10, 10, 50)
    s = P.from_processed(x, x + 0j, 100.0)
    with pytest.raises(ValueError):
        P.op_twopoint_bg(s, x1=1.0, y1=0.0, x2=1.0, y2=5.0)


class TestWurstCorrect:
    """A4: the amplitude half of swept-pulse physics (the phase half --
    autophase's p2 -- has been in for a while)."""

    def test_profile_shape(self):
        from larmor.processing import wurst_profile

        x = np.linspace(-500.0, 500.0, 2001)
        w = wurst_profile(x, 216.0, 0.0, 2.0 * 216.0 * 500.0 / 1e3,  # exact span
                          n=80.0, floor=0.05)
        centre = w[np.abs(x) < 100]
        assert np.all(centre > 0.999)              # flat over the middle
        assert w[0] == 0.05 and w[-1] == 0.05      # clamped at the edges
        # monotone roll-off from centre to edge above the floor
        half = w[x >= 0]
        assert np.all(np.diff(half) < 1e-9)

    def test_divide_out_restores_flat_intensity(self):
        from larmor.processing import Spectrum1D, apply, wurst_profile

        x = np.linspace(-400.0, 400.0, 1601)
        sweep_khz = 216.0 * 800.0 / 1e3            # sweep spans the window
        w = wurst_profile(x, 216.0, 0.0, sweep_khz, n=80.0, floor=0.10)
        true = np.ones_like(x)
        s = Spectrum1D(x_ppm=x, y=true * w, sfo1_MHz=216.0,
                       sw_Hz=216.0 * 800.0, domain="freq")
        out = apply(s, [{"op": "wurst_correct", "centre_ppm": 0.0,
                         "sweep_kHz": sweep_khz, "n": 80.0, "floor": 0.10}])
        above = w > 0.10 + 1e-12                   # only where not clamped
        assert np.allclose(out.y[above].real, 1.0, atol=1e-9)

    def test_requires_frequency_domain(self):
        import pytest as _pytest

        from larmor.processing import Spectrum1D, op_wurst_correct

        s = Spectrum1D(x_ppm=None, y=np.zeros(8, complex), sfo1_MHz=216.0,
                       sw_Hz=1000.0, domain="time")
        with _pytest.raises(ValueError, match="frequency-domain"):
            op_wurst_correct(s, centre_ppm=0.0, sweep_kHz=100.0)
