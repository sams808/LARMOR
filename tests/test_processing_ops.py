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
