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
    for name in ("subtract_avg", "scale_sw", "ift", "real", "imag", "conj"):
        assert name in P.OPS
