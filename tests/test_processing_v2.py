"""Linear prediction, whole-echo, region tools, algebra and peak picking."""
import numpy as np
import pytest

from larmor import processing


def _fid(n=1024, sw=10000.0, freqs=(1000.0,), t2=0.02):
    t = np.arange(n) / sw
    y = np.zeros(n, complex)
    for f in freqs:
        y += np.exp(2j * np.pi * f * t) * np.exp(-t / t2)
    return processing.Spectrum1D(x_ppm=None, y=y, sfo1_MHz=100.0, sw_Hz=sw,
                                 domain="time")


# ---------------------------------------------------------------- LP
def test_lp_forward_extends_and_keeps_frequency():
    full = _fid(n=1024)
    truncated = _fid(n=256)
    processing.op_lp(truncated, n_predict=768, n_coeff=16, mode="forward")
    assert truncated.y.size == 1024
    # the predicted tail must follow the true decay, not diverge
    assert np.abs(truncated.y[-1]) < np.abs(truncated.y[0])
    err = np.abs(truncated.y - full.y).max() / np.abs(full.y).max()
    assert err < 0.1

    # and the spectrum keeps the peak where it belongs
    processing.apply(truncated, [{"op": "zf"}, {"op": "ft"}])
    peak = truncated.x_ppm[np.argmax(np.abs(truncated.y))]
    assert abs(abs(peak) - 10.0) < 0.2


def test_lp_forward_sharpens_a_truncated_fid():
    """The point of forward LP: truncation broadens; LP restores width."""
    def fwhm_of(spec):
        y = np.abs(spec.y)
        half = y.max() / 2
        idx = np.where(y > half)[0]
        return abs(spec.x_ppm[idx[-1]] - spec.x_ppm[idx[0]])

    trunc = _fid(n=200)
    processing.apply(trunc, [{"op": "zf", "factor": 8}, {"op": "ft"}])
    lp = _fid(n=200)
    processing.apply(lp, [{"op": "lp", "n_predict": 824, "n_coeff": 16},
                          {"op": "zf", "factor": 2}, {"op": "ft"}])
    assert fwhm_of(lp) < fwhm_of(trunc)


def test_lp_backward_repairs_first_points():
    s = _fid(n=512)
    truth = s.y.copy()
    s.y[:4] = 0.0                       # simulate dead-time corruption
    processing.op_lp(s, n_coeff=16, mode="backward", n_replace=4)
    assert s.y.size == 512
    err = np.abs(s.y[:4] - truth[:4]).max() / np.abs(truth).max()
    assert err < 0.15


def test_lp_rejects_short_fid():
    s = _fid(n=20)
    with pytest.raises(ValueError, match="too short"):
        processing.op_lp(s, n_predict=10, n_coeff=16)


# ---------------------------------------------------------------- whole echo
def test_swap_echo_and_symmetric_apodization():
    n, sw = 512, 10000.0
    t = np.arange(n) / sw
    top = 128
    echo = np.exp(-np.abs(np.arange(n) - top) / 40.0).astype(complex)
    s = processing.Spectrum1D(x_ppm=None, y=echo, sfo1_MHz=100.0, sw_Hz=sw,
                              domain="time")
    processing.op_swap_echo(s, top)
    assert s.whole_echo is True
    assert np.abs(s.y[0]) == pytest.approx(np.abs(echo).max())  # top first

    before = np.abs(s.y[s.y.size // 2])
    processing.op_echo_apodize(s, lb_hz=50.0)
    # symmetric window: the middle (farthest from both ends) is damped most
    assert np.abs(s.y[s.y.size // 2]) < before
    assert np.abs(s.y[0]) == pytest.approx(np.abs(echo).max(), rel=1e-6)


# ---------------------------------------------------------------- regions
def test_extract_region():
    s = processing.apply(_fid(), [{"op": "ft"}])
    n0 = s.y.size
    processing.op_extract(s, 20.0, 0.0)
    assert s.y.size < n0
    assert s.x_ppm.min() >= -1e-9 and s.x_ppm.max() <= 20.0 + 1e-9


def test_normalize_window():
    """normalize scales so the WINDOW's maximum is 1 -- that is the point:
    comparing spectra on a feature you choose, not on their global maximum."""
    s = processing.apply(_fid(), [{"op": "ft"}, {"op": "magnitude"}])
    peak_ppm = float(s.x_ppm[np.argmax(np.abs(s.y))])
    hi, lo = peak_ppm + 5.0, peak_ppm - 5.0
    processing.op_normalize(s, hi, lo)
    sel = (s.x_ppm >= lo) & (s.x_ppm <= hi)
    assert np.abs(s.y[sel].real).max() == pytest.approx(1.0, rel=1e-6)


def test_algebra_subtract_background():
    x = np.linspace(-50, 50, 1000)
    from larmor.engine import gauss_lor

    sample = processing.from_processed(x, gauss_lor(x, 0, 5, 10, 1.0) + 2.0, 100.0)
    background = processing.from_processed(x, np.full_like(x, 2.0), 100.0)
    out = processing.combine(sample, background, "subtract")
    assert out.y.real.max() == pytest.approx(10.0, rel=0.02)
    assert abs(out.y.real[0]) < 0.1        # background gone at the edges


def test_algebra_interpolates_a_different_grid():
    from larmor.engine import gauss_lor

    xa = np.linspace(-50, 50, 1000)
    xb = np.linspace(-50, 50, 373)          # deliberately different sampling
    a = processing.from_processed(xa, gauss_lor(xa, 0, 5, 10, 1.0), 100.0)
    b = processing.from_processed(xb, gauss_lor(xb, 0, 5, 10, 1.0), 100.0)
    out = processing.combine(a, b, "subtract")
    assert np.abs(out.y.real).max() < 0.05   # same peak cancels itself


def test_align_finds_a_known_shift():
    from larmor.engine import gauss_lor

    x = np.linspace(-50, 50, 2000)
    a = processing.from_processed(x, gauss_lor(x, 0.0, 4.0, 10.0, 1.0), 100.0)
    b = processing.from_processed(x, gauss_lor(x, 3.0, 4.0, 10.0, 1.0), 100.0)
    shift = processing.align(a, b)
    assert shift == pytest.approx(-3.0, abs=0.15)


# ---------------------------------------------------------------- peaks
def test_pick_peaks():
    from larmor.engine import gauss_lor

    x = np.linspace(-60, 60, 4000)
    y = (gauss_lor(x, 20.0, 4.0, 100.0, 1.0)
         + gauss_lor(x, -10.0, 6.0, 40.0, 1.0)
         + gauss_lor(x, 40.0, 3.0, 2.0, 1.0))      # below threshold
    peaks = processing.pick_peaks(x, y, threshold_frac=0.05, min_sep_ppm=1.0)
    assert len(peaks) == 2
    assert peaks[0]["ppm"] == pytest.approx(20.0, abs=0.1)
    assert peaks[0]["height"] == pytest.approx(100.0, rel=0.02)
    assert peaks[0]["fwhm_ppm"] == pytest.approx(4.0, abs=0.6)
    assert peaks[1]["ppm"] == pytest.approx(-10.0, abs=0.1)


def test_pick_peaks_min_separation():
    from larmor.engine import gauss_lor

    x = np.linspace(-20, 20, 2000)
    y = gauss_lor(x, 0.0, 3.0, 10.0, 1.0) + gauss_lor(x, 1.0, 3.0, 9.0, 1.0)
    wide = processing.pick_peaks(x, y, min_sep_ppm=5.0)
    assert len(wide) == 1
