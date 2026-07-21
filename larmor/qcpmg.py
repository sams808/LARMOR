"""QCPMG processing (Larsen, Jakobsen & Nielsen, J. Phys. Chem. A 1997/1998).

A QCPMG acquisition is a train of echoes recorded after a single excitation.
Two standard ways to turn it into a spectrum:

  * **spikelet** -- Fourier-transform the whole echo train. The result is a
    manifold of sharp "spikelets" separated by 1/τ_echo (the echo period), whose
    intensities trace the underlying quadrupolar/CSA powder pattern. This is the
    familiar QCPMG display.

  * **coadded envelope** -- add the individual echoes together (aligning their
    tops) into a single echo and transform that. The spikelets vanish and you
    recover the continuous powder lineshape at much higher S/N -- the form you
    fit.

The echo period is set in the pulse program (a fixed spikelet separation, e.g.
cnst = 2000 Hz here); we recover it robustly from the magnitude FID's
autocorrelation so no pulse-program bookkeeping is needed.
"""
from __future__ import annotations

import numpy as np


def detect_period(fid: np.ndarray, min_period: int = 8) -> int:
    """Echo period in points, from the first autocorrelation maximum of the
    magnitude FID (the echo train is quasi-periodic at the echo spacing)."""
    mag = np.abs(np.asarray(fid))
    m = mag - mag.mean()
    ac = np.correlate(m, m, mode="full")[len(m) - 1:]
    ac[:max(1, min_period)] = 0.0
    half = max(2, len(ac) // 2)
    return int(np.argmax(ac[:half]))


def _axis_ppm(nfft: int, sw_Hz: float, sfo_MHz: float,
              carrier_ppm: float) -> np.ndarray:
    freq = np.fft.fftshift(np.fft.fftfreq(nfft, d=1.0 / sw_Hz))
    return carrier_ppm + freq / (sfo_MHz or 1.0)


def spikelet_spectrum(fid: np.ndarray, sw_Hz: float, sfo_MHz: float,
                      carrier_ppm: float = 0.0, lb_Hz: float = 50.0,
                      zf: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """FT the whole echo train -> the spikelet spectrum (complex)."""
    fid = np.asarray(fid, complex)
    n = fid.size
    t = np.arange(n) / sw_Hz
    w = np.exp(-np.pi * lb_Hz * t)
    nfft = int(2 ** np.ceil(np.log2(n * max(1, zf))))
    spec = np.fft.fftshift(np.fft.fft(fid * w, n=nfft))
    ppm = _axis_ppm(nfft, sw_Hz, sfo_MHz, carrier_ppm)
    order = np.argsort(ppm)
    return ppm[order], spec[order]


def coadd_echoes(fid: np.ndarray, period: int,
                 drop_first: int = 1) -> np.ndarray:
    """Sum the echoes into one, aligning each echo's top (magnitude peak).

    The first partial echo after excitation is dropped by default.
    """
    fid = np.asarray(fid, complex)
    if period < 4:
        return fid
    n_full = fid.size // period
    if n_full < 2:
        return fid
    block = fid[:n_full * period].reshape(n_full, period)
    block = block[max(0, drop_first):]
    # reference top from a mid-train echo (clean, past any startup transient)
    ref_top = int(np.argmax(np.abs(block[block.shape[0] // 2])))
    out = np.zeros(period, complex)
    for echo in block:
        shift = ref_top - int(np.argmax(np.abs(echo)))
        out += np.roll(echo, shift)
    return out


def coadd_spectrum(fid: np.ndarray, period: int, sw_Hz: float, sfo_MHz: float,
                   carrier_ppm: float = 0.0, lb_Hz: float = 100.0, zf: int = 16,
                   drop_first: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Coadd the echoes, then FT the single coadded echo -> the continuous
    powder envelope (magnitude, spikelet-free)."""
    echo = coadd_echoes(fid, period, drop_first=drop_first)
    top = int(np.argmax(np.abs(echo)))
    echo = np.roll(echo, -top)                 # echo top -> t = 0
    m = echo.size
    t = np.arange(m) / sw_Hz
    w = np.exp(-np.pi * lb_Hz * t)
    nfft = int(2 ** np.ceil(np.log2(max(m, 2) * max(1, zf))))
    spec = np.fft.fftshift(np.fft.fft(echo * w, n=nfft))
    ppm = _axis_ppm(nfft, sw_Hz, sfo_MHz, carrier_ppm)
    order = np.argsort(ppm)
    return ppm[order], np.abs(spec[order])


def spikelet_spacing_ppm(period: int, sw_Hz: float, sfo_MHz: float) -> float:
    return (sw_Hz / period) / (sfo_MHz or 1.0) if period else 0.0
