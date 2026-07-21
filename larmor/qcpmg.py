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


# --------------------------------------------------------------------------
# Full ssNake-style "sum echo" workflow: split -> (T2 fit / weight) -> sum ->
# whole-echo processing (swap the echo top to t=0 -> a clean absorption
# lineshape you FIT, instead of the spikelet comb).

def split_echoes(fid: np.ndarray, period: int) -> np.ndarray:
    """Reshape the train into (n_echoes, period), dropping trailing all-zero
    echoes (appended zeros / fully relaxed slots)."""
    fid = np.asarray(fid, complex)
    n = fid.size // period
    if n < 1:
        return fid[None, :]
    block = fid[:n * period].reshape(n, period)
    e = np.abs(block).sum(axis=1)
    nz = e > 1e-9 * (e.max() or 1.0)
    if nz.any():
        block = block[:int(np.max(np.where(nz))) + 1]
    return block


def echo_top_point(echoes: np.ndarray) -> int:
    """Point index of the echo maximum (from a clean mid-train echo)."""
    ref = echoes[echoes.shape[0] // 2] if echoes.shape[0] > 2 else echoes[0]
    return int(np.argmax(np.abs(ref)))


def echo_decay(echoes: np.ndarray, top: int) -> np.ndarray:
    """Echo-top intensity vs echo number -- the transverse (T2') decay."""
    return np.abs(echoes[:, top])


def fit_t2(tau_s: float, decay: np.ndarray):
    """Mono-exponential T2 from the echo-top decay. tau_s = echo spacing (s).
    Returns (T2_seconds, model_callable(times_s))."""
    from scipy.optimize import curve_fit

    d = np.asarray(decay, float)
    d = d / (d.max() or 1.0)
    t = np.arange(d.size) * tau_s
    try:
        popt, _ = curve_fit(
            lambda tt, a, T2: a * np.exp(-tt / T2), t, d,
            p0=[1.0, tau_s * max(5, d.size / 3)],
            bounds=([0.0, tau_s], [np.inf, tau_s * d.size * 20]), maxfev=10000)
        a, T2 = float(popt[0]), float(popt[1])
    except Exception:
        a, T2 = 1.0, tau_s * max(1, d.size / 2)
    return T2, (lambda tt, a=a, T2=T2: a * np.exp(-np.asarray(tt, float) / T2))


def sum_echoes(echoes: np.ndarray, tau_s: float,
               t2_weight_s: float | None = None) -> np.ndarray:
    """Coherently add the echoes. With ``t2_weight_s`` set, weight echo k by
    exp(-k·tau/T2) -- the matched filter that maximises S/N (ssNake's T2
    weighting via a Lorentzian LB = 1/(πT2) along the echo dimension)."""
    n = echoes.shape[0]
    w = np.ones(n)
    if t2_weight_s:
        w = np.exp(-np.arange(n) * tau_s / t2_weight_s)
    return (echoes * w[:, None]).sum(axis=0)


def _gaussian_apod(n: int, sw_Hz: float, gb_Hz: float) -> np.ndarray:
    if not gb_Hz:
        return np.ones(n)
    t = np.arange(n) / sw_Hz
    return np.exp(-((np.pi * gb_Hz * t) ** 2) / (4.0 * np.log(2.0)))


def sum_echo_spectrum(fid: np.ndarray, period: int, sw_Hz: float, sfo_MHz: float,
                      carrier_ppm: float = 0.0, top: int | None = None,
                      t2_weight_s: float | None = None, p0_deg: float = 0.0,
                      p1_deg: float = 0.0, gb_Hz: float = 0.0, zf: int = 16,
                      ) -> tuple[np.ndarray, np.ndarray]:
    """The fittable QCPMG spectrum: sum the echoes, whole-echo process (swap the
    top to t=0), Gaussian-apodize, FT, and apply a p0/p1 phase. Returns
    (ppm, complex spectrum); take .real for the absorption lineshape to fit."""
    echoes = split_echoes(fid, period)
    if top is None:
        top = echo_top_point(echoes)
    tau = period / sw_Hz
    summed = sum_echoes(echoes, tau, t2_weight_s)
    echo = np.roll(summed, -int(top))                 # whole echo: top -> t=0
    m = echo.size
    echo = echo * _gaussian_apod(m, sw_Hz, gb_Hz)
    nfft = int(2 ** np.ceil(np.log2(max(m, 2) * max(1, zf))))
    spec = np.fft.fftshift(np.fft.fft(echo, n=nfft))
    freq = np.fft.fftshift(np.fft.fftfreq(nfft, d=1.0 / sw_Hz))
    ppm = carrier_ppm + freq / (sfo_MHz or 1.0)
    ramp = np.arange(nfft) / nfft
    spec = spec * np.exp(-1j * (np.deg2rad(p0_deg) + np.deg2rad(p1_deg) * ramp))
    order = np.argsort(ppm)
    return ppm[order], spec[order]


def autophase0(spec: np.ndarray) -> float:
    """Zero-order phase (deg) that maximises the real integral of a spectrum."""
    ph = np.linspace(-np.pi, np.pi, 361)
    scores = [np.real(spec * np.exp(1j * p)).sum() for p in ph]
    return float(np.degrees(ph[int(np.argmax(scores))]))
