"""Spectral processing: a small, ordered pipeline of operations.

Operations act on a Spectrum1D (complex or real) and are described as JSON
dicts so a processing chain can be stored in a recipe and replayed:

    [{"op": "em", "lb_hz": 100},
     {"op": "zf", "factor": 2},
     {"op": "ft"},
     {"op": "autophase"},
     {"op": "phase", "p0": 12.0, "p1": 0.0},
     {"op": "baseline", "order": 3}]

Time-domain ops (em, zf, ft) apply only when starting from a raw fid;
frequency-domain ops (phase, autophase, baseline) work on any spectrum,
including TopSpin-processed 1r data.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Spectrum1D:
    x_ppm: np.ndarray | None      # None while still in the time domain
    y: np.ndarray                 # complex during processing, real at the end
    sfo1_MHz: float
    sw_Hz: float
    domain: str = "freq"          # "time" | "freq"


# --------------------------------------------------------------------------

def from_bruker_fid(expno_path: str) -> Spectrum1D:
    """Load a raw Bruker fid (read-only) ready for time-domain processing."""
    import nmrglue as ng

    dic, fid = ng.bruker.read(str(expno_path))
    fid = ng.bruker.remove_digital_filter(dic, fid)
    acqus = dic["acqus"]
    return Spectrum1D(x_ppm=None, y=fid.astype(complex),
                      sfo1_MHz=float(acqus["SFO1"]),
                      sw_Hz=float(acqus["SW_h"]), domain="time")


def from_processed(x_ppm: np.ndarray, y: np.ndarray, sfo1_MHz: float,
                   sw_Hz: float = 0.0) -> Spectrum1D:
    return Spectrum1D(x_ppm=np.asarray(x_ppm, float),
                      y=np.asarray(y).astype(complex),
                      sfo1_MHz=sfo1_MHz, sw_Hz=sw_Hz, domain="freq")


# --------------------------------------------------------------------------

def _taxis(s: Spectrum1D) -> np.ndarray:
    return np.arange(s.y.size) / s.sw_Hz


def _need_time(s: Spectrum1D, name: str):
    if s.domain != "time":
        raise ValueError(f"{name} needs time-domain data (load the raw fid)")


def op_em(s: Spectrum1D, lb_hz: float = 0.0) -> Spectrum1D:
    """TopSpin EM: exponential line broadening, w(t) = exp(-pi*LB*t)."""
    _need_time(s, "em")
    s.y = s.y * np.exp(-np.pi * lb_hz * _taxis(s))
    return s


def op_gm(s: Spectrum1D, lb_hz: float = -10.0, gb: float = 0.1) -> Spectrum1D:
    """TopSpin GM (Lorentz-to-Gauss): w(t) = exp(-a*t - b*t^2),
    a = pi*LB (LB is negative), b = -a / (2*GB*AQ)."""
    _need_time(s, "gm")
    t = _taxis(s)
    aq = t[-1] if t[-1] > 0 else 1.0
    a = np.pi * lb_hz
    b = -a / (2.0 * max(gb, 1e-6) * aq)
    s.y = s.y * np.exp(-a * t - b * t * t)
    return s


def op_sine(s: Spectrum1D, ssb: float = 2.0, power: int = 1) -> Spectrum1D:
    """TopSpin SINE/QSINE bell. ssb >= 2 shifts the start phase by pi/ssb
    (2 = cosine bell); ssb 0/1 = pure sine bell. power 2 = QSINE."""
    _need_time(s, "sine")
    t = _taxis(s)
    aq = t[-1] if t[-1] > 0 else 1.0
    phi = np.pi / ssb if ssb >= 2 else 0.0
    w = np.sin(phi + (np.pi - phi) * (t / aq))
    s.y = s.y * (w ** max(1, int(power)))
    return s


def op_traf(s: Spectrum1D, lb_hz: float = 10.0) -> Spectrum1D:
    """Traficante-Ziessow window: resolution enhancement preserving S/N."""
    _need_time(s, "traf")
    t = _taxis(s)
    aq = t[-1] if t[-1] > 0 else 1.0
    tau = 1.0 / (np.pi * max(lb_hz, 1e-6))
    E = np.exp(-t / tau)
    F = np.exp(-(aq - t) / tau)
    s.y = s.y * (E / (E * E + F * F))
    return s


def op_tdeff(s: Spectrum1D, points: int) -> Spectrum1D:
    """TopSpin TDeff: use only the first `points` of the fid."""
    _need_time(s, "tdeff")
    if 0 < points < s.y.size:
        s.y = s.y[:points].copy()
    return s


def op_shift_fid(s: Spectrum1D, points: int = 0) -> Spectrum1D:
    """Left-shift the fid (drop leading points, e.g. before an echo top)."""
    _need_time(s, "shift_fid")
    if points > 0:
        s.y = s.y[points:].copy()
    elif points < 0:
        s.y = np.pad(s.y, (-points, 0))
    return s


def op_fcor(s: Spectrum1D, factor: float = 0.5) -> Spectrum1D:
    """TopSpin FCOR: scale the first fid point (0.5 removes the DC ridge)."""
    _need_time(s, "fcor")
    s.y[0] = s.y[0] * factor
    return s


def op_zf(s: Spectrum1D, factor: int = 2, si: int = 0) -> Spectrum1D:
    """Zero-fill: either by a power-of-two factor, or to an absolute SI."""
    _need_time(s, "zf")
    if si and si > s.y.size:
        n = int(si)
    else:
        n = int(2 ** np.ceil(np.log2(s.y.size * max(1, int(factor)))))
    s.y = np.pad(s.y, (0, max(0, n - s.y.size)))
    return s


def op_ft(s: Spectrum1D, offset_ppm: float = 0.0) -> Spectrum1D:
    if s.domain != "time":
        raise ValueError("data is already in the frequency domain")
    spec = np.fft.fftshift(np.fft.fft(s.y))
    n = spec.size
    freq_hz = np.linspace(s.sw_Hz / 2, -s.sw_Hz / 2, n, endpoint=False)
    s.x_ppm = freq_hz / s.sfo1_MHz + offset_ppm
    s.y = spec
    s.domain = "freq"
    return s


def op_phase(s: Spectrum1D, p0: float = 0.0, p1: float = 0.0,
             pivot_frac: float = 0.5) -> Spectrum1D:
    """Zero- and first-order phase (degrees); p1 pivots at pivot_frac."""
    if s.domain != "freq":
        raise ValueError("phase correction needs frequency-domain data")
    n = s.y.size
    idx = np.arange(n) / max(n - 1, 1)
    ph = np.deg2rad(p0 + p1 * (idx - pivot_frac))
    s.y = s.y * np.exp(1j * ph)
    return s


def op_autophase(s: Spectrum1D, method: str = "scan") -> Spectrum1D:
    """Automatic phasing.

    "scan" (default): fine p0 sweep maximizing positive real signal with a
    negativity penalty, then a Nelder-Mead (p0, p1) refinement -- robust on
    wide solid-state lines. "acme": nmrglue's entropy minimization.
    """
    if s.domain != "freq":
        raise ValueError("autophase needs frequency-domain data")
    if method == "acme":
        import nmrglue as ng

        s.y = ng.process.proc_autophase.autops(s.y, "acme", disp=False)
        return s

    y = s.y
    scale = np.abs(y).max() or 1.0

    def score(p0, p1):
        n = y.size
        ph = np.exp(1j * (p0 + p1 * (np.arange(n) / max(n - 1, 1) - 0.5)))
        r = (y * ph).real / scale
        return r.sum() - 4.0 * np.abs(r[r < 0]).sum()

    phis = np.linspace(-np.pi, np.pi, 1441)
    best = phis[int(np.argmax([score(p, 0.0) for p in phis]))]
    from scipy.optimize import minimize

    res = minimize(lambda v: -score(v[0], v[1]), x0=[best, 0.0],
                   method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-6})
    p0, p1 = res.x
    n = y.size
    s.y = y * np.exp(1j * (p0 + p1 * (np.arange(n) / max(n - 1, 1) - 0.5)))
    return s


def op_baseline(s: Spectrum1D, order: int = 3, k_clip: float = 1.5,
                iterations: int = 10) -> Spectrum1D:
    """Polynomial baseline by iterative asymmetric clipping.

    Fit a polynomial to all points, discard points sticking up more than
    k_clip*sigma above it (i.e. the peaks), refit, repeat until stable --
    the standard automatic baseline used by most NMR software.
    """
    if s.domain != "freq":
        raise ValueError("baseline correction needs frequency-domain data")
    y = s.y.real
    t = np.linspace(-1.0, 1.0, y.size)          # conditioned abscissa
    mask = np.ones(y.size, dtype=bool)
    base = np.zeros_like(y)
    for _ in range(iterations):
        coeffs = np.polynomial.polynomial.polyfit(t[mask], y[mask], order)
        base = np.polynomial.polynomial.polyval(t, coeffs)
        r = y - base
        sigma = r[mask].std() or 1.0
        new_mask = (r < k_clip * sigma) & (r > -4.0 * sigma)
        if new_mask.sum() < (order + 1) * 3 or (new_mask == mask).all():
            break
        mask = new_mask
    s.y = (y - base) + 1j * s.y.imag
    return s


def op_sr(s: Spectrum1D, sr_hz: float = 0.0) -> Spectrum1D:
    """TopSpin SR (spectral reference): shift the ppm axis by SR/SFO1."""
    if s.domain != "freq":
        raise ValueError("sr applies to the frequency domain")
    if s.sfo1_MHz:
        s.x_ppm = s.x_ppm + sr_hz / s.sfo1_MHz
    return s


def op_magnitude(s: Spectrum1D) -> Spectrum1D:
    """Magnitude spectrum (phase-insensitive)."""
    if s.domain != "freq":
        raise ValueError("magnitude applies to the frequency domain")
    s.y = np.abs(s.y) + 0j
    return s


def op_hilbert(s: Spectrum1D) -> Spectrum1D:
    """Rebuild the imaginary part from a real-only spectrum (e.g. TopSpin 1r)
    so that phase correction becomes possible (ssNake's Hilbert)."""
    if s.domain != "freq":
        raise ValueError("hilbert applies to the frequency domain")
    from scipy.signal import hilbert as _hilbert

    real = s.y.real
    # a DC offset (uncorrected first fid point) has no dispersive partner and
    # corrupts the reconstruction: remove it from H, keep it in the real part
    dc = float(np.median(np.concatenate([real[:real.size // 20 or 1],
                                         real[-(real.size // 20 or 1):]])))
    analytic = _hilbert(real - dc)
    s.y = (analytic.conj() + dc)   # real preserved, imag = -H(real - dc)
    return s


OPS = {
    # time domain
    "em": op_em,
    "gm": op_gm,
    "sine": op_sine,
    "traf": op_traf,
    "tdeff": op_tdeff,
    "shift_fid": op_shift_fid,
    "fcor": op_fcor,
    "zf": op_zf,
    "ft": op_ft,
    # frequency domain
    "phase": op_phase,
    "autophase": op_autophase,
    "baseline": op_baseline,
    "sr": op_sr,
    "magnitude": op_magnitude,
    "hilbert": op_hilbert,
}


def apply(s: Spectrum1D, ops: list[dict]) -> Spectrum1D:
    """Apply an ordered list of {"op": name, ...kwargs} steps."""
    for step in ops:
        step = dict(step)
        name = step.pop("op")
        if name not in OPS:
            raise ValueError(f"unknown processing op {name!r} "
                             f"(valid: {sorted(OPS)})")
        s = OPS[name](s, **step)
    return s
