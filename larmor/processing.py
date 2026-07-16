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

def op_em(s: Spectrum1D, lb_hz: float = 0.0) -> Spectrum1D:
    if s.domain != "time":
        raise ValueError("em (exponential apodization) needs time-domain data")
    t = np.arange(s.y.size) / s.sw_Hz
    s.y = s.y * np.exp(-np.pi * lb_hz * t)
    return s


def op_zf(s: Spectrum1D, factor: int = 2) -> Spectrum1D:
    if s.domain != "time":
        raise ValueError("zf (zero filling) needs time-domain data")
    n = int(2 ** np.ceil(np.log2(s.y.size * max(1, int(factor)))))
    s.y = np.pad(s.y, (0, n - s.y.size))
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


def op_autophase(s: Spectrum1D) -> Spectrum1D:
    """Automatic p0/p1 via nmrglue's ACME minimization."""
    if s.domain != "freq":
        raise ValueError("autophase needs frequency-domain data")
    import nmrglue as ng

    s.y = ng.process.proc_autophase.autops(s.y, "acme", disp=False)
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


OPS = {
    "em": op_em,
    "zf": op_zf,
    "ft": op_ft,
    "phase": op_phase,
    "autophase": op_autophase,
    "baseline": op_baseline,
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
