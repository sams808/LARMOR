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
    whole_echo: bool = False      # set by swap_echo; magnitude is then usual


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
    # frequency axis matching fftshift(fft): ascending, so a +f Hz component
    # lands at +f (verified against real Bruker 1r data: a raw-fid FT peaks at
    # the same ppm as TopSpin's own processed spectrum, up to the SR offset).
    freq_hz = np.fft.fftshift(np.fft.fftfreq(n, d=1.0 / s.sw_Hz))
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


def op_lp(s: Spectrum1D, n_predict: int = 0, n_coeff: int = 16,
          mode: str = "forward", n_replace: int = 0) -> Spectrum1D:
    """Linear prediction (Burg-style autoregression on the analytic fid).

    mode "forward": extend the fid by `n_predict` points -- recovers
        resolution lost to truncation (ssNake lpsvd / TopSpin LPfr).
    mode "backward": rebuild the FIRST `n_replace` points -- repairs
        receiver dead-time distortion (the usual cause of a rolling
        baseline that no polynomial can fix).
    """
    _need_time(s, "lp")
    y = s.y
    if y.size <= n_coeff * 2:
        raise ValueError("fid too short for the requested LP order")

    def ar_coeffs(sig: np.ndarray, order: int) -> np.ndarray:
        # least-squares AR: sig[n] = sum_k a[k] * sig[n-1-k]
        rows = sig.size - order
        A = np.empty((rows, order), dtype=complex)
        for k in range(order):
            A[:, k] = sig[order - 1 - k: order - 1 - k + rows]
        b = sig[order:order + rows]
        a, *_ = np.linalg.lstsq(A, b, rcond=None)
        return a

    if mode == "forward":
        if n_predict <= 0:
            return s
        a = ar_coeffs(y, n_coeff)
        out = list(y)
        for _ in range(n_predict):
            nxt = np.dot(a, np.array(out[-1:-n_coeff - 1:-1]))
            out.append(nxt)
        s.y = np.array(out)
        return s

    if mode == "backward":
        if n_replace <= 0:
            return s
        # predict forward on the time-reversed GOOD part: that extrapolates
        # backwards in real time, into the dead-time-corrupted first points
        good = y[n_replace:]
        if good.size <= n_coeff * 2:
            raise ValueError("not enough good points left for backward LP")
        rev = good[::-1]
        a = ar_coeffs(rev, n_coeff)
        out = list(rev)
        for _ in range(n_replace):
            out.append(np.dot(a, np.array(out[-1:-n_coeff - 1:-1])))
        rebuilt_head = np.array(out[-n_replace:])[::-1]
        s.y = np.concatenate([rebuilt_head, good])
        return s

    raise ValueError(f"unknown lp mode {mode!r} (forward|backward)")


def op_swap_echo(s: Spectrum1D, point: int) -> Spectrum1D:
    """Rotate the fid so the echo top becomes the first point (ssNake
    swapEcho) -- the standard whole-echo preparation."""
    _need_time(s, "swap_echo")
    p = int(point)
    if not (0 < p < s.y.size):
        raise ValueError("echo top must be inside the fid")
    s.y = np.concatenate([s.y[p:], s.y[:p]])
    s.whole_echo = True
    return s


def op_echo_apodize(s: Spectrum1D, lb_hz: float = 0.0) -> Spectrum1D:
    """Symmetric apodization about the echo top for whole-echo data: the
    window decays away from BOTH ends (ssNake wholeEcho)."""
    _need_time(s, "echo_apodize")
    n = s.y.size
    t = np.arange(n) / s.sw_Hz
    t_sym = np.minimum(t, t[::-1])
    s.y = s.y * np.exp(-np.pi * lb_hz * t_sym)
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


def op_extract(s: Spectrum1D, hi_ppm: float, lo_ppm: float) -> Spectrum1D:
    """Keep only a ppm region (ssNake extract)."""
    if s.domain != "freq":
        raise ValueError("extract applies to the frequency domain")
    sel = (s.x_ppm >= min(hi_ppm, lo_ppm)) & (s.x_ppm <= max(hi_ppm, lo_ppm))
    if sel.sum() < 2:
        raise ValueError("extract region contains no data")
    s.x_ppm, s.y = s.x_ppm[sel], s.y[sel]
    return s


def op_scale(s: Spectrum1D, factor: float = 1.0) -> Spectrum1D:
    s.y = s.y * factor
    return s


def op_subtract_avg(s: Spectrum1D, hi_ppm: float | None = None,
                    lo_ppm: float | None = None) -> Spectrum1D:
    """Subtract the mean of a (signal-free) region — a DC/offset correction
    (ssNake 'Subtract Averages'). Defaults to the outer 10% edges."""
    if s.domain != "freq":
        raise ValueError("subtract-averages applies to the frequency domain")
    if hi_ppm is not None and lo_ppm is not None:
        sel = (s.x_ppm >= min(hi_ppm, lo_ppm)) & (s.x_ppm <= max(hi_ppm, lo_ppm))
    else:
        e = max(3, s.y.size // 10)
        sel = np.zeros(s.y.size, bool); sel[:e] = True; sel[-e:] = True
    s.y = s.y - np.mean(s.y[sel].real)
    return s


def op_scale_sw(s: Spectrum1D, factor: float = 1.0) -> Spectrum1D:
    """Scale the spectral width / ppm axis about its centre (ssNake 'Scale SW')
    — a stretch used to correct a mis-set SW or overlay mismatched axes."""
    if s.domain != "freq" or s.x_ppm is None:
        raise ValueError("scale-SW applies to a frequency-domain spectrum")
    c = 0.5 * (float(s.x_ppm[0]) + float(s.x_ppm[-1]))
    s.x_ppm = c + (s.x_ppm - c) * factor
    return s


def op_ift(s: Spectrum1D) -> Spectrum1D:
    """Inverse Fourier transform back to the time domain (ssNake Toggle
    Time/Frequency) so you can re-apodize / reprocess."""
    if s.domain != "freq":
        raise ValueError("inverse FT needs frequency-domain data")
    y = np.asarray(s.y, complex)
    s.y = np.fft.ifft(np.fft.ifftshift(y))          # inverse of fftshift(fft)
    s.x_ppm = None
    s.domain = "time"
    return s


def op_real(s: Spectrum1D) -> Spectrum1D:
    s.y = s.y.real + 0j
    return s


def op_imag(s: Spectrum1D) -> Spectrum1D:
    s.y = s.y.imag + 0j
    return s


def op_conj(s: Spectrum1D) -> Spectrum1D:
    """Complex conjugate — reverses the spectral sense (ssNake)."""
    s.y = np.conj(s.y)
    return s


def op_offset(s: Spectrum1D, value: float = 0.0) -> Spectrum1D:
    s.y = s.y + value
    return s


def op_normalize(s: Spectrum1D, hi_ppm: float | None = None,
                 lo_ppm: float | None = None) -> Spectrum1D:
    """Peak-normalize to 1, optionally inside a window (NMRVEW norm_0_to_1)."""
    if s.domain != "freq":
        raise ValueError("normalize applies to the frequency domain")
    if hi_ppm is not None and lo_ppm is not None:
        sel = (s.x_ppm >= min(hi_ppm, lo_ppm)) & (s.x_ppm <= max(hi_ppm, lo_ppm))
    else:
        sel = np.ones(s.y.shape, bool)
    peak = np.abs(s.y[sel].real).max() or 1.0
    s.y = s.y / peak
    return s


def combine(a: Spectrum1D, b: Spectrum1D, op: str = "subtract",
            scale: float = 1.0) -> Spectrum1D:
    """Spectra algebra on a common axis (dmfit Dual / background removal).

    b is interpolated onto a's ppm axis, so the two need not share a grid.
    """
    if a.domain != "freq" or b.domain != "freq":
        raise ValueError("algebra needs two frequency-domain spectra")
    bi = np.interp(a.x_ppm, b.x_ppm, b.y.real, left=0.0, right=0.0) + \
        1j * np.interp(a.x_ppm, b.x_ppm, b.y.imag, left=0.0, right=0.0)
    bi = bi * scale
    if op == "subtract":
        y = a.y - bi
    elif op == "add":
        y = a.y + bi
    elif op == "multiply":
        y = a.y * bi
    elif op == "divide":
        y = np.divide(a.y, bi, out=np.zeros_like(a.y), where=np.abs(bi) > 1e-12)
    else:
        raise ValueError(f"unknown algebra op {op!r}")
    return Spectrum1D(x_ppm=a.x_ppm.copy(), y=y, sfo1_MHz=a.sfo1_MHz,
                      sw_Hz=a.sw_Hz, domain="freq")


def align(a: Spectrum1D, b: Spectrum1D, hi_ppm: float | None = None,
          lo_ppm: float | None = None) -> float:
    """ppm shift TO APPLY TO b so it lands on a (ssNake align).

    Cross-correlation of the real parts. Apply the result with
    op_sr(b, shift * b.sfo1_MHz), or add it to b.x_ppm.
    """
    x = a.x_ppm
    ya = a.y.real
    yb = np.interp(x, b.x_ppm, b.y.real, left=0.0, right=0.0)
    if hi_ppm is not None and lo_ppm is not None:
        sel = (x >= min(hi_ppm, lo_ppm)) & (x <= max(hi_ppm, lo_ppm))
        x, ya, yb = x[sel], ya[sel], yb[sel]
    ya = ya - ya.mean()
    yb = yb - yb.mean()
    corr = np.correlate(ya, yb, mode="full")
    # lag maximizing sum(ya[n] * yb[n - lag]) is (peak_a - peak_b) in points,
    # which is exactly the displacement b must undergo to reach a
    lag = int(np.argmax(corr)) - (len(yb) - 1)
    dppm = float(np.mean(np.diff(x)))
    return lag * dppm


def pick_peaks(x_ppm: np.ndarray, y: np.ndarray, threshold_frac: float = 0.05,
               min_sep_ppm: float = 0.0) -> list[dict]:
    """Peak picking with parabolic sub-point interpolation.

    Returns [{"ppm":…, "height":…, "fwhm_ppm":…}] sorted by descending height
    -- directly usable to seed one line per peak.
    """
    y = np.asarray(y, float)
    order = np.argsort(x_ppm)
    x, yy = np.asarray(x_ppm)[order], y[order]
    thr = threshold_frac * float(np.abs(yy).max() or 1.0)
    peaks = []
    for i in range(1, len(yy) - 1):
        if yy[i] < thr or not (yy[i] >= yy[i - 1] and yy[i] >= yy[i + 1]):
            continue
        # parabolic vertex through the three points
        d = yy[i - 1] - 2 * yy[i] + yy[i + 1]
        delta = 0.5 * (yy[i - 1] - yy[i + 1]) / d if d else 0.0
        dx = float(np.mean(np.diff(x)))
        ppm = float(x[i] + delta * dx)
        height = float(yy[i] - 0.25 * (yy[i - 1] - yy[i + 1]) * delta)
        # local FWHM by walking down to half height
        half = height / 2.0
        li = i
        while li > 0 and yy[li] > half:
            li -= 1
        ri = i
        while ri < len(yy) - 1 and yy[ri] > half:
            ri += 1
        peaks.append({"ppm": ppm, "height": height,
                      "fwhm_ppm": float(abs(x[ri] - x[li])) or abs(dx)})
    peaks.sort(key=lambda p: -p["height"])
    if min_sep_ppm > 0:
        kept: list[dict] = []
        for p in peaks:
            if all(abs(p["ppm"] - q["ppm"]) >= min_sep_ppm for q in kept):
                kept.append(p)
        peaks = kept
    return peaks


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
    "lp": op_lp,
    "swap_echo": op_swap_echo,
    "echo_apodize": op_echo_apodize,
    "ft": op_ft,
    # frequency domain
    "phase": op_phase,
    "autophase": op_autophase,
    "baseline": op_baseline,
    "sr": op_sr,
    "magnitude": op_magnitude,
    "hilbert": op_hilbert,
    "extract": op_extract,
    "scale": op_scale,
    "offset": op_offset,
    "normalize": op_normalize,
    "subtract_avg": op_subtract_avg,
    "scale_sw": op_scale_sw,
    "ift": op_ift,
    "real": op_real,
    "imag": op_imag,
    "conj": op_conj,
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
