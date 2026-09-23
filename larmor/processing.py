"""Spectral processing: a small, ordered pipeline of operations.

Operations act on a Spectrum1D (complex or real) and are described as JSON
dicts so a processing chain can be stored in a recipe and replayed:

    [{"op": "em", "lb_hz": 100},
     {"op": "zf", "factor": 2},
     {"op": "ft"},
     {"op": "autophase"},
     {"op": "phase", "p0": 12.0, "p1": 0.0},
     {"op": "baseline", "order": 3}]

Time-domain ops (em, zf, ft) apply from a raw fid OR after an `ift` step,
so a chain may start in the frequency domain (`[hilbert, ift, em, ft]`
re-apodizes a TopSpin-processed 1r); frequency-domain ops (phase, autophase,
baseline) work on any spectrum. `chain_start_domain` says which kind of input
a recorded chain needs; TIME_DOMAIN_OPS lists the ops that refuse
frequency-domain data.
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
    #: the frequency axis parked by op_ift so op_ft can restore it exactly
    #: (a processed spectrum has no SW of its own; see op_ift / op_ft)
    x_ppm_hold: np.ndarray | None = None


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


def time_axis_s(s: Spectrum1D) -> np.ndarray:
    """Acquisition time of every point of a time-domain Spectrum1D, in s
    (t = k / SW). Public name of the abscissa every window function uses."""
    return _taxis(s)


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
    # never write into the caller's array: fourier.ft1d wraps a fid with
    # np.asarray (no copy for complex input), so an in-place s.y[0] *= f halved
    # the Open-FID dialog's instrument data again on EVERY preview -- the
    # spectrum drifted by a DC step per control change and no longer matched
    # the chain that claimed to have produced it
    y = np.array(s.y, dtype=complex, copy=True)
    y[0] = y[0] * factor
    s.y = y
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
    hold = s.x_ppm_hold
    if hold is not None and hold.size == n:
        # back from an ift of a processed spectrum with the same point count:
        # the original axis is restored exactly (a ppm grid rebuilt from a
        # derived SW would be off by float noise and would lose any SR/offset
        # already applied to the axis)
        s.x_ppm = hold + offset_ppm
    elif hold is not None and hold.size:
        # zero-filled (or truncated) since the ift: fftshift puts the 0 Hz bin
        # at index n//2 on BOTH grids, so anchoring the new grid on the held
        # zero bin keeps every peak at its ppm. The mid-span value
        # 0.5*(x[0]+x[-1]) would be half a bin off for even n.
        s.x_ppm = freq_hz / s.sfo1_MHz + hold[hold.size // 2] + offset_ppm
    else:
        s.x_ppm = freq_hz / s.sfo1_MHz + offset_ppm
    s.x_ppm_hold = None
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


def op_iterbaseline(s: Spectrum1D, dead_time_pts: int = 0,
                    smoothness: float = 1.0, threshold_factor: float = 1.0,
                    max_iter: int = 60) -> Spectrum1D:
    """Iterative baseline correction for dead-time-truncated spectra.

    Adaptation of Yon, Fayon, Massiot & Sarou-Kanian, Solid State Nucl. Magn.
    Reson. 110, 101699 (2020) (github.com/maximeYon/Baseline_Corrector): an
    iterative histogram-thresholded smoothing baseline, optionally restricted to
    broad (dead-time) components in the time domain. Use it when a polynomial
    baseline fails on a rolling baseline from receiver dead time. Set
    `dead_time_pts` ~ 2*DE/DW to enable the time-domain restriction.
    """
    if s.domain != "freq":
        raise ValueError("baseline correction needs frequency-domain data")
    from larmor.baseline import iterative_baseline

    res = iterative_baseline(
        s.y.real, x=s.x_ppm, dead_time_pts=int(dead_time_pts),
        smoothness=float(smoothness), threshold_factor=float(threshold_factor),
        max_iter=int(max_iter))
    s.y = res.corrected + 1j * s.y.imag
    return s


def op_flat_baseline(s: Spectrum1D, edge_frac: float = 0.05) -> Spectrum1D:
    """Flat baseline: subtract the median of the two spectrum edges (a quick
    DC-offset correction when the baseline is genuinely flat, not rolling).
    `edge_frac` of the points at each end are pooled for the median."""
    if s.domain != "freq":
        raise ValueError("baseline correction needs frequency-domain data")
    y = s.y.real
    n = max(3, int(y.size * edge_frac))
    level = float(np.median(np.concatenate([y[:n], y[-n:]])))
    s.y = (y - level) + 1j * s.y.imag
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


def wurst_profile(ppm: "np.ndarray", sfo1_MHz: float, centre_ppm: float,
                  sweep_kHz: float, n: float = 80.0,
                  floor: float = 0.10) -> "np.ndarray":
    """Amplitude weighting a WURST-N sweep imprints across its band.

    A linear chirp visits each frequency at one moment of the pulse, so the
    frequency-domain weighting IS the pulse's amplitude envelope read at
    that moment: WURST-N has A(t) = 1 - |cos(pi t / tau_p)|^N, t in
    [0, tau_p], which maps to W(nu) = 1 - |cos(pi f)|^N with f the
    fractional position of nu across the sweep. Outside the sweep, and
    wherever W falls below ``floor`` (relative), the returned profile is
    clamped to ``floor`` -- dividing by less would amplify noise, not signal.
    For large N (WURST-80 is typical) the profile is near-flat over ~90 % of
    the band and only the outer edges are corrected.
    """
    x = np.asarray(ppm, float)
    if sfo1_MHz <= 0 or sweep_kHz <= 0:
        return np.ones_like(x)
    half_ppm = (sweep_kHz * 1e3 / 2.0) / sfo1_MHz
    f = (x - (centre_ppm - half_ppm)) / (2.0 * half_ppm)   # 0..1 across sweep
    w = 1.0 - np.abs(np.cos(np.pi * np.clip(f, 0.0, 1.0))) ** float(n)
    w[(f < 0.0) | (f > 1.0)] = 0.0
    floor = float(np.clip(floor, 1e-3, 1.0))
    return np.maximum(w, floor)


def op_wurst_correct(s: Spectrum1D, centre_ppm: float, sweep_kHz: float,
                     n: float = 80.0, floor: float = 0.10) -> Spectrum1D:
    """Divide out the WURST/chirp excitation profile (see wurst_profile), so
    intensities across a swept-excitation wideline pattern (WCPMG and
    friends) are comparable and quantification means something. LARMOR
    already fits the QUADRATIC PHASE such a sweep imprints (autophase p2);
    this is the amplitude half of the same physics. Frequency domain only."""
    if s.domain != "freq" or s.x_ppm is None:
        raise ValueError("wurst_correct applies to a frequency-domain "
                         "spectrum")
    w = wurst_profile(s.x_ppm, s.sfo1_MHz, centre_ppm, sweep_kHz,
                      n=n, floor=floor)
    s.y = s.y / w
    return s


def op_scale(s: Spectrum1D, factor: float = 1.0) -> Spectrum1D:
    s.y = s.y * factor
    return s


def op_twopoint_bg(s: Spectrum1D, x1: float, y1: float, x2: float, y2: float
                   ) -> Spectrum1D:
    """Subtract the straight line through two user-picked points — a flat or
    tilted 2-point linear background, extrapolated across the whole spectrum
    (dmfit-style; the manual pick tool in the app and the batch dialog)."""
    if s.domain != "freq" or s.x_ppm is None:
        raise ValueError("2-point background needs frequency-domain data")
    if abs(x2 - x1) < 1e-12:
        raise ValueError("the two background points must be at different positions")
    m = (y2 - y1) / (x2 - x1)
    base = m * (s.x_ppm - x1) + y1
    s.y = s.y - base
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
    Time/Frequency) so you can re-apodize / reprocess.

    Works on a processed spectrum too (TopSpin 1r, CSV): a descending axis is
    reversed first (ifftshift assumes the ascending grid op_ft produces), the
    spectral width is derived from the axis when the spectrum carries none
    (`sw_Hz = dx * n * sfo1`, which needs the Larmor frequency), and the axis
    is parked in `x_ppm_hold` so the next op_ft restores it exactly instead
    of rebuilding a grid about 0 ppm. A real-only spectrum should go through
    `hilbert` first: its ift is two-sided (hermitian), and a one-sided window
    on it gives a dispersive line.
    """
    if s.domain != "freq":
        raise ValueError("inverse FT needs frequency-domain data")
    y = np.asarray(s.y, complex)
    x = None if s.x_ppm is None else np.asarray(s.x_ppm, float)
    if x is not None and x.size > 1 and x[-1] < x[0]:
        x, y = x[::-1].copy(), y[::-1]
    if not s.sw_Hz > 0:
        if x is None or x.size < 2:
            raise ValueError("inverse FT needs a spectral width or a ppm axis")
        if not s.sfo1_MHz > 0:
            raise ValueError(
                "inverse FT of a processed spectrum needs the Larmor frequency "
                "— set it in Process ▸ Experiment parameters")
        n = x.size
        dx = abs(float(x[-1]) - float(x[0])) / max(n - 1, 1)
        s.sw_Hz = dx * n * s.sfo1_MHz
    s.x_ppm_hold = None if x is None else x.copy()
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
    "wurst_correct": op_wurst_correct,
    "autophase": op_autophase,
    "baseline": op_baseline,
    "iterbaseline": op_iterbaseline,
    "flat_baseline": op_flat_baseline,
    "twopoint_bg": op_twopoint_bg,
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

#: the ops that refuse frequency-domain input: every op guarded by _need_time,
#: plus ft itself. The single truth for "this chain needs a time-domain start"
#: (loader.apply_processing, the desktop's processing panel). Must stay a
#: subset of OPS -- pinned by tests/test_processing_ops.py.
TIME_DOMAIN_OPS = frozenset({
    "em", "gm", "sine", "traf", "tdeff", "shift_fid", "fcor", "zf", "lp",
    "swap_echo", "echo_apodize", "ft",
})

#: ops that work in either domain and therefore say nothing about where a
#: chain starts
DOMAIN_AGNOSTIC_OPS = frozenset({"scale", "offset", "real", "imag", "conj"})

#: display channels of a complex spectrum or fid (see channel_view)
CHANNELS = ("real", "imag", "magnitude")


def chain_start_domain(ops: list[dict] | None) -> str:
    """Which input a recorded chain needs: "time" (a raw fid) or "freq".

    Walks the chain in order; the first op that restricts the domain decides
    -- a time-domain op (TIME_DOMAIN_OPS) means the raw fid, any other
    registered, non-agnostic op (ift, hilbert, phase, baseline, ...) means a
    processed spectrum, so `[hilbert, ift, em, ft]` replays from the processed
    arrays while `[em, ft]` needs the fid. Empty or all-agnostic: "freq".
    """
    for step in ops or []:
        name = step.get("op")
        if name in TIME_DOMAIN_OPS:
            return "time"
        if name in OPS and name not in DOMAIN_AGNOSTIC_OPS:
            return "freq"
    return "freq"


def channel_view(y: np.ndarray, channel: str) -> np.ndarray:
    """One real-valued display channel of complex data: "real", "imag" or
    "magnitude" (|y|). A projection for display -- unlike op_real / op_imag /
    op_magnitude it changes nothing in the pipeline."""
    y = np.asarray(y)
    if channel == "real":
        return np.real(y).astype(float, copy=False)
    if channel == "imag":
        return np.imag(y).astype(float, copy=False)
    if channel == "magnitude":
        return np.abs(y).astype(float, copy=False)
    raise ValueError(f"unknown display channel {channel!r} (valid: {CHANNELS})")


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
