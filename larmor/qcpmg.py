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

The echo period is set in the pulse program (a fixed spikelet separation), so
it is READ from ``acqus`` (CNST7/CNST8) whenever the pulse program recorded it
-- ``detect_period`` is only a fallback for datasets that did not.

Two time axes, and why it matters
---------------------------------
ssNake's ``Matrix -> Split`` copies the ACQUISITION sweep width onto the new
echo dimension, so a T2 fitted there is expressed in "one dwell per echo"
units, not seconds: ``T2_physical = T2_ssNake * points_per_echo``. Both are
self-consistent, but only the physical value transfers to another program.
:class:`T2Fit` therefore carries both, and the matched filter for each
(``lb_Hz`` and ``lb_ssnake_Hz``), so a value copied out of ssNake can be
reconciled instead of silently misread by a factor of the echo length.

Reference: `QCPMG processing in ssNake` (2019) and Larsen et al.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def detect_period(fid: np.ndarray, min_period: int = 8) -> int:
    """Echo period in points from the magnitude FID's autocorrelation.

    FALLBACK only -- prefer :func:`echo_period_from_meta`, which reads the
    period from the pulse program exactly. Returns 0 when no periodicity is
    found, so a caller can say "not detected" instead of using a plausible
    wrong integer.

    The central autocorrelation lobe is blanked by its OWN measured width, not
    by a fixed ``min_period``: an echo wider than that constant leaves the lobe
    partially intact and it wins the argmax, which is why a real 35Cl train
    (true period 293) used to come back as 8.
    """
    mag = np.abs(np.asarray(fid, complex))
    m = mag - mag.mean()
    if m.size < 16 or not np.any(m):
        return 0
    ac = np.correlate(m, m, mode="full")[m.size - 1:]
    if ac[0] <= 0:
        return 0
    ac = ac / ac[0]
    # width of the central lobe: first rise after the initial monotonic fall
    fall = np.where(np.diff(ac) > 0)[0]
    lobe = int(fall[0]) + 1 if fall.size else max(1, min_period)
    start = max(int(min_period), lobe)
    half = ac.size // 2
    if start + 2 > half:                # nothing left to search
        return 0
    seg = ac[start:half]
    k = int(np.argmax(seg))
    # require a real peak, not just the biggest noise sample
    if seg[k] < 0.15:
        return 0
    return start + k


def period_correlation(fid: np.ndarray, period: int) -> float:
    """How well the train repeats at ``period``: normalised correlation of
    two adjacent period-blocks anchored on the strongest echo. 1.0 = the next
    block is a scaled copy (the right period); near 0 = the split cuts
    through echoes.

    Anchoring on the echo matters: blocks taken blindly from the start of
    the train are dominated by the smooth pre-echo ramp, which correlates at
    ANY small lag. A VERIFIER, not a detector: a smooth echo is coherent over
    a few samples, so tiny periods always score high -- only compare
    candidates at least as long as the echo itself (which is what
    :func:`find_period_by_correlation` arranges). On a real 81Br train this
    scored 0.98 at the true 475 points and 0.01 at a plausible-looking wrong
    candidate."""
    x = np.asarray(fid, complex)
    p = int(period)
    n = x.size
    if p < 4 or 2 * p > n:
        return 0.0
    k = int(np.argmax(np.abs(x)))
    i0 = min(max(k - p // 3, 0), n - 2 * p)
    a, b = x[i0:i0 + p], x[i0 + p:i0 + 2 * p]
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na == 0.0 or nb == 0.0:
        return 0.0
    return float(abs(np.vdot(a, b)) / (na * nb))


def find_period_by_correlation(fid: np.ndarray, min_period: int = 8
                               ) -> tuple[int, float]:
    """The period MEASURED from the data: autocorrelation proposes
    candidates (with harmonics, divisors and ±1 jitter), the echo-anchored
    :func:`period_correlation` verifies them, and the SMALLEST candidate
    within 90 % of the best score wins — sub-multiples of the true period cut
    echoes apart and collapse the score, while multiples always score well,
    so the smallest good one is the fundamental.

    Returns (0, best_score) when nothing repeats convincingly (score < 0.5),
    so a caller can say "not found" instead of guessing."""
    x = np.asarray(fid, complex)
    n = x.size
    d = detect_period(x, min_period)
    if not d:
        return 0, 0.0
    cands: set[int] = set()
    for f in (1, 2, 3, 4):
        cands.add(d * f)
        if d % f == 0:
            cands.add(d // f)
    for c in list(cands):
        cands.update({c - 1, c + 1})
    scored = [(p, period_correlation(x, p)) for p in sorted(cands)
              if min_period <= p <= n // 2]
    if not scored:
        return 0, 0.0
    best = max(s for _, s in scored)
    if best < 0.5:
        return 0, best
    pick = next(p for p, s in scored if s >= 0.9 * best)

    # refine ±2: adjacent-block correlation barely separates p from p±1, but
    # a one-point error drifts k points over k periods, so add a three-
    # periods-apart block to the score (a real 81Br train: 475 vs 474)
    def long_score(p: int) -> float:
        s1 = period_correlation(x, p)
        if n < 4 * p:
            return s1
        k = int(np.argmax(np.abs(x)))
        i0 = min(max(k - p // 3, 0), n - 4 * p)
        a, b = x[i0:i0 + p], x[i0 + 3 * p:i0 + 4 * p]
        na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
        if na == 0.0 or nb == 0.0:
            return s1
        return s1 + float(abs(np.vdot(a, b)) / (na * nb))

    lo = max(int(min_period), pick - 2)
    hi2 = min(n // 2, pick + 2)
    pick = max(range(lo, hi2 + 1), key=long_score)
    return pick, period_correlation(x, pick)


#: widening of a MAGNITUDE line relative to its absorption part, for a CAUSAL
#: (one-sided) FID, where absorption and dispersion are Hilbert partners:
#: |L| = 1/sqrt(1+x^2) halves at x = sqrt(3), the absorption 1/(1+x^2) at
#: x = 1. It does NOT apply to whole-echo processing -- see
#: :func:`magnitude_spectrum`.
MAGNITUDE_LORENTZ_WIDENING = 3.0 ** 0.5


def noise_floor(y: np.ndarray, edge_frac: float = 0.10) -> float:
    """Baseline level estimated from the two EDGES of a spectrum, which are
    signal-free in any sanely-referenced acquisition.

    The median of the WHOLE trace is not a floor once the pattern occupies a
    large fraction of the window -- exactly the wide-line case that matters
    here. On a real 81Br pattern filling ~59 % of the axis the global median
    came out 6.2x the true floor, and subtracting it moved delta_CG by 44 ppm.
    """
    y = np.asarray(y, float)
    k = max(1, int(round(float(np.clip(edge_frac, 0.01, 0.49)) * y.size)))
    return float(np.median(np.concatenate([y[:k], y[-k:]])))


def magnitude_spectrum(spec: np.ndarray, *, subtract_floor: bool = True,
                       edge_frac: float = 0.10) -> np.ndarray:
    """|spectrum| -- TopSpin's ``mc`` -- with its rectified noise floor
    removed (see :func:`noise_floor`).

    Magnitude is phase-independent by construction, which makes it the
    fallback whenever a phase error cannot be written as p0/p1 (a swept
    WURST/chirp refocusing pulse imprints a QUADRATIC phase, for instance).

    The usual objection -- that magnitude costs a factor sqrt(3) in width --
    does NOT apply to a whole-echo spectrum. That factor comes from the
    dispersion of a CAUSAL FID adding in quadrature; a whole echo is
    symmetric about t = 0, so its ideal spectrum is REAL and |spec| simply
    recovers the absorption lineshape (measured widening 1.000 through
    :func:`whole_echo_ft`, against 1.70 for the same linewidth from a
    one-sided FID). What magnitude really costs is the sign: rectification
    turns zero-mean noise into a positive pedestal, and folds any genuinely
    negative feature upward.
    """
    m = np.abs(np.asarray(spec))
    if subtract_floor and m.size:
        m = m - noise_floor(m, edge_frac)
    return m


def carrier_ppm(meta: dict) -> tuple[float, bool]:
    """Transmitter offset in ppm, correctly referenced, and whether it is.

    Returns ``(ppm, referenced)``. The right answer is ``(SFO1-SF)*1e6/SF``
    using SF from ``pdata/1/procs`` -- i.e. the same reference TopSpin puts on
    the axis. Falling back to ``O1/BF1`` ignores the spectrometer referencing
    (SR) and was off by 50.8 ppm on a real, correctly referenced 35Cl dataset;
    ``referenced=False`` tells the UI to say so rather than quietly shifting
    every chemical shift the user reads off.
    """
    sfo = float(meta.get("larmor_MHz", 0.0) or 0.0)
    sf = float(meta.get("sf_MHz", 0.0) or 0.0)
    if sf > 0 and sfo > 0:
        return (sfo - sf) * 1e6 / sf, True
    bf1 = float(meta.get("bf1_MHz", 0.0) or 0.0)
    return (float(meta.get("o1_Hz", 0.0)) / bf1 if bf1 else 0.0), False


# --------------------------------------------------------------------------
# Model-independent observables. For a broad, distribution-dominated pattern
# at a single field, delta_CG and the central-band width are the defensible
# numbers -- a single-field lineshape fit cannot separate delta_iso from the
# second-order quadrupolar shift (both distributed). See help/qcpmg.md.

def cg_window(ppm: np.ndarray, y: np.ndarray, *, rel_floor: float = 0.05
              ) -> tuple[float, float]:
    """(hi_ppm, lo_ppm) bracketing the central band, cut at the first intensity
    MINIMA either side of the tallest peak -- the reproducible window
    definition used for delta_CG. Falls back to a relative-height cut when no
    clean minimum exists, and is only ever a SEED: the caller should let the
    user drag it and must report the sensitivity (see :func:`centre_of_gravity`).
    """
    ppm = np.asarray(ppm, float); y = np.asarray(y, float)
    if ppm.size < 5:
        return (float(np.max(ppm)) if ppm.size else 0.0,
                float(np.min(ppm)) if ppm.size else 0.0)
    order = np.argsort(ppm)
    x, v = ppm[order], y[order]
    k = int(np.argmax(v))
    base = float(np.median(v)) + rel_floor * (float(v[k]) - float(np.median(v)))

    def edge(step: int) -> int:
        i = k
        while 0 < i < v.size - 1:
            j = i + step
            if not 0 <= j < v.size:
                break
            if v[j] > v[i] and v[i] <= base:      # a genuine minimum, low enough
                return i
            if v[j] <= base and abs(j - k) > 3 and v[j] < v[i]:
                i = j
                continue
            i = j
        return int(np.clip(i, 0, v.size - 1))

    lo_i, hi_i = edge(-1), edge(+1)
    if abs(hi_i - lo_i) < 5:
        # the peak sits on an edge (or the trace is flat/inverted) and the walk
        # never moved: a zero-width window would report sigma = 0, i.e. maximum
        # confidence, on the most degenerate input there is. Fall back to a
        # relative-height cut about the peak, as the docstring promises.
        above = np.where(v >= base)[0]
        if above.size >= 5:
            lo_i, hi_i = int(above[0]), int(above[-1])
        else:
            return float("nan"), float("nan")
    return float(max(x[lo_i], x[hi_i])), float(min(x[lo_i], x[hi_i]))


def centre_of_gravity(ppm: np.ndarray, y: np.ndarray,
                      window: tuple[float, float] | None = None, *,
                      jitter_frac: float = 0.10
                      ) -> tuple[float, float]:
    """Intensity-weighted centre of gravity (ppm) and its window sensitivity.

    Returns ``(delta_cg, sigma)`` where sigma is the spread of delta_CG when
    each window edge is jittered by ``jitter_frac`` of the window width -- a
    deterministic, reproducible replacement for "integrate it three times by
    hand and quote the standard deviation". sigma is also a QUALITY FLAG: a
    few ppm means the window is well defined; tens of ppm means the edges are
    running down a tail and should be placed by hand.
    """
    ppm = np.asarray(ppm, float); y = np.asarray(y, float)
    if window is None:
        window = cg_window(ppm, y)
    if not np.all(np.isfinite(window)) or max(window) <= min(window):
        return float("nan"), float("nan")
    hi, lo = max(window), min(window)
    width = hi - lo

    def cg(a: float, b: float) -> float:
        m = (ppm >= min(a, b)) & (ppm <= max(a, b))
        w = y[m]
        if w.size < 3:
            return float("nan")
        den = float(w.sum())
        # a window straddling equal +/- lobes has den ~ 0 and the ratio blows
        # up to values far outside the window itself -- refuse rather than
        # report a centre of gravity the data does not support
        if abs(den) <= 1e-3 * float(np.abs(w).sum()):
            return float("nan")
        val = float((ppm[m] * w).sum() / den)
        return val if min(a, b) <= val <= max(a, b) else float("nan")

    base = cg(lo, hi)
    d = jitter_frac * width
    trials = [cg(lo + sa * d, hi + sb * d)
              for sa in (-1, 0, 1) for sb in (-1, 0, 1)]
    trials = [t for t in trials if np.isfinite(t)]
    return base, (float(np.std(trials)) if len(trials) > 1 else 0.0)


def fwhm_hz(ppm: np.ndarray, y: np.ndarray, sfo_MHz: float,
            window: tuple[float, float] | None = None) -> float:
    """Full width at half maximum of the central band, in Hz."""
    ppm = np.asarray(ppm, float); y = np.asarray(y, float)
    if window is not None:
        m = (ppm >= min(window)) & (ppm <= max(window))
        ppm, y = ppm[m], y[m]
    if ppm.size < 3:
        return 0.0
    mx = float(np.max(y))
    if not np.isfinite(mx) or mx <= 0:
        return 0.0
    half = 0.5 * mx
    above = np.where(y >= half)[0]
    if above.size < 2:
        return 0.0
    # the OUTERMOST half-max crossings, by convention: a two-horned pattern
    # whose saddle dips below half still spans both horns. No noise-spike
    # pruning here -- any sample-count rule is resolution-dependent (the same
    # lineshape got a different FWHM at different zero-fills), and a spike
    # inside the window is the WINDOW's problem: it is draggable, and delta_CG
    # sigma already flags a badly placed one.
    return abs(float(ppm[above[-1]] - ppm[above[0]])) * (sfo_MHz or 0.0)


def overlay_pair(y_a: np.ndarray, y_b: np.ndarray, mode: str = "max"
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Scale two traces onto a common display scale so the sum-echo envelope
    and the spikelet manifold can be read against each other."""
    a = np.asarray(y_a, float); b = np.asarray(y_b, float)
    if mode == "area":
        na, nb = float(np.abs(a).sum()), float(np.abs(b).sum())
    elif mode == "raw":
        na = nb = 1.0
    else:
        na, nb = float(np.max(np.abs(a))), float(np.max(np.abs(b)))
    return a / (na or 1.0), b / (nb or 1.0)


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
    powder envelope (magnitude, spikelet-free).

    Delegates the transform to :func:`whole_echo_ft` so there is exactly ONE
    whole-echo convention in this module. It previously had its own copy,
    which kept the pre-fix one-sided window and end-appended zero-fill and so
    disagreed with ``sum_echo_spectrum`` (measured +53 % on the width of a
    known Lorentzian, and a 6 Hz peak shift)."""
    echo = coadd_echoes(fid, period, drop_first=drop_first)
    top = int(np.argmax(np.abs(echo)))
    ppm, spec = whole_echo_ft(echo, top, sw_Hz, sfo_MHz, carrier_ppm,
                              lb_Hz=lb_Hz, zf=zf)
    return ppm, np.abs(spec)


def spikelet_spacing_ppm(period: int, sw_Hz: float, sfo_MHz: float) -> float:
    return (sw_Hz / period) / (sfo_MHz or 1.0) if period else 0.0


# --------------------------------------------------------------------------
# Full ssNake-style "sum echo" workflow: split -> (T2 fit / weight) -> sum ->
# whole-echo processing (swap the echo top to t=0 -> a clean absorption
# lineshape you FIT, instead of the spikelet comb).

def echo_period_from_meta(meta: dict, n_points: int = 0) -> tuple[float, str]:
    """Echo period in POINTS, read from the pulse program -- exact, not guessed.

QCPMG pulse programs record the echo period, but not all in the same
    constant. Bruker's own use CNST7 (spikelet spacing, Hz) or CNST8 (points);
    the widely-circulated NMRFAM/Perras ``qcpmg.av4.nmrfam`` sequence writes
    CNST11 (spikelet spacing, Hz), CNST14 (points per echo) and CNST15 (echo
    period, us) instead, leaving CNST7 at its default of 1. Returns
    ``(period_points, source)`` naming the constant it used, so the UI can say
    whether the number was READ or GUESSED. ``period`` may be fractional
    (e.g. 292.9688) -- callers round, but see ``sum_echoes(realign=...)``.

    ``n_points`` (the train length) enables a sanity range. It is needed:
    CNST defaults to 1.0 on a Bruker dataset, so an untouched CNST7 would
    otherwise "read" a period of sw/1 = the whole sweep width. Anything
    outside 8..n_points/2 is rejected as not-a-period.
    """
    cnst = list(meta.get("cnst") or [])
    sw = float(meta.get("sw_Hz", 0.0) or 0.0)
    hi = (n_points / 2.0) if n_points else float("inf")

    def ok(p: float) -> bool:
        return bool(np.isfinite(p) and 8.0 <= p <= hi)

    if sw > 0 and len(cnst) > 7 and cnst[7] > 0:
        p = sw / float(cnst[7])
        if ok(p):
            return p, "CNST7"
    if len(cnst) > 8 and ok(float(cnst[8])):
        return float(cnst[8]), "CNST8"
    # NMRFAM / Perras qcpmg.av4.nmrfam: the same three facts, different slots
    if sw > 0 and len(cnst) > 11 and cnst[11] > 0:
        p = sw / float(cnst[11])
        if ok(p):
            return p, "CNST11"
    if sw > 0 and len(cnst) > 15 and cnst[15] > 0:
        p = float(cnst[15]) * 1e-6 * sw            # echo period in us
        if ok(p):
            return p, "CNST15"
    if len(cnst) > 14 and ok(float(cnst[14])):
        return float(cnst[14]), "CNST14"
    masr = float(meta.get("masr_Hz") or 0.0)
    if sw > 0 and masr > 0:                 # rotor-synchronised echo, 1 rotor
        p = sw / masr
        if ok(p):
            return p, "MASR"
    return 0.0, "none"


def centre_offset(fid: np.ndarray, period: int, *, tol_frac: float = 0.25
                  ) -> int:
    """Points to skip so each block holds ONE whole, CENTRED echo.

    Acquisition does not always begin half an echo before the first top. The
    NMRFAM sequence starts recording AT a top, so the natural blocks each
    hold the right half of one echo and the left half of the NEXT -- two
    different echoes, of different amplitude, glued into a fake one. The
    consequences are not subtle: on a real 35Cl train it put T2 at 4.0 ms
    instead of 10.3, demanded p1 = 493 deg and p2 = 331 deg to phase, and
    inflated the FWHM by 24 %.

    Returns 0 when the top is already within ``tol_frac`` of the centre, so a
    train that was acquired conventionally is left exactly as it was.
    """
    p = int(period)
    if p < 4 or 2 * p > np.asarray(fid).size:
        return 0
    top = echo_top_point(split_echoes(fid, p))
    centre = p // 2
    if abs(top - centre) <= tol_frac * p:
        return 0
    return int((top - centre) % p)


def split_echoes(fid: np.ndarray, period: int, *, first: int = 0,
                 n_echoes: int | None = None, drop_first: int = 0
                 ) -> np.ndarray:
    """Reshape the train into (n_echoes, period).

    ``first`` skips leading points, ``drop_first`` discards whole leading
    echoes (a partial first echo), ``n_echoes`` caps the count -- ssNake's
    "trim to an exact multiple of the echo length" step, made explicit
    instead of inferred from an all-zero test that never fires on real
    (noisy) data.
    """
    fid = np.asarray(fid, complex)
    period = int(period)
    if period < 2:
        raise ValueError(f"echo period must be >= 2 points (got {period})")
    fid = fid[int(first):]
    n = fid.size // period
    if n < 1:
        raise ValueError(
            f"echo period {period} exceeds the {fid.size}-point train")
    block = fid[:n * period].reshape(n, period)
    if drop_first:
        block = block[int(drop_first):]
        if block.shape[0] == 0:
            raise ValueError(
                f"drop_first={drop_first} leaves no echoes (train has {n})")
    if n_echoes is not None:
        block = block[:max(1, int(n_echoes))]
    return block


def n_usable_echoes(echoes: np.ndarray, top: int | None = None,
                    snr_floor: float = 3.0) -> int:
    """How many echoes rise above the noise floor -- seeds the "echoes to use"
    control so the user is not summing pure noise into the spectrum."""
    e = np.asarray(echoes)
    if e.ndim != 2 or e.shape[0] < 2:
        return 1
    if top is None:
        top = echo_top_point(e)
    amp = np.abs(e[:, int(top)])
    # Noise level = MEDIAN |.| of the last echo. Most points in an echo block
    # are baseline, so the median tracks the noise whether or not the train has
    # fully decayed (and if it hasn't, it over-estimates, which errs toward
    # keeping fewer echoes -- the safe direction). The smallest quartile of a
    # magnitude is near zero by construction and badly under-estimates it.
    sigma = float(np.median(np.abs(e[-1]))) or 1e-12
    good = np.where(amp > snr_floor * sigma)[0]
    return int(good[-1]) + 1 if good.size else 1


def echo_top_point(echoes: np.ndarray, method: str = "mean") -> int:
    """Point index of the echo top.

    ``method="mean"`` (default) takes the argmax of the COHERENT average,
    ``|mean(echoes)|`` -- the echoes add in phase at the top, so the average
    sharpens it. Averaging the MAGNITUDES instead adds the noise floor too and
    can land a point off the true top; one point matters enormously here
    (moving the top by 1 changed a measured T2 by up to 7400 % across a real
    12-sample set), so this is not a cosmetic choice.
    ``method="magnitude"`` and ``"mid"`` keep the older behaviours."""
    e = np.asarray(echoes, complex)
    if e.ndim == 1:
        return int(np.argmax(np.abs(e)))
    if method == "mid":
        ref = e[e.shape[0] // 2] if e.shape[0] > 2 else e[0]
        return int(np.argmax(np.abs(ref)))
    if method == "magnitude":
        return int(np.argmax(np.abs(e).mean(axis=0)))
    return int(np.argmax(np.abs(e.mean(axis=0))))


def echo_top_candidates(echoes: np.ndarray) -> dict[str, int]:
    """The independent echo-top estimates. When they disagree the top is
    genuinely ambiguous and the UI should say so rather than pick silently."""
    e = np.asarray(echoes, complex)
    return {"coherent": echo_top_point(e, "mean"),
            "magnitude": echo_top_point(e, "magnitude"),
            "first": int(np.argmax(np.abs(e[0])))}


def echo_centre(period: int) -> int:
    """Index of the echo block's centre. A whole echo needs its top HERE; the
    dialog shows both so the condition is visible rather than assumed."""
    return int(period) // 2


def first_last_echo(echoes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(first, last) echo magnitudes -- the ssNake split validation: their
    features (and the flat tail) must line up, or the period is wrong."""
    e = np.asarray(echoes)
    return np.abs(e[0]), np.abs(e[-1])


def split_alignment(echoes: np.ndarray) -> float:
    """0..1 score of how well the echoes line up (mean |correlation| of each
    echo with the mean echo shape). A wrong period drops this sharply, so the
    UI can flag a bad split instead of asking the user to eyeball it."""
    e = np.abs(np.asarray(echoes, complex))
    if e.ndim != 2 or e.shape[0] < 2:
        return 0.0
    ref = e.mean(axis=0)
    ref = ref - ref.mean()
    den_r = float(np.sqrt((ref * ref).sum())) or 1e-12
    scores = []
    for row in e:
        r = row - row.mean()
        den = float(np.sqrt((r * r).sum())) or 1e-12
        scores.append(abs(float((r * ref).sum()) / (den * den_r)))
    return float(np.mean(scores))


def echo_decay(echoes: np.ndarray, top: int, mode: str = "real") -> np.ndarray:
    """Echo-top intensity vs echo number -- the transverse (T2') decay.

    ``mode="real"`` (default) is the SIGNED real part, which is what ssNake
    samples and what reproduces published T2 values; ``"magnitude"`` has a
    rectified noise floor that biases the tail upward. The signed decay
    legitimately runs negative when the echo is not perfectly phased.
    """
    e = np.asarray(echoes)
    top = int(top)
    if not 0 <= top < e.shape[-1]:
        raise ValueError(f"decay point {top} outside the echo "
                         f"(0..{e.shape[-1] - 1})")
    col = e[:, top]
    if mode == "magnitude":
        return np.abs(col)
    if mode == "complex":
        return col
    return np.real(col)


@dataclass(frozen=True)
class T2Fit:
    """A T2 measurement, on BOTH time axes (see the module docstring)."""
    T2_s: float                # physical (tau_echo per echo)
    T2_err_s: float
    const: float               # ssNake "Cst"
    coeff: float               # ssNake "Coeff"
    r2: float
    n_used: int
    ok: bool
    pinned: bool
    T2_ssnake_s: float         # on ssNake's D1 pseudo-axis (1 dwell per echo)
    lb_Hz: float               # matched filter, physical  = 1/(pi*T2)
    lb_ssnake_Hz: float        # the number you type into ssNake

    def model(self, t_s):
        """The fitted curve on the physical time axis."""
        t = np.asarray(t_s, float)
        return self.const + self.coeff * np.exp(-t / self.T2_s)


def autophase_best(spec: np.ndarray, *, gain: float = 0.25):
    """(p0, p1, p2) using the LOWEST order that actually phases the spectrum.

    Fits p0/p1, then p0/p1/p2, and keeps the quadratic term only when it
    reduces the residual negative area by more than ``gain`` (25 %) -- so an
    ordinary echo is phased exactly as before and reports p2 = 0, while a
    frequency-swept (WURST/chirp) dataset gets the term it genuinely needs.
    On a real WCPMG 81Br sample this took the negative dips from -47 % to
    -3.8 %, and brought the phased delta_CG into agreement with the
    magnitude-mode value (-314 vs -310 ppm) where p0/p1 alone gave -94.
    """
    def neg_frac(y: np.ndarray) -> float:
        mx = float(np.abs(y).max()) or 1.0
        return float(np.abs(y[y < 0]).sum()) / (y.size * mx)

    p0, p1 = autophase(spec, order=1)
    base = neg_frac(phase_spectrum(spec, p0, p1).real)
    q0, q1, q2 = autophase(spec, order=2)
    quad = neg_frac(phase_spectrum(spec, q0, q1, p2_deg=q2).real)
    if base > 0 and quad < (1.0 - float(gain)) * base:
        return float(q0), float(q1), float(q2)
    return float(p0), float(p1), 0.0


def matched_lb_Hz(T2_s: float) -> float:
    """The matched-filter Lorentzian broadening, LB = 1/(pi*T2)."""
    return 1.0 / (np.pi * T2_s) if T2_s and T2_s > 0 else 0.0


def fit_t2(tau_s: float, decay: np.ndarray, *, offset: bool = True,
           period: int | None = None, n_fit: int | None = None,
           t_s: np.ndarray | None = None) -> T2Fit:
    """Fit ``C + B*exp(-t/T2)`` to the echo-top decay (ssNake's model).

    ``tau_s`` is the echo spacing. ``offset=True`` fits the constant C, which
    matters: on real data B/C is ~18 (median over 12 samples) and dropping the
    constant moves T2 by ~26 % (median). ``period`` (points per echo) fills in
    the ssNake pseudo-axis fields.

    ``t_s`` gives the time of each supplied point EXPLICITLY. Pass it whenever
    ``decay`` is not a contiguous run of echoes -- e.g. after excluding an
    outlier -- otherwise the remaining points are silently re-timed onto
    ``0, tau, 2*tau, ...`` and every later echo is credited with an earlier
    time (excluding the second echo of a real train shifted T2 by -36 %).

    Never fabricates a number: a fit that failed, pinned a bound, or is simply
    not supported by the data comes back ``ok=False`` so the UI can refuse to
    offer a matched filter derived from it.
    """
    from scipy.optimize import curve_fit

    d = np.asarray(decay, float)
    t_in = None if t_s is None else np.asarray(t_s, float)
    if n_fit:
        d = d[:max(4, int(n_fit))]
        if t_in is not None:
            t_in = t_in[:max(4, int(n_fit))]
    n = d.size
    if t_in is not None and t_in.size != n:
        raise ValueError(f"t_s has {t_in.size} points but decay has {n}")
    bad = T2Fit(0.0, float("nan"), 0.0, 0.0, 0.0, n, False, False,
                0.0, 0.0, 0.0)
    if n < 4 or not np.all(np.isfinite(d)):
        return bad
    scale = float(np.max(np.abs(d))) or 1.0
    # No dynamic range = no decay to measure. Without this an all-zero or flat
    # decay fits perfectly (R2 = 1) and reports a confident, meaningless T2.
    if float(np.ptp(d)) <= 1e-9 * scale:
        return bad
    y = d / scale
    t = np.arange(n) * tau_s if t_in is None else t_in
    lo, hi = tau_s, tau_s * n * 20.0

    def model(tt, C, B, T2):
        return C + B * np.exp(-tt / T2)

    p0 = [y[-1] if offset else 0.0, y[0] - y[-1], tau_s * max(3.0, n / 3.0)]
    bounds = ([-np.inf, -np.inf, lo], [np.inf, np.inf, hi])
    if not offset:
        bounds[0][0], bounds[1][0] = -1e-12, 1e-12
        p0[0] = 0.0
    try:
        popt, cov = curve_fit(model, t, y, p0=p0, bounds=bounds, maxfev=20000)
    except Exception:
        return bad
    C, B, T2 = (float(v) for v in popt)
    if not np.isfinite(T2) or T2 <= 0:
        return bad
    err = float(np.sqrt(np.diag(cov))[2]) if np.all(np.isfinite(cov)) else float("nan")
    pinned = bool(T2 <= lo * 1.001 or T2 >= hi * 0.9)
    resid = y - model(t, C, B, T2)
    ss_tot = float(((y - y.mean()) ** 2).sum()) or 1e-30
    r2 = float(1.0 - (resid ** 2).sum() / ss_tot)
    per = int(period) if period else 0
    t2_ss = T2 / per if per else 0.0
    # "ok" means USABLE, not merely "curve_fit returned". Pure noise converges
    # happily (R2 ~ 0.02) and would otherwise hand the user a matched filter
    # computed from nothing; r2 and pinned stay on the dataclass so the UI can
    # still distinguish "did not converge" from "converged but meaningless".
    ok = bool((not pinned) and np.isfinite(err) and r2 > 0.5
              and err < 0.5 * T2)
    return T2Fit(T2_s=T2, T2_err_s=err, const=C * scale, coeff=B * scale,
                 r2=r2, n_used=n, ok=ok, pinned=pinned, T2_ssnake_s=t2_ss,
                 lb_Hz=matched_lb_Hz(T2), lb_ssnake_Hz=matched_lb_Hz(t2_ss))


def apodize_echoes(echoes: np.ndarray, tau_s: float, t2_s: float | None,
                   *, normalise: bool = False) -> np.ndarray:
    """The T2-weighted echo MATRIX -- ssNake's "Apodised echoes" figure.
    ``sum_echoes`` collapses this immediately, so without it the weighting
    can never be plotted."""
    e = np.asarray(echoes, complex)
    w = echo_weights(e.shape[0], tau_s, t2_s, normalise=normalise)
    return e * w[:, None]


def echo_weights(n: int, tau_s: float, t2_s: float | None,
                 *, normalise: bool = False) -> np.ndarray:
    """The per-echo matched-filter weights exp(-k*tau/T2) (ssNake's "Apodised
    D1" curve). Unit weights when no T2 is given."""
    w = np.ones(int(n))
    if t2_s and t2_s > 0:
        w = np.exp(-np.arange(int(n)) * tau_s / t2_s)
    if normalise:
        w = w / (w.sum() or 1.0)
    return w


def sum_echoes(echoes: np.ndarray, tau_s: float,
               t2_weight_s: float | None = None, *, normalise: bool = True,
               realign: bool = False, top: int | None = None) -> np.ndarray:
    """Coherently add the echoes. With ``t2_weight_s`` set, weight echo k by
    exp(-k·tau/T2) -- the matched filter that maximises S/N (ssNake's T2
    weighting via a Lorentzian LB = 1/(πT2) along the echo dimension).

    ``normalise`` divides by the weight sum so absolute intensity does not
    jump when weighting is toggled (it changed by 0.56x before). ``realign``
    rolls each echo onto the reference top first -- insurance against a
    fractional echo period drifting by a point across a long train.
    """
    e = np.asarray(echoes, complex)
    n = e.shape[0]
    if realign:
        ref = int(top) if top is not None else echo_top_point(e)
        e = np.stack([np.roll(row, ref - int(np.argmax(np.abs(row))))
                      for row in e])
    w = echo_weights(n, tau_s, t2_weight_s)
    out = (e * w[:, None]).sum(axis=0)
    return out / (w.sum() or 1.0) if normalise else out


def _circ_index(n: int, npos: int | None) -> np.ndarray:
    """|t| in points about index 0, wrapping at ``npos`` (the number of
    non-negative-time samples). Samples from ``npos`` on are negative times,
    so their |t| counts back from n. ``npos=None`` assumes the top sits at
    the block centre, i.e. the plain circular distance min(k, n-k)."""
    k = np.arange(n)
    if npos is None:
        return np.minimum(k, n - k)
    npos = int(np.clip(npos, 1, n))
    return np.where(k < npos, k, n - k)


def _gaussian_apod(n: int, sw_Hz: float, gb_Hz: float,
                   whole_echo: bool = False, npos: int | None = None) -> np.ndarray:
    """Gaussian window. With ``whole_echo`` the window is symmetric about
    index 0 in the CIRCULAR sense (``min(k, n-k)``).

    That symmetry is not cosmetic. After the whole-echo swap the echo top sits
    at index 0 and the echo's LEFT half is wrapped to the END of the array; a
    one-sided window decaying from index 0 multiplies that half by ~0, i.e.
    throws away half the echo (measured: 44 % of |signal| at gb=2 kHz on a
    512-point echo) and destroys the very symmetry whole-echo processing
    exists to exploit."""
    if not gb_Hz:
        return np.ones(n)
    k = _circ_index(n, npos) if whole_echo else np.arange(n)
    t = k / sw_Hz
    return np.exp(-((np.pi * gb_Hz * t) ** 2) / (4.0 * np.log(2.0)))


def _lorentz_apod(n: int, sw_Hz: float, lb_Hz: float,
                  whole_echo: bool = False, npos: int | None = None) -> np.ndarray:
    """Lorentzian (exponential) window; ``whole_echo`` as in _gaussian_apod."""
    if not lb_Hz:
        return np.ones(n)
    k = _circ_index(n, npos) if whole_echo else np.arange(n)
    return np.exp(-np.pi * lb_Hz * (k / sw_Hz))


def _zerofill_whole_echo(sig: np.ndarray, nfft: int,
                         npos: int | None = None) -> np.ndarray:
    """Zero-fill a whole echo whose top is already at index 0.

    The zeros MUST go in the middle (at maximum |t|), not at the end: the
    samples at the end of ``sig`` are NEGATIVE times, and appending zeros
    after them re-interprets them as large positive times. That single
    mistake turns a pure absorption spectrum into a 60 %-dispersive one
    (measured |imag|max/|real|max 0.603 -> 0.000 once fixed)."""
    sig = np.asarray(sig, complex)
    m = sig.size
    if nfft <= m:
        return sig
    # after roll(-top) there are exactly (m - top) non-negative-time samples,
    # NOT m//2; the two coincide only when the top sits at the block centre
    half = int(np.clip(m // 2 if npos is None else npos, 1, m))
    out = np.zeros(nfft, complex)
    out[:half] = sig[:half]
    out[nfft - (m - half):] = sig[half:]        # negative times stay at the end
    return out


def whole_echo_ft(echo: np.ndarray, top: int, sw_Hz: float, sfo_MHz: float,
                  carrier_ppm: float = 0.0, *, lb_Hz: float = 0.0,
                  gb_Hz: float = 0.0, zf: int = 16,
                  ) -> tuple[np.ndarray, np.ndarray]:
    """Whole-echo process ONE echo -> (ppm, complex spectrum).

    The single place the swap/apodize/zero-fill convention lives: roll the
    echo top to t=0, apodize SYMMETRICALLY about it, zero-fill MID-ARRAY, FT.
    Done correctly the result is essentially pure absorption, so the caller
    normally needs only a zero-order phase."""
    echo = np.asarray(echo, complex)
    m = echo.size
    if not 0 <= int(top) < m:
        raise ValueError(f"echo top {top} outside the echo (0..{m - 1})")
    npos = m - int(top)                # samples at t >= 0 after the roll
    echo = np.roll(echo, -int(top))                   # whole echo: top -> t=0
    echo = echo * _lorentz_apod(m, sw_Hz, lb_Hz, whole_echo=True, npos=npos)
    echo = echo * _gaussian_apod(m, sw_Hz, gb_Hz, whole_echo=True, npos=npos)
    nfft = int(2 ** np.ceil(np.log2(max(m, 2) * max(1, zf))))
    spec = np.fft.fftshift(np.fft.fft(_zerofill_whole_echo(echo, nfft, npos)))
    ppm = _axis_ppm(nfft, sw_Hz, sfo_MHz, carrier_ppm)
    order = np.argsort(ppm)
    return ppm[order], spec[order]


def sum_echo_spectrum(fid: np.ndarray, period: int, sw_Hz: float, sfo_MHz: float,
                      carrier_ppm: float = 0.0, top: int | None = None,
                      t2_weight_s: float | None = None, p0_deg: float = 0.0,
                      p1_deg: float = 0.0, gb_Hz: float = 0.0, zf: int = 16,
                      *, first: int = 0, n_echoes: int | None = None,
                      drop_first: int = 0, lb_Hz: float = 0.0,
                      normalise: bool = True, realign: bool = False,
                      ) -> tuple[np.ndarray, np.ndarray]:
    """The fittable QCPMG spectrum: sum the echoes, whole-echo process (swap the
    top to t=0), apodize, FT, and apply a p0/p1 phase. Returns
    (ppm, complex spectrum); take .real for the absorption lineshape to fit."""
    echoes = split_echoes(fid, period, first=first, n_echoes=n_echoes,
                          drop_first=drop_first)
    if top is None:
        top = echo_top_point(echoes)
    tau = period / sw_Hz
    summed = sum_echoes(echoes, tau, t2_weight_s, normalise=normalise,
                        realign=realign, top=top)
    ppm, spec = whole_echo_ft(summed, int(top), sw_Hz, sfo_MHz, carrier_ppm,
                              lb_Hz=lb_Hz, gb_Hz=gb_Hz, zf=zf)
    return ppm, phase_spectrum(spec, p0_deg, p1_deg)


def phase_spectrum(spec: np.ndarray, p0_deg: float = 0.0, p1_deg: float = 0.0,
                   pivot_frac: float = 0.5, p2_deg: float = 0.0) -> np.ndarray:
    """Apply a p0/p1/p2 phase, pivoted at ``pivot_frac`` of the spectrum
    (0.5 = centre, matching larmor.processing.op_phase). Pivoting on the left
    edge instead makes the orders fight each other.

    ``p2_deg`` is the SECOND-order (quadratic) term, in degrees at the ends of
    the axis. A frequency-swept refocusing pulse (WURST/chirp, as in
    WURST-CPMG) imprints exactly such a quadratic phase across the swept band,
    and no amount of p0/p1 removes it: on a real WCPMG 81Br dataset the
    residual negative dips fall from -47 % to about -4 % once p2 is allowed.
    """
    n = np.asarray(spec).size
    ramp = np.arange(n) / max(n - 1, 1) - float(pivot_frac)
    phi = (np.deg2rad(p0_deg) + np.deg2rad(p1_deg) * ramp
           + np.deg2rad(p2_deg) * (2.0 * ramp) ** 2)
    return spec * np.exp(-1j * phi)


def autophase0(spec: np.ndarray) -> float:
    """Zero-order phase (deg) that maximises the real integral.

    Scored with ``exp(-1j*p)`` -- the SAME sign convention every consumer
    applies it with. Scoring the opposite sign returned the negated angle, so
    applying the result actually INVERTED the spectrum."""
    ph = np.linspace(-np.pi, np.pi, 721)
    scores = [np.real(spec * np.exp(-1j * p)).sum() for p in ph]
    return float(np.degrees(ph[int(np.argmax(scores))]))


def autophase(spec: np.ndarray, order: int = 1):
    """Phase (deg) for a mostly-absorptive lineshape.

    Minimises the area of the negative part of the real spectrum (a robust
    criterion for a powder pattern that should be all-positive), with a light
    penalty on p1 to avoid runaway first-order twists.

    ``order=1`` returns (p0, p1) -- the default, unchanged. ``order=2``
    returns (p0, p1, p2) and is what a frequency-swept (WURST/chirp)
    refocusing pulse needs: its quadratic phase is invisible to p0/p1.
    """
    from scipy.optimize import minimize

    spec = np.asarray(spec)
    n = spec.size
    ramp = np.arange(n) / max(n - 1, 1) - 0.5          # pivot at the centre
    norm = np.abs(spec).sum() or 1.0

    quad = (2.0 * ramp) ** 2

    def penalty(ph):
        p0, p1 = ph[0], ph[1]
        p2 = ph[2] if len(ph) > 2 else 0.0
        r = np.real(spec * np.exp(-1j * (np.deg2rad(p0) + np.deg2rad(p1) * ramp
                                         + np.deg2rad(p2) * quad)))
        # maximise (positive area - 4x negative area): minimising the negative
        # area ALONE is degenerate on a non-negative pattern -- it is flat over
        # a wide window of p0, so the optimiser could stop anywhere in it.
        score = (r.sum() - 4.0 * np.abs(r[r < 0]).sum()) / norm
        # keep the p1 penalty far below the score scale: at 1e-6 a line
        # sitting at ~25% of the axis (legitimate p1 ~ 200 deg) was paying
        # 0.04 -- enough to drag p1 visibly toward 0. 1e-8 only breaks ties.
        return -score + 1e-8 * (p1 ** 2)

    p0_0 = autophase0(spec)
    best = None
    # a few p1 (and p2) starts so we don't fall into a local twist
    p2_starts = (0.0,) if order < 2 else (0.0, 180.0, -180.0, 360.0, -360.0)
    for p1_0 in (0.0, 180.0, -180.0, 360.0, -360.0):
        for p2_0 in p2_starts:
            x0 = [p0_0, p1_0] if order < 2 else [p0_0, p1_0, p2_0]
            res = minimize(penalty, x0, method="Nelder-Mead",
                           options={"xatol": 0.5, "fatol": 1e-4,
                                    "maxiter": 800 * (2 if order > 1 else 1)})
            if best is None or res.fun < best.fun:
                best = res
    p0 = ((best.x[0] + 180) % 360) - 180
    if order < 2:
        return float(p0), float(best.x[1])
    return float(p0), float(best.x[1]), float(best.x[2])
