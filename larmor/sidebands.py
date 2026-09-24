"""Spinning-sideband detection from the data itself.

Under magic-angle spinning every anisotropic pattern breaks into a centreband
and a manifold of sidebands at delta_iso +- k nu_rot / nu_0 (ppm, k > 0 towards
higher shift -- the convention of larmor.herzfeld_berger). The user used to
have to SPOT that repeat and then open "Add a copy of this spectrum..." twice.
This module finds it: the spectrum is correlated with a shifted copy of
itself, the lag at which it repeats is the sideband spacing, and the comb of
peaks is then verified on the trace before anything is offered.

Method
------
* The trace is sorted, block-averaged / resampled onto a uniform grid of at
  most MAX_POINTS, the edge-median floor (:func:`larmor.qcpmg.noise_floor`,
  never the global median) is subtracted and negatives are clipped.
* One zero-padded rfft/irfft gives the LINEAR autocorrelation at every lag;
  each lag is normalised by the energies of the two overlapping segments, so
  c[0] == 1, |c| <= 1 and a self-similar repeat scores ~1 however much of the
  manifold falls off the axis. The profile is cut at half the span: a
  manifold with +-1 on both sides needs an axis of at least twice its spacing.
* With a recorded rate the spacing is the INTERIOR maximum of c inside
  spacing (1 +- RATE_TOLERANCE); a maximum on the window boundary is a slope,
  not a repeat (the lesson of ``qcpmg.detect_period``). Without a usable rate
  (``scan=True``) every interior maximum between the lags of SCAN_HZ is a
  candidate, scored by its prominence over the local median of c, and the
  SMALLEST candidate within 90 % of the best prominence wins -- harmonics of
  the true spacing always score well, sub-multiples collapse, so the smallest
  good one is the fundamental (``qcpmg.find_period_by_correlation``).
* The spacing is then VERIFIED on the spectrum: scipy.signal.find_peaks with a
  prominence floor (not a height threshold -- noise wiggles on the body of a
  broad static pattern are local maxima above any height fraction), the
  tallest peak is the centreband, and orders k = 1, 2, ... are matched on each
  side. The autocorrelation is an even function and cannot tell a manifold
  from two unrelated lines or a two-horn static pattern; the symmetric
  criterion -- order +1 AND order -1 both matched -- can, and is required.

Everything returns a :class:`SidebandDetection` carrying ``ok`` and a
``message`` (the shape of ``staticct.StaticReading``), so a caller can say WHY
nothing was offered instead of staying silent. Nothing here raises on
degenerate input.

Qt-free; the desktop banner (larmor/desktop/sideband_offer.py) and the main
window are consumers.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks, peak_widths

from larmor.qcpmg import noise_floor

#: half-width of the search window around a recorded rate (fraction)
RATE_TOLERANCE = 0.02
#: autocorrelation prominence that counts as full confidence
PROMINENCE_FULL = 0.15
#: below this confidence a detection is reported but not ``ok``
MIN_CONFIDENCE = 0.5
#: rate range searched when the recorded rate is unknown or wrong (Hz)
SCAN_HZ = (1000.0, 80000.0)
#: the trace is reduced to at most this many points before correlating
MAX_POINTS = 16384
#: highest sideband order followed on each side
MAX_ORDER = 12
#: a peak must rise above this fraction of the tallest one to be a tooth
PEAK_FRAC = 0.02
#: a tooth must be RESOLVED: prominence >= this fraction of its height (a
#: bump riding on the body of a broad pattern is not)
RESOLVED_FRAC = 0.3
#: measured vs recorded rate: beyond this the banner warns / the model line
#: writes the rate
RATE_MISMATCH = 0.005
#: confidence penalty when the rate had to be scanned rather than confirmed
SCAN_PENALTY = 0.85

__all__ = [
    "RATE_TOLERANCE", "PROMINENCE_FULL", "MIN_CONFIDENCE", "SCAN_HZ",
    "MAX_POINTS", "MAX_ORDER", "PEAK_FRAC", "RESOLVED_FRAC", "RATE_MISMATCH",
    "SCAN_PENALTY", "SidebandOrder", "SidebandDetection",
    "correlation_profile", "shift_correlation", "detect", "seed_ratio",
    "describe", "format_hz",
]


@dataclass(frozen=True)
class SidebandOrder:
    #: signed order: +1 is the first sideband towards higher shift
    k: int
    #: expected position, centre + k * spacing (ppm)
    ppm: float
    #: floor-subtracted height of the matched peak, else the trace at ``ppm``
    height: float
    #: a resolved peak sits within tolerance of the expected position
    matched: bool


@dataclass
class SidebandDetection:
    ok: bool
    message: str
    nu_rot_Hz: float = 0.0
    spacing_ppm: float = 0.0
    centre_ppm: float = 0.0
    centre_fwhm_ppm: float = 0.0
    centre_height: float = 0.0
    #: both signs, |k| ascending
    orders: tuple = ()
    #: highest |k| matched on either side
    n_orders: int = 0
    #: normalised autocorrelation at the spacing
    score: float = 0.0
    #: score above the local median of the profile
    prominence: float = 0.0
    confidence: float = 0.0
    #: the rate was found by scanning, not confirmed around a recorded one
    scanned: bool = False
    #: the rate was confirmed within RATE_TOLERANCE of the recorded one
    from_meta: bool = False

    def matched(self, sign: int | None = None) -> list:
        """The matched orders, optionally of one sign (+1 / -1)."""
        return [o for o in self.orders
                if o.matched and (sign is None or (o.k > 0) == (sign > 0))]

    def sides(self) -> tuple:
        """The signs whose first order is matched: (1, -1), (1,), (-1,) or ()."""
        return tuple(s for s in (1, -1)
                     if any(o.matched and o.k == s for o in self.orders))


def linked_position_expr(parent: int, k: int, spacing_ppm: float) -> str:
    """The position constraint of a linked sideband copy of order ``k``:
    ``s<parent>.isotropic_chemical_shift_ppm ± |k|·spacing`` -- the form the
    table shows as ``A+124.6`` / ``A-124.6`` (cellparse.format_link)."""
    off = float(k) * float(spacing_ppm)
    return (f"s{int(parent)}.isotropic_chemical_shift_ppm "
            f"{'+' if off >= 0 else '-'} {abs(off):.6g}")


def refresh_linked(sites: list, larmor_MHz: float, spin_rate_Hz: float) -> int:
    """Move every linked sideband copy to the CURRENT spin rate.

    Recipe-dict sites carrying ``site["sideband"] = {"parent": i, "k": k}``
    get their position constraint and value recomputed from
    ``spin_rate_Hz / larmor_MHz``; the offset used to be a constant frozen
    at creation, so changing νrot afterwards (Experiment dialog, the
    detector's **Use**, a fit kept onto a new spectrum) left the copies on
    the old comb. A copy whose position the user has since unlinked or
    re-linked to another line loses its marker and is left alone; sites
    without a marker are never touched. Rate or field unknown (≤ 0):
    nothing moves. Returns the number of copies moved."""
    import re

    lar, nu = float(larmor_MHz or 0.0), float(spin_rate_Hz or 0.0)
    if lar <= 0 or nu <= 0 or not sites:
        return 0
    spacing = nu / lar
    moved = 0
    for i, s in enumerate(sites):
        mark = s.get("sideband") if isinstance(s, dict) else None
        if not mark:
            continue
        try:
            parent, k = int(mark["parent"]), int(mark["k"])
        except (KeyError, TypeError, ValueError):
            s.pop("sideband", None)
            continue
        pos = (s.get("params") or {}).get("isotropic_chemical_shift_ppm")
        expr = pos.get("expr") if isinstance(pos, dict) else None
        still_linked = (
            parent != i and 0 <= parent < len(sites) and expr
            and re.fullmatch(
                rf"\s*s{parent}\.isotropic_chemical_shift_ppm\s*[-+]\s*[0-9.eE+-]+\s*",
                expr) is not None)
        if not still_linked:
            s.pop("sideband", None)                 # unlinked by hand: theirs now
            continue
        new_expr = linked_position_expr(parent, k, spacing)
        p_pos = (sites[parent].get("params") or {}).get(
            "isotropic_chemical_shift_ppm", {})
        if new_expr != expr:
            moved += 1
        pos["expr"] = new_expr
        pos["stderr"] = None
        if isinstance(p_pos, dict) and p_pos.get("value") is not None:
            pos["value"] = float(p_pos["value"]) + k * spacing
    return moved


def format_hz(value: float) -> str:
    """20000 -> '20 000' (thousands separated by a space, no decimals)."""
    return f"{float(value):,.0f}".replace(",", " ")


# --------------------------------------------------------------- preparation
def _uniform(x, y, max_points: int):
    """Sorted, finite, block-averaged to <= max_points and resampled onto a
    uniform grid when the axis is not one. Returns (xu, yu, dx)."""
    x = np.asarray(x, float).ravel()
    y = np.real(np.asarray(y)).astype(float).ravel()
    n = min(x.size, y.size)
    x, y = x[:n], y[:n]
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    n = x.size
    if n > max_points:                     # average blocks, do not point-sample
        f = int(np.ceil(n / max_points))
        m = (n // f) * f
        x = x[:m].reshape(-1, f).mean(axis=1)
        y = y[:m].reshape(-1, f).mean(axis=1)
        n = x.size
    d = np.diff(x)
    if n < 3 or x[-1] <= x[0]:
        return x, y, float(d.mean()) if d.size else 0.0
    if d.min() > 0 and (d.max() - d.min()) <= 1e-3 * d.mean():
        return x, y, float(d.mean())
    xu = np.linspace(x[0], x[-1], n)
    return xu, np.interp(xu, x, y), float((x[-1] - x[0]) / (n - 1))


def _noise_rms(y: np.ndarray) -> float:
    """RMS of the quiet outer 10 % after median removal (the S/N convention)."""
    edge = max(5, y.size // 10)
    region = np.concatenate([y[:edge], y[-edge:]])
    return float(np.std(region - np.median(region)))


def _prepared(x, y, max_points: int):
    """(xu, floor-subtracted trace, clipped trace, dx, noise_rms)."""
    xu, yu, dx = _uniform(x, y, max_points)
    if xu.size < 3:
        return xu, yu, yu, dx, 0.0
    z = yu - noise_floor(yu)
    return xu, z, np.clip(z, 0.0, None), dx, _noise_rms(yu)


def _profile(z: np.ndarray, dx: float):
    """Overlap-normalised linear autocorrelation of the clipped trace at lags
    0 .. n//2 - 1 (half the span). Returns (lags_ppm, c)."""
    n = z.size
    n_lag = max(2, n // 2)
    lags = np.arange(n_lag) * dx
    e = z * z
    cs = np.cumsum(e)
    total = float(cs[-1]) if cs.size else 0.0
    if total <= 0.0:
        return lags, np.zeros(n_lag)
    m = 1 << int(np.ceil(np.log2(2 * n)))
    spec = np.fft.rfft(z, m)
    ac = np.fft.irfft(np.abs(spec) ** 2, m)[:n_lag]
    L = np.arange(n_lag)
    e_head = cs[n - 1 - L]                            # energy of z[0 .. n-1-L]
    cs_before = np.concatenate([[0.0], cs[:-1]])      # cs[L-1], cs[-1] := 0
    e_tail = total - cs_before[L]                     # energy of z[L .. n-1]
    denom = np.sqrt(e_head * e_tail)
    valid = (e_head > 1e-6 * total) & (e_tail > 1e-6 * total) & (denom > 0)
    c = np.zeros(n_lag)
    c[valid] = ac[valid] / denom[valid]
    return lags, np.clip(c, -1.0, 1.0)


def correlation_profile(x_ppm, y, *, max_points: int = MAX_POINTS):
    """(lags_ppm ascending from 0, c): the overlap-normalised linear
    autocorrelation of the floor-subtracted, clipped spectrum, c[0] == 1 and
    |c| <= 1, zero where an overlap energy is below 1e-6 of the total, cut at
    half the axis span."""
    _xu, _z, zc, dx, _noise = _prepared(x_ppm, y, max_points)
    return _profile(zc, dx)


def shift_correlation(x_ppm, y, shift_ppm: float) -> float:
    """Direct (np.interp) reference implementation of one lag of the same
    quantity as :func:`correlation_profile`: the spectrum correlated with a
    copy of itself displaced by ``shift_ppm``, over their overlap."""
    xu, _z, zc, _dx, _noise = _prepared(x_ppm, y, MAX_POINTS)
    if xu.size < 3:
        return 0.0
    xs = xu + float(shift_ppm)
    inside = (xs >= xu[0]) & (xs <= xu[-1])
    if inside.sum() < 2:
        return 0.0
    a = zc[inside]
    b = np.interp(xs[inside], xu, zc)
    denom = float(np.sqrt(np.sum(a * a) * np.sum(b * b)))
    return float(np.sum(a * b) / denom) if denom > 0 else 0.0


# ------------------------------------------------------------- lag search
def _refine(c: np.ndarray, i: int) -> float:
    """Parabolic vertex through c[i-1], c[i], c[i+1] as a fractional index
    (the sub-point formula processing.pick_peaks uses)."""
    if i <= 0 or i >= c.size - 1:
        return float(i)
    d = c[i - 1] - 2.0 * c[i] + c[i + 1]
    delta = 0.5 * (c[i - 1] - c[i + 1]) / d if d != 0 else 0.0
    return float(i + np.clip(delta, -1.0, 1.0))


def _vertex_height(c: np.ndarray, i: int) -> float:
    """Height of the same parabola at its vertex."""
    if i <= 0 or i >= c.size - 1:
        return float(c[i])
    d = c[i - 1] - 2.0 * c[i] + c[i + 1]
    delta = 0.5 * (c[i - 1] - c[i + 1]) / d if d != 0 else 0.0
    delta = float(np.clip(delta, -1.0, 1.0))
    return float(c[i] - 0.25 * (c[i - 1] - c[i + 1]) * delta)


def _interior_max(c: np.ndarray, lo: int, hi: int):
    """Index of the maximum of c[lo..hi] if it is strictly inside the window
    and rises above the window's floor; None for a slope or a flat."""
    lo, hi = max(int(lo), 0), min(int(hi), c.size - 1)
    if hi - lo < 2:
        return None
    seg = c[lo:hi + 1]
    i = lo + int(np.argmax(seg))
    if i == lo or i == hi:
        return None
    if seg.max() - seg.min() < 1e-6:
        return None
    return i


def _prominence(lags: np.ndarray, c: np.ndarray, i: int) -> float:
    """c[i] above the median of c over [0.5, 1.5] x lag, excluding +-10 %
    around the lag itself -- a pedestal (rectified noise) or a broad lobe
    contributes nothing, a genuine repeat all of its height."""
    lag = lags[i]
    sel = ((lags >= 0.5 * lag) & (lags <= 1.5 * lag)
           & (np.abs(lags - lag) > 0.1 * lag))
    base = float(np.median(c[sel])) if sel.sum() >= 4 else 0.0
    return float(c[i] - base)


def _scan(lags: np.ndarray, c: np.ndarray, lo_ppm: float, hi_ppm: float):
    """(index, prominence) of the fundamental repeat between the two lags, or
    None. Candidates are the interior local maxima of the profile; a
    candidate within +-20 % of a more prominent one is the same peak (a noise
    wiggle on its flank) and is dropped; of the survivors within 90 % of the
    best prominence the SMALLEST lag wins."""
    if c.size < 3:
        return None
    i = np.arange(1, c.size - 1)
    is_max = (c[i] > c[i - 1]) & (c[i] >= c[i + 1])
    cand = i[is_max]
    cand = cand[(lags[cand] >= lo_ppm) & (lags[cand] <= hi_ppm)]
    if cand.size == 0:
        return None
    prom = np.array([_prominence(lags, c, int(j)) for j in cand])
    keep_idx, keep_prom = [], []
    for j in np.argsort(-prom):                       # strongest first
        lag = lags[cand[j]]
        if any(0.8 * lags[k] <= lag <= 1.25 * lags[k] for k in keep_idx):
            continue
        keep_idx.append(int(cand[j]))
        keep_prom.append(float(prom[j]))
    best = max(keep_prom)
    if best <= 1e-3:                                   # ripples, not a repeat
        return None
    good = [(k, p) for k, p in zip(keep_idx, keep_prom) if p >= 0.9 * best]
    k, p = min(good)                                   # smallest lag
    return k, p


# ---------------------------------------------------------------- the comb
def _comb(xu, z, dx, spacing, centre_ppm, noise):
    """Verify the spacing on the floor-subtracted trace. Returns
    (centre_ppm, centre_fwhm, centre_height, orders) or None when no peak
    stands out of the noise."""
    top = float(z.max())
    if top <= 0.0:
        return None
    floor = max(PEAK_FRAC * top, 3.0 * noise, 1e-12 * top)
    dist = max(1, int(round(0.3 * spacing / dx)))
    peaks, props = find_peaks(z, prominence=floor, distance=dist)
    if peaks.size == 0:
        return None
    pos = np.array([xu[0] + _refine(z, int(i)) * dx for i in peaks])
    hts = np.array([_vertex_height(z, int(i)) for i in peaks])
    # a tooth is a RESOLVED peak: its prominence is a fair share of its
    # height. Noise on the body of a broad pattern rides high but shallow.
    resolved = props["prominences"] >= RESOLVED_FRAC * np.maximum(hts, 1e-300)
    resolved |= peaks == peaks[int(np.argmax(hts))]     # the top always is
    peaks, pos, hts = peaks[resolved], pos[resolved], hts[resolved]
    if peaks.size == 0:
        return None
    if centre_ppm is None:
        ci = int(np.argmax(hts))
    else:
        ci = int(np.argmin(np.abs(pos - float(centre_ppm))))
    centre, c_height = float(pos[ci]), float(hts[ci])
    widths = peak_widths(z, np.array([peaks[ci]]), rel_height=0.5)[0]
    c_fwhm = float(widths[0]) * dx if widths.size else dx
    tol = max(0.05 * spacing, 0.5 * c_fwhm, 3.0 * dx)
    tol = min(tol, 0.45 * spacing)          # one peak matches one order
    orders = []
    for sign in (1, -1):
        misses = 0
        for k in range(1, MAX_ORDER + 1):
            expected = centre + sign * k * spacing
            if expected < xu[0] or expected > xu[-1]:
                break
            j = int(np.argmin(np.abs(pos - expected)))
            if abs(pos[j] - expected) <= tol:
                orders.append(SidebandOrder(sign * k, float(expected),
                                            float(hts[j]), True))
                misses = 0
            else:
                orders.append(SidebandOrder(
                    sign * k, float(expected),
                    float(np.interp(expected, xu, z)), False))
                misses += 1
                if misses >= 2:
                    break
    orders.sort(key=lambda o: (abs(o.k), -o.k))
    return centre, c_fwhm, c_height, tuple(orders)


# ------------------------------------------------------------------ detect
def detect(x_ppm, y, larmor_MHz: float, nu_rot_Hz: float = 0.0, *,
           scan: bool = False, centre_ppm: float | None = None,
           tolerance: float = RATE_TOLERANCE,
           min_confidence: float = MIN_CONFIDENCE) -> SidebandDetection:
    """Find the +-nu_rot repeat in (x_ppm, y).

    ``nu_rot_Hz > 0``: the spacing is searched within ``tolerance`` of the
    recorded rate (``from_meta``). When that fails and ``scan`` is set, or
    when the rate is unknown (<= 0) and ``scan`` is set, SCAN_HZ is scanned
    for the fundamental repeat (``scanned``). An unknown rate without
    ``scan`` is a static experiment: nothing is detected. ``centre_ppm``
    anchors the comb on the peak nearest to it instead of the tallest one.

    Never raises: degenerate input (< 64 points, no Larmor frequency, a flat
    trace, a spacing beyond half the axis) comes back as ``ok=False`` with
    the reason in ``message``.
    """
    def fail(msg, **kw):
        return SidebandDetection(False, msg, **kw)

    try:
        x = np.asarray(x_ppm, float).ravel()
        yy = np.real(np.asarray(y)).astype(float).ravel()
    except (TypeError, ValueError):
        return fail("data is not numeric")
    n = min(x.size, yy.size)
    if n < 64:
        return fail(f"too few points ({n}) to look for a repeat")
    lar = float(larmor_MHz or 0.0)
    if lar <= 0:
        return fail("Larmor frequency unknown — set it in the experiment "
                    "parameters")
    nu = float(nu_rot_Hz or 0.0)
    if nu <= 0 and not scan:
        return fail("static (νrot = 0) — nothing to detect")

    xu, z, zc, dx, noise = _prepared(x, yy, MAX_POINTS)
    if xu.size < 64 or dx <= 0:
        return fail(f"too few finite points ({xu.size}) to look for a repeat")
    if not np.any(zc > 0) or np.ptp(z) <= 0:
        return fail("flat trace — nothing to detect")
    lags, c = _profile(zc, dx)
    half = float(lags[-1])

    idx = None
    from_meta = scanned = False
    if nu > 0:
        spacing0 = nu / lar
        if spacing0 > half:
            if not scan:
                return fail(f"νrot spacing ({spacing0:.3g} ppm) exceeds half "
                            f"the axis ({half:.3g} ppm)")
        else:
            i0 = spacing0 / dx
            lo = min(int(np.floor(spacing0 * (1 - tolerance) / dx)),
                     int(round(i0)) - 2)
            hi = max(int(np.ceil(spacing0 * (1 + tolerance) / dx)),
                     int(round(i0)) + 2)
            idx = _interior_max(c, lo, hi)
            if idx is not None and _prominence(lags, c, idx) <= 1e-3:
                idx = None                       # a numerical ripple, not a repeat
            if idx is not None:
                from_meta = True
            elif not scan:
                return fail(f"no repeat within ±{tolerance * 100:g} % of "
                            f"{format_hz(nu)} Hz")
    if idx is None:
        lo_ppm = SCAN_HZ[0] / lar
        hi_ppm = min(SCAN_HZ[1] / lar, half)
        if hi_ppm <= lo_ppm:
            return fail(f"axis too narrow to scan ({half * lar / 1000:.3g} kHz "
                        f"half-span, needs > {SCAN_HZ[0] / 1000:g} kHz)")
        found = _scan(lags, c, lo_ppm, hi_ppm)
        if found is None:
            return fail(f"no repeat between {SCAN_HZ[0] / 1000:g} and "
                        f"{hi_ppm * lar / 1000:.3g} kHz")
        idx = found[0]
        scanned = True

    spacing = _refine(c, idx) * dx
    nu_meas = spacing * lar
    score = float(c[idx])
    prom = _prominence(lags, c, idx)
    confidence = float(np.clip(prom / PROMINENCE_FULL, 0.0, 1.0))
    if scanned:
        confidence *= SCAN_PENALTY
    base = dict(nu_rot_Hz=float(nu_meas), spacing_ppm=float(spacing),
                score=score, prominence=float(prom), confidence=confidence,
                scanned=scanned, from_meta=from_meta)

    comb = _comb(xu, z, dx, spacing, centre_ppm, noise)
    if comb is None:
        return fail("no peak stands out of the noise", **base)
    centre, c_fwhm, c_height, orders = comb
    matched_k = [o.k for o in orders if o.matched]
    det = SidebandDetection(
        False, "", centre_ppm=float(centre), centre_fwhm_ppm=float(c_fwhm),
        centre_height=float(c_height), orders=orders,
        n_orders=max((abs(k) for k in matched_k), default=0), **base)
    plus, minus = 1 in matched_k, -1 in matched_k
    if not (plus and minus):
        where = ("either side" if not (plus or minus)
                 else ("the −1 side" if plus else "the +1 side"))
        det.message = (f"no symmetric ±1 pair about the tallest band at "
                       f"{centre:.3g} ppm (missing on {where})")
        return det
    if confidence < min_confidence:
        det.message = (f"repeat at {format_hz(nu_meas)} Hz too weak "
                       f"(confidence {confidence:.2f} < {min_confidence:g})")
        return det
    det.ok = True
    det.message = (f"spinning sidebands at ±k·{format_hz(nu_meas)} Hz, "
                   f"{det.n_orders} order(s)")
    return det


# ----------------------------------------------------------------- seeding
def seed_ratio(det: SidebandDetection) -> float:
    """Median height ratio between successive matched orders (the centre is
    order 0), clipped to [0.02, 0.98] -- the ``sidebands`` model's ssb_ratio
    seed. 0.3 (the model default) when nothing can be measured."""
    heights = {0: float(det.centre_height)}
    for o in det.orders:
        if o.matched:
            heights[o.k] = float(o.height)
    ratios = []
    for k, h in heights.items():
        if k == 0:
            continue
        inner = k - 1 if k > 0 else k + 1
        h0 = heights.get(inner)
        if h0 is not None and h0 > 0 and h > 0:
            ratios.append(h / h0)
    if not ratios:
        return 0.3
    return float(np.clip(np.median(ratios), 0.02, 0.98))


def describe(det: SidebandDetection, recipe_rate_Hz: float = 0.0) -> str:
    """The banner sentence: 'Spinning sidebands: νrot 20 000 Hz (124.6 ppm at
    160.46 MHz) · 2 orders each side · confidence 0.90', plus
    ' · ⚠ acquisition says 20 300 Hz' when the recorded rate differs by more
    than RATE_MISMATCH."""
    lar = det.nu_rot_Hz / det.spacing_ppm if det.spacing_ppm else 0.0
    n_plus = max((o.k for o in det.orders if o.matched and o.k > 0), default=0)
    n_minus = max((-o.k for o in det.orders if o.matched and o.k < 0), default=0)
    if n_plus == n_minus:
        orders = f"{n_plus} order{'s' if n_plus != 1 else ''} each side"
    else:
        orders = f"orders +{n_plus} / −{n_minus}"
    text = (f"Spinning sidebands: νrot {format_hz(det.nu_rot_Hz)} Hz "
            f"({det.spacing_ppm:.1f} ppm at {lar:.2f} MHz) · {orders} · "
            f"confidence {det.confidence:.2f}")
    rec = float(recipe_rate_Hz or 0.0)
    if rec > 0 and abs(det.nu_rot_Hz - rec) > RATE_MISMATCH * rec:
        text += f" · ⚠ acquisition says {format_hz(rec)} Hz"
    return text
