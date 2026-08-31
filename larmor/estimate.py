"""Starting values MEASURED from the spectrum, not guessed from a table.

A new site used to start from the model's static defaults (Cq 5 MHz, width
5 ppm), optionally nudged by a per-nucleus table. For a nucleus not in that
table -- 81Br, 127I, 209Bi, anything heavy -- the result is a needle on a
several-thousand-ppm axis: the fit has no gradient to follow and does not
move. Since LARMOR is holding the data anyway, measure the starting values
from it.

Qt-free and testable.
"""
from __future__ import annotations

import numpy as np

#: the parameter that carries the BREADTH of each model, and whether that
#: breadth is a quadrupolar coupling (which has to be solved for) or simply a
#: line width in ppm (which can be set directly)
#: height fraction used when matching a MODEL's breadth to the data's. Low on
#: purpose: see band_width_ppm.
CALIB_FRAC = 0.10

_WIDTH_KEY = {
    "czjzek": ("sigma_Cq_MHz", True),
    "ext_czjzek": ("Cq_MHz", True),
    "quad_ct": ("Cq_MHz", True),
    "quad_csa": ("Cq_MHz", True),
    "quad_first": ("Cq_MHz", True),
    "amorphous": ("Cq_MHz", True),
    "gauss_lor": ("shift_fwhm_ppm", False),
    "gl_norm": ("shift_fwhm_ppm", False),
    "jmultiplet": ("shift_fwhm_ppm", False),
    "sidebands": ("shift_fwhm_ppm", False),
    "csa_mas": ("shift_fwhm_ppm", False),
    "csa_czjzek": ("shift_fwhm_ppm", False),
    "voigt": ("gauss_fwhm_ppm", False),
}


def band_width_ppm(ppm: np.ndarray, amp: np.ndarray,
                   centre_ppm: float | None = None, frac: float = 0.5
                   ) -> tuple[float, float]:
    """(centre, width) in ppm of the band containing ``centre_ppm`` -- or of
    the tallest feature when no centre is given -- measured at ``frac`` of
    the peak height above the baseline (edges of the axis).

    ``frac=0.5`` is the FWHM a user would read off. For CALIBRATING a
    quadrupolar model use a lower fraction: a second-order powder pattern is
    a pair of sharp horns, so its half-height width is the width of ONE horn
    and says almost nothing about the breadth that Cq controls."""
    x = np.asarray(ppm, float)
    y = np.asarray(amp, float)
    if x.size < 3:
        return (float(centre_ppm or 0.0), 0.0)
    order = np.argsort(x)
    x, y = x[order], y[order]
    k = max(1, x.size // 10)
    y = y - float(np.median(np.concatenate([y[:k], y[-k:]])))
    if centre_ppm is None:
        i = int(np.argmax(y))
    else:
        i = int(np.argmin(np.abs(x - float(centre_ppm))))
        while 0 < i < y.size - 1 and (y[i + 1] > y[i] or y[i - 1] > y[i]):
            i = i + 1 if y[i + 1] >= y[i - 1] else i - 1
    peak = float(y[i])
    if peak <= 0:
        return (float(x[i]), 0.0)
    half = float(np.clip(frac, 0.01, 0.99)) * peak
    lo, hi = i, i
    while lo > 0 and y[lo] > half:
        lo -= 1
    while hi < y.size - 1 and y[hi] > half:
        hi += 1
    return (float(x[i]), abs(float(x[hi] - x[lo])))


def _model_width(model: str, key: str, val: float, nucleus: str,
                 larmor_MHz: float, eta: float, span_ppm: float) -> float:
    """FWHM (ppm) of ``model`` when ``key`` is set to ``val``. 0.0 when the
    pattern does not fit inside ``span_ppm`` (so a caller can widen)."""
    from larmor import engine
    from larmor.models.base import REGISTRY
    from larmor.recipe import Param, Recipe, SiteModel

    m = REGISTRY.get(model)
    if m is None:
        return 0.0
    params = {p.name: Param(p.default) for p in m.params}
    params["isotropic_chemical_shift_ppm"] = Param(0.0)
    params["amplitude"] = Param(1.0)
    if "eta" in params:
        params["eta"] = Param(float(eta))
    if "shift_fwhm_ppm" in params:
        params["shift_fwhm_ppm"] = Param(max(span_ppm * 0.002, 0.5))
    params[key] = Param(float(val))
    rec = Recipe(nucleus=nucleus, larmor_frequency_MHz=larmor_MHz,
                 spin_rate_Hz=0.0,
                 sites=[SiteModel(model=model, label="r", params=params)])
    x = np.linspace(-span_ppm, span_ppm, 3001)
    try:
        # simulate() returns the axis it CHOSE (the Czjzek kernel axis for a
        # kernel model, not the one passed in) -- measure on that one
        gx, y, _ = engine.simulate(rec, exp_ppm=x)
    except Exception:                                     # noqa: BLE001
        return 0.0
    if not np.any(y > 0):
        return 0.0
    # a pattern touching the edge of its own axis is not measurable here
    if max(float(y[0]), float(y[-1])) > 0.05 * float(y.max()):
        return 0.0
    return band_width_ppm(gx, y, frac=CALIB_FRAC)[1]


def cq_for_width(fwhm_ppm: float, nucleus: str, larmor_MHz: float,
                 eta: float = 0.6, *, model: str = "quad_ct",
                 key: str = "Cq_MHz") -> float:
    """The value of ``key`` that makes ``model`` a pattern ``fwhm_ppm`` wide.

    Solved by bisection on the model itself rather than an analytic formula:
    the breadth is monotone in the coupling for every model here, but a
    Czjzek sigma and a discrete Cq are not related by any fixed factor, and
    the convolved shift distribution rides along with both.
    """
    if fwhm_ppm <= 0 or larmor_MHz <= 0:
        return 0.0
    span = max(8.0 * fwhm_ppm, 500.0)
    lo, hi = 0.05, 200.0

    def w(v):
        return _model_width(model, key, v, nucleus, larmor_MHz, eta, span)

    w_lo, w_hi = w(lo), w(hi)
    if w_lo <= 0:
        return 0.0
    if w_lo >= fwhm_ppm:                       # already too broad at the floor
        return lo
    while w_hi <= 0 or w_hi < fwhm_ppm:        # shrink the ceiling until usable
        hi *= 0.6
        if hi <= lo * 1.01:
            return 0.0
        w_hi = w(hi)
    for _ in range(10):
        mid = 0.5 * (lo + hi)
        wm = w(mid)
        if wm <= 0:                            # outgrew the window: too big
            hi = mid
            continue
        if abs(wm - fwhm_ppm) <= 0.05 * fwhm_ppm:
            return float(mid)
        if wm < fwhm_ppm:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi))


def start_values(model: str, ppm, amp, nucleus: str, larmor_MHz: float,
                 centre_ppm: float | None = None) -> dict:
    """Starting parameter values for a new ``model`` site, measured from the
    spectrum. Returns only the keys the data justifies; the caller keeps the
    model's own defaults for everything else."""
    centre, fwhm = band_width_ppm(ppm, amp, centre_ppm)
    if fwhm <= 0:
        return {}
    entry = _WIDTH_KEY.get(model)
    if entry is None:
        return {}
    key, is_cq = entry
    if not is_cq:
        out = {key: float(fwhm)}
        if model == "voigt":                   # split the width between both
            out["gauss_fwhm_ppm"] = float(fwhm * 0.7)
            out["lorentz_fwhm_ppm"] = float(fwhm * 0.4)
        return out
    _, breadth = band_width_ppm(ppm, amp, centre_ppm, frac=CALIB_FRAC)
    # probe with quad_ct always: it is the plain second-order CT pattern all
    # of these models are built on, it needs no kernel (0.1 s rather than
    # several seconds), and it is stable. The seed only has to land in the
    # right ballpark -- the fit refines it from there.
    _, breadth = band_width_ppm(ppm, amp, centre_ppm, frac=CALIB_FRAC)
    cq = cq_for_width(max(breadth, fwhm), nucleus, larmor_MHz,
                      model="quad_ct", key="Cq_MHz")
    if cq <= 0:
        return {"shift_fwhm_ppm": float(fwhm)}
    if key == "sigma_Cq_MHz":
        # a Czjzek sigma spreads Cq over a distribution, so the same breadth
        # comes from a different number -- and not by a fixed factor, since
        # the distribution saturates at large sigma. One measured rescale of
        # the real model gets close enough to start from.
        # measured on real cases: within about a factor 1.5 in sigma, which
        # is a perfectly good STARTING point -- refining it here cost seconds
        # per click and did not converge more reliably than the fit itself.
        cq *= 0.5
    return {key: float(cq), "shift_fwhm_ppm": float(max(fwhm * 0.10, 1.0))}
