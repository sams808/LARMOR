"""Baseline correction through the pybaselines library (Qt-free core).

A curated registry of pybaselines algorithms chosen for solid-state NMR
spectra -- the Whittaker smoothers (arPLS, asLS, airPLS, IarPLS, drPLS,
asPLS), the iterative polynomials (ModPoly, IModPoly, penalized poly) and
the window methods (SNIP, rolling ball, Mor, MorMol) -- with, for each one,
the parameters worth exposing, their ranges and defaults that follow the
spectrum length. ``compute`` is the single entry point: the processing op
``pybaseline`` (``larmor.processing.op_pybaseline``) and the desktop dialog
both go through it, so what the preview shows is what the recipe records
and what a saved fit replays.

Reference (cite it when you publish a corrected spectrum):
    D. Erb, "pybaselines: A Python library of algorithms for the baseline
    correction of experimental data", Zenodo, doi:10.5281/zenodo.5608581;
    documentation: https://pybaselines.readthedocs.io

The parameter names below are pybaselines 1.2's own (read from the
installed package's docstrings), so a recorded step is also a valid direct
call of ``pybaselines.Baseline``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "CITATION", "DOCS_URL", "FAMILIES", "LAM_MAX", "LAM_MIN", "LAM_REF_POINTS",
    "METHODS", "MethodSpec", "ParamSpec", "available", "compute",
    "decimate", "default_params", "preview_params", "resolved_max",
    "suggest_half_window", "suggest_lam", "validate_params",
]

CITATION = ('D. Erb, "pybaselines: A Python library of algorithms for the '
            'baseline correction of experimental data", Zenodo, '
            "doi:10.5281/zenodo.5608581")
DOCS_URL = "https://pybaselines.readthedocs.io"
_INSTALL_HINT = ("pybaselines is not installed -- install it with "
                 "`pip install pybaselines` in LARMOR's environment "
                 "(`conda activate larmor` first)")

#: the decade range of the Whittaker smoothing parameter
LAM_MIN, LAM_MAX = 1e2, 1e12
#: the spectrum length the registry's λ defaults are written for
LAM_REF_POINTS = 2048

#: family key -> the heading the dialog groups the methods under
FAMILIES = {
    "whittaker": "Whittaker smoothers (penalized least squares)",
    "polynomial": "Iterative polynomials",
    "window": "Peak clipping and morphology",
}


@dataclass(frozen=True)
class ParamSpec:
    """One exposed parameter of a method.

    ``kind`` is ``float`` / ``int`` / ``bool`` / ``choice``; ``log`` marks a
    parameter that lives on a log scale (the dialog adds a slider in
    decades); ``scale`` says how the default follows the spectrum length and
    how a decimated preview rescales the value: ``"lam"`` (Whittaker λ, see
    ``suggest_lam``) or ``"window"`` (a half-window in points). ``max`` None
    on a window means half the spectrum length (``resolved_max``)."""
    key: str
    label: str
    kind: str
    default: object
    min: float | None = None
    max: float | None = None
    step: float | None = None
    log: bool = False
    tooltip: str = ""
    scale: str = ""
    choices: tuple = ()


@dataclass(frozen=True)
class MethodSpec:
    key: str            # also the pybaselines.Baseline method name
    label: str
    family: str
    description: str    # one line: which baseline shape it suits
    params: tuple[ParamSpec, ...]

    def param(self, key: str) -> ParamSpec | None:
        for p in self.params:
            if p.key == key:
                return p
        return None


# ------------------------------------------------------------ parameters
def _lam(default: float) -> ParamSpec:
    return ParamSpec(
        "lam", "λ (smoothness)", "float", default, LAM_MIN, LAM_MAX, None,
        log=True, scale="lam",
        tooltip="Whittaker smoothing parameter: larger = stiffer, smoother "
                "baseline. Decades matter, not digits. The default follows the "
                "spectrum length (×16 per doubling of the points at difference "
                "order 2).")


def _diff(min_order: int = 1) -> ParamSpec:
    return ParamSpec(
        "diff_order", "difference order", "int", 2, min_order, 3, 1,
        tooltip="Order of the finite-difference penalty: 2 penalises "
                "curvature (the usual choice), 1 slope, 3 the change of "
                "curvature.")


def _max_iter(default: int) -> ParamSpec:
    return ParamSpec(
        "max_iter", "max iterations", "int", default, 1, 1000, 10,
        tooltip="Upper bound on the reweighting iterations; the method stops "
                "earlier once the weights settle.")


def _window(key: str, label: str, tooltip: str) -> ParamSpec:
    return ParamSpec(key, label, "int", 64, 1, None, 1, scale="window",
                     tooltip=tooltip)


_P = ParamSpec(
    "p", "asymmetry p", "float", 0.01, 1e-4, 0.5, None, log=True,
    tooltip="Weight of the points above the baseline (1 − p for those below). "
            "0.001–0.05 keeps the baseline under positive peaks.")
_ETA = ParamSpec(
    "eta", "η", "float", 0.5, 0.0, 1.0, 0.05,
    tooltip="drPLS balance: 0 gives the smoothest baseline, higher values let "
            "it follow curvature under the peaks more closely.")
_ASYM = ParamSpec(
    "asymmetric_coef", "asymmetry coefficient", "float", 0.5, 0.05, 10.0, 0.05,
    tooltip="asPLS weighting steepness: larger = a more step-like weight, "
            "sharper separation of peak and baseline points.")
_NORMW = ParamSpec(
    "normalize_weights", "normalise weights", "bool", False,
    tooltip="Rescale the weights to 0–1 at each pass (steadier on very strong "
            "peaks); off reproduces the original airPLS.")
_POLY = ParamSpec(
    "poly_order", "polynomial order", "int", 3, 1, 12, 1,
    tooltip="Order of the baseline polynomial: 2–4 for a gently curved "
            "baseline; higher orders start to follow the peaks.")
_MASK = ParamSpec(
    "mask_initial_peaks", "mask initial peaks", "bool", False,
    tooltip="Exclude from the first fit every point above the initial "
            "baseline plus one residual standard deviation.")
_NUMSTD = ParamSpec(
    "num_std", "noise threshold (σ)", "float", 1.0, 0.0, 5.0, 0.1,
    tooltip="How many standard deviations of the residual count as noise "
            "rather than peak when the polynomial is re-fitted.")
_COST = ParamSpec(
    "cost_function", "cost function", "choice", "asymmetric_truncated_quadratic",
    choices=("asymmetric_truncated_quadratic", "symmetric_truncated_quadratic",
             "asymmetric_huber", "symmetric_huber",
             "asymmetric_indec", "symmetric_indec"),
    tooltip="Robust loss of the polynomial fit: 'asymmetric' variants ignore "
            "positive residuals (peaks pointing up), 'symmetric' ones treat "
            "both signs alike.")
_SMOOTH = ParamSpec(
    "smooth_half_window", "smoothing half-window", "int", 0, 0, None, 1,
    scale="window",
    tooltip="Moving-average half-window in points; 0 = the method's default "
            "(no smoothing for SNIP and MorMol, the half-window itself for "
            "the rolling ball).")
_DECR = ParamSpec(
    "decreasing", "decreasing windows", "bool", False,
    tooltip="Clip from the largest window down to 1 instead of up from 1 -- "
            "a smoother baseline.")
_FILTER = ParamSpec(
    "filter_order", "filter order", "choice", 2, choices=(2, 4, 6, 8),
    tooltip="2 approximates a linear baseline between the clipping points; "
            "4–8 follow a more complicated curved baseline.")
_HALF = _window(
    "half_window", "half-window (points)",
    "Half-width in points of the structuring window; make it wider than the "
    "widest peak's half-width or the baseline climbs into the peaks. The "
    "default is 1/32 of the spectrum.")
_SNIP_HALF = _window(
    "max_half_window", "max half-window (points)",
    "Largest clipping half-window in points, about half the width of the "
    "widest peak; the default is 1/32 of the spectrum.")

# ---------------------------------------------------------------- methods
_M = (
    MethodSpec("arpls", "arPLS", "whittaker",
               "A smooth, slowly curving baseline under peaks of any width; "
               "the weights adapt to the noise, so only λ is left to set. The "
               "general first choice.",
               (_lam(1e5), _diff(), _max_iter(50))),
    MethodSpec("asls", "asLS", "whittaker",
               "The classic Whittaker smoother with a fixed asymmetry p; a "
               "smooth baseline when every peak points the same way.",
               (_lam(1e6), _P, _diff(), _max_iter(50))),
    MethodSpec("airpls", "airPLS", "whittaker",
               "Adaptive reweighting that discounts the peak regions more "
               "strongly at each pass; a drifting baseline under dense, "
               "overlapping peaks.",
               (_lam(1e6), _diff(), _max_iter(50), _NORMW)),
    MethodSpec("iarpls", "IarPLS", "whittaker",
               "arPLS with a gentler weighting for weak signals; keeps low, "
               "broad features from being taken as baseline.",
               (_lam(1e5), _diff(), _max_iter(50))),
    MethodSpec("drpls", "drPLS", "whittaker",
               "Doubly reweighted: a smooth baseline that still follows a "
               "strong curvature; η balances smoothness against fit.",
               (_lam(1e5), _ETA, _diff(2), _max_iter(50))),
    MethodSpec("aspls", "asPLS", "whittaker",
               "Adaptive smoothness: λ varies along the axis, so a baseline "
               "flat in one region and curved in another; the slowest of the "
               "family.",
               (_lam(1e5), _ASYM, _diff(), _max_iter(100))),
    MethodSpec("modpoly", "ModPoly", "polynomial",
               "Iterative polynomial (Lieber & Mahadevan 2003): a gently "
               "curved baseline that a low-order polynomial describes.",
               (_POLY, _max_iter(250), _MASK)),
    MethodSpec("imodpoly", "IModPoly", "polynomial",
               "ModPoly with a noise-aware threshold (Zhao et al. 2007): the "
               "same polynomial baseline on a noisy spectrum.",
               (_POLY, _NUMSTD, _max_iter(250))),
    MethodSpec("penalized_poly", "Penalized poly", "polynomial",
               "Polynomial fit under a robust, non-quadratic cost (Mazet et "
               "al. 2005): the peaks count for little without iterating on a "
               "mask.",
               (_POLY, _COST, _max_iter(250))),
    MethodSpec("snip", "SNIP", "window",
               "Statistics-sensitive peak clipping (Ryan et al. 1988): a "
               "baseline under narrow lines; set the half-window wider than "
               "the widest peak.",
               (_SNIP_HALF, _DECR, _SMOOTH, _FILTER)),
    MethodSpec("rolling_ball", "Rolling ball", "window",
               "A ball of the given radius rolled under the spectrum (Kneen & "
               "Annegarn 1996): a wide rolling baseline under narrow lines.",
               (_HALF, _SMOOTH)),
    MethodSpec("mor", "Mor", "window",
               "Morphological opening: a fast baseline under narrow lines; "
               "the half-window must exceed the widest peak's half-width.",
               (_HALF,)),
    MethodSpec("mormol", "MorMol", "window",
               "Morphological opening, smoothed and iterated (Koch et al. "
               "2017): the baseline of Mor without its flat steps; slower.",
               (_HALF, _SMOOTH, _max_iter(250))),
)
#: ordered registry: key -> MethodSpec
METHODS: dict[str, MethodSpec] = {m.key: m for m in _M}


# ---------------------------------------------------------------- rules
def suggest_lam(n_points: int, diff_order: int = 2, lam_ref: float = 1e5) -> float:
    """The Whittaker λ that keeps the baseline stiffness of ``lam_ref`` on a
    ``LAM_REF_POINTS``-point spectrum when the spectrum has ``n_points``.

    The penalty λ·Σ(Δᵈz)² is written in points: sampling the same curve ``f``
    times more densely makes every d-th difference ``f``ᵈ times smaller, so
    the same curve costs λ·f^(2d) -- λ must grow as (n / n_ref)^(2·diff_order)
    (×16 per doubling at order 2). Rounded to a decade, clipped to
    [LAM_MIN, LAM_MAX].
    """
    f = max(int(n_points), 2) / LAM_REF_POINTS
    lam = float(lam_ref) * f ** (2 * max(int(diff_order), 1))
    lam = 10.0 ** round(math.log10(lam))
    return float(min(max(lam, LAM_MIN), LAM_MAX))


def suggest_half_window(n_points: int) -> int:
    """Default half-window of the window methods: 1/32 of the spectrum (at
    least 2 points) -- wider than a MAS line, narrower than a static pattern;
    raise it when the baseline climbs into a broad peak."""
    return max(2, int(n_points) // 32)


def resolved_max(p: ParamSpec, n_points: int) -> float | None:
    """The upper bound of ``p`` on an ``n_points`` spectrum: a window's None
    means half the spectrum."""
    if p.max is None and p.scale == "window":
        return max(2, int(n_points) // 2)
    return p.max


def default_params(method: str, n_points: int | None = None) -> dict:
    """Every parameter of ``method`` at its default, the λ and window
    defaults scaled to ``n_points`` when given (see ``suggest_lam`` /
    ``suggest_half_window``)."""
    spec = _spec(method)
    out: dict = {}
    for p in spec.params:
        v = p.default
        if n_points:
            if p.scale == "lam":
                d = spec.param("diff_order")
                v = suggest_lam(n_points, int(d.default) if d else 2, float(p.default))
            elif p.scale == "window" and p.key != "smooth_half_window":
                v = min(suggest_half_window(n_points), resolved_max(p, n_points))
        out[p.key] = v
    return out


def preview_params(method: str, params: dict, factor: int) -> dict:
    """``params`` rescaled for a spectrum decimated by ``factor`` (every
    factor-th point), so a preview on fewer points shows the baseline the
    full-length call will give: λ ÷ factor^(2·diff_order) (the inverse of
    ``suggest_lam``'s rule), windows ÷ factor (0 = 'default' stays 0)."""
    spec = _spec(method)
    factor = max(int(factor), 1)
    out = dict(params)
    if factor == 1:
        return out
    d = spec.param("diff_order")
    diff_order = int(out.get("diff_order", d.default if d else 2))
    for p in spec.params:
        if p.key not in out:
            continue
        if p.scale == "lam":
            out[p.key] = max(float(out[p.key]) / factor ** (2 * diff_order), 1e-6)
        elif p.scale == "window" and int(out[p.key]) > 0:
            out[p.key] = max(1, int(round(int(out[p.key]) / factor)))
    return out


def decimate(x, y, max_points: int = 8192):
    """``(x, y, factor)`` with every ``factor``-th point kept so that at most
    ``max_points`` remain (factor 1 = unchanged). A stride, not a block mean:
    the noise statistics the threshold methods read stay those of the data."""
    y = np.asarray(y)
    n = y.size
    factor = max(1, int(math.ceil(n / max(int(max_points), 8))))
    if factor == 1:
        return x, y, 1
    xd = None if x is None else np.asarray(x)[::factor]
    return xd, y[::factor], factor


# -------------------------------------------------------------- compute
def available() -> bool:
    """Is pybaselines importable?"""
    try:
        _baseline_class()
    except ImportError:
        return False
    return True


def _baseline_class():
    try:
        from pybaselines import Baseline
    except ImportError as exc:                      # pragma: no cover - env
        raise ImportError(_INSTALL_HINT) from exc
    return Baseline


def _spec(method: str) -> MethodSpec:
    spec = METHODS.get(str(method))
    if spec is None:
        raise ValueError(f"unknown pybaselines method {method!r} "
                         f"(valid: {', '.join(METHODS)})")
    return spec


def _cast(spec: MethodSpec, p: ParamSpec, v, n_points: int):
    """``v`` as the Python value ``p`` takes, range-checked; a readable
    ValueError names the method and the parameter."""
    who = f"{spec.label}: {p.label}"
    if p.kind == "bool":
        if isinstance(v, (bool, np.bool_)) or v in (0, 1):
            return bool(v)
        raise ValueError(f"{who} must be true or false (got {v!r})")
    if p.kind == "choice":
        ctype = type(p.choices[0])
        try:
            cv = ctype(v)
        except (TypeError, ValueError):
            cv = v
        if cv not in p.choices:
            raise ValueError(f"{who} must be one of "
                             f"{', '.join(str(c) for c in p.choices)} (got {v!r})")
        return cv
    try:
        fv = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"{who} must be a number (got {v!r})") from None
    if not math.isfinite(fv):
        raise ValueError(f"{who} must be finite (got {v!r})")
    if p.kind == "int":
        if abs(fv - round(fv)) > 1e-9:
            raise ValueError(f"{who} must be a whole number (got {v!r})")
        fv = int(round(fv))
    lo, hi = p.min, resolved_max(p, n_points)
    if (lo is not None and fv < lo) or (hi is not None and fv > hi):
        rng = (f"between {lo:g} and {hi:g}" if lo is not None and hi is not None
               else f"at least {lo:g}" if lo is not None else f"at most {hi:g}")
        raise ValueError(f"{who} must be {rng} (got {fv:g})")
    return int(fv) if p.kind == "int" else float(fv)


def validate_params(method: str, params: dict, n_points: int) -> dict:
    """Complete, typed and range-checked parameters of ``method`` for an
    ``n_points`` spectrum: missing ones take their (length-scaled) default,
    an unknown name or an out-of-range value raises a readable ValueError."""
    spec = _spec(method)
    known = {p.key for p in spec.params}
    extra = sorted(set(params) - known)
    if extra:
        raise ValueError(f"{spec.label} has no parameter "
                         f"{', '.join(repr(k) for k in extra)} "
                         f"(valid: {', '.join(sorted(known))})")
    defaults = default_params(method, n_points)
    out = {}
    for p in spec.params:
        v = params.get(p.key, defaults[p.key])
        out[p.key] = _cast(spec, p, v, n_points)
    return out


def _call_kwargs(params: dict) -> dict:
    """The validated parameters as pybaselines takes them: a 0 smoothing
    half-window means 'the method's default' (None)."""
    kw = dict(params)
    if kw.get("smooth_half_window") == 0:
        kw["smooth_half_window"] = None
    return kw


def compute(x, y, method: str = "arpls", **params) -> np.ndarray:
    """The baseline of ``y`` (real, any axis direction) by ``method``.

    Returns a float array of ``y``'s length; NaN / inf points are bridged by
    linear interpolation before the call, so subtracting the result keeps
    them NaN and nothing else. Raises ValueError with a readable message on
    an unknown method or a bad parameter, ImportError with the install hint
    when pybaselines is missing.
    """
    spec = _spec(method)
    y = np.asarray(y, float).ravel()
    n = y.size
    if n < 8:
        raise ValueError("baseline correction needs at least 8 points")
    if x is None:
        x = np.arange(n, dtype=float)
    else:
        x = np.asarray(x, float).ravel()
        if x.size != n:
            raise ValueError(f"x has {x.size} points, y has {n}")
        if not np.isfinite(x).all() or np.unique(x).size != n:
            x = np.arange(n, dtype=float)      # a broken axis: use the index
    kw = _call_kwargs(validate_params(method, params, n))
    finite = np.isfinite(y)
    if finite.sum() < 8:
        raise ValueError("the spectrum has fewer than 8 finite points")
    if not finite.all():
        idx = np.arange(n, dtype=float)
        y = y.copy()
        y[~finite] = np.interp(idx[~finite], idx[finite], y[finite])
    Baseline = _baseline_class()
    fn = getattr(Baseline(x_data=x), spec.key)
    try:
        base, _info = fn(y, **kw)
    except ValueError as exc:
        raise ValueError(f"{spec.label}: {exc}") from exc
    base = np.asarray(base, float).ravel()
    if base.shape != y.shape:
        raise ValueError(f"{spec.label} returned {base.size} points for {n}")
    if not np.isfinite(base).all():
        raise ValueError(f"{spec.label} did not reach a finite baseline with "
                         "these parameters -- change λ or the window")
    return base
