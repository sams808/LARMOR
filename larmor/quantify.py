"""Quantification: per-site integrals and relative populations.

This is dmfit's results table -- the thing that actually gets published for a
glass: each site's integrated intensity as a percentage of the total, with an
uncertainty.

Error model: for a fixed lineshape, a site's integral is proportional to its
amplitude parameter, so the fractional error on the integral is taken from the
fractional error on the amplitude (first-order approximation; shape-parameter
covariance is neglected and this is stated in the output).
"""
from __future__ import annotations

import numpy as np

from larmor import families
from larmor.engine import make_context, simulate_site
from larmor.recipe import Recipe

#: share (%) of a line's simulated area outside the integration window above
#: which the fit-health strip flags the line (its population is biased low)
TAIL_LIMIT_PCT = 2.0
#: share of every line's area the widened window keeps (0.25 % per side), so
#: the widened window leaves at most 0.5 % outside and the chip clears
TAIL_COVERAGE = 0.995
#: the two data-backed background models engine.simulate_site special-cases;
#: a background has no tail notion
_NO_TAIL_MODELS = frozenset({"spectrum", "function"})


def _cumulative_area(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(x ascending, cumulative trapezoid of |y| from the low end)."""
    order = np.argsort(x)
    xs = np.asarray(x, float)[order]
    a = np.abs(np.asarray(y, float))[order]
    if a.size < 2:
        return xs, np.zeros_like(xs)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (a[1:] + a[:-1]) * np.diff(xs))])
    return xs, cum


def _tail_outside_pct(y: np.ndarray, x: np.ndarray, lo: float, hi: float):
    """Share (%) of |y|'s area outside [lo, hi], with the window edges
    interpolated on the cumulative area (exact to the grid, not to the
    nearest grid point); None when the line has no area."""
    xs, cum = _cumulative_area(y, x)
    full = float(cum[-1]) if cum.size else 0.0
    if full <= 0:
        return None
    inside = (float(np.interp(min(hi, xs[-1]), xs, cum))
              - float(np.interp(max(lo, xs[0]), xs, cum)))
    return float(np.clip(100.0 * (1.0 - inside / full), 0.0, 100.0))


def _integrate(recipe: Recipe, window_ppm):
    """Simulate every site and integrate it over the window: ``(ctx, hi, lo,
    [y per site], integrals)`` -- the one loop quantify() and
    site_integrals() share, so the Report and the Monte-Carlo worker
    integrate identically."""
    ctx = make_context(recipe)
    window = window_ppm or recipe.fit_window_ppm or \
        (float(ctx.x_ppm.max()), float(ctx.x_ppm.min()))
    hi, lo = max(window), min(window)
    sel = (ctx.x_ppm >= lo) & (ctx.x_ppm <= hi)
    ys = [simulate_site(site, ctx) for site in recipe.sites]
    integrals = np.asarray([float(np.trapezoid(y[sel], ctx.x_ppm[sel]))
                            for y in ys], float)
    return ctx, hi, lo, ys, integrals


def site_integrals(recipe: Recipe, window_ppm: tuple[float, float] | None = None
                   ) -> tuple[np.ndarray, tuple[float, float]]:
    """Per-site window integrals (signed, same trapezoid as quantify()) and
    the ``(hi, lo)`` window actually integrated -- what a Monte-Carlo trial
    records so family sums can be re-integrated per trial (larmor.families).
    """
    _ctx, hi, lo, _ys, integrals = _integrate(recipe, window_ppm)
    return integrals, (hi, lo)


def _same_window(a, b) -> bool:
    """Two (hi, lo) windows agree to 1e-6 of their span (float view ranges)."""
    if a is None or b is None or len(a) != 2 or len(b) != 2:
        return False
    tol = 1e-6 * max(1.0, abs(float(b[0]) - float(b[1])))
    return (abs(max(a) - max(b)) <= tol) and (abs(min(a) - min(b)) <= tol)


def quantify(recipe: Recipe, window_ppm: tuple[float, float] | None = None,
             *, uvars: dict | None = None, mc=None) -> dict:
    """Integrate every site over the window. Returns a JSON-friendly table.

    Each row also carries ``tail_outside_pct``: the share of the line's
    simulated |area| that lies outside the window, measured on the
    simulation axis (``axis_ppm``) — a lower bound when the line extends
    beyond that axis; None for the background models.

    Sites carrying a ``family`` tag are additionally summed into families
    with the named ratios of the nucleus (``families`` / ``ratios`` /
    ``untagged`` / ``family_basis`` / ``family_note``, see larmor.families).
    ``uvars`` (``fit.amplitude_uvars`` of the fit that produced these values)
    selects the covariance basis; ``mc`` (a ``MonteCarloResult`` carrying
    ``site_integrals`` over THIS window) the Monte-Carlo basis; otherwise
    the family block is flagged ``independent``. The per-site rows are
    identical whatever the basis. An untagged recipe adds empty lists.
    """
    ctx, hi, lo, ys, integrals = _integrate(recipe, window_ppm)

    rows = []
    for i, site in enumerate(recipe.sites):
        y = ys[i]
        integral = float(integrals[i])
        tail = (None if site.model in _NO_TAIL_MODELS
                else _tail_outside_pct(y, ctx.x_ppm, lo, hi))
        amp = site.params["amplitude"]
        # an ill-conditioned covariance (e.g. an amplitude at/near a bound, or
        # degenerate with another free parameter) makes lmfit report stderr as
        # NaN rather than None -- and NaN is truthy, so a plain "if amp.stderr"
        # would let it through as if it were a real error. Require it finite.
        amp_err = amp.stderr if (amp.stderr is not None
                                 and np.isfinite(amp.stderr)) else None
        rel_err = (amp_err / amp.value) if (amp_err and amp.value) else None
        pos = site.params["isotropic_chemical_shift_ppm"]
        pos_err = (pos.stderr if (pos.stderr is not None
                                  and np.isfinite(pos.stderr)) else None)
        rows.append({
            "site": f"s{i}",
            "label": site.label or site.model,
            "model": site.model,
            "position_ppm": pos.value,
            "position_err": pos_err,
            "integral": integral,
            "integral_err": abs(integral) * rel_err if rel_err is not None else None,
            "tail_outside_pct": tail,
        })

    total = sum(abs(r["integral"]) for r in rows) or 1.0
    for r in rows:
        r["fraction_pct"] = 100.0 * abs(r["integral"]) / total
        r["fraction_err_pct"] = (
            100.0 * r["integral_err"] / total
            if r["integral_err"] is not None else None)

    # ---- family sums and named ratios (only when a line is tagged)
    samples = None
    fam_note_extra = ""
    if mc is not None:
        mc_ints = getattr(mc, "site_integrals", None)
        mc_win = getattr(mc, "window_ppm", None)
        if mc_ints is not None and _same_window(mc_win, (hi, lo)):
            samples = mc_ints
        elif mc_ints is not None and mc_win is not None:
            fam_note_extra = (f"Monte-Carlo integrals not used: they cover "
                              f"{max(mc_win):g}..{min(mc_win):g} ppm, the "
                              f"table {hi:g}..{lo:g} ppm")
    fam = families.summarize(
        integrals, [s.params["amplitude"].value for s in recipe.sites],
        [getattr(s, "family", "") for s in recipe.sites], recipe.nucleus,
        sigma=[r["integral_err"] for r in rows], uvars=uvars, samples=samples)
    fam_note = fam["note"]
    if fam_note_extra and fam["families"]:
        fam_note = fam_note_extra + ("; " + fam_note if fam_note else "")

    return {
        "window_ppm": [hi, lo],
        "axis_ppm": [float(ctx.x_ppm.max()), float(ctx.x_ppm.min())],
        "rows": rows,
        "note": "fraction errors are first-order (amplitude covariance only; "
                "lineshape-parameter covariance neglected); tail_outside_pct "
                "is a lower bound measured on the simulation axis",
        "families": fam["families"],
        "ratios": fam["ratios"],
        "untagged": fam["untagged"],
        "family_basis": fam["basis"],
        "family_note": fam_note,
    }


def window_containing_tails(recipe: Recipe, coverage: float = TAIL_COVERAGE,
                            exp_ppm=None) -> tuple[float, float, dict]:
    """The narrowest (hi, lo) ppm window holding ``coverage`` of every
    non-background line's |area| (the two tails split the remainder equally),
    from the cumulative area of each line on the simulation axis. With
    ``exp_ppm`` the window is clipped to the acquired axis and the third
    element names, per line, the share (%) clipped away ({label: pct}).
    """
    ctx = make_context(recipe)
    x = None
    q_lo, q_hi = (1.0 - coverage) / 2.0, 1.0 - (1.0 - coverage) / 2.0
    highs, lows, cums = [], [], []
    for site in recipe.sites:
        if site.model in _NO_TAIL_MODELS:
            continue
        x, cum = _cumulative_area(simulate_site(site, ctx), ctx.x_ppm)
        if cum.size < 2 or cum[-1] <= 0:
            continue
        cum = cum / cum[-1]
        lows.append(float(np.interp(q_lo, cum, x)))
        highs.append(float(np.interp(q_hi, cum, x)))
        cums.append((site.label or site.model, cum))
    if not highs:
        return float(ctx.x_ppm.max()), float(ctx.x_ppm.min()), {}
    hi, lo = max(highs), min(lows)
    clipped: dict = {}
    if exp_ppm is not None and len(exp_ppm):
        e = np.asarray(exp_ppm, float)
        ehi, elo = float(e.max()), float(e.min())
        if hi > ehi or lo < elo:
            hi, lo = min(hi, ehi), max(lo, elo)
            for label, cum in cums:
                inside = (float(np.interp(hi, x, cum))
                          - float(np.interp(lo, x, cum)))
                pct = 100.0 * max(0.0, 1.0 - inside)
                if pct > 0.01:
                    clipped[label] = pct
    return hi, lo, clipped
