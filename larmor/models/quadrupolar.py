"""Quadrupolar models: Czjzek distributions and discrete second-order CT sites.

All share the same physics engine (mrsimulator BlochDecayCTSpectrum):
  - czjzek reweights a precomputed (Cq, eta) kernel -- fast in fits
  - czjzek_d: the same kernel with Czjzek's dimensionality d as a parameter
    (d = 5 == czjzek == Le Caer & Brand's Gaussian Isotropic Model; dmfit
    CzSimple's <d>)
  - czjzek_corr: d = 5 weights, each C_Q group translated by a linear
    (delta_iso, C_Q) correlation about the mean C_Q (slope 0 == czjzek)
  - ext_czjzek / amorphous: other reweightings of the same kernel
  - quad_ct etc. simulate one site on demand with an LRU cache -- exact in Cq/eta
"""
from __future__ import annotations

import numpy as np

from larmor.models.base import Model, ParamDef, SimContext, register
#: upper bound offered for a discrete Cq. 40 MHz was arbitrary and too low
#: for the heavy quadrupolar halides this program is used on -- a real 81Br
#: glass needs ~34 MHz, so the old ceiling sat right on top of the answer.
CQ_MAX_MHZ = 120.0

from larmor.models.analytic import FWHM_TO_SIGMA


def _clamp_fwhm(x: np.ndarray, fwhm_ppm: float) -> float:
    """A Gaussian broadening is never wider than the axis it is applied on.

    The optimiser's trial steps are not physical: scipy's trust-region
    (least_squares "trf") opens its first radius at the largest parameter
    magnitude -- the AMPLITUDE, ~3e6 on a Bruker spectrum -- so the first
    step can move a width by that much. Measured on the 27Al example with a
    czjzek_corr line and five linked sideband copies: the second residual
    evaluation asked for dCS = 3.15e6 ppm, a 2.4-million-point Gaussian on a
    2560-point axis, and gaussian_filter1d (cost ~ npts x 8 sigma) ran for
    more than 40 minutes on ONE evaluation -- the "long to compute" report,
    and why Stop could not get a word in. A broadening wider than the axis
    is already a flat smear that _finish renormalises to a plateau; clamping
    it to the axis span leaves that outcome and every physical width (a
    fraction of the span) untouched, and bounds one evaluation to
    ~4 npts^2 operations (~20 ms at 2560 points)."""
    if x.size < 2:
        return fwhm_ppm
    return min(float(fwhm_ppm), abs(float(x[-1]) - float(x[0])))


def _broaden_shift(x: np.ndarray, y: np.ndarray, pos_ppm: float,
                   fwhm_ppm: float) -> np.ndarray:
    """Translate a delta_iso=0 lineshape to pos and apply Gaussian broadening."""
    y = np.interp(x - pos_ppm, x, y, left=0.0, right=0.0)
    dppm = abs(x[1] - x[0])
    sigma_pts = _clamp_fwhm(x, fwhm_ppm) * FWHM_TO_SIGMA / dppm
    if sigma_pts > 0.05:
        from scipy.ndimage import gaussian_filter1d   # deferred: startup cost
        y = gaussian_filter1d(y, sigma_pts, mode="constant")
    return y


def _czjzek_fwhm(v: dict) -> float:
    """Total 1D Gaussian broadening: the isotropic-shift distribution (dmfit's
    dCS) and the round point broadening (dmfit's wid) both blur the single MAS
    dimension, so they add in quadrature.  In 2D they differ (diagonal vs
    round) -- see larmor.twod.simulate_site_2d."""
    cs = float(v.get("shift_fwhm_ppm", 0.0))
    line = float(v.get("line_fwhm_ppm", 0.0))
    return float(np.hypot(cs, line)) if line > 0.0 else cs


def _lorentz_convolve(y: np.ndarray, fwhm_pts: float) -> np.ndarray:
    """Convolve with a normalised Lorentzian of the given FWHM (in points).

    ``np.convolve(a, v, mode="same")`` returns length ``max(len(a), len(v))``,
    NOT ``len(a)`` -- if the kernel `k` (sized off `fwhm_pts`, unbounded) ever
    comes out longer than `y`, the result silently grows past `y`'s length,
    which breaks every caller's assumption that broadening preserves the
    array length (crashing downstream, e.g. `_render_amorphous`'s final
    `np.interp(ctx.x_ppm, kernel.x_ppm, y, ...)` with "fp and xp are not of
    the same length"). Real trigger: an lmfit errorbar-rescue retry step can
    push a poorly-determined FWHM parameter far outside typical values while
    probing the Jacobian -- caught fitting real 11B glass data (LARMOR
    validation pass, 2026-08). Clamp the kernel to never exceed `y`'s length;
    a Lorentzian that wide relative to the data is already a bad fit the
    optimiser should reject via the residual, not something worth crashing
    over."""
    hwhm = fwhm_pts / 2.0
    if hwhm <= 0.0 or y.size == 0:
        return y
    half = int(np.ceil(hwhm * 20.0))
    half = min(half, max(0, (y.size - 1) // 2))
    t = np.arange(-half, half + 1)
    k = 1.0 / (1.0 + (t / hwhm) ** 2)
    k /= k.sum()
    return np.convolve(y, k, mode="same")


def _broaden_shift_pv(x: np.ndarray, y: np.ndarray, pos_ppm: float,
                      gauss_fwhm_ppm: float, lor_fwhm_ppm: float,
                      gl: float) -> np.ndarray:
    """Translate to pos, then apply a Gaussian (shift distribution) and a
    pseudo-Voigt line broadening.  `gl` is dmfit's Gaus/Lor mix: gl=1 makes the
    `lor_fwhm` broadening Gaussian, gl=0 makes it Lorentzian (the Amorphous
    default).  The Gaussian shift distribution and the Gaussian part of the line
    add in quadrature; the Lorentzian part is convolved separately."""
    y = np.interp(x - pos_ppm, x, y, left=0.0, right=0.0)
    dppm = abs(x[1] - x[0])
    gl = float(np.clip(gl, 0.0, 1.0))
    lor_fwhm_ppm = max(lor_fwhm_ppm, 0.0)          # negative lb (dmfit resolution
    #                                                enhancement) is not a convolution
    g_fwhm = _clamp_fwhm(x, float(np.hypot(max(gauss_fwhm_ppm, 0.0),
                                           gl * lor_fwhm_ppm)))
    sigma_pts = g_fwhm * FWHM_TO_SIGMA / dppm
    if sigma_pts > 0.05:
        from scipy.ndimage import gaussian_filter1d   # deferred: startup cost
        y = gaussian_filter1d(y, sigma_pts, mode="constant")
    l_fwhm = (1.0 - gl) * lor_fwhm_ppm
    if l_fwhm > 0.5 * dppm:
        y = _lorentz_convolve(y, l_fwhm / dppm)
    return y


def _gaussian_weight(grid: np.ndarray, mean: float, fwhm: float) -> np.ndarray:
    """Weights of a Gaussian(mean, fwhm) sampled on `grid`.  When the FWHM is
    below one grid step the distribution is a delta: put the weight on the two
    grid points bracketing `mean` (linear interpolation) so an off-grid mean is
    still placed accurately."""
    grid = np.asarray(grid, float)
    step = float(grid[1] - grid[0]) if grid.size > 1 else 1.0
    w = np.zeros_like(grid)
    if fwhm is None or fwhm < step:                # delta -> nearest-two interp
        mean = float(np.clip(mean, grid[0], grid[-1]))
        j = int(np.clip(np.searchsorted(grid, mean) - 1, 0, grid.size - 2))
        frac = (mean - grid[j]) / (grid[j + 1] - grid[j]) if grid.size > 1 else 0.0
        w[j] = 1.0 - frac
        w[j + 1] += frac
        return w
    sig = fwhm * FWHM_TO_SIGMA
    w = np.exp(-0.5 * ((grid - mean) / sig) ** 2)
    return w


# --------------------------------------------------------------------------
# Czjzek distribution (kernel-reweighting; kernel built once in larmor.engine)

#: how far, in units of the stored sigma, the (Cq, eta) kernel grid must reach
#: for a Czjzek family model. The stored sigma is mrsimulator's, whose
#: distribution has sigma_Cz = 2 sigma: the |Cq| marginal peaks at 3.73 sigma
#: and 21.3 % of the d = 5 mass lies beyond 5 sigma (4.5e-5 beyond 10 sigma).
#: The old 5-sigma request dropped -- and renormalised away -- a fifth of the
#: distribution whenever 5 sigma sat just under a ladder step. Shared by
#: czjzek, czjzek_d and czjzek_corr so the three ask for ONE cached kernel and
#: their identities (d = 5, slope = 0) hold exactly.
CZJZEK_KERNEL_HEADROOM = 10.0


def _kernel_for(ctx: SimContext, needed_cq_MHz: float):
    """The cached (Cq, eta) CT basis covering ``ctx`` up to ``needed_cq_MHz``
    (quantised to engine.CQ_MAX_LADDER), at the engine's resolution rules."""
    from larmor import engine

    sw, ref = engine.kernel_window(ctx.x_ppm, ctx.larmor_MHz)
    cq_max = engine.kernel_cq_max(float(needed_cq_MHz))
    return engine.build_kernel(
        ctx.nucleus, ctx.larmor_MHz, ctx.spin_rate_Hz, sw_Hz=sw,
        npts=min(int(engine.KERNEL_SETTINGS["npts"]
                     * max(1.0, sw / engine.KERNEL_MIN_SW_HZ)), 16384),
        ref_offset_ppm=ref, cq_max_MHz=cq_max,
        n_cq=max(engine.KERNEL_SETTINGS["n_cq"],
                 int(engine.KERNEL_SETTINGS["n_cq"] * cq_max / 25.0)),
        n_eta=int(engine.KERNEL_SETTINGS["n_eta"]))


#: per-kernel memo of the (sigma[, d])-dependent reweighting, the part of a
#: Czjzek-family render that does NOT depend on position or amplitude. A
#: fit with linked spinning-sideband copies renders the same shape five or
#: six times per residual evaluation at different positions; the copies now
#: share one w @ K (czjzek, czjzek_d) or one eta-summed basis (czjzek_corr)
#: per evaluation instead of each redoing the 880 x npts multiply. The memo
#: lives on the kernel object (it dies with it; no id() reuse hazard) and
#: holds the last few entries -- a Jacobian probe alternates a handful of
#: parameter sets. Results are bit-identical: the cached array IS the array
#: the uncached render would have computed.
_REWEIGHT_MEMO_SIZE = 8


def _reweight(kernel, key: tuple, compute):
    """``compute()`` memoised on ``kernel`` under ``key``."""
    memo = kernel.__dict__.get("_reweight_memo")
    if memo is None:
        memo = kernel.__dict__["_reweight_memo"] = {}
    hit = memo.get(key)
    if hit is None:
        hit = compute()
        if len(memo) >= _REWEIGHT_MEMO_SIZE:
            memo.pop(next(iter(memo)))           # oldest insertion
        memo[key] = hit
    return hit


def _finish(kernel, y: np.ndarray, v: dict, ctx: SimContext) -> np.ndarray:
    """Peak-normalise to the site amplitude and move from the kernel axis to
    the context axis (no interpolation when they coincide)."""
    peak = y.max()
    y = v["amplitude"] * (y / peak) if peak > 0 else y
    if kernel.x_ppm.shape == ctx.x_ppm.shape and \
            np.allclose(kernel.x_ppm, ctx.x_ppm):
        return y
    return np.interp(ctx.x_ppm, kernel.x_ppm, y, left=0.0, right=0.0)


def _render_czjzek(v: dict, ctx: SimContext) -> np.ndarray:
    kernel = _kernel_for(
        ctx, CZJZEK_KERNEL_HEADROOM * float(v.get("sigma_Cq_MHz", 2.0)))
    sigma = float(v["sigma_Cq_MHz"])
    y = _reweight(kernel, ("czjzek", sigma),
                  lambda: kernel.weights(sigma) @ kernel.K)
    y = _broaden_shift(kernel.x_ppm, y, v["isotropic_chemical_shift_ppm"],
                       _czjzek_fwhm(v))
    return _finish(kernel, y, v, ctx)


register(Model(
    name="czjzek",
    label="Czjzek (quad. distribution)",
    description="Czjzek distribution of quadrupolar tensors for disordered "
                "materials (dmfit's CzSimple). sigma is HALF of dmfit's sCZ_CQ.",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        # max: the kernel's Cq ladder tops out at 400 MHz and the render
        # requests 10*sigma of headroom (CZJZEK_KERNEL_HEADROOM: 21 % of the
        # mass lies beyond 5 sigma), so sigma beyond 40 MHz cannot be
        # represented -- it saturated into a plain Gaussian that LOOKED
        # converged. At the bound the fit's at-bounds diagnosis fires
        # instead. (The largest published glass sigmas are ~20 MHz.)
        ParamDef("sigma_Cq_MHz", "sigma", 2.0, "MHz",
                 "Czjzek width parameter (mrsimulator sigma; sigma_Cz = "
                 "dmfit sCZ_CQ = 2 sigma; mode of |Cq| = 3.7 sigma)",
                 min=0.05, max=40.0),
        ParamDef("shift_fwhm_ppm", "dCS", 10.0, "ppm",
                 "isotropic-shift distribution FWHM (dmfit dCS; diagonal in 2D)",
                 min=0.1),
        ParamDef("line_fwhm_ppm", "line", 0.0, "ppm",
                 "round point/line broadening (dmfit wid; isotropic in 2D) -- "
                 "held at 0 unless freed, like dmfit's greyed Lb",
                 min=0.0, vary=False, default_fixed=True),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_czjzek,
))


# --------------------------------------------------------------------------
# extended Czjzek: perturbation of a dominant tensor (same kernel grid)

def _render_ext_czjzek(v: dict, ctx: SimContext) -> np.ndarray:
    from mrsimulator.models import ExtCzjzekDistribution

    kernel = _kernel_for(ctx, 2.5 * float(v.get("Cq_MHz", 5.0))
                         * (1.0 + float(v.get("eps", 0.3))))
    # the dominant tensor must share the pdf grid's unit system (MHz here)
    dominant = {"Cq": v["Cq_MHz"], "eta": v["eta"]}
    res = ExtCzjzekDistribution(dominant, eps=max(v["eps"], 1e-3)).pdf(
        pos=[kernel.cq_grid_MHz, kernel.eta_grid])
    amp = np.asarray(res[-1] if isinstance(res, (tuple, list)) else res)
    w = amp.ravel()
    s = w.sum()
    if s > 0:
        w = w / s
    y = w @ kernel.K
    y = _broaden_shift(kernel.x_ppm, y, v["isotropic_chemical_shift_ppm"],
                       _czjzek_fwhm(v))
    return _finish(kernel, y, v, ctx)


register(Model(
    name="ext_czjzek",
    label="ext. Czjzek",
    description="Extended Czjzek: random perturbation (eps) around a dominant "
                "quadrupolar tensor -- partially ordered environments.",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("Cq_MHz", "cq", 5.0, "MHz", "dominant quadrupolar coupling",
                 min=0.05, max=CQ_MAX_MHZ),
        ParamDef("eta", "eta", 0.2, "", "dominant asymmetry", min=0.0, max=1.0),
        ParamDef("eps", "eps", 0.3, "", "perturbation fraction",
                 min=0.01, max=3.0),
        ParamDef("shift_fwhm_ppm", "dCS", 5.0, "ppm",
                 "isotropic-shift distribution FWHM (dmfit dCS; diagonal in 2D)",
                 min=0.1),
        ParamDef("line_fwhm_ppm", "line", 0.0, "ppm",
                 "round point/line broadening (dmfit wid; isotropic in 2D) -- "
                 "held at 0 unless freed, like dmfit's greyed Lb",
                 min=0.0, vary=False, default_fixed=True),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_ext_czjzek,
))


# --------------------------------------------------------------------------
# Czjzek, general d: Czjzek's dimensionality as a parameter (dmfit CzSimple <d>)
#
# Le Caer & Brand's Gaussian Isotropic Model (five i.i.d. Gaussian components
# of the traceless symmetric EFG tensor) IS the d = 5 Czjzek distribution that
# mrsimulator's CzjzekDistribution implements, so a separate "GIM" model would
# be byte-identical to `czjzek`; their extended GIM is `ext_czjzek`. What
# mrsimulator does not offer is Czjzek's exponent d - 1 as a free parameter,
# which dmfit exposes as CzSimple's <d> box. This model adds it on the shared
# kernel through czjzek_dist.czjzek_weights, which reproduces
# CzjzekKernel.weights at d = 5 -- so czjzek_d at d = 5 renders exactly as
# czjzek on the same cached kernel. The parameter is named czjzek_d (not d):
# the `function` model already owns a parameter literally called d, and
# table columns, labels and multi-field sharing are keyed by parameter name.

def _render_czjzek_d(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor import czjzek_dist

    sigma = float(v.get("sigma_Cq_MHz", 2.0))
    d = float(v.get("czjzek_d", 5.0))
    kernel = _kernel_for(ctx, CZJZEK_KERNEL_HEADROOM * sigma)

    def _basis():
        w = czjzek_dist.czjzek_weights(sigma, d, kernel.cq_grid_MHz,
                                       kernel.eta_grid)
        return w @ kernel.K

    y = _reweight(kernel, ("czjzek_d", sigma, d), _basis)
    y = _broaden_shift(kernel.x_ppm, y, v["isotropic_chemical_shift_ppm"],
                       _czjzek_fwhm(v))
    return _finish(kernel, y, v, ctx)


register(Model(
    name="czjzek_d",
    label="Czjzek, general d  (GIM at d = 5)",
    description="Czjzek distribution with Czjzek's dimensionality d as a "
                "parameter (dmfit CzSimple <d>). d = 5 is the fully isotropic "
                "Gaussian Isotropic Model == the plain czjzek model; d < 5 is "
                "an empirical, more constrained family.",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("sigma_Cq_MHz", "sigma", 2.0, "MHz",
                 "Czjzek width parameter (mrsimulator sigma; sigma_Cz = "
                 "dmfit sCZ_CQ = 2 sigma)",
                 min=0.05, max=40.0),
        ParamDef("czjzek_d", "d", 5.0, "",
                 "Czjzek dimensionality d (5 = GIM/standard Czjzek; set 3-4 "
                 "before freeing it)", min=1.0, max=5.0, vary=False),
        ParamDef("shift_fwhm_ppm", "dCS", 10.0, "ppm",
                 "isotropic-shift distribution FWHM (dmfit dCS)", min=0.1),
        ParamDef("line_fwhm_ppm", "line", 0.0, "ppm",
                 "round point/line broadening (dmfit wid) -- held at 0 unless "
                 "freed, like dmfit's greyed Lb",
                 min=0.0, vary=False, default_fixed=True),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_czjzek_d,
))


# --------------------------------------------------------------------------
# Czjzek + correlated isotropic shift: delta_iso = pos + slope (C_Q - <C_Q>)
#
# Real glasses correlate the isotropic shift with the quadrupolar coupling
# (17O, 23Na, 27Al MQMAS: Vermillion 1998, Clark 2004, Vasconcelos 2013). The
# kernel rows are eta-major, so one eta-summed subspectrum per C_Q group is
# formed (80 rows on the 25 MHz step, not 880), each translated by the
# correlated shift and summed; the residual dCS and line broadening follow.

def _shift_sum(x: np.ndarray, Y: np.ndarray, shifts: np.ndarray,
               weights: np.ndarray | None = None, tol: float = 1e-6
               ) -> np.ndarray:
    """Sum the rows of ``Y`` (each sampled on the ascending axis ``x``) after
    translating row q by ``shifts[q]`` along x -- linear interpolation, zero
    outside the axis, the primitive _broaden_shift uses. When ``weights`` is
    given, rows below ``tol`` times the largest weight are skipped (they
    carry no visible intensity)."""
    x = np.asarray(x, float)
    Y = np.asarray(Y, float)
    shifts = np.asarray(shifts, float)
    keep = np.ones(Y.shape[0], dtype=bool)
    if weights is not None:
        weights = np.asarray(weights, float)
        keep = weights >= tol * (weights.max() if weights.size else 0.0)
    # a per-row np.interp loop: a vectorised searchsorted + fancy-indexing
    # version reproduced it bit for bit but ran 6x SLOWER (31 vs 4.7 ms for
    # 80 rows x 2048 points), so the plain loop stays
    out = np.zeros(x.shape[0])
    for q in np.flatnonzero(keep):
        out += np.interp(x - shifts[q], x, Y[q], left=0.0, right=0.0)
    return out


def _render_czjzek_corr(v: dict, ctx: SimContext) -> np.ndarray:
    sigma = float(v.get("sigma_Cq_MHz", 2.0))
    kernel = _kernel_for(ctx, CZJZEK_KERNEL_HEADROOM * sigma)
    n_eta, n_cq = kernel.eta_grid.size, kernel.cq_grid_MHz.size

    def _basis():
        # kernel rows follow np.meshgrid(cq, eta, indexing='xy').ravel():
        # eta-major. One eta-summed subspectrum per C_Q: the shift
        # correlates with C_Q only.
        W = kernel.weights(sigma).reshape(n_eta, n_cq)
        Y = np.einsum("eq,eqx->qx", W, kernel.K.reshape(n_eta, n_cq, -1))
        wq = W.sum(axis=0)
        return Y, wq, float(wq @ kernel.cq_grid_MHz)

    Y, wq, cq_mean = _reweight(kernel, ("czjzek_corr", sigma), _basis)
    pos = float(v["isotropic_chemical_shift_ppm"])
    slope = float(v.get("shift_slope_ppm_per_MHz", 0.0))
    # pivot at <C_Q>: pos stays the ensemble-MEAN shift for any slope and the
    # centre of gravity is slope-invariant (sum_q wq (C_Q - <C_Q>) = 0), so
    # pos and slope are decorrelated in the fit
    shifts = pos + slope * (kernel.cq_grid_MHz - cq_mean)
    y = _shift_sum(kernel.x_ppm, Y, shifts, weights=wq)
    y = _broaden_shift(kernel.x_ppm, y, 0.0, _czjzek_fwhm(v))
    return _finish(kernel, y, v, ctx)


register(Model(
    name="czjzek_corr",
    label="Czjzek + δiso–C_Q correlation",
    description="Czjzek (d = 5) distribution whose isotropic shift depends "
                "linearly on C_Q about the ensemble mean: δiso = pos + "
                "slope·(C_Q − ⟨C_Q⟩), plus a residual Gaussian shift "
                "distribution dCS. slope = 0 is the plain Czjzek model.",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "ensemble-mean isotropic chemical shift"),
        ParamDef("sigma_Cq_MHz", "sigma", 2.0, "MHz",
                 "Czjzek width parameter (mrsimulator sigma; sigma_Cz = "
                 "dmfit sCZ_CQ = 2 sigma)",
                 min=0.05, max=40.0),
        ParamDef("shift_slope_ppm_per_MHz", "slope", 0.0, "ppm/MHz",
                 "dδiso/dC_Q about the mean C_Q (0 = plain Czjzek)",
                 min=-50.0, max=50.0),
        ParamDef("shift_fwhm_ppm", "dCS", 5.0, "ppm",
                 "residual shift-distribution FWHM at fixed C_Q", min=0.1),
        ParamDef("line_fwhm_ppm", "line", 0.0, "ppm",
                 "round point/line broadening (dmfit wid) -- held at 0 unless "
                 "freed, like dmfit's greyed Lb",
                 min=0.0, vary=False, default_fixed=True),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_czjzek_corr,
))


# --------------------------------------------------------------------------
# Amorphous: independent Gaussian distributions of Cq and eta (dmfit "Amorphous")
#
# Unlike Czjzek (a coupled d=5 distribution whose |Cq| mode is 2*sigma and which
# always includes Cq -> 0), the Amorphous model spreads a *non-zero mean* Cq and
# eta by INDEPENDENT Gaussians.  That is the right statistics for a well-defined
# quadrupolar site with modest disorder -- e.g. trigonal BO3 boron in borate /
# borosilicate glasses (11B Cq ~ 2.5-2.7 MHz, eta ~ 0.1-0.2), which Czjzek fits
# poorly.  We reuse the validated (Cq, eta) CT kernel (larmor.engine.build_kernel)
# on a finer Cq grid and reweight it by Gaussian(Cq) x Gaussian(eta).

#: finer kernel for the (narrow) Amorphous Cq distribution; cached per process
AMORPH_CQ_MAX = 6.0
AMORPH_N_CQ = 120
AMORPH_N_ETA = 11


def _amorphous_weights(kernel, cq_MHz: float, eta: float,
                       cq_fwhm_MHz: float, eta_fwhm: float) -> np.ndarray:
    """Row weights for kernel.K: independent Gaussians on the Cq and eta grids,
    ordered to match np.meshgrid(cq_grid, eta_grid, indexing='xy').ravel()."""
    wq = _gaussian_weight(kernel.cq_grid_MHz, cq_MHz, cq_fwhm_MHz)
    we = _gaussian_weight(kernel.eta_grid, eta, eta_fwhm)
    W = np.outer(we, wq)                    # (n_eta, n_cq) == meshgrid xy layout
    w = W.ravel()
    s = w.sum()
    return w / s if s > 0 else w


def _render_amorphous(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor import engine

    sw, ref = engine.kernel_window(ctx.x_ppm, ctx.larmor_MHz)
    kernel = engine.build_kernel(
        ctx.nucleus, ctx.larmor_MHz, ctx.spin_rate_Hz, sw_Hz=sw,
        npts=min(int(engine.KERNEL_SETTINGS["npts"]
                     * max(1.0, sw / engine.KERNEL_MIN_SW_HZ)), 16384),
        ref_offset_ppm=ref,
        cq_max_MHz=AMORPH_CQ_MAX, n_cq=AMORPH_N_CQ, n_eta=AMORPH_N_ETA)
    w = _amorphous_weights(kernel, v["Cq_MHz"], v.get("eta", 0.0),
                           v.get("Cq_fwhm_MHz", 0.0), v.get("eta_fwhm", 0.0))
    y = w @ kernel.K
    y = _broaden_shift_pv(kernel.x_ppm, y, v["isotropic_chemical_shift_ppm"],
                          v.get("shift_fwhm_ppm", 0.0),
                          v.get("line_fwhm_ppm", 0.0), v.get("gl", 0.0))
    return _finish(kernel, y, v, ctx)


register(Model(
    name="amorphous",
    label="Amorphous (Gaussian Cq/eta dist.)",
    description="dmfit's 'Amorphous': a second-order quadrupolar CT lineshape "
                "with INDEPENDENT Gaussian distributions of Cq and eta, a "
                "chemical-shift distribution (dCS) and a pseudo-Voigt line "
                "broadening (lb/gl). For BO3 in 11B and other well-defined "
                "quadrupolar sites with modest disorder (unlike Czjzek).",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 15.0, "ppm",
                 "isotropic chemical shift (dmfit pos)"),
        ParamDef("Cq_MHz", "cq", 2.6, "MHz",
                 "mean quadrupolar coupling (dmfit CQ/1000; nuQ = 3Cq/[2I(2I-1)])",
                 min=0.05, max=AMORPH_CQ_MAX),
        ParamDef("eta", "eta", 0.1, "", "mean quadrupolar asymmetry (dmfit etaQ)",
                 min=0.0, max=1.0),
        ParamDef("Cq_fwhm_MHz", "dcq", 0.3, "MHz",
                 "Gaussian FWHM of the Cq distribution (dmfit FWHM_CQ/1000)",
                 min=0.0, max=3.0),
        ParamDef("eta_fwhm", "deta", 0.0, "",
                 "Gaussian FWHM of the eta distribution (dmfit FWHM_etaQ)",
                 min=0.0, max=1.0),
        ParamDef("shift_fwhm_ppm", "dCS", 3.0, "ppm",
                 "isotropic-shift distribution FWHM (dmfit dCS, Gaussian)",
                 min=0.0),
        ParamDef("line_fwhm_ppm", "lb", 0.5, "ppm",
                 "line broadening (dmfit lb; Lorentzian when gl=0) -- held "
                 "unless freed, like dmfit's greyed Lb",
                 min=0.0, vary=False, default_fixed=True),
        ParamDef("gl", "gl", 0.0, "",
                 "Gaus/Lor mix of the line broadening (1=Gaussian, 0=Lorentzian)",
                 min=0.0, max=1.0, vary=False),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_amorphous,
))


# --------------------------------------------------------------------------
# discrete second-order quadrupolar CT lineshape (crystalline sites)

def _render_quad_ct(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor.models._singlesite import render_single_site

    return render_single_site(v, ctx, cq_key="Cq_MHz", eta_q_key="eta",
                              ct_only=True, n_ssb=8)


register(Model(
    name="quad_ct",
    label="Quad CT (2nd order)",
    description="Second-order quadrupolar central-transition lineshape for a "
                "single crystalline site (MAS or static via spin rate).",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("Cq_MHz", "cq", 3.0, "MHz", "quadrupolar coupling constant",
                 min=0.01, max=CQ_MAX_MHZ),
        ParamDef("eta", "eta", 0.2, "", "quadrupolar asymmetry", min=0.0, max=1.0),
        ParamDef("shift_fwhm_ppm", "fwhm", 2.0, "ppm", "Gaussian broadening",
                 min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_quad_ct,
))


# --------------------------------------------------------------------------
# first-order quadrupolar: full satellite manifold with spinning sidebands

def _render_quad_first(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor.models._singlesite import render_single_site

    return render_single_site(v, ctx, cq_key="Cq_MHz", eta_q_key="eta",
                              ct_only=False, n_ssb=64)


register(Model(
    name="quad_first",
    label="Quad 1st order (satellites)",
    description="Full quadrupolar pattern including satellite transitions and "
                "their spinning-sideband manifold (dmfit's 'quad 1st order').",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("Cq_MHz", "cq", 1.0, "MHz", "quadrupolar coupling constant",
                 min=0.001, max=CQ_MAX_MHZ),
        ParamDef("eta", "eta", 0.1, "", "quadrupolar asymmetry", min=0.0, max=1.0),
        ParamDef("shift_fwhm_ppm", "fwhm", 1.0, "ppm", "Gaussian broadening",
                 min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_quad_first,
))


# --------------------------------------------------------------------------
# combined second-order quad CT + CSA on the same site

def _render_quad_csa(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor.models._singlesite import render_single_site

    return render_single_site(v, ctx, cq_key="Cq_MHz", eta_q_key="eta_q",
                              zeta_key="zeta_ppm", eta_cs_key="eta_cs",
                              ct_only=True, n_ssb=16)


register(Model(
    name="quad_csa",
    label="Quad CT + CSA",
    description="Central transition with BOTH second-order quadrupolar and "
                "shielding-anisotropy interactions on the same site.",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("Cq_MHz", "cq", 3.0, "MHz", "quadrupolar coupling constant",
                 min=0.01, max=CQ_MAX_MHZ),
        ParamDef("eta_q", "etaq", 0.2, "", "quadrupolar asymmetry",
                 min=0.0, max=1.0),
        ParamDef("zeta_ppm", "zeta", 50.0, "ppm", "shielding anisotropy",
                 min=-1000.0, max=1000.0),
        ParamDef("eta_cs", "etacs", 0.3, "", "shielding asymmetry",
                 min=0.0, max=1.0),
        ParamDef("shift_fwhm_ppm", "fwhm", 2.0, "ppm", "Gaussian broadening",
                 min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_quad_csa,
))
