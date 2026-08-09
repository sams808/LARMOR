"""Quadrupolar models: Czjzek distribution and discrete second-order CT sites.

Both share the same physics engine (mrsimulator BlochDecayCTSpectrum):
  - czjzek reweights a precomputed (Cq, eta) kernel -- fast in fits
  - quad_ct simulates one site on demand with an LRU cache -- exact in Cq/eta
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d

from larmor.models.base import Model, ParamDef, SimContext, register
from larmor.models.analytic import FWHM_TO_SIGMA


def _broaden_shift(x: np.ndarray, y: np.ndarray, pos_ppm: float,
                   fwhm_ppm: float) -> np.ndarray:
    """Translate a delta_iso=0 lineshape to pos and apply Gaussian broadening."""
    y = np.interp(x - pos_ppm, x, y, left=0.0, right=0.0)
    dppm = abs(x[1] - x[0])
    sigma_pts = fwhm_ppm * FWHM_TO_SIGMA / dppm
    if sigma_pts > 0.05:
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
    g_fwhm = float(np.hypot(max(gauss_fwhm_ppm, 0.0), gl * lor_fwhm_ppm))
    sigma_pts = g_fwhm * FWHM_TO_SIGMA / dppm
    if sigma_pts > 0.05:
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

def _render_czjzek(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor import engine

    kernel = engine.build_kernel(ctx.nucleus, ctx.larmor_MHz, ctx.spin_rate_Hz)
    y = kernel.weights(v["sigma_Cq_MHz"]) @ kernel.K
    y = _broaden_shift(kernel.x_ppm, y, v["isotropic_chemical_shift_ppm"],
                       _czjzek_fwhm(v))
    peak = y.max()
    y = v["amplitude"] * (y / peak) if peak > 0 else y
    if kernel.x_ppm.shape == ctx.x_ppm.shape and \
            np.allclose(kernel.x_ppm, ctx.x_ppm):
        return y
    return np.interp(ctx.x_ppm, kernel.x_ppm, y, left=0.0, right=0.0)


register(Model(
    name="czjzek",
    label="Czjzek (quad. distribution)",
    description="Czjzek distribution of quadrupolar tensors for disordered "
                "materials (dmfit's CzSimple). sigma is HALF of dmfit's sCZ_CQ.",
    needs_quadrupolar=True,
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("sigma_Cq_MHz", "sigma", 2.0, "MHz",
                 "Czjzek width parameter (mode of |Cq| = 2 sigma)", min=0.05),
        ParamDef("shift_fwhm_ppm", "dCS", 10.0, "ppm",
                 "isotropic-shift distribution FWHM (dmfit dCS; diagonal in 2D)",
                 min=0.1),
        ParamDef("line_fwhm_ppm", "line", 0.0, "ppm",
                 "round point/line broadening (dmfit wid; isotropic in 2D)",
                 min=0.0),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_czjzek,
))


# --------------------------------------------------------------------------
# extended Czjzek: perturbation of a dominant tensor (same kernel grid)

def _render_ext_czjzek(v: dict, ctx: SimContext) -> np.ndarray:
    from mrsimulator.models import ExtCzjzekDistribution

    from larmor import engine

    kernel = engine.build_kernel(ctx.nucleus, ctx.larmor_MHz, ctx.spin_rate_Hz)
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
    peak = y.max()
    y = v["amplitude"] * (y / peak) if peak > 0 else y
    if kernel.x_ppm.shape == ctx.x_ppm.shape and \
            np.allclose(kernel.x_ppm, ctx.x_ppm):
        return y
    return np.interp(ctx.x_ppm, kernel.x_ppm, y, left=0.0, right=0.0)


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
                 min=0.05, max=40.0),
        ParamDef("eta", "eta", 0.2, "", "dominant asymmetry", min=0.0, max=1.0),
        ParamDef("eps", "eps", 0.3, "", "perturbation fraction",
                 min=0.01, max=3.0),
        ParamDef("shift_fwhm_ppm", "dCS", 5.0, "ppm",
                 "isotropic-shift distribution FWHM (dmfit dCS; diagonal in 2D)",
                 min=0.1),
        ParamDef("line_fwhm_ppm", "line", 0.0, "ppm",
                 "round point/line broadening (dmfit wid; isotropic in 2D)",
                 min=0.0),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_ext_czjzek,
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

    kernel = engine.build_kernel(
        ctx.nucleus, ctx.larmor_MHz, ctx.spin_rate_Hz,
        cq_max_MHz=AMORPH_CQ_MAX, n_cq=AMORPH_N_CQ, n_eta=AMORPH_N_ETA)
    w = _amorphous_weights(kernel, v["Cq_MHz"], v.get("eta", 0.0),
                           v.get("Cq_fwhm_MHz", 0.0), v.get("eta_fwhm", 0.0))
    y = w @ kernel.K
    y = _broaden_shift_pv(kernel.x_ppm, y, v["isotropic_chemical_shift_ppm"],
                          v.get("shift_fwhm_ppm", 0.0),
                          v.get("line_fwhm_ppm", 0.0), v.get("gl", 0.0))
    peak = y.max()
    y = v["amplitude"] * (y / peak) if peak > 0 else y
    if kernel.x_ppm.shape == ctx.x_ppm.shape and \
            np.allclose(kernel.x_ppm, ctx.x_ppm):
        return y
    return np.interp(ctx.x_ppm, kernel.x_ppm, y, left=0.0, right=0.0)


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
                 "line broadening (dmfit lb; Lorentzian when gl=0)", min=0.0),
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
                 min=0.01, max=40.0),
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
                 min=0.001, max=40.0),
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
                 min=0.01, max=40.0),
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
