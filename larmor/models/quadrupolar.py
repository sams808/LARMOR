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
