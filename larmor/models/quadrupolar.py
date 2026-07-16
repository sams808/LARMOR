"""Quadrupolar models: Czjzek distribution and discrete second-order CT sites.

Both share the same physics engine (mrsimulator BlochDecayCTSpectrum):
  - czjzek reweights a precomputed (Cq, eta) kernel -- fast in fits
  - quad_ct simulates one site on demand with an LRU cache -- exact in Cq/eta
"""
from __future__ import annotations

from functools import lru_cache

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


# --------------------------------------------------------------------------
# Czjzek distribution (kernel-reweighting; kernel built once in larmor.engine)

def _render_czjzek(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor import engine

    kernel = engine.build_kernel(ctx.nucleus, ctx.larmor_MHz, ctx.spin_rate_Hz)
    y = kernel.weights(v["sigma_Cq_MHz"]) @ kernel.K
    y = _broaden_shift(kernel.x_ppm, y, v["isotropic_chemical_shift_ppm"],
                       v["shift_fwhm_ppm"])
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
        ParamDef("shift_fwhm_ppm", "fwhm", 10.0, "ppm",
                 "isotropic shift distribution width", min=0.1),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_czjzek,
))


# --------------------------------------------------------------------------
# discrete second-order quadrupolar CT lineshape (crystalline sites)

@lru_cache(maxsize=256)
def _quad_ct_shape(nucleus: str, larmor_khz: int, mas_hz: int,
                   cq_khz: int, eta_milli: int,
                   x0: float, x1: float, npts: int) -> tuple:
    """Simulate one quadrupolar site at delta_iso=0 (cached on rounded params)."""
    from mrsimulator import Simulator, Site, SpinSystem
    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.method import SpectralDimension
    from mrsimulator.spin_system.isotope import Isotope
    from mrsimulator.spin_system.tensors import SymmetricTensor

    larmor_MHz = larmor_khz / 1000.0
    B0 = larmor_MHz / abs(Isotope(symbol=nucleus).gyromagnetic_ratio)
    site = Site(isotope=nucleus, isotropic_chemical_shift=0.0,
                quadrupolar=SymmetricTensor(Cq=cq_khz * 1e3, eta=eta_milli / 1000.0))
    sw = abs(x1 - x0) * larmor_MHz
    method = BlochDecayCTSpectrum(
        channels=[nucleus], magnetic_flux_density=B0, rotor_frequency=mas_hz,
        spectral_dimensions=[SpectralDimension(
            count=npts, spectral_width=sw,
            reference_offset=(x0 + x1) / 2.0 * larmor_MHz)],
    )
    sim = Simulator(spin_systems=[SpinSystem(sites=[site])], methods=[method])
    sim.config.number_of_sidebands = 8
    sim.run()
    ds = sim.methods[0].simulation
    coords = ds.x[0].coordinates
    x = coords.value if str(coords.unit) == "ppm" else coords.to("Hz").value / larmor_MHz
    y = np.asarray(ds.y[0].components[0].real, dtype=float)
    order = np.argsort(x)
    return tuple(x[order]), tuple(y[order])


def _render_quad_ct(v: dict, ctx: SimContext) -> np.ndarray:
    x0, x1 = float(ctx.x_ppm[0]), float(ctx.x_ppm[-1])
    npts = min(len(ctx.x_ppm), 2048)
    xs, ys = _quad_ct_shape(
        ctx.nucleus, int(round(ctx.larmor_MHz * 1000)),
        int(round(ctx.spin_rate_Hz)),
        int(round(v["Cq_MHz"] * 1000)), int(round(v["eta"] * 1000)),
        round(x0, 2), round(x1, 2), npts)
    xs, ys = np.array(xs), np.array(ys)
    y = _broaden_shift(xs, ys, v["isotropic_chemical_shift_ppm"],
                       v["shift_fwhm_ppm"])
    peak = y.max()
    y = v["amplitude"] * (y / peak) if peak > 0 else y
    return np.interp(ctx.x_ppm, xs, y, left=0.0, right=0.0)


register(Model(
    name="quad_ct",
    label="Quadrupolar CT (discrete)",
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
