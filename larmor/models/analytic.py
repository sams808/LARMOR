"""Analytic lineshapes: Gauss/Lorentz (pseudo-Voigt), as in dmfit's Gaus/Lor."""
from __future__ import annotations

import numpy as np

from larmor.models.base import Model, ParamDef, SimContext, register

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def gauss_lor(x_ppm: np.ndarray, pos_ppm: float, fwhm_ppm: float,
              amplitude: float, gl: float) -> np.ndarray:
    """Peak-normalized pseudo-Voigt: y = gl*Gaussian + (1-gl)*Lorentzian."""
    dx = x_ppm - pos_ppm
    sig = max(fwhm_ppm, 1e-6) * FWHM_TO_SIGMA
    g = np.exp(-0.5 * (dx / sig) ** 2)
    hwhm = max(fwhm_ppm, 1e-6) / 2.0
    l = 1.0 / (1.0 + (dx / hwhm) ** 2)
    return amplitude * (gl * g + (1.0 - gl) * l)


def _render(v: dict, ctx: SimContext) -> np.ndarray:
    return gauss_lor(ctx.x_ppm, v["isotropic_chemical_shift_ppm"],
                     v["shift_fwhm_ppm"], v["amplitude"], v.get("gl", 1.0))


register(Model(
    name="gauss_lor",
    label="Gauss/Lorentz",
    description="Pseudo-Voigt peak (dmfit's Gaus/Lor). gl=1 pure Gaussian, "
                "gl=0 pure Lorentzian.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm", "peak position"),
        ParamDef("shift_fwhm_ppm", "fwhm", 5.0, "ppm", "full width at half maximum",
                 min=0.1),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
        ParamDef("gl", "gl", 1.0, "", "Gaussian fraction (1=G, 0=L)",
                 min=0.0, max=1.0, vary=False),
    ),
    render=_render,
))
