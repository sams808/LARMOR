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


def gl_unit_area(x_ppm: np.ndarray, pos_ppm: float, fwhm_ppm: float,
                 gl: float) -> np.ndarray:
    """Unit-AREA pseudo-Voigt (∫ = 1): area-normalized Gauss/Lorentz mix."""
    dx = x_ppm - pos_ppm
    w = max(fwhm_ppm, 1e-9)
    sig = w * FWHM_TO_SIGMA
    g = np.exp(-0.5 * (dx / sig) ** 2) / (sig * np.sqrt(2.0 * np.pi))
    hw = w / 2.0
    lo = (hw / np.pi) / (dx ** 2 + hw ** 2)
    return gl * g + (1.0 - gl) * lo


def _render_gl_norm(v: dict, ctx: SimContext) -> np.ndarray:
    return v["amplitude"] * gl_unit_area(
        ctx.x_ppm, v["isotropic_chemical_shift_ppm"], v["shift_fwhm_ppm"],
        v.get("gl", 1.0))


register(Model(
    name="gl_norm",
    label="Gauss/Lorentz (area)",
    description="Area-normalized Gauss/Lorentz (dmfit's GL Norm): amplitude is "
                "the integral, so amplitudes read directly as populations for "
                "quantification.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm", "peak position"),
        ParamDef("shift_fwhm_ppm", "fwhm", 5.0, "ppm", "FWHM", min=0.1),
        ParamDef("amplitude", "amp", 1.0, "", "area (integral)", min=0.0),
        ParamDef("gl", "gl", 1.0, "", "Gaussian fraction (1=G, 0=L)",
                 min=0.0, max=1.0, vary=False),
    ),
    render=_render_gl_norm,
))


def _render_jmultiplet(v: dict, ctx: SimContext) -> np.ndarray:
    """n equivalent scalar couplings → binomial (n+1)-line multiplet."""
    from math import comb

    n = int(round(v.get("n_j", 1)))
    j_hz = v.get("j_hz", 0.0)
    larmor = ctx.larmor_MHz or 1.0
    d_ppm = j_hz / larmor
    pos = v["isotropic_chemical_shift_ppm"]
    fwhm = v["shift_fwhm_ppm"]
    amp = v["amplitude"]
    gl = v.get("gl", 1.0)
    total = float(2 ** n)
    y = np.zeros_like(ctx.x_ppm)
    for i in range(n + 1):
        w = comb(n, i) / total
        centre = pos + (i - n / 2.0) * d_ppm
        y += amp * w * gauss_lor(ctx.x_ppm, centre, fwhm, 1.0, gl)
    return y


register(Model(
    name="jmultiplet",
    label="J-multiplet",
    description="Scalar-coupling multiplet: n equivalent spin-½ couplings give a "
                "binomial (n+1)-line pattern split by J (Hz). Each component is a "
                "Gauss/Lorentz line.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "multiplet centre"),
        ParamDef("j_hz", "j", 100.0, "Hz", "scalar coupling J", min=0.0),
        ParamDef("n_j", "nj", 1.0, "", "number of equivalent couplings (integer)",
                 min=0.0, max=12.0, vary=False),
        ParamDef("shift_fwhm_ppm", "fwhm", 1.0, "ppm", "component FWHM", min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "total intensity", min=0.0),
        ParamDef("gl", "gl", 1.0, "", "Gaussian fraction", min=0.0, max=1.0,
                 vary=False),
    ),
    render=_render_jmultiplet,
))


def voigt(x_ppm: np.ndarray, pos_ppm: float, gauss_fwhm_ppm: float,
          lorentz_fwhm_ppm: float, amplitude: float) -> np.ndarray:
    """Peak-normalized TRUE Voigt: the convolution of a Gaussian and a
    Lorentzian (not the pseudo-Voigt sum). Independent Gaussian and Lorentzian
    widths, e.g. Gaussian from disorder + Lorentzian from T2."""
    from scipy.special import voigt_profile

    sigma = max(gauss_fwhm_ppm, 1e-9) * FWHM_TO_SIGMA
    gamma = max(lorentz_fwhm_ppm, 0.0) / 2.0
    y = voigt_profile(x_ppm - pos_ppm, sigma, gamma)
    peak = voigt_profile(0.0, sigma, gamma) or 1.0
    return amplitude * y / peak


def _render_voigt(v: dict, ctx: SimContext) -> np.ndarray:
    return voigt(ctx.x_ppm, v["isotropic_chemical_shift_ppm"],
                 v["gauss_fwhm_ppm"], v["lorentz_fwhm_ppm"], v["amplitude"])


register(Model(
    name="voigt",
    label="Voigt (true)",
    description="True Voigt profile: a Gaussian convolved with a Lorentzian, "
                "with independent widths (Gaussian ← disorder, Lorentzian ← T2). "
                "Unlike Gauss/Lorentz this is a genuine convolution.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm", "peak position"),
        ParamDef("gauss_fwhm_ppm", "gfwhm", 3.0, "ppm", "Gaussian FWHM", min=0.0),
        ParamDef("lorentz_fwhm_ppm", "lfwhm", 3.0, "ppm", "Lorentzian FWHM",
                 min=0.0),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_voigt,
))
