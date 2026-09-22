"""Analytic lineshapes: Gauss/Lorentz (pseudo-Voigt) as in dmfit's Gaus/Lor,
the true Voigt, J-multiplets, empirical sidebands, and the two-site
chemical-exchange lineshape (Gutowsky-Holm / McConnell, closed form)."""
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


def _render_sidebands(v: dict, ctx: SimContext) -> np.ndarray:
    """Empirical spinning-sideband manifold (dmfit 'ss band'): a centre band
    plus sidebands at ±k·νrot with a geometric intensity ratio — for modelling
    an observed sideband pattern without a full CSA fit."""
    pos = v["isotropic_chemical_shift_ppm"]
    fwhm = v["shift_fwhm_ppm"]; amp = v["amplitude"]; gl = v.get("gl", 1.0)
    r = float(v.get("ssb_ratio", 0.3))
    n = int(round(v.get("n_ssb", 4)))
    nu = ctx.spin_rate_Hz; lar = ctx.larmor_MHz or 1.0
    y = amp * gauss_lor(ctx.x_ppm, pos, fwhm, 1.0, gl)
    if nu > 0:
        spacing = nu / lar                          # ppm between sidebands
        for k in range(1, n + 1):
            w = amp * r ** k
            y += w * gauss_lor(ctx.x_ppm, pos + k * spacing, fwhm, 1.0, gl)
            y += w * gauss_lor(ctx.x_ppm, pos - k * spacing, fwhm, 1.0, gl)
    return y


register(Model(
    name="sidebands",
    label="Spinning sidebands",
    description="Empirical MAS sideband manifold (dmfit 'ss band'): a centre band "
                "plus sidebands at ±k·νrot with a geometric intensity ratio. Set "
                "the MAS rate in the experiment parameters.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm", "centre band"),
        ParamDef("shift_fwhm_ppm", "fwhm", 2.0, "ppm", "band FWHM", min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "centre-band height", min=0.0),
        ParamDef("ssb_ratio", "r", 0.3, "", "intensity ratio between successive "
                 "sidebands", min=0.0, max=1.0),
        ParamDef("n_ssb", "nssb", 4.0, "", "sidebands each side (integer)",
                 min=0.0, max=32.0, vary=False),
        ParamDef("gl", "gl", 1.0, "", "Gaussian fraction", min=0.0, max=1.0,
                 vary=False),
    ),
    render=_render_sidebands,
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


def two_site_exchange(nu_hz: np.ndarray, nu_a_hz: float, nu_b_hz: float,
                      p_a: float, k_ex: float, fwhm_hz: float) -> np.ndarray:
    """Steady-state Bloch-McConnell absorption for two sites A <-> B in
    chemical exchange (Gutowsky & Holm 1956; McConnell 1958), unnormalised.

    Populations p_A and p_B = 1 - p_A, one intrinsic Lorentzian FWHM
    ``fwhm_hz`` for both sites (1/(pi T2)) and the exchange rate
    k_ex = k_AB + k_BA = 1/tau (s^-1) under detailed balance
    (k_AB = p_B k_ex, k_BA = p_A k_ex). Solving the 2x2 steady-state system
    and adding the two magnetisations gives the cancellation-free form

        M(nu) = [p_A beta_B + p_B beta_A + k_ex]
                / [beta_A beta_B + k_ex (p_A beta_A + p_B beta_B)]
        beta_X = pi w + 2 pi i (nu - nu_X)

    and the absorption is Re M >= 0. Limits: k_ex -> 0 gives
    p_A/beta_A + p_B/beta_B (two Lorentzians, areas p_A : p_B); k_ex -> inf
    gives 1/(p_A beta_A + p_B beta_B), one Lorentzian of FWHM w at the
    population-weighted mean frequency. Equal populations coalesce at
    k_ex = sqrt(2) pi |nu_A - nu_B|; the fast-exchange residual FWHM is
    w + 4 pi p_A p_B (nu_A - nu_B)^2 / k_ex. The area does not depend on k_ex.
    """
    nu = np.asarray(nu_hz, float)
    p_a = float(p_a)
    p_b = 1.0 - p_a
    w = max(float(fwhm_hz), 1e-9)
    k = max(float(k_ex), 0.0)
    beta_a = np.pi * w + 2j * np.pi * (nu - float(nu_a_hz))
    beta_b = np.pi * w + 2j * np.pi * (nu - float(nu_b_hz))
    num = p_a * beta_b + p_b * beta_a + k
    den = beta_a * beta_b + k * (p_a * beta_a + p_b * beta_b)
    return np.real(num / den)


def _render_exchange2(v: dict, ctx: SimContext) -> np.ndarray:
    lar = ctx.larmor_MHz or 1.0                 # ppm -> Hz (jmultiplet precedent)
    pos = float(v["isotropic_chemical_shift_ppm"])
    split = float(v.get("split_ppm", 5.0))
    pa = float(np.clip(v.get("pop_a", 0.5), 0.0, 1.0))
    # pos is the population-weighted mean: p_A d_A + p_B d_B = pos, and
    # d_A - d_B = split with A the higher-ppm site
    d_a = pos + (1.0 - pa) * split
    d_b = pos - pa * split
    y = two_site_exchange(ctx.x_ppm * lar, d_a * lar, d_b * lar, pa,
                          float(v.get("k_ex_hz", 100.0)),
                          float(v.get("lorentz_fwhm_ppm", 1.0)) * lar)
    peak = float(y.max())
    return v["amplitude"] * y / peak if peak > 0 else y


register(Model(
    name="exchange2",
    label="Two-site exchange  (Bloch–McConnell)",
    description="Two isotropic sites A/B in chemical exchange (Gutowsky–Holm "
                "1956 / McConnell 1958): populations p_A, 1−p_A, separation "
                "Δδ = δ_A − δ_B, exchange rate k_ex = k_AB + k_BA (s⁻¹) and one "
                "intrinsic (T2) Lorentzian width. pos is the population-weighted "
                "mean shift. k→0: two lines; k→∞: one line at pos.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "population-weighted mean shift (the fast-exchange position)"),
        ParamDef("split_ppm", "dd", 5.0, "ppm",
                 "δ_A − δ_B (A is the higher-ppm site)", min=0.0),
        ParamDef("pop_a", "pa", 0.5, "", "population of A (p_B = 1 − p_A)",
                 min=0.01, max=0.99),
        ParamDef("k_ex_hz", "kex", 100.0, "Hz",
                 "exchange rate k_ex = k_AB + k_BA = 1/τ (s⁻¹)",
                 min=0.0, max=1e8),
        ParamDef("lorentz_fwhm_ppm", "lfwhm", 1.0, "ppm",
                 "intrinsic Lorentzian FWHM of each site (1/πT2)", min=0.01),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_exchange2,
))
