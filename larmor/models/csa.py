"""Chemical-shift-anisotropy model: spin-1/2 powder pattern with physical
spinning sidebands (replaces dmfit's ad-hoc 'ss band' lines)."""
from __future__ import annotations

import numpy as np

from larmor.models.base import Model, ParamDef, SimContext, register


def _n_ssb_for(zeta_ppm: float, ctx: SimContext, floor: int = 8,
               ceiling: int = 256) -> int:
    """Sidebands mrsimulator must carry so none is truncated: sized from the
    static span |zeta| nu0 / nu_rot (Herzfeld-Berger), rounded up to a power
    of two (the simulator's efficient sizes). A fixed 32 silently clipped the
    manifold of a slowly spun, strongly anisotropic site."""
    from larmor.herzfeld_berger import n_sidebands_needed

    need = n_sidebands_needed(zeta_ppm, ctx.larmor_MHz, ctx.spin_rate_Hz)
    if need <= 0:
        return floor
    n = 1 << int(np.ceil(np.log2(max(2 * need, floor))))
    return int(min(max(n, floor), ceiling))


def _render_csa(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor.models._singlesite import render_single_site

    return render_single_site(v, ctx, zeta_key="zeta_ppm", eta_cs_key="eta",
                              ct_only=False,
                              n_ssb=_n_ssb_for(v["zeta_ppm"], ctx))


def _render_csa_czjzek(v: dict, ctx: SimContext) -> np.ndarray:
    """Disordered CSA: the shielding-anisotropy powder pattern averaged over a
    Gaussian distribution of ζ (width σζ) — the CSA analogue of a Czjzek quad
    distribution, for glassy/amorphous shielding."""
    from larmor.models._singlesite import render_single_site

    z0 = v["isotropic_chemical_shift_ppm"] * 0 + v["zeta_ppm"]
    sig = float(v.get("sigma_zeta_ppm", 20.0))
    base = dict(v); base["amplitude"] = 1.0

    def one(zeta):
        return render_single_site({**base, "zeta_ppm": zeta}, ctx,
                                  zeta_key="zeta_ppm", eta_cs_key="eta",
                                  ct_only=False, n_ssb=16)

    if sig <= 0.5:
        y = one(z0)
    else:
        zs = np.linspace(z0 - 2 * sig, z0 + 2 * sig, 5)
        w = np.exp(-0.5 * ((zs - z0) / sig) ** 2); w /= w.sum()
        y = sum(wi * one(z) for z, wi in zip(zs, w))
    peak = float(np.max(y))
    return v["amplitude"] * (y / peak) if peak > 0 else y


register(Model(
    name="csa_czjzek",
    label="CSA distribution (disordered)",
    description="Disordered shielding anisotropy: the CSA powder pattern averaged "
                "over a Gaussian distribution of ζ (width σζ) — the CSA analogue "
                "of the quadrupolar Czjzek, for glasses.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("zeta_ppm", "zeta", 50.0, "ppm", "mean shielding anisotropy",
                 min=-1000.0, max=1000.0),
        ParamDef("sigma_zeta_ppm", "sz", 20.0, "ppm", "ζ distribution width",
                 min=0.0, max=500.0),
        ParamDef("eta", "eta", 0.3, "", "shielding asymmetry", min=0.0, max=1.0),
        ParamDef("shift_fwhm_ppm", "fwhm", 2.0, "ppm", "Gaussian broadening",
                 min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_csa_czjzek,
))


register(Model(
    name="csa_mas",
    label="CSA powder (MAS/static)",
    description="Shielding-anisotropy powder pattern with physical spinning "
                "sidebands whose intensities follow Herzfeld & Berger (1980) "
                "from ζ and η (spin rate 0 = static). Replaces manual 'ss "
                "band' lines; Tools > Herzfeld–Berger sideband analysis reads "
                "ζ, η off the measured manifold to seed it.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("zeta_ppm", "zeta", 50.0, "ppm",
                 "shielding anisotropy (Haeberlen zeta)", min=-1000.0, max=1000.0),
        ParamDef("eta", "eta", 0.3, "", "shielding asymmetry", min=0.0, max=1.0),
        ParamDef("shift_fwhm_ppm", "fwhm", 2.0, "ppm", "Gaussian broadening",
                 min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_csa,
))
