"""Chemical-shift-anisotropy model: spin-1/2 powder pattern with physical
spinning sidebands (replaces dmfit's ad-hoc 'ss band' lines)."""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from larmor.models.base import Model, ParamDef, SimContext, register
from larmor.models.quadrupolar import _broaden_shift


@lru_cache(maxsize=256)
def _csa_shape(nucleus: str, larmor_khz: int, mas_hz: int,
               zeta_cppm: int, eta_milli: int,
               x0: float, x1: float, npts: int) -> tuple:
    from mrsimulator import Simulator, Site, SpinSystem
    from mrsimulator.method.lib import BlochDecaySpectrum
    from mrsimulator.method import SpectralDimension
    from mrsimulator.spin_system.isotope import Isotope
    from mrsimulator.spin_system.tensors import SymmetricTensor

    larmor_MHz = larmor_khz / 1000.0
    B0 = larmor_MHz / abs(Isotope(symbol=nucleus).gyromagnetic_ratio)
    site = Site(isotope=nucleus, isotropic_chemical_shift=0.0,
                shielding_symmetric=SymmetricTensor(
                    zeta=zeta_cppm / 100.0, eta=eta_milli / 1000.0))
    sw = abs(x1 - x0) * larmor_MHz
    method = BlochDecaySpectrum(
        channels=[nucleus], magnetic_flux_density=B0, rotor_frequency=mas_hz,
        spectral_dimensions=[SpectralDimension(
            count=npts, spectral_width=sw,
            reference_offset=(x0 + x1) / 2.0 * larmor_MHz)],
    )
    sim = Simulator(spin_systems=[SpinSystem(sites=[site])], methods=[method])
    sim.config.number_of_sidebands = 32
    sim.run()
    ds = sim.methods[0].simulation
    coords = ds.x[0].coordinates
    x = coords.value if str(coords.unit) == "ppm" else coords.to("Hz").value / larmor_MHz
    y = np.asarray(ds.y[0].components[0].real, dtype=float)
    order = np.argsort(x)
    return tuple(x[order]), tuple(y[order])


def _render_csa(v: dict, ctx: SimContext) -> np.ndarray:
    x0, x1 = float(ctx.x_ppm[0]), float(ctx.x_ppm[-1])
    npts = min(len(ctx.x_ppm), 2048)
    xs, ys = _csa_shape(
        ctx.nucleus, int(round(ctx.larmor_MHz * 1000)),
        int(round(ctx.spin_rate_Hz)),
        int(round(v["zeta_ppm"] * 100)), int(round(v["eta"] * 1000)),
        round(x0, 2), round(x1, 2), npts)
    xs, ys = np.array(xs), np.array(ys)
    y = _broaden_shift(xs, ys, v["isotropic_chemical_shift_ppm"],
                       v["shift_fwhm_ppm"])
    peak = y.max()
    y = v["amplitude"] * (y / peak) if peak > 0 else y
    return np.interp(ctx.x_ppm, xs, y, left=0.0, right=0.0)


register(Model(
    name="csa_mas",
    label="CSA powder (MAS/static)",
    description="Shielding-anisotropy powder pattern with physical spinning "
                "sidebands (spin rate 0 = static). Replaces manual 'ss band' "
                "lines.",
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
