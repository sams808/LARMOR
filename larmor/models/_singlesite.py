"""Shared on-demand single-site simulator with an LRU cache.

quad_ct, quad_first, csa_mas and quad_csa all reduce to: simulate ONE site at
delta_iso = 0 on a window, cache by rounded parameters, then translate /
broaden / scale. Only the Site tensors and the method differ.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np


@lru_cache(maxsize=512)
def simulate_single_site(nucleus: str, larmor_khz: int, mas_hz: int,
                         cq_khz: int, eta_q_milli: int,
                         zeta_cppm: int, eta_cs_milli: int,
                         ct_only: bool, n_ssb: int,
                         x0: float, x1: float, npts: int) -> tuple:
    """Returns (x_ppm tuple, y tuple) for one site at delta_iso = 0.

    cq_khz = 0 disables the quadrupolar tensor; zeta_cppm = 0 disables the
    shielding tensor (zeta is in centi-ppm for a hashable integer key).
    """
    from mrsimulator import Simulator, Site, SpinSystem
    from mrsimulator.method.lib import BlochDecayCTSpectrum, BlochDecaySpectrum
    from mrsimulator.method import SpectralDimension
    from mrsimulator.spin_system.isotope import Isotope
    from mrsimulator.spin_system.tensors import SymmetricTensor

    larmor_MHz = larmor_khz / 1000.0
    B0 = larmor_MHz / abs(Isotope(symbol=nucleus).gyromagnetic_ratio)
    kwargs: dict = {"isotope": nucleus, "isotropic_chemical_shift": 0.0}
    if cq_khz:
        kwargs["quadrupolar"] = SymmetricTensor(
            Cq=cq_khz * 1e3, eta=eta_q_milli / 1000.0)
    if zeta_cppm:
        kwargs["shielding_symmetric"] = SymmetricTensor(
            zeta=zeta_cppm / 100.0, eta=eta_cs_milli / 1000.0)
    site = Site(**kwargs)

    method_cls = BlochDecayCTSpectrum if ct_only else BlochDecaySpectrum
    sw = abs(x1 - x0) * larmor_MHz
    method = method_cls(
        channels=[nucleus], magnetic_flux_density=B0, rotor_frequency=mas_hz,
        spectral_dimensions=[SpectralDimension(
            count=npts, spectral_width=sw,
            reference_offset=(x0 + x1) / 2.0 * larmor_MHz)],
    )
    sim = Simulator(spin_systems=[SpinSystem(sites=[site])], methods=[method])
    sim.config.number_of_sidebands = n_ssb
    sim.run()
    ds = sim.methods[0].simulation
    coords = ds.x[0].coordinates
    x = coords.value if str(coords.unit) == "ppm" \
        else coords.to("Hz").value / larmor_MHz
    y = np.asarray(ds.y[0].components[0].real, dtype=float)
    order = np.argsort(x)
    return tuple(np.asarray(x)[order]), tuple(y[order])


def render_single_site(v: dict, ctx, *, cq_key: str | None = None,
                       eta_q_key: str | None = None,
                       zeta_key: str | None = None,
                       eta_cs_key: str | None = None,
                       ct_only: bool = True, n_ssb: int = 16) -> np.ndarray:
    """Common render path: simulate (cached), shift, broaden, scale."""
    from larmor.models.quadrupolar import _broaden_shift

    x0, x1 = float(ctx.x_ppm[0]), float(ctx.x_ppm[-1])
    npts = min(len(ctx.x_ppm), 2048)
    xs, ys = simulate_single_site(
        ctx.nucleus, int(round(ctx.larmor_MHz * 1000)),
        int(round(ctx.spin_rate_Hz)),
        int(round(v[cq_key] * 1000)) if cq_key else 0,
        int(round(v.get(eta_q_key, 0.0) * 1000)) if eta_q_key else 0,
        int(round(v[zeta_key] * 100)) if zeta_key else 0,
        int(round(v.get(eta_cs_key, 0.0) * 1000)) if eta_cs_key else 0,
        ct_only, n_ssb, round(x0, 2), round(x1, 2), npts)
    xs, ys = np.array(xs), np.array(ys)
    y = _broaden_shift(xs, ys, v["isotropic_chemical_shift_ppm"],
                       v["shift_fwhm_ppm"])
    peak = y.max()
    y = v["amplitude"] * (y / peak) if peak > 0 else y
    return np.interp(ctx.x_ppm, xs, y, left=0.0, right=0.0)
