"""Simulation engine: shared kernels + registry-dispatched site rendering.

The expensive part of a Czjzek fit is simulating the quadrupolar lineshape for
every (Cq, eta) grid point. That basis does not depend on the fit parameters,
so it is simulated ONCE per (nucleus, field, spin rate, window) and cached;
every fit iteration afterwards is a cheap reweighting. Discrete models
(quad_ct, csa_mas) simulate on demand with parameter-level LRU caches instead.

Site rendering itself is dispatched through larmor.models.REGISTRY, so new
models plug in without touching this module or the fit engine.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from larmor import models as model_registry
from larmor.models.base import SimContext
from larmor.models.analytic import gauss_lor  # noqa: F401  (back-compat export)
from larmor.recipe import Recipe, SiteModel

_KERNEL_CACHE: dict[tuple, "CzjzekKernel"] = {}

#: user-tunable 1D Czjzek kernel resolution (dmfit's Computing parameters).
#: Edited via the Computing-parameters dialog; the cache is cleared on change.
KERNEL_SETTINGS = {"npts": 2048, "cq_max_MHz": 25.0, "n_cq": 80, "n_eta": 11}


def clear_kernel_cache():
    _KERNEL_CACHE.clear()


@dataclass
class Axis:
    """Bare ppm axis for recipes with no kernel-based site."""

    x_ppm: np.ndarray


@dataclass
class CzjzekKernel:
    x_ppm: np.ndarray            # ascending ppm axis, shape (npts,)
    K: np.ndarray                # basis subspectra, shape (ngrid, npts)
    cq_grid_MHz: np.ndarray
    eta_grid: np.ndarray

    def weights(self, sigma_MHz: float) -> np.ndarray:
        from mrsimulator.models import CzjzekDistribution

        res = CzjzekDistribution(sigma=sigma_MHz).pdf(
            pos=[self.cq_grid_MHz, self.eta_grid])
        amp = np.asarray(res[-1] if isinstance(res, (tuple, list)) else res)
        w = amp.ravel()
        return w / w.sum()


def build_kernel(nucleus: str, larmor_MHz: float, spin_rate_Hz: float,
                 sw_Hz: float = 150000.0, npts: int = 2048,
                 ref_offset_ppm: float = 30.0,
                 cq_max_MHz: float = 25.0, n_cq: int = 80, n_eta: int = 11,
                 ) -> CzjzekKernel:
    """Simulate the (Cq, eta) basis once with mrsimulator (cached per process)."""
    key = (nucleus, round(larmor_MHz, 3), round(spin_rate_Hz), round(sw_Hz),
           npts, round(ref_offset_ppm, 1), round(cq_max_MHz, 1), n_cq, n_eta)
    if key in _KERNEL_CACHE:
        return _KERNEL_CACHE[key]

    from mrsimulator import Simulator
    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.method import SpectralDimension
    from mrsimulator.spin_system.isotope import Isotope
    from mrsimulator.utils.collection import single_site_system_generator

    B0 = larmor_MHz / abs(Isotope(symbol=nucleus).gyromagnetic_ratio)
    cq_grid = np.linspace(0.05, cq_max_MHz, n_cq)
    eta_grid = np.linspace(0, 1, n_eta)
    CQ, ETA = np.meshgrid(cq_grid, eta_grid, indexing="xy")
    n = CQ.size

    systems = single_site_system_generator(
        isotope=nucleus,
        isotropic_chemical_shift=0.0,
        quadrupolar={"Cq": (CQ * 1e6).ravel(), "eta": ETA.ravel()},
        abundance=np.full(n, 100.0 / n),
    )
    method = BlochDecayCTSpectrum(
        channels=[nucleus],
        magnetic_flux_density=B0,
        rotor_frequency=spin_rate_Hz,
        spectral_dimensions=[SpectralDimension(
            count=npts, spectral_width=sw_Hz,
            reference_offset=ref_offset_ppm * larmor_MHz)],
    )
    sim = Simulator(spin_systems=systems, methods=[method])
    sim.config.decompose_spectrum = "spin_system"
    sim.config.number_of_sidebands = 4
    sim.run()

    ds = sim.methods[0].simulation
    coords = ds.x[0].coordinates
    x = coords.value if str(coords.unit) == "ppm" else coords.to("Hz").value / larmor_MHz
    K = np.array([np.asarray(dv.components[0].real, dtype=float) for dv in ds.y])
    order = np.argsort(x)
    kernel = CzjzekKernel(x_ppm=np.asarray(x)[order], K=K[:, order],
                          cq_grid_MHz=cq_grid, eta_grid=eta_grid)
    _KERNEL_CACHE[key] = kernel
    return kernel


# --------------------------------------------------------------------------

def needs_kernel(recipe: Recipe) -> bool:
    return any(s.model == "czjzek" for s in recipe.sites)


def make_context(recipe: Recipe, exp_ppm: np.ndarray | None = None) -> SimContext:
    """Build the simulation context; picks the axis a recipe should render on."""
    if needs_kernel(recipe):
        kernel = build_kernel(recipe.nucleus, recipe.larmor_frequency_MHz,
                              recipe.spin_rate_Hz,
                              npts=KERNEL_SETTINGS["npts"],
                              cq_max_MHz=KERNEL_SETTINGS["cq_max_MHz"],
                              n_cq=KERNEL_SETTINGS["n_cq"],
                              n_eta=KERNEL_SETTINGS["n_eta"])
        x = kernel.x_ppm
    elif exp_ppm is not None:
        x = np.asarray(exp_ppm)[np.argsort(exp_ppm)]
    else:
        x = np.linspace(-300, 300, 2048)
    return SimContext(nucleus=recipe.nucleus,
                      larmor_MHz=recipe.larmor_frequency_MHz,
                      spin_rate_Hz=recipe.spin_rate_Hz, x_ppm=x)


def simulate_site(site: SiteModel, ctx) -> np.ndarray:
    """Render one site on the context axis. Accepts a SimContext (preferred)
    or, for backward compatibility, a CzjzekKernel/Axis."""
    if isinstance(ctx, (CzjzekKernel, Axis)):
        ctx = SimContext(nucleus="27Al", larmor_MHz=0.0, spin_rate_Hz=0.0,
                         x_ppm=ctx.x_ppm) if isinstance(ctx, Axis) else _ctx_from_kernel(ctx)
    if site.model == "spectrum":
        return _render_spectrum(site, ctx)
    values = {k: v.value for k, v in site.params.items()}
    return model_registry.get(site.model).render(values, ctx)


def _render_spectrum(site, ctx) -> np.ndarray:
    """Render an external-spectrum component: its reference trace (unit peak),
    interpolated onto the fit axis, rigidly shifted, and scaled by amplitude."""
    ref = getattr(site, "ref", None) or {}
    rp = np.asarray(ref.get("ppm", []), float)
    ra = np.asarray(ref.get("amp", []), float)
    if rp.size < 2 or ra.size != rp.size:
        return np.zeros_like(ctx.x_ppm)
    amp = site.params["amplitude"].value
    shift = site.params["shift_ppm"].value if "shift_ppm" in site.params else 0.0
    order = np.argsort(rp)
    y = np.interp(ctx.x_ppm, rp[order] + shift, ra[order], left=0.0, right=0.0)
    return amp * y


def _ctx_from_kernel(kernel: CzjzekKernel) -> SimContext:
    # legacy path: infer nothing, just carry the axis; czjzek render rebuilds
    # its kernel from the cache so this only needs the axis to be right
    return SimContext(nucleus="27Al", larmor_MHz=0.0, spin_rate_Hz=0.0,
                      x_ppm=kernel.x_ppm)


def simulate(recipe: Recipe, kernel=None, exp_ppm: np.ndarray | None = None,
             ) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Simulate a recipe. Returns (x_ppm, total, per_site)."""
    if kernel is not None and isinstance(kernel, (CzjzekKernel, Axis)):
        ctx = SimContext(nucleus=recipe.nucleus,
                         larmor_MHz=recipe.larmor_frequency_MHz,
                         spin_rate_Hz=recipe.spin_rate_Hz, x_ppm=kernel.x_ppm)
    else:
        ctx = make_context(recipe, exp_ppm=exp_ppm)
    per_site = [simulate_site(s, ctx) for s in recipe.sites]
    total = np.sum(per_site, axis=0) if per_site else np.zeros_like(ctx.x_ppm)
    return ctx.x_ppm, total, per_site
