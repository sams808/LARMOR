"""Fast lineshape engine: precomputed Czjzek kernel + cheap reweighting.

The expensive part of a Czjzek fit is simulating the quadrupolar lineshape for
every (Cq, eta) grid point. That basis does not depend on the fit parameters:
sigma only reweights the grid, the isotropic shift only translates the
spectrum, and the shift-distribution width is a convolution. So the kernel is
simulated ONCE per (nucleus, field, spin rate, spectral window) via
mrsimulator, and every fit iteration afterwards costs milliseconds -- the same
factorization mrinversion uses.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d

from larmor.recipe import Recipe, SiteModel

_KERNEL_CACHE: dict[tuple, "CzjzekKernel"] = {}

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


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

    def spectrum(self, sigma_MHz: float, pos_ppm: float,
                 shift_fwhm_ppm: float, amplitude: float) -> np.ndarray:
        """One Czjzek site's lineshape on self.x_ppm, peak-normalized basis."""
        y = self.weights(sigma_MHz) @ self.K
        # kernel is simulated at delta_iso = 0: translate by pos_ppm
        y = np.interp(self.x_ppm - pos_ppm, self.x_ppm, y, left=0.0, right=0.0)
        dppm = abs(self.x_ppm[1] - self.x_ppm[0])
        sigma_pts = shift_fwhm_ppm * FWHM_TO_SIGMA / dppm
        if sigma_pts > 0.05:
            y = gaussian_filter1d(y, sigma_pts, mode="constant")
        peak = y.max()
        return amplitude * (y / peak) if peak > 0 else y


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
def gauss_lor(x_ppm: np.ndarray, pos_ppm: float, fwhm_ppm: float,
              amplitude: float, gl: float) -> np.ndarray:
    """Analytic pseudo-Voigt, dmfit-style: y = gl*Gaussian + (1-gl)*Lorentzian.

    Both components are peak-normalized so `amplitude` is the peak height.
    """
    dx = x_ppm - pos_ppm
    sig = max(fwhm_ppm, 1e-6) * FWHM_TO_SIGMA
    g = np.exp(-0.5 * (dx / sig) ** 2)
    hwhm = max(fwhm_ppm, 1e-6) / 2.0
    l = 1.0 / (1.0 + (dx / hwhm) ** 2)
    return amplitude * (gl * g + (1.0 - gl) * l)


def simulate_site(site: SiteModel, kernel: CzjzekKernel) -> np.ndarray:
    p = {k: v.value for k, v in site.params.items()}
    if site.model == "czjzek":
        return kernel.spectrum(
            sigma_MHz=p["sigma_Cq_MHz"],
            pos_ppm=p["isotropic_chemical_shift_ppm"],
            shift_fwhm_ppm=p["shift_fwhm_ppm"],
            amplitude=p["amplitude"],
        )
    if site.model == "gauss_lor":
        return gauss_lor(
            kernel.x_ppm, p["isotropic_chemical_shift_ppm"],
            p["shift_fwhm_ppm"], p["amplitude"], p.get("gl", 1.0),
        )
    raise ValueError(f"unknown site model {site.model!r}")


def simulate(recipe: Recipe, kernel: CzjzekKernel | None = None,
             ) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Simulate a recipe. Returns (x_ppm, total, per_site)."""
    if kernel is None:
        kernel = build_kernel(recipe.nucleus, recipe.larmor_frequency_MHz,
                              recipe.spin_rate_Hz)
    per_site = [simulate_site(s, kernel) for s in recipe.sites]
    total = np.sum(per_site, axis=0) if per_site else np.zeros_like(kernel.x_ppm)
    return kernel.x_ppm, total, per_site
