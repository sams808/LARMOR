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

#: the Czjzek kernel's own simulated window. It used to be a FIXED 150 kHz,
#: which is a different width in ppm for every nucleus -- 1152 ppm for 27Al at
#: 130 MHz (fine for aluminosilicates, which is why it went unnoticed) but only
#: 696 ppm for 81Br at 216 MHz. Everything outside it is interpolated to ZERO,
#: so a wide-line site simply vanished: a real 81Br site at 617 ppm rendered as
#: all zeros, and ext_czjzek could not produce a pattern wider than ~694 ppm no
#: matter what Cq was asked for.
KERNEL_MIN_SW_HZ = 150000.0
#: how much wider than the data the kernel is built, so a pattern can spill
#: past the fit window without being clipped
KERNEL_SPAN_MARGIN = 1.25
#: the kernel's (Cq, eta) GRID ceiling. Fixing the spectral window alone was
#: not enough: the grid stopped at 25 MHz, so a Czjzek distribution whose
#: weight lies above that could not be represented and the pattern saturated
#: (81Br: 1013 ppm however large sigma got, against 1488 ppm of real data).
#: Quantised to a ladder so the cache holds a handful of kernels rather than
#: one per optimiser step.
CQ_MAX_LADDER = (25.0, 50.0, 100.0, 200.0, 400.0)


def kernel_cq_max(needed_MHz: float) -> float:
    """Smallest ladder step covering ``needed_MHz`` (the largest Cq the model
    puts real weight on)."""
    for step in CQ_MAX_LADDER:
        if needed_MHz <= step:
            return step
    return CQ_MAX_LADDER[-1]


def kernel_window(x_ppm, larmor_MHz: float) -> tuple[float, float]:
    """(sw_Hz, ref_offset_ppm) for a kernel that COVERS ``x_ppm``.

    Never narrower than the historical 150 kHz, so every axis that already
    fitted keeps its kernel; only a wider request extends it.
    """
    if x_ppm is None or larmor_MHz <= 0:
        return KERNEL_MIN_SW_HZ, 30.0
    x = np.asarray(x_ppm, float)
    if x.size < 2:
        return KERNEL_MIN_SW_HZ, 30.0
    lo, hi = float(np.min(x)), float(np.max(x))
    need = (hi - lo) * KERNEL_SPAN_MARGIN * larmor_MHz
    if need <= KERNEL_MIN_SW_HZ:
        return KERNEL_MIN_SW_HZ, 30.0      # unchanged: the historical window
    return float(need), 0.5 * (lo + hi)


def needs_kernel(recipe: Recipe) -> bool:
    return any(s.model == "czjzek" for s in recipe.sites)


#: models safe to simulate on a grid RESTRICTED to the fit window (+ margin)
#: rather than the full experimental axis. Deliberately an ALLOWLIST, not a
#: denylist: a new/unaudited model defaults to the full grid (always correct,
#: just not optimized) until someone verifies it belongs here. A model
#: qualifies if either (a) its render() is a closed-form formula evaluated
#: pointwise at each x (gauss_lor, voigt, ...) -- restricting the grid can't
#: change any value, since points don't interact -- or (b) it simulates on its
#: OWN independent grid (the Czjzek kernel family, engine.build_kernel) and
#: only interpolates the result onto ctx.x_ppm at the very end -- restricting
#: ctx.x_ppm only shrinks that final, always-safe downsampling step.
#:
#: Explicitly EXCLUDED: quad_ct/quad_first/quad_csa/csa_mas/csa_czjzek (all via
#: models/_singlesite.py's render_single_site) build their OWN cached
#: simulation directly from ctx.x_ppm's first/last VALUE, and broaden it with a
#: real convolution (gaussian_filter1d, mode="constant" == zero-padded edges) —
#: restricting the grid there would truncate the sideband/satellite manifold
#: and change the result, not just its cost. Also excluded: "function" (an
#: arbitrary user expression -- unknowable in general).
_GRID_RESTRICTABLE = frozenset({
    "gauss_lor", "gl_norm", "jmultiplet", "sidebands", "voigt",   # pointwise
    "czjzek", "ext_czjzek", "amorphous",                          # own kernel grid
    "spectrum",                                                    # interpolates a reference trace pointwise
})


def grid_restrictable(recipe: Recipe) -> bool:
    """True if every site's model tolerates simulating on a window-restricted
    grid instead of the full experimental axis (see _GRID_RESTRICTABLE)."""
    return all(s.model in _GRID_RESTRICTABLE for s in recipe.sites)


#: parameter names treated as a lineshape "width" when estimating how far a
#: site's visible extent reaches beyond its center (shared by figures.py's
#: data-less-fit preview and fit.py's windowed simulation grid, so the two
#: never drift apart).
_WIDTH_PARAM_NAMES = ("shift_fwhm_ppm", "line_fwhm_ppm", "lorentz_fwhm_ppm",
                     "gauss_fwhm_ppm", "fwhm")


def site_width_margin(sites, default: float = 10.0, factor: float = 6.0) -> float:
    """A generous margin (ppm) beyond a site's center that a lineshape needs
    room to be simulated in: factor x the widest declared width among the
    sites (dmfit/ssNake lineshapes are negligible beyond a few widths)."""
    spans = [default]
    for s in sites:
        for wn in _WIDTH_PARAM_NAMES:
            if wn in s.params:
                spans.append(abs(float(s.params[wn].value)))
    return max(spans) * factor


def make_context(recipe: Recipe, exp_ppm: np.ndarray | None = None) -> SimContext:
    """Build the simulation context; picks the axis a recipe should render on."""
    if needs_kernel(recipe):
        sw, ref = kernel_window(exp_ppm, recipe.larmor_frequency_MHz)
        # keep the RESOLUTION when the window is widened, or a broad-line
        # dataset would be simulated on a coarser grid than its own data
        npts = int(KERNEL_SETTINGS["npts"]
                   * max(1.0, sw / KERNEL_MIN_SW_HZ))
        kernel = build_kernel(recipe.nucleus, recipe.larmor_frequency_MHz,
                              recipe.spin_rate_Hz, sw_Hz=sw,
                              npts=min(npts, 16384), ref_offset_ppm=ref,
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
    if site.model == "function":
        return _render_function(site, ctx)
    values = {k: v.value for k, v in site.params.items()}
    return model_registry.get(site.model).render(values, ctx)


_SAFE_FUNCS = ("sin", "cos", "tan", "exp", "log", "log10", "sqrt", "abs",
               "tanh", "arctan", "sign", "sinc")


def _render_function(site, ctx) -> np.ndarray:
    """Evaluate a user y(x; a,b,c,d) expression (ssNake Function fit) on the ppm
    axis, scaled by amplitude. Restricted namespace (numpy funcs + the params)."""
    expr = getattr(site, "func", None)
    if not expr:
        return np.zeros_like(ctx.x_ppm)
    ns = {"x": ctx.x_ppm, "pi": np.pi, "np": np}
    ns.update({fn: getattr(np, fn) for fn in _SAFE_FUNCS})
    ns.update({k: v.value for k, v in site.params.items()})
    try:
        y = eval(expr, {"__builtins__": {}}, ns)  # noqa: S307 (trusted local user)
    except Exception:
        return np.zeros_like(ctx.x_ppm)
    amp = site.params["amplitude"].value if "amplitude" in site.params else 1.0
    return amp * (np.asarray(y, float) * np.ones_like(ctx.x_ppm))


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
