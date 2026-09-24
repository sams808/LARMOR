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

import threading
from dataclasses import dataclass

from collections import OrderedDict

import numpy as np

from larmor import models as model_registry
from larmor.models.base import SimContext
from larmor.models.analytic import gauss_lor  # back-compat export, see __all__
from larmor.recipe import Recipe, SiteModel

__all__ = [
    "Axis",
    "CQ_MAX_LADDER",
    "CzjzekKernel",
    "KERNEL_CACHE_BUDGET_MB",
    "KERNEL_DISK_CACHE_MB",
    "KERNEL_MIN_SW_HZ",
    "KERNEL_SETTINGS",
    "KERNEL_SPAN_MARGIN",
    "KernelBuildCancelled",
    "build_kernel",
    "cancel_registered",
    "check_cancel",
    "clear_kernel_cache",
    "gauss_lor",
    "grid_restrictable",
    "kernel_axis_ppm",
    "kernel_build_feedback",
    "kernel_cache_info",
    "kernel_cq_max",
    "kernel_window",
    "make_context",
    "needs_kernel",
    "simulate",
    "simulate_site",
    "site_width_margin",
]

#: kernels in most-recently-used order (kept bounded; see _cache_put). A
#: plain dict here only ever grew: one 81Br wideline kernel is ~58 MB
#: (float32) and the Cq ladder makes several per dataset, so a long session
#: quietly held hundreds of MB it would never look at again.
_KERNEL_CACHE: "OrderedDict[tuple, CzjzekKernel]" = OrderedDict()

#: user-tunable 1D Czjzek kernel resolution (dmfit's Computing parameters).
#: Edited via the Computing-parameters dialog; the cache is cleared on change.
#: NOTE: no cq_max here -- the 1D ceiling is AUTOMATIC (kernel_cq_max's
#: ladder follows each model's own width request). A user "Cq max" knob
#: existed for months and did nothing: its only reader built a kernel that
#: make_context immediately discarded.
KERNEL_SETTINGS = {"npts": 2048, "n_cq": 80, "n_eta": 11}

#: kernel-cache byte budget. Roomy enough for a wideline dataset's whole Cq
#: ladder; small enough that a day of mixed datasets cannot pin the machine.
KERNEL_CACHE_BUDGET_MB = 1024.0
#: never evict the N most recent whatever the budget: a fit whose sigma
#: straddles a ladder step alternates between two kernels every iteration,
#: and evicting one of those would mean a full rebuild per iteration.
_CACHE_KEEP_MIN = 4


class KernelBuildCancelled(RuntimeError):
    """Raised out of build_kernel when the registered cancel check fired."""


#: process-wide kernel-build feedback: {"cb": progress(done, total) | None,
#: "cancel": () -> bool | None}. A 23 s wideline build inside a worker
#: thread was indistinguishable from a hang; the build now runs in chunks
#: and reports between them. Set via kernel_build_feedback().
_BUILD_FEEDBACK = {"cb": None, "cancel": None}


#: systems simulated per chunk when feedback hooks are registered: ~7 ms
#: per (Cq, eta) system on the 27Al example, so a chunk is ~0.7 s and a Stop
#: pressed during a build is honoured within about that. The old fixed 8
#: chunks made a 14 080-system build (the 400 MHz ladder step) check the
#: flag every ~12 s.
KERNEL_CHUNK_SYSTEMS = 96

#: one in-flight build per kernel key: the KernelWarmWorker started at load
#: and the SimWorker of the first line dropped seconds later used to build
#: the SAME kernel side by side (the race the warm-up docstring accepted),
#: doubling the CPU time of exactly the wait the user is watching. A second
#: caller now waits on the first build (polling, so a Stop still cancels it)
#: and takes the cached result.
_BUILD_LOCKS: dict = {}
_BUILD_LOCKS_GUARD = threading.Lock()


def _build_lock(key: tuple) -> threading.Lock:
    with _BUILD_LOCKS_GUARD:
        lock = _BUILD_LOCKS.get(key)
        if lock is None:
            lock = _BUILD_LOCKS[key] = threading.Lock()
        return lock


def _n_chunks(n_systems: int) -> int:
    """Chunks for ``n_systems`` with feedback registered: at least the
    historical 8 (progress ticks on small builds), never more systems per
    chunk than KERNEL_CHUNK_SYSTEMS."""
    if n_systems <= 0:
        return 1
    return max(min(8, n_systems), -(-n_systems // KERNEL_CHUNK_SYSTEMS))


def cancel_registered() -> bool:
    """True while a cancel hook is registered on this thread's process."""
    return _BUILD_FEEDBACK.get("cancel") is not None


def check_cancel(where: str = "") -> None:
    """Raise KernelBuildCancelled if the registered cancel hook has fired.

    The kernel build calls it between chunks; larmor.fit calls it between
    site renders so a Stop is honoured within one render even when no kernel
    is being built (a residual evaluation with six wide sites is otherwise
    the unit of latency)."""
    cancel = _BUILD_FEEDBACK.get("cancel")
    if cancel is not None and cancel():
        raise KernelBuildCancelled("cancelled" + (f" {where}" if where else ""))


class kernel_build_feedback:
    """Context manager registering progress/cancel hooks for kernel builds
    on this thread's process (workers each have their own module copy)."""

    def __init__(self, cb=None, cancel=None):
        self._new = {"cb": cb, "cancel": cancel}

    def __enter__(self):
        self._old = dict(_BUILD_FEEDBACK)
        _BUILD_FEEDBACK.update(self._new)
        return self

    def __exit__(self, *exc):
        _BUILD_FEEDBACK.update(self._old)
        return False


#: on-disk kernel cache. A wideline kernel costs ~23 s to simulate and its
#: key is a clean tuple of physical parameters, so re-opening a dataset
#: tomorrow -- or in a fresh WORKER PROCESS (Windows spawn starts every
#: Monte-Carlo/batch worker with an empty in-memory cache) -- should load it
#: in milliseconds instead. Files are versioned by mrsimulator release; the
#: whole folder is bounded by mtime. LARMOR_NO_SESSION disables it (tests).
KERNEL_DISK_CACHE_MB = 2048.0


def _disk_cache_dir():
    import os
    from pathlib import Path

    if os.environ.get("LARMOR_NO_SESSION"):
        return None
    base = os.environ.get("LOCALAPPDATA") or str(Path.home())
    d = Path(base) / "LARMOR" / "kernels"
    try:
        d.mkdir(parents=True, exist_ok=True)
    except OSError:
        return None
    return d


def _disk_key_file(key: tuple):
    import hashlib

    d = _disk_cache_dir()
    if d is None:
        return None
    import mrsimulator

    stamp = f"{key!r}|mrsim={mrsimulator.__version__}"
    return d / (hashlib.sha1(stamp.encode()).hexdigest() + ".npz")


def _disk_load(key: tuple):
    f = _disk_key_file(key)
    if f is None or not f.exists():
        return None
    try:
        with np.load(f, allow_pickle=False) as z:
            kernel = CzjzekKernel(x_ppm=z["x_ppm"], K=z["K"],
                                  cq_grid_MHz=z["cq"], eta_grid=z["eta"])
        f.touch()                                # refresh mtime for the LRU
        return kernel
    except Exception:                            # corrupt/truncated: rebuild
        try:
            f.unlink()
        except OSError:
            pass
        return None


def _disk_store(key: tuple, kernel: "CzjzekKernel") -> None:
    f = _disk_key_file(key)
    if f is None:
        return
    try:
        tmp = f.with_suffix(".tmp.npz")
        np.savez(tmp, x_ppm=kernel.x_ppm, K=kernel.K,
                 cq=kernel.cq_grid_MHz, eta=kernel.eta_grid)
        tmp.replace(f)                           # atomic; racing writers OK
        # bound the folder: drop oldest files past the byte budget
        d = f.parent
        files = sorted(d.glob("*.npz"), key=lambda x: x.stat().st_mtime)
        total = sum(x.stat().st_size for x in files)
        while files and total > KERNEL_DISK_CACHE_MB * 1e6:
            victim = files.pop(0)
            total -= victim.stat().st_size
            victim.unlink()
    except OSError:
        pass


def clear_kernel_cache():
    _KERNEL_CACHE.clear()


def kernel_cache_info() -> dict:
    """{'entries': n, 'mb': total} -- shown in Computing parameters."""
    total = sum(k.K.nbytes + k.x_ppm.nbytes for k in _KERNEL_CACHE.values())
    return {"entries": len(_KERNEL_CACHE), "mb": total / 1e6}


def _cache_put(key: tuple, kernel: "CzjzekKernel") -> None:
    _KERNEL_CACHE[key] = kernel
    _KERNEL_CACHE.move_to_end(key)
    budget = KERNEL_CACHE_BUDGET_MB * 1e6
    while (len(_KERNEL_CACHE) > _CACHE_KEEP_MIN
           and sum(k.K.nbytes for k in _KERNEL_CACHE.values()) > budget):
        _KERNEL_CACHE.popitem(last=False)          # least recently used


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
        _KERNEL_CACHE.move_to_end(key)             # LRU freshness
        return _KERNEL_CACHE[key]
    lock = _build_lock(key)
    while not lock.acquire(timeout=0.2):
        check_cancel("while waiting for a kernel build")   # a Stop still works
    try:
        if key in _KERNEL_CACHE:                   # built while we waited
            _KERNEL_CACHE.move_to_end(key)
            return _KERNEL_CACHE[key]
        kernel = _disk_load(key)
        if kernel is not None:
            _cache_put(key, kernel)
            return kernel
        return _build_kernel_uncached(key, nucleus, larmor_MHz, spin_rate_Hz,
                                      sw_Hz, npts, ref_offset_ppm, cq_max_MHz,
                                      n_cq, n_eta)
    finally:
        lock.release()


def _build_kernel_uncached(key, nucleus, larmor_MHz, spin_rate_Hz, sw_Hz, npts,
                           ref_offset_ppm, cq_max_MHz, n_cq, n_eta) -> CzjzekKernel:
    """The simulation itself (build_kernel holds the key's lock)."""
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
    # the (Cq, eta) systems are simulated in CHUNKS so a 23 s wideline build
    # can report progress and be cancelled between chunks (each system is
    # independent; per-chunk outputs concatenate in grid order). Verified
    # bit-identical to the old monolithic run in tests/test_engine_fit.py.
    cb = _BUILD_FEEDBACK.get("cb")
    cancel = _BUILD_FEEDBACK.get("cancel")
    n_chunks = _n_chunks(len(systems)) if (cb or cancel) else 1
    bounds = np.linspace(0, len(systems), n_chunks + 1).astype(int)
    rows, x = [], None
    for ci in range(n_chunks):
        if cancel is not None and cancel():
            raise KernelBuildCancelled(
                f"kernel build cancelled at chunk {ci}/{n_chunks}")
        chunk = systems[bounds[ci]:bounds[ci + 1]]
        if not chunk:
            continue
        sim = Simulator(spin_systems=chunk, methods=[method])
        sim.config.decompose_spectrum = "spin_system"
        sim.config.number_of_sidebands = 4
        sim.run()
        ds = sim.methods[0].simulation
        coords = ds.x[0].coordinates
        x = (coords.value if str(coords.unit) == "ppm"
             else coords.to("Hz").value / larmor_MHz)
        # float32: the basis is a simulation on a coarse (Cq, eta) grid whose
        # own quantisation dwarfs 1e-7 relative precision; halves both memory
        # AND the weights @ K multiply in every residual evaluation.
        # Example-fit RMSDs verified unchanged to <0.1 % against float64.
        rows.extend(np.asarray(dv.components[0].real, dtype=np.float32)
                    for dv in ds.y)
        if cb is not None:
            cb(ci + 1, n_chunks)
    K = np.array(rows)
    order = np.argsort(x)
    kernel = CzjzekKernel(x_ppm=np.asarray(x)[order], K=K[:, order],
                          cq_grid_MHz=cq_grid, eta_grid=eta_grid)
    _cache_put(key, kernel)
    _disk_store(key, kernel)
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


def kernel_axis_ppm(nucleus: str, larmor_MHz: float, sw_Hz: float, npts: int,
                    ref_offset_ppm: float) -> np.ndarray:
    """The EXACT ppm axis ``build_kernel`` would produce -- without simulating.

    make_context used to build a full (Cq, eta) kernel just to read its
    ``x_ppm`` (23 s and ~115 MB thrown away on a cold wideline dataset, the
    documented "two kernel builds" defect). The axis is pure arithmetic:
    mrsimulator's spectral dimension is ``(arange(n) - n//2) * (sw/n) + ro``
    in Hz, converted to ppm against the isotope's REFERENCE frequency
    (``Isotope.B0_to_ref_freq(B0)``, csdmpy's origin_offset) -- which differs
    from the bare gamma*B0 by the IUPAC reference ratio (0.08 % for 27Al,
    ~0.45 ppm at 550 ppm), so dividing by ``larmor_MHz`` would NOT reproduce
    the kernel axis. Verified bin-exact against a real build in
    tests/test_engine_fit.py.
    """
    from mrsimulator.spin_system.isotope import Isotope

    iso = Isotope(symbol=nucleus)
    B0 = larmor_MHz / abs(iso.gyromagnetic_ratio)
    origin_Hz = float(iso.B0_to_ref_freq(B0))
    hz = ((np.arange(int(npts)) - int(npts) // 2) * (sw_Hz / int(npts))
          + ref_offset_ppm * larmor_MHz)
    x = hz / origin_Hz * 1e6
    return np.sort(x)


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


#: the models that render on the Czjzek kernel at the czjzek axis/resolution
#: rules: czjzek and the two kernel models that must equal it exactly at
#: d = 5 / slope = 0 (one interpolation, not two). ext_czjzek and amorphous
#: staying out is historical, not a design rule.
_KERNEL_AXIS_MODELS = ("czjzek", "czjzek_d", "czjzek_corr")


def needs_kernel(recipe: Recipe) -> bool:
    return any(s.model in _KERNEL_AXIS_MODELS for s in recipe.sites)


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
    "exchange2",                                                   # pointwise (closed-form Bloch-McConnell)
    "czjzek", "czjzek_d", "czjzek_corr", "ext_czjzek", "amorphous",  # own kernel grid
    "spectrum",                                                    # interpolates a reference trace pointwise
})

#: the audited complement: models that MUST simulate on the full experimental
#: grid (render_single_site derives its cached simulation from the axis
#: endpoints and convolves with zero-padded edges), plus "function" (an
#: arbitrary user expression, unknowable). Together with _GRID_RESTRICTABLE
#: this must partition the registry -- tests/test_model_tables.py enforces it.
_GRID_FULL_REQUIRED = frozenset({
    "quad_ct", "quad_first", "quad_csa", "csa_mas", "csa_czjzek", "function",
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
        # the axis alone -- the kernel itself is built (and cached) by the
        # czjzek render when it is actually needed, at the Cq ceiling the
        # model's own sigma asks for
        x = kernel_axis_ppm(recipe.nucleus, recipe.larmor_frequency_MHz,
                            sw_Hz=sw, npts=min(npts, 16384),
                            ref_offset_ppm=ref)
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
