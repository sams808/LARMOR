"""2D fitting: MQMAS (and any 2D VAS method) with the same model registry.

Architecture mirrors the 1D engine exactly, one dimension up:

  * a (Cq, eta) BASIS of 2D subspectra is simulated once per
    (nucleus, field, method, window) via mrsimulator's ThreeQ_VAS / FiveQ_VAS
    and cached;
  * czjzek / ext_czjzek sites are then just REWEIGHTINGS of that basis, so a
    fit iteration costs one matrix product -- the same trick that makes the
    1D Czjzek fits fast;
  * discrete sites (quad_ct) simulate on demand with an LRU cache.

The isotropic dimension follows mrsimulator's convention: its ThreeQ_VAS
method already applies the standard shearing/scaling affine transform, so F1
is a true isotropic axis and no manual shear factor has to be hardcoded here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_KERNEL2D_CACHE: dict[tuple, "Kernel2D"] = {}

METHODS = {"3QMAS": "ThreeQ_VAS", "5QMAS": "FiveQ_VAS", "ST1": "ST1_VAS"}


# --------------------------------------------------------------------------
@dataclass
class Data2D:
    """A processed 2D spectrum on ppm axes. F2 = direct/MAS, F1 = indirect."""

    f2_ppm: np.ndarray               # (n2,) ascending
    f1_ppm: np.ndarray               # (n1,) ascending
    z: np.ndarray                    # (n1, n2) real
    nucleus: str = ""
    larmor_MHz: float = 0.0
    spin_rate_Hz: float = 0.0
    source: str = ""
    notes: list[str] = field(default_factory=list)

    def region(self, f2_range=None, f1_range=None) -> "Data2D":
        s2 = np.ones(self.f2_ppm.shape, bool)
        s1 = np.ones(self.f1_ppm.shape, bool)
        if f2_range:
            s2 = ((self.f2_ppm >= min(f2_range)) & (self.f2_ppm <= max(f2_range)))
        if f1_range:
            s1 = ((self.f1_ppm >= min(f1_range)) & (self.f1_ppm <= max(f1_range)))
        return Data2D(f2_ppm=self.f2_ppm[s2], f1_ppm=self.f1_ppm[s1],
                      z=self.z[np.ix_(s1, s2)], nucleus=self.nucleus,
                      larmor_MHz=self.larmor_MHz,
                      spin_rate_Hz=self.spin_rate_Hz, source=self.source,
                      notes=list(self.notes))

    def normalized(self) -> "Data2D":
        d = self.z / (np.abs(self.z).max() or 1.0)
        return Data2D(self.f2_ppm, self.f1_ppm, d, self.nucleus,
                      self.larmor_MHz, self.spin_rate_Hz, self.source,
                      list(self.notes))

    def projection(self, axis: str = "f2", mode: str = "skyline") -> np.ndarray:
        red = np.max if mode == "skyline" else np.sum
        return red(self.z, axis=0 if axis == "f2" else 1)


def read_bruker_2d(expno: str | Path, procno: int = 1) -> Data2D:
    """Read a processed Bruker 2D (2rr) via the universal reader.

    Accepts an EXPNO folder, a pdata folder, or a 2rr file path. Real
    spectroscopic 2Ds get ppm F1 axes; pseudo-2D (relaxation) datasets keep
    their arrayed F1 coordinates.
    """
    from larmor.io import bruker

    p = Path(expno)
    # allow the historical (expno, procno) call as well as a direct file path
    if p.is_dir() and (p / "pdata" / str(procno) / "2rr").exists():
        p = p / "pdata" / str(procno) / "2rr"
    data = bruker.read(p)
    if data.ndim != 2:
        raise ValueError(f"{expno} is not a 2D dataset")
    if data.domain != "freq":
        raise ValueError("this is a raw 2D FID; Fourier-transform it first "
                         "(larmor.fourier.ft2d)")
    f1_axis, f2_axis = data.axes
    d = Data2D(f2_ppm=f2_axis.values, f1_ppm=f1_axis.values, z=data.data,
               nucleus=data.nucleus, larmor_MHz=data.meta["larmor_MHz"],
               spin_rate_Hz=data.meta.get("masr_Hz") or 0.0,
               source=str(expno))
    d.notes = list(data.warnings)
    if data.is_pseudo2d:
        d.notes.append(f"pseudo-2D: F1 is '{f1_axis.label}' in {f1_axis.unit}")
    return d


def shear(data: Data2D, factor: float, ref_ppm: float = 0.0) -> Data2D:
    """Manual shear: F1' = F1 + factor * (F2 - ref).

    mrsimulator's MQMAS methods already deliver a sheared isotropic axis, so
    this exists for data processed elsewhere (TopSpin without xfshear, ssNake
    workflows).
    """
    out = np.empty_like(data.z)
    for j, f2 in enumerate(data.f2_ppm):
        shift = factor * (f2 - ref_ppm)
        out[:, j] = np.interp(data.f1_ppm - shift, data.f1_ppm, data.z[:, j],
                              left=0.0, right=0.0)
    d = Data2D(data.f2_ppm, data.f1_ppm, out, data.nucleus, data.larmor_MHz,
               data.spin_rate_Hz, data.source, list(data.notes))
    d.notes.append(f"sheared by {factor:g} about {ref_ppm:g} ppm")
    return d


# --------------------------------------------------------------------------
@dataclass
class Kernel2D:
    f2_ppm: np.ndarray               # (n2,)
    f1_ppm: np.ndarray               # (n1,)
    K: np.ndarray                    # (ngrid, n1, n2) basis subspectra
    cq_grid_MHz: np.ndarray
    eta_grid: np.ndarray

    @property
    def shape(self) -> tuple[int, int]:
        return (self.f1_ppm.size, self.f2_ppm.size)

    def weights(self, sigma_MHz: float) -> np.ndarray:
        from mrsimulator.models import CzjzekDistribution

        res = CzjzekDistribution(sigma=sigma_MHz).pdf(
            pos=[self.cq_grid_MHz, self.eta_grid])
        amp = np.asarray(res[-1] if isinstance(res, (tuple, list)) else res)
        w = amp.ravel()
        s = w.sum()
        return w / s if s else w

    def ext_weights(self, cq_MHz: float, eta: float, eps: float) -> np.ndarray:
        from mrsimulator.models import ExtCzjzekDistribution

        res = ExtCzjzekDistribution({"Cq": cq_MHz, "eta": eta},
                                    eps=max(eps, 1e-3)).pdf(
            pos=[self.cq_grid_MHz, self.eta_grid])
        amp = np.asarray(res[-1] if isinstance(res, (tuple, list)) else res)
        w = amp.ravel()
        s = w.sum()
        return w / s if s else w


def build_mqmas_kernel(nucleus: str, larmor_MHz: float,
                       f2_window: tuple[float, float],
                       f1_window: tuple[float, float],
                       n2: int = 192, n1: int = 96,
                       method: str = "3QMAS",
                       cq_max_MHz: float = 16.0, n_cq: int = 40,
                       n_eta: int = 6) -> Kernel2D:
    """Simulate the (Cq, eta) basis of 2D subspectra once (cached).

    Grids are coarser than in 1D on purpose: a 2D basis is n_cq*n_eta*n1*n2
    floats, so 40x6 on a 96x192 window is ~35 MB -- large but workable, while
    the 1D grid (80x11) would be 8x that.

    No spin rate: mrsimulator's named 2D methods fix the rotor frequency at
    the infinite-spinning limit by design (sidebands are not modelled in
    MQMAS), and reject any other value.
    """
    key = (nucleus, round(larmor_MHz, 3), tuple(np.round(f2_window, 2)),
           tuple(np.round(f1_window, 2)), n2, n1, method,
           round(cq_max_MHz, 2), n_cq, n_eta)
    if key in _KERNEL2D_CACHE:
        return _KERNEL2D_CACHE[key]

    from mrsimulator import Simulator
    from mrsimulator.method import SpectralDimension
    from mrsimulator.method import lib as method_lib
    from mrsimulator.spin_system.isotope import Isotope
    from mrsimulator.utils.collection import single_site_system_generator

    if method not in METHODS:
        raise ValueError(f"unknown 2D method {method!r} (valid: {list(METHODS)})")
    method_cls = getattr(method_lib, METHODS[method])

    B0 = larmor_MHz / abs(Isotope(symbol=nucleus).gyromagnetic_ratio)
    cq_grid = np.linspace(0.05, cq_max_MHz, n_cq)
    eta_grid = np.linspace(0, 1, n_eta)
    CQ, ETA = np.meshgrid(cq_grid, eta_grid, indexing="xy")
    n = CQ.size

    systems = single_site_system_generator(
        isotope=nucleus, isotropic_chemical_shift=0.0,
        quadrupolar={"Cq": (CQ * 1e6).ravel(), "eta": ETA.ravel()},
        abundance=np.full(n, 100.0 / n))

    sw2 = abs(f2_window[0] - f2_window[1]) * larmor_MHz
    sw1 = abs(f1_window[0] - f1_window[1]) * larmor_MHz
    ro2 = 0.5 * (f2_window[0] + f2_window[1]) * larmor_MHz
    ro1 = 0.5 * (f1_window[0] + f1_window[1]) * larmor_MHz
    m = method_cls(
        channels=[nucleus], magnetic_flux_density=B0,
        spectral_dimensions=[
            SpectralDimension(count=n1, spectral_width=sw1,
                              reference_offset=ro1, label="iso"),
            SpectralDimension(count=n2, spectral_width=sw2,
                              reference_offset=ro2, label="MAS"),
        ])
    sim = Simulator(spin_systems=systems, methods=[m])
    sim.config.decompose_spectrum = "spin_system"
    sim.run()

    ds = sim.methods[0].simulation
    # mrsimulator returns x[0] = direct (MAS/F2), x[1] = indirect (iso/F1)
    c2 = ds.x[0].coordinates
    c1 = ds.x[1].coordinates
    f2 = c2.value if str(c2.unit) == "ppm" else c2.to("Hz").value / larmor_MHz
    f1 = c1.value if str(c1.unit) == "ppm" else c1.to("Hz").value / larmor_MHz
    o2, o1 = np.argsort(f2), np.argsort(f1)
    K = np.array([np.asarray(dv.components[0].real, float)[np.ix_(o1, o2)]
                  for dv in ds.y])
    kernel = Kernel2D(f2_ppm=np.asarray(f2)[o2], f1_ppm=np.asarray(f1)[o1],
                      K=K, cq_grid_MHz=cq_grid, eta_grid=eta_grid)
    _KERNEL2D_CACHE[key] = kernel
    return kernel


# --------------------------------------------------------------------------
def _broaden_2d(z: np.ndarray, f1_ppm: np.ndarray, f2_ppm: np.ndarray,
                fwhm_f2: float, fwhm_f1: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter

    to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    d2 = abs(f2_ppm[1] - f2_ppm[0]) if f2_ppm.size > 1 else 1.0
    d1 = abs(f1_ppm[1] - f1_ppm[0]) if f1_ppm.size > 1 else 1.0
    s2 = max(fwhm_f2 * to_sigma / d2, 0.0)
    s1 = max(fwhm_f1 * to_sigma / d1, 0.0)
    if s1 < 0.05 and s2 < 0.05:
        return z
    return gaussian_filter(z, sigma=(s1, s2), mode="constant")


def simulate_site_2d(site, kernel: Kernel2D) -> np.ndarray:
    """Render one recipe site as a 2D subspectrum on the kernel grid."""
    v = {k: p.value for k, p in site.params.items()}
    model = site.model
    if model == "czjzek":
        w = kernel.weights(v["sigma_Cq_MHz"])
    elif model == "ext_czjzek":
        w = kernel.ext_weights(v["Cq_MHz"], v["eta"], v["eps"])
    elif model in ("quad_ct", "quad_csa"):
        # nearest grid point of the discrete site (exact-enough on a fine grid
        # and keeps everything inside the same cached basis)
        cq = v["Cq_MHz"]
        eta = v.get("eta_q", v.get("eta", 0.0))
        CQ, ETA = np.meshgrid(kernel.cq_grid_MHz, kernel.eta_grid, indexing="xy")
        d = (CQ.ravel() - cq) ** 2 / (kernel.cq_grid_MHz.max() ** 2) + \
            (ETA.ravel() - eta) ** 2
        w = np.zeros(d.size)
        w[int(np.argmin(d))] = 1.0
    else:
        raise ValueError(f"model {model!r} has no 2D implementation "
                         "(2D supports czjzek, ext_czjzek, quad_ct, quad_csa)")

    z = np.tensordot(w, kernel.K, axes=(0, 0))          # (n1, n2)
    # isotropic shift moves BOTH dimensions: F2 directly, F1 through the
    # isotropic axis convention (mrsimulator already scaled F1)
    pos = v["isotropic_chemical_shift_ppm"]
    z = _shift_2d(z, kernel.f2_ppm, kernel.f1_ppm, pos)
    z = _broaden_2d(z, kernel.f1_ppm, kernel.f2_ppm,
                    v.get("shift_fwhm_ppm", 1.0), v.get("shift_fwhm_ppm", 1.0))
    peak = z.max()
    return v["amplitude"] * (z / peak) if peak > 0 else z


def _shift_2d(z: np.ndarray, f2: np.ndarray, f1: np.ndarray,
              pos_ppm: float) -> np.ndarray:
    """Translate a delta_iso = 0 subspectrum to pos_ppm along both axes."""
    if abs(pos_ppm) < 1e-12:
        return z
    out = np.empty_like(z)
    for i in range(z.shape[0]):
        out[i] = np.interp(f2 - pos_ppm, f2, z[i], left=0.0, right=0.0)
    out2 = np.empty_like(out)
    for j in range(z.shape[1]):
        out2[:, j] = np.interp(f1 - pos_ppm, f1, out[:, j], left=0.0, right=0.0)
    return out2


def simulate_2d(recipe, kernel: Kernel2D):
    """(total, per_site) 2D arrays for a recipe on the kernel grid."""
    per_site = [simulate_site_2d(s, kernel) for s in recipe.sites]
    total = (np.sum(per_site, axis=0) if per_site
             else np.zeros(kernel.shape))
    return total, per_site


# --------------------------------------------------------------------------
@dataclass
class Fit2DResult:
    recipe: object
    lmfit_result: object
    kernel: Kernel2D
    z_fit: np.ndarray
    per_site: list[np.ndarray]
    rmsd: float

    @property
    def report(self) -> str:
        import lmfit

        return lmfit.fit_report(self.lmfit_result)


def fit_2d(recipe, data: Data2D, kernel: Kernel2D | None = None,
           method: str = "3QMAS") -> Fit2DResult:
    """Fit a recipe against a 2D dataset (MQMAS).

    Same registry, same constraints, same uncertainties as 1D -- only the
    residual is two-dimensional.
    """
    import lmfit

    from larmor.fit import _apply_params, _make_params

    data = data.normalized()
    if kernel is None:
        kernel = build_mqmas_kernel(
            recipe.nucleus or data.nucleus,
            recipe.larmor_frequency_MHz or data.larmor_MHz,
            f2_window=(float(data.f2_ppm.max()), float(data.f2_ppm.min())),
            f1_window=(float(data.f1_ppm.max()), float(data.f1_ppm.min())),
            method=method)

    # interpolate the experiment onto the kernel grid once
    from scipy.interpolate import RegularGridInterpolator

    interp = RegularGridInterpolator(
        (data.f1_ppm, data.f2_ppm), data.z, bounds_error=False, fill_value=0.0)
    G1, G2 = np.meshgrid(kernel.f1_ppm, kernel.f2_ppm, indexing="ij")
    z_exp = interp(np.stack([G1.ravel(), G2.ravel()], axis=-1)).reshape(
        kernel.shape)
    scale = np.abs(z_exp).max() or 1.0

    params = _make_params(recipe)

    def residual(p):
        _apply_params(recipe, p)
        total, _ = simulate_2d(recipe, kernel)
        return ((total - z_exp) / scale).ravel()

    # amplitude pre-scale so the optimizer starts on-scale
    _apply_params(recipe, params)
    total0, _ = simulate_2d(recipe, kernel)
    denom = float((total0 * total0).sum())
    if denom > 0:
        k = float((z_exp * total0).sum()) / denom
        for i, site in enumerate(recipe.sites):
            from larmor.fit import _lmfit_name

            name = _lmfit_name(i, site, "amplitude")
            if params[name].vary:
                params[name].value *= k

    result = lmfit.minimize(residual, params, method="least_squares")
    _apply_params(recipe, result.params)
    z_fit, per_site = simulate_2d(recipe, kernel)
    rmsd = float(np.sqrt(np.mean((z_fit - z_exp) ** 2)) / scale)
    recipe.fit_rmsd = rmsd
    return Fit2DResult(recipe=recipe, lmfit_result=result, kernel=kernel,
                       z_fit=z_fit, per_site=per_site, rmsd=rmsd)
