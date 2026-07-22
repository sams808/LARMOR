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

#: user-tunable MQMAS kernel resolution (dmfit Computing parameters, nuQ block).
MQMAS_SETTINGS = {"n2": 192, "n1": 96, "n_cq": 40, "n_eta": 6, "cq_max_MHz": 16.0}


def clear_kernel_cache():
    _KERNEL2D_CACHE.clear()


# --------------------------------------------------------------------------
@dataclass
class Data2D:
    """A processed 2D spectrum on ppm axes. F2 = direct/MAS, F1 = indirect.

    ``z`` is the real-real quadrant. When the imaginary quadrants ``ri`` (imag
    along F2), ``ir`` (imag along F1) and ``ii`` are present, phase correction
    is exact (true hypercomplex) instead of Hilbert-reconstructed.
    """

    f2_ppm: np.ndarray               # (n2,) ascending
    f1_ppm: np.ndarray               # (n1,) ascending
    z: np.ndarray                    # (n1, n2) real (rr)
    nucleus: str = ""
    larmor_MHz: float = 0.0
    spin_rate_Hz: float = 0.0
    source: str = ""
    notes: list[str] = field(default_factory=list)
    ri: np.ndarray | None = None     # imag along F2 (direct)
    ir: np.ndarray | None = None     # imag along F1 (indirect)
    ii: np.ndarray | None = None     # imag-imag

    @property
    def has_hyper(self) -> bool:
        return self.ri is not None and self.ir is not None and self.ii is not None

    def _like(self, z, ri=None, ir=None, ii=None) -> "Data2D":
        return Data2D(self.f2_ppm, self.f1_ppm, z, self.nucleus,
                      self.larmor_MHz, self.spin_rate_Hz, self.source,
                      list(self.notes), ri, ir, ii)

    def region(self, f2_range=None, f1_range=None) -> "Data2D":
        s2 = np.ones(self.f2_ppm.shape, bool)
        s1 = np.ones(self.f1_ppm.shape, bool)
        if f2_range:
            s2 = ((self.f2_ppm >= min(f2_range)) & (self.f2_ppm <= max(f2_range)))
        if f1_range:
            s1 = ((self.f1_ppm >= min(f1_range)) & (self.f1_ppm <= max(f1_range)))
        ix = np.ix_(s1, s2)
        cut = (lambda a: a[ix] if a is not None else None)
        return Data2D(f2_ppm=self.f2_ppm[s2], f1_ppm=self.f1_ppm[s1],
                      z=self.z[ix], nucleus=self.nucleus,
                      larmor_MHz=self.larmor_MHz,
                      spin_rate_Hz=self.spin_rate_Hz, source=self.source,
                      notes=list(self.notes),
                      ri=cut(self.ri), ir=cut(self.ir), ii=cut(self.ii))

    def normalized(self) -> "Data2D":
        k = np.abs(self.z).max() or 1.0
        sc = (lambda a: a / k if a is not None else None)
        return self._like(self.z / k, sc(self.ri), sc(self.ir), sc(self.ii))

    def projection(self, axis: str = "f2", mode: str = "skyline") -> np.ndarray:
        red = np.max if mode == "skyline" else np.sum
        return red(self.z, axis=0 if axis == "f2" else 1)

    def transposed(self) -> "Data2D":
        """Swap F1 and F2 (dmfit Transpose RR)."""
        t = (lambda a: a.T if a is not None else None)
        return Data2D(self.f1_ppm, self.f2_ppm, self.z.T, self.nucleus,
                      self.larmor_MHz, self.spin_rate_Hz, self.source,
                      list(self.notes),
                      ri=t(self.ir), ir=t(self.ri), ii=t(self.ii))

    def reversed_axis(self, axis: str) -> "Data2D":
        """Mirror the intensity along one axis (dmfit Reverse F1/F2)."""
        ax = 1 if axis == "f2" else 0
        flip = (lambda a: np.flip(a, ax) if a is not None else None)
        return self._like(np.flip(self.z, ax), flip(self.ri), flip(self.ir),
                          flip(self.ii))

    def symmetrized(self) -> "Data2D":
        """Symmetrize about the F1=F2 diagonal (dmfit 2D 'Symmetric'): average z
        with its transpose on a common square grid. Useful to clean MQMAS
        auto-correlation ridges."""
        from scipy.interpolate import RegularGridInterpolator

        lo = max(self.f1_ppm.min(), self.f2_ppm.min())
        hi = min(self.f1_ppm.max(), self.f2_ppm.max())
        n = max(self.z.shape)
        g = np.linspace(lo, hi, n)
        interp = RegularGridInterpolator((self.f1_ppm, self.f2_ppm), self.z,
                                         bounds_error=False, fill_value=0.0)
        G1, G2 = np.meshgrid(g, g, indexing="ij")
        zg = interp(np.stack([G1.ravel(), G2.ravel()], -1)).reshape(n, n)
        zs = 0.5 * (zg + zg.T)
        return Data2D(g, g, zs, self.nucleus, self.larmor_MHz,
                      self.spin_rate_Hz, self.source, list(self.notes) +
                      ["symmetrized about the diagonal"])

    def diagonal(self) -> tuple[np.ndarray, np.ndarray]:
        """The trace along F1 = F2 (dmfit Extract Diag). Returns (ppm, amp)."""
        from scipy.interpolate import RegularGridInterpolator

        lo = max(self.f1_ppm.min(), self.f2_ppm.min())
        hi = min(self.f1_ppm.max(), self.f2_ppm.max())
        g = np.linspace(lo, hi, max(self.z.shape))
        interp = RegularGridInterpolator((self.f1_ppm, self.f2_ppm), self.z,
                                         bounds_error=False, fill_value=0.0)
        return g, interp(np.stack([g, g], axis=-1))

    def phased(self, axis: str, p0_deg: float, p1_deg: float,
               pivot_ppm: float | None = None) -> "Data2D":
        """Zero/first-order phase along one axis, TopSpin-style.

        axis 'f2' phases the direct dimension (rows), 'f1' the indirect one
        (columns). A processed 2rr keeps only the real part, so the dispersive
        companion is reconstructed by a Hilbert transform along the axis
        (Kramers-Kronig) -- the same trick ssNake uses to re-phase an
        already-transformed spectrum. p1 pivots about ``pivot_ppm`` (default:
        the axis centre).
        """
        if p0_deg == 0.0 and p1_deg == 0.0:
            return self
        ax = 1 if axis == "f2" else 0
        coords = self.f2_ppm if axis == "f2" else self.f1_ppm
        n = self.z.shape[ax]
        if n < 2:
            return self
        piv = _pivot_index(coords, pivot_ppm)
        ramp = (np.arange(n) - piv) / max(n - 1, 1)
        phase = np.deg2rad(p0_deg + p1_deg * ramp)
        shape = [1, 1]; shape[ax] = n
        if self.has_hyper:
            # exact hypercomplex rotation of all four quadrants
            c = np.cos(phase).reshape(tuple(shape))
            s = np.sin(phase).reshape(tuple(shape))
            rr, ri, ir, ii = self.z, self.ri, self.ir, self.ii
            if axis == "f2":
                return self._like(rr * c + ri * s, -rr * s + ri * c,
                                  ir * c + ii * s, -ir * s + ii * c)
            return self._like(rr * c + ir * s, ri * c + ii * s,
                              -rr * s + ir * c, -ri * s + ii * c)
        # fallback: reconstruct the dispersive part from the real spectrum
        from scipy.signal import hilbert

        analytic = hilbert(np.asarray(self.z, float), axis=ax)
        out = np.real(analytic * np.exp(-1j * phase).reshape(tuple(shape)))
        return self._like(out)

    def phase_line(self, axis: str, idx: int, p0_deg: float, p1_deg: float,
                   pivot_ppm: float | None = None) -> np.ndarray:
        """Phased 1D reference row/column for the live phasing preview, exact
        when the imaginary quadrant is available."""
        if axis == "f2":
            rr = self.z[idx]
            im = self.ri[idx] if self.ri is not None else None
            coords = self.f2_ppm
        else:
            rr = self.z[:, idx]
            im = self.ir[:, idx] if self.ir is not None else None
            coords = self.f1_ppm
        if im is None:
            return phase_1d(rr, coords, p0_deg, p1_deg, pivot_ppm)
        n = rr.size
        piv = _pivot_index(coords, pivot_ppm)
        ramp = (np.arange(n) - piv) / max(n - 1, 1)
        phase = np.deg2rad(p0_deg + p1_deg * ramp)
        return rr * np.cos(phase) + im * np.sin(phase)


def _pivot_index(coords: np.ndarray, pivot_ppm: float | None) -> int:
    n = len(coords)
    if pivot_ppm is None:
        return n // 2
    return int(np.argmin(np.abs(np.asarray(coords) - pivot_ppm)))


def phase_1d(y: np.ndarray, coords: np.ndarray, p0_deg: float, p1_deg: float,
             pivot_ppm: float | None = None) -> np.ndarray:
    """Phase one real trace exactly as :meth:`Data2D.phased` phases each line,
    so a live single-row/column preview matches the applied 2D result."""
    if p0_deg == 0.0 and p1_deg == 0.0:
        return np.asarray(y, float)
    from scipy.signal import hilbert

    y = np.asarray(y, float)
    n = y.size
    piv = _pivot_index(coords, pivot_ppm)
    ramp = (np.arange(n) - piv) / max(n - 1, 1)
    return np.real(hilbert(y) * np.exp(-1j * np.deg2rad(p0_deg + p1_deg * ramp)))


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
    h = data.hyper or {}
    d = Data2D(f2_ppm=f2_axis.values, f1_ppm=f1_axis.values, z=data.data,
               nucleus=data.nucleus, larmor_MHz=data.meta["larmor_MHz"],
               spin_rate_Hz=data.meta.get("masr_Hz") or 0.0,
               source=str(expno),
               ri=h.get("ri"), ir=h.get("ir"), ii=h.get("ii"))
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
_CS_SCALE_CACHE: dict = {}


def f1_cs_scale(nucleus: str, larmor_MHz: float, method: str = "3QMAS") -> float:
    """The factor c by which an isotropic chemical shift moves the F1 (indirect)
    axis relative to F2 in mrsimulator's named MQMAS method.

    A pure-CS site at δiso lands at F2 = δiso but F1 = c·δiso — for ²⁷Al 3QMAS
    c = -17/31 ≈ -0.548, NOT +1. This is the shear/scaling convention of the
    isotropic dimension and is spin- and coherence-order dependent; measuring it
    once from a reference simulation keeps it correct for every nucleus/method.
    """
    key = (nucleus, round(larmor_MHz, 3), method)
    if key in _CS_SCALE_CACHE:
        return _CS_SCALE_CACHE[key]
    from mrsimulator import Simulator
    from mrsimulator.method import SpectralDimension
    from mrsimulator.method import lib as method_lib
    from mrsimulator.spin_system.isotope import Isotope
    from mrsimulator.utils.collection import single_site_system_generator

    B0 = larmor_MHz / abs(Isotope(symbol=nucleus).gyromagnetic_ratio)
    ref = 50.0                                   # a pure-CS probe shift (ppm)
    sw = 400.0 * larmor_MHz                       # wide window so it never clips
    systems = single_site_system_generator(
        isotope=nucleus, isotropic_chemical_shift=[ref],
        quadrupolar={"Cq": [2e4], "eta": [0.0]}, abundance=[100.0])
    m = getattr(method_lib, METHODS[method])(
        channels=[nucleus], magnetic_flux_density=B0,
        spectral_dimensions=[SpectralDimension(count=256, spectral_width=sw),
                             SpectralDimension(count=256, spectral_width=sw)])
    sim = Simulator(spin_systems=systems, methods=[m]); sim.run()
    ds = sim.methods[0].simulation
    c1 = ds.x[1].coordinates
    f1 = c1.value if str(c1.unit) == "ppm" else c1.to("Hz").value / larmor_MHz
    z = np.asarray(ds.y[0].components[0].real, float)   # (x1, x0)
    proj = np.clip(z.sum(axis=1), 0, None)              # onto F1
    c = float((proj * f1).sum() / proj.sum() / ref) if proj.sum() else 1.0
    _CS_SCALE_CACHE[key] = c
    return c


@dataclass
class Kernel2D:
    f2_ppm: np.ndarray               # (n2,)
    f1_ppm: np.ndarray               # (n1,)
    K: np.ndarray                    # (ngrid, n1, n2) basis subspectra
    cq_grid_MHz: np.ndarray
    eta_grid: np.ndarray
    f1_cs_scale: float = 1.0         # δiso moves F1 by c·δiso (see f1_cs_scale())

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

    # f1_window is the desired δ1-ISOTROPIC range (CS on the diagonal, ≈ δiso).
    # mrsimulator returns its sheared F1 scaled by c (= f1_cs_scale), so simulate
    # over the native window c·(iso window) and rescale the result back to iso.
    cs = f1_cs_scale(nucleus, larmor_MHz, method)
    nat = sorted((cs * f1_window[0], cs * f1_window[1]))
    f1n_lo, f1n_hi = nat[0], nat[1]

    sw2 = abs(f2_window[0] - f2_window[1]) * larmor_MHz
    sw1 = abs(f1n_hi - f1n_lo) * larmor_MHz
    ro2 = 0.5 * (f2_window[0] + f2_window[1]) * larmor_MHz
    ro1 = 0.5 * (f1n_hi + f1n_lo) * larmor_MHz
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
    f1 = np.asarray(f1) / cs                       # sheared native -> δ1-isotropic
    o2, o1 = np.argsort(f2), np.argsort(f1)        # (cs < 0 flips F1 -> re-sort)
    K = np.array([np.asarray(dv.components[0].real, float)[np.ix_(o1, o2)]
                  for dv in ds.y])
    kernel = Kernel2D(f2_ppm=np.asarray(f2)[o2], f1_ppm=np.asarray(f1)[o1],
                      K=K, cq_grid_MHz=cq_grid, eta_grid=eta_grid, f1_cs_scale=cs)
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
    pos = v["isotropic_chemical_shift_ppm"]
    # The kernel's F1 axis is the δ1-isotropic convention (rescaled at build
    # time, see build_mqmas_kernel), so a pure-CS site lies on the diagonal:
    # δiso moves F2 and F1 together (slope 1).
    if model in ("czjzek", "ext_czjzek"):
        # Glass: the isotropic-shift *distribution* (dmfit's dCS) elongates the
        # peak ALONG THE DIAGONAL (both dims move with δiso), while the round
        # point broadening (dmfit's wid) is isotropic.
        z = _shift_smear_2d(z, kernel.f2_ppm, kernel.f1_ppm, pos,
                            v.get("shift_fwhm_ppm", 0.0))
        line = v.get("line_fwhm_ppm", 0.0)
        if line > 0.05:
            z = _broaden_2d(z, kernel.f1_ppm, kernel.f2_ppm, line, line)
    else:
        z = _shift_2d(z, kernel.f2_ppm, kernel.f1_ppm, pos, pos)
        fw = v.get("shift_fwhm_ppm", 1.0)
        z = _broaden_2d(z, kernel.f1_ppm, kernel.f2_ppm, fw, fw)
    peak = z.max()
    return v["amplitude"] * (z / peak) if peak > 0 else z


def _shift_smear_2d(z: np.ndarray, f2: np.ndarray, f1: np.ndarray,
                    pos_ppm: float, cs_fwhm_ppm: float) -> np.ndarray:
    """Place at δiso = pos_ppm on the diagonal (F2 and F1 both +pos) and, if
    cs_fwhm_ppm > 0, smear the chemical-shift *distribution* along the diagonal
    — the along-CS elongation of a disordered/glassy site, distinct from the
    round point broadening. (The kernel F1 axis is the δ1-isotropic convention,
    so the CS axis is the diagonal.)"""
    if cs_fwhm_ppm <= 0.1:
        return _shift_2d(z, f2, f1, pos_ppm, pos_ppm)
    sigma = cs_fwhm_ppm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    offs = np.linspace(-2.0 * sigma, 2.0 * sigma, 5)
    w = np.exp(-0.5 * (offs / sigma) ** 2)
    w /= w.sum()
    out = np.zeros_like(z)
    for o, wi in zip(offs, w):
        d = pos_ppm + float(o)
        out += wi * _shift_2d(z, f2, f1, d, d)
    return out


def _shift_2d(z: np.ndarray, f2: np.ndarray, f1: np.ndarray,
              pos_f2: float, pos_f1: float) -> np.ndarray:
    """Translate a δiso = 0 subspectrum by pos_f2 along F2 and pos_f1 along F1."""
    if abs(pos_f2) < 1e-12 and abs(pos_f1) < 1e-12:
        return z
    out = np.empty_like(z)
    for i in range(z.shape[0]):
        out[i] = np.interp(f2 - pos_f2, f2, z[i], left=0.0, right=0.0)
    if abs(pos_f1) >= 1e-12:
        out2 = np.empty_like(out)
        for j in range(z.shape[1]):
            out2[:, j] = np.interp(f1 - pos_f1, f1, out[:, j], left=0.0, right=0.0)
        return out2
    return out


def mqmas_f1_axis(kernel: Kernel2D, recipe) -> np.ndarray:
    """The kernel's δ1-isotropic F1 axis mapped to the experiment's F1 axis by
    the fitted reference offset β (Recipe.mqmas_f1_ref_ppm)."""
    b = float(getattr(recipe, "mqmas_f1_ref_ppm", 0.0))
    return np.asarray(kernel.f1_ppm) + b


def simulate_2d(recipe, kernel: Kernel2D):
    """(total, per_site) 2D model arrays on the kernel's NATIVE grid.

    The isotropic-axis affine (mqmas_f1_scale/ref) is a property of how the
    experiment's F1 was referenced, not of the model — apply it to the F1 axis
    (mqmas_f1_axis) when comparing to or displaying over the experiment.
    """
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
        kernel = _kernel_for(recipe, data, method)

    # The kernel F1 axis is the δ1-isotropic convention (CS on the diagonal); an
    # experimental F1 axis differs only by a small referencing offset β, so
    # sample the experiment at (f1_iso + β) and fit β (the isotropic-axis
    # reference). δiso sets the diagonal position, so β is independent.
    from scipy.interpolate import RegularGridInterpolator

    interp = RegularGridInterpolator(
        (data.f1_ppm, data.f2_ppm), data.z, bounds_error=False, fill_value=0.0)
    G1n, G2 = np.meshgrid(kernel.f1_ppm, kernel.f2_ppm, indexing="ij")

    def sample_exp(b):
        pts = np.stack([(G1n + b).ravel(), G2.ravel()], axis=-1)
        return interp(pts).reshape(kernel.shape)

    scale = float(np.abs(data.z).max()) or 1.0
    params = _make_params(recipe)
    _apply_params(recipe, params)
    vary_ref = getattr(recipe, "mqmas_f1_ref_vary", True)
    # β is the isotropic-axis (F1) referencing offset. Peaks that start far apart
    # in F1 give the local optimiser no gradient, so coarse-search β once for the
    # best model/experiment overlap at the initial params (place δiso near the
    # peaks first — δiso sets the diagonal and is anchored by F2).
    if vary_ref:
        m0, _ = simulate_2d(recipe, kernel)
        m0f = m0.ravel(); m0n = float(np.sqrt((m0f * m0f).sum())) or 1.0
        best = (-1.0, 0.0)
        for b in np.linspace(-40.0, 40.0, 81):
            e = sample_exp(b).ravel()
            den = np.sqrt((e * e).sum()) * m0n
            corr = float((e * m0f).sum()) / den if den > 0 else 0.0
            if corr > best[0]:
                best = (corr, float(b))
        b0 = best[1]
    else:
        b0 = float(getattr(recipe, "mqmas_f1_ref_ppm", 0.0))
    params.add("mqmas_f1_ref_ppm", value=b0, min=-40.0, max=40.0, vary=vary_ref)

    def residual(p):
        _apply_params(recipe, p)
        recipe.mqmas_f1_ref_ppm = float(p["mqmas_f1_ref_ppm"].value)
        total, _ = simulate_2d(recipe, kernel)
        return ((total - sample_exp(recipe.mqmas_f1_ref_ppm)) / scale).ravel()

    # amplitude pre-scale so the optimizer starts on-scale
    recipe.mqmas_f1_ref_ppm = b0
    total0, _ = simulate_2d(recipe, kernel)
    e0 = sample_exp(b0)
    denom = float((total0 * total0).sum())
    if denom > 0:
        kk = float((e0 * total0).sum()) / denom
        for i, site in enumerate(recipe.sites):
            from larmor.fit import _lmfit_name

            name = _lmfit_name(i, site, "amplitude")
            if params[name].vary:
                params[name].value *= kk

    result = lmfit.minimize(residual, params, method="least_squares")
    _apply_params(recipe, result.params)
    recipe.mqmas_f1_ref_ppm = float(result.params["mqmas_f1_ref_ppm"].value)
    z_fit, per_site = simulate_2d(recipe, kernel)
    z_exp = sample_exp(recipe.mqmas_f1_ref_ppm)
    rmsd = float(np.sqrt(np.mean((z_fit - z_exp) ** 2)) / scale)
    recipe.fit_rmsd = rmsd
    return Fit2DResult(recipe=recipe, lmfit_result=result, kernel=kernel,
                       z_fit=z_fit, per_site=per_site, rmsd=rmsd)


def _kernel_for(recipe, data: "Data2D", method: str = "3QMAS") -> Kernel2D:
    """Build (or fetch) the MQMAS kernel for a dataset. The kernel F1 grid is the
    δ1-isotropic convention; a Bruker MQMAS F1 axis is already ≈ that convention,
    so use the experiment's F1 range (padded) as the iso window — the sites'
    diagonal placements (≈ δiso) then land on-grid."""
    f2w = (float(data.f2_ppm.max()), float(data.f2_ppm.min()))
    f1lo, f1hi = float(data.f1_ppm.min()), float(data.f1_ppm.max())
    pad = 0.15 * (f1hi - f1lo)
    return build_mqmas_kernel(
        recipe.nucleus or data.nucleus,
        recipe.larmor_frequency_MHz or data.larmor_MHz,
        f2_window=f2w, f1_window=(f1hi + pad, f1lo - pad),
        method=method, n2=MQMAS_SETTINGS["n2"], n1=MQMAS_SETTINGS["n1"],
        n_cq=MQMAS_SETTINGS["n_cq"], n_eta=MQMAS_SETTINGS["n_eta"],
        cq_max_MHz=MQMAS_SETTINGS["cq_max_MHz"])
