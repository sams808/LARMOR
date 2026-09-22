"""Herzfeld-Berger analysis of MAS spinning-sideband intensities.

Under magic-angle spinning a chemical-shift-anisotropy powder pattern breaks
into a centreband and sidebands at +-N nu_rot whose INTENSITIES map the
shielding tensor (Herzfeld & Berger, J. Chem. Phys. 73, 6021 (1980)). This
module computes those intensities exactly for a given (zeta, eta), measures
them off a spectrum, and inverts the measured set to (zeta, eta) -- the
classical sideband analysis, which is a measurement rather than a lineshape
fit and therefore works where the individual sidebands are too broad or too
overlapped for a full csa_mas fit to be trusted.

Conventions match the rest of LARMOR: ``zeta_ppm`` is the Haeberlen SHIELDING
anisotropy exactly as the csa_mas model passes it to mrsimulator (the shift
anisotropy is -zeta), ``eta`` in [0, 1], sideband N sits at
delta_iso + N nu_rot / nu_0 (ppm, N > 0 towards higher shift). The
intensities are computed from the exact time-dependent frequency of each
crystallite: the FID of one orientation is exp(i Phi(t)) with Phi periodic in
the rotor period, so its Fourier coefficients are the sideband amplitudes and
the powder average of their squares is the intensity (gamma-averaging over
the rotor phase is exact, so only (beta, gamma) are sampled). No sideband
count is truncated: the time grid is sized from |zeta| nu_0 / nu_rot.

Qt-free; the desktop dialog (larmor/desktop/herzfeld_berger_dialog.py) is a
consumer.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: the magic angle
MAGIC_RAD = float(np.arccos(1.0 / np.sqrt(3.0)))

__all__ = [
    "MAGIC_RAD", "HBFit", "sideband_intensities", "measure_manifold",
    "fit_sidebands", "sideband_positions_ppm", "n_sidebands_needed",
]


def n_sidebands_needed(zeta_ppm: float, larmor_MHz: float,
                       spin_rate_Hz: float) -> int:
    """How many sidebands each side carry intensity: the static pattern
    spans ~1.5|zeta| (eta = 1 edge to edge), so orders up to that span in
    units of nu_rot, plus a margin for the tails."""
    if spin_rate_Hz <= 0 or larmor_MHz <= 0:
        return 0
    span_Hz = 1.5 * abs(float(zeta_ppm)) * float(larmor_MHz)
    return int(np.ceil(span_Hz / float(spin_rate_Hz))) + 3


def _shift_tensor_aniso(zeta_ppm: float, eta: float) -> np.ndarray:
    """Anisotropic part of the SHIFT tensor in its PAS (ppm), from the
    shielding-convention (zeta, eta) LARMOR stores: delta_zz - delta_iso =
    -zeta, delta_xx - delta_iso = zeta (1 + eta) / 2, delta_yy - delta_iso =
    zeta (1 - eta) / 2 (traceless)."""
    z, e = float(zeta_ppm), float(eta)
    return np.diag([z * (1.0 + e) / 2.0, z * (1.0 - e) / 2.0, -z])


def _rotations(gamma: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """PAS -> rotor frame, R = Ry(beta) Rz(gamma), for arrays of angles ->
    (K, 3, 3). The third Euler angle (a rotation about the rotor axis) only
    shifts the time origin of the periodic trajectory and leaves every |c_N|
    unchanged, so it is not sampled; gamma, the azimuth of the tensor about
    its own z axis, DOES change the trajectory and must be averaged."""
    cg, sg = np.cos(gamma), np.sin(gamma)
    cb, sb = np.cos(beta), np.sin(beta)
    zero, one = np.zeros_like(cg), np.ones_like(cg)
    rz = np.stack([np.stack([cg, -sg, zero], -1),
                   np.stack([sg, cg, zero], -1),
                   np.stack([zero, zero, one], -1)], -2)
    ry = np.stack([np.stack([cb, zero, sb], -1),
                   np.stack([zero, one, zero], -1),
                   np.stack([-sb, zero, cb], -1)], -2)
    return ry @ rz


def _powder(n_beta: int, n_gamma: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gauss-Legendre in cos(beta) on [-1, 1] x uniform gamma: (gamma, beta,
    weights) flattened, weights summing to 1."""
    xb, wb = np.polynomial.legendre.leggauss(int(n_beta))
    beta = np.arccos(xb)
    gamma = np.linspace(0.0, 2.0 * np.pi, int(n_gamma), endpoint=False)
    G, B = np.meshgrid(gamma, beta, indexing="ij")
    W = np.broadcast_to(wb[None, :], G.shape)
    w = W.ravel() / W.sum()
    return G.ravel(), B.ravel(), w


def sideband_intensities(zeta_ppm: float, eta: float, larmor_MHz: float,
                         spin_rate_Hz: float, n_max: int | None = None,
                         n_beta: int = 24, n_gamma: int = 48) -> dict[int, float]:
    """Relative intensities {N: I_N} of the MAS sidebands of a CSA site,
    normalised to sum to 1 over ALL orders (Parseval: |exp(i Phi)| = 1), then
    returned for |N| <= n_max (default: every order that carries intensity).

    A static experiment (spin_rate_Hz <= 0) or a vanishing anisotropy has no
    sidebands: {0: 1.0}.
    """
    nu0 = float(larmor_MHz)
    nur = float(spin_rate_Hz)
    if nur <= 0 or nu0 <= 0 or abs(float(zeta_ppm)) < 1e-9:
        return {0: 1.0}
    need = n_sidebands_needed(zeta_ppm, nu0, nur)
    if n_max is None:
        n_max = need
    # time grid over one rotor period: enough points that no order aliases
    m = 1 << int(np.ceil(np.log2(max(4 * (need + 2), 32))))
    m = min(m, 1 << 14)
    gamma, beta, w = _powder(n_beta, n_gamma)
    rot = _rotations(gamma, beta)                                  # (K,3,3)
    sig = rot @ _shift_tensor_aniso(zeta_ppm, eta) @ np.swapaxes(rot, 1, 2)
    phase = 2.0 * np.pi * np.arange(m) / m                          # omega_r t
    st, ct = np.sin(MAGIC_RAD), np.cos(MAGIC_RAD)
    b = np.stack([st * np.cos(phase), st * np.sin(phase),
                  np.full(m, ct)], axis=-1)                         # (M,3)
    # frequency (Hz) of every crystallite along the rotor period
    freq = np.einsum("kij,ti,tj->kt", sig, b, b) * nu0
    # integrate the harmonics analytically: Phi(t) = 2 pi Int freq dt
    W = np.fft.fft(freq, axis=1) / m
    k = np.fft.fftfreq(m, d=1.0 / m)                                # harmonic index
    with np.errstate(divide="ignore", invalid="ignore"):
        coef = np.where(k != 0, W / (1j * k * nur), 0.0)
    Phi = np.fft.ifft(coef, axis=1).real * m                        # radians
    s = np.exp(1j * Phi)
    c = np.fft.fft(s, axis=1) / m                                   # sideband amplitudes
    I = np.einsum("k,kn->n", w, np.abs(c) ** 2)                     # powder average
    orders = np.fft.fftfreq(m, d=1.0 / m).astype(int)
    total = float(I.sum())
    out = {}
    for n, val in zip(orders, I):
        if abs(int(n)) <= int(n_max):
            out[int(n)] = float(val / total)
    return dict(sorted(out.items()))


def sideband_positions_ppm(delta_iso_ppm: float, larmor_MHz: float,
                           spin_rate_Hz: float, orders) -> dict[int, float]:
    """ppm position of each sideband order: delta_iso + N nu_rot / nu_0."""
    step = float(spin_rate_Hz) / float(larmor_MHz)
    return {int(n): float(delta_iso_ppm) + int(n) * step for n in orders}


def measure_manifold(x_ppm, y, centre_ppm: float, larmor_MHz: float,
                     spin_rate_Hz: float, n_each_side: int,
                     half_width_ppm: float | None = None,
                     floor: float | None = None) -> dict[int, float]:
    """Integrate the spectrum in a window of +-half_width_ppm around each
    sideband position (centre + N nu_rot/nu_0) for N in -n..n, after
    subtracting a flat floor (default: the median of the outer 10 % of the
    trace). Negative integrals are clipped to 0. Returns {N: integral}; the
    caller normalises."""
    x = np.asarray(x_ppm, float)
    yy = np.asarray(y, float)
    order = np.argsort(x)
    x, yy = x[order], yy[order]
    step = float(spin_rate_Hz) / float(larmor_MHz)
    if half_width_ppm is None:
        half_width_ppm = 0.35 * step
    if floor is None:
        k = max(1, int(0.05 * x.size))
        floor = float(np.median(np.concatenate([yy[:k], yy[-k:]])))
    yy = yy - floor
    out = {}
    for n in range(-int(n_each_side), int(n_each_side) + 1):
        c = float(centre_ppm) + n * step
        m = (x >= c - half_width_ppm) & (x <= c + half_width_ppm)
        out[n] = float(max(_trapz(yy[m], x[m]), 0.0)) if m.sum() > 1 else 0.0
    return out


_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")


@dataclass
class HBFit:
    zeta_ppm: float
    eta: float
    rms: float                       # RMS misfit of normalised intensities
    orders: list = field(default_factory=list)
    measured: dict = field(default_factory=dict)     # normalised over `orders`
    fitted: dict = field(default_factory=dict)       # normalised over `orders`
    ok: bool = True
    message: str = ""


def fit_sidebands(measured: dict[int, float], larmor_MHz: float,
                  spin_rate_Hz: float, zeta_max_ppm: float | None = None) -> HBFit:
    """Invert measured sideband intensities to (zeta, eta).

    Both the measurement and the model are normalised over the SAME set of
    orders, so a truncated manifold is compared fairly. A coarse (zeta, eta)
    grid (both signs of zeta -- the N/-N asymmetry fixes the sign) seeds a
    bounded least-squares refinement. Needs at least three orders with
    intensity, two of them sidebands.
    """
    from scipy.optimize import least_squares

    orders = sorted(int(n) for n in measured)
    vals = np.array([max(float(measured[n]), 0.0) for n in orders])
    if len(orders) < 3 or vals.sum() <= 0 or np.count_nonzero(vals) < 3:
        return HBFit(0.0, 0.0, np.inf, orders, {}, {}, ok=False,
                     message="need a centreband and at least two sidebands "
                             "with intensity")
    meas = vals / vals.sum()
    step_ppm = float(spin_rate_Hz) / float(larmor_MHz)
    if zeta_max_ppm is None:
        # the static span is ~1.5 |zeta|; a manifold measured to order n_max
        # cannot come from a tensor much wider than n_max spacings
        zeta_max_ppm = max(2.0 * max(abs(n) for n in orders) * step_ppm,
                           2.0 * step_ppm)

    def model(z, e):
        I = sideband_intensities(z, e, larmor_MHz, spin_rate_Hz,
                                 n_max=max(abs(n) for n in orders),
                                 n_beta=16, n_gamma=32)
        v = np.array([I.get(n, 0.0) for n in orders])
        s = v.sum()
        return v / s if s > 0 else v

    def resid(p):
        return model(p[0], p[1]) - meas

    # coarse grid: the cost landscape is narrow in zeta (its scale is the
    # sideband spacing), so step at a fraction of nu_rot/nu_0 and refine from
    # the best point of EACH sign -- a coarse grid once handed the refinement
    # a wrong-sign start it could not leave
    n_z = max(60, int(np.ceil(2 * zeta_max_ppm / (0.2 * step_ppm))))
    n_z = min(n_z, 400)
    grid = []
    for z in np.linspace(-zeta_max_ppm, zeta_max_ppm, n_z):
        if abs(z) < 0.05 * step_ppm:
            continue
        for e in np.linspace(0.0, 1.0, 11):
            grid.append((float(np.sum(resid((z, e)) ** 2)), z, e))
    starts = []
    for sign in (1, -1):
        cands = [g for g in grid if np.sign(g[1]) == sign]
        if cands:
            starts.append(min(cands)[1:])
    sols = []
    for start in starts:
        sol = least_squares(resid, start, bounds=([-zeta_max_ppm, 0.0],
                                                  [zeta_max_ppm, 1.0]),
                            diff_step=(1e-3, 1e-3))
        sols.append((float(np.sum(sol.fun ** 2)), sol))
    sol = min(sols, key=lambda t: t[0])[1]
    z, e = float(sol.x[0]), float(np.clip(sol.x[1], 0.0, 1.0))
    fit = model(z, e)
    rms = float(np.sqrt(np.mean((fit - meas) ** 2)))
    return HBFit(z, e, rms, orders,
                 {n: float(v) for n, v in zip(orders, meas)},
                 {n: float(v) for n, v in zip(orders, fit)},
                 ok=True,
                 message="" if abs(z) < 0.95 * zeta_max_ppm else
                 "zeta at the search limit: include more sidebands")
