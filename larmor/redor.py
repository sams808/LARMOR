"""REDOR analysis: dipolar couplings and distances from S0/S dephasing.

Physics
-------
REDOR reintroduces the heteronuclear dipolar coupling that MAS averages out.
The dephasing curve ΔS/S0 = (S0 - S)/S0 vs recoupling time Ntr depends only
on the dipolar coupling D (Gullion & Schaefer, J. Magn. Reson. 81, 196, 1989).

Two analysis levels, both provided:

1. **Universal short-time parabola** (model-free, robust):
       ΔS/S0 ≈ (4/15) (D Ntr)^2 ... valid while ΔS/S0 < ~0.2
   This is the standard way to extract D without assuming a spin geometry,
   and the only defensible one when the number of neighbours is unknown
   (Bertmer & Eckert, Solid State Nucl. Magn. Reson. 15, 139, 1999).

2. **Full isolated-pair curve** (Bessel expansion): the exact powder-averaged
   S/S0 for one I-S pair, fitted over the whole curve.

For a heterogeneous sample the meaningful quantity is the second moment
M2 = (4/15)·D²·... -- LARMOR reports D, M2 and the equivalent pair distance,
and states which regime was used.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: mu0/(4 pi) * hbar / (2 pi)  in  Hz * m^3 * T^-2 ... assembled below
_MU0_OVER_4PI = 1e-7          # T m / A
_HBAR = 1.054571817e-34       # J s

#: gyromagnetic ratios (MHz/T) for the common REDOR partners
GAMMA_MHZ_T = {
    "1H": 42.577478, "19F": 40.077, "13C": 10.7084, "15N": -4.3173,
    "31P": 17.235, "27Al": 11.1031, "23Na": 11.2688, "29Si": -8.4654,
    "11B": 13.663, "7Li": 16.546, "25Mg": -2.6083, "17O": -5.7743,
    "35Cl": 4.1765, "45Sc": 10.3591, "51V": 11.2133, "71Ga": 13.0208,
    "93Nb": 10.4523, "133Cs": 5.6234, "207Pb": 8.8073,
}


def dipolar_constant_hz(iso1: str, iso2: str, r_angstrom: float) -> float:
    """D = (mu0/4pi) * gamma1 * gamma2 * hbar / r^3   [Hz]."""
    g1 = GAMMA_MHZ_T[iso1] * 1e6 * 2 * np.pi      # rad/(s T)
    g2 = GAMMA_MHZ_T[iso2] * 1e6 * 2 * np.pi
    r = r_angstrom * 1e-10
    return abs(_MU0_OVER_4PI * g1 * g2 * _HBAR / (2 * np.pi) / r ** 3)


def distance_angstrom(iso1: str, iso2: str, d_hz: float) -> float:
    """Inverse of dipolar_constant_hz."""
    g1 = GAMMA_MHZ_T[iso1] * 1e6 * 2 * np.pi
    g2 = GAMMA_MHZ_T[iso2] * 1e6 * 2 * np.pi
    r3 = abs(_MU0_OVER_4PI * g1 * g2 * _HBAR / (2 * np.pi) / d_hz)
    return float(r3 ** (1 / 3) * 1e10)


#: Dephasing phase after N rotor periods for an isolated pair, derived by
#: integrating the MAS dipolar frequency over a rotor cycle with a pi pulse at
#: Tr/2 (Gullion & Schaefer):
#:
#:   omega_D(t)/2pi = -D [ sqrt2 sin2b cos(a + w_r t) + sin^2 b cos(2a + 2w_r t) ]
#:   dphi_1 = int_0^{Tr/2} - int_{Tr/2}^{Tr}  ->  the cos(2a+2w_r t) term
#:            integrates to zero over each half rotor period
#:   dphi   = N * 4 sqrt2 * D * Tr * sin(2b) * sin(a)
#:
#: so with lambda = D * N * Tr  (D in Hz, N*Tr in s):
#:   dphi = GEOM_PREFACTOR * lambda * sin(2b) * sin(a)
GEOM_PREFACTOR = 4.0 * np.sqrt(2.0)

#: small-lambda limit of the powder average (see test_redor_parabola_is_the_
#: small_lambda_limit): 1 - <cos(dphi)> ~ <dphi^2>/2 = (64/15) lambda^2
#: because <sin^2 a> = 1/2 and <sin^2 2b>_powder = 8/15.
SHORT_TIME_COEFF = 64.0 / 15.0


def _powder_grid(na: int = 200, nb: int = 100):
    alpha = (np.arange(na) + 0.5) * 2 * np.pi / na
    cosb = (np.arange(nb) + 0.5) / nb              # uniform in cos(beta)
    A, C = np.meshgrid(alpha, cosb, indexing="ij")
    beta = np.arccos(C)
    return GEOM_PREFACTOR * np.sin(2 * beta) * np.sin(A)


def redor_pair_curve(d_hz: float, ntr_s: np.ndarray) -> np.ndarray:
    """Powder-averaged ΔS/S0 for an isolated spin pair.

    ΔS/S0 = 1 - <cos(Δφ)>_powder, by direct numerical powder averaging over
    (α, cos β) -- exact within the quadrature and independently checkable,
    unlike a memorized Bessel prefactor.
    """
    lam = np.atleast_1d(np.asarray(ntr_s, float)) * d_hz
    geom = _powder_grid()
    out = np.empty(lam.shape)
    for i, l in enumerate(lam):
        out[i] = 1.0 - float(np.mean(np.cos(l * geom)))
    return out


def short_time_curve(d_hz: float, ntr_s: np.ndarray) -> np.ndarray:
    """Universal short-time parabola: ΔS/S0 = (64/15)(D·N·Tr)^2.

    Valid while ΔS/S0 < ~0.2. Model-free: it assumes only the second-order
    expansion of the powder average, no spin geometry and no neighbour count.
    """
    lam = np.asarray(ntr_s, float) * d_hz
    return SHORT_TIME_COEFF * lam ** 2


@dataclass
class RedorResult:
    ntr_s: np.ndarray
    ds_s0: np.ndarray
    d_hz: float
    d_err: float | None
    m2: float                       # second moment, rad^2 s^-2
    regime: str                     # "short-time parabola" | "isolated pair"
    n_used: int
    distance_A: float | None = None
    pair: tuple[str, str] | None = None
    curve: object = None
    notes: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        err = f" ± {self.d_err:.1f}" if self.d_err else ""
        s = f"D = {self.d_hz:.1f}{err} Hz  ({self.regime}, {self.n_used} pts)"
        if self.distance_A:
            s += f"   r({self.pair[0]}–{self.pair[1]}) = {self.distance_A:.2f} Å"
        return s


def read_redor_txt(path: str | Path) -> tuple[np.ndarray, np.ndarray, float]:
    """Parse TopSpin's redor.txt: returns (n_rotor_cycles, ΔS/S0, masr_Hz)."""
    text = Path(path).read_text().splitlines()
    masr = None
    rows = []
    for ln in text:
        parts = ln.split()
        if ln.strip().startswith("Spinning speed"):
            masr = float(parts[-1])
        if len(parts) == 5 and parts[0].isdigit():
            rows.append([float(v) for v in parts])
    if not rows:
        raise ValueError(f"no S0/S rows found in {path}")
    arr = np.array(rows)
    n, s0, s = arr[:, 0], arr[:, 1], arr[:, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        ds = 1.0 - np.where(s0 != 0, s / s0, np.nan)
    return n, ds, (masr or 0.0)


def analyze(ntr_s: np.ndarray, ds_s0: np.ndarray,
            pair: tuple[str, str] | None = None,
            regime: str = "auto", max_ds: float = 0.2) -> RedorResult:
    """Fit a REDOR dephasing curve for the dipolar coupling D.

    regime "short" uses only points below `max_ds` with the universal
    parabola; "pair" fits the full isolated-pair curve; "auto" picks "short"
    when at least 3 points lie in the parabolic regime (the safe default for
    unknown coordination), else "pair".
    """
    import lmfit

    ntr = np.asarray(ntr_s, float)
    ds = np.asarray(ds_s0, float)
    ok = np.isfinite(ntr) & np.isfinite(ds)
    ntr, ds = ntr[ok], ds[ok]
    notes = []

    low = ds < max_ds
    if regime == "auto":
        regime = "short" if low.sum() >= 3 else "pair"
        notes.append(f"regime auto-selected: {regime} "
                     f"({int(low.sum())} points below ΔS/S0 = {max_ds})")

    if regime == "short":
        x, y = ntr[low], ds[low]
        if x.size < 2:
            raise ValueError("not enough points in the short-time regime; "
                             "use regime='pair'")
        params = lmfit.Parameters()
        params.add("d", value=max(1.0, 1.0 / (x.max() or 1.0)), min=0.0)
        out = lmfit.minimize(
            lambda p: short_time_curve(p["d"].value, x) - y, params,
            method="leastsq")
        d = float(out.params["d"].value)
        derr = (float(out.params["d"].stderr)
                if out.params["d"].stderr else None)
        n_used = int(x.size)
        label = "short-time parabola"
        curve = (lambda tt, dd=d: short_time_curve(dd, tt))
        notes.append("model-free: assumes only the short-time expansion, no "
                     "spin geometry")
    else:
        params = lmfit.Parameters()
        params.add("d", value=max(1.0, 1.0 / (ntr.max() or 1.0)), min=0.0)
        out = lmfit.minimize(
            lambda p: redor_pair_curve(p["d"].value, ntr) - ds, params,
            method="leastsq")
        d = float(out.params["d"].value)
        derr = (float(out.params["d"].stderr)
                if out.params["d"].stderr else None)
        n_used = int(ntr.size)
        label = "isolated pair"
        curve = (lambda tt, dd=d: redor_pair_curve(dd, np.asarray(tt)))
        notes.append("assumes ONE isolated I–S pair; invalid for multi-spin "
                     "environments (D would be overestimated)")

    # van Vleck second moment of the heteronuclear pair, rad^2 s^-2
    m2 = (4.0 / 15.0) * (2 * np.pi * d) ** 2
    dist = None
    if pair:
        try:
            dist = distance_angstrom(pair[0], pair[1], d)
        except KeyError as exc:
            notes.append(f"no gyromagnetic ratio for {exc}; distance skipped")
    return RedorResult(ntr_s=ntr, ds_s0=ds, d_hz=d, d_err=derr, m2=m2,
                       regime=label, n_used=n_used, distance_A=dist,
                       pair=pair, curve=curve, notes=notes)


def analyze_expno(expno: str | Path, pair: tuple[str, str] | None = None,
                  procno: int = 1, **kwargs) -> RedorResult:
    """Analyze TopSpin's redor.txt inside an EXPNO (read-only)."""
    p = Path(expno) / "pdata" / str(procno) / "redor.txt"
    n, ds, masr = read_redor_txt(p)
    if not masr:
        raise ValueError("no spinning speed in redor.txt; cannot convert "
                         "rotor cycles to seconds")
    return analyze(n / masr, ds, pair=pair, **kwargs)
