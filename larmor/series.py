"""Relaxation / series suite: one engine for every arrayed experiment.

Kinds (auto-detected from the pulse program when possible):
  satrec  I(t) = I0 (1 - exp(-(t/T1)^b))          saturation recovery
  invrec  I(t) = I0 (1 - 2f exp(-(t/T1)^b))       inversion recovery
  cpmg    I(t) = I0 exp(-(t/T2)^b)                T2 / echo-train decay
  t1rho   I(t) = I0 exp(-(t/T1rho)^b)             spin-lock decay
  decay   generic mono-exponential decay

Two analysis modes:
  * window integration (fast, dmfit/TopSpin style) -- `analyze`
  * **per-site fitting** -- `analyze_per_site`: the 1D model is fitted to
    every slice with shared lineshapes and per-slice amplitudes, so
    OVERLAPPING sites get separate relaxation times. Window integration
    cannot do this: it lumps everything under the window into one number.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from larmor.satrec import process_slices, read_vdlist

KINDS = ("satrec", "invrec", "cpmg", "t1rho", "decay")

#: pulse-program name fragments -> series kind
_PULPROG_HINTS = (
    ("satrec", "satrec"), ("sr_", "satrec"), ("satur", "satrec"),
    ("t1ir", "invrec"), ("invrec", "invrec"), ("ir_", "invrec"),
    ("cpmg", "cpmg"), ("t2", "cpmg"), ("echotrain", "cpmg"),
    ("t1rho", "t1rho"), ("splock", "t1rho"), ("spinlock", "t1rho"),
)


def detect_kind(expno: str | Path) -> str | None:
    """Guess the series kind from acqus PULPROG. Returns None if unsure --
    LARMOR asks rather than guessing wrong."""
    try:
        import nmrglue as ng

        dic, _ = ng.bruker.read_acqus_file(str(expno))
        pulprog = str(dic.get("PULPROG", "")).lower()
    except Exception:
        return None
    for frag, kind in _PULPROG_HINTS:
        if frag in pulprog:
            return kind
    return None


def read_delays(expno: str | Path) -> tuple[np.ndarray, str]:
    """vdlist (delays, s) or vclist (loop counts). Returns (values, source)."""
    expno = Path(expno)
    if (expno / "vdlist").exists():
        return read_vdlist(expno), "vdlist"
    if (expno / "vclist").exists():
        vals = [float(v) for v in (expno / "vclist").read_text().split() if v.strip()]
        return np.array(vals), "vclist"
    raise FileNotFoundError(f"no vdlist or vclist in {expno}")


# --------------------------------------------------------------------------
def model_curve(kind: str, p, t: np.ndarray, stretched: bool) -> np.ndarray:
    b = p["beta"].value if stretched else 1.0
    tau = p["tau"].value
    if kind == "satrec":
        return p["i0"] * (1.0 - np.exp(-((t / tau) ** b)))
    if kind == "invrec":
        return p["i0"] * (1.0 - 2.0 * p["f"] * np.exp(-((t / tau) ** b)))
    if kind in ("cpmg", "t1rho", "decay"):
        return p["i0"] * np.exp(-((t / tau) ** b))
    raise ValueError(f"unknown series kind {kind!r} (valid: {KINDS})")


def fit_series(t: np.ndarray, y: np.ndarray, kind: str = "satrec",
               stretched: bool = False):
    """Fit one build-up/decay curve. Returns (lmfit result, model callable)."""
    import lmfit

    t = np.asarray(t, float)
    y = np.asarray(y, float)
    scale = np.abs(y).max() or 1.0
    yn = y / scale

    params = lmfit.Parameters()
    params.add("i0", value=float(np.abs(yn).max()), min=0)
    pos = t > 0
    if kind in ("satrec", "invrec"):
        target = 0.63 * np.abs(yn).max()
        idx = int(np.argmin(np.abs(np.abs(yn) - target)))
    else:
        idx = int(np.argmin(np.abs(np.abs(yn) - 0.37 * np.abs(yn).max())))
    tau0 = float(t[idx]) or float(t[pos].min() if pos.any() else 1.0)
    params.add("tau", value=max(tau0, 1e-6), min=1e-9)
    if kind == "invrec":
        params.add("f", value=1.0, min=0.3, max=1.2)
    if stretched:
        params.add("beta", value=0.9, min=0.2, max=2.0)

    out = lmfit.minimize(
        lambda p: model_curve(kind, p, t, stretched) - yn,
        params, method="leastsq")
    return out, (lambda tt, r=out: model_curve(kind, r.params, np.asarray(tt),
                                               stretched) * scale)


@dataclass
class SeriesResult:
    kind: str
    x: np.ndarray                    # delays (s) or loop counts
    y: np.ndarray                    # integrals (normalized)
    tau: float                       # T1 / T2 / T1rho
    tau_err: float | None
    beta: float
    beta_err: float | None
    curve: object
    label: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def tau_name(self) -> str:
        return {"satrec": "T1", "invrec": "T1", "cpmg": "T2",
                "t1rho": "T1ρ"}.get(self.kind, "τ")

    @property
    def summary(self) -> str:
        err = f" ± {self.tau_err:.3g}" if self.tau_err else ""
        s = f"{self.tau_name} = {self.tau:.4g}{err} s"
        if self.beta != 1.0:
            b = f" ± {self.beta_err:.2f}" if self.beta_err else ""
            s += f", β = {self.beta:.2f}{b}"
        return (f"{self.label}: " if self.label else "") + s


def analyze(expno: str | Path, kind: str | None = None,
            window_ppm: tuple[float, float] | None = None,
            lb_hz: float = 100.0, stretched: bool = False,
            mode: str = "phase") -> SeriesResult:
    """Window-integration analysis of an arrayed EXPNO."""
    expno = Path(expno)
    kind = kind or detect_kind(expno) or "satrec"
    x, src = read_delays(expno)
    ppm, slices = process_slices(expno, lb_hz=lb_hz, mode=mode)
    n = min(len(x), slices.shape[0])
    x, slices = x[:n], slices[:n]

    notes = [f"delays from {src}", f"kind: {kind}"]
    if window_ppm is None:
        ref = slices[np.argmax(np.abs(slices).max(axis=1))]
        i = int(np.argmax(np.abs(ref)))
        half = np.abs(ref) > 0.5 * abs(ref[i])
        idx = np.where(half)[0]
        w = max(idx.max() - idx.min(), 5)
        lo_i, hi_i = max(0, idx.min() - w), min(len(ppm) - 1, idx.max() + w)
        window_ppm = (float(ppm[hi_i]), float(ppm[lo_i]))
        notes.append(f"auto window {window_ppm[0]:.1f} … {window_ppm[1]:.1f} ppm")
    hi, lo = max(window_ppm), min(window_ppm)
    sel = (ppm >= lo) & (ppm <= hi)
    y = slices[:, sel].sum(axis=1)
    ref_i = int(np.argmax(np.abs(y)))
    if y[ref_i] < 0:
        y = -y
        notes.append("integrals sign-flipped")
    y = y / (np.abs(y).max() or 1.0)

    out, curve = fit_series(x, y, kind=kind, stretched=stretched)
    p = out.params
    return SeriesResult(
        kind=kind, x=x, y=y,
        tau=float(p["tau"].value),
        tau_err=float(p["tau"].stderr) if p["tau"].stderr else None,
        beta=float(p["beta"].value) if stretched else 1.0,
        beta_err=(float(p["beta"].stderr)
                  if stretched and p["beta"].stderr else None),
        curve=curve, notes=notes)


# --------------------------------------------------------------------------
def fit_slice_amplitudes(recipe, ppm: np.ndarray, slices: np.ndarray,
                         window_ppm: tuple[float, float] | None = None,
                         ) -> np.ndarray:
    """Amplitude of every site in every slice, lineshapes held fixed.

    The model is linear in the amplitudes, so each slice is an ordinary
    non-negative least-squares problem on the fixed basis -- fast and free of
    local minima.
    """
    from scipy.optimize import nnls

    from larmor.engine import make_context, simulate_site

    ctx = make_context(recipe, exp_ppm=ppm)
    basis = []
    for site in recipe.sites:
        p = site.params["amplitude"]
        keep = p.value
        p.value = 1.0                      # unit-amplitude basis function
        basis.append(np.interp(ppm, ctx.x_ppm, simulate_site(site, ctx)))
        p.value = keep
    A = np.vstack(basis).T                 # (npts, nsites)

    if window_ppm:
        hi, lo = max(window_ppm), min(window_ppm)
        sel = (ppm >= lo) & (ppm <= hi)
    else:
        sel = np.ones(ppm.shape, bool)

    amps = np.empty((slices.shape[0], len(recipe.sites)))
    for k in range(slices.shape[0]):
        coef, _ = nnls(A[sel], slices[k][sel])
        amps[k] = coef
    return amps


def analyze_per_site(expno: str | Path, recipe, kind: str | None = None,
                     lb_hz: float = 100.0, stretched: bool = False,
                     mode: str = "phase",
                     window_ppm: tuple[float, float] | None = None,
                     ) -> list[SeriesResult]:
    """Relaxation time PER SITE: decompose every slice on the fitted
    lineshapes, then fit each site's amplitude series.

    This is what window integration cannot do -- two overlapping sites with
    different T1 come out separated.
    """
    expno = Path(expno)
    kind = kind or detect_kind(expno) or "satrec"
    x, src = read_delays(expno)
    ppm, slices = process_slices(expno, lb_hz=lb_hz, mode=mode)
    n = min(len(x), slices.shape[0])
    x, slices = x[:n], slices[:n]

    amps = fit_slice_amplitudes(recipe, ppm, slices, window_ppm=window_ppm)
    results = []
    for i, site in enumerate(recipe.sites):
        y = amps[:, i]
        norm = np.abs(y).max() or 1.0
        out, curve = fit_series(x, y / norm, kind=kind, stretched=stretched)
        p = out.params
        results.append(SeriesResult(
            kind=kind, x=x, y=y / norm,
            tau=float(p["tau"].value),
            tau_err=float(p["tau"].stderr) if p["tau"].stderr else None,
            beta=float(p["beta"].value) if stretched else 1.0,
            beta_err=(float(p["beta"].stderr)
                      if stretched and p["beta"].stderr else None),
            curve=curve, label=site.label or f"s{i}",
            notes=[f"per-site decomposition (NNLS on fixed lineshapes)",
                   f"delays from {src}"]))
    return results
