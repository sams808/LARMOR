"""Saturation-recovery T1 analysis from Bruker pseudo-2D acquisitions.

Reads ser + vdlist (read-only), processes every slice (EM -> ZF -> FT, phased
on the longest-delay slice or magnitude mode), integrates a ppm window, and
fits I(t) = I0 * (1 - exp(-(t/T1)^beta)) -- beta fixed to 1 unless stretched.

The per-slice integrals can be cross-checked against TopSpin's own t1ints.txt
when present (the test suite does exactly that on real data). ``read_ct1t2``
reads TopSpin's own per-integral T1 result (both the legacy and the SIMFIT
layout) for the quantitativity chip.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class SatrecResult:
    delays_s: np.ndarray
    integrals: np.ndarray            # normalized to max = 1
    t1_s: float
    t1_err: float | None
    beta: float
    beta_err: float | None
    i0: float
    window_ppm: tuple[float, float]
    x_ppm: np.ndarray                # slice axis (for display)
    slices: np.ndarray               # processed real spectra, shape (nvd, npts)
    curve: object = None             # callable model(t)
    notes: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        err = f" ± {self.t1_err:.3g}" if self.t1_err else ""
        b = f", β = {self.beta:.2f}" + (f" ± {self.beta_err:.2f}"
                                        if self.beta_err else "")
        return f"T1 = {self.t1_s:.4g}{err} s" + (b if self.beta != 1.0 else "")


def _seconds(tok: str) -> float:
    """A TopSpin duration token in seconds: 'm' → 1e-3, 'u' → 1e-6, 'n' →
    1e-9 (TopSpin writes a zero delay as '0.000n'), 's' or a bare number →
    seconds. Raises ValueError on anything else (callers decide)."""
    tok = tok.strip().lower()
    mult = 1.0
    if tok.endswith("m"):
        mult, tok = 1e-3, tok[:-1]
    elif tok.endswith("u"):
        mult, tok = 1e-6, tok[:-1]
    elif tok.endswith("n"):
        mult, tok = 1e-9, tok[:-1]
    elif tok.endswith("s"):
        tok = tok[:-1]
    return float(tok) * mult


def read_vdlist(expno: str | Path) -> np.ndarray:
    """Parse Bruker vdlist (supports m/s/u suffixes)."""
    text = (Path(expno) / "vdlist").read_text().split()
    out = []
    for tok in text:
        tok = tok.strip().lower()
        if not tok:
            continue
        out.append(_seconds(tok))
    return np.array(out)


@dataclass(frozen=True)
class T1Region:
    """One integral region (or peak point) of a TopSpin t1/t2 analysis and
    the T1 it fitted there. ``hi_ppm == lo_ppm`` for a 'Peak Point' block;
    ``components_s`` keeps every SIMFIT component, ``t1_s`` the longest."""

    index: int
    hi_ppm: float
    lo_ppm: float
    t1_s: float
    components_s: tuple = ()
    sd: float | None = None
    npoints: int | None = None


_CT1T2_BLOCK_RE = re.compile(
    r"(\d+) points for (?:Integral|Peak) (\d+),\s+"
    r"(?:Integral Region from (\S+) to (\S+) ppm|Peak Point at (\S+) ppm)")
_CT1T2_T1_RE = re.compile(r"^\s*T1\s*=\s*(.+)$")
_CT1T2_SD_RE = re.compile(r"^\s*SD\s*=\s*(\S+)")


def read_ct1t2(pdata_dir: str | Path) -> list[T1Region]:
    """TopSpin's per-integral T1 result from ``pdata/N/ct1t2.txt``.

    Understands the legacy layout ('AREA fit : I[t]=I[0]+P*exp(-t/T1)', one
    'T1 = 4.643s' line per '7 points for Integral 1, Integral Region from
    27.382 to 5.846 ppm' block), the SIMFIT layout ('T1 = 60.764m 1.759s
    4.539s 4.466s' — every component parsed, the longest governs recovery)
    and the INTENSITY 'Peak Point at -68.561 ppm' variant. Returns [] for a
    missing file, a file without a T1 line or garbage; never raises. Accepts
    the pdata folder or the file itself.
    """
    p = Path(pdata_dir)
    if p.is_dir():
        p = p / "ct1t2.txt"
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except (OSError, ValueError):
        return []
    regions: list[T1Region] = []
    cur: dict | None = None

    def _flush():
        if cur is None or not cur.get("components"):
            return
        comps = tuple(cur["components"])
        regions.append(T1Region(
            index=cur["index"], hi_ppm=cur["hi"], lo_ppm=cur["lo"],
            t1_s=max(comps), components_s=comps, sd=cur.get("sd"),
            npoints=cur.get("npoints")))

    try:
        for line in text.splitlines():
            m = _CT1T2_BLOCK_RE.search(line)
            if m:
                _flush()
                npts, idx, a, b, pk = m.groups()
                try:
                    if pk is not None:
                        hi = lo = float(pk)
                    else:
                        hi, lo = max(float(a), float(b)), min(float(a), float(b))
                except ValueError:
                    cur = None
                    continue
                cur = {"index": int(idx), "npoints": int(npts), "hi": hi,
                       "lo": lo, "components": [], "sd": None}
                continue
            if cur is None:
                continue
            m = _CT1T2_T1_RE.match(line)
            if m:
                comps = []
                for tok in m.group(1).split():
                    try:
                        comps.append(_seconds(tok))
                    except ValueError:
                        continue
                if comps and not cur["components"]:
                    cur["components"] = comps
                continue
            m = _CT1T2_SD_RE.match(line)
            if m and cur.get("sd") is None:
                try:
                    cur["sd"] = float(m.group(1))
                except ValueError:
                    pass
        _flush()
    except Exception:
        return []
    return regions


def process_slices(expno: str | Path, lb_hz: float = 100.0,
                   zf_factor: int = 2, mode: str = "phase",
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Read ser (read-only), process every row identically.

    mode "phase": EM, ZF, FT, then apply the p0 that maximizes the real
    integral of the LAST (fully relaxed) slice to all slices.
    mode "magnitude": phase-insensitive |S|.
    Returns (x_ppm ascending, slices_real) with slices shape (nvd, npts).

    The ser is read through the universal reader, so it is correctly reshaped
    into rows (nmrglue leaves it flat) and the digital filter is removed.
    """
    from larmor.io import bruker

    # read the RAW arrayed data (ser/fid), never the processed 2rr: an EXPNO
    # folder resolves to its processed pdata by default, so point at the ser
    p = Path(expno)
    if p.is_dir():
        if (p / "ser").exists():
            p = p / "ser"
        elif (p / "fid").exists():
            p = p / "fid"
    data = bruker.read(p)              # 2D time-domain NMRData, rows = slices
    if data.domain != "time":
        raise ValueError("process_slices needs the raw ser/fid, but got a "
                         "processed spectrum")
    if data.ndim == 1:
        ser = data.data[None, :]
    else:
        ser = data.data
    sw = data.axes[-1].sw_Hz or data.meta.get("sw_Hz", 0.0)
    sfo1 = data.meta["larmor_MHz"]

    ser = np.asarray(ser, complex).copy()
    ser[:, 0] *= 0.5                     # FCOR: kill the DC ridge / edge spike
    n = ser.shape[1]
    t = np.arange(n) / sw
    window = np.exp(-np.pi * lb_hz * t)
    nfft = int(2 ** np.ceil(np.log2(n * max(1, zf_factor))))
    spec = np.fft.fftshift(np.fft.fft(ser * window, n=nfft, axis=1), axes=1)
    # ascending frequency axis (fftshift(fftfreq), matches the rest of LARMOR)
    freq = np.fft.fftshift(np.fft.fftfreq(nfft, d=1.0 / sw))
    x_ppm = freq / sfo1
    order = np.argsort(x_ppm)
    x_ppm, spec = x_ppm[order], spec[:, order]

    if mode == "magnitude":
        out = np.abs(spec)
    else:
        # phase on the MAX-SIGNAL slice, not blindly the last row. A Bruker ser
        # often carries a trailing dummy/setup acquisition (acqu2s TD1 exceeds
        # len(vdlist)) whose phase is unrelated; using it as the reference would
        # corrupt the phase -- and therefore the build-up -- of every slice.
        ref = spec[int(np.argmax(np.abs(spec).sum(axis=1)))]
        phis = np.linspace(-np.pi, np.pi, 721)
        scores = [(np.real(ref * np.exp(1j * p))).sum() for p in phis]
        p0 = float(phis[int(np.argmax(scores))])
        out = np.real(spec * np.exp(1j * p0))
    return x_ppm, out


def fit_t1(delays: np.ndarray, integrals: np.ndarray,
           stretched: bool = False):
    """Fit the saturation-recovery build-up. Returns an lmfit result and the
    model callable."""
    import lmfit

    y = integrals / np.abs(integrals).max()

    def model(p, t):
        beta = p["beta"].value if stretched else 1.0
        return p["i0"] * (1.0 - np.exp(-((t / p["t1"]) ** beta)))

    params = lmfit.Parameters()
    params.add("i0", value=float(y.max()), min=0)
    half = delays[np.searchsorted(y, 0.63 * y.max())] if y.max() > 0 else 1.0
    params.add("t1", value=float(max(half, delays[delays > 0].min())), min=1e-6)
    if stretched:
        params.add("beta", value=0.8, min=0.2, max=2.0)

    sel = delays >= 0
    out = lmfit.minimize(lambda p: model(p, delays[sel]) - y[sel], params,
                         method="leastsq")
    return out, (lambda tt, r=out: (
        r.params["i0"].value *
        (1 - np.exp(-((tt / r.params["t1"].value) **
                      (r.params["beta"].value if stretched else 1.0))))))


def analyze(expno: str | Path, window_ppm: tuple[float, float] | None = None,
            lb_hz: float = 100.0, stretched: bool = False,
            mode: str = "phase") -> SatrecResult:
    """Full pipeline: vdlist + ser -> per-slice integrals -> T1."""
    expno = Path(expno)
    delays = read_vdlist(expno)
    x, slices = process_slices(expno, lb_hz=lb_hz, mode=mode)
    if slices.shape[0] != delays.size:
        n = min(slices.shape[0], delays.size)
        delays, slices = delays[:n], slices[:n]

    notes = []
    if window_ppm is None:
        # auto-window: FWHM region of the fully relaxed slice, widened 3x
        ref = slices[-1]
        imax = int(np.argmax(np.abs(ref)))
        half = np.abs(ref) > 0.5 * np.abs(ref[imax])
        idx = np.where(half)[0]
        lo_i, hi_i = idx.min(), idx.max()
        width = max(hi_i - lo_i, 5)
        lo_i = max(0, lo_i - width)
        hi_i = min(len(x) - 1, hi_i + width)
        window_ppm = (float(x[hi_i]), float(x[lo_i]))
        notes.append(f"auto integration window {window_ppm[0]:.1f} … "
                     f"{window_ppm[1]:.1f} ppm")
    hi, lo = max(window_ppm), min(window_ppm)
    sel = (x >= lo) & (x <= hi)
    integrals = slices[:, sel].sum(axis=1)
    if integrals[-1] < 0:            # sign convention: relaxed slice positive
        integrals = -integrals
        notes.append("integrals sign-flipped (negative phasing)")
    norm = np.abs(integrals).max() or 1.0
    integrals = integrals / norm

    out, curve = fit_t1(delays, integrals, stretched=stretched)
    p = out.params
    return SatrecResult(
        delays_s=delays, integrals=integrals,
        t1_s=float(p["t1"].value),
        t1_err=float(p["t1"].stderr) if p["t1"].stderr else None,
        beta=float(p["beta"].value) if stretched else 1.0,
        beta_err=(float(p["beta"].stderr)
                  if stretched and p["beta"].stderr else None),
        i0=float(p["i0"].value),
        window_ppm=(hi, lo), x_ppm=x, slices=slices, curve=curve, notes=notes)
