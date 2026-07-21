"""Plain-text spectrum I/O — a simple (ppm, intensity) CSV that carries the
experiment metadata in a comment header, so a spectrum LARMOR produces (e.g. a
background subtraction) round-trips straight back into the fit workbench.

Format::

    # LARMOR spectrum
    # nucleus=27Al
    # larmor_MHz=195.483000
    # spin_rate_Hz=20000.0
    # sample=CaAlGlass minus empty rotor
    ppm,intensity
    150.0,0.0123
    ...

Comment lines (``#``) are optional; the data may be comma- or whitespace-
separated, so ordinary two-column exports load too.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

_FLOAT_KEYS = {"larmor_mhz": "larmor_MHz", "spin_rate_hz": "spin_rate_Hz"}


def write_csv(path: str | Path, ppm: np.ndarray, amp: np.ndarray,
              meta: dict | None = None) -> str:
    """Write a (ppm, intensity) CSV with a metadata header. Returns the path."""
    meta = meta or {}
    ppm = np.asarray(ppm, float)
    amp = np.asarray(amp, float)
    order = np.argsort(ppm)[::-1]            # descending ppm, NMR convention
    lines = ["# LARMOR spectrum"]
    for key in ("nucleus", "larmor_MHz", "spin_rate_Hz", "sample"):
        if meta.get(key) not in (None, "", 0) or key == "nucleus":
            lines.append(f"# {key}={meta.get(key, '')}")
    lines.append("ppm,intensity")
    for x, y in zip(ppm[order], amp[order]):
        lines.append(f"{x:.6f},{y:.8g}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def read_csv(path: str | Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Read a (ppm, intensity) text/CSV spectrum. Returns (ppm, amp, meta)."""
    meta: dict = {}
    xs, ys = [], []
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            body = line.lstrip("#").strip()
            if "=" in body:
                k, val = body.split("=", 1)
                k = k.strip().lower()
                key = _FLOAT_KEYS.get(k, k)
                if k in _FLOAT_KEYS:
                    try:
                        meta[key] = float(val)
                    except ValueError:
                        pass
                else:
                    meta[key] = val.strip()
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            continue                          # header row like "ppm intensity"
        xs.append(x); ys.append(y)
    if len(xs) < 2:
        raise ValueError(f"no 2-column numeric data found in {path}")
    ppm = np.asarray(xs); amp = np.asarray(ys)
    order = np.argsort(ppm)
    return ppm[order], amp[order], meta


def subtract(sample_ppm: np.ndarray, sample_amp: np.ndarray,
             bg_ppm: np.ndarray, bg_amp: np.ndarray,
             scale: float = 1.0, shift_ppm: float = 0.0) -> np.ndarray:
    """sample − scale·background, with the background interpolated onto the
    sample axis (and optionally shifted). Returns the difference amplitude."""
    bp = np.asarray(bg_ppm, float) + shift_ppm
    order = np.argsort(bp)
    bg_on = np.interp(np.asarray(sample_ppm, float), bp[order],
                      np.asarray(bg_amp, float)[order], left=0.0, right=0.0)
    return np.asarray(sample_amp, float) - scale * bg_on


def best_scale(sample_ppm, sample_amp, bg_ppm, bg_amp, shift_ppm: float = 0.0,
               window=None) -> float:
    """Least-squares background scale: argmin ||sample − s·bg||² over a window
    (default: the whole overlap). Clamped to be non-negative."""
    bp = np.asarray(bg_ppm, float) + shift_ppm
    order = np.argsort(bp)
    bg_on = np.interp(np.asarray(sample_ppm, float), bp[order],
                      np.asarray(bg_amp, float)[order], left=0.0, right=0.0)
    s = np.asarray(sample_amp, float)
    sel = np.ones(s.shape, bool)
    if window:
        hi, lo = max(window), min(window)
        sel = (np.asarray(sample_ppm) >= lo) & (np.asarray(sample_ppm) <= hi)
    denom = float(bg_on[sel] @ bg_on[sel])
    if denom <= 0:
        return 1.0
    return max(0.0, float(s[sel] @ bg_on[sel]) / denom)
