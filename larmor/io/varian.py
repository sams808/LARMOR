"""Read Varian / Agilent (VnmrJ) data — a ``.fid`` directory with ``procpar`` +
``fid`` (and ``log`` / ``text``).

The binary FID + parameter table are parsed with nmrglue; this module adds the
LARMOR-facing pieces: the observe nucleus (VnmrJ writes "Al27" → we return
"27Al"), the VnmrJ-referenced ppm axis (validated against real ²⁷Al glass data —
the tetrahedral-Al peak lands at ~50 ppm), and a default EM+FT+phase so opening a
dataset shows a spectrum immediately. Returns the same ``NMRData`` the Bruker
reader does, so the rest of the app is vendor-agnostic.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from larmor.io.bruker import Axis, NMRData


def _fid_dir(path) -> Path:
    p = Path(path)
    return p if p.is_dir() else p.parent


def is_varian(path) -> bool:
    """True for a Varian ``.fid`` directory (or a file inside one)."""
    d = _fid_dir(path)
    return (d / "procpar").exists() and (d / "fid").exists()


def varian_nucleus(tn: str) -> str:
    """VnmrJ isotope naming → LARMOR: 'Al27' → '27Al', 'Na23' → '23Na'."""
    m = re.match(r"\s*([A-Za-z]+)\s*(\d+)\s*", tn or "")
    return f"{m.group(2)}{m.group(1)}" if m else (tn or "")


def _pv(pp, key, default=None):
    try:
        return pp[key]["values"][0]
    except Exception:
        return default


def varian_ppm_axis(sfrq, sw, rfl, rfp, n) -> np.ndarray:
    """VnmrJ-referenced chemical-shift axis (ppm, high→low).

    The reference peak ``rfp`` (Hz) sits ``rfl`` Hz from the right edge, so the
    spectrum centre is at ``rfp + (sw/2 − rfl)`` Hz; the axis then spans ±sw/2
    about it, divided by the observe frequency. Validated on real ²⁷Al data.
    """
    sfrq, sw, rfl, rfp, n = float(sfrq), float(sw), float(rfl), float(rfp), int(n)
    centre_hz = rfp + (sw / 2.0 - rfl)
    hz = centre_hz + sw / 2.0 - np.arange(n) * (sw / max(n - 1, 1))
    return hz / (sfrq or 1.0)


def _autophase0(spec: np.ndarray) -> np.ndarray:
    """Coarse zeroth-order autophase: pick φ maximising the real signal over the
    peak region, return the real part (a first-look absorption spectrum)."""
    mag = np.abs(spec)
    mask = mag > 0.2 * (mag.max() or 1.0)
    if not mask.any():
        return spec.real
    phis = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    best = max(phis, key=lambda p: float(np.sum((spec * np.exp(1j * p)).real[mask])))
    return (spec * np.exp(1j * best)).real


def read(path) -> NMRData:
    """Read a Varian dataset as an ``NMRData`` (time-domain FID; 1D or arrayed)."""
    import nmrglue as ng

    d = _fid_dir(path)
    dic, data = ng.varian.read(str(d))
    pp = dic["procpar"]
    nucleus = varian_nucleus(_pv(pp, "tn", ""))
    sfrq = float(_pv(pp, "sfrq", 0.0) or 0.0)
    reffrq = float(_pv(pp, "reffrq", sfrq) or sfrq)
    sw = float(_pv(pp, "sw", 0.0) or 0.0)
    data = np.asarray(data)
    meta = {
        "nucleus": nucleus, "larmor_MHz": sfrq, "reffrq_MHz": reffrq,
        "sw_Hz": sw, "sr_hz": (sfrq - reffrq) * 1e6,
        "rfl": float(_pv(pp, "rfl", 0.0) or 0.0),
        "rfp": float(_pv(pp, "rfp", 0.0) or 0.0),
        "pulse_program": _pv(pp, "seqfil", ""),
        "title": (_pv(pp, "comment", "") or d.name),
        "vendor": "varian",
    }
    if data.ndim == 2 and data.shape[0] > 1:
        axes = [Axis(label="index", unit="point"),
                Axis(label=nucleus, unit="point", obs_MHz=sfrq, sw_Hz=sw)]
        return NMRData(ndim=2, domain="time", data=data.astype(complex),
                       axes=axes, meta=meta, is_pseudo2d=True, source=str(d),
                       warnings=["Varian arrayed data — open with the series tools"])
    fid = data.ravel().astype(complex)
    axes = [Axis(label=nucleus, unit="point", values=None, obs_MHz=sfrq, sw_Hz=sw)]
    return NMRData(ndim=1, domain="time", data=fid, axes=axes, meta=meta,
                   source=str(d))


def read_spectrum(path, lb_hz: float = 30.0):
    """Default processing (EM + FT + coarse phase) → ``(ppm, real_amp, meta)`` so a
    Varian dataset opens straight to a spectrum. Re-process/re-phase as needed."""
    nd = read(path)
    if nd.ndim != 1:
        raise ValueError("Varian dataset is arrayed/2D — use the series/2D tools")
    fid = np.asarray(nd.data, complex)
    sw, sfrq = nd.meta["sw_Hz"], nd.meta["larmor_MHz"]
    n = fid.size
    t = np.arange(n) / (sw or 1.0)
    fid = fid * np.exp(-np.pi * lb_hz * t)          # exponential apodisation
    fid = fid.copy(); fid[0] *= 0.5                 # first-point correction
    spec = np.fft.fftshift(np.fft.fft(fid))
    ppm = varian_ppm_axis(sfrq, sw, nd.meta["rfl"], nd.meta["rfp"], n)
    real = _autophase0(spec)
    order = np.argsort(ppm)
    return ppm[order], real[order], nd.meta
