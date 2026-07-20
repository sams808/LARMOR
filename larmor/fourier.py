"""Fourier processing for raw data, 1D and 2D, ssNake-style.

1D time-domain processing (apodize / zero-fill / FT / phase / baseline) already
lives in ``larmor.processing``. This module adds:

  * ``ft1d``   -- convenience: run a processing pipeline on a 1D FID.
  * the 2D machinery ssNake exposes and TopSpin's xfb hides: hypercomplex
    recombination of the indirect dimension according to the acquisition mode
    (States / States-TPPI / Echo-Antiecho / TPPI / QF), then an F2 then F1
    transform with independent apodization and phasing per dimension.

The indirect-dimension quadrature is the part beginners get wrong, so the mode
is explicit and defaults are taken from the Bruker FnMODE when a dataset is
read through ``larmor.io.bruker``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from larmor import processing as proc


# --------------------------------------------------------------------------
def ft1d(fid: np.ndarray, sw_Hz: float, sfo1_MHz: float,
         ops: list[dict] | None = None, offset_ppm: float = 0.0,
         ) -> tuple[np.ndarray, np.ndarray]:
    """Process one FID to a spectrum. Returns (ppm, complex spectrum).

    ops is a processing pipeline (see larmor.processing.OPS); if it contains no
    'ft' step one is appended, so passing e.g. [{"op":"em","lb_hz":100}] just
    works.
    """
    s = proc.Spectrum1D(x_ppm=None, y=np.asarray(fid, complex),
                        sfo1_MHz=sfo1_MHz, sw_Hz=sw_Hz, domain="time")
    ops = list(ops or [])
    if not any(o.get("op") == "ft" for o in ops):
        ops = ops + [{"op": "ft", "offset_ppm": offset_ppm}]
    s = proc.apply(s, ops)
    return s.x_ppm, s.y


# --------------------------------------------------------------------------
def states_recombine(ser: np.ndarray, mode: str) -> np.ndarray:
    """Recombine the indirect quadrature into a hypercomplex t1 series.

    Bruker stores the two indirect components on ALTERNATE rows. This returns a
    complex array whose real/imag along t1 are the cosine/sine-modulated data,
    ready for a complex F1 transform.

      States / States-TPPI : rows (cos, sin, cos, sin, ...) -> cos + i sin
      Echo-Antiecho        : rows (echo, antiecho) -> (echo+anti)/2, (echo-anti)/2i
      QF / undefined       : single component; nothing to recombine
    """
    mode = (mode or "QF").lower()
    if mode in ("qf", "undefined", "qseq", ""):
        return ser
    n1 = ser.shape[0] // 2
    even, odd = ser[0:2 * n1:2], ser[1:2 * n1:2]
    if mode in ("states", "states-tppi"):
        return even + 1j * odd
    if mode == "echo-antiecho":
        a = 0.5 * (even + odd)
        b = 0.5 * (even - odd) / 1j
        return a + 1j * b
    if mode == "tppi":
        # TPPI is a real F1 series; keep as-is, the F1 FT is a real transform
        return ser
    return ser


def _apodize_axis(data: np.ndarray, sw_Hz: float, ops: list[dict],
                  axis: int) -> np.ndarray:
    """Apply time-domain apodization ops along one axis via processing.py."""
    if not ops:
        return data
    out = np.empty_like(data)
    it = np.moveaxis(data, axis, 0)
    res = np.moveaxis(out, axis, 0)
    n = it.shape[0]
    dummy = proc.Spectrum1D(x_ppm=None, y=np.zeros(n, complex),
                            sfo1_MHz=1.0, sw_Hz=sw_Hz, domain="time")
    for idx in np.ndindex(it.shape[1:]):
        dummy.y = it[(slice(None),) + idx].astype(complex)
        dummy.domain = "time"
        proc.apply(dummy, ops)
        res[(slice(None),) + idx] = dummy.y
    return out


@dataclass
class FT2DParams:
    f2_ops: list = field(default_factory=list)     # apodization on the direct dim
    f1_ops: list = field(default_factory=list)     # apodization on the indirect dim
    f2_zf: int = 1
    f1_zf: int = 1
    mode: str = "States"                           # indirect quadrature mode
    f2_p0: float = 0.0
    f2_p1: float = 0.0
    f1_p0: float = 0.0
    f1_p1: float = 0.0


def ft2d(ser: np.ndarray, sw2_Hz: float, sw1_Hz: float,
         sfo2_MHz: float, sfo1_MHz: float, params: FT2DParams | None = None,
         offset2_ppm: float = 0.0, offset1_ppm: float = 0.0,
         ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Full 2D transform of a raw ser. Returns (f2_ppm, f1_ppm, complex Z).

    Order: recombine indirect quadrature -> apodize+ZF+FT the direct (F2)
    dimension -> apodize+ZF+FT the indirect (F1) dimension -> phase each.
    """
    p = params or FT2DParams()
    ser = np.asarray(ser, complex)

    # 1) indirect quadrature -> hypercomplex t1
    hyper = states_recombine(ser, p.mode)

    # 2) direct dimension (F2): apodize, zero-fill, FT
    z = _apodize_axis(hyper, sw2_Hz, p.f2_ops, axis=1)
    n2 = int(2 ** np.ceil(np.log2(z.shape[1] * max(1, p.f2_zf))))
    z = np.pad(z, ((0, 0), (0, n2 - z.shape[1])))
    z = np.fft.fftshift(np.fft.fft(z, axis=1), axes=1)
    if p.f2_p0 or p.f2_p1:
        z = _phase_axis(z, p.f2_p0, p.f2_p1, axis=1)

    # 3) indirect dimension (F1): apodize, zero-fill, FT
    z = _apodize_axis(z, sw1_Hz, p.f1_ops, axis=0)
    n1 = int(2 ** np.ceil(np.log2(z.shape[0] * max(1, p.f1_zf))))
    z = np.pad(z, ((0, n1 - z.shape[0]), (0, 0)))
    z = np.fft.fftshift(np.fft.fft(z, axis=0), axes=0)
    if p.f1_p0 or p.f1_p1:
        z = _phase_axis(z, p.f1_p0, p.f1_p1, axis=0)

    f2 = _axis_ppm(sw2_Hz, sfo2_MHz, z.shape[1], offset2_ppm)
    f1 = _axis_ppm(sw1_Hz, sfo1_MHz, z.shape[0], offset1_ppm)
    o2, o1 = np.argsort(f2), np.argsort(f1)
    return f2[o2], f1[o1], z[np.ix_(o1, o2)]


def _phase_axis(z: np.ndarray, p0: float, p1: float, axis: int) -> np.ndarray:
    n = z.shape[axis]
    idx = np.arange(n) / max(n - 1, 1)
    ph = np.exp(1j * np.deg2rad(p0 + p1 * (idx - 0.5)))
    shape = [1, 1]
    shape[axis] = n
    return z * ph.reshape(shape)


def _axis_ppm(sw_Hz: float, sfo_MHz: float, n: int,
              offset_ppm: float) -> np.ndarray:
    """Frequency axis matching fftshift(fft): ascending, so a component at
    +f Hz lands at +f. Display inverts the ppm axis separately."""
    if sfo_MHz <= 0:
        return np.arange(n, dtype=float)
    freq = np.fft.fftshift(np.fft.fftfreq(n, d=1.0 / sw_Hz))
    return freq / sfo_MHz + offset_ppm


def ft2d_from_nmrdata(data, params: FT2DParams | None = None):
    """Transform a raw 2D NMRData (from larmor.io.bruker) into a Data2D."""
    from larmor.twod import Data2D

    if data.ndim != 2 or data.domain != "time":
        raise ValueError("ft2d_from_nmrdata needs a raw 2D FID")
    f1_axis, f2_axis = data.axes
    p = params or FT2DParams(mode=data.meta.get("fnmode", "States"))
    if data.is_pseudo2d:
        # arrayed experiment: only the direct dimension is spectroscopic
        rows = []
        for k in range(data.data.shape[0]):
            ppm, spec = ft1d(data.data[k], f2_axis.sw_Hz,
                             f2_axis.obs_MHz or data.meta["larmor_MHz"],
                             ops=p.f2_ops)
            rows.append(spec.real)
        z = np.array(rows)
        f1_vals = (f1_axis.values if f1_axis.values is not None
                   else np.arange(z.shape[0], dtype=float))
        return Data2D(f2_ppm=ppm, f1_ppm=f1_vals, z=z,
                      nucleus=data.nucleus, larmor_MHz=data.meta["larmor_MHz"],
                      source=data.source,
                      notes=["pseudo-2D: F1 kept as the arrayed axis"] +
                            list(data.warnings))
    f2, f1, z = ft2d(
        data.data, f2_axis.sw_Hz, f1_axis.sw_Hz,
        f2_axis.obs_MHz or data.meta["larmor_MHz"],
        f1_axis.obs_MHz or data.meta["larmor_MHz"], params=p)
    return Data2D(f2_ppm=f2, f1_ppm=f1, z=z.real, nucleus=data.nucleus,
                  larmor_MHz=data.meta["larmor_MHz"], source=data.source,
                  notes=list(data.warnings))
