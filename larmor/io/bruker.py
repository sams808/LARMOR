"""Universal read-only Bruker TopSpin reader built on nmrglue.

Accepts anything a spectroscopist would drag in -- a processed file
(``pdata/1/1r`` or ``2rr``), a raw ``fid``/``ser``, a ``pdata/N`` folder, or an
EXPNO folder -- and figures out for itself what it is: 1D or 2D, raw
(time-domain) or processed (frequency-domain), a genuine spectroscopic 2D or a
pseudo-2D arrayed experiment (relaxation, REDOR) whose indirect dimension is a
list of delays rather than a chemical-shift axis.

Design guarantees kept from the first version:
  - never writes into instrument folders; ``snapshot`` / ``verify_untouched``
    make it checkable and ``read_expno(..., verify=True)`` enforces it;
  - surfaces metadata conflicts (acqus MASR vs the operator-typed title).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import nmrglue as ng

#: Bruker FnMODE codes for the indirect dimension
FNMODE = {0: "undefined", 1: "QF", 2: "QSEQ", 3: "TPPI", 4: "States",
          5: "States-TPPI", 6: "Echo-Antiecho"}


# --------------------------------------------------------------------------
def snapshot(root: Path) -> dict[str, tuple[int, int]]:
    return {
        str(p): (p.stat().st_mtime_ns, p.stat().st_size)
        for p in sorted(Path(root).rglob("*")) if p.is_file()
    }


def verify_untouched(root: Path, before: dict) -> None:
    after = snapshot(root)
    if after != before:
        changed = {k for k in set(before) | set(after)
                   if before.get(k) != after.get(k)}
        raise RuntimeError(f"read-only violation in {root}: {changed}")


# --------------------------------------------------------------------------
@dataclass
class Axis:
    """One dimension of an NMR dataset."""

    label: str                        # isotope ("27Al") or role ("delay")
    unit: str                         # "ppm" | "Hz" | "s" | "point"
    values: np.ndarray | None = None  # coordinates; None for un-FT'd time data
    obs_MHz: float = 0.0              # observe frequency (SFO/SF)
    sw_Hz: float = 0.0


@dataclass
class NMRData:
    """A read Bruker dataset, 1D or 2D, time- or frequency-domain."""

    ndim: int
    domain: str                       # "time" | "freq"
    data: np.ndarray                  # 1D (n,) or 2D (n1, n2); complex in time
    axes: list                        # innermost (direct/F2) LAST, like numpy
    meta: dict = field(default_factory=dict)
    is_pseudo2d: bool = False
    warnings: list = field(default_factory=list)
    source: str = ""
    #: for a processed 2D, the imaginary quadrants {"ri","ir","ii"} oriented
    #: like ``data`` (rr), enabling exact hypercomplex phase correction
    hyper: dict | None = None

    @property
    def nucleus(self) -> str:
        return self.meta.get("nucleus", "")

    @property
    def summary(self) -> str:
        kind = {1: "1D", 2: "pseudo-2D" if self.is_pseudo2d else "2D"}[self.ndim]
        dom = "FID" if self.domain == "time" else "spectrum"
        bits = [f"{kind} {dom}", self.meta.get("nucleus", "?")]
        if self.meta.get("larmor_MHz"):
            bits.append(f"{self.meta['larmor_MHz']:.3f} MHz")
        if self.meta.get("pulse_program"):
            bits.append(self.meta["pulse_program"])
        return " · ".join(bits)


# --------------------------------------------------------------------------
@dataclass
class BrukerRef:
    """Where a path points, once resolved."""

    expno: Path
    procno: int
    target: str                       # "1r" | "2rr" | "fid" | "ser"
    ndim: int


def is_expno(path: str | Path) -> bool:
    p = Path(path)
    return p.is_dir() and (p / "acqus").exists()


def _find_expno(start: Path) -> Path:
    """Walk up until a folder containing acqus is found."""
    p = start if start.is_dir() else start.parent
    for cand in [p, *p.parents]:
        if (cand / "acqus").exists():
            return cand
    raise FileNotFoundError(f"no EXPNO (acqus) above {start}")


def resolve(path: str | Path) -> BrukerRef:
    """Turn any Bruker-ish path into (expno, procno, what-to-read, ndim)."""
    p = Path(path)
    is2d_expno = lambda e: (e / "acqu2s").exists() or (e / "ser").exists()

    # a processed data file, e.g. .../pdata/1/1r or .../pdata/1/2rr
    if p.is_file() and p.name in ("1r", "1i", "2rr", "2ri", "2ir", "2ii"):
        procno = int(p.parent.name)
        expno = _find_expno(p.parent)
        target = "2rr" if p.name.startswith("2") else "1r"
        return BrukerRef(expno, procno, target, 2 if target == "2rr" else 1)

    # a raw data file
    if p.is_file() and p.name in ("fid", "ser"):
        expno = _find_expno(p)
        return BrukerRef(expno, 1, p.name, 2 if p.name == "ser" else 1)

    # a pdata/N folder
    if p.is_dir() and p.parent.name == "pdata":
        procno = int(p.name)
        expno = _find_expno(p)
        target = "2rr" if (p / "2rr").exists() else "1r"
        return BrukerRef(expno, procno, target, 2 if target == "2rr" else 1)

    # an EXPNO folder: prefer processed pdata/1 if present, else raw
    if is_expno(p):
        if (p / "pdata" / "1" / "2rr").exists():
            return BrukerRef(p, 1, "2rr", 2)
        if (p / "pdata" / "1" / "1r").exists():
            return BrukerRef(p, 1, "1r", 1)
        if (p / "ser").exists():
            return BrukerRef(p, 1, "ser", 2)
        if (p / "fid").exists():
            return BrukerRef(p, 1, "fid", 1)
    raise ValueError(f"not a recognizable Bruker path: {path}")


# --------------------------------------------------------------------------
def _ppm_axis(procs: dict, npts: int) -> np.ndarray:
    """TopSpin processed axis: ppm[i] = OFFSET - i * SW_p / SF / SI.

    (SW_p in Hz, SF in MHz -> SW_p/SF is the span in ppm.)
    """
    offset = float(procs["OFFSET"])
    sw_p, sf, si = float(procs["SW_p"]), float(procs["SF"]), int(procs["SI"])
    return offset - np.arange(npts) * (sw_p / sf / si)


def _isotope(dic: dict, key: str = "acqus") -> str:
    return str(dic.get(key, {}).get("NUC1", "")).strip()


def _read_delays(expno: Path) -> tuple[np.ndarray | None, str]:
    from larmor.satrec import read_vdlist

    if (expno / "vdlist").exists():
        return read_vdlist(expno), "vdlist (s)"
    if (expno / "vclist").exists():
        vals = [float(v) for v in (expno / "vclist").read_text().split()
                if v.strip()]
        return np.array(vals), "vclist"
    return None, ""


def _pseudo2d(expno: Path, proc2s: dict) -> bool:
    """A real spectroscopic F1 has a sane ppm OFFSET and no delay list."""
    if (expno / "vdlist").exists() or (expno / "vclist").exists():
        return True
    try:
        offset = abs(float(proc2s.get("OFFSET", 0.0)))
    except (TypeError, ValueError):
        return True
    return offset > 100000.0          # absurd F1 offset => not FT'd as ppm


def read(path: str | Path, verify: bool = True) -> NMRData:
    """Read any Bruker path into a unified NMRData (read-only)."""
    ref = resolve(path)
    before = snapshot(ref.expno) if verify else None
    try:
        if ref.target == "1r":
            out = _read_1r(ref)
        elif ref.target == "2rr":
            out = _read_2rr(ref)
        elif ref.target == "fid":
            out = _read_fid(ref)
        else:
            out = _read_ser(ref)
    finally:
        if verify:
            verify_untouched(ref.expno, before)
    out.source = str(path)
    return out


def _read_1r(ref: BrukerRef) -> NMRData:
    pdata = ref.expno / "pdata" / str(ref.procno)
    dic, real = ng.bruker.read_pdata(str(pdata))
    procs = dic["procs"]
    ppm = _ppm_axis(procs, real.size)
    acqus = dic.get("acqus", {})
    order = np.argsort(ppm)
    title = _read_title(pdata)
    meta = _meta_1d(acqus, title, ref.expno)
    ax = Axis(meta["nucleus"], "ppm", ppm[order],
              obs_MHz=float(procs.get("SF", 0.0)), sw_Hz=meta.get("sw_Hz", 0.0))
    return NMRData(ndim=1, domain="freq", data=np.asarray(real, float)[order],
                   axes=[ax], meta=meta, warnings=_conflicts(meta, title))


def _read_2rr(ref: BrukerRef) -> NMRData:
    pdata = ref.expno / "pdata" / str(ref.procno)
    hyper_raw = None
    z = None
    dic = None
    # read all four quadrants when they genuinely exist, so phasing can be exact
    try:
        dic2, comps = ng.bruker.read_pdata(str(pdata), all_components=True)
        if (isinstance(comps, (list, tuple)) and len(comps) == 4
                and all(np.ndim(c) == 2 for c in comps)):
            dic = dic2
            z = np.asarray(comps[0], float)
            hyper_raw = [np.asarray(comps[1], float),   # ri
                         np.asarray(comps[2], float),   # ir
                         np.asarray(comps[3], float)]   # ii
    except Exception:
        pass
    if z is None:                          # not hypercomplex: plain real read
        dic, z = ng.bruker.read_pdata(str(pdata))
        z = np.asarray(z, float)
    procs, proc2s = dic["procs"], dic["proc2s"]
    f2 = _ppm_axis(procs, z.shape[1])
    pseudo = _pseudo2d(ref.expno, proc2s)
    warnings = []
    if pseudo:
        delays, src = _read_delays(ref.expno)
        if delays is not None and len(delays) >= z.shape[0]:
            f1_vals, f1_unit, f1_label = delays[:z.shape[0]], "s", "delay"
        else:
            f1_vals, f1_unit, f1_label = np.arange(z.shape[0]), "point", "index"
        warnings.append("indirect dimension is an arrayed/relaxation axis, "
                        "not a chemical shift")
    else:
        f1_vals = _ppm_axis(proc2s, z.shape[0])
        f1_unit, f1_label = "ppm", _isotope(dic, "acqu2s") or "F1"
    o2, o1 = np.argsort(f2), np.argsort(f1_vals)
    ix = np.ix_(o1, o2)
    z = np.asarray(z, float)[ix]
    hyper = None
    if hyper_raw is not None and not pseudo:
        hyper = {"ri": hyper_raw[0][ix], "ir": hyper_raw[1][ix],
                 "ii": hyper_raw[2][ix]}
    acqus = dic.get("acqus", {})
    meta = _meta_1d(acqus, _read_title(pdata), ref.expno)
    meta["fnmode"] = FNMODE.get(int(proc2s.get("MC2", 0) or 0), "?")
    axes = [Axis(f1_label, f1_unit, f1_vals[o1],
                 obs_MHz=float(proc2s.get("SF", 0.0))),
            Axis(meta["nucleus"], "ppm", f2[o2],
                 obs_MHz=float(procs.get("SF", 0.0)))]
    return NMRData(ndim=2, domain="freq", data=z, axes=axes, meta=meta,
                   is_pseudo2d=pseudo, warnings=warnings, hyper=hyper)


def _read_fid(ref: BrukerRef) -> NMRData:
    dic, fid = ng.bruker.read(str(ref.expno))
    fid = ng.bruker.remove_digital_filter(dic, fid)
    acqus = dic["acqus"]
    meta = _meta_1d(acqus, _read_title(ref.expno / "pdata" / "1"), ref.expno)
    ax = Axis(meta["nucleus"], "s", None, obs_MHz=meta["larmor_MHz"],
              sw_Hz=meta["sw_Hz"])
    return NMRData(ndim=1, domain="time", data=fid.astype(complex), axes=[ax],
                   meta=meta, warnings=["digital filter removed (GRPDLY)"])


def _read_ser(ref: BrukerRef) -> NMRData:
    dic, ser = ng.bruker.read(str(ref.expno))
    acqus, acqu2s = dic["acqus"], dic.get("acqu2s", {})
    td1 = int(acqu2s.get("TD", 0)) or 1
    if ser.ndim == 1:
        # nmrglue could not reshape: do it from the direct-dim TD
        td2 = int(acqus.get("TD", 0))
        ncomplex = td2 // 2
        # round the row length up to the 256-word Bruker block boundary
        block = ((ncomplex + 127) // 128) * 128
        if block and ser.size % block == 0:
            ser = ser.reshape(-1, block)[:, :ncomplex]
        elif ncomplex and ser.size % ncomplex == 0:
            ser = ser.reshape(-1, ncomplex)
        else:
            ser = ser.reshape(td1, -1)
    ser = ng.bruker.remove_digital_filter(dic, ser)
    # drop trailing all-zero rows: a ser is often stored with one extra blank
    # acquisition beyond the real array length
    nonzero = np.abs(ser).sum(axis=1) > 0
    if nonzero.any():
        last = int(np.max(np.where(nonzero))) + 1
        ser = ser[:last]
    meta = _meta_1d(acqus, _read_title(ref.expno / "pdata" / "1"), ref.expno)
    pseudo = _pseudo2d(ref.expno, dic.get("proc2s", {}))
    warnings = ["digital filter removed (GRPDLY)"]
    delays, src = _read_delays(ref.expno)
    if pseudo and delays is not None:
        f1 = Axis("delay", "s", delays[:ser.shape[0]])
        warnings.append(f"indirect dimension = {src}")
    else:
        f1 = Axis(_isotope(dic, "acqu2s") or "F1", "point", None,
                  sw_Hz=float(acqu2s.get("SW_h", 0.0)))
    meta["fnmode"] = FNMODE.get(int(acqu2s.get("FnMODE", 0) or 0), "?")
    axes = [f1, Axis(meta["nucleus"], "s", None, obs_MHz=meta["larmor_MHz"],
                     sw_Hz=meta["sw_Hz"])]
    return NMRData(ndim=2, domain="time", data=ser.astype(complex), axes=axes,
                   meta=meta, is_pseudo2d=pseudo, warnings=warnings)


# --------------------------------------------------------------------------
def _read_title(pdata: Path) -> str:
    tp = pdata / "title"
    return tp.read_text(encoding="utf-8", errors="replace") if tp.exists() else ""


def _meta_1d(acqus: dict, title: str, expno: Path) -> dict:
    return {
        "nucleus": str(acqus.get("NUC1", "")).strip(),
        "larmor_MHz": float(acqus.get("SFO1", 0.0)),
        "sw_Hz": float(acqus.get("SW_h", 0.0)),
        "td": int(acqus.get("TD", 0)),
        "pulse_program": str(acqus.get("PULPROG", "")).strip(),
        "masr_Hz": (float(acqus["MASR"]) if acqus.get("MASR") is not None
                    else None),
        "title": title,
        "grpdly": acqus.get("GRPDLY"),
        "expno": str(expno),
    }


def _conflicts(meta: dict, title: str) -> list[str]:
    out = []
    m = re.search(r"MASR?\s*[=:]?\s*([\d.]+)\s*kHz", title, re.IGNORECASE)
    if m and meta.get("masr_Hz"):
        title_hz = float(m.group(1)) * 1000.0
        if abs(title_hz - meta["masr_Hz"]) > 0.02 * max(title_hz, meta["masr_Hz"]):
            out.append(f"MAS rate: acqus says {meta['masr_Hz']:.0f} Hz but the "
                       f"title says {title_hz:.0f} Hz -- confirm before fitting")
    return out


# --------------------------------------------------------------------------
# Backward-compatible EXPNO API (used by satrec, tests, older callers)

@dataclass
class BrukerExperiment:
    path: str
    nucleus: str
    sfo1_MHz: float
    pulse_program: str
    td: int
    sw_Hz: float
    masr_Hz: float | None
    title: str
    fid: np.ndarray | None
    processed: np.ndarray | None
    processed_ppm: np.ndarray | None
    conflicts: list = field(default_factory=list)

    @property
    def summary(self) -> str:
        lines = [
            f"EXPNO: {self.path}",
            f"nucleus: {self.nucleus}   SFO1: {self.sfo1_MHz:.4f} MHz",
            f"pulse program: {self.pulse_program}   TD: {self.td}   "
            f"SW: {self.sw_Hz:.0f} Hz",
            f"MASR (acqus): {self.masr_Hz} Hz",
            f"title: {self.title.splitlines()[0] if self.title else '(none)'}",
        ]
        for c in self.conflicts:
            lines.append(f"CONFLICT: {c}")
        return "\n".join(lines)


def read_expno(path: str | Path, procno: int = 1,
               verify: bool = True) -> BrukerExperiment:
    """Read one EXPNO folder (kept for satrec and existing callers)."""
    path = Path(path)
    before = snapshot(path) if verify else None
    dic, fid = ng.bruker.read(str(path))
    acqus = dic["acqus"]

    processed = processed_ppm = None
    pdata_dir = path / "pdata" / str(procno)
    if (pdata_dir / "procs").exists() and (pdata_dir / "1r").exists():
        dic_p, processed = ng.bruker.read_pdata(str(pdata_dir))
        processed_ppm = _ppm_axis(dic_p["procs"], processed.size)

    title = _read_title(pdata_dir)
    meta = _meta_1d(acqus, title, path)
    exp = BrukerExperiment(
        path=str(path), nucleus=meta["nucleus"], sfo1_MHz=meta["larmor_MHz"],
        pulse_program=meta["pulse_program"], td=meta["td"],
        sw_Hz=meta["sw_Hz"], masr_Hz=meta["masr_Hz"], title=title, fid=fid,
        processed=processed, processed_ppm=processed_ppm,
        conflicts=_conflicts(meta, title))
    if verify:
        verify_untouched(path, before)
    return exp
