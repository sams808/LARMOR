"""Referencing audit: is every spectrum of a session referenced the way its
1H adamantane spectrum says it should be?

Bruker stores the reference as SF (procs); SR = (SF - BF1) * 1e6 Hz. Indirect
referencing (TopSpin ``xiref``, the IUPAC unified scale) sets

    SF_X = SF_1H * Xi_X / Xi_1H

where SF_1H is the spectrometer frequency of the 1H spectrum once its
adamantane line has been put at ADAMANTANE_1H_PPM, and Xi is the IUPAC
reference-frequency ratio of the nucleus (mrsimulator's isotope table
carries every Xi through ``Isotope.B0_to_ref_freq``). A session where xiref
was forgotten keeps SF = BF1 (SR = 0) or a value carried over from another
day; both are hundreds of Hz off, far above the magnet's drift within a
session, so the check is unambiguous.

Scope rule: a session is one month folder of the NMRFAM tree
(``DATA/2026-05/<sample>/<EXPNO>``); references are never borrowed from
another month. The audit is meant to run ONCE on a new dataset: it writes a
TopSpin-ready list of ``sr`` values and appends every old and new value to
an append-only log, so any correction typed at the spectrometer can be
reversed later. Instrument folders are never written to.

Qt-free; the desktop window (larmor/desktop/referencing_dialog.py) and the
``larmor srcheck`` command are consumers.
"""
from __future__ import annotations

import csv
import datetime as _dt
import json
import os
import re
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

#: adamantane 1H line used for referencing (the group's convention; some
#: sources quote 1.85 ppm -- the difference is 0.03 ppm, below the tolerance)
ADAMANTANE_1H_PPM = 1.82
#: |stored - expected| above this (ppm of the X nucleus) is flagged
DEFAULT_TOL_PPM = 0.05
#: an SR closer to zero than this is "never referenced" (SF = BF1)
UNREFERENCED_HZ = 0.5
#: a 1H peak this close to delta_ada confirms the reference spectrum
PEAK_TOL_PPM = 0.15
_REF_TITLE_RE = re.compile(r"referenc|adamant|\bada\b|\bkbr\b", re.IGNORECASE)

__all__ = [
    "ADAMANTANE_1H_PPM", "DEFAULT_TOL_PPM", "Acquisition", "AuditRow",
    "ReferenceCheck", "xi_ratio", "expected_sf_MHz", "expected_sr_hz",
    "axis_shift_ppm", "session_root", "read_acquisition", "scan_session",
    "find_references", "pick_reference", "check_reference", "sf_from_peak",
    "audit", "consistency_notes", "summary", "topspin_list", "to_csv",
    "log_path", "append_log", "previous_audits",
]


# ---------------------------------------------------------------- physics
@lru_cache(maxsize=None)
def xi_ratio(nucleus: str) -> float:
    """Xi_X / Xi_1H: the reference frequency of ``nucleus`` divided by that
    of 1H (TMS) at the same field -- field independent."""
    from mrsimulator.spin_system.isotope import Isotope

    b0 = 9.4
    return float(Isotope(symbol=nucleus).B0_to_ref_freq(b0)
                 / Isotope(symbol="1H").B0_to_ref_freq(b0))


def expected_sf_MHz(sf_h_MHz: float, nucleus: str) -> float:
    return float(sf_h_MHz) * xi_ratio(nucleus)


def expected_sr_hz(sf_h_MHz: float, nucleus: str, bf1_MHz: float) -> float:
    return (expected_sf_MHz(sf_h_MHz, nucleus) - float(bf1_MHz)) * 1e6


def axis_shift_ppm(old_sr_hz: float, new_sr_hz: float, sf_MHz: float) -> float:
    """How the ppm axis moves when SR goes from old to new (TopSpin's
    direction): delta = (nu - SF) / SF, so a LARGER SR (larger SF) moves every
    peak to LOWER ppm by (SR_new - SR_old) / SF."""
    return -(float(new_sr_hz) - float(old_sr_hz)) / float(sf_MHz) if sf_MHz else 0.0


def sf_from_peak(sf_MHz: float, peak_ppm: float, ada_ppm: float) -> float:
    """The SF that puts a line now at ``peak_ppm`` exactly at ``ada_ppm``."""
    return float(sf_MHz) * (1.0 + float(peak_ppm) * 1e-6) / (1.0 + float(ada_ppm) * 1e-6)


# ---------------------------------------------------------------- reading
@dataclass
class Acquisition:
    path: str
    sample: str
    expno: int
    nucleus: str
    bf1_MHz: float
    sfo1_MHz: float
    o1_Hz: float
    sf_MHz: float | None            # procs SF (None: no processed data)
    sr_hz: float | None             # (SF - BF1) * 1e6
    date: float                     # acqus DATE, epoch seconds (0 if absent)
    probe: str = ""
    pulprog: str = ""
    title: str = ""
    procno: int = 1
    ns: int = 0                     # acqus NS (0 if absent)
    d1_s: float = 0.0               # acqus D[1], the recycle delay in s (0 if absent)

    @property
    def is_1h(self) -> bool:
        return self.nucleus == "1H"

    @property
    def referenced(self) -> bool:
        return self.sr_hz is not None and abs(self.sr_hz) > UNREFERENCED_HZ

    @property
    def label(self) -> str:
        return f"{self.sample}/{self.expno}"

    @property
    def date_iso(self) -> str:
        if not self.date:
            return ""
        return _dt.datetime.fromtimestamp(self.date).strftime("%Y-%m-%d %H:%M")

    @property
    def magnet_1h_MHz(self) -> float:
        """The 1H basic frequency this acquisition's BF1 implies -- the
        magnet's signature, so two spectrometers never share a reference."""
        return self.bf1_MHz / xi_ratio(self.nucleus) if self.nucleus else 0.0


def _jcamp(path: Path) -> dict | None:
    try:
        import nmrglue as ng

        return ng.fileio.bruker.read_jcamp(str(path))
    except Exception:
        return None


def _clean(v) -> str:
    return str(v).strip().strip("<>").strip() if v is not None else ""


def read_acquisition(expno_dir, procno: int = 1) -> Acquisition | None:
    """One EXPNO's referencing-relevant parameters (acqus + pdata/procno/procs
    + title); None when there is no readable acqus."""
    p = Path(expno_dir)
    acqus = _jcamp(p / "acqus")
    if not acqus:
        return None
    try:
        bf1 = float(acqus.get("BF1", 0.0) or 0.0)
        sfo1 = float(acqus.get("SFO1", 0.0) or 0.0)
        o1 = float(acqus.get("O1", 0.0) or 0.0)
        date = float(acqus.get("DATE", 0) or 0)
        ns = int(float(acqus.get("NS", 0) or 0))
        # nmrglue's read_jcamp returns the "(0..63)" delay array as a list
        d = acqus.get("D")
        d1 = float(d[1]) if isinstance(d, (list, tuple)) and len(d) > 1 else 0.0
    except (TypeError, ValueError):
        return None
    sf = sr = None
    procs = _jcamp(p / "pdata" / str(procno) / "procs")
    if procs and procs.get("SF") not in (None, 0, 0.0):
        try:
            sf = float(procs["SF"])
            sr = (sf - bf1) * 1e6 if bf1 > 0 else float(procs.get("SR", 0.0) or 0.0)
        except (TypeError, ValueError):
            sf = sr = None
    title = ""
    tp = p / "pdata" / str(procno) / "title"
    if tp.exists():
        try:
            title = tp.read_text(encoding="utf-8", errors="replace").strip()
        except OSError:
            title = ""
    try:
        expno = int(p.name)
    except ValueError:
        expno = 0
    return Acquisition(path=str(p), sample=p.parent.name, expno=expno,
                       nucleus=_clean(acqus.get("NUC1")), bf1_MHz=bf1,
                       sfo1_MHz=sfo1, o1_Hz=o1, sf_MHz=sf, sr_hz=sr,
                       date=date, probe=_clean(acqus.get("PROBHD")),
                       pulprog=_clean(acqus.get("PULPROG")), title=title,
                       procno=procno, ns=ns, d1_s=d1)


def session_root(path) -> Path:
    """The month/session folder a path belongs to: an EXPNO -> its month;
    a sample folder -> its month; anything else is taken as the session."""
    p = Path(path)
    if p.is_file():
        p = p.parent
    if (p / "acqus").exists():                  # EXPNO
        return p.parent.parent
    if p.parent.name == "pdata":                # pdata/<procno>
        return p.parent.parent.parent.parent
    if p.is_dir() and any((c / "acqus").exists() for c in p.iterdir()
                          if c.is_dir()):      # sample folder
        return p.parent
    return p


def scan_session(month_dir, progress=None) -> list[Acquisition]:
    """Every EXPNO under ``<month>/<sample>/<EXPNO>`` (and directly under a
    sample folder), oldest first."""
    root = Path(month_dir)
    if not root.is_dir():
        return []
    out = []
    samples = sorted(c for c in root.iterdir() if c.is_dir())
    for k, sample in enumerate(samples):
        if progress:
            progress(k, len(samples), sample.name)
        for e in sorted(sample.iterdir(), key=lambda c: (not c.name.isdigit(),
                                                         int(c.name) if c.name.isdigit() else 0)):
            if e.is_dir() and (e / "acqus").exists():
                a = read_acquisition(e)
                if a is not None:
                    out.append(a)
    # an EXPNO directly under the given folder (a sample folder was passed)
    for e in samples:
        if (e / "acqus").exists():
            a = read_acquisition(e)
            if a is not None and a.path not in {o.path for o in out}:
                out.append(a)
    out.sort(key=lambda a: (a.date, a.sample, a.expno))
    return out


# ---------------------------------------------------------------- references
def find_references(acqs: list[Acquisition]) -> list[Acquisition]:
    """The 1H acquisitions that carry a reference (SR set), oldest first;
    ones whose title says so are listed first within a date."""
    refs = [a for a in acqs if a.is_1h and a.referenced]
    refs.sort(key=lambda a: (a.date, 0 if _REF_TITLE_RE.search(a.title or "") else 1))
    return refs


def pick_reference(refs: list[Acquisition], acq: Acquisition) -> Acquisition | None:
    """The reference nearest in time on the SAME magnet (BF1 implies the 1H
    basic frequency; 599.77 MHz and 1100.35 MHz never mix)."""
    same = [r for r in refs
            if abs(r.bf1_MHz - acq.magnet_1h_MHz) < 0.002 * r.bf1_MHz]
    if not same:
        return None
    return min(same, key=lambda r: abs(r.date - acq.date))


@dataclass
class ReferenceCheck:
    acq: Acquisition
    peak_ppm: float | None
    status: str
    ok: bool
    ppm: np.ndarray | None = None
    amp: np.ndarray | None = None


def check_reference(acq: Acquisition, ada_ppm: float = ADAMANTANE_1H_PPM,
                    load=None) -> ReferenceCheck:
    """Read the 1H reference spectrum and confirm its tallest line sits at
    delta_ada: the referencing itself can be wrong (adamantane put at 0 ppm,
    or the wrong line picked)."""
    if load is None:
        from larmor.loader import load_any

        def load(path):
            ppm, amp, _r, _m, _w = load_any(path)
            return np.asarray(ppm, float), np.asarray(amp, float)
    try:
        ppm, amp = load(acq.path)
    except Exception as exc:
        return ReferenceCheck(acq, None, f"spectrum not readable: {exc}", False)
    if ppm.size == 0:
        return ReferenceCheck(acq, None, "empty spectrum", False)
    peak = float(ppm[int(np.argmax(np.abs(amp)))])
    if abs(peak - ada_ppm) <= PEAK_TOL_PPM:
        return ReferenceCheck(acq, peak, f"adamantane line at {peak:.2f} ppm "
                              f"(expected {ada_ppm:g}) — reference confirmed",
                              True, ppm, amp)
    if abs(peak) <= PEAK_TOL_PPM:
        return ReferenceCheck(acq, peak, f"tallest line at {peak:.2f} ppm: the 1H "
                              "seems referenced to 0 ppm, not to adamantane "
                              f"at {ada_ppm:g}", False, ppm, amp)
    return ReferenceCheck(acq, peak, f"tallest line at {peak:.2f} ppm — not "
                          f"adamantane at {ada_ppm:g}; check the reference",
                          False, ppm, amp)


# ---------------------------------------------------------------- the audit
@dataclass
class AuditRow:
    acq: Acquisition
    ref: Acquisition | None
    expected_sf_MHz: float | None
    expected_sr_hz: float | None
    delta_ppm: float | None          # (stored SF - expected SF) in ppm of X
    delta_hz: float | None           # stored SR - expected SR
    verdict: str                     # ok | off | unreferenced | no reference | unprocessed | 1H reference | 1H
    note: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def flagged(self) -> bool:
        return self.verdict in ("off", "unreferenced")


def audit(acqs: list[Acquisition], *, tol_ppm: float = DEFAULT_TOL_PPM,
          ada_ppm: float = ADAMANTANE_1H_PPM,
          ref_sf: dict[str, float] | None = None) -> list[AuditRow]:
    """Judge every acquisition of a session against the nearest 1H reference
    on its magnet. ``ref_sf`` optionally overrides the SF of a reference
    (by path) -- e.g. one re-derived from its adamantane peak."""
    refs = find_references(acqs)
    rows = []
    for a in acqs:
        if a.is_1h:
            rows.append(AuditRow(a, None, None, None, None, None,
                                 "1H reference" if a.referenced else "1H",
                                 "" if a.referenced else "SR = 0: not a usable reference"))
            continue
        ref = pick_reference(refs, a)
        if ref is None:
            rows.append(AuditRow(a, None, None, None, None, None, "no reference",
                                 "no referenced 1H spectrum on this magnet in the session"))
            continue
        sf_h = (ref_sf or {}).get(ref.path, ref.sf_MHz)
        exp_sf = expected_sf_MHz(sf_h, a.nucleus)
        exp_sr = (exp_sf - a.bf1_MHz) * 1e6
        if a.sf_MHz is None:
            rows.append(AuditRow(a, ref, exp_sf, exp_sr, None, None, "unprocessed",
                                 "no processed data (procs) yet"))
            continue
        d_ppm = (a.sf_MHz - exp_sf) / exp_sf * 1e6
        d_hz = float(a.sr_hz) - exp_sr
        if not a.referenced:
            verdict, note = "unreferenced", "SR = 0 (SF = BF1): xiref was never applied"
        elif abs(d_ppm) <= tol_ppm:
            verdict, note = "ok", ""
        else:
            verdict, note = "off", f"stored SR is {d_ppm:+.3f} ppm from the 1H reference"
        rows.append(AuditRow(a, ref, exp_sf, exp_sr, d_ppm, d_hz, verdict, note))
    consistency_notes(rows, tol_ppm)
    return rows


def consistency_notes(rows: list[AuditRow], tol_ppm: float = DEFAULT_TOL_PPM) -> None:
    """Flag a nucleus whose stored SF differs between EXPNOs of the session
    (a partly referenced day) -- independent of any 1H reference."""
    by = {}
    for r in rows:
        if not r.acq.is_1h and r.acq.sf_MHz:
            by.setdefault(r.acq.nucleus, []).append(r)
    for nuc, group in by.items():
        sfs = np.array([r.acq.sf_MHz for r in group])
        if sfs.size > 1 and (sfs.max() - sfs.min()) / sfs.mean() * 1e6 > tol_ppm:
            distinct = len({round(s, 6) for s in sfs})
            for r in group:
                r.notes.append(f"{nuc}: {distinct} different SF values within the session")


def summary(rows: list[AuditRow]) -> dict:
    counts = {}
    for r in rows:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    return counts


# ---------------------------------------------------------------- outputs
def topspin_list(rows: list[AuditRow], only_flagged: bool = True,
                 session: str = "") -> str:
    """A paste-ready list: for each EXPNO the ``sr`` value to type in TopSpin,
    with the value it replaces so the change is reversible by hand too."""
    lines = [f"# LARMOR referencing audit — {session or 'session'} — "
             f"{_dt.datetime.now().strftime('%Y-%m-%d %H:%M')}",
             "# sample/EXPNO  procno  nucleus   sr <new Hz>      (was <old Hz>; offset ppm)"]
    for r in rows:
        if r.acq.is_1h or r.expected_sr_hz is None:
            continue
        if only_flagged and not r.flagged:
            continue
        old = f"{r.acq.sr_hz:.2f}" if r.acq.sr_hz is not None else "none"
        off = f"{r.delta_ppm:+.3f} ppm" if r.delta_ppm is not None else "unprocessed"
        lines.append(f"{r.acq.label:<28} {r.acq.procno:<3} {r.acq.nucleus:<6} "
                     f"sr {r.expected_sr_hz:10.2f}   (was {old}; {off})")
    if len(lines) == 2:
        lines.append("# nothing to correct")
    return "\n".join(lines) + "\n"


def to_csv(rows: list[AuditRow], path) -> Path:
    p = Path(path)
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["path", "sample", "expno", "nucleus", "date", "stored_sf_MHz",
                    "stored_sr_hz", "reference", "expected_sf_MHz",
                    "expected_sr_hz", "delta_ppm", "delta_hz", "verdict", "note"])
        for r in rows:
            w.writerow([r.acq.path, r.acq.sample, r.acq.expno, r.acq.nucleus,
                        r.acq.date_iso,
                        f"{r.acq.sf_MHz:.6f}" if r.acq.sf_MHz else "",
                        f"{r.acq.sr_hz:.2f}" if r.acq.sr_hz is not None else "",
                        r.ref.label if r.ref else "",
                        f"{r.expected_sf_MHz:.6f}" if r.expected_sf_MHz else "",
                        f"{r.expected_sr_hz:.2f}" if r.expected_sr_hz is not None else "",
                        f"{r.delta_ppm:+.4f}" if r.delta_ppm is not None else "",
                        f"{r.delta_hz:+.2f}" if r.delta_hz is not None else "",
                        r.verdict, "; ".join([r.note] + r.notes).strip("; ")])
    return p


def log_path() -> Path:
    """Append-only record of every audit: %LOCALAPPDATA%/LARMOR/referencing_log.jsonl
    (LARMOR_REF_LOG overrides). Never truncated by the application."""
    env = os.environ.get("LARMOR_REF_LOG")
    if env:
        return Path(env)
    base = os.environ.get("LOCALAPPDATA") or str(Path.home())
    return Path(base) / "LARMOR" / "referencing_log.jsonl"


def append_log(rows: list[AuditRow], session, *, ada_ppm: float, tol_ppm: float,
               action: str, path=None) -> Path:
    """Append one record holding the OLD and NEW SR/SF of every flagged EXPNO
    (and the reference used), so a correction can be reversed years later."""
    from larmor import __version__

    p = Path(path) if path else log_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    entries = []
    for r in rows:
        if r.acq.is_1h or r.expected_sr_hz is None or not r.flagged:
            continue
        entries.append({
            "path": r.acq.path, "sample": r.acq.sample, "expno": r.acq.expno,
            "procno": r.acq.procno, "nucleus": r.acq.nucleus,
            "old_sf_MHz": r.acq.sf_MHz, "old_sr_hz": r.acq.sr_hz,
            "new_sf_MHz": r.expected_sf_MHz, "new_sr_hz": r.expected_sr_hz,
            "delta_ppm": r.delta_ppm, "verdict": r.verdict,
            "reference": r.ref.path if r.ref else None,
            "reference_sf_MHz": r.ref.sf_MHz if r.ref else None,
        })
    record = {"time": _dt.datetime.now().isoformat(timespec="seconds"),
              "larmor": __version__, "session": str(Path(session)),
              "action": action, "ada_ppm": ada_ppm, "tol_ppm": tol_ppm,
              "counts": summary(rows), "entries": entries}
    with open(p, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
    return p


def previous_audits(session, path=None) -> list[dict]:
    """Earlier log records for this session (the audit is meant to run once)."""
    p = Path(path) if path else log_path()
    if not p.exists():
        return []
    want = str(Path(session)).lower()
    out = []
    for line in p.read_text(encoding="utf-8").splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if str(rec.get("session", "")).lower() == want:
            out.append(rec)
    return out


def row_dict(r: AuditRow) -> dict:
    """JSON-friendly view of a row (for the CLI and tests)."""
    d = asdict(r.acq)
    d.update(reference=r.ref.label if r.ref else None,
             expected_sf_MHz=r.expected_sf_MHz, expected_sr_hz=r.expected_sr_hz,
             delta_ppm=r.delta_ppm, delta_hz=r.delta_hz, verdict=r.verdict,
             note=r.note, notes=list(r.notes))
    return d
