"""Comparability of a spectrum series: were its members acquired and
processed the same way?

A batch fit shares one lineshape model across a series, so every member is
implicitly assumed to have been *measured and processed alike*. TopSpin never
says otherwise: one glass of a composition series apodised with LB = 100 Hz
next to one at 0 Hz, a fid cut at TDeff 384 of 9590 points, a repeat run at
D1 = 90 s instead of 300 s -- all look like ordinary 1r files. This module
reads what the instrument wrote (``acqus``, ``pdata/N/procs``, the operator's
``title`` and TopSpin's own command trail ``auditp.txt``), compares every
member of a series against the series *majority*, and maps the TopSpin
processing parameters onto LARMOR's own processing ops so a series can be
rebuilt from its fids with one common pipeline.

Vocabulary (aligned with :mod:`larmor.fithealth`):

- ``bad``   -- a shared model is meaningless (mixed nuclei, different fields);
- ``check`` -- a processing or acquisition parameter differs from the
  majority and biases a shared-width model (window, LB, TDeff, SI, phase
  mode, baseline mode, referencing, pulse program, pulse, power, D1, ...);
- ``info``  -- reported, never flagged (NS, RG, date, the phase values, the
  audit command line, the title): legitimately per spectrum.

Majority rule: the most frequent value (floats within a tolerance); an exact
tie takes the first spectrum's value and says so; a key where every member
differs (three or more) reports once as *varies min-max* with no
per-spectrum flag.

TopSpin conventions pinned on real data (tests/test_comparability.py):
``TDeff`` counts real points while the fid is complex, so the LARMOR
``tdeff`` step takes ``TDeff // 2`` points; the stored phase replays as
``phase(p0=-PHC0, p1=PHC1, pivot_frac=1.0)`` on ``op_ft``'s ascending axis.
``abs`` (ABSG) and ``BC_mod`` are reported but never replayed -- the batch
dialog's *Fit baseline...* applies one identical method to every spectrum.

Read-only by construction: the only file access is nmrglue's JCAMP reader
(through :func:`larmor.referencing._jcamp`) and ``Path.read_text``. Qt-free;
the desktop widgets live in :mod:`larmor.desktop.comparability_dialog` and
the command line is ``larmor compare``.
"""
from __future__ import annotations

import csv
import datetime as _dt
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from larmor.io import bruker
from larmor.qcpmg import carrier_ppm
from larmor.referencing import DEFAULT_TOL_PPM, _clean, _jcamp

__all__ = [
    "TOPSPIN_WDW", "TOPSPIN_PH_MOD", "LEVELS", "PROC_SHAPE_KEYS", "KeyRule",
    "KEY_RULES", "SpectrumParams", "read_params", "parse_auditp", "KeyReport",
    "Comparison", "compare", "majority_template", "template_from", "ops_for",
    "describe_template", "reprocess", "as_reprocessed", "caveat_note",
]

#: TopSpin WDW codes -> window names
TOPSPIN_WDW = {0: "none", 1: "EM", 2: "GM", 3: "SINE", 4: "QSINE", 5: "TRAP",
               6: "USER", 7: "SINC", 8: "QSINC", 9: "TRAF", 10: "TRAFS"}
#: TopSpin PH_mod codes -> phase-correction modes
TOPSPIN_PH_MOD = {0: "no", 1: "pk", 2: "mc", 3: "ps"}
#: severity vocabulary (same words as larmor.fithealth.LEVELS)
LEVELS = ("bad", "check", "info")
#: the processing keys a reprocess-from-fid makes identical by construction
PROC_SHAPE_KEYS = ("WDW", "LB", "GB", "SSB", "TDeff", "SI", "FCOR")
#: the Larmor-frequency spread above which a series mixes fields
#: (batch.homogeneity's rule)
FIELD_TOL_REL = 0.05
#: recycle delays within this relative spread count as the same
D1_TOL_REL = 0.05


# ------------------------------------------------------------------ formatting
def _g(v) -> str:
    return format(float(v), "g")


def _si_fmt(v, _p=None) -> str:
    v = int(v)
    return f"{v // 1024}k" if v and v % 1024 == 0 else str(v)


def _tdeff_fmt(v, p=None) -> str:
    v = int(v)
    td = int((p.acq.get("TD") if p is not None else 0) or 0)
    return "all" if v <= 0 or (td and v >= td) else str(v)


def _plw_fmt(v, p=None) -> str:
    unit = (p.acq.get("PLW1_unit") if p is not None else "") or ""
    return (_g(v) + (f" {unit}" if unit else "")).strip()


def _date_fmt(v, _p=None) -> str:
    try:
        return _dt.datetime.fromtimestamp(float(v)).strftime("%Y-%m-%d %H:%M")
    except (OverflowError, OSError, ValueError):
        return str(v)


def _title_fmt(v, _p=None) -> str:
    return str(v).splitlines()[0] if str(v).strip() else ""


def _cmds_fmt(v, _p=None) -> str:
    return " · ".join(v) if isinstance(v, (list, tuple)) else str(v)


# ------------------------------------------------------------------ rules
@dataclass(frozen=True)
class KeyRule:
    """How one parameter is read, displayed and judged."""

    key: str
    label: str                      # table / CSV header (carries the unit)
    group: str                      # "processing" | "referencing" | "acquisition"
    level: str                      # "bad" | "check" | "info"
    source: str                     # "proc" | "acq" | "derived"
    unit: str = ""                  # prose unit ("Hz", "s"); "" when fmt carries it
    tol: float | None = None        # relative tolerance for floats (None: default)
    fmt: object = None              # callable(value, params) -> str
    when: object = None             # callable(params) -> bool: row applies to it
    mention: bool = False           # info key whose differences the summary names

    def display(self, v, params=None) -> str:
        if v is None:
            return "—"
        if self.fmt is not None:
            return self.fmt(v, params)
        if isinstance(v, float):
            return _g(v)
        return str(v)


def _wdw_is(*codes):
    return lambda p: int(p.proc.get("WDW") or 0) in codes


KEY_RULES: tuple[KeyRule, ...] = (
    # ---- processing (pdata/N/procs, the audit trail)
    KeyRule("WDW", "window", "processing", "check", "proc",
            fmt=lambda v, p: TOPSPIN_WDW.get(int(v), str(v))),
    KeyRule("LB", "LB (Hz)", "processing", "check", "proc", unit="Hz"),
    KeyRule("GB", "GB", "processing", "check", "proc", when=_wdw_is(2)),
    KeyRule("SSB", "SSB", "processing", "check", "proc", when=_wdw_is(3, 4)),
    KeyRule("TDeff", "TDeff (points)", "processing", "check", "proc", fmt=_tdeff_fmt),
    KeyRule("SI", "SI", "processing", "check", "proc", fmt=_si_fmt),
    KeyRule("FCOR", "FCOR", "processing", "check", "proc"),
    KeyRule("PH_mod", "phase mode", "processing", "check", "proc",
            fmt=lambda v, p: TOPSPIN_PH_MOD.get(int(v), str(v))),
    KeyRule("PHC0", "PHC0 (°)", "processing", "info", "proc",
            fmt=lambda v, p: f"{float(v):.2f}"),
    KeyRule("PHC1", "PHC1 (°)", "processing", "info", "proc",
            fmt=lambda v, p: f"{float(v):.2f}"),
    KeyRule("ABSG", "ABSG (abs order)", "processing", "check", "proc"),
    KeyRule("BC_mod", "BC_mod", "processing", "check", "proc"),
    KeyRule("FT_mod", "FT_mod", "processing", "info", "proc"),
    KeyRule("commands", "TopSpin commands (auditp)", "processing", "info", "derived",
            fmt=_cmds_fmt),
    # ---- referencing
    KeyRule("SF", "SF (MHz)", "referencing", "check", "proc", unit="MHz",
            tol=DEFAULT_TOL_PPM * 1e-6, fmt=lambda v, p: f"{float(v):.6f}"),
    # ---- acquisition (acqus)
    KeyRule("NUC1", "nucleus", "acquisition", "bad", "acq"),
    KeyRule("SFO1", "SFO1 (MHz)", "acquisition", "bad", "acq", unit="MHz",
            tol=FIELD_TOL_REL, fmt=lambda v, p: f"{float(v):.4f}"),
    KeyRule("PULPROG", "pulse program", "acquisition", "check", "acq"),
    KeyRule("SW_h", "SW (Hz)", "acquisition", "check", "acq", unit="Hz"),
    KeyRule("TD", "TD (points)", "acquisition", "check", "acq"),
    KeyRule("D1", "D1 (s)", "acquisition", "check", "acq", unit="s", tol=D1_TOL_REL),
    KeyRule("P1", "P1 (µs)", "acquisition", "check", "acq", unit="µs"),
    KeyRule("PLW1", "PLW1", "acquisition", "check", "acq", fmt=_plw_fmt),
    KeyRule("NS", "NS", "acquisition", "info", "acq", mention=True),
    KeyRule("RG", "RG", "acquisition", "info", "acq", mention=True),
    KeyRule("PROBHD", "probe", "acquisition", "check", "acq"),
    KeyRule("DATE", "date", "acquisition", "info", "acq", fmt=_date_fmt),
    KeyRule("title", "title", "acquisition", "info", "derived", fmt=_title_fmt),
)
_RULES = {r.key: r for r in KEY_RULES}


# ------------------------------------------------------------------ reading
@dataclass
class SpectrumParams:
    """One Bruker spectrum's acquisition / processing parameters, flattened."""

    path: str                       # the path given (1r, pdata/N, EXPNO, fid)
    expno: str                      # the EXPNO folder
    procno: int
    sample: str                     # the sample folder's name
    acq: dict = field(default_factory=dict)
    proc: dict = field(default_factory=dict)
    commands: list = field(default_factory=list)
    title: str = ""
    has_fid: bool = False

    def value(self, key: str):
        rule = _RULES.get(key)
        if rule is None:
            return None
        if rule.when is not None and not rule.when(self):
            return None
        if key == "commands":
            return list(self.commands) if self.commands else None
        if key == "title":
            return self.title or None
        src = self.proc if rule.source == "proc" else self.acq
        return src.get(key)


def _f(d: dict, key: str):
    v = d.get(key)
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _i(d: dict, key: str):
    v = _f(d, key)
    return None if v is None else int(v)


def _idx(v, i: int):
    lst = bruker._num_list(v)
    return lst[i] if len(lst) > i else None


def read_params(path, procno: int | None = None) -> SpectrumParams | None:
    """Read ``path``'s acqus + pdata/<procno>/procs + title + auditp.txt.

    ``path`` is anything :func:`larmor.io.bruker.resolve` accepts (a 1r file,
    a pdata/N folder, an EXPNO, a fid); ``procno`` defaults to the path's
    own pdata/N (else 1). None for a CSV / fxmla / recipe / anything without
    a readable acqus. Nothing is written."""
    p = Path(str(path))
    try:
        ref = bruker.resolve(p)
        expno, ref_procno = ref.expno, int(ref.procno or 1)
    except (ValueError, FileNotFoundError, OSError, NotADirectoryError):
        # an EXPNO holding acqus/procs but no 1r/fid yet (resolve wants data)
        if not bruker.is_expno(p):
            return None
        expno, ref_procno = p, 1
    acqus = _jcamp(expno / "acqus")
    if not acqus:
        return None
    if procno is None:
        procno = ref_procno
    pdata = expno / "pdata" / str(procno)
    procs = _jcamp(pdata / "procs") or {}

    acq = {
        "NUC1": _clean(acqus.get("NUC1")),
        "PULPROG": _clean(acqus.get("PULPROG")),
        "NS": _i(acqus, "NS"), "TD": _i(acqus, "TD"), "RG": _f(acqus, "RG"),
        "SW_h": _f(acqus, "SW_h"), "SFO1": _f(acqus, "SFO1"),
        "BF1": _f(acqus, "BF1"), "O1": _f(acqus, "O1"), "DATE": _f(acqus, "DATE"),
        "PROBHD": _clean(acqus.get("PROBHD")),
        "D1": _idx(acqus.get("D"), 1), "P1": _idx(acqus.get("P"), 1),
    }
    plw1 = _idx(acqus.get("PLW"), 1)
    if plw1 is not None:
        acq["PLW1"], acq["PLW1_unit"] = plw1, "W"
    else:
        pl1 = _idx(acqus.get("PL"), 1)
        acq["PLW1"], acq["PLW1_unit"] = pl1, ("dB" if pl1 is not None else "")

    proc = {}
    if procs:
        proc = {
            "WDW": _i(procs, "WDW"), "LB": _f(procs, "LB"), "GB": _f(procs, "GB"),
            "SSB": _f(procs, "SSB"), "TDeff": _i(procs, "TDeff"), "SI": _i(procs, "SI"),
            "FCOR": _f(procs, "FCOR"), "PH_mod": _i(procs, "PH_mod"),
            "PHC0": _f(procs, "PHC0"), "PHC1": _f(procs, "PHC1"),
            "ABSG": _i(procs, "ABSG"), "BC_mod": _i(procs, "BC_mod"),
            "FT_mod": _i(procs, "FT_mod"), "SF": _f(procs, "SF"),
            "OFFSET": _f(procs, "OFFSET"), "SW_p": _f(procs, "SW_p"),
        }
        pknl = procs.get("PKNL")
        if pknl is not None:
            proc["PKNL"] = bool(pknl) if not isinstance(pknl, str) else \
                pknl.strip().lower() in ("yes", "true", "1")
        sf = proc.get("SF")
        if sf is not None and sf <= 0:      # an unprocessed procs writes SF = 0
            proc["SF"] = None
    title = _clean_title(bruker._read_title(pdata))
    return SpectrumParams(
        path=str(path), expno=str(expno), procno=procno, sample=expno.parent.name,
        acq=acq, proc=proc, commands=parse_auditp(pdata), title=title,
        has_fid=(expno / "fid").is_file())


def _clean_title(t: str) -> str:
    return (t or "").replace("\r", "").strip()


#: one audit-trail record: (NUMBER,<WHEN>,<WHO>,<WHERE>,<PROCESS>,<VERSION>,<WHAT>)
_AUDIT_REC = re.compile(
    r"\(\s*\d+\s*,<[^>]*>,<[^>]*>,<[^>]*>,<([^>]*)>,<[^>]*>,\s*<(.*?)>\s*\)",
    re.DOTALL)
_HEX_LINE = re.compile(r"^(?:[0-9A-Fa-f]{2}\s*)+$")


def parse_auditp(path) -> list[str]:
    """The TopSpin processing commands behind the current spectrum, from
    ``auditp.txt`` (a file or its pdata folder): the command line of every
    ``proc1d`` record since the last *Start of raw data processing* -- e.g.
    ``['em LB = 0 TDeff = 384 SI = 32K', 'ft FT_mod = 6 PKNL = 1 SI = 32K',
    'apk', 'abs n ABSG = 5']``. Acquisition records (``go4``, *created by
    zg*), parameter edits (``cprserver``) and user comments (``audit``) are
    excluded. ``[]`` when the file is absent or unparsable."""
    p = Path(str(path))
    if p.is_dir():
        p = p / "auditp.txt"
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except (OSError, ValueError):
        return []
    cmds: list[str] = []
    try:
        for m in _AUDIT_REC.finditer(text):
            process, what = m.group(1).strip(), m.group(2)
            if process != "proc1d":
                continue
            lines = [ln.strip() for ln in what.splitlines() if ln.strip()]
            if lines and lines[0].lower().startswith("start of raw data processing"):
                cmds = []            # a fresh chain from the fid supersedes the old
                lines = lines[1:]
            for ln in lines:
                low = ln.lower()
                if ("data hash" in low or "hash md5" in low or _HEX_LINE.match(ln)
                        or low.startswith("configuration hash")):
                    continue
                cmds.append(ln)
                break
    except Exception:                    # noqa: BLE001 -- a truncated trail
        return cmds
    return cmds


# ------------------------------------------------------------------ comparing
@dataclass
class KeyReport:
    """One parameter across the series."""

    key: str
    label: str
    group: str
    level: str
    values: list                    # raw per spectrum (None: not available)
    display: list                   # str per spectrum ("—" when None)
    majority: object                # the reference value (None when nothing)
    majority_display: str
    deviants: list                  # indices that differ from the majority
    varies: bool = False            # every member differs (>= 3): no per-spectrum flag
    tie: bool = False               # the majority was a tie: first spectrum taken
    message: str = ""               # one prose line ("LB 0 / 100 Hz (EM)")
    mention: bool = False           # an info key whose differences are named

    @property
    def differs(self) -> bool:
        return bool(self.deviants) or self.varies

    @property
    def flagged(self) -> bool:
        return self.differs and self.level != "info"


def _same(a, b, tol: float | None) -> bool:
    if isinstance(a, (int, float, np.integer, np.floating)) and \
            isinstance(b, (int, float, np.integer, np.floating)) and \
            not isinstance(a, bool) and not isinstance(b, bool):
        fa, fb = float(a), float(b)
        if fa == fb:
            return True
        t = 1e-6 if tol is None else float(tol)
        return abs(fa - fb) <= t * max(abs(fa), abs(fb))
    return a == b


def _clusters(values: list, tol: float | None) -> list[list[int]]:
    """Group indices of equal (within tol) non-None values, in first-seen order."""
    out: list[list[int]] = []
    reps: list = []
    for i, v in enumerate(values):
        if v is None:
            continue
        for c, r in zip(out, reps):
            if _same(v, r, tol):
                c.append(i)
                break
        else:
            out.append([i])
            reps.append(v)
    return out


def _sorted_distinct(reps: list) -> list:
    try:
        return sorted(reps)
    except TypeError:
        return list(reps)


def _report(rule: KeyRule, params: list) -> KeyReport | None:
    values = [p.value(rule.key) if p is not None else None for p in params]
    if all(v is None for v in values):
        return None
    display = [rule.display(v, p) for v, p in zip(values, params)]
    clusters = _clusters(values, rule.tol)
    n_present = sum(1 for v in values if v is not None)
    best = max(len(c) for c in clusters)
    top = [c for c in clusters if len(c) == best]
    maj_cluster = min(top, key=lambda c: c[0])
    k0 = maj_cluster[0]
    majority, maj_disp = values[k0], display[k0]
    tie = len(top) > 1
    varies = best == 1 and n_present >= 3
    deviants = [] if (varies or len(clusters) == 1) else \
        sorted(i for c in clusters if c is not maj_cluster for i in c)
    rep = KeyReport(rule.key, rule.label, rule.group, rule.level, values, display,
                    majority, maj_disp, deviants, varies=varies, tie=tie,
                    mention=rule.mention)
    if rep.differs:
        rep.message = _message(rule, rep, params)
    return rep


def _message(rule: KeyRule, rep: KeyReport, params: list) -> str:
    reps = [rep.values[c[0]] for c in _clusters(rep.values, rule.tol)]
    p0 = next((p for p in params if p is not None), None)
    disp = [rule.display(v, p0) for v in _sorted_distinct(reps)]
    unit = f" {rule.unit}" if rule.unit else ""
    if rep.varies:
        txt = f"{rule.key} varies {disp[0]}–{disp[-1]}{unit}"
    elif len(disp) > 4:
        txt = f"{rule.key} {disp[0]}–{disp[-1]}{unit} ({len(disp)} values)"
    else:
        txt = f"{rule.key} {' / '.join(disp)}{unit}"
    if rule.key == "LB":
        wins = sorted({TOPSPIN_WDW.get(int(p.proc.get("WDW") or 0), "?")
                       for p in params if p is not None and p.proc.get("WDW") is not None})
        if wins:
            txt += f" ({'/'.join(wins)})"
    if rule.key == "SF" and len(reps) >= 2:
        lo, hi = min(reps), max(reps)
        txt += (f" — {(hi - lo) * 1e6 / hi:.2f} ppm apart; see Process ▸ "
                "Reference ▸ Referencing audit")
    if rep.tie and not rep.varies:
        txt += " (tie: first spectrum taken as reference)"
    return txt


def _bad_messages(reports: dict) -> None:
    """Give the NUC1 / SFO1 rows batch.homogeneity's own wording (the batch
    report and this check must never disagree) and demote an SFO1 row whose
    spread is within the field tolerance to info."""
    nuc, sfo = reports.get("NUC1"), reports.get("SFO1")
    if sfo is not None and not sfo.differs:
        sfo.level = "info"
    if not ((nuc is not None and nuc.differs) or (sfo is not None and sfo.differs)):
        return
    from larmor import batch

    entries = []
    n = len((nuc or sfo).values)
    for i in range(n):
        entries.append(SimpleNamespace(
            nucleus=(nuc.values[i] if nuc is not None else "") or "",
            larmor_MHz=float((sfo.values[i] if sfo is not None else 0.0) or 0.0),
            models=()))
    for note in batch.homogeneity(entries):
        if "mixed nuclei" in note and nuc is not None and nuc.differs:
            nuc.message = note
        elif "Larmor frequencies" in note and sfo is not None and sfo.differs:
            sfo.message = note


@dataclass
class Comparison:
    """A series compared against its majority; the object the bar, the
    Details dialog, the CLI and the recipes' caveat notes all read."""

    params: list                    # SpectrumParams | None per spectrum
    labels: list
    reports: list = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.params)

    @property
    def n_bruker(self) -> int:
        return sum(1 for p in self.params if p is not None)

    def flagged(self) -> list:
        return [r for r in self.reports if r.flagged]

    def noted(self) -> list:
        """Info rows the summary still names when they differ (NS, RG)."""
        return [r for r in self.reports if r.mention and r.differs
                and r.level == "info"]

    @property
    def any_flagged(self) -> bool:
        return bool(self.flagged())

    @property
    def level(self) -> str:
        if self.n_bruker < 2 or not self.reports:
            return "none"
        fl = self.flagged()
        if any(r.level == "bad" for r in fl):
            return "bad"
        return "check" if fl else "ok"

    def report(self, key: str) -> KeyReport | None:
        return next((r for r in self.reports if r.key == key), None)

    # ---- prose
    def summary(self) -> str:
        lvl = self.level
        if lvl == "none":
            return ""
        if lvl == "ok":
            bits = []
            for key in ("WDW", "TDeff", "SI", "NS", "D1"):
                r = self.report(key)
                if r is None or r.majority is None:
                    continue
                rule = _RULES[key]
                unit = f" {rule.unit}" if rule.unit else ""
                if key == "WDW":
                    lb = self.report("LB")
                    txt = r.majority_display if r.majority else "no window"
                    if lb is not None and lb.majority is not None and r.majority:
                        txt += f" {lb.majority_display} Hz"
                elif r.differs:
                    reps = _sorted_distinct(
                        [r.values[c[0]] for c in _clusters(r.values, rule.tol)])
                    p0 = next((p for p in self.params if p is not None), None)
                    txt = (f"{key} {rule.display(reps[0], p0)}–"
                           f"{rule.display(reps[-1], p0)}{unit}")
                else:
                    txt = f"{key} {r.majority_display}{unit}"
                bits.append(txt)
            return (f"✓ {self.n_bruker} spectra acquired and processed alike"
                    + (" — " + " · ".join(bits) if bits else ""))
        parts = []
        heads = {"processing": "processed differently",
                 "referencing": "referenced differently",
                 "acquisition": "acquired differently"}
        for group in ("processing", "referencing", "acquisition"):
            msgs = [r.message for r in self.flagged() if r.group == group and r.message]
            if group == "acquisition":
                msgs += [r.message for r in self.noted() if r.message]
            if msgs:
                parts.append(f"{heads[group]}: " + " · ".join(msgs))
        tail = (" — a shared model cannot describe these spectra together"
                if lvl == "bad" else
                " — populations from one shared model are not strictly comparable")
        return "⚠ " + " · ".join(parts) + tail

    def _deviant_reports(self, k: int) -> list:
        return [r for r in self.flagged() if k in r.deviants]

    def level_of(self, k: int) -> str:
        """'bad' / 'check' / '' for spectrum k (the worst flagged row it deviates on)."""
        lv = [r.level for r in self._deviant_reports(k)]
        return "bad" if "bad" in lv else ("check" if lv else "")

    def signature(self, k: int) -> str:
        """'LB 100 Hz; TDeff 768' for spectrum k, '' when it is not a deviant."""
        out = []
        for r in self._deviant_reports(k):
            rule = _RULES[r.key]
            unit = f" {rule.unit}" if rule.unit else ""
            out.append(f"{r.key} {r.display[k]}{unit}")
        return "; ".join(out)

    def details(self, k: int) -> list:
        """Tooltip lines: 'LB 100 Hz — series majority 0 Hz'."""
        out = []
        for r in self._deviant_reports(k):
            rule = _RULES[r.key]
            unit = f" {rule.unit}" if rule.unit else ""
            out.append(f"{r.key} {r.display[k]}{unit} — series majority "
                       f"{r.majority_display}{unit}")
        return out

    def without(self, keys) -> "Comparison":
        keys = set(keys)
        return replace(self, reports=[r for r in self.reports if r.key not in keys])

    # ---- tables
    def rows(self, only_differences: bool = False) -> list:
        head = ["group", "key", "label", "level", "majority", *self.labels]
        body = []
        for r in self.reports:
            if only_differences and not (r.flagged or (r.mention and r.differs)):
                continue
            body.append([r.group, r.key, r.label, r.level, r.majority_display,
                         *r.display])
        return [head, *body]

    def to_text(self, only_differences: bool = True) -> str:
        rows = self.rows(only_differences)
        if len(rows) == 1:
            return ("no Bruker parameters to compare" if self.level == "none"
                    else "no differences from the series majority")
        marks = []
        for r in self.reports:
            if only_differences and not (r.flagged or (r.mention and r.differs)):
                continue
            marks.append(r)
        out = [rows[0]]
        for rep, row in zip(marks, rows[1:]):
            cells = list(row[:5])
            for k, d in enumerate(row[5:]):
                cells.append(f"{d} *" if k in rep.deviants else d)
            out.append(cells)
        widths = [max(len(str(r[c])) for r in out) for c in range(len(out[0]))]
        lines = ["  ".join(str(c).ljust(w) for c, w in zip(r, widths)).rstrip()
                 for r in out]
        if any(r.deviants for r in marks):
            lines.append("* differs from the series majority")
        return "\n".join(lines)

    def to_csv(self, path) -> None:
        with open(path, "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows(self.rows(False))


def compare(params: list, labels: list | None = None) -> Comparison:
    """Compare a series (SpectrumParams or None per spectrum) against its
    majority. Fewer than two Bruker rows -> no reports, level 'none'."""
    params = list(params or [])
    if labels is None:
        labels = [p.sample if p is not None else "" for p in params]
    labels = [str(x) for x in labels]
    if sum(1 for p in params if p is not None) < 2:
        return Comparison(params, labels, [])
    reports = []
    for rule in KEY_RULES:
        rep = _report(rule, params)
        if rep is not None:
            reports.append(rep)
    _bad_messages({r.key: r for r in reports})
    return Comparison(params, labels, reports)


# ------------------------------------------------------------------ the pipeline
def template_from(params: SpectrumParams) -> dict:
    """The processing-shape template one spectrum's procs describe."""
    pr = params.proc
    return {"wdw": int(pr.get("WDW") or 0), "lb_hz": float(pr.get("LB") or 0.0),
            "gb": float(pr.get("GB") or 0.0), "ssb": float(pr.get("SSB") or 0.0),
            "tdeff": int(pr.get("TDeff") or 0), "si": int(pr.get("SI") or 0),
            "fcor": float(pr.get("FCOR") if pr.get("FCOR") is not None else 1.0)}


def majority_template(cmp: Comparison) -> dict:
    """The template built from the series' majority values."""
    def maj(key, default):
        r = cmp.report(key)
        return default if r is None or r.majority is None else r.majority
    return {"wdw": int(maj("WDW", 0)), "lb_hz": float(maj("LB", 0.0)),
            "gb": float(maj("GB", 0.0)), "ssb": float(maj("SSB", 0.0)),
            "tdeff": int(maj("TDeff", 0)), "si": int(maj("SI", 0)),
            "fcor": float(maj("FCOR", 1.0))}


def ops_for(params: SpectrumParams, template: dict | None = None,
            phase: str = "own") -> list:
    """The LARMOR processing chain that rebuilds ``params``' spectrum from its
    fid: ``template`` (shape keys wdw / lb_hz / gb / ssb / tdeff / si / fcor,
    defaults = the spectrum's own procs) then ``ft`` referenced by the
    spectrum's OWN SF, then its own TopSpin phase (``phase='own'``),
    ``autophase`` or nothing. Only existing :data:`larmor.processing.OPS`."""
    t = template_from(params)
    if template:
        t.update({k: v for k, v in template.items() if v is not None})
    ops: list = []
    td = int(params.acq.get("TD") or 0)
    tdeff = int(t.get("tdeff") or 0)
    if 0 < tdeff and (not td or tdeff < td):
        ops.append({"op": "tdeff", "points": tdeff // 2})   # TopSpin counts real points
    fcor = float(t.get("fcor", 1.0))
    if fcor != 1.0:
        ops.append({"op": "fcor", "factor": fcor})
    wdw, lb = int(t.get("wdw") or 0), float(t.get("lb_hz") or 0.0)
    if wdw == 0:
        pass
    elif wdw == 1:
        if lb:
            ops.append({"op": "em", "lb_hz": lb})
    elif wdw == 2:
        ops.append({"op": "gm", "lb_hz": -abs(lb), "gb": float(t.get("gb") or 0.1)})
    elif wdw in (3, 4):
        ops.append({"op": "sine", "ssb": float(t.get("ssb") or 2.0),
                    "power": 1 if wdw == 3 else 2})
    elif wdw == 9:
        ops.append({"op": "traf", "lb_hz": lb or 10.0})
    else:
        raise ValueError(f"WDW {wdw} ({TOPSPIN_WDW.get(wdw, '?')}) is not supported "
                         "-- choose none, EM, GM or QSINE")
    si = int(t.get("si") or 0)
    if si:
        ops.append({"op": "zf", "si": si})
    offset, _referenced = carrier_ppm({
        "larmor_MHz": params.acq.get("SFO1") or 0.0,
        "sf_MHz": params.proc.get("SF") or 0.0,
        "bf1_MHz": params.acq.get("BF1") or 0.0,
        "o1_Hz": params.acq.get("O1") or 0.0})
    ops.append({"op": "ft", "offset_ppm": float(offset)})
    if phase == "own":
        ph_mod = int(params.proc.get("PH_mod") or 0)
        if ph_mod == 1:
            ops.append({"op": "phase", "p0": -float(params.proc.get("PHC0") or 0.0),
                        "p1": float(params.proc.get("PHC1") or 0.0),
                        "pivot_frac": 1.0})
        elif ph_mod == 2:
            ops.append({"op": "magnitude"})
    elif phase == "auto":
        ops.append({"op": "autophase"})
    elif phase not in ("none", "", None):
        raise ValueError(f"unknown phase mode {phase!r} (own / auto / none)")
    return ops


def describe_template(t: dict, phase: str = "own") -> str:
    """'TDeff 1024 · EM 0 Hz · SI 32k · phase: own PHC0/PHC1'."""
    wdw, lb = int(t.get("wdw") or 0), float(t.get("lb_hz") or 0.0)
    if wdw == 0:
        win = "no window"
    elif wdw == 1:
        win = f"EM {_g(lb)} Hz"
    elif wdw == 2:
        win = f"GM {_g(lb)} Hz / GB {_g(t.get('gb') or 0.1)}"
    elif wdw in (3, 4):
        win = f"{TOPSPIN_WDW[wdw]} ssb {_g(t.get('ssb') or 2)}"
    else:
        win = TOPSPIN_WDW.get(wdw, f"WDW {wdw}")
    tdeff = int(t.get("tdeff") or 0)
    bits = [f"TDeff {tdeff if tdeff > 0 else 'all'}", win]
    if t.get("si"):
        bits.append(f"SI {_si_fmt(t['si'])}")
    if float(t.get("fcor", 1.0)) != 1.0:
        bits.append(f"FCOR {_g(t['fcor'])}")
    bits.append({"own": "phase: own PHC0/PHC1", "auto": "autophase"}.get(phase, "no phase"))
    return " · ".join(bits)


def reprocess(params: SpectrumParams, ops: list):
    """Rebuild one spectrum from its fid through :func:`larmor.loader.
    apply_processing` -- the very path a saved recipe replays on reopen and
    in the Plotting studio, so what the dialog shows is what reopening
    reproduces. Returns ``(ppm, amp, notes)``; raises when there is no fid."""
    if not params.has_fid:
        raise ValueError(f"{params.sample}/{Path(params.expno).name}: no raw fid "
                         f"in {params.expno}")
    from larmor.loader import apply_processing
    from larmor.recipe import Recipe

    rec = Recipe.from_dict({
        "processing": [dict(o) for o in ops], "processing_from_raw": True,
        "larmor_frequency_MHz": float(params.acq.get("SFO1") or 0.0),
        "nucleus": params.acq.get("NUC1") or ""})
    empty = np.zeros(0)
    ppm, amp, notes = apply_processing(rec, empty, empty, source_path=params.expno)
    return np.asarray(ppm, float), np.asarray(amp, float), list(notes)


#: the procs keys a reprocess-from-fid rewrites: the shape keys become the
#: template's, and no `abs` / baseline command is replayed
_REPROCESSED_PROC = {"WDW": "wdw", "LB": "lb_hz", "GB": "gb", "SSB": "ssb",
                     "TDeff": "tdeff", "SI": "si", "FCOR": "fcor"}


def as_reprocessed(params: SpectrumParams, template: dict) -> SpectrumParams:
    """``params`` as they read AFTER a reprocess with ``template``: the shape
    keys take the template's values and ABSG / BC_mod drop to 0, so comparing
    a partly reprocessed series flags exactly the members that still carry
    their TopSpin processing (no fid) and nothing else."""
    proc = dict(params.proc)
    for key, tkey in _REPROCESSED_PROC.items():
        if tkey in template and template[tkey] is not None:
            proc[key] = template[tkey]
    proc["ABSG"] = 0
    proc["BC_mod"] = 0
    return replace(params, proc=proc)


def caveat_note(cmp: Comparison, k: int) -> str:
    """The note a recipe carries when its spectrum was NOT reprocessed and
    differs from the series majority; '' otherwise."""
    reps = cmp._deviant_reports(k)
    if not reps:
        return ""
    bits = []
    for r in reps:
        rule = _RULES[r.key]
        unit = f" {rule.unit}" if rule.unit else ""
        bits.append(f"{r.key} {r.display[k]}{unit} (majority {r.majority_display}{unit})")
    return ("acquired/processed differently from the series majority: "
            + "; ".join(bits))
