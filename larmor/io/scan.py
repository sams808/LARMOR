"""Scan a Bruker sample folder and auto-identify every spectrum in it.

A "sample" in TopSpin is a folder holding many numbered EXPNOs. This reads the
acqus of each (read-only, no data load) and reports what it is: nucleus,
1D/2D, the experiment kind inferred from the pulse program, and which data
files are available (fid / ser / 1r / 2rr). Powers the left explorer panel and
the "Open sample" action.

The same module owns the **sample identity** rule (``sample_name`` and
friends): the name a spectrum is known by everywhere -- plot title, workspace
row, default save names, batch scopes -- derived from the sample folder
(``MMDDYYYY_<rotor>_<sample>_SS_ALP`` in the group's instrument tree) or, for
EXPNO-per-sample layouts, from the title's ``Sample <name>`` line; never from
a pulse note such as "11B with short tip angle".
"""
from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "ExperimentInfo", "TreeEntry", "SampleName", "read_experiment",
    "is_sample_folder", "scan_sample", "list_dir", "name_parts", "sample_key",
    "sample_name", "sample_label", "disambiguate", "DATE_PREFIX_RE",
    "ROTOR_PREFIX_RE", "SUFFIXES", "SESSION_DIR_RE", "TITLE_SAMPLE_RE",
]

#: pulse-program fragments -> a human experiment kind (checked in order)
_KIND_HINTS = (
    ("satrec", "Saturation recovery (T1)"),
    ("satrect", "Saturation recovery (T1)"),
    ("t1ir", "Inversion recovery (T1)"),
    ("invrec", "Inversion recovery (T1)"),
    ("t1rho", "T1ρ"),
    ("cpmg", "CPMG (T2)"),
    ("hahnecho", "Hahn echo"),
    ("spinecho", "Spin echo"),
    ("echo", "Echo"),
    ("mp3q", "MQMAS"),
    ("mqmas", "MQMAS"),
    ("3qmas", "MQMAS"),
    ("mqvas", "MQMAS"),
    ("5qmas", "5Q-MAS"),
    ("stmas", "ST-MAS"),
    ("redor", "REDOR"),
    ("cpmas", "CP-MAS"),
    ("hetcor", "HETCOR"),
    ("cp", "Cross-polarization"),
    ("mas", "MAS"),
    ("hpdec", "Single pulse (dec.)"),
    ("onepulse", "Single pulse"),
    ("zg", "Single pulse"),
    ("bp", "Bloch decay"),
)


@dataclass
class ExperimentInfo:
    expno: str                       # folder name (e.g. "2702")
    path: str                        # absolute EXPNO path
    nucleus: str                     # NUC1, e.g. "27Al"
    ndim: int                        # 1 or 2
    kind: str                        # human experiment kind
    pulse_program: str
    title: str
    has_fid: bool
    has_ser: bool
    has_1r: bool
    has_2rr: bool
    #: pdata/1/ct1t2.txt exists: TopSpin's t1/t2 analysis was run here
    has_t1_fit: bool = False
    # acquisition facts the session inventory reads (all defaulted so the
    # positional constructors above keep working)
    ns: int = 0                      # acqus NS (0 if absent)
    d1_s: float = 0.0                # acqus D[1], the recycle delay in s
    date: float = 0.0                # acqus DATE, epoch seconds (0 if absent)
    n_procs: int = 0                 # numeric pdata/<procno> folders
    has_popt: bool = False           # popt.array present (parameter optimisation)
    title_full: str = ""             # the whole title text, stripped

    @property
    def date_iso(self) -> str:
        if not self.date:
            return ""
        return _dt.datetime.fromtimestamp(self.date).strftime("%Y-%m-%d %H:%M")

    @property
    def openable(self) -> str | None:
        """The best default target path for opening this experiment."""
        p = Path(self.path)
        if self.has_2rr:
            return str(p / "pdata" / "1" / "2rr")
        if self.has_1r:
            return str(p / "pdata" / "1" / "1r")
        if self.has_ser:
            return str(p / "ser")
        if self.has_fid:
            return str(p / "fid")
        return None

    @property
    def label(self) -> str:
        dim = "2D" if self.ndim == 2 else "1D"
        return f"{self.expno} · {self.nucleus} {dim} · {self.kind}"

    @property
    def has_raw(self) -> bool:
        return self.has_fid or self.has_ser


def _classify(pulprog: str) -> str:
    p = (pulprog or "").lower()
    for frag, kind in _KIND_HINTS:
        if frag in p:
            return kind
    return pulprog or "unknown"


def _read_acqus_min(expno: Path) -> dict:
    """Parse just the few acqus keys we need, without loading data: NUC1,
    PULPROG, PARMODE, NS (int), DATE (float) and D (the ``##$D= (0..63)``
    delay array -- its values follow on the next lines until the next ``##``
    key; D[1] is the recycle delay). Still a pure text parse, no nmrglue."""
    keys = {"NUC1": "", "PULPROG": "", "PARMODE": None, "NS": 0, "DATE": 0.0,
            "D": []}
    text = (expno / "acqus").read_text(errors="replace")
    in_d = False
    for line in text.splitlines():
        if in_d:
            if line.startswith("##"):
                in_d = False
            else:
                for tok in line.split():
                    try:
                        keys["D"].append(float(tok))
                    except ValueError:
                        pass
                continue
        if line.startswith("##$NUC1="):
            keys["NUC1"] = line.split("=", 1)[1].strip().strip("<>")
        elif line.startswith("##$PULPROG="):
            keys["PULPROG"] = line.split("=", 1)[1].strip().strip("<>")
        elif line.startswith("##$PARMODE="):
            try:
                keys["PARMODE"] = int(line.split("=", 1)[1])
            except ValueError:
                pass
        elif line.startswith("##$NS="):
            try:
                keys["NS"] = int(float(line.split("=", 1)[1]))
            except ValueError:
                pass
        elif line.startswith("##$DATE="):
            try:
                keys["DATE"] = float(line.split("=", 1)[1])
            except ValueError:
                pass
        elif line.startswith("##$D= ("):
            in_d = True
    return keys


def read_experiment(expno: str | Path) -> ExperimentInfo:
    """Identify one EXPNO folder (read-only, metadata only)."""
    p = Path(expno)
    acq = _read_acqus_min(p)
    is2d = (p / "acqu2s").exists() or (p / "ser").exists() \
        or (p / "pdata" / "1" / "2rr").exists() or (acq["PARMODE"] or 0) >= 1
    title = title_full = ""
    tp = p / "pdata" / "1" / "title"
    if tp.exists():
        text = tp.read_text(errors="replace")
        first = text.splitlines()
        title = first[0][:80] if first else ""
        title_full = text.strip()
    pdata = p / "pdata"
    n_procs = 0
    if pdata.is_dir():
        try:
            n_procs = sum(1 for c in pdata.iterdir() if c.is_dir() and c.name.isdigit())
        except OSError:
            n_procs = 0
    d = acq["D"]
    return ExperimentInfo(
        expno=p.name, path=str(p),
        nucleus=acq["NUC1"], ndim=2 if is2d else 1,
        kind=_classify(acq["PULPROG"]), pulse_program=acq["PULPROG"],
        title=title,
        has_fid=(p / "fid").exists(), has_ser=(p / "ser").exists(),
        has_1r=(p / "pdata" / "1" / "1r").exists(),
        has_2rr=(p / "pdata" / "1" / "2rr").exists(),
        has_t1_fit=(p / "pdata" / "1" / "ct1t2.txt").exists(),
        ns=int(acq["NS"] or 0), d1_s=float(d[1]) if len(d) > 1 else 0.0,
        date=float(acq["DATE"] or 0.0), n_procs=n_procs,
        has_popt=(p / "popt.array").exists(), title_full=title_full)


def is_sample_folder(folder: str | Path) -> bool:
    """True if the folder directly contains at least one EXPNO."""
    p = Path(folder)
    if not p.is_dir():
        return False
    for child in p.iterdir():
        if child.is_dir() and (child / "acqus").exists():
            return True
    return False


def scan_sample(folder: str | Path) -> list[ExperimentInfo]:
    """Every EXPNO in a sample folder, sorted by number then name."""
    p = Path(folder)
    out = []
    for child in sorted(p.iterdir(), key=_expno_sort_key):
        if child.is_dir() and (child / "acqus").exists():
            try:
                out.append(read_experiment(child))
            except Exception:
                continue
    return out


#: experiment kinds that measure T1 (the sibling the quantitativity chip reads)
T1_KINDS = ("Saturation recovery (T1)", "Inversion recovery (T1)")


def find_sibling_t1(expno_dir: str | Path, nucleus: str, *,
                    kinds: tuple = T1_KINDS) -> list[ExperimentInfo]:
    """The T1 measurements of ``nucleus`` in the sample folder that holds
    ``expno_dir`` (the EXPNO itself excluded), best first: those with a
    TopSpin fit (``has_t1_fit``) before those without, then the nearest EXPNO
    number, then the number itself; non-numeric EXPNO names sort last. []
    when the parent is not a sample folder."""
    p = Path(expno_dir)
    parent = p.parent
    if not nucleus or not is_sample_folder(parent):
        return []
    try:
        target = int(p.name)
    except ValueError:
        target = None

    def _key(info: ExperimentInfo):
        try:
            n = int(info.expno)
        except ValueError:
            return (not info.has_t1_fit, 1, float("inf"), info.expno.lower())
        dist = abs(n - target) if target is not None else 0
        return (not info.has_t1_fit, 0, dist, n)

    out = [info for info in scan_sample(parent)
           if info.nucleus == nucleus and info.kind in kinds
           and Path(info.path).resolve() != p.resolve()]
    return sorted(out, key=_key)


def _expno_sort_key(p: Path):
    return (0, int(p.name)) if p.name.isdigit() else (1, p.name.lower())


# ---------------------------------------------------------------- sample identity
#: the group's folder convention is ``MMDDYYYY_<rotor>_<sample>_SS_ALP``:
#: an acquisition date (one tree carries a 9-digit "031172026_" typo), an
#: optional rotor token (SR31648, RS40175, SR31602LS), the sample, and the
#: operators' initials as one suffix (_SS_ALP, _ALP_SS, _SS or _ALP)
DATE_PREFIX_RE = re.compile(r"^\d{8,9}_")
ROTOR_PREFIX_RE = re.compile(r"^(?:SR|RS)\d{4,}[A-Z]*_")
SUFFIXES = ("_SS_ALP", "_ALP_SS", "_SS", "_ALP")
#: a session (month) folder: ``2026-05``, ``2026-06_35Cl``
SESSION_DIR_RE = re.compile(r"^\d{4}-\d{2}")
#: the MagLab convention puts the sample on its own title line, "Sample X"
TITLE_SAMPLE_RE = re.compile(r"^\s*sample\b[\s:=-]*(\S.*)$", re.IGNORECASE)


@dataclass(frozen=True)
class SampleName:
    """A sample's identity and where it came from.

    ``source`` says which rule produced ``key``: ``folder`` (the tokenised
    folder name), ``month-folder`` (a plain folder directly under a session
    month), ``title-sample`` (a ``Sample <name>`` title line), ``title`` (the
    title's first line) or ``fallback`` (the raw folder name)."""
    key: str
    folder: str
    date: str = ""
    rotor: str = ""
    suffix: str = ""
    title_first: str = ""
    source: str = "folder"

    @property
    def from_folder(self) -> bool:
        """True when the folder name carried at least one convention token."""
        return bool(self.date or self.rotor or self.suffix)


def name_parts(folder) -> SampleName:
    """Split a sample folder name into date, rotor, key and suffix.

    ``01192026_SR31649_Base0Ca_SS_ALP`` -> key ``Base0Ca``; a name without any
    of the tokens keeps itself as the key. Accepts a path or a bare name."""
    raw = Path(str(folder)).name
    rest = raw
    date = rotor = suffix = ""
    m = DATE_PREFIX_RE.match(rest)
    if m:
        date = m.group(0)[:-1]
        rest = rest[m.end():]
    m = ROTOR_PREFIX_RE.match(rest)
    if m:
        rotor = m.group(0)[:-1]
        rest = rest[m.end():]
    for suf in SUFFIXES:                       # ONE suffix, longest first
        if rest.endswith(suf) and len(rest) > len(suf):
            suffix = suf
            rest = rest[:-len(suf)]
            break
    return SampleName(key=rest or raw, folder=raw, date=date, rotor=rotor,
                      suffix=suffix, source="folder")


def sample_key(folder) -> str:
    """The sample key of a folder name (see ``name_parts``)."""
    return name_parts(folder).key


def sample_name(expno_dir, title: str = "") -> SampleName:
    """The sample an EXPNO belongs to, by a layered rule:

    (a) the EXPNO's parent folder follows the group convention (a date, rotor
        or ``_SS_ALP`` token was stripped) -> the remaining key;
    (b) else the parent's parent looks like a session month -> the folder
        name as is (``2025-12/Cryolite/16`` -> ``Cryolite``);
    (c) else a title line reads ``Sample <name>`` -> that name (the MagLab
        EXPNO-per-sample layout);
    (d) else the title's first non-empty line, else the folder name.

    A display name given to the sample folder in LARMOR (``larmor.aliases``,
    the Explorer's Rename…) comes before every rule, with ``source ==
    "alias"``. Rules (a) and (b) never consult the title; its first line is
    recorded in ``title_first`` for provenance in every case."""
    from larmor import aliases

    p = Path(expno_dir)
    folder = p.parent
    lines = [ln.strip() for ln in (title or "").splitlines()]
    first = next((ln for ln in lines if ln), "")
    parts = name_parts(folder.name)
    alias = aliases.alias_for(folder)
    if alias:
        return SampleName(alias, folder.name, parts.date, parts.rotor,
                          parts.suffix, first, "alias")
    if parts.from_folder:
        return SampleName(parts.key, folder.name, parts.date, parts.rotor,
                          parts.suffix, first, "folder")
    if SESSION_DIR_RE.match(folder.parent.name or ""):
        return SampleName(folder.name, folder.name, "", "", "", first, "month-folder")
    for ln in lines:
        m = TITLE_SAMPLE_RE.match(ln)
        if m and m.group(1).strip():
            return SampleName(m.group(1).strip(), folder.name, "", "", "", first,
                              "title-sample")
    if first:
        return SampleName(first, folder.name, "", "", "", first, "title")
    return SampleName(folder.name, folder.name, "", "", "", "", "fallback")


_DATA_FILES = ("1r", "2rr", "fid", "ser", "1i", "2ii")


def _expno_dir_of(path) -> Path | None:
    """The first ancestor (or the path itself) that holds an ``acqus``."""
    p = Path(path)
    try:
        candidates = [p] + list(p.parents)
    except (TypeError, ValueError):
        return None
    for c in candidates[:6]:
        try:
            if (c / "acqus").exists():
                return c
        except OSError:
            return None
    return None


def _read_title(expno_dir: Path) -> str:
    tp = expno_dir / "pdata" / "1" / "title"
    try:
        return tp.read_text(errors="replace") if tp.exists() else ""
    except OSError:
        return ""


def sample_label(path, rec) -> str:
    """A meaningful sample name for a spectrum: the recipe's sample if it is not
    just the nucleus, else the sample derived from the path -- through
    ``sample_name`` when the EXPNO is on disk, else from the first path
    segment that is not a proc/EXPNO number, tokenised by ``name_parts`` -- so
    a title of "31P" becomes the real sample name and a Bruker file is never
    labelled "1r". A display name given in LARMOR to the EXPNO, else to its
    sample folder (``larmor.aliases``), comes before the recipe's sample."""
    from larmor import aliases

    rec = rec or {}
    nucleus = (rec.get("nucleus") or "").strip()
    name = (rec.get("sample") or "").strip()
    expno_dir = _expno_dir_of(path)
    if expno_dir is not None:
        alias = aliases.alias_for(expno_dir) or aliases.alias_for(expno_dir.parent)
        if alias:
            return alias
    if name and name.lower() != nucleus.lower():
        return name
    if expno_dir is not None:
        return sample_name(expno_dir, _read_title(expno_dir)).key
    parts = Path(str(path)).parts
    for k, seg in enumerate(reversed(parts)):
        low = seg.lower()
        if low.endswith(".fid"):                 # a Varian dataset folder
            return seg[:-4]
        if seg == "pdata" or seg.isdigit() or low in _DATA_FILES:
            continue
        if k == 0 and "." in seg and not seg.startswith("."):
            stem = Path(seg).stem                # a file (csv, fxmla, recipe)
            return stem[:-7] if stem.lower().endswith(".recipe") else stem
        return name_parts(seg).key
    return name or Path(str(path)).stem


def _folder_and_expno(path) -> tuple[str, str]:
    """(sample folder name, EXPNO) of a Bruker-ish path, from disk when the
    EXPNO exists, else from the path's segments."""
    expno_dir = _expno_dir_of(path)
    if expno_dir is not None:
        return expno_dir.parent.name, expno_dir.name
    parts = Path(str(path)).parts
    if "pdata" in parts:
        i = parts.index("pdata")
        return (parts[i - 2] if i >= 2 else ""), (parts[i - 1] if i >= 1 else "")
    if parts and parts[-1].lower() in _DATA_FILES and len(parts) >= 3:
        return parts[-3], parts[-2]
    if parts and parts[-1].isdigit() and len(parts) >= 2:
        return parts[-2], parts[-1]
    return (parts[-2] if len(parts) >= 2 else ""), (Path(str(path)).stem if parts else "")


def disambiguate(labels: list[str], paths: list[str]) -> list[str]:
    """Make equal labels distinct: the same name from different sample folders
    gets the folder's date token (else the folder name) appended as
    " (04272026)"; the same name from one folder gets " · <EXPNO>". Unique
    labels are returned unchanged."""
    labels = list(labels)
    if len(labels) != len(paths):
        return labels
    where = [_folder_and_expno(p) for p in paths]
    groups: dict[str, list[int]] = {}
    for i, lab in enumerate(labels):
        groups.setdefault(lab, []).append(i)
    for lab, idx in groups.items():
        if len(idx) < 2:
            continue
        folders = sorted({where[i][0] for i in idx})
        if len(folders) > 1:
            tags = {f: (name_parts(f).date or f) for f in folders}
            if len(set(tags.values())) < len(folders):   # same day twice
                tags = {f: f for f in folders}
            for i in idx:
                tag = tags[where[i][0]]
                if tag:
                    labels[i] = f"{lab} ({tag})"
    # still equal (same folder, or folders without a tag) -> the EXPNO
    groups = {}
    for i, lab in enumerate(labels):
        groups.setdefault(lab, []).append(i)
    for lab, idx in groups.items():
        if len(idx) < 2:
            continue
        for i in idx:
            if where[i][1]:
                labels[i] = f"{lab} · {where[i][1]}"
    return labels


@dataclass
class TreeEntry:
    name: str
    path: str
    is_dir: bool
    is_sample: bool = False
    is_expno: bool = False
    info: ExperimentInfo | None = None
    children: list = field(default_factory=list)


def list_dir(folder: str | Path) -> list[TreeEntry]:
    """One level of the filesystem for the explorer: subfolders (flagged if
    they are samples), and EXPNOs identified inline."""
    p = Path(folder)
    entries = []
    try:
        children = sorted(p.iterdir(), key=_expno_sort_key)
    except (PermissionError, OSError):
        return entries
    for child in children:
        if not child.is_dir():
            continue
        if (child / "acqus").exists():
            try:
                info = read_experiment(child)
            except Exception:
                info = None
            entries.append(TreeEntry(child.name, str(child), True,
                                     is_expno=True, info=info))
        else:
            entries.append(TreeEntry(child.name, str(child), True,
                                     is_sample=is_sample_folder(child)))
    return entries
