"""Scan a Bruker sample folder and auto-identify every spectrum in it.

A "sample" in TopSpin is a folder holding many numbered EXPNOs. This reads the
acqus of each (read-only, no data load) and reports what it is: nucleus,
1D/2D, the experiment kind inferred from the pulse program, and which data
files are available (fid / ser / 1r / 2rr). Powers the left explorer panel and
the "Open sample" action.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

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
    """Parse just the few acqus keys we need, without loading data."""
    keys = {"NUC1": "", "PULPROG": "", "PARMODE": None}
    text = (expno / "acqus").read_text(errors="replace")
    for line in text.splitlines():
        if line.startswith("##$NUC1="):
            keys["NUC1"] = line.split("=", 1)[1].strip().strip("<>")
        elif line.startswith("##$PULPROG="):
            keys["PULPROG"] = line.split("=", 1)[1].strip().strip("<>")
        elif line.startswith("##$PARMODE="):
            try:
                keys["PARMODE"] = int(line.split("=", 1)[1])
            except ValueError:
                pass
    return keys


def read_experiment(expno: str | Path) -> ExperimentInfo:
    """Identify one EXPNO folder (read-only, metadata only)."""
    p = Path(expno)
    acq = _read_acqus_min(p)
    is2d = (p / "acqu2s").exists() or (p / "ser").exists() \
        or (p / "pdata" / "1" / "2rr").exists() or (acq["PARMODE"] or 0) >= 1
    title = ""
    tp = p / "pdata" / "1" / "title"
    if tp.exists():
        first = tp.read_text(errors="replace").splitlines()
        title = first[0][:80] if first else ""
    return ExperimentInfo(
        expno=p.name, path=str(p),
        nucleus=acq["NUC1"], ndim=2 if is2d else 1,
        kind=_classify(acq["PULPROG"]), pulse_program=acq["PULPROG"],
        title=title,
        has_fid=(p / "fid").exists(), has_ser=(p / "ser").exists(),
        has_1r=(p / "pdata" / "1" / "1r").exists(),
        has_2rr=(p / "pdata" / "1" / "2rr").exists())


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


def _expno_sort_key(p: Path):
    return (0, int(p.name)) if p.name.isdigit() else (1, p.name.lower())


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
