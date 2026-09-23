"""Session inventory: every EXPNO of a month folder as a sample x nucleus grid
with the production spectrum pre-picked per block.

A session is one month folder of the instrument tree
(``DATA/2026-05/<sample>/<EXPNO>``); an EXPNO-per-sample set such as the
MagLab ``35Cl_2025-12/<EXPNO>`` layout is a session too. Every EXPNO becomes
a ``Row`` -- the ``larmor.io.scan`` identification (nucleus, 1D/2D, kind, NS,
D1, date, has-1r, proc count) plus its sample identity
(``scan.sample_name``) -- and rows are grouped in **blocks** keyed
``(folder, sample key, nucleus)``.

The production pick follows the operator's habit: the highest EXPNO of the
block that has a ``pdata/1/1r``, demoted when its NS is a small fraction of
the block's maximum (a last quick shot) or when the title says the shot was
a power check / optimisation / test / failure. Demotion is two-tiered: soft
words (``SOFT_TITLE_RE``), a ``popt.array`` or a 2D pulse program acquired
as 1D only steer the pick when an unflagged candidate exists (copy-pasted
titles are common); ``failed`` / ``trash`` / ``abort`` / ``stopped`` always
demote. Roles: production / candidate / short / setup / failed / arrayed /
2D / reference / unprocessed.

Title-vs-folder **flags** are advisory and never move a pick: a rotor ID in
the title that contradicts the folder's, a first title line naming another
sample of the session, a leading nucleus token that is not NUC1, ``zg`` in
the title of a non-zg pulse program, two folders sharing one sample name,
and an unprocessed EXPNO with at least the pick's NS.

Qt-free; the desktop window (larmor/desktop/inventory_dialog.py) and the
``larmor inventory`` command are consumers. Instrument folders are never
written to.
"""
from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from larmor import referencing
from larmor.io import scan

__all__ = [
    "DEFAULT_NS_FRAC", "HARD_TITLE_RE", "SOFT_TITLE_RE", "ROTOR_TOKEN_RE",
    "TITLE_NUCLEUS_RE", "SETUP_KINDS_1D", "RELAX_KINDS", "ROLES", "Row",
    "Inventory", "inventory_root", "classify_title", "assign_roles",
    "flag_titles", "build", "join_sr", "to_csv", "picks_text", "grid_text",
    "text_report",
]

#: a pick whose NS is below this fraction of the block's maximum is a short shot
DEFAULT_NS_FRAC = 0.25
#: title words that always demote (tier 0)
HARD_TITLE_RE = re.compile(r"\bfail(?:ed|ure)?\b|\btrash\b|\babort|\bstopped\b",
                           re.IGNORECASE)
#: title words that demote only when an unflagged candidate exists (tier 2)
SOFT_TITLE_RE = re.compile(
    r"power check|power optimi[sz]ation|\bopt\b|\btest\b|set(?:ting)?[ -]?up|"
    r"\binitial\b|check empty rotor|non-?quantitat|\bdraft\b|\bbackground\b|referenc",
    re.IGNORECASE)
#: "Rotor SR31648" in a title (the folder carries the same token)
ROTOR_TOKEN_RE = re.compile(r"\bRotor[ :#]*((?:SR|RS)\d{4,}[A-Z]*)", re.IGNORECASE)
#: a title that starts with a nucleus token ("27Al zg power check")
TITLE_NUCLEUS_RE = re.compile(r"^\s*(\d{1,3}[A-Z][a-z]?)(?![A-Za-z0-9])")
#: 2D pulse programs sometimes run as 1D to set them up -- a soft demotion
SETUP_KINDS_1D = ("MQMAS", "5Q-MAS", "ST-MAS")
#: 'arrayed' when ndim == 2 only (a QCPMG 1D is a production spectrum)
RELAX_KINDS = ("Saturation recovery (T1)", "Inversion recovery (T1)", "T1ρ",
               "CPMG (T2)", "REDOR")
ROLES = ("production", "candidate", "short", "setup", "failed", "arrayed", "2D",
         "reference", "unprocessed")
_ZG_RE = re.compile(r"\bzg\b", re.IGNORECASE)


# ---------------------------------------------------------------- rows
@dataclass
class Row:
    info: scan.ExperimentInfo
    folder: str                      # the EXPNO's parent folder name
    name: scan.SampleName
    role: str = "candidate"
    reasons: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    pick: bool = False
    sr_verdict: str = ""
    sr_note: str = ""

    @property
    def expno(self) -> int:
        try:
            return int(self.info.expno)
        except ValueError:
            return 0

    @property
    def nucleus(self) -> str:
        return self.info.nucleus

    @property
    def path(self) -> str:
        return self.info.path

    @property
    def openable(self) -> str | None:
        return self.info.openable

    @property
    def ns(self) -> int:
        return self.info.ns

    @property
    def d1_s(self) -> float:
        return self.info.d1_s

    @property
    def date_iso(self) -> str:
        return self.info.date_iso

    @property
    def title_line(self) -> str:
        for ln in (self.info.title_full or self.info.title or "").splitlines():
            if ln.strip():
                return ln.strip()
        return ""

    @property
    def sample(self) -> str:
        return self.name.key

    @property
    def sample_id(self) -> tuple[str, str]:
        return (self.folder, self.name.key)

    @property
    def block_id(self) -> tuple[str, str, str]:
        return (self.folder, self.name.key, self.info.nucleus)

    @property
    def label(self) -> str:
        return f"{self.folder}/{self.info.expno}"


def _natural(s: str) -> list:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s or "")]


def _mass(nucleus: str) -> int:
    m = re.match(r"(\d+)", nucleus or "")
    return int(m.group(1)) if m else 10**6


# ---------------------------------------------------------------- the inventory
@dataclass
class Inventory:
    root: Path
    rows: list[Row]
    ns_frac: float = DEFAULT_NS_FRAC

    def blocks(self) -> dict[tuple[str, str, str], list[Row]]:
        """(folder, sample key, nucleus) -> rows in EXPNO order."""
        out: dict[tuple[str, str, str], list[Row]] = {}
        for r in self.rows:
            out.setdefault(r.block_id, []).append(r)
        for rows in out.values():
            rows.sort(key=lambda r: (r.expno, r.info.expno))
        return out

    def samples(self) -> list[tuple[str, str]]:
        """(folder, key) pairs, natural-sorted by key, then folder date, then folder."""
        seen = {r.sample_id for r in self.rows}
        return sorted(seen, key=lambda s: (_natural(s[1]),
                                           scan.name_parts(s[0]).date, s[0].lower()))

    def labels(self) -> dict[tuple[str, str], str]:
        """Display label per sample: the key, with the folder's date token (else
        the folder name) appended when two folders of the session share a key."""
        folders: dict[str, list[str]] = {}
        for folder, key in self.samples():
            folders.setdefault(key, []).append(folder)
        out = {}
        for folder, key in self.samples():
            if len(folders[key]) > 1:
                tags = [scan.name_parts(f).date or f for f in folders[key]]
                if len(set(tags)) < len(tags):          # same day twice: the folder
                    tags = list(folders[key])
                out[(folder, key)] = f"{key} ({tags[folders[key].index(folder)]})"
            else:
                out[(folder, key)] = key
        return out

    def cell(self, sample: tuple[str, str], nucleus: str) -> Row | None:
        """What the grid shows for (sample, nucleus): the pick, else -- for
        1H, which is never picked -- the highest reference spectrum, else None."""
        rows = self.candidates(sample, nucleus)
        pick = next((r for r in rows if r.pick), None)
        if pick is not None or nucleus != "1H":
            return pick
        refs = [r for r in rows if r.role == "reference"]
        return max(refs, key=lambda r: r.expno) if refs else None

    def nuclei(self) -> list[str]:
        """The 1D nuclei of the session, mass-number order, 1H last."""
        nucs = {r.nucleus for r in self.rows if r.info.ndim == 1 and r.nucleus}
        return sorted(nucs, key=lambda n: (n == "1H", _mass(n), n))

    def candidates(self, sample: tuple[str, str], nucleus: str) -> list[Row]:
        return self.blocks().get((sample[0], sample[1], nucleus), [])

    def grid(self) -> dict[tuple[tuple[str, str], str], Row | None]:
        """The pick per (sample, nucleus); None where the block has none."""
        blocks = self.blocks()
        out = {}
        for s in self.samples():
            for n in self.nuclei():
                rows = blocks.get((s[0], s[1], n), [])
                out[(s, n)] = next((r for r in rows if r.pick), None)
        return out

    def picks(self, nucleus: str | None = None) -> list[Row]:
        """The picked rows in grid order (samples, then nuclei)."""
        g = self.grid()
        out = []
        for s in self.samples():
            for n in self.nuclei():
                if nucleus and n != nucleus:
                    continue
                r = g.get((s, n))
                if r is not None:
                    out.append(r)
        return out

    def set_pick(self, sample: tuple[str, str], nucleus: str, expno: int | None) -> None:
        """Manual override: make ``expno`` the block's only pick (reason
        'chosen by user'); ``None`` clears the block's pick."""
        rows = self.candidates(sample, nucleus)
        target = None
        if expno is not None:
            target = next((r for r in rows if r.expno == int(expno)), None)
            if target is None:
                raise ValueError(f"EXPNO {expno} is not in block {sample} / {nucleus}")
        for r in rows:
            if r.pick:
                r.pick = False
                if r.role == "production":
                    r.role = "candidate"
        if target is not None:
            target.pick = True
            target.role = "production"
            if "chosen by user" not in target.reasons:
                target.reasons.append("chosen by user")

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for r in self.rows:
            out[r.role] = out.get(r.role, 0) + 1
        return out

    def n_flags(self) -> int:
        """Rows carrying at least one title-vs-folder flag."""
        return sum(1 for r in self.rows if r.flags)

    def n_folders(self) -> int:
        return len({r.folder for r in self.rows})


# ---------------------------------------------------------------- the rules
def inventory_root(path) -> Path:
    """The session folder a path belongs to: an EXPNO, a pdata folder or an
    NMRFAM sample folder -> its month; an EXPNO-per-sample set (a folder of
    EXPNOs whose name carries none of the convention tokens and whose parent
    is not a month) -> itself; anything else is taken as the session."""
    p = Path(path)
    if p.is_file():
        p = p.parent
    if (p / "acqus").exists():                        # EXPNO -> its folder
        p = p.parent
    elif p.name == "pdata":
        p = p.parent.parent
    elif p.parent.name == "pdata":                    # pdata/<procno>
        p = p.parent.parent.parent
    if scan.is_sample_folder(p):
        if scan.name_parts(p.name).from_folder or \
                scan.SESSION_DIR_RE.match(p.parent.name or ""):
            return p.parent
        return p
    return p


def _first_line(title: str) -> str:
    for ln in (title or "").splitlines():
        if ln.strip():
            return ln.strip()
    return ""


def classify_title(title: str) -> str:
    """'failed' for a hard word anywhere in the title, 'setup' for a soft word
    on its FIRST line (the operator's note lines say things like "SR is set
    for IUPAC referencing" under a production spectrum), '' otherwise."""
    t = title or ""
    if HARD_TITLE_RE.search(t):
        return "failed"
    if SOFT_TITLE_RE.search(_first_line(t)):
        return "setup"
    return ""


def _pct(ns: int, ns_max: int) -> str:
    if not ns_max:
        return "0"
    pct = 100.0 * ns / ns_max
    return f"{pct:.1f}" if 0 < pct < 1 else f"{pct:.0f}"


def assign_roles(rows: list[Row], ns_frac: float = DEFAULT_NS_FRAC) -> None:
    """Roles and the production pick for every block of ``rows``."""
    blocks: dict[tuple[str, str, str], list[Row]] = {}
    for r in rows:
        r.role, r.reasons, r.pick = "candidate", [], False
        blocks.setdefault(r.block_id, []).append(r)
    for block in blocks.values():
        block.sort(key=lambda r: (r.expno, r.info.expno))
        tier1, tier2 = [], []
        # popt.array tells rows apart only when some rows of the block lack
        # it (one spectrometer leaves it in every EXPNO of a session)
        oned = [r for r in block if r.info.ndim == 1 and r.info.has_1r]
        popt_discriminates = any(r.info.has_popt for r in oned) and \
            not all(r.info.has_popt for r in oned)
        for r in block:
            info = r.info
            if info.ndim == 2:
                r.role = "arrayed" if info.kind in RELAX_KINDS else "2D"
                continue
            if not info.has_1r:
                r.role = "unprocessed"
                continue
            if r.nucleus == "1H":
                r.role = "reference"
                continue
            m = HARD_TITLE_RE.search(info.title_full or info.title or "")
            if m:
                r.role = "failed"
                r.reasons.append(f"title says '{m.group(0)}'")
                continue
            soft = []
            m = SOFT_TITLE_RE.search(_first_line(info.title_full or info.title))
            if m:
                soft.append(f"title says '{m.group(0)}'")
            if info.has_popt and popt_discriminates:
                soft.append("popt.array present (parameter optimisation)")
            if info.kind in SETUP_KINDS_1D:
                soft.append("2D pulse program acquired as 1D")
            if soft:
                r.role = "setup"
                r.reasons.extend(soft)
                tier2.append(r)
            else:
                tier1.append(r)
        tier = tier1 or tier2
        if not tier:
            continue
        ns_max = max(r.ns for r in tier)
        eligible = []
        for r in tier:
            if ns_max > 0 and r.ns < ns_frac * ns_max:
                # a tier-2 row keeps its telling 'setup' role; the NS reason is added
                if r.role != "setup":
                    r.role = "short"
                r.reasons.append(f"NS {r.ns} is {_pct(r.ns, ns_max)} % of the block's {ns_max}")
            else:
                eligible.append(r)
        if not eligible:
            continue
        prod = max(eligible, key=lambda r: r.expno)
        prod.role, prod.pick = "production", True
        if tier is tier2:
            prod.reasons.append("no unflagged spectrum in the block -- picked among "
                                "the flagged ones")
        # say why a higher EXPNO with a 1r was passed over
        for r in block:
            if r.expno > prod.expno and r.info.ndim == 1 and r.info.has_1r \
                    and r.role != "reference":
                if r.role == "short":
                    prod.reasons.append(
                        f"highest EXPNO with pdata/1/1r is {r.expno} but NS {r.ns} is "
                        f"{_pct(r.ns, ns_max)} % of the block's {ns_max} -> demoted; "
                        f"picked {prod.expno}")
                else:
                    why = r.reasons[0] if r.reasons else r.role
                    prod.reasons.append(f"{r.expno} skipped ({r.role}: {why})")


def flag_titles(inv: Inventory) -> None:
    """Advisory title-vs-folder flags on every row (never move a pick)."""
    for r in inv.rows:
        r.flags = []
    samples = inv.samples()
    keys = {key for _, key in samples}
    folders_by_key: dict[str, list[str]] = {}
    for folder, key in samples:
        folders_by_key.setdefault(key, []).append(folder)
    nuclei = {r.nucleus for r in inv.rows if r.nucleus}
    sibling_res = {
        key: re.compile(r"(?<![A-Za-z0-9_])" + re.escape(key) + r"(?![A-Za-z0-9_])",
                        re.IGNORECASE)
        for key in keys if len(key) >= 3 and not SOFT_TITLE_RE.fullmatch(key)
        and not HARD_TITLE_RE.fullmatch(key)}
    for r in inv.rows:
        full = r.info.title_full or r.info.title or ""
        first = r.title_line
        m = ROTOR_TOKEN_RE.search(full)
        if m and r.name.rotor and m.group(1).upper() != r.name.rotor.upper():
            r.flags.append(f"title says rotor {m.group(1)} but the folder is {r.name.rotor}")
        for key, rx in sibling_res.items():
            if key.lower() == r.sample.lower():
                continue
            if rx.search(first):
                r.flags.append(f"title names another sample of the session: {key}")
        m = TITLE_NUCLEUS_RE.match(first)
        if m and m.group(1) in nuclei and m.group(1) != r.nucleus:
            r.flags.append(f"title starts with {m.group(1)} but NUC1 is {r.nucleus}")
        pp = (r.info.pulse_program or "").lower()
        if _ZG_RE.search(first) and "zg" not in pp:
            r.flags.append(f"title says zg but the pulse program is "
                           f"{r.info.pulse_program or 'unknown'}")
    for key, folders in folders_by_key.items():
        if len(folders) < 2:
            continue
        for r in inv.rows:
            if r.pick and r.sample == key:
                others = ", ".join(f for f in folders if f != r.folder)
                r.flags.append(f"sample name '{key}' is shared with folder {others}")
    for block in inv.blocks().values():
        pick = next((r for r in block if r.pick), None)
        if pick is None:
            continue
        unproc = [r for r in block if r.role == "unprocessed" and r.ns > 0
                  and r.ns >= pick.ns]
        if unproc:
            if len({r.ns for r in unproc}) == 1:
                who = ", ".join(str(r.expno) for r in unproc)
                pick.flags.append(f"{who} have NS {unproc[0].ns} but no pdata/1/1r "
                                  "-- process them in TopSpin")
            else:
                who = ", ".join(f"{r.expno} (NS {r.ns})" for r in unproc)
                pick.flags.append(f"{who} have at least the pick's NS but no pdata/1/1r "
                                  "-- process them in TopSpin")


def _row(info: scan.ExperimentInfo, folder: str) -> Row:
    return Row(info=info, folder=folder,
               name=scan.sample_name(info.path, info.title_full or info.title))


def build(path, *, ns_frac: float = DEFAULT_NS_FRAC, sr_audit: bool = True,
          progress=None) -> Inventory:
    """Inventory the session ``path`` belongs to. Folders starting with '~'
    (TopSpin's ~TEMP) are skipped; EXPNOs directly under the root are rows
    too (EXPNO-per-sample sets). ``progress(k, n, name)`` is called per
    folder. With ``sr_audit`` the referencing audit's verdict joins each row
    (a second pass through nmrglue; needs a referenced 1H in the session)."""
    root = inventory_root(path)
    rows: list[Row] = []
    if root.is_dir():
        try:
            children = sorted((c for c in root.iterdir()
                               if c.is_dir() and not c.name.startswith("~")),
                              key=lambda c: (not c.name.isdigit(),
                                             int(c.name) if c.name.isdigit() else 0,
                                             c.name.lower()))
        except OSError:
            children = []
        for k, child in enumerate(children):
            if progress:
                progress(k, len(children), child.name)
            if (child / "acqus").exists():
                try:
                    rows.append(_row(scan.read_experiment(child), root.name))
                except Exception:
                    continue
            else:
                for info in scan.scan_sample(child):
                    rows.append(_row(info, child.name))
    inv = Inventory(root=root, rows=rows, ns_frac=ns_frac)
    assign_roles(inv.rows, ns_frac)
    flag_titles(inv)
    if sr_audit and rows:
        acqs = referencing.scan_session(root)
        if acqs:
            join_sr(inv, referencing.audit(acqs))
    return inv


def _norm(p) -> str:
    return os.path.normcase(os.path.normpath(str(p)))


def join_sr(inv: Inventory, rows: list[referencing.AuditRow]) -> None:
    """Attach the referencing audit's verdict and note to each row (by path)."""
    by_path = {_norm(a.acq.path): a for a in rows}
    for r in inv.rows:
        a = by_path.get(_norm(r.path))
        if a is not None:
            r.sr_verdict = a.verdict
            r.sr_note = "; ".join([a.note] + list(a.notes)).strip("; ")


# ---------------------------------------------------------------- outputs
CSV_HEADER = ["path", "session", "folder", "sample", "label", "expno", "nucleus",
              "ndim", "kind", "pulse_program", "ns", "d1_s", "has_1r", "has_2rr",
              "n_procs", "has_popt", "date", "sr_verdict", "sr_note", "role",
              "reasons", "flags", "pick", "title"]


def to_csv(inv: Inventory, path) -> Path:
    p = Path(path)
    labels = inv.labels()
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(CSV_HEADER)
        for r in sorted(inv.rows, key=lambda r: (r.folder.lower(), r.expno, r.info.expno)):
            i = r.info
            w.writerow([r.path, inv.root.name, r.folder, r.sample,
                        labels.get(r.sample_id, r.sample), i.expno, i.nucleus, i.ndim,
                        i.kind, i.pulse_program, i.ns, f"{i.d1_s:g}", int(i.has_1r),
                        int(i.has_2rr), i.n_procs, int(i.has_popt), i.date_iso,
                        r.sr_verdict, r.sr_note, r.role, "; ".join(r.reasons),
                        "; ".join(r.flags), int(r.pick), r.title_line])
    return p


def picks_text(inv: Inventory, nucleus: str | None = None) -> str:
    """One line per pick: ``sample<TAB>nucleus<TAB>EXPNO<TAB>path``."""
    labels = inv.labels()
    lines = ["# sample\tnucleus\tEXPNO\tpath"]
    for r in inv.picks(nucleus):
        lines.append(f"{labels.get(r.sample_id, r.sample)}\t{r.nucleus}\t{r.expno}\t"
                     f"{r.openable or r.path}")
    return "\n".join(lines) + "\n"


def grid_text(inv: Inventory) -> str:
    """The sample x nucleus grid as fixed-width text: the pick's EXPNO per cell,
    '*' when the block was demoted or flagged, '-' when it has no 1D
    production spectrum."""
    labels = inv.labels()
    samples, nucs = inv.samples(), inv.nuclei()

    def cell(r: Row | None) -> str:
        if r is None:
            return "-"
        if r.role == "reference":
            return f"{r.expno} ref"
        return str(r.expno) + ("*" if (r.reasons or r.flags) else "")

    w0 = max([len("sample")] + [len(labels[s]) for s in samples])
    table = {n: [cell(inv.cell(s, n)) for s in samples] for n in nucs}
    widths = {n: max([len(n)] + [len(c) for c in table[n]]) for n in nucs}
    lines = ["  ".join([f"{'sample':<{w0}}"] + [f"{n:>{widths[n]}}" for n in nucs])]
    for k, s in enumerate(samples):
        lines.append("  ".join([f"{labels[s]:<{w0}}"]
                               + [f"{table[n][k]:>{widths[n]}}" for n in nucs]))
    lines.append("* demoted or flagged block (see the report)   - no 1D production spectrum"
                 "   ref = 1H reference (never picked)")
    return "\n".join(lines) + "\n"


def text_report(inv: Inventory, all_rows: bool = False) -> str:
    """One line per row -- the demoted / flagged ones, or every row."""
    lines = []
    for r in sorted(inv.rows, key=lambda r: (r.folder.lower(), r.expno, r.info.expno)):
        interesting = (r.role in ("short", "setup", "failed") or r.flags
                       or (r.pick and r.reasons))
        if not all_rows and not interesting:
            continue
        mark = "*" if r.pick else " "
        d1 = f"{r.d1_s:g} s" if r.d1_s else "-"
        why = "; ".join(r.reasons + r.flags)
        sr = f"  SR {r.sr_verdict}" if r.sr_verdict else ""
        lines.append(f"{mark} {r.label:<40} {r.nucleus:<6} {r.role:<11} NS {r.ns:>7}  "
                     f"D1 {d1:>8}  {r.title_line[:40]:<40}{sr}"
                     + (f"  | {why}" if why else ""))
    if not lines:
        lines.append("nothing demoted or flagged")
    return "\n".join(lines) + "\n"
