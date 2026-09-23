"""MAS rate from three sources -- acqus MASR, the operator title, the NMRFAM
booking sidecar -- with a majority rule and a per-session confirmation store.

Why three sources. On the NMRFAM spectrometer ``acqus`` records a constant
``MASR`` (4200 Hz in 681 of 694 EXPNOs, 5000 in one session) that has nothing
to do with the rotor: the controller value is a leftover. The operator types
the rate into the title (``MASR 20 kHz``) and the booking system writes the
reserved rate into ``<EXPNO>/experiment_addenda.xml``
(``<magic_angle_spinning_rate>22</...>``, in kHz). Two independent human
entries agreeing within 2 % settle the rate; the third source is reported as
outvoted. All three disagreeing keeps the documented ``max()`` with the red
"MAS -- check!" badge, so the user decides.

Why ``<magic_angle_spinning>false</...>`` is not static evidence. The booking
form's spinning box was left unticked in 214 of 684 sidecars, several of them
MQMAS and 19F rotor-synchronised echoes titled ``MASR 35.714 kHz``. None of
the ``false`` files carries a rate. The flag is informational; the sidecar
contributes a positive rate or nothing.

Why a title rate is bounded. ``MASR 35741 kHz`` (a typo for 35 741 Hz) and
``MASR 35.714 Hz`` (a typo for 35.714 kHz) exist in the tree. A parsed rate
outside 1-150 kHz is re-read in the other unit and marked repaired; a
repaired value never settles the rate alone -- it is flagged unless another
source agrees with it within 2 %.

Why the store lives here and not in the recipe. A confirmation applies to
every later spectrum of the same session folder, rotor and nucleus that
shows the same three source values, in the desktop AND in the loader (Batch
fit, Sequential fit, the CLI). The store is an append-only jsonl under
%LOCALAPPDATA%/LARMOR (``LARMOR_MAS_LOG`` overrides; tests point it at a
temporary file), never inside an instrument folder. The recipe carries the
outcome in ``provenance["mas_rate"]`` so a saved fit says where its rate
came from. The match requires the exact source triplet because one 2026-05
rotor was booked at both 20 and 22 kHz under a constant title of 20 kHz: a
key-only match would cross-apply a confirmation.

Qt-free: stdlib only, plus a lazy ``larmor.referencing.session_root``.
``larmor.io.bruker`` imports this module, never the reverse.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import re
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path

__all__ = [
    "MAS_FALLBACK_HZ", "PLAUSIBLE_HZ", "AGREE_FRACTION", "ADDENDA_FILE",
    "TitleRate", "parse_title_rate", "read_booking_rate",
    "MasEvidence", "MasResolution", "resolve", "agree",
    "rotor_id", "confirmation_key", "describe_key", "format_hz",
    "log_path", "remember", "forget", "lookup", "resolve_for_load",
]

#: default MAS rate when none can be found anywhere (flagged uncertain so the
#: user is warned in the app). 35714 Hz is a common 27Al fast-MAS rate here.
MAS_FALLBACK_HZ = 35714.0

#: a title rate outside this window (Hz) is a unit typo or garbage
PLAUSIBLE_HZ = (1000.0, 150000.0)

#: two rates within this fraction of the larger one agree
AGREE_FRACTION = 0.02

#: the NMRFAM booking sidecar, next to acqus in the EXPNO folder
ADDENDA_FILE = "experiment_addenda.xml"

#: 'MASR 20 kHz', 'MAS 12.5 kHz', 'MASR = 20kHz', 'MASR 35714 Hz'; the unit is
#: mandatory (no unit-less rate exists in either data tree), the prefix too
#: ('1 kHz spikelet spacing' must not read as a rate)
_TITLE_RE = re.compile(r"MASR?\s*[=:]?\s*([\d.]+)\s*(k?Hz)", re.IGNORECASE)

#: 'Rotor RS2427418', 'Rotor SR31602LS - NaCl'; the ID shape rules out
#: 'rotor period', 'Rotor x' and 'Rotor RS' placeholders
_ROTOR_TITLE_RE = re.compile(
    r"\bRotor\s*[:=]?\s*([A-Za-z]{1,2}\d{4,}[A-Za-z]{0,2})\b", re.IGNORECASE)
_ROTOR_FOLDER_RE = re.compile(r"(?:^|_)((?:SR|RS)\d{4,})(?=_|$)", re.IGNORECASE)
_ROTOR_SHAPE_RE = re.compile(r"[A-Za-z]{1,2}\d{4,}[A-Za-z]{0,2}")


def format_hz(value: float) -> str:
    """20000 -> '20 000' (same rendering as ``sidebands.format_hz``; kept
    local so this module stays stdlib-only)."""
    return f"{float(value):,.0f}".replace(",", " ")


def agree(a: float, b: float) -> bool:
    """True when two positive rates are within AGREE_FRACTION of the larger."""
    return abs(a - b) <= AGREE_FRACTION * max(a, b)


# ------------------------------------------------------------------- title
@dataclass
class TitleRate:
    hz: float | None
    static_zero: bool = False
    note: str = ""
    raw: str = ""
    repaired: bool = False


def parse_title_rate(title: str) -> TitleRate:
    """The operator-typed rate: 'MASR 35.714 kHz' -> 35714 Hz. '0 kHz' is a
    static declaration (no candidate). A value outside PLAUSIBLE_HZ is
    re-read in the other unit and marked repaired ('35741 kHz' -> 35741 Hz,
    '35.714 Hz' -> 35714 Hz); implausible both ways -> None with a note."""
    m = _TITLE_RE.search(title or "")
    if not m:
        return TitleRate(None)
    raw = m.group(0).strip()
    try:
        value = float(m.group(1))
    except ValueError:                      # '.' or '1.2.3'
        return TitleRate(None)
    kilo = m.group(2).lower() == "khz"
    hz = value * 1000.0 if kilo else value
    if hz == 0.0:
        return TitleRate(None, static_zero=True, raw=raw)
    lo, hi = PLAUSIBLE_HZ
    if lo <= hz <= hi:
        return TitleRate(hz, raw=raw)
    other = value if kilo else value * 1000.0
    if lo <= other <= hi:
        note = (f'title says "{raw}", read as {other:.0f} Hz (unit typo?)')
        return TitleRate(other, note=note, raw=raw, repaired=True)
    return TitleRate(None, note=f'title says "{raw}" -- implausible, ignored',
                     raw=raw)


# ----------------------------------------------------------------- booking
def read_booking_rate(expno_dir) -> tuple[float | None, bool | None]:
    """(rate_Hz, spinning_flag) from <EXPNO>/experiment_addenda.xml.

    The rate is written in kHz. Returns (None, None) when the file is
    missing or malformed, (None, False) when the booking form's spinning box
    was left unticked (never a rate in that case), (None, True) when the box
    is ticked but the rate is empty or non-numeric."""
    p = Path(expno_dir) / ADDENDA_FILE
    if not p.is_file():
        return None, None
    try:
        root = ET.parse(str(p)).getroot()
    except (ET.ParseError, OSError, ValueError):
        return None, None
    flag_el = root.find("magic_angle_spinning")
    flag: bool | None = None
    if flag_el is not None and flag_el.text is not None:
        t = flag_el.text.strip().lower()
        flag = True if t == "true" else False if t == "false" else None
    if flag is False:
        return None, False
    rate_el = root.find("magic_angle_spinning_rate")
    if rate_el is None or rate_el.text is None:
        return None, flag
    try:
        khz = float(rate_el.text.strip())
    except ValueError:
        return None, flag
    hz = khz * 1000.0
    if not (PLAUSIBLE_HZ[0] <= hz <= PLAUSIBLE_HZ[1]):
        return None, flag
    return hz, flag


# ---------------------------------------------------------------- resolver
@dataclass
class MasEvidence:
    """The three parsed source values (Hz or None) with their annotations."""
    acqus_Hz: float | None = None
    title_Hz: float | None = None
    booking_Hz: float | None = None
    title_note: str = ""
    title_repaired: bool = False
    title_raw: str = ""
    booking_flag: bool | None = None

    def values(self) -> tuple:
        """(acqus, title, booking) rounded to 1 Hz -- the store's match key."""
        return tuple(None if v is None else float(round(v))
                     for v in (self.acqus_Hz, self.title_Hz, self.booking_Hz))

    def as_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "MasEvidence":
        d = d or {}
        return cls(acqus_Hz=d.get("acqus_Hz"), title_Hz=d.get("title_Hz"),
                   booking_Hz=d.get("booking_Hz"),
                   title_note=str(d.get("title_note") or ""),
                   title_repaired=bool(d.get("title_repaired", False)),
                   title_raw=str(d.get("title_raw") or ""),
                   booking_flag=d.get("booking_flag"))


@dataclass
class MasResolution:
    rate_Hz: float
    uncertain: bool
    #: 'acqus' | 'title' | 'booking' | 'acqus+title' | 'acqus+booking' |
    #: 'title+booking' | 'all' | 'static' | 'fallback' | 'highest'
    source: str
    note: str = ""


_NAMES = {"acqus": "acqus MASR", "title": "title", "booking": "booking sidecar"}


def _pick(ev: MasEvidence, names: list[str]) -> float:
    """The value to use inside an agreeing set: a non-repaired title, else
    the booking, else acqus (a repaired title never supplies the number)."""
    if "title" in names and not ev.title_repaired:
        return float(ev.title_Hz)
    if "booking" in names:
        return float(ev.booking_Hz)
    if "acqus" in names:
        return float(ev.acqus_Hz)
    return float(ev.title_Hz)


def resolve(ev: MasEvidence, *, masr_zero: bool = False,
            static_title: bool = False, static_pp: bool = False) -> MasResolution:
    """Majority rule over the present candidates; the static matrix of
    ``bruker._resolve_mas`` is kept exactly (a positive rate always wins over
    static hints but any static hint flags it)."""
    cands = [(n, float(v)) for n, v in (("acqus", ev.acqus_Hz),
                                        ("title", ev.title_Hz),
                                        ("booking", ev.booking_Hz))
             if v is not None and v > 0]
    forced = masr_zero or static_title or static_pp
    # a title rate that was present but unreadable (implausible both ways)
    # is a disagreement the user must see, whatever the other sources say
    title_dropped = ev.title_Hz is None and bool(ev.title_note)

    if not cands:
        if masr_zero:
            return MasResolution(0.0, not (static_title or static_pp), "static",
                                 "acqus MASR = 0" +
                                 ("" if (static_title or static_pp)
                                  else " with no second source"))
        if static_title:
            return MasResolution(0.0, False, "static", "title declares static")
        if static_pp:
            return MasResolution(0.0, True, "static",
                                 "static pulse program, nothing else found")
        return MasResolution(MAS_FALLBACK_HZ, True, "fallback",
                             f"no source found -- {format_hz(MAS_FALLBACK_HZ)} "
                             "Hz assumed")

    names = [n for n, _ in cands]
    if len(cands) == 1:
        n, v = cands[0]
        unc = forced or title_dropped or (n == "title" and ev.title_repaired)
        note = f"{_NAMES[n]} only"
        if n == "title" and ev.title_repaired:
            note += " (repaired unit typo, unconfirmed)"
        return MasResolution(v, unc, n, note)

    if all(agree(a, b) for i, (_, a) in enumerate(cands)
           for _, b in cands[i + 1:]):
        rate = _pick(ev, names)
        src = "all" if len(cands) == 3 else "+".join(names)
        note = ("all sources agree" if len(cands) == 3 else
                f"{_NAMES[names[0]]} and {_NAMES[names[1]]} agree")
        return MasResolution(rate, forced or title_dropped, src, note)

    for i, (na, a) in enumerate(cands):
        for nb, b in cands[i + 1:]:
            if agree(a, b):
                pair = [na, nb]
                rate = _pick(ev, pair)
                (out_name, out_v), = [c for c in cands if c[0] not in pair]
                note = (f"{_NAMES[na]} and {_NAMES[nb]} agree; "
                        f"{_NAMES[out_name]} ({format_hz(out_v)} Hz) is outvoted")
                return MasResolution(rate, forced or title_dropped,
                                     "+".join(pair), note)

    rate = max(v for _, v in cands)
    parts = [f"{_NAMES[n]} {format_hz(v)} Hz" for n, v in cands]
    note = (", ".join(parts[:-1]) + " and " + parts[-1] +
            " disagree; the highest was taken")
    return MasResolution(rate, True, "highest", note)


# ------------------------------------------------------------- rotor / key
def rotor_id(title: str, sample_folder: str) -> str:
    """The rotor actually in the probe: the title's 'Rotor <ID>' line, else
    the SR/RS token of the sample folder, else the folder name itself."""
    m = _ROTOR_TITLE_RE.search(title or "")
    if m:
        return m.group(1).upper()
    m = _ROTOR_FOLDER_RE.search(sample_folder or "")
    if m:
        return m.group(1).upper()
    return sample_folder or ""


def _expno_dir(path) -> Path:
    """Normalise a 1r/fid file or a pdata/N folder to its EXPNO folder."""
    p = Path(path)
    if p.is_file():
        p = p.parent
    if p.parent.name == "pdata":
        p = p.parent.parent
    return p


def confirmation_key(expno_dir, nucleus: str, title: str) -> tuple | None:
    """(session month folder, rotor ID, nucleus) or None when the nucleus is
    empty or the path does not exist."""
    nuc = str(nucleus or "").strip().strip("<>")
    if not nuc or not expno_dir:
        return None
    p = _expno_dir(expno_dir)
    if not p.exists():
        return None
    from larmor.referencing import session_root

    return (str(session_root(p)), rotor_id(title, p.parent.name), nuc)


def describe_key(key) -> str:
    """'2026-05 · rotor RS2427418 · 31P' ('folder <name>' when no rotor ID
    could be read and the sample folder stood in)."""
    session, rotor, nucleus = key
    what = "rotor" if _ROTOR_SHAPE_RE.fullmatch(rotor or "") else "folder"
    return f"{Path(session).name} · {what} {rotor} · {nucleus}"


# ------------------------------------------------------------------- store
def log_path() -> Path:
    """Append-only store: %LOCALAPPDATA%/LARMOR/mas_confirmations.jsonl
    (LARMOR_MAS_LOG overrides). Never truncated by the application."""
    env = os.environ.get("LARMOR_MAS_LOG")
    if env:
        return Path(env)
    base = os.environ.get("LOCALAPPDATA") or str(Path.home())
    return Path(base) / "LARMOR" / "mas_confirmations.jsonl"


def _append(record: dict, path) -> Path:
    from larmor import __version__

    p = Path(path) if path else log_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    record = {"time": _dt.datetime.now().isoformat(timespec="seconds"),
              "larmor": __version__, **record}
    with open(p, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
    return p


def remember(key, ev: MasEvidence, rate_Hz: float, *, note: str = "",
             path=None) -> Path:
    """Record a confirmed rate for this key and this exact source triplet."""
    session, rotor, nucleus = key
    return _append({"action": "remember", "session": str(session),
                    "rotor": rotor, "nucleus": nucleus,
                    "evidence": list(ev.values()), "rate_Hz": float(rate_Hz),
                    "note": note}, path)


def forget(key, ev: MasEvidence, path=None) -> Path:
    """Append a reversal: later lookups of this key + triplet find nothing."""
    session, rotor, nucleus = key
    return _append({"action": "forget", "session": str(session),
                    "rotor": rotor, "nucleus": nucleus,
                    "evidence": list(ev.values())}, path)


def _same_key(rec: dict, key, ev: MasEvidence) -> bool:
    session, rotor, nucleus = key
    if str(rec.get("session", "")).lower() != str(session).lower():
        return False
    if str(rec.get("rotor", "")).lower() != str(rotor).lower():
        return False
    if str(rec.get("nucleus", "")) != str(nucleus):
        return False
    return list(rec.get("evidence") or []) == list(ev.values())


def lookup(key, ev: MasEvidence, path=None) -> dict | None:
    """The latest record for this key and triplet, or None when there is none
    or the latest one is a 'forget'. Corrupt lines are skipped."""
    p = Path(path) if path else log_path()
    if key is None or not p.exists():
        return None
    latest = None
    for line in p.read_text(encoding="utf-8").splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(rec, dict) or not _same_key(rec, key, ev):
            continue
        latest = rec
    if latest is None or latest.get("action") != "remember":
        return None
    if latest.get("rate_Hz") is None:
        return None
    return latest


# -------------------------------------------------------------- the loader
_BLOCK_KEYS = ("acqus_Hz", "title_Hz", "title_note", "title_raw", "booking_Hz",
               "booking_flag", "source", "uncertain", "session", "rotor",
               "nucleus", "confirmed", "confirmed_note")


def resolve_for_load(meta: dict, expno_dir) -> dict:
    """Apply a stored confirmation to a freshly read Bruker ``meta`` and build
    the recipe's ``provenance["mas_rate"]`` block.

    Returns {'spin_rate_Hz', 'mas_uncertain', 'note', 'provenance'}. The store
    is consulted only when the resolution is uncertain; a certain rate is
    never overridden from memory."""
    src = meta.get("mas_sources") or {}
    ev = MasEvidence.from_dict(src)
    spin = float(meta.get("spin_rate_Hz") or meta.get("masr_Hz") or 0.0)
    uncertain = bool(meta.get("mas_uncertain", False))
    nucleus = str(meta.get("nucleus") or "").strip().strip("<>")
    key = confirmation_key(expno_dir, nucleus, meta.get("title", ""))
    block = {
        "acqus_Hz": ev.acqus_Hz, "title_Hz": ev.title_Hz,
        "title_note": ev.title_note, "title_raw": ev.title_raw,
        "booking_Hz": ev.booking_Hz, "booking_flag": ev.booking_flag,
        "source": str(meta.get("mas_source") or ""), "uncertain": uncertain,
        "session": key[0] if key else None,
        "rotor": key[1] if key else (meta.get("rotor") or None),
        "nucleus": nucleus, "confirmed": None, "confirmed_note": "",
    }
    note = ""
    if uncertain and key is not None:
        rec = lookup(key, ev)
        if rec is not None:
            spin = float(rec["rate_Hz"])
            uncertain = False
            when = str(rec.get("time") or "")
            block.update(source="confirmed", uncertain=False, confirmed=when,
                         confirmed_note=str(rec.get("note") or ""))
            trip = " / ".join(
                f"{n} {'-' if v is None else format_hz(v)}"
                for n, v in zip(("acqus", "title", "booking"), ev.values()))
            note = (f"MAS {format_hz(spin)} Hz — confirmed on {when[:10]} for "
                    f"{describe_key(key)} ({trip})")
    assert set(block) == set(_BLOCK_KEYS)
    return {"spin_rate_Hz": spin, "mas_uncertain": uncertain, "note": note,
            "provenance": block}
