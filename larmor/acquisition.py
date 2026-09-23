"""The acquisition record of a Bruker EXPNO and the Experimental paragraph /
Table S1 a paper needs from it.

``read_block`` turns one EXPNO (acqus, ``pdata/<procno>/procs``, the
operator title, ``uxnmr.info``, the NMRFAM booking sidecar) into a FLAT dict
of JSON scalars (``BLOCK_KEYS``) that the loader stores in
``recipe.acquisition``. Everything the block says the recipe already relies
on is read through the same functions, so block and recipe never disagree:
the MAS candidates and their resolution come from ``bruker._meta_1d`` +
``masrate.resolve_for_load`` (the per-session confirmation store included),
the booking rate from ``masrate.read_booking_rate``, the title's ``P1(90)=``
/ tip angle from ``quantitativity.parse_title``, the sample from
``scan.sample_name``, the rotor from ``masrate``'s title / folder rule, B0
from ``uxnmr.info``'s 1H frequency (else BF1 / Xi) over ``nuclei.GAMMA_1H``.

``sentences`` / ``paragraph`` write the Experimental text under hard rules a
student can trust -- nothing is stated that the files do not support:

* the MAS rate is the recipe's confirmed value or a bracketed candidate list
  ``[confirm MAS rate: acqus 4.2 kHz, title 20 kHz, booking 22 kHz]``;
* the flip angle is printed only from title evidence (a stated tip, or
  computed from the title's ``P1(90)=``) and carries its source; otherwise
  the pulse length and the title's own words (``short tip angle``);
* decoupling and 'ambient temperature, no VT gas flow' only from title
  words; a temperature never from acqus TE (0 everywhere);
* referencing to adamantane only with evidence (the recipe's audit block or
  an ``ok`` verdict of the session's referencing audit), else the SR value
  and ``[state the reference standard]``;
* a parameter that varies across a series becomes a range
  ``recycle delays of 12.5–36 s (Table S1)``, never a single value.

``table`` / ``table_csv`` / ``table_latex`` / ``table_markdown`` build Table
S1 with every column whose value differs across the series marked; the
LaTeX form is pure ASCII/TeX so pdflatex compiles it.

Qt-free, read-only (nothing is written into an instrument folder);
``larmor.io.bruker``, ``larmor.referencing`` and ``larmor.masrate`` are
imported inside functions so the module imports without mrsimulator.
"""
from __future__ import annotations

import csv
import datetime as _dt
import io
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

__all__ = [
    "BLOCK_KEYS", "read_block", "parse_title", "flip_angle", "block_for_recipe",
    "folder_expnos",
    "summary_lines", "OP_PHRASES", "describe_processing", "sentences",
    "AcqTable", "ALWAYS_VARIES", "COLUMNS", "table", "table_csv", "table_latex",
    "table_markdown", "paragraph", "referencing_evidence", "mas_candidates",
]

#: the flat record (every value a JSON scalar or None)
BLOCK_KEYS: tuple[str, ...] = (
    "expno_path", "expno", "procno", "data_file", "sample", "sample_folder",
    "title", "nucleus", "bf1_MHz", "sfo1_MHz", "sf_MHz", "sr_hz",
    "magnet_1h_MHz", "b0_T", "spectrometer", "topspin", "probe", "pulprog",
    "ns", "ds", "rg", "d1_s", "aq_s", "p1_us", "plw1_W", "p90_us", "flip_deg",
    "flip_source", "tip_words", "rotor", "vt_note", "decoupling_note", "sw_Hz",
    "td", "date_utc", "mas_acqus_Hz", "mas_title_Hz", "mas_booking_Hz",
    "spin_rate_Hz", "mas_uncertain", "mas_source", "wdw", "lb_hz", "gb", "ssb",
    "si", "tdeff", "ph_mod", "phc0", "phc1", "absg",
)

#: TopSpin codes (same tables as larmor.comparability)
_WDW = {0: "none", 1: "EM", 2: "GM", 3: "SINE", 4: "QSINE", 5: "TRAP", 6: "USER",
        7: "SINC", 8: "QSINC", 9: "TRAF", 10: "TRAFS"}
_PH_MOD = {0: "no", 1: "pk", 2: "mc", 3: "ps"}

#: title conventions (verified over the NMRFAM tree: 'with 19F decoupling' x5,
#: 'without 19F decoupling' x5, 'VT ambient - no flow' x672, 'VT at ambient,
#: no flow' x13, 'short tip angle' x16)
_DECOUPLING_RE = re.compile(
    r"\b(with|without)\s+(\d{1,3}\s*[A-Za-z]{1,2})\s+(?:[A-Za-z]+\s+)?decoupling\b",
    re.IGNORECASE)
_VT_AMBIENT_RE = re.compile(r"\bVT\b[^\n]{0,12}?\bambient\b[^\n]{0,20}?\bno\s*flow\b",
                            re.IGNORECASE)
_TIP_WORDS_RE = re.compile(r"\bshort\s+tip\s+angle\b", re.IGNORECASE)
_TOPSPIN_RE = re.compile(r"TopSpin\s+([\d.]+[A-Za-z0-9.]*)", re.IGNORECASE)

# ---------------------------------------------------------------- readers


def _read_text(p: Path) -> str:
    try:
        return p.read_text(encoding="utf-8", errors="replace")
    except (OSError, ValueError):
        return ""


@lru_cache(maxsize=256)
def _uxnmr_info_cached(path_str: str, _mtime_ns: int) -> dict:
    d: dict = {}
    for line in _read_text(Path(path_str)).splitlines()[:40]:
        if ":" not in line:
            continue
        key, _, val = line.partition(":")
        key, val = key.strip(), val.strip()
        if key and val and key not in d:
            d[key] = val
    return d


def _uxnmr_info(expno) -> dict:
    """The header block of ``<EXPNO>/uxnmr.info`` as {key: value} ('System',
    '1H-frequency', 'Release', ...); {} when absent. Cached by path + mtime;
    latin-1 bytes are tolerated."""
    p = Path(expno) / "uxnmr.info"
    try:
        st = p.stat()
    except OSError:
        return {}
    return _uxnmr_info_cached(str(p), st.st_mtime_ns)


def _topspin_version(expno: Path, info: dict) -> str:
    """'3.6.2' from the acqus header line (``##TITLE= Parameter file, TopSpin
    3.6.2``), else from uxnmr.info's Release line, else ''."""
    try:
        with open(expno / "acqus", "r", encoding="utf-8", errors="replace") as f:
            head = f.readline() + f.readline()
    except OSError:
        head = ""
    m = _TOPSPIN_RE.search(head) or _TOPSPIN_RE.search(info.get("Release", ""))
    return m.group(1).rstrip(".") if m else ""


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if x == x else None


def _i(v):
    x = _f(v)
    return None if x is None else int(x)


def _clean(v) -> str:
    return str(v).strip().strip("<>").strip() if v is not None else ""


def _first_line(title: str) -> str:
    for ln in (title or "").splitlines():
        if ln.strip():
            return ln.strip()
    return ""


def parse_title(title: str) -> dict:
    """The operator's claims in a TopSpin title: ``p90_us`` and
    ``flip_title_deg`` (``quantitativity.parse_title``: 'P1(90)=3.750; 30 deg
    tip', '11 degree tip angle', '30 degree pulse' -- a stated angle above
    90 is not a flip angle), ``tip_words`` ('short tip angle' | ''),
    ``rotor`` (the 'Rotor <ID>' token), ``vt_note`` ('ambient temperature,
    no VT gas flow' only from 'VT ambient - no flow' / 'VT at ambient, no
    flow'), ``decoupling_note`` ('with 19F decoupling' | 'without 19F
    decoupling' | ''), ``mas_title_Hz`` (``masrate.parse_title_rate``)."""
    from larmor import masrate, quantitativity

    text = title or ""
    q = quantitativity.parse_title(text)
    flip = q.get("flip_deg")
    if flip is not None and not (0.0 < flip <= 90.0):
        flip = None
    m = masrate._ROTOR_TITLE_RE.search(text)
    rotor = m.group(1).upper() if m else ""
    m = _DECOUPLING_RE.search(text)
    dec = (f"{m.group(1).lower()} {m.group(2).replace(' ', '')} decoupling"
           if m else "")
    return {"p90_us": q.get("p90_us"), "flip_title_deg": flip,
            "tip_words": "short tip angle" if _TIP_WORDS_RE.search(text) else "",
            "rotor": rotor,
            "vt_note": ("ambient temperature, no VT gas flow"
                        if _VT_AMBIENT_RE.search(text) else ""),
            "decoupling_note": dec,
            "mas_title_Hz": masrate.parse_title_rate(text).hz}


def flip_angle(p1_us, p90_us, title: str) -> tuple[float | None, str]:
    """(flip angle in degrees, source): a tip angle stated in the title
    ('30 deg tip', '11 degree tip angle', '30 degree pulse'; <= 90) wins as
    ``'title'``; else 90 * P1 / P1(90) as ``'computed'`` when the 90-degree
    pulse is known; else ``(None, '')`` -- the paragraph then prints the pulse
    length with the title's own words and no degrees."""
    claims = parse_title(title)
    if claims["flip_title_deg"] is not None:
        return float(claims["flip_title_deg"]), "title"
    p90 = _f(p90_us) if p90_us is not None else claims["p90_us"]
    p1 = _f(p1_us)
    if p1 is not None and p90 and p90 > 0 and p1 > 0:
        return 90.0 * p1 / p90, "computed"
    return None, ""


def _rotor(title: str, sample_folder: str) -> str:
    """The rotor in the probe from the title's 'Rotor <ID>' line or the
    folder's SR/RS token (``masrate.rotor_id``); '' when neither carries an
    ID (rotor_id then falls back to the folder name, which is not a rotor)."""
    from larmor import masrate

    r = masrate.rotor_id(title or "", sample_folder or "")
    return "" if not r or r == (sample_folder or "") else r


def _date_utc(epoch) -> str:
    t = _f(epoch)
    if not t:
        return ""
    try:
        return _dt.datetime.fromtimestamp(t, _dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    except (OverflowError, OSError, ValueError):
        return ""


def read_block(path, procno: int | None = None) -> dict:
    """The acquisition + processing record of one EXPNO as a flat dict
    (``BLOCK_KEYS``). ``path`` is anything ``bruker.resolve`` accepts (an
    EXPNO, a pdata/N folder, a 1r or fid file) or a bare EXPNO holding only
    acqus; ``procno`` defaults to the path's own pdata/N, else 1. {} when
    there is no readable acqus. Never raises; nothing is written."""
    try:
        return _read_block(path, procno)
    except Exception:                                    # noqa: BLE001
        return {}


def _read_block(path, procno) -> dict:
    from larmor import masrate
    from larmor.io import bruker, scan
    from larmor.referencing import _jcamp, xi_ratio

    p = Path(str(path))
    try:
        ref = bruker.resolve(p)
        expno, ref_procno = Path(ref.expno), int(ref.procno or 1)
    except (ValueError, FileNotFoundError, OSError, NotADirectoryError):
        if not bruker.is_expno(p):
            return {}
        expno, ref_procno = p, 1
    acqus = _jcamp(expno / "acqus")
    if not acqus:
        return {}
    if procno is None:
        procno = ref_procno
    procno = int(procno)
    pdata = expno / "pdata" / str(procno)
    title = bruker._read_title(pdata) or bruker._read_title(expno / "pdata" / "1")
    # the recipe's own MAS resolution: the three sources, the majority rule
    # and the per-session confirmation (never a second rule here)
    meta = bruker._meta_1d(acqus, title, expno)
    mas = masrate.resolve_for_load(meta, expno)
    src = mas["provenance"]
    procs = _jcamp(pdata / "procs") or {}
    info = _uxnmr_info(expno)
    name = scan.sample_name(expno, title)
    claims = parse_title(title)

    nucleus = _clean(acqus.get("NUC1"))
    bf1 = _f(acqus.get("BF1"))
    sfo1 = _f(acqus.get("SFO1"))
    sf = _f(procs.get("SF")) if procs else None
    if sf is not None and sf <= 0:
        sf = None
    sr = (sf - bf1) * 1e6 if (sf is not None and bf1) else None

    magnet = None
    m = re.search(r"([\d.]+)\s*MHz", info.get("1H-frequency", ""))
    if m:
        magnet = _f(m.group(1))
    if magnet is None and bf1 and nucleus:
        try:
            magnet = bf1 / xi_ratio(nucleus)
        except Exception:                                # noqa: BLE001
            magnet = None
    from larmor.nuclei import GAMMA_1H
    b0 = magnet / GAMMA_1H if magnet else None

    system = info.get("System", "")
    spectrometer = re.sub(r"\s*NMR spectrometer\s*$", "", system).strip()
    if not spectrometer:
        instrum = _clean(acqus.get("INSTRUM"))
        spectrometer = "" if instrum.lower() in ("", "spect") else instrum

    p_us = bruker._num_list(acqus.get("P"))
    plw = bruker._num_list(acqus.get("PLW"))
    p1 = p_us[1] if len(p_us) > 1 else None
    plw1 = plw[1] if len(plw) > 1 else None
    flip, flip_src = flip_angle(p1, claims["p90_us"], title)
    td = _i(acqus.get("TD"))
    sw = _f(acqus.get("SW_h"))

    data_file = ""
    if (pdata / "1r").is_file():
        data_file = f"pdata/{procno}/1r"
    elif (expno / "fid").is_file():
        data_file = "fid"

    block = {
        "expno_path": str(expno), "expno": expno.name, "procno": procno,
        "data_file": data_file, "sample": name.key, "sample_folder": name.folder,
        "title": _first_line(title), "nucleus": nucleus, "bf1_MHz": bf1,
        "sfo1_MHz": sfo1, "sf_MHz": sf, "sr_hz": sr, "magnet_1h_MHz": magnet,
        "b0_T": b0, "spectrometer": spectrometer,
        "topspin": _topspin_version(expno, info), "probe": _clean(acqus.get("PROBHD")),
        "pulprog": _clean(acqus.get("PULPROG")), "ns": _i(acqus.get("NS")),
        "ds": _i(acqus.get("DS")), "rg": _f(acqus.get("RG")),
        "d1_s": meta.get("d1_s"), "aq_s": meta.get("aq_s"), "p1_us": p1,
        "plw1_W": plw1, "p90_us": claims["p90_us"], "flip_deg": flip,
        "flip_source": flip_src, "tip_words": claims["tip_words"],
        "rotor": _rotor(title, expno.parent.name), "vt_note": claims["vt_note"],
        "decoupling_note": claims["decoupling_note"], "sw_Hz": sw, "td": td,
        "date_utc": _date_utc(acqus.get("DATE")),
        "mas_acqus_Hz": src.get("acqus_Hz"), "mas_title_Hz": src.get("title_Hz"),
        "mas_booking_Hz": src.get("booking_Hz"),
        "spin_rate_Hz": float(mas["spin_rate_Hz"]),
        "mas_uncertain": bool(mas["mas_uncertain"]), "mas_source": str(src.get("source") or ""),
        "wdw": _WDW.get(_i(procs.get("WDW")) or 0, str(procs.get("WDW"))) if procs else "",
        "lb_hz": _f(procs.get("LB")) if procs else None,
        "gb": _f(procs.get("GB")) if procs else None,
        "ssb": _f(procs.get("SSB")) if procs else None,
        "si": _i(procs.get("SI")) if procs else None,
        "tdeff": _i(procs.get("TDeff")) if procs else None,
        "ph_mod": _PH_MOD.get(_i(procs.get("PH_mod")) or 0, str(procs.get("PH_mod"))) if procs else "",
        "phc0": _f(procs.get("PHC0")) if procs else None,
        "phc1": _f(procs.get("PHC1")) if procs else None,
        "absg": _i(procs.get("ABSG")) if procs else None,
    }
    assert set(block) == set(BLOCK_KEYS)
    return block


def folder_expnos(folder) -> list[str]:
    """The EXPNO folders a dropped folder stands for: a sample folder -> its
    own EXPNO children (numeric order); a month / session folder -> every
    EXPNO ``referencing.scan_session`` finds; an EXPNO or pdata folder ->
    itself. [] when nothing holds an acqus."""
    from larmor import referencing
    from larmor.io import bruker

    p = Path(str(folder))
    if not p.is_dir():
        return []
    if bruker.is_expno(p) or p.parent.name == "pdata":
        return [str(p)]
    kids = [c for c in p.iterdir() if c.is_dir() and (c / "acqus").exists()]
    if kids:
        kids.sort(key=lambda c: (not c.name.isdigit(), int(c.name) if c.name.isdigit() else 0,
                                 c.name))
        return [str(c) for c in kids]
    return [a.path for a in referencing.scan_session(referencing.session_root(p))]


def block_for_recipe(recipe_dict: dict) -> dict:
    """The stored ``acquisition`` block of a recipe dict; else, for a Bruker
    source that still exists, ``read_block(source_path)``; {} for CSV / dmfit
    / Varian sources or when nothing is readable."""
    rec = recipe_dict or {}
    stored = rec.get("acquisition")
    if stored:
        return dict(stored)
    kind = str(rec.get("source_kind") or "")
    src = str(rec.get("source_path") or "")
    if not src or kind in ("csv", "fxmla", "varian"):
        return {}
    return read_block(src)


# ------------------------------------------------------------- MAS state
def mas_candidates(block_or_mas: dict) -> list[tuple[str, float]]:
    """[('acqus', 4200.0), ('title', 20000.0), ('booking', 22000.0)] from an
    acquisition block (``mas_*_Hz``) or a recipe ``provenance['mas_rate']``
    block (``acqus_Hz`` ...); only positive rates."""
    d = block_or_mas or {}
    out = []
    for name in ("acqus", "title", "booking"):
        v = d.get(f"mas_{name}_Hz", d.get(f"{name}_Hz"))
        v = _f(v)
        if v is not None and v > 0:
            out.append((name, v))
    return out


def _khz(v: float) -> str:
    """4200 -> '4.2 kHz', 20000 -> '20 kHz', 35714 -> '35.7 kHz'."""
    k = float(v) / 1000.0
    s = f"{k:.3g}" if k < 100 else f"{k:.0f}"
    return f"{s} kHz"


def _mas_state(block: dict, spin_rate_Hz, mas_uncertain, mas_rate: dict | None):
    """(rate_Hz or Range, uncertain, candidates): the recipe's confirmed state
    when given, else the block's own resolution; the candidate list comes
    from the recipe's ``provenance['mas_rate']`` when given, else the
    block."""
    rate = spin_rate_Hz if spin_rate_Hz is not None else block.get("spin_rate_Hz")
    unc = mas_uncertain if mas_uncertain is not None else bool(block.get("mas_uncertain"))
    if mas_rate:
        if mas_rate.get("confirmed"):
            unc = False
        cands = mas_candidates(mas_rate)
    else:
        cands = mas_candidates(block)
    return rate, bool(unc), cands


# --------------------------------------------------------- summary lines
def _g(v, unit: str = "") -> str:
    if v is None or v == "":
        return "—"
    if isinstance(v, _Range):
        return str(v) + (f" {unit}" if unit else "")
    if isinstance(v, float):
        s = format(v, "g")
    else:
        s = str(v)
    return s + (f" {unit}" if unit else "")


def summary_lines(block: dict, software: dict | None = None, spin_rate_Hz=None,
                  mas_uncertain=None) -> list[str]:
    """The compact read-only summary of an acquisition block for the
    Experiment dialog, the experiment strip's tooltip and ``larmor info``:
    acqus facts, probe / spectrometer / date, the MAS sources and verdict,
    the procs processing keys, and -- when a software block is given -- the
    fitting versions. [] for an empty block."""
    b = block or {}
    if not b:
        return []
    out = []
    flip = ""
    if b.get("flip_deg") is not None:
        how = {"title": "from the title", "computed": "from the title's P1(90)"}.get(
            b.get("flip_source"), "")
        flip = f" ({b['flip_deg']:.0f}°{', ' + how if how else ''})"
    elif b.get("tip_words"):
        flip = f" ({b['tip_words']})"
    acq = [f"acqus: {b.get('pulprog') or '?'}", f"NS {_g(b.get('ns'))}"]
    if b.get("ds"):
        acq.append(f"DS {b['ds']}")
    acq += [f"D1 {_g(b.get('d1_s'), 's')}"]
    if b.get("aq_s") is not None:
        acq.append(f"AQ {b['aq_s'] * 1000:.3g} ms")
    acq += [f"P1 {_g(b.get('p1_us'), 'µs')}{flip}", f"PLW1 {_g(b.get('plw1_W'), 'W')}"]
    if b.get("sw_Hz"):
        acq.append(f"SW {b['sw_Hz'] / 1000:g} kHz / TD {_g(b.get('td'))}")
    if b.get("decoupling_note"):
        acq.append(b["decoupling_note"])
    out.append(" · ".join(acq))

    inst = []
    if b.get("probe"):
        inst.append(f"probe {b['probe']}")
    if b.get("rotor"):
        inst.append(f"rotor {b['rotor']}")
    field_txt = ""
    if b.get("magnet_1h_MHz"):
        field_txt = f"{b['magnet_1h_MHz']:.2f} MHz"
        if b.get("b0_T"):
            field_txt += f", {b['b0_T']:.1f} T"
    if b.get("spectrometer"):
        inst.append(b["spectrometer"] + (f" ({field_txt})" if field_txt else ""))
    elif field_txt:
        inst.append(field_txt)
    if b.get("vt_note"):
        inst.append(b["vt_note"])
    if b.get("date_utc"):
        inst.append(b["date_utc"].replace("T", " ").replace("Z", " UTC"))
    if inst:
        out.append(" · ".join(inst))

    rate, unc, cands = _mas_state(b, spin_rate_Hz, mas_uncertain, None)
    if cands or rate:
        parts = [f"{n} {v:.0f} Hz" if n == "acqus" else f"{n} {_khz(v)}" for n, v in cands]
        verdict = ("static" if not rate else f"using {float(rate):.0f} Hz")
        if unc:
            verdict += " (confirm)"
        out.append("MAS: " + (" · ".join(parts) + " → " if parts else "") + verdict)

    if b.get("wdw") or b.get("si"):
        proc = [f"procs (pdata/{b.get('procno', 1)}): {b.get('wdw') or '?'}"]
        if b.get("lb_hz") is not None:
            proc[0] += f" LB {b['lb_hz']:g} Hz"
        if b.get("wdw") == "GM" and b.get("gb") is not None:
            proc.append(f"GB {b['gb']:g}")
        if b.get("wdw") in ("SINE", "QSINE") and b.get("ssb") is not None:
            proc.append(f"SSB {b['ssb']:g}")
        if b.get("tdeff"):
            proc.append(f"TDeff {b['tdeff']}")
        if b.get("si"):
            proc.append(f"SI {b['si']}")
        if b.get("ph_mod"):
            proc.append(b["ph_mod"])
        if b.get("absg") is not None:
            proc.append(f"ABSG {b['absg']}")
        if b.get("sf_MHz"):
            sr = b.get("sr_hz")
            proc.append(f"SF {b['sf_MHz']:.8f} MHz"
                        + (f" (SR {sr:.2f} Hz)" if sr is not None else ""))
        out.append(" · ".join(proc))

    if software:
        sw = software
        v = f"fitted with LARMOR {sw.get('larmor') or '?'}"
        if sw.get("git_commit"):
            v += f" (commit {str(sw['git_commit'])[:7]})"
        for name in ("mrsimulator", "lmfit", "numpy"):
            if sw.get(name):
                v += f" · {name} {sw[name]}"
        if sw.get("fitted"):
            v += f" · {str(sw['fitted']).replace('T', ' ')}"
        out.append(v)
    return out


# ------------------------------------------------------- LARMOR processing
#: op -> phrase; ``{name}`` placeholders are the op's own parameters
OP_PHRASES: dict[str, str] = {
    "em": "exponential apodization (LB {lb_hz:g} Hz)",
    "gm": "Gaussian apodization (LB {lb_hz:g} Hz, GB {gb:g})",
    "sine": "sine-bell apodization (SSB {ssb:g})",
    "traf": "Traficante apodization (LB {lb_hz:g} Hz)",
    "tdeff": "truncation of the fid to {points} points",
    "zf": "zero-filling",
    "ft": "Fourier transformation",
    "phase": "phase correction (p0 {p0:g}°, p1 {p1:g}°)",
    "autophase": "automatic phase correction",
    "baseline": "polynomial baseline (order {order})",
    "iterbaseline": "iterative polynomial baseline",
    "flat_baseline": "flat baseline (edge offset)",
    "twopoint_bg": "two-point linear baseline",
    "subtract_avg": "constant-offset subtraction",
    "sr": "re-referencing (SR {sr_hz:g} Hz)",
    "extract": "extraction of {hi_ppm:g} to {lo_ppm:g} ppm",
    "wurst_correct": "WURST excitation-profile correction",
    "hilbert": "Hilbert transform (imaginary part restored)",
    "ift": "inverse Fourier transformation",
    "magnitude": "magnitude calculation",
    "scale": "scaling (× {factor:g})",
    "normalize": "normalisation",
    "fcor": "first-point correction (× {factor:g})",
    "lp": "linear prediction",
    "swap_echo": "echo top moved to the first point",
    "echo_apodize": "symmetric echo apodization (LB {lb_hz:g} Hz)",
    "shift_fid": "fid shift by {points} points",
    "offset": "constant offset",
    "scale_sw": "spectral-width rescaling",
    "real": "real part", "imag": "imaginary part", "conj": "complex conjugate",
}


def _phrase(op: dict) -> str:
    name = str(op.get("op", ""))
    tmpl = OP_PHRASES.get(name)
    if tmpl is not None:
        try:
            return tmpl.format(**{k: v for k, v in op.items() if k != "op"})
        except (KeyError, ValueError, TypeError, IndexError):
            pass
    kw = ", ".join(f"{k}={v}" for k, v in op.items() if k != "op")
    return name + (f" ({kw})" if kw else "")


def describe_processing(ops) -> str:
    """'exponential apodization (LB 100 Hz), two-point linear baseline and
    re-referencing (SR 120 Hz)' for a recipe's processing list; an unknown
    op prints as ``op (k=v, …)`` exactly as *Process > Processing steps*
    lists it. '' for no ops."""
    phrases = [_phrase(o) for o in (ops or []) if isinstance(o, dict) and o.get("op")]
    if not phrases:
        return ""
    if len(phrases) == 1:
        return phrases[0]
    return ", ".join(phrases[:-1]) + " and " + phrases[-1]


# --------------------------------------------------------------- ranges
class _Range:
    """A parameter that varies across a series: rendered 'lo–hi' for numbers,
    'a / b' for strings, always followed by '(Table S1)' in prose."""

    def __init__(self, values):
        self.values = [v for v in values if v is not None and v != ""]
        nums = [v for v in self.values if isinstance(v, (int, float)) and not isinstance(v, bool)]
        self.numeric = bool(nums) and len(nums) == len(self.values)
        self.lo = min(nums) if self.numeric else None
        self.hi = max(nums) if self.numeric else None

    def __str__(self) -> str:
        if self.numeric:
            return f"{format(float(self.lo), 'g')}–{format(float(self.hi), 'g')}"
        seen = []
        for v in self.values:
            if v not in seen:
                seen.append(v)
        return " / ".join(str(v) for v in seen)

    def fmt(self, f) -> str:
        """lo–hi through a per-value formatter (kHz, 1 decimal ...)."""
        if self.numeric:
            return f"{f(self.lo)}–{f(self.hi)}"
        return str(self)


def _is_range(v) -> bool:
    return isinstance(v, _Range)


def _val(v, f=None, unit: str = "") -> str:
    """One value or a range in prose; a range carries '(Table S1)'."""
    if _is_range(v):
        s = v.fmt(f) if f else str(v)
        return f"{s}{' ' + unit if unit else ''} (Table S1)"
    s = f(v) if f else (format(float(v), "g") if isinstance(v, float) else str(v))
    return f"{s}{' ' + unit if unit else ''}"


def _plural(v) -> bool:
    return _is_range(v)


# ------------------------------------------------------------- sentences
def _fmt_sr(sr) -> str:
    return f"{float(sr):.2f}".replace("-", "−")


def sentences(block: dict, *, spin_rate_Hz=None, mas_uncertain=None,
              mas_rate: dict | None = None, sr_hz=None,
              referencing_text: str | None = None, larmor_ops=None,
              n_spectra: int = 1) -> list[str]:
    """The acquisition / processing sentences of the Experimental paragraph
    for one block (or a merged series view whose varying values are ranges).
    ``spin_rate_Hz`` / ``mas_uncertain`` are the recipe's confirmed state
    (the block's own resolution otherwise); ``mas_rate`` is the recipe's
    ``provenance['mas_rate']`` (its three candidates feed the bracket);
    ``sr_hz`` overrides the block's SR; ``referencing_text`` is the evidence
    sentence fragment (``referencing_evidence``) or None for a bracket;
    ``larmor_ops`` is the recipe's processing list."""
    b = dict(block or {})
    if not b:
        return []
    many = n_spectra > 1
    nuc = b.get("nucleus") or "the"
    were, spectra = ("were", "spectra") if many else ("was", "spectrum")
    rate, unc, cands = _mas_state(b, spin_rate_Hz, mas_uncertain, mas_rate)
    static = (not _is_range(rate)) and (rate is None or float(rate or 0) == 0.0) and not unc
    out: list[str] = []

    # 1 -- instrument and field
    nu = b.get("sf_MHz") or b.get("sfo1_MHz")
    field_bits = []
    if b.get("b0_T"):
        field_bits.append(f"B₀ = {float(b['b0_T']):.1f} T")
    if b.get("magnet_1h_MHz"):
        field_bits.append(f"{float(b['magnet_1h_MHz']):.2f} MHz for 1H")
    if nu:
        field_bits.append(f"{float(nu):.2f} MHz for {nuc}")
    field_txt = f" ({'; '.join(field_bits)})" if field_bits else ""
    kind = "static" if static else "MAS"
    s = f"The {nuc} {kind} NMR {spectra} {were} acquired"
    if b.get("spectrometer"):
        s += f" on a Bruker {_val(b['spectrometer'])} spectrometer{field_txt}"
    elif field_bits:
        s += " on a spectrometer operating at " + field_bits[0]
        if field_bits[1:]:
            s += f" ({'; '.join(field_bits[1:])})"
    if b.get("probe"):
        s += f" with a {_val(b['probe'])} probe"
        if b.get("rotor"):
            s += f" (rotor{'s' if _plural(b['rotor']) else ''} {_val(b['rotor'])})"
    elif b.get("rotor"):
        s += f" (rotor{'s' if _plural(b['rotor']) else ''} {_val(b['rotor'])})"
    out.append(s + ".")

    # 2 -- MAS
    if static:
        out.append(f"The {spectra} {were} acquired static (no sample spinning).")
    elif unc:
        if cands:
            lst = ", ".join(f"{n} {_khz(v)}" for n, v in cands)
            pp = b.get("pulprog") or ""
            hint = (f"; static pulse program {pp}" if (
                pp and any(k in str(pp).lower() for k in ("wcpmg", "wurst", "static"))) else "")
            out.append(f"[confirm MAS rate: {lst}{hint}]")
        else:
            out.append("[confirm MAS rate: no source found in acqus, title or booking sidecar]")
    else:
        out.append("The magic-angle spinning rate was "
                   f"{_val(rate, lambda v: f'{float(v) / 1000:.1f}', 'kHz')}.")

    # 3 -- pulse, flip, power, recycle, transients
    from larmor.io import scan

    pp = b.get("pulprog")
    pp_first = pp.values[0] if _is_range(pp) else pp
    kind_txt = scan._classify(str(pp_first or ""))
    if pp:
        if kind_txt and kind_txt not in ("unknown", str(pp_first)):
            # 'Single pulse' -> 'a single pulse acquisition'; an acronym
            # ('CPMG (T2)', 'MQMAS') keeps its case
            low = (kind_txt[0].lower() + kind_txt[1:]
                   if len(kind_txt) > 1 and kind_txt[1].islower() else kind_txt)
            art = "An" if low[:1].lower() in "aeiou" else "A"
            head = f"{kind_txt} acquisitions" if many else f"{art} {low} acquisition"
            head += f" (pulse program {_val(pp)})"
        else:
            head = f"The pulse program {_val(pp)}"
    else:
        head = f"The {spectra}"
    parts = []
    if b.get("p1_us") is not None:
        pulse = f"a {_val(b['p1_us'], lambda v: format(float(v), 'g'), 'µs')} pulse"
        flip, src = b.get("flip_deg"), b.get("flip_source")
        p90 = b.get("p90_us")
        p90_txt = format(float(p90), "g") + " µs" if p90 and not _is_range(p90) else ""
        if flip is not None and not _is_range(flip):
            if src == "title":
                pulse += (f" ({float(flip):.0f}° flip angle, as stated in the acquisition title"
                          + (f"; P1(90) = {p90_txt}" if p90_txt else "") + ")")
            elif src == "computed":
                pulse += (f" (≈{float(flip):.1f}° flip angle, from the title's 90° pulse"
                          + (f" of {p90_txt}" if p90_txt else "") + ")")
        elif flip is not None:
            pulse += f" ({_val(flip, lambda v: f'{float(v):.0f}', '°')} flip angle)"
        elif b.get("tip_words"):
            pulse += f" ({_val(b['tip_words'])})"
        if b.get("plw1_W") is not None:
            pulse += f" at {_val(b['plw1_W'], lambda v: format(float(v), 'g'), 'W')}"
        parts.append(pulse)
    if b.get("decoupling_note"):
        parts.append(str(_val(b["decoupling_note"])))
    if b.get("d1_s") is not None:
        d1 = b["d1_s"]
        parts.append(("recycle delays of " if _plural(d1) else "a recycle delay of ")
                     + _val(d1, lambda v: format(float(v), "g"), "s"))
    if b.get("ns") is not None:
        ns = b["ns"]
        t = f"{_val(ns, lambda v: str(int(v)))} transients"
        ds = b.get("ds")
        if ds and not _is_range(ds):
            t += f" ({int(ds)} dummy scan{'s' if int(ds) != 1 else ''})"
        parts.append(t)
    if parts:
        joined = parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + " and " + parts[-1]
        out.append(f"{head} used {joined}.")
    elif pp:
        out.append(head + ".")

    # 4 -- spectral width, points, temperature words
    if b.get("sw_Hz"):
        s = f"The spectral width was {_val(b['sw_Hz'], lambda v: format(float(v) / 1000, 'g'), 'kHz')}"
        if b.get("td"):
            s += f" (TD = {_val(b['td'], lambda v: str(int(v)))} points)"
        if b.get("vt_note"):
            s += ", at ambient temperature without VT gas flow"
        out.append(s + ".")
    elif b.get("vt_note"):
        out.append(f"The {spectra} {were} acquired at ambient temperature without VT gas flow.")

    # 5 -- referencing (evidence or a bracket)
    sr = sr_hz if sr_hz is not None else b.get("sr_hz")
    sr_txt = ""
    if sr is not None and not _is_range(sr):
        sr_txt = f"SR = {_fmt_sr(sr)} Hz"
    elif _is_range(sr):
        sr_txt = f"SR = {sr.fmt(_fmt_sr)} Hz (Table S1)"
    if referencing_text:
        out.append(f"Chemical shifts were {referencing_text}"
                   + (f" ({sr_txt})" if sr_txt else "") + ".")
    elif sr_txt and not (not _is_range(sr) and abs(float(sr)) < 0.5):
        out.append(f"The spectral reference was {sr_txt} [state the reference standard].")
    else:
        out.append("[state the reference standard; the processed data carry no SR]")

    # 6 -- TopSpin processing, then LARMOR's own steps
    proc = []
    wdw = b.get("wdw")
    wdw_first = wdw.values[0] if _is_range(wdw) else wdw
    if wdw_first:
        lb = b.get("lb_hz")
        lb_txt = f" (LB = {_val(lb, lambda v: format(float(v), 'g'), 'Hz')})" if lb is not None else ""
        if _is_range(wdw):
            proc.append(f"{_val(wdw)} apodization{lb_txt}")
        elif wdw == "EM":
            proc.append(f"exponential apodization{lb_txt}")
        elif wdw == "GM":
            gb = b.get("gb")
            proc.append("Gaussian apodization" + (
                f" (LB = {_val(lb, lambda v: format(float(v), 'g'), 'Hz')}, GB = {_val(gb)})"
                if lb is not None and gb is not None else lb_txt))
        elif wdw in ("SINE", "QSINE"):
            ssb = b.get("ssb")
            proc.append(("a squared " if wdw == "QSINE" else "a ") + "sine-bell window"
                        + (f" (SSB = {_val(ssb)})" if ssb is not None else ""))
        elif wdw == "none":
            proc.append("no apodization")
        else:
            proc.append(f"a {wdw} window{lb_txt}")
    tdeff, td = b.get("tdeff"), b.get("td")
    if tdeff and (not td or _is_range(tdeff) or _is_range(td) or int(tdeff) < int(td)):
        proc.append(f"the first {_val(tdeff, lambda v: str(int(v)))} points of the fid (TDeff)")
    if b.get("si"):
        proc.append(f"zero-filling to {_val(b['si'], lambda v: str(int(v)))} points (SI)")
    ph = b.get("ph_mod")
    ph_first = ph.values[0] if _is_range(ph) else ph
    if ph_first:
        ph_txt = {"pk": "manual phase correction", "mc": "magnitude calculation",
                  "ps": "power-spectrum calculation", "no": "no phase correction"}
        proc.append(_val(ph) if _is_range(ph) else ph_txt.get(str(ph), f"phase mode {ph}"))
    if proc:
        ts = f" {b['topspin']}" if b.get("topspin") and not _is_range(b["topspin"]) else ""
        s = f"Spectra were processed in TopSpin{ts} with " + (
            proc[0] if len(proc) == 1 else ", ".join(proc[:-1]) + " and " + proc[-1])
        lar = describe_processing(larmor_ops)
        if lar:
            s += f"; in LARMOR, {lar} preceded the fit"
        out.append(s + ".")
    else:
        lar = describe_processing(larmor_ops)
        if lar:
            out.append(f"In LARMOR, {lar} preceded the fit.")
    return out


# ---------------------------------------------------------------- table
@dataclass
class AcqTable:
    """Table S1: ``columns`` (key, header, unit, fmt, tex_header), one row
    dict per spectrum (every BLOCK_KEY + the confirmed ``spin_rate_Hz`` /
    ``mas_uncertain``), the keys that vary across the rows and the constant
    values of the rest."""

    columns: list[tuple] = field(default_factory=list)
    rows: list[dict] = field(default_factory=list)
    varying: set = field(default_factory=set)
    constants: dict = field(default_factory=dict)

    def header(self, key: str) -> str:
        for k, h, *_ in self.columns:
            if k == key:
                return h
        return key

    def varying_headers(self) -> list[str]:
        return [h for k, h, *_ in self.columns if k in self.varying]


#: keys that legitimately differ per spectrum and never count as variation
ALWAYS_VARIES: frozenset = frozenset({
    "sample", "sample_folder", "expno", "procno", "date_utc", "expno_path",
    "data_file", "title", "mas_uncertain", "mas_source", "rg", "phc0", "phc1",
    "sr_hz", "sf_MHz", "ds", "aq_s"})


def _fmt_khz(v):
    return f"{float(v) / 1000:.1f}"


def _fmt_g(v):
    return format(float(v), "g")


def _fmt_int(v):
    return str(int(v))


def _fmt_date(v):
    return str(v)[:10]


def _fmt_flip(v):
    return f"{float(v):.1f}"


def _fmt_sr_cell(v):
    return f"{float(v):.2f}"


def _fmt_mas(v):
    return "check" if v else ""


#: (key, header, unit, fmt, LaTeX header) -- unicode headers for the screen /
#: CSV / Markdown, TeX forms for pdflatex
COLUMNS: tuple[tuple, ...] = (
    ("sample", "Sample", "", str, "Sample"),
    ("expno", "EXPNO", "", str, "EXPNO"),
    ("nucleus", "Nucleus", "", str, "Nucleus"),
    ("sfo1_MHz", "ν₀ (MHz)", "MHz", lambda v: f"{float(v):.3f}", r"$\nu_0$ (MHz)"),
    ("b0_T", "B₀ (T)", "T", lambda v: f"{float(v):.2f}", r"$B_0$ (T)"),
    ("spectrometer", "Spectrometer", "", str, "Spectrometer"),
    ("probe", "Probe", "", str, "Probe"),
    ("rotor", "Rotor", "", str, "Rotor"),
    ("spin_rate_Hz", "νrot (kHz)", "kHz", _fmt_khz, r"$\nu_\mathrm{rot}$ (kHz)"),
    ("mas_uncertain", "MAS", "", _fmt_mas, "MAS"),
    ("pulprog", "Pulse program", "", str, "Pulse program"),
    ("p1_us", "P1 (µs)", "µs", _fmt_g, "P1 (us)"),
    ("flip_deg", "Flip (°)", "°", _fmt_flip, "Flip (deg)"),
    ("plw1_W", "PLW1 (W)", "W", _fmt_g, "PLW1 (W)"),
    ("d1_s", "D1 (s)", "s", _fmt_g, "D1 (s)"),
    ("ns", "NS", "", _fmt_int, "NS"),
    ("sw_Hz", "SW (kHz)", "kHz", lambda v: format(float(v) / 1000, "g"), "SW (kHz)"),
    ("td", "TD", "", _fmt_int, "TD"),
    ("tdeff", "TDeff", "", _fmt_int, "TDeff"),
    ("si", "SI", "", _fmt_int, "SI"),
    ("wdw", "Window", "", str, "Window"),
    ("lb_hz", "LB (Hz)", "Hz", _fmt_g, "LB (Hz)"),
    ("gb", "GB", "", _fmt_g, "GB"),
    ("ssb", "SSB", "", _fmt_g, "SSB"),
    ("ph_mod", "Phase", "", str, "Phase"),
    ("absg", "ABSG", "", _fmt_int, "ABSG"),
    ("sr_hz", "SR (Hz)", "Hz", _fmt_sr_cell, "SR (Hz)"),
    ("date_utc", "Date", "", _fmt_date, "Date"),
)


def _same(a, b) -> bool:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and \
            not isinstance(a, bool) and not isinstance(b, bool):
        fa, fb = float(a), float(b)
        return fa == fb or abs(fa - fb) <= 1e-6 * max(abs(fa), abs(fb))
    return a == b


def _cell(t_or_cols, key: str, v) -> str:
    """A cell through the column's formatter ('' for None)."""
    if v is None or v == "":
        return ""
    cols = t_or_cols.columns if isinstance(t_or_cols, AcqTable) else t_or_cols
    for k, _h, _u, fmt, *_ in cols:
        if k == key:
            try:
                return fmt(v)
            except (TypeError, ValueError):
                return str(v)
    return str(v)


def table(blocks, spin_rates: dict | None = None) -> AcqTable:
    """Table S1 over acquisition blocks: rows sorted (nucleus, date), the
    confirmed rate per spectrum from ``spin_rates`` ({expno_path: (rate_Hz,
    uncertain)}, e.g. from the fitted recipes) overriding the block's own,
    ``varying`` = the column keys whose values differ across rows (never the
    ALWAYS_VARIES keys), ``constants`` = value of every other column key
    that is present. GB / SSB columns appear only when a Gaussian / sine
    window is in use."""
    rows = []
    for b in blocks or []:
        if not b:
            continue
        r = dict(b)
        sp = (spin_rates or {}).get(r.get("expno_path"))
        if sp is not None:
            rate, unc = sp
            if rate is not None:
                r["spin_rate_Hz"] = float(rate)
                r["mas_uncertain"] = bool(unc)
        rows.append(r)
    rows.sort(key=lambda r: (str(r.get("nucleus") or ""), str(r.get("date_utc") or ""),
                             str(r.get("sample") or ""), str(r.get("expno") or "")))
    cols = []
    for c in COLUMNS:
        k = c[0]
        if k == "gb" and not any(r.get("wdw") == "GM" for r in rows):
            continue
        if k == "ssb" and not any(r.get("wdw") in ("SINE", "QSINE") for r in rows):
            continue
        cols.append(c)
    varying: set = set()
    constants: dict = {}
    for k, *_ in cols:
        vals = [r.get(k) for r in rows]
        present = [v for v in vals if v is not None and v != ""]
        if k in ALWAYS_VARIES:
            continue
        if not present:
            continue
        first = vals[0]
        if len(present) != len(vals) or any(not _same(v, first) for v in vals[1:]):
            if len(rows) > 1:
                varying.add(k)
        else:
            constants[k] = first
    return AcqTable(columns=cols, rows=rows, varying=varying, constants=constants)


def table_csv(t: AcqTable) -> str:
    """Every column, one row per spectrum, then a '# varies' row with 'yes'
    under every column that differs across the series."""
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow([h for _k, h, *_ in t.columns])
    for r in t.rows:
        w.writerow([_cell(t, k, r.get(k)) for k, *_ in t.columns])
    w.writerow(["# varies"] + ["yes" if k in t.varying else "" for k, *_ in t.columns[1:]])
    return buf.getvalue()


_TEX_MAP = (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"),
            ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"),
            ("^", r"\textasciicircum{}"), ("µ", "u"), ("°", " deg"), ("±", r"$\pm$"),
            ("–", "--"), ("−", "-"), ("→", r"$\rightarrow$"), ("×", r"$\times$"),
            ("δ", r"$\delta$"), ("Ξ", r"$\Xi$"), ("ν", r"$\nu$"), ("₀", "0"), ("⚠", ""),
            ("≈", r"$\approx$"), ("·", "."), ("¹", "1"), ("²", "2"), ("³", "3"), ("σ", r"$\sigma$"))


def _tex(s) -> str:
    """Plain text -> TeX: escapes plus the unicode the headers and cells use;
    anything else above U+007F degrades to '?' so pdflatex never chokes."""
    out = str(s)
    for a, b in _TEX_MAP:
        out = out.replace(a, b)
    return out.encode("ascii", "replace").decode("ascii")


def table_latex(t: AcqTable, caption: str = "", label: str = "tab:acq",
                compact: bool = True) -> str:
    """A booktabs ``table`` of the acquisition parameters (ASCII/TeX only).
    ``compact`` keeps Sample / EXPNO / Nucleus plus the varying columns and
    lists the constant parameters in the caption; otherwise every column is
    printed. Varying cells are set in ``\\textbf{}`` with a footnote."""
    keep = [c for c in t.columns
            if not compact or c[0] in ("sample", "expno", "nucleus") or c[0] in t.varying]
    cap = _tex(caption) if caption else "Acquisition and processing parameters"
    if compact:
        consts = []
        for k, h, _u, _f, tex_h in t.columns:
            if k in t.constants and k not in ("sample", "expno", "nucleus"):
                consts.append(f"{tex_h} = {_tex(_cell(t, k, t.constants[k]))}")
        if consts:
            cap += " Constant across the series: " + "; ".join(consts) + "."
    lines = [r"% LARMOR acquisition table (Table S1)", r"\begin{table}[htbp]", r"\centering",
             r"\caption{" + cap + "}", r"\label{" + label + "}",
             r"\begin{tabular}{" + "l" * len(keep) + "}", r"\toprule"]
    lines.append(" & ".join(c[4] for c in keep) + r" \\")
    lines.append(r"\midrule")
    for r in t.rows:
        cells = []
        for k, *_ in keep:
            txt = _tex(_cell(t, k, r.get(k)))
            cells.append(r"\textbf{" + txt + "}" if k in t.varying and txt else txt)
        lines.append(" & ".join(cells) + r" \\")
    lines.append(r"\bottomrule")
    if t.varying:
        lines.append(r"\multicolumn{" + str(len(keep)) + r"}{l}{\footnotesize "
                     r"Columns in bold vary across the series.} \\")
    lines += [r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def table_markdown(t: AcqTable) -> str:
    """Every column as a Markdown table; varying headers and cells in bold."""
    heads = [f"**{h}**" if k in t.varying else h for k, h, *_ in t.columns]
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(t.columns)]
    for r in t.rows:
        cells = []
        for k, *_ in t.columns:
            txt = _cell(t, k, r.get(k)).replace("|", "\\|")
            cells.append(f"**{txt}**" if k in t.varying and txt else txt)
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


# ------------------------------------------------------------- paragraph
def _merged_view(t: AcqTable) -> dict:
    """One block-like dict for a nucleus group: constant keys keep their
    value, varying keys become ranges."""
    view: dict = {}
    if not t.rows:
        return view
    keys = set(BLOCK_KEYS) | {"spin_rate_Hz", "mas_uncertain"}
    for k in keys:
        vals = [r.get(k) for r in t.rows]
        present = [v for v in vals if v is not None and v != ""]
        if k in t.varying:
            view[k] = _Range(present) if present else None
            continue
        if k in ALWAYS_VARIES and k not in ("mas_uncertain",):
            # per-spectrum keys: a shared value stays, a differing one is a
            # range too (SR / SF differ only when referencing differs)
            if present and all(_same(v, present[0]) for v in present):
                view[k] = present[0]
            elif present and k in ("sr_hz", "sf_MHz"):
                view[k] = _Range(present)
            else:
                view[k] = present[0] if present else None
            continue
        view[k] = vals[0] if vals else None
    view["mas_uncertain"] = any(bool(r.get("mas_uncertain")) for r in t.rows)
    return view


def paragraph(blocks, *, spin_rates: dict | None = None,
              referencing_text: str | None = None, larmor_ops=None) -> str:
    """The acquisition / processing part of the Experimental section for a
    set of blocks: one paragraph per nucleus (blank line between), every
    varying parameter as a range with '(Table S1)', the MAS rate confirmed
    per spectrum through ``spin_rates`` or bracketed. '' for no blocks."""
    blocks = [b for b in (blocks or []) if b]
    if not blocks:
        return ""
    groups: dict[str, list[dict]] = {}
    for b in sorted(blocks, key=lambda b: (str(b.get("nucleus") or ""),
                                           str(b.get("date_utc") or ""))):
        groups.setdefault(str(b.get("nucleus") or ""), []).append(b)
    paras = []
    for _nuc, grp in groups.items():
        t = table(grp, spin_rates)
        view = _merged_view(t)
        rate = view.get("spin_rate_Hz")
        unc = bool(view.get("mas_uncertain"))
        sents = sentences(view, spin_rate_Hz=rate, mas_uncertain=unc,
                          sr_hz=view.get("sr_hz"), referencing_text=referencing_text,
                          larmor_ops=larmor_ops, n_spectra=len(grp))
        paras.append(" ".join(sents))
    return "\n\n".join(paras)


# ------------------------------------------------------------ referencing
@lru_cache(maxsize=16)
def _audit_rows(root_str: str) -> dict:
    """{EXPNO path: AuditRow} of a session's referencing audit, once per
    root -- the evidence that a spectrum's SR is the one its 1H adamantane
    reference gives by indirect (Xi) referencing."""
    from larmor import referencing as R

    try:
        acqs = R.scan_session(Path(root_str))
        rows = R.audit(acqs)
    except Exception:                                    # noqa: BLE001
        return {}
    return {str(Path(r.acq.path)): r for r in rows}


def _ada_ppm_from(note: str) -> float | None:
    m = re.search(r"adamantane at\s*([-\d.]+)\s*ppm", note or "")
    return _f(m.group(1)) if m else None


def referencing_evidence(blocks, session_root=None, recipe: dict | None = None) -> str | None:
    """The referencing clause ONLY when the files prove it:

    * the recipe's ``provenance['referencing']`` (the audit re-referenced this
      fit to the session's 1H adamantane reference) -> its adamantane shift;
    * else every block's EXPNO gets an ``ok`` verdict in the session's
      referencing audit (``referencing.audit`` over ``session_root`` or the
      first block's month folder: a referenced 1H spectrum on the same
      magnet in the same month whose Xi prediction matches the stored SR).

    Returns 'referenced indirectly (IUPAC Ξ) to the ¹H resonance of
    adamantane at 1.82 ppm' or None (the paragraph then prints a bracket)."""
    from larmor.referencing import ADAMANTANE_1H_PPM, session_root as _root

    ref = ((recipe or {}).get("provenance") or {}).get("referencing") or {}
    if ref:
        ppm = _ada_ppm_from(str(ref.get("note") or "")) or ADAMANTANE_1H_PPM
        return (f"referenced indirectly (IUPAC Ξ) to the ¹H resonance of adamantane "
                f"at {ppm:g} ppm")
    blocks = [b for b in (blocks or []) if b and b.get("expno_path")]
    if not blocks:
        return None
    try:
        root = Path(session_root) if session_root else _root(blocks[0]["expno_path"])
        root = _root(root)
    except Exception:                                    # noqa: BLE001
        return None
    rows = _audit_rows(str(root))
    if not rows:
        return None
    for b in blocks:
        r = rows.get(str(Path(b["expno_path"])))
        if r is None or r.verdict != "ok":
            return None
    return (f"referenced indirectly (IUPAC Ξ) to the ¹H resonance of adamantane "
            f"at {ADAMANTANE_1H_PPM:g} ppm")
