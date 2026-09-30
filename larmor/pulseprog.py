"""Bruker pulse programs: find, parse, resolve against acqus, and lay out
for a TopSpin-like timing diagram.

Every Bruker EXPNO carries ``pulseprogram`` -- the sequence TopSpin compiled
for that acquisition (the cpp-preprocessed source: ``# 1 "..."`` line
markers, expanded ``#include``s and the ``mc`` macro spelled out as
MCWRK / MCREST delays). Newer TopSpin versions may write
``pulseprogram.precomp`` instead. This module turns that text into:

* :func:`parse` -> :class:`Program`: the header comments, the includes, the
  quoted definitions (``"p2=p1*2"``), the ``define delay`` /
  ``define loopcounter`` names, the body as a list of :class:`Statement`
  (delays, pulses, shaped pulses, decoupling on / off, acquisition, loops,
  labels, gotos, power settings, composite ``(center ...)`` groups), the
  phase programs (``ph1=(12) 0 2 4``, ``{0}*24``, continuation lines) and
  the legend comments (``;p1 : f1 channel - 90 degree pulse``);
* :func:`resolve` -> the value of every symbol the body uses, from the acqus
  arrays (P in us, D in s, PLW / SPW in W, L, CNST, IN in s, INF in us,
  PCPD in us ...) and from the definitions when they only use known
  symbols (a small safe arithmetic evaluator: numbers with u / m / s
  suffixes, ``+ - * /``, parentheses, ``trunc``);
* :func:`timeline` -> the rows of a NOT-to-scale diagram (one row per
  channel in use plus the acquisition row; every element a fixed nominal
  width; loop brackets with their counts) for the desktop widget to paint;
* :func:`describe` -> a one-line summary of the sequence.

Qt-free; the drawing lives in ``larmor.desktop.pulseprog_dialog``. Nothing
is ever written into an EXPNO folder.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "Statement", "Program", "Value", "Event", "Loop", "Timeline", "Record",
    "find_pulseprogram", "clean_text", "parse", "resolve", "timeline",
    "describe", "channel_nuclei", "load", "category", "evaluate",
]

# ----------------------------------------------------------------- syntax
#: kinds a :class:`Statement` can have
KINDS = ("delay", "pulse", "shaped_pulse", "decouple_on", "decouple_off",
         "acquire", "loop", "label", "goto", "set_power", "misc", "comment",
         "group")

_MARKER = re.compile(r'^#\s+\d+\s+"(?P<file>[^"]*)"(?P<flags>(?:\s+\d+)*)\s*$')
_INCLUDE = re.compile(r'^#\s*include\s*[<"]([^>"]+)[>"]')
_DEFINE = re.compile(r'^define\s+(delay|loopcounter|pulse|list)\b\s*'
                     r'(?:<[^>]*>\s*)?([A-Za-z_]\w*)')
_QUOTED_DEF = re.compile(r'^"\s*([A-Za-z_]\w*)\s*=\s*(.*?)\s*"\s*$')
_PHASE_LINE = re.compile(r'^(ph\d+)\s*=\s*(.*)$')
_PHASE_CONT = re.compile(r'^\s+[\d\s{}*^()]+$')
_LEGEND = re.compile(r'^;\s*([A-Za-z_]\w*)\s*:\s*(.*)$')
_META = re.compile(r'^;\s*\$(\w+)\s*=\s*(.*)$')
_DIMENSION = re.compile(r'dimension\s+(\S+)\s*;\s*AQ_mode\s*(.*)$')
_LEGEND_KEY = re.compile(
    r'^(?:p|d|pl|plw|sp|spw|spnam|spoffs|l|cnst|in|inf|pcpd|cpdprg|gp|gpnam|'
    r'gpz|ph|o|sfo|nuc|td|vd|vp|vc|va)\d*$'
    r'|^(?:ns|ds|aq|de|dw|vdlist|vplist|vclist|valist|zgoptns|fnmode|masr|'
    r'acqt0|nbl|digmod|aunm|tdav|rg)$', re.I)
_LITERAL = re.compile(r'^(\d+(?:\.\d+)?)(u|m|s)$')
_NUMBER = re.compile(r'^\d+(?:\.\d+)?$')
_DELAY_SYM = re.compile(r'^d\d+$')
_PULSE_SYM = re.compile(r'^p\d+$')
_PHASE_TOK = re.compile(r'^(ph\d+)(?::r)?$')
_POWER_TOK = re.compile(r'^(plw?\d+|spw?\d+|gp\d+)$')
_DECOUPLE_ON = re.compile(r'^(cpds?\d*|cw|cwlp)$')
_ACQUIRE = re.compile(r'^(go(?:=\d+)?|gosc|goscnp|gonp|adc)$')
_RCYC = re.compile(r'^rcyc(?:=(\d+))?$')
_MISC_WORDS = frozenset({
    "ze", "zd", "wr", "rf", "if", "st", "st0", "ivd", "ivc", "ivp", "iva",
    "rvd", "rvc", "rvp", "UNBLKGRAD", "BLKGRAD", "UNBLKGRAMP", "BLKGRAMP",
    "dccorr", "aqseq", "reset", "rpp", "syrec", "eoscnp", "eosc",
})
_INCREMENT = re.compile(r'^(?:ip|rp|id|dd|iu|ru|ipp|rpp|dp)\d+(?:\*\d+)?$')
_CHANNEL_SUFFIX = re.compile(r'^(.*?):(f\d)$')
_SYMBOL = re.compile(r'^([A-Za-z]+?)(\d+)$')
_SYMBOLS_IN_EXPR = re.compile(r'[A-Za-z_]\w*')
_MC_DELAYS = ("MCWRK", "MCREST")


# ------------------------------------------------------------- dataclasses
@dataclass
class Statement:
    """One element of the body.

    ``kind`` is one of :data:`KINDS`; ``channel`` is ``f1``.. (``grad`` for a
    gradient), ``name`` the symbol (``p1``, ``d1``, ``MCWRK``, ``10u``,
    ``cpds2``, ``go=2``), ``phase`` the phase program (``ph1``), ``power``
    the power level (``pl1`` / ``sp1`` / ``gp1``), ``duration_expr`` the
    delay / pulse length as written (``d1``, ``MCWRK*2``), ``times_expr``
    the loop count (``l20``, ``td1``, ``ns``), ``target_label`` a loop /
    goto target, ``children`` the members of a ``group`` (a simultaneous or
    ``center`` composite, or one channel's element list inside it),
    ``extras`` the housekeeping tokens that rode on the same line
    (``wr #0``, ``zd``, ``ivd``), ``cond`` the ``#ifdef`` branch the line
    belongs to, ``condition`` the test of an ``if "..." goto``."""

    kind: str
    raw: str = ""
    channel: str = ""
    name: str = ""
    phase: str = ""
    power: str = ""
    duration_expr: str = ""
    times_expr: str = ""
    target_label: str = ""
    children: list["Statement"] = field(default_factory=list)
    line_no: int = 0
    extras: list[str] = field(default_factory=list)
    cond: str = ""
    condition: str = ""

    def walk(self):
        yield self
        for c in self.children:
            yield from c.walk()


@dataclass
class Program:
    """A parsed pulse program (see the module docstring)."""

    name: str = ""
    header: list[str] = field(default_factory=list)
    includes: list[str] = field(default_factory=list)
    definitions: dict[str, str] = field(default_factory=dict)
    defined_delays: list[str] = field(default_factory=list)
    loopcounters: list[str] = field(default_factory=list)
    defined_pulses: list[str] = field(default_factory=list)
    statements: list[Statement] = field(default_factory=list)
    phase_programs: dict[str, list[int]] = field(default_factory=dict)
    phase_base: dict[str, int] = field(default_factory=dict)
    legend: dict[str, str] = field(default_factory=dict)
    meta: dict[str, str] = field(default_factory=dict)
    prosol: str = ""
    dimension: str = ""
    aq_mode: str = ""
    text: str = ""
    warnings: list[str] = field(default_factory=list)

    # ---- queries
    def walk(self):
        for s in self.statements:
            yield from s.walk()

    def channels(self) -> list[str]:
        """The rf channels the body drives, ``f1`` first."""
        seen = []
        for s in self.walk():
            if s.channel and s.kind in ("pulse", "shaped_pulse", "decouple_on",
                                        "decouple_off", "set_power") \
                    and s.channel not in seen:
                seen.append(s.channel)
        return sorted(seen, key=lambda c: (c == "grad", c))

    def legend_for(self, symbol: str) -> str:
        """The legend comment of ``symbol``; case-insensitive, and a power
        level is looked up under its watt-named twin too (the body says
        ``sp1`` / ``pl11``, the legend often documents ``spw1`` / ``plw11``)."""
        cands = [symbol]
        m = re.match(r"^(pl|plw|sp|spw)(\d+)$", symbol, re.I)
        if m:
            twin = {"pl": "plw", "plw": "pl", "sp": "spw", "spw": "sp"}[m.group(1).lower()]
            cands.append(f"{twin}{m.group(2)}")
        for cand in cands:
            if cand in self.legend:
                return self.legend[cand]
            low = cand.lower()
            for k, v in self.legend.items():
                if k.lower() == low:
                    return v
        return ""

    def phase_cycle_length(self) -> int:
        return max((len(v) for v in self.phase_programs.values()), default=0)

    def symbols_used(self) -> list[str]:
        """Every parameter symbol the body refers to, in order of first use
        (durations, powers, loop counts, decoupling programs, phases)."""
        out: list[str] = []

        def add(sym: str):
            if sym and sym not in out and not _NUMBER.match(sym) \
                    and not _LITERAL.match(sym):
                out.append(sym)

        for s in self.walk():
            if s.kind in ("delay", "pulse", "shaped_pulse"):
                for part in re.split(r"[*/+\-()]", s.duration_expr or s.name):
                    add(part.strip())
            if s.power:
                add(s.power)
                m = re.match(r"^sp(\d+)$", s.power)
                if m:
                    add(f"spnam{m.group(1)}")
            if s.kind in ("decouple_on",):
                m = re.match(r"^cpds?(\d+)$", s.name)
                if m:
                    add(f"cpdprg{m.group(1)}")
                    add(f"pcpd{m.group(1)}")
            if s.kind == "acquire":
                add("ns")
                add("aq")
            if s.kind == "loop" and s.times_expr:
                add(s.times_expr)
        for s in self.walk():
            if s.phase:
                add(s.phase)
        return out


@dataclass
class Value:
    """A resolved symbol: ``value`` in the display ``unit`` (us / ms / s /
    W / count / % / ''), ``si`` the same number in canonical units (seconds,
    watts, counts) for arithmetic, ``source`` where it came from (acqus,
    definition, literal, list, acqu2s, derived), ``note`` a short aside
    (the dB equivalent, the shape file, the definition), ``text`` the
    ready-to-print number with its unit."""

    symbol: str
    value: float | None = None
    unit: str = ""
    si: float | None = None
    source: str = ""
    note: str = ""
    text: str = ""
    kind: str = ""          # time | power | count | number | string | list | phase


@dataclass
class Event:
    """One drawn element: ``start`` / ``width`` in nominal units."""

    kind: str
    channel: str
    start: float
    width: float
    label: str = ""
    sub: str = ""
    value: str = ""
    note: str = ""
    cond: str = ""
    stmt: Statement | None = None

    @property
    def end(self) -> float:
        return self.start + self.width


@dataclass
class Loop:
    """A bracket under the rows: ``times`` as written (``l20`` / ``ns`` /
    ``2``), ``label`` the target label, ``level`` the bracket row (0 nearest
    the channels), ``text`` the printed count (``× ns = 512``)."""

    start: float
    end: float
    times: str
    label: str = ""
    level: int = 0
    text: str = ""


@dataclass
class Timeline:
    channels: list[str]
    events: list[Event]
    loops: list[Loop]
    total: float
    notes: list[str] = field(default_factory=list)
    labels: dict[str, float] = field(default_factory=dict)


@dataclass
class Record:
    """What :func:`load` returns for an EXPNO."""

    expno: Path
    path: Path | None
    program: Program | None
    values: dict[str, Value]
    nuclei: dict[str, str]
    acqus: dict
    pulprog: str
    warnings: list[str] = field(default_factory=list)


# ------------------------------------------------------------------ files
def find_pulseprogram(expno) -> Path | None:
    """The pulse program file of an EXPNO (``pulseprogram``, else
    ``pulseprogram.precomp``); a fid / ser / pdata path is walked up to
    the folder holding acqus. None when there is none."""
    p = Path(str(expno))
    if p.is_file():
        p = p.parent
    for cand in [p, *list(p.parents)[:4]]:
        if (cand / "acqus").is_file():
            p = cand
            break
    for name in ("pulseprogram", "pulseprogram.precomp"):
        f = p / name
        if f.is_file():
            return f
    return None


def clean_text(text: str) -> str:
    """The program text without the cpp line markers (``# 1 "..."``)."""
    out = []
    for line in text.splitlines():
        if _MARKER.match(line.strip()):
            continue
        out.append(line.rstrip())
    # collapse runs of blank lines the markers left behind
    res, blank = [], 0
    for line in out:
        if line.strip():
            blank = 0
        else:
            blank += 1
            if blank > 1:
                continue
        res.append(line)
    return "\n".join(res).strip("\n") + "\n"


# ---------------------------------------------------------------- parsing
def _strip_comment(line: str) -> str:
    out, inq = [], False
    for ch in line:
        if ch == '"':
            inq = not inq
        if ch == ";" and not inq:
            break
        out.append(ch)
    return "".join(out).strip()


def _paren_balance(s: str) -> int:
    bal, inq = 0, False
    for ch in s:
        if ch == '"':
            inq = not inq
        elif not inq:
            if ch == "(":
                bal += 1
            elif ch == ")":
                bal -= 1
    return bal


def _split_label(s: str) -> tuple[str, str]:
    m = re.match(r"^(\d+)(?:\s+(.*))?$", s)
    if m:
        return m.group(1), (m.group(2) or "").strip()
    m = re.match(r"^([A-Za-z_]\w*)\s*,\s*(.*)$", s)
    if m:
        return m.group(1), m.group(2).strip()
    return "", s


def _tokens(s: str) -> list[str]:
    """Split a body line into tokens; quoted strings and balanced
    parenthesised groups (with a glued ``:fN`` suffix) stay whole, and a
    ``(`` glued to an identifier (``F0(zd)``, ``trunc(td1/2)``) is part of
    that identifier's token."""
    out: list[str] = []
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c.isspace():
            i += 1
            continue
        if c == '"':
            j = s.find('"', i + 1)
            j = n if j < 0 else j + 1
            out.append(s[i:j])
            i = j
            continue
        if c == "(":
            depth, j = 0, i
            while j < n:
                if s[j] == "(":
                    depth += 1
                elif s[j] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            j = min(j + 1, n)
            m = re.match(r":f\d", s[j:j + 3])
            if m:
                j += 3
            out.append(s[i:j])
            i = j
            continue
        if c == ")":
            i += 1
            continue
        j = i
        while j < n and not s[j].isspace() and s[j] != '"':
            if s[j] == "(":
                depth = 0
                while j < n:
                    if s[j] == "(":
                        depth += 1
                    elif s[j] == ")":
                        depth -= 1
                        if depth == 0:
                            break
                    j += 1
                j = min(j + 1, n)
                continue
            if s[j] == ")":
                break
            j += 1
        out.append(s[i:j])
        i = j
    return out


def _split_suffix(tok: str) -> tuple[str, str]:
    m = _CHANNEL_SUFFIX.match(tok)
    return (m.group(1), m.group(2)) if m else (tok, "")


class _Parser:
    def __init__(self, program: Program):
        self.p = program
        self.delay_names = set(_MC_DELAYS)
        self.pulse_names: set[str] = set()

    # ---- predicates
    def is_delay(self, base: str) -> bool:
        head = base.split("*", 1)[0]
        return bool(_DELAY_SYM.match(head) or _LITERAL.match(head)
                    or head in ("vd", "aq") or head in self.delay_names)

    def is_pulse(self, base: str) -> bool:
        head = base.split(":", 1)[0].split("*", 1)[0]
        return bool(_PULSE_SYM.match(head) or head == "vp"
                    or head in self.pulse_names)

    # ---- elements
    def elements(self, tokens: list[str], channel: str, line_no: int, raw: str,
                 cond: str, top: bool) -> list[Statement]:
        stmts: list[Statement] = []
        cur: Statement | None = None
        pending_groups: list[Statement] = []

        def flush_groups():
            nonlocal cur
            if not pending_groups:
                return
            if len(pending_groups) == 1:
                g = pending_groups[0]
                if g.name == "center":
                    stmts.append(g)
                else:
                    stmts.extend(g.children)
            else:
                stmts.append(Statement("group", raw, name="sim", line_no=line_no,
                                       children=list(pending_groups), cond=cond))
            pending_groups.clear()
            cur = None

        def new(st: Statement) -> Statement:
            st.line_no = st.line_no or line_no
            st.raw = st.raw or raw
            st.cond = cond
            stmts.append(st)
            return st

        i, n = 0, len(tokens)
        while i < n:
            t = tokens[i]
            if t.startswith("("):
                pending_groups.append(self.group(t, line_no, raw, cond, channel))
                i += 1
                continue
            flush_groups()
            base, chan = _split_suffix(t)
            ch = chan or channel
            # multiplier that follows a delay / pulse: "MCWRK * 2"
            if t == "*" and cur is not None and cur.kind in ("delay", "pulse", "shaped_pulse") \
                    and i + 1 < n:
                cur.duration_expr = f"{cur.duration_expr}*{tokens[i + 1]}"
                i += 2
                continue
            if t == "lo":
                target = tokens[i + 2] if i + 2 < n and tokens[i + 1] == "to" else ""
                times = tokens[i + 4] if i + 4 < n and tokens[i + 3] == "times" else ""
                cur = new(Statement("loop", name="lo", target_label=target, times_expr=times))
                i += 5
                continue
            if t == "goto":
                cur = new(Statement("goto", name="goto",
                                    target_label=tokens[i + 1] if i + 1 < n else ""))
                i += 2
                continue
            if t == "if" and i + 1 < n and tokens[i + 1].startswith('"'):
                condition = tokens[i + 1].strip('"')
                target = tokens[i + 3] if i + 3 < n and tokens[i + 2] == "goto" else ""
                cur = new(Statement("goto", name="if", target_label=target, condition=condition))
                i += 4
                continue
            if t == "print":
                msg = tokens[i + 1] if i + 1 < n else ""
                new(Statement("misc", name="print", extras=[msg.strip('"')]))
                i += 2
                continue
            if t == "mc":
                rest = tokens[i:]
                target = ""
                if "to" in rest:
                    k = rest.index("to")
                    target = rest[k + 1] if k + 1 < len(rest) else ""
                times = "td1" if any(x.startswith(("F1", "F2", "F3")) for x in rest) else "td0"
                if cur is not None and cur.kind == "delay":
                    cur.extras.append(" ".join(rest))
                new(Statement("loop", name="mc", target_label=target, times_expr=times,
                              extras=[" ".join(rest)]))
                i = n
                continue
            m = _RCYC.match(base)
            if m:
                new(Statement("loop", name=base, target_label=m.group(1) or "",
                              times_expr="ns"))
                cur = None
                i += 1
                continue
            if _ACQUIRE.match(base):
                target = base.split("=", 1)[1] if "=" in base else ""
                cur = new(Statement("acquire", name=base, target_label=target, channel=""))
                i += 1
                continue
            m = _PHASE_TOK.match(base)
            if m:
                if cur is not None and cur.kind in ("pulse", "shaped_pulse", "acquire"):
                    cur.phase = m.group(1)
                else:
                    new(Statement("misc", name=t))
                i += 1
                continue
            if base == "do":
                new(Statement("decouple_off", name="do", channel=ch or "f2"))
                i += 1
                continue
            if _DECOUPLE_ON.match(base):
                new(Statement("decouple_on", name=base, channel=ch or "f2"))
                i += 1
                continue
            m = _POWER_TOK.match(base)
            if m:
                if chan or cur is None or cur.kind not in ("pulse", "shaped_pulse"):
                    new(Statement("set_power", name=base, power=base, channel=ch or "f1"))
                else:
                    cur.power = base
                i += 1
                continue
            if self.is_delay(base):
                head, _, mult = base.partition("*")
                cur = new(Statement("delay", name=head, channel=chan if chan else ("" if top else channel),
                                    duration_expr=base))
                i += 1
                continue
            if self.is_pulse(base):
                name, _, shape = base.partition(":")
                head, _, mult = name.partition("*")
                kind = "pulse"
                power = ""
                pch = ch or "f1"
                if shape.startswith("sp"):
                    kind, power = "shaped_pulse", shape
                elif shape.startswith("gp"):
                    power, pch = shape, "grad"
                cur = new(Statement(kind, name=head, channel=pch, power=power,
                                    duration_expr=name))
                i += 1
                continue
            if t in ("#0", "#1", "#2") and cur is not None and cur.extras:
                cur.extras[-1] = f"{cur.extras[-1]} {t}"
                i += 1
                continue
            if base in ("ze", "zd") and cur is not None and cur.kind == "delay":
                cur.extras.append(base)
                i += 1
                continue
            if base in ("ze", "zd", "exit") or base in _MISC_WORDS or _INCREMENT.match(base) \
                    or base.startswith(("fq", "trig", "setnmr", "setrtp", "#")):
                if cur is not None and cur.kind == "delay":
                    cur.extras.append(t)
                else:
                    new(Statement("misc", name=t))
                i += 1
                continue
            # anything else: keep it visible rather than dropping it
            if cur is not None and cur.kind == "delay":
                cur.extras.append(t)
            else:
                new(Statement("misc", name=t))
            i += 1
        flush_groups()
        return stmts

    def group(self, tok: str, line_no: int, raw: str, cond: str, outer_channel: str) -> Statement:
        """``(p1 pl1 ph1):f1`` -> a group of that channel's elements;
        ``(center (...):f2 (...):f1)`` -> a centred composite of groups."""
        body, chan = _split_suffix(tok)
        inner = body[1:-1] if body.startswith("(") and body.endswith(")") else body.lstrip("(")
        toks = _tokens(inner)
        if toks and toks[0] in ("center", "align"):
            kids = []
            loose = [t for t in toks[1:] if not t.startswith("(")]
            for t in toks[1:]:
                if t.startswith("("):
                    kids.append(self.group(t, line_no, raw, cond, outer_channel))
            if loose:
                kids.append(Statement("group", raw, channel=chan or outer_channel or "f1",
                                      line_no=line_no, cond=cond,
                                      children=self.elements(loose, chan or outer_channel or "f1",
                                                             line_no, raw, cond, False)))
            return Statement("group", raw, name=toks[0], line_no=line_no, cond=cond,
                             children=kids)
        channel = chan or outer_channel or "f1"
        kids = self.elements(toks, channel, line_no, raw, cond, False)
        return Statement("group", raw, channel=channel, line_no=line_no, cond=cond,
                         children=kids)


def _phase_list(spec: str) -> tuple[list[int], int]:
    base = 4
    m = re.match(r"^\((\d+)\)\s*(.*)$", spec.strip())
    if m:
        base = int(m.group(1))
        spec = m.group(2)
    out: list[int] = []
    for t in re.findall(r"\{[^}]*\}(?:\*\d+)?(?:\^\d+)*|\d+", spec):
        if t.startswith("{"):
            m2 = re.match(r"^\{([^}]*)\}(?:\*(\d+))?((?:\^\d+)*)$", t)
            if not m2:
                continue
            vals = [int(x) for x in re.findall(r"\d+", m2.group(1))]
            seq = vals * int(m2.group(2) or 1)
            full = list(seq)
            for a in (int(x) for x in re.findall(r"\d+", m2.group(3) or "")):
                full += [(v + a) % base for v in seq]
            out += full
        else:
            out.append(int(t))
    return out, base


def parse(text: str, name: str = "") -> Program:
    """Parse a pulse program (compiled ``pulseprogram`` or plain source)
    into a :class:`Program`. Never raises on odd input: a line that does
    not parse becomes a ``misc`` statement and a warning."""
    prog = Program(name=name, text=clean_text(text))
    prs = _Parser(prog)
    # pre-pass: the defined names (their legend lines may precede them)
    for m in re.finditer(r'^\s*define\s+(delay|loopcounter|pulse|list)\b\s*(?:<[^>]*>\s*)?'
                         r'([A-Za-z_]\w*)', text, flags=re.M):
        kind, sym = m.group(1), m.group(2)
        if kind == "delay":
            prs.delay_names.add(sym)
        elif kind == "pulse":
            prs.pulse_names.add(sym)
    defined_names = set(prs.delay_names) | set(prs.pulse_names) | {
        m.group(1) for m in re.finditer(r'^\s*define\s+loopcounter\s+([A-Za-z_]\w*)',
                                        text, flags=re.M)}

    state = "head"
    main_file = ""
    cond_stack: list[str] = []
    lines = text.splitlines()
    i = 0
    pending: tuple[str, int] | None = None     # multi-line parenthesised statement
    last_phase = ""
    while i < len(lines):
        raw_line = lines[i]
        line_no = i + 1
        i += 1
        stripped = raw_line.strip()
        if not stripped:
            continue
        m = _MARKER.match(stripped)
        if m:
            f = m.group("file")
            flags = m.group("flags").split()
            if not main_file and not f.startswith("<"):
                main_file = f
                if not prog.name:
                    prog.name = Path(f).name
            elif flags and flags[0] == "1" and not f.startswith(("mc_line", "<")):
                inc = Path(f).name
                if inc not in prog.includes:
                    prog.includes.append(inc)
            continue
        # comments
        if stripped.startswith(";"):
            mm = _META.match(stripped)
            if mm:
                prog.meta[mm.group(1)] = mm.group(2).strip()
            md = _DIMENSION.search(stripped)
            if md:
                prog.dimension, prog.aq_mode = md.group(1), md.group(2).strip()
            ml = _LEGEND.match(stripped)
            if ml and (_LEGEND_KEY.match(ml.group(1)) or ml.group(1) in defined_names):
                prog.legend.setdefault(ml.group(1), ml.group(2).strip())
            if state == "head":
                prog.header.append(stripped[1:].rstrip())
            elif state == "body":
                prog.statements.append(Statement("comment", raw_line.rstrip(), name=stripped[1:].strip(),
                                                 line_no=line_no, cond=" & ".join(cond_stack)))
            continue
        # preprocessor conditionals (plain sources; compiled files have none)
        if stripped.startswith("#"):
            mi = _INCLUDE.match(stripped)
            if mi:
                if mi.group(1) not in prog.includes:
                    prog.includes.append(mi.group(1))
                continue
            low = stripped.split()
            if low[0] in ("#ifdef", "#ifndef") and len(low) > 1:
                cond_stack.append(("" if low[0] == "#ifdef" else "!") + low[1])
            elif low[0] == "#else" and cond_stack:
                c = cond_stack.pop()
                cond_stack.append(c[1:] if c.startswith("!") else "!" + c)
            elif low[0] == "#endif" and cond_stack:
                cond_stack.pop()
            elif low[0] in ("#define",):
                pass
            continue
        code = _strip_comment(raw_line)
        if not code:
            continue
        # declarations (anywhere before the body; tolerated inside it too)
        md = _DEFINE.match(code)
        if md:
            kind, sym = md.group(1), md.group(2)
            {"delay": prog.defined_delays, "loopcounter": prog.loopcounters,
             "pulse": prog.defined_pulses}.get(kind, prog.defined_delays).append(sym)
            continue
        mq = _QUOTED_DEF.match(code)
        if mq:
            prog.definitions[mq.group(1)] = mq.group(2)
            continue
        if code.startswith("prosol"):
            prog.prosol = code.split("=", 1)[1].strip("<> ") if "=" in code else code
            continue
        if state == "tail":
            mp = _PHASE_LINE.match(code)
            if mp:
                vals, base = _phase_list(mp.group(2))
                prog.phase_programs[mp.group(1)] = vals
                prog.phase_base[mp.group(1)] = base
                last_phase = mp.group(1)
            elif last_phase and _PHASE_CONT.match(raw_line.split(";")[0].rstrip()):
                vals, _ = _phase_list(code)
                prog.phase_programs[last_phase].extend(vals)
            else:
                last_phase = ""
            continue
        if state == "head":
            if code in ("dccorr",):
                prog.statements.append(Statement("misc", raw_line.rstrip(), name=code, line_no=line_no))
                continue
            state = "body"
        # ---- body
        raw = raw_line.rstrip()
        if pending is not None:
            code = f"{pending[0]} {code}"
            line_no = pending[1]
            raw = code                       # a multi-line statement keeps its joined text
            pending = None
        if _paren_balance(code) > 0:
            pending = (code, line_no)
            continue
        if code == "exit":
            prog.statements.append(Statement("misc", raw, name="exit", line_no=line_no))
            state = "tail"
            continue
        cond = " & ".join(cond_stack)
        label, rest = _split_label(code)
        if label:
            prog.statements.append(Statement("label", raw, name=label, line_no=line_no, cond=cond))
        if not rest:
            continue
        try:
            prog.statements.extend(prs.elements(_tokens(rest), "", line_no, raw, cond, True))
        except Exception as exc:                           # noqa: BLE001 -- never lose a line
            prog.warnings.append(f"line {line_no}: could not parse '{rest}' ({exc})")
            prog.statements.append(Statement("misc", raw, name=rest, line_no=line_no, cond=cond))
    if pending is not None:
        prog.warnings.append(f"line {pending[1]}: unbalanced parentheses")
        prog.statements.append(Statement("misc", pending[0], name=pending[0], line_no=pending[1]))
    if not prog.name and prog.header:
        first = prog.header[0].strip()          # ";zg" -- the conventional first line
        if first and not first.startswith("$") and " " not in first:
            prog.name = first
    return prog


# --------------------------------------------------------------- resolver
_UNITS = {"u": 1e-6, "m": 1e-3, "s": 1.0, "n": 1e-9}
_FUNCS = {"trunc": math.trunc, "int": int, "abs": abs, "sqrt": math.sqrt,
          "floor": math.floor, "ceil": math.ceil, "round": round,
          "log": math.log, "log10": math.log10, "exp": math.exp,
          "max": max, "min": min}


class _Unknown(Exception):
    pass


def _expr_tokens(expr: str) -> list[tuple[str, object]]:
    toks: list[tuple[str, object]] = []
    i, n = 0, len(expr)
    while i < n:
        c = expr[i]
        if c.isspace():
            i += 1
            continue
        if c.isdigit() or (c == "." and i + 1 < n and expr[i + 1].isdigit()):
            j = i
            while j < n and (expr[j].isdigit() or expr[j] == "."):
                j += 1
            if j < n and expr[j] in "eE" and j + 1 < n and (expr[j + 1].isdigit() or expr[j + 1] in "+-"):
                k = j + 2
                while k < n and expr[k].isdigit():
                    k += 1
                j = k
            val = float(expr[i:j])
            if j < n and expr[j] in _UNITS and (j + 1 == n or not (expr[j + 1].isalnum() or expr[j + 1] == "_")):
                val *= _UNITS[expr[j]]
                j += 1
            toks.append(("num", val))
            i = j
            continue
        if c.isalpha() or c == "_":
            j = i
            while j < n and (expr[j].isalnum() or expr[j] == "_"):
                j += 1
            toks.append(("sym", expr[i:j]))
            i = j
            continue
        if c in "+-*/(),":
            toks.append(("op", c))
            i += 1
            continue
        raise ValueError(f"bad character {c!r} in {expr!r}")
    return toks


class _Eval:
    def __init__(self, toks, env):
        self.t, self.env, self.i = toks, env, 0

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else (None, None)

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def expr(self):
        v = self.term()
        while self.peek() == ("op", "+") or self.peek() == ("op", "-"):
            op = self.take()[1]
            w = self.term()
            v = v + w if op == "+" else v - w
        return v

    def term(self):
        v = self.factor()
        while self.peek() == ("op", "*") or self.peek() == ("op", "/"):
            op = self.take()[1]
            w = self.factor()
            v = v * w if op == "*" else v / w
        return v

    def factor(self):
        kind, val = self.take()
        if kind == "op" and val in "+-":
            v = self.factor()
            return -v if val == "-" else v
        if kind == "num":
            return val
        if kind == "sym":
            if self.peek() == ("op", "("):
                self.take()
                args = [self.expr()]
                while self.peek() == ("op", ","):
                    self.take()
                    args.append(self.expr())
                if self.take() != ("op", ")"):
                    raise ValueError("missing )")
                fn = _FUNCS.get(val)
                if fn is None:
                    raise _Unknown(val)
                return float(fn(*args))
            if val not in self.env:
                raise _Unknown(val)
            return float(self.env[val])
        if kind == "op" and val == "(":
            v = self.expr()
            if self.take() != ("op", ")"):
                raise ValueError("missing )")
            return v
        raise ValueError("unexpected token")


def evaluate(expr: str, env: dict) -> float | None:
    """Evaluate a pulse-program arithmetic expression (``1s*l1/cnst31``,
    ``trunc(td1 / 2)``, ``-p1*2/3.1416``) with ``env`` holding symbol
    values in canonical units (seconds, watts, counts). None when a symbol
    is unknown or the expression does not evaluate."""
    try:
        toks = _expr_tokens(expr)
        if not toks:
            return None
        ev = _Eval(toks, env)
        v = ev.expr()
        if ev.i != len(toks):
            return None
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return None
        return float(v)
    except (_Unknown, ValueError, ZeroDivisionError, TypeError, OverflowError, IndexError):
        return None


def _fmt_time(seconds: float) -> tuple[float, str, str]:
    """(value, unit, text) with the unit chosen by magnitude."""
    a = abs(seconds)
    if a == 0:
        return 0.0, "s", "0 s"
    if a >= 1:
        v, u = seconds, "s"
    elif a >= 1e-3:
        v, u = seconds * 1e3, "ms"
    else:
        v, u = seconds * 1e6, "µs"
    v = float(f"{v:.10g}")                  # 2.1e-6 * 1e6 is 2.0999999999999996
    return v, u, f"{_g(v)} {u}"


def _g(v: float, digits: int = 4) -> str:
    if v == int(v) and abs(v) < 1e6:
        return str(int(v))
    s = f"{v:.{digits}g}"
    if "e" in s:
        s = f"{v:.{digits}f}".rstrip("0").rstrip(".")
    return s


def _arr(acqus: dict, key: str, idx: int):
    v = acqus.get(key)
    if isinstance(v, (list, tuple)):
        return v[idx] if 0 <= idx < len(v) else None
    return v if idx == 0 else None


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


_ARRAYS = {
    "p": ("P", 1e-6, "time"), "d": ("D", 1.0, "time"), "pl": ("PLW", 1.0, "power"),
    "plw": ("PLW", 1.0, "power"), "sp": ("SPW", 1.0, "power"), "spw": ("SPW", 1.0, "power"),
    "l": ("L", 1.0, "count"), "cnst": ("CNST", 1.0, "number"), "in": ("IN", 1.0, "time"),
    "inf": ("INF", 1e-6, "time"), "pcpd": ("PCPD", 1e-6, "time"), "gp": ("GPZ", 1.0, "percent"),
    "gpz": ("GPZ", 1.0, "percent"), "spoffs": ("SPOFFS", 1.0, "number"),
}
_STRINGS = {"spnam": "SPNAM", "cpdprg": "CPDPRG", "gpnam": "GPNAM"}
_SCALARS = {"ns": ("NS", 1.0, "count"), "ds": ("DS", 1.0, "count"), "td": ("TD", 1.0, "count"),
            "td0": ("TD0", 1.0, "count"), "de": ("DE", 1e-6, "time"), "tdav": ("TDav", 1.0, "count"),
            "nbl": ("NBL", 1.0, "count"), "rg": ("RG", 1.0, "number")}
_LISTS = {"vd": ("VDLIST", "delay list"), "vp": ("VPLIST", "pulse list"),
          "vc": ("VCLIST", "counter list"), "va": ("VALIST", "amplitude list")}


def category(symbol: str) -> str:
    """pulse | delay | power | shape | loop | constant | increment |
    decoupling | acquisition | list | phase | other -- for grouping."""
    s = symbol.lower()
    if re.match(r"^p\d+$", s) or s == "vp":
        return "pulse"
    if re.match(r"^d\d+$", s) or s in ("vd", "aq", "de", "dw"):
        return "delay"
    if re.match(r"^(pl|plw|sp|spw)\d+$", s):
        return "power"
    if re.match(r"^(spnam|spoffs|gpnam)\d+$", s):
        return "shape"
    if re.match(r"^l\d+$", s) or s in ("ns", "ds", "td", "td0", "td1", "tdav", "nbl", "vc"):
        return "loop"
    if re.match(r"^cnst\d+$", s):
        return "constant"
    if re.match(r"^(in|inf)\d+$", s):
        return "increment"
    if re.match(r"^(cpdprg|pcpd)\d+$", s):
        return "decoupling"
    if re.match(r"^ph\d+$", s):
        return "phase"
    if re.match(r"^gp\d+$", s):
        return "gradient"
    return "other"


def _make(symbol: str, si: float | None, kind: str, source: str, note: str = "") -> Value:
    v = Value(symbol, si=si, source=source, note=note, kind=kind)
    if si is None:
        return v
    if kind == "time":
        v.value, v.unit, v.text = _fmt_time(si)
    elif kind == "power":
        v.value, v.unit, v.text = si, "W", f"{_g(si, 3)} W"
        if si > 0:
            db = -10.0 * math.log10(si)
            v.note = (f"{db:+.1f} dB" + (f" · {note}" if note else ""))
        elif not note:
            v.note = "0 W"
    elif kind == "count":
        v.value, v.unit, v.text = si, "count", _g(si)
    elif kind == "percent":
        v.value, v.unit, v.text = si, "%", f"{_g(si)} %"
    else:
        v.value, v.unit, v.text = si, "", _g(si, 6)
    return v


def _from_acqus(symbol: str, acqus: dict, acqu2s: dict | None) -> Value | None:
    low = symbol.lower()
    if low in _SCALARS:
        key, scale, kind = _SCALARS[low]
        x = _num(acqus.get(key))
        return _make(symbol, x * scale if x is not None else None, kind, "acqus")
    if low == "td1":
        x = _num((acqu2s or {}).get("TD"))
        return _make(symbol, x, "count", "acqu2s") if x is not None else None
    if low == "aq":
        td, sw = _num(acqus.get("TD")), _num(acqus.get("SW_h"))
        if td and sw:
            return _make(symbol, td / (2.0 * sw), "time", "derived", "TD / (2·SW_h)")
        return None
    if low == "dw":
        sw = _num(acqus.get("SW_h"))
        return _make(symbol, 1.0 / (2.0 * sw), "time", "derived", "1 / (2·SW_h)") if sw else None
    if low in _LISTS:
        key, what = _LISTS[low]
        name = str(acqus.get(key, "") or "").strip("<> ")
        v = Value(symbol, source="list", kind="list", unit="list", text=key,
                  note=what + (f" {name}" if name else ""))
        return v
    m = _SYMBOL.match(low)
    if not m:
        return None
    prefix, idx = m.group(1), int(m.group(2))
    if prefix in _STRINGS:
        s = _arr(acqus, _STRINGS[prefix], idx)
        if s is None:
            return None
        s = str(s).strip("<> ")
        v = Value(symbol, source="acqus", kind="string", text=s or "(none)")
        return v
    if prefix in _ARRAYS:
        key, scale, kind = _ARRAYS[prefix]
        x = _num(_arr(acqus, key, idx))
        if x is None:
            if prefix in ("pl", "plw", "sp", "spw"):
                # old TopSpin (dB only): PL / SP arrays
                x = _num(_arr(acqus, "PL" if prefix.startswith("pl") else "SP", idx))
                if x is None:
                    return None
                v = Value(symbol, value=x, unit="dB", si=None, source="acqus", kind="power",
                          text=f"{_g(x, 4)} dB")
                return v
            return None
        note = ""
        if prefix in ("sp", "spw"):
            shape = _arr(acqus, "SPNAM", idx)
            if shape:
                note = f"shape {str(shape).strip('<> ')}"
        return _make(symbol, x * scale, kind, "acqus", note)
    return None


def _kind_for(symbol: str, program: Program) -> str:
    s = symbol.lower()
    if re.match(r"^(p|d|in|inf|pcpd)\d+$", s) or s in ("acqt0", "aq", "de", "dw") \
            or symbol in program.defined_delays or symbol in program.defined_pulses \
            or symbol in _MC_DELAYS:
        return "time"
    if re.match(r"^(pl|plw|sp|spw)\d+$", s):
        return "power"
    if re.match(r"^(l|td)\d*$", s) or symbol in program.loopcounters or s in ("ns", "ds", "tdav"):
        return "count"
    return "number"


def _expr_symbols(expr: str) -> list[str]:
    out = []
    for m in _SYMBOLS_IN_EXPR.finditer(expr):
        s = m.group(0)
        # a unit suffix glued to a number ("1s", "30m") is not a symbol
        start = m.start()
        if start > 0 and (expr[start - 1].isdigit() or expr[start - 1] == "."):
            continue
        if s not in _FUNCS and s not in out:
            out.append(s)
    return out


def resolve(program: Program, acqus: dict | None, acqu2s: dict | None = None) -> dict[str, Value]:
    """The value of every symbol the body uses (plus the symbols the
    definitions of those need), from the acqus arrays and the quoted
    definitions. A definition that evaluates wins over the acqus entry of
    the same symbol (TopSpin recomputes it at run time); the acqus number
    is kept in the note when the two differ."""
    acqus = acqus or {}
    wanted = list(program.symbols_used())
    if "acqt0" in program.definitions and "acqt0" not in wanted:
        wanted.append("acqt0")
    k = 0                                   # closure over the definitions' operands
    while k < len(wanted):
        s = wanted[k]
        k += 1
        if s in program.definitions:
            for dep in _expr_symbols(program.definitions[s]):
                if dep not in wanted:
                    wanted.append(dep)
    values: dict[str, Value] = {}
    env: dict[str, float] = {}
    for sym in wanted:
        if _LITERAL.match(sym):
            m = _LITERAL.match(sym)
            si = float(m.group(1)) * _UNITS[m.group(2)]
            values[sym] = _make(sym, si, "time", "literal")
            env[sym] = si
            continue
        if sym.startswith("ph"):
            vals = program.phase_programs.get(sym)
            base = program.phase_base.get(sym, 4)
            v = Value(sym, kind="phase", source="program" if vals is not None else "",
                      unit=f"× {360 // base if base else 90}°" if vals is not None else "")
            v.text = " ".join(str(x) for x in vals) if vals is not None else ""
            v.note = f"{len(vals)} steps" if vals else "no phase program"
            values[sym] = v
            continue
        v = _from_acqus(sym, acqus, acqu2s)
        if v is not None:
            values[sym] = v
            if v.si is not None:
                env[sym] = v.si
    # definitions, iterated so that a definition may use another
    pending = {s: e for s, e in program.definitions.items() if s in wanted}
    for _ in range(12):
        progress = False
        for sym, expr in list(pending.items()):
            si = evaluate(expr, env)
            if si is None:
                continue
            kind = _kind_for(sym, program)
            prev = values.get(sym)
            v = _make(sym, si, kind, "definition", "")
            note = f"= {expr}"
            if prev is not None and prev.si is not None and prev.source == "acqus":
                if not math.isclose(prev.si, si, rel_tol=5e-3, abs_tol=1e-12):
                    note += f" (acqus has {prev.text})"
            if kind == "power" and v.note:
                note = f"{v.note} · {note}"
            v.note = note
            values[sym] = v
            env[sym] = si
            del pending[sym]
            progress = True
        if not progress:
            break
    for sym, expr in pending.items():
        v = values.get(sym)
        if v is None:
            v = Value(sym, kind=_kind_for(sym, program), source="definition")
            values[sym] = v
        missing = [s for s in _expr_symbols(expr) if s not in env]
        v.note = f"= {expr}" + (f" (needs {', '.join(missing)})" if missing else " (could not evaluate)")
    return values


# ---------------------------------------------------------------- timeline
W_PULSE, W_SHAPED, W_DELAY, W_SHORT, W_ACQ, W_POWER = 1.0, 2.0, 2.0, 1.0, 4.0, 0.7


def _delay_width(stmt: Statement) -> float:
    if stmt.name in _MC_DELAYS or _LITERAL.match(stmt.name):
        return W_SHORT
    return W_DELAY


def _value_text(values: dict[str, Value], sym: str) -> str:
    v = values.get(sym)
    return v.text if v is not None and v.text else ""


def _duration_text(values: dict[str, Value], expr: str) -> str:
    """The resolved text of ``d1`` or ``MCWRK*2``."""
    if not expr:
        return ""
    m = _LITERAL.match(expr)
    if m:
        return _fmt_time(float(m.group(1)) * _UNITS[m.group(2)])[2]
    if "*" in expr or "/" in expr:
        env = {k: v.si for k, v in values.items() if v.si is not None}
        si = evaluate(expr, env)
        return _fmt_time(si)[2] if si is not None else ""
    return _value_text(values, expr)


def timeline(program: Program, values: dict[str, Value] | None = None) -> Timeline:
    """Lay the body out left to right in nominal width units."""
    values = values or {}
    events: list[Event] = []
    loops: list[Loop] = []
    notes: list[str] = []
    labels: dict[str, float] = {}
    open_dec: dict[str, tuple[float, Statement]] = {}
    unresolved: list[str] = []

    def dur(stmt: Statement) -> str:
        t = _duration_text(values, stmt.duration_expr or stmt.name)
        if not t and stmt.kind in ("delay", "pulse", "shaped_pulse") \
                and stmt.name not in ("vd", "vp") and not _LITERAL.match(stmt.name):
            if stmt.name not in unresolved:
                unresolved.append(stmt.name)
        return t

    def measure(stmts: list[Statement]) -> float:
        """Dry-run width of a sequential list."""
        x = 0.0
        for s in stmts:
            x += width_of(s)
        return x

    def width_of(s: Statement) -> float:
        if s.kind == "delay":
            return _delay_width(s)
        if s.kind == "pulse":
            return W_PULSE
        if s.kind == "shaped_pulse":
            return W_SHAPED
        if s.kind == "acquire":
            return W_ACQ
        if s.kind == "set_power":
            return W_POWER
        if s.kind == "group":
            if s.name in ("sim", "center", "align"):
                return max((measure(c.children) if c.kind == "group" else width_of(c)
                            for c in s.children), default=0.0)
            return measure(s.children)
        return 0.0

    def lay(stmts: list[Statement], x: float, channel: str) -> float:
        i = 0
        while i < len(stmts):
            s = stmts[i]
            i += 1
            ch = s.channel or channel
            if s.kind == "label":
                labels[s.name] = x
            elif s.kind == "delay":
                if s.name == "aq" and events and events[-1].kind == "acquire" \
                        and abs(events[-1].end - x) < 1e-9:
                    events[-1].value = f"aq {dur(s)}".strip()
                    continue
                w = _delay_width(s)
                mult = s.duration_expr.split("*", 1)[1] if "*" in s.duration_expr else ""
                ev = Event("delay", ch, x, w, label=s.name + (f" ×{mult}" if mult else ""),
                           value=dur(s), cond=s.cond, stmt=s)
                if s.extras:
                    ev.note = " ".join(s.extras)
                events.append(ev)
                x += w
            elif s.kind in ("pulse", "shaped_pulse"):
                w = W_SHAPED if s.kind == "shaped_pulse" else W_PULSE
                mult = s.duration_expr.split("*", 1)[1] if "*" in s.duration_expr else ""
                label = s.name + (f" ×{mult}" if mult else "")
                if s.kind == "shaped_pulse" and s.power:
                    label = f"{s.name}:{s.power}"
                note = ""
                if s.power:
                    pw = values.get(s.power)
                    note = s.power + (f" {pw.text}" if pw is not None and pw.text else "")
                    if s.kind == "shaped_pulse":
                        m = re.match(r"^sp(\d+)$", s.power)
                        shape = values.get(f"spnam{m.group(1)}") if m else None
                        if shape is not None and shape.text and shape.text != "(none)":
                            note = f"{shape.text} · {note}"
                events.append(Event(s.kind, ch or "f1", x, w, label=label, sub=s.phase,
                                    value=dur(s), note=note, cond=s.cond, stmt=s))
                x += w
            elif s.kind == "acquire":
                ev = Event("acquire", "acq", x, W_ACQ, label=s.name, sub=s.phase,
                           value=_value_text(values, "aq") and f"aq {_value_text(values, 'aq')}",
                           cond=s.cond, stmt=s)
                events.append(ev)
                x += W_ACQ
                if s.target_label:
                    loops.append(Loop(labels.get(s.target_label, 0.0), x, "ns",
                                      label=s.target_label))
            elif s.kind == "set_power":
                events.append(Event("set_power", ch or "f1", x, W_POWER, label=s.power,
                                    value=_value_text(values, s.power), cond=s.cond, stmt=s))
                x += W_POWER
            elif s.kind == "decouple_on":
                x0 = x
                # `go=2 ph31 cpd2:f2`: decoupling during the acquisition it rides on
                if events and events[-1].kind == "acquire" and events[-1].stmt is not None \
                        and events[-1].stmt.line_no == s.line_no:
                    x0 = events[-1].start
                open_dec[ch or "f2"] = (x0, s)
            elif s.kind == "decouple_off":
                c = ch or "f2"
                if c in open_dec:
                    x0, st = open_dec.pop(c)
                    events.append(_decouple_event(st, c, x0, max(x - x0, W_SHORT), values))
            elif s.kind == "loop":
                start = labels.get(s.target_label)
                if start is None:
                    notes.append(f"loop to unknown label {s.target_label!r}")
                    start = 0.0
                loops.append(Loop(start, x, s.times_expr or "?", label=s.target_label))
            elif s.kind == "goto":
                notes.append(f"{s.name} {s.condition + ' ' if s.condition else ''}"
                             f"→ {s.target_label}".replace("if  →", "goto →"))
            elif s.kind == "group":
                if s.name in ("sim", "center", "align"):
                    gw = width_of(s)
                    end = x
                    for c in s.children:
                        cw = measure(c.children) if c.kind == "group" else width_of(c)
                        x0 = x + (gw - cw) / 2 if s.name in ("center", "align") else x
                        kids = c.children if c.kind == "group" else [c]
                        end = max(end, lay(kids, x0, c.channel or channel))
                    x = max(end, x + gw)
                else:
                    x = lay(s.children, x, s.channel or channel)
            # label / misc / comment: no width
        return x

    total = lay(program.statements, 0.0, "")
    for c, (x0, st) in open_dec.items():
        events.append(_decouple_event(st, c, x0, max(total - x0, W_SHORT), values))
        notes.append(f"{st.name} on {c} is never switched off (no do:{c})")
    channels = program.channels()
    for ev in events:
        if ev.channel and ev.channel != "acq" and ev.channel not in channels:
            channels.append(ev.channel)
    channels = sorted(channels, key=lambda c: (c == "grad", c)) or ["f1"]
    if any(ev.kind == "acquire" for ev in events):
        channels.append("acq")
    # bracket rows: the shortest span sits nearest the channel rows (level 0);
    # a longer loop that overlaps it takes the next free level (TopSpin's
    # mc loops start at staggered labels, so they overlap without nesting)
    placed: list[Loop] = []
    for lp in sorted(loops, key=lambda o: (o.end - o.start, o.start)):
        used = {o.level for o in placed if o.start < lp.end and lp.start < o.end}
        lvl = 0
        while lvl in used:
            lvl += 1
        lp.level = lvl
        placed.append(lp)
        cnt = values.get(lp.times)
        lp.text = f"× {lp.times}"
        if cnt is not None and cnt.text and cnt.kind == "count" and not _NUMBER.match(lp.times):
            lp.text += f" = {cnt.text}"
    # a delay that an id<N> / dd<N> increments is the t1 (or tau) axis
    incs: dict[str, str] = {}
    for s in program.walk():
        for tok in ([s.name] if s.kind == "misc" else []) + list(s.extras):
            m = re.match(r"^(id|dd)(\d+)$", tok)
            if m:
                incs.setdefault(f"d{m.group(2)}",
                                f"{'+' if m.group(1) == 'id' else '−'}in{m.group(2)} per increment")
            if tok == "ivd":
                incs.setdefault("vd", "next list entry per slice")
    for ev in events:
        if ev.kind == "delay" and ev.stmt is not None and ev.stmt.name in incs:
            ev.note = (f"{ev.note} · " if ev.note else "") + incs[ev.stmt.name]
    for sym in unresolved:
        v = values.get(sym)
        why = f" ({v.note})" if v is not None and v.note else ""
        notes.append(f"{sym}: duration not resolved{why}")
    return Timeline(channels, events, loops, max(total, 1.0), notes, labels)


def _decouple_event(st: Statement, channel: str, x0: float, w: float, values: dict[str, Value]) -> Event:
    m = re.match(r"^cpds?(\d+)$", st.name)
    label = st.name
    value = ""
    if m:
        prg = values.get(f"cpdprg{m.group(1)}")
        if prg is not None and prg.text and prg.text != "(none)":
            label = f"{st.name} · {prg.text}"
        pc = values.get(f"pcpd{m.group(1)}")
        if pc is not None and pc.text:
            value = f"pcpd{m.group(1)} {pc.text}"
    return Event("decouple", channel, x0, w, label=label, value=value, cond=st.cond, stmt=st)


# ---------------------------------------------------------------- describe
def _elem_text(s: Statement, values: dict[str, Value], channel: str = "",
               in_group: bool = False) -> str:
    ch = s.channel or channel
    suffix = "" if in_group else (f" {ch}" if ch and ch not in ("", "f1") else "")
    if s.kind == "delay":
        if s.name in _MC_DELAYS or _LITERAL.match(s.name):
            return ""
        return s.name
    if s.kind == "pulse":
        return f"{s.name}" + (f" ({s.phase})" if s.phase else "") + suffix
    if s.kind == "shaped_pulse":
        return f"{s.name}:{s.power} shaped" + (f" ({s.phase})" if s.phase else "") + suffix
    if s.kind == "decouple_on":
        return f"{s.name} on {s.channel}"
    if s.kind == "decouple_off":
        return f"do {s.channel}"
    if s.kind == "acquire":
        return f"acquire {s.name}" + (f" ({s.phase})" if s.phase else "")
    if s.kind == "group":
        if s.name in ("sim", "center", "align"):
            parts = []
            for c in s.children:
                kids = c.children if c.kind == "group" else [c]
                inner = [t for t in (_elem_text(k, values, c.channel, True) for k in kids) if t]
                parts.append(f"{c.channel}: " + " – ".join(inner) if c.channel else " – ".join(inner))
            return ("centre [" if s.name in ("center", "align") else "[") + " | ".join(parts) + "]"
        return " – ".join(t for t in (_elem_text(k, values, s.channel) for k in s.children) if t)
    return ""


#: loop counts that belong to the acquisition machinery, not to the sequence
_MC_TIMES = frozenset({"td0", "td1", "td2", "td3", "tdav", "ns"})


def describe(program: Program, values: dict[str, Value] | None = None) -> str:
    """``zg: d1 – p1 (ph1) – acquire go=2 (ph31); phase cycle of 8 steps``;
    a local loop is folded into brackets (``[d20 – p1 (ph4)] ×l20``), the
    2D / averaging loops TopSpin's mc macro adds are listed afterwards."""
    values = values or {}
    parts: list[str] = []
    marks: dict[str, int] = {}
    folded: set[int] = set()
    for s in program.statements:
        if s.kind == "label":
            marks[s.name] = len(parts)
        elif s.kind == "loop" and s.name == "lo":
            pos = marks.get(s.target_label)
            if pos is not None and pos < len(parts) and not s.target_label.startswith("LBL") \
                    and s.times_expr not in _MC_TIMES:
                inner = " – ".join(parts[pos:])
                del parts[pos:]
                parts.append(f"[{inner}] ×{s.times_expr}")
                folded.add(id(s))
        else:
            t = _elem_text(s, values)
            if t:
                parts.append(t)
    out = f"{program.name or 'pulse program'}: " + " – ".join(parts)
    loops = [s for s in program.statements if s.kind == "loop" and id(s) not in folded
             and s.times_expr not in ("td0",) and s.name != "rcyc"]
    if loops:
        out += "; loops: " + ", ".join(
            f"×{s.times_expr} ({s.target_label})" if s.target_label else f"×{s.times_expr}"
            for s in loops)
    n = program.phase_cycle_length()
    if n:
        out += f"; phase cycle of {n} steps"
    return out


# -------------------------------------------------------------------- load
def channel_nuclei(acqus: dict | None) -> dict[str, str]:
    """{'f1': '27Al', 'f2': '1H'} from NUC1..NUC4 (``off`` dropped)."""
    out = {}
    for i in range(1, 5):
        nuc = str((acqus or {}).get(f"NUC{i}", "") or "").strip("<> ")
        if nuc and nuc.lower() != "off":
            out[f"f{i}"] = nuc
    return out


def load(expno) -> Record:
    """Everything the viewer needs for an EXPNO: the file, the parsed
    program, the resolved values, the channel nuclei. ``program`` is None
    when there is no pulse program file; never raises."""
    from larmor.referencing import _jcamp

    p = Path(str(expno))
    if p.is_file():
        p = p.parent
    for cand in [p, *list(p.parents)[:4]]:
        if (cand / "acqus").is_file():
            p = cand
            break
    acqus = _jcamp(p / "acqus") or {}
    acqu2s = _jcamp(p / "acqu2s") if (p / "acqu2s").is_file() else None
    pulprog = str(acqus.get("PULPROG", "") or "").strip("<> ")
    path = find_pulseprogram(p)
    warnings: list[str] = []
    if path is None:
        return Record(p, None, None, {}, channel_nuclei(acqus), acqus, pulprog,
                      ["no pulse program file in this EXPNO"])
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        return Record(p, path, None, {}, channel_nuclei(acqus), acqus, pulprog,
                      [f"could not read {path.name}: {exc}"])
    program = parse(text, name=pulprog)
    warnings.extend(program.warnings)
    try:
        values = resolve(program, acqus, acqu2s)
    except Exception as exc:                                  # noqa: BLE001
        values = {}
        warnings.append(f"values not resolved: {exc}")
    return Record(p, path, program, values, channel_nuclei(acqus), acqus, pulprog, warnings)
