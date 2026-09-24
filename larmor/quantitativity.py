"""Can the populations of a fitted spectrum be trusted? The acquisition side.

``quantify.py`` measures how much of each line the integration window cuts
off. This module judges the two acquisition conditions the glass-fitting
protocol (help/glass-fitting.md section 6) asks for before an integral is
read as a population:

* **recycle delay** — D1 + AQ of the fitted EXPNO against the T1 of each
  site, taken from TopSpin's ``pdata/1/ct1t2.txt`` of the same-nucleus
  saturation/inversion-recovery EXPNO in the same sample folder
  (``scan.find_sibling_t1`` + ``satrec.read_ct1t2``), matched to the site by
  its δiso falling inside the TopSpin integral region (nearest region
  otherwise, marked). Recovery uses the steady-state formula
  (1 − E)/(1 − E·cos θ_eff), E = exp(−(D1+AQ)/T1), with θ_eff = (I+½)·θ for
  half-integer quadrupolar nuclei (CT-selective limit, capped at 180°) and
  θ = 90° when no flip angle is known or the program is not a single pulse
  (both conservative);
* **flip angle** — 90·P1/P1(90) from acqus and the title's ``P1(90)=`` (or
  the 90° pulse typed in Experiment parameters), else the title's stated
  ``n deg tip``, against the short-pulse limit 30°/(I+½) (Edén Eq. 38) for
  single-pulse programs on half-integer quadrupolar nuclei.

Everything here reads files and returns plain scalars (``Check``); every
reader fails soft to None / [] and nothing is ever written into an
instrument folder. ``fithealth`` turns a ``Check`` into chips — this module
must not import it. The 90 % / 30° thresholds live here; ``fithealth``
imports them (its 'no new judgement' rule).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

#: steady-state recovery below which the recycle-delay chip is amber
#: (≈ D1 ≥ 4.6 T1 at 90°; the group's '5*T1' rule gives 99.3 %)
RECOVERY_MIN = 0.99
#: the short-pulse linear regime: flip ≤ 30°/(I+½) (Edén Eq. 38)
FLIP_LIMIT_DEG = 30.0
#: a TopSpin T1 longer than this many times the longest vdlist delay is a
#: diverged fit (Base2Ca/12: 3.5e6 s against a 256 s delay), not a T1
T1_IMPLAUSIBLE_FACTOR = 10.0
#: scan._classify kinds whose excitation uniformity the flip rule judges
#: (the T1 kinds the sibling search accepts are scan.T1_KINDS)
SINGLE_PULSE_KINDS = ("Single pulse", "Single pulse (dec.)", "Bloch decay")
#: the two background models: no δiso meaning, skipped for recovery
_BACKGROUND_MODELS = frozenset({"spectrum", "function"})
#: echo-family kinds named as such in the passing line
_ECHO_KINDS = ("Hahn echo", "Spin echo", "Echo", "CPMG (T2)")
#: fithealth.TEXT_LIMIT (chip wording bound); kept literal so that fithealth
#: can import this module without a cycle
_TEXT_LIMIT = 60

TITLE_P90_RE = re.compile(r"P1\(90\)\s*=\s*([0-9.]+)", re.IGNORECASE)
#: requires the tip/pulse/flip word, so 'VT 25 degrees' never matches
TITLE_FLIP_RE = re.compile(
    r"([0-9.]+)\s*(?:deg(?:ree)?s?|°)\s*(?:tip|pulse|flip)", re.IGNORECASE)
TITLE_T1_MULTIPLE_RE = re.compile(r"([0-9.]+)\s*[x×*]\s*T1", re.IGNORECASE)
#: operator T1 notes in a relaxation title: 'satrec 4.64 3.74s',
#: 'T1 approximately 53s', 'Longest T1 at 27.6 ms', 'Satrec 993m'
T1_NOTE_RE = re.compile(
    r"(?:T1|satrec)\b[^0-9]{0,40}?([0-9.]+(?:\s+[0-9.]+)*)\s*(ms|m|s)\b",
    re.IGNORECASE)
_UNIT_S = {"s": 1.0, "m": 1e-3, "ms": 1e-3}


# ------------------------------------------------------------------ facts
@dataclass
class Acquisition:
    """What acqus and the title say about the fitted EXPNO."""

    expno: str
    nucleus: str
    pulprog: str
    kind: str
    ns: int | None
    d1_s: float | None
    aq_s: float | None
    p1_us: float | None
    plw1_w: float | None
    probhd: str
    title: str
    p90_us_title: float | None
    flip_deg_title: float | None
    t1_multiple_claimed: float | None

    @property
    def recycle_s(self) -> float:
        return (self.d1_s or 0.0) + (self.aq_s or 0.0)

    @property
    def single_pulse(self) -> bool:
        return self.kind in SINGLE_PULSE_KINDS


@dataclass
class T1Source:
    """Where the T1 per site comes from: TopSpin's ct1t2 of a sibling EXPNO
    ('ct1t2'), a typed value ('user') or LARMOR's per-site analysis
    ('per_site')."""

    expno: str
    kind: str
    regions: list = field(default_factory=list)      # satrec.T1Region
    vdlist_max_s: float | None = None
    plausible: bool = True
    note: str = ""
    title_notes_s: list = field(default_factory=list)

    @property
    def name(self) -> str:
        return Path(self.expno).name if self.expno else ""

    def t1_for_ppm(self, delta_ppm: float):
        """(t1_s, region, exact): the region holding ``delta_ppm`` (exact
        True), else the nearest region by edge distance (exact False);
        (None, None, False) without regions."""
        if not self.regions:
            return None, None, False
        for r in self.regions:
            if r.lo_ppm <= delta_ppm <= r.hi_ppm:
                return r.t1_s, r, True
        near = min(self.regions, key=lambda r: min(abs(delta_ppm - r.lo_ppm),
                                                   abs(delta_ppm - r.hi_ppm)))
        return near.t1_s, near, False


@dataclass
class SiteRecovery:
    index: int
    label: str
    ppm: float
    t1_s: float
    exact: bool
    ratio: float
    recovery: float
    #: 'ct1t2' (TopSpin region) | 'user' (typed T1) | 'per_site' (typed or
    #: measured per site)
    source: str = "ct1t2"


@dataclass
class AcqFacts:
    """Everything read from disk for one fitted spectrum (cached per source
    path by the desktop)."""

    acquisition: Acquisition | None
    t1: T1Source | None
    #: the nearest same-nucleus T1 EXPNO path (click-through even without a fit)
    sibling_expno: str | None
    spin: float | None
    folder: str
    #: 'ok' | 'missing' (sibling without a usable ct1t2) | 'implausible'
    #: (TopSpin fit diverged) | 'none' (no T1 EXPNO of this nucleus)
    t1_status: str = "none"
    t1_note: str = ""


# ------------------------------------------------------------- arithmetic
def spin_of(nucleus: str) -> float | None:
    """Nuclear spin of '27Al' → 2.5; None for an unknown symbol."""
    if not nucleus:
        return None
    try:
        from larmor.nuclei import all_isotopes

        return next((float(i.spin) for i in all_isotopes()
                     if i.symbol == nucleus), None)
    except Exception:
        return None


def _half_integer_quadrupolar(spin) -> bool:
    return (spin is not None and spin > 0.5
            and abs(spin * 2 - round(spin * 2)) < 1e-9
            and int(round(spin * 2)) % 2 == 1)


def flip_limit_deg(spin) -> float | None:
    """30°/(I+½) for a half-integer quadrupolar nucleus; None for spin-½,
    integer spins or an unknown nucleus."""
    if not _half_integer_quadrupolar(spin):
        return None
    return FLIP_LIMIT_DEG / (spin + 0.5)


def effective_flip_deg(flip_deg: float, spin) -> float:
    """The central-transition nutation angle (I+½)·θ of a CT-selective pulse
    on a half-integer quadrupolar nucleus (θ otherwise), capped at 180°."""
    theta = float(flip_deg)
    if _half_integer_quadrupolar(spin):
        theta *= (spin + 0.5)
    return min(theta, 180.0)


def steady_state_recovery(recycle_s: float, t1_s: float, flip_deg=None,
                          spin=None) -> float:
    """Steady-state magnetisation before each pulse, as a fraction of the
    equilibrium value: (1 − E)/(1 − E·cos θ_eff), E = exp(−recycle/T1).
    ``flip_deg`` None (or 90°) gives the saturation-recovery form 1 − E."""
    if t1_s is None or t1_s <= 0:
        return 1.0
    e = math.exp(-max(float(recycle_s), 0.0) / float(t1_s))
    if flip_deg is None:
        return 1.0 - e
    c = math.cos(math.radians(effective_flip_deg(flip_deg, spin)))
    denom = 1.0 - e * c
    return (1.0 - e) / denom if denom > 0 else 1.0


def spin_text(spin) -> str:
    """1.5 → '3/2', 0.5 → '1/2', 1.0 → '1'."""
    if spin is None:
        return "?"
    n = int(round(spin * 2))
    return f"{n // 2}" if n % 2 == 0 else f"{n}/2"


def _deg(v: float) -> str:
    """An angle for a chip: '90', '10.8', '7.5' (one decimal unless whole)."""
    return f"{v:.0f}" if abs(v - round(v)) < 0.05 else f"{v:.1f}"


# ----------------------------------------------------------------- titles
def _float(tok) -> float | None:
    try:
        v = float(tok)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def parse_title(title: str) -> dict:
    """The operator's claims in a TopSpin title: {'p90_us', 'flip_deg',
    't1_multiple', 't1_notes_s'} (None / [] when absent)."""
    text = title or ""
    out = {"p90_us": None, "flip_deg": None, "t1_multiple": None,
           "t1_notes_s": []}
    m = TITLE_P90_RE.search(text)
    if m:
        out["p90_us"] = _float(m.group(1))
    m = TITLE_FLIP_RE.search(text)
    if m:
        out["flip_deg"] = _float(m.group(1))
    m = TITLE_T1_MULTIPLE_RE.search(text)
    if m:
        out["t1_multiple"] = _float(m.group(1))
    for m in T1_NOTE_RE.finditer(text):
        mult = _UNIT_S.get(m.group(2).lower(), 1.0)
        for tok in m.group(1).split():
            v = _float(tok)
            if v is not None:
                out["t1_notes_s"].append(v * mult)
    return out


# ---------------------------------------------------------------- readers
def read_acquisition(expno_dir) -> Acquisition | None:
    """acqus + title of one EXPNO as an ``Acquisition``; None when the folder
    has no readable acqus."""
    from larmor.io import bruker, scan

    p = Path(expno_dir)
    meta = bruker.read_acqus_meta(p)
    if not meta:
        return None
    title = str(meta.get("title") or "")
    claims = parse_title(title)
    p_us = meta.get("p_us") or []
    plw = meta.get("plw_w") or []
    pulprog = str(meta.get("pulse_program") or "")
    return Acquisition(
        expno=str(p), nucleus=str(meta.get("nucleus") or ""),
        pulprog=pulprog, kind=scan._classify(pulprog),
        ns=meta.get("ns"), d1_s=meta.get("d1_s"), aq_s=meta.get("aq_s"),
        p1_us=(p_us[1] if len(p_us) > 1 else None),
        plw1_w=(plw[1] if len(plw) > 1 else None),
        probhd=str(meta.get("probhd") or ""), title=title,
        p90_us_title=claims["p90_us"], flip_deg_title=claims["flip_deg"],
        t1_multiple_claimed=claims["t1_multiple"])


def read_t1_source(expno_dir) -> T1Source | None:
    """TopSpin's per-region T1 of a relaxation EXPNO (pdata/1/ct1t2.txt) with
    a plausibility check against the vdlist; regions [] when the file is
    missing or carries no T1 line (the caller reads that as 'missing').
    None only when the folder does not exist."""
    from larmor import satrec
    from larmor.io import bruker

    p = Path(expno_dir)
    if not p.is_dir():
        return None
    regions = satrec.read_ct1t2(p / "pdata" / "1")
    vdmax = None
    try:
        vd = satrec.read_vdlist(p)
        vdmax = float(vd.max()) if vd.size else None
    except Exception:
        vdmax = None
    try:
        title = bruker._read_title(p / "pdata" / "1")
    except Exception:
        title = ""
    notes = parse_title(title)["t1_notes_s"]
    src = T1Source(expno=str(p), kind="ct1t2", regions=list(regions),
                   vdlist_max_s=vdmax, plausible=True, note="",
                   title_notes_s=notes)
    if not regions:
        src.note = f"EXPNO {p.name}: no T1 result in pdata/1/ct1t2.txt"
        return src
    bad = [r for r in regions if r.t1_s <= 0
           or (vdmax and r.t1_s > T1_IMPLAUSIBLE_FACTOR * vdmax)]
    if bad:
        src.plausible = False
        worst = max(bad, key=lambda r: r.t1_s)
        vd_txt = f" vs a {vdmax:g} s longest delay" if vdmax else ""
        src.note = (f"EXPNO {p.name}: TopSpin fit did not converge "
                    f"(T1 {worst.t1_s:.2g} s{vd_txt})")
    return src


def facts_for(source_path) -> AcqFacts | None:
    """The acquisition facts behind any Bruker source path (an EXPNO folder,
    a pdata folder, a 1r/fid file); None for anything ``bruker.resolve``
    cannot handle (a synthetic 'src', an fxmla or CSV file)."""
    if not source_path:
        return None
    from larmor.io import bruker, scan

    try:
        ref = bruker.resolve(source_path)
        expno = Path(ref.expno)
    except Exception:
        return None
    acq = read_acquisition(expno)
    if acq is None:
        return None
    facts = AcqFacts(acquisition=acq, t1=None, sibling_expno=None,
                     spin=spin_of(acq.nucleus), folder=str(expno.parent))
    try:
        sibs = scan.find_sibling_t1(expno, acq.nucleus)
    except Exception:
        sibs = []
    if not sibs:
        facts.t1_note = (f"no T1 EXPNO of {acq.nucleus or 'this nucleus'} "
                         f"in {expno.parent.name or expno.parent}")
        return facts

    def _dist(info):
        try:
            return abs(int(info.expno) - int(expno.name))
        except ValueError:
            return float("inf")

    facts.sibling_expno = min(sibs, key=_dist).path
    facts.t1_status = "missing"
    notes = []
    for info in sibs:
        src = read_t1_source(info.path)
        if src is None:
            continue
        if src.regions and src.plausible:
            facts.t1, facts.t1_status, facts.t1_note = src, "ok", ""
            return facts
        if not src.plausible and facts.t1_status != "implausible":
            facts.t1_status = "implausible"
        if src.note:
            notes.append(src.note)
    facts.t1_note = "; ".join(notes) or "no usable T1 result in the sample folder"
    return facts


# ------------------------------------------------------------------ check
def _site_view(sites) -> list:
    """(index, label, model, δiso) of Recipe.sites or site dicts."""
    out = []
    for i, s in enumerate(sites or []):
        if isinstance(s, dict):
            model = str(s.get("model", ""))
            label = str(s.get("label") or model)
            p = (s.get("params") or {}).get("isotropic_chemical_shift_ppm")
            v = p.get("value") if isinstance(p, dict) else getattr(p, "value", p)
        else:
            model = str(getattr(s, "model", ""))
            label = str(getattr(s, "label", "") or model)
            p = getattr(s, "params", {}).get("isotropic_chemical_shift_ppm")
            v = getattr(p, "value", None)
        out.append((i, label, model, _float(v)))
    return out


@dataclass
class Check:
    """The judgement fithealth renders: plain scalars, strings and lists only
    (``Health.__eq__`` compares it)."""

    facts: AcqFacts
    flip_deg: float | None
    flip_source: str                   # 'P1(90)' | 'title' | 'user' | 'echo(90)' | ''
    flip_assumed_90: bool
    flip_limit_deg: float | None
    flip_effective_deg: float | None
    excitation_judged: bool
    recoveries: list                   # SiteRecovery
    recovery_min: float | None
    recovery_max: float | None
    recovery_spread: float | None
    t1_status: str                     # 'ok' | 'missing' | 'implausible' | 'none'

    # ------------------------------------------------------------ helpers
    @property
    def acquisition(self) -> Acquisition | None:
        return self.facts.acquisition if self.facts is not None else None

    @property
    def spin(self):
        return self.facts.spin if self.facts is not None else None

    def _ratio_text(self) -> str:
        rs = [r.ratio for r in self.recoveries]
        lo, hi = f"{min(rs):.1f}", f"{max(rs):.1f}"
        return lo if lo == hi else f"{lo}–{hi}"

    def _recovery_pct_text(self) -> str:
        lo = f"{100 * self.recovery_min:.0f}"
        hi = f"{100 * self.recovery_max:.0f}"
        return lo if lo == hi else f"{lo}–{hi}"

    def _angle_text(self) -> str:
        if self.flip_assumed_90 or self.flip_deg is None:
            return ""
        return f" at {_deg(self.flip_deg)}°"

    def _flip_source_text(self) -> str:
        acq = self.acquisition
        if self.flip_source == "P1(90)":
            return (f"P1(90)={acq.p90_us_title:g} in the title"
                    if acq and acq.p90_us_title else "P1(90) in the title")
        if self.flip_source == "title":
            return "the tip angle stated in the title"
        if self.flip_source == "user":
            return "typed in Experiment parameters"
        if self.flip_source == "echo(90)":
            return "a 90° excitation assumed for this pulse program"
        return "unknown"

    def recovery_ok(self) -> bool:
        return self.recovery_min is not None and self.recovery_min >= RECOVERY_MIN

    def excitation_over(self) -> bool:
        return (self.excitation_judged and self.flip_deg is not None
                and self.flip_limit_deg is not None
                and self.flip_deg > self.flip_limit_deg)

    # -------------------------------------------------------------- texts
    def recovery_text(self) -> str | None:
        """'D1 = 3.0–3.6 T1 → 95–97 % (90° assumed)' / 'D1 = 2.1 T1 → 88 %'
        / 'D1 = 2.1–4.7 T1 at 90° → 88–99 %'; None without recoveries."""
        if not self.recoveries:
            return None
        text = (f"D1 = {self._ratio_text()} T1{self._angle_text()} → "
                f"{self._recovery_pct_text()} %")
        if self.flip_assumed_90:
            text += " (90° assumed)"
        return text[:_TEXT_LIMIT]

    def recovery_detail(self) -> str:
        acq = self.acquisition
        parts = []
        if acq is not None:
            d1 = acq.d1_s if acq.d1_s is not None else 0.0
            aq = acq.aq_s if acq.aq_s is not None else 0.0
            parts.append(f"recycle D1 + AQ = {acq.recycle_s:.4g} s "
                         f"({d1:g} s + {aq:.3g} s)")
        src = self.facts.t1 if self.facts is not None else None
        if src is not None and src.kind == "ct1t2":
            parts.append(f"T1 per integral region from EXPNO {src.name} "
                         "(TopSpin ct1t2.txt)")
        elif self.recoveries:
            parts.append("T1 typed in Experiment parameters / measured per site")
        for r in self.recoveries:
            reg = ""
            if r.source == "ct1t2" and src is not None and src.regions:
                _, region, _ = src.t1_for_ppm(r.ppm)
                if region is not None:
                    reg = (f" ({region.hi_ppm:.1f} … {region.lo_ppm:.1f} ppm region"
                           + ("" if r.exact else
                              " — nearest; δiso outside every TopSpin integral")
                           + ")")
            elif r.source != "ct1t2":
                reg = " (typed / per-site T1)"
            parts.append(f"{r.label} at {r.ppm:g} ppm{reg}: T1 {r.t1_s:.3g} s, "
                         f"{r.ratio:.1f} T1, {100 * r.recovery:.1f} %")
        if self.recovery_spread is not None and len(self.recoveries) > 1 \
                and self.recovery_max:
            bias = 100.0 * (1.0 - self.recovery_min / self.recovery_max)
            parts.append(f"site ratios biased by up to {bias:.1f} %")
        theta = self.flip_effective_deg
        if self.flip_assumed_90:
            angle = ("θ = 90° assumed — enter the 90° pulse in Process ▸ "
                     "Experiment parameters… for the steady-state value of "
                     "the short pulse")
        elif self.flip_source == "echo(90)":
            angle = f"θ = 90° ({self._flip_source_text()})"
        elif theta is not None and self.flip_deg is not None \
                and abs(theta - self.flip_deg) > 1e-9:
            angle = (f"θ_eff = (I+½)·{_deg(self.flip_deg)}° = {_deg(theta)}° "
                     f"({self._flip_source_text()})")
        elif theta is not None:
            angle = f"θ = {_deg(theta)}° ({self._flip_source_text()})"
        else:
            angle = "θ = 90°"
        parts.append("steady state (1−E)/(1−E·cos θ), E = exp(−(D1+AQ)/T1); "
                     + angle)
        if acq is not None and acq.t1_multiple_claimed and self.recoveries:
            parts.append(f"title says {acq.t1_multiple_claimed:g}×T1; measured "
                         f"{self._ratio_text()} T1")
        if src is not None and src.title_notes_s:
            vals = " ".join(f"{v:g}" for v in src.title_notes_s)
            parts.append(f"note in EXPNO {src.name}: T1 {vals} s")
        parts.append(f"limit {100 * RECOVERY_MIN:.0f} % (≈ 4.6 T1 at 90°); "
                     "click for Tools ▸ Relaxation on that EXPNO")
        return "; ".join(parts)

    def recovery_unknown_text(self) -> str:
        acq = self.acquisition
        val = acq.d1_s if (acq is not None and acq.d1_s is not None) else (
            acq.recycle_s if acq is not None else 0.0)
        return f"recycle {val:g} s — T1 unknown"

    def recovery_unknown_detail(self) -> str:
        acq = self.acquisition
        parts = []
        if acq is not None:
            parts.append(f"D1 = {acq.d1_s if acq.d1_s is not None else 0:g} s "
                         f"of EXPNO {Path(acq.expno).name}")
        note = self.facts.t1_note if self.facts is not None else ""
        if self.t1_status == "none":
            parts.append(note or "no T1 EXPNO of this nucleus in the sample folder")
        elif self.t1_status == "implausible":
            parts.append(note or "the TopSpin T1 fit did not converge")
        else:
            parts.append(note or "the sibling T1 EXPNO carries no TopSpin "
                         "result (pdata/1/ct1t2.txt)")
        if self.facts is not None and self.facts.sibling_expno:
            parts.append("click for Tools ▸ Relaxation on that EXPNO, or F7 ▸ "
                         "'Measure T1 per site' with this fit")
        else:
            parts.append("click for Tools ▸ Relaxation, or type a T1 in Process ▸ "
                         "Experiment parameters…")
        return "; ".join(parts)

    def excitation_text(self) -> str | None:
        """'flip 18° > 15° limit (I = 3/2)'; None unless the flip is over."""
        if not self.excitation_over():
            return None
        return (f"flip {_deg(self.flip_deg)}° > {_deg(self.flip_limit_deg)}° limit "
                f"(I = {spin_text(self.spin)})")[:_TEXT_LIMIT]

    def excitation_detail(self) -> str:
        acq = self.acquisition
        parts = []
        if acq is not None and acq.p1_us is not None:
            power = f" at {acq.plw1_w:g} W" if acq.plw1_w else ""
            parts.append(f"P1 {acq.p1_us:g} µs{power}")
        if self.flip_deg is not None:
            parts.append(f"flip {_deg(self.flip_deg)}° from {self._flip_source_text()}")
        if self.flip_limit_deg is not None and self.flip_deg is not None:
            excess = (self.flip_deg - self.flip_limit_deg) / self.flip_limit_deg
            if 0 < excess < 0.25:
                parts.append(f"just above the limit ({100 * excess:.0f} %)")
        if self.flip_limit_deg is not None:
            parts.append(f"linear regime τ ≤ 1/[12(I+½)ν₁] ⇔ θ ≤ 30°/(I+½) = "
                         f"{_deg(self.flip_limit_deg)}° for I = {spin_text(self.spin)} "
                         "(Edén Eq. 38)")
        else:
            parts.append("spin-½: any flip angle is quantitative")
        parts.append("sites with different C_Q were excited with different "
                     "efficiency — no fit repairs this; state it in the Methods")
        parts.append("click to enter the 90° pulse in Process ▸ Experiment "
                     "parameters…")
        return "; ".join(parts)

    def excitation_unknown_text(self) -> str:
        return f"flip angle unknown (I = {spin_text(self.spin)})"

    def excitation_unknown_detail(self) -> str:
        acq = self.acquisition
        parts = []
        if acq is not None:
            parts.append(f"no P1(90)= or 'n deg tip' in the title of EXPNO "
                         f"{Path(acq.expno).name}")
            if acq.p1_us is not None:
                power = f" at {acq.plw1_w:g} W" if acq.plw1_w else ""
                parts.append(f"P1 = {acq.p1_us:g} µs{power}")
        if self.flip_limit_deg is not None:
            parts.append(f"limit 30°/(I+½) = {_deg(self.flip_limit_deg)}° for "
                         f"I = {spin_text(self.spin)}")
        parts.append("click to type the 90° pulse at this power in Process ▸ "
                     "Experiment parameters… (remembered per nucleus, probe "
                     "and power)")
        return "; ".join(parts)

    def unchecked_lines(self) -> list:
        """What could NOT be judged — no usable T1, no flip angle — as neutral
        tooltip lines (never chips: an unknown fact is not a problem). Each
        names where the fact can be supplied."""
        out = []
        acq = self.acquisition
        if acq is None:
            return out
        if not self.recoveries and self.t1_status in ("missing", "implausible",
                                                      "none"):
            note = self.facts.t1_note if self.facts is not None else ""
            if self.t1_status == "none":
                why = note or "no T1 EXPNO of this nucleus in the sample folder"
            elif self.t1_status == "implausible":
                why = note or "the TopSpin T1 fit did not converge"
            else:
                why = note or "the sibling T1 EXPNO carries no TopSpin result"
            out.append(f"{self.recovery_unknown_text()} — not judged ({why}); "
                       "Tools ▸ Relaxation measures it, or type a T1 in Process "
                       "▸ Experiment parameters…")
        if self.excitation_judged and self.flip_deg is None:
            out.append(f"{self.excitation_unknown_text()} — not judged; the 90° "
                       "pulse can be typed in Process ▸ Experiment parameters…")
        return out

    def passing_lines(self, exclude=()) -> list:
        """What passed, for the pill tooltip; ``exclude`` names the flag kinds
        ('recovery' / 'excitation') already shown as chips."""
        out = []
        acq = self.acquisition
        if acq is None:
            return out
        if "recovery" not in exclude:
            if self.recoveries and self.recovery_ok():
                src = self.facts.t1
                where = f" (EXPNO {src.name})" if (src is not None and src.name) else ""
                out.append(f"recycle {acq.recycle_s:.3g} s = {self._ratio_text()} "
                           f"T1 → {self._recovery_pct_text()} %{where}")
            elif not self.recoveries and self.t1_status == "ok":
                out.append(f"recycle {acq.recycle_s:.3g} s — T1 not checked")
        if "excitation" not in exclude:
            if self.spin is not None and self.spin == 0.5:
                out.append("spin-½: any flip angle is quantitative")
            elif not acq.single_pulse:
                what = "echo/CPMG" if acq.kind in _ECHO_KINDS else acq.kind
                out.append(f"{what}: excitation uniformity not judged")
            elif (self.excitation_judged and self.flip_deg is not None
                  and not self.excitation_over()):
                out.append(f"flip {_deg(self.flip_deg)}° within the linear regime "
                           f"(≤ {_deg(self.flip_limit_deg)}°)")
        return out


def _override_t1(override: dict | None, label: str):
    if not override:
        return None
    by_site = override.get("t1_by_site") or {}
    v = _float(by_site.get(label)) if isinstance(by_site, dict) else None
    if v is not None and v > 0:
        return v, "per_site"
    v = _float(override.get("t1_s"))
    if v is not None and v > 0:
        return v, "user"
    return None


def check(facts: AcqFacts, sites, override: dict | None = None) -> Check:
    """Judge the acquisition against the fitted sites.

    ``sites``: Recipe.sites or site dicts; ``override``: the recipe's
    provenance['quantitativity'] ({'p90_us', 'flip_deg', 't1_s',
    't1_by_site': {label: s}}). Flip precedence: typed flip > 90·P1/P1(90)
    (typed 90° pulse before the title's) > the title's stated tip > unknown;
    a non-single-pulse program counts as a 90° excitation ('echo(90)'). T1
    precedence: t1_by_site > t1_s > the sibling's ct1t2.
    """
    override = override or {}
    acq = facts.acquisition if facts is not None else None
    spin = facts.spin if facts is not None else None
    limit = flip_limit_deg(spin)

    flip = None
    source = ""
    assumed = False
    # the nominal short pulse the (I+½) scaling applies to; None keeps the
    # saturation-recovery form 1 − E (a 90° central-transition excitation:
    # unknown angle, or an echo / other program)
    nominal = None
    if acq is not None:
        if not acq.single_pulse:
            flip, source = 90.0, "echo(90)"
        else:
            user_flip = _float(override.get("flip_deg"))
            user_p90 = _float(override.get("p90_us"))
            if user_flip is not None and user_flip > 0:
                flip, source = user_flip, "user"
            elif acq.p1_us and user_p90 and user_p90 > 0:
                flip, source = 90.0 * acq.p1_us / user_p90, "user"
            elif acq.p1_us and acq.p90_us_title and acq.p90_us_title > 0:
                flip, source = 90.0 * acq.p1_us / acq.p90_us_title, "P1(90)"
            elif acq.flip_deg_title is not None and acq.flip_deg_title > 0:
                flip, source = acq.flip_deg_title, "title"
            else:
                assumed = True
            nominal = flip
    theta = None
    if acq is not None:
        theta = (effective_flip_deg(nominal, spin) if nominal is not None
                 else 90.0)
    judged = bool(acq is not None and acq.single_pulse and limit is not None)

    recoveries = []
    t1_status = facts.t1_status if facts is not None else "none"
    if acq is not None and acq.recycle_s > 0:
        t1src = facts.t1
        for i, label, model, ppm in _site_view(sites):
            if model in _BACKGROUND_MODELS:
                continue
            ov = _override_t1(override, label)
            if ov is not None:
                t1_s, exact, src_kind = ov[0], True, ov[1]
            elif t1src is not None and t1src.regions and ppm is not None:
                t1_s, _region, exact = t1src.t1_for_ppm(ppm)
                src_kind = "ct1t2"
            else:
                continue
            if not t1_s or t1_s <= 0:
                continue
            ratio = acq.recycle_s / t1_s
            rec = steady_state_recovery(acq.recycle_s, t1_s, nominal, spin)
            recoveries.append(SiteRecovery(index=i, label=label,
                                           ppm=ppm if ppm is not None else 0.0,
                                           t1_s=float(t1_s), exact=bool(exact),
                                           ratio=float(ratio),
                                           recovery=float(rec), source=src_kind))
    if recoveries:
        t1_status = "ok"
        rmin = min(r.recovery for r in recoveries)
        rmax = max(r.recovery for r in recoveries)
        spread = rmax - rmin
    else:
        rmin = rmax = spread = None

    return Check(facts=facts, flip_deg=(float(flip) if flip is not None else None),
                 flip_source=source, flip_assumed_90=assumed,
                 flip_limit_deg=limit, flip_effective_deg=theta,
                 excitation_judged=judged, recoveries=recoveries,
                 recovery_min=rmin, recovery_max=rmax, recovery_spread=spread,
                 t1_status=t1_status)
