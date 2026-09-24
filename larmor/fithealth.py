"""One verdict for a fit: sanity, identifiability and residual diagnostics
gathered into a single Qt-free ``Health`` object.

``sanity.py`` (physical values), ``identifiability.py`` (degenerate parameter
pairs) and ``diagnostics.py`` (residual structure) each judge one aspect of a
fit. This module runs them together, adds the residual/noise ratio and the
population rule of the glass-fitting protocol, and turns the outcome into
``Flag`` objects with fixed wording, a click-through target and the table
cells they concern. The desktop strip (``larmor/desktop/fithealth_strip.py``),
the Report header and the status bar all render a ``Health``; nothing here
imports Qt, so the same verdict is available headless.

Levels are classes of flag, not grades of the fit. ``bad`` ("Fit: not
physical", red) is reserved for values that cannot be physical at all -- a
negative amplitude, η outside [0, 1], a width of zero, a C_Q above the
registry's maximum (the sanity rules). Everything a physically sound fit can
still show -- a structured residual, a parameter at a bound, a degenerate
pair, a population the data do not support, δiso outside the fit window, a
too-short recycle delay or a too-long pulse -- is ``check`` ("Fit: check",
amber): something to look at, never an alarm. ``ok`` ("Fit OK") is a fitted
model with no flag and ``none`` a model that has not been fitted yet. An
UNKNOWN acquisition fact (no sibling T1, no 90° pulse) raises no flag: only a
measured problem does, and the unknowns stay in the pill tooltip and in
Process ▸ Experiment parameters. Every threshold is the one its source module
already applies (1.5x edge noise, |z| > 3 / lag-1 > 0.4, |r| >= 0.95, a
population error of 100 %, the 2 % tail, 99 % recovery and 30°/(I+½) rules of
quantify.py / quantitativity.py); no new judgement is introduced.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, replace

import numpy as np

from larmor import diagnostics, quantitativity, sanity
from larmor.identifiability import IDENTIFIABILITY_THRESHOLD, unidentifiable_pairs
from larmor.quantify import TAIL_LIMIT_PCT
from larmor.quantitativity import RECOVERY_MIN
from larmor.recipe import Recipe

LEVELS = ("none", "ok", "check", "bad")
LEVEL_RANK = {lvl: i for i, lvl in enumerate(LEVELS)}
#: residual RMS in the signal region / edge noise above which the residual is
#: flagged (the |z| > 3 / lag-1 > 0.4 thresholds stay in diagnostics.py, the
#: 0.95 identifiability threshold in identifiability.py)
NOISE_RATIO_LIMIT = 1.5
#: a population whose relative error reaches this fraction is not supported
#: by the data (help/glass-fitting.md section 5)
POPULATION_REL_ERR_LIMIT = 1.0
#: chip wording is bounded by construction; the strip elides at 48 characters
TEXT_LIMIT = 60

NO_FIT_TEXT = "no fit yet — F5 fits the lines"
STALE_SUFFIX = " · edited since fit"
STALE_PREFIX = "from the last fit — F5 to refresh: "
GLYPH = {"bad": "✗", "check": "⚠", "info": "·"}
#: the pill wording per level (``Fit`` becomes ``Model`` for a live model)
PILL = {"ok": "✓ Fit OK", "check": "⚠ Fit: check", "bad": "✗ Fit: not physical"}

#: display order of the flags: what makes the numbers unreadable first, then
#: the things to check, then information
KIND_ORDER = ("physical", "degenerate", "nocov", "at_bounds", "outside",
              "residual", "population", "tail", "recovery", "excitation", "frozen")
#: the flags derived from the quantify rows and the acquisition facts, the
#: ones ``with_quantification`` rebuilds without a new fit
QUANT_KINDS = ("population", "tail", "recovery", "excitation")
_WINDOW_MSG = "is outside the fit window"
_BOUND_RE = re.compile(r"^s(\d+)\.(.+)$")
#: the residual chip's tooltip: what usually leaves structure behind
RESIDUAL_CAUSES = ("usual causes: spinning sidebands not modelled (Decomposition "
                   "▸ Add spinning sidebands…), a phasing or baseline error, one "
                   "component too few — the fitted values may still be sound")
#: short parameter names for the chips, the table's own headers where the
#: registry key differs from them ('line' → lb); the model's registry key for
#: everything else (fwhm, dCS, cq, eta, …)
_SHORT = {"isotropic_chemical_shift_ppm": "pos", "sigma_Cq_MHz": "σCq",
          "amplitude": "amp", "line_fwhm_ppm": "lb"}


@dataclass
class Flag:
    """One thing the verdict says about the fit.

    ``kind`` names the check (see ``KIND_ORDER``), ``level`` is ``bad`` /
    ``check`` / ``info``, ``text`` is the chip wording, ``detail`` the one-line
    explanation (tooltip / details menu), ``target`` what a click opens
    (``residual`` / ``param`` / ``correlations`` / ``errors`` / ``report`` /
    ``widen`` (widen the integration window) / ``relaxation`` (Tools ▸
    Relaxation on the sibling T1 EXPNO) / ``flip`` (Experiment parameters,
    the 90° pulse) or empty) and ``params`` the ``(site index, parameter
    name)`` cells the flag concerns, the first being the focus target.
    """

    kind: str
    level: str
    text: str
    detail: str
    target: str
    params: list = field(default_factory=list)
    #: True on flags carried forward from the last fit by ``reassess_live``
    stale: bool = False
    #: degenerate / at_bounds / nocov / population / tail: derived from the
    #: FIT (its covariance or its quantify rows), so a live pass cannot
    #: recompute them
    covariance_based: bool = False


@dataclass
class Health:
    """The verdict on a model: the raw measurements plus the flags built from
    them, with every rendering the desktop app needs."""

    fitted: bool
    stale: bool
    rmsd: float | None
    redchi: float | None
    noise_ratio: float | None
    runs_z: float | None
    lag1: float | None
    warns: list                 # every sanity warning, suppressed ones included
    pairs: list                 # (name_i, name_j, r) with |r| >= threshold
    at_bounds: list
    frozen: list
    flags: list
    recipe_sig: tuple | None = None
    data_sig: tuple | None = None
    #: the quantitativity.Check behind the recovery / excitation flags (plain
    #: scalars; in memory only, like the rest of the snapshot); None when the
    #: source is not a Bruker dataset
    acquisition: object = None
    #: the quantify rows of the fit carried a tail measurement
    tail_checked: bool = False
    #: the residual failed the runs / lag-1 test (behind the residual chip)
    structured: bool = False
    #: acquisition facts that could NOT be judged (no sibling T1, no 90°
    #: pulse): tooltip lines, never chips
    unchecked: list = field(default_factory=list)

    # ------------------------------------------------------------ verdict
    @property
    def level(self) -> str:
        """``bad`` / ``check`` from the flags (info does not count); ``ok``
        for a fitted model with none; ``none`` for an unfitted one."""
        levels = {f.level for f in self.flags}
        if "bad" in levels:
            return "bad"
        if "check" in levels:
            return "check"
        return "ok" if self.fitted else "none"

    def kinds(self) -> set:
        return {f.kind for f in self.flags}

    def _word(self) -> str:
        return "Fit" if (self.fitted and not self.stale) else "Model"

    def pill_text(self, stale_hint: bool = False) -> str:
        """'✓ Fit OK' / '⚠ Fit: check' / '✗ Fit: not physical'; 'Model' in
        place of 'Fit' (and ' · edited since fit' after it) for a live model
        edited since the fit -- ``stale_hint`` forces that wording during a
        paddle drag, before the full pass has run."""
        lvl = self.level
        if lvl == "none":
            return NO_FIT_TEXT
        live = self.stale or stale_hint
        word = "Fit" if (self.fitted and not live) else "Model"
        suffix = STALE_SUFFIX if (self.fitted and live) else ""
        core = PILL[lvl]
        if word == "Model":
            core = core.replace("Fit", "Model", 1)
            if lvl == "ok":
                core = core[2:]                      # no tick on a live model
        return core + suffix

    def status_suffix(self) -> str:
        """The verdict appended to the status-bar 'fit done' message."""
        t = self.pill_text()
        for g in ("✓ ", "⚠ ", "✗ "):
            if t.startswith(g):
                t = t[len(g):]
                break
        return "  ·  " + t

    def chi_text(self) -> str:
        """'RMSD 0.0462 · χ²ᵣ 1.03' — the Fit-parameters footer text."""
        if self.rmsd is None:
            return ""
        t = f"RMSD {self.rmsd:.4f}"
        if self.redchi is not None and math.isfinite(self.redchi):
            t += f" · χ²ᵣ {self.redchi:.2f}"
        return t

    def passing(self) -> list:
        """The checks that ran and passed, for the pill tooltip."""
        kinds = self.kinds()
        out = []
        if "residual" not in kinds:
            if self.noise_ratio is not None:
                out.append(f"residual ≈ noise ({self.noise_ratio:.1f}×)")
            if self.runs_z is not None:
                out.append("residual unstructured")
        if "physical" not in kinds:
            out.append("values physical")
        if self.fitted and not self.stale:
            if "nocov" not in kinds:
                out.append("errors from covariance")
                if "degenerate" not in kinds:
                    out.append("parameters separable")
            if "at_bounds" not in kinds:
                out.append("none at a bound")
            if self.tail_checked and "tail" not in kinds:
                out.append(f"tails inside the window (≤ {TAIL_LIMIT_PCT:g} %)")
        if self.acquisition is not None:
            try:
                out.extend(self.acquisition.passing_lines(
                    exclude=kinds & {"recovery", "excitation"}))
            except Exception:
                pass
        return out

    def tooltip(self) -> str:
        if self.fitted and not self.stale:
            chi = self.chi_text()
            head = "Fit health — last fit" + (f" ({chi})" if chi else "")
        elif self.fitted:
            head = ("Model health — edited since the last fit; covariance-based "
                    "flags are from the last fit (F5 to refresh)")
        else:
            head = "Model health — not fitted yet (F5 to fit)"
        lines = [head]
        passing = self.passing()
        if passing:
            lines.append(" · ".join(passing))
        if self.fitted and self.acquisition is None:
            lines.append("acquisition not checked (not a Bruker dataset)")
        for f in self.flags:
            lines.append(f"{GLYPH[f.level]} {f.detail}"
                         + (" — from the last fit" if f.stale else ""))
        for u in self.unchecked:
            lines.append(f"· {u}")
        if self.flags:
            lines.append("click a chip for its detail · F7 lists everything")
        else:
            lines.append("F7 lists the checks and the analysis tools")
        return "\n".join(lines)

    # ---------------------------------------------- legacy Report header
    def summary(self) -> str:
        """The Report-dock header sentence, in the wording the app has always
        shown there (bits joined by '   ·   ')."""
        kinds = self.kinds()
        bits = []
        chi = self.chi_text()
        if chi:
            bits.append(chi)
        if self.noise_ratio is not None:
            bits.append("residual within noise"
                        if self.noise_ratio < NOISE_RATIO_LIMIT
                        else f"⚠ residual {self.noise_ratio:.1f}× noise "
                             "(structure left)")
        if self.structured:
            bits.append("⚠ structured residual")
        if self.frozen:
            bits.append("frozen: " + ", ".join(self.frozen))
        if self.at_bounds:
            bits.append("⚠ at bounds: " + ", ".join(self.at_bounds))
        if self.warns:
            n = len(self.warns)
            bits.append(f"⚠ {n} physical warning{'s' if n != 1 else ''}")
        if self.pairs:
            n = len(self.pairs)
            bits.append(f"⚠ {n} unidentifiable pair{'s' if n != 1 else ''} "
                        "(see Correlations)")
        if "nocov" in kinds:
            bits.append("⚠ no covariance (no error bars)")
        for f in self.flags:
            if f.kind in QUANT_KINDS and f.level == "check":
                bits.append("⚠ " + f.text)
        return "   ·   ".join(bits)

    def summary_tooltip(self) -> str:
        """The Report-header tooltip: structure message, sanity summary and the
        first eight unidentifiable pairs (lmfit names)."""
        lines = []
        struct = next((f for f in self.flags if f.kind == "residual"), None)
        if struct is not None:
            lines.append(struct.detail)
        if sanity.summarize(self.warns):
            lines.append(sanity.summarize(self.warns))
        if self.pairs:
            lines.append("unidentifiable: " + ", ".join(
                f"{a}↔{b} ({r:+.2f})" for a, b, r in self.pairs[:8]))
        return "\n".join(lines)


# ---------------------------------------------------------------- helpers
def residual_noise_ratio(y_exp, y_fit) -> float | None:
    """RMS of the residual in the signal region ÷ the baseline noise (RMS of
    the quiet edges). ≈1 means the model captures the data down to the noise;
    ≫1 means there is unmodelled structure left in the residual. None when
    the arrays are too short (< 40 points), mismatched, or the edges carry
    no noise at all."""
    try:
        y = np.asarray(y_exp, float)
        f = np.asarray(y_fit, float)
        r = y - f
        n = r.size
        if n < 40:
            return None
        edge = max(5, n // 12)
        noise = float(np.std(np.concatenate([r[:edge], r[-edge:]])))
        if noise <= 0:
            return None
        sig = np.abs(f) > 0.05 * (np.abs(f).max() or 1.0)
        r_sig = r[sig] if sig.any() else r
        return float(np.sqrt(np.mean(r_sig ** 2)) / noise)
    except Exception:
        return None


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return repr(v)
    return "nan" if math.isnan(x) else x


def _site_signature(site) -> tuple:
    if isinstance(site, dict):
        items = []
        for pn, p in (site.get("params") or {}).items():
            if isinstance(p, dict):
                v, vary, expr = p.get("value"), p.get("vary", True), p.get("expr")
            else:
                v = getattr(p, "value", None)
                vary, expr = getattr(p, "vary", True), getattr(p, "expr", None)
            items.append((str(pn), _num(v), bool(vary), expr or ""))
        return (str(site.get("model", "")), tuple(sorted(items)))
    items = [(str(pn), _num(p.value), bool(p.vary), p.expr or "")
             for pn, p in site.params.items()]
    return (str(site.model), tuple(sorted(items)))


def recipe_signature(recipe) -> tuple:
    """What 'edited since the fit' compares: per site the model and every
    parameter's (name, value, vary, expr), plus the fit zones. Labels, bounds,
    stderr, reference traces and hidden flags do not count. O(n params), no
    serialisation; accepts a recipe dict or a Recipe."""
    if recipe is None:
        return ()
    if isinstance(recipe, Recipe):
        sites, zones = recipe.sites, recipe.fit_zones or []
    else:
        sites, zones = recipe.get("sites") or [], recipe.get("fit_zones") or []
    zsig = []
    for z in zones:
        try:
            zsig.append(tuple(_num(v) for v in z))
        except TypeError:
            zsig.append((repr(z),))
    return tuple(_site_signature(s) for s in sites) + tuple(zsig)


def data_signature(ppm, amp) -> tuple:
    """An O(1) value-based stamp of the experimental arrays (size, both ends,
    the middle sample) — survives the copies workspace snapshots make, where
    an identity check would not."""
    try:
        n = len(amp)
    except TypeError:
        return ()
    if not n:
        return ()
    a = amp
    sig = [n, float(a[0]), float(a[-1]), float(a[n // 2])]
    if ppm is not None:
        try:
            if len(ppm) == n:
                sig += [float(ppm[0]), float(ppm[-1])]
        except TypeError:
            pass
    return tuple(sig)


def short_name(lmfit_name: str) -> str:
    """The short form the correlation dialog shows: 's0_amplitude' → 's0.amp'."""
    return (str(lmfit_name).replace("isotropic_chemical_shift_ppm", "pos")
            .replace("sigma_Cq_MHz", "σCq").replace("shift_fwhm_ppm", "dCS")
            .replace("amplitude", "amp").replace("_", "."))


def parse_bound_name(name: str):
    """'s0.eta' (fit.py's at-bounds naming) → (0, 'eta'); None if it does not
    parse."""
    m = _BOUND_RE.match(str(name))
    if not m:
        return None
    return int(m.group(1)), m.group(2)


def short_param(model: str | None, pname: str) -> str:
    """'lb' for line_fwhm_ppm, 'σCq' for sigma_Cq_MHz, …; the registry key
    of the model's ParamDef for anything the short table does not know."""
    if pname in _SHORT:
        return _SHORT[pname]
    try:
        from larmor import models as model_registry

        return model_registry.get(model).key_of(pname)
    except Exception:
        return str(pname)


def param_label(rec: Recipe, i: int, pname: str) -> str:
    """'A lb' — the table letter of site ``i`` and the short parameter name,
    the way a chip names a cell."""
    from larmor.cellparse import index_to_letter

    model = rec.sites[i].model if 0 <= i < len(rec.sites) else None
    return f"{index_to_letter(i)} {short_param(model, pname)}"


def _bound_words(rec: Recipe, at_bounds: list) -> tuple:
    """The at-bounds names of a fit ('s0.line_fwhm_ppm') as cells and words:
    ``(cells, short, long)`` — ``short`` 'A lb at its lower bound 0' per
    parameter (for the chip), ``long`` the same with the value (for the
    tooltip). A name whose parameter is no longer free (pinned or linked since
    the fit, or held at its model default) is dropped: a held value is not
    'at a bound', whatever the fit said."""
    from larmor.paramstatus import effective_bounds, param_status

    cells, short, long = [], [], []
    for name in at_bounds:
        parsed = parse_bound_name(name)
        if parsed is None:
            short.append(str(name))
            long.append(str(name))
            continue
        i, pname = parsed
        p = rec.sites[i].params.get(pname) if 0 <= i < len(rec.sites) else None
        model = rec.sites[i].model if p is not None else None
        if p is not None and (p.expr or not p.vary):
            continue
        label = param_label(rec, i, pname) if p is not None else name
        st = param_status(model, pname, p) if p is not None else None
        if st is not None and st.kind == "at_bound":
            side, bound = st.side, st.bound
        elif p is not None:
            lo, hi = effective_bounds(model, pname, p)
            side, bound = ("min", lo) if lo is not None else ("max", hi)
        else:
            side, bound = "", None
        word = {"min": "lower", "max": "upper"}.get(side, "")
        b = f" {bound:.4g}" if bound is not None else ""
        cells.append((i, pname))
        short.append(f"{label} at its {word} bound{b}".replace("  ", " "))
        value = f" (value {float(p.value):.4g})" if p is not None else ""
        long.append(f"{label} at its {word} bound{b}{value}")
    return cells, short, long


def _bounded(prefix: str, names: list) -> str:
    text = f"{prefix}: " + ", ".join(names)
    if len(text) > TEXT_LIMIT:
        text = f"{prefix} ×{len(names)}: {names[0]}, …"
        if len(text) > TEXT_LIMIT:
            text = f"{prefix} ×{len(names)}"
    return text


def _frozen_indices(rec: Recipe, frozen) -> dict:
    """Map each frozen name (fit.py: label or 'site-<i>') to its site index."""
    out = {}
    for name in frozen:
        for i, s in enumerate(rec.sites):
            if (s.label and s.label == name) or name == f"site-{i}":
                out[name] = i
                break
    return out


def _aligned(y_exp, y_fit, ppm, x_fit):
    """The experimental and model arrays on one axis (model interpolated from
    its own axis when that differs — kernel models simulate on a wider grid),
    or (None, None) when no residual can be formed."""
    if y_exp is None or y_fit is None:
        return None, None
    try:
        y = np.asarray(y_exp, float).ravel()
        f = np.asarray(y_fit, float).ravel()
        if x_fit is not None and ppm is not None:
            xf = np.asarray(x_fit, float).ravel()
            xp = np.asarray(ppm, float).ravel()
            if xf.size == f.size and xp.size == y.size and (
                    xf.size != xp.size or not np.array_equal(xf, xp)):
                order = np.argsort(xf)
                f = np.interp(xp, xf[order], f[order])
    except Exception:
        return None, None
    if y.size != f.size or y.size < 40:
        return None, None
    return y, f


def _unsupported_populations(rows) -> list:
    out = []
    for r in rows or []:
        err = r.get("fraction_err_pct")
        if err is None:
            continue
        frac = r.get("fraction_pct") or 0.0
        try:
            if not math.isfinite(float(err)):
                continue
            if frac <= 0 or float(err) >= POPULATION_REL_ERR_LIMIT * float(frac):
                out.append(str(r.get("label") or r.get("site") or "?"))
        except (TypeError, ValueError):
            continue
    return out


def _tail_rows(rows) -> list:
    """[(label, tail_outside_pct)] of the quantify rows whose line lies more
    than TAIL_LIMIT_PCT outside the integration window; rows without the
    key → []."""
    out = []
    for r in rows or []:
        pct = r.get("tail_outside_pct")
        if pct is None:
            continue
        try:
            pct = float(pct)
        except (TypeError, ValueError):
            continue
        if math.isfinite(pct) and pct > TAIL_LIMIT_PCT:
            idx = None
            m = re.match(r"^s(\d+)$", str(r.get("site", "")))
            if m:
                idx = int(m.group(1))
            out.append((str(r.get("label") or r.get("site") or "?"), pct, idx))
    return out


def _tails_checked(rows) -> bool:
    return any(r.get("tail_outside_pct") is not None for r in rows or [])


def _quant_flags(rec: Recipe, rows, facts, fitted: bool) -> tuple:
    """The flags derived from the quantify rows (population, tail — a FIT's
    numbers, so fitted only) and from the acquisition facts (recovery,
    excitation — live, and only when the fact IS known: a missing T1 or an
    unknown flip angle is a tooltip line, not a flag). Returns (flags,
    quantitativity.Check or None)."""
    flags = []
    if fitted:
        unsupported = _unsupported_populations(rows)
        if unsupported:
            flags.append(Flag(
                kind="population", level="check",
                text=_bounded("population ±≥100 %", unsupported),
                detail=("relative uncertainty ≥ 100 % on the population of "
                        + ", ".join(unsupported)
                        + " — the data do not support that component (report "
                        "it, or remove it)"),
                target="report", covariance_based=True))
        tails = _tail_rows(rows)
        if tails:
            flags.append(Flag(
                kind="tail", level="check",
                text=_bounded("tail outside window",
                              [f"{label} {pct:.0f} %" for label, pct, _ in tails]),
                detail=("; ".join(f"{pct:.0f} % of {label}'s simulated area lies "
                                  "outside the integration window"
                                  for label, pct, _ in tails)
                        + f" (limit {TAIL_LIMIT_PCT:g} %) — its population is "
                        "biased low; click to widen the window to 99.5 % of "
                        "every line, then F6 re-integrates and F5 refits over it"),
                target="widen",
                params=[(i, "amplitude") for _, _, i in tails if i is not None],
                covariance_based=True))
    chk = None
    if facts is not None:
        prov = getattr(rec, "provenance", None)
        override = prov.get("quantitativity") if isinstance(prov, dict) else None
        try:
            chk = quantitativity.check(facts, rec.sites, override)
        except Exception:
            chk = None
    if chk is not None and chk.acquisition is not None:
        if chk.recovery_min is not None and chk.recovery_min < RECOVERY_MIN:
            flags.append(Flag(
                kind="recovery", level="check", text=chk.recovery_text(),
                detail=chk.recovery_detail(), target="relaxation"))
        if chk.excitation_judged and chk.excitation_over():
            flags.append(Flag(
                kind="excitation", level="check", text=chk.excitation_text(),
                detail=chk.excitation_detail(), target="flip"))
    return flags, chk


def _unchecked(chk) -> list:
    """The acquisition facts a Check could not judge, as tooltip lines."""
    if chk is None:
        return []
    try:
        return list(chk.unchecked_lines())
    except Exception:
        return []


def residual_flag(noise_ratio, struct) -> Flag | None:
    """The one amber residual chip: a residual at or above NOISE_RATIO_LIMIT
    times the edge noise and/or a structured residual (``struct`` is
    ``diagnostics.residual_structure``'s dict, or None). Its tooltip gives
    the measurement and the usual causes; None when neither applies."""
    structured = bool(struct is not None and struct.get("structured"))
    noisy = noise_ratio is not None and noise_ratio >= NOISE_RATIO_LIMIT
    if not (structured or noisy):
        return None
    if noisy:
        text = f"residual {noise_ratio:.1f}× noise"
        measured = (f"residual RMS in the signal region is {noise_ratio:.1f}× "
                    f"the edge noise (limit {NOISE_RATIO_LIMIT:g}×)")
        if structured:
            text += ", structured"
            measured += "; " + str(struct.get("message", ""))
    else:
        text = "residual: structure left"
        measured = str(struct.get("message", ""))
    return Flag(kind="residual", level="check", text=text,
                detail=(f"the residual carries structure the model does not "
                        f"describe — {measured}; {RESIDUAL_CAUSES}"),
                target="residual")


def _sort_flags(flags: list) -> list:
    return sorted(flags, key=lambda f: KIND_ORDER.index(f.kind)
                  if f.kind in KIND_ORDER else len(KIND_ORDER))


# ------------------------------------------------------------------ assess
def assess(recipe, y_exp, y_fit, *, lmfit_result=None, at_bounds=(), frozen=(),
           window=None, rmsd=None, quant_rows=None, fitted=True, ppm=None,
           x_fit=None, acquisition=None) -> Health:
    """Judge a model against the data.

    ``recipe`` is a Recipe or a recipe dict (converted with Recipe.from_dict,
    never mutated); ``window`` defaults to the recipe's fit window. Residual
    flags need ``y_exp`` / ``y_fit`` on the same axis with at least 40 points
    (``x_fit`` + ``ppm`` interpolate a model simulated on its own grid);
    covariance flags (degenerate, no covariance, at bounds, population, tail)
    are only emitted for a fitted model. ``acquisition`` is the
    ``quantitativity.AcqFacts`` of the source EXPNO (None for a non-Bruker
    source): the recycle-delay and flip-angle flags are judged from it and
    the current sites, fitted or not, and only when T1 / the flip angle are
    known. Sanity's δiso-outside-window warning is a ``check`` flag of its
    own (``outside``), not repeated for a site the fit froze for exactly that
    reason (it stays in ``warns``; the frozen chip reports it once); the other
    sanity warnings are the one ``bad`` flag.
    """
    rec = recipe if isinstance(recipe, Recipe) else Recipe.from_dict(recipe)
    if window is None:
        window = rec.fit_window_ppm
    fitted = bool(fitted)
    # at-bounds, like every covariance-derived flag, describes a FIT
    at_bounds = [str(n) for n in (at_bounds or [])] if fitted else []
    frozen = [str(n) for n in (frozen or [])]

    redchi = getattr(lmfit_result, "redchi", None) if lmfit_result is not None else None
    try:
        redchi = float(redchi) if redchi is not None and math.isfinite(float(redchi)) else None
    except (TypeError, ValueError):
        redchi = None
    try:
        rmsd = float(rmsd) if rmsd is not None and math.isfinite(float(rmsd)) else None
    except (TypeError, ValueError):
        rmsd = None

    # residual: within the noise? white?
    noise_ratio = runs_z = lag1 = None
    struct = None
    y, f = _aligned(y_exp, y_fit, ppm, x_fit)
    if y is not None:
        noise_ratio = residual_noise_ratio(y, f)
        struct = diagnostics.residual_structure(y - f)
        runs_z, lag1 = float(struct["runs_z"]), float(struct["lag1"])

    # physical values (sanity): the impossible ones are the red flag, a δiso
    # outside the fit window is a thing to check (minus the frozen sites,
    # whose chip already says so)
    warns = sanity.check_recipe(rec, window)
    frozen_idx = _frozen_indices(rec, frozen)
    frozen_sites = set(frozen_idx.values())
    outside = [w for w in warns
               if w["param"] == "isotropic_chemical_shift_ppm"
               and _WINDOW_MSG in w["message"]]
    impossible = [w for w in warns if w not in outside]
    outside = [w for w in outside if w["site"] not in frozen_sites]

    pairs = unidentifiable_pairs(lmfit_result) if fitted else []
    names = list(getattr(lmfit_result, "var_names", []) or []) if lmfit_result is not None else []
    nocov = fitted and (lmfit_result is None or not names
                        or getattr(lmfit_result, "covar", None) is None)
    quant, chk = _quant_flags(rec, quant_rows, acquisition, fitted)

    flags = []
    if impossible:
        flags.append(Flag(
            kind="physical", level="bad", text=f"unphysical ×{len(impossible)}",
            detail="; ".join(f"{w['label']}: {w['message']}" for w in impossible),
            target="param",
            params=[(int(w["site"]), str(w["param"])) for w in impossible]))
    if pairs:
        n = len(pairs)
        a, b, r = pairs[0]
        text = f"degenerate pair: {short_name(a)}↔{short_name(b)} ({r:+.2f})"
        if n > 1:
            text = f"degenerate ×{n}: {short_name(a)}↔{short_name(b)} ({r:+.2f})"
        if len(text) > TEXT_LIMIT:
            text = f"degenerate ×{n}"
        listed = ", ".join(f"{short_name(a)}↔{short_name(b)} ({r:+.2f})"
                           for a, b, r in pairs[:8])
        if n > 8:
            listed += f" … +{n - 8} more"
        flags.append(Flag(
            kind="degenerate", level="check", text=text,
            detail=(f"the data cannot separate these parameters (|r| ≥ "
                    f"{IDENTIFIABILITY_THRESHOLD:.2f}): {listed} — their "
                    "values are fine to read, their error bars are not; fix "
                    "or link one, or add a constraint from chemistry"),
            target="correlations", covariance_based=True))
    if nocov:
        flags.append(Flag(
            kind="nocov", level="check", text="no error bars (no covariance)",
            detail=("the fit returned no covariance matrix, so no standard "
                    "error is available — a parameter at a bound or a "
                    "degenerate pair usually causes this; Errors Analysis "
                    "(χ² profile) still gives confidence intervals"),
            target="errors", covariance_based=True))
    if at_bounds:
        cells, short, long = _bound_words(rec, at_bounds)
        if short:
            flags.append(Flag(
                kind="at_bounds", level="check",
                text=_bounded("at a bound", short),
                detail=("finished at a bound — its uncertainty is conditional "
                        "on the bound: loosen it if a larger range is "
                        "physical, or hold the parameter on purpose: "
                        + "; ".join(long)),
                target="param", params=cells, covariance_based=True))
    if outside:
        flags.append(Flag(
            kind="outside", level="check",
            text=_bounded("δiso outside the fit window",
                          [str(w["label"]) for w in outside]),
            detail=("; ".join(f"{w['label']}: {w['message']}" for w in outside)
                    + " — widen the window or move the line; the fit cannot "
                    "constrain a centre it does not see"),
            target="param",
            params=[(int(w["site"]), str(w["param"])) for w in outside]))
    structured = bool(struct is not None and struct["structured"])
    res = residual_flag(noise_ratio, struct)
    if res is not None:
        flags.append(res)
    flags.extend(quant)
    if frozen:
        flags.append(Flag(
            kind="frozen", level="info", text=_bounded("frozen", frozen),
            detail=("held by the fit — centre outside the fit window, "
                    "amplitude set to 0 and parameters fixed: "
                    + ", ".join(frozen)),
            target="param",
            params=[(frozen_idx[n], "isotropic_chemical_shift_ppm")
                    for n in frozen if n in frozen_idx]))

    return Health(fitted=fitted, stale=False, rmsd=rmsd, redchi=redchi,
                  noise_ratio=noise_ratio, runs_z=runs_z, lag1=lag1,
                  warns=warns, pairs=list(pairs), at_bounds=at_bounds,
                  frozen=frozen, flags=_sort_flags(flags), acquisition=chk,
                  tail_checked=bool(fitted and _tails_checked(quant_rows)),
                  structured=structured, unchecked=_unchecked(chk))


def reassess_live(prev, recipe_dict, y_exp, y_model_on_exp, *, ppm=None,
                  window=None, acquisition=None) -> Health:
    """Re-judge the live model after an edit, cheaply.

    With a previous fit (``prev`` fitted): when the recipe and data signatures
    still match the fit, ``prev`` is returned unchanged (a Compute / theme
    re-simulation keeps the fit verdict). Otherwise the residual and physical
    checks are recomputed from ``y_model_on_exp`` (the live model on the
    experimental axis) and the covariance-based flags of the fit are carried
    forward marked stale. Before any fit, only the physical checks run — a
    hand-placed model is always far from the data, so its residual says
    nothing yet. The acquisition facts (``acquisition``, else the ones behind
    ``prev``) are re-judged against the current sites on every pass: the
    recycle-delay and flip-angle chips are live, the tail chip is carried.
    """
    facts = acquisition if acquisition is not None else getattr(
        getattr(prev, "acquisition", None), "facts", None)
    fitted = prev is not None and prev.fitted
    if fitted:
        if (recipe_signature(recipe_dict) == prev.recipe_sig
                and data_signature(ppm, y_exp) == prev.data_sig):
            return prev
        h = assess(recipe_dict, y_exp, y_model_on_exp, frozen=prev.frozen,
                   window=window, fitted=False, ppm=ppm, acquisition=facts)
        live = [replace(f, stale=True) if f.kind == "frozen" else f
                for f in h.flags if not f.covariance_based]
        carried = [replace(f, stale=True) for f in prev.flags if f.covariance_based]
        h.flags = _sort_flags(live + carried)
        h.fitted, h.stale = True, True
        h.rmsd, h.redchi = prev.rmsd, prev.redchi
        h.pairs, h.at_bounds = list(prev.pairs), list(prev.at_bounds)
        h.recipe_sig, h.data_sig = prev.recipe_sig, prev.data_sig
        h.tail_checked = prev.tail_checked
        return h
    return assess(recipe_dict, None, None, window=window, fitted=False, ppm=ppm,
                  acquisition=facts)


def with_quantification(prev: Health, recipe, quant_rows, acquisition=None) -> Health:
    """``prev`` with its quantification-derived flags (population, tail,
    recovery, excitation) rebuilt from fresh quantify rows and acquisition
    facts — a re-integration over a wider window or a typed 90° pulse / T1
    re-judged without a new fit. Everything else (residual flags, signatures,
    stale, rmsd) is kept; ``acquisition`` defaults to the facts behind
    ``prev``."""
    rec = recipe if isinstance(recipe, Recipe) else Recipe.from_dict(recipe)
    facts = acquisition if acquisition is not None else getattr(
        getattr(prev, "acquisition", None), "facts", None)
    keep = [f for f in prev.flags if f.kind not in QUANT_KINDS]
    quant, chk = _quant_flags(rec, quant_rows, facts, prev.fitted)
    return replace(prev, flags=_sort_flags(keep + quant), acquisition=chk,
                   tail_checked=bool(prev.fitted and _tails_checked(quant_rows)),
                   unchecked=_unchecked(chk))
