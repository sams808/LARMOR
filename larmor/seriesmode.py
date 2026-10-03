"""Series mode of the workbench (Qt-free).

A spectral series -- a composition or temperature run -- is fitted one
spectrum at a time in the MAIN window: every member is an ordinary workspace
of the Workspaces dock, so every editing tool (add / remove lines, bounds and
links, paddles, the processing panel, Auto fit, undo) applies per member, and
a thin series bar above the plot walks the series, carries a model from one
neighbour to the next and runs the forward / backward sweeps. This module
holds what that mode needs without Qt:

* :class:`SeriesSpec` / :class:`SeriesMember` -- the series record: an id,
  the ordered members (name, replicate group, sample folder, title line,
  proc number, source path, the acqus / procs parameters read for the
  comparability check) and the options (seed on move, the carried
  parameters, passes / start / smooth of the auto sweep). ``tag_for`` is the
  small dict each member workspace carries (``ws["series"]``), persisted in
  the project bundle and rebuilt by :meth:`SeriesSpec.from_tags`.
* :func:`carry_into` -- the carry rules between two members: an EMPTY target
  takes a copy of the source's lines (structure, labels, families, links --
  valid because the structure is identical) with every amplitude scaled by
  the ratio of the two spectra's maxima (a display seed; the fit's own
  analytic pre-scale does the rest); a target WITH lines keeps its own
  structure (the user's add / remove per spectrum is respected) and receives
  :func:`larmor.seqfit.seed_from` values for the carried parameters on the
  lines that pair by LABEL (:func:`larmor.components.pair_sites`; by index
  only where a label cannot identify its line), clipped to its own bounds,
  links untouched. When the two carry happens is the window's business
  (``mw_series``): on a move only an empty member takes anything, so a fit
  already made is never disturbed by walking the series; the value carry is
  an explicit action. A member marked ``locked`` ("keep this fit") is left
  alone by every carry and by the sweep, which still starts its neighbours
  from it.
* :func:`member_status` / :func:`rmsd_of` -- what the bar's status dots say.
* :func:`entries_for_sweep` / :func:`apply_sweep_result` -- the bridge to
  :func:`larmor.seqfit.run_sequential` (the auto sweep; ``larmor seqfit``
  runs the same function from the command line) and back, by series order.

Nothing here imports the desktop; ``tests/test_seriesmode.py`` covers it.
"""
from __future__ import annotations

import copy
import datetime as _dt
import hashlib
import secrets
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from larmor import fithealth
from larmor.components import label_of, pair_sites, usable_labels
from larmor.recipe import Recipe

__all__ = [
    "SeriesMember", "SeriesOptions", "SeriesSpec", "PASS_CHOICES",
    "SMOOTH_CHOICES", "START_CHOICES", "STATUSES", "proc_number",
    "new_series_id", "signature_key", "amp_max", "carry_candidates",
    "carry_into", "member_status", "rmsd_of", "entries_for_sweep",
    "apply_sweep_result", "slug", "auto_name", "member_records",
]

#: the auto-sweep choices the bar offers
PASS_CHOICES = (1, 2, 3, 4, 5, 6)
START_CHOICES = ("first", "last")
SMOOTH_CHOICES = (0, 3, 5)
#: member statuses, in the order the legend reads them
STATUSES = ("unfitted", "fitted", "edited", "failed")

#: parameter names that are a line's position / width, for the carry note
_POSITION_NAMES = frozenset({"isotropic_chemical_shift_ppm", "shift_ppm"})


# ------------------------------------------------------------------ helpers
def proc_number(path) -> str:
    """The ``pdata/<N>`` number of a Bruker path ('' for anything else) --
    the Qt-free twin of the batch dialog's helper."""
    parts = Path(str(path)).parts
    if "pdata" in parts:
        i = parts.index("pdata")
        if i + 1 < len(parts):
            return parts[i + 1]
    return ""


def new_series_id() -> str:
    """A fresh series id: the start time plus two random bytes, so two series
    started within a second (or saved in two projects) never collide."""
    return ("series-" + _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
            + "-" + secrets.token_hex(2))


def signature_key(recipe) -> str:
    """A short, JSON-friendly stamp of ``fithealth.recipe_signature`` (per
    site the model and every parameter's value / vary / expr, plus the fit
    zones): what 'edited since the fit' compares across workspaces and
    project reloads. '' for no recipe / no lines."""
    if not recipe:
        return ""
    sites = recipe.sites if isinstance(recipe, Recipe) else recipe.get("sites")
    if not sites:
        return ""
    sig = fithealth.recipe_signature(recipe)
    return hashlib.sha1(repr(sig).encode("utf-8")).hexdigest()[:16]


def amp_max(amp) -> float:
    """The largest |value| of a spectrum (0.0 for an empty / unusable array)."""
    try:
        a = np.asarray(amp, float)
    except (TypeError, ValueError):
        return 0.0
    if a.size == 0:
        return 0.0
    m = float(np.nanmax(np.abs(a)))
    return m if np.isfinite(m) else 0.0


def slug(s: str) -> str:
    """A file-name-safe stem (the batch dialog's rule)."""
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in (s or ""))[:80]


def auto_name(recipe, name: str = "", when=None) -> str:
    """'sample_nucleus_seq_YYYYMMDD_HHMM' -- the automatic file stem of a
    member's saved fit; ``name`` (the member's series name) comes before
    the recipe's sample."""
    rec = recipe.to_dict() if isinstance(recipe, Recipe) else (recipe or {})
    when = when or _dt.datetime.now()
    parts = [name or rec.get("sample") or "fit", rec.get("nucleus") or "", "seq",
             when.strftime("%Y%m%d_%H%M")]
    return slug("_".join(p for p in parts if p))


def _sites(recipe) -> list:
    if recipe is None:
        return []
    if isinstance(recipe, Recipe):
        return [s.__dict__ for s in recipe.sites]
    return list(recipe.get("sites") or [])


def carry_candidates(recipes) -> tuple:
    """Every parameter name present in the members' lines except amplitude,
    sorted -- the checklist of the Carry menu and the default carried set
    (``batchfit.all_but_amplitude`` on dicts, no Recipe conversion)."""
    names: set = set()
    for rec in recipes or []:
        for s in _sites(rec):
            params = s.get("params") if isinstance(s, dict) else s.params
            names.update((params or {}).keys())
    names.discard("amplitude")
    return tuple(sorted(names))


def _describe_carried(names) -> str:
    cats = []
    names = list(names)
    if any(n in _POSITION_NAMES for n in names):
        cats.append("positions")
    if any("fwhm" in n for n in names):
        cats.append("widths")
    if any(n not in _POSITION_NAMES and "fwhm" not in n for n in names):
        cats.append("shape parameters" if cats else "parameters")
    if not cats:
        return "nothing"
    if len(cats) == 1:
        return cats[0]
    return ", ".join(cats[:-1]) + " and " + cats[-1]


# ------------------------------------------------------------- the record
@dataclass
class SeriesMember:
    """One spectrum of the series, by workspace ``key``."""
    key: str
    name: str
    group: str = ""
    folder: str = ""
    title: str = ""
    proc: str = ""
    source_path: str = ""
    #: ``comparability.SpectrumParams`` of the source (None for CSV / dmfit);
    #: in memory only, re-read when a project is reopened
    params: object = None
    #: the fit the status dot refers to: ``signature_key`` of the recipe the
    #: last successful fit left, '' when never fitted; ``failed`` after a fit
    #: that raised
    fit_sig: str = ""
    failed: bool = False
    #: "keep this fit": no carry reaches it and a sweep never refits it (it
    #: still seeds its neighbours)
    locked: bool = False


@dataclass
class SeriesOptions:
    seed_on_move: bool = True
    #: parameter names unticked in the Carry menu (the carried set is every
    #: name present minus these, so a model added later carries by default)
    carry_off: tuple = ()
    passes: int = 2
    start: str = "first"
    smooth: int = 0

    def to_dict(self) -> dict:
        return {"seed_on_move": bool(self.seed_on_move),
                "carry_off": list(self.carry_off),
                "passes": int(self.passes), "start": str(self.start),
                "smooth": int(self.smooth)}

    @classmethod
    def from_dict(cls, d) -> "SeriesOptions":
        d = d or {}
        passes = int(d.get("passes", 2) or 2)
        smooth = int(d.get("smooth", 0) or 0)
        start = str(d.get("start", "first") or "first")
        return cls(seed_on_move=bool(d.get("seed_on_move", True)),
                   carry_off=tuple(str(n) for n in (d.get("carry_off") or [])),
                   passes=passes if passes in PASS_CHOICES else 2,
                   start=start if start in START_CHOICES else "first",
                   smooth=smooth if smooth in SMOOTH_CHOICES else 0)


def member_records(paths, recipes, names=None) -> list:
    """The members of a series from its source paths and loaded recipe dicts
    (what the old dialog's ``_load`` collected): the name is
    ``io.scan.sample_label`` (or the given name) made distinct by
    ``io.scan.disambiguate``, the group the label before disambiguation
    (replicates share it), folder and title from the loader's provenance,
    the proc number from the path and the acqus / procs parameters from
    ``comparability.read_params`` (None where there are none)."""
    from larmor import comparability
    from larmor.io.scan import disambiguate, sample_label

    paths = [str(p) for p in paths]
    recipes = list(recipes) + [None] * (len(paths) - len(recipes))
    given = list(names or []) + [None] * (len(paths) - len(names or []))
    base = []
    for p, rec, nm in zip(paths, recipes, given):
        nm = (nm or "").strip() if isinstance(nm, str) else ""
        base.append(nm or sample_label(p, rec or {}) or Path(p).stem)
    labels = disambiguate(list(base), paths)
    out = []
    for k, (p, rec) in enumerate(zip(paths, recipes)):
        prov = (rec or {}).get("provenance") or {}
        try:
            params = comparability.read_params(p)
        except Exception:                                 # noqa: BLE001
            params = None
        out.append(SeriesMember(
            key=f"m{k}", name=labels[k], group=base[k],
            folder=str(prov.get("sample_folder", "") or ""),
            title=str(prov.get("title", "") or ""), proc=proc_number(p),
            source_path=p, params=params))
    return out


@dataclass
class SeriesSpec:
    """The active series: id, ordered members, options, the comparability
    verdict and the Series table (names / groups / order / composition
    columns) the Series table dialog edits."""
    id: str
    members: list = field(default_factory=list)
    options: SeriesOptions = field(default_factory=SeriesOptions)
    comparison: object = None
    table: object = None

    # ---- construction
    @classmethod
    def build(cls, paths, recipes, names=None, series_id: str | None = None,
              options: SeriesOptions | None = None) -> "SeriesSpec":
        sid = series_id or new_series_id()
        members = member_records(paths, recipes, names)
        for k, m in enumerate(members):
            m.key = f"{sid}:{k}"
        spec = cls(id=sid, members=members, options=options or SeriesOptions())
        spec.refresh_comparison()
        return spec

    @classmethod
    def from_tags(cls, tags: list, params: list | None = None) -> "SeriesSpec":
        """Rebuild a series from member tags (``ws["series"]`` dicts, any
        order; all of one id). ``params`` are the re-read comparability
        parameters aligned with ``tags``."""
        tags = [t for t in tags if isinstance(t, dict)]
        if not tags:
            raise ValueError("no series tags")
        order = sorted(range(len(tags)), key=lambda i: int(tags[i].get("index", i)))
        sid = str(tags[order[0]].get("id") or new_series_id())
        members = []
        table = None
        for i in order:
            t = tags[i]
            members.append(SeriesMember(
                key=str(t.get("key") or f"{sid}:{len(members)}"),
                name=str(t.get("name") or f"spectrum {len(members) + 1}"),
                group=str(t.get("group") or t.get("name") or ""),
                folder=str(t.get("folder") or ""), title=str(t.get("title") or ""),
                proc=str(t.get("proc") or ""),
                source_path=str(t.get("source_path") or ""),
                params=(params[i] if params is not None and i < len(params) else None),
                fit_sig=str(t.get("fit_sig") or ""), failed=bool(t.get("failed")),
                locked=bool(t.get("locked"))))
            if table is None and isinstance(t.get("table"), dict):
                table = t["table"]
        opts = SeriesOptions.from_dict(tags[order[0]].get("options"))
        spec = cls(id=sid, members=members, options=opts)
        if table is not None:
            try:
                from larmor.series_table import SeriesTable
                spec.table = SeriesTable.from_dict(table)
            except Exception:                             # noqa: BLE001
                spec.table = None
        spec.refresh_comparison()
        return spec

    # ---- read
    @property
    def n(self) -> int:
        return len(self.members)

    def keys(self) -> list:
        return [m.key for m in self.members]

    def names(self) -> list:
        return [m.name for m in self.members]

    def paths(self) -> list:
        return [m.source_path for m in self.members]

    def index_of(self, key) -> int | None:
        for k, m in enumerate(self.members):
            if m.key == key:
                return k
        return None

    def member(self, key) -> SeriesMember | None:
        k = self.index_of(key)
        return None if k is None else self.members[k]

    def tag_for(self, k: int) -> dict:
        """The dict member ``k``'s workspace carries (JSON-safe). The Series
        table rides on the first member only."""
        m = self.members[k]
        tag = {"id": self.id, "key": m.key, "index": int(k), "name": m.name,
               "group": m.group, "folder": m.folder, "title": m.title,
               "proc": m.proc, "source_path": m.source_path,
               "options": self.options.to_dict(), "fit_sig": m.fit_sig,
               "failed": bool(m.failed), "locked": bool(m.locked)}
        if k == 0 and self.table is not None and hasattr(self.table, "to_dict"):
            tag["table"] = self.table.to_dict()
        return tag

    def carry_for(self, recipes) -> tuple:
        """The parameters carried between neighbours: every name present in
        ``recipes`` except amplitude, minus the names unticked in Carry."""
        off = set(self.options.carry_off)
        return tuple(n for n in carry_candidates(recipes) if n not in off)

    def set_carry(self, names_on, candidates) -> None:
        """Record the Carry checklist: ``names_on`` ticked among
        ``candidates``; an unticked name that is not a candidate right now
        stays remembered."""
        cand = set(candidates)
        on = set(names_on)
        kept = {n for n in self.options.carry_off if n not in cand}
        self.options.carry_off = tuple(sorted(kept | (cand - on)))

    def refresh_comparison(self) -> None:
        from larmor import comparability
        try:
            self.comparison = comparability.compare(
                [m.params for m in self.members], self.names())
        except Exception:                                 # noqa: BLE001
            self.comparison = None

    def comparability_flag(self, k: int) -> str:
        """'LB 100 Hz; TDeff 768' when member ``k`` was acquired / processed
        unlike the series majority, '' otherwise."""
        cmp = self.comparison
        if cmp is None:
            return ""
        try:
            return cmp.signature(k)
        except Exception:                                 # noqa: BLE001
            return ""

    # ---- the Series table
    def series_table(self):
        """A ``SeriesTable`` over the members in order: the edited one when
        it still pairs with the members (composition columns kept), else a
        fresh one from the member records."""
        from larmor.series_table import SeriesTable
        fresh = SeriesTable.from_paths(
            self.paths(), names=self.names(), groups=[m.group for m in self.members],
            folders=[m.folder for m in self.members],
            titles=[m.title for m in self.members])
        if self.table is None:
            return fresh
        try:
            return self.table.aligned_to(self.paths(), fresh)
        except Exception:                                 # noqa: BLE001
            return fresh

    def apply_table(self, table, perm) -> bool:
        """Adopt an edited Series table: ``perm`` (new position -> current
        index) reorders the members, the rows' names and groups are copied
        onto them, the table is kept for the plot's x axis. Returns whether
        the order changed."""
        n = self.n
        perm = [int(k) for k in perm]
        if sorted(perm) != list(range(n)):
            raise ValueError("the permutation does not cover the series")
        moved = perm != list(range(n))
        self.members = [self.members[k] for k in perm]
        rows = list(getattr(table, "rows", []) or [])
        for m, row in zip(self.members, rows):
            m.name = str(row.display_name)
            m.group = str(row.group or row.display_name)
        self.table = table
        self.refresh_comparison()
        return moved

    def remove(self, key) -> bool:
        k = self.index_of(key)
        if k is None:
            return False
        del self.members[k]
        self.refresh_comparison()
        return True


# ------------------------------------------------------------- carry rules
def carry_into(dst_recipe, src_recipe, carry, dst_amp_max: float,
               src_amp_max: float, src_name: str = "") -> tuple:
    """Carry ``src_recipe``'s model into ``dst_recipe`` (dicts; neither is
    mutated). Returns ``(new_dst_recipe, note)``.

    * ``dst`` has NO lines: the source's lines are copied -- structure,
      values, labels, families and the constraints by letter stay valid
      because the structure is identical -- with every free amplitude
      scaled by ``dst_amp_max / src_amp_max`` (1 when either is unknown);
      ``carry`` plays no part (a copy is whole), the destination's fit
      window is taken from the source when it has none, and the copy is
      unfitted (no ``fit_rmsd``).
    * ``dst`` HAS lines: its structure is kept and :func:`seqfit.seed_from`
      semantics apply for ``carry`` (``None`` = every parameter except
      amplitude) on the lines that pair (:func:`larmor.components.pair_sites`:
      by label, by index only where a label cannot identify its line),
      clipped to the destination's own bounds; linked parameters follow
      their master. A destination line without a counterpart -- a component
      the source does not have -- is left exactly as it is, and the note
      names it.
    """
    dst = copy.deepcopy(dst_recipe) if dst_recipe else {}
    src = src_recipe or {}
    src_sites = list(src.get("sites") or [])
    who = f" from {src_name}" if src_name else ""
    if not src_sites:
        return dst, "nothing to carry — the source spectrum has no lines"
    if not dst.get("sites"):
        sites = copy.deepcopy(src_sites)
        scale = 1.0
        try:
            d, s = float(dst_amp_max), float(src_amp_max)
            if d > 0 and s > 0 and np.isfinite(d) and np.isfinite(s):
                scale = d / s
        except (TypeError, ValueError):
            scale = 1.0
        for site in sites:
            for pn, p in (site.get("params") or {}).items():
                if isinstance(p, dict):
                    p["stderr"] = None
                    if pn == "amplitude" and not p.get("expr") and scale != 1.0:
                        v = float(p.get("value", 0.0) or 0.0) * scale
                        lo, hi = p.get("min"), p.get("max")
                        if lo is not None and np.isfinite(lo):
                            v = max(v, float(lo))
                        if hi is not None and np.isfinite(hi):
                            v = min(v, float(hi))
                        p["value"] = v
        dst["sites"] = sites
        if not dst.get("fit_window_ppm") and src.get("fit_window_ppm"):
            dst["fit_window_ppm"] = list(src["fit_window_ppm"])
        dst.pop("fit_rmsd", None)
        n = len(sites)
        note = f"{n} line{'s' if n != 1 else ''} copied{who}"
        if abs(scale - 1.0) > 1e-9:
            note += f", amplitudes scaled ×{scale:.2f}"
        return dst, note

    names = (carry_candidates([src]) if carry is None
             else tuple(str(n) for n in carry))
    pairs = pair_sites(dst["sites"], src_sites)
    matched_sites = 0
    matched_names: set = set()
    for i, site in enumerate(dst["sites"]):
        j = pairs[i]
        if j is None:
            continue
        sparams = src_sites[j].get("params") or {}
        hit = False
        for pn, p in (site.get("params") or {}).items():
            if not isinstance(p, dict) or p.get("expr"):
                continue                               # linked: follows its master
            if pn not in names:
                continue
            sp = sparams.get(pn)
            if not isinstance(sp, dict) or sp.get("value") is None:
                continue
            v = float(sp["value"])
            lo, hi = p.get("min"), p.get("max")
            if lo is not None and np.isfinite(lo):
                v = max(v, float(lo))
            if hi is not None and np.isfinite(hi):
                v = min(v, float(hi))
            p["value"] = v
            p["stderr"] = None
            matched_names.add(pn)
            hit = True
        matched_sites += int(hit)
    n_dst = len(dst["sites"])
    own = f"this spectrum keeps its own {n_dst} line{'s' if n_dst != 1 else ''}"
    if not matched_sites:
        if all(j is None for j in pairs) and all(usable_labels(dst["sites"])) \
                and all(usable_labels(src_sites)):
            mine = ", ".join(label_of(s) for s in dst["sites"])
            theirs = ", ".join(label_of(s) for s in src_sites)
            return dst, (f"no line here shares a name with {src_name or 'the source'}'s "
                         f"({mine} vs {theirs}) — nothing seeded; {own}. Name the "
                         "same component alike on both (the label column) to "
                         "seed between them")
        return dst, f"no matching parameters to seed{who}; {own}"
    note = (f"{_describe_carried(matched_names)} of {matched_sites} matching "
            f"line{'s' if matched_sites != 1 else ''} seeded{who}; {own}")
    alone = [label_of(s) or f"line {i + 1}" for i, (s, j) in enumerate(zip(dst["sites"], pairs))
             if j is None]
    if alone:
        note += (f" ({', '.join(alone)} {'has' if len(alone) == 1 else 'have'} no "
                 "counterpart there)")
    return dst, note


# ------------------------------------------------------------- statuses
def rmsd_of(recipe):
    """The recipe's ``fit_rmsd`` (LARMOR's normalised definition) as a float,
    None when it was never fitted."""
    if recipe is None:
        return None
    v = recipe.fit_rmsd if isinstance(recipe, Recipe) else recipe.get("fit_rmsd")
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


def member_status(recipe, health_stale: bool = False, failed: bool = False) -> str:
    """'unfitted' (no lines, or lines never fitted) · 'fitted' · 'edited'
    (fitted, then a value / vary / link / zone changed -- ``health_stale`` is
    the fit-health machinery's verdict) · 'failed' (the last fit raised)."""
    if failed:
        return "failed"
    if not _sites(recipe) or rmsd_of(recipe) is None:
        return "unfitted"
    return "edited" if health_stale else "fitted"


# ------------------------------------------------------------- the sweep
def entries_for_sweep(members: list, default_window=None) -> list:
    """``[(Recipe, ppm, amp, window), …]`` for :func:`seqfit.run_sequential`
    from member dicts ``{"recipe", "ppm", "amp", "name"}`` in series order.
    Each member fits in its own ``fit_window_ppm`` (else ``default_window``,
    else the whole spectrum). Every member needs at least one line: the
    error names the ones without."""
    out = []
    missing = []
    for k, m in enumerate(members):
        rec = m.get("recipe") or {}
        if not rec.get("sites"):
            missing.append(str(m.get("name") or f"spectrum {k + 1}"))
            continue
        window = rec.get("fit_window_ppm") or default_window
        if window is not None:
            window = (float(max(window)), float(min(window)))
        out.append((Recipe.from_dict(rec), np.asarray(m.get("ppm"), float),
                    np.asarray(m.get("amp"), float), window))
    if missing:
        raise ValueError(
            "every spectrum of the series needs at least one line before a "
            "sweep — none on " + ", ".join(missing)
            + " (Carry ▾ Copy this model to every spectrum)")
    return out


def _finite(v) -> bool:
    try:
        return bool(np.isfinite(float(v)))
    except (TypeError, ValueError):
        return False


def apply_sweep_result(members: list, result) -> list:
    """Which member gets which fitted recipe, by series order: one entry per
    member -- ``{"recipe": dict, "rmsd": float, "x", "y_fit"}`` for a member
    the sweep fitted in at least one pass, None for one it never reached (a
    Stop keeps what was fitted; a member the sweep never touched keeps its
    own model and status) and None for a kept member (``result.fixed``: it
    seeded its neighbours and was never refitted, so its own fit, errors and
    verdict stand)."""
    n = len(members)
    recs = list(getattr(result, "recipes", []) or [])
    if len(recs) != n:
        raise ValueError(f"{len(recs)} fitted recipes for {n} members")
    hist = list(getattr(result, "history", []) or [])
    rmsd = list(getattr(result, "rmsd", []) or [])
    per = list(getattr(result, "per_dataset", []) or [])
    fixed = set(getattr(result, "fixed", ()) or ())
    out = []
    for k in range(n):
        if k in fixed:
            out.append(None)
            continue
        if hist:
            fitted = any(k < len(h.get("rmsd") or []) and _finite(h["rmsd"][k])
                         for h in hist)
        else:
            fitted = True
        if not fitted:
            out.append(None)
            continue
        rec = recs[k]
        d = rec.to_dict() if hasattr(rec, "to_dict") else dict(rec)
        pd = per[k] if k < len(per) and isinstance(per[k], dict) else {}
        out.append({"recipe": d,
                    "rmsd": float(rmsd[k]) if k < len(rmsd) and _finite(rmsd[k])
                    else float("nan"),
                    "x": pd.get("x"), "y_fit": pd.get("y_fit")})
    return out
