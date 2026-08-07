"""Resolve a batch of already-fitted spectra into panels for a publication
grid: point it at a batch_table*.csv (from the batch-fit dialog's Save
table.../Export CSV...), a folder of saved fits, or an explicit path list, and
it finds the .recipe.json files, auto-matching a CSV's rows to sibling saved
fits by sample name -- so the Plotting studio's "Batch grid" kind can list
what it found, let the user pick/reorder/hide panels, and render each one's
experiment + fit + components without asking the user to re-locate anything.

The actual per-panel data/model resolution at render time reuses figures.py's
existing load_trace({"recipe": path, "part": ...}) machinery (which already
follows a recipe's own source_path, kernel or not, data-less or not) -- this
module's only job is finding and listing which recipes to feed it. Qt-free
and testable; the desktop dialog is only a front-end.

Newer batch-fit CSV exports carry each row's own "source_path"/"model"
columns (see batchfit.shared_table/error_table) -- a per-scope hint this
module uses two ways: as an extra folder to search for sibling .recipe.json
fits (the CSV may live somewhere else entirely), and, failing that, as a
raw-data-only fallback so the sample still shows up as *something* rather
than silently vanishing. Older CSVs (no such columns) fall back to the
original filename/sample matching unchanged -- read_csv_hints() simply
returns {} for them. A scope neither method resolves comes back as a Panel
with needs_manual=True, for the studio to ask the user about directly
(one popup per sample, "locate the data for g3") instead of just warning.
"""
from __future__ import annotations

import csv
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

from larmor.recipe import Param, Recipe, SiteModel

#: fit files this module recognizes when scanning a folder
_FIT_GLOBS = ("*.recipe.json", "*.fxmla", "*.fxml")

#: where CSV-reconstructed recipes get written (one per session; the studio
#: only needs these to exist for as long as it's open) -- lazily created
_recon_dir: Path | None = None


def _reconstructed_recipe_path(scope: str) -> Path:
    global _recon_dir
    if _recon_dir is None:
        _recon_dir = Path(tempfile.mkdtemp(prefix="larmor_csv_recipes_"))
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in scope)
    return _recon_dir / f"{safe}.recipe.json"


@dataclass
class Panel:
    """One resolved spectrum: enough to list in the studio and to reference
    in a saved batch_grid figure spec (which stores `path`, not raw arrays,
    so the spec stays diffable/reloadable like every other LARMOR figure).

    `path` is "" when no .recipe.json was found for this sample -- either
    `data_path` names a raw spectrum with no fit yet (drawn experiment-only),
    or, if that's empty too, `needs_manual` is True and the studio should
    prompt the user to locate this sample's data itself."""
    path: str
    sample: str
    nucleus: str
    models: tuple[str, ...]         # each site's model, e.g. ("gauss_lor",)*7
    has_data: bool                  # False = data-less fit (model-only preview)
    n_sites: int = 0
    data_path: str = ""             # raw spectrum location when there's no fit
    needs_manual: bool = False      # True: ask the user to locate the data
    reconstructed: bool = False     # True: `path` is a fit rebuilt from the CSV's own rows, not a saved file


def read_csv_hints(csv_path: str | Path) -> dict[str, dict]:
    """Per-scope hints from a batch CSV's own "source_path"/"model" columns
    (added to batch-fit exports alongside this feature) -- {} for an older,
    column-less CSV, so callers degrade to filename/sample matching alone.
    A scope's source_path is its LAST non-empty value across its rows (every
    per-spectrum row carries the same one); models are collected in the order
    first seen, one per site."""
    hints: dict[str, dict] = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "source_path" not in reader.fieldnames:
            return {}
        for row in reader:
            scope = (row.get("scope") or "").strip()
            if not scope or scope.lower() == "shared":
                continue
            h = hints.setdefault(scope, {"source_path": "", "models": []})
            sp = (row.get("source_path") or "").strip()
            if sp:
                h["source_path"] = sp
            m = (row.get("model") or "").strip()
            if m and m not in h["models"]:
                h["models"].append(m)
    return {k: {"source_path": v["source_path"], "models": tuple(v["models"])}
            for k, v in hints.items()}


def recipe_from_csv_rows(shared_rows: list[dict], scope_rows: list[dict],
                         source_path: str = "") -> Recipe:
    """Rebuild one spectrum's fitted Recipe purely from a batch CSV's own
    rows -- a shared-ladder batch fit spreads one site's parameters across
    two row groups (params held fixed across every spectrum land under
    scope "shared"; params released or always-free, like amplitude, land
    under the spectrum's own scope), so a full per-site Param set needs
    both merged. Needs each row's `model` (only present in CSVs exported
    after this feature); raises ValueError otherwise so callers fall back
    to a real saved .recipe.json or a data-only panel instead of a bad guess.

    `source_path`, when it exists on disk, seeds nucleus/larmor_frequency_MHz
    /spin_rate_Hz from the real spectrum (via loader.load_any) so the
    reconstructed recipe is a proper, axis-correct Recipe, not a units-blind
    shell -- exactly what render_batch_grid needs to draw experiment+fit
    together like any other saved fit.
    """
    by_site: dict[str, dict] = {}
    for r in list(shared_rows) + list(scope_rows):
        model = (r.get("model") or "").strip()
        if not model:
            raise ValueError(
                "this CSV has no 'model' column (an older export) -- can't "
                "rebuild a fit from it; use a folder of saved .recipe.json "
                "fits instead")
        site = by_site.setdefault(r["site"], {"label": r.get("label", ""),
                                              "model": model, "params": {}})
        try:
            value = float(r["value"])
        except (TypeError, ValueError):
            continue
        stderr = None
        raw_stderr = r.get("stderr")
        if raw_stderr not in (None, ""):
            try:
                stderr = float(raw_stderr)
            except ValueError:
                stderr = None
        site["params"][r["param"]] = Param(value, stderr=stderr)

    if not by_site:
        raise ValueError("no fittable rows for this scope")

    from larmor import models

    def _site_key(k: str):
        return int(k[1:]) if k[1:].isdigit() else k

    sites = []
    for k in sorted(by_site, key=_site_key):
        s = by_site[k]
        try:
            needed = models.get(s["model"]).param_names
        except KeyError:
            raise ValueError(f"unknown model {s['model']!r} for site {k}") from None
        missing = [n for n in needed if n not in s["params"]]
        if missing:
            raise ValueError(
                f"site {k} ({s['label'] or k}) is missing {', '.join(missing)} "
                "in this CSV -- can't rebuild a complete fit from it")
        sites.append(SiteModel(model=s["model"], label=s["label"], params=s["params"]))

    nucleus, larmor_MHz, spin_rate_Hz = "", 0.0, 0.0
    if source_path and Path(source_path).exists():
        try:
            from larmor.loader import load_any
            _, _, rec_dict, _, _ = load_any(source_path)
            nucleus = rec_dict.get("nucleus", "")
            larmor_MHz = float(rec_dict.get("larmor_frequency_MHz", 0.0) or 0.0)
            spin_rate_Hz = float(rec_dict.get("spin_rate_Hz", 0.0) or 0.0)
        except Exception:
            pass    # still usable data-less/unit-less; render falls back too

    return Recipe(nucleus=nucleus, larmor_frequency_MHz=larmor_MHz,
                  spin_rate_Hz=spin_rate_Hz, source_path=source_path, sites=sites)


def find_recipes_near_csv(csv_path: str | Path) -> list[str]:
    """Match a batch_table*.csv's rows to sibling .recipe.json fits by sample
    name -- the natural pairing when "Save individual fits..." and "Save
    table.../Export CSV..." came from the same batch session. Searches the
    CSV's own folder plus, when the CSV carries source_path hints, each
    hinted spectrum's folder too (the CSV and the fits needn't live in the
    same place). Returns what it could match, in the CSV's own row order;
    scopes without a matching file are simply absent."""
    folder = Path(csv_path).parent
    scopes: list[str] = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            s = (row.get("scope") or "").strip()
            if s and s.lower() != "shared" and s not in scopes:
                scopes.append(s)

    hints = read_csv_hints(csv_path)
    search_dirs = {folder}
    for h in hints.values():
        sp = h.get("source_path")
        if sp:
            search_dirs.add(Path(sp).parent)

    candidates: set[Path] = set()
    for d in search_dirs:
        candidates |= set(d.glob("*.recipe.json")) | set(d.glob("*/*.recipe.json"))
    candidates = sorted(candidates)

    matched: list[str] = []
    by_sample: dict[str, Path] | None = None    # built lazily, only if needed
    for scope in scopes:
        hit = next((c for c in candidates if scope in c.stem), None)
        if hit is None:
            if by_sample is None:
                by_sample = {}
                for c in candidates:
                    try:
                        by_sample[Recipe.load(c).sample] = c
                    except Exception:
                        continue
            hit = by_sample.get(scope)
        if hit is not None and str(hit) not in matched:
            matched.append(str(hit))
    return matched


def resolve_paths(source) -> tuple[list[str], list[str]]:
    """`source` is a batch CSV path, a folder, or a list of fit paths --
    normalize to an ordered list of fit-file paths + any warnings."""
    warnings: list[str] = []
    if isinstance(source, (list, tuple)):
        return [str(p) for p in source], warnings

    p = Path(source)
    if p.is_dir():
        paths: list[str] = []
        for pat in _FIT_GLOBS:
            paths += sorted(str(f) for f in p.glob(pat))
        if not paths:
            warnings.append(f"no fit files (.recipe.json/.fxmla) found in {p}")
        return paths, warnings

    if p.suffix.lower() == ".csv":
        paths = find_recipes_near_csv(p)
        if not paths:
            warnings.append(
                f"no .recipe.json fits found next to {p.name} matching its "
                "rows -- population-vs-composition plots can still use the "
                "CSV's own values directly, but a deconvolution grid needs "
                "the saved fits (batch fit's \"Save individual fits…\")")
        return paths, warnings

    return [str(p)], warnings


def load_panels(source) -> tuple[list[Panel], list[str]]:
    """Resolve `source` to Panels for the studio's spectrum list (checkable,
    reorderable) -- cheap (reads each recipe's header, not its data).

    For a CSV source, every scope it mentions gets a Panel, even when no
    .recipe.json matched, in decreasing order of fidelity: (1) a real saved
    .recipe.json (sample-name or source_path-folder matched); (2) failing
    that, a full fit REBUILT from the CSV's own rows (recipe_from_csv_rows,
    needs the `model` column -- most batch-CSV loads land here, since
    "Save individual fits…" is easy to skip); (3) failing that (an older
    CSV with no model column), a data-only panel from its source_path hint
    (raw spectrum, no fit); (4) failing even that, needs_manual=True rather
    than being dropped, so the studio can ask the user to locate it directly
    instead of the plot silently missing a sample."""
    paths, warnings = resolve_paths(source)
    panels: list[Panel] = []
    matched_samples: set[str] = set()
    for p in paths:
        try:
            rec = Recipe.load(p)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{Path(p).name}: could not load ({exc}) — skipped")
            continue
        sample = rec.sample or Path(p).stem
        panels.append(Panel(
            path=str(p), sample=sample, nucleus=rec.nucleus,
            models=tuple(s.model for s in rec.sites),
            has_data=bool(rec.source_path and Path(rec.source_path).exists()),
            n_sites=len(rec.sites)))
        matched_samples.add(sample)

    src_path = Path(source) if isinstance(source, (str, Path)) else None
    if src_path is not None and src_path.suffix.lower() == ".csv" and src_path.exists():
        hints = read_csv_hints(src_path)
        rows_by_scope = csv_rows_by_scope(src_path)
        shared_rows = rows_by_scope.get("shared", [])
        scopes = [s for s in rows_by_scope if s.lower() != "shared"]
        unresolved: list[str] = []
        for scope in scopes:
            if not scope or scope in matched_samples:
                continue
            hint = hints.get(scope, {})
            sp = hint.get("source_path", "")
            sp = sp if sp and Path(sp).exists() else ""
            try:
                rec = recipe_from_csv_rows(shared_rows, rows_by_scope[scope], sp)
                tmp = _reconstructed_recipe_path(scope)
                rec.sample = scope
                rec.save(tmp)
                panels.append(Panel(path=str(tmp), sample=scope, nucleus=rec.nucleus,
                                    models=tuple(s.model for s in rec.sites),
                                    has_data=bool(sp), n_sites=len(rec.sites),
                                    data_path=sp, reconstructed=True))
                continue
            except ValueError:
                pass    # no model column (an older CSV) -- fall through below
            if sp:
                panels.append(Panel(path="", sample=scope, nucleus="",
                                    models=hint.get("models", ()), has_data=True,
                                    data_path=sp))
            else:
                panels.append(Panel(path="", sample=scope, nucleus="",
                                    models=hint.get("models", ()), has_data=False,
                                    needs_manual=True))
                unresolved.append(scope)
        if unresolved:
            warnings.append(
                "couldn't locate data for: " + ", ".join(unresolved) +
                " — locate them manually to include in the plot")
    return panels, warnings


def resolve_manual(panel: Panel, data_path: str | Path) -> Panel:
    """Fold a user's manually-picked file (the "successive popups" fallback,
    one per unresolved sample) back into its Panel -- has_data flips True and
    needs_manual clears, so the studio's list can drop the "locate…" prompt
    for that row and the panel renders experiment-only, same as a CSV-hinted
    data-only match."""
    return replace(panel, data_path=str(data_path), has_data=True, needs_manual=False)


def csv_rows_by_scope(csv_path: str | Path) -> dict[str, list[dict]]:
    """The raw rows of a batch_table*.csv grouped by "scope" (sample) --
    for the composition-trend/species-distribution templates, which plot the
    CSV's own values directly and don't need the saved fits at all."""
    out: dict[str, list[dict]] = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            out.setdefault(row.get("scope", ""), []).append(row)
    return out
