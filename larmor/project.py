"""The ``.larproj.json`` project bundle: schema version, migrations, entries.

A project bundle is the "save the whole session" file: every row of the
Workspaces dock, in order, plus which one was active. Version 1 carried only
1D workspaces (spectrum + processing + fit + overlays by reference). Version 2
tags every entry with a ``kind`` and adds three more:

* ``2d`` -- a 2D map BY REFERENCE (its source path plus a project-relative
  path) together with the replayable list of processing operations the
  contour view recorded (phase, shear, transpose / reverse / symmetrize,
  calibrate), its contour settings and its 1D projection overlays (also by
  reference). Never the arrays: a 2rr is 10^5-10^6 points, four quadrants
  when hypercomplex, and the rule everywhere in this codebase is "never embed
  what you can reload".
* ``figure`` -- a Plotting-studio / Figure-studio spec (the ``.figure.json``
  dict itself).
* ``batch`` -- a batch-fit session: spectrum paths, model, release / baseline
  / error settings and the fitted result (recipes, RMSDs, error detail). The
  per-spectrum model curves are recomputed on reopen and the results table is
  rebuilt from the result (``larmor.batchfit.shared_table`` / ``error_table``),
  so neither is stored.

Qt-free: the desktop app builds the entries from its plain workspace snapshot
dicts, ``larmor.batch`` reads bundles through :func:`load_bundle`, and the
migration chain runs through the same ``recipe.run_migrations`` loop the
recipe schema uses, so both tables migrate through one tested path.
"""
from __future__ import annotations

import json
import os
from collections import Counter
from pathlib import Path
from typing import Callable

import numpy as np

from larmor import recipe as _recipe

#: bundle schema version written by :func:`build_bundle`. v1: 1D workspaces
#: only, no ``kind`` on the entries. v2: every entry carries a ``kind`` and
#: 2D maps, figures and batch-fit sessions are included.
PROJECT_BUNDLE_VERSION = 2

#: entry kinds that are documents (they own the central view when active)
DOC_KINDS = ("1d", "2d")
#: entry kinds that are sessions reopened in their own dialog
SESSION_KINDS = ("figure", "batch")


def _v1_to_v2(d: dict) -> dict:
    """v1 entries were all 1D workspaces and carried no ``kind``. Only an
    existing list is rewritten: a JSON object that is not a bundle (a recipe
    handed to the batch reader) must come out without a phantom empty
    ``workspaces`` list, or it would read as an empty project."""
    d = dict(d)
    if isinstance(d.get("workspaces"), list):
        d["workspaces"] = [{**ws, "kind": ws.get("kind", "1d")}
                           for ws in d["workspaces"]]
    return d


#: {from_version: transform(dict) -> dict}, walked by recipe.run_migrations
_MIGRATIONS: dict[int, Callable[[dict], dict]] = {1: _v1_to_v2}


def migrate_bundle(d: dict) -> tuple[dict, list[str]]:
    """Bring a bundle dict up to :data:`PROJECT_BUNDLE_VERSION`.

    A missing or non-integer version means v1 (the first release wrote the
    field but never read it). A bundle written by a NEWER LARMOR is opened
    anyway with a note; whatever this build does not understand is ignored by
    the reader, never refused. Returns ``(bundle, notes)``."""
    v = d.get("larmor_project_version", 1)
    if not isinstance(v, int) or isinstance(v, bool):
        v = 1
    notes: list[str] = []
    if v < PROJECT_BUNDLE_VERSION:
        d, notes = _recipe.run_migrations(d, v, _MIGRATIONS,
                                          PROJECT_BUNDLE_VERSION, "project")
        d["larmor_project_version"] = v + len(notes)
    elif v > PROJECT_BUNDLE_VERSION:
        notes.append(
            f"this project was saved by a newer LARMOR (bundle v{v}, this "
            f"build reads v{PROJECT_BUNDLE_VERSION}) — unknown content ignored")
    return d, notes


def load_bundle(path: str | Path) -> tuple[dict, list[str]]:
    """Read and migrate a ``.larproj.json``. Raises on an unreadable file or
    on JSON that is not a bundle -- no object, or no ``workspaces`` list (a
    recipe file, say); the caller reports it or falls back."""
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(d, dict) or not isinstance(d.get("workspaces"), list):
        raise ValueError("not a LARMOR project bundle (no workspaces list)")
    return migrate_bundle(d)


# ------------------------------------------------------------------ paths
def relpath_or_none(path: str | None, project_dir: str | Path) -> str | None:
    """``path`` relative to the project folder with forward slashes, or None
    when no relative path exists (another Windows drive) or ``path`` is
    empty."""
    if not path:
        return None
    try:
        rel = os.path.relpath(str(path), str(project_dir))
    except ValueError:
        return None
    return rel.replace("\\", "/")


def relocate(path: str | None, project_dir: str | Path,
             rel: str | None = None) -> str | None:
    """The first existing of ``path`` (as saved) and ``project_dir/rel`` (the
    project moved together with its data), else None."""
    if path and Path(path).exists():
        return str(path)
    if rel:
        cand = Path(project_dir) / rel
        if cand.exists():
            return str(cand)
    return None


def json_safe(obj):
    """Recursively turn numpy scalars / arrays, sets, tuples and paths into
    plain JSON types (floats, ints, lists, strings)."""
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (set, frozenset)):
        try:
            items = sorted(obj)
        except TypeError:
            items = list(obj)
        return [json_safe(v) for v in items]
    if isinstance(obj, np.ndarray):
        return json_safe(obj.tolist())
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, Path):
        return str(obj)
    return obj


# ------------------------------------------------------------------ entries
def view2d_persisted(state: dict) -> dict:
    """The JSON subset of a ``Contour2DView.get_state()`` dict (live or a
    stored workspace snapshot): the recorded op log, contour settings, display
    mode / colormap / axis flips, and the 1D projection overlays whose source
    path is known (colour and scale ride along; the arrays are reloaded)."""
    hs = state.get("hmqc_source") or {}
    hc = state.get("hmqc_color") or {}
    sc = state.get("hmqc_scale") or {}
    projections = {}
    for axis in ("f2", "f1"):
        src = hs.get(axis) or ""
        if src:
            projections[axis] = {"source": str(src),
                                 "color": str(hc.get(axis) or ""),
                                 "scale": float(sc.get(axis, 1.0))}
    flip = {"f2": True, "f1": False, **(state.get("flip") or {})}
    return {
        "ops": json_safe(list(state.get("ops") or [])),
        "sign": str(state.get("sign", "positive")),
        "nlevels": float(state.get("nlevels", 12)),
        "floor": float(state.get("floor", 8.0)),
        "title": str(state.get("title", "")),
        "disp": str(state.get("disp", "contour")),
        "cmap": str(state.get("cmap", "viridis")),
        "flip": {"f2": bool(flip["f2"]), "f1": bool(flip["f1"])},
        "projections": projections,
    }


def entry_1d(ws: dict) -> dict | None:
    """A 1D workspace: spectrum arrays, processing recipe, fit, and its
    overlays by reference. None when the snapshot holds no spectrum."""
    snap = ws.get("snap") or {}
    if snap.get("exp_ppm") is None:
        return None
    # overlays are saved by REFERENCE (label/color/visible/source), not their
    # embedded arrays -- same "never embed what you can reload" principle as
    # everything else in this codebase (Recipe, batch CSVs, ...); an overlay
    # with no source (e.g. one hand-typed from an in-memory array with nothing
    # on disk) can't be restored and is dropped rather than erroring.
    overlays = [{"label": o.get("label", ""), "color": o.get("color", ""),
                 "visible": bool(o.get("visible", True)),
                 "source": o.get("source", "")}
                for o in snap.get("overlays", []) if o.get("source")]
    return {
        "kind": "1d",
        "title": ws.get("title", ""),
        "source_path": snap.get("source_path"),
        "recipe": snap.get("recipe"),
        "hidden": sorted(snap.get("hidden", set())),
        "exp_ppm": np.asarray(snap["exp_ppm"], float).tolist(),
        "exp_amp": np.asarray(snap["exp_amp"], float).tolist(),
        "overlays": overlays,
    }


def entry_2d(ws: dict, project_dir: str | Path) -> dict | None:
    """A 2D workspace by reference: source path (+ project-relative path),
    recipe, hidden sites, fittability and the persisted view (op log, contour
    settings, projection overlays). None when the source cannot be reloaded
    (no path, or the file is gone) -- the arrays are never embedded."""
    snap = ws.get("snap") or {}
    src = snap.get("source_path") or getattr(snap.get("data2d"), "source", "")
    if not src or not Path(src).exists():
        return None
    return {
        "kind": "2d",
        "title": ws.get("title", ""),
        "source_path": str(src),
        "source_rel": relpath_or_none(str(src), project_dir),
        "recipe": snap.get("recipe"),
        "hidden": sorted(snap.get("hidden", set())),
        "fittable": bool(snap.get("fittable", False)),
        "view": view2d_persisted(snap.get("view2d") or {}),
    }


def entry_figure(ws: dict) -> dict:
    """A figure kept from the Plotting / Figure studio: its spec verbatim."""
    snap = ws.get("snap") or {}
    return {"kind": "figure", "title": ws.get("title", ""),
            "spec": json_safe(snap.get("spec") or {})}


def entry_batch(ws: dict, project_dir: str | Path) -> dict:
    """A batch-fit session: the dialog's ``session_state()`` made JSON-safe,
    plus project-relative spectrum paths so the session follows a project
    moved together with its data."""
    snap = ws.get("snap") or {}
    state = json_safe(snap.get("state") or {})
    state["paths_rel"] = [relpath_or_none(p, project_dir)
                          for p in (state.get("paths") or [])]
    return {"kind": "batch", "title": ws.get("title", ""), "state": state}


def relocate_batch_state(state: dict, project_dir: str | Path
                         ) -> tuple[dict, list[str]]:
    """Relocate a batch session's spectrum paths (absolute first, then the
    project-relative copy) and rewrite every reference to a moved path -- the
    per-spectrum settings keyed by path and each restored recipe's
    ``source_path``, which is what aligns results to reloaded spectra. Paths
    found nowhere are kept as saved (the dialog skips what it cannot read)
    and returned as ``missing``."""
    st = dict(state)
    paths = list(st.get("paths") or [])
    rels = list(st.get("paths_rel") or [])
    mapping: dict[str, str] = {}
    missing: list[str] = []
    new_paths = []
    for i, p in enumerate(paths):
        rel = rels[i] if i < len(rels) else None
        new = relocate(p, project_dir, rel)
        if new is None:
            missing.append(p)
            new_paths.append(p)
            continue
        if new != p:
            mapping[p] = new
        new_paths.append(new)
    st["paths"] = new_paths
    if mapping:
        per = st.get("per_spectrum") or {}
        st["per_spectrum"] = {mapping.get(k, k): v for k, v in per.items()}
        res = st.get("result")
        if isinstance(res, dict):
            res = dict(res)
            recs = []
            for r in res.get("recipes") or []:
                r = dict(r)
                sp = r.get("source_path", "")
                if sp in mapping:
                    r["source_path"] = mapping[sp]
                recs.append(r)
            res["recipes"] = recs
            st["result"] = res
        ser = st.get("series")
        if isinstance(ser, dict):
            # the Series table's rows pair by source_path too (SeriesTable.aligned_to)
            ser = dict(ser)
            ser["rows"] = [{**r, "source_path": mapping.get(r.get("source_path", ""),
                                                            r.get("source_path", ""))}
                           for r in (ser.get("rows") or [])]
            st["series"] = ser
    return st, missing


def build_bundle(workspaces: list[dict], active: int | None,
                 project_dir: str | Path) -> tuple[dict, int]:
    """The bundle dict for a list of workspace entries (the app's
    ``MainWindow.workspaces``) and how many were dropped -- a 1D without a
    spectrum or a 2D without a reloadable source. ``active`` is remapped onto
    the saved list so it survives the drops."""
    entries: list[dict] = []
    remap: dict[int, int] = {}
    dropped = 0
    for i, ws in enumerate(workspaces):
        kind = ws.get("kind", "1d")
        if kind == "1d":
            e = entry_1d(ws)
        elif kind == "2d":
            e = entry_2d(ws, project_dir)
        elif kind == "figure":
            e = entry_figure(ws)
        elif kind == "batch":
            e = entry_batch(ws, project_dir)
        else:
            e = None
        if e is None:
            dropped += 1
            continue
        remap[i] = len(entries)
        entries.append(e)
    act = remap.get(active) if active is not None else None
    from larmor.provenance import software_stamp

    # the software that wrote the bundle (the recipes inside carry the
    # software of their own fits); load_bundle needs only the workspaces list
    return {"larmor_project_version": PROJECT_BUNDLE_VERSION,
            "active": act, "workspaces": entries,
            "software": software_stamp()}, dropped


def summary(bundle: dict) -> str:
    """'3 spectra, 1 2D map, 2 figures, 1 batch session' -- only the kinds
    present, singular / plural."""
    counts = Counter(w.get("kind", "1d") for w in bundle.get("workspaces", []))
    names = (("1d", "spectrum", "spectra"), ("2d", "2D map", "2D maps"),
             ("figure", "figure", "figures"),
             ("batch", "batch session", "batch sessions"))
    parts = [f"{counts[k]} {one if counts[k] == 1 else many}"
             for k, one, many in names if counts.get(k)]
    return ", ".join(parts) or "no workspaces"
