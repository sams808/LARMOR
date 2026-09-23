"""Series publication bundle: everything about a Batch fit / Sequential fit
in one folder, written once, ready for a notebook or a paper's SI.

Per spectrum ``NN_<sample>_curves.csv`` (ppm, the experiment EXACTLY as
fitted -- after the recipe's processing and any per-spectrum baseline --
model, residual and one column per component on the experimental axis) and
``NN_<sample>.recipe.json`` (a copy with source_path and source_sha256
filled); for the series ``manifest.csv`` (source file + SHA-256, EXPNO,
NS, D1, SF/SR, processing, fit window, a recomputed normalised RMSD, error
method, versions), the long tables (``<kind>_table.csv`` and the error
table when one was computed) and ``README.txt`` (software versions, the
Methods paragraph, conventions, a file glossary).

Distinct from :mod:`larmor.project`'s ``.larproj`` bundle (a reopenable
session) and from the Report dock's single-fit publication bundle
(``mw_fitting.export_publication_bundle``: methods.txt, table.tex,
figure.*, report.md) -- this one writes README.txt and uniquely prefixed
files, so both can share a folder.

Qt-free: the two dialogs and the CLI only call :func:`write_bundle`.
"""
from __future__ import annotations

import copy
import csv
import datetime as _dt
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from larmor import batchfit
from larmor.recipe import Recipe, sha256_of

#: manifest.csv's fixed column contract, one row per spectrum in series order
MANIFEST_COLUMNS: tuple[str, ...] = (
    "index", "stem", "name", "recipe_file", "curves_file", "source_path",
    "data_file", "source_kind", "source_sha256", "expno", "procno", "nucleus",
    "larmor_MHz", "spin_rate_Hz", "bf1_MHz", "sf_MHz", "sr_hz", "ns", "d1_s",
    "pulprog", "date", "title", "processing", "fit_window_hi_ppm",
    "fit_window_lo_ppm", "n_points_window", "rmsd", "n_sites", "excluded_sites",
    "released", "error_method", "larmor_version", "mrsimulator_version", "note")

#: the acquisition facts read from a Bruker EXPNO ('' for anything else)
ACQUISITION_KEYS: tuple[str, ...] = ("expno", "procno", "bf1_MHz", "sf_MHz", "sr_hz",
                                     "ns", "d1_s", "pulprog", "date", "title")
#: everything source_facts() returns
SOURCE_KEYS: tuple[str, ...] = ("data_file", "source_kind", "source_sha256",
                                *ACQUISITION_KEYS)

_TEXT_SUFFIXES = (".csv", ".txt", ".dat")
_FXMLA_SUFFIXES = (".fxmla", ".fxml")


# ---------------------------------------------------------------- names
def _slug(s: str) -> str:
    """The batch dialog / CLI file-name rule: letters, digits, '-' and '_'
    kept, anything else '_', capped at 60 characters, 'fit' when empty."""
    return ("".join(c if c.isalnum() or c in "-_" else "_" for c in (s or ""))[:60]
            or "fit")


def stems(labels: list[str], width: int | None = None) -> list[str]:
    """``'01_Base0Ca'``: the zero-padded series position (``width`` digits,
    default ``max(2, len(str(n)))``) + the sample slug -- directory order
    equals series order and identically titled spectra stay distinct.
    Always unique (case-insensitively, for Windows): a collision gets
    ``_2``, ``_3`` …"""
    n = len(labels)
    width = width or max(2, len(str(n)))
    out: list[str] = []
    seen: set[str] = set()
    for k, lab in enumerate(labels):
        base = f"{k + 1:0{width}d}_{_slug(str(lab))}"
        s, j = base, 2
        while s.lower() in seen:
            s = f"{base}_{j}"
            j += 1
        seen.add(s.lower())
        out.append(s)
    return out


# ---------------------------------------------------------------- sources
def _deref(source_path: str) -> str:
    """A ``.recipe.json`` source stands for its own source_path (one level)."""
    p = Path(str(source_path))
    if p.suffix.lower() == ".json" and p.is_file():
        try:
            return Recipe.load(p).source_path or ""
        except Exception:
            return ""
    return str(source_path)


def data_file(source_path: str) -> tuple[Path | None, str]:
    """``(the file actually read, source_kind)`` for a spectrum's source:
    any Bruker path (EXPNO dir, pdata/N, 1r…) -> ``expno/pdata/<procno>/
    <target>`` and ``'bruker'``; a .csv/.txt/.dat file itself and ``'csv'``;
    .fxmla and ``'fxmla'``; a .recipe.json -> its own source (once);
    ``(None, '')`` when unresolvable or missing. Never raises."""
    try:
        sp = _deref(source_path) if source_path else ""
        if not sp:
            return None, ""
        p = Path(sp)
        suf = p.suffix.lower()
        if suf in _TEXT_SUFFIXES:
            return (p, "csv") if p.is_file() else (None, "")
        if suf in _FXMLA_SUFFIXES:
            return (p, "fxmla") if p.is_file() else (None, "")
        if suf == ".json":                      # a recipe pointing at a recipe
            return None, ""
        from larmor.io import bruker

        try:
            ref = bruker.resolve(p)
        except (ValueError, FileNotFoundError):
            ref = None
        if ref is not None:
            if ref.target in ("1r", "2rr"):
                f = ref.expno / "pdata" / str(ref.procno) / ref.target
            else:
                f = ref.expno / ref.target
            return (f, "bruker") if f.is_file() else (None, "")
        if p.is_file():
            return p, (suf.lstrip(".") or "file")
        return None, ""
    except Exception:
        return None, ""


def acquisition_fields(source_path: str) -> dict:
    """expno, procno, bf1_MHz, sf_MHz, sr_hz, ns, d1_s, pulprog, date (ISO),
    title (first line) -- ``referencing.read_acquisition`` on the EXPNO
    ``bruker.resolve`` finds. Every key present; ``''`` for a non-Bruker
    source or on any failure. Read-only."""
    blank = {k: "" for k in ACQUISITION_KEYS}
    try:
        sp = _deref(source_path) if source_path else ""
        if not sp:
            return blank
        from larmor import referencing
        from larmor.io import bruker

        ref = bruker.resolve(Path(sp))
        acq = referencing.read_acquisition(ref.expno, ref.procno)
        if acq is None:
            return blank
        return {"expno": str(acq.expno or ref.expno.name),
                "procno": str(acq.procno),
                "bf1_MHz": acq.bf1_MHz,
                "sf_MHz": "" if acq.sf_MHz is None else acq.sf_MHz,
                "sr_hz": "" if acq.sr_hz is None else acq.sr_hz,
                "ns": acq.ns, "d1_s": acq.d1_s,
                "pulprog": acq.pulprog, "date": acq.date_iso,
                "title": acq.title.splitlines()[0].strip() if acq.title else ""}
    except Exception:
        return blank


def source_facts(source_path: str) -> dict:
    """data_file + source_kind + source_sha256 (of the file actually read) +
    :func:`acquisition_fields`. Blanks on failure, never raises."""
    out = {k: "" for k in SOURCE_KEYS}
    f, kind = data_file(source_path)
    if f is None:
        return out
    out["data_file"] = str(f)
    out["source_kind"] = kind
    try:
        out["source_sha256"] = sha256_of(f)
    except Exception:
        pass
    if kind == "bruker":
        out.update(acquisition_fields(source_path))
    return out


# ---------------------------------------------------------------- numbers
def curve_rmsd(ppm, experiment, model, window) -> tuple[float, int]:
    """``fit.fit``'s definition: sqrt(mean((model - experiment)^2)) over the
    fit window (inclusive bounds; the whole axis when ``window`` is None or
    selects nothing), divided by the experiment's maximum in the window (or
    1). Returns ``(rmsd, n_points_window)``."""
    ppm = np.asarray(ppm, float)
    y = np.asarray(experiment, float)
    m = np.asarray(model, float)
    if window:
        hi, lo = max(window), min(window)
        sel = (ppm >= lo) & (ppm <= hi)
        if not sel.any():
            sel = np.ones(ppm.shape, dtype=bool)
    else:
        sel = np.ones(ppm.shape, dtype=bool)
    yw, mw = y[sel], m[sel]
    rmsd = float(np.sqrt(np.mean((mw - yw) ** 2)) / (yw.max() or 1.0))
    return rmsd, int(sel.sum())


# ---------------------------------------------------------------- result
@dataclass
class BundleResult:
    outdir: str
    files: list[str] = field(default_factory=list)
    manifest: list[dict] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    kind: str = "batch"

    @property
    def summary(self) -> str:
        n = len(self.manifest)
        curves = sum(1 for f in self.files if f.endswith("_curves.csv"))
        recipes = sum(1 for f in self.files if f.endswith(".recipe.json"))
        rest = [f for f in self.files
                if not (f.endswith("_curves.csv") or f.endswith(".recipe.json"))]
        parts = [f"{curves} curves", f"{recipes} recipes", *rest]
        return (f"publication bundle: {n} spectra → " + " + ".join(parts)
                + f" in {self.outdir}")


# ---------------------------------------------------------------- README
def readme_text(result, kind: str, files_glossary: list[str],
                error_method: str, warnings: list[str] | None = None) -> str:
    """The README body: title, timestamp, the software block
    (``methods.software_versions``), the Methods paragraph
    (``methods.methods_sentence`` on the first recipe), the conventions
    (RMSD, populations, number precision, excluded sites, how to read the
    files) and the file glossary."""
    from larmor import methods

    recs = list(result.recipes)
    n = len(recs)
    what = "sequential" if kind == "seq" else "batch"
    v = methods.software_versions()
    lines = [f"LARMOR publication bundle — {what} fit of {n} spectra",
             f"written {_dt.datetime.now().strftime('%Y-%m-%d %H:%M')}", ""]
    lines += ["SOFTWARE",
              f"  LARMOR {v['larmor']}" + (f" (git {v['git_commit']})" if v["git_commit"] else ""),
              f"  mrsimulator {v['mrsimulator']}",
              f"  lmfit {v['lmfit']}",
              f"  numpy {v['numpy']}",
              f"  Python {v['python']}", ""]
    lines += ["METHODS"]
    if recs:
        lines.append(methods.methods_sentence(recs[0].to_dict(), error_method))
    shared = tuple(getattr(result, "shared", ()) or ())
    released = tuple(getattr(result, "released", ()) or ())
    if kind == "seq":
        lines.append("Each spectrum was fitted independently, warm-started from its "
                     "neighbour in the series (sequential forward–backward fit); "
                     "every parameter is per spectrum.")
    else:
        frac = float(getattr(result, "release_frac", 0.0) or 0.0)
        s_txt = ", ".join(shared) if shared else "none"
        r_txt = (", ".join(released) + f" (±{frac:.0%} around the shared value)"
                 if released else "none")
        lines.append(f"One model was fitted to every spectrum with the amplitudes "
                     f"free per spectrum. Parameters shared across the series: "
                     f"{s_txt}. Released per spectrum: {r_txt}.")
    lines.append("")
    lines += ["CONVENTIONS",
              "  RMSD = sqrt(mean((model - experiment)^2)) over the fit window, divided by "
              "the maximum of the experiment in the window (LARMOR's definition, fit.fit); "
              "manifest.csv recomputes it from the columns of each _curves.csv.",
              "  population_pct rows in the tables are each component's integral over the "
              "fit window, in percent of the sum over the components present.",
              "  In _curves.csv the ppm, experiment and experiment_raw columns are written "
              "with Python repr() (bit-exact, the arrays the fit saw); model, residual and "
              "the component columns with 9 significant digits.",
              "  experiment is the array as fitted (after the recipe's processing steps and "
              "any per-spectrum baseline); experiment_raw, when present, is the trace "
              "before that baseline. The processing column of manifest.csv (and the "
              "'# processing=' header line) is the exact step list, replayable with "
              "larmor.loader.apply_processing.",
              "  An excluded component (amplitude locked at zero for that spectrum) stays "
              "as a zero column in _curves.csv, is named in that file's '# excluded=' "
              "header and in manifest.csv (excluded_sites), and is left out of the tables.",
              "  Read _curves.csv and manifest.csv with pandas.read_csv(path, comment='#'); "
              "a _curves.csv also opens in LARMOR (File > Open) as the fitted spectrum.",
              ""]
    lines += ["FILES", *[f"  {g}" for g in files_glossary], ""]
    if warnings:
        lines += ["WARNINGS", *[f"  {w}" for w in warnings], ""]
    return "\n".join(lines)


# ---------------------------------------------------------------- writer
def write_bundle(result, data: list[tuple], outdir, *, kind: str = "batch",
                 names: list[str] | None = None,
                 source_paths: list[str] | None = None,
                 raw: list | None = None, error_method: str | None = None,
                 write_recipes: bool = True, write_tables: bool = True,
                 ) -> BundleResult:
    """Write the bundle into ``outdir`` (created if needed).

    ``result``: a ``BatchFitResult`` (or any object with recipes / labels /
    rmsd / shared / released / error_method / error_detail); ``data``:
    ``[(ppm, amp, window), …]`` aligned with ``result.recipes`` -- the arrays
    the fit saw, the list ``batchfit.batch_error_analysis`` takes.
    ``names``: the stems to use (default :func:`stems` of the labels);
    ``source_paths``: per-spectrum source to record and to fill into the
    recipe copies when the recipe lacks one; ``raw``: per-spectrum
    pre-baseline arrays (``experiment_raw`` is written where ``raw[k]`` is
    given and differs from ``amp``); ``error_method``: the error table to
    write as ``<kind>_fit_<method>.csv`` when that method is present in
    ``result.error_detail`` (default ``result.error_method``); never
    triggers an error computation. ``write_recipes=False`` references the
    ``<stem>.recipe.json`` a caller already wrote (the CLI);
    ``write_tables=False`` skips ``<kind>_table.csv`` likewise.

    One bad spectrum never aborts the bundle: its manifest row carries the
    message in ``note`` and the warning is collected. Source facts never
    raise (blank cells + a warning)."""
    from larmor import methods
    from larmor.io import export

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    recs: list[Recipe] = list(result.recipes)
    n = len(recs)
    labels = list(getattr(result, "labels", []) or [])
    if len(labels) != n:
        labels = [r.sample or f"spectrum{k + 1}" for k, r in enumerate(recs)]
    names = list(names) if names else stems(labels)
    if len(names) != n:
        raise ValueError(f"{len(names)} names for {n} spectra")
    if len(data) != n:
        raise ValueError(f"{len(data)} data entries for {n} spectra")

    versions = methods.software_versions()
    detail = getattr(result, "error_detail", {}) or {}
    em_wanted = error_method if error_method is not None else \
        getattr(result, "error_method", None)
    em_file = em_wanted if em_wanted and em_wanted in detail else None
    em_manifest = em_file or (getattr(result, "error_method", "") or "")
    released = ",".join(getattr(result, "released", ()) or ())

    res = BundleResult(outdir=str(outdir), kind=kind)
    copies: list[Recipe] = []
    for k, rec in enumerate(recs):
        ppm, amp, window = data[k]
        ppm = np.asarray(ppm, float)
        amp = np.asarray(amp, float)
        sp = ""
        if source_paths and k < len(source_paths) and source_paths[k]:
            sp = str(source_paths[k])
        sp = sp or (rec.source_path or "")
        facts = source_facts(sp) if sp else {key: "" for key in SOURCE_KEYS}
        if sp and not facts["data_file"]:
            res.warnings.append(f"{names[k]}: source not found or unreadable, "
                                f"no hash recorded: {sp}")

        # the recipe copy: source_path / source_sha256 filled, the in-memory
        # result untouched
        rc = Recipe.from_dict(rec.to_dict())
        rc.source_path = sp
        if facts["source_sha256"]:
            rc.source_sha256 = facts["source_sha256"]
        if facts["source_kind"] and not rc.source_kind:
            rc.source_kind = facts["source_kind"]
        copies.append(rc)
        recipe_file = f"{names[k]}.recipe.json"
        if write_recipes:
            rc.save(outdir / recipe_file)
            res.files.append(recipe_file)
        elif not (outdir / recipe_file).exists():
            recipe_file = ""

        excluded = [f"s{i}_{s.label or f'site{i}'}" for i, s in enumerate(rec.sites)
                    if batchfit.is_zeroed_out(s.params.get("amplitude"))]
        win = rec.fit_window_ppm or window
        hi = lo = ""
        if win:
            hi, lo = max(win), min(win)
        processing = json.dumps(list(rec.processing or []))
        header = {"source_path": sp}
        if win:
            header["fit_window_ppm"] = f"{hi},{lo}"
        header["processing"] = processing
        if excluded:
            header["excluded"] = ",".join(excluded)
        raw_k = None
        if raw and k < len(raw) and raw[k] is not None:
            raw_k = np.asarray(raw[k], float)
            if raw_k.shape != amp.shape or np.array_equal(raw_k, amp):
                raw_k = None

        curves_file = f"{names[k]}_curves.csv"
        note = "" if rec.fit_rmsd is not None else "not fitted"
        rmsd: float | str = ""
        npts: int | str = ""
        try:
            out = export.export_curves_csv(rec, ppm, amp, outdir / curves_file,
                                           raw=raw_k, header=header)
            res.files.append(curves_file)
            r, npts = curve_rmsd(out["ppm"], out["experiment"], out["model"], win)
            rmsd = format(r, ".10g")
        except Exception as exc:
            curves_file = ""
            note = f"curves failed: {exc}"
            res.warnings.append(f"{names[k]}: {note}")

        res.manifest.append({
            "index": k + 1, "stem": names[k], "name": labels[k],
            "recipe_file": recipe_file, "curves_file": curves_file,
            "source_path": sp, "data_file": facts["data_file"],
            "source_kind": facts["source_kind"] or rec.source_kind or "",
            "source_sha256": facts["source_sha256"],
            "expno": facts["expno"], "procno": facts["procno"],
            "nucleus": rec.nucleus, "larmor_MHz": rec.larmor_frequency_MHz,
            "spin_rate_Hz": rec.spin_rate_Hz,
            "bf1_MHz": facts["bf1_MHz"], "sf_MHz": facts["sf_MHz"],
            "sr_hz": facts["sr_hz"], "ns": facts["ns"], "d1_s": facts["d1_s"],
            "pulprog": facts["pulprog"], "date": facts["date"],
            "title": facts["title"], "processing": processing,
            "fit_window_hi_ppm": hi, "fit_window_lo_ppm": lo,
            "n_points_window": npts, "rmsd": rmsd, "n_sites": len(rec.sites),
            "excluded_sites": ",".join(excluded), "released": released,
            "error_method": em_manifest,
            "larmor_version": versions["larmor"],
            "mrsimulator_version": versions["mrsimulator"], "note": note})

    # the tables, over the copies so their source_path column is filled
    tres = copy.copy(result)
    tres.recipes = copies
    if write_tables:
        batchfit.write_shared_csv(tres, outdir / f"{kind}_table.csv")
        res.files.append(f"{kind}_table.csv")
    if em_file:
        batchfit.write_error_csv(tres, outdir / f"{kind}_fit_{em_file}.csv", em_file)
        res.files.append(f"{kind}_fit_{em_file}.csv")

    with open(outdir / "manifest.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(MANIFEST_COLUMNS), lineterminator="\n")
        w.writeheader()
        for row in res.manifest:
            w.writerow(row)
    res.files.append("manifest.csv")

    glossary = ["README.txt — this file",
                "manifest.csv — one row per spectrum: stem, source file and its SHA-256, "
                "EXPNO/procno, NS, D1, BF1/SF/SR, processing steps, fit window, "
                "recomputed RMSD, excluded sites, error method, versions"]
    if write_tables or not write_recipes:
        glossary.append(f"{kind}_table.csv — the long table: shared parameters once, "
                        "then per spectrum the amplitudes, released parameters and "
                        "population_pct (value, stderr, model, source_path)")
    if em_file:
        glossary.append(f"{kind}_fit_{em_file}.csv — the same rows with the {em_file} "
                        "errors (stderr, sigma_pct, ci68_lo/hi)")
    glossary.append("<stem>.recipe.json — the fitted LARMOR recipe of that spectrum "
                    "(source_path and source_sha256 filled); opens in LARMOR")
    glossary.append("<stem>_curves.csv — ppm, experiment (as fitted), [experiment_raw,] "
                    "model, residual, s<i>_<label> per component, on the experimental "
                    "ppm axis")
    (outdir / "README.txt").write_text(
        readme_text(tres, kind, glossary, em_manifest, res.warnings), encoding="utf-8")
    res.files.append("README.txt")
    return res
