"""Batch publication report from a set of already-made fits.

Point it at several saved fits (LARMOR .recipe.json, dmfit .fxmla, or .larproj)
— ideally the same nucleus / acquisition — and it re-loads each with its data,
re-fits for fresh errors (covariance or Monte-Carlo), quantifies the site
populations, and writes a **publication table** (CSV + LaTeX + Markdown), a
per-fit overlay figure, and a Markdown report to a chosen folder.

The heavy lifting reuses the rest of LARMOR: `loader.load_any` (any fit format
+ its data), `fit.fit` (covariance errors), `autofit.monte_carlo_errors`
(bootstrap errors), `quantify.quantify` (populations), and `engine.simulate`
(the overlay curves). This module is Qt-free and testable; the desktop dialog is
only a front-end.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from larmor import autofit
from larmor import fit as fitmod
from larmor import quantify as quantmod
from larmor.recipe import Recipe


# --------------------------------------------------------------------------
@dataclass
class FitEntry:
    path: str
    sample: str
    nucleus: str
    larmor_MHz: float
    spin_rate_Hz: float
    recipe: dict
    ppm: np.ndarray
    amp: np.ndarray
    window: tuple | None
    warnings: list[str] = field(default_factory=list)

    @property
    def models(self) -> tuple[str, ...]:
        return tuple(s.get("model", "?") for s in self.recipe.get("sites", []))


def _fit_entry(path: str, rec: dict, ppm, amp, warns=()) -> FitEntry:
    sample = rec.get("sample") or Path(path).stem
    return FitEntry(
        path=path, sample=sample, nucleus=rec.get("nucleus", ""),
        larmor_MHz=float(rec.get("larmor_frequency_MHz", 0.0) or 0.0),
        spin_rate_Hz=float(rec.get("spin_rate_Hz", 0.0) or 0.0),
        recipe=rec, ppm=np.asarray(ppm, float), amp=np.asarray(amp, float),
        window=tuple(rec["fit_window_ppm"]) if rec.get("fit_window_ppm")
        else None, warnings=list(warns or []))


def _larproj_workspaces(path: str) -> list[dict] | None:
    """The saved workspaces of a LARMOR project bundle (app.py's
    save_project/open_project -- ppm/amp/recipe embedded directly, no
    separate data load needed), or None if `path` isn't one -- callers fall
    back to the normal single-fit load_any path for anything else."""
    import json

    try:
        d = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(d, dict) or "workspaces" not in d:
        return None
    return d.get("workspaces") or []


def load_entries(paths) -> tuple[list[FitEntry], list[str]]:
    """Load each fit file with its source data. Returns (entries, warnings).

    A LARMOR project bundle (`.larproj.json`, several spectra saved as one
    file) expands into one entry per workspace that carries a fit -- it
    isn't itself a single Recipe, so load_any's plain `.json` branch can't
    read it (previously: silently skipped with "could not load", the project
    file's own source_path-less shape never matching a Recipe well enough to
    find real data -- this is what the docs already claimed worked)."""
    from larmor.loader import load_any

    entries, warnings = [], []
    for p in paths:
        p = str(p)
        workspaces = _larproj_workspaces(p)
        if workspaces is not None:
            if not workspaces:
                warnings.append(f"{Path(p).name}: no spectra in this project — skipped")
                continue
            for w in workspaces:
                rec = w.get("recipe")
                ppm = w.get("exp_ppm")
                if not rec or not rec.get("sites") or not ppm:
                    continue
                entries.append(_fit_entry(p, rec, ppm, w.get("exp_amp", [])))
            continue
        try:
            ppm, amp, rec, meta, warns = load_any(p)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{Path(p).name}: could not load ({exc}) — skipped")
            continue
        if not rec.get("sites"):
            warnings.append(f"{Path(p).name}: no fitted sites — skipped")
            continue
        entries.append(_fit_entry(p, rec, ppm, amp, warns))
    return entries, warnings


def homogeneity(entries: list[FitEntry]) -> list[str]:
    """Warn about the mixes the user asked to be told about."""
    notes = []
    nuclei = {e.nucleus for e in entries}
    if len(nuclei) > 1:
        notes.append(f"⚠ mixed nuclei: {', '.join(sorted(nuclei))} — the table "
                     "columns and populations are only comparable within a nucleus")
    lar = [e.larmor_MHz for e in entries if e.larmor_MHz]
    if lar and (max(lar) - min(lar)) / (np.mean(lar) or 1) > 0.05:
        notes.append(f"⚠ Larmor frequencies span {min(lar):.1f}–{max(lar):.1f} MHz "
                     "(different fields) — δ₂-dependent widths are not directly "
                     "comparable")
    modelsets = {e.models for e in entries}
    if len(modelsets) > 1:
        notes.append("⚠ different model sets across fits — some table cells will "
                     "be blank where a column does not apply")
    return notes


# --------------------------------------------------------------------------
# the publication columns (model-aware; each cell carries an error)

def _p(site, name):
    p = site.get("params", {}).get(name)
    return p["value"] if p else None


def _site_columns(site: dict, errors: dict) -> list[tuple[str, float, float | None]]:
    """Ordered (header, value, error) for one site. `errors` maps param name ->
    error (covariance stderr or MC σ). Derived quantities propagate the error."""
    P = site.get("params", {})
    out: list[tuple[str, float, float | None]] = []

    pos = _p(site, "isotropic_chemical_shift_ppm")
    if pos is not None:
        out.append(("δiso (ppm)", pos, errors.get("isotropic_chemical_shift_ppm")))

    if "sigma_Cq_MHz" in P:                       # Czjzek family
        s = P["sigma_Cq_MHz"]["value"]
        se = errors.get("sigma_Cq_MHz")
        from larmor.czjzek_dist import rms_pq
        out.append(("σ (MHz)", s, se))
        out.append(("C_Q=2σ (MHz)", 2.0 * s, (2.0 * se) if se else None))
        out.append(("√⟨P_Q²⟩ (MHz)", rms_pq(s),
                    (5.0 ** 0.5 * se) if se else None))
    if "Cq_MHz" in P:                             # discrete / amorphous
        out.append(("C_Q (MHz)", P["Cq_MHz"]["value"], errors.get("Cq_MHz")))
    for ek in ("eta", "etaQ", "eta_q"):
        if ek in P:
            out.append(("η", P[ek]["value"], errors.get(ek)))
            break
    if "eps" in P:
        out.append(("ε", P["eps"]["value"], errors.get("eps")))
    if "Cq_fwhm_MHz" in P:
        out.append(("ΔC_Q (MHz)", P["Cq_fwhm_MHz"]["value"],
                    errors.get("Cq_fwhm_MHz")))
    for wk, wl in (("shift_fwhm_ppm", "dCS/FWHM (ppm)"),
                   ("gauss_fwhm_ppm", "Gauss FWHM (ppm)"),
                   ("lorentz_fwhm_ppm", "Lorentz FWHM (ppm)")):
        if wk in P:
            out.append((wl, P[wk]["value"], errors.get(wk)))
    return out


# --------------------------------------------------------------------------
@dataclass
class BatchTable:
    headers: list[str]                 # column order (union across all rows)
    rows: list[dict]                   # each: sample, site, label, model, cells{header:(v,e)}, pop
    notes: list[str] = field(default_factory=list)
    error_method: str = "covariance"


def _errors_for(entry: FitEntry, method: str, n_mc: int, seed: int,
                progress=None) -> tuple[Recipe, dict, dict]:
    """Return (fitted recipe, {(site,param):err}, quant table). Always does one
    covariance fit (best fit + populations); for 'montecarlo' the parameter
    errors are replaced by the bootstrap σ."""
    rec = Recipe.from_dict(entry.recipe)
    fitmod.fit(rec, entry.ppm, entry.amp, window_ppm=entry.window)
    errors = {(i, pn): p.stderr
              for i, s in enumerate(rec.sites) for pn, p in s.params.items()}
    if method == "montecarlo":
        mc = autofit.monte_carlo_errors(
            rec, entry.ppm, entry.amp, window_ppm=entry.window,
            n_trials=n_mc, seed=seed, progress=progress, parallel=True)
        for mp in mc.params:
            errors[(mp.site, mp.param)] = mp.std
    quant = quantmod.quantify(rec, window_ppm=entry.window)
    return rec, errors, quant


def build_table(entries: list[FitEntry], error_method: str = "covariance",
                n_mc: int = 200, seed: int = 0, progress=None,
                should_stop=None) -> BatchTable:
    """Refit every entry, compute errors + populations, and assemble the table."""
    rows: list[dict] = []
    headers: list[str] = []
    for k, e in enumerate(entries):
        if should_stop is not None and should_stop():
            break
        rec, errors, quant = _errors_for(e, error_method, n_mc, seed)
        qrows = {r["site"]: r for r in quant["rows"]}
        for i, s in enumerate(rec.sites):
            sdict = {"model": s.model, "label": s.label,
                     "params": {pn: {"value": p.value, "stderr": p.stderr}
                                for pn, p in s.params.items()}}
            perr = {pn: errors.get((i, pn)) for pn in sdict["params"]}
            cells = {h: (v, err) for h, v, err in _site_columns(sdict, perr)}
            for h in cells:
                if h not in headers:
                    headers.append(h)
            qr = qrows.get(f"s{i}", {})
            rows.append({
                "sample": e.sample, "site": f"s{i}",
                "label": s.label or s.model, "model": s.model,
                "cells": cells,
                "pop": (qr.get("fraction_pct"), qr.get("fraction_err_pct")),
            })
        if progress:
            progress(k + 1, len(entries))
    notes = homogeneity(entries)
    return BatchTable(headers=headers, rows=rows, notes=notes,
                      error_method=error_method)


# --------------------------------------------------------------------------
# formatting helpers

def _fmt(v, e) -> str:
    if v is None:
        return ""
    if e is not None and np.isfinite(e) and e > 0:
        return f"{v:.4g} ± {e:.2g}"
    return f"{v:.4g}"


def _slug(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", s).strip("_") or "fit"


def _all_columns(t: BatchTable) -> list[str]:
    return t.headers + ["pop (%)"]


def _cell_text(row: dict, h: str) -> str:
    if h == "pop (%)":
        return _fmt(*row["pop"])
    v, e = row["cells"].get(h, (None, None))
    return _fmt(v, e)


def write_csv(t: BatchTable, path: Path) -> None:
    cols = _all_columns(t)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sample", "site", "label", "model", *cols,
                    *[c + " err" for c in cols]])
        for r in t.rows:
            vals, errs = [], []
            for h in cols:
                if h == "pop (%)":
                    v, e = r["pop"]
                else:
                    v, e = r["cells"].get(h, (None, None))
                vals.append("" if v is None else f"{v:.6g}")
                errs.append("" if e is None else f"{e:.6g}")
            w.writerow([r["sample"], r["site"], r["label"], r["model"],
                        *vals, *errs])


def write_latex(t: BatchTable, path: Path) -> None:
    cols = _all_columns(t)
    lines = [r"% LARMOR batch publication table",
             r"\begin{table}[htbp]\centering",
             r"\caption{Fitted parameters (%s errors).}" % t.error_method,
             r"\begin{tabular}{ll" + "r" * len(cols) + "}",
             r"\hline"]
    header = " & ".join(["Sample", "Site"] + [c.replace("_", r"\_")
                                              .replace("²", r"$^2$")
                                              .replace("√", r"$\sqrt{~}$")
                                              for c in cols]) + r" \\"
    lines += [header, r"\hline"]
    for r in t.rows:
        cells = [_cell_text(r, h).replace("±", r"$\pm$") for h in cols]
        lines.append(" & ".join([r["sample"].replace("_", r"\_"),
                                 r["label"].replace("_", r"\_")] + cells) + r" \\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def _markdown_table(t: BatchTable) -> str:
    cols = _all_columns(t)
    head = "| Sample | Site | Model | " + " | ".join(cols) + " |"
    sep = "|" + "---|" * (3 + len(cols))
    out = [head, sep]
    for r in t.rows:
        cells = [_cell_text(r, h) for h in cols]
        out.append(f"| {r['sample']} | {r['label']} | {r['model']} | "
                   + " | ".join(cells) + " |")
    return "\n".join(out)


# --------------------------------------------------------------------------
# per-fit overlay figure (data + model + components + residual)

def render_overlay(entry: FitEntry, recipe: Recipe, path: Path) -> bool:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from larmor import engine

    try:
        x, total, per = engine.simulate(recipe, exp_ppm=entry.ppm)
    except Exception:
        return False
    total = np.asarray(total, float)
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.6, 4.2), sharex=True,
                                  gridspec_kw={"height_ratios": [4, 1],
                                               "hspace": 0.05})
    ax.plot(entry.ppm, entry.amp, color="#222", lw=1.0, label="experiment")
    ax.plot(entry.ppm, total, color="#c0392b", lw=1.3, ls=(0, (4, 2)),
            label="model")
    from larmor.figures import site_color
    for i, ys in enumerate(per):
        ax.fill_between(entry.ppm, np.asarray(ys, float),
                        color=site_color(i), alpha=0.18)
    ax.set_xlim(entry.ppm.max(), entry.ppm.min())
    ax.set_yticks([])
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title(f"{entry.sample}  ·  {entry.nucleus}", fontsize=9)
    axr.plot(entry.ppm, entry.amp - total, color="#888", lw=0.7)
    axr.axhline(0, color="#ccc", lw=0.6)
    axr.set_xlim(entry.ppm.max(), entry.ppm.min())
    axr.set_yticks([])
    axr.set_xlabel(f"{entry.nucleus} shift (ppm)")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
@dataclass
class BatchResult:
    outdir: str
    files: list[str]
    n_fits: int
    n_sites: int
    error_method: str
    warnings: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        return (f"{self.n_fits} fits · {self.n_sites} sites · {self.error_method} "
                f"errors → {len(self.files)} files in {self.outdir}")


def run_batch(paths, outdir, *, error_method: str = "covariance", n_mc: int = 200,
              seed: int = 0, make_plots: bool = True,
              formats=("csv", "latex", "markdown"), progress=None,
              should_stop=None) -> BatchResult:
    """Full pipeline. `progress(stage, k, n)` is called for 'load'/'fit'/'plot'."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    files: list[str] = []

    entries, warnings = load_entries(paths)
    if not entries:
        raise ValueError("no usable fits: " + "; ".join(warnings))

    table = build_table(
        entries, error_method=error_method, n_mc=n_mc, seed=seed,
        progress=(lambda k, n: progress("fit", k, n)) if progress else None,
        should_stop=should_stop)

    # per-fit refit recipes (rebuild once more for the plots, cheap vs the fit)
    figdir = outdir / "figures"
    plot_links: dict[str, str] = {}
    if make_plots:
        figdir.mkdir(exist_ok=True)
        for k, e in enumerate(entries):
            rec = Recipe.from_dict(e.recipe)
            try:
                fitmod.fit(rec, e.ppm, e.amp, window_ppm=e.window)
            except Exception:
                pass
            fp = figdir / f"{_slug(e.sample)}.png"
            if render_overlay(e, rec, fp):
                plot_links[e.sample] = f"figures/{fp.name}"
                files.append(str(fp))
            if progress:
                progress("plot", k + 1, len(entries))

    if "csv" in formats:
        p = outdir / "table.csv"; write_csv(table, p); files.append(str(p))
    if "latex" in formats:
        p = outdir / "table.tex"; write_latex(table, p); files.append(str(p))
    if "markdown" in formats:
        p = outdir / "report.md"
        _write_report(p, entries, table, plot_links, warnings)
        files.append(str(p))

    n_sites = len(table.rows)
    return BatchResult(outdir=str(outdir), files=files, n_fits=len(entries),
                       n_sites=n_sites, error_method=error_method,
                       warnings=warnings + table.notes)


def _write_report(path: Path, entries, table: BatchTable, plot_links, warnings):
    from datetime import date

    method = {"covariance": "covariance matrix (lmfit)",
              "montecarlo": "Monte-Carlo (parametric bootstrap)"}.get(
        table.error_method, table.error_method)
    lines = [
        f"# LARMOR batch fit report",
        "",
        f"*{len(entries)} fits · {len(table.rows)} sites · generated "
        f"{date.today().isoformat()}*",
        "",
        "## Fits", "",
    ]
    nuclei = sorted({e.nucleus for e in entries})
    fields_ = sorted({round(e.larmor_MHz, 1) for e in entries if e.larmor_MHz})
    lines.append(f"- **Nucleus:** {', '.join(nuclei) or '—'}")
    lines.append(f"- **Larmor:** {', '.join(f'{f:g}' for f in fields_)} MHz")
    lines.append(f"- **Errors:** {method}")
    lines.append("")
    for n in table.notes:
        lines.append(f"> {n}")
    if table.notes:
        lines.append("")
    lines += ["## Table", "", _markdown_table(table), "",
              "Errors are ± one standard error. Population % is the integrated "
              "area over the fit window (first-order amplitude error). "
              "`table.csv` / `table.tex` hold the same data.", ""]
    if plot_links:
        lines += ["## Fits (overlays)", ""]
        for e in entries:
            if e.sample in plot_links:
                lines.append(f"### {e.sample}")
                lines.append("")
                lines.append(f"![{e.sample}]({plot_links[e.sample]})")
                lines.append("")
    if warnings:
        lines += ["## Load notes", ""] + [f"- {w}" for w in warnings] + [""]
    path.write_text("\n".join(lines), encoding="utf-8")
