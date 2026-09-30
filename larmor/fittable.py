"""One parameter table from several fits (Qt-free).

The desktop's Fit parameter table dialog, the Workspaces dock's right-click
and the Explorer's fit rows all end here: a list of ``FitEntry`` (a name, a
recipe dict and the fit's facts) becomes either

* the **wide** layout -- one row per (fit, site), one column per parameter
  in ``PARAM_COLUMNS`` order (only the parameters some row uses, unknown keys
  appended with a header derived from the model's own ``ParamDef``), the
  cells ``value ± err`` with the ``larmor.paramstatus`` glyph († fixed,
  ‡ at a bound, § linked), a held value printed without an error bar, an
  empty cell when the site's model has no such parameter; or
* the **long** layout -- one row per parameter with value, stderr, the
  status as a word (the CSV vocabulary), bounds and link expression.

``PARAM_COLUMNS`` -- the curated column order and short headers every
parameter table shares -- lives here so the desktop's lines table
(``larmor.desktop.table``) and this module cannot drift apart; the desktop
imports it back. Statuses are never re-derived: ``paramstatus.param_status``
is the one rule.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

from larmor import cellparse, paramstatus

__all__ = [
    "PARAM_COLUMNS", "LONG_HEADERS", "FitEntry", "auto_header", "header_text",
    "fit_name", "load_fit_file", "entry_from_recipe", "format_value",
    "build_wide", "build_long", "cell_text", "to_csv", "to_tsv",
]

#: column order and short headers, dmfit-style
#: curated column order + headers for every parameter the registry's models
#: use today, roughly dmfit's panel order: amplitude/position/width first,
#: then quadrupolar, then CSA, then couplings/sidebands. This is a DISPLAY
#: preference only -- every table appends an automatic column for any key
#: this list does not know (header from the model's own ParamDef through
#: ``auto_header``), so a model can never again have a fitted-but-invisible
#: parameter, as the Amorphous ΔCq FWHM was for weeks.
PARAM_COLUMNS = [
    ("amplitude", "Amplitude"),
    ("isotropic_chemical_shift_ppm", "Position\n(ppm)"),
    ("shift_ppm", "Shift\n(ppm)"),                 # external-spectrum offset
    ("shift_fwhm_ppm", "Width\n(ppm)"),
    ("gauss_fwhm_ppm", "Gauss w\n(ppm)"),          # true Voigt
    ("lorentz_fwhm_ppm", "Lorentz w\n(ppm)"),      # Voigt, two-site exchange
    ("split_ppm", "Δδ A−B\n(ppm)"),                # two-site exchange
    ("pop_a", "p(A)"),
    ("k_ex_hz", "k_ex\n(s⁻¹)"),
    ("gl", "xG/(1-x)L"),
    ("sigma_Cq_MHz", "σ(Cq)\n(MHz)"),
    ("czjzek_d", "d\n(Czjzek)"),                   # general-d Czjzek
    ("shift_slope_ppm_per_MHz", "dδ/dC_Q\n(ppm/MHz)"),   # correlated Czjzek
    ("Cq_MHz", "Cq\n(MHz)"),
    ("Cq_fwhm_MHz", "ΔCq FWHM\n(MHz)"),
    ("eta", "η"),
    ("eta_q", "η(Q)"),
    ("eta_fwhm", "Δη FWHM"),
    ("eps", "ε"),                                  # extended Czjzek
    ("line_fwhm_ppm", "lb\n(ppm)"),
    ("zeta_ppm", "ζ CSA\n(ppm)"),
    ("sigma_zeta_ppm", "σ(ζ)\n(ppm)"),
    ("eta_cs", "η(CSA)"),
    ("j_hz", "J\n(Hz)"),
    ("n_j", "n(J)"),
    ("ssb_ratio", "SSB ratio"),
    ("n_ssb", "n(SSB)"),
]

#: the long layout's columns: one row per parameter, the status as the CSV
#: word ('fixed' | 'linked' | 'at_min' | 'at_max' | ''), never a glyph
LONG_HEADERS = ("fit", "site", "model", "param", "value", "stderr", "status",
                "min", "max", "expr")

#: the wide layout's leading and trailing (non-parameter) columns
_LEAD = ("fit", "site", "model")
_TAIL = ("RMSD", "nucleus", "ν0 (MHz)", "νrot (Hz)", "source")

#: file suffixes a fit file may carry, longest first (``fit_name`` strips one)
_FIT_SUFFIXES = (".recipe.json", ".json", ".fxmla", ".fxml")


def auto_header(model: str, key: str) -> str:
    """Header for a parameter the curated list does not know, taken from the
    model's own ParamDef (short key + unit) so it is never blank."""
    try:
        from larmor import models

        for pd in models.get(model).params:
            if pd.name == key:
                return f"{pd.key}\n({pd.unit})" if pd.unit else pd.key
    except Exception:
        pass
    return key


def header_text(header: str) -> str:
    """A curated header on one line: the lines table stacks 'Position' over
    '(ppm)'; a spreadsheet column reads 'Position (ppm)'."""
    return " ".join(part.strip() for part in str(header).split("\n") if part.strip())


@dataclass
class FitEntry:
    """One fit in the table: a display name, where it came from (the fit
    file, or the data path of an open workspace), the recipe dict and the
    fit's facts read from it."""

    name: str
    source: str = ""
    recipe: dict = field(default_factory=dict)
    rmsd: float | None = None
    nucleus: str = ""
    larmor_MHz: float = 0.0
    spin_rate_Hz: float = 0.0

    @property
    def sites(self) -> list[dict]:
        return [s for s in (self.recipe.get("sites") or []) if isinstance(s, dict)]


def fit_name(path) -> str:
    """The file name without its fit suffix: 'pCABS2-4_27Al.recipe.json' ->
    'pCABS2-4_27Al', '1r.fxml' -> '1r'."""
    name = Path(str(path)).name
    low = name.lower()
    for suf in _FIT_SUFFIXES:
        if low.endswith(suf) and len(name) > len(suf):
            return name[:-len(suf)]
    return name


def _num_or_none(v) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def entry_from_recipe(name: str, recipe: dict, source: str | None = None
                      ) -> FitEntry:
    """A ``FitEntry`` for a recipe dict already in memory (an open workspace):
    the fit RMSD, nucleus, Larmor frequency and spin rate come from the
    recipe; ``source`` defaults to the recipe's ``source_path``."""
    recipe = recipe if isinstance(recipe, dict) else {}
    return FitEntry(
        name=str(name or "") or "fit",
        source=str(source if source is not None else recipe.get("source_path") or ""),
        recipe=recipe,
        rmsd=_num_or_none(recipe.get("fit_rmsd")),
        nucleus=str(recipe.get("nucleus") or ""),
        larmor_MHz=float(_num_or_none(recipe.get("larmor_frequency_MHz")) or 0.0),
        spin_rate_Hz=float(_num_or_none(recipe.get("spin_rate_Hz")) or 0.0))


def load_fit_file(path) -> FitEntry:
    """Read a saved fit without its data: a LARMOR ``.recipe.json`` (or any
    recipe ``.json``) through ``Recipe.load``, a dmfit ``.fxml`` / ``.fxmla``
    through ``larmor.io.fxmla`` -- the same two readers the workbench's
    'Apply recipe' and 'Compare with a saved fit' use. Raises on a file
    neither reader accepts; the entry's source is the fit file itself."""
    from larmor.recipe import Recipe

    p = Path(str(path))
    low = p.name.lower()
    if low.endswith((".fxml", ".fxmla")):
        from larmor.io import fxmla

        recipe, _warnings = fxmla.to_recipe(fxmla.read(p))
        d = recipe.to_dict()
    else:
        d = Recipe.load(p).to_dict()
    return entry_from_recipe(fit_name(p), d, source=str(p))


# --------------------------------------------------------------- formatting
def _fmt_num(v) -> str:
    f = _num_or_none(v)
    if f is None:
        return "" if v is None else str(v)
    return f"{f:.5g}"


def _site_name(site: dict, index: int) -> str:
    return str(site.get("label") or "").strip() or cellparse.index_to_letter(index)


def format_value(p: dict, status: paramstatus.ParamStatus | None = None) -> str:
    """``value ± err`` plus the status glyph -- the wide layout's cell. A HELD
    value (fixed, or at its model default) prints without its error bar, like
    every other LARMOR table (a saved recipe carries stderr 0.0 on a fixed
    parameter, which would read as a fitted ``1.00 ± 0.00``); a value that
    finished at a bound has no covariance error after the fit's retry, so ‡
    is its explanation."""
    if not isinstance(p, dict) or p.get("value") is None:
        return ""
    txt = _fmt_num(p.get("value"))
    err = _num_or_none(p.get("stderr"))
    held = status is not None and status.held
    if err is not None and err > 0 and not held:
        txt += f" ± {err:.2g}"
    if status is not None and status.marker:
        txt += " " + status.marker
    return txt


def cell_text(c) -> str:
    """One CSV / TSV cell: '' for None, %.8g for a float, str otherwise."""
    if c is None:
        return ""
    if isinstance(c, float):
        return f"{c:.8g}" if math.isfinite(c) else ""
    return str(c)


# ------------------------------------------------------------------ layouts
def _param_columns(entries: list[FitEntry]) -> list[tuple[str, str]]:
    """The parameter columns in PARAM_COLUMNS order, only those some site
    uses, then every key the curated list does not know (in order of first
    appearance) with its ParamDef-derived header."""
    used_keys: set[str] = set()
    for e in entries:
        for s in e.sites:
            used_keys.update((s.get("params") or {}).keys())
    cols = [(k, h) for k, h in PARAM_COLUMNS if k in used_keys]
    known = {k for k, _ in PARAM_COLUMNS}
    for e in entries:
        for s in e.sites:
            for k in (s.get("params") or {}):
                if k not in known:
                    known.add(k)
                    cols.append((k, auto_header(str(s.get("model") or ""), k)))
    return cols


def build_wide(entries: list[FitEntry]) -> tuple[list[str], list[list[str]]]:
    """(headers, rows): one row per (fit, site) -- fit | site | model | the
    parameter columns | RMSD | nucleus | ν0 (MHz) | νrot (Hz) | source. Rows
    keep the entries' order, sites their recipe order. Every cell is text."""
    cols = _param_columns(entries)
    headers = list(_LEAD) + [header_text(h) for _, h in cols] + list(_TAIL)
    rows: list[list[str]] = []
    for e in entries:
        statuses = paramstatus.recipe_statuses(e.recipe)
        tail = [_fmt_num(e.rmsd) if e.rmsd is not None else "",
                e.nucleus,
                f"{e.larmor_MHz:.3f}" if e.larmor_MHz else "",
                f"{e.spin_rate_Hz:.0f}" if e.spin_rate_Hz else "",
                e.source]
        for i, site in enumerate(e.sites):
            params = site.get("params") or {}
            row = [e.name, _site_name(site, i), str(site.get("model") or "")]
            for key, _h in cols:
                p = params.get(key)
                row.append(format_value(p, statuses.get((i, key)))
                           if isinstance(p, dict) else "")
            rows.append(row + tail)
    return headers, rows


def build_long(entries: list[FitEntry]) -> tuple[tuple[str, ...], list[list]]:
    """(LONG_HEADERS, rows): one row per parameter -- fit, site, model,
    param, value, stderr (None when the value was held), status word
    (``ParamStatus.csv_flag``: 'fixed' | 'linked' | 'at_min' | 'at_max' |
    ''), min, max, expr. Numbers stay numbers; ``to_csv`` / ``to_tsv``
    format them."""
    rows: list[list] = []
    for e in entries:
        statuses = paramstatus.recipe_statuses(e.recipe)
        for i, site in enumerate(e.sites):
            model = str(site.get("model") or "")
            for key, p in (site.get("params") or {}).items():
                if not isinstance(p, dict):
                    continue
                st = statuses.get((i, key)) or paramstatus.param_status(model, key, p)
                rows.append([e.name, _site_name(site, i), model, key,
                             _num_or_none(p.get("value")),
                             None if st.held else _num_or_none(p.get("stderr")),
                             st.csv_flag,
                             _num_or_none(p.get("min")), _num_or_none(p.get("max")),
                             str(p.get("expr") or "")])
    return LONG_HEADERS, rows


# ------------------------------------------------------------------ outputs
def to_csv(headers, rows, path) -> None:
    """Write headers + rows as UTF-8 CSV (the csv module's quoting, '\\n'
    line ends); None cells are empty, floats %.8g."""
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow([str(h) for h in headers])
        for r in rows:
            w.writerow([cell_text(c) for c in r])


def to_tsv(headers, rows) -> str:
    """The table as tab-separated text for the clipboard (one header line,
    one line per row; a spreadsheet pastes it into cells)."""
    lines = ["\t".join(str(h) for h in headers)]
    for r in rows:
        lines.append("\t".join(cell_text(c).replace("\t", " ") for c in r))
    return "\n".join(lines)
