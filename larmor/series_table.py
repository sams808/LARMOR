"""Series table: the identity, order and metadata of a batch / sequential series.

A batch fit (or a sequential sweep) is a list of spectra. What each spectrum
*is* -- its name, which glass it belongs to, where it sits along the series,
its composition -- lives here, Qt-free, so the batch dialog, the sequential
dialog, the Series evolution plot and the exports all read one object.

Identity comes from :mod:`larmor.io.scan` (``sample_label`` /
``sample_name`` / ``disambiguate``): the sample folder without its date, rotor
and operator tokens, made distinct across folders by the date token
(``P5-Bi8-12 (04272026)``). This module never re-derives a name; it keeps
what the loader produced and adds:

* **unique names** -- a name typed by hand that collides with another row
  gets `` (2)``, `` (3)``… (:func:`unique_names`), because names become CSV
  scopes, recipe stems and Plotting-studio panel keys;
* a **group** per row, defaulting to the base sample name (the label before
  disambiguation), which is what replicate averaging collapses;
* an explicit **order** (``move`` / ``apply_order`` / ``sort_perm`` -- natural
  sort, so ``Base0Ca, Base1Ca, …, Base10Ca``);
* typed or CSV-joined **numeric columns** (:class:`SeriesColumn`, with an
  optional ± partner paired once by suffix -- :data:`ERR_SUFFIXES`),
  behind a proposed-then-confirmed row mapping (:func:`propose_mapping`);
* **replicate statistics** (:func:`replicate_stats`: mean, sample std with
  ddof = 1, n) and an **OLS** line (:func:`ols`: polyfit with covariance when
  n >= 3, Pearson r);
* the wide ``<stem>_series.csv`` companion of the long batch CSV
  (:func:`write_series_csv`) and a **DUST**-shaped composition CSV
  (:func:`dust_rows`: ``Sample``, canonical oxide columns as DUST's importer
  recognises them -- :data:`DUST_OXIDES` mirrors ``core/oxides.py`` of the
  DUST repository -- plus ``N4_measured`` / ``N4_measured_err``, which DUST
  maps to *Ignore* on import and which are meant for a join against its
  Results CSV);
* JSON persistence (``.series.json``).

Named ``series_table`` because :mod:`larmor.series` is the relaxation-series
engine.
"""
from __future__ import annotations

import csv
import json
import re
import unicodedata
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from larmor.io import scan

__all__ = [
    "DUST_OXIDES", "ERR_SUFFIXES", "SeriesColumn", "SeriesRow", "SeriesTable",
    "Match", "natural_key", "norm_key", "unique_names", "as_float", "read_wide_csv",
    "is_long_batch_csv", "read_batch_csv_wide", "numeric_headers",
    "pair_error_columns", "propose_mapping", "apply_join", "replicate_stats",
    "ols", "ols_label", "write_series_csv", "canonical_oxide", "dust_rows",
    "species_bar_spec", "ols_trace", "correlation_traces",
]

#: DUST's canonical oxide list (``core/oxides.py`` of github.com/sams808/DUST),
#: the plain names its CSV importer maps a header to after normalisation
DUST_OXIDES: tuple[str, ...] = (
    "SiO2", "B2O3", "Al2O3", "P2O5", "GeO2", "TeO2", "As2O5", "Sb2O3",
    "Fe2O3", "Cr2O3", "TiO2", "ZrO2", "HfO2", "Nb2O5", "Ta2O5",
    "WO3", "MoO3", "V2O5", "SO3", "SnO2",
    "Na2O", "Li2O", "K2O", "Cs2O", "Rb2O",
    "CaO", "MgO", "SrO", "BaO", "BeO",
    "ZnO", "PbO", "CdO", "MnO", "NiO", "CoO", "CuO",
    "La2O3", "Y2O3", "Bi2O3", "CeO2", "Nd2O3", "Sm2O3", "Gd2O3",
    "In2O3", "Ga2O3", "Sc2O3",
    "ThO2", "UO3", "NpO2", "PuO2",
)

#: how a CSV names the uncertainty of a value column: ``P2O5_mol_sd``,
#: ``Vm_u``, ``x_err``, ``"A population_pct ±"`` -- paired once at join time
ERR_SUFFIXES: tuple[str, ...] = ("_sd", "_err", "_u", "_std", "_stderr",
                                 "_sigma", "_unc", " ±", " +/-", "_±", "±")

#: candidate sample-name headers of a wide CSV, in order of preference
_SAMPLE_HEADERS = ("sample", "name", "scope", "label", "id", "glass")

_NUM_RE = re.compile(r"(\d+)")
_NORM_RE = re.compile(r"[^a-z0-9]+")
_SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


# ---------------------------------------------------------------- names
def natural_key(s: str) -> tuple:
    """Sort key that orders embedded numbers numerically (``9Ca`` < ``10Ca``)."""
    return tuple((0, int(tok)) if tok.isdigit() else (1, tok.casefold())
                 for tok in _NUM_RE.split(str(s)) if tok != "")


def norm_key(s: str) -> str:
    """A loose identity for matching names across files: casefolded,
    ``[a-z0-9]`` only (``P5-Bi8-12`` == ``P5Bi8-12`` == ``p5 bi8 12``)."""
    return _NORM_RE.sub("", str(s).casefold())


def unique_names(names: list[str]) -> list[str]:
    """Make equal names distinct with `` (2)``, `` (3)``… (first occurrence
    unchanged). Blank names become ``spectrum <k>``."""
    out: list[str] = []
    seen: set[str] = set()
    for k, raw in enumerate(names):
        name = (str(raw) if raw is not None else "").strip() or f"spectrum {k + 1}"
        cand, i = name, 1
        while cand in seen:
            i += 1
            cand = f"{name} ({i})"
        seen.add(cand)
        out.append(cand)
    return out


def _folder_of(path: str) -> str:
    """The sample folder of a Bruker-shaped path ('' for a file with an
    extension -- CSV, fxmla, recipe -- which has no sample folder)."""
    p = Path(str(path))
    if p.suffix:
        return ""
    return scan._folder_and_expno(p)[0]


# ---------------------------------------------------------------- data model
@dataclass
class SeriesColumn:
    """One numeric metadata column: ``key`` indexes ``SeriesRow.values``,
    ``label`` is what the tables and axes show (the analysed / nominal tag is
    a prefix of the label), ``err_key`` names the ± partner in ``values``."""
    key: str
    label: str
    err_key: str | None = None
    tag: str = ""          # 'analysed' | 'nominal' | ''
    source: str = ""       # the file the column was joined from


@dataclass
class SeriesRow:
    source_path: str
    display_name: str
    group: str
    folder: str = ""
    title: str = ""
    values: dict = field(default_factory=dict)


@dataclass
class SeriesTable:
    rows: list
    columns: list = field(default_factory=list)

    # ---- construction
    @classmethod
    def from_paths(cls, paths, names=None, groups=None, folders=None,
                   titles=None) -> "SeriesTable":
        """Rows in the given order. ``names`` default to the loader's rule
        (``scan.sample_label`` made distinct by ``scan.disambiguate``);
        ``groups`` default to the base label (before disambiguation);
        ``folders`` to the sample folder read off the path."""
        paths = [str(p) for p in paths]
        n = len(paths)
        base = None
        if names is None:
            base = [scan.sample_label(p, {}) for p in paths]
            names = scan.disambiguate(base, paths)
        names = unique_names(list(names))
        if groups is None:
            groups = base if base is not None else list(names)
        if folders is None:
            folders = [_folder_of(p) for p in paths]
        if titles is None:
            titles = [""] * n
        rows = [SeriesRow(source_path=paths[k], display_name=names[k],
                          group=(str(groups[k]) if groups[k] else names[k]),
                          folder=str(folders[k] or ""), title=str(titles[k] or ""))
                for k in range(n)]
        return cls(rows=rows, columns=[])

    @classmethod
    def from_spectra(cls, data: list) -> "SeriesTable":
        """From the batch / sequential dialogs' per-spectrum dicts (keys
        ``path``, ``sample``, and optionally ``group``, ``folder``, ``title``)."""
        return cls.from_paths(
            [d["path"] for d in data], names=[d.get("sample", "") for d in data],
            groups=[d.get("group") or d.get("sample", "") for d in data],
            folders=[d.get("folder", "") for d in data],
            titles=[d.get("title", "") for d in data])

    # ---- read
    def labels(self) -> list[str]:
        return [r.display_name for r in self.rows]

    def paths(self) -> list[str]:
        return [r.source_path for r in self.rows]

    def groups(self) -> dict:
        """group -> row indices, in first-seen order."""
        out: dict = {}
        for k, r in enumerate(self.rows):
            out.setdefault(r.group or r.display_name, []).append(k)
        return out

    def has_replicates(self) -> bool:
        return any(len(idx) > 1 for idx in self.groups().values())

    def numeric_columns(self) -> list:
        return list(self.columns)

    def column(self, key: str):
        return next((c for c in self.columns if c.key == key), None)

    def x_values(self, key: str) -> tuple:
        """(x, xerr) aligned to the rows, NaN where missing."""
        col = self.column(key)
        x = np.array([_as_float(r.values.get(key)) for r in self.rows], float)
        if col is not None and col.err_key:
            xerr = np.array([_as_float(r.values.get(col.err_key)) for r in self.rows],
                            float)
        else:
            xerr = np.full(len(self.rows), np.nan)
        return x, xerr

    # ---- edit
    def rename(self, i: int, name: str) -> str:
        """Rename row ``i``; the applied name is made unique against the other
        rows. A row whose group equalled its old name moves its group too."""
        others = [r.display_name for k, r in enumerate(self.rows) if k != i]
        row = self.rows[i]
        new = unique_names(others + [name])[-1]
        if row.group == row.display_name or not row.group:
            row.group = new
        row.display_name = new
        return new

    def set_group(self, i: int, group: str) -> None:
        self.rows[i].group = (str(group) or "").strip() or self.rows[i].display_name

    def move(self, src: int, dst: int) -> None:
        row = self.rows.pop(src)
        self.rows.insert(dst, row)

    def apply_order(self, perm: list) -> None:
        """Reorder in place: ``perm[new_pos] = old index``."""
        self.rows = [self.rows[k] for k in perm]

    def permuted(self, perm: list) -> "SeriesTable":
        import copy
        return SeriesTable(rows=[copy.deepcopy(self.rows[k]) for k in perm],
                           columns=copy.deepcopy(self.columns))

    def sort_perm(self, key=None, descending: bool = False) -> list:
        """The permutation sorting the rows: ``key`` None = natural name
        order, else a numeric column (stable; rows without a value last)."""
        n = len(self.rows)
        if key is None:
            perm = sorted(range(n), key=lambda k: natural_key(self.rows[k].display_name))
            return perm[::-1] if descending else perm
        x, _ = self.x_values(key)
        have = [k for k in range(n) if np.isfinite(x[k])]
        missing = [k for k in range(n) if not np.isfinite(x[k])]
        have.sort(key=lambda k: x[k], reverse=descending)
        return have + missing

    def add_column(self, col: SeriesColumn, values=None) -> None:
        """Add (or replace, by key) a column; ``values`` maps row index ->
        value (missing rows None)."""
        self.remove_column(col.key)
        self.columns.append(col)
        for k, r in enumerate(self.rows):
            r.values[col.key] = _as_float((values or {}).get(k)) if values else None
            if col.err_key and col.err_key not in r.values:
                r.values[col.err_key] = None

    def remove_column(self, key: str) -> None:
        col = self.column(key)
        if col is None:
            return
        self.columns = [c for c in self.columns if c.key != key]
        for r in self.rows:
            r.values.pop(key, None)
            if col.err_key:
                r.values.pop(col.err_key, None)

    # ---- pairing
    def aligned_to(self, paths, fresh: "SeriesTable | None" = None) -> "SeriesTable":
        """A table over ``paths`` (in that order): rows paired by source_path
        (pop-first, the batchfit.align_result rule); a path with no row takes
        the matching row of ``fresh`` (the dialog's default table) or a bare
        row. Names are re-uniquified."""
        import copy
        slots: dict = {}
        for k, r in enumerate(self.rows):
            slots.setdefault(r.source_path, []).append(k)
        fresh_by_path: dict = {}
        if fresh is not None:
            for r in fresh.rows:
                fresh_by_path.setdefault(r.source_path, []).append(r)
        rows = []
        for p in paths:
            p = str(p)
            idx = slots.get(p)
            if idx:
                rows.append(copy.deepcopy(self.rows[idx.pop(0)]))
                continue
            fr = fresh_by_path.get(p)
            if fr:
                r = copy.deepcopy(fr.pop(0))
            else:
                r = SeriesRow(source_path=p, display_name=Path(p).stem,
                              group=Path(p).stem, folder=_folder_of(p))
            for c in self.columns:
                r.values.setdefault(c.key, None)
                if c.err_key:
                    r.values.setdefault(c.err_key, None)
            rows.append(r)
        names = unique_names([r.display_name for r in rows])
        for r, nm in zip(rows, names):
            if r.group == r.display_name:
                r.group = nm
            r.display_name = nm
        return SeriesTable(rows=rows, columns=copy.deepcopy(self.columns))

    def merge_from(self, other: "SeriesTable") -> list:
        """Take names, groups and columns from a loaded table: rows pair by
        source_path, else by exact display name, else by normalised name.
        Returns, per row, the index of the ``other`` row it took (None when
        unmatched) -- the saved order, for a caller that wants to follow it."""
        by_path = {r.source_path: k for k, r in enumerate(other.rows) if r.source_path}
        by_name = {r.display_name: k for k, r in enumerate(other.rows)}
        by_norm: dict = {}
        for k, r in enumerate(other.rows):
            by_norm.setdefault(norm_key(r.display_name), k)
        pairs: list = []
        for r in self.rows:
            j = by_path.get(r.source_path)
            if j is None:
                j = by_name.get(r.display_name)
            if j is None:
                j = by_norm.get(norm_key(r.display_name))
            pairs.append(j)
        for c in other.columns:
            self.remove_column(c.key)
        self.columns.extend(SeriesColumn(**asdict(c)) for c in other.columns)
        for r, j in zip(self.rows, pairs):
            for c in other.columns:
                r.values[c.key] = None
                if c.err_key:
                    r.values[c.err_key] = None
            if j is None:
                continue
            o = other.rows[j]
            r.group = o.group or o.display_name
            r.display_name = o.display_name
            for c in other.columns:
                r.values[c.key] = _as_float(o.values.get(c.key))
                if c.err_key:
                    r.values[c.err_key] = _as_float(o.values.get(c.err_key))
        names = unique_names([r.display_name for r in self.rows])
        for r, nm in zip(self.rows, names):
            if r.group == r.display_name:
                r.group = nm
            r.display_name = nm
        return pairs

    # ---- persistence
    def to_dict(self) -> dict:
        return {"version": 1,
                "rows": [asdict(r) for r in self.rows],
                "columns": [asdict(c) for c in self.columns]}

    @classmethod
    def from_dict(cls, d: dict) -> "SeriesTable":
        rows = []
        for r in d.get("rows") or []:
            r = dict(r)
            r.pop("raw_name", None)
            rows.append(SeriesRow(
                source_path=str(r.get("source_path", "")),
                display_name=str(r.get("display_name", "")),
                group=str(r.get("group") or r.get("display_name", "")),
                folder=str(r.get("folder", "") or ""),
                title=str(r.get("title", "") or ""),
                values={k: _as_float(v) for k, v in (r.get("values") or {}).items()}))
        cols = [SeriesColumn(key=str(c.get("key")), label=str(c.get("label", c.get("key"))),
                             err_key=c.get("err_key") or None, tag=str(c.get("tag", "") or ""),
                             source=str(c.get("source", "") or ""))
                for c in (d.get("columns") or [])]
        return cls(rows=rows, columns=cols)

    def save(self, path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path) -> "SeriesTable":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def _as_float(v):
    """float(v) or None (blank, None, non-numeric, non-finite)."""
    if v is None:
        return None
    if isinstance(v, str):
        s = v.strip().replace(",", ".") if v.count(",") == 1 and "." not in v else v.strip()
        if not s:
            return None
        try:
            f = float(s)
        except ValueError:
            return None
    else:
        try:
            f = float(v)
        except (TypeError, ValueError):
            return None
    return f if np.isfinite(f) else None


#: the public spelling, for the dialogs' cell parsing
as_float = _as_float


# ---------------------------------------------------------------- CSV join
def _sniff(text: str) -> str:
    head = "\n".join(text.splitlines()[:20])
    try:
        return csv.Sniffer().sniff(head, delimiters=",;\t").delimiter
    except csv.Error:
        first = text.splitlines()[0] if text else ""
        counts = {d: first.count(d) for d in (",", ";", "\t")}
        return max(counts, key=counts.get) if any(counts.values()) else ","


def read_wide_csv(path) -> tuple:
    """(headers, rows, sample_header) of a wide CSV -- one row per sample.
    utf-8-sig (Excel's BOM), delimiter sniffed among ``, ; \\t`` with a
    comma fallback; the sample header is the first of sample / name / scope /
    label / id (case-insensitive) else the first non-numeric column."""
    text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
    delim = _sniff(text)
    reader = csv.reader(text.splitlines(), delimiter=delim)
    headers: list = []
    rows: list = []
    for rec in reader:
        if not any(c.strip() for c in rec):
            continue
        if not headers:
            headers = [h.strip() for h in rec]
            continue
        rows.append({h: (rec[i].strip() if i < len(rec) else "")
                     for i, h in enumerate(headers)})
    return headers, rows, _sample_header(headers, rows)


def _sample_header(headers: list, rows: list) -> str:
    low = {h.strip().casefold(): h for h in headers}
    for cand in _SAMPLE_HEADERS:
        if cand in low:
            return low[cand]
    nums = set(numeric_headers(headers, rows))
    for h in headers:
        if h not in nums:
            return h
    return headers[0] if headers else ""


def is_long_batch_csv(headers) -> bool:
    """A LARMOR long table (batch_table*.csv / seq_table.csv / the error
    CSV): scope, site, param, value columns."""
    return {"scope", "site", "param", "value"} <= {str(h).strip() for h in headers}


def read_batch_csv_wide(path) -> tuple:
    """A LARMOR long CSV pivoted to one row per scope (sample) -- columns
    ``<label or site> <param>`` and ``<…> ±`` from stderr -- so another
    nucleus' batch results join a series like a composition table would.
    The ``shared`` scope is dropped. Returns (headers, rows, 'scope')."""
    from larmor import series_grid

    by_scope = series_grid.csv_rows_by_scope(path)
    headers: list = ["scope"]
    rows: list = []
    for scope, recs in by_scope.items():
        if scope == "shared" or not scope:
            continue
        row = {"scope": scope}
        for r in recs:
            site = (r.get("label") or "").strip() or (r.get("site") or "").strip()
            param = (r.get("param") or "").strip()
            key = f"{site} {param}".strip()
            if not key:
                continue
            row[key] = (r.get("value") or "").strip()
            if key not in headers:
                headers.append(key)
            err = (r.get("stderr") or "").strip()
            if err:
                ek = f"{key} ±"
                row[ek] = err
                if ek not in headers:
                    headers.append(ek)
        rows.append(row)
    for r in rows:
        for h in headers:
            r.setdefault(h, "")
    return headers, rows, "scope"


def numeric_headers(headers, rows) -> list:
    """Headers whose non-blank cells all parse as numbers (and at least one
    does) -- 'method' / 'structure' text columns are excluded."""
    out = []
    for h in headers:
        cells = [(r.get(h) or "").strip() for r in rows]
        filled = [c for c in cells if c]
        if not filled:
            continue
        if all(_as_float(c) is not None for c in filled):
            out.append(h)
    return out


def pair_error_columns(headers) -> dict:
    """value header -> its ± header, by :data:`ERR_SUFFIXES` (``P2O5_mol`` ->
    ``P2O5_mol_sd``, ``Vm`` -> ``Vm_u``, ``A population_pct`` -> ``… ±``)."""
    hs = [str(h) for h in headers]
    have = set(hs)
    out: dict = {}
    for h in hs:
        for suf in ERR_SUFFIXES:
            if h.endswith(suf) and len(h) > len(suf):
                base = h[:-len(suf)].rstrip()
                if base in have and base not in out:
                    out[base] = h
                break
    return out


@dataclass
class Match:
    row: int
    csv_row: int | None
    how: str                       # 'exact' | 'normalised' | 'none' | 'ambiguous'
    candidates: list = field(default_factory=list)


def propose_mapping(table: SeriesTable, csv_rows: list, sample_header: str) -> list:
    """Pair series rows with CSV rows: exact on the display name, then exact
    on the group, then normalised-key equality (name, then group), then the
    CSV name read as a sample FOLDER (``04272026_P5-Bi8-12_SS_ALP`` -> key
    ``P5-Bi8-12``, as a batch CSV written before the folder-derived names
    carries; ``how`` = ``folder``). More than one candidate at a stage ->
    ``ambiguous`` with the candidates listed and no default -- the person
    decides."""
    names = [(r.get(sample_header) or "").strip() for r in csv_rows]
    exact: dict = {}
    normed: dict = {}
    folder: dict = {}
    for i, nm in enumerate(names):
        exact.setdefault(nm, []).append(i)
        normed.setdefault(norm_key(nm), []).append(i)
        parts = scan.name_parts(nm)
        if parts.from_folder:
            folder.setdefault(norm_key(parts.key), []).append(i)
    out = []
    for k, row in enumerate(table.rows):
        found = None
        for cand_name, how in ((row.display_name, "exact"), (row.group, "exact"),
                               (row.display_name, "normalised"), (row.group, "normalised"),
                               (row.group, "folder"), (row.display_name, "folder")):
            if not cand_name:
                continue
            idx = (exact.get(cand_name) if how == "exact"
                   else normed.get(norm_key(cand_name)) if how == "normalised"
                   else folder.get(norm_key(cand_name)))
            if not idx:
                continue
            if len(idx) == 1:
                found = Match(k, idx[0], how)
            else:
                found = Match(k, None, "ambiguous", list(idx))
            break
        out.append(found or Match(k, None, "none"))
    return out


def apply_join(table: SeriesTable, csv_rows: list, sample_header: str, matches: list,
               columns: list, tag: str = "", source: str = "") -> list:
    """Write the chosen CSV columns into the table as numeric columns
    (``float`` or None), labelled ``'<tag> <header>'``; a ± partner found by
    :func:`pair_error_columns` rides along as ``err_key``. Returns the
    columns added."""
    headers = list(csv_rows[0].keys()) if csv_rows else []
    pairs = pair_error_columns(headers)
    where = {m.row: m.csv_row for m in matches if m.csv_row is not None}
    added = []
    for h in columns:
        if h == sample_header:
            continue
        err = pairs.get(h)
        col = SeriesColumn(key=h, label=f"{tag} {h}".strip(), err_key=err,
                           tag=tag, source=source)
        table.remove_column(h)
        table.columns.append(col)
        for k, r in enumerate(table.rows):
            j = where.get(k)
            r.values[h] = _as_float(csv_rows[j].get(h)) if j is not None else None
            if err:
                r.values[err] = _as_float(csv_rows[j].get(err)) if j is not None else None
        added.append(col)
    return added


# ---------------------------------------------------------------- statistics
def replicate_stats(groups: dict, y, yerr=None, x=None, xerr=None) -> dict:
    """Collapse replicates: per group the mean of ``y``; its spread is the
    sample std (ddof = 1) when two or more finite members exist, else the
    single member's own error. ``x`` is averaged the same way; its bar is
    the sample std when the members' x differ, else their RMS uncertainty
    (replicates of one composition keep the composition's own ± ).
    Returns {'labels', 'y', 'yerr', 'x', 'xerr', 'n'} (x / xerr None when no x)."""
    y = np.asarray(y, float)
    yerr = np.asarray(yerr, float) if yerr is not None else np.full(y.shape, np.nan)
    labels = list(groups)
    ym, ys, ns = [], [], []
    xm, xs = [], []
    for g in labels:
        idx = [k for k in groups[g] if k < len(y)]
        m, s, n = _mean_spread(y[idx], yerr[idx])
        ym.append(m); ys.append(s); ns.append(n)
        if x is not None:
            xx = np.asarray(x, float)[idx]
            xe = (np.asarray(xerr, float)[idx] if xerr is not None
                  else np.full(len(idx), np.nan))
            fin = np.isfinite(xx)
            if fin.sum() == 0:
                xm.append(np.nan); xs.append(np.nan)
                continue
            mean = float(np.mean(xx[fin]))
            std = float(np.std(xx[fin], ddof=1)) if fin.sum() >= 2 else 0.0
            if std > 0:
                spread = std
            else:
                e = xe[fin & np.isfinite(xe)]
                spread = float(np.sqrt(np.mean(e ** 2))) if e.size else np.nan
            xm.append(mean); xs.append(spread)
    out = {"labels": labels, "y": np.array(ym, float), "yerr": np.array(ys, float),
           "n": np.array(ns, int), "x": None, "xerr": None}
    if x is not None:
        out["x"] = np.array(xm, float)
        out["xerr"] = np.array(xs, float)
    return out


def _mean_spread(vals: np.ndarray, errs: np.ndarray) -> tuple:
    fin = np.isfinite(vals)
    n = int(fin.sum())
    if n == 0:
        return np.nan, np.nan, 0
    if n == 1:
        e = errs[fin][0]
        return float(vals[fin][0]), (float(e) if np.isfinite(e) else np.nan), 1
    return float(np.mean(vals[fin])), float(np.std(vals[fin], ddof=1)), n


def ols(x, y) -> dict:
    """Ordinary least squares y = slope x + intercept on the finite pairs:
    np.polyfit(deg 1) with the covariance errors when n >= 3 (NaN at n = 2),
    Pearson r from np.corrcoef; n < 2 -> every field NaN. Never raises."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    out = {"slope": np.nan, "intercept": np.nan, "slope_err": np.nan,
           "intercept_err": np.nan, "r": np.nan, "n": n}
    if n < 2:
        return out
    xx, yy = x[m], y[m]
    try:
        if n >= 3:
            p, cov = np.polyfit(xx, yy, 1, cov=True)
            out["slope_err"] = float(np.sqrt(abs(cov[0, 0])))
            out["intercept_err"] = float(np.sqrt(abs(cov[1, 1])))
        else:
            p = np.polyfit(xx, yy, 1)
        out["slope"], out["intercept"] = float(p[0]), float(p[1])
    except (np.linalg.LinAlgError, ValueError, TypeError):
        try:
            p = np.polyfit(xx, yy, 1)
            out["slope"], out["intercept"] = float(p[0]), float(p[1])
        except Exception:
            return out
    with np.errstate(invalid="ignore", divide="ignore"):
        try:
            r = np.corrcoef(xx, yy)[0, 1]
        except Exception:
            r = np.nan
    out["r"] = float(r) if np.isfinite(r) else np.nan
    return out


def ols_label(fit: dict) -> str:
    """'slope a ± e · r = … · n = …' for a legend."""
    s, e, r, n = fit["slope"], fit["slope_err"], fit["r"], fit["n"]
    err = f"{e:.2g}" if np.isfinite(e) else "—"
    rr = f"{r:.3f}" if np.isfinite(r) else "—"
    return f"slope {s:.3g} ± {err} · r = {rr} · n = {n}"


# ---------------------------------------------------------------- exports
def write_series_csv(table: SeriesTable, path, extra: dict | None = None) -> str:
    """The wide companion of the long batch CSV: position, name, group,
    folder, title, source_path, every series column (value then ±), then the
    ``extra`` columns (e.g. RMSD, aligned with the rows). Returns the path."""
    header = ["position", "name", "group", "folder", "title", "source_path"]
    for c in table.columns:
        header.append(c.label)
        if c.err_key:
            header.append(f"{c.label} ±")
    extra = extra or {}
    header += list(extra)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for k, r in enumerate(table.rows):
            cells = [k + 1, r.display_name, r.group, r.folder, r.title, r.source_path]
            for c in table.columns:
                cells.append(_fmt(r.values.get(c.key)))
                if c.err_key:
                    cells.append(_fmt(r.values.get(c.err_key)))
            for name, vals in extra.items():
                cells.append(_fmt(vals[k]) if k < len(vals) else "")
            w.writerow(cells)
    return str(path)


def _fmt(v) -> str:
    f = _as_float(v)
    return "" if f is None else f"{f:.6g}"


def _nan(v) -> float:
    f = _as_float(v)
    return np.nan if f is None else f


def _oxide_norm(name: str) -> str:
    s = unicodedata.normalize("NFKC", str(name)).translate(_SUBSCRIPTS).strip()
    s = re.sub(r"[\[(].*?[\])]", "", s)                     # (mol%), [wt%]
    s = re.sub(r"(?i)(mol|wt|percent|pct)", "", s)
    s = re.sub(r"[\s._\-%]+", "", s)
    return s.casefold()


_OXIDE_LOOKUP = {_oxide_norm(ox): ox for ox in DUST_OXIDES}


def canonical_oxide(header: str) -> str | None:
    """The DUST canonical oxide a header denotes (``na2o_wt%``, ``Bi2O3
    (mol%)``, ``P2O5_mol`` -> ``Na2O``, ``Bi2O3``, ``P2O5``), None for a
    non-oxide column (``Vm``, ``n_matrix``, ``R'``)."""
    return _OXIDE_LOOKUP.get(_oxide_norm(header))


def _oxide_of_column(col: SeriesColumn) -> str | None:
    ox = canonical_oxide(col.key)
    if ox:
        return ox
    label = col.label
    if col.tag and label.startswith(col.tag):
        label = label[len(col.tag):].strip()
    return canonical_oxide(label)


def dust_rows(table: SeriesTable | None, n4, n4_err=None, average: bool = False,
              names=None) -> tuple:
    """Rows of a DUST-importable composition CSV: header ``['Sample',
    <canonical oxides present>, 'N4_measured', 'N4_measured_err']``, one row
    per series row -- or per group when ``average`` (mean composition, the
    replicate mean / spread of N4 through :func:`replicate_stats`).
    ``n4`` / ``n4_err`` are fractions aligned with the rows. Returns
    (header, rows, skipped non-oxide labels). ``names`` label the rows when
    ``table`` is None (no series table: Sample + N4 only)."""
    n4 = np.asarray(n4, float)
    n4_err = (np.asarray(n4_err, float) if n4_err is not None
              else np.full(n4.shape, np.nan))
    ox_cols: list = []
    skipped: list = []
    seen: set = set()
    if table is not None:
        for col in table.columns:
            ox = _oxide_of_column(col)
            if ox and ox not in seen:
                ox_cols.append((ox, col)); seen.add(ox)
            else:
                skipped.append(col.label)
        groups = table.groups() if average else {
            r.display_name: [k] for k, r in enumerate(table.rows)}
    else:
        labels = list(names) if names is not None else [f"spectrum {k + 1}"
                                                       for k in range(len(n4))]
        groups = {lab: [k] for k, lab in enumerate(labels)}
    header = ["Sample"] + [ox for ox, _ in ox_cols] + ["N4_measured", "N4_measured_err"]
    rows = []
    for g, idx in groups.items():
        cells: list = [g]
        for _ox, col in ox_cols:
            vals = np.array([_nan(table.rows[k].values.get(col.key)) for k in idx], float)
            cells.append(_fmt(float(np.mean(vals[np.isfinite(vals)]))
                              if np.isfinite(vals).any() else None))
        st = replicate_stats({g: idx}, n4, n4_err)
        cells.append(_fmt(st["y"][0]))
        cells.append(_fmt(st["yerr"][0]))
        rows.append(cells)
    return header, rows, skipped


# ---------------------------------------------------------------- figure specs
def species_bar_spec(categories, site_labels, values, xlabel=None, title=None,
                     colors=None) -> dict:
    """A ``species_bar`` figure spec (figures.render_species_bar) from an
    (n_categories × n_sites) array of populations: one series per site, in
    site order, categories in the given (series or x) order; NaN -> 0."""
    vals = np.nan_to_num(np.asarray(values, float))
    series = []
    for j, lab in enumerate(site_labels):
        s = {"label": str(lab), "values": [float(v) for v in vals[:, j]]}
        if colors is not None and j < len(colors) and colors[j]:
            s["color"] = colors[j]
        series.append(s)
    spec = {"kind": "species_bar", "categories": [str(c) for c in categories],
            "series": series, "ylabel": "population (%)",
            "xtick_rotation": 45 if len(categories) > 4 else 0}
    if xlabel:
        spec["xlabel"] = xlabel
    if title:
        spec["title"] = title
    return spec


def ols_trace(x, fit: dict, label: str = "", color=None) -> dict:
    """The dashed OLS line over [min x, max x] as a 1d inline trace, labelled
    by :func:`ols_label` (prefixed with ``label``)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    lo, hi = (float(x.min()), float(x.max())) if x.size else (0.0, 1.0)
    xx = [lo, hi]
    yy = [fit["slope"] * v + fit["intercept"] for v in xx]
    t = {"data": {"x": xx, "y": yy},
         "label": (f"{label} " if label else "") + ols_label(fit),
         "linestyle": "--", "marker": None, "linewidth": 1.2}
    if color:
        t["color"] = color
    return t


def correlation_traces(x, y, xerr, yerr, fit: dict, label: str, color=None) -> list:
    """The marker trace (with data.xerr / data.yerr) plus the OLS line."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    data = {"x": [float(v) for v in x], "y": [float(v) for v in y]}
    if xerr is not None and np.isfinite(np.asarray(xerr, float)).any():
        data["xerr"] = [float(e) if np.isfinite(e) else 0.0 for e in np.asarray(xerr, float)]
    if yerr is not None and np.isfinite(np.asarray(yerr, float)).any():
        data["yerr"] = [float(e) if np.isfinite(e) else 0.0 for e in np.asarray(yerr, float)]
    pts = {"data": data, "label": label, "marker": "o", "linestyle": "none"}
    if color:
        pts["color"] = color
    return [pts, ols_trace(x, fit, label, color)]
