"""Provenance of a fit: the software that produced it, the exact data file
it was fitted on, and what changed since.

Three things a published fit must be able to answer years later:

* **which software** -- ``software_stamp()`` is written into every fitted
  recipe by ``fit.fit`` (``recipe.software``): LARMOR version and git commit
  when run from a checkout, mrsimulator, lmfit, numpy, scipy, nmrglue,
  Python, platform and the fitting time. The versions come from
  :func:`larmor.methods.software_versions` (one reader for the bundle README
  and the stamp) extended with the two libraries the lineshape and the
  Fourier path also depend on;
* **which bytes** -- ``source_sha256`` hashes the data file a recipe was
  fitted on (``pdata/<procno>/1r``, the ``fid`` behind a File > Open FID
  fit, a CSV, a Varian ``fid``), through the one hasher ``recipe.sha256_of``;
  ``verify_source`` re-hashes on reopen and compares the stored SR with the
  file's, so a TopSpin reprocess or a later ``sr`` shows up as a status-bar
  line, never as a modal;
* **which read-out** -- ``READOUT_CHANGES`` records the LARMOR versions at
  which a stored number started to be *read* differently (the stored
  parameters never change, the derived text does); ``version_notes`` voices
  the entries that lie between the fitting version and the current one and
  apply to the recipe's models. A mere version difference is not voiced:
  a warning on every reopen trains students to ignore warnings.

``carry_source`` is the per-spectrum copy the batch / sequential / CLI
recipe builders use so a saved individual fit keeps the source kind, hash
and acquisition block of the spectrum it was made from (the 32 published
Final2 recipes were written with an empty ``source_kind`` and
``source_sha256`` while the README promised path + hash).

Qt-free, stdlib + ``larmor.recipe``; ``larmor.methods`` and
``larmor.io.bruker`` are imported inside functions. The git commit is read
from the ``.git`` files (a worktree's ``gitdir:`` pointer and ``commondir``
are followed) -- no subprocess from a worker thread or the frozen exe, which
has no ``.git`` and records ``""``.
"""
from __future__ import annotations

import datetime as _dt
import re
from functools import lru_cache
from pathlib import Path

__all__ = [
    "SOURCE_KEYS", "carry_source", "software_stamp", "git_commit",
    "parse_version", "data_file_for", "source_sha256", "SOURCE_CHANGED_PREFIX",
    "REREFERENCED_PREFIX", "verify_source", "CZJZEK_FAMILY", "READOUT_CHANGES",
    "version_notes",
]

#: the source-identity fields a per-spectrum recipe copies from the loaded
#: record (source_path is set by the caller: EXPNO for a fid replay, else path)
SOURCE_KEYS: tuple[str, ...] = ("source_kind", "source_sha256", "acquisition")

_DEFAULTS = {"source_kind": "", "source_sha256": "", "acquisition": {}}

#: status-bar prefixes (tests and the manuals pin them)
SOURCE_CHANGED_PREFIX = "source data changed"
REREFERENCED_PREFIX = "re-referenced since the fit"

#: |stored SR - file SR| above this (Hz) means the ppm axis moved
SR_TOL_HZ = 0.5

#: models whose read-outs the 0.13.0 entry concerns
CZJZEK_FAMILY: tuple[str, ...] = ("czjzek", "czjzek_d", "czjzek_corr",
                                  "ext_czjzek", "csa_czjzek")

_HEX40 = re.compile(r"^[0-9a-f]{40}$")


# ------------------------------------------------------------ carry-over
def carry_source(rec: dict, *, raw_expno=None) -> dict:
    """``{source_kind, source_sha256, acquisition}`` of a loaded recipe dict
    (defaults for missing or empty values), ready to spread into the
    per-spectrum recipe a batch / sequential / CLI builder assembles. With
    ``raw_expno`` (a batch reprocessed from its fids) the hash is that of
    ``<EXPNO>/fid`` and the block's ``data_file`` says so."""
    rec = rec or {}
    out = {}
    for k in SOURCE_KEYS:
        v = rec.get(k)
        out[k] = v if v else (dict(_DEFAULTS[k]) if isinstance(_DEFAULTS[k], dict)
                              else _DEFAULTS[k])
    if raw_expno:
        fid = Path(str(raw_expno)) / "fid"
        if fid.is_file():
            out["source_sha256"] = source_sha256(fid)
            if out["acquisition"]:
                out["acquisition"] = {**out["acquisition"], "data_file": "fid"}
    return out


# ------------------------------------------------------------- versions
@lru_cache(maxsize=1)
def _versions() -> dict:
    """The version block, read once per process: ``methods.software_versions``
    (LARMOR, mrsimulator, lmfit, numpy, python, git_commit) plus scipy,
    nmrglue and the platform. A git commit the subprocess reader could not
    give (no git binary, a frozen exe) is read from the ``.git`` files."""
    import platform

    from larmor import methods

    v = dict(methods.software_versions())
    v["scipy"] = methods._dist_version("scipy")
    v["nmrglue"] = methods._dist_version("nmrglue")
    v["platform"] = platform.platform()
    if not v.get("git_commit"):
        v["git_commit"] = git_commit()
    return v


def software_stamp() -> dict:
    """The block ``fit.fit`` writes into ``recipe.software``: the cached
    versions plus ``fitted`` (local time, ISO seconds). Never raises."""
    try:
        out = dict(_versions())
    except Exception:                                    # noqa: BLE001
        out = {"larmor": "", "git_commit": "", "mrsimulator": "", "lmfit": "",
               "numpy": "", "scipy": "", "nmrglue": "", "python": "", "platform": ""}
        try:
            from larmor import __version__
            out["larmor"] = __version__
        except Exception:                                # noqa: BLE001
            pass
    out["fitted"] = _dt.datetime.now().isoformat(timespec="seconds")
    return out


def _read(p: Path) -> str:
    try:
        return p.read_text(encoding="utf-8", errors="replace").strip()
    except OSError:
        return ""


def _git_dir(root: Path) -> tuple[Path | None, Path | None]:
    """(gitdir, commondir) of a checkout root: ``.git`` may be a directory or
    a worktree's ``gitdir: <path>`` file; a worktree gitdir names the shared
    object store in its ``commondir`` file."""
    dot = root / ".git"
    if dot.is_dir():
        return dot, dot
    if dot.is_file():
        txt = _read(dot)
        if txt.lower().startswith("gitdir:"):
            g = Path(txt.split(":", 1)[1].strip())
            if not g.is_absolute():
                g = (root / g).resolve()
            if g.is_dir():
                common = g
                c = _read(g / "commondir")
                if c:
                    cp = Path(c)
                    common = (g / cp).resolve() if not cp.is_absolute() else cp
                return g, common
    return None, None


def _resolve_ref(ref: str, dirs: list[Path]) -> str:
    """The 40-hex commit ``refs/heads/x`` points at: the loose ref file in
    any of ``dirs``, else a ``packed-refs`` line."""
    for d in dirs:
        txt = _read(d / ref)
        if _HEX40.match(txt):
            return txt
    for d in dirs:
        packed = _read(d / "packed-refs")
        for line in packed.splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1] == ref and _HEX40.match(parts[0]):
                return parts[0]
    return ""


def git_commit(package_dir=None) -> str:
    """The checked-out commit (40 hex) of the checkout ``larmor/`` lives in,
    by file reads only: ``.git/HEAD`` -> ``ref: refs/heads/x`` -> that loose
    ref or a ``packed-refs`` line; a detached HEAD is the hash itself; a
    ``.git`` FILE (worktree) is followed through ``gitdir:`` and
    ``commondir``. ``""`` when anything is missing -- an installed copy, the
    frozen exe. ``package_dir`` is the folder holding ``.git`` (default: the
    parent of the ``larmor`` package)."""
    try:
        if package_dir is None:
            import larmor

            root = Path(larmor.__file__).resolve().parents[1]
        else:
            root = Path(package_dir)
        gitdir, common = _git_dir(root)
        if gitdir is None:
            return ""
        head = _read(gitdir / "HEAD")
        if not head:
            return ""
        if _HEX40.match(head):
            return head
        if head.startswith("ref:"):
            ref = head.split(":", 1)[1].strip()
            dirs = [gitdir] + ([common] if common is not None and common != gitdir else [])
            return _resolve_ref(ref, dirs)
        return ""
    except Exception:                                    # noqa: BLE001
        return ""


def parse_version(s) -> tuple[int, ...]:
    """'0.13.0' -> (0, 13, 0); non-numeric tails are dropped ('0.14.0rc1' ->
    (0, 14, 0)); an empty or unparsable string -> ()."""
    out = []
    for part in str(s or "").split("."):
        m = re.match(r"\d+", part)
        if not m:
            break
        out.append(int(m.group(0)))
    return tuple(out)


# ------------------------------------------------------------- the bytes
def data_file_for(recipe_dict: dict, from_raw: bool = False) -> Path | None:
    """The data file a recipe was (or would be) fitted on: the acquisition
    block's ``data_file`` under ``source_path`` when present; else any Bruker
    ``source_path`` -> ``expno/pdata/<procno>/1r`` (``expno/fid`` when the
    recipe replays from the raw fid or ``from_raw`` is set); a CSV / txt /
    dat / fxmla file -> itself; a Varian ``.fid`` folder -> its ``fid``.
    None when unresolvable. Never raises."""
    try:
        rec = recipe_dict or {}
        src = str(rec.get("source_path") or "")
        if not src:
            return None
        p = Path(src)
        acq = rec.get("acquisition") or {}
        raw = from_raw or bool(rec.get("processing_from_raw"))
        from larmor.io import bruker

        if bruker.is_expno(p):
            if raw:
                f = p / "fid"
                return f if f.is_file() else None
            df = str(acq.get("data_file") or "")
            if df:
                f = p / Path(df)
                if f.is_file():
                    return f
            procno = int(acq.get("procno") or 1)
            f = p / "pdata" / str(procno) / "1r"
            return f if f.is_file() else None
        suf = p.suffix.lower()
        if p.is_file() and suf in (".csv", ".txt", ".dat", ".fxmla", ".fxml"):
            return p
        if p.is_dir() and (p / "fid").is_file() and (p / "procpar").is_file():
            return p / "fid"                                 # Varian
        try:
            ref = bruker.resolve(p)
        except (ValueError, FileNotFoundError, OSError):
            return p if p.is_file() else None
        if raw or ref.target in ("fid", "ser"):
            f = ref.expno / ("ser" if ref.target == "ser" else "fid")
        else:
            f = ref.expno / "pdata" / str(ref.procno) / ref.target
        return f if f.is_file() else None
    except Exception:                                    # noqa: BLE001
        return None


def source_sha256(path_or_recipe, from_raw: bool = False) -> str:
    """SHA-256 of a data file (a path) or of the data file behind a recipe
    dict (:func:`data_file_for`); ``""`` on any failure so loading never
    breaks. The one hasher is :func:`larmor.recipe.sha256_of`."""
    from larmor.recipe import sha256_of

    try:
        if isinstance(path_or_recipe, dict):
            f = data_file_for(path_or_recipe, from_raw)
        else:
            f = Path(str(path_or_recipe))
        if f is None or not Path(f).is_file():
            return ""
        return sha256_of(f)
    except Exception:                                    # noqa: BLE001
        return ""


def _sr_explained(stored: dict, fresh_sr: float) -> bool:
    """True when the recipe itself explains why its SR differs from the
    file's: the referencing audit re-referenced it (``provenance
    ['referencing']`` keeps the old value, which is the file's) or an ``sr``
    processing step moves the axis on replay."""
    prov = (stored.get("provenance") or {}).get("referencing") or {}
    try:
        old = float(prov.get("old_sr_hz"))
        if abs(old - fresh_sr) <= SR_TOL_HZ:
            return True
    except (TypeError, ValueError):
        pass
    return any((op or {}).get("op") == "sr" for op in (stored.get("processing") or []))


def verify_source(stored: dict, fresh: dict | None) -> list[str]:
    """Reload checks of a saved recipe against the files it points at:

    * the data file's SHA-256 against ``stored['source_sha256']`` (silent when
      the recipe carries no hash -- every pre-0.14 recipe -- or the file is
      gone: the loader reports a missing source itself);
    * ``stored['sr_hz']`` against the freshly read ``fresh['sr_hz']`` (the
      ppm axis moved -- a TopSpin ``sr`` after the fit leaves the 1r bytes
      unchanged, which a hash alone misses), unless the recipe explains the
      difference itself.

    Returns status-bar lines; [] when nothing changed. Never raises."""
    out: list[str] = []
    try:
        if not stored:
            return out
        want = str(stored.get("source_sha256") or "")
        if want:
            f = data_file_for(stored)
            if f is not None:
                have = source_sha256(f)
                if have and have != want:
                    out.append(
                        f"{SOURCE_CHANGED_PREFIX} (SHA-256 differs): {Path(f).name} was "
                        "reprocessed or replaced after this fit — refit before publishing")
        if fresh:
            try:
                s_sr = float(stored.get("sr_hz") or 0.0)
                f_sr = float(fresh.get("sr_hz") or 0.0)
            except (TypeError, ValueError):
                return out
            if abs(s_sr - f_sr) > SR_TOL_HZ and not _sr_explained(stored, f_sr):
                out.append(f"{REREFERENCED_PREFIX} (SR {s_sr:.2f} → {f_sr:.2f} Hz): "
                           "the ppm axis moved")
    except Exception:                                    # noqa: BLE001
        return out
    return out


# --------------------------------------------------------- read-out changes
def _has_czjzek(rec: dict) -> bool:
    return any((s or {}).get("model") in CZJZEK_FAMILY
               for s in (rec.get("sites") or []))


#: (version at which the read-out changed, applies(recipe_dict) -> bool, note)
#: -- wording from docs/workplan-status.md (the 0.13.0 release notes)
READOUT_CHANGES: list[tuple[str, object, str]] = [
    ("0.13.0", _has_czjzek,
     "Czjzek read-outs (P(C_Q) dialog, √⟨P_Q²⟩, Methods text) changed at 0.13.0; "
     "the stored σ is unchanged — re-copy tables and Methods text"),
]


def version_notes(recipe_dict: dict, current: str | None = None) -> list[str]:
    """'fitted with LARMOR <v>: <note>' for every READOUT_CHANGES entry with
    ``fitted < change <= current`` whose predicate holds for the recipe; []
    when the recipe carries no software block, was fitted with the current
    version, or no entry applies."""
    rec = recipe_dict or {}
    sw = rec.get("software") or {}
    fitted = parse_version(sw.get("larmor"))
    if not fitted:
        return []
    if current is None:
        from larmor import __version__ as current
    cur = parse_version(current)
    if not cur or fitted >= cur:
        return []
    out = []
    for ver, applies, note in READOUT_CHANGES:
        v = parse_version(ver)
        if fitted < v <= cur:
            try:
                ok = bool(applies(rec))
            except Exception:                            # noqa: BLE001
                ok = False
            if ok:
                out.append(f"fitted with LARMOR {sw.get('larmor')}: {note}")
    return out
