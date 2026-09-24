"""Display names for sample folders and EXPNOs, and folder renames on disk.

A Bruker dataset is named by its folder (``01192026_SR31649_Base0Ca_SS_ALP``,
``2702``) and that name is read everywhere -- the Explorer tree, the window
title, the Datasets dock, the recipe's ``sample`` (through
``larmor.io.scan.sample_name`` / ``sample_label``), batch and sequential
labels. The Explorer's right-click **Rename…** offers two ways to change it:

* an **alias** -- a display name kept by LARMOR only, in
  ``%LOCALAPPDATA%/LARMOR/aliases.json`` keyed by the folder path
  (``LARMOR_ALIASES`` overrides, tests point it at a temporary file). The
  desktop writes it; the Qt-free core reads it, so a batch fit run from the
  CLI labels the sample the same way the window does;
* a **rename on disk** -- ``os.rename`` of the sample folder or EXPNO folder
  after the checks in :func:`check_rename` (the target must not exist, no
  file inside may be open in LARMOR, an EXPNO stays a number so the TopSpin
  layout holds), recorded in the append-only
  ``%LOCALAPPDATA%/LARMOR/rename_log.jsonl`` (``LARMOR_RENAME_LOG``
  overrides) the way the referencing audit records its corrections, so a
  rename is traceable years later. Aliases keyed under the old path follow
  the folder.

Everything here is plain Python: no Qt.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
from pathlib import Path

__all__ = [
    "aliases_path", "rename_log_path", "load", "alias_for", "set_alias",
    "display_name", "window_label", "is_expno", "RenameError",
    "check_rename", "rename_folder", "append_rename_log",
]


# ------------------------------------------------------------------ paths
def _app_dir() -> Path:
    base = os.environ.get("LOCALAPPDATA") or str(Path.home())
    return Path(base) / "LARMOR"


def aliases_path() -> Path:
    """``%LOCALAPPDATA%/LARMOR/aliases.json`` (``LARMOR_ALIASES`` overrides)."""
    env = os.environ.get("LARMOR_ALIASES")
    return Path(env) if env else _app_dir() / "aliases.json"


def rename_log_path() -> Path:
    """Append-only record of every rename on disk:
    ``%LOCALAPPDATA%/LARMOR/rename_log.jsonl`` (``LARMOR_RENAME_LOG``
    overrides). Never truncated by the application."""
    env = os.environ.get("LARMOR_RENAME_LOG")
    return Path(env) if env else _app_dir() / "rename_log.jsonl"


def _key(path) -> str:
    """The store key of a folder: absolute, normalised, case-folded the way
    the platform compares paths (Windows is case-insensitive)."""
    return os.path.normcase(os.path.normpath(os.path.abspath(str(path))))


# ------------------------------------------------------------------ store
_cache: dict = {"path": None, "mtime": None, "data": {}}


def load() -> dict[str, str]:
    """``{folder key: alias}`` -- re-read only when the file changed, so the
    lookups in ``scan.sample_name`` stay cheap for a batch of hundreds."""
    p = aliases_path()
    try:
        mtime = p.stat().st_mtime_ns
    except OSError:
        _cache.update(path=str(p), mtime=None, data={})
        return {}
    if _cache["path"] == str(p) and _cache["mtime"] == mtime:
        return _cache["data"]
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
        data = {str(k): str(v) for k, v in raw.items()
                if isinstance(v, str) and v.strip()} if isinstance(raw, dict) else {}
    except (OSError, ValueError):
        data = {}
    _cache.update(path=str(p), mtime=mtime, data=data)
    return data


def _write(data: dict[str, str]) -> None:
    p = aliases_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=1, sort_keys=True),
                   encoding="utf-8")
    os.replace(tmp, p)
    _cache.update(path=None, mtime=None, data={})      # force a re-read


def alias_for(path) -> str:
    """The alias of a folder, ``""`` when it has none."""
    if not path:
        return ""
    return load().get(_key(path), "")


def set_alias(path, name: str | None) -> None:
    """Give a folder a display name; an empty name, or the folder's own
    name, removes the alias."""
    name = (name or "").strip()
    data = dict(load())
    k = _key(path)
    if name and name != Path(str(path)).name:
        data[k] = name
    else:
        data.pop(k, None)
    _write(data)


def display_name(path) -> str:
    """The alias when there is one, else the folder's own name."""
    return alias_for(path) or Path(str(path)).name


def _expno_dir_of(path) -> Path | None:
    """The first ancestor (or the path itself) that holds an ``acqus``."""
    p = Path(str(path))
    for c in [p] + list(p.parents)[:5]:
        try:
            if (c / "acqus").exists():
                return c
        except OSError:
            return None
    return None


def window_label(path) -> str:
    """What the title bar shows for a data path: the EXPNO's alias, else
    ``<sample alias> · <expno>``, else -- unchanged from before aliases
    existed -- the file's own name."""
    p = Path(str(path))
    expno = _expno_dir_of(p)
    if expno is None:
        return p.name
    a_exp = alias_for(expno)
    if a_exp:
        return a_exp
    a_smp = alias_for(expno.parent)
    if a_smp:
        return f"{a_smp} · {expno.name}"
    return p.name


def is_expno(folder) -> bool:
    """A TopSpin experiment folder: it holds an ``acqus``."""
    try:
        return (Path(str(folder)) / "acqus").exists()
    except OSError:
        return False


# ----------------------------------------------------------------- rename
class RenameError(ValueError):
    """A rename on disk that must not happen; the message says why."""


def _within(path: str, folder_key: str) -> bool:
    k = _key(path)
    return k == folder_key or k.startswith(folder_key.rstrip(os.sep) + os.sep)


def check_rename(old, new_name: str, open_paths=()) -> Path:
    """Validate a rename of the folder ``old`` to ``new_name`` (a bare name
    in the same parent) and return the target path. Raises
    :class:`RenameError` when the source is missing, the name is empty or
    not a single path component, an EXPNO would stop being a number, the
    target exists, or a file inside the folder is open in LARMOR."""
    src = Path(str(old))
    if not src.is_dir():
        raise RenameError(f"{src} is not a folder")
    name = (new_name or "").strip()
    if not name:
        raise RenameError("the new name is empty")
    if name in (".", "..") or any(ch in name for ch in "/\\") or name != Path(name).name:
        raise RenameError(f"{name!r} is not a plain folder name")
    if is_expno(src) and not name.isdigit():
        raise RenameError(f"an EXPNO must stay a number (TopSpin layout); "
                          f"{name!r} is not -- use a display name instead")
    if name == src.name:
        raise RenameError("the name is unchanged")
    target = src.with_name(name)
    if target.exists():
        raise RenameError(f"{target} already exists")
    old_key = _key(src)
    busy = [p for p in open_paths if p and _within(str(p), old_key)]
    if busy:
        raise RenameError(f"{Path(busy[0]).name} inside {src.name} is open in "
                          "LARMOR -- close that workspace first")
    return target


def append_rename_log(old, new, kind: str, *, alias_moved: bool = False,
                      path=None) -> Path:
    """Append one record of a rename on disk (old and new absolute paths, the
    folder kind, whether an alias followed), so it can be traced later."""
    from larmor import __version__

    p = Path(path) if path else rename_log_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    record = {"time": _dt.datetime.now().isoformat(timespec="seconds"),
              "larmor": __version__, "action": "rename", "kind": kind,
              "old": os.path.abspath(str(old)), "new": os.path.abspath(str(new)),
              "alias_moved": bool(alias_moved)}
    with open(p, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
    return p


def rename_folder(old, new_name: str, *, open_paths=()) -> Path:
    """Rename a sample folder or an EXPNO on disk after :func:`check_rename`,
    move the aliases keyed under it (an alias equal to the new name is
    dropped -- the folder now says it itself) and log the rename. Returns
    the new path."""
    src = Path(str(old))
    target = check_rename(src, new_name, open_paths)
    kind = "expno" if is_expno(src) else "sample"
    old_key = _key(src)
    data = dict(load())
    moved = {}
    for k, v in list(data.items()):
        if _within(k, old_key):
            moved[k] = data.pop(k)
    os.rename(src, target)
    new_key = _key(target)
    for k, v in moved.items():
        nk = new_key + k[len(old_key):]
        if nk == new_key and v == target.name:
            continue                    # renamed to its alias: alias redundant
        data[nk] = v
    if moved:
        _write(data)
    append_rename_log(src, target, kind, alias_moved=bool(moved))
    return target
