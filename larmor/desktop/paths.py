"""Where a save dialog points, and what a save may never replace.

Save dialogs used to open on whatever folder the last file dialog visited,
which after switching datasets meant the PREVIOUS dataset's folder (opening
from the Explorer never goes through a file dialog, so the remembered folder
went stale). They now seed from the data being saved: its own folder, so a
fit lands next to its spectrum by default — the same place dmfit keeps its
fits and the Explorer lists them.

Saving is allowed anywhere, including inside EXPNO folders. The only thing
that stays forbidden is REPLACING one of the acquired files themselves:
those are the measurement, and no fit output is worth clobbering them.
"""
from __future__ import annotations

from pathlib import Path

#: the files that ARE the measurement / its processing, by exact name
_INSTRUMENT_FILES = {
    "fid", "ser", "acqu", "acqus", "acqu2", "acqu2s", "acqu3", "acqu3s",
    "pulseprogram", "pdata", "1r", "1i", "2rr", "2ri", "2ir", "2ii", "3rrr",
    "proc", "procs", "proc2", "proc2s", "proc3", "proc3s", "title",
}


def is_instrument_file(target) -> bool:
    """True when writing to ``target`` would REPLACE an existing acquired
    file (fid, ser, 1r, acqus, ...). New files anywhere are fine."""
    p = Path(str(target))
    return p.exists() and p.name.lower() in _INSTRUMENT_FILES


def suggest_save_dir(data_path=None, fallback: str = "") -> str:
    """Directory to seed a save dialog with, for data loaded from data_path:
    the data's own folder when it exists, else ``fallback``."""
    if not data_path:
        return fallback
    p = Path(str(data_path))
    if not p.exists():
        return fallback
    return str(p if p.is_dir() else p.parent)
