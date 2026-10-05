"""Saved-preference reads that survive a corrupt or foreign value.

``QSettings("LARMOR", "app")`` is the Windows registry. Only LARMOR writes
it, but a value can still come back unusable: an older version's layout, a
registry cleaner, a half-written key, a type the platform plugin coerced to
a string. Every read that is parsed -- JSON libraries, recent-file lists --
goes through here, so the worst case is "the preference is lost", never
"Add line raises ValueError" or "the window does not open" (the start-up
audit for 0.16 found three ``json.loads`` of registry values with no guard
and a font size parsed with ``int()`` before the main window existed).
"""
from __future__ import annotations

import copy
import json

from PySide6.QtCore import QSettings

__all__ = ["settings", "json_setting", "list_setting", "str_setting", "font_pt_setting"]

ORG, APP = "LARMOR", "app"


def settings() -> QSettings:
    return QSettings(ORG, APP)


def font_pt_setting(default: int = 9) -> int:
    """The saved UI font size (View ▸ Text size), or ``default`` for anything
    unreadable or outside 6-24 pt. Read before the main window exists
    (app.main) and again when its menu is built."""
    try:
        v = settings().value("fontPt", default)
        pt = int(float(v)) if v not in (None, "") else default
    except (TypeError, ValueError):
        return default
    return pt if 6 <= pt <= 24 else default


def json_setting(key: str, default):
    """The JSON stored under ``key`` as a Python value of the same kind as
    ``default`` (a dict for a dict, a list for a list); a copy of
    ``default`` when the value is missing, unparsable or of another kind."""
    raw = settings().value(key, None)
    if raw in (None, ""):
        return copy.deepcopy(default)
    if isinstance(raw, (dict, list)):
        value = raw
    else:
        try:
            value = json.loads(str(raw))
        except (TypeError, ValueError):
            return copy.deepcopy(default)
    if not isinstance(value, type(default)):
        return copy.deepcopy(default)
    return value


def list_setting(key: str) -> list:
    """A list of strings stored under ``key`` -- a single string becomes a
    one-item list, anything else an empty one (QSettings hands a one-item
    list back as a bare string on some platforms)."""
    raw = settings().value(key, [])
    if raw is None or raw == "":
        return []
    if isinstance(raw, (list, tuple)):
        return [str(x) for x in raw if x not in (None, "")]
    if isinstance(raw, str):
        return [raw]
    return []


def str_setting(key: str, default: str, allowed=None) -> str:
    """A string preference, ``default`` when missing or (with ``allowed``)
    not one of the known values."""
    raw = settings().value(key, default)
    text = str(raw) if raw is not None else default
    if allowed is not None and text not in allowed:
        return default
    return text
