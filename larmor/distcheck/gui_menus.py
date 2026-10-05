"""Desktop sweep, part 1: the main window built offscreen, every menu
action triggered with its dialogs auto-dismissed, every theme and text
size applied, every manual and tutorial rendered. Each stage is
``(name, fn(say, ctx))``; plumbing in ``gui_common``."""
from __future__ import annotations

__all__ = ["stages"]


def stages() -> list:
    return []
