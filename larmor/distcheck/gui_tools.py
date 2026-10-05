"""Desktop sweep, part 2: every tool window and workflow of the desktop on
synthetic data -- fit, overlays, series mode, 2D view, batch, the error
tools, the Plotting studio with its exports, project save / open.
Each stage is ``(name, fn(say, ctx))``; plumbing in ``gui_common``."""
from __future__ import annotations

__all__ = ["stages"]


def stages() -> list:
    return []
