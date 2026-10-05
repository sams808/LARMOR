"""Engine stages of the distribution check: every Qt-free capability on
synthetic data (readers, processing, fitting, 2D, batch, series, error
tools, exports, provenance, help rendering ...). Each stage is
``(name, fn(say, ctx))``; see ``larmor.distcheck`` for the contract."""
from __future__ import annotations

__all__ = ["stages"]


def stages() -> list:
    return []
