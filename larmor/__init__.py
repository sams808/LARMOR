"""LARMOR: open successor to dmfit for solid-state NMR lineshape fitting.

Reuse-first design: the physics comes from mrsimulator + lmfit; LARMOR adds
ingestion (Bruker/TopSpin, legacy dmfit .fxmla), a diffable JSON recipe format,
and a fast kernel-based fitting engine with real parameter uncertainties.
"""

__version__ = "0.1.0"

from larmor.recipe import Param, SiteModel, Recipe  # noqa: F401
