"""External-spectrum component: fit a measured trace (e.g. a background, an
impurity, or a reference lineshape) as a scalable, shiftable basis function.

The reference data lives on the site (SiteModel.ref); the registry render here
is only a placeholder so the model self-describes its parameters to the UI and
the fit machinery. The real rendering is done in engine.simulate_site, which
has access to the site and its reference trace.
"""
from __future__ import annotations

import numpy as np

from larmor.models.base import Model, ParamDef, register


def _placeholder(values: dict, ctx) -> np.ndarray:
    return np.zeros_like(ctx.x_ppm)


register(Model(
    name="function",
    label="Function fit",
    description="a user y(x; a,b,c,d) expression (ssNake Function fit), evaluated "
                "with numpy — x is the ppm axis, amplitude scales the result. "
                "Set the expression when adding the line.",
    params=(
        ParamDef("amplitude", "amp", 1.0, "", "overall scale", min=0.0),
        ParamDef("a", "a", 0.0, "", "free parameter a"),
        ParamDef("b", "b", 0.0, "", "free parameter b"),
        ParamDef("c", "c", 0.0, "", "free parameter c"),
        ParamDef("d", "d", 0.0, "", "free parameter d"),
    ),
    render=_placeholder,
))


register(Model(
    name="spectrum",
    label="Spectrum (background)",
    description="an external measured spectrum used as a fit component "
                "(background, impurity, or reference), scaled by amplitude and "
                "rigidly shiftable in ppm",
    params=(
        ParamDef("amplitude", "amp", 1.0, "", "scale factor", min=0.0),
        ParamDef("shift_ppm", "sh", 0.0, "ppm", "rigid ppm shift"),
    ),
    render=_placeholder,
))
