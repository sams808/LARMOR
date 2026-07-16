"""Chemical-shift-anisotropy model: spin-1/2 powder pattern with physical
spinning sidebands (replaces dmfit's ad-hoc 'ss band' lines)."""
from __future__ import annotations

import numpy as np

from larmor.models.base import Model, ParamDef, SimContext, register


def _render_csa(v: dict, ctx: SimContext) -> np.ndarray:
    from larmor.models._singlesite import render_single_site

    return render_single_site(v, ctx, zeta_key="zeta_ppm", eta_cs_key="eta",
                              ct_only=False, n_ssb=32)


register(Model(
    name="csa_mas",
    label="CSA powder (MAS/static)",
    description="Shielding-anisotropy powder pattern with physical spinning "
                "sidebands (spin rate 0 = static). Replaces manual 'ss band' "
                "lines.",
    params=(
        ParamDef("isotropic_chemical_shift_ppm", "pos", 0.0, "ppm",
                 "isotropic chemical shift"),
        ParamDef("zeta_ppm", "zeta", 50.0, "ppm",
                 "shielding anisotropy (Haeberlen zeta)", min=-1000.0, max=1000.0),
        ParamDef("eta", "eta", 0.3, "", "shielding asymmetry", min=0.0, max=1.0),
        ParamDef("shift_fwhm_ppm", "fwhm", 2.0, "ppm", "Gaussian broadening",
                 min=0.05),
        ParamDef("amplitude", "amp", 1.0, "", "peak height", min=0.0),
    ),
    render=_render_csa,
))
