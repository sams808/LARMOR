"""Multi-dataset fitting: one physical model, several spectra.

The flagship use case is the multi-field 1D fit: the same sites measured at
two or more B0 fields. Quadrupolar widths scale as 1/B0 while chemical-shift
dispersion is constant in ppm, so a simultaneous fit lifts the Cq/delta_iso
degeneracy that a single field cannot resolve. The same machinery accepts any
list of (recipe, data) pairs -- an MQMAS projection alongside a 1D, spectra of
a composition series with shared widths, etc.

Sharing model: dataset 0 is the master. For every parameter name in `share`,
dataset k>0's parameter is expression-linked to dataset 0's (full error
propagation, exactly like intra-recipe links).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import lmfit

from larmor import models as model_registry
from larmor.engine import make_context, simulate_site
from larmor.fit import ConstraintError, _EXPR_REF
from larmor.recipe import Recipe

#: parameter names it makes physical sense to share across fields
DEFAULT_SHARE = ("isotropic_chemical_shift_ppm", "sigma_Cq_MHz", "Cq_MHz",
                 "eta", "eta_q", "eps", "zeta_ppm", "eta_cs", "shift_fwhm_ppm")


@dataclass
class MultiFitResult:
    recipes: list[Recipe]
    lmfit_result: lmfit.minimizer.MinimizerResult
    rmsd: list[float]
    per_dataset: list[dict]          # {"x":..., "y_fit":..., "per_site":...}

    @property
    def report(self) -> str:
        return lmfit.fit_report(self.lmfit_result)


def _key(site, pname):
    return model_registry.get(site.model).key_of(pname)


def _pname(k: int, i: int, site, pname: str) -> str:
    return f"d{k}_s{i}_{_key(site, pname)}"


def fit_multi(entries: list[tuple[Recipe, np.ndarray, np.ndarray]],
              share: tuple[str, ...] = DEFAULT_SHARE,
              windows: list[tuple[float, float] | None] | None = None,
              ) -> MultiFitResult:
    """Simultaneously fit several (recipe, ppm, amp) datasets.

    All recipes must have the same number of sites with matching models
    (site j of every dataset is the same physical site).
    """
    if len(entries) < 2:
        raise ValueError("multi-fit needs at least two datasets")
    n_sites = len(entries[0][0].sites)
    for r, _, _ in entries[1:]:
        if len(r.sites) != n_sites:
            raise ValueError("all recipes must have the same number of sites")
        for j in range(n_sites):
            if r.sites[j].model != entries[0][0].sites[j].model:
                raise ValueError(
                    f"site {j} model differs between datasets "
                    f"({r.sites[j].model} vs {entries[0][0].sites[j].model})")

    ctxs = [make_context(r, exp_ppm=ppm) for r, ppm, _ in entries]

    # selections per dataset (zones > window > full range)
    sels, weights = [], []
    for k, (r, ppm, amp) in enumerate(entries):
        zones = [z for z in (r.fit_zones or []) if z and len(z) == 2]
        if zones:
            sel = np.zeros(ppm.shape, bool)
            for zhi, zlo in zones:
                sel |= (ppm >= min(zhi, zlo)) & (ppm <= max(zhi, zlo))
        elif windows and windows[k]:
            hi, lo = max(windows[k]), min(windows[k])
            sel = (ppm >= lo) & (ppm <= hi)
        else:
            sel = np.ones(ppm.shape, bool)
        sels.append(sel)
        weights.append(1.0 / (np.abs(amp[sel]).max() or 1.0))

    # ---------------- parameters
    params = lmfit.Parameters()
    for k, (r, _, _) in enumerate(entries):
        for i, site in enumerate(r.sites):
            for pname, p in site.params.items():
                params.add(_pname(k, i, site, pname), value=p.value,
                           vary=p.vary,
                           min=p.min if p.min is not None else -np.inf,
                           max=p.max if p.max is not None else np.inf)
    # intra-recipe links (translated with the dataset prefix)
    for k, (r, _, _) in enumerate(entries):
        for i, site in enumerate(r.sites):
            for pname, p in site.params.items():
                if p.expr:
                    trans = _EXPR_REF.sub(
                        lambda m: _pname(k, int(m.group(1)),
                                         r.sites[int(m.group(1))], m.group(2)),
                        p.expr)
                    params[_pname(k, i, site, pname)].expr = trans
    # cross-dataset sharing: dataset 0 is the master
    for k, (r, _, _) in enumerate(entries):
        if k == 0:
            continue
        for i, site in enumerate(r.sites):
            for pname in site.params:
                if pname in share and not site.params[pname].expr:
                    params[_pname(k, i, site, pname)].expr = \
                        _pname(0, i, entries[0][0].sites[i], pname)
    try:
        params.update_constraints()
    except Exception as exc:
        raise ConstraintError(f"invalid multi-fit constraints: {exc}") from exc

    # ---------------- residual
    def apply_params(p):
        for k, (r, _, _) in enumerate(entries):
            for i, site in enumerate(r.sites):
                for pname in site.params:
                    lp = p[_pname(k, i, site, pname)]
                    site.params[pname].value = float(lp.value)
                    site.params[pname].stderr = (float(lp.stderr)
                                                 if lp.stderr is not None else None)

    def residual(p):
        apply_params(p)
        chunks = []
        for k, (r, ppm, amp) in enumerate(entries):
            per_site = [simulate_site(s, ctxs[k]) for s in r.sites]
            total = np.sum(per_site, axis=0)
            yi = np.interp(ppm[sels[k]], ctxs[k].x_ppm, total)
            chunks.append((yi - amp[sels[k]]) * weights[k])
        return np.concatenate(chunks)

    result = lmfit.minimize(residual, params, method="least_squares")
    if not result.errorbars:
        retry = lmfit.minimize(residual, result.params.copy(), method="leastsq")
        if retry.errorbars:
            result = retry
    apply_params(result.params)

    # ---------------- outputs
    rmsds, per_dataset = [], []
    for k, (r, ppm, amp) in enumerate(entries):
        per_site = [simulate_site(s, ctxs[k]) for s in r.sites]
        total = np.sum(per_site, axis=0)
        yi = np.interp(ppm[sels[k]], ctxs[k].x_ppm, total)
        yw = amp[sels[k]]
        rmsds.append(float(np.sqrt(np.mean((yi - yw) ** 2)) / (np.abs(yw).max() or 1)))
        r.fit_rmsd = rmsds[-1]
        per_dataset.append({"x": ctxs[k].x_ppm, "y_fit": total,
                            "per_site": per_site})
    return MultiFitResult(recipes=[e[0] for e in entries],
                          lmfit_result=result, rmsd=rmsds,
                          per_dataset=per_dataset)
