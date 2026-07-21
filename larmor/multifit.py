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
    """Simultaneously fit several 1D (recipe, ppm, amp) datasets.

    Backward-compatible wrapper around :func:`fit_cofit`.
    """
    conv = [(r, (np.asarray(ppm), np.asarray(amp))) for r, ppm, amp in entries]
    return fit_cofit(conv, share=share, windows=windows)


def _prepare_1d(rec: Recipe, ppm, amp, window):
    ctx = make_context(rec, exp_ppm=ppm)
    zones = [z for z in (rec.fit_zones or []) if z and len(z) == 2]
    if zones:
        sel = np.zeros(ppm.shape, bool)
        for zhi, zlo in zones:
            sel |= (ppm >= min(zhi, zlo)) & (ppm <= max(zhi, zlo))
    elif window:
        hi, lo = max(window), min(window)
        sel = (ppm >= lo) & (ppm <= hi)
    else:
        sel = np.ones(ppm.shape, bool)
    exp = amp[sel]
    return {"kind": "1d", "r": rec, "ppm": ppm, "ctx": ctx, "sel": sel,
            "exp": exp, "w": 1.0 / (np.abs(exp).max() or 1.0)}


def _prepare_2d(rec: Recipe, data2d, method: str):
    from scipy.interpolate import RegularGridInterpolator

    from larmor import twod

    d = data2d.normalized()
    kernel = twod.build_mqmas_kernel(
        rec.nucleus or d.nucleus, rec.larmor_frequency_MHz or d.larmor_MHz,
        f2_window=(float(d.f2_ppm.max()), float(d.f2_ppm.min())),
        f1_window=(float(d.f1_ppm.max()), float(d.f1_ppm.min())), method=method)
    interp = RegularGridInterpolator((d.f1_ppm, d.f2_ppm), d.z,
                                     bounds_error=False, fill_value=0.0)
    G1, G2 = np.meshgrid(kernel.f1_ppm, kernel.f2_ppm, indexing="ij")
    z_exp = interp(np.stack([G1.ravel(), G2.ravel()], axis=-1)).reshape(
        kernel.shape)
    return {"kind": "2d", "r": rec, "kernel": kernel, "z_exp": z_exp,
            "exp": z_exp.ravel(), "w": 1.0 / (np.abs(z_exp).max() or 1.0)}


def fit_cofit(entries: list[tuple], share: tuple[str, ...] = DEFAULT_SHARE,
              windows: list | None = None, method: str = "3QMAS",
              ) -> MultiFitResult:
    """Co-fit a mix of 1D and 2D (MQMAS) datasets with shared physical model.

    Each entry is ``(recipe, spec)`` where ``spec`` is a ``(ppm, amp)`` pair for
    a 1D dataset or a :class:`larmor.twod.Data2D` for a 2D one. Site j of every
    recipe is the same physical site; parameters in ``share`` are tied to
    dataset 0, so e.g. an MQMAS 2D and a 1D MAS spectrum of the same sample are
    fit together with a common Cq/eta/delta while amplitudes stay independent.
    """
    from larmor import twod

    if len(entries) < 2:
        raise ValueError("co-fit needs at least two datasets")
    recipes = [rec for rec, _ in entries]
    n_sites = len(recipes[0].sites)
    for r in recipes[1:]:
        if len(r.sites) != n_sites:
            raise ValueError("all recipes must have the same number of sites")
        for j in range(n_sites):
            if r.sites[j].model != recipes[0].sites[j].model:
                raise ValueError(
                    f"site {j} model differs between datasets "
                    f"({r.sites[j].model} vs {recipes[0].sites[j].model})")

    prep = []
    for k, (rec, spec) in enumerate(entries):
        win = windows[k] if windows and k < len(windows) else None
        if hasattr(spec, "z"):                       # a Data2D
            prep.append(_prepare_2d(rec, spec, method))
        else:
            ppm, amp = spec
            prep.append(_prepare_1d(rec, np.asarray(ppm, float),
                                    np.asarray(amp, float), win))

    # ---------------- parameters
    params = lmfit.Parameters()
    for k, e in enumerate(prep):
        for i, site in enumerate(e["r"].sites):
            for pname, p in site.params.items():
                params.add(_pname(k, i, site, pname), value=p.value, vary=p.vary,
                           min=p.min if p.min is not None else -np.inf,
                           max=p.max if p.max is not None else np.inf)
    for k, e in enumerate(prep):
        for i, site in enumerate(e["r"].sites):
            for pname, p in site.params.items():
                if p.expr:
                    trans = _EXPR_REF.sub(
                        lambda m: _pname(k, int(m.group(1)),
                                         e["r"].sites[int(m.group(1))],
                                         m.group(2)), p.expr)
                    params[_pname(k, i, site, pname)].expr = trans
    for k, e in enumerate(prep):
        if k == 0:
            continue
        for i, site in enumerate(e["r"].sites):
            for pname in site.params:
                if pname in share and not site.params[pname].expr:
                    params[_pname(k, i, site, pname)].expr = \
                        _pname(0, i, recipes[0].sites[i], pname)
    try:
        params.update_constraints()
    except Exception as exc:
        raise ConstraintError(f"invalid co-fit constraints: {exc}") from exc

    def apply_params(p):
        for k, e in enumerate(prep):
            for i, site in enumerate(e["r"].sites):
                for pname in site.params:
                    lp = p[_pname(k, i, site, pname)]
                    site.params[pname].value = float(lp.value)
                    site.params[pname].stderr = (float(lp.stderr)
                                                 if lp.stderr is not None else None)

    def model_vec(e):
        r = e["r"]
        if e["kind"] == "1d":
            total = np.sum([simulate_site(s, e["ctx"]) for s in r.sites], axis=0)
            return np.interp(e["ppm"][e["sel"]], e["ctx"].x_ppm, total)
        total, _ = twod.simulate_2d(r, e["kernel"])
        return total.ravel()

    def residual(p):
        apply_params(p)
        return np.concatenate([(model_vec(e) - e["exp"]) * e["w"] for e in prep])

    result = lmfit.minimize(residual, params, method="least_squares")
    if not result.errorbars:
        retry = lmfit.minimize(residual, result.params.copy(), method="leastsq")
        if retry.errorbars:
            result = retry
    apply_params(result.params)

    rmsds, per_dataset = [], []
    for e in prep:
        r = e["r"]
        yi = model_vec(e)
        yw = e["exp"]
        rmsds.append(float(np.sqrt(np.mean((yi - yw) ** 2)) / (np.abs(yw).max() or 1)))
        r.fit_rmsd = rmsds[-1]
        if e["kind"] == "1d":
            total = np.sum([simulate_site(s, e["ctx"]) for s in r.sites], axis=0)
            per_dataset.append({"kind": "1d", "x": e["ctx"].x_ppm,
                                "y_fit": total})
        else:
            total, per = twod.simulate_2d(r, e["kernel"])
            per_dataset.append({"kind": "2d", "f2": e["kernel"].f2_ppm,
                                "f1": e["kernel"].f1_ppm, "z_fit": total,
                                "per_site": per})
    return MultiFitResult(recipes=recipes, lmfit_result=result, rmsd=rmsds,
                          per_dataset=per_dataset)
