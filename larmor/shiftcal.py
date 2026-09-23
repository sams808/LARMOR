"""Shielding-to-shift calibration: delta = a*sigma + b over reference compounds.

A GIPAW shielding sigma is an absolute number on a scale set by the code,
functional, cutoffs and pseudopotential of the calculation; a measured
chemical shift delta is relative to a reference compound. The two are
related linearly, but neither the slope nor the intercept is universal: in
practice a != -1 (typically -0.7 to -1.1 for 19F, 29Si, 31P), so a single
sigma_ref (delta = sigma_ref - sigma, slope -1) is a special case, not the
rule. The line has to be FITTED, per DFT setup, over reference compounds
whose shifts were measured under the same referencing as the sample.

What this module does that a two-point OLS does not:
  * weights each reference by its measured shift error and propagates the
    covariance of (a, b) into every predicted delta (delta_err);
  * says when the line has no residual degrees of freedom -- with two
    references the line passes through both points BY CONSTRUCTION, the
    residuals are zero and the +- reflects the reference-shift errors only,
    never the quality of the calibration;
  * refuses references computed with different [calculation] settings
    (code / version / functional / cutoffs / the observed element's
    pseudopotential), unless the mismatch is explicitly overridden, in which
    case it is recorded;
  * serialises to <name>.shiftcal.json and writes the provenance record and
    the Methods clause a paper needs.

The weighted line is qcpmg_fields.weighted_line -- the same closed form the
QCPMG two-field extrapolation uses. Qt-free.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from larmor import dft
from larmor.qcpmg_fields import weighted_line

__all__ = [
    "SHIFTCAL_VERSION", "FLAG_ZERO_DOF", "FLAG_UNWEIGHTED", "FLAG_NO_COV",
    "CalibrationPoint", "ShiftCalibration", "fit_calibration",
    "from_sigma_ref", "compatible", "mismatches", "point_from_magres",
]

SHIFTCAL_VERSION = 1

FLAG_ZERO_DOF = ("zero residual degrees of freedom: the line is exact by "
                 "construction; the ± reflects the reference-shift errors only")
FLAG_UNWEIGHTED = "unweighted: some references carry no uncertainty"
FLAG_NO_COV = ("no covariance: an unweighted two-point line carries no "
               "uncertainty estimate")

#: the fingerprint keys compatible() compares, in reporting order
_FINGERPRINT_KEYS = ("calc_code", "calc_code_version", "calc_xcfunctional",
                     "calc_cutoffenergy", "calc_cutoffenergy_rho")


@dataclass
class CalibrationPoint:
    """One reference compound: its computed isotropic shielding (the mean
    over the equivalent atoms) and its measured shift."""
    name: str
    sigma_iso: float
    delta_exp: float | None = None
    delta_err: float | None = None
    source: str = ""                  # 'recipe: CaF2.recipe.json' / 'typed'
    magres: str = ""
    n_atoms: int = 1
    calc: dft.MagresCalc | None = None

    @property
    def ready(self) -> bool:
        return (self.delta_exp is not None
                and math.isfinite(float(self.delta_exp)))

    def to_dict(self) -> dict:
        return {"name": self.name, "sigma_iso": float(self.sigma_iso),
                "delta_exp": None if self.delta_exp is None
                else float(self.delta_exp),
                "delta_err": None if self.delta_err is None
                else float(self.delta_err),
                "source": self.source, "magres": self.magres,
                "n_atoms": int(self.n_atoms),
                "calc": self.calc.to_dict() if self.calc is not None else None}

    @classmethod
    def from_dict(cls, d: dict) -> "CalibrationPoint":
        calc = d.get("calc")
        return cls(name=str(d.get("name", "")),
                   sigma_iso=float(d.get("sigma_iso", 0.0)),
                   delta_exp=dft._float_or_none(d.get("delta_exp")),
                   delta_err=dft._float_or_none(d.get("delta_err")),
                   source=str(d.get("source", "") or ""),
                   magres=str(d.get("magres", "") or ""),
                   n_atoms=int(d.get("n_atoms", 1) or 1),
                   calc=dft.MagresCalc.from_dict(calc) if calc else None)


@dataclass
class ShiftCalibration:
    """The line delta = slope*sigma + intercept with its statistics.

    ``kind`` is 'fit' (over reference compounds) or 'sigma_ref' (slope -1,
    intercept sigma_ref, no statistics)."""
    isotope: str
    slope: float
    intercept: float
    cov: list[list[float]] | None = None      # [[var_a, cov_ab], [cov_ab, var_b]]
    n: int = 0
    dof: int = 0
    rmse_ppm: float = 0.0
    residuals: list[float] = field(default_factory=list)
    points: list[CalibrationPoint] = field(default_factory=list)
    calc: dft.MagresCalc | None = None
    kind: str = "fit"
    flags: list[str] = field(default_factory=list)
    chi2_red: float | None = None
    file: str = ""

    # ---- use --------------------------------------------------------
    def delta(self, sigma: float) -> float:
        return float(self.slope * float(sigma) + self.intercept)

    def delta_err(self, sigma: float) -> float | None:
        """1σ uncertainty of delta(sigma) from cov(a, b); None without a
        covariance."""
        if self.cov is None:
            return None
        c = np.asarray(self.cov, float)
        if c.shape != (2, 2) or not np.all(np.isfinite(c)):
            return None
        v = np.array([float(sigma), 1.0])
        var = float(v @ c @ v)
        return math.sqrt(var) if var > 0 else 0.0

    @property
    def zeta_scale(self) -> float:
        """Factor applied to a computed shielding anisotropy zeta_sigma to
        seed LARMOR's csa_mas ``zeta_ppm`` (the SHIELDING anisotropy passed
        to mrsimulator, delta_aniso = -zeta): the shift tensor is
        a*sigma + b, its anisotropy a*zeta_sigma, and the shielding-
        convention zeta of that shift tensor is -a*zeta_sigma. For
        sigma_ref (a = -1) this is +1, today's behaviour."""
        return -float(self.slope)

    @property
    def slope_err(self) -> float | None:
        return None if self.cov is None else math.sqrt(max(self.cov[0][0], 0.0))

    @property
    def intercept_err(self) -> float | None:
        return None if self.cov is None else math.sqrt(max(self.cov[1][1], 0.0))

    # ---- text -------------------------------------------------------
    def describe(self) -> str:
        if self.kind == "sigma_ref":
            return f"δ = σ_ref − σ   σ_ref = {self.intercept:.1f} ppm"
        sa, sb = self.slope_err, self.intercept_err
        a_txt = (f"a = {self.slope:.4f}" if sa is None
                 else f"a = {self.slope:.4f} ± {sa:.4f}")
        b_txt = (f"b = {self.intercept:.2f} ppm" if sb is None
                 else f"b = {self.intercept:.2f} ± {sb:.2f} ppm")
        cov_txt = ("" if self.cov is None
                   else f"   cov(a,b) = {self.cov[0][1]:.3g}")
        return (f"δ = a·σ + b   {a_txt}   {b_txt}{cov_txt}   "
                f"n = {self.n}, dof = {self.dof}, RMSE {self.rmse_ppm:.2f} ppm")

    def methods_clause(self) -> str:
        """The conversion half of the Methods sentence (no leading space)."""
        if self.kind == "sigma_ref":
            return (f"shieldings were converted to chemical shifts with "
                    f"σ_ref = {self.intercept:.1f} ppm (δ = σ_ref − σ).")
        names = ", ".join(p.name for p in self.points if p.name) or "—"
        sa, sb = self.slope_err, self.intercept_err
        a_txt = (f"a = {self.slope:.3f}" if sa is None
                 else f"a = {self.slope:.3f} ± {sa:.3f}")
        b_txt = (f"b = {self.intercept:.1f}" if sb is None
                 else f"b = {self.intercept:.1f} ± {sb:.1f}")
        n = self.n or len(self.points)
        if self.dof <= 0:
            tail = ("; two references: the line is exact by construction, "
                    "no residual degrees of freedom)")
        else:
            tail = f", RMSE {self.rmse_ppm:.1f} ppm)"
        return (f"shieldings were converted to chemical shifts through "
                f"δ = a·σ + b fitted over {n} reference compound"
                f"{'s' if n != 1 else ''} ({names}; {a_txt}, {b_txt} ppm"
                + tail)

    # ---- persistence -----------------------------------------------
    def to_dict(self) -> dict:
        return {
            "larmor_shiftcal_version": SHIFTCAL_VERSION,
            "isotope": self.isotope, "kind": self.kind,
            "slope": float(self.slope), "intercept": float(self.intercept),
            "cov": None if self.cov is None
            else [[float(v) for v in row] for row in self.cov],
            "n": int(self.n), "dof": int(self.dof),
            "rmse_ppm": float(self.rmse_ppm),
            "chi2_red": None if self.chi2_red is None else float(self.chi2_red),
            "residuals": [float(r) for r in self.residuals],
            "flags": list(self.flags),
            "points": [p.to_dict() for p in self.points],
            "calc": self.calc.to_dict() if self.calc is not None else None,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ShiftCalibration":
        d = dict(d or {})
        cov = d.get("cov")
        pts = []
        for p in d.get("points") or []:
            if isinstance(p, dict):
                pts.append(CalibrationPoint.from_dict(p))
            else:                       # a bare name (hand-written record)
                pts.append(CalibrationPoint(name=str(p), sigma_iso=0.0))
        calc = d.get("calc")
        return cls(
            isotope=str(d.get("isotope", "") or ""),
            slope=float(d.get("slope", -1.0)),
            intercept=float(d.get("intercept", 0.0)),
            cov=None if cov is None else [[float(v) for v in row] for row in cov],
            n=int(d.get("n", len(pts)) or 0), dof=int(d.get("dof", 0) or 0),
            rmse_ppm=float(d.get("rmse_ppm", 0.0) or 0.0),
            residuals=[float(r) for r in (d.get("residuals") or [])],
            points=pts,
            calc=dft.MagresCalc.from_dict(calc) if calc else None,
            kind=str(d.get("kind", "fit") or "fit"),
            flags=list(d.get("flags") or []),
            chi2_red=dft._float_or_none(d.get("chi2_red")),
            file=str(d.get("file", "") or ""))

    def save(self, path: str | Path) -> None:
        p = Path(path)
        p.write_text(json.dumps(self.to_dict(), indent=2, ensure_ascii=False),
                     encoding="utf-8")
        self.file = str(p)

    @classmethod
    def load(cls, path: str | Path) -> "ShiftCalibration":
        p = Path(path)
        d = json.loads(p.read_text(encoding="utf-8"))
        if "larmor_shiftcal_version" not in d or "slope" not in d:
            raise ValueError(f"{p.name} is not a LARMOR shift calibration "
                             "(.shiftcal.json)")
        cal = cls.from_dict(d)
        cal.file = str(p)
        return cal


# ------------------------------------------------------------- building
def from_sigma_ref(sigma_ref: float, isotope: str = "") -> ShiftCalibration:
    """The one-number conversion delta = sigma_ref - sigma (slope -1)."""
    return ShiftCalibration(isotope=isotope, slope=-1.0,
                            intercept=float(sigma_ref), kind="sigma_ref")


def _same(a, b) -> bool:
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    if isinstance(a, (int, float)) or isinstance(b, (int, float)):
        try:
            return abs(float(a) - float(b)) <= 1e-6
        except (TypeError, ValueError):
            return False
    return str(a).strip() == str(b).strip()


def compatible(a: dft.MagresCalc | None, b: dft.MagresCalc | None,
               element: str) -> list[str]:
    """The fingerprint keys on which two [calculation] headers differ
    (empty = compatible). Either side None (no header recorded, or a
    sigma_ref line) compares as compatible: there is nothing to check."""
    if a is None or b is None:
        return []
    fa, fb = a.fingerprint(element), b.fingerprint(element)
    return [k for k in fa if not _same(fa[k], fb[k])]


def mismatches(points: list[CalibrationPoint], element: str) -> list[str]:
    """'CaF2: calc_cutoffenergy' for every point whose header differs from
    the first point's header (points without a header are skipped)."""
    ref = next((p.calc for p in points if p.calc is not None), None)
    out = []
    if ref is None:
        return out
    for p in points:
        if p.calc is None:
            continue
        diff = compatible(ref, p.calc, element)
        if diff:
            out.append(f"{p.name}: " + ", ".join(diff))
    return out


def point_from_magres(mf: dft.MagresFile, isotope: str, name: str = "",
                      delta_exp: float | None = None,
                      delta_err: float | None = None, source: str = "",
                      ) -> CalibrationPoint:
    """One reference point from a magres file: the isotropic shielding
    averaged over every atom of ``isotope`` (a compound with more than one
    inequivalent site of that isotope gets a note in ``source`` -- give the
    sites separately through the dialog when their shifts are resolved)."""
    if not any(s.isotope for s in mf.sites):
        dft.assign_isotopes(mf.sites)
    sites = dft.sites_for_isotope(mf.sites, isotope)
    if not sites:
        raise ValueError(f"{mf.name} has no {isotope} site")
    groups = dft.group_equivalent(sites)
    vals = [s.shielding()["iso_ppm"] for s in sites if s.shielding()]
    if not vals:
        raise ValueError(f"{mf.name}: no shielding tensor for {isotope}")
    note = source
    if len(groups) > 1:
        note = (source + ("; " if source else "")
                + f"mean over {len(groups)} inequivalent sites")
    return CalibrationPoint(
        name=name or mf.calc.prefix or Path(mf.path).stem.split(".")[0],
        sigma_iso=float(np.mean(vals)), delta_exp=delta_exp,
        delta_err=delta_err, source=note, magres=mf.path,
        n_atoms=len(sites), calc=mf.calc)


def fit_calibration(points: list[CalibrationPoint], isotope: str, *,
                    allow_mismatch: bool = False) -> ShiftCalibration:
    """Weighted least-squares line over the READY points (those with a
    measured shift).

    Weights 1/err² when every point has an error, else unweighted and
    flagged. The covariance is the propagated one from the shift errors,
    inflated by the reduced χ² when dof > 0 and χ²_red > 1 (never deflated:
    over-reported reference errors do not make the line better). With two
    points the line is exact -- the propagated covariance is kept and the
    zero-dof flag set; an unweighted two-point line has no covariance at
    all. References whose [calculation] settings differ raise ValueError
    unless ``allow_mismatch`` (then the mismatch is a recorded flag).
    """
    pts = [p for p in points if p.ready]
    if len(pts) < 2:
        raise ValueError("need at least two reference compounds with a "
                         "measured shift")
    element = dft.element_of(isotope)
    mism = mismatches(pts, element)
    if mism and not allow_mismatch:
        raise ValueError("references computed with different [calculation] "
                         "settings: " + "; ".join(mism))
    x = np.array([float(p.sigma_iso) for p in pts])
    y = np.array([float(p.delta_exp) for p in pts])
    weighted = all(p.delta_err is not None and float(p.delta_err) > 0
                   for p in pts)
    err = (np.array([float(p.delta_err) for p in pts]) if weighted
           else np.ones(len(pts)))
    try:
        slope, intercept, cov = weighted_line(x, y, err)
    except ValueError:
        raise ValueError("the reference shieldings are identical: no line "
                         "can be fitted") from None
    res = y - (slope * x + intercept)
    n = len(pts)
    dof = n - 2
    rmse = float(np.sqrt(np.mean(res ** 2)))
    flags: list[str] = []
    chi2_red: float | None = None
    if mism:
        flags.append("mixed [calculation] settings (override recorded): "
                     + "; ".join(mism))
    cov_out: list[list[float]] | None
    if weighted:
        if dof > 0:
            chi2_red = float(np.sum((res / err) ** 2)) / dof
            if chi2_red > 1.0:
                cov = cov * chi2_red
        else:
            flags.append(FLAG_ZERO_DOF)
        cov_out = cov.tolist()
    else:
        flags.append(FLAG_UNWEIGHTED)
        if dof > 0:
            # unit weights are fictitious: the residual variance IS the scale
            chi2_red = float(np.sum(res ** 2)) / dof
            cov_out = (cov * chi2_red).tolist()
        else:
            flags.append(FLAG_ZERO_DOF)
            flags.append(FLAG_NO_COV)
            cov_out = None
    calc = next((p.calc for p in pts if p.calc is not None), None)
    return ShiftCalibration(
        isotope=isotope, slope=slope, intercept=intercept, cov=cov_out,
        n=n, dof=dof, rmse_ppm=rmse, residuals=[float(r) for r in res],
        points=list(pts), calc=calc, kind="fit", flags=flags,
        chi2_red=chi2_red)
