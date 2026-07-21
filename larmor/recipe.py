"""The LARMOR recipe: a small, diffable JSON description of a fit.

This is the file format both Guided and Expert modes share. It references the
source data by path + SHA-256 rather than embedding it (the core fix over
dmfit's .fxmla, which inlines the full spectrum into every saved fit).
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path

RECIPE_VERSION = 1


@dataclass
class Param:
    """One fittable quantity: value, uncertainty, and fit behavior.

    Constraints (ssNake-inspired, but algebraic):
      - vary=False           freezes the parameter at `value`
      - min / max            box bounds enforced during the fit
      - expr                 links this parameter to others by an algebraic
                             expression, e.g. "0.5 * s0.amplitude" or
                             "s0.shift_fwhm_ppm" (shared linewidth) or
                             "s0.isotropic_chemical_shift_ppm + 30".
                             Site parameters are addressed as s<index>.<name>.
                             A linked parameter is not varied independently;
                             its value and stderr are derived (with full error
                             propagation via lmfit).
    """

    value: float
    stderr: float | None = None
    vary: bool = True
    min: float | None = None
    max: float | None = None
    expr: str | None = None


@dataclass
class SiteModel:
    """One spectral component.

    model:
      - "czjzek": Czjzek distribution of quadrupolar tensors
          params: isotropic_chemical_shift_ppm, sigma_Cq_MHz, shift_fwhm_ppm, amplitude
      - "gauss_lor": analytic pseudo-Voigt, y = gl*Gaussian + (1-gl)*Lorentzian
          params: isotropic_chemical_shift_ppm, shift_fwhm_ppm, amplitude, gl
    """

    model: str
    label: str = ""
    params: dict[str, Param] = field(default_factory=dict)
    #: for the "spectrum" (background) model: the reference trace it renders,
    #: {"ppm": [...], "amp": [...]} (amp pre-normalized to unit peak)
    ref: dict | None = None


@dataclass
class Recipe:
    sample: str = ""
    source_kind: str = ""        # "fxmla" | "bruker"
    source_path: str = ""
    source_sha256: str = ""
    nucleus: str = ""
    larmor_frequency_MHz: float = 0.0
    spin_rate_Hz: float = 0.0
    engine: str = "czjzek-kernel+lmfit"
    sites: list[SiteModel] = field(default_factory=list)
    #: ordered processing pipeline ({"op": name, ...}) applied to the source
    #: data before fitting; replayed on load so a recipe carries its
    #: processing, making the fit reproducible end to end
    processing: list = field(default_factory=list)
    #: True when `processing` starts from the raw fid rather than pdata
    processing_from_raw: bool = False
    fit_window_ppm: tuple[float, float] | None = None
    #: optional list of [hi_ppm, lo_ppm] regions (dmfit's "Zones"); when set,
    #: the fit residual is evaluated only inside their union
    fit_zones: list = field(default_factory=list)
    fit_rmsd: float | None = None
    notes: list[str] = field(default_factory=list)

    # ---------- serialization ----------
    def to_dict(self) -> dict:
        d = asdict(self)
        d["larmor_recipe_version"] = RECIPE_VERSION
        return d

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def from_dict(cls, d: dict) -> "Recipe":
        d = dict(d)
        d.pop("larmor_recipe_version", None)
        sites = []
        for s in d.pop("sites", []):
            s = dict(s)   # never mutate the caller's dicts
            params = {k: Param(**p) for k, p in s.pop("params", {}).items()}
            sites.append(SiteModel(params=params, **s))
        window = d.pop("fit_window_ppm", None)
        # forward compatibility: a recipe written by a NEWER LARMOR may carry
        # fields this version does not know. Drop them with a note instead of
        # refusing to open the file.
        known = {f.name for f in fields(cls)} - {"sites", "fit_window_ppm"}
        unknown = [k for k in d if k not in known]
        extra_note = None
        if unknown:
            extra_note = ("ignored unknown recipe fields (written by a newer "
                          "LARMOR?): " + ", ".join(sorted(unknown)))
            for k in unknown:
                d.pop(k)
        recipe = cls(sites=sites, **d)
        recipe.fit_window_ppm = tuple(window) if window else None
        if extra_note:
            recipe.notes.append(extra_note)
        return recipe

    @classmethod
    def load(cls, path: str | Path) -> "Recipe":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def sha256_of(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
