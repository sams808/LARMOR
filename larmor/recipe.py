"""The LARMOR recipe: a small, diffable JSON description of a fit.

This is the file format both Guided and Expert modes share. It references the
source data by path + SHA-256 rather than embedding it (the core fix over
dmfit's .fxmla, which inlines the full spectrum into every saved fit).
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path

RECIPE_VERSION = 1


@dataclass
class Param:
    """One fittable quantity: value, uncertainty, and fit behavior."""

    value: float
    stderr: float | None = None
    vary: bool = True
    min: float | None = None
    max: float | None = None


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
    fit_window_ppm: tuple[float, float] | None = None
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
            params = {k: Param(**p) for k, p in s.pop("params", {}).items()}
            sites.append(SiteModel(params=params, **s))
        window = d.pop("fit_window_ppm", None)
        recipe = cls(sites=sites, **d)
        recipe.fit_window_ppm = tuple(window) if window else None
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
