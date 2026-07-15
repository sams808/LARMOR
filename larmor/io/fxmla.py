"""Reader for dmfit .fxmla fit files, and conversion to LARMOR recipes.

A .fxmla file is XML holding:
  - <FitParameters>: per-dimension list of <line> elements, each one lineshape
    model (CzSimple = simple Czjzek, "Gaus/Lor", "ss band", ...) with parameters
    carrying optional Unit= and Fix= attributes,
  - <ExpData>: the experimental spectrum inlined as a dmfit "SIMP" ASCII block.

Czjzek width convention (established empirically in Phase 0 against
CaAlGlass.fxmla, RMSD minimum sharply at exactly 1/2):

    mrsimulator CzjzekDistribution sigma [MHz] = dmfit sCZ_CQ [kHz] / 2000

dmfit also stores CQ = 2 * sCZ_CQ (the mode of the |Cq| distribution) and
CQ_max as a derived duplicate of CQ.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from larmor.recipe import Param, Recipe, SiteModel, sha256_of


@dataclass
class DmfitParam:
    value: float
    unit: str | None = None
    fix_flag: bool = False  # dmfit's Fix="*" attribute, stored verbatim


@dataclass
class DmfitLine:
    model_name: str
    model_nb: int | None
    name: str
    params: dict[str, DmfitParam] = field(default_factory=dict)


@dataclass
class DmfitDimension:
    label: str               # "F2", "F1"
    nucleus: str
    frequency_MHz: float
    spin_rate: float | None
    lines: list[DmfitLine] = field(default_factory=list)


@dataclass
class DmfitSpectrum:
    header: dict[str, float]
    amplitude: np.ndarray        # real part; imaginary column is stored separately
    imaginary: np.ndarray

    @property
    def ppm(self) -> np.ndarray:
        """Frequency axis in ppm, dmfit SIMP convention: (X0 + i*dX - Sr)/Sf."""
        h = self.header
        n = int(h["NP"])
        if "dX" in h:
            freq = h["X0"] + np.arange(n) * h["dX"]
        else:  # older SIMP blocks: X0 is the left edge, SW spans the axis
            freq = h["X0"] - np.arange(n) * (h["SW"] / n)
        return (freq - h.get("Sr", 0.0)) / h["Sf"]


@dataclass
class DmfitFile:
    path: str
    version: str
    fit_mode: str                    # e.g. "Fit 1D", "MQMAS"
    dimensions: list[DmfitDimension]
    spectrum: DmfitSpectrum | None
    comment: str = ""

    @property
    def is_2d(self) -> bool:
        return "NI" in (self.spectrum.header if self.spectrum else {})


# --------------------------------------------------------------------------
def read(path: str | Path) -> DmfitFile:
    """Parse a dmfit .fxmla file. Read-only; never modifies the source."""
    path = Path(path)
    root = ET.fromstring(path.read_text(encoding="utf-8", errors="replace"))

    fitparams = root.find("FitParameters")
    version = fitparams.findtext("DMFitVersion", default="")
    fit_mode = fitparams.findtext("FitModeAsc", default="")

    dimensions = []
    for dim_el in fitparams.findall("Dimension"):
        label = (dim_el.text or "").strip()
        spin_rate_txt = dim_el.findtext("spinrate")
        dim = DmfitDimension(
            label=label,
            nucleus=dim_el.findtext("nucleus", default=""),
            frequency_MHz=float(dim_el.findtext("frequency", default="0")),
            spin_rate=float(spin_rate_txt) if spin_rate_txt else None,
        )
        for line_el in dim_el.findall("line"):
            model_nb_txt = line_el.findtext("ModelNb")
            line = DmfitLine(
                model_name=line_el.findtext("ModelName", default=""),
                model_nb=int(model_nb_txt) if model_nb_txt else None,
                name=line_el.findtext("Name", default="") or "",
            )
            # every leaf element with a float-parsable text is a parameter
            for group in line_el:
                if group.tag in ("ModelName", "ModelNb", "Name"):
                    continue
                for p in group:
                    try:
                        value = float((p.text or "").strip())
                    except ValueError:
                        continue
                    line.params[p.tag] = DmfitParam(
                        value=value,
                        unit=p.get("Unit"),
                        fix_flag=p.get("Fix") == "*",
                    )
            dim.lines.append(line)
        dimensions.append(dim)

    spectrum = _parse_simp_block(root)
    comment = spectrum.header.pop("_comment", "") if spectrum else ""

    return DmfitFile(
        path=str(path), version=version, fit_mode=fit_mode,
        dimensions=dimensions, spectrum=spectrum, comment=str(comment),
    )


def _parse_simp_block(root: ET.Element) -> DmfitSpectrum | None:
    data_el = root.find("ExpData/Data")
    if data_el is None or not data_el.text:
        return None
    header: dict[str, float] = {}
    rows: list[tuple[float, float]] = []
    in_data = False
    for raw in data_el.text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line == "DATA":
            in_data = True
            continue
        if not in_data:
            if "=" in line:
                key, _, val = line.partition("=")
                key, val = key.strip(), val.strip()
                if key.lower() == "comment":
                    header["_comment"] = val  # type: ignore[assignment]
                else:
                    try:
                        header[key] = float(val)
                    except ValueError:
                        pass
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    arr = np.array(rows)
    return DmfitSpectrum(header=header, amplitude=arr[:, 0], imaginary=arr[:, 1])


# --------------------------------------------------------------------------
# dmfit -> LARMOR recipe conversion

#: The Phase 0 empirical convention: mrsimulator sigma = dmfit sCZ_CQ / 2.
SCZ_TO_SIGMA = 0.5


def to_recipe(dm: DmfitFile, dimension: int = 0) -> tuple[Recipe, list[str]]:
    """Convert one dimension of a parsed dmfit file to a LARMOR recipe.

    Returns (recipe, warnings). Lines with models LARMOR does not yet fit
    (spinning sidebands, exotic models) are reported in warnings, not dropped
    silently.
    """
    dim = dm.dimensions[dimension]
    warnings: list[str] = []
    recipe = Recipe(
        sample=dm.comment,
        source_kind="fxmla",
        source_path=dm.path,
        source_sha256=sha256_of(dm.path),
        nucleus=dim.nucleus,
        larmor_frequency_MHz=dim.frequency_MHz,
        spin_rate_Hz=_spin_rate_hz(dim.spin_rate),
    )
    if dm.fit_mode not in ("Fit 1D", ""):
        warnings.append(
            f"fit mode {dm.fit_mode!r}: only the 1D lineshape content is converted; "
            "2D methods (MQMAS, ...) are a Phase 2 feature"
        )

    for i, line in enumerate(dim.lines):
        if line.model_name == "CzSimple":
            scz_khz = line.params["sCZ_CQ"].value
            site = SiteModel(
                model="czjzek",
                label=line.name or f"CzSimple-{i}",
                params={
                    "isotropic_chemical_shift_ppm": Param(line.params["pos"].value),
                    "sigma_Cq_MHz": Param(scz_khz * SCZ_TO_SIGMA / 1000.0, min=0.05),
                    "shift_fwhm_ppm": Param(
                        abs(line.params["dCS"].value) if "dCS" in line.params else 10.0,
                        min=0.1,
                    ),
                    "amplitude": Param(line.params["amp"].value, min=0.0),
                },
            )
            recipe.sites.append(site)
        elif line.model_name == "Gaus/Lor":
            site = SiteModel(
                model="gauss_lor",
                label=line.name or f"GausLor-{i}",
                params={
                    "isotropic_chemical_shift_ppm": Param(line.params["pos"].value),
                    "shift_fwhm_ppm": Param(abs(line.params["wid"].value), min=0.1),
                    "amplitude": Param(abs(line.params["amp"].value), min=0.0),
                    "gl": Param(line.params.get("gl", DmfitParam(1.0)).value,
                                vary=False, min=0.0, max=1.0),
                },
            )
            recipe.sites.append(site)
        elif line.model_name == "ss band":
            warnings.append(
                f"line {i} ('ss band' at {line.params.get('pos', DmfitParam(0)).value:.1f} ppm) "
                "skipped: explicit sideband lines are handled by the simulation itself in LARMOR"
            )
        else:
            warnings.append(f"line {i} (model {line.model_name!r}) not yet supported, skipped")

    recipe.notes.append(
        f"imported from dmfit {dm.version} ({dm.fit_mode}); "
        f"Czjzek sigma = sCZ_CQ x {SCZ_TO_SIGMA} (Phase 0 convention)"
    )
    recipe.notes.extend(warnings)
    return recipe, warnings


def _spin_rate_hz(value: float | None) -> float:
    """dmfit stores spinrate sometimes in Hz, sometimes in kHz. Disambiguate."""
    if value is None:
        return 0.0
    return value * 1000.0 if value < 200 else value
