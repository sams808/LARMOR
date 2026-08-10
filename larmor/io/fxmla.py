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
    """Parse a dmfit .fxmla file. Read-only; never modifies the source.

    Handles both generations of the format: the classic all-XML file (SIMP
    ASCII spectrum in ``<ExpData><Data>``), and the newer dmfit export
    (seen from dmfit #20230120) whose ``<ExpData><Data>`` holds a **CSDM
    JSON** block instead — base64 float32 spectrum plus the full Bruker
    metadata verbatim. That metadata (audit trails, acqus) is full of raw
    unescaped ``<...>`` tokens, so the newer files are NOT well-formed XML
    and a strict parse raises ParseError; ``_read_lenient`` then extracts
    the (well-formed) ``<FitParameters>`` island and the JSON payload
    separately, so a fit dmfit itself wrote today still imports."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return _read_lenient(path, text)

    fitparams = root.find("FitParameters")
    version, fit_mode, dimensions = _parse_fitparameters(fitparams)

    spectrum = _parse_simp_block(root)
    comment = spectrum.header.pop("_comment", "") if spectrum else ""

    return DmfitFile(
        path=str(path), version=version, fit_mode=fit_mode,
        dimensions=dimensions, spectrum=spectrum, comment=str(comment),
    )


def _parse_fitparameters(fitparams: ET.Element):
    """(version, fit_mode, dimensions) from a parsed <FitParameters> element —
    shared by the strict and lenient read paths."""
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
            # every leaf element with a float-parsable text is a parameter.
            # dmfit .fxmla nests params one level under groups; the flatter
            # .fxml (fit-parameters-only export) lists them directly under
            # <line> -- handle both by treating any leaf (childless) element as
            # a parameter and descending into any element that has children.
            def _collect(el):
                if el.tag in ("ModelName", "ModelNb", "Name"):
                    return
                kids = list(el)
                if kids:
                    for kid in kids:
                        _collect(kid)
                    return
                try:
                    value = float((el.text or "").strip())
                except (ValueError, TypeError):
                    return
                # .fxmla marks fitted params Fix="*"; the flat .fxml uses Fit="*"
                line.params[el.tag] = DmfitParam(
                    value=value, unit=el.get("Unit"),
                    fix_flag=el.get("Fix") == "*" or el.get("Fit") == "*",
                )

            for group in line_el:
                _collect(group)
            dim.lines.append(line)
        dimensions.append(dim)

    return version, fit_mode, dimensions


def _read_lenient(path: Path, text: str) -> DmfitFile:
    """Fallback for a dmfit file that is not well-formed XML as a whole (the
    newer CSDM-JSON ``<ExpData>`` embeds raw Bruker audit text with bare
    ``<``/``>``): parse the two islands separately — ``<FitParameters>`` is
    always clean XML, and the ``<Data>`` payload is valid JSON (JSON strings
    are perfectly happy holding '<'; only the XML wrapper ever chokes)."""
    import json as _json
    import re as _re

    m = _re.search(r"<FitParameters>.*?</FitParameters>", text, _re.DOTALL)
    if not m:
        raise ValueError(f"{path.name}: no <FitParameters> block found")
    version, fit_mode, dimensions = _parse_fitparameters(ET.fromstring(m.group(0)))

    spectrum, comment = None, ""
    md = _re.search(r"<Data>\s*(\{.*\})\s*</Data>", text, _re.DOTALL)
    if md:
        try:
            spectrum, comment = _spectrum_from_csdm(_json.loads(md.group(1)))
        except Exception:
            spectrum = None                # fit params still usable without it

    return DmfitFile(
        path=str(path), version=version, fit_mode=fit_mode,
        dimensions=dimensions, spectrum=spectrum, comment=comment,
    )


def _spectrum_from_csdm(payload: dict) -> tuple[DmfitSpectrum, str]:
    """A DmfitSpectrum from dmfit's CSDM-JSON <Data> payload (1D linear
    frequency dimension, base64 float32 components). The ppm axis maps onto
    the existing SIMP header convention ((X0 + i*dX - Sr)/Sf) by storing
    X0 = coordinates_offset [Hz], dX = increment [Hz], Sf = origin_offset
    [MHz], Sr = 0 — Hz/MHz = ppm relative to the origin offset, exactly how
    TopSpin's own axis is defined."""
    import base64

    cs = payload["csdm"]
    dim = cs["dimensions"][0]
    dv = cs["dependent_variables"][0]

    def _qty(s: str, unit: str) -> float:
        v, u = s.split()
        if u.lower() != unit.lower():
            raise ValueError(f"unexpected unit {u!r} (wanted {unit})")
        return float(v)

    n = int(dim["count"])
    x0 = _qty(dim["coordinates_offset"], "Hz")
    dx = _qty(dim["increment"], "Hz")
    sf = _qty(dim["origin_offset"], "MHz")
    if dv.get("encoding") != "base64":
        raise ValueError(f"unsupported CSDM encoding {dv.get('encoding')!r}")
    dtype = {"float32": "<f4", "float64": "<f8"}[dv.get("numeric_type", "float32")]
    comp = dv["components"][0]
    amp = np.frombuffer(base64.b64decode(comp), dtype=dtype).astype(float)
    if amp.size != n:
        raise ValueError(f"CSDM component length {amp.size} != count {n}")
    header = {"NP": float(n), "X0": x0, "dX": dx, "Sf": sf, "Sr": 0.0}
    return (DmfitSpectrum(header=header, amplitude=amp,
                          imaginary=np.zeros(n)),
            str(cs.get("description", "")))


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

#: dmfit's "Amorphous" amp scales the integrated AREA (its Gaus/Lor amp is the
#: peak). To make an imported Amorphous line overlay at the correct height
#: relative to the Gaus/Lor lines, convert the area amp to a peak amp:
#:     peak_amp = dmfit_amp * _DMFIT_AMORPHOUS_AMP_FACTOR / (unit-peak area, ppm)
#: The factor is dmfit's area<->peak convention constant, calibrated on the
#: Piepel0 11B BO3 fit (recovers the measured BO3:BO4 height ratio to ~2 %); like
#: the Czjzek export factor it is a dmfit-interop calibration, not a LARMOR
#: physical quantity, so verify it if you rely on absolute imported amplitudes.
_DMFIT_AMORPHOUS_AMP_FACTOR = 3.55


def _rescale_amorphous_amplitudes(recipe: Recipe, warnings: list[str]) -> None:
    """Convert dmfit Amorphous *area* amplitudes to LARMOR *peak* amplitudes so an
    imported fit overlays with the right relative heights (see the factor above).

    Renders each Amorphous site once at unit peak to measure its area in ppm; if
    the render is unavailable (no mrsimulator) the amplitudes are left as dmfit's
    area values and a warning is added."""
    sites = [s for s in recipe.sites if s.model == "amorphous"]
    if not sites:
        return
    try:
        import numpy as np

        from larmor import engine

        pos = [s.params["isotropic_chemical_shift_ppm"].value for s in sites]
        x = np.linspace(min(pos) - 150.0, max(pos) + 150.0, 4000)
        for s in sites:
            probe = SiteModel(
                model="amorphous", label="_probe",
                params={k: Param(p.value) for k, p in s.params.items()})
            probe.params["amplitude"] = Param(1.0)          # unit peak
            r = Recipe(nucleus=recipe.nucleus,
                       larmor_frequency_MHz=recipe.larmor_frequency_MHz,
                       spin_rate_Hz=recipe.spin_rate_Hz, sites=[probe])
            _, _, per = engine.simulate(r, exp_ppm=x)
            area = float(np.trapezoid(np.clip(per[0], 0.0, None), x))
            if area > 0:
                s.params["amplitude"].value *= _DMFIT_AMORPHOUS_AMP_FACTOR / area
    except Exception as exc:  # noqa: BLE001 - import must not fail on render issues
        warnings.append(
            f"Amorphous area->peak amplitude conversion skipped ({exc}); imported "
            "amplitudes are dmfit's area values -- refit to your data")


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
        elif line.model_name == "Amorphous":
            p = line.params
            lb = p.get("lb", DmfitParam(0.0))
            site = SiteModel(
                model="amorphous",
                label=line.name or f"Amorphous-{i}",
                params={
                    "isotropic_chemical_shift_ppm": Param(p["pos"].value),
                    # dmfit stores CQ and FWHM_CQ in kHz; LARMOR uses MHz
                    "Cq_MHz": Param(min(p["CQ"].value / 1000.0, 6.0),
                                    min=0.05, max=6.0),
                    "eta": Param(p.get("etaQ", DmfitParam(0.0)).value,
                                 min=0.0, max=1.0),
                    "Cq_fwhm_MHz": Param(
                        p.get("FWHM_CQ", DmfitParam(0.0)).value / 1000.0,
                        min=0.0, max=3.0),
                    "eta_fwhm": Param(p.get("FWHM_etaQ", DmfitParam(0.0)).value,
                                      min=0.0, max=1.0),
                    "shift_fwhm_ppm": Param(
                        abs(p.get("dCS", DmfitParam(0.0)).value), min=0.0),
                    "line_fwhm_ppm": Param(abs(lb.value), min=0.0),
                    "gl": Param(p.get("gl", DmfitParam(0.0)).value,
                                vary=False, min=0.0, max=1.0),
                    "amplitude": Param(abs(p["amp"].value), min=0.0),
                },
            )
            recipe.sites.append(site)
            if lb.value < 0:
                warnings.append(
                    f"line {i} 'Amorphous' at {p['pos'].value:.1f} ppm: dmfit lb "
                    "was negative (resolution enhancement); imported as 0 -- refit "
                    "the line broadening if needed."
                )
        elif line.model_name == "ss band":
            warnings.append(
                f"line {i} ('ss band' at {line.params.get('pos', DmfitParam(0)).value:.1f} ppm) "
                "skipped: explicit sideband lines are handled by the simulation itself in LARMOR"
            )
        else:
            warnings.append(f"line {i} (model {line.model_name!r}) not yet supported, skipped")

    # dmfit's Amorphous amp is an AREA; convert to a peak amp (rendering each line
    # once) so the imported fit overlays with the right relative heights.
    n_amorph = sum(1 for s in recipe.sites if s.model == "amorphous")
    if n_amorph:
        _rescale_amorphous_amplitudes(recipe, warnings)
        recipe.notes.append(
            f"{n_amorph} Amorphous line(s): dmfit area amplitudes converted to "
            f"peak amplitudes (x{_DMFIT_AMORPHOUS_AMP_FACTOR}/area) for correct "
            "relative heights; absolute scale still needs a fit."
        )

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
