"""DFT / MD tensor import: turn computed tensors into fittable LARMOR sites.

Supported inputs
  * CASTEP / QE `.magres` (the community standard, magres 1.0 text format):
    per-atom shielding (`ms`) and electric field gradient (`efg`) tensors.
  * Plain 3x3 tensors from any code (Gaussian/ORCA log parsing is left to the
    user's tooling; the conversion below is the part that is easy to get wrong).

The conversion math comes from mrsimulator (Haeberlen convention, EFG -> Cq
with the isotope's quadrupole moment) -- never re-derived here.

LARMOR is a CONSUMER of DFT output: it never runs a DFT code.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class ComputedSite:
    label: str                       # e.g. "Al1"
    isotope: str                     # e.g. "27Al"
    index: int
    ms_tensor: np.ndarray | None = None      # 3x3 shielding, ppm
    efg_tensor: np.ndarray | None = None     # 3x3 EFG, atomic units
    notes: list[str] = field(default_factory=list)

    # ---- derived NMR parameters -------------------------------------
    def shielding(self) -> dict | None:
        """{iso_ppm, zeta_ppm, eta} in the Haeberlen convention."""
        if self.ms_tensor is None:
            return None
        sym = 0.5 * (self.ms_tensor + self.ms_tensor.T)
        eig = np.sort(np.linalg.eigvalsh(sym))
        iso = float(eig.mean())
        # Haeberlen ordering: |zz-iso| >= |xx-iso| >= |yy-iso|
        dev = eig - iso
        order = np.argsort(np.abs(dev))          # yy, xx, zz
        yy, xx, zz = dev[order]
        zeta = float(zz)
        eta = float((yy - xx) / zz) if zz else 0.0
        return {"iso_ppm": iso, "zeta_ppm": zeta, "eta": abs(eta)}

    def quadrupolar(self) -> dict | None:
        """{Cq_MHz, eta} using the isotope's quadrupole moment."""
        if self.efg_tensor is None:
            return None
        from mrsimulator.spin_system.isotope import Isotope

        sym = 0.5 * (self.efg_tensor + self.efg_tensor.T)
        eig = np.linalg.eigvalsh(sym)
        order = np.argsort(np.abs(eig))          # |Vxx| <= |Vyy| <= |Vzz|
        vxx, vyy, vzz = eig[order]
        eta = float(abs((vxx - vyy) / vzz)) if vzz else 0.0
        iso = Isotope(symbol=self.isotope)
        # mrsimulator exposes the atomic-unit EFG -> Cq (Hz) conversion
        cq_hz = float(vzz * iso.efg_to_Cq)
        return {"Cq_MHz": cq_hz / 1e6, "eta": min(eta, 1.0)}

    def to_site_dict(self, model: str = "quad_ct",
                     reference_ppm: float | None = None) -> dict:
        """A LARMOR recipe site seeded from the computed tensors.

        reference_ppm: sigma_ref for shielding -> chemical shift
        (delta = sigma_ref - sigma_iso). Without it the ISOTROPIC SHIELDING is
        written as-is and flagged -- computed shieldings are not chemical
        shifts and pretending otherwise is a classic silent error.
        """
        from larmor import models as model_registry

        m = model_registry.get(model)
        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        notes = list(self.notes)

        sh = self.shielding()
        if sh:
            if reference_ppm is not None:
                params["isotropic_chemical_shift_ppm"]["value"] = \
                    reference_ppm - sh["iso_ppm"]
            else:
                params["isotropic_chemical_shift_ppm"]["value"] = sh["iso_ppm"]
                notes.append("value is the ISOTROPIC SHIELDING, not a chemical "
                             "shift: give reference_ppm to convert")
            if "zeta_ppm" in params:
                params["zeta_ppm"]["value"] = sh["zeta_ppm"]
            if "eta_cs" in params:
                params["eta_cs"]["value"] = sh["eta"]

        q = self.quadrupolar()
        if q:
            if "Cq_MHz" in params:
                params["Cq_MHz"]["value"] = q["Cq_MHz"]
            key = "eta_q" if "eta_q" in params else "eta"
            if key in params:
                params[key]["value"] = q["eta"]

        return {"model": model, "label": self.label, "params": params,
                "notes": notes}


# --------------------------------------------------------------------------
_MAGRES_BLOCK = re.compile(r"\[(/?)(\w+)\]")


def read_magres(path: str | Path) -> list[ComputedSite]:
    """Parse a CASTEP/QE .magres file (magres 1.0 text format).

    Reads the `ms` and `efg` records of the [magres] block; atom labels come
    from the records themselves, so no separate structure file is needed.
    """
    text = Path(path).read_text(errors="replace")
    sites: dict[tuple[str, int], ComputedSite] = {}

    def _tensor(vals: list[str]) -> np.ndarray:
        return np.array([float(v) for v in vals[:9]]).reshape(3, 3)

    in_magres = False
    for raw in text.splitlines():
        line = raw.strip()
        m = _MAGRES_BLOCK.fullmatch(line)
        if m:
            in_magres = (m.group(2) == "magres") and not m.group(1)
            continue
        if not in_magres or not line or line.startswith("#"):
            continue
        parts = line.split()
        if not parts:
            continue
        kind = parts[0]
        if kind not in ("ms", "efg", "efg_local", "efg_nonlocal"):
            continue
        # format: ms <label> <index> <9 tensor components>
        try:
            label, index = parts[1], int(parts[2])
            tensor = _tensor(parts[3:])
        except (ValueError, IndexError):
            continue
        key = (label, index)
        site = sites.get(key)
        if site is None:
            site = ComputedSite(label=f"{label}{index}", isotope="",
                                index=index)
            sites[key] = site
        if kind == "ms":
            site.ms_tensor = tensor
        elif kind == "efg":
            site.efg_tensor = tensor
    if not sites:
        raise ValueError(f"no ms/efg records found in {path} "
                         "(is this a magres file?)")
    return list(sites.values())


#: element -> the isotope solid-state NMR actually measures
DEFAULT_ISOTOPES = {
    "H": "1H", "Li": "7Li", "B": "11B", "C": "13C", "N": "15N", "O": "17O",
    "F": "19F", "Na": "23Na", "Mg": "25Mg", "Al": "27Al", "Si": "29Si",
    "P": "31P", "S": "33S", "Cl": "35Cl", "K": "39K", "Ca": "43Ca",
    "Sc": "45Sc", "Ti": "47Ti", "V": "51V", "Ga": "71Ga", "Ge": "73Ge",
    "Rb": "87Rb", "Sr": "87Sr", "Y": "89Y", "Nb": "93Nb", "Sn": "119Sn",
    "Cs": "133Cs", "Ba": "137Ba", "La": "139La", "Pb": "207Pb",
}


def assign_isotopes(sites: list[ComputedSite],
                    overrides: dict[str, str] | None = None) -> list[str]:
    """Fill in `isotope` from each site's element symbol. Returns warnings."""
    warnings = []
    for s in sites:
        element = re.match(r"([A-Za-z]+)", s.label)
        el = element.group(1) if element else ""
        iso = (overrides or {}).get(el) or DEFAULT_ISOTOPES.get(el)
        if iso:
            s.isotope = iso
        else:
            warnings.append(f"no default isotope for element {el!r} "
                            f"(site {s.label})")
    return warnings


def sites_for_isotope(sites: list[ComputedSite], isotope: str,
                      ) -> list[ComputedSite]:
    return [s for s in sites if s.isotope == isotope]
