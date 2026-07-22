"""Isotope database + periodic-table layout for the NMR table utility.

Numeric data (spin, γ, Q, abundance) come from mrsimulator's ISOTOPE_DATA; this
module adds element names and a standard periodic-table grid, and the Larmor-
frequency algebra (set B0 in T, or set your magnet's ¹H frequency and read every
nucleus).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

GAMMA_1H = 42.57747920984721        # MHz/T (mrsimulator)

_ELEMENTS = [
    (1, "H", "Hydrogen"), (2, "He", "Helium"), (3, "Li", "Lithium"),
    (4, "Be", "Beryllium"), (5, "B", "Boron"), (6, "C", "Carbon"),
    (7, "N", "Nitrogen"), (8, "O", "Oxygen"), (9, "F", "Fluorine"),
    (10, "Ne", "Neon"), (11, "Na", "Sodium"), (12, "Mg", "Magnesium"),
    (13, "Al", "Aluminium"), (14, "Si", "Silicon"), (15, "P", "Phosphorus"),
    (16, "S", "Sulfur"), (17, "Cl", "Chlorine"), (18, "Ar", "Argon"),
    (19, "K", "Potassium"), (20, "Ca", "Calcium"), (21, "Sc", "Scandium"),
    (22, "Ti", "Titanium"), (23, "V", "Vanadium"), (24, "Cr", "Chromium"),
    (25, "Mn", "Manganese"), (26, "Fe", "Iron"), (27, "Co", "Cobalt"),
    (28, "Ni", "Nickel"), (29, "Cu", "Copper"), (30, "Zn", "Zinc"),
    (31, "Ga", "Gallium"), (32, "Ge", "Germanium"), (33, "As", "Arsenic"),
    (34, "Se", "Selenium"), (35, "Br", "Bromine"), (36, "Kr", "Krypton"),
    (37, "Rb", "Rubidium"), (38, "Sr", "Strontium"), (39, "Y", "Yttrium"),
    (40, "Zr", "Zirconium"), (41, "Nb", "Niobium"), (42, "Mo", "Molybdenum"),
    (43, "Tc", "Technetium"), (44, "Ru", "Ruthenium"), (45, "Rh", "Rhodium"),
    (46, "Pd", "Palladium"), (47, "Ag", "Silver"), (48, "Cd", "Cadmium"),
    (49, "In", "Indium"), (50, "Sn", "Tin"), (51, "Sb", "Antimony"),
    (52, "Te", "Tellurium"), (53, "I", "Iodine"), (54, "Xe", "Xenon"),
    (55, "Cs", "Caesium"), (56, "Ba", "Barium"), (57, "La", "Lanthanum"),
    (58, "Ce", "Cerium"), (59, "Pr", "Praseodymium"), (60, "Nd", "Neodymium"),
    (61, "Pm", "Promethium"), (62, "Sm", "Samarium"), (63, "Eu", "Europium"),
    (64, "Gd", "Gadolinium"), (65, "Tb", "Terbium"), (66, "Dy", "Dysprosium"),
    (67, "Ho", "Holmium"), (68, "Er", "Erbium"), (69, "Tm", "Thulium"),
    (70, "Yb", "Ytterbium"), (71, "Lu", "Lutetium"), (72, "Hf", "Hafnium"),
    (73, "Ta", "Tantalum"), (74, "W", "Tungsten"), (75, "Re", "Rhenium"),
    (76, "Os", "Osmium"), (77, "Ir", "Iridium"), (78, "Pt", "Platinum"),
    (79, "Au", "Gold"), (80, "Hg", "Mercury"), (81, "Tl", "Thallium"),
    (82, "Pb", "Lead"), (83, "Bi", "Bismuth"), (84, "Po", "Polonium"),
    (85, "At", "Astatine"), (86, "Rn", "Radon"), (87, "Fr", "Francium"),
    (88, "Ra", "Radium"), (89, "Ac", "Actinium"), (90, "Th", "Thorium"),
    (91, "Pa", "Protactinium"), (92, "U", "Uranium"), (93, "Np", "Neptunium"),
    (94, "Pu", "Plutonium"), (95, "Am", "Americium"), (96, "Cm", "Curium"),
]
ELEMENT_NAME = {sym: name for _, sym, name in _ELEMENTS}
SYMBOL_OF_Z = {z: sym for z, sym, _ in _ELEMENTS}

#: standard periodic-table grid (row-major), None = gap; f-block on rows 8-9
PERIODIC_ROWS = [
    ["H", None, None, None, None, None, None, None, None, None, None, None,
     None, None, None, None, None, "He"],
    ["Li", "Be", None, None, None, None, None, None, None, None, None, None,
     "B", "C", "N", "O", "F", "Ne"],
    ["Na", "Mg", None, None, None, None, None, None, None, None, None, None,
     "Al", "Si", "P", "S", "Cl", "Ar"],
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
     "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
     "In", "Sn", "Sb", "Te", "I", "Xe"],
    ["Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
     "Tl", "Pb", "Bi", "Po", "At", "Rn"],
    ["Fr", "Ra", "Ac", None, None, None, None, None, None, None, None, None,
     None, None, None, None, None, None],
    [None, None, None, "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
     "Ho", "Er", "Tm", "Yb", "Lu", None],
    [None, None, None, "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", None, None,
     None, None, None, None, None, None],
]


@dataclass(frozen=True)
class Isotope:
    symbol: str                     # e.g. "27Al"
    element: str                    # "Al"
    mass: int                       # 27
    atomic_number: int
    spin: float                     # 2.5
    gamma_MHz_T: float              # signed, MHz/T
    quad_moment_barn: float
    abundance: float                # %

    @property
    def name(self) -> str:
        return ELEMENT_NAME.get(self.element, self.element)

    def larmor_MHz(self, B0_T: float) -> float:
        return abs(self.gamma_MHz_T) * B0_T

    def larmor_from_1H(self, h1_MHz: float) -> float:
        return abs(self.gamma_MHz_T) / GAMMA_1H * h1_MHz

    @property
    def receptivity_1H(self) -> float:
        """Receptivity relative to ¹H = (γ/γH)³ · abundance · I(I+1)/(0.75)."""
        rel = (abs(self.gamma_MHz_T) / GAMMA_1H) ** 3
        return rel * (self.abundance / 100.0) * (self.spin * (self.spin + 1.0) / 0.75)


def _split_symbol(sym: str) -> tuple[int, str]:
    i = 0
    while i < len(sym) and sym[i].isdigit():
        i += 1
    return int(sym[:i] or 0), sym[i:]


@lru_cache(maxsize=1)
def all_isotopes() -> list[Isotope]:
    from mrsimulator.spin_system.isotope import ISOTOPE_DATA

    out = []
    for sym, d in ISOTOPE_DATA.items():
        mass, el = _split_symbol(sym)
        if el not in ELEMENT_NAME:
            continue
        spin = (d["spin_multiplicity"] - 1) / 2.0
        out.append(Isotope(
            symbol=sym, element=el, mass=mass,
            atomic_number=int(d["atomic_number"]), spin=spin,
            gamma_MHz_T=float(d["gyromagnetic_ratio"]),
            quad_moment_barn=float(d["quadrupole_moment"]),
            abundance=float(d["natural_abundance"])))
    return out


@lru_cache(maxsize=1)
def isotopes_by_element() -> dict[str, list[Isotope]]:
    d: dict[str, list[Isotope]] = {}
    for iso in all_isotopes():
        if iso.spin > 0:                       # NMR-active only
            d.setdefault(iso.element, []).append(iso)
    for el in d:
        d[el].sort(key=lambda i: i.receptivity_1H, reverse=True)
    return d


def primary_isotope(element: str) -> Isotope | None:
    """Most receptive NMR-active isotope of an element (for the table cell)."""
    lst = isotopes_by_element().get(element)
    return lst[0] if lst else None


def b0_from_1H(h1_MHz: float) -> float:
    return h1_MHz / GAMMA_1H
