"""DFT / MD tensor import: turn computed tensors into fittable LARMOR sites.

Supported inputs
  * CASTEP / QE `.magres` (the community standard, magres 1.0 text format):
    per-atom shielding (`ms`) and electric field gradient (`efg`) tensors,
    plus the `[calculation]` header (code, functional, cutoffs,
    pseudopotentials) that decides whether two files are comparable.
  * Plain 3x3 tensors from any code (Gaussian/ORCA log parsing is left to the
    user's tooling; the conversion below is the part that is easy to get wrong).

The conversion math comes from mrsimulator (Haeberlen convention, EFG -> Cq
with the isotope's quadrupole moment) -- never re-derived here.

Seeding is ROUTED BY ROLE through ``_SEED_KEYS``: each model declares which
of its parameters receive the shielding anisotropy, the shielding asymmetry,
C_Q and the quadrupolar asymmetry. Two models call different quantities
``eta`` (csa_mas: shielding asymmetry; quad_ct: quadrupolar asymmetry), so
matching on the shared NAME once wrote the quadrupolar eta into a CSA site
and handed a spin-1/2 nucleus a 3 MHz C_Q. ``_SEED_KEYS`` / ``_NOT_SEEDABLE``
partition the registry (tests/test_model_tables.py), so a new model must say
what it seeds before the import offers it.

Shielding is not shift: the isotropic shielding sigma is converted through a
calibration line delta = a*sigma + b (larmor.shiftcal) or a single sigma_ref
(a = -1). The seeded CSA anisotropy follows the same line: LARMOR's
``zeta_ppm`` is the SHIELDING anisotropy as mrsimulator takes it
(delta_aniso = -zeta), so with slope a the seeded value is -a*zeta_sigma
(identical to zeta_sigma for sigma_ref).

LARMOR is a CONSUMER of DFT output: it never runs a DFT code.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = [
    "ComputedSite", "MagresCalc", "MagresFile", "DEFAULT_ISOTOPES",
    "read_magres", "read_calculation", "read_magres_file", "assign_isotopes",
    "sites_for_isotope", "spin_of", "element_of", "group_equivalent",
    "seedable_models", "default_model", "multiplicity_links",
    "sites_to_recipe_dicts", "import_provenance", "model_choice_label",
]

#: 1 Ry in eV (CODATA 2018), for CASTEP headers that quote cutoffs in eV
_EV_PER_RY = 13.605693122994


def spin_of(isotope: str) -> float | None:
    """Nuclear spin of an isotope symbol ('19F' -> 0.5, '27Al' -> 2.5) from
    mrsimulator's table; None when unknown or when mrsimulator is missing."""
    if not isotope:
        return None
    try:
        from mrsimulator.spin_system.isotope import Isotope

        return float(Isotope(symbol=str(isotope)).spin)
    except Exception:
        return None


def element_of(isotope: str) -> str:
    """'19F' -> 'F', '27Al' -> 'Al' (the mass number stripped)."""
    return re.sub(r"^\d+", "", str(isotope or "")).strip()


# --------------------------------------------------------------------------
#: which parameter of each model receives (shielding zeta, shielding eta,
#: C_Q, quadrupolar eta) -- None = the model has no such parameter. With
#: ``_NOT_SEEDABLE`` this partitions models.REGISTRY (the fifth tested
#: partition table); order = the order the import dialog lists them.
_SEED_KEYS: dict[str, tuple[str | None, str | None, str | None, str | None]] = {
    "quad_ct": (None, None, "Cq_MHz", "eta"),
    "quad_csa": ("zeta_ppm", "eta_cs", "Cq_MHz", "eta_q"),
    "quad_first": (None, None, "Cq_MHz", "eta"),
    "ext_czjzek": (None, None, "Cq_MHz", "eta"),
    "amorphous": (None, None, "Cq_MHz", "eta"),
    "csa_mas": ("zeta_ppm", "eta", None, None),
    "csa_czjzek": ("zeta_ppm", "eta", None, None),
    "czjzek": (None, None, None, None),          # position only
    "czjzek_d": (None, None, None, None),
    "czjzek_corr": (None, None, None, None),
    "gl_norm": (None, None, None, None),
    "gauss_lor": (None, None, None, None),
    "voigt": (None, None, None, None),
}

#: models a computed tensor cannot seed: a multiplet needs J, a sideband
#: comb its ratio, exchange its rate, a function its expression and a
#: spectrum component its reference trace
_NOT_SEEDABLE = frozenset({"jmultiplet", "sidebands", "exchange2", "function",
                           "spectrum"})

#: models whose ``amplitude`` is an AREA, so a multiplicity ratio between
#: amplitudes is exactly a population ratio (height models only approximate it)
_AREA_AMPLITUDE = frozenset({"gl_norm"})


@dataclass
class ComputedSite:
    label: str                       # e.g. "Al1"
    isotope: str                     # e.g. "27Al"
    index: int
    ms_tensor: np.ndarray | None = None      # 3x3 shielding, ppm
    efg_tensor: np.ndarray | None = None     # 3x3 EFG, atomic units
    notes: list[str] = field(default_factory=list)
    #: how many symmetry-equivalent atoms this site stands for, and which
    multiplicity: int = 1
    members: list[str] = field(default_factory=list)

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
        """{Cq_MHz, eta} using the isotope's quadrupole moment.

        None for a spin-1/2 nucleus even when the file carries an EFG record:
        a spin-1/2 nucleus has no quadrupole moment, so the record is a
        property of the crystal, not of the spectrum."""
        if self.efg_tensor is None:
            return None
        spin = spin_of(self.isotope)
        if spin is not None and spin < 1.0:
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

    def to_site_dict(self, model: str = "quad_ct", *, calibration=None,
                     reference_ppm: float | None = None) -> dict:
        """A LARMOR recipe site seeded from the computed tensors.

        calibration: a ``larmor.shiftcal.ShiftCalibration`` (delta = a*sigma
        + b); reference_ppm: the one-number alternative, sigma_ref
        (delta = sigma_ref - sigma_iso, slope -1). Without either the
        ISOTROPIC SHIELDING is written as-is and flagged -- computed
        shieldings are not chemical shifts and pretending otherwise is a
        classic silent error.

        Tensor components go to the parameters ``_SEED_KEYS`` names for the
        model; a quadrupolar model on a spin-1/2 nucleus and a model in
        ``_NOT_SEEDABLE`` raise ValueError.
        """
        from larmor import models as model_registry

        if model in _NOT_SEEDABLE:
            raise ValueError(f"{model} cannot be seeded from a computed tensor")
        if model not in _SEED_KEYS:
            model_registry.get(model)          # raises the registry's message
            raise ValueError(f"{model} has no entry in dft._SEED_KEYS")
        m = model_registry.get(model)
        spin = spin_of(self.isotope)
        if m.needs_quadrupolar and spin is not None and spin <= 0.5:
            raise ValueError(f"{model} needs a quadrupolar nucleus; "
                             f"{self.isotope} is spin-1/2")
        if calibration is None and reference_ppm is not None:
            from larmor import shiftcal

            calibration = shiftcal.from_sigma_ref(reference_ppm, self.isotope)
        zeta_key, eta_cs_key, cq_key, eta_q_key = _SEED_KEYS[model]

        params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                           "min": p.min, "max": p.max, "expr": None}
                  for p in m.params}
        notes = list(self.notes)

        sh = self.shielding()
        if sh:
            if calibration is not None:
                params["isotropic_chemical_shift_ppm"]["value"] = \
                    float(calibration.delta(sh["iso_ppm"]))
                zeta_scale = float(calibration.zeta_scale)
            else:
                params["isotropic_chemical_shift_ppm"]["value"] = sh["iso_ppm"]
                notes.append("value is the ISOTROPIC SHIELDING, not a chemical "
                             "shift: give reference_ppm to convert")
                zeta_scale = 1.0
            if zeta_key is not None:
                params[zeta_key]["value"] = zeta_scale * sh["zeta_ppm"]
            if eta_cs_key is not None:
                params[eta_cs_key]["value"] = sh["eta"]

        q = self.quadrupolar()
        if q and cq_key is not None:
            params[cq_key]["value"] = q["Cq_MHz"]
            if eta_q_key is not None:
                params[eta_q_key]["value"] = q["eta"]
        elif self.efg_tensor is not None and cq_key is None:
            notes.append("EFG tensor ignored (model has no quadrupolar "
                         "parameters)")

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


# ---------------------------------------------------------- [calculation]
@dataclass
class MagresCalc:
    """The `[calculation]` header of a magres file: what decides whether two
    computed shieldings sit on the same absolute scale."""
    code: str = ""
    version: str = ""
    xc: str = ""
    cutoff_wfc_Ry: float | None = None
    cutoff_rho_Ry: float | None = None
    pspots: dict[str, str] = field(default_factory=dict)
    kgrid: tuple[int, ...] = ()
    prefix: str = ""
    extra: dict = field(default_factory=dict)

    def describe(self) -> str:
        """'QE-GIPAW 7.5 · PBE · 60/480 Ry' (whatever the header has)."""
        bits = []
        if self.code or self.version:
            bits.append(" ".join(b for b in (self.code, self.version) if b))
        if self.xc:
            bits.append(self.xc)
        if self.cutoff_wfc_Ry is not None:
            cut = f"{self.cutoff_wfc_Ry:g}"
            if self.cutoff_rho_Ry is not None:
                cut += f"/{self.cutoff_rho_Ry:g}"
            bits.append(cut + " Ry")
        return " · ".join(bits) if bits else "no calculation header"

    def fingerprint(self, element: str) -> dict:
        """The six settings that must agree for two files to share a
        shielding scale: code, version, functional, both cutoffs and the
        pseudopotential of the OBSERVED element. The k-grid and the other
        elements' pseudopotentials are cell-specific and are recorded, not
        compared."""
        el = element_of(element)
        return {"calc_code": self.code, "calc_code_version": self.version,
                "calc_xcfunctional": self.xc,
                "calc_cutoffenergy": self.cutoff_wfc_Ry,
                "calc_cutoffenergy_rho": self.cutoff_rho_Ry,
                f"calc_pspot[{el}]": self.pspots.get(el, "")}

    def to_dict(self) -> dict:
        return {"code": self.code, "version": self.version, "xc": self.xc,
                "cutoff_wfc_Ry": self.cutoff_wfc_Ry,
                "cutoff_rho_Ry": self.cutoff_rho_Ry,
                "pspots": dict(self.pspots), "kgrid": list(self.kgrid),
                "prefix": self.prefix, "extra": dict(self.extra)}

    @classmethod
    def from_dict(cls, d: dict | None) -> "MagresCalc":
        d = dict(d or {})
        return cls(code=str(d.get("code", "") or ""),
                   version=str(d.get("version", "") or ""),
                   xc=str(d.get("xc", "") or ""),
                   cutoff_wfc_Ry=_float_or_none(d.get("cutoff_wfc_Ry")),
                   cutoff_rho_Ry=_float_or_none(d.get("cutoff_rho_Ry")),
                   pspots=dict(d.get("pspots") or {}),
                   kgrid=tuple(int(k) for k in (d.get("kgrid") or ())),
                   prefix=str(d.get("prefix", "") or ""),
                   extra=dict(d.get("extra") or {}))


def _float_or_none(v) -> float | None:
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def _cutoff_Ry(tokens: list[str]) -> float | None:
    """'60.00 Ry' / '816.34 eV' / '60' -> Ry."""
    if not tokens:
        return None
    try:
        val = float(tokens[0])
    except ValueError:
        return None
    unit = tokens[1].lower() if len(tokens) > 1 else "ry"
    if unit.startswith("ev"):
        return val / _EV_PER_RY
    if unit.startswith("ha"):
        return val * 2.0
    return val


def read_calculation(path: str | Path) -> MagresCalc:
    """Parse the `[calculation]` block (tolerant: a CASTEP header without
    the QE cutoff keys, or no block at all, gives a partly empty record)."""
    text = Path(path).read_text(errors="replace")
    calc = MagresCalc()
    in_block = False
    for raw in text.splitlines():
        line = raw.strip()
        m = _MAGRES_BLOCK.fullmatch(line)
        if m:
            in_block = (m.group(2) == "calculation") and not m.group(1)
            continue
        if not in_block or not line or line.startswith("#"):
            continue
        parts = line.split()
        key, vals = parts[0], parts[1:]
        value = " ".join(vals)
        if key == "calc_code":
            calc.code = value
        elif key == "calc_code_version":
            calc.version = value
        elif key == "calc_xcfunctional":
            calc.xc = value
        elif key == "calc_cutoffenergy":
            calc.cutoff_wfc_Ry = _cutoff_Ry(vals)
        elif key == "calc_cutoffenergy_rho":
            calc.cutoff_rho_Ry = _cutoff_Ry(vals)
        elif key == "calc_pspot":
            if not vals:
                continue
            # QE: 'F.pbe-n-kjpaw_psl.1.0.0.UPF'; CASTEP: 'F 2|1.4|...'
            if len(vals) > 1 and re.fullmatch(r"[A-Z][a-z]?", vals[0]):
                calc.pspots[vals[0]] = " ".join(vals[1:])
            else:
                el = re.split(r"[._\-]", vals[0], maxsplit=1)[0]
                calc.pspots[el] = value
        elif key == "calc_kpoint_mp_grid":
            try:
                calc.kgrid = tuple(int(v) for v in vals)
            except ValueError:
                calc.extra[key] = value
        elif key == "calc_prefix":
            calc.prefix = value
        elif value:
            calc.extra[key] = value
    return calc


@dataclass
class MagresFile:
    path: str
    sites: list[ComputedSite]
    calc: MagresCalc
    sha256: str

    @property
    def name(self) -> str:
        return Path(self.path).name


def read_magres_file(path: str | Path) -> MagresFile:
    """read_magres + read_calculation + the file's SHA-256 (for provenance)."""
    from larmor.recipe import sha256_of

    p = Path(path)
    return MagresFile(path=str(p), sites=read_magres(p),
                      calc=read_calculation(p), sha256=sha256_of(p))


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


# ------------------------------------------------------- equivalent atoms
def _features(site: ComputedSite) -> tuple:
    sh = site.shielding()
    q = site.quadrupolar()
    return (sh["iso_ppm"] if sh else None, sh["zeta_ppm"] if sh else None,
            sh["eta"] if sh else None, q["Cq_MHz"] if q else None,
            q["eta"] if q else None)


def _close(a, b, tol) -> bool:
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    return abs(a - b) <= tol


def _diag_mean(tensors: list[np.ndarray]) -> np.ndarray:
    """Orientation-free average: the mean of the sorted principal values as a
    diagonal tensor. Symmetry-equivalent atoms carry the SAME tensor in
    different orientations, so averaging the raw 3x3 would shrink the
    anisotropy; the invariants (iso, zeta, eta, C_Q) are what agree."""
    eig = [np.sort(np.linalg.eigvalsh(0.5 * (t + t.T))) for t in tensors]
    return np.diag(np.mean(eig, axis=0))


def group_equivalent(sites: list[ComputedSite], tol_ppm: float = 0.05,
                     tol_eta: float = 0.02, tol_cq_MHz: float = 0.01,
                     ) -> list[ComputedSite]:
    """Collapse symmetry-equivalent atoms into one site per crystallographic
    position: greedy clustering in file order on (sigma_iso, zeta, eta_CS,
    C_Q, eta_Q) within the tolerances. The representative keeps the first
    member's label, carries orientation-averaged tensors, ``multiplicity``
    = the number of atoms and ``members`` = their labels."""
    groups: list[list[ComputedSite]] = []
    feats: list[tuple] = []
    for s in sites:
        f = _features(s)
        for g, fg in zip(groups, feats):
            if g[0].isotope != s.isotope:
                continue
            if (_close(f[0], fg[0], tol_ppm) and _close(f[1], fg[1], tol_ppm)
                    and _close(f[2], fg[2], tol_eta)
                    and _close(f[3], fg[3], tol_cq_MHz)
                    and _close(f[4], fg[4], tol_eta)):
                g.append(s)
                break
        else:
            groups.append([s])
            feats.append(f)
    out: list[ComputedSite] = []
    for g in groups:
        first = g[0]
        ms = [s.ms_tensor for s in g if s.ms_tensor is not None]
        efg = [s.efg_tensor for s in g if s.efg_tensor is not None]
        rep = ComputedSite(
            label=first.label, isotope=first.isotope, index=first.index,
            ms_tensor=_diag_mean(ms) if ms else None,
            efg_tensor=_diag_mean(efg) if efg else None,
            notes=list(first.notes), multiplicity=len(g),
            members=[s.label for s in g])
        if len(g) > 1:
            rep.notes.append(f"{len(g)} equivalent atoms: "
                             + ", ".join(rep.members))
        out.append(rep)
    return out


# ------------------------------------------------------------ model choice
def seedable_models(spin: float | None) -> list[str]:
    """The models a computed tensor can seed for a nucleus of this spin,
    recommended first (gl_norm for spin-1/2, quad_ct otherwise), then the
    tensor-seeded models, then the position-only ones. A spin-1/2 nucleus
    never sees a ``needs_quadrupolar`` model."""
    from larmor import models as model_registry

    names = [n for n in _SEED_KEYS if n in model_registry.REGISTRY]
    if spin is not None and spin <= 0.5:
        names = [n for n in names
                 if not model_registry.get(n).needs_quadrupolar]
    rec = "gl_norm" if (spin is not None and spin <= 0.5) else "quad_ct"
    if rec in names:
        names.remove(rec)
        names.insert(0, rec)
    return names


def default_model(spin: float | None, spin_rate_Hz: float | None = None) -> str:
    """quad_ct for a quadrupolar (or unknown) nucleus; for spin-1/2, csa_mas
    when the spectrum is static (spin_rate_Hz == 0.0 exactly -- a static
    spin-1/2 spectrum IS a CSA pattern), else gl_norm (area amplitudes, so
    multiplicity locks are exact populations)."""
    if spin is None or spin > 0.5:
        return "quad_ct"
    if spin_rate_Hz is not None and float(spin_rate_Hz) == 0.0:
        return "csa_mas"
    return "gl_norm"


def model_choice_label(model: str) -> str:
    """Combo text: 'gl_norm — Gauss/Lorentz (area): exact multiplicity locks'."""
    from larmor import models as model_registry

    m = model_registry.get(model)
    zeta_key, _eta, cq_key, _etaq = _SEED_KEYS.get(model, (None,) * 4)
    if model in _AREA_AMPLITUDE:
        tail = "exact multiplicity locks"
    elif zeta_key and cq_key:
        tail = "ζ, η_CS and C_Q, η_Q seeded from the tensors"
    elif zeta_key:
        tail = "ζ, η_CS seeded from the tensor"
    elif cq_key:
        tail = "C_Q, η_Q seeded from the EFG"
    else:
        tail = "position only"
    return f"{model} — {m.label}: {tail}"


# ------------------------------------------------------------ links
def _ratio_str(r: float) -> str:
    return f"{r:g}"


def multiplicity_links(site_dicts: list[dict], multiplicities: list[int],
                       model: str, *, lock_amplitude: bool = True,
                       share_width: bool = True, first_index: int = 0,
                       ) -> None:
    """Link every site after the first to the first (the master, recipe
    index ``first_index``): amplitude = (m_k / m_0) * master amplitude (the
    coefficient omitted when 1, the form cellparse.format_link shows as
    '0.5A'), and one shared breadth parameter (from estimate._WIDTH_KEY when
    it is a width, never a C_Q). Mutates the dicts in place."""
    if not site_dicts:
        return
    from larmor.estimate import _WIDTH_KEY

    wkey, is_cq = _WIDTH_KEY.get(model, (None, True))
    m0 = float(multiplicities[0]) if multiplicities else 1.0
    master = f"s{first_index}"
    for k, sd in enumerate(site_dicts):
        if k == 0:
            continue
        params = sd["params"]
        if lock_amplitude and "amplitude" in params:
            ratio = float(multiplicities[k]) / m0 if m0 else 1.0
            params["amplitude"]["expr"] = (
                f"{master}.amplitude" if ratio == 1.0
                else f"{_ratio_str(ratio)}*{master}.amplitude")
            params["amplitude"]["vary"] = False
        if (share_width and wkey is not None and not is_cq
                and wkey in params):
            params[wkey]["expr"] = f"{master}.{wkey}"
            params[wkey]["vary"] = False


def sites_to_recipe_dicts(groups: list[ComputedSite], model: str, *,
                          calibration=None, reference_ppm: float | None = None,
                          first_index: int = 0, lock_amplitude: bool = True,
                          share_width: bool = True,
                          ) -> tuple[list[dict], list[str]]:
    """The single path the dialog and the CLI share: one recipe site dict
    (model / label / params only -- Recipe.from_dict must never see another
    key) per grouped site, plus the human notes ('s3 F5: ×8 (F5, F6, …),
    δ_pred −108.8 ± 0.4 ppm, amplitude locked to s3')."""
    if calibration is None and reference_ppm is not None:
        from larmor import shiftcal

        calibration = shiftcal.from_sigma_ref(reference_ppm, "")
    dicts: list[dict] = []
    notes: list[str] = []
    extra: list[str] = []
    for k, g in enumerate(groups):
        sd = g.to_site_dict(model, calibration=calibration)
        for n in sd.get("notes", []):
            if n not in extra and not n.startswith(f"{g.multiplicity} equivalent"):
                extra.append(n)
        dicts.append({"model": sd["model"], "label": sd["label"],
                      "params": sd["params"]})
    mults = [max(int(g.multiplicity), 1) for g in groups]
    multiplicity_links(dicts, mults, model, lock_amplitude=lock_amplitude,
                       share_width=share_width, first_index=first_index)
    for k, (g, sd) in enumerate(zip(groups, dicts)):
        idx = first_index + k
        sh = g.shielding()
        bits = [f"s{idx} {g.label}: ×{g.multiplicity}"]
        if g.members and g.multiplicity > 1:
            bits[0] += " (" + ", ".join(g.members) + ")"
        if sh:
            if calibration is not None:
                d = float(calibration.delta(sh["iso_ppm"]))
                e = calibration.delta_err(sh["iso_ppm"])
                bits.append(f"δ_pred {d:.1f} ppm" if e is None
                            else f"δ_pred {d:.1f} ± {e:.1f} ppm")
            else:
                bits.append(f"σ_iso {sh['iso_ppm']:.1f} ppm (shielding, "
                            "unconverted)")
        amp = sd["params"].get("amplitude", {})
        if amp.get("expr"):
            bits.append(f"amplitude locked to s{first_index} "
                        f"({amp['expr']})")
        notes.append(", ".join(bits))
    notes.extend(extra)
    return dicts, notes


def import_provenance(mf: MagresFile, isotope: str, model: str,
                      groups: list[ComputedSite], calibration=None, *,
                      calc_override: list[str] | None = None,
                      lock_amplitude: bool = True, share_width: bool = True,
                      ) -> dict:
    """The ``provenance['dft_import']`` record: file + hash, isotope, model,
    the [calculation] header, every grouped site with its members and
    predicted shift, the conversion (a calibration record or sigma_ref) and
    whether a settings mismatch was overridden."""
    sites = []
    for g in groups:
        sh = g.shielding() or {}
        q = g.quadrupolar() or {}
        d_pred = d_err = None
        if calibration is not None and sh:
            d_pred = float(calibration.delta(sh["iso_ppm"]))
            d_err = calibration.delta_err(sh["iso_ppm"])
        sites.append({
            "label": g.label, "multiplicity": int(g.multiplicity),
            "members": list(g.members) or [g.label],
            "sigma_iso_ppm": sh.get("iso_ppm"),
            "zeta_sigma_ppm": sh.get("zeta_ppm"), "eta_cs": sh.get("eta"),
            "Cq_MHz": q.get("Cq_MHz"), "eta_q": q.get("eta"),
            "delta_pred_ppm": d_pred,
            "delta_err_ppm": None if d_err is None else float(d_err)})
    return {
        "file": mf.path, "sha256": mf.sha256, "isotope": isotope,
        "model": model, "calculation": mf.calc.to_dict(), "sites": sites,
        "calibration": calibration.to_dict() if calibration is not None
        else None,
        "calc_override": list(calc_override or []),
        "lock_amplitude": bool(lock_amplitude),
        "share_width": bool(share_width),
    }
