"""SIMPSON bridge: export a LARMOR spin system, run it, read the result back.

Why: mrsimulator is fast because it assumes weak coupling and no rotational
resonance. Recoupling experiments (REDOR, RFDR, ...) deliberately violate
that -- they exist to reintroduce a coupling MAS removes. For those, an exact
density-matrix simulation is required, and SIMPSON (Bak, Rasmussen & Nielsen,
J. Magn. Reson. 147, 296, 2000) is the reference tool.

LARMOR does NOT bundle SIMPSON. It writes the input, runs the binary if it is
on PATH, and parses the output. Everything here works without SIMPSON
installed except `run`.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from larmor.recipe import Recipe


def simpson_available() -> str | None:
    """Path to the simpson binary, or None."""
    return shutil.which("simpson")


def spinsys_block(recipe: Recipe, site_index: int = 0,
                  partner: dict | None = None) -> str:
    """SIMPSON `spinsys` section from a LARMOR site.

    partner: optional {"isotope": "19F", "dipolar_hz": -1234.0} to add a
    coupled spin (the REDOR case).
    """
    site = recipe.sites[site_index]
    p = {k: v.value for k, v in site.params.items()}
    iso = p.get("isotope", recipe.nucleus)
    lines = ["spinsys {"]
    channels = [iso] + ([partner["isotope"]] if partner else [])
    lines.append("  channels " + " ".join(dict.fromkeys(channels)))
    lines.append("  nuclei " + " ".join(channels))
    # shift in ppm -> SIMPSON wants p (ppm) for shift
    lines.append(f"  shift 1 {p.get('isotropic_chemical_shift_ppm', 0.0)}p "
                 f"{p.get('zeta_ppm', 0.0)}p {p.get('eta_cs', p.get('eta', 0.0)) if 'zeta_ppm' in p else 0.0} 0 0 0")
    cq = p.get("Cq_MHz")
    if cq:
        eta_q = p.get("eta_q", p.get("eta", 0.0))
        lines.append(f"  quadrupole 1 2 {cq * 1e6:.6g} {eta_q} 0 0 0")
    if partner:
        lines.append(f"  dipole 1 2 {partner['dipolar_hz']:.6g} 0 0 0")
    lines.append("}")
    return "\n".join(lines)


REDOR_TEMPLATE = """{spinsys}

par {{
  method           direct
  proton_frequency {proton_freq:.6g}
  spin_rate        {spin_rate:.6g}
  crystal_file     {crystal}
  gamma_angles     {gamma_angles}
  start_operator   I1x
  detect_operator  I1p
  np               {np_points}
  sw               {sw:.6g}
  variable tr      {tr:.10g}
  verbose          0000
}}

proc pulseq {{}} {{
  global par
  maxdt 1.0
  set tr2 [expr $par(tr)/2.0]
  for {{set i 0}} {{$i < $par(np)}} {{incr i}} {{
    acq
    delay $tr2
    pulse_shaped ... ;# placeholder: recoupling pulses inserted per experiment
    delay $tr2
  }}
}}

proc main {{}} {{
  global par
  set f [fsimpson]
  fsave $f $par(name).fid
}}
"""


def redor_input(recipe: Recipe, partner_isotope: str, dipolar_hz: float,
                spin_rate_hz: float, n_points: int = 16,
                proton_freq_hz: float = 400e6,
                crystal: str = "rep168", gamma_angles: int = 16) -> str:
    """A REDOR SIMPSON input for the current spin system.

    The pulse sequence body is intentionally a template stub: recoupling
    schemes differ per lab (xy-8, xy-16, ...), and silently guessing one
    would produce authoritative-looking but wrong curves. Fill it in, or use
    the analytic `larmor.redor` module which needs no simulation at all.
    """
    ss = spinsys_block(recipe, 0,
                       {"isotope": partner_isotope, "dipolar_hz": dipolar_hz})
    tr = 1.0 / spin_rate_hz
    return REDOR_TEMPLATE.format(
        spinsys=ss, proton_freq=proton_freq_hz, spin_rate=spin_rate_hz,
        crystal=crystal, gamma_angles=gamma_angles, np_points=n_points,
        sw=spin_rate_hz, tr=tr)


@dataclass
class SimpsonResult:
    x: np.ndarray
    y: np.ndarray
    stdout: str = ""


def parse_fid(path: str | Path) -> SimpsonResult:
    """Read a SIMPSON .fid (SIMP text format) into arrays."""
    text = Path(path).read_text().splitlines()
    header = {}
    data_start = None
    for i, ln in enumerate(text):
        s = ln.strip()
        if s == "DATA":
            data_start = i + 1
            break
        if "=" in s:
            k, _, v = s.partition("=")
            header[k.strip()] = v.strip()
    if data_start is None:
        raise ValueError(f"{path} is not a SIMPSON SIMP file (no DATA line)")
    rows = []
    for ln in text[data_start:]:
        s = ln.strip()
        if s in ("END", ""):
            continue
        parts = s.split()
        if len(parts) >= 2:
            rows.append((float(parts[0]), float(parts[1])))
    arr = np.array(rows)
    n = len(arr)
    sw = float(header.get("SW", n))
    x = np.arange(n) / sw
    return SimpsonResult(x=x, y=arr[:, 0] + 1j * arr[:, 1])


def run(input_text: str, workdir: str | Path | None = None,
        timeout: int = 600) -> SimpsonResult:
    """Run SIMPSON on `input_text` and return the parsed fid.

    Raises RuntimeError with a clear message when SIMPSON is not installed --
    LARMOR never pretends to have simulated something it did not.
    """
    exe = simpson_available()
    if not exe:
        raise RuntimeError(
            "SIMPSON is not on PATH. Install it from "
            "https://inano.au.dk/about/research-centers/nmr/software/simpson "
            "or use the analytic larmor.redor module instead.")
    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="larmor_simpson_"))
    tmp.mkdir(parents=True, exist_ok=True)
    inp = tmp / "sim.in"
    inp.write_text(input_text)
    proc = subprocess.run([exe, str(inp)], cwd=str(tmp), capture_output=True,
                          text=True, timeout=timeout)
    fid = tmp / "sim.fid"
    if not fid.exists():
        raise RuntimeError(f"SIMPSON produced no output.\n"
                           f"stdout: {proc.stdout[-2000:]}\n"
                           f"stderr: {proc.stderr[-2000:]}")
    out = parse_fid(fid)
    out.stdout = proc.stdout
    return out
