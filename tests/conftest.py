"""Shared test fixtures, and the real-data layer's root.

Real local datasets used by the integration/acceptance tests. Tests that need
them skip cleanly when they are missing -- WHICH SILENTLY REMOVES THE
STRONGEST REGRESSION LAYER: on the development machine the suite is ~663
passed / 15 skipped, on a machine with none of this data it is still green at
~639 passed / 39 skipped with the published-value acceptance tests gone. The
terminal summary below says loudly which case a run was.

Point ``LARMOR_TEST_DATA`` at a folder holding the same RELATIVE layout to run
the real-data layer elsewhere (the default root is the original development
machine's home, so existing runs are unchanged).
"""
import os
import tempfile
from pathlib import Path

import pytest

#: root of the real-data tree; every dataset path is relative to this
DATA_ROOT = Path(os.environ.get("LARMOR_TEST_DATA", r"C:\Users\samso"))

# No background kernel pre-build in test windows: every MainWindow that
# received a quadrupolar spectrum used to start a KernelWarmWorker thread
# that outlived the fixture's window, and a QThread destroyed while running
# aborts the interpreter -- runs that printed all their dots and died, or
# hung, right after the tests that load 11B/27Al data. The one test of the
# warm-up itself deletes this variable and stubs the worker.
os.environ.setdefault("LARMOR_NO_KERNEL_WARM", "1")

# The MAS-rate confirmation store (larmor.masrate) lives under
# %LOCALAPPDATA%/LARMOR on a real machine. A confirmation recorded there
# would silently change what the real-data cases assert (2702 must load
# FLAGGED), and a test writing there would pollute the developer's store: every
# run gets its own empty store, the way LARMOR_REF_LOG isolates the audit log.
os.environ.setdefault(
    "LARMOR_MAS_LOG",
    str(Path(tempfile.mkdtemp(prefix="larmor_mas_")) / "mas_confirmations.jsonl"))


def _data(rel: str) -> Path:
    return DATA_ROOT / rel


# moved into Desktop/larmor_tests at some point -- the old Desktop paths
# made these acceptance tests skip SILENTLY for who knows how long; the
# terminal-summary manifest below is what finally surfaced it
CAALGLASS = _data("Desktop/larmor_tests/CaAlGlass.fxmla")
CAALGLASS_MQ = _data("Desktop/larmor_tests/CaAlGlassMQ.fxmla")
EXPNO_1903 = _data("Desktop/WSU_work/NMR/NMRFAM/DATA/2026-07"
                   "/07062026_SR31648_0Ca-9F-Al_SS_ALP/1903")
EXPNO_1901 = _data("Desktop/WSU_work/NMR/NMRFAM/DATA/2026-07"
                   "/07062026_SR31648_0Ca-9F-Al_SS_ALP/1901")
NMRVEW_2D = _data("Desktop/SAVE_PC_PRO/THESE/NMRVEW/data/B15Na20-17O/1181")

# Universal-reader fixtures: a 1D processed spectrum, a raw fid, a pseudo-2D
# (saturation-recovery) 2rr + ser, and a real MQMAS 2rr.
_D = _data("Desktop/WSU_work/NMR/NMRFAM/DATA")
BRUKER_1R = _D / "2026-05/04272026_P1-Bi1-12_SS_ALP/2702/pdata/1/1r"
BRUKER_FID = _D / "2026-05/04272026_P1-Bi1-12_SS_ALP/2702/fid"
BRUKER_2RR_PSEUDO = _D / "2026-01/01202026_SR31649_Base4Ca_SS_ALP/32/pdata/1/2rr"
BRUKER_SER = _D / "2026-01/01202026_SR31649_Base4Ca_SS_ALP/32/ser"
BRUKER_2RR_MQMAS = _D / "2025-12/neg8A2O-0F/35/pdata/1/2rr"
# a real DQ/SQ 2D correlation (indirect dim is double-quantum, not a delay)
BRUKER_2RR_DQSQ = _data("Desktop/SAVE_PC_PRO/THESE/RMN/Data/ISG-2GPa"
                        "/1181/pdata/1")
# the published MagLab 35Cl QCPMG acceptance set (12 samples, ssNake T2s)
MAGLAB_35CL = _data("Desktop/WSU_work/NMR/MagLab/DATA/35Cl_2025-12")
# Tutorial 7: the static 81Br WURST-CPMG set (EXPNOs 30-34 = 0-4 Ca glasses;
# acqus records a stale MASR=14000 on a static probe)
MAGLAB_81BR = _data("Desktop/WSU_work/NMR/MagLab/DATA/81Br_2026-08")
# Tutorial 4: the five-glass 11B composition series (EXPNO 24 of each sample)
LAW_CA_11B = tuple(_D / f"2026-01/{s}/24" for s in (
    "01192026_SR31649_Base0Ca_SS_ALP", "01202026_SR31649_Base1Ca_SS_ALP",
    "01202026_SR31649_Base2Ca_SS_ALP", "01202026_SR31648_Base3Ca_SS_ALP",
    "01202026_SR31649_Base4Ca_SS_ALP"))
# DFT tensor import: the QE-GIPAW 7.5 / PBE / 60-480 Ry magres files of the
# two fluoride standards and their 19F MAS spectra (2026-03, 564.27 MHz;
# the acqus MASR = 4200 is stale, the title says 35.714 kHz)
_QE = _data("Desktop/WSU_work/Python/NEW/DFT/QE_work")
CAF2_MAGRES = _QE / "CaF2/CaF2.nmr.magres"
NAF_MAGRES = _QE / "NaF/NaF.nmr.magres"
F19_STD_CAF2 = _D / "2026-03/03132026_SR31649_CaF2_SS_ALP/4"
F19_STD_NAF = _D / "2026-03/03132026_SR31648_NaF_SS_ALP/4"

#: the full manifest, so the summary can report what a MACHINE is missing
#: (not only what this particular selection of tests happened to request)
ALL_DATASETS = {
    "CaAlGlass.fxmla": CAALGLASS, "CaAlGlassMQ.fxmla": CAALGLASS_MQ,
    "19F EXPNO 1903": EXPNO_1903, "19F EXPNO 1901": EXPNO_1901,
    "17O 2D (NMRVEW)": NMRVEW_2D, "27Al 1r": BRUKER_1R, "27Al fid": BRUKER_FID,
    "pseudo-2D 2rr": BRUKER_2RR_PSEUDO, "ser": BRUKER_SER,
    "MQMAS 2rr": BRUKER_2RR_MQMAS, "DQ/SQ 2rr": BRUKER_2RR_DQSQ,
    "MagLab 35Cl QCPMG set": MAGLAB_35CL,
    "81Br WCPMG set (Tutorial 7)": MAGLAB_81BR / "30",
    "LAW Ca 11B series (Tutorial 4)": LAW_CA_11B[0],
    "CaF2 magres (GIPAW)": CAF2_MAGRES, "NaF magres (GIPAW)": NAF_MAGRES,
    "19F CaF2 standard 2026-03": F19_STD_CAF2,
    "19F NaF standard 2026-03": F19_STD_NAF,
}


# ------------------------------------------------- synthetic magres files
def synthetic_magres(records: list[tuple[str, int, list[float]]], *,
                     prefix: str = "test", code: str = "QE-GIPAW",
                     version: str = "7.5", xc: str = "PBE",
                     cutoff_Ry: float | None = 60.0,
                     cutoff_rho_Ry: float | None = 480.0,
                     pspots: tuple[str, ...] = ("Ca.pbe-spn-kjpaw_psl.1.0.0.UPF",
                                                "F.pbe-n-kjpaw_psl.1.0.0.UPF"),
                     kgrid: tuple[int, int, int] = (4, 4, 4),
                     efg: list[tuple[str, int, list[float]]] = ()) -> str:
    """A magres 1.0 text with a QE-GIPAW-style [calculation] header.
    ``records`` are (element, index, 9 ms components); ``efg`` likewise."""
    head = ["#$magres-abinitio-v1.0", "[calculation]", f"calc_code {code}"]
    if version:
        head.append(f"calc_code_version {version}")
    head += [f"calc_prefix {prefix}", f"calc_xcfunctional {xc}"]
    if cutoff_Ry is not None:
        head.append(f"calc_cutoffenergy {cutoff_Ry:.2f} Ry")
    if cutoff_rho_Ry is not None:
        head.append(f"calc_cutoffenergy_rho {cutoff_rho_Ry:.2f} Ry")
    head += [f"calc_pspot {p}" for p in pspots]
    head.append("calc_kpoint_mp_grid " + " ".join(str(k) for k in kgrid))
    head += ["[/calculation]", "[atoms]", "units lattice Angstrom",
             "units atom Angstrom"]
    head += [f"atom {el} {el} {i} 0.0 0.0 0.0" for el, i, _ in records]
    head += ["[/atoms]", "[magres]", "units ms ppm"]
    body = [f"ms {el} {i} " + " ".join(f"{v:.6f}" for v in vals)
            for el, i, vals in records]
    body += [f"efg {el} {i} " + " ".join(f"{v:.6f}" for v in vals)
             for el, i, vals in efg]
    return "\n".join(head + body + ["[/magres]", ""])


def diag9(a: float, b: float, c: float) -> list[float]:
    return [a, 0.0, 0.0, 0.0, b, 0.0, 0.0, 0.0, c]


def caf2_like_records():
    """Two Ca, eight F identical to 1e-3 ppm (CaF2's values) and one F
    5 ppm away -- the grouping fixture."""
    recs = [("Ca", 1, diag9(1153.8779, 1153.8779, 1153.8779)),
            ("Ca", 2, diag9(1153.8690, 1153.8690, 1153.8688))]
    for k in range(8):
        d = 0.0006 * ((-1) ** k)
        recs.append(("F", 3 + k, [232.9226, d, d, d, 232.9226, d, d, d,
                                  232.9226]))
    recs.append(("F", 11, diag9(237.9226, 237.9226, 237.9226)))
    return recs


def kogarkoite_like_records():
    """Six PAIRS of F (twelve atoms) at distinct shieldings, pairwise equal
    to 1e-4 ppm, plus two Na -- the multiplicity fixture."""
    recs = [("Na", 1, diag9(520.0, 521.0, 522.0)),
            ("Na", 2, diag9(520.0, 521.0, 522.0))]
    sig = [250.0, 262.0, 275.0, 291.0, 304.0, 318.0]
    idx = 3
    for s in sig:
        for rep in range(2):
            eps = 1e-4 * rep
            recs.append(("F", idx, diag9(s - 12.0 + eps, s + eps, s + 12.0 + eps)))
            idx += 1
    return recs


def require(path: Path):
    if not path.exists():
        pytest.skip(f"local test data not present: {path}")
    return path


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Say LOUDLY whether the real-data acceptance layer actually ran.

    A green bar with this layer missing validates far less than a green bar
    with it present -- and the difference used to be invisible."""
    tr = terminalreporter
    missing = {name: p for name, p in ALL_DATASETS.items() if not p.exists()}
    if missing:
        tr.write_sep("=", "LARMOR real-data layer: INCOMPLETE", red=True)
        tr.write_line(
            f"{len(missing)}/{len(ALL_DATASETS)} dataset(s) missing under "
            f"LARMOR_TEST_DATA={DATA_ROOT} -- their published-value "
            "acceptance tests SKIP, so this green bar proves less than a "
            "full run:")
        for name, p in sorted(missing.items()):
            tr.write_line(f"  - {name}: {p}")
    else:
        tr.write_sep("-", "LARMOR real-data layer: complete "
                          f"(all {len(ALL_DATASETS)} datasets present)",
                     green=True)



# --------------------------------------------------------------------------
# Simulated 35Cl central-transition patterns for the QCPMG multi-field tests
# (direct mrsimulator: every site carries the same CT area, so the whole-
# manifold centre of gravity obeys the first-moment theorem exactly).

def simulate_ct_single(nucleus: str, larmor_MHz: float, rotor_Hz: float,
                       cq_MHz: float, eta: float, delta_iso_ppm: float,
                       span_ppm=(-600.0, 400.0), npts: int = 8000,
                       shift_fwhm_ppm: float = 3.0, n_ssb: int = 8):
    """(ppm, amp) of one crystalline site's CT pattern (quad_ct's own
    mrsimulator route, ``larmor.models._singlesite.simulate_single_site``),
    shifted to ``delta_iso_ppm`` and Gaussian-broadened."""
    import numpy as np
    from larmor.models._singlesite import simulate_single_site

    xs, ys = simulate_single_site(
        nucleus, int(round(larmor_MHz * 1000)), int(round(rotor_Hz)),
        int(round(cq_MHz * 1000)), int(round(eta * 1000)), 0, 0, True, n_ssb,
        float(span_ppm[0]), float(span_ppm[1]), int(npts))
    x = np.asarray(xs, float) + delta_iso_ppm
    y = np.asarray(ys, float)
    return x, _gauss_broaden(x, y, shift_fwhm_ppm)


def simulate_ct_czjzek(nucleus: str, larmor_MHz: float, rotor_Hz: float,
                       sigma_MHz: float, delta_iso_ppm: float, n: int = 200,
                       seed: int = 7, span_ppm=(-3000.0, 2000.0),
                       npts: int = 16384, shift_fwhm_ppm: float = 10.0,
                       n_ssb: int = 16):
    """(ppm, amp, rms_PQ_of_sample) for a Czjzek distribution of ``n`` sites
    (mrsimulator's CzjzekDistribution in LARMOR's sigma convention), each
    at ``delta_iso_ppm``; the returned rms P_Q is that of the DRAWN sample,
    which is what the first-moment theorem predicts for this spectrum."""
    import numpy as np
    from mrsimulator import Simulator, Site, SpinSystem
    from mrsimulator.method import SpectralDimension
    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.models import CzjzekDistribution
    from mrsimulator.spin_system.isotope import Isotope

    state = np.random.get_state()
    try:
        np.random.seed(seed)
        cq, eta = CzjzekDistribution(sigma=sigma_MHz).rvs(size=n)
    finally:
        np.random.set_state(state)
    b0 = larmor_MHz / abs(Isotope(symbol=nucleus).gyromagnetic_ratio)
    systems = [SpinSystem(sites=[Site(
        isotope=nucleus, isotropic_chemical_shift=delta_iso_ppm,
        quadrupolar={"Cq": float(abs(c)) * 1e6, "eta": float(e)})])
        for c, e in zip(cq, eta)]
    method = BlochDecayCTSpectrum(
        channels=[nucleus], magnetic_flux_density=b0, rotor_frequency=rotor_Hz,
        spectral_dimensions=[SpectralDimension(
            count=npts, spectral_width=(span_ppm[1] - span_ppm[0]) * larmor_MHz,
            reference_offset=0.5 * (span_ppm[0] + span_ppm[1]) * larmor_MHz)])
    sim = Simulator(spin_systems=systems, methods=[method])
    sim.config.number_of_sidebands = n_ssb
    sim.run()
    ds = sim.methods[0].simulation
    x = np.asarray(ds.x[0].coordinates.value, float)
    y = np.asarray(ds.y[0].components[0].real, float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    y = _gauss_broaden(x, y, shift_fwhm_ppm)
    rms = float(np.sqrt(np.mean(cq ** 2 * (1.0 + eta ** 2 / 3.0))))
    return x, y / y.max(), rms


def _gauss_broaden(x, y, fwhm_ppm: float):
    import numpy as np
    if fwhm_ppm <= 0:
        return y
    dx = float(abs(x[1] - x[0]))
    s = fwhm_ppm / 2.3548 / dx
    k = np.arange(-int(5 * s), int(5 * s) + 1)
    g = np.exp(-0.5 * (k / s) ** 2)
    return np.convolve(y, g / g.sum(), "same")

def fake_spectrum_params(sample: str, **over):
    """An in-memory larmor.comparability.SpectrumParams (the 11B series'
    acqus/procs values) for tests that stub comparability.read_params. Any
    acqus / procs key can be overridden by name (lb=100 -> proc LB), plus
    ``path`` / ``expno`` / ``has_fid`` / ``commands``."""
    from larmor.comparability import SpectrumParams

    acq = {"NUC1": "11B", "PULPROG": "zg", "NS": 256, "TD": 7988, "RG": 194.07,
           "SW_h": 100000.0, "SFO1": 192.430693, "BF1": 192.430693, "O1": 0.0,
           "DATE": 1768800000.0, "PROBHD": "16_Solenoid (PMAS16)", "D1": 12.5,
           "P1": 0.425, "PLW1": 100.0, "PLW1_unit": "W"}
    proc = {"WDW": 1, "LB": 0.0, "GB": 0.0, "SSB": 0.0, "TDeff": 1024, "SI": 32768,
            "FCOR": 1.0, "PH_mod": 1, "PHC0": -160.0, "PHC1": 29.9084, "ABSG": 0,
            "BC_mod": 0, "FT_mod": 6, "PKNL": True, "SF": 192.431528,
            "OFFSET": 100.0, "SW_p": 100000.0}
    meta = {k: over.pop(k) for k in ("path", "expno", "has_fid", "commands", "title")
            if k in over}
    for k, v in over.items():                 # case-insensitive: lb -> LB, tdeff -> TDeff
        hits = [(d, kk) for d in (acq, proc) for kk in d if kk.lower() == k.lower()]
        if hits:
            hits[0][0][hits[0][1]] = v
        else:
            proc[k] = v
    return SpectrumParams(
        path=meta.get("path", f"/fake/{sample}/24/pdata/1/1r"),
        expno=meta.get("expno", f"/fake/{sample}/24"), procno=1, sample=sample,
        acq=acq, proc=proc, commands=meta.get("commands", ["efp LB = 0 SI = 32K"]),
        title=meta.get("title", sample), has_fid=bool(meta.get("has_fid", False)))
