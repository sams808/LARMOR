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
