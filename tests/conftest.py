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
}


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
