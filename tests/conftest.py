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
