from pathlib import Path

import pytest

# Real local data used by the integration tests. Tests that need these skip
# cleanly on machines that don't have them (CI, collaborators).
CAALGLASS = Path("C:/Users/samso/Desktop/CaAlGlass.fxmla")
CAALGLASS_MQ = Path("C:/Users/samso/Desktop/CaAlGlassMQ.fxmla")
EXPNO_1903 = Path(
    "C:/Users/samso/Desktop/WSU_work/NMR/NMRFAM/DATA/2026-07"
    "/07062026_SR31648_0Ca-9F-Al_SS_ALP/1903")
EXPNO_1901 = Path(
    "C:/Users/samso/Desktop/WSU_work/NMR/NMRFAM/DATA/2026-07"
    "/07062026_SR31648_0Ca-9F-Al_SS_ALP/1901")
NMRVEW_2D = Path(
    "C:/Users/samso/Desktop/SAVE_PC_PRO/THESE/NMRVEW/data/B15Na20-17O/1181")

# Universal-reader fixtures: a 1D processed spectrum, a raw fid, a pseudo-2D
# (saturation-recovery) 2rr + ser, and a real MQMAS 2rr.
_D = Path("C:/Users/samso/Desktop/WSU_work/NMR/NMRFAM/DATA")
BRUKER_1R = _D / "2026-05/04272026_P1-Bi1-12_SS_ALP/2702/pdata/1/1r"
BRUKER_FID = _D / "2026-05/04272026_P1-Bi1-12_SS_ALP/2702/fid"
BRUKER_2RR_PSEUDO = _D / "2026-01/01202026_SR31649_Base4Ca_SS_ALP/32/pdata/1/2rr"
BRUKER_SER = _D / "2026-01/01202026_SR31649_Base4Ca_SS_ALP/32/ser"
BRUKER_2RR_MQMAS = _D / "2025-12/neg8A2O-0F/35/pdata/1/2rr"


def require(path: Path):
    if not path.exists():
        pytest.skip(f"local test data not present: {path}")
    return path
