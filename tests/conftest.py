from pathlib import Path

import pytest

# Real local data used by the integration tests. Tests that need these skip
# cleanly on machines that don't have them (CI, collaborators).
CAALGLASS = Path("C:/Users/samso/Desktop/CaAlGlass.fxmla")
CAALGLASS_MQ = Path("C:/Users/samso/Desktop/CaAlGlassMQ.fxmla")
EXPNO_1903 = Path(
    "C:/Users/samso/Desktop/WSU_work/NMR/NMRFAM/DATA/2026-07"
    "/07062026_SR31648_0Ca-9F-Al_SS_ALP/1903")


def require(path: Path):
    if not path.exists():
        pytest.skip(f"local test data not present: {path}")
    return path
