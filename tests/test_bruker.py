import pytest

from larmor.io import bruker

from conftest import EXPNO_1903, require


def test_read_expno_1903_readonly():
    path = require(EXPNO_1903)
    before = bruker.snapshot(path)
    exp = bruker.read_expno(path, verify=True)  # raises on any modification
    bruker.verify_untouched(path, before)

    assert exp.nucleus == "19F"
    assert exp.sfo1_MHz == pytest.approx(564.27052307)
    assert exp.pulse_program == "hahnecho.nmrfam"
    assert exp.fid is not None and exp.fid.size > 0
    if (EXPNO_1903 / "pdata" / "1" / "1r").exists():
        assert exp.processed is not None
        assert exp.processed_ppm is not None
        assert exp.processed_ppm.size == exp.processed.size
    else:  # the user reorganized this dataset; raw-fid path must still work
        assert exp.processed is None


def test_metadata_conflict_detected():
    """acqus MASR (4200 Hz) disagrees with the title (35.714 kHz) in this EXPNO."""
    exp = bruker.read_expno(require(EXPNO_1903))
    assert any("MAS rate" in c for c in exp.conflicts)


def test_is_expno():
    assert bruker.is_expno(require(EXPNO_1903))
    assert not bruker.is_expno(EXPNO_1903 / "pdata")
