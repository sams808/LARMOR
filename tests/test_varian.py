"""Varian/Agilent (VnmrJ) reader: nucleus naming + referenced ppm axis."""
import numpy as np
import pytest

from larmor.io import varian


def test_nucleus_naming():
    assert varian.varian_nucleus("Al27") == "27Al"
    assert varian.varian_nucleus("Na23") == "23Na"
    assert varian.varian_nucleus("H1") == "1H"
    assert varian.varian_nucleus("C13") == "13C"
    assert varian.varian_nucleus("") == ""


def test_ppm_axis_reference_and_direction():
    # ref peak rfp=0 sits rfl Hz from the right edge; sw/sfrq = 100 ppm wide
    ppm = varian.varian_ppm_axis(sfrq=100.0, sw=10000.0, rfl=5000.0, rfp=0.0, n=1001)
    assert ppm[0] == pytest.approx(50.0)          # high ppm on the left
    assert ppm[-1] == pytest.approx(-50.0)        # low ppm on the right
    assert ppm[500] == pytest.approx(0.0, abs=1e-6)   # reference at the centre
    assert ppm[0] > ppm[-1]                        # descending (NMR convention)


def test_ppm_axis_offset_reference():
    # rfl closer to the right edge shifts the whole window up
    ppm = varian.varian_ppm_axis(sfrq=100.0, sw=10000.0, rfl=2500.0, rfp=0.0, n=1001)
    assert ppm.max() == pytest.approx(75.0)
    assert ppm.min() == pytest.approx(-25.0)


def test_is_varian_detects_fid_dir(tmp_path):
    d = tmp_path / "sample.fid"
    d.mkdir()
    assert not varian.is_varian(d)                 # needs procpar + fid
    (d / "procpar").write_text("x")
    (d / "fid").write_bytes(b"\0")
    assert varian.is_varian(d)
    assert varian.is_varian(d / "fid")             # a file inside works too
