"""Hz-or-ppm value entry in the fit-parameters table (P0 width rework)."""
import pytest
from larmor import cellparse as C


def _w(text, unit="ppm", larmor=100.0):
    return C.parse_cell(text, param_name="shift_fwhm_ppm", param_unit=unit,
                        this_index=1, n_sites=2, larmor_MHz=larmor)


def test_hz_ppm_khz_entry():
    assert _w("300Hz").value == pytest.approx(3.0)     # 300 Hz / 100 MHz
    assert _w("2ppm").value == pytest.approx(2.0)
    assert _w("1.5kHz").value == pytest.approx(15.0)
    assert _w("62.6").value == pytest.approx(62.6)      # plain number unchanged


def test_unit_on_mhz_param():
    r = C.parse_cell("3000kHz", param_name="Cq_MHz", param_unit="MHz",
                     this_index=0, n_sites=1, larmor_MHz=100.0)
    assert r.value == pytest.approx(3.0)                # 3000 kHz = 3 MHz


def test_ppm_needs_larmor():
    r = _w("300Hz", larmor=0.0)
    assert r.error and "Larmor" in r.error
