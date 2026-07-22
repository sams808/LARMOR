"""NMR conversion calculators and the isotope/NMR-table backend."""
import numpy as np
import pytest

from larmor import convert as C
from larmor import nuclei as N


def test_chemical_shift():
    assert C.ppm_to_Hz(10.0, 156.0) == pytest.approx(1560.0)
    assert C.Hz_to_ppm(1560.0, 156.0) == pytest.approx(10.0)


def test_quadrupole():
    assert C.pq_from_cq_eta(5.0, 0.5) == pytest.approx(5.0 * np.sqrt(1 + 0.25 / 3))
    assert C.cq_from_pq_eta(C.pq_from_cq_eta(5.0, 0.7), 0.7) == pytest.approx(5.0)
    assert C.nu_q(5.0, 2.5) == pytest.approx(0.75)          # 3·5/(2·2.5·4)
    assert C.cq_from_nu_q(0.75, 2.5) == pytest.approx(5.0)
    # central-transition second-order shift is negative and scales as (PQ/ν0)²
    s1 = C.ct_second_order_shift_ppm(5.0, 2.5, 156.0)
    s2 = C.ct_second_order_shift_ppm(10.0, 2.5, 156.0)
    assert s1 < 0 and s2 == pytest.approx(4 * s1, rel=1e-6)


def test_dipolar_roundtrip_and_known_values():
    gH = N.GAMMA_1H
    d = C.dipolar_Hz(gH, gH, 1.5)
    assert d / 1e3 == pytest.approx(35.5, abs=0.5)          # ¹H-¹H @ 1.5 Å
    assert C.distance_from_dipolar(gH, gH, d) == pytest.approx(1.5, rel=1e-6)
    # scales as 1/r³
    assert C.dipolar_Hz(gH, gH, 3.0) == pytest.approx(d / 8.0, rel=1e-9)


def test_isotope_data_and_larmor():
    al = next(i for i in N.all_isotopes() if i.symbol == "27Al")
    assert al.spin == 2.5 and al.element == "Al"
    assert al.larmor_MHz(9.4) == pytest.approx(104.37, abs=0.05)
    assert al.larmor_from_1H(600.0) == pytest.approx(156.4, abs=0.2)
    assert N.b0_from_1H(500.0) == pytest.approx(11.74, abs=0.02)


def test_nmr_table_layout():
    # every symbol in the grid is a known element
    grid = {s for row in N.PERIODIC_ROWS for s in row if s}
    assert grid <= set(N.ELEMENT_NAME)
    assert {"H", "Al", "Cl", "Pb"} <= grid
    # the primary isotope of chlorine is the most receptive one, 35Cl
    assert N.primary_isotope("Cl").symbol == "35Cl"
    assert N.primary_isotope("H").symbol == "1H"


def test_measure_region():
    from larmor import measure as M
    x = np.linspace(-50, 50, 4000)
    g = lambda c, w, a: a * np.exp(-4 * np.log(2) * ((x - c) / w) ** 2)
    y = g(0, 6, 1.0) + g(20, 6, 0.5)
    assert M.fwhm(x, y, (10, -10)) == pytest.approx(6.0, abs=0.05)
    assert M.centre_of_mass(x, y, (10, -10)) == pytest.approx(0.0, abs=0.1)
    rows = M.integrate_regions(x, y, [(10, -10), (30, 10)])
    assert rows[0]["percent"] == pytest.approx(66.7, abs=0.5)
    assert rows[1]["percent"] == pytest.approx(33.3, abs=0.5)
