"""Infinite-field δiso extrapolation from the CT centre of gravity at two
fields (Sandland 2004 Eq.1 / Baasner 2014 Fig.6)."""
import numpy as np
import pytest

from larmor.qcpmg_fields import (
    FieldPoint, cq_from_slope, dcg_at_field, infinite_field_diso,
    centre_of_gravity,
)


def test_forward_inverse_roundtrip():
    I, diso, cq, eta = 1.5, -50.0, 3.0, 0.7        # 35Cl
    nu = [58.726, 81.599]                           # 14.1 T, 19.6 T
    pts = [FieldPoint(n, dcg_at_field(diso, cq, n, I, eta)) for n in nu]
    res = infinite_field_diso(pts, spin=I, eta=eta)
    assert res.delta_iso_ppm == pytest.approx(diso, abs=1e-6)
    assert res.cq_MHz == pytest.approx(cq, abs=1e-6)
    assert res.pq_MHz == pytest.approx(cq * np.sqrt(1 + eta ** 2 / 3), abs=1e-6)


def test_second_order_shift_sign_and_scale():
    # the CT second-order shift is negative and larger at lower field
    lo = dcg_at_field(0.0, 3.0, 58.726, 1.5, 0.7)
    hi = dcg_at_field(0.0, 3.0, 81.599, 1.5, 0.7)
    assert lo < hi < 0                              # both negative, low field lower
    assert lo == pytest.approx(-75.9, abs=1.0)      # realistic 35Cl magnitude


def test_uncertainties_match_paper_scale():
    I = 1.5
    pts = [FieldPoint(58.726, -125.9, 5.0), FieldPoint(81.599, -89.31, 5.0)]
    res = infinite_field_diso(pts, spin=I, eta=0.7)
    # Baasner quotes ~±16 ppm δiso, ~±0.3 MHz Cq for two fields at ±5 ppm δcg
    assert 8 < res.delta_iso_err_ppm < 20
    assert 0.15 < res.cq_err_MHz < 0.5


def test_cq_from_slope_matches_spin_factor():
    # I=5/2 (27Al) inversion is self-consistent
    for I in (1.5, 2.5, 3.5):
        cq = 4.2
        d0 = dcg_at_field(10.0, cq, 100.0, I, 0.6)
        d1 = dcg_at_field(10.0, cq, 130.0, I, 0.6)
        slope = (d1 - d0) / (1 / 130.0 ** 2 - 1 / 100.0 ** 2)
        assert cq_from_slope(slope, I, 0.6) == pytest.approx(cq, rel=1e-6)


def test_centre_of_gravity_window():
    x = np.linspace(-200, 100, 600)
    y = np.exp(-0.5 * ((x + 120) / 8) ** 2)         # peak at -120
    assert centre_of_gravity(x, y) == pytest.approx(-120, abs=0.5)
    # windowing excludes a second peak
    y2 = y + 0.4 * np.exp(-0.5 * ((x - 50) / 8) ** 2)
    assert centre_of_gravity(x, y2, -160, -80) == pytest.approx(-120, abs=1.0)


def test_dialog_computes(qapp=None):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication, QTableWidgetItem
    QApplication.instance() or QApplication([])
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    dlg = QcpmgFieldsDialog(None, "35Cl", None)
    dlg.table.setItem(0, 0, QTableWidgetItem("58.726"))
    dlg.table.setItem(0, 1, QTableWidgetItem("-125.9"))
    dlg.table.setItem(1, 0, QTableWidgetItem("81.599"))
    dlg.table.setItem(1, 1, QTableWidgetItem("-89.31"))
    dlg._compute()
    assert "iso" in dlg.result.text()
    assert dlg.spin.value() == 1.5                  # 35Cl
    dlg.close()
