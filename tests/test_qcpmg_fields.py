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


def test_two_field_width_split_roundtrip():
    from larmor.qcpmg_fields import two_field_widths
    n1, n2, wq1, wcsd = 58.726, 81.599, 60.0, 20.0
    wq2 = wq1 * (n1 / n2) ** 2
    f1 = np.hypot(wq1, wcsd); f2 = np.hypot(wq2, wcsd)
    ws = two_field_widths(n1, f1, n2, f2)
    assert ws.ok
    assert ws.wq_lo_ppm == pytest.approx(wq1)
    assert ws.wcsd_ppm == pytest.approx(wcsd)
    # order does not matter, and an unphysical case is flagged
    assert two_field_widths(n2, f2, n1, f1).wcsd_ppm == pytest.approx(wcsd)
    assert not two_field_widths(n1, 20.0, n2, 60.0).ok


def test_dialog_computes(qapp=None):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication, QTableWidgetItem
    QApplication.instance() or QApplication([])
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    dlg = QcpmgFieldsDialog(None, "35Cl", None)
    dlg.table.setItem(0, 0, QTableWidgetItem("58.726"))
    dlg.table.setItem(0, 1, QTableWidgetItem("-125.9"))
    dlg.table.setItem(0, 3, QTableWidgetItem("63.2"))     # FWHM (ppm)
    dlg.table.setItem(1, 0, QTableWidgetItem("81.599"))
    dlg.table.setItem(1, 1, QTableWidgetItem("-89.31"))
    dlg.table.setItem(1, 3, QTableWidgetItem("37.0"))
    dlg._compute()
    dlg._compute_widths()
    assert "iso" in dlg.result.text()
    assert "csd" in dlg.wresult.text().lower()
    assert dlg.spin.value() == 1.5                  # 35Cl
    dlg.close()


@pytest.fixture(scope="module")
def qapp():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("LARMOR_NO_SESSION", "1")
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_add_dataset_spectrum_fills_row_and_supervises(qapp):
    """'Add from datasets…': δcg ± σ and FWHM are read off automatically over
    an auto window; selecting the row shows the spectrum with a draggable
    band, and moving the band recomputes the row."""
    import numpy as np

    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    d = QcpmgFieldsDialog(None, "35Cl")
    x = np.linspace(-300.0, 100.0, 4001)
    y = np.exp(-(((x + 105.0) / 25.0) ** 2))         # CT band at -105 ppm
    r = d.add_dataset_spectrum(78.354, x, y)
    assert float(d.table.item(r, 0).text()) == pytest.approx(78.354, rel=1e-6)
    assert float(d.table.item(r, 1).text()) == pytest.approx(-105.0, abs=2.0)
    fwhm_true = 2.0 * 25.0 * np.sqrt(np.log(2.0))     # ~41.6 ppm
    assert float(d.table.item(r, 3).text()) == pytest.approx(fwhm_true, rel=0.1)
    # the row was auto-selected -> supervision view with a draggable band
    assert d._region is not None and d._cg_line is not None
    d._region.setRegion((-170.0, -40.0))              # user drags the band
    ds_id = d._row_ds_id(r)
    d._region_moved(r, ds_id)
    assert d._ds[ds_id]["window"] == (-170.0, -40.0)
    assert float(d.table.item(r, 1).text()) == pytest.approx(-105.0, abs=2.0)
    # two dataset rows + Compute = the whole two-field workflow
    d.add_dataset_spectrum(160.0, x, np.exp(-(((x + 101.0) / 15.0) ** 2)))
    d._compute()
    assert "δiso" in d.result.text()
    d.close()


# ---------------------------------------------------------------- fix 1
NU_LO, NU_HI = 78.354, 107.811                     # the tutorial's 35Cl pair


def test_non_negative_slope_is_an_upper_bound_not_a_zero():
    """A slope >= 0 or within 2 sigma of zero used to print C_Q = 0.000 +-
    0.000 MHz with maximal confidence; it is a 2-sigma upper bound."""
    from larmor.qcpmg_fields import report_text

    for (d_lo, d_hi), upper in (((-90.0, -92.0), 2.331), ((-93.0, -92.0), 2.603)):
        pts = [FieldPoint(NU_LO, d_lo, 5.0), FieldPoint(NU_HI, d_hi, 5.0)]
        res = infinite_field_diso(pts, spin=1.5, eta=0.7)
        assert res.note != ""
        assert res.cq_MHz == 0.0 and res.cq_err_MHz == 0.0     # format-safe
        assert res.cq_upper_2sigma_MHz == pytest.approx(upper, abs=1e-3)
        assert np.isfinite(res.slope_err) and res.slope_err > 0
        txt = report_text({"s": res}, 1.5, 0.7, "35Cl")
        assert "0.000 +- 0.000" not in txt
        assert "<=" in txt and "upper bound" in txt
        slope_line = [ln for ln in txt.splitlines() if ln.strip().startswith("slope")][0]
        assert "+-" in slope_line
    # the positive-slope case is named as such
    res = infinite_field_diso([FieldPoint(NU_LO, -80.0, 0.5),
                               FieldPoint(NU_HI, -92.0, 0.5)], 1.5, 0.7)
    assert "POSITIVE" in res.note and res.cq_upper_2sigma_MHz == 0.0


def test_well_determined_case_is_algebraically_unchanged():
    pts = [FieldPoint(NU_LO, -113.1, 1.5), FieldPoint(NU_HI, -92.1, 1.1)]
    res = infinite_field_diso(pts, spin=1.5, eta=0.7)
    assert res.delta_iso_ppm == pytest.approx(-68.5899, abs=1e-3)
    assert res.delta_iso_err_ppm == pytest.approx(2.8733, abs=1e-3)
    assert res.cq_MHz == pytest.approx(3.0653, abs=1e-3)
    assert res.cq_err_MHz == pytest.approx(0.1358, abs=1e-3)
    assert res.note == "" and not res.cq_is_bound


# ---------------------------------------------------------------- fix 2
def test_missing_sigma_is_not_rewritten_as_a_default():
    """err = 0 used to become 1 ppm silently (identical result to err = 1);
    now a missing sigma means an unweighted fit that says so."""
    from larmor.qcpmg_fields import report_text

    base = (-112.76, -95.78)
    r0 = infinite_field_diso([FieldPoint(NU_LO, base[0], 0.0),
                              FieldPoint(NU_HI, base[1], 5.0)], 1.5, 0.7)
    r1 = infinite_field_diso([FieldPoint(NU_LO, base[0], 1.0),
                              FieldPoint(NU_HI, base[1], 5.0)], 1.5, 0.7)
    assert r1.weighted and r1.delta_iso_err_ppm == pytest.approx(10.657, abs=1e-2)
    assert not r0.weighted and np.isnan(r0.delta_iso_err_ppm)     # DIFFERENT
    # both missing: the line is exact, nothing to propagate, and it is said
    r00 = infinite_field_diso([FieldPoint(NU_LO, base[0], 0.0),
                               FieldPoint(NU_HI, base[1], 0.0)], 1.5, 0.7)
    assert r00.delta_iso_ppm == pytest.approx(-76.77, abs=0.01)
    assert np.isnan(r00.delta_iso_err_ppm) and np.isnan(r00.cq_err_MHz)
    assert "not propagated" in r00.note
    txt = report_text({"s": r00}, 1.5, 0.7)
    assert "n/a" in txt and "+- --" in txt and "0.00 +- 0.00" not in txt
    # mixed valid/missing -> unweighted
    assert not infinite_field_diso([FieldPoint(NU_LO, base[0], 0.0),
                                    FieldPoint(NU_HI, base[1], 2.0)],
                                   1.5, 0.7).weighted
    # negative or NaN sigma are "missing", never squared into a weight, and
    # never leak a NaN without a note
    for bad in (-1.5, float("nan")):
        r = infinite_field_diso([FieldPoint(NU_LO, base[0], bad),
                                 FieldPoint(NU_HI, base[1], 5.0)], 1.5, 0.7)
        assert np.isfinite(r.delta_iso_ppm) and np.isfinite(r.cq_MHz)
        assert not r.weighted and r.note
    # three points without sigma: +- from the scatter about the line
    pts = [FieldPoint(n, dcg_at_field(-70.0, 3.0, n, 1.5, 0.7) + dy, 0.0)
           for n, dy in ((58.79, 0.5), (78.354, -0.5), (107.811, 0.3))]
    r3 = infinite_field_diso(pts, 1.5, 0.7)
    assert np.isfinite(r3.delta_iso_err_ppm) and r3.delta_iso_err_ppm > 0
    assert "scatter" in r3.note


def test_static_chain_sigma_near_zero_does_not_cancel_the_denominator():
    """sigma ~1e-14 gave w ~1e28 and sw*sxx - sx*sx == 0.0 exactly: 'the two
    fields are too close' on fields 30 MHz apart."""
    res = infinite_field_diso([FieldPoint(NU_LO, -119.4, 0.0),
                               FieldPoint(NU_HI, -96.1, 1e-14)], 1.5, 0.7)
    assert np.isfinite(res.delta_iso_ppm)
    res = infinite_field_diso([FieldPoint(NU_LO, -119.4, 1e-14),
                               FieldPoint(NU_HI, -96.1, 1e-14)], 1.5, 0.7)
    assert res.weighted and np.isfinite(res.delta_iso_err_ppm)     # floored
    # the tutorial pair is untouched by the floor
    res = infinite_field_diso([FieldPoint(NU_LO, -113.1, 1.5),
                               FieldPoint(NU_HI, -92.1, 1.1)], 1.5, 0.7)
    assert res.delta_iso_err_ppm == pytest.approx(2.873, abs=1e-3)


# ---------------------------------------------------------------- fix 3
@pytest.mark.filterwarnings("error")
def test_near_identical_fields_raise_and_short_lever_arm_warns():
    for nu2 in (78.354, 78.3541):
        with pytest.raises(ValueError, match="lever arm"):
            infinite_field_diso([FieldPoint(78.354, -113.1, 1.5),
                                 FieldPoint(nu2, -113.0, 1.1)], 1.5, 0.7)
    for lo, hi in ((58.8, 83.3), (83.3, 107.8)):
        res = infinite_field_diso([FieldPoint(lo, -113.1, 1.5),
                                   FieldPoint(hi, -100.0, 1.1)], 1.5, 0.7)
        assert res.warning == "" and res.lever_arm > 0.4
    res = infinite_field_diso([FieldPoint(78.354, -113.1, 1.5),
                               FieldPoint(81.0, -110.0, 1.1)], 1.5, 0.7)
    assert res.lever_arm == pytest.approx(0.064, abs=0.002)
    assert "lever arm" in res.warning
    assert np.isfinite(res.delta_iso_err_ppm)      # no sqrt(negative) NaN


# ---------------------------------------------------------------- fix 4
def test_zero_larmor_is_refused_with_a_named_reason():
    with pytest.raises(ValueError) as ei:
        infinite_field_diso([FieldPoint(0.0, -90.0, 1.0),
                             FieldPoint(107.8, -92.0, 1.0)], 1.5, 0.7)
    msg = str(ei.value)
    assert "Larmor" in msg and "MHz" in msg and "division" not in msg


def test_implied_field_flags_the_1h_frequency_typed_for_35cl():
    from larmor.qcpmg_fields import field_plausibility_warning, implied_B0_T

    assert implied_B0_T(78.354, "35Cl") == pytest.approx(18.76, abs=0.01)
    assert implied_B0_T(78.354, "") is None
    assert implied_B0_T(78.354, "not-a-nucleus") is None
    assert field_plausibility_warning(78.354, "35Cl") == ""
    assert field_plausibility_warning(107.811, "35Cl") == ""
    w = field_plausibility_warning(850.0, "35Cl")
    assert "203.5 T" in w and "1H" in w


def test_fields_dialog_warns_on_an_implausible_field(qapp):
    from PySide6.QtWidgets import QTableWidgetItem

    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    dlg = QcpmgFieldsDialog(None, "35Cl", None)
    for r, (nu, dcg) in enumerate(((800.0, -90.0), (1100.0, -80.0))):
        dlg.table.setItem(r, 0, QTableWidgetItem(f"{nu:g}"))
        dlg.table.setItem(r, 1, QTableWidgetItem(f"{dcg:g}"))
        dlg.table.setItem(r, 2, QTableWidgetItem("1"))
    dlg._compute()
    assert " T" in dlg.result.text() and "1H" in dlg.result.text()
    for r, nu in enumerate((78.354, 107.811)):
        dlg.table.setItem(r, 0, QTableWidgetItem(f"{nu:g}"))
    dlg._compute()
    assert " T " not in dlg.result.text() and "δiso" in dlg.result.text()
    # an empty +- cell is "no sigma", never 5 ppm
    assert dlg.table.item(0, 2).text() == "1"
    dlg._add_row(58.7)
    assert dlg.table.item(dlg.table.rowCount() - 1, 2).text() == ""
    dlg.table.setItem(dlg.table.rowCount() - 1, 1, QTableWidgetItem("-120"))
    pts = dlg._points()
    assert pts[-1].dcg_err_ppm == 0.0 and not pts[-1].has_err
    dlg._compute()
    assert "not propagated" in dlg.result.text()
    dlg.close()
