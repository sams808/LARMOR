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


# ---------------------------------------------------------------- fix 5
def test_shared_dialog_reconciles_the_nucleus(qapp):
    """The shared dialog kept the FIRST nucleus it was created with: a 35Cl
    send into a dialog opened for 27Al was fitted at I = 5/2. An empty
    table adopts the new nucleus; a table with rows keeps its own and warns."""
    from larmor.desktop import qcpmg_fields_dialog as qfd

    qfd._shared = None
    try:
        assert qfd.QcpmgFieldsDialog(None, "27Al").spin.value() == 2.5
        d = qfd.shared_fields_dialog(None, "27Al")
        assert d._nucleus == "27Al" and d.spin.value() == 2.5
        x = np.linspace(-300.0, 100.0, 2001)
        cur = (78.354, x, np.exp(-(((x + 105.0) / 25.0) ** 2)))
        d2 = qfd.shared_fields_dialog(None, "35Cl", cur)
        assert d2 is d and d._nucleus == "35Cl" and d.spin.value() == 1.5
        assert "35Cl" in d.lblNuc.text()
        # a dataset row is present -> a different nucleus is refused
        r = d.add_dataset_spectrum(78.354, x, cur[2], nucleus="35Cl")
        assert r >= 0
        r2 = d.add_dataset_spectrum(130.3, x, cur[2], nucleus="27Al")
        assert r2 == -1 and "27Al" in d.wresult.text()
        assert d._nucleus == "35Cl" and d.spin.value() == 1.5
        qfd.shared_fields_dialog(None, "27Al")
        assert d._nucleus == "35Cl"                       # unchanged
        # an anonymous dialog is not silently relabelled 35Cl
        e = qfd.QcpmgFieldsDialog(None, "")
        assert e._nucleus == "" and "—" in e.lblNuc.text()
        e.close(); d.close()
    finally:
        qfd._shared = None


# ---------------------------------------------------------------- fix 6
def _comb_spectrum(asymmetric=True):
    """A CT-like envelope sampled by a spikelet comb whose spacing is ~1/25
    of the band -- what a TopSpin 1r of a QCPMG EXPNO looks like."""
    x = np.linspace(-300.0, 100.0, 4001)
    if asymmetric:
        env = np.where(x < -105.0, np.exp(-(((x + 105.0) / 80.0) ** 2)),
                       np.exp(-(((x + 105.0) / 30.0) ** 2)))
    else:
        env = np.exp(-(((x + 105.0) / 60.0) ** 2))
    comb = np.zeros_like(x)
    for c in np.arange(-300.0, 100.0, 4.0):
        comb += np.exp(-((x - c) / 0.3) ** 2)
    return x, env, env * comb


def test_seed_window_detects_a_spikelet_comb_and_seeds_from_the_envelope():
    from larmor import qcpmg

    x, env, y = _comb_spectrum()
    seed = qcpmg.seed_window(x, y, {"pulse_program": "qcpmg_dec.ih"})
    assert seed.comb and seed.period_ppm == pytest.approx(4.0, abs=0.2)
    assert seed.n_spikelets >= 5 and "QCPMG spikelet" in seed.note
    # the seeded window spans the envelope, the raw first-minima one does not
    hi0, lo0 = qcpmg.cg_window(x, y)
    assert (seed.hi_ppm - seed.lo_ppm) > 5 * (hi0 - lo0)
    cg_env = qcpmg.centre_of_gravity(x, env, seed.window)[0]
    cg_comb = qcpmg.centre_of_gravity(x, y, seed.window)[0]
    cg_old = qcpmg.centre_of_gravity(x, y)[0]
    assert abs(cg_comb - cg_env) < 1.0                 # tracks the envelope
    assert abs(cg_old - cg_env) > 8.0                  # the bug this pins
    # a smooth band, with noise at any S/N, is NOT a comb
    rng = np.random.default_rng(0)
    for sn in (100.0, 20.0, 8.0):
        ys = env + rng.standard_normal(x.size) / sn
        assert not qcpmg.seed_window(x, ys).comb
        assert qcpmg.detect_comb(x, ys)[0] == 0


def test_add_dataset_spectrum_flags_a_comb_and_measures_its_envelope(qapp):
    from larmor import qcpmg
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    x, env, y = _comb_spectrum()
    d = QcpmgFieldsDialog(None, "35Cl")
    r = d.add_dataset_spectrum(78.354, x, y, magnitude=None,
                               source="Bruker 1r test/1/pdata/1")
    ds = d._ds[d._row_ds_id(r)]
    assert ds["comb"] and "spikelet" in ds["seed_note"]
    assert "spikelet" in d.table.item(r, 1).toolTip()
    assert "spikelet" in d.wresult.text()
    cg_env = qcpmg.centre_of_gravity(x, env, ds["window"])[0]
    assert float(d.table.item(r, 1).text()) == pytest.approx(cg_env, abs=1.0)
    assert "mode unknown" in d.table.item(r, 0).toolTip()      # magnitude=None
    assert "source: Bruker 1r" in d.table.item(r, 0).toolTip()
    d.close()


def test_pick_datasets_accepts_a_larmor_csv(qapp, tmp_path, monkeypatch):
    """The picker promised 'saved datasets' but only opened Bruker 1r files:
    the LAW*.csv sum-echo datasets were 'skipped'."""
    from PySide6.QtWidgets import QFileDialog, QMessageBox

    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog
    from larmor.io import spectra

    x = np.linspace(-300.0, 100.0, 3001)
    y = np.exp(-(((x + 105.0) / 25.0) ** 2))
    p = tmp_path / "LAW0Ca-3Cl_850_MHz.csv"
    spectra.write_csv(p, x, y, {"nucleus": "35Cl", "larmor_MHz": 78.354,
                                "sample": "x · QCPMG sum echo (LB 75 Hz, magnitude)"})
    monkeypatch.setattr(QFileDialog, "getOpenFileNames",
                        staticmethod(lambda *a, **k: ([str(p)], "")))
    warned = []
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: warned.append(a)))
    d = QcpmgFieldsDialog(None, "35Cl")
    n0 = d.table.rowCount()
    d._pick_datasets()
    assert not warned
    assert d.table.rowCount() == n0 + 1
    r = n0
    assert float(d.table.item(r, 0).text()) == pytest.approx(78.354)
    assert float(d.table.item(r, 1).text()) == pytest.approx(-105.0, abs=1.0)
    ds = d._ds[d._row_ds_id(r)]
    assert ds["magnitude"] is True                   # legacy ', magnitude)'
    assert "LAW0Ca-3Cl_850_MHz.csv" in ds["source"]
    d.close()


def test_real_qcpmg_1r_is_a_comb_read_as_magnitude():
    """The MagLab 35Cl EXPNO 1 1r (PH_mod = 2): its first-minima window is one
    spikelet gap (dcg -105.1) while the sum-echo CSV gives -112.8; seeded
    from the envelope the 1r tracks the CSV to ~1.3 ppm."""
    from larmor.io import spectra
    from larmor.qcpmg_fields import read_field_spectrum
    from larmor import qcpmg
    from tests.conftest import MAGLAB_35CL, require

    one_r = require(MAGLAB_35CL / "1" / "pdata" / "1" / "1r")
    csv = require(MAGLAB_35CL.parent / "LAW0Ca-3Cl_850_MHz.csv")
    fs = read_field_spectrum(str(one_r))
    assert fs["magnitude"] is True and fs["meta"].get("ph_mod") == 2
    assert fs["nucleus"] == "35Cl"
    assert fs["larmor"] == pytest.approx(78.354, abs=1e-3)
    seed = fs["seed"]
    assert seed.comb and seed.period_ppm == pytest.approx(6.8, abs=0.2)
    assert seed.lo_ppm == pytest.approx(-208.1, abs=3.0)
    assert seed.hi_ppm == pytest.approx(-38.8, abs=4.0)
    cg_1r = qcpmg.centre_of_gravity(fs["ppm"], fs["amp"], seed.window)[0]
    px, py, _ = spectra.read_csv(csv)
    cg_csv = qcpmg.centre_of_gravity(px, py)[0]
    assert cg_csv == pytest.approx(-112.76, abs=0.05)
    assert abs(cg_1r - cg_csv) < 2.0
    # and the raw first-minima window really was the 7.7 ppm bug
    assert qcpmg.centre_of_gravity(fs["ppm"], fs["amp"])[0] == pytest.approx(-105.1, abs=0.3)


# ---------------------------------------------------------------- fix 7 / 8
def _three_peaks(step_ppm=200.0, centre=-100.0, side=0.4):
    """A centreband with +-1 sidebands at +-step -- the MAS manifold of a
    narrow pattern, for the window/sideband checks."""
    x = np.linspace(-700.0, 500.0, 6001)
    y = np.exp(-(((x - centre) / 15.0) ** 2))
    for k in (-1, 1):
        y += side * np.exp(-(((x - centre - k * step_ppm) / 15.0) ** 2))
    return x, y


def test_sideband_ticks_and_the_one_sided_window_flag():
    from larmor import qcpmg

    nu, rot = 78.354, 16000.0                        # 204.2 ppm per sideband
    x, y = _three_peaks(step_ppm=rot / nu)
    peak, ticks = qcpmg.sideband_ticks(x, y, (-160.0, -40.0), rot, nu)
    assert peak == pytest.approx(-100.0, abs=0.5)
    # k = +-1, +-2 fall inside the -700..500 ppm axis, +-3 do not
    assert len(ticks) == 4 and ticks[2] == pytest.approx(-100.0 + rot / nu, abs=0.5)
    assert not qcpmg.one_sided_sideband((-160.0, -40.0), peak, rot, nu)
    # centreband + the +1 sideband, not the -1: one-sided -> flag
    m = qcpmg.measure_cg(x, y, window=(-160.0, 160.0), rotor_Hz=rot, larmor_MHz=nu)
    assert "! window catches one sideband only" in m.flags
    assert m.ticks_ppm and "wider than nu_r" in m.note
    # symmetric manifold window: no flag, and the CG is the centreband's
    m2 = qcpmg.measure_cg(x, y, window=(-360.0, 160.0), rotor_Hz=rot, larmor_MHz=nu)
    assert "! window catches one sideband only" not in m2.flags
    assert m2.cg_ppm == pytest.approx(-100.0, abs=0.5)
    # static / unknown rate: nothing to tick, nothing to flag
    m3 = qcpmg.measure_cg(x, y, window=(-160.0, 160.0))
    assert m3.ticks_ppm == [] and not any("sideband" in f for f in m3.flags)


def test_cg_convergence_sees_a_cut_tail_and_a_pedestal():
    """The jitter sigma is a local sensitivity; the drift |CG(2w) - CG(w)|
    is what a window that cuts a tail shows. A raw magnitude pedestal is
    subtracted (out-of-window median) and flagged, never mistaken for
    convergence."""
    from larmor import qcpmg

    x = np.linspace(-600.0, 300.0, 9001)
    # an asymmetric pattern with a long low-frequency tail
    y = np.where(x < -100.0, np.exp(-(((x + 100.0) / 120.0) ** 2)),
                 np.exp(-(((x + 100.0) / 30.0) ** 2)))
    true_cg = float((x * y).sum() / y.sum())
    m = qcpmg.measure_cg(x, y, window=(-180.0, -20.0))     # cuts the tail
    conv = m.convergence
    assert conv.drift_ppm > 5.0 and conv.drift_lo_ppm > conv.drift_hi_ppm
    assert "! CG not converged -- window cuts the pattern" in m.flags
    assert m.sigma_ppm == pytest.approx(max(m.jitter_ppm, m.drift_ppm))
    assert m.sigma_ppm >= m.drift_ppm > m.jitter_ppm
    assert "CG(w, 1.5w, 2w, 3w) =" in conv.sequence()
    wide = qcpmg.measure_cg(x, y, window=(-500.0, 50.0))
    assert wide.drift_ppm < 0.5 and wide.cg_ppm == pytest.approx(true_cg, abs=0.5)
    assert not any("converged" in f for f in wide.flags)
    # a positive pedestal (raw |spectrum| noise floor) on a clean line
    rng = np.random.default_rng(1)
    xg = np.linspace(-400.0, 200.0, 6001)
    yg = np.exp(-(((xg + 110.0) / 30.0) ** 2))
    ped = np.abs(yg + 0.05 * rng.standard_normal(xg.size))
    mp = qcpmg.measure_cg(xg, ped, window=(-260.0, 40.0))
    assert mp.convergence.floor_frac > 0.02
    assert any("floor" in f for f in mp.flags)
    assert mp.cg_ppm == pytest.approx(-110.0, abs=1.0)
    # ... and the signed (unclipped) CG is the one used
    ys = yg.copy(); ys[xg > -60.0] -= 0.3 * np.exp(-(((xg[xg > -60.0] + 30.0) / 20.0) ** 2))
    ms = qcpmg.measure_cg(xg, ys, window=(-260.0, 40.0))
    assert ms.cg_ppm == pytest.approx(
        qcpmg.centre_of_gravity(xg, ys, (40.0, -260.0))[0], abs=1e-9)


def test_manifold_and_centreband_modes_on_simulated_patterns():
    """The whole-manifold CG recovers a Czjzek glass's delta_iso and rms P_Q
    at the tutorial's MAS rates where the first-minima window is one-sided
    and biased; the centreband mode is accepted only for a pattern narrower
    than nu_r."""
    from larmor import qcpmg
    from larmor.convert import ct_second_order_shift_ppm, pq_from_cq_eta
    from larmor.qcpmg_fields import FieldPoint, infinite_field_diso
    from tests.conftest import simulate_ct_czjzek, simulate_ct_single

    fields = ((78.354, 16000.0), (107.811, 20000.0))
    man, mini = [], []
    for nu, rot in fields:
        x, y, rms = simulate_ct_czjzek("35Cl", nu, rot, 1.0, -70.0)
        m = qcpmg.measure_cg(x, y, mode="manifold", rotor_Hz=rot, larmor_MHz=nu)
        man.append(FieldPoint(nu, m.cg_ppm, max(m.sigma_ppm, 0.1)))
        mm = qcpmg.measure_cg(x, y, mode="minima", rotor_Hz=rot, larmor_MHz=nu)
        mini.append(mm)
        assert mm.window[1] - mm.window[0] > rot / nu          # wider than nu_r
        assert "wider than nu_r" in mm.note
    res = infinite_field_diso(man, spin=1.5, eta=0.7)
    assert res.delta_iso_ppm == pytest.approx(-70.0, abs=1.0)
    assert res.pq_MHz == pytest.approx(rms, rel=0.03)
    # the first-minima window at 78 MHz cuts the pattern: drift > 5 ppm
    assert mini[0].drift_ppm > 5.0
    assert "! CG not converged -- window cuts the pattern" in mini[0].flags
    # centreband: refused for the glass, accepted for one crystalline site
    x, y, _ = simulate_ct_czjzek("35Cl", 78.354, 16000.0, 1.0, -70.0)
    with pytest.raises(ValueError, match="centreband window refused"):
        qcpmg.measure_cg(x, y, mode="centreband", rotor_Hz=16000.0, larmor_MHz=78.354)
    x, y, _ = simulate_ct_czjzek("35Cl", 107.811, 20000.0, 1.5, -70.0)
    with pytest.raises(ValueError, match="centreband window refused"):
        qcpmg.measure_cg(x, y, mode="centreband", rotor_Hz=20000.0, larmor_MHz=107.811)
    for nu, rot in fields:
        xs, ys = simulate_ct_single("35Cl", nu, rot, 3.0, 0.7, -70.0)
        m = qcpmg.measure_cg(xs, ys, mode="centreband", rotor_Hz=rot, larmor_MHz=nu)
        expect = -70.0 + ct_second_order_shift_ppm(pq_from_cq_eta(3.0, 0.7), 1.5, nu)
        assert m.cg_ppm == pytest.approx(expect, abs=0.5)
        assert m.mode == "centreband" and m.window[1] - m.window[0] == pytest.approx(rot / nu)
    # magnitude + whole manifold is flagged as the biased combination
    mg = qcpmg.measure_cg(x, np.abs(y), mode="manifold", magnitude=True)
    assert any("magnitude + whole manifold" in f for f in mg.flags)


def test_static_czjzek_full_axis_is_converged_and_recovers_rms_pq():
    from larmor import qcpmg
    from larmor.qcpmg_fields import FieldPoint, infinite_field_diso
    from tests.conftest import simulate_ct_czjzek

    pts = []
    for nu in (78.354, 107.811):
        x, y, rms = simulate_ct_czjzek("35Cl", nu, 0.0, 1.0, -70.0)
        m = qcpmg.measure_cg(x, y, window=(float(x.min()), float(x.max())))
        assert m.drift_ppm < 0.5 and not any("converged" in f for f in m.flags)
        pts.append(FieldPoint(nu, m.cg_ppm, max(m.sigma_ppm, 0.1)))
    res = infinite_field_diso(pts, spin=1.5, eta=0.7)
    assert res.delta_iso_ppm == pytest.approx(-70.0, abs=1.0)
    assert res.pq_MHz == pytest.approx(rms, rel=0.03)


def test_field_point_from_measurement_carries_the_provenance():
    from larmor import qcpmg
    from larmor.qcpmg_fields import FieldPoint, report_text, infinite_field_diso

    x, y = _three_peaks(step_ppm=204.2)
    m = qcpmg.measure_cg(x, y, window=(-160.0, 160.0), rotor_Hz=16000.0,
                         larmor_MHz=78.354)
    p = FieldPoint.from_measurement(78.354, m, magnitude=True, source="LAW.csv",
                                    rotor_Hz=16000.0)
    assert p.window == m.window and p.flags == tuple(m.flags)
    assert p.dcg_err_ppm == pytest.approx(max(m.sigma_ppm, 0.1))
    assert "CG(w, 1.5w, 2w, 3w)" in p.cg_sequence and p.rotor_Hz == 16000.0
    q = FieldPoint(107.811, -80.0, 1.0)
    res = infinite_field_diso([p, q], 1.5, 0.7)
    txt = report_text({"s": res}, 1.5, 0.7, "35Cl")
    assert "! window catches one sideband only" in txt
    assert "CG(w, 1.5w, 2w, 3w) =" in txt


def test_real_law_series_one_sided_flag_fires_only_on_law4ca():
    """LAW4Ca at 78 MHz: the auto window (-391..+129 ppm) catches the
    centreband and ONE sideband of the ~-80 ppm centreband; LAW2/3Ca do not
    trip the check although their windows are wider than nu_r."""
    from larmor import qcpmg
    from larmor.io import spectra
    from tests.conftest import MAGLAB_35CL, require

    root = require(MAGLAB_35CL).parent
    flags = {}
    for k in (2, 3, 4):
        x, y, _ = spectra.read_csv(require(root / f"LAW{k}Ca-3Cl_850_MHz.csv"))
        m = qcpmg.measure_cg(x, y, rotor_Hz=16000.0, larmor_MHz=78.354)
        flags[k] = m.flags
        assert m.window[1] - m.window[0] > 16000.0 / 78.354
    assert "! window catches one sideband only" in flags[4]
    assert "! window sensitive" in flags[4]
    for k in (2, 3):
        assert "! window catches one sideband only" not in flags[k]
    # LAW0Ca: a genuine centreband window; the drift, not the jitter, sets sigma
    x, y, _ = spectra.read_csv(require(root / "LAW0Ca-3Cl_850_MHz.csv"))
    m = qcpmg.measure_cg(x, y, rotor_Hz=16000.0, larmor_MHz=78.354)
    assert m.cg_ppm == pytest.approx(-112.76, abs=0.05)
    assert m.jitter_ppm == pytest.approx(1.28, abs=0.05)
    assert m.drift_ppm > m.jitter_ppm and m.sigma_ppm == pytest.approx(m.drift_ppm)
    assert not any("sideband" in f for f in m.flags)


# ---------------------------------------------------------------- fix 9
def test_module_centre_of_gravity_is_the_signed_estimator():
    """qcpmg_fields.centre_of_gravity clipped negatives (a positive noise
    pedestal pulling the CG to the window centre); it is now a thin wrapper
    on the signed qcpmg.centre_of_gravity, so the two routes agree."""
    from larmor import qcpmg

    rng = np.random.default_rng(0)
    x = np.linspace(-300.0, -20.0, 2801)
    y = np.exp(-(((x + 110.0) / 30.0) ** 2)) + 0.05 * rng.standard_normal(x.size)
    assert (y < 0).any()
    a = centre_of_gravity(x, y, -300.0, -20.0)
    b = qcpmg.centre_of_gravity(x, y, (-20.0, -300.0))[0]
    assert a == pytest.approx(b, abs=0.05)
    # the clipped estimator this replaces was biased by several ppm here
    yc = np.clip(y, 0.0, None)
    clipped = float((x * yc).sum() / yc.sum())
    assert abs(clipped + 110.0) > 1.0 and abs(a + 110.0) < 1.0


def test_from_current_goes_through_the_dataset_route(qapp):
    """'dcg from open spectrum (visible range)' used to write a clipped
    whole-visible-range CG into the first empty cell with the default 5 ppm
    sigma and no window record. It now adds a dataset row: signed CG,
    data-derived sigma, FWHM, stored window, supervision band."""
    import pyqtgraph as pg
    from PySide6.QtWidgets import QWidget

    from larmor import qcpmg
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog
    from larmor.qcpmg_fields import report_text

    class _Parent(QWidget):
        def __init__(self):
            super().__init__()
            self.view = pg.PlotWidget()
            self.recipe = {"qcpmg_magnitude": False}

    rng = np.random.default_rng(3)
    x = np.linspace(-400.0, 200.0, 6001)
    y = np.exp(-(((x + 110.0) / 30.0) ** 2)) + 0.05 * rng.standard_normal(x.size)
    parent = _Parent()
    parent.view.setXRange(-300.0, 20.0, padding=0)       # wider than the line
    d = QcpmgFieldsDialog(parent, "35Cl", (78.354, x, y))
    n0 = d.table.rowCount()
    d._from_current()
    assert d.table.rowCount() == n0 + 1
    r = n0
    ds_id = d._row_ds_id(r)
    assert ds_id is not None                              # a dataset row
    ds = d._ds[ds_id]
    x0, x1 = parent.view.getPlotItem().getViewBox().viewRange()[0]
    assert ds["window"] == pytest.approx((min(x0, x1), max(x0, x1)), abs=1.0)
    assert ds["magnitude"] is False and ds["source"] == "workspace spectrum"
    cg_cell = float(d.table.item(r, 1).text())
    err_cell = d.table.item(r, 2).text()
    ref = qcpmg.centre_of_gravity(x, y, tuple(ds["window"]))
    assert cg_cell == pytest.approx(ref[0], abs=max(ref[1], 0.05))
    assert err_cell != "5" and float(err_cell) >= 0.1       # data-derived
    assert d._region is not None                             # supervision band
    # the same spectrum through add_dataset_spectrum with the same window
    r2 = d.add_dataset_spectrum(107.811, x, y, window=tuple(ds["window"]))
    assert float(d.table.item(r2, 1).text()) == pytest.approx(cg_cell, abs=1e-6)
    d._compute()
    txt = report_text(d._result_map(), 1.5, 0.7, "35Cl")
    assert "workspace spectrum" in txt or "window" in txt   # provenance travels
    # a window straddling equal +/- lobes never writes 'nan' into a cell
    y2 = np.exp(-(((x + 110.0) / 20.0) ** 2)) - np.exp(-(((x + 40.0) / 20.0) ** 2))
    d2 = QcpmgFieldsDialog(None, "35Cl", (78.354, x, y2))
    r3 = d2.add_dataset_spectrum(78.354, x, y2, window=(-200.0, 50.0))
    assert "nan" not in d2.table.item(r3, 1).text().lower()
    d.close(); d2.close(); parent.close()


# ---------------------------------------------------------------- fix 10 / 11
def test_report_prints_the_mode_and_a_declared_selectivity():
    from larmor.qcpmg_fields import mixed_modes, report_text

    pts = [FieldPoint(NU_LO, -112.76, 1.28, None, magnitude=True),
           FieldPoint(NU_HI, -95.78, 0.51, True, magnitude=True)]
    res = infinite_field_diso(pts, 1.5, 0.7)
    txt = report_text({"s": res}, 1.5, 0.7, "35Cl")
    assert "CT-selective (declared)" in txt
    rows = [ln for ln in txt.splitlines() if ln.strip().startswith(("78.", "107."))]
    assert rows[0].rstrip().endswith("?") and rows[1].rstrip().endswith("yes")
    assert "magnitude" in rows[0]
    assert "all points measured on magnitude (mc) spectra" in txt
    assert "operator's declaration" in txt
    assert not mixed_modes(pts)
    # one magnitude + one absorption: not the same observable
    pts[1].magnitude = False
    assert mixed_modes(pts)
    txt = report_text({"s": infinite_field_diso(pts, 1.5, 0.7)}, 1.5, 0.7)
    assert "NOT COMPARABLE" in txt and "absorption" in txt
    # unknown modes are neither mixed nor 'all magnitude'
    pts[1].magnitude = None
    txt = report_text({"s": infinite_field_diso(pts, 1.5, 0.7)}, 1.5, 0.7)
    assert "NOT COMPARABLE" not in txt and "all points measured" not in txt
    assert "  ?" in txt


def test_ct_selective_box_starts_unknown(qapp):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QCheckBox, QTableWidgetItem

    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog
    from larmor.qcpmg_fields import report_text

    d = QcpmgFieldsDialog(None, "35Cl")
    x = np.linspace(-300.0, 100.0, 2001)
    r = d.add_dataset_spectrum(78.354, x, np.exp(-(((x + 105.0) / 25.0) ** 2)))
    chk = d.table.cellWidget(r, 4).findChild(QCheckBox)
    assert chk.isTristate() and chk.checkState() == Qt.PartiallyChecked
    d.table.setItem(1, 0, QTableWidgetItem("107.811"))
    d.table.setItem(1, 1, QTableWidgetItem("-92.1"))
    d.table.cellWidget(1, 4).findChild(QCheckBox).setCheckState(Qt.Checked)
    d.table.setItem(0, 0, QTableWidgetItem("58.7"))
    d.table.setItem(0, 1, QTableWidgetItem("-130.0"))
    d.table.cellWidget(0, 4).findChild(QCheckBox).setCheckState(Qt.Unchecked)
    pts = d._points()
    by_nu = {round(p.larmor_MHz, 1): p.ct_selective for p in pts}
    assert by_nu[78.4] is None and by_nu[107.8] is True and by_nu[58.7] is False
    txt = report_text(d._result_map(), 1.5, 0.7, "35Cl")
    rows = {ln.split()[0]: ln.split()[-1] for ln in txt.splitlines()
            if ln.strip().startswith(("58.", "78.", "107."))}
    assert rows["58.7000"] == "no" and rows["78.3540"] == "?" and rows["107.8110"] == "yes"
    assert "(declared)" in d.table.horizontalHeaderItem(4).text()
    d.close()


# ---------------------------------------------------------------- fix 12
def test_multi_field_widths_uses_every_field():
    """The dialogs passed only the first two FWHM rows to the split with
    3+ fields loaded; the third field was silently ignored."""
    from larmor.qcpmg_fields import multi_field_widths, two_field_widths

    nus = (58.726, 78.354, 107.811)
    wq_ref, wcsd = 60.0, 20.0
    f = [np.hypot(wq_ref * (nus[0] / n) ** 2, wcsd) for n in nus]
    rng = np.random.default_rng(0)
    pts = [(n, fi * (1 + 0.02 * rng.standard_normal())) for n, fi in zip(nus, f)]
    ws = multi_field_widths(pts)
    assert ws.ok and ws.n_fields == 3 and ws.fields_MHz == pytest.approx(nus)
    assert ws.wcsd_ppm == pytest.approx(wcsd, abs=1.0)
    assert ws.wq_lo_ppm == pytest.approx(wq_ref, rel=0.05)
    assert np.isfinite(ws.chi2) and len(ws.residuals_ppm2) == 3
    # perturbing the THIRD field changes the answer (it is no longer ignored)
    pts2 = list(pts); pts2[2] = (pts2[2][0], pts2[2][1] * 1.3)
    ws2 = multi_field_widths(pts2)
    assert ws2.wcsd_ppm != pytest.approx(ws.wcsd_ppm, abs=0.5)
    # the row order of the input is irrelevant (sorted by field inside)
    assert multi_field_widths(pts[::-1]).wcsd_ppm == pytest.approx(ws.wcsd_ppm)
    # n == 2 equals the closed form exactly
    n1, n2, wq1 = 58.726, 81.599, 60.0
    f1, f2 = np.hypot(wq1, wcsd), np.hypot(wq1 * (n1 / n2) ** 2, wcsd)
    a = two_field_widths(n1, f1, n2, f2); b = multi_field_widths([(n1, f1), (n2, f2)])
    assert a.wcsd_ppm == pytest.approx(b.wcsd_ppm, abs=1e-9)
    assert a.wq_lo_ppm == pytest.approx(wq1, abs=1e-9) and a.ok
    # a negative fitted intercept: ok False, never NaN
    bad = multi_field_widths([(n1, 20.0), (n2, 60.0)])
    assert not bad.ok and np.isfinite(bad.wcsd_ppm) and np.isfinite(bad.wq_lo_ppm)
    bad3 = multi_field_widths([(58.726, 100.0), (78.354, 95.0), (107.811, 99.0)])
    assert np.isfinite(bad3.wcsd_ppm) and np.isfinite(bad3.wq_lo_ppm)
    # sigma on every point -> 1/sigma^2 weights and parameter errors
    wpts = [(n, fi, 0.5) for n, fi in pts]
    ww = multi_field_widths(wpts)
    assert ww.weighted and np.isfinite(ww.wcsd_err_ppm) and ww.wcsd_err_ppm > 0
    # a processing LB is removed in quadrature (in ppm) at each field
    lb = 300.0
    pts_lb = [(n, np.hypot(fi, lb / n)) for n, fi in zip(nus, f)]
    assert multi_field_widths(pts_lb, lb_Hz=lb).wcsd_ppm == pytest.approx(wcsd, abs=1e-6)


def test_width_split_gate_flags_the_pure_quadrupolar_regime():
    from larmor.qcpmg_fields import multi_field_widths

    n1, n2 = 78.354, 107.811
    wq = 100.0
    pure = multi_field_widths([(n1, wq), (n2, wq * (n1 / n2) ** 2 * 1.02)])
    assert pure.gate and "not resolved" in pure.gate
    small = multi_field_widths([(n1, np.hypot(wq, 15.0)), (n2, np.hypot(wq * (n1 / n2) ** 2, 15.0))])
    assert small.ok and small.wcsd_ppm == pytest.approx(15.0, abs=1e-6)
    # W_csd = 0.15 W_q: for two fields the ratio gate always covers the
    # 0.3 W_q regime (c = 0.3 a gives a 10 % ratio deviation), so either gate
    # may speak -- what matters is that one does, with the right W_csd
    assert small.gate and ("0.3 W_q" in small.gate or "not resolved" in small.gate)
    mid = multi_field_widths([(n1, np.hypot(wq, 28.0)), (n2, np.hypot(wq * (n1 / n2) ** 2, 28.0))])
    assert mid.ok and mid.gate and mid.wcsd_ppm == pytest.approx(28.0, abs=1e-6)
    fine = multi_field_widths([(n1, np.hypot(wq, 60.0)), (n2, np.hypot(wq * (n1 / n2) ** 2, 60.0))])
    assert fine.ok and fine.gate == ""


def test_report_labels_wcsd_as_shift_distribution_plus_csa():
    from larmor.qcpmg_fields import multi_field_widths, report_text

    pts = [FieldPoint(NU_LO, -113.1, 1.5), FieldPoint(NU_HI, -92.1, 1.1)]
    res = infinite_field_diso(pts, 1.5, 0.7)
    ws = multi_field_widths([(NU_LO, 59.4), (NU_HI, 43.0)])
    txt = report_text({"s": res}, 1.5, 0.7, "35Cl", {"s": ws})
    assert "W_csd" in txt and "CSA" in txt.split("W_csd", 1)[1].splitlines()[0]
    assert "width split over 2 fields" in txt


def test_second_moment_split_recovers_an_injected_shift_distribution():
    """Variances add under convolution, so the second-moment split needs no
    Gaussian assumption: a 30 ppm Gaussian shift distribution injected on a
    static CT pattern at two fields comes back within 1 ppm, where the FWHM
    split on the same pure-CT pattern returns a spurious W_csd and trips
    the gate."""
    from larmor import qcpmg
    from larmor.qcpmg_fields import multi_field_widths, second_moment_split
    from tests.conftest import simulate_ct_single

    fields = (78.354, 107.811)
    sm, fw = [], []
    for nu in fields:
        x, y = simulate_ct_single("35Cl", nu, 0.0, 3.0, 0.7, -50.0,
                                  span_ppm=(-700.0, 400.0), npts=11000,
                                  shift_fwhm_ppm=30.0)
        win = (float(x.min()), float(x.max()))
        sm.append((nu, qcpmg.second_moment_ppm(x, y, win)))
        hi, lo = qcpmg.cg_window(x, y)
        fw.append((nu, qcpmg.fwhm_hz(x, y, 1.0, (hi, lo))))
    ws = second_moment_split(sm)
    assert ws.ok and ws.kind == "second_moment"
    assert ws.wcsd_ppm == pytest.approx(30.0, abs=1.0)
    # the pure-CT FWHM split: shift_fwhm 0.3 -> a spurious W_csd, gated
    fw0 = []
    for nu in fields:
        x, y = simulate_ct_single("35Cl", nu, 0.0, 3.0, 0.7, -50.0,
                                  span_ppm=(-700.0, 400.0), npts=11000,
                                  shift_fwhm_ppm=0.3)
        hi, lo = qcpmg.cg_window(x, y)
        fw0.append((nu, qcpmg.fwhm_hz(x, y, 1.0, (hi, lo))))
    spurious = multi_field_widths(fw0)
    assert spurious.wcsd_ppm > 3.0                    # not the injected 0.3
    assert spurious.gate != ""
    sm0 = []
    for nu in fields:
        x, y = simulate_ct_single("35Cl", nu, 0.0, 3.0, 0.7, -50.0,
                                  span_ppm=(-700.0, 400.0), npts=11000,
                                  shift_fwhm_ppm=0.3)
        sm0.append((nu, qcpmg.second_moment_ppm(x, y, (float(x.min()), float(x.max())))))
    assert second_moment_split(sm0).wcsd_ppm < 3.0


def test_fields_dialog_width_split_takes_all_rows_sorted(qapp):
    from PySide6.QtWidgets import QTableWidgetItem

    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    dlg = QcpmgFieldsDialog(None, "35Cl", None)
    nus = (107.811, 58.726, 78.354)                     # unsorted on purpose
    wq_ref, wcsd = 60.0, 20.0
    for r, n in enumerate(nus):
        if r >= dlg.table.rowCount():
            dlg._add_row()
        dlg.table.setItem(r, 0, QTableWidgetItem(f"{n}"))
        dlg.table.setItem(r, 1, QTableWidgetItem("-100"))
        dlg.table.setItem(r, 3, QTableWidgetItem(
            f"{np.hypot(wq_ref * (58.726 / n) ** 2, wcsd):.4f}"))
    assert [n for n, _ in dlg._fields_fwhm()] == sorted(nus)
    dlg._compute_widths()
    assert "csd" in dlg.wresult.text().lower()          # existing assertion
    assert "3 fields" in dlg.wresult.text() and "CSA" in dlg.wresult.text()
    assert "20.0 ppm" in dlg.wresult.text()
    dlg.close()


# ---------------------------------------------------------------- fix 13
def test_three_field_misfit_is_measured_and_scales_the_errors():
    """Fields 58.79 / 78.354 / 107.811 MHz from delta_iso -70, C_Q 3.0 with
    the middle point displaced +8 ppm, all sigma 1 ppm: the a-priori +-
    (1.36 ppm) put -70 at 2.8 sigma; chi2 = 41.8 for 1 dof and the scaled
    +- (8.8 ppm) covers it."""
    from larmor.qcpmg_fields import report_text

    nus = (58.79, 78.354, 107.811)
    pts = [FieldPoint(n, dcg_at_field(-70.0, 3.0, n, 1.5, 0.7) + (8.0 if i == 1 else 0.0), 1.0)
           for i, n in enumerate(nus)]
    res = infinite_field_diso(pts, 1.5, 0.7)
    assert res.delta_iso_ppm == pytest.approx(-66.21, abs=0.01)
    assert res.delta_iso_err_ppm == pytest.approx(1.364, abs=0.01)
    assert res.chi2 == pytest.approx(41.8, rel=1e-2) and res.dof == 1
    assert res.chi2_red == pytest.approx(41.8, rel=1e-2)
    assert res.delta_iso_err_scaled_ppm == pytest.approx(8.8, abs=0.2)
    assert abs(res.delta_iso_ppm + 70.0) < res.delta_iso_err_scaled_ppm
    assert res.cq_err_scaled_MHz == pytest.approx(res.cq_err_MHz * np.sqrt(res.chi2_red))
    assert res.residuals_ppm == pytest.approx((-1.98, 5.23, -3.25), abs=0.01)
    assert res.p_value < 1e-6 and res.misfit and res.scaled
    w = 1.0 / np.array([p.dcg_err_ppm for p in pts]) ** 2
    assert abs((w * np.array(res.residuals_ppm)).sum()) < 1e-9     # weighted residuals sum to 0
    txt = report_text({"s": res}, 1.5, 0.7)
    assert "chi2/dof = 41.8/1" in txt and "p = " in txt
    assert "(a priori) / +- 8.8" in txt
    assert "resid" in txt and "5.23" in txt
    assert "improbable" in res.warning
    # perfectly linear three points: chi2 ~ 0, +- unscaled
    lin = [FieldPoint(n, dcg_at_field(-70.0, 3.0, n, 1.5, 0.7), 1.0) for n in nus]
    r2 = infinite_field_diso(lin, 1.5, 0.7)
    assert r2.chi2 < 1e-9 and not r2.scaled and not r2.misfit
    assert r2.delta_iso_err_scaled_ppm == pytest.approx(r2.delta_iso_err_ppm)
    # two points: exact, chi2 0, chi2_red NaN, said so
    r3 = infinite_field_diso(pts[1:], 1.5, 0.7)
    assert r3.chi2 == pytest.approx(0.0, abs=1e-12) and r3.dof == 0
    assert np.isnan(r3.chi2_red) and np.isnan(r3.p_value)
    assert "exact (2 points, no redundancy)" in report_text({"s": r3}, 1.5, 0.7)
    # no input errors: the scatter IS the estimate, flagged as such
    r4 = infinite_field_diso([FieldPoint(p.larmor_MHz, p.dcg_ppm, 0.0) for p in pts], 1.5, 0.7)
    assert r4.errors_from_scatter and np.isfinite(r4.delta_iso_err_ppm)
    assert r4.delta_iso_err_ppm == pytest.approx(8.8, abs=0.2)     # same scale
    assert "+- from the scatter" in report_text({"s": r4}, 1.5, 0.7)


def test_dialog_label_shows_the_chi2_line(qapp):
    from PySide6.QtWidgets import QTableWidgetItem

    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    dlg = QcpmgFieldsDialog(None, "35Cl", None)
    nus = (58.79, 78.354, 107.811)
    for r, n in enumerate(nus):
        if r >= dlg.table.rowCount():
            dlg._add_row()
        dcg = dcg_at_field(-70.0, 3.0, n, 1.5, 0.7) + (8.0 if r == 1 else 0.0)
        dlg.table.setItem(r, 0, QTableWidgetItem(f"{n}"))
        dlg.table.setItem(r, 1, QTableWidgetItem(f"{dcg:.3f}"))
        dlg.table.setItem(r, 2, QTableWidgetItem("1"))
    dlg._compute()
    assert "chi2/dof = 41.8/1" in dlg.result.text()
    assert "#c0392b" in dlg.result.text()                 # misfit colour
    assert "scaled by sqrt(chi2/dof)" in dlg.result.text()
    dlg.close()


# ---------------------------------------------------------------- fix 15
def test_report_carries_per_field_provenance_and_checks_referencing():
    from larmor.qcpmg_fields import (REF_TOLERANCE_PPM, reference_deviation_ppm,
                                     referencing_checks, report_text)

    lo = FieldPoint(NU_LO, -112.76, 1.28, magnitude=True, window=(-206.8, -35.4),
                    window_mode="minima", source="LAW0Ca-3Cl_850_MHz.csv",
                    lb_Hz=75.0, sr_hz=3982.9, referenced=True, rotor_Hz=16000.0)
    hi = FieldPoint(NU_HI, -95.78, 0.51, magnitude=True, window=(-169.0, -34.0),
                    window_mode="minima", source="LAW0Ca-3Cl_1p1GHz.csv",
                    lb_Hz=51.0, sr_hz=263.6, referenced=True, rotor_Hz=20000.0)
    res = infinite_field_diso([lo, hi], 1.5, 0.7)
    txt = report_text({"LAW0Ca": res}, 1.5, 0.7, "35Cl")
    assert "-206.8 ... -35.4 ppm" in txt and "magnitude (mc)" in txt
    assert "LB 75 Hz" in txt and "SR +3982.9 Hz" in txt and "MAS 16000 Hz" in txt
    assert "LAW0Ca-3Cl_850_MHz.csv" in txt and "LAW0Ca-3Cl_1p1GHz.csv" in txt
    assert "not checked against a 1H reference" in txt
    assert "unreferenced" not in txt
    # SR = 0 (the SF = BF1 EXPNO): a hard warning
    hi0 = FieldPoint(NU_HI, -95.78, 0.51, sr_hz=0.0)
    lines = referencing_checks([lo, hi0])
    assert any("unreferenced" in ln and "107.8110" in ln for ln in lines)
    assert not any("unreferenced" in ln and "78.3540" in ln for ln in lines)
    # +268 Hz: referenced, informational only
    assert not any("!" in ln for ln in referencing_checks([FieldPoint(NU_HI, -95.0, 0.5, sr_hz=268.0)]))
    # SF vs the session 1H reference: -0.043 ppm passes, -2.49 ppm fails
    sf_h = 1100.35
    from larmor.referencing import expected_sf_MHz
    exp = expected_sf_MHz(sf_h, "35Cl")
    good = FieldPoint(NU_HI, -95.0, 0.5, sr_hz=263.6, sf_MHz=exp * (1 - 0.043e-6))
    bad = FieldPoint(NU_HI, -95.0, 0.5, sr_hz=263.6, sf_MHz=exp * (1 - 2.49e-6))
    for pnt in (good, bad):
        pnt.ref_dev_ppm = reference_deviation_ppm(pnt.sf_MHz, "35Cl", sf_h)
    assert good.ref_dev_ppm == pytest.approx(-0.043, abs=1e-3)
    assert bad.ref_dev_ppm == pytest.approx(-2.49, abs=1e-2)
    assert abs(good.ref_dev_ppm) < REF_TOLERANCE_PPM < abs(bad.ref_dev_ppm)
    lines = referencing_checks([good, bad])
    assert not lines[0].startswith("!") and lines[1].startswith("!")
    assert "-2.49 ppm" in lines[1]


def test_fill_referencing_from_a_dataset_header_and_bruker_meta():
    from larmor import qcpmg

    x = np.linspace(-300.0, 100.0, 2001)
    y = np.exp(-(((x + 105.0) / 25.0) ** 2))
    m = qcpmg.measure_cg(x, y)
    hdr = {"lb_Hz": 75.0, "sf_MHz": 78.3621718681, "sr_hz": 3982.87,
           "referenced": True, "spectrum_mode": "magnitude(mc)"}
    p = FieldPoint.from_measurement(78.354, m, magnitude=True, meta=hdr,
                                    source="x.csv")
    assert p.lb_Hz == 75.0 and p.sr_hz == pytest.approx(3982.87)
    assert p.referenced is True and p.sf_MHz == pytest.approx(78.36217, abs=1e-5)
    # a Bruker meta with SR 0 and no 'referenced' key derives False
    q = FieldPoint.from_measurement(107.811, m, meta={"sf_MHz": 107.811292, "sr_hz": 0.0})
    assert q.referenced is False
    assert "SR +0.0 Hz" in q.provenance()


def test_figure_export_writes_a_json_record_with_the_windows(tmp_path):
    import json

    import matplotlib
    matplotlib.use("Agg")
    from larmor import figures
    from larmor.qcpmg_fields import point_provenance

    pts = [FieldPoint(NU_LO, -112.76, 1.28, window=(-206.8, -35.4), window_mode="minima",
                      magnitude=True, source="a.csv"),
           FieldPoint(NU_HI, -95.78, 0.51, window=(-169.0, -34.0), window_mode="manual",
                      magnitude=True, source="b.csv")]
    spec = {"kind": "infinite_field", "style": "article", "nucleus": "35Cl",
            "spin": 1.5, "eta": 0.7,
            "samples": [{"label": "LAW0Ca",
                         "points": [[p.larmor_MHz, p.dcg_ppm, p.dcg_err_ppm] for p in pts],
                         "provenance": [point_provenance(p) for p in pts]}]}
    written = figures.export(spec, tmp_path / "inf", formats=("png", "json"))
    assert (tmp_path / "inf.json").exists() and (tmp_path / "inf.png").exists()
    rec = json.loads((tmp_path / "inf.json").read_text(encoding="utf-8"))
    prov = rec["samples"][0]["provenance"]
    assert prov[0]["window"] == [-206.8, -35.4] and prov[1]["window_mode"] == "manual"
    assert prov[0]["source"] == "a.csv" and prov[0]["magnitude"] is True
    assert len(written) == 2
    # render still ignores the extra key
    fig = figures.render(spec)
    assert len(fig.axes) == 1
