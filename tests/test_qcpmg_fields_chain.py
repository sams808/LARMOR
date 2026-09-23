"""The whole QCPMG multi-field chain on SIMULATED 35Cl central-transition
patterns: window -> centre of gravity -> sigma -> extrapolation, at the
tutorial's two fields (78.354 and 107.811 MHz), static and under MAS.

tests/test_qcpmg_fields.py validates Eq. 1 against exact points and
tests/test_physics_validation.py the full-range centroid of one pattern; this
module runs cg_window / measure_cg / fit_samples on mrsimulator patterns so
window truncation, sideband and sigma behaviour are exercised end to end.
The dataflow paths that had no coverage (1r picker, report export, the two
spin sources) are here as well.
"""
import os

import pytest

from tests.conftest import simulate_ct_single

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

FIELDS = (78.354, 107.811)          # the tutorial's 35Cl pair
SPIN, ETA, DISO = 1.5, 0.7, -50.0


def _chain(cq, rotor_Hz, *, window="minima", span=(-700.0, 400.0), npts=11000,
           n_ssb=8):
    """(result, [measurement per field]) through the same functions the
    dialogs use. n_ssb 8 is quad_ct's own setting (agrees with 64 to 0.2
    ppm for C_Q <= 4.5); a C_Q 5.5 MHz manifold at 78 MHz needs more."""
    from larmor import qcpmg
    from larmor.qcpmg_fields import FieldPoint, fit_samples

    rows, meas = [], []
    for nu in FIELDS:
        x, y = simulate_ct_single("35Cl", nu, rotor_Hz, cq, ETA, DISO,
                                  span_ppm=span, npts=npts, n_ssb=n_ssb)
        if window == "full":
            m = qcpmg.measure_cg(x, y, window=(float(x.min()), float(x.max())),
                                 rotor_Hz=rotor_Hz, larmor_MHz=nu)
        else:
            m = qcpmg.measure_cg(x, y, mode=window, rotor_Hz=rotor_Hz, larmor_MHz=nu)
        meas.append(m)
        rows.append(("sim", FieldPoint.from_measurement(nu, m, rotor_Hz=rotor_Hz)))
    return fit_samples(rows, spin=SPIN, eta=ETA)["sim"], meas


@pytest.mark.parametrize("cq", (3.0, 4.5))
def test_mas_22khz_first_minima_chain_recovers_diso_and_cq(cq):
    """At 22 kHz the auto window is a genuine centreband for C_Q <= 4.5 MHz:
    delta_iso within 1 ppm, C_Q within 0.05 MHz."""
    res, meas = _chain(cq, 22000.0)
    assert abs(res.delta_iso_ppm - DISO) < 1.0, (cq, res.delta_iso_ppm)
    assert abs(res.cq_MHz - cq) < 0.05, (cq, res.cq_MHz)
    assert res.note == ""
    for m in meas:
        assert not any("sideband" in f for f in m.flags)


@pytest.mark.parametrize("rotor_Hz", (22000.0, 16000.0))
def test_mas_cq_5p5_is_flagged_not_silently_wrong(rotor_Hz):
    """C_Q 5.5 MHz: the pattern exceeds nu_r, the first-minima window is
    sensitive (sigma > 3 ppm) and, at 16 kHz, catches one sideband only."""
    res, meas = _chain(5.5, rotor_Hz)
    assert max(m.sigma_ppm for m in meas) > 3.0
    flagged = [f for m in meas for f in m.flags]
    assert flagged, "a wrong window must not pass silently"
    # the seeded window stops between the centreband and the +-1 sidebands
    # at both rates, so it is the CONVERGENCE flag (drift 20 ppm at 78 MHz,
    # 9-10 at 108 MHz) that speaks, not the one-sided-sideband one
    assert "! CG not converged -- window cuts the pattern" in flagged
    assert any(m.drift_ppm > 5.0 for m in meas)
    assert abs(res.delta_iso_ppm - DISO) > 1.0          # and it IS wrong
    # the whole manifold recovers what the window could not (all sidebands
    # simulated: n_ssb 64)
    res_m, _ = _chain(5.5, rotor_Hz, window="manifold", span=(-2500.0, 1200.0),
                      npts=16000, n_ssb=64)
    assert abs(res_m.delta_iso_ppm - DISO) < 1.5
    assert abs(res_m.cq_MHz - 5.5) < 0.1


@pytest.mark.parametrize("cq", (3.0, 4.5, 5.5))
def test_static_full_range_window_is_exact_and_sigma_zero_does_not_raise(cq):
    """Static patterns with an explicit full-range window: delta_iso within
    0.1 ppm, and the fit does not raise although the jitter sigma is ~0
    (the fitter floors it -- the old exact-zero denominator test raised
    'too close' here)."""
    res, meas = _chain(cq, 0.0, window="full", span=(-1400.0, 700.0), npts=21000)
    assert abs(res.delta_iso_ppm - DISO) < 0.1, (cq, res.delta_iso_ppm)
    assert abs(res.cq_MHz - cq) < 0.02, (cq, res.cq_MHz)
    assert all(m.sigma_ppm < 0.2 for m in meas)


def test_static_auto_window_on_a_sharp_pattern_is_unsafe_and_documented():
    """For C_Q >= 4.5 the first-minima walk on a STATIC pattern stops at the
    dip between the two horns, catching one horn (delta_iso off by tens of
    ppm) -- the manual says so."""
    from pathlib import Path

    # an acquisition-like axis (-700..400 ppm holds the whole pattern: the
    # full-window CG equals the one on a 2x wider axis)
    res, meas = _chain(5.5, 0.0)
    res_full, _ = _chain(5.5, 0.0, window="full")
    res_wide, _ = _chain(5.5, 0.0, window="full", span=(-1400.0, 700.0), npts=21000)
    assert abs(res_full.delta_iso_ppm - res_wide.delta_iso_ppm) < 0.2
    assert abs(res.delta_iso_ppm - res_full.delta_iso_ppm) > 5.0
    assert any(m.flags for m in meas)                    # not silent either
    manual = Path(__file__).resolve().parents[1] / "larmor/help/qcpmg.md"
    text = manual.read_text(encoding="utf-8")
    assert "static" in text and "horn" in text and "first-minima" in text


def test_spin_from_the_nucleus_in_both_dialogs():
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.qcpmg_batch_dialog import _spin_of
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    assert QcpmgFieldsDialog(None, "27Al").spin.value() == 2.5
    assert QcpmgFieldsDialog(None, "35Cl").spin.value() == 1.5
    assert _spin_of("27Al") == 2.5 and _spin_of("93Nb") == 4.5
    assert _spin_of("bogus") == 1.5 and _spin_of("") == 1.5


def test_pick_datasets_reads_a_real_1r_and_skips_a_fid(monkeypatch):
    """The 1r route end to end on the MagLab 35Cl EXPNO 1: Larmor from the
    meta, PH_mod 2 -> magnitude, the comb seeded from its envelope; a fid
    in the same selection lands in the skipped list."""
    from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox
    QApplication.instance() or QApplication([])
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog
    from tests.conftest import MAGLAB_35CL, require

    one_r = require(MAGLAB_35CL / "1" / "pdata" / "1" / "1r")
    fid = require(MAGLAB_35CL / "1" / "fid")
    monkeypatch.setattr(QFileDialog, "getOpenFileNames",
                        staticmethod(lambda *a, **k: ([str(one_r), str(fid)], "")))
    skipped = []
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: skipped.append(a[2])))
    d = QcpmgFieldsDialog(None, "35Cl")
    n0 = d.table.rowCount()
    d._pick_datasets()
    assert d.table.rowCount() == n0 + 1
    assert skipped and "fid" in skipped[0] and "1r" not in skipped[0].split("\n")[0].split("fid")[0]
    r = n0
    assert float(d.table.item(r, 0).text()) == pytest.approx(78.354, abs=1e-3)
    ds = d._ds[d._row_ds_id(r)]
    assert ds["magnitude"] is True and ds["comb"] is True
    assert ds["rotor_Hz"] == pytest.approx(16000.0)
    assert float(d.table.item(r, 1).text()) == pytest.approx(-114.0, abs=2.0)
    d.close()


def test_export_report_contains_every_input_and_the_assumptions(tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication, QFileDialog, QTableWidgetItem
    QApplication.instance() or QApplication([])
    from larmor.desktop.qcpmg_fields_dialog import QcpmgFieldsDialog

    d = QcpmgFieldsDialog(None, "35Cl", None)
    rows = ((78.354, -113.1, 1.5), (107.811, -92.1, 1.1))
    for r, (nu, dcg, err) in enumerate(rows):
        d.table.setItem(r, 0, QTableWidgetItem(f"{nu}"))
        d.table.setItem(r, 1, QTableWidgetItem(f"{dcg}"))
        d.table.setItem(r, 2, QTableWidgetItem(f"{err}"))
    d.eta.setValue(0.7)
    out = tmp_path / "report.txt"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(out), "")))
    d._export_report()
    txt = out.read_text(encoding="utf-8")
    for nu, dcg, err in rows:
        assert f"{nu:10.4f}" in txt and f"{dcg:11.2f}" in txt and f"{err:7.2f}" in txt
    assert "nucleus            : 35Cl" in txt
    assert "spin I             : 1.5" in txt
    assert "eta (ASSUMED)      : 0.7" in txt
    assert "P_Q" in txt and "+-" in txt.split("P_Q", 1)[1].splitlines()[0]
    assert "delta_iso" in txt and "-68.59" in txt
    d.close()
