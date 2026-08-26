"""Smoke tests for the guided QCPMG dialog (offscreen).

These pin the WIRING, not the physics (that is tests/test_qcpmg.py): that the
six stages exist, the interactive echo-top pick stays in sync with its
spinbox, the matched-LB button writes the field it advertises, both spectra
can be overlaid at once, and the result carries its full provenance.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _synthetic_qcpmg(tmp_path, period=128, n_echoes=24, top=64, t2_echoes=6.0):
    """A Bruker-less stand-in: the dialog only needs a readable time-domain
    source, so drive it through the same public API the loader produces."""
    t = (np.arange(period) - top)
    echo = np.exp(-np.abs(t) / 9.0) * np.exp(2j * np.pi * 0.06 * t)
    fid = np.concatenate([echo * np.exp(-k / t2_echoes) for k in range(n_echoes)])
    return fid.astype(complex)


def _dialog_with(qapp, fid, period=128, top=64, sw=50000.0):
    from larmor.desktop.qcpmg_dialog import QcpmgDialog

    d = QcpmgDialog(None, None)
    d._loading = True
    d.fid = fid
    d.meta = {"sw_Hz": sw, "larmor_MHz": 78.0, "nucleus": "35Cl",
              "title": "synthetic", "expno": "1"}
    d._carrier, d._referenced = 0.0, True
    d._period_src = "manual"
    d.period.setValue(period)
    d.periodHz.setValue(sw / period)
    d.nEch.setMaximum(fid.size // period)
    d.nEch.setValue(fid.size // period)
    d.top.setMaximum(period - 1)
    d.top.setValue(top)
    d._loading = False
    d._recompute()
    for w in (d.btnCsv, d.btnSend):
        w.setEnabled(True)
    return d


def test_dialog_constructs_with_no_source(qapp):
    from larmor.desktop.qcpmg_dialog import QcpmgDialog

    d = QcpmgDialog(None, None)
    assert d.tabs.count() == 6
    titles = [d.tabs.tabText(i) for i in range(6)]
    assert titles[0].startswith("1") and titles[-1].startswith("6")
    # nothing to send before anything is loaded
    assert not d.btnSend.isEnabled() and not d.btnCsv.isEnabled()
    d.close()


def test_echo_top_line_and_spinbox_stay_in_sync(qapp, tmp_path):
    """The interactive pick replaces ssNake's 'type the number' step -- the
    draggable line and the numeric field are one value, and neither may
    recurse into an endless recompute."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.topLine.setValue(70.0)
    d._on_top_line()
    assert d.top.value() == 70
    d.top.setValue(55)
    assert int(round(float(d.topLine.value()))) == 55
    # a drag past the end of the echo is clamped, not wrapped
    d.topLine.setValue(9999.0)
    d._on_top_line()
    assert d.top.value() == d.period.value() - 1
    d.close()


def test_matched_lb_button_writes_the_lb_field(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    assert d.t2 is not None and d.t2.ok
    assert f"{d.t2.lb_Hz:,.1f}" in d.btnMatched.text()   # the value is on the button
    d.lb.setValue(0.0)
    d._use_matched()
    assert d.lb.value() == pytest.approx(round(d.t2.lb_Hz, 1), rel=1e-6)
    d.close()


def test_matched_lb_button_disabled_when_the_fit_fails(qapp, tmp_path):
    """A failed T2 must not offer a matched filter derived from nothing."""
    d = _dialog_with(qapp, np.zeros(128 * 12, complex))
    assert d.t2 is not None and not d.t2.ok
    assert not d.btnMatched.isEnabled()
    assert "did not converge" in d.lbl3.text()
    d.close()


def test_sum_echo_and_spikelets_can_be_shown_together(qapp, tmp_path):
    """ssNake's 'Superimposed with spikelets' figure: the old dialog had a
    mutually-exclusive mode combo, so this was structurally impossible."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d.showSum.setChecked(True)
    d.showSpk.setChecked(True)
    assert len(d.p_spec.getPlotItem().listDataItems()) == 2
    assert "⊕" in d.lbl5.text()
    d.showSpk.setChecked(False)
    assert len(d.p_spec.getPlotItem().listDataItems()) == 1
    d.close()


def test_clicking_a_decay_point_excludes_it_and_refits(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    assert d.keep is not None and d.keep.all()
    before = d.t2.T2_s

    class _Pt:
        def data(self):
            return 3
    d._decay_clicked(None, [_Pt()])
    assert not d.keep[3]                       # excluded
    assert d.t2 is not None                    # and refitted
    assert d.t2.T2_s != before or d.t2.n_used < d.keep.size
    d.close()


def test_send_carries_full_provenance(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    got = {}
    d.accepted_1d.connect(lambda p, a, m: got.update(ppm=p, amp=a, meta=m))
    d._send()
    m = got["meta"]
    assert got["ppm"].size == got["amp"].size > 0
    for key in ("qcpmg_period_pts", "qcpmg_echo_top", "qcpmg_n_echoes",
                "qcpmg_lb_Hz", "qcpmg_gb_Hz", "qcpmg_p0_deg", "qcpmg_p1_deg",
                "qcpmg_carrier_ppm", "qcpmg_T2_s", "qcpmg_matched_lb_Hz"):
        assert key in m, key
    # a QCPMG sum echo is not spinning-sideband modulated: never inherit a
    # stale MAS rate from acqus
    assert m["spin_rate_Hz"] == 0.0 and m["mas_uncertain"] is False
    d.close()


def test_period_in_points_and_hz_stay_consistent(qapp, tmp_path):
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path), period=128, sw=50000.0)
    d.periodHz.setValue(250.0)                 # 50000/250 = 200 points
    assert d.period.value() == 200
    d.period.setValue(100)                     # -> 500 Hz
    assert d.periodHz.value() == pytest.approx(500.0, rel=1e-6)
    d.close()


def test_copy_csv_includes_both_time_axes(qapp, tmp_path):
    """The CSV must never say a bare 'T2' -- the ssNake pseudo-axis value and
    the physical one differ by the echo length and are both in the wild."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    d._copy_csv()
    text = QApplication.clipboard().text()
    assert "echo,t_s,intensity,kept" in text
    assert "T2_physical_s=" in text and "T2_ssnake_D1_s=" in text
    assert "LB_physical_Hz=" in text and "LB_ssnake_D1_Hz=" in text
    assert "period_pts=" in text and "echo_top=" in text
    d.close()
