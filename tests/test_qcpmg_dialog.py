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


def test_send_flushes_a_pending_debounce(qapp, tmp_path):
    """Send inside the 150 ms debounce window used to ship the STALE spectrum
    while the meta recorded the NEW settings -- provenance that lies about how
    the data was processed."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    stale = d._spec.copy()                     # spectrum at LB = 0
    d.lb.setValue(800.0)                       # queues a recompute...
    assert d._timer.isActive()                 # ...which has NOT run yet
    got = {}
    d.accepted_1d.connect(lambda p, a, m: got.update(amp=a, meta=m))
    d._send()                                  # must flush, then emit
    assert not d._timer.isActive()
    assert got["meta"]["qcpmg_lb_Hz"] == 800.0
    assert not np.array_equal(got["amp"], stale)
    d.close()


def test_failed_second_load_keeps_nothing_sendable(qapp, tmp_path, monkeypatch):
    """A second _load that reads fine but cannot be processed (e.g. an aborted
    3-point fid) used to leave the OLD spectrum sendable under the NEW file's
    name; now nothing is sendable and the state is fully reset."""
    from larmor.io import bruker

    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))
    assert d.btnSend.isEnabled()

    class _Tiny:
        domain, ndim = "time", 1
        data = np.ones(3, complex)
        meta = {"sw_Hz": 50000.0, "title": "aborted", "expno": "2"}
    monkeypatch.setattr(bruker, "read", lambda _s: _Tiny())
    d._load(str(tmp_path / "2" / "fid"))
    assert "cannot process" in d.res.text()
    assert not d.btnSend.isEnabled() and not d.btnCsv.isEnabled()
    assert d._spec is None and d._spec_raw is None
    got = {}
    d.accepted_1d.connect(lambda p, a, m: got.update(meta=m))
    d._send()
    assert not got                             # nothing emitted
    d._copy_csv()                              # must not raise either
    d.close()


def test_period_change_resets_echo_bookkeeping(qapp, tmp_path):
    """Exclusions and the echo count are indices INTO one particular split;
    after a period change they would point at different stretches of data, so
    they are reset (excluded={3} at 128 pts is not echo 3 at 64 pts)."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path))

    class _Pt:
        def data(self):
            return 3
    d._decay_clicked(None, [_Pt()])
    assert d.excluded == {3}
    d.dropFirst.setValue(1)
    d.period.setValue(64)
    assert d.excluded == set()
    assert d.nEch.value() == d.fid.size // 64
    assert d.dropFirst.value() == 0
    d.close()


def test_spinboxes_do_not_react_per_keystroke(qapp):
    """Typing '293' digit by digit must not act on the intermediate '29' --
    _on_period_pts would clamp the echo top to 28 and the recompute would
    silently run from the wrong point."""
    from PySide6.QtWidgets import QDoubleSpinBox, QSpinBox

    from larmor.desktop.qcpmg_dialog import QcpmgDialog
    d = QcpmgDialog(None, None)
    boxes = d.findChildren(QSpinBox) + d.findChildren(QDoubleSpinBox)
    assert boxes
    assert all(not b.keyboardTracking() for b in boxes)
    d.close()


def test_period_hz_field_always_shows_a_realisable_spacing(qapp, tmp_path):
    """Typing 392 Hz quantises to the same 128-point period; the field must
    snap back to the 390.625 Hz that period actually realises, not keep
    showing a spacing the data cannot have."""
    d = _dialog_with(qapp, _synthetic_qcpmg(tmp_path), period=128, sw=50000.0)
    d.periodHz.setValue(392.0)
    assert d.period.value() == 128
    # the field shows 2 decimals, so the snap lands on 390.62
    assert d.periodHz.value() == pytest.approx(390.625, abs=0.011)
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
