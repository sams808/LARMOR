"""F6 -- FID <-> spectrum toggle and real / imag / |S| display channels.

Display projections of the absolute processing chain: the ProcessingPanel
owns the (domain, channel) state, the SpectrumView draws either the real
frequency-domain spectrum, another channel of the complex pipeline result, or
the windowed FID captured before the pipeline's last `ft`; the workbench
replays the RECORDED chain (not the panel's leftovers) on the first toggle.
Kept out of test_desktop.py's conflict surface. Real-data tests skip cleanly.
"""
import os

import numpy as np
import pytest

from conftest import BRUKER_1R, BRUKER_FID, require

pyside = pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QPointF, Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")  # never inherit a real session
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _write_csv_spectrum(path, centres=((0.0, 3.0, 1.0), (15.0, 5.0, 0.45))):
    """A tiny CSV spectrum with the metadata header io/spectra reads."""
    ppm = np.linspace(-60.0, 100.0, 512)
    amp = np.zeros_like(ppm)
    for c, w, a in centres:
        amp += a * np.exp(-4 * np.log(2) * ((ppm - c) / w) ** 2)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# LARMOR spectrum\n# nucleus=11B\n"
                "# larmor_MHz=160.46\n# spin_rate_Hz=20000.0\n")
        f.write("ppm,intensity\n")
        for x, y in zip(ppm, amp):
            f.write(f"{x},{y}\n")
    return ppm, amp


def _flush(qapp, win):
    """Fire the panel's live debounce now and drain the event queue."""
    win.proc_panel._live_timer.stop()
    win.proc_panel._live_timer.timeout.emit()
    qapp.processEvents()


def _names(ops):
    return [o["op"] for o in ops]


# ============================================================ ProcessingPanel
def test_panel_view_state_isolation_and_reapodize_chain(qapp):
    from larmor.desktop.panels import ProcessingPanel

    p = ProcessingPanel()
    views = []
    p.view_changed.connect(lambda d, c: views.append((d, c)))
    assert p.view_state() == ("freq", "real")

    # the channel radios must NOT share Qt's implicit sibling group with the
    # source radios: checking 'imag' would otherwise flip the pipeline to raw
    p.rb_imag.setChecked(True)
    assert p.rb_pdata.isChecked() and not p.rb_raw.isChecked()
    assert views == [("freq", "imag")]
    p.btnDomain.setChecked(True)
    assert views[-1] == ("time", "imag")
    n = len(views)
    p.reset_view()
    assert p.view_state() == ("freq", "real") and len(views) == n   # silent

    assert not p.chain_has_ft()
    p.arm_reapodize()
    assert p.chkReapod.isChecked()
    assert p.chkHilbert.isChecked() and not p.chkHilbert.isEnabled()
    assert p.adv_toggle.isChecked() and not p._adv.isHidden()
    assert not p._live_timer.isActive()
    assert p.chain_has_ft()

    got = []
    p.apply_requested.connect(lambda ops, raw: got.append((ops, raw)))
    p.btnApply.click()
    ops, raw = got[-1]
    names = _names(ops)
    assert names[:2] == ["hilbert", "ift"]
    assert names.count("em") == 1 and names.count("hilbert") == 1
    em = next(o for o in ops if o["op"] == "em")
    assert em["lb_hz"] == 50                                 # panel default
    ft = next(o for o in ops if o["op"] == "ft")
    assert "offset_ppm" not in ft
    assert "fcor" not in names
    assert names.index("ft") > names.index("em")
    assert raw is False

    # unticking restores the user's own Hilbert choice (it was unchecked)
    p.chkReapod.setChecked(False)
    assert p.chkHilbert.isEnabled() and not p.chkHilbert.isChecked()

    # raw mode: re-apodize is irrelevant, the raw chain is untouched
    p.rb_raw.setChecked(True)
    assert not p.chkReapod.isEnabled()
    p.btnApply.click()
    ops, raw = got[-1]
    assert _names(ops) == ["fcor", "em", "zf", "ft"] and raw is True
    assert ops[-1] == {"op": "ft", "offset_ppm": 0.0}


def test_panel_sync_from_ops_round_trips_and_carries_unknown_ops(qapp):
    from larmor.desktop.panels import ProcessingPanel

    p = ProcessingPanel()
    got = []
    p.apply_requested.connect(lambda ops, raw: got.append((ops, raw)))

    raw_chain = [{"op": "tdeff", "points": 4096}, {"op": "fcor", "factor": 0.5},
                 {"op": "gm", "lb_hz": -30.0, "gb": 0.2}, {"op": "zf", "factor": 4},
                 {"op": "ft", "offset_ppm": 1.5},
                 {"op": "phase", "p0": 33.0, "p1": -10.0, "pivot_frac": 0.3},
                 {"op": "magnitude"}, {"op": "sr", "sr_hz": 120.0}]
    assert p.sync_from_ops(raw_chain, True) is True
    assert p.rb_raw.isChecked() and p.wdw.currentText() == "GM"
    assert p.lb.value() == 30 and p.gb.value() == pytest.approx(0.2)
    assert p.tdeff.value() == 4096 and p.zf.value() == 4
    assert p.off.value() == 1.5 and p.p0v.value() == 33 and p.p1v.value() == -10
    assert p.chkMag.isChecked() and p.sr.value() == 120
    assert p.adv_toggle.isChecked()
    assert got == [] and not p._live_timer.isActive()        # silent
    p._emit([])
    ops, raw = got[-1]
    expect = [dict(o) for o in raw_chain]
    expect[5].pop("pivot_frac")                     # the live pivot governs
    assert raw is True and ops == expect

    # a re-apodized pdata chain with steps the widgets cannot express
    pdata_chain = [{"op": "hilbert"}, {"op": "ift"}, {"op": "em", "lb_hz": 80.0},
                   {"op": "ft"}, {"op": "phase", "p0": 12.0, "p1": 0.0},
                   {"op": "baseline", "order": 3}, {"op": "subtract_avg"}]
    assert p.sync_from_ops(pdata_chain, False) is True
    assert p.rb_pdata.isChecked() and p.chkReapod.isChecked()
    assert p.chkHilbert.isChecked() and not p.chkHilbert.isEnabled()
    assert p.wdw.currentText() == "EM" and p.lb.value() == 80 and p.zf.value() == 1
    assert p._carried_ops == [{"op": "baseline", "order": 3}, {"op": "subtract_avg"}]
    p._emit([])
    ops, raw = got[-1]
    assert raw is False and ops == pdata_chain

    # an empty chain is all-neutral and emits nothing
    assert p.sync_from_ops([], False) is True
    assert not p.chkReapod.isChecked() and p.chkHilbert.isEnabled()
    assert not p.chkHilbert.isChecked() and p.p0v.value() == 0
    p._emit([])
    assert got[-1] == ([], False)

    # a chain the panel cannot drive is refused without touching a widget
    p.lb.setValue(77.0)
    bad = [{"op": "swap_echo", "point": 147}, {"op": "echo_apodize", "lb_hz": 20},
           {"op": "ft"}]
    assert p.sync_from_ops(bad, True) is False
    assert p.lb.value() == 77.0 and p.rb_pdata.isChecked()
    assert p.sync_from_ops([{"op": "zf", "si": 8192}, {"op": "ft"}], True) is False
    assert p.sync_from_ops([{"op": "em", "lb_hz": 5}], True) is False  # no ft

    # Reset clears the carried steps
    p.sync_from_ops(pdata_chain, False)
    assert p._carried_ops
    p.btnReset.click()
    assert p._carried_ops == []


# ================================================================ SpectrumView
def _two_lines():
    x = np.linspace(-20.0, 40.0, 400)
    y = (np.exp(-4 * np.log(2) * ((x - 10.0) / 2.0) ** 2)
         + 0.3 * np.exp(-4 * np.log(2) * ((x - 25.0) / 3.0) ** 2))
    return x, y


def _fid_trace():
    t = np.linspace(0.0, 50.0, 256)
    return t, np.exp(-t / 10.0)


def test_spectrum_view_fid_mode_hides_ppm_items_and_restores(qapp):
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    fired = []
    v.experiment_set.connect(lambda: fired.append(1))
    x, y = _two_lines()
    v.set_experiment(x, y)
    assert fired == [1]
    v.show_phase_pivot(True)
    v.set_model(x, y, [y], ["a"], set(), x, y)
    v.set_paddles([(0, 10.0, 1.0, 2.0, True)])
    v.setXRange(0, 30, padding=0)
    frac0 = v.phase_pivot_frac()
    assert v.domain == "freq"

    t, f = _fid_trace()
    v.set_fid(t, f, "FID (real)")
    vb = v.getPlotItem().getViewBox()
    assert v.domain == "time" and not vb.xInverted()
    assert v.getPlotItem().getAxis("bottom")._factor == 1.0
    for it in (v._pivot, v._model, v._components[0], v._paddles[0]):
        assert it.isVisible() is False
    assert v.phase_pivot_frac() == frac0             # frequency-based pivot
    assert v.current_xrange() == (30.0, 0.0)         # the saved ppm window
    assert np.array_equal(v._exp.xData, t)
    assert v._legend.getLabel(v._exp).text == "FID (real)"
    assert "time" in v.getPlotItem().getAxis("bottom").labelText

    # items created WHILE the FID is displayed must not appear on the ms axis
    v.set_model(x, y, [y, y], ["a", "b"], set(), x, y)
    v.set_paddles([(0, 10.0, 1.0, 2.0, True), (1, 25.0, 0.3, 3.0, True)])
    v.show_paddles(True)
    v.set_zones([[30.0, 0.0]])
    assert all(not c.isVisible() for c in v._components)
    assert all(not p.isVisible() for p in v._paddles)
    assert not v._zones[0].isVisible()

    v.set_experiment(x, y)
    assert v.domain == "freq" and vb.xInverted()
    assert v._pivot.isVisible() and v._model.isVisible()
    assert all(c.isVisible() for c in v._components)
    assert all(p.isVisible() for p in v._paddles) and v._zones[0].isVisible()
    xr = vb.viewRange()[0]
    assert (min(xr), max(xr)) == pytest.approx((0.0, 30.0))
    assert "chemical shift" in v.getPlotItem().getAxis("bottom").labelText
    assert v._legend.getLabel(v._exp).text == "experiment"
    assert fired == [1, 1]


def test_spectrum_view_channel_trace_keeps_the_real_trace_for_the_pivot(qapp):
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    x, y = _two_lines()
    v.set_experiment(x, y)
    fired = []
    v.experiment_set.connect(lambda: fired.append(1))
    v.set_channel_trace(-y, "experiment (imag)")
    assert np.array_equal(v._exp.yData, -y)
    assert v._legend.getLabel(v._exp).text == "experiment (imag)"
    assert fired == []                                # not a new spectrum
    v.show_phase_pivot(True)
    top = float(x[int(np.argmax(y))])
    # the tallest REAL peak, not the displayed channel's maximum
    assert v._pivot.value() == pytest.approx(top)
    assert v._snap_peak(top + 0.4) == pytest.approx(top)


class _FakeClick:
    def __init__(self, button, pos):
        self._b, self._p = button, pos

    def button(self):
        return self._b

    def scenePos(self):
        return self._p

    def accept(self):
        pass


def test_spectrum_view_ignores_left_clicks_in_fid_mode(qapp):
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    v.resize(600, 400)
    v.show()
    qapp.processEvents()
    x, y = _two_lines()
    v.set_experiment(x, y)
    v.set_add_mode("gauss_lor")
    got = []
    v.add_requested.connect(lambda px, py: got.append((px, py)))
    rect = v.sceneBoundingRect()
    assert not rect.isEmpty()
    pos = QPointF(rect.center())
    t, f = _fid_trace()
    v.set_fid(t, f, "FID (real)")
    v._on_click(_FakeClick(Qt.LeftButton, pos))
    assert got == []                                  # inert on the ms axis
    v.set_experiment(x, y)
    v._on_click(_FakeClick(Qt.LeftButton, pos))
    assert len(got) == 1
    v.close()


# ================================================================= MainWindow
def _load_csv(win, qapp, tmp_path, name="s.csv", **kw):
    path = tmp_path / name
    ppm, amp = _write_csv_spectrum(path, **kw)
    win.load_source(str(path), keep_fit=False)
    qapp.processEvents()
    assert win.central_stack.currentWidget() is win.view
    assert win.exp_ppm.size == 512
    return path, ppm, amp


def _fwhm_at(x, y, x0):
    """Half-height width of the line whose maximum is nearest x0, the two
    crossings linearly interpolated (the grid step is 0.31 ppm)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    i = int(np.argmin(np.abs(x - x0)))
    while 0 < i < x.size - 1 and (y[i + 1] > y[i] or y[i - 1] > y[i]):
        i += 1 if y[i + 1] > y[i] else -1
    half = y[i] / 2.0
    li = i
    while li > 0 and y[li - 1] > half:
        li -= 1
    ri = i
    while ri < x.size - 1 and y[ri + 1] > half:
        ri += 1

    def cross(a, b):
        if y[a] == y[b]:
            return x[a]
        return x[a] + (half - y[a]) * (x[b] - x[a]) / (y[b] - y[a])

    left = cross(li, li - 1) if li > 0 else x[li]
    right = cross(ri, ri + 1) if ri < x.size - 1 else x[ri]
    return abs(right - left)


def _legend_text(win):
    return win.view._legend.getLabel(win.view._exp).text


def test_imag_channel_shows_the_imaginary_part_and_keeps_the_fit_data_real(
        win, qapp, tmp_path):
    _load_csv(win, qapp, tmp_path)
    amp0 = win.exp_amp.copy()
    pp = win.proc_panel
    pp.rb_imag.setChecked(True)
    qapp.processEvents()
    assert pp.chkHilbert.isChecked()          # armed: a 1r/CSV has no imag channel
    assert win._proc_spec_valid()
    assert np.any(win._proc_spec.y.imag)
    assert np.allclose(win.view._exp.yData, win._proc_spec.y.imag)
    # op_hilbert keeps the real part (to FFT round-off): the fit data is unchanged
    assert np.allclose(win.exp_amp, amp0, atol=1e-12)
    assert _legend_text(win) == "experiment (imag)"
    assert win.actChannel["imag"].isChecked()
    assert win.recipe["processing"] == [{"op": "hilbert"}]
    assert win.recipe["processing_from_raw"] is False
    assert "imaginary" in win.statusBar().currentMessage()

    win._cycle_channel()
    qapp.processEvents()
    assert pp.view_state() == ("freq", "magnitude")
    assert _legend_text(win) == "|experiment|"
    assert np.allclose(win.view._exp.yData, np.abs(win._proc_spec.y))
    assert win.actChannel["magnitude"].isChecked()

    win._set_channel("real")
    qapp.processEvents()
    assert np.array_equal(win.view._exp.yData, win.exp_amp)
    assert _legend_text(win) == "experiment"
    assert win.actChannel["real"].isChecked()


def test_fid_toggle_on_pdata_arms_reapodize_shows_a_one_sided_fid_and_returns(
        win, qapp, tmp_path):
    from larmor import processing as proc

    _load_csv(win, qapp, tmp_path)
    x0 = win.exp_ppm.copy()
    # one line, so run_fit gets past its 'add a line first' check
    name = next(iter(win._model_actions))
    win._set_add_mode(name)
    win.add_site_at(0.0, 1.0)
    win._set_add_mode(None)
    qapp.processEvents()
    pp = win.proc_panel

    win._toggle_time_domain(True)
    qapp.processEvents()
    assert not win.proc_dock.isHidden()
    assert pp.chkReapod.isChecked() and pp.adv_toggle.isChecked()
    assert pp.chkHilbert.isChecked() and not pp.chkHilbert.isEnabled()
    assert pp.view_state() == ("time", "real")
    assert win.view.domain == "time"
    assert win.actTimeDomain.isChecked() and win.sbFid.isChecked()
    names = _names(win.recipe["processing"])
    assert names[:2] == ["hilbert", "ift"] and "ft" in names
    assert win.recipe["processing_from_raw"] is False
    fid = win._proc_fid
    assert fid is not None and fid.domain == "time" and fid.sw_Hz > 0
    n = fid.y.size
    # one-sided: the Hilbert-reconstructed fid decays, it is not hermitian
    assert np.abs(fid.y[-n // 8:]).sum() < 0.05 * np.abs(fid.y[:n // 8]).sum()
    assert np.allclose(win.exp_ppm, x0)              # the axis op_ft restored
    assert np.allclose(win.view._exp.xData, proc.time_axis_s(fid) * 1e3)
    assert _legend_text(win) == "FID (real)"
    assert win._format_x(1.5).endswith("ms")
    assert "time domain" in win.statusBar().currentMessage()
    win.run_fit()
    assert win._fit_worker is None
    assert "Ctrl+T" in win.statusBar().currentMessage()

    win._toggle_time_domain(False)
    qapp.processEvents()
    assert win.view.domain == "freq"
    assert pp.view_state() == ("freq", "real")
    assert np.allclose(win.view._exp.yData, win.exp_amp)
    assert _legend_text(win) == "experiment"
    assert win._format_x(100.0) == "100.00 ppm"
    assert not win.actTimeDomain.isChecked() and not win.sbFid.isChecked()


def test_changing_lb_in_fid_view_replays_the_chain_live(win, qapp, tmp_path):
    _load_csv(win, qapp, tmp_path)
    w0 = _fwhm_at(win.exp_ppm, win.exp_amp, 0.0)
    assert w0 == pytest.approx(3.0, abs=0.15)
    win._toggle_time_domain(True)
    qapp.processEvents()
    pp = win.proc_panel
    assert pp.wdw.currentText() == "EM" and pp.lb.value() == 0   # LB is live at once
    n = win._proc_fid.y.size
    tail0 = np.abs(win._proc_fid.y[n // 4:]).sum()

    pp.lb.setValue(400)
    _flush(qapp, win)
    assert win.view.domain == "time"                 # still the FID, redrawn
    assert np.abs(win._proc_fid.y[n // 4:]).sum() < tail0
    em = next(o for o in win.recipe["processing"] if o["op"] == "em")
    assert em["lb_hz"] == 400
    assert win.exp_amp.size == 512
    assert np.allclose(win.view._exp.yData, win._proc_fid.y.real)

    win._toggle_time_domain(False)
    qapp.processEvents()
    w1 = _fwhm_at(win.exp_ppm, win.exp_amp, 0.0)
    # a 3 ppm Gaussian convolved with the 400 Hz / 160.46 MHz Lorentzian EM
    # adds: Voigt FWHM (Olivero & Longbothum) 0.5346 fL + sqrt(0.2166 fL^2 + fG^2)
    f_l = 400.0 / 160.46
    expect = 0.5346 * f_l + np.sqrt(0.2166 * f_l ** 2 + w0 ** 2)
    assert w1 == pytest.approx(expect, rel=0.1)
    assert w1 > w0 + 1.0


def test_reapodized_recipe_reopens_to_the_same_spectrum(win, qapp, tmp_path):
    """The reproducibility contract for a re-apodized pdata spectrum."""
    from larmor import loader
    from larmor.recipe import Recipe

    csv, _, _ = _load_csv(win, qapp, tmp_path)
    win._toggle_time_domain(True)
    qapp.processEvents()
    win.proc_panel.lb.setValue(400)
    _flush(qapp, win)
    win._toggle_time_domain(False)
    qapp.processEvents()
    assert _names(win.recipe["processing"]) == ["hilbert", "ift", "em", "ft"]
    r = Recipe.from_dict(win.recipe)
    r.source_path = str(csv)
    path = tmp_path / "x.recipe.json"
    r.save(path)
    ppm, amp, rd, meta, warnings = loader.load_any(path)
    assert np.allclose(ppm, win.exp_ppm)
    assert np.allclose(amp, win.exp_amp, atol=1e-9 * win.exp_amp.max())
    assert any("replayed" in w for w in warnings)
    assert rd["processing"] == win.recipe["processing"]


def test_reopened_recipe_seeds_the_unprocessed_base_and_the_toggle_does_not_compound(
        win, qapp, tmp_path):
    from larmor.recipe import Recipe

    csv = tmp_path / "raw.csv"
    _, raw_amp = _write_csv_spectrum(csv)
    r = Recipe(nucleus="11B", larmor_frequency_MHz=160.46, source_kind="csv",
               source_path=str(csv),
               processing=[{"op": "hilbert"},
                           {"op": "phase", "p0": 40.0, "p1": 0.0}])
    path = tmp_path / "phased.recipe.json"
    r.save(path)
    win.load_source(str(path), keep_fit=False)
    qapp.processEvents()
    assert win._proc_base is not None
    assert np.allclose(win._proc_base[1], raw_amp)          # unprocessed base
    assert not np.allclose(win.exp_amp, raw_amp)            # shown phased by 40°
    amp_replayed = win.exp_amp.copy()

    pp = win.proc_panel
    pp.rb_imag.setChecked(True)                             # forces the synced Apply
    qapp.processEvents()
    assert pp.p0v.value() == 40 and pp.chkHilbert.isChecked()
    assert pp.view_state() == ("freq", "imag")
    assert np.allclose(win.exp_amp, amp_replayed,
                       atol=1e-9 * np.abs(amp_replayed).max())   # NOT phased twice
    assert _names(win.recipe["processing"]) == ["hilbert", "phase"]
    assert win.recipe["processing"][1]["p0"] == 40.0


def test_new_data_and_bypassing_edits_reset_the_display_to_real(win, qapp, tmp_path):
    _load_csv(win, qapp, tmp_path, name="a.csv")
    pp = win.proc_panel
    pp.rb_imag.setChecked(True)
    qapp.processEvents()
    assert pp.view_state() == ("freq", "imag") and win._proc_spec_valid()

    # new data: the display returns to the real spectrum, controls mirrored
    _load_csv(win, qapp, tmp_path, name="b.csv", centres=((-20.0, 4.0, 1.0),))
    assert pp.view_state() == ("freq", "real")
    assert win.view.domain == "freq" and win.actChannel["real"].isChecked()
    assert win._proc_spec is None
    assert _legend_text(win) == "experiment"

    pp.rb_imag.setChecked(True)
    qapp.processEvents()
    assert pp.view_state() == ("freq", "imag") and win._proc_spec_valid()
    # a bypassing mutation of the exp arrays, as apply_twopoint_bg does
    win.exp_amp = win.exp_amp * 2.0
    win.view.set_experiment(win.exp_ppm, win.exp_amp)
    assert not win._proc_spec_valid()
    assert pp.view_state() == ("freq", "real")          # at once, via experiment_set
    assert win.actChannel["real"].isChecked() and not win.actTimeDomain.isChecked()

    pp.rb_mag.setChecked(True)                          # the guard rebuilt via Apply
    qapp.processEvents()
    assert pp.view_state() == ("freq", "magnitude") and win._proc_spec_valid()
    assert np.allclose(win.view._exp.yData, np.abs(win._proc_spec.y))


def test_pivot_fraction_is_frequency_based_while_the_fid_is_shown(win, qapp, tmp_path):
    _load_csv(win, qapp, tmp_path)
    win.view.show_phase_pivot(True)
    win.view._pivot.setValue(15.0)
    f0 = win.view.phase_pivot_frac()
    assert 0.0 < f0 < 1.0 and f0 != 0.5
    win._toggle_time_domain(True)
    qapp.processEvents()
    assert win.view.domain == "time"
    assert win.view.phase_pivot_frac() == f0
    assert win.view._pivot.isVisible() is False
    win.proc_panel.p1v.setValue(90)
    _flush(qapp, win)
    assert win.view.domain == "time"
    ph = next(o for o in win.recipe["processing"] if o["op"] == "phase")
    assert ph["p1"] == 90 and ph["pivot_frac"] == pytest.approx(f0)
    win._toggle_time_domain(False)
    qapp.processEvents()
    assert win.view._pivot.isVisible() and win.view._pivot.value() == 15.0


def test_process_menu_and_sidebar_expose_toggle_and_channels(win, qapp, tmp_path):
    from PySide6.QtGui import QKeySequence
    from PySide6.QtWidgets import QMenu

    menus = win.menuBar().findChildren(QMenu)
    acts = [a for m in menus for a in m.actions()]

    def plain(a):
        return a.text().replace("&", "")       # drop the mnemonic marker

    fid = next(a for a in acts if "FID ⇄ spectrum" in plain(a))
    assert fid is win.actTimeDomain
    assert fid.shortcut() == QKeySequence("Ctrl+T") and fid.isCheckable()
    chan = next(m for m in menus if m.title() == "Display &channel")
    checkable = [a for a in chan.actions() if a.isCheckable()]
    assert len(checkable) == 3
    assert [a.isChecked() for a in checkable] == [True, False, False]
    assert checkable[0].actionGroup() is not None
    assert checkable[0].actionGroup().isExclusive()
    cycle = next(a for a in chan.actions() if "Cycle channel" in plain(a))
    assert cycle.shortcut() == QKeySequence("Ctrl+I")
    assert win.sbFid.isCheckable() and not win.sbFid.isChecked()

    # with NO data loaded everything stays at (freq, real) and unchecked
    win.actTimeDomain.trigger()
    qapp.processEvents()
    assert win.proc_panel.view_state() == ("freq", "real")
    assert not win.actTimeDomain.isChecked() and not win.sbFid.isChecked()
    win.actChannel["imag"].trigger()
    qapp.processEvents()
    assert win.proc_panel.view_state() == ("freq", "real")
    assert win.actChannel["real"].isChecked()

    _load_csv(win, qapp, tmp_path)
    win.actChannel["imag"].trigger()
    qapp.processEvents()
    assert win.proc_panel.rb_imag.isChecked()
    assert _legend_text(win) == "experiment (imag)"
    win.sbFid.trigger()
    qapp.processEvents()
    assert win.view.domain == "time"
    assert win.actTimeDomain.isChecked() and win.sbFid.isChecked()
    assert win.proc_panel.view_state() == ("time", "imag")
    assert _legend_text(win) == "FID (imag)"
    win.sbFid.trigger()
    qapp.processEvents()
    assert win.view.domain == "freq" and not win.actTimeDomain.isChecked()
    assert win.proc_panel.view_state() == ("freq", "imag")
    assert _legend_text(win) == "experiment (imag)"


# ================================================================== real data
def test_raw_fid_pipeline_captures_the_windowed_fid_and_resolves_a_fid_file_path(
        win, qapp):
    win.load_source(str(require(BRUKER_1R)), keep_fit=False)
    qapp.processEvents()
    assert win.central_stack.currentWidget() is win.view
    pp = win.proc_panel
    pp.rb_raw.setChecked(True)
    pp._live_timer.stop()
    pp.btnApply.click()
    qapp.processEvents()
    assert "processing failed" not in win.statusBar().currentMessage()
    fid = win._proc_fid
    assert fid is not None and fid.domain == "time" and fid.sw_Hz > 0
    n = fid.y.size
    assert n & (n - 1) == 0                           # zero-filled to 2^k (ZF 2)
    assert win.recipe["processing_from_raw"] is True

    win._toggle_time_domain(True)
    qapp.processEvents()
    assert win.view.domain == "time"
    assert win.view._exp.xData[-1] == pytest.approx((n - 1) / fid.sw_Hz * 1e3)

    # a dropped `fid` FILE leaves source_path at the file (the preview does):
    # raw processing must resolve it to the EXPNO instead of refusing
    win.source_path = str(require(BRUKER_FID))
    win._proc_fid = None
    pp.btnApply.click()
    qapp.processEvents()
    assert "processing failed" not in win.statusBar().currentMessage()
    assert win._proc_fid is not None and win.view.domain == "time"


def test_fid_preview_records_its_chain_and_toggles_to_the_true_fid(win, qapp):
    win.load_source(str(require(BRUKER_FID)), keep_fit=False)
    qapp.processEvents()
    assert win.central_stack.currentWidget() is win.view
    assert _names(win.recipe["processing"]) == ["fcor", "em", "ft", "magnitude"]
    assert win.recipe["processing_from_raw"] is True
    amp_preview = win.exp_amp.copy()

    pp = win.proc_panel
    win._toggle_time_domain(True)
    qapp.processEvents()
    assert pp.rb_raw.isChecked() and pp.wdw.currentText() == "EM"
    assert pp.lb.value() == 100 and pp.chkMag.isChecked()
    assert pp.zf.value() == 1 and pp.fcor.value() == 0.5
    assert win.view.domain == "time"
    assert win._proc_fid.y.size == win.exp_amp.size          # no zero-filling
    assert np.allclose(win.exp_amp, amp_preview)   # the synced replay IS the preview
    assert win.recipe["processing"] == [
        {"op": "fcor", "factor": 0.5}, {"op": "em", "lb_hz": 100},
        {"op": "ft", "offset_ppm": 0.0}, {"op": "magnitude"}]


def test_fid_dialog_records_its_processing_chain(qapp):
    from larmor import loader
    from larmor.desktop.fid_dialog import FidDialog
    from larmor.recipe import Recipe

    dlg = FidDialog(None, str(require(BRUKER_FID)))
    assert dlg.data is not None and dlg.data.ndim == 1
    got = []
    dlg.accepted_1d.connect(lambda ppm, amp, meta: got.append((ppm, amp, meta)))
    dlg.wdw.setCurrentText("GM")
    dlg.lb.setValue(30)
    dlg.zf.setValue(4)
    dlg.p0.setValue(12)
    assert "transform failed" not in dlg.status.text()  # p0 != 0 used to break it
    dlg._accept()
    ppm, amp, meta = got[-1]
    assert meta["processing"] == [
        {"op": "fcor", "factor": 0.5}, {"op": "gm", "lb_hz": -30.0, "gb": 0.1},
        {"op": "zf", "factor": 4}, {"op": "ft", "offset_ppm": 0.0},
        {"op": "phase", "p0": 12.0, "p1": 0.0}]
    assert meta["processing_from_raw"] is True
    assert not np.iscomplexobj(amp)
    assert "processing" not in dlg.data.meta            # the reader's meta untouched

    def replay(m):
        r = Recipe(nucleus=m.get("nucleus", ""), larmor_frequency_MHz=m["larmor_MHz"],
                   processing=m["processing"], processing_from_raw=True)
        return loader.apply_processing(r, None, None, source_path=m["expno"])

    rp, ra, _ = replay(meta)
    order = np.argsort(ppm)
    assert np.allclose(rp, np.asarray(ppm)[order])
    assert np.allclose(ra, np.asarray(amp)[order], atol=1e-9 * np.abs(amp).max())

    dlg._autophase()
    dlg._accept()
    ppm2, amp2, meta2 = got[-1]
    names2 = _names(meta2["processing"])
    assert names2[-1] == "autophase" and "phase" not in names2
    _, ra2, _ = replay(meta2)
    assert np.allclose(ra2, np.asarray(amp2)[np.argsort(ppm2)],
                       atol=1e-6 * np.abs(amp2).max())
    dlg.close()
