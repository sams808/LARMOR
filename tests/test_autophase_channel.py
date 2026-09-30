"""Phasing a processed spectrum acts on a COMPLEX spectrum.

A Bruker 1r was loaded real-only, so 'Autophase' and a typed p0 rotated a
spectrum whose imaginary part was identically zero: y * exp(i*phi) taken
real is y * cos(phi) -- Autophase 'did nothing' and p0 only scaled the
line. Now the 1i TopSpin writes next to the 1r is read and used, and a
source without one (CSV, dmfit, magnitude data) gets its imaginary channel
reconstructed by a Hilbert transform before the first phase step.
"""
import os
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO = Path(__file__).resolve().parents[1]
EXPNO_27AL = REPO / "examples" / "pCABS2-4" / "3616"
ONE_R = EXPNO_27AL / "pdata" / "1" / "1r"


# ------------------------------------------------------------- Qt-free parts
def test_read_imag_gives_the_1i_on_the_readers_ppm_order():
    import nmrglue as ng

    from larmor.io import bruker

    imag = bruker.read_imag(str(ONE_R))
    assert imag is not None and imag.shape == (2048,)
    data = bruker.read(str(ONE_R))
    assert data.data.shape == imag.shape
    # the raw 1i, scaled and sorted like the reader sorts the 1r
    dic, raw = ng.bruker.read_pdata(str(ONE_R.parent), bin_files=["1i"])
    ppm = bruker._ppm_axis(dic["procs"], raw.size)
    assert np.allclose(imag, np.asarray(raw, float)[np.argsort(ppm)])
    assert np.abs(imag).max() > 0
    # an EXPNO folder resolves to its 1r as well
    assert bruker.read_imag(str(EXPNO_27AL)) is not None
    # not a processed 1D Bruker spectrum: None, never an exception
    assert bruker.read_imag(str(REPO / "examples" / "README.md")) is None
    assert bruker.read_imag(str(EXPNO_27AL / "fid")) is None


def test_complex_for_phase_inserts_hilbert_only_when_needed():
    from larmor import processing as proc
    from larmor.desktop.mw_processing import _complex_for_phase

    x = np.linspace(-10, 10, 64)
    real_only = proc.from_processed(x, np.exp(-x ** 2), 100.0)
    ops, added = _complex_for_phase(real_only, [{"op": "autophase"}])
    assert added and ops[0] == {"op": "hilbert"} and ops[1] == {"op": "autophase"}
    ops, added = _complex_for_phase(real_only, [{"op": "baseline", "order": 3},
                                                {"op": "phase", "p0": 30.0}])
    assert added and [o["op"] for o in ops] == ["hilbert", "baseline", "phase"]
    # nothing to phase: untouched
    ops, added = _complex_for_phase(real_only, [{"op": "baseline", "order": 3}])
    assert not added and ops == [{"op": "baseline", "order": 3}]
    # the chain already makes a complex signal before the phase step
    for maker in ("hilbert", "ift"):
        ops, added = _complex_for_phase(real_only, [{"op": maker}, {"op": "phase"}])
        assert not added
    # a true imaginary channel: untouched
    cplx = proc.from_processed(x, np.exp(-x ** 2) * (1 + 0.3j), 100.0)
    ops, added = _complex_for_phase(cplx, [{"op": "autophase"}])
    assert not added and ops == [{"op": "autophase"}]
    # time-domain input (a raw fid) is complex by nature
    fid = proc.Spectrum1D(x_ppm=None, y=np.ones(64, complex), sfo1_MHz=100.0,
                          sw_Hz=1e4, domain="time")
    ops, added = _complex_for_phase(fid, [{"op": "ft"}, {"op": "autophase"}])
    assert not added


# ------------------------------------------------------------- the window
pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    monkeypatch.setenv("LARMOR_NO_KERNEL_WARM", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def _write_csv_spectrum(path):
    ppm = np.linspace(-60.0, 100.0, 512)
    amp = np.exp(-4 * np.log(2) * ((ppm - 10.0) / 4.0) ** 2)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# LARMOR spectrum\n# nucleus=11B\n# larmor_MHz=160.46\n"
                "# spin_rate_Hz=20000.0\nppm,intensity\n")
        for x, y in zip(ppm, amp):
            f.write(f"{x},{y}\n")
    return ppm, amp


def _neg_fraction(y, ppm=None, window=(100.0, -40.0)):
    """|negative| / |all| inside the line's window (the shipped 27Al 1r has
    a rolling baseline far from the line that would swamp a whole-trace
    figure)."""
    if ppm is not None:
        sel = (ppm < window[0]) & (ppm > window[1])
        y = y[sel]
    return float(np.abs(y[y < 0]).sum() / np.abs(y).sum())


def _corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def test_bruker_1r_phases_with_its_own_1i(win, qapp):
    win.load_source(str(ONE_R), keep_fit=False)
    qapp.processEvents()
    assert win._exp_imag is not None and win._exp_imag.shape == win.exp_amp.shape
    real0, imag0 = win.exp_amp.copy(), win._exp_imag.copy()

    # a 90 degree zero-order rotation: the real channel becomes -imag
    win.apply_processing([{"op": "phase", "p0": 90.0, "p1": 0.0}], False)
    qapp.processEvents()
    assert np.allclose(win.exp_amp, -imag0, atol=1e-6 * np.abs(real0).max())
    assert not win.proc_panel.chkHilbert.isChecked()    # the true 1i, no stand-in
    assert win.recipe["processing"] == [{"op": "phase", "p0": 90.0, "p1": 0.0,
                                         "pivot_frac": win.view.phase_pivot_frac()}]
    assert "Hilbert" not in win.statusBar().currentMessage()

    # dephase by 40 degrees, then Autophase brings the line back upright.
    # Measured: dephased vs TopSpin's 1r correlates 0.77 and has 21 % of its
    # area negative around the line (3.5 % before); autophased 0.95 and
    # 0.0 % -- the 'scan' criterion zeroes the negative lobes, TopSpin's
    # operator left a little; both are 'phased', the sign of the change is
    # what this test pins
    ppm = win.exp_ppm
    win.apply_processing([{"op": "phase", "p0": 40.0, "p1": 0.0}], False)
    qapp.processEvents()
    dephased = win.exp_amp.copy()
    assert _neg_fraction(dephased, ppm) > 3 * _neg_fraction(real0, ppm)
    assert _corr(dephased, real0) < 0.85
    win.apply_processing([{"op": "phase", "p0": 40.0, "p1": 0.0},
                          {"op": "autophase"}], False)
    qapp.processEvents()
    assert _corr(win.exp_amp, real0) > 0.93
    assert _corr(win.exp_amp, real0) > _corr(dephased, real0) + 0.1
    assert _neg_fraction(win.exp_amp, ppm) < _neg_fraction(dephased, ppm) / 3
    assert np.abs(win.exp_amp).max() == pytest.approx(np.abs(real0).max(), rel=0.1)


def test_csv_without_imag_gets_a_hilbert_reconstruction_before_phasing(
        win, qapp, tmp_path):
    data = tmp_path / "line.csv"
    ppm, amp = _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    assert win._exp_imag is None
    assert not win.proc_panel.chkHilbert.isChecked()

    win.apply_processing([{"op": "phase", "p0": 90.0, "p1": 0.0}], False)
    qapp.processEvents()
    names = [o["op"] for o in win.recipe["processing"]]
    assert names == ["hilbert", "phase"]
    assert win.proc_panel.chkHilbert.isChecked()
    assert "Hilbert" in win.statusBar().currentMessage()
    # a real rotation, not a cos(90) = 0 scaling: the dispersive shape has
    # an antisymmetric lobe pair around the line
    y = win.exp_amp
    assert np.abs(y).max() > 0.3 * amp.max()
    assert y[np.argmax(ppm > 12.0)] * y[np.argmax(ppm > 8.0)] < 0

    # Autophase from the dephased state lands back on the absorption line
    # (0.98 measured: the Hilbert reconstruction of a 512-point CSV is not
    # exact at the edges; the dephased state correlates 0.5)
    win.apply_processing([{"op": "phase", "p0": 60.0}], False)
    qapp.processEvents()
    dephased = win.exp_amp.copy()
    win.apply_processing([{"op": "phase", "p0": 60.0}, {"op": "autophase"}], False)
    qapp.processEvents()
    assert [o["op"] for o in win.recipe["processing"]] == ["hilbert", "phase",
                                                            "autophase"]
    assert _corr(win.exp_amp, amp) > 0.95
    assert _corr(win.exp_amp, amp) > _corr(dephased, amp) + 0.2


def test_imag_survives_a_workspace_round_trip(win, qapp, tmp_path):
    win.load_source(str(ONE_R), keep_fit=False)
    qapp.processEvents()
    imag = win._exp_imag.copy()
    data = tmp_path / "line.csv"
    _write_csv_spectrum(data)
    win._ws_mode = "new"          # an unfitted 1D workspace is otherwise reused
    win.load_source(str(data), keep_fit=False)         # a second workspace
    qapp.processEvents()
    assert len(win.workspaces) == 2 and win.active_ws == 1
    assert win._exp_imag is None
    win.switch_workspace(0)
    qapp.processEvents()
    assert win._exp_imag is not None and np.array_equal(win._exp_imag, imag)
    win.switch_workspace(1)
    qapp.processEvents()
    assert win._exp_imag is None
