"""Autophase resolves to explicit p0 / p1 that land in the Phase controls.

An {"op": "autophase"} step used to stay opaque in the recorded chain and
invisible to the processing panel: the panel emits the absolute chain from
its widgets on every live tick, so nudging p0 by one degree after Autophase
re-applied p0 = 1 to the RAW spectrum and the line snapped back -- 'autophase
does not stick'. The menu entries (Autophase, Polynomial, Subtract averages,
Iterative) replaced the recorded chain outright, dropping an earlier phase.
Now the angles Autophase finds are folded into a plain phase step (the way
TopSpin's apk fills PHC0 / PHC1), the panel is synced to the recorded chain
after every apply, and the menu entries append to the chain.
"""
import os
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------- Qt-free parts
def _dephased_line(p0_deg=37.0, p1_deg=-20.0, n=1024):
    """A complex Lorentzian on a ppm axis, rotated by (p0, p1) about the
    centre -- what a badly phased spectrum looks like."""
    from larmor import processing as proc

    x = np.linspace(-50.0, 50.0, n)
    y = 1.0 / (1.0 + ((x - 5.0) / 3.0) ** 2) + 0.0j
    y = y + 0.6 / (1.0 + ((x + 20.0) / 4.0) ** 2)
    s = proc.from_processed(x, y, 100.0)
    s.y = s.y * np.exp(1j * np.deg2rad(p0_deg + p1_deg * (np.arange(n) / (n - 1) - 0.5)))
    return s


def _copy(s):
    import dataclasses
    return dataclasses.replace(s, y=np.array(s.y, copy=True),
                               x_ppm=np.array(s.x_ppm, copy=True))


def test_autophase_angles_reproduce_op_autophase_and_find_the_rotation():
    from larmor import processing as proc

    s = _dephased_line()
    p0, p1 = proc.autophase_angles(_copy(s))
    # the explicit step is the same operation as the opaque one
    a = proc.op_phase(_copy(s), p0, p1, pivot_frac=0.5)
    b = proc.op_autophase(_copy(s))
    assert np.allclose(a.y, b.y, atol=1e-9 * np.abs(b.y).max())
    # and it undoes the rotation the line was given (p0 wrapped to +-180)
    assert -180.0 < p0 <= 180.0
    assert p0 == pytest.approx(-37.0, abs=3.0)
    assert p1 == pytest.approx(20.0, abs=8.0)
    # the result is upright: hardly any negative area
    yr = b.y.real
    assert np.abs(yr[yr < 0]).sum() < 0.03 * np.abs(yr).sum()


def test_autophase_angles_acme_matches_nmrglue():
    ng = pytest.importorskip("nmrglue")
    from larmor import processing as proc

    s = _dephased_line(p0_deg=25.0, p1_deg=0.0)
    ref = ng.process.proc_autophase.autops(np.array(s.y, copy=True), "acme",
                                           disp=False)
    p0, p1 = proc.autophase_angles(_copy(s), method="acme")
    out = proc.op_phase(_copy(s), p0, p1, pivot_frac=0.5)
    assert np.allclose(out.y, ref, atol=1e-6 * np.abs(ref).max())
    assert np.allclose(proc.op_autophase(_copy(s), method="acme").y, ref,
                       atol=1e-6 * np.abs(ref).max())


def test_resolve_autophase_replaces_the_step_in_place_and_folds_into_a_phase():
    from larmor import processing as proc

    s = _dephased_line()
    f = 0.3                                              # a peak pivot

    # alone: one explicit phase step about the requested pivot
    ops, folds = proc.resolve_autophase(_copy(s), [{"op": "autophase"}], pivot_frac=f)
    assert len(folds) == 1 and len(ops) == 1
    assert ops[0]["op"] == "phase" and ops[0]["pivot_frac"] == f
    assert np.allclose(proc.apply(_copy(s), ops).y,
                       proc.apply(_copy(s), [{"op": "autophase"}]).y,
                       atol=1e-9 * np.abs(s.y).max())
    assert -180.0 < ops[0]["p0"] <= 180.0

    # after a phase step about the same pivot: folded into that step
    chain = [{"op": "phase", "p0": 10.0, "p1": 4.0, "pivot_frac": f},
             {"op": "autophase"}]
    ops2, folds2 = proc.resolve_autophase(_copy(s), chain, pivot_frac=f)
    assert len(ops2) == 1 and ops2[0]["op"] == "phase"
    assert ops2[0]["p0"] == pytest.approx(10.0 + folds2[0][0])
    assert ops2[0]["p1"] == pytest.approx(4.0 + folds2[0][1])
    assert np.allclose(proc.apply(_copy(s), ops2).y, proc.apply(_copy(s), chain).y,
                       atol=1e-9 * np.abs(s.y).max())

    # after a baseline step: the phase goes where the controls keep it --
    # AHEAD of the steps the panel merely carries -- so the panel's next
    # live tick re-emits the same chain; the angles are the ones measured
    # after the baseline, so the result agrees to the baseline's own
    # sensitivity to the phase, not bit for bit
    chain3 = [{"op": "baseline", "order": 1}, {"op": "autophase"},
              {"op": "baseline", "order": 1}]
    ops3, folds3 = proc.resolve_autophase(_copy(s), chain3, pivot_frac=f)
    assert [o["op"] for o in ops3] == ["phase", "baseline", "baseline"]
    assert len(folds3) == 1
    direct = proc.apply(_copy(s), chain3).y.real
    folded = proc.apply(_copy(s), ops3).y.real
    assert np.corrcoef(direct, folded)[0, 1] > 0.999

    # ... and it folds into an earlier phase step across a carried baseline
    chain3b = [{"op": "hilbert"}, {"op": "phase", "p0": 30.0, "p1": 0.0, "pivot_frac": f},
               {"op": "baseline", "order": 1}, {"op": "autophase"}]
    ops3b, folds3b = proc.resolve_autophase(_copy(s), chain3b, pivot_frac=f)
    assert [o["op"] for o in ops3b] == ["hilbert", "phase", "baseline"]
    assert ops3b[1]["p0"] == pytest.approx(30.0 + folds3b[0][0]) or \
        abs(abs(ops3b[1]["p0"] - (30.0 + folds3b[0][0])) - 360.0) < 1e-9
    # a transform block ahead: the phase is inserted right after it
    chain3c = [{"op": "hilbert"}, {"op": "baseline", "order": 1}, {"op": "autophase"}]
    ops3c, _ = proc.resolve_autophase(_copy(s), chain3c, pivot_frac=f)
    assert [o["op"] for o in ops3c] == ["hilbert", "phase", "baseline"]

    # a phase step about ANOTHER pivot is not folded (different pivot, kept)
    chain4 = [{"op": "phase", "p0": 10.0, "p1": 4.0, "pivot_frac": 0.5},
              {"op": "autophase"}]
    ops4, _ = proc.resolve_autophase(_copy(s), chain4, pivot_frac=f)
    assert [o["op"] for o in ops4] == ["phase", "phase"]
    assert ops4[1]["pivot_frac"] == f

    # nothing to resolve: the chain comes back as it went in
    plain = [{"op": "baseline", "order": 2}]
    ops5, folds5 = proc.resolve_autophase(_copy(s), plain, pivot_frac=f)
    assert ops5 == plain and folds5 == []


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


def _flush(qapp, win):
    win.proc_panel._live_timer.stop()
    win.proc_panel._live_timer.timeout.emit()
    qapp.processEvents()


def _names(ops):
    return [o["op"] for o in ops]


def test_panel_autophase_fills_the_phase_controls_and_sticks(win, qapp, tmp_path):
    data = tmp_path / "line.csv"
    ppm, amp = _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    pp = win.proc_panel

    # dephase with the panel first (a typed p0), then press Autophase
    pp.p0v.setValue(60.0)
    _flush(qapp, win)
    assert _names(win.recipe["processing"]) == ["hilbert", "phase"]
    dephased = win.exp_amp.copy()

    pp.btnAuto.click()
    qapp.processEvents()
    chain = win.recipe["processing"]
    assert _names(chain) == ["hilbert", "phase"]          # no opaque step left
    step = chain[1]
    assert pp.p0v.value() == pytest.approx(step["p0"], abs=0.01)
    assert pp.p1v.value() == pytest.approx(step["p1"], abs=0.01)
    assert "autophase" in win.statusBar().currentMessage()
    assert "p0" in win.statusBar().currentMessage()
    after_auto = win.exp_amp.copy()
    assert np.corrcoef(after_auto, amp)[0, 1] > np.corrcoef(dephased, amp)[0, 1] + 0.2

    # a nudge of the p0 control now moves the spectrum by ONE degree from the
    # autophased state, not back to the raw spectrum plus one degree
    pp.p0v.setValue(pp.p0v.value() + 1.0)
    _flush(qapp, win)
    assert _names(win.recipe["processing"]) == ["hilbert", "phase"]
    assert win.recipe["processing"][1]["p0"] == pytest.approx(step["p0"] + 1.0, abs=0.01)
    assert np.corrcoef(win.exp_amp, after_auto)[0, 1] > 0.999


def test_menu_entries_append_to_the_recorded_chain(win, qapp, tmp_path):
    data = tmp_path / "line.csv"
    _write_csv_spectrum(data)
    win.load_source(str(data), keep_fit=False)
    qapp.processEvents()
    pp = win.proc_panel
    pp.p0v.setValue(30.0)
    _flush(qapp, win)
    assert _names(win.recipe["processing"]) == ["hilbert", "phase"]

    # Process > Baseline > Polynomial keeps the phase in front of it
    assert win.append_processing_step({"op": "baseline", "order": 3})
    assert _names(win.recipe["processing"]) == ["hilbert", "phase", "baseline"]
    assert win.recipe["processing"][1]["p0"] == pytest.approx(30.0)
    # ... and the panel carries the baseline through its next live tick
    pp.p0v.setValue(31.0)
    _flush(qapp, win)
    assert _names(win.recipe["processing"]) == ["hilbert", "phase", "baseline"]
    assert win.recipe["processing"][1]["p0"] == pytest.approx(31.0)

    # Process > Phase > Autophase folds into the panel's phase step
    assert win.append_processing_step({"op": "autophase"})
    chain = win.recipe["processing"]
    assert _names(chain) == ["hilbert", "phase", "baseline"]
    assert pp.p0v.value() == pytest.approx(chain[1]["p0"], abs=0.01)
    assert pp.p1v.value() == pytest.approx(chain[1]["p1"], abs=0.01)

    # the panel's own Correct button is carried too (it used to be dropped
    # by the very next live tick)
    pp.btnBaseline.click()
    qapp.processEvents()
    assert _names(win.recipe["processing"]).count("baseline") == 2
    pp.p0v.setValue(pp.p0v.value() + 1.0)
    _flush(qapp, win)
    assert _names(win.recipe["processing"]).count("baseline") == 2


def test_autophase_on_a_bruker_1r_keeps_topspins_channel(win, qapp):
    one_r = REPO / "examples" / "pCABS2-4" / "3616" / "pdata" / "1" / "1r"
    win.load_source(str(one_r), keep_fit=False)
    qapp.processEvents()
    pp = win.proc_panel
    real0 = win.exp_amp.copy()
    pp.btnAuto.click()
    qapp.processEvents()
    chain = win.recipe["processing"]
    assert _names(chain) == ["phase"]                     # no Hilbert: the 1i is real
    assert not pp.chkHilbert.isChecked()
    assert pp.p0v.value() == pytest.approx(chain[0]["p0"], abs=0.01)
    # TopSpin's operator phased this one: the autophase only trims it
    assert np.corrcoef(win.exp_amp, real0)[0, 1] > 0.9
