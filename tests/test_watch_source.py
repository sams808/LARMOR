"""File > Watch the source file: the spectrum reloads when its file is
rewritten (spectrometer-side use), keeping the fit model."""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _write_csv(path, amp_scale):
    ppm = np.linspace(-60.0, 100.0, 512)
    amp = amp_scale * np.exp(-4 * np.log(2) * ((ppm - 10.0) / 4.0) ** 2)
    with open(path, "w", encoding="utf-8") as f:
        f.write("# LARMOR spectrum\n# nucleus=11B\n"
                "# larmor_MHz=160.46\n# spin_rate_Hz=20000.0\n")
        f.write("ppm,intensity\n")
        for x, y in zip(ppm, amp):
            f.write(f"{x},{y}\n")


def test_watch_reloads_on_change_and_keeps_the_fit(qapp, tmp_path):
    from PySide6.QtTest import QTest
    from larmor.desktop.app import MainWindow

    data = tmp_path / "growing.csv"
    _write_csv(data, 1.0)
    win = MainWindow()
    try:
        win.load_source(str(data), keep_fit=False)
        qapp.processEvents()
        assert np.max(win.exp_amp) == pytest.approx(1.0, rel=1e-2)   # grid misses the apex slightly
        from larmor import models
        m = models.get("gauss_lor")                # a model to keep
        win.recipe["sites"].append({"model": "gauss_lor", "label": "L", "params": {
            p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                     "min": p.min, "max": p.max, "expr": None} for p in m.params}})
        win.on_structure_changed()
        assert len(win.recipe["sites"]) == 1

        win.actWatch.setChecked(True)
        win._toggle_watch(True)
        assert win._watch_path == str(data)
        assert str(data) in win._watcher.files()

        _write_csv(data, 3.0)                      # the spectrometer rewrote it
        win._watched_changed()                     # (the OS signal, delivered)
        QTest.qWait(900)                           # debounce + reload
        assert np.max(win.exp_amp) == pytest.approx(3.0, rel=1e-2)
        assert len(win.recipe["sites"]) == 1       # keep_fit honoured
        assert win._watch_path == str(data)        # re-armed after the reload

        win.actWatch.setChecked(False)
        win._toggle_watch(False)
        assert win._watch_path is None and not win._watcher.files()
    finally:
        win.close()


def test_watch_without_a_source_is_refused(qapp):
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    try:
        win.actWatch.setChecked(True)
        win._toggle_watch(True)
        assert not win.actWatch.isChecked() and win._watch_path is None
    finally:
        win.close()
