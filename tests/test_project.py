"""Project files: save/reopen all open 1D spectra + fits (idea #8)."""
import os
import tempfile

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")
import pytest

pytest.importorskip("PySide6")


def test_project_roundtrip():
    from PySide6.QtWidgets import QApplication, QFileDialog
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    w = MainWindow()

    def mkws(sample, pos):
        x = np.linspace(-40, 120, 300)
        y = np.exp(-0.5 * ((x - pos) / 8) ** 2)
        w.exp_ppm, w.exp_amp = x, y
        w.recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
                    "sample": sample, "sites": [
                        {"model": "gauss_lor", "label": "A", "params": {
                            "isotropic_chemical_shift_ppm": {"value": pos},
                            "shift_fwhm_ppm": {"value": 8},
                            "gl": {"value": 0.5},
                            "amplitude": {"value": 1.0}}}]}
        w.source_path = "src_" + sample
        w.hidden = set()
        w.view.set_experiment(x, y)
        w.lines_table.rebuild(w.recipe, w.hidden)
        w._ws_mode = "new"
        w._register_ws("1d")

    mkws("glassA", 60)
    mkws("glassB", 30)
    assert len(w.workspaces) == 2

    proj = tempfile.mktemp(suffix=".larproj.json")
    QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.save_project()
    assert os.path.exists(proj)

    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.open_project()
    assert len(w.workspaces) == 2
    samples = [ws["snap"]["recipe"]["sample"] for ws in w.workspaces]
    assert samples == ["glassA", "glassB"]
    pos = [ws["snap"]["recipe"]["sites"][0]["params"]
           ["isotropic_chemical_shift_ppm"]["value"] for ws in w.workspaces]
    assert pos == [60, 30]
    w.close()
