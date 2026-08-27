"""Project files: save/reopen all open 1D spectra + fits."""
import os
import tempfile
from pathlib import Path

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


def _write_csv_spectrum(path, nucleus, x, y):
    path.write_text("# nucleus = " + nucleus + "\n# larmor_MHz = 130.3\n"
                    + "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, y)))


def test_add_overlay_dialog_reads_a_recipe_json_source(tmp_path):
    """Regression: add_overlay_dialog used to unpack load_any()'s return
    tuple in the wrong order (recipe, ppm, amp) instead of the real
    (ppm, amp, recipe, ...), so recipe.get(...) always raised and EVERY
    overlay source _load_any actually supports (recipe.json, fxmla,
    csv/txt/dat) silently fell through to the Bruker-only fallback and
    failed there too -- caught while wiring overlay round-trip into
    save_project/open_project."""
    from PySide6.QtWidgets import QApplication, QFileDialog
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    x = np.linspace(-40, 120, 300)
    y = np.exp(-0.5 * ((x - 20) / 8) ** 2)
    csv_path = tmp_path / "overlay_source.csv"
    _write_csv_spectrum(csv_path, "27Al", x, y)

    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (str(csv_path), ""))
    n_before = len(w._overlays)
    w.add_overlay_dialog()
    assert len(w._overlays) == n_before + 1
    ov = w._overlays[-1]
    assert ov["source"] == str(csv_path)
    assert ov["label"] == "overlay_source"        # from the csv's own stem
    # abs=5e-4: the csv round-trips through "%.4f" text, not full precision
    assert np.allclose(ov["ppm"], x, atol=5e-4)
    assert np.allclose(ov["amp"], y, atol=5e-4)
    w.close()


def test_project_roundtrip_includes_overlays(tmp_path):
    """save_project captured overlays in the live snapshot but never wrote
    them into the saved file, so every overlay silently vanished on
    reopen -- fixed alongside the add_overlay_dialog bug above."""
    from PySide6.QtWidgets import QApplication, QFileDialog
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    x = np.linspace(-40, 120, 300)
    y = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    w.exp_ppm, w.exp_amp = x, y
    w.recipe = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3,
               "sample": "glassA", "sites": []}
    w.source_path = str(tmp_path / "main.csv")
    _write_csv_spectrum(Path(w.source_path), "27Al", x, y)
    w.hidden = set()
    w.view.set_experiment(x, y)
    w.lines_table.rebuild(w.recipe, w.hidden)

    ov_path = tmp_path / "compare.csv"
    ov_y = np.exp(-0.5 * ((x - 30) / 8) ** 2)
    _write_csv_spectrum(ov_path, "27Al", x, ov_y)
    w._add_overlay("compare", x, ov_y, str(ov_path))
    w._overlays[-1]["visible"] = False

    w._ws_mode = "new"
    w._register_ws("1d")

    proj = tempfile.mktemp(suffix=".larproj.json")
    QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.save_project()

    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (proj, ""))
    w.open_project()

    assert len(w.workspaces) == 1
    overlays = w.workspaces[0]["snap"]["overlays"]
    assert len(overlays) == 1
    assert overlays[0]["label"] == "compare"
    assert overlays[0]["source"] == str(ov_path)
    assert overlays[0]["visible"] is False
    # peak position/height, not a full-array compare: a Gaussian's tail spans
    # many orders of magnitude and the csv's "%.4f" text rounds far-tail
    # values to exactly 0 -- checking the peak is what actually matters here
    # (that the right FILE came back), not bit-for-bit tail precision
    restored_amp = np.asarray(overlays[0]["amp"])
    restored_ppm = np.asarray(overlays[0]["ppm"])
    assert restored_amp.max() == pytest.approx(ov_y.max(), abs=1e-3)
    assert restored_ppm[restored_amp.argmax()] == pytest.approx(30.0, abs=0.5)
    w.close()
