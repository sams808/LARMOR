"""Offscreen smoke tests for the desktop app (no window is shown)."""
import json
import os

import numpy as np
import pytest

from conftest import CAALGLASS, require

pyside = pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture()
def win(qapp):
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


@pytest.mark.slow
def test_load_fit_quantify_undo(qapp, win):
    require(CAALGLASS)
    win.load_source(str(CAALGLASS))
    qapp.processEvents()
    assert win.recipe is not None and len(win.recipe["sites"]) == 5
    assert win.exp_ppm.size == 8192

    # from_dict must never mutate the window's recipe (regression:
    # shallow-copy bug found by this very scenario)
    from larmor.recipe import Recipe

    before = json.dumps(win.recipe)
    Recipe.from_dict(win.recipe)
    assert json.dumps(win.recipe) == before

    from larmor import fit as fitmod

    result = fitmod.fit(Recipe.from_dict(win.recipe),
                        win.exp_ppm, win.exp_amp, window_ppm=(150.0, -80.0))
    win._fit_done(result)
    qapp.processEvents()
    assert result.rmsd < 0.01
    assert win.qtable.rowCount() == 5

    # undo restores the site count
    n = len(win.recipe["sites"])
    win.snapshot()
    win.recipe["sites"].pop()
    win._after_structural_change()
    win.undo()
    assert len(win.recipe["sites"]) == n


def test_add_site_and_markers(qapp, win):
    require(CAALGLASS)
    win.load_source(str(CAALGLASS))
    qapp.processEvents()
    n = len(win.recipe["sites"])
    win._model_actions["gauss_lor"].setChecked(True)
    win.add_site_at(42.0, 123.0)
    assert len(win.recipe["sites"]) == n + 1
    new = win.recipe["sites"][-1]
    assert new["model"] == "gauss_lor"
    assert new["params"]["isotropic_chemical_shift_ppm"]["value"] == 42.0
    assert new["params"]["amplitude"]["value"] == 123.0
    # a draggable marker exists for it
    assert any(getattr(m, "site_index", None) == n
               for m in win.view._markers)

    # marker drag updates the recipe
    win.marker_moved(n, 55.5)
    assert win.recipe["sites"][n]["params"][
        "isotropic_chemical_shift_ppm"]["value"] == 55.5
