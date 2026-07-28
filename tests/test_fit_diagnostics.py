"""Residual diagnostics (idea #4) and per-nucleus seeding (idea #17)."""
import os
import types

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")
import pytest

pytest.importorskip("PySide6")


def _win():
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow
    return MainWindow()


def test_residual_noise_ratio():
    from larmor.desktop.app import MainWindow
    x = np.linspace(-50, 150, 400)
    peak = np.exp(-0.5 * ((x - 60) / 8) ** 2)
    noise = np.random.RandomState(0).normal(0, 0.02, 400)
    good = types.SimpleNamespace(y_exp=peak + noise, y_fit=peak)
    bad = types.SimpleNamespace(y_exp=peak + noise, y_fit=0.7 * peak)
    assert MainWindow._residual_noise_ratio(good) < 1.5      # within noise
    assert MainWindow._residual_noise_ratio(bad) > 3         # structure left


def test_per_nucleus_seed():
    w = _win()
    w.recipe = {"nucleus": "27Al", "sites": []}
    p = {"sigma_Cq_MHz": {"value": 0}, "shift_fwhm_ppm": {"value": 0},
         "eta": {"value": 0}}
    w._seed_nucleus_defaults("czjzek", p)
    assert p["sigma_Cq_MHz"]["value"] == 1.5
    assert p["shift_fwhm_ppm"]["value"] == 12.0
    # a different nucleus seeds different values
    w.recipe["nucleus"] = "11B"
    p2 = {"sigma_Cq_MHz": {"value": 0}}
    w._seed_nucleus_defaults("czjzek", p2)
    assert p2["sigma_Cq_MHz"]["value"] == 1.0
    # non-quadrupolar models are left alone
    p3 = {"shift_fwhm_ppm": {"value": 5.0}}
    w._seed_nucleus_defaults("gauss_lor", p3)
    assert p3["shift_fwhm_ppm"]["value"] == 5.0
    # unknown nucleus -> no change
    w.recipe["nucleus"] = "207Pb"
    p4 = {"sigma_Cq_MHz": {"value": 0.3}}
    w._seed_nucleus_defaults("czjzek", p4)
    assert p4["sigma_Cq_MHz"]["value"] == 0.3
    w.close()
