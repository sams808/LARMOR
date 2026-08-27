"""Parameter correlation heat-map."""
import os
import types

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest

pytest.importorskip("PySide6")


def test_corr_color_endpoints():
    from larmor.desktop.correlation_dialog import _corr_color
    r = _corr_color(1.0); b = _corr_color(-1.0); w = _corr_color(0.0)
    assert (r.red(), r.blue()) == (255, 0)          # +1 red
    assert (b.blue(), b.red()) == (255, 0)          # −1 blue
    assert (w.red(), w.green(), w.blue()) == (255, 255, 255)   # 0 white


def test_dialog_builds_and_handles_empty():
    from PySide6.QtWidgets import QApplication, QTableWidget
    QApplication.instance() or QApplication([])
    from larmor.desktop.correlation_dialog import CorrelationDialog

    # off-diagonal chosen so corr[0,2] = -0.95/√(1·2)·√2 = -0.95 exactly
    c02 = -0.95 * np.sqrt(1.0 * 2.0)
    cov = np.array([[1.0, 0.0, c02], [0.0, 4.0, 0.1], [c02, 0.1, 2.0]])
    lm = types.SimpleNamespace(
        var_names=["s0_sigma_Cq_MHz", "s0_amplitude", "s1_amplitude"], covar=cov)
    d = CorrelationDialog(None, lm)
    t = d.findChild(QTableWidget)
    assert t is not None and t.rowCount() == 3
    assert t.item(0, 0).text() == "+1.00"           # diagonal
    assert t.item(0, 2).text() == "-0.95"           # strong anti-correlation
    d.close()
    # missing covariance is handled gracefully (no table, no crash)
    d2 = CorrelationDialog(None, types.SimpleNamespace(var_names=[], covar=None))
    d2.close()
