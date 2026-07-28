"""The Czjzek C_Q distribution implied by a fitted σ (idea #2)."""
import numpy as np
import pytest

from larmor import czjzek_dist as cd


def test_marginal_peaks_near_2sigma_and_normalised():
    sigma = 1.8
    cq = cd.suggested_cq_axis(sigma, 800)
    m = cd.marginal_cq(sigma, cq)
    assert m.sum() == pytest.approx(1.0)
    mode = cq[int(np.argmax(m))]
    assert 1.7 * sigma < mode < 2.1 * sigma        # peaks a touch below 2σ


def test_rms_pq_is_sqrt5_sigma():
    for sigma in (0.5, 1.8, 3.0):
        assert cd.rms_pq(sigma) == pytest.approx(np.sqrt(5) * sigma)


def test_joint_pdf_normalised_and_nonneg():
    cq = cd.suggested_cq_axis(2.0, 200)
    eta = np.linspace(0, 1, 60)
    p = cd.czjzek_pdf(2.0, cq, eta)
    assert p.shape == (60, 200)
    assert (p >= 0).all()
    assert p.sum() == pytest.approx(1.0)


def test_dialog_builds():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.czjzek_dist_dialog import CzjzekDistDialog

    r = {"nucleus": "27Al", "sites": [
        {"model": "czjzek", "label": "AlIV",
         "params": {"sigma_Cq_MHz": {"value": 1.8}}}]}
    d = CzjzekDistDialog(None, r)
    d.close()
