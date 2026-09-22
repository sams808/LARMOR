"""The Czjzek C_Q distribution implied by a fitted σ."""
import numpy as np
import pytest

from larmor import czjzek_dist as cd


def test_marginal_mode_and_rms_in_kernel_convention():
    """The stored σ is mrsimulator's (σ_Cz = 2σ): the |C_Q| marginal peaks at
    ≈ 3.73σ (1.87 σ_Cz, ≈ dmfit's displayed CQ = 4σ), NOT at 2σ -- the old
    read-out described a distribution half as wide as the one the fit used."""
    sigma = 1.8
    cq = cd.suggested_cq_axis(sigma, 2000)
    m = cd.marginal_cq(sigma, cq)
    assert m.sum() == pytest.approx(1.0)
    mode = cq[int(np.argmax(m))]
    assert 3.5 * sigma < mode < 4.0 * sigma        # ≈ 3.73σ


def test_rms_pq_is_two_sqrt_d_sigma():
    """P_Q is chi-distributed (d dof, scale σ_Cz = 2σ): √⟨P_Q²⟩ = 2√d·σ."""
    for sigma in (0.5, 1.8, 3.0):
        assert cd.rms_pq(sigma) == pytest.approx(2.0 * np.sqrt(5) * sigma)
        for d in (2, 3, 5):
            assert cd.rms_pq(sigma, d) == pytest.approx(2.0 * np.sqrt(d) * sigma)
    assert cd.mode_pq(1.0) == pytest.approx(4.0)               # 2σ√(d−1), d = 5
    assert cd.mode_pq(1.0, 1.0) == 0.0


def test_joint_pdf_normalised_and_nonneg():
    cq = cd.suggested_cq_axis(2.0, 200)
    eta = np.linspace(0, 1, 60)
    p = cd.czjzek_pdf(2.0, cq, eta)
    assert p.shape == (60, 200)
    assert (p >= 0).all()
    assert p.sum() == pytest.approx(1.0)
    # the general-d family keeps the same contract
    p3 = cd.czjzek_pdf(2.0, cq, eta, d=3)
    assert p3.shape == (60, 200) and (p3 >= 0).all()
    assert p3.sum() == pytest.approx(1.0)


def test_suggested_axis_covers_the_distribution():
    """10σ = 5σ_Cz holds all but 4.5e-5 of the d = 5 mass; 5σ dropped 21 %."""
    cq = cd.suggested_cq_axis(1.0, 4000)
    assert cq[-1] == pytest.approx(10.0)
    m = cd.marginal_cq(1.0, cq)
    assert m[cq > 5.0].sum() == pytest.approx(0.213, abs=0.01)


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
