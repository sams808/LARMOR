"""Variable-temperature Arrhenius / VFT fits."""
import numpy as np
import pytest

from larmor import vt


def test_arrhenius_recovers_activation_energy():
    R = 8.314462618
    Ea, A = 40000.0, 1e12
    T = np.linspace(250, 350, 12)
    k = A * np.exp(-Ea / (R * T))
    fit = vt.fit_arrhenius(T, k)
    assert fit["Ea_kJmol"] == pytest.approx(40.0, abs=0.1)
    assert fit["A"] == pytest.approx(1e12, rel=0.05)


def test_vft_recovers_parameters():
    T = np.linspace(250, 350, 15)
    k = 1e13 * np.exp(-800.0 / (T - 180.0))
    fit = vt.fit_vft(T, k)
    assert fit["B_K"] == pytest.approx(800.0, rel=0.05)
    assert fit["T0_K"] == pytest.approx(180.0, abs=3.0)
