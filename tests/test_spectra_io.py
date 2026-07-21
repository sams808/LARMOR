"""Background subtraction and the reusable CSV spectrum format."""
import numpy as np
import pytest

from larmor.io import spectra
from larmor.loader import load_any


def test_csv_roundtrip_carries_metadata(tmp_path):
    x = np.linspace(-100, 100, 200)
    y = np.exp(-(x / 20) ** 2)
    p = tmp_path / "s.csv"
    spectra.write_csv(p, x, y, {"nucleus": "27Al", "larmor_MHz": 195.483,
                                "spin_rate_Hz": 20000.0, "sample": "glass"})
    ppm, amp, meta = spectra.read_csv(p)
    assert meta["nucleus"] == "27Al"
    assert meta["larmor_MHz"] == pytest.approx(195.483)
    assert meta["spin_rate_Hz"] == pytest.approx(20000.0)
    assert np.allclose(np.sort(amp), np.sort(y), atol=1e-6)


def test_load_any_reads_csv(tmp_path):
    x = np.linspace(0, 50, 100); y = np.sin(x)
    p = tmp_path / "spec.csv"
    spectra.write_csv(p, x, y, {"nucleus": "11B", "larmor_MHz": 160.0})
    ppm, amp, recipe, meta, warn = load_any(str(p))
    assert recipe["nucleus"] == "11B"
    assert recipe["larmor_frequency_MHz"] == pytest.approx(160.0)
    assert ppm.size == 100


def test_plain_two_column_txt_loads(tmp_path):
    p = tmp_path / "plain.txt"
    p.write_text("10 1.0\n5 2.0\n0 3.0\n")          # no header, whitespace
    ppm, amp, meta = spectra.read_csv(p)
    assert ppm.size == 3 and amp.size == 3


def test_subtract_removes_background():
    x = np.linspace(-100, 100, 300)
    peak = np.exp(-((x - 10) / 5) ** 2)
    bg = 0.4 * np.exp(-((x + 30) / 8) ** 2)
    diff = spectra.subtract(x, peak + bg, x, bg, scale=1.0)
    assert abs(diff[np.argmin(np.abs(x + 30))]) < 0.02      # background gone
    assert diff[np.argmin(np.abs(x - 10))] == pytest.approx(1.0, abs=0.02)


def test_best_scale_recovers_amount():
    x = np.linspace(-50, 50, 400)
    bg = np.exp(-((x + 10) / 6) ** 2)
    sample = np.exp(-((x - 5) / 4) ** 2) + 2.5 * bg
    assert spectra.best_scale(x, sample, x, bg) == pytest.approx(2.5, rel=0.02)
