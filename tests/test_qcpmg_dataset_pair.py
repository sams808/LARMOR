"""*Save as dataset…* writes the sum-echo envelope AND the spikelet spectrum.

One Save dialog, two LARMOR .csv files next to each other
(``<base>_sumecho.csv`` + ``<base>_spikelets.csv``), each on its own ppm axis
and -- by default -- normalised to unit maximum, with the factor divided out
kept in the header as ``intensity_scale`` so the raw intensity stays
recoverable. Both must open through ``load_any`` and the multi-field readers.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from larmor import qcpmg  # noqa: E402
from larmor.io import spectra  # noqa: E402
from larmor.loader import load_any  # noqa: E402


# ------------------------------------------------------------ Qt-free core
def test_dataset_pair_paths_from_any_spelling(tmp_path):
    want = {"sumecho": tmp_path / "qcpmg_7_sumecho.csv",
            "spikelets": tmp_path / "qcpmg_7_spikelets.csv"}
    for chosen in ("qcpmg_7_sumecho.csv", "qcpmg_7_spikelets.csv",
                   "qcpmg_7.csv", "qcpmg_7", "qcpmg_7_SumEcho.csv",
                   "qcpmg_7 spikelet.csv"):
        assert qcpmg.dataset_pair_paths(tmp_path / chosen) == want, chosen
    # a name that IS only the token keeps it as the base rather than vanishing
    got = qcpmg.dataset_pair_paths(tmp_path / "spikelets.csv")
    assert got["sumecho"].name == "spikelets_sumecho.csv"
    # a dotted sample name survives (only a trailing .csv is a suffix)
    got = qcpmg.dataset_pair_paths(tmp_path / "LAW3CL0CA_18.8T")
    assert got["spikelets"].name == "LAW3CL0CA_18.8T_spikelets.csv"


def test_normalise_to_max_returns_the_factor_divided_out():
    y = np.array([0.5, -4.0, 2.0])
    n, s = qcpmg.normalise_to_max(y)
    assert s == 4.0 and np.allclose(n * s, y)
    assert np.max(np.abs(n)) == 1.0
    # an all-zero trace is returned unchanged with a unit scale (no 0/0)
    z, s0 = qcpmg.normalise_to_max(np.zeros(5))
    assert s0 == 1.0 and np.all(z == 0.0)
    # a NaN gap does not poison the maximum
    n2, s2 = qcpmg.normalise_to_max(np.array([1.0, np.nan, 3.0]))
    assert s2 == 3.0 and n2[-1] == 1.0


def _twin_inputs():
    x1 = np.linspace(-300.0, 100.0, 801)
    env = 7.5 * np.exp(-((x1 + 100.0) / 40.0) ** 2)
    x2 = np.linspace(-300.0, 100.0, 3201)
    comb = 40.0 * np.exp(-((x2 + 100.0) / 40.0) ** 2) * (
        1.0 + np.cos(2 * np.pi * x2 / 5.0)) / 2.0
    base = {"nucleus": "35Cl", "larmor_MHz": 78.354, "spin_rate_Hz": 0.0,
            "qcpmg_rotor_Hz": 16000.0, "spikelet_spacing_Hz": 391.9,
            "echo_period_pts": 293, "split_offset_pts": 0, "n_echoes": 38,
            "carrier_ppm": -102.8, "referenced": True}
    m_sum = dict(base, sample="LAW · QCPMG sum echo (LB 75 Hz)",
                 spectrum_mode="absorption", lb_Hz=75.0)
    m_spk = dict(base, sample="LAW · QCPMG spikelets (LB 75 Hz, magnitude)",
                 spectrum_mode="magnitude", lb_Hz=75.0)
    return x1, env, x2, comb, m_sum, m_spk


def test_write_dataset_pair_normalises_and_cross_references(tmp_path):
    x1, env, x2, comb, m_sum, m_spk = _twin_inputs()
    out = qcpmg.write_dataset_pair(tmp_path / "LAW.csv", x1, env, x2, comb,
                                   m_sum, m_spk)
    assert set(out) == {"sumecho", "spikelets"}
    assert os.path.basename(out["sumecho"]) == "LAW_sumecho.csv"
    assert os.path.basename(out["spikelets"]) == "LAW_spikelets.csv"

    ppm_s, amp_s, meta_s = spectra.read_csv(out["sumecho"])
    ppm_k, amp_k, meta_k = spectra.read_csv(out["spikelets"])
    # both normalised to unit maximum, each on its OWN axis
    assert amp_s.max() == pytest.approx(1.0) and amp_k.max() == pytest.approx(1.0)
    assert ppm_s.size == x1.size and ppm_k.size == x2.size
    # the factor divided out is the raw maximum, so raw = intensity * scale
    assert meta_s["intensity_scale"] == pytest.approx(env.max())
    assert meta_k["intensity_scale"] == pytest.approx(comb.max())
    assert np.allclose(np.sort(amp_s * meta_s["intensity_scale"]), np.sort(env))
    # each header says what it is and names its twin
    assert meta_s["spectrum_kind"] == "sumecho"
    assert meta_k["spectrum_kind"] == "spikelets"
    assert meta_s["twin_file"] == "LAW_spikelets.csv"
    assert meta_k["twin_file"] == "LAW_sumecho.csv"
    # the shared QCPMG provenance travels in both, typed on the way back
    for meta in (meta_s, meta_k):
        assert meta["qcpmg_rotor_Hz"] == pytest.approx(16000.0)
        assert meta["spikelet_spacing_Hz"] == pytest.approx(391.9)
        assert meta["echo_period_pts"] == 293 and isinstance(meta["echo_period_pts"], int)
        assert meta["n_echoes"] == 38 and meta["split_offset_pts"] == 0
        assert meta["referenced"] is True
    assert meta_s["spectrum_mode"] == "absorption"
    assert meta_k["spectrum_mode"] == "magnitude"
    # the inputs were not normalised in place
    assert env.max() == pytest.approx(7.5) and comb.max() == pytest.approx(40.0)


def test_write_dataset_pair_raw_keeps_the_data_and_writes_unit_scale(tmp_path):
    x1, env, x2, comb, m_sum, m_spk = _twin_inputs()
    out = qcpmg.write_dataset_pair(tmp_path / "raw", x1, env, x2, comb,
                                   m_sum, m_spk, normalise=False)
    _, amp_s, meta_s = spectra.read_csv(out["sumecho"])
    _, amp_k, meta_k = spectra.read_csv(out["spikelets"])
    assert amp_s.max() == pytest.approx(7.5) and amp_k.max() == pytest.approx(40.0)
    assert meta_s["intensity_scale"] == 1.0 and meta_k["intensity_scale"] == 1.0


def test_both_files_load_through_load_any_and_the_field_reader(tmp_path):
    from larmor.qcpmg_fields import read_field_spectrum

    x1, env, x2, comb, m_sum, m_spk = _twin_inputs()
    out = qcpmg.write_dataset_pair(tmp_path / "LAW", x1, env, x2, comb,
                                   m_sum, m_spk)
    for kind, path in out.items():
        ppm, amp, recipe, summary, warnings = load_any(path)
        assert recipe["nucleus"] == "35Cl"
        assert recipe["larmor_frequency_MHz"] == pytest.approx(78.354)
        assert recipe["spin_rate_Hz"] == 0.0          # nothing to model
        assert amp.max() == pytest.approx(1.0)
        assert recipe["source_kind"] == "csv" and not warnings
    fs = read_field_spectrum(out["sumecho"])
    assert fs["larmor"] == pytest.approx(78.354) and fs["nucleus"] == "35Cl"
    assert fs["magnitude"] is False and fs["seed"].comb is False
    assert fs["rotor_Hz"] == pytest.approx(16000.0)
    fk = read_field_spectrum(out["spikelets"])
    assert fk["larmor"] == pytest.approx(78.354) and fk["magnitude"] is True
    # the spikelet twin is recognised as the comb it is, by name
    assert fk["seed"].comb is True
    assert "QCPMG spikelet spectrum" in fk["seed"].note
    assert fk["seed"].period_ppm == pytest.approx(5.0, rel=0.05)


def test_sample_label_strips_the_spikelet_suffix_too():
    from larmor.qcpmg_fields import sample_label

    assert sample_label({"sample": "Sample LAW3CL0CA · QCPMG spikelets (LB 1 Hz, magnitude)"},
                        "x/qcpmg_1_spikelets.csv") == "Sample LAW3CL0CA"
    assert sample_label({"sample": "Sample LAW3CL0CA · QCPMG sum echo (LB 75 Hz)"},
                        "x/qcpmg_1_sumecho.csv") == "Sample LAW3CL0CA"


# ------------------------------------------------------------ the dialog
pytest.importorskip("PySide6")


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _synthetic_dialog(qapp, tmp_path):
    from test_qcpmg_dialog import _dialog_with, _synthetic_qcpmg
    return _dialog_with(qapp, _synthetic_qcpmg(tmp_path))


def _pick(monkeypatch, path):
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getSaveFileName",
        staticmethod(lambda *a, **k: (str(path), "")))


def test_save_as_dataset_writes_both_files_normalised(qapp, tmp_path, monkeypatch):
    d = _synthetic_dialog(qapp, tmp_path)
    assert d.normSave.isChecked()                      # default: normalised
    _pick(monkeypatch, tmp_path / "qcpmg_1_sumecho.csv")
    d._save_dataset()
    sum_p = tmp_path / "qcpmg_1_sumecho.csv"
    spk_p = tmp_path / "qcpmg_1_spikelets.csv"
    assert sum_p.exists() and spk_p.exists()
    # both names are reported
    assert "qcpmg_1_sumecho.csv" in d.res.text() and "qcpmg_1_spikelets.csv" in d.res.text()

    ppm_s, amp_s, meta_s = spectra.read_csv(sum_p)
    ppm_k, amp_k, meta_k = spectra.read_csv(spk_p)
    assert amp_s.max() == pytest.approx(1.0) and amp_k.max() == pytest.approx(1.0)
    assert meta_s["intensity_scale"] == pytest.approx(float(np.abs(d._spec).max()))
    assert meta_k["intensity_scale"] == pytest.approx(float(np.abs(d._spk).max()))
    assert np.allclose(np.sort(ppm_s), np.sort(d._ppm))
    assert np.allclose(np.sort(ppm_k), np.sort(d._spk_ppm))
    assert meta_s["twin_file"] == spk_p.name and meta_k["twin_file"] == sum_p.name
    assert meta_s["spectrum_kind"] == "sumecho" and meta_k["spectrum_kind"] == "spikelets"
    # the spikelet spacing (1/tau_echo) and the echo-train bookkeeping
    for meta in (meta_s, meta_k):
        assert meta["spikelet_spacing_Hz"] == pytest.approx(d.periodHz.value())
        assert meta["echo_period_pts"] == d.period.value()
        assert meta["n_echoes"] == d.nEch.value()
        assert meta["split_offset_pts"] == d.offset.value()
        assert meta["nucleus"] == "35Cl" and meta["larmor_MHz"] == pytest.approx(78.0)
    assert meta_s["spectrum_mode"] == "absorption"
    assert meta_k["spectrum_mode"] == "magnitude"        # |FT(train)|
    assert "QCPMG spikelets" in meta_k["sample"]
    d.close()


def test_save_as_dataset_raw_when_unticked(qapp, tmp_path, monkeypatch):
    d = _synthetic_dialog(qapp, tmp_path)
    d.normSave.setChecked(False)
    _pick(monkeypatch, tmp_path / "raw.csv")
    d._save_dataset()
    _, amp_s, meta_s = spectra.read_csv(tmp_path / "raw_sumecho.csv")
    _, amp_k, meta_k = spectra.read_csv(tmp_path / "raw_spikelets.csv")
    assert meta_s["intensity_scale"] == 1.0 and meta_k["intensity_scale"] == 1.0
    assert amp_s.max() == pytest.approx(float(d._spec.max()))
    assert amp_k.max() == pytest.approx(float(d._spk.max()))
    d.close()


def test_saved_twins_open_in_the_multi_field_readers(qapp, tmp_path, monkeypatch):
    from larmor.qcpmg_fields import read_field_spectrum

    d = _synthetic_dialog(qapp, tmp_path)
    d.meta["spin_rate_Hz"] = 16000.0
    _pick(monkeypatch, tmp_path / "LAW.csv")
    d._save_dataset()
    fs = read_field_spectrum(str(tmp_path / "LAW_sumecho.csv"))
    assert fs["larmor"] == pytest.approx(78.0) and fs["nucleus"] == "35Cl"
    assert fs["rotor_Hz"] == pytest.approx(16000.0)
    assert fs["magnitude"] is False and fs["seed"].comb is False
    fk = read_field_spectrum(str(tmp_path / "LAW_spikelets.csv"))
    assert fk["larmor"] == pytest.approx(78.0)
    assert fk["magnitude"] is True and fk["seed"].comb is True
    assert "QCPMG spikelet spectrum" in fk["seed"].note
    d.close()


def test_dataset_meta_rejects_an_unknown_kind(qapp, tmp_path):
    d = _synthetic_dialog(qapp, tmp_path)
    with pytest.raises(ValueError):
        d.dataset_meta("envelope")
    assert d.dataset_meta()["spectrum_mode"] == "absorption"     # the default
    d.close()


# ------------------------------------------------------------ real data
def test_real_35cl_expno_writes_twin_datasets_that_round_trip(qapp, tmp_path):
    """The MagLab 35Cl EXPNO 1 (LAW3CL0CA, CNST7 period 293, 38 echoes):
    both files written, both maxima 1.0, intensity_scale = the raw maximum,
    the acquisition provenance in both headers, and both round-trip through
    load_any and the field reader (the spikelet twin flagged as a comb)."""
    from conftest import MAGLAB_35CL, require
    from larmor.desktop.qcpmg_dialog import QcpmgDialog
    from larmor.qcpmg_fields import read_field_spectrum

    fid = require(MAGLAB_35CL / "1" / "fid")
    d = QcpmgDialog(None, str(fid))
    assert d._spec is not None and d._spk is not None
    out = d.write_datasets(tmp_path / "qcpmg_1_sumecho.csv")
    assert os.path.basename(out["sumecho"]) == "qcpmg_1_sumecho.csv"
    assert os.path.basename(out["spikelets"]) == "qcpmg_1_spikelets.csv"
    raw_max = {"sumecho": float(np.abs(d._spec).max()),
               "spikelets": float(np.abs(d._spk).max())}
    for kind, path in out.items():
        ppm, amp, meta = spectra.read_csv(path)
        assert amp.max() == pytest.approx(1.0), kind
        assert meta["intensity_scale"] == pytest.approx(raw_max[kind], rel=1e-9), kind
        assert meta["spectrum_kind"] == kind
        assert meta["nucleus"] == "35Cl"
        assert meta["larmor_MHz"] == pytest.approx(78.354, abs=1e-3)
        assert meta["qcpmg_rotor_Hz"] == pytest.approx(16000.0)
        assert meta["echo_period_pts"] == 293 and meta["n_echoes"] == 38
        assert meta["spikelet_spacing_Hz"] == pytest.approx(533.3, abs=0.1)
        assert meta["referenced"] is True
        assert meta["carrier_ppm"] == pytest.approx(-102.82, abs=0.01)
        assert meta["sample"].startswith("Sample LAW3CL0CA")
        assert meta["twin_file"] == os.path.basename(
            out["spikelets" if kind == "sumecho" else "sumecho"])
        _, amp2, recipe, _, _ = load_any(path)
        assert recipe["nucleus"] == "35Cl"
        assert recipe["larmor_frequency_MHz"] == pytest.approx(78.354, abs=1e-3)
        assert recipe["spin_rate_Hz"] == 0.0
        assert amp2.max() == pytest.approx(1.0)
    fs = read_field_spectrum(out["sumecho"])
    assert fs["rotor_Hz"] == pytest.approx(16000.0) and fs["seed"].comb is False
    hi, lo = fs["seed"].hi_ppm, fs["seed"].lo_ppm
    assert lo < -113.0 < hi                    # the CT band straddles delta_CG
    fk = read_field_spectrum(out["spikelets"])
    assert fk["seed"].comb is True and fk["magnitude"] is True
    assert fk["seed"].period_ppm == pytest.approx(533.3 / 78.354, rel=0.05)
    d.close()
