"""Sample scanner (auto-identify) and multi-format fit export."""
import numpy as np
import pytest

from larmor.io import export, fxmla, scan
from larmor.recipe import Param, Recipe, SiteModel

from conftest import BRUKER_1R, require

SAMPLE = None
try:
    from conftest import BRUKER_1R as _b
    SAMPLE = _b.parent.parent.parent      # .../2702/pdata/1/1r → the sample
except Exception:
    pass


# ---------------------------------------------------------------- scanner
def test_scan_sample_identifies_experiments():
    require(BRUKER_1R)
    sample = BRUKER_1R.parents[3]   # the dated sample folder
    exps = scan.scan_sample(sample)
    assert len(exps) >= 2
    nuclei = {e.nucleus for e in exps}
    assert "27Al" in nuclei
    # each has a best openable target and a human kind
    for e in exps:
        assert e.openable is None or e.openable.endswith(
            ("1r", "2rr", "fid", "ser"))
        assert e.kind and e.ndim in (1, 2)


def test_classify_pulse_programs():
    assert scan._classify("zg") == "Single pulse"
    assert scan._classify("hahnecho.nmrfam") == "Hahn echo"
    assert scan._classify("satrect1") == "Saturation recovery (T1)"
    assert scan._classify("mp3q") == "MQMAS"
    assert scan._classify("cpmg_something") == "CPMG (T2)"
    assert scan._classify("weird_pp") == "weird_pp"


def test_is_sample_folder():
    require(BRUKER_1R)
    sample = BRUKER_1R.parents[3]
    assert scan.is_sample_folder(sample)
    # a month folder is NOT a sample (its children are samples, not EXPNOs)
    month = sample.parent
    assert not scan.is_sample_folder(month)


def test_list_dir_flags_samples_and_expnos():
    require(BRUKER_1R)
    sample = BRUKER_1R.parents[3]
    entries = scan.list_dir(sample)
    assert any(e.is_expno for e in entries)
    month = scan.list_dir(sample.parent)
    assert any(e.is_sample for e in month)


# ---------------------------------------------------------------- export
def _recipe():
    return Recipe(sample="glass", nucleus="27Al", larmor_frequency_MHz=195.483,
                  spin_rate_Hz=20000,
                  sites=[SiteModel(model="czjzek", label="AlO4", params={
                      "isotropic_chemical_shift_ppm": Param(60.0),
                      "sigma_Cq_MHz": Param(2.0, stderr=0.1),
                      "shift_fwhm_ppm": Param(10.0),
                      "amplitude": Param(1000.0, min=0)}),
                         SiteModel(model="gauss_lor", label="imp", params={
                      "isotropic_chemical_shift_ppm": Param(0.0),
                      "shift_fwhm_ppm": Param(5.0),
                      "amplitude": Param(100.0),
                      "gl": Param(1.0, vary=False)})])


def test_export_text_columns(tmp_path):
    r = _recipe()
    ppm = np.linspace(-50, 150, 800)
    amp = np.exp(-((ppm - 60) / 12) ** 2) * 900
    txt = export.export_text(r, ppm, amp, tmp_path / "o.txt")
    head = txt.splitlines()[0]
    assert head.startswith("# ppm") and "experiment" in head and "model" in head
    assert "AlO4" in head and "imp" in head
    # a data row has one column per label
    ncol = len(head.lstrip("# ").split("\t"))
    assert len(txt.splitlines()[1].split("\t")) == ncol


def test_export_csv_params(tmp_path):
    csv = export.export_csv_params(_recipe(), tmp_path / "o.csv")
    assert csv.splitlines()[0] == "line,model,parameter,value,stderr,min,max,link"
    # components are lettered A, B…
    assert any(row.startswith("A,czjzek,") for row in csv.splitlines())
    assert any(row.startswith("B,gauss_lor,") for row in csv.splitlines())


def test_export_fxmla_roundtrips_through_our_parser(tmp_path):
    r = _recipe()
    ppm = np.linspace(-50, 150, 1200)
    amp = np.exp(-((ppm - 60) / 12) ** 2) * 900
    fx = tmp_path / "out.fxmla"
    export.export_fxmla(r, ppm, amp, fx)

    dm = fxmla.read(fx)
    assert dm.fit_mode == "Fit 1D"
    assert len(dm.dimensions[0].lines) == 2
    models = [ln.model_name for ln in dm.dimensions[0].lines]
    assert "CzSimple" in models and "Gaus/Lor" in models

    back, warns = fxmla.to_recipe(dm)
    cz = [s for s in back.sites if s.model == "czjzek"][0]
    # σ = sCZ_CQ/2 survives the write→read round trip
    assert cz.params["sigma_Cq_MHz"].value == pytest.approx(2.0, rel=1e-4)
    assert cz.params["isotropic_chemical_shift_ppm"].value == pytest.approx(60.0)
    # the embedded spectrum peaks where it should
    assert dm.spectrum is not None
    assert dm.spectrum.ppm[dm.spectrum.amplitude.argmax()] == pytest.approx(
        60.0, abs=1.0)


def test_export_fxmla_writes_unknown_models_as_gausslor(tmp_path):
    r = Recipe(sample="s", nucleus="27Al", larmor_frequency_MHz=195.5,
               sites=[SiteModel(model="quad_ct", label="q", params={
                   "isotropic_chemical_shift_ppm": Param(30.0),
                   "Cq_MHz": Param(3.0), "eta": Param(0.2),
                   "shift_fwhm_ppm": Param(2.0), "amplitude": Param(1.0)})])
    ppm = np.linspace(-50, 100, 500)
    fx = tmp_path / "q.fxmla"
    text = export.export_fxmla(r, ppm, np.zeros_like(ppm), fx)
    assert "Gaus/Lor" in text and "exported as Gaus/Lor" in text
    dm = fxmla.read(fx)                      # still parses
    assert len(dm.dimensions[0].lines) == 1


def test_formats_registry():
    assert set(export.FORMATS) >= {"text (.txt)", "parameters (.csv)",
                                   "LARMOR recipe (.json)", "dmfit (.fxmla)"}
