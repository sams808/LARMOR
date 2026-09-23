"""Sample scanner (auto-identify) and multi-format fit export."""
from pathlib import Path

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
    # the acquisition facts the session inventory reads come off the same
    # text parse of acqus (NS, DATE, the D array) and the pdata listing
    assert all(e.ns > 0 and e.date > 0 for e in exps)
    assert any(e.d1_s > 0 for e in exps)
    assert all(e.n_procs >= 1 for e in exps)
    assert all(e.date_iso.startswith("20") for e in exps)


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


def _export_text_reference(recipe, exp_ppm, exp_amp):
    """export_text's body before curve_table was factored out of it -- the
    golden for the byte-identity check, portable across mrsimulator builds."""
    from larmor import engine
    x, total, per_site = engine.simulate(recipe, exp_ppm=exp_ppm)
    exp_on_x = np.interp(x, exp_ppm, exp_amp)
    cols = [x, exp_on_x, total, exp_on_x - total, *per_site]
    labels = ["ppm", "experiment", "model", "residual"] + \
        [s.label or f"site{i}" for i, s in enumerate(recipe.sites)]
    lines = ["# " + "\t".join(labels)]
    for row in zip(*cols):
        lines.append("\t".join(f"{v:.6g}" for v in row))
    return "\n".join(lines) + "\n"


def test_curve_table_experiment_axis_is_exact_and_export_text_is_unchanged(tmp_path):
    r = _recipe()
    ppm_asc = np.linspace(-50, 150, 800)
    amp_asc = np.exp(-((ppm_asc - 60) / 12) ** 2) * 900
    ppm, amp = ppm_asc[::-1].copy(), amp_asc[::-1].copy()        # DESCENDING input

    labels, cols = export.curve_table(r, ppm, amp, on_experiment_axis=True,
                                      site_prefix=True)
    assert labels == ["ppm", "experiment", "model", "residual", "s0_AlO4", "s1_imp"]
    x, experiment, model, residual = cols[:4]
    assert np.all(np.diff(x) > 0)
    order = np.argsort(ppm, kind="stable")
    assert np.array_equal(x, ppm[order]) and np.array_equal(experiment, amp[order])
    assert np.allclose(residual, experiment - model)
    assert np.allclose(model, cols[4] + cols[5], atol=1e-9 * model.max())
    assert x.size == ppm.size                       # never the model's kernel axis

    # export_text (model axis, .6g) is byte-identical to its former body
    txt = export.export_text(r, ppm_asc, amp_asc, tmp_path / "o.txt")
    assert txt == _export_text_reference(r, ppm_asc, amp_asc)
    head = txt.splitlines()[0]
    assert "AlO4" in head and "imp" in head


def _read_curves(path):
    """A curves CSV as a structured array: the '# key=value' lines dropped
    first, because genfromtxt(names=True) would otherwise take the first
    commented line ('# LARMOR curves') as the header. pandas reads the
    same file with comment="#"."""
    rows = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines()
            if ln and not ln.startswith("#")]
    return np.genfromtxt(rows, delimiter=",", names=True)


def test_export_curves_csv_roundtrips_exactly_and_reopens_in_larmor(tmp_path):
    from larmor.loader import load_any

    r = _recipe()
    ppm = np.linspace(-50, 150, 300)[::-1].copy()
    amp = (np.exp(-((ppm - 60) / 12) ** 2) * 900
           + np.random.default_rng(1).normal(0, 3.0, ppm.size))
    p = tmp_path / "glass_curves.csv"
    out = export.export_curves_csv(r, ppm, amp, p, header={"fit_window_ppm": "150.0,-50.0"})
    lines = p.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# LARMOR curves"
    assert "# sample=glass" in lines and "# nucleus=27Al" in lines
    assert "# larmor_MHz=195.483" in lines and "# fit_window_ppm=150.0,-50.0" in lines
    hdr = next(ln for ln in lines if not ln.startswith("#"))
    assert hdr == "ppm,experiment,model,residual,s0_AlO4,s1_imp"
    assert out["columns"] == hdr.split(",")

    tab = _read_curves(p)
    order = np.argsort(ppm, kind="stable")
    assert np.array_equal(tab["ppm"], ppm[order])              # repr round-trip
    assert np.array_equal(tab["experiment"], amp[order])
    assert np.allclose(tab["model"], out["model"], rtol=1e-8)
    assert np.allclose(tab["model"], tab["s0_AlO4"] + tab["s1_imp"],
                       atol=1e-6 * tab["model"].max())

    # File > Open reopens the curves file as the fitted spectrum
    ppm2, amp2, rec2, _meta, _w = load_any(p)
    assert np.array_equal(ppm2, ppm[order]) and np.array_equal(amp2, amp[order])
    assert rec2["nucleus"] == "27Al" and rec2["sample"] == "glass"
    assert rec2["larmor_frequency_MHz"] == pytest.approx(195.483)

    # a raw (pre-baseline) trace goes right after the experiment
    raw = amp + 5.0
    p2 = tmp_path / "glass_raw_curves.csv"
    out2 = export.export_curves_csv(r, ppm, amp, p2, raw=raw)
    hdr2 = next(ln for ln in p2.read_text(encoding="utf-8").splitlines()
                if not ln.startswith("#"))
    assert hdr2.split(",")[:3] == ["ppm", "experiment", "experiment_raw"]
    tab2 = _read_curves(p2)
    assert np.array_equal(tab2["experiment_raw"], raw[order])
    assert np.array_equal(out2["experiment_raw"], raw[order])
