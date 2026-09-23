"""Series evolution plot, per-fit baseline and the reusable export options."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _result(positions=(15.0, 15.3, 14.7), ampsA=(100, 80, 120), ampsB=(50, 70, 40)):
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor.batchfit import BatchFitResult

    recs = []
    for k, pos in enumerate(positions):
        recs.append(Recipe(nucleus="11B", larmor_frequency_MHz=160.0,
                           spin_rate_Hz=0.0, sample=f"g{k}", sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(pos),
                "shift_fwhm_ppm": Param(6.0), "amplitude": Param(ampsA[k]),
                "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="B", params={
                "isotropic_chemical_shift_ppm": Param(2.0),
                "shift_fwhm_ppm": Param(3.0), "amplitude": Param(ampsB[k]),
                "gl": Param(1.0, vary=False)})]))
    return BatchFitResult(recipes=recs, labels=[f"g{k}" for k in range(len(recs))],
                          rmsd=[0.0] * len(recs), per_dataset=[], shared=(),
                          released=())


def test_series_options_lists_params_and_popfrac():
    from larmor.desktop.series_plot import series_options
    opts = series_options(_result())
    params = {(o["site"], o["param"], o["kind"]) for o in opts}
    assert (0, "isotropic_chemical_shift_ppm", "param") in params
    assert (0, "amplitude", "popfrac") in params
    assert not any(o["param"] == "gl" for o in opts)     # gl excluded


def test_series_values_track_positions_and_populations():
    from larmor.desktop.series_plot import series_values
    res = _result()
    pos = series_values(res, {"site": 0, "param": "isotropic_chemical_shift_ppm",
                              "kind": "param"})[0]
    assert list(pos) == pytest.approx([15.0, 15.3, 14.7])
    frac = series_values(res, {"site": 0, "param": "amplitude",
                               "kind": "popfrac"})[0]
    assert frac[0] == pytest.approx(100.0 / 150.0 * 100.0, abs=1e-6)


def test_series_dialog_builds(qapp):
    from larmor.desktop.series_plot import SeriesPlotDialog
    dlg = SeriesPlotDialog(None, _result())
    assert dlg.list.count() > 0
    assert len(dlg._selected()) == 1                     # first option preselected


def test_series_send_to_studio_spec_is_correct(qapp):
    # right-click a subplot → studio must get THIS parameter's values, an upright
    # axis and the sample names as x-ticks (not a mislabelled ppm spectrum)
    from larmor.desktop.series_plot import SeriesPlotDialog
    dlg = SeriesPlotDialog(None, _result())
    pop_spec = next(s for s in dlg._params if s["kind"] == "pop_integral")
    spec = dlg._studio_spec_for(pop_spec)
    assert spec["x_is_ppm"] is False and spec["hide_yaxis"] is False
    assert [lab for _, lab in spec["xticks"]] == dlg._labels   # sample names
    assert spec["ylabel"] == pop_spec["label"]
    assert spec["traces"] and len(spec["traces"][0]["data"]["y"]) == len(dlg._labels)


def _result_with_mc():
    """A batch result carrying a Monte-Carlo error detail on the amplitudes."""
    from larmor.batchfit import ParamError
    res = _result()
    # mark amplitude free so it counts as a fitted parameter, and attach MC errors
    detail = []
    for k, rec in enumerate(res.recipes):
        rec.sites[0].params["amplitude"].vary = True
        detail.append({(0, "amplitude"): ParamError(
            0, "amplitude", "s0.amplitude",
            float(rec.sites[0].params["amplitude"].value), 1.5 + k, (None, None),
            2.0)})
    res.error_detail = {"montecarlo": detail}
    res.error_method = "montecarlo"
    return res


def test_series_values_use_selected_error_method():
    from larmor.desktop.series_plot import series_values, error_methods
    res = _result_with_mc()
    assert error_methods(res) == ["none", "covariance", "montecarlo"]
    opt = {"site": 0, "param": "amplitude", "kind": "param"}
    _, mc = series_values(res, opt, "montecarlo")
    assert list(mc) == pytest.approx([1.5, 2.5, 3.5])      # the MC σ per spectrum
    _, none = series_values(res, opt, "none")
    assert np.isnan(none).all()                            # 'none' → no error bars


def test_population_integral_carries_an_error():
    from larmor.desktop.series_plot import population_integral, series_values
    res = _result_with_mc()          # MC error stored for site 0 only
    vals, errs = population_integral(res, "montecarlo")
    assert np.isfinite(vals).all()
    assert np.isfinite(errs[:, 0]).all() and (errs[:, 0] > 0).all()
    assert np.isnan(errs[:, 1]).all()          # site 1 has no MC error -> NaN
    # 'none' gives no error at all; the same call routed through series_values agrees
    _, e_none = population_integral(res, "none")
    assert np.isnan(e_none).all()
    v2, e2 = series_values(
        res, {"site": 0, "param": "population_pct", "kind": "pop_integral"},
        "montecarlo")
    assert v2 == pytest.approx(vals[:, 0]) and e2 == pytest.approx(errs[:, 0])


def test_series_values_for_family_and_ratio_kinds_and_dialog_items(qapp):
    """N3: Σ family and named-ratio series beside the sites -- values from
    quantify per spectrum, the error the stored exact per-spectrum value
    when an error run wrote one, the dialog lists them and adds a ratio
    subplot only when a ratio is defined; a SeqFitResult still works."""
    from larmor.batchfit import ParamError
    from larmor.desktop.series_plot import (SeriesPlotDialog, _param_specs,
                                            series_options, series_values)
    res = _result()
    for rec in res.recipes:
        rec.sites[0].family = "BO3"
        rec.sites[1].family = "BO4"
    opts = series_options(res)
    kinds = {(o["kind"], o.get("family"), o.get("ratio")) for o in opts}
    assert {("family", "BO3", None), ("family", "BO4", None), ("ratio", None, "N4")} <= kinds
    assert any(o["text"] == "Σ BO4: population % (integral)" for o in opts)
    fam_opt = {"kind": "family", "family": "BO4", "param": "family_pct", "site": None}
    rat_opt = {"kind": "ratio", "ratio": "N4", "param": "ratio", "site": None}
    pop = series_values(res, {"site": 1, "param": "population_pct",
                              "kind": "pop_integral"}, "none")[0]
    fam, e_none = series_values(res, fam_opt, "none")
    assert list(fam) == pytest.approx(list(pop)) and np.isnan(e_none).all()
    rat, _ = series_values(res, rat_opt, "none")
    assert list(rat) == pytest.approx(list(pop / 100.0))
    # a stored Monte-Carlo family error is used verbatim; without a stored
    # entry the fallback is quantify's own (None here: no stderr on the fit)
    res.error_detail = {"montecarlo": [
        {(-1, "family:BO4"): ParamError(-1, "family:BO4", "family.BO4", float(v),
                                        0.5 + k, (None, None), None)}
        for k, v in enumerate(fam)]}
    res.error_method = "montecarlo"
    _, e_mc = series_values(res, fam_opt, "montecarlo")
    assert list(e_mc) == pytest.approx([0.5, 1.5, 2.5])
    _, e_cov = series_values(res, fam_opt, "covariance")
    assert np.isnan(e_cov).all()
    _, e_rat = series_values(res, rat_opt, "montecarlo")     # no stored ratio error
    assert np.isnan(e_rat).all()

    dlg = SeriesPlotDialog(None, res)
    texts = [dlg.list.item(i).text() for i in range(dlg.list.count())]
    assert texts[-3:] == ["Σ BO3", "Σ BO4", "N4"]
    assert "ratio" in dlg._subplots and _param_specs(res)[-1]["kind"] == "ratio"
    dlg.list.clearSelection()
    dlg.list.item(dlg.list.count() - 2).setSelected(True)          # Σ BO4
    assert dlg._selected() == ["f:BO4"]
    pop_spec = next(s for s in dlg._params if s["kind"] == "pop_integral")
    spec = dlg._studio_spec_for(pop_spec)
    assert [tr["label"] for tr in spec["traces"]] == ["Σ BO4"]
    assert spec["traces"][0]["data"]["y"] == pytest.approx(list(fam))
    assert "yerr" in spec["traces"][0]["data"]                   # the MC errors
    amp_spec = next(s for s in dlg._params if s["param"] == "amplitude")
    assert dlg._studio_spec_for(amp_spec)["traces"] == []        # families: population only
    dlg.list.item(dlg.list.count() - 1).setSelected(True)          # + N4
    ratio_spec = next(s for s in dlg._params if s["kind"] == "ratio")
    assert [tr["label"] for tr in dlg._studio_spec_for(ratio_spec)["traces"]] == ["N4"]
    dlg._draw()
    assert dlg._sel_key("f:BO4", pop_spec) == "family:BO4"
    assert dlg._sel_key("r:N4", ratio_spec) == "ratio:N4"
    assert dlg._sel_key(0, amp_spec) == "s0:amplitude"
    # a SeqFitResult (no error_detail / shared / released) still works
    from larmor.seqfit import SeqFitResult
    seq = SeqFitResult(recipes=res.recipes, labels=res.labels, rmsd=res.rmsd,
                       per_dataset=[], history=[], passes=1, propagated=())
    v, e = series_values(seq, fam_opt, "montecarlo")
    assert list(v) == pytest.approx(list(fam)) and np.isnan(e).all()
    dlg_seq = SeriesPlotDialog(None, seq)
    assert "ratio" in dlg_seq._subplots
    # untagged: no Σ items, no ratio subplot, options as before
    plain = SeriesPlotDialog(None, _result())
    assert not any(plain.list.item(i).text().startswith("Σ") for i in range(plain.list.count()))
    assert "ratio" not in plain._subplots
    assert not any(o["kind"] in ("family", "ratio") for o in series_options(_result()))


def test_series_dialog_error_selector_and_studio_yerr(qapp):
    from larmor.desktop.series_plot import SeriesPlotDialog
    dlg = SeriesPlotDialog(None, _result_with_mc())
    methods = [dlg.errSel.itemData(i) for i in range(dlg.errSel.count())]
    assert methods == ["none", "covariance", "montecarlo"]
    assert dlg._error_method() == "montecarlo"             # defaults to computed
    amp = next(s for s in dlg._params if s["param"] == "amplitude")
    spec = dlg._studio_spec_for(amp)
    # the studio spec carries the selected error as per-point yerr
    assert "yerr" in spec["traces"][0]["data"]
    assert len(spec["traces"][0]["data"]["yerr"]) == len(dlg._labels)


def test_estimate_baseline_recovers_a_slope():
    from larmor.desktop.batchfit_dialog import estimate_baseline
    x = np.linspace(-20, 60, 600)
    slope = 0.5 * x + 3.0
    peak = 80.0 * np.exp(-0.5 * ((x - 15.0) / 4.0) ** 2)
    y = slope + peak
    base = estimate_baseline(x, y, "Polynomial", order=1)
    # away from the peak the polynomial baseline tracks the true slope
    off = np.abs(x - 15.0) > 20
    assert np.max(np.abs(base[off] - slope[off])) < 3.0


def test_batch_snr_helper():
    from larmor.desktop.batchfit_dialog import _snr
    x = np.linspace(-20, 60, 600)
    clean = 100.0 * np.exp(-0.5 * ((x - 15) / 4) ** 2)
    noisy = clean + np.random.default_rng(0).normal(0, 20, x.size)
    assert _snr(clean) > _snr(noisy)               # more noise -> lower S/N
    assert _snr(clean) > 50


def test_export_options_values(qapp):
    from larmor.desktop.export_dialog import ExportOptions, CM_PER_IN
    dlg = ExportOptions(None, ["PNG", "PDF", "SVG"], dpi=600, width_cm=10, height_cm=8)
    v = dlg.values()
    assert v["format"] == "PNG" and v["dpi"] == 600
    assert v["width_cm"] == 10 and v["height_cm"] == 8
    assert CM_PER_IN == pytest.approx(2.54)


def test_export_matplotlib_writes_file(qapp, tmp_path, monkeypatch):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from larmor.desktop import export_dialog

    fig = plt.figure(); fig.add_subplot(111).plot([0, 1], [0, 1])
    out = tmp_path / "fig.png"
    monkeypatch.setattr(export_dialog, "choose",
                        lambda *a, **k: {"format": "PNG", "dpi": 150,
                                         "width_cm": 12, "height_cm": 9})
    monkeypatch.setattr(export_dialog, "_ask_path", lambda *a, **k: str(out))
    path = export_dialog.export_matplotlib(None, fig, "fig")
    assert path and out.exists() and out.stat().st_size > 0
