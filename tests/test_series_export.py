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


def _result(positions=(15.0, 15.3, 14.7), ampsA=(100, 80, 120), ampsB=(50, 70, 40),
            samples=None):
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor.batchfit import BatchFitResult

    samples = list(samples) if samples else [f"g{k}" for k in range(len(positions))]
    recs = []
    for k, pos in enumerate(positions):
        recs.append(Recipe(nucleus="11B", larmor_frequency_MHz=160.0,
                           spin_rate_Hz=0.0, sample=samples[k], sites=[
            SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(pos),
                "shift_fwhm_ppm": Param(6.0), "amplitude": Param(ampsA[k]),
                "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="B", params={
                "isotropic_chemical_shift_ppm": Param(2.0),
                "shift_fwhm_ppm": Param(3.0), "amplitude": Param(ampsB[k]),
                "gl": Param(1.0, vary=False)})]))
    return BatchFitResult(recipes=recs, labels=samples,
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


# ---------------------------------------------------------------- series table (N6)
def _series(names, groups=None, column=None):
    """A SeriesTable over fake paths; ``column`` = (key, label, values, errs)."""
    from larmor.series_table import SeriesColumn, SeriesTable
    t = SeriesTable.from_paths([f"C:/d/{n}.csv" for n in names], names=list(names),
                               groups=groups)
    if column:
        key, label, vals, errs = column
        t.add_column(SeriesColumn(key, label, f"{key}_sd" if errs else None),
                     {k: v for k, v in enumerate(vals)})
        if errs:
            for r, e in zip(t.rows, errs):
                r.values[f"{key}_sd"] = e
    return t


def test_series_dialog_without_series_keeps_todays_behaviour(qapp):
    from larmor.desktop.series_plot import SeriesPlotDialog
    dlg = SeriesPlotDialog(None, _result())                  # no series kwarg
    assert dlg.xSel.count() == 1 and dlg.xSel.currentData() is None
    assert dlg.chkAvg.isHidden() and not dlg.chkOls.isChecked()
    pop = next(s for s in dlg._params if s["kind"] == "pop_integral")
    pts = dlg._points(pop, 0)
    assert pts["names"] == dlg._labels and pts["x"].tolist() == [1.0, 2.0, 3.0]
    assert pts["xerr"] is None and pts["n"].tolist() == [1, 1, 1]
    spec = dlg._studio_spec_for(pop)
    assert [lab for _, lab in spec["xticks"]] == dlg._labels
    assert spec["xlabel"] == "sample" and "xerr" not in spec["traces"][0]["data"]
    assert len(spec["traces"]) == 1
    assert dlg.msg.text() == ""


def test_series_dialog_numeric_x_axis_gives_xerr_and_column_label_in_studio_spec(
        qapp, monkeypatch):
    from larmor.desktop import export_dialog
    from larmor.desktop.series_plot import SeriesPlotDialog
    ser = _series(["g0", "g1", "g2"],
                  column=("P2O5", "analysed P2O5 (mol%)", [5.0, 1.0, 8.0], [.1, .1, .2]))
    dlg = SeriesPlotDialog(None, _result(), series=ser)
    assert [dlg.xSel.itemText(i) for i in range(dlg.xSel.count())] == \
        ["series order (names)", "analysed P2O5 (mol%)"]
    assert dlg.chkAvg.isHidden()                            # no replicate group
    dlg.xSel.setCurrentIndex(1)
    amp = next(s for s in dlg._params if s["param"] == "amplitude")
    spec = dlg._studio_spec_for(amp)
    assert spec["xlabel"] == "analysed P2O5 (mol%)" and "xticks" not in spec
    d = spec["traces"][0]["data"]
    assert d["x"] == [1.0, 5.0, 8.0]                         # sorted by x
    assert d["y"] == [80.0, 100.0, 120.0]                    # g1, g0, g2 amplitudes
    assert d["xerr"] == [.1, .1, .2]
    assert len(spec["traces"]) == 1
    dlg.chkOls.setChecked(True)
    spec = dlg._studio_spec_for(amp)
    assert len(spec["traces"]) == 2
    line = spec["traces"][1]
    assert line["linestyle"] == "--" and "slope" in line["label"] and "r =" in line["label"]
    assert line["data"]["x"] == [1.0, 8.0]
    pw = dlg._subplots["amplitude"]
    assert pw.getAxis("bottom").labelText == "analysed P2O5 (mol%)"
    # the figure export carries the numeric axis and the x bars
    got = {}
    monkeypatch.setattr(export_dialog, "export_matplotlib",
                        lambda parent, fig, name: got.setdefault("fig", fig))
    dlg._export_fig()
    from matplotlib.container import ErrorbarContainer
    ax0 = got["fig"].axes[0]
    conts = [c for c in ax0.containers if isinstance(c, ErrorbarContainer)]
    assert ax0.get_xlabel() == "analysed P2O5 (mol%)" and conts and conts[0].has_xerr
    dlg.xSel.setCurrentIndex(0)
    spec = dlg._studio_spec_for(amp)
    assert [lab for _, lab in spec["xticks"]] == ["g0", "g1", "g2"]
    assert spec["traces"][0]["data"]["y"] == [100.0, 80.0, 120.0]
    # a table that does not pair with the result is ignored, not half-applied
    bad = SeriesPlotDialog(None, _result(), series=_series(["a", "b"]))
    assert bad._series is None and bad.xSel.count() == 1


def test_series_dialog_averages_replicates_with_spread_when_ticked(qapp, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog
    from larmor.desktop.series_plot import SeriesPlotDialog, population_integral
    res = _result(samples=["g0", "g0 (2)", "g1"])
    ser = _series(["g0", "g0 (2)", "g1"], groups=["g0", "g0", "g1"],
                  column=("Na2O", "Na2O (mol%)", [20.0, 20.0, 30.0], [.5, .5, .5]))
    dlg = SeriesPlotDialog(None, res, series=ser)
    assert not dlg.chkAvg.isHidden() and dlg.chkAvg.isChecked()
    pop = next(s for s in dlg._params if s["kind"] == "pop_integral")
    frac = population_integral(res)[0][:, 0]
    pts = dlg._points(pop, 0)
    assert pts["names"] == ["g0", "g1"] and len(pts["y"]) == 2
    assert pts["y"][0] == pytest.approx(np.mean(frac[:2]))
    assert pts["yerr"][0] == pytest.approx(np.std(frac[:2], ddof=1))
    assert pts["y"][1] == pytest.approx(frac[2]) and np.isnan(pts["yerr"][1])
    assert pts["n"].tolist() == [2, 1]
    dlg.chkAvg.setChecked(False)
    assert len(dlg._points(pop, 0)["y"]) == 3
    dlg.chkAvg.setChecked(True)
    out = tmp_path / "series.csv"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(out), "")))
    dlg._export_csv()
    lines = out.read_text(encoding="utf-8").splitlines()
    head = lines[0].split(",")
    assert head[:2] == ["spectrum", "n"] and len(lines) == 3
    assert lines[1].split(",")[:2] == ["g0", "2"] and lines[2].split(",")[:2] == ["g1", "1"]
    assert any(h.endswith("±(sd; covariance when n=1)") for h in head)
    dlg.xSel.setCurrentIndex(1)                             # the Na2O column
    dlg._export_csv()
    head = out.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert head[:4] == ["spectrum", "x (Na2O (mol%))", "x ±", "n"]
    row = out.read_text(encoding="utf-8").splitlines()[1].split(",")
    assert row[:4] == ["g0", "20", "0.5", "2"]               # identical x -> RMS ±


def test_series_dialog_species_bar_and_dust_export(qapp, tmp_path, monkeypatch):
    import csv
    from PySide6.QtWidgets import QFileDialog
    from larmor.desktop.series_plot import SeriesPlotDialog, population_integral
    res = _result()
    ser = _series(["g0", "g1", "g2"], column=("Na2O_mol", "analysed Na2O_mol",
                                              [20.0, 25.0, 30.0], None))
    dlg = SeriesPlotDialog(None, res, series=ser)
    pop = population_integral(res)[0]
    spec = dlg._species_bar_spec()
    assert spec["kind"] == "species_bar" and spec["categories"] == ["g0", "g1", "g2"]
    assert [s["label"] for s in spec["series"]] == ["s0 A", "s1 B"]
    for j, s in enumerate(spec["series"]):
        assert s["values"] == pytest.approx(pop[:, j].tolist())
    assert "xlabel" not in spec
    dlg.xSel.setCurrentIndex(1)
    spec = dlg._species_bar_spec()
    assert spec["xlabel"] == "in order of analysed Na2O_mol"
    dlg.xSel.setCurrentIndex(0)
    out = tmp_path / "dust.csv"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(out), "")))
    assert dlg._selected() == [0]
    dlg._export_dust()
    with open(out, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["Sample", "Na2O", "N4_measured", "N4_measured_err"]
    assert [r[0] for r in rows[1:]] == ["g0", "g1", "g2"]
    assert [float(r[1]) for r in rows[1:]] == [20.0, 25.0, 30.0]
    for k, r in enumerate(rows[1:]):
        assert float(r[2]) == pytest.approx(pop[k, 0] / 100.0, abs=1e-6)   # .6g cells
        assert r[3] == ""                                    # no error computed
    msg = dlg.msg.text()
    assert "dust.csv" in msg and "N4_measured = s0" in msg and "Ignore" in msg
    assert "Na2O" in msg and "skipped" not in msg
    # both sites selected: N4 = their sum; no series table: Sample + N4 only
    dlg.list.item(1).setSelected(True)
    dlg._export_dust(str(out))
    with open(out, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert float(rows[1][2]) == pytest.approx((pop[0, 0] + pop[0, 1]) / 100.0, abs=1e-6)
    assert "N4_measured = s0 + s1" in dlg.msg.text()
    bare = SeriesPlotDialog(None, res)
    bare._export_dust(str(out))
    with open(out, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["Sample", "N4_measured", "N4_measured_err"]
    assert [r[0] for r in rows[1:]] == ["g0", "g1", "g2"]
    assert "none — join a composition CSV" in bare.msg.text()
    bare.list.clearSelection()
    bare._export_dust(str(out))
    assert "select the site" in bare.msg.text()


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
