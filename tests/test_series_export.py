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
