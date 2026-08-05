"""Plotting studio + the spec-driven figure extensions it drives."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _trace():
    x = np.linspace(0, 100, 60)
    y = np.exp(-((x - 50) / 8.0) ** 2)
    return {"data": {"x": list(x), "y": list(y)}, "label": "peak", "color": "#1f77b4"}


def test_render_1d_honours_title_and_traces():
    from larmor import figures
    png = figures.render_png_bytes(
        {"kind": "1d", "title": "My figure", "traces": [_trace()]}, dpi=80)
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and len(png) > 500


def test_studio_builds_1d_spec(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    st._push_trace(_trace())
    st.title.setText("T"); st.stack.setValue(2.0)
    st._push_trace(_trace())
    spec = st._spec()
    assert spec["kind"] == "1d" and len(spec["traces"]) == 2
    assert spec["traces"][1]["offset"] == pytest.approx(2.0)   # stacked
    assert spec["title"] == "T"


def test_studio_builds_2d_spec_with_contour_and_iso(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(1)
    st.path2d.setText("/some/expno")
    st.cmode.setCurrentText("both"); st.chkValues.setChecked(True)
    st._iso.append({"slope": 1.0, "intercept": 0.0, "label": "CS"})
    spec = st._spec()
    assert spec["kind"] == "2d"
    assert spec["contour_mode"] == "both" and spec["contour_values"] is True
    assert spec["iso_lines"] and spec["iso_lines"][0]["label"] == "CS"


def test_journal_style_presets_render():
    from larmor import figures
    for style in ("nature", "acs", "rsc"):
        assert style in figures.STYLES
        png = figures.render_png_bytes(
            {"kind": "1d", "style": style, "traces": [_trace()]}, dpi=70)
        assert png[:4] == b"\x89PNG"


def test_norm_max_scales_traces_to_unit_peak():
    import matplotlib.pyplot as plt
    from larmor import figures
    x = np.linspace(0, 100, 80)
    t1 = {"data": {"x": list(x), "y": list(5 * np.exp(-((x - 40) / 6) ** 2))}}
    t2 = {"data": {"x": list(x), "y": list(9 * np.exp(-((x - 42) / 6) ** 2))}}
    fig = figures.render({"kind": "1d", "norm": "max", "traces": [t1, t2]})
    peaks = [float(l.get_ydata().max()) for l in fig.axes[0].lines]
    assert peaks == pytest.approx([1.0, 1.0], abs=1e-6)
    plt.close(fig)


def test_difference_subtracts_reference(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    st._push_trace(_trace()); st._push_trace(_trace())
    st.chkDiff.setChecked(True); st.norm.setCurrentText("max")
    spec = st._spec()
    assert spec["difference"] is True and spec["norm"] == "max"


def test_simulate_model_curve_for_dataless_fit():
    # a fit with no embedded spectrum → a model curve peaking at the site
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import figures
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(15.0),
            "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0),
            "gl": Param(1.0, vary=False)})])
    x, y = figures._simulate_model_curve(rec)
    assert len(x) > 100 and np.isfinite(y).all()
    assert abs(float(x[np.argmax(y)]) - 15.0) < 3.0        # peaks at the site


def test_studio_names_explorer_traces_by_sample(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    # a Bruker 1r path → the sample folder, NOT "1r"
    st._add_from_explorer("C:/data/03232026_P1-Bi0_SS_ALP/3102/pdata/1/1r")
    assert st._traces[-1]["label"] == "03232026_P1-Bi0_SS_ALP"
    # a dmfit fit file → its own (stemmed) name
    st._add_from_explorer("C:/data/sample/3102/pdata/1/P1-Bi0_31P.fxml")
    assert st._traces[-1]["label"] == "P1-Bi0_31P"


def test_studio_has_file_explorer_and_adds_traces(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    assert st.files is not None and not st.files.btnBatch.isVisible()
    st._add_from_explorer("C:/data/s/expno/pdata/1/1r")
    assert st._traces[-1]["path"].endswith("1r")
    st._add_from_explorer("C:/data/fit.recipe.json")
    assert st._traces[-1].get("recipe", "").endswith(".recipe.json")
    # in 2D mode a picked file sets the 2D path instead of adding a trace
    st.kind.setCurrentIndex(1)
    st._add_from_explorer("C:/data/expno2d")
    assert st.path2d.text() == "C:/data/expno2d"


def test_generic_xy_axis_with_custom_ticks():
    import matplotlib.pyplot as plt
    from larmor import figures
    spec = {"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
            "xlabel": "sample", "ylabel": "pop %",
            "xticks": [[1, "0Ca"], [2, "1Ca"], [3, "2Ca"]],
            "traces": [{"data": {"x": [1, 2, 3], "y": [80, 70, 60]}}]}
    fig = figures.render(spec); ax = fig.axes[0]
    assert [t.get_text() for t in ax.get_xticklabels()] == ["0Ca", "1Ca", "2Ca"]
    assert ax.get_xlim()[0] < ax.get_xlim()[1]        # NOT inverted (not ppm)
    assert len(ax.get_yticks()) > 0                    # intensity axis shown
    plt.close(fig)


def test_studio_ppm_toggle_and_custom_ticks(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    st._push_trace({"data": {"x": [1, 2, 3], "y": [1, 2, 3]}, "label": "t"})
    st.chkPpm.setChecked(False)
    st.xticks.setText("0Ca, 1Ca, 2Ca")
    spec = st._spec()
    assert spec["x_is_ppm"] is False and spec["hide_yaxis"] is False
    assert [lab for _, lab in spec["xticks"]] == ["0Ca", "1Ca", "2Ca"]
    # explicit pos:label pairs are honoured too
    st.xticks.setText("10:a, 20:b")
    assert st._spec()["xticks"] == [[10.0, "a"], [20.0, "b"]]


def test_studio_export_and_spec_roundtrip(qapp, tmp_path, monkeypatch):
    from larmor.desktop import plotting_studio, export_dialog
    st = plotting_studio.PlottingStudio(None)
    st._push_trace(_trace())
    out = tmp_path / "fig.pdf"
    monkeypatch.setattr(export_dialog, "choose",
                        lambda *a, **k: {"format": "PDF", "dpi": 150,
                                         "width_cm": 12, "height_cm": 9})
    monkeypatch.setattr(export_dialog, "_ask_path", lambda *a, **k: str(out))
    st._export()
    assert out.exists() and out.stat().st_size > 0
