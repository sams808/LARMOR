"""Plotting studio + the spec-driven figure extensions it drives."""
import os
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtWidgets import QApplication, QTableWidgetItem  # noqa: E402


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


def test_studio_2d_spec_carries_nucleus_and_fit_overlay(qapp, tmp_path):
    from larmor.desktop.plotting_studio import PlottingStudio
    from larmor.recipe import Recipe, SiteModel, Param
    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=195.483, sites=[
        SiteModel(model="czjzek", label="Al", params={
            "isotropic_chemical_shift_ppm": Param(60.0),
            "sigma_Cq_MHz": Param(2.0), "shift_fwhm_ppm": Param(6.0),
            "amplitude": Param(1.0)})])
    p = tmp_path / "fit.recipe.json"
    rec.save(p)

    st = PlottingStudio(None)
    st.kind.setCurrentIndex(1)
    st.path2d.setText("/some/expno")
    st.nuc2d.setText("27Al")
    st.fitRecipe2d.setText(str(p))
    st.mqmasMethod.setCurrentText("5QMAS")
    spec = st._spec()
    assert spec["nucleus"] == "27Al"
    assert spec["fit_recipe"] == str(p)
    assert spec["mqmas_method"] == "5QMAS"

    st2 = PlottingStudio(None)
    st2._apply_spec(spec)
    assert st2.nuc2d.text() == "27Al"
    assert st2.fitRecipe2d.text() == str(p)
    assert st2.mqmasMethod.currentText() == "5QMAS"


def test_studio_reference_line_dialog_computes_cs_and_qis_slopes(qapp):
    from larmor.desktop.plotting_studio import _ReferenceLineDialog
    dlg = _ReferenceLineDialog(None, "27Al", 195.483)
    dlg.kind.setCurrentIndex(1)                # CS axis
    dlg._compute()
    cs_slope = dlg.slope.value()
    assert cs_slope != 1.0 and dlg.intercept.value() == pytest.approx(0.0)
    assert dlg.labelEdit.text() == "CS axis"

    dlg2 = _ReferenceLineDialog(None, "27Al", 195.483)
    dlg2.kind.setCurrentIndex(2)                # QIS axis
    dlg2.anchor.setValue(60.0)
    dlg2._compute()
    assert dlg2.slope.value() != 1.0
    # the QIS line passes through the anchor's OWN point on the CS axis
    # (loose tolerance: slope/intercept round-trip through 4/3-decimal
    # QDoubleSpinBox display precision, not exact float storage)
    f1_at_anchor = cs_slope * 60.0
    assert dlg2.slope.value() * 60.0 + dlg2.intercept.value() == pytest.approx(
        f1_at_anchor, abs=0.01)

    dlg3 = _ReferenceLineDialog(None)
    dlg3.kind.setCurrentIndex(0)                # Manual: fields disabled, no crash
    assert not dlg3.nucleus.isEnabled()
    dlg3._compute()                             # a no-op for Manual
    assert dlg3.slope.value() == 1.0


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


def test_render_1d_draws_error_bars_from_trace_yerr():
    import matplotlib.pyplot as plt
    from larmor import figures
    spec = {"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
            "traces": [{"data": {"x": [1, 2, 3], "y": [10, 12, 9],
                                 "yerr": [0.5, 0.8, 0.4]}, "label": "s0"}]}
    fig = figures.render(spec)
    # an ErrorbarContainer is present on the axes
    assert fig.axes[0].containers, "no error bars drawn"
    plt.close(fig)


def test_render_1d_axis_customisation_options():
    import matplotlib.pyplot as plt
    from larmor import figures
    base = {"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
            "traces": [{"data": {"x": [1, 2, 3], "y": [1, 2, 3]}, "label": "a"}]}
    # y-limits (order-independent), legend hidden, grid on
    fig = figures.render({**base, "ylim": (5, 0), "legend_loc": "none",
                          "grid": True})
    ax = fig.axes[0]
    assert ax.get_ylim() == pytest.approx((0.0, 5.0))
    assert ax.get_legend() is None                      # 'none' hides it
    plt.close(fig)
    # tick step + minor ticks + font-size override all render
    fig2 = figures.render({**base, "ytick_step": 0.5, "minor_ticks": True,
                           "tick_direction": "out", "font_size": 6})
    assert len(fig2.axes[0].get_yticks()) > 3
    plt.close(fig2)


def test_studio_spec_carries_axis_customisation(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    st._push_trace({"data": {"x": [1, 2, 3], "y": [1, 2, 3]}, "label": "t"})
    st.chkPpm.setChecked(False)
    st.chkYlim.setChecked(True); st.ylo.setValue(0); st.yhi.setValue(10)
    st.legloc.setCurrentText("upper left"); st.legncol.setValue(2)
    st.chkGrid.setChecked(True); st.chkMinor.setChecked(True)
    st.xstep.setValue(1.0); st.fontsz.setValue(9)
    spec = st._spec()
    assert spec["ylim"] == (0.0, 10.0)
    assert spec["legend_loc"] == "upper left" and spec["legend_ncol"] == 2
    assert spec["grid"] is True and spec["minor_ticks"] is True
    assert spec["xtick_step"] == pytest.approx(1.0) and spec["font_size"] == 9
    # round-trips through _apply_spec
    st2 = PlottingStudio(None); st2._apply_spec(spec)
    assert st2.legloc.currentText() == "upper left"
    assert st2.chkYlim.isChecked() and st2.chkGrid.isChecked()


def test_trace_editor_exposes_error_bar_style_only_when_data_has_yerr(qapp):
    from larmor.desktop.plotting_studio import _TraceEditor

    no_err = _TraceEditor(None, {"data": {"x": [1, 2], "y": [1, 2]}})
    assert no_err.errVisible is None
    assert "err_visible" not in no_err.values()

    with_err = _TraceEditor(
        None, {"data": {"x": [1, 2], "y": [1, 2], "yerr": [0.1, 0.1]}})
    assert with_err.errVisible is not None
    with_err.errVisible.setChecked(False)
    with_err.errWidth.setValue(2.0)
    with_err.errCap.setValue(5.0)
    with_err._err_color = "#ff0000"
    v = with_err.values()
    assert v["err_visible"] is False
    assert v["err_width"] == pytest.approx(2.0)
    assert v["err_capsize"] == pytest.approx(5.0)
    assert v["err_color"] == "#ff0000"


def test_render_1d_error_bar_style_is_controllable():
    import matplotlib.pyplot as plt
    from larmor import figures
    trace = {"data": {"x": [1, 2, 3], "y": [20, 25, 30],
                      "yerr": [0.1, 0.2, 0.15]}, "label": "pop%"}
    shown = figures.render({"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
                            "traces": [{**trace, "err_visible": True,
                                       "err_color": "#ff0000", "err_width": 2.5,
                                       "err_capsize": 6.0}]})
    assert len(shown.axes[0].containers) == 1
    plt.close(shown)
    hidden = figures.render({"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
                             "traces": [{**trace, "err_visible": False}]})
    assert len(hidden.axes[0].containers) == 0        # explicitly hidden
    plt.close(hidden)


def _fit(tmp_path, sample, pos=10.0, amp=80.0, with_data=True):
    from larmor import engine
    from larmor.recipe import Recipe, SiteModel, Param
    sites = [SiteModel(model="gauss_lor", label="A", params={
        "isotropic_chemical_shift_ppm": Param(pos), "shift_fwhm_ppm": Param(5.0),
        "amplitude": Param(amp), "gl": Param(1.0, vary=False)})]
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=sample, sites=sites)
    if with_data:
        x = np.linspace(-30, 30, 200)
        _, y, _ = engine.simulate(rec, exp_ppm=x)
        p = tmp_path / f"{sample}_raw.csv"
        p.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                     "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, y)))
        rec.source_path = str(p)
    path = tmp_path / f"{sample}.recipe.json"
    rec.save(path)
    return path


def test_studio_batch_grid_loads_folder_and_builds_spec(qapp, tmp_path):
    from larmor.desktop.plotting_studio import PlottingStudio
    _fit(tmp_path, "g0"); _fit(tmp_path, "g1")
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(3)
    st._grid_load(str(tmp_path))
    assert len(st._panels) == 2 and st.gridList.count() == 2
    spec = st._spec()
    assert spec["kind"] == "batch_grid"
    assert {Path(p["recipe"]).name.split(".")[0] for p in spec["panels"]} == {"g0", "g1"}
    # unchecking a panel excludes it from the spec without removing it
    item = st.gridList.item(0)
    item.setCheckState(Qt.Unchecked)
    spec2 = st._spec()
    assert len(spec2["panels"]) == 1 and len(st._panels) == 2


def _fit2(tmp_path, sample):
    """A two-site fit, for exercising per-component color/legend/hide UI."""
    from larmor.recipe import Recipe, SiteModel, Param
    sites = [SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(10.0), "shift_fwhm_ppm": Param(5.0),
                "amplitude": Param(80.0), "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="B", params={
                "isotropic_chemical_shift_ppm": Param(-8.0), "shift_fwhm_ppm": Param(3.0),
                "amplitude": Param(40.0), "gl": Param(1.0, vary=False)})]
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=sample, sites=sites)
    path = tmp_path / f"{sample}.recipe.json"
    rec.save(path)
    return path


def test_studio_grid_component_dialog_sets_colors_and_legend_visibility(qapp, tmp_path):
    from larmor.desktop.plotting_studio import PlottingStudio
    _fit2(tmp_path, "g0")
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(3)
    st._grid_load(str(tmp_path))

    sites = st._grid_detect_sites()
    assert [label for _, label in sites] == ["A", "B"]

    # drive the state the dialog would set (offscreen-safe: exercising the
    # real exec() loop needs a live event loop / mocked QColorDialog, which
    # buys nothing over checking the state it commits and _spec() reads)
    st._component_colors[1] = "#ff00ff"
    st._legend_hide.add(1)
    spec = st._spec()
    assert spec["component_colors"] == {1: "#ff00ff"}
    assert spec["legend_hide"] == [1]

    st2 = PlottingStudio(None)
    st2._apply_spec(spec)
    assert st2._component_colors == {1: "#ff00ff"}
    assert st2._legend_hide == {1}


def test_studio_grid_hide_field_roundtrips(qapp, tmp_path):
    from larmor.desktop.plotting_studio import PlottingStudio
    _fit2(tmp_path, "g0")
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(3)
    st._grid_load(str(tmp_path))
    st.gridHide.setText("1")
    spec = st._spec()
    assert spec["hide_components"] == [1]

    st2 = PlottingStudio(None)
    st2._apply_spec(spec)
    assert st2.gridHide.text() == "1"


def test_studio_batch_grid_reorder_and_remove(qapp, tmp_path):
    from larmor.desktop.plotting_studio import PlottingStudio
    _fit(tmp_path, "g0"); _fit(tmp_path, "g1")
    st = PlottingStudio(None)
    st._grid_load(str(tmp_path))
    order = [p["sample"] for p in st._panels]
    st.gridList.setCurrentRow(0)
    st._grid_move(1)
    assert [p["sample"] for p in st._panels] == list(reversed(order))
    st._grid_remove()
    assert len(st._panels) == 1 and st.gridList.count() == 1


def test_studio_batch_grid_asks_for_data_when_csv_cant_auto_match(qapp, tmp_path, monkeypatch):
    """The "successive popups" fallback: an unresolved scope must trigger a
    file-locate dialog, and a path picked there is folded back into the panel
    (has_data flips True, needs_manual clears) rather than being dropped."""
    from larmor.desktop.plotting_studio import PlottingStudio
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text("scope,site,label,param,value,stderr\ng0,s0,A,amplitude,100,2\n")
    found = tmp_path / "g0_manual.csv"
    found.write_text("# nucleus = 11B\n# larmor_MHz = 160\n-10 1\n0 2\n10 1\n")

    monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getOpenFileName",
                        staticmethod(lambda *a, **k: (str(found), "")))
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(3)
    st._grid_load(str(csv_path))
    assert len(st._panels) == 1
    p = st._panels[0]
    assert p["sample"] == "g0" and p["has_data"] and not p["needs_manual"]
    assert p["data_path"] == str(found)
    spec = st._spec()
    assert spec["panels"][0]["data_path"] == str(found)
    assert "recipe" not in spec["panels"][0]


def test_studio_batch_grid_double_click_relocates_unresolved_panel(qapp, tmp_path, monkeypatch):
    """Cancelling the initial popups leaves a panel needing data; a later
    double-click on it re-opens the locate flow rather than being a dead end."""
    from larmor.desktop.plotting_studio import PlottingStudio
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text("scope,site,label,param,value,stderr\ng0,s0,A,amplitude,100,2\n")

    monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getOpenFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getExistingDirectory",
                        staticmethod(lambda *a, **k: ""))
    st = PlottingStudio(None)
    st._grid_load(str(csv_path))
    assert st._panels[0]["needs_manual"]
    assert "⚠" in st.gridList.item(0).text()

    found = tmp_path / "g0_manual.csv"
    found.write_text("# nucleus = 11B\n# larmor_MHz = 160\n-10 1\n0 2\n10 1\n")
    monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getOpenFileName",
                        staticmethod(lambda *a, **k: (str(found), "")))
    st._grid_item_double_clicked(st.gridList.item(0))
    assert not st._panels[0]["needs_manual"] and st._panels[0]["has_data"]
    assert "⚠" not in st.gridList.item(0).text()


def test_studio_grid_arrow_buttons_are_native_qt_arrows(qapp):
    """Regression: text-glyph arrows ("▲"/"▼") silently render blank on some
    Windows font setups -- native QToolButton arrows always draw."""
    from PySide6.QtCore import Qt as _Qt
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    assert st.gridBtnUp.arrowType() == _Qt.UpArrow
    assert st.gridBtnDown.arrowType() == _Qt.DownArrow


def test_studio_grid_drag_reorder_updates_panels(qapp, tmp_path):
    """The list's native InternalMove drag-drop reorders the *visual* items;
    _grid_rows_reordered must fold that back into self._panels by each row's
    stamped original index, not just leave the two out of sync."""
    from larmor.desktop.plotting_studio import PlottingStudio
    _fit(tmp_path, "g0"); _fit(tmp_path, "g1"); _fit(tmp_path, "g2")
    st = PlottingStudio(None)
    st._grid_load(str(tmp_path))
    order = [p["sample"] for p in st._panels]
    # simulate what a drag-drop leaves behind: the visual items reordered,
    # each still carrying its original Qt.UserRole stamp
    item0 = st.gridList.takeItem(0)
    st.gridList.insertItem(st.gridList.count(), item0)
    st._grid_rows_reordered()
    assert [p["sample"] for p in st._panels] == order[1:] + order[:1]
    # re-stamped for the next reorder
    assert [st.gridList.item(i).data(Qt.UserRole) for i in range(3)] == [0, 1, 2]


def test_studio_grid_reconstructs_fit_from_csv_and_marks_it(qapp, tmp_path):
    """The headline ask: a batch CSV with model+source_path columns, and NO
    saved .recipe.json anywhere, still gets a full fit panel (not "data
    only") -- listed as rebuilt, spec references the rebuilt recipe file."""
    from larmor.desktop.plotting_studio import PlottingStudio
    x = np.linspace(-30, 30, 100)
    raw = tmp_path / "g0_raw.csv"
    raw.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                   "\n".join(f"{xi:.4f} {xi:.4f}" for xi in x))
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr,model,source_path\n"
        "shared,s0,A,gl,1,,gauss_lor,\n"
        f"g0,s0,A,isotropic_chemical_shift_ppm,10.0,0.1,gauss_lor,{raw}\n"
        f"g0,s0,A,shift_fwhm_ppm,5.0,0.2,gauss_lor,{raw}\n"
        f"g0,s0,A,amplitude,100.0,3.0,gauss_lor,{raw}\n")
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(3)
    st._grid_load(str(csv_path))
    assert len(st._panels) == 1
    p = st._panels[0]
    assert p["reconstructed"] and p["path"].endswith(".recipe.json")
    assert "(rebuilt from CSV)" in st.gridList.item(0).text()
    spec = st._spec()
    assert spec["panels"][0]["recipe"] == p["path"]
    assert "data_path" not in spec["panels"][0]


def test_studio_grid_data_roots_setting_seeds_locate_dialog(qapp, tmp_path, monkeypatch):
    from PySide6.QtCore import QSettings
    from larmor.desktop.plotting_studio import PlottingStudio
    QSettings("LARMOR", "app").remove("plottingStudio/dataRoots")
    try:
        root = tmp_path / "raw_data"
        root.mkdir()
        st = PlottingStudio(None)
        monkeypatch.setattr("PySide6.QtWidgets.QInputDialog.getText",
                            staticmethod(lambda *a, **k: (str(root), True)))
        st._grid_set_data_roots()
        assert st._grid_data_roots() == [str(root)]

        seen = {}
        def fake_open(*a, **k):
            seen["dir"] = a[2] if len(a) > 2 else k.get("dir", "")
            return ("", "")
        monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getOpenFileName",
                            staticmethod(fake_open))
        monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getExistingDirectory",
                            staticmethod(lambda *a, **k: ""))
        st._grid_ask_manual_path("g0")
        assert seen["dir"] == str(root)
    finally:
        QSettings("LARMOR", "app").remove("plottingStudio/dataRoots")


def test_studio_batch_grid_spec_roundtrips_through_apply_spec(qapp, tmp_path):
    from larmor.desktop.plotting_studio import PlottingStudio
    _fit(tmp_path, "g0")
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(3)
    st._grid_load(str(tmp_path))
    st.gridCols.setValue(2); st.gridComp.setCurrentText("dashed")
    st.gridShade.setText("0,1"); st.gridLabels.setCurrentText("position")
    spec = st._spec()
    assert spec["cols"] == 2 and spec["component_mode"] == "dashed"
    assert spec["shade_only"] == [0, 1] and spec["peak_labels"] == "position"

    st2 = PlottingStudio(None)
    st2._apply_spec(spec)
    assert st2.kind.currentIndex() == 3
    assert st2.gridCols.value() == 2 and st2.gridComp.currentText() == "dashed"
    assert st2.gridShade.text() == "0,1"
    assert len(st2._panels) == 1 and st2._panels[0]["sample"] == "g0"


def test_studio_species_bar_table_builds_spec(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(4)
    st.barTable.setRowCount(2); st.barTable.setColumnCount(3)
    st.barTable.setHorizontalHeaderLabels(["category", "Q2", "Q3"])
    for r, (cat, q2, q3) in enumerate([("P-5", 30, 70), ("P-10", 60, 40)]):
        st.barTable.setItem(r, 0, QTableWidgetItem(cat))
        st.barTable.setItem(r, 1, QTableWidgetItem(str(q2)))
        st.barTable.setItem(r, 2, QTableWidgetItem(str(q3)))
    spec = st._spec()
    assert spec["kind"] == "species_bar"
    assert spec["categories"] == ["P-5", "P-10"]
    assert {s["label"] for s in spec["series"]} == {"Q2", "Q3"}
    q2 = next(s for s in spec["series"] if s["label"] == "Q2")
    assert q2["values"] == pytest.approx([30.0, 60.0])
    import matplotlib.pyplot as plt
    from larmor import figures
    fig = figures.render(spec)
    plt.close(fig)


def test_studio_species_bar_roundtrips_through_apply_spec(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    spec = {"kind": "species_bar", "categories": ["a", "b"],
            "series": [{"label": "Q2", "values": [1, 2]},
                      {"label": "Q3", "values": [3, 4]}]}
    st = PlottingStudio(None)
    st._apply_spec(spec)
    assert st.kind.currentIndex() == 4
    assert st.barTable.item(0, 0).text() == "a"
    assert st.barTable.item(1, 2).text() == "4"
    spec2 = st._spec()
    assert spec2["categories"] == ["a", "b"]


def test_studio_species_bar_load_from_csv_pivots_a_parameter(qapp, monkeypatch, tmp_path):
    from larmor.desktop.plotting_studio import PlottingStudio
    csv_path = tmp_path / "batch_table.csv"
    csv_path.write_text(
        "scope,site,label,param,value,stderr\n"
        "P-5,s0,Q2,amplitude,30,1\nP-5,s1,Q3,amplitude,70,2\n"
        "P-10,s0,Q2,amplitude,60,1\nP-10,s1,Q3,amplitude,40,2\n"
        "shared,s0,Q2,shift_fwhm_ppm,4.0,\n")
    monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getOpenFileName",
                        staticmethod(lambda *a, **k: (str(csv_path), "")))
    monkeypatch.setattr("PySide6.QtWidgets.QInputDialog.getItem",
                        staticmethod(lambda *a, **k: ("amplitude", True)))
    st = PlottingStudio(None)
    st.kind.setCurrentIndex(4)
    st._bar_load_csv()
    spec = st._spec()
    assert spec["categories"] == ["P-5", "P-10"]
    by_label = {s["label"]: s["values"] for s in spec["series"]}
    assert by_label["Q2"] == pytest.approx([30.0, 60.0])
    assert by_label["Q3"] == pytest.approx([70.0, 40.0])


def test_studio_template_applies_kind_and_spec_defaults(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    from larmor import figures
    st = PlottingStudio(None)
    name = next(n for n, t in figures.TEMPLATES.items() if t["kind"] == "batch_grid")
    idx = st.template.findData(name)
    assert idx >= 0
    st.template.setCurrentIndex(idx)
    assert st.kind.currentIndex() == 3
    assert st.templateDesc.text() == figures.TEMPLATES[name]["description"]


def test_studio_trace_defaults_from_template_apply_to_new_traces(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    st._trace_defaults = {"end_label": True}
    st._push_trace({"data": {"x": [1, 2], "y": [1, 2]}, "label": "s0"})
    assert st._traces[-1]["end_label"] is True
    assert st._traces[-1]["label"] == "s0"     # trace-specific fields still win


def test_studio_auto_update_defaults_off_and_schedule_is_a_noop(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    assert not st.chkAuto.isChecked()
    st._debounce.stop()
    st._schedule()
    assert not st._debounce.isActive()   # off -> _schedule() does nothing
    st.chkAuto.setChecked(True)
    st._schedule()
    assert st._debounce.isActive()       # on -> the debounce timer arms


def test_studio_preview_button_renders_regardless_of_auto_update(qapp):
    from larmor.desktop.plotting_studio import PlottingStudio
    st = PlottingStudio(None)
    assert not st.chkAuto.isChecked()
    st._push_trace(_trace())
    before = st._canvas.figure
    st.title.setText("new title")     # a control change; auto update is off
    assert st._canvas.figure is before     # nothing rendered yet
    st._refresh()                          # what the Preview button calls
    assert st._canvas.figure.axes[0].get_title() == "new title"


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
