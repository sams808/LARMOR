from pathlib import Path

import numpy as np
import pytest

from larmor import figures

from conftest import BRUKER_2RR_MQMAS, CAALGLASS, EXPNO_1901, NMRVEW_2D, require


def test_styles_complete():
    for name, s in figures.STYLES.items():
        assert "figsize" in s and "rc" in s, name


def test_nucleus_xlabel():
    assert figures.nucleus_xlabel("27Al") == r"$^{27}$Al NMR shift (ppm)"
    assert figures.nucleus_xlabel("19F") == r"$^{19}$F NMR shift (ppm)"


def test_render_1d_inline():
    x = np.linspace(-50, 50, 500)
    spec = {
        "kind": "1d",
        "traces": [
            {"data": {"x": x.tolist(), "y": np.exp(-(x / 8.0) ** 2).tolist()},
             "label": "a", "normalize": True},
            {"data": {"x": x.tolist(), "y": (2 * np.exp(-(x / 4.0) ** 2)).tolist()},
             "label": "b", "offset": 1.2, "scale": 0.5, "normalize": True},
        ],
        "xlim": [40, -40],
        "annotations": [{"x": 30, "y": 1.5, "text": "test"}],
    }
    fig = figures.render(spec)
    ax = fig.axes[0]
    assert len(ax.lines) == 2
    assert ax.get_xlim()[0] > ax.get_xlim()[1]  # ppm axis reversed
    # normalize+scale+offset arithmetic: trace b peaks at 0.5*1 + 1.2
    assert ax.lines[1].get_ydata().max() == pytest.approx(1.7, abs=1e-6)


def test_render_1d_from_fxmla_and_recipe():
    require(CAALGLASS)
    recipe_path = Path(__file__).resolve().parents[1] / "examples" / "CaAlGlass.recipe.json"
    if not recipe_path.exists():
        pytest.skip("example recipe not present")
    spec = {
        "kind": "1d",
        "traces": [
            {"path": str(CAALGLASS), "label": "exp"},
            {"recipe": str(recipe_path), "part": "total", "label": "fit"},
            {"recipe": str(recipe_path), "part": "site", "site": 0},
            {"recipe": str(recipe_path), "part": "residual"},
        ],
        "xlim": [150, -80],
    }
    png = figures.render_png_bytes(spec)
    assert len(png) > 10_000


def test_render_2d_nmrvew():
    require(NMRVEW_2D)
    spec = {"kind": "2d", "path": str(NMRVEW_2D), "style": "article",
            "levels": {"n": 8}, "slopes": [{"slope": 1.0, "intercept": 0.0}]}
    png = figures.render_png_bytes(spec)
    assert len(png) > 10_000


def test_render_2d_nucleus_labels_when_given():
    import matplotlib.pyplot as plt
    require(NMRVEW_2D)
    fig = figures.render({"kind": "2d", "path": str(NMRVEW_2D), "nucleus": "17O"})
    ax = fig.axes[0]
    assert "17" in ax.get_xlabel() and "O" in ax.get_xlabel()
    assert r"\delta_1" in ax.get_ylabel()
    plt.close(fig)
    # unset -> the original generic defaults (no regression for existing specs)
    fig2 = figures.render({"kind": "2d", "path": str(NMRVEW_2D)})
    assert fig2.axes[0].get_xlabel() == "F2 shift (ppm)"
    plt.close(fig2)


@pytest.mark.slow
def test_render_2d_overlays_a_real_mqmas_fit(tmp_path):
    """B1: a fitted 2D recipe overlays as a dashed contour on the real
    experimental map it was fit against -- the missing piece render_2d had
    for a genuinely publication-ready MQMAS figure (previously: raw contour
    only, no way to show the fit)."""
    require(BRUKER_2RR_MQMAS)
    from larmor import twod
    from larmor.recipe import Recipe, SiteModel, Param

    data = twod.read_bruker_2d(str(BRUKER_2RR_MQMAS)).region(
        f2_range=(-50, 150), f1_range=(-20, 150)).normalized()
    k = twod.build_mqmas_kernel(
        data.nucleus, data.larmor_MHz, f2_window=(150.0, -50.0),
        f1_window=(150.0, -20.0), n2=96, n1=64, n_cq=12, n_eta=4,
        cq_max_MHz=14.0)
    rec = Recipe(nucleus=data.nucleus, larmor_frequency_MHz=data.larmor_MHz,
                spin_rate_Hz=data.spin_rate_Hz, sites=[
        SiteModel(model="czjzek", label="Al", params={
            "isotropic_chemical_shift_ppm": Param(65.0, min=30, max=100),
            "sigma_Cq_MHz": Param(3.0, min=0.2, max=8.0),
            "shift_fwhm_ppm": Param(6.0, min=1.0, max=25.0),
            "amplitude": Param(1.0, min=0.0)})])
    twod.fit_2d(rec, data, kernel=k)
    recipe_path = tmp_path / "mqmas_fit.recipe.json"
    rec.save(recipe_path)

    expno = BRUKER_2RR_MQMAS.parents[2]     # .../35/pdata/1/2rr -> .../35
    import matplotlib.pyplot as plt
    fig_plain = figures.render({"kind": "2d", "path": str(expno)})
    fig_fit = figures.render({"kind": "2d", "path": str(expno),
                              "fit_recipe": str(recipe_path)})
    ax_plain, ax_fit = fig_plain.axes[0], fig_fit.axes[0]
    # the overlay adds contour line collections beyond the plain experimental figure
    assert len(ax_fit.collections) > len(ax_plain.collections)
    assert "27" in ax_fit.get_xlabel() and "Al" in ax_fit.get_xlabel()
    plt.close(fig_plain); plt.close(fig_fit)


def test_render_series_satrec_and_redor():
    require(EXPNO_1901)
    satrec = figures.load_series({"path": str(EXPNO_1901), "mode": "satrec"})
    assert satrec["x"].size >= 5
    assert satrec["x"][0] == 0.0 and satrec["x"][-1] > 10  # delays in seconds
    assert np.all(satrec["y"] >= 0) and satrec["y"].max() == pytest.approx(1.0)

    png = figures.render_png_bytes({
        "kind": "series", "mode": "satrec", "path": str(EXPNO_1901),
        "stretched": True})
    assert len(png) > 10_000

    redor = figures.load_series({"path": str(EXPNO_1901), "mode": "redor"})
    assert redor["x"].size >= 3
    png = figures.render_png_bytes({
        "kind": "series", "mode": "redor", "path": str(EXPNO_1901)})
    assert len(png) > 5_000


def test_unknown_kind_rejected():
    with pytest.raises(ValueError, match="unknown figure kind"):
        figures.render({"kind": "3d-hologram"})


def _saved_fit(tmp_path, sample, pos, amp, seed, with_data=True):
    from larmor import engine
    from larmor.recipe import Recipe, SiteModel, Param
    sites = [SiteModel(model="gauss_lor", label=chr(65 + i), params={
        "isotropic_chemical_shift_ppm": Param(p), "shift_fwhm_ppm": Param(w),
        "amplitude": Param(a), "gl": Param(1.0, vary=False)})
            for i, (p, w, a) in enumerate(zip(pos, [5.0] * len(pos), amp))]
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample=sample, sites=sites)
    if with_data:
        x = np.linspace(-40, 40, 400)
        _, y, _ = engine.simulate(rec, exp_ppm=x)
        data = y + np.random.default_rng(seed).normal(0, 1.0, x.size)
        csv_path = tmp_path / f"{sample}_raw.csv"
        csv_path.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                            "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, data)))
        rec.source_path = str(csv_path)
    p = tmp_path / f"{sample}.recipe.json"
    rec.save(p)
    return str(p)


def test_load_trace_dataless_recipe_falls_back_to_a_framed_model(tmp_path):
    """A non-kernel recipe with NO source_path (or an unreachable one) must
    fall back to a sites-framed model grid, not raise -- the same graceful
    behaviour a data-less dmfit .fxmla already gets."""
    p = _saved_fit(tmp_path, "g0", [10.0], [50.0], 0, with_data=False)
    x, y, meta = figures.load_trace({"recipe": p, "part": "total"})
    assert len(x) > 100 and np.isfinite(y).all()
    assert abs(float(x[np.argmax(y)]) - 10.0) < 3.0
    with pytest.raises(ValueError, match="source_path"):
        figures.load_trace({"recipe": p, "part": "residual"})


def test_site_colors_are_the_fixed_categorical_palette():
    assert len(figures.SITE_COLORS) >= 8
    assert figures.site_color(0) == figures.SITE_COLORS[0]
    assert figures.site_color(len(figures.SITE_COLORS)) == figures.SITE_COLORS[0]  # wraps


def test_render_batch_grid_deconvolution(tmp_path):
    paths = [_saved_fit(tmp_path, f"g{k}", [10.0, -5.0], [80.0, 40.0], k)
            for k in range(3)]
    spec = {"kind": "batch_grid", "panels": [{"recipe": p} for p in paths],
            "component_mode": "fill", "peak_labels": "position", "cols": 2}
    fig = figures.render(spec)
    assert len(fig.axes) == 4                       # 3 panels + 1 hidden filler
    visible = [a for a in fig.axes if a.get_visible()]
    assert len(visible) == 3
    assert {a.get_title() for a in visible} == {"g0", "g1", "g2"}
    # a shared legend names the two components once, not per panel
    assert fig.legends and {t.get_text() for t in fig.legends[0].get_texts()} == {"A", "B"}
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_render_batch_grid_shade_only_hides_the_rest():
    """Category-3 style: only the flagged component(s) draw; everything else
    is invisible except the total fit line."""
    import tempfile
    from pathlib import Path as _P
    tmp_path = _P(tempfile.mkdtemp())
    paths = [_saved_fit(tmp_path, f"g{k}", [10.0, -5.0], [80.0, 40.0], k)
            for k in range(2)]
    spec = {"kind": "batch_grid", "panels": [{"recipe": p} for p in paths],
            "component_mode": "fill", "shade_only": [1],
            "peak_labels": "position+pct"}
    fig = figures.render(spec)
    # only site B (index 1) gets a legend entry
    assert {t.get_text() for t in fig.legends[0].get_texts()} == {"B"}
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_render_batch_grid_hidden_components_shows_no_legend():
    import tempfile
    from pathlib import Path as _P
    tmp_path = _P(tempfile.mkdtemp())
    p = _saved_fit(tmp_path, "g0", [10.0, -5.0], [80.0, 40.0], 0)
    fig = figures.render({"kind": "batch_grid",
                          "panels": [{"recipe": p}],
                          "component_mode": "hidden"})
    assert not fig.legends
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_render_batch_grid_component_colors_override_the_palette(tmp_path):
    import matplotlib.pyplot as plt
    p = _saved_fit(tmp_path, "g0", [10.0, -5.0], [80.0, 40.0], 0)
    fig = figures.render({"kind": "batch_grid", "panels": [{"recipe": p}],
                          "component_colors": {1: "#ff00ff"}})
    ax = fig.axes[0]
    colors = {tuple(ln.get_color() if isinstance(ln.get_color(), tuple)
                    else plt.matplotlib.colors.to_rgb(ln.get_color()))
             for ln in ax.lines if ln.get_label() not in ("_experiment", "_fit")}
    assert plt.matplotlib.colors.to_rgb("#ff00ff") in colors
    assert figures.site_color(0) not in ("#ff00ff",)   # site 0 kept its default
    plt.close(fig)


def test_render_batch_grid_hide_components_drops_line_and_legend(tmp_path):
    import matplotlib.pyplot as plt
    p = _saved_fit(tmp_path, "g0", [10.0, -5.0], [80.0, 40.0], 0)
    fig = figures.render({"kind": "batch_grid", "panels": [{"recipe": p}],
                          "hide_components": [1]})
    assert {t.get_text() for t in fig.legends[0].get_texts()} == {"A"}
    plt.close(fig)


def test_render_batch_grid_legend_hide_keeps_line_drops_label_only(tmp_path):
    import matplotlib.pyplot as plt
    p = _saved_fit(tmp_path, "g0", [10.0, -5.0], [80.0, 40.0], 0)
    fig = figures.render({"kind": "batch_grid", "panels": [{"recipe": p}],
                          "legend_hide": [1]})
    ax = fig.axes[0]
    # both components still drew (2 component lines + total + experiment)
    assert len(ax.lines) == 4
    assert {t.get_text() for t in fig.legends[0].get_texts()} == {"A"}
    plt.close(fig)


def test_render_batch_grid_explicit_legend_loc_is_honoured(tmp_path):
    import matplotlib.pyplot as plt
    from matplotlib.legend import Legend
    p = _saved_fit(tmp_path, "g0", [10.0], [80.0], 0)
    fig = figures.render({"kind": "batch_grid", "panels": [{"recipe": p}],
                          "legend_loc": "upper right"})
    assert fig.legends[0]._loc == Legend.codes["upper right"]
    plt.close(fig)


def test_render_batch_grid_excluded_site_never_drawn_even_in_a_real_recipe(tmp_path):
    """A saved recipe can carry a site explicitly zero-locked (batch fit's
    "Exclude component") -- it must never draw a line or a legend entry,
    even though the recipe file itself still lists it (full fidelity)."""
    import matplotlib.pyplot as plt
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor import engine

    x = np.linspace(-30, 30, 300)
    sites = [SiteModel(model="gauss_lor", label="A", params={
                "isotropic_chemical_shift_ppm": Param(10.0), "shift_fwhm_ppm": Param(5.0),
                "amplitude": Param(80.0), "gl": Param(1.0, vary=False)}),
            SiteModel(model="gauss_lor", label="B", params={
                "isotropic_chemical_shift_ppm": Param(-8.0), "shift_fwhm_ppm": Param(3.0),
                "amplitude": Param(0.0, vary=False, min=0.0, max=0.0),
                "gl": Param(1.0, vary=False)})]
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sample="g0", sites=sites)
    _, y, _ = engine.simulate(rec, exp_ppm=x)
    raw = tmp_path / "g0_raw.csv"
    raw.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                   "\n".join(f"{xi:.4f} {yi:.4f}" for xi, yi in zip(x, y)))
    rec.source_path = str(raw)
    p = tmp_path / "g0.recipe.json"
    rec.save(p)

    fig = figures.render({"kind": "batch_grid", "panels": [{"recipe": str(p)}]})
    assert fig.legends and {t.get_text() for t in fig.legends[0].get_texts()} == {"A"}
    plt.close(fig)


def test_render_batch_grid_requires_at_least_one_panel():
    with pytest.raises(ValueError, match="at least one panel"):
        figures.render({"kind": "batch_grid", "panels": []})


def test_render_batch_grid_accepts_data_only_and_empty_panels(tmp_path):
    """A batch CSV can name a sample series_grid couldn't pair with a saved
    fit -- panels without "recipe" must still render (experiment-only, or a
    placeholder), never break the rest of the grid."""
    import matplotlib.pyplot as plt

    fitted = _saved_fit(tmp_path, "g0", [10.0], [50.0], 0)
    data_only = tmp_path / "g1_raw.csv"
    x = np.linspace(-30, 30, 100)
    data_only.write_text("# nucleus = 11B\n# larmor_MHz = 160\n" +
                         "\n".join(f"{xi:.4f} {np.exp(-(xi/6)**2):.4f}" for xi in x))

    fig = figures.render({"kind": "batch_grid", "panels": [
        {"recipe": fitted},
        {"data_path": str(data_only), "title": "g1"},
        {"title": "g2 (unresolved)"},
    ], "cols": 3})
    titles = [a.get_title() for a in fig.axes if a.get_visible()]
    assert titles == ["g0", "g1", "g2 (unresolved)"]
    # the data-only panel drew a line; the fully-empty one drew none
    axes = [a for a in fig.axes if a.get_visible()]
    assert len(axes[1].lines) == 1
    assert len(axes[2].lines) == 0
    assert any(t.get_text() == "no data located" for t in axes[2].texts)
    plt.close(fig)


def test_render_species_bar_normalizes_and_stacks():
    spec = {"kind": "species_bar", "categories": ["P-5", "P-10"],
            "series": [{"label": "Q0", "values": [10, 20]},
                      {"label": "Q1", "values": [30, 20]},
                      {"label": "Q2", "values": [60, 60]}]}
    fig = figures.render(spec)
    ax = fig.axes[0]
    assert ax.get_ylim() == pytest.approx((0.0, 100.0))
    assert [t.get_text() for t in ax.get_xticklabels()] == ["P-5", "P-10"]
    bars = [c for c in ax.containers]
    assert len(bars) == 3                            # one per series
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_render_species_bar_needs_categories_and_series():
    with pytest.raises(ValueError):
        figures.render({"kind": "species_bar", "categories": [], "series": []})


def test_end_label_places_text_at_the_traces_own_displayed_edge():
    import matplotlib.pyplot as plt
    spec = {"kind": "1d", "x_is_ppm": True, "hide_yaxis": True,
            "legend_loc": "none",
            "traces": [{"data": {"x": [1, 2, 3], "y": [1, 2, 1]},
                       "label": "sample-A", "end_label": True}]}
    fig = figures.render(spec)
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert texts == ["sample-A"]
    assert not fig.axes[0].get_legend()              # legend_loc none honoured
    plt.close(fig)


def test_templates_are_well_formed_and_generic():
    assert len(figures.TEMPLATES) >= 6
    for name, tpl in figures.TEMPLATES.items():
        assert tpl["kind"] in figures.RENDERERS
        assert "description" in tpl and "spec" in tpl
        # generic by construction: no nucleus name leaks into a template name
        for nuc in ("1H", "13C", "27Al", "29Si", "31P", "11B", "23Na"):
            assert nuc not in name
