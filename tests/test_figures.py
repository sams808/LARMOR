from pathlib import Path

import numpy as np
import pytest

from larmor import figures

from conftest import CAALGLASS, EXPNO_1901, NMRVEW_2D, require


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
