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
