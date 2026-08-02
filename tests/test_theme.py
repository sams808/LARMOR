"""Guards for the colour themes.

The whole point of the presets is that everything stays visible and readable, so
these tests compute WCAG contrast for every text-on-background and
mark-on-plot-background pair and fail if any theme drops below the floor. If a
future colour tweak breaks readability, this catches it.
"""
import pytest

from larmor.desktop import theme as T


def test_ten_presets_each_with_a_full_series():
    assert len(T.THEMES) == 10
    for name, t in T.THEMES.items():
        assert len(t.series) >= 10, name
        assert t.name == name


@pytest.mark.parametrize("name", list(T.THEMES))
def test_contrast_floors(name):
    t = T.THEMES[name]
    C = T.contrast
    # text must be clearly readable on every surface it sits on (WCAG AA 4.5)
    for role, bg in (("window", t.window), ("base", t.base), ("alt_base", t.alt_base),
                     ("hover", t.hover), ("surface", t.surface)):
        assert C(t.text, bg) >= 4.5, (name, "text/" + role, round(C(t.text, bg), 2))
    assert C(t.text, t.header) >= 4.0, (name, "text/header")
    # secondary ink (headers, axes) at least the 3:1 large-text/graphics floor
    assert C(t.text_dim, t.window) >= 3.0, (name, "text_dim/window")
    # accent text on the accent (default buttons, selections)
    assert C(t.accent_text, t.accent) >= 4.5, (name, "accent_text/accent")
    assert C(t.accent, t.base) >= 3.0, (name, "accent/base")
    # plot marks vs the plot background (3:1 graphical floor)
    assert C(t.experiment, t.plot_bg) >= 4.5, (name, "experiment/plot")
    for role in ("model", "baseline", "pivot", "measure"):
        v = C(getattr(t, role), t.plot_bg)
        assert v >= 3.0, (name, role + "/plot", round(v, 2))
    # every categorical series colour distinguishable from the plot background
    for i, s in enumerate(t.series):
        v = C(s, t.plot_bg)
        assert v >= 3.0, (name, f"series[{i}]={s}", round(v, 2))


def test_qss_is_complete_and_parametric():
    t = T.get("Dark")
    css = T.qss(t)
    assert t.window in css and t.text in css and t.accent in css
    for widget in ("QMenu", "QPushButton", "QTableWidget", "QHeaderView",
                   "QToolTip", "QTabBar", "QLineEdit"):
        assert widget in css
    # no leftover hardcoded light colour from the old stylesheet
    assert "#f0f2f0" not in css or t.window == "#f0f2f0"


def test_set_active_and_get_roundtrip():
    T.set_active("Nord")
    assert T.active().name == "Nord"
    T.set_active(T.DEFAULT)
    assert T.active().name == T.DEFAULT
    assert T.get("does-not-exist").name == T.DEFAULT   # falls back, never raises


def test_no_hardcoded_light_surfaces_in_desktop():
    """Guard the theme rollout: no hardcoded LIGHT background stays in the
    desktop layer, or it would be an unreadable white block on the dark themes.
    (The help viewer renders manuals as a deliberate light 'document' — allowed.)"""
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[1] / "larmor" / "desktop"
    allow = {"theme.py", "help_dialog.py", "mdrender.py"}
    bad_literals = [
        # light surfaces (an unreadable white block on dark themes)
        "#fcfdfc", "#f3f5f3", "#f6f8f6", "#eef1ee", "#e2f0f0",
        "background: #ff", "background:#ff", 'background="w"', "background='w'",
        # old fixed 'dim/secondary/accent' ink that goes low-contrast on dark —
        # these must come from theme roles (text_dim / accent), not literals
        "#93a0a8", "#5a6871", "#4a5560", "#a0aec0", "#0a5a62",
    ]
    offenders = []
    for f in sorted(root.glob("*.py")):
        if f.name in allow:
            continue
        text = f.read_text(encoding="utf-8")
        for lit in bad_literals:
            if lit in text:
                offenders.append(f"{f.name}: {lit}")
    assert not offenders, "hardcoded light surfaces (break dark themes): " + \
        "; ".join(offenders)


def test_palette_builds_for_every_theme():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])   # noqa: F841
    for name in T.names():
        pal = T.palette(T.get(name))
        assert pal is not None
