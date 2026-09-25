"""The themes are what they say, and can be told apart.

Sam: "make sure the themes are actually what they say: the ones in 'More
themes' are not discernible enough and don't match their names". So, beyond
the contrast floor of tests/test_theme.py, every theme is (a) pinned to the
canonical palette of the scheme it is named after (Solarized, Nord) or to the
identity it was given (Ocean, Slate, Sepia, Y2K, Vaporwave, ...), (b) held a
perceptual distance from every other theme in CIE Lab, and (c) followed by the
plot -- background, axes, grid, pens and site colours all read from the theme.
"""
import itertools
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from larmor.desktop import theme as T  # noqa: E402

ALL = {**T.THEMES, **T.AESTHETIC_THEMES}

# Solarized -- Ethan Schoonover's sixteen
SOL = dict(base03="#002b36", base02="#073642", base01="#586e75", base00="#657b83",
           base0="#839496", base1="#93a1a1", base2="#eee8d5", base3="#fdf6e3",
           yellow="#b58900", orange="#cb4b16", red="#dc322f", magenta="#d33682",
           violet="#6c71c4", blue="#268bd2", cyan="#2aa198", green="#859900")
# Nord -- polar night, snow storm, frost, aurora
NORD = ["#2e3440", "#3b4252", "#434c5e", "#4c566a", "#d8dee9", "#e5e9f0", "#eceff4",
        "#8fbcbb", "#88c0d0", "#81a1c1", "#5e81ac", "#bf616a", "#d08770", "#ebcb8b",
        "#a3be8c", "#b48ead"]


# ------------------------------------------------------------ Lab helpers
def test_lab_and_delta_e_reference_values():
    L, a, b = T.lab("#ffffff")
    assert L == pytest.approx(100.0, abs=0.01) and abs(a) < 0.01 and abs(b) < 0.01
    L, a, b = T.lab("#000000")
    assert L == pytest.approx(0.0, abs=1e-6)
    # sRGB red and the mid grey, against the textbook conversions
    L, a, b = T.lab("#ff0000")
    assert (L, a, b) == pytest.approx((53.24, 80.09, 67.20), abs=0.1)
    L, a, b = T.lab("#808080")
    assert L == pytest.approx(53.59, abs=0.05) and abs(a) < 0.02 and abs(b) < 0.02
    assert T.delta_e("#123456", "#123456") == 0.0
    assert T.delta_e("#000000", "#ffffff") == pytest.approx(100.0, abs=0.01)
    # hue: red ~40°, yellow ~103°, green ~136°, blue ~306° in Lab; grey has no chroma
    h, c = T.hue_chroma("#ff0000")
    assert h == pytest.approx(40.0, abs=1.0) and c > 100
    h, c = T.hue_chroma("#0000ff")
    assert h == pytest.approx(306.3, abs=1.0)
    _, c = T.hue_chroma("#777777")
    assert c < 0.1
    assert T.hue_difference("#ff0000", "#0000ff") == pytest.approx(360 - 306.3 + 40.0, abs=1.5)
    assert T.hue_difference("#2dd4bf", "#2dd4bf") == 0.0


# ------------------------------------------------------------ (a) identity
def test_solarized_themes_use_schoonovers_palette():
    lt, dk = T.THEMES["Solarized Light"], T.THEMES["Solarized Dark"]
    assert (lt.window, lt.base) == (SOL["base2"], SOL["base3"])
    assert (dk.window, dk.base) == (SOL["base02"], SOL["base03"])
    assert lt.accent == dk.accent == SOL["blue"]
    assert lt.plot_bg == SOL["base3"] and dk.plot_bg == SOL["base03"]
    # inks: the darkest / lightest canonical tones (base01 and base1, the
    # scheme's emphasised text, miss the 4.5 floor -- see the comments)
    assert lt.text == SOL["base02"] and lt.text_dim == SOL["base01"]
    assert dk.text == SOL["base2"] and dk.text_dim == SOL["base1"]
    assert lt.border == SOL["base1"] and dk.border == SOL["base01"]
    # the dark site palette IS the eight accents plus base1 / base0
    assert set(dk.series) == {SOL[k] for k in ("blue", "orange", "green", "magenta",
                                               "yellow", "cyan", "violet", "red",
                                               "base1", "base0")}
    # the light one keeps every accent that clears 3:1 on base3 verbatim and
    # moves the three that do not (green, yellow, cyan) by ΔE < 2
    canonical = {SOL[k] for k in SOL}
    moved = [s for s in lt.series if s not in canonical]
    assert len(moved) == 3
    for s, ref in zip(moved, (SOL["green"], SOL["yellow"], SOL["cyan"])):
        assert T.delta_e(s, ref) < 2.0, (s, ref)
        assert T.contrast(s, lt.plot_bg) >= 3.0
    assert "Solarized" in lt.scheme and "Solarized" in dk.scheme


def test_nord_theme_uses_the_nord_palette():
    t = T.THEMES["Nord"]
    for role in (t.window, t.base, t.text, t.text_dim, t.border, t.accent, t.plot_bg):
        assert role in NORD, role
    assert t.window == "#3b4252" and t.base == t.plot_bg == "#2e3440"   # nord1 / nord0
    assert t.accent == "#88c0d0"                                        # nord8
    assert t.text == "#eceff4" and t.text_dim == "#81a1c1"              # nord6 / nord9
    assert set(t.series) <= set(NORD) and len(t.series) == 10
    assert "Nord" in t.scheme


def test_named_identities_read_as_their_names():
    """Themes without a published scheme: check the identity they were given
    in Lab terms (hue of the window and accent), not just a hex string."""
    def hue(h):
        return T.hue_chroma(h)[0]

    def chroma(h):
        return T.hue_chroma(h)[1]

    ocean = T.THEMES["Ocean"]
    assert ocean.is_dark and 200 < hue(ocean.window) < 290 and chroma(ocean.window) > 12   # blue-teal water
    assert T.lab(ocean.accent)[0] > 88 and 150 < hue(ocean.accent) < 220              # sea foam: near-white aqua
    slate = T.THEMES["Slate"]
    assert slate.is_dark and 230 < hue(slate.window) < 300 and 5 < chroma(slate.window) < 15   # cool blue-grey
    assert 60 < hue(slate.accent) < 100 and chroma(slate.accent) > 40                    # brass
    sepia = T.THEMES["Sepia (paper)"]
    assert not sepia.is_dark and 70 < hue(sepia.window) < 100 and chroma(sepia.window) > 12  # aged paper (warm)
    assert 40 < hue(sepia.accent) < 80 and T.lab(sepia.accent)[0] < 45                    # brown ink
    assert 40 < hue(sepia.text) < 90                                                      # brown-black ink
    dark = T.THEMES["Dark"]
    assert dark.is_dark and chroma(dark.window) < 5                                      # neutral graphite
    assert 165 < hue(dark.accent) < 200 and chroma(dark.accent) > 30                      # brand teal
    light = T.THEMES["Light"]
    assert not light.is_dark and chroma(light.window) < 3 and T.lab(light.window)[0] > 93
    assert 180 < hue(light.accent) < 215                                                  # brand teal (light)
    y2k = T.AESTHETIC_THEMES["Y2K"]
    assert not y2k.is_dark and chroma(y2k.window) < 8 and 78 < T.lab(y2k.window)[0] < 90  # brushed silver
    assert 125 < hue(y2k.accent) < 150 and chroma(y2k.accent) > 50                        # lime
    assert "#9be32a" in y2k.flourish_qss                                                  # the gloss carries the bright lime
    vapor = T.AESTHETIC_THEMES["Vaporwave"]
    assert vapor.is_dark and 290 < hue(vapor.window) < 340                                # indigo
    assert 320 < hue(vapor.accent) or hue(vapor.accent) < 10                              # hot pink
    assert chroma(vapor.accent) > 70
    dream = T.AESTHETIC_THEMES["Dreamcore"]
    assert not dream.is_dark and 290 < hue(dream.window) < 340                            # lilac
    assert 290 < hue(dream.accent) < 320 and chroma(dream.accent) > 40                    # violet
    soft = T.AESTHETIC_THEMES["Gen X Soft Club"]
    assert not soft.is_dark and 140 < hue(soft.window) < 185 and chroma(soft.window) > 5   # mint
    assert (hue(soft.accent) < 30 or hue(soft.accent) > 340) and chroma(soft.accent) > 30  # rose
    for t in ALL.values():
        assert t.scheme, t.name                          # every theme says what it is


def test_high_contrast_themes_are_pure_black_and_white():
    hl, hd = T.THEMES["High Contrast Light"], T.THEMES["High Contrast Dark"]
    assert hl.window == hl.base == hl.plot_bg == "#ffffff" and hl.text == "#000000"
    assert hd.window == hd.plot_bg == "#000000" and hd.text == "#ffffff"
    assert T.contrast(hl.text, hl.window) == T.contrast(hd.text, hd.window) == 21.0


# ------------------------------------------------------------ (b) discernible
WINDOW_MIN, ACCENT_MIN = 6.0, 25.0
SAME_HUE_DEG, SAME_HUE_WINDOW_MIN = 20.0, 8.0
TUPLE_MIN = 40.0


@pytest.mark.parametrize("a,b", list(itertools.combinations(ALL, 2)))
def test_every_pair_of_themes_is_discernible(a, b):
    ta, tb = ALL[a], ALL[b]
    d_win = T.delta_e(ta.window, tb.window)
    d_acc = T.delta_e(ta.accent, tb.accent)
    # the two things a glance registers: the big surface or the highlight
    assert d_win >= WINDOW_MIN or d_acc >= ACCENT_MIN, (
        f"{a} vs {b}: window ΔE {d_win:.1f} and accent ΔE {d_acc:.1f} -- "
        "two themes that look the same")
    if ta.is_dark != tb.is_dark:
        return                          # a light and a dark theme never blur
    # same polarity: the (window, base, accent, text) tuple must be well apart...
    d_tuple = (d_win + T.delta_e(ta.base, tb.base) + d_acc
               + T.delta_e(ta.text, tb.text))
    assert d_tuple >= TUPLE_MIN, (a, b, round(d_tuple, 1))
    # ...and two same-hue accents (both teal, both blue...) may only coexist
    # on windows that are clearly different shades
    ha, ca = T.hue_chroma(ta.accent)
    hb, cb = T.hue_chroma(tb.accent)
    same_hue = T.hue_difference(ta.accent, tb.accent) < SAME_HUE_DEG and min(ca, cb) > 15
    if same_hue:
        assert d_win >= SAME_HUE_WINDOW_MIN, (
            f"{a} vs {b}: both accents are the same hue "
            f"({ha:.0f}° / {hb:.0f}°) on windows only ΔE {d_win:.1f} apart")


def test_dark_theme_accents_are_not_all_teal():
    """The complaint in one number: among the dark presets, at most two
    accents may share a hue family, and Nord / Ocean / Solarized Dark --
    the three that used to be teal-on-dark-blue -- must all differ."""
    darks = [t for t in T.THEMES.values() if t.is_dark]
    families = {}
    for t in darks:
        h, c = T.hue_chroma(t.accent)
        families.setdefault(int(h // 30), []).append(t.name)
    assert max(len(v) for v in families.values()) <= 2, families
    trio = ("Nord", "Ocean", "Solarized Dark")
    for a, b in itertools.combinations(trio, 2):
        assert T.delta_e(T.THEMES[a].accent, T.THEMES[b].accent) >= 20, (a, b)
        assert T.delta_e(T.THEMES[a].window, T.THEMES[b].window) >= 10, (a, b)


def test_light_windows_are_distinct_shades():
    """Sepia vs Solarized Light and Light vs Gen X Soft Club used to sit ΔE
    2-3 apart -- indistinguishable side by side."""
    assert T.delta_e(T.THEMES["Sepia (paper)"].window,
                     T.THEMES["Solarized Light"].window) >= 8
    assert T.delta_e(T.THEMES["Light"].window,
                     T.AESTHETIC_THEMES["Gen X Soft Club"].window) >= 7
    assert T.delta_e(T.THEMES["Dark"].window, T.THEMES["Slate"].window) >= 8


def test_contrast_floor_still_holds_with_margin():
    """The identity rework moved colours; the readability floor of
    test_theme.py must hold for every theme, including the per-theme site
    palettes (Solarized, Nord) on their own plot backgrounds."""
    C = T.contrast
    for name, t in ALL.items():
        assert C(t.text, t.base) >= 4.5 and C(t.text, t.window) >= 4.5, name
        assert C(t.text, t.hover) >= 4.5, name
        assert C(t.text_dim, t.window) >= 3.0, name
        assert C(t.accent_text, t.accent) >= 4.5, name
        assert C(t.accent, t.base) >= 3.0, name
        assert len(t.series) >= 10, name
        for s in t.series:
            assert C(s, t.plot_bg) >= 3.0, (name, s)
        for role in ("model", "baseline", "pivot", "measure"):
            assert C(getattr(t, role), t.plot_bg) >= 3.0, (name, role)


# ------------------------------------------------------------ (c) the plot
@pytest.fixture(scope="module")
def qapp():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _qcolor_hex(c) -> str:
    return c.name().lower()


@pytest.mark.parametrize("name", ["Nord", "Solarized Dark", "Ocean", "Sepia (paper)",
                                  "Y2K", "Vaporwave"])
def test_plot_follows_the_theme(qapp, name):
    """A Nord plot on a Nord background: background, axis pens, grid alpha,
    the experiment / model / residual / baseline pens and the site colours
    all come from the active theme."""
    import pyqtgraph as pg
    from larmor.desktop import plot as P

    previous = T.active().name
    try:
        # app=None: set the active theme and pyqtgraph's defaults WITHOUT
        # restyling the QApplication -- QApplication.setStyleSheet re-polishes
        # every live widget, and late in the full suite (thousands of widgets
        # from earlier tests still alive) 14 themes x 2 calls took > 30 min
        # and looked like a hang (faulthandler dumps at theme.apply). The
        # view built next reads theme.active(), which is what is under test.
        t = T.apply(None, name)
        view = P.SpectrumView()
        pi = view.getPlotItem()
        assert _qcolor_hex(view.backgroundBrush().color()) == t.plot_bg
        assert _qcolor_hex(pi.getAxis("bottom").pen().color()) == t.axis
        assert _qcolor_hex(pi.getAxis("left").pen().color()) == t.axis
        assert _qcolor_hex(pi.getAxis("bottom").textPen().color()) == t.axis
        assert _qcolor_hex(pi.getAxis("top").pen().color()) == t.axis_minor
        assert _qcolor_hex(view._exp.opts["pen"].color()) == t.experiment == t.text
        assert _qcolor_hex(view._model.opts["pen"].color()) == t.model
        assert _qcolor_hex(view._resid.opts["pen"].color()) == t.resid
        assert _qcolor_hex(view._bl_curve.opts["pen"].color()) == t.baseline
        assert pi.ctrl.gridAlphaSlider.value() == int(t.grid_alpha * 255)
        # pyqtgraph's global defaults for plots created later follow too
        assert pg.getConfigOption("background") == t.plot_bg
        assert pg.getConfigOption("foreground") == t.axis
        # site colours come from the theme's own series, readable on its plot
        for i in range(10):
            assert P.site_color(i) == t.series[i]
            assert T.contrast(P.site_color(i), t.plot_bg) >= 3.0
        # a baseline anchor's fill follows the theme too (no white dot on Nord)
        view.set_baseline_mode(True)
        view._add_baseline_anchor(0.0, 0.0)
        anchor = view._bl_anchors[-1]
        assert _qcolor_hex(anchor.brush.color()) == t.base
        view.clear_baseline()
        view.close()
    finally:
        T.apply(None, previous)             # symmetric: the app was never restyled


def test_swatch_icon_shows_window_accent_and_plot(qapp):
    from PySide6.QtCore import Qt

    for t in ALL.values():
        icon = T.swatch_icon(t, 30, 14)
        pm = icon.pixmap(30, 14)
        img = pm.toImage()
        assert img.width() == 30 and img.height() == 14
        # sample the middle of each band, inside the hairline border
        assert img.pixelColor(5, 7).name().lower() == t.window, t.name
        assert img.pixelColor(15, 7).name().lower() == t.accent, t.name
        assert img.pixelColor(25, 7).name().lower() == t.plot_bg, t.name
        assert img.pixelColor(0, 0).name().lower() == t.border, t.name
        assert img.pixelColor(0, 0).alpha() == 255
        assert img.pixelFormat().alphaUsage() == img.pixelFormat().alphaUsage()  # sanity
    d = T.describe(T.THEMES["Nord"])
    assert "Nord" in d and "#3b4252" in d and "#88c0d0" in d and "#2e3440" in d
    assert Qt is not None


# ------------------------------------------------------------ the chooser
def test_theme_menu_previews_every_theme_and_states_the_restart(qapp, monkeypatch):
    """View ▸ Theme: every entry (the ten presets and the four styles) carries
    a swatch icon and a tooltip naming its scheme; the More styles… chooser
    states, in a section heading and in every tooltip, that a style is
    applied now and needs a restart for every dialog to follow."""
    from PySide6.QtWidgets import QMenu

    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    try:
        for act in win._theme_group.actions():
            t = T.get(act.text())
            assert not act.icon().isNull(), act.text()
            assert t.scheme in act.toolTip() and t.accent in act.toolTip(), act.text()
        more = None
        for menu in win.findChildren(QMenu):
            if menu.title() == "More styles…":
                more = menu
        assert more is not None and more.toolTipsVisible()
        sections = [a.text() for a in more.actions() if a.isSeparator() and a.text()]
        assert any("restart" in s.lower() for s in sections), sections
        assert any("swatch" in s.lower() for s in sections), sections
        for act in win._aesthetic_group.actions():
            if act.text() == "Normal":
                assert act.icon().isNull()
                continue
            t = T.get(act.text())
            assert not act.icon().isNull(), act.text()
            assert t.scheme in act.toolTip() and "restart" in act.toolTip().lower()
        # the menu TREE is unchanged: sections are separators, not entries
        entries = [a.text() for a in more.actions() if not a.isSeparator()]
        assert entries == ["Normal"] + T.aesthetic_names()
    finally:
        win.close()
