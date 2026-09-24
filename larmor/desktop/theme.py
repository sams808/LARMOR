"""Colour themes for the LARMOR desktop app.

Ten presets the user can switch between (View ▸ Theme), plus four hidden
"aesthetic" styles (View ▸ Theme ▸ More styles…). Each theme is defined by a
handful of *primitive* colours; every other role (surfaces, borders, hovers,
disabled text, plot axes/curves) is **derived** by blending, so a theme is
always internally consistent. Site colours and the signal marks (model /
baseline / pivot / measure) come from contrast-checked palettes — shared light
and dark ones, and a scheme's own accents where it has them (Solarized, Nord)
— so lines stay distinguishable and readable whatever the theme.

A theme named after a published scheme uses that scheme's canonical values
(Solarized, Nord); a theme without one is given an identity it cannot be
mistaken for (Ocean = deep blue-teal water with a foam accent, Slate = cool
blue-grey with brass, Sepia = paper and ink, Y2K = glossy silver and lime,
Vaporwave = pink and cyan on indigo). The comments on each theme say where
a canonical value had to move for the contrast floor, and by how much.

`tests/test_theme.py` checks that every text-on-background and
mark-on-plot-background pair clears a WCAG contrast floor, so "keep everything
visible and readable" is enforced, not eyeballed; `tests/test_theme_identity.py`
pins each theme to its scheme and holds every pair of themes a perceptual
distance apart (CIE Lab ΔE, see `delta_e`), so "the themes are what they say
and can be told apart" is enforced too.
"""
from __future__ import annotations

from dataclasses import dataclass, field


# ---- colour maths ---------------------------------------------------------
def _rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _hex(rgb: tuple[float, float, float]) -> str:
    return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(c)))) for c in rgb)


def blend(a: str, b: str, t: float) -> str:
    """(1-t)*a + t*b — t=0 gives a, t=1 gives b."""
    ra, ga, ba = _rgb(a)
    rb, gb, bb = _rgb(b)
    return _hex((ra + (rb - ra) * t, ga + (gb - ga) * t, ba + (bb - ba) * t))


def _lin(c: int) -> float:
    s = c / 255.0
    return s / 12.92 if s <= 0.03928 else ((s + 0.055) / 1.055) ** 2.4


def luminance(h: str) -> float:
    r, g, b = _rgb(h)
    return 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b)


def contrast(a: str, b: str) -> float:
    la, lb = luminance(a), luminance(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def best_text_on(bg: str) -> str:
    """Black or white, whichever reads better on `bg`."""
    return "#111111" if contrast("#111111", bg) >= contrast("#ffffff", bg) else "#ffffff"


# ---- perceptual distance (CIE L*a*b*, D65) ---------------------------------
# WCAG contrast says whether text READS; it says nothing about whether two
# themes LOOK different. That is a perceptual question, answered in CIE Lab:
# a ΔE*ab (CIE76) of ~2 is a just-noticeable difference side by side, ~6 is
# a clearly different shade at a glance, ~25 is a different colour.
# tests/test_theme_identity.py holds every pair of themes to such distances.
def lab(h: str) -> tuple[float, float, float]:
    """sRGB hex -> CIE L*a*b* (D65 white, 2° observer)."""
    r, g, b = (_lin(c) for c in _rgb(h))
    x = (0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / 0.95047
    y = (0.2126729 * r + 0.7151522 * g + 0.0721750 * b)
    z = (0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / 1.08883

    def f(t: float) -> float:
        return t ** (1.0 / 3.0) if t > 216.0 / 24389.0 else (841.0 / 108.0) * t + 4.0 / 29.0

    fx, fy, fz = f(x), f(y), f(z)
    return 116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)


def delta_e(a: str, b: str) -> float:
    """CIE76 colour difference ΔE*ab between two hex colours."""
    la, lb = lab(a), lab(b)
    return sum((p - q) ** 2 for p, q in zip(la, lb)) ** 0.5


def hue_chroma(h: str) -> tuple[float, float]:
    """(hue angle in degrees 0-360, chroma) of a hex colour in CIE Lab --
    the two numbers that say "these accents are both teal"."""
    import math

    _, a, b = lab(h)
    return math.degrees(math.atan2(b, a)) % 360.0, math.hypot(a, b)


def hue_difference(a: str, b: str) -> float:
    """Smallest angle (degrees) between the Lab hues of two colours."""
    ha, _ = hue_chroma(a)
    hb, _ = hue_chroma(b)
    d = abs(ha - hb) % 360.0
    return min(d, 360.0 - d)


# ---- shared, contrast-checked palettes ------------------------------------
# categorical site colours: deep on light backgrounds, bright on dark ones
LIGHT_SERIES = ["#0072b2", "#d55e00", "#009e73", "#b0568c", "#8f6e00",
                "#0e7c86", "#7f3fbf", "#117733", "#882255", "#5a5a5a"]
DARK_SERIES = ["#4fa3e0", "#ff8c42", "#3fd07f", "#e58fbf", "#f2c14e",
               "#4bd6df", "#b48ee8", "#5fd3a0", "#e06c9f", "#b7c0c8"]
# signal marks (model / baseline / pivot / measure)
LIGHT_SIGNAL = dict(model="#c0392b", baseline="#96650f", pivot="#8e44ad",
                    measure="#0e7c86")
DARK_SIGNAL = dict(model="#ff6b6b", baseline="#f2b134", pivot="#c58cff",
                   measure="#3fd0d8")


@dataclass(frozen=True)
class Theme:
    name: str
    is_dark: bool
    window: str          # main background
    base: str            # inputs / tables / menus
    text: str            # primary ink
    text_dim: str        # secondary ink (headers, axes)
    border: str
    accent: str          # highlight / selection / default button
    plot_bg: str
    series: list[str] = field(default_factory=list)
    #: extra QSS appended after the normal parametric stylesheet -- empty for
    #: every normal theme (zero behaviour change); a hidden "aesthetic" theme
    #: (see AESTHETIC_THEMES) uses this for a decorative flourish (gradients,
    #: glow) beyond what the flat primitive-colour system expresses
    flourish_qss: str = ""
    #: the reference scheme the colours are taken from (or the identity a
    #: theme without one is given) -- shown in the theme menu's tooltips and
    #: pinned by tests/test_theme_identity.py
    scheme: str = ""

    # ---- derived roles (consistent by construction) ----
    @property
    def accent_text(self) -> str:
        return best_text_on(self.accent)

    @property
    def surface(self) -> str:          # menubar / toolbar / statusbar
        return blend(self.window, self.base, 0.35)

    @property
    def header(self) -> str:           # table headers / dock titles / tabs
        return blend(self.window, self.text, 0.06)

    @property
    def alt_base(self) -> str:         # alternate table rows
        return blend(self.base, self.text, 0.045)

    @property
    def border_soft(self) -> str:
        return blend(self.border, self.base, 0.5)

    @property
    def hover(self) -> str:            # menu/selection/hover wash
        return blend(self.base, self.accent, 0.20)

    @property
    def disabled_text(self) -> str:
        return blend(self.text, self.window, 0.5)

    # ---- plot roles ----
    @property
    def experiment(self) -> str:
        return self.text

    @property
    def resid(self) -> str:
        return blend(self.text, self.plot_bg, 0.55)

    @property
    def resid_zero(self) -> str:
        return blend(self.text, self.plot_bg, 0.78)

    @property
    def axis(self) -> str:
        return self.text_dim

    @property
    def axis_minor(self) -> str:
        return blend(self.text_dim, self.plot_bg, 0.55)

    @property
    def grid_alpha(self) -> float:
        return 0.16 if self.is_dark else 0.09

    @property
    def legend_bg(self) -> str:
        return self.base

    @property
    def _signal(self) -> dict:
        return DARK_SIGNAL if self.is_dark else LIGHT_SIGNAL

    @property
    def model(self) -> str:
        return self._signal["model"]

    @property
    def baseline(self) -> str:
        return self._signal["baseline"]

    @property
    def pivot(self) -> str:
        return self._signal["pivot"]

    @property
    def measure(self) -> str:
        return self._signal["measure"]


def _theme(name, is_dark, window, base, text, text_dim, border, accent, plot_bg,
          flourish_qss="", series=None, scheme=""):
    return Theme(name=name, is_dark=is_dark, window=window, base=base, text=text,
                 text_dim=text_dim, border=border, accent=accent, plot_bg=plot_bg,
                 series=list(series or (DARK_SERIES if is_dark else LIGHT_SERIES)),
                 flourish_qss=flourish_qss, scheme=scheme)


# ---- canonical schemes ------------------------------------------------------
# Solarized (Ethan Schoonover, 2011): base03 #002b36, base02 #073642,
# base01 #586e75, base00 #657b83, base0 #839496, base1 #93a1a1,
# base2 #eee8d5, base3 #fdf6e3; yellow #b58900, orange #cb4b16, red #dc322f,
# magenta #d33682, violet #6c71c4, blue #268bd2, cyan #2aa198, green #859900.
# The site palette of the dark mode is the eight accents plus base1/base0,
# every one of which clears the 3:1 floor on base03.
SOLARIZED_DARK_SERIES = ["#268bd2", "#cb4b16", "#859900", "#d33682", "#b58900",
                         "#2aa198", "#6c71c4", "#dc322f", "#93a1a1", "#839496"]
# On base3 (#fdf6e3) three canonical accents miss the 3:1 floor by a hair
# (green 2.97, yellow 2.98, cyan 2.93): they are darkened by 2-3 %, the
# smallest step that clears 3.05 (ΔE ≈ 1.5-1.9, invisible side by side):
# green #859900 -> #829600, yellow #b58900 -> #b18600, cyan #2aa198 -> #299c93.
# The greys are base01 / base00, the scheme's own emphasised / body ink.
SOLARIZED_LIGHT_SERIES = ["#268bd2", "#cb4b16", "#829600", "#d33682", "#b18600",
                          "#299c93", "#6c71c4", "#dc322f", "#586e75", "#657b83"]
# Nord (Arctic Ice Studio): polar night nord0-3 #2e3440 #3b4252 #434c5e
# #4c566a, snow storm nord4-6 #d8dee9 #e5e9f0 #eceff4, frost nord7-10 #8fbcbb
# #88c0d0 #81a1c1 #5e81ac, aurora nord11-15 #bf616a #d08770 #ebcb8b #a3be8c
# #b48ead. The site palette is aurora + frost + nord4 on a nord0 plot; the
# darkest two (nord11 3.05, nord10 3.10) sit just above the 3:1 floor.
NORD_SERIES = ["#88c0d0", "#d08770", "#a3be8c", "#b48ead", "#ebcb8b",
               "#bf616a", "#81a1c1", "#8fbcbb", "#5e81ac", "#d8dee9"]


# ---- the ten presets ------------------------------------------------------
# Every theme is meant to be recognisable at a glance AND to be what its name
# says. tests/test_theme_identity.py pins a few key colours of each to the
# scheme it is named after and holds every pair of themes apart in CIE Lab
# (window ΔE ≥ 6 or accent ΔE ≥ 25; no two themes of the same polarity with
# both a near-identical window and a same-hue accent). Where a canonical
# value had to move for the contrast floor, the comment says so and by how
# much.
THEMES: dict[str, Theme] = {t.name: t for t in [
    # LARMOR's own neutral light: near-white surfaces, dark blue-grey ink and
    # the brand teal of the logo as the accent
    _theme("Light", False, "#f0f2f0", "#ffffff", "#16202a", "#37424a",
           "#cfd6d1", "#0e7c86", "#fcfdfc",
           scheme="LARMOR neutral light — brand teal accent"),
    # warm paper and ink: aged-paper chrome (#e6d7b8) around parchment pages
    # (#f7f1e3 -- the paper must stay this light for the shared site palette
    # to clear 3:1 on it), dark-brown ink, faded-ink secondary text and a
    # saddle-brown accent
    _theme("Sepia (paper)", False, "#e6d7b8", "#f7f1e3", "#3b2a1a", "#6b5236",
           "#c9b48d", "#8b4513", "#f7f1e3",
           scheme="Sepia — aged paper, brown ink, saddle-brown accent"),
    # Solarized Light: base2 chrome around base3 pages, base1 borders, the
    # Solarized blue as accent (black label text on it, 5.1:1). Primary ink
    # is base02 rather than the scheme's base01/base00 body text: base01
    # reads 4.3:1 on base2 and base00 4.0:1 on base3, both under the 4.5
    # floor, so the darkest canonical tone carries the text and base01
    # becomes the secondary ink
    _theme("Solarized Light", False, "#eee8d5", "#fdf6e3", "#073642", "#586e75",
           "#93a1a1", "#268bd2", "#fdf6e3", series=SOLARIZED_LIGHT_SERIES,
           scheme="Solarized Light (Schoonover) — base3/base2 paper, blue accent"),
    _theme("High Contrast Light", False, "#ffffff", "#ffffff", "#000000",
           "#1a1a1a", "#3a3a3a", "#0033cc", "#ffffff",
           scheme="High contrast — pure white, black ink, cobalt accent"),
    # LARMOR's own neutral dark: graphite surfaces with the brand teal
    _theme("Dark", True, "#262b2e", "#1e2225", "#e6ebe8", "#aab4ad",
           "#3a4247", "#2dd4bf", "#1b1f22",
           scheme="LARMOR neutral dark — graphite, brand teal accent"),
    # cool blue-grey slate (a definite blue cast, unlike the neutral Dark)
    # with a brass accent -- the one warm highlight among the dark themes,
    # so it is never mistaken for Nord or Ocean
    _theme("Slate", True, "#263241", "#1c2634", "#d7e0ea", "#9fb0c3",
           "#46586c", "#e3b341", "#1c2634",
           scheme="Slate — cool blue-grey, brass accent"),
    # Nord: nord1 chrome around nord0 pages and plot, nord6 ink, nord9
    # (frost) secondary ink -- nord3, the scheme's comment colour, reads
    # 1.7:1 and cannot serve as text -- nord3 borders and nord8, the
    # scheme's primary accent
    _theme("Nord", True, "#3b4252", "#2e3440", "#eceff4", "#81a1c1",
           "#4c566a", "#88c0d0", "#2e3440", series=NORD_SERIES,
           scheme="Nord (Arctic Ice Studio) — polar night, snow storm, frost"),
    # deep ocean blue (bluer and brighter than Solarized's base02, so the two
    # never blur), foam-white ink, sea-glass secondary ink and a foam accent
    _theme("Ocean", True, "#0b3d5c", "#072b42", "#e3f6f5", "#8fc7cf",
           "#1c5a7a", "#bff7ef", "#072b42",
           scheme="Ocean — deep blue-teal water, sea-foam accent"),
    # Solarized Dark: base02 chrome around base03 pages, base1 secondary
    # ink, base01 borders, the Solarized blue as accent. Primary ink is
    # base2: base1, the scheme's emphasised text, reads 4.1:1 on the
    # selection wash (base03 blended 20 % toward the accent), under the 4.5
    # floor, so the next canonical tone up carries the text
    _theme("Solarized Dark", True, "#073642", "#002b36", "#eee8d5", "#93a1a1",
           "#586e75", "#268bd2", "#002b36", series=SOLARIZED_DARK_SERIES,
           scheme="Solarized Dark (Schoonover) — base03/base02, blue accent"),
    _theme("High Contrast Dark", True, "#000000", "#0a0a0a", "#ffffff",
           "#e6e6e6", "#8a8a8a", "#ffd400", "#000000",
           scheme="High contrast — pure black, white ink, yellow accent"),
]}

DEFAULT = "Light"
_ACTIVE = THEMES[DEFAULT]


# ---- hidden "aesthetic" themes ---------------------------------------------
# A just-for-fun easter egg (View ▸ Theme ▸ More styles…, restart required —
# see app.py) — deliberately NOT in THEMES/names() so the normal Theme menu
# stays exactly the 10 presets above. Same Theme dataclass, same derived
# roles, same contrast floor as every other theme (see test_theme.py) — a
# "fun" theme that's illegible isn't worth shipping — plus a `flourish_qss`
# snippet for the bit of personality flat primitive colours can't express.
# Y2K: the glossy brushed-silver chrome of 2000-era UIs with a lime accent.
# The bright lime itself (#9be32a) reads 1.6:1 on the silver page, far under
# the 3:1 accent floor, so it lives in the flourish -- hover borders and the
# button gloss -- while the accent proper is the darkest lime that still
# passes (#339400: 3.45:1 on the page, black label text 4.8:1).
_Y2K_QSS = """
QPushButton { border: 1px solid #9aa5b3; border-radius: 8px;
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #ffffff, stop:0.45 #e6eaef, stop:0.5 #d4dae2, stop:1 #c2c9d2); }
QPushButton:hover { border: 2px solid #9be32a; }
QPushButton:pressed { background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #c2c9d2, stop:1 #e6eaef); }
QPushButton:default { background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #7fd12a, stop:0.5 #4fae12, stop:1 #339400); border-color: #339400; }
QMenuBar, QToolBar { background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #eef1f5, stop:1 #c9d0d9); }
QGroupBox, QTabWidget::pane { border: 1px solid #9aa5b3; border-radius: 8px; }
QTabBar::tab:selected { border-bottom: 2px solid #9be32a; }
"""
# Vaporwave: hot pink and cyan on deep indigo (this one already read right)
_VAPORWAVE_QSS = """
QPushButton { border: 1px solid #ff2fd6; border-radius: 4px;
  background: qlineargradient(x1:0,y1:0,x2:1,y2:0, stop:0 #241238, stop:1 #2e1648); }
QPushButton:hover { border-color: #4be8ff; }
QMenuBar { border-bottom: 1px solid #ff2fd6; }
"""
# Gen X Soft Club: pastel mint, soft-focus rounded buttons, a rose accent
_SOFT_QSS = """
QPushButton { border: 1px solid #b8d9cc; border-radius: 10px;
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #f7fcf9, stop:1 #e6f4ed); }
QPushButton:hover { border-color: #c2617f; }
"""
# Dreamcore: a lilac haze, rounded buttons, a dream-violet accent
_DREAMCORE_QSS = """
QPushButton { border: 1px solid #d3c4e6; border-radius: 10px;
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #faf7fd, stop:1 #efe7f7); }
QPushButton:hover { border-color: #7d4fc9; }
"""

AESTHETIC_THEMES: dict[str, Theme] = {t.name: t for t in [
    # a LIGHT theme now: brushed-silver chrome, near-white pages, ink-dark
    # text, a lime accent -- the old dark-graphite-and-neon Y2K was a Matrix
    # look, not a 2000-era UI
    _theme("Y2K", False, "#cdd3da", "#eef1f5", "#10161d", "#4a5666",
           "#9aa5b3", "#339400", "#eef1f5", flourish_qss=_Y2K_QSS,
           scheme="Y2K — glossy brushed silver, lime accent"),
    # lilac haze: lavender chrome, pale-violet pages, deep-violet ink
    _theme("Dreamcore", False, "#e9e0f3", "#f6f1fb", "#3a2a4e", "#77678f",
           "#d3c4e6", "#7d4fc9", "#f8f4fc", flourish_qss=_DREAMCORE_QSS,
           scheme="Dreamcore — lilac haze, dream-violet accent"),
    # pastel mint chrome (deep enough to be told from Light at a glance),
    # mint-white pages, deep-green ink and a soft rose accent
    _theme("Gen X Soft Club", False, "#dcefe6", "#f3fbf7", "#23443a",
           "#5f8a7c", "#b8d9cc", "#c2617f", "#f4fbf7", flourish_qss=_SOFT_QSS,
           scheme="Gen X Soft Club — pastel mint, soft rose accent"),
    _theme("Vaporwave", True, "#1a0b2e", "#241238", "#f2e6ff", "#c9a0e8",
           "#3d1f57", "#ff2fd6", "#160a26", flourish_qss=_VAPORWAVE_QSS,
           scheme="Vaporwave — hot pink and cyan on indigo"),
]}


def names() -> list[str]:
    return list(THEMES)


def aesthetic_names() -> list[str]:
    """The hidden, just-for-fun styles (View ▸ Theme ▸ More styles…) —
    kept out of names()/THEMES so the normal Theme menu is unaffected."""
    return list(AESTHETIC_THEMES)


def active() -> Theme:
    return _ACTIVE


def get(name: str) -> Theme:
    return THEMES.get(name) or AESTHETIC_THEMES.get(name) or THEMES[DEFAULT]


def set_active(name: str) -> Theme:
    global _ACTIVE
    _ACTIVE = get(name)
    return _ACTIVE


# ---- Qt application styling ------------------------------------------------
def qss(t: Theme) -> str:
    """A complete stylesheet built from the theme's roles (parametric APP_STYLE)."""
    return f"""
QMainWindow {{ background: {t.window}; }}
QWidget {{ color: {t.text}; }}
QMenuBar {{ background: {t.surface}; color: {t.text}; border-bottom: 1px solid {t.border}; }}
QMenuBar::item {{ padding: 4px 10px; background: transparent; }}
QMenuBar::item:selected {{ background: {t.hover}; border-radius: 4px; }}
QMenu {{ background: {t.base}; color: {t.text}; border: 1px solid {t.border}; }}
QMenu::item {{ padding: 4px 26px 4px 18px; }}
QMenu::item:selected {{ background: {t.hover}; color: {t.text}; }}
QMenu::item:checked {{ font-weight: 600; }}
QMenu::separator {{ height: 1px; background: {t.border_soft}; margin: 4px 8px; }}
QToolBar {{ background: {t.surface}; border-bottom: 1px solid {t.border}; spacing: 3px; padding: 3px; }}
QToolBar#sidebar {{ border-right: 1px solid {t.border}; border-bottom: none; padding: 3px 2px; }}
QToolButton {{ padding: 4px 9px; border-radius: 4px; color: {t.text}; border: 1px solid transparent; }}
QToolButton:hover {{ background: {t.hover}; border-color: {t.border_soft}; }}
QToolButton:checked {{ background: {t.accent}; color: {t.accent_text}; }}
QDockWidget {{ color: {t.text}; }}
QDockWidget::title {{ background: {t.header}; color: {t.text}; padding: 4px 8px; border-top: 1px solid {t.border}; }}
QTableWidget {{ background: {t.base}; color: {t.text}; gridline-color: {t.border_soft};
               alternate-background-color: {t.alt_base}; selection-background-color: {t.accent};
               selection-color: {t.accent_text}; }}
QHeaderView::section {{ background: {t.header}; color: {t.text_dim}; font-weight: 600;
                       border: none; border-right: 1px solid {t.border_soft};
                       border-bottom: 1px solid {t.border}; padding: 3px 6px; }}
QTableCornerButton::section {{ background: {t.header}; border: none; }}
QDoubleSpinBox, QSpinBox, QLineEdit {{ color: {t.text}; background: {t.base};
                                      border: 1px solid {t.border_soft}; border-radius: 3px;
                                      selection-background-color: {t.accent};
                                      selection-color: {t.accent_text}; }}
QComboBox {{ color: {t.text}; background: {t.base}; border: 1px solid {t.border_soft};
            border-radius: 3px; padding: 2px 6px; }}
QComboBox QAbstractItemView {{ background: {t.base}; color: {t.text};
                              selection-background-color: {t.accent}; selection-color: {t.accent_text}; }}
QPushButton {{ color: {t.text}; background: {t.base}; border: 1px solid {t.border};
              border-radius: 4px; padding: 4px 14px; }}
QPushButton:hover {{ background: {t.hover}; }}
QPushButton:disabled {{ color: {t.disabled_text}; }}
QPushButton:default {{ background: {t.accent}; color: {t.accent_text}; border-color: {t.accent}; }}
QCheckBox {{ color: {t.text}; }}
QRadioButton {{ color: {t.text}; }}
QLabel {{ color: {t.text}; background: transparent; }}
QGroupBox {{ color: {t.text}; border: 1px solid {t.border_soft}; border-radius: 4px; margin-top: 6px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 8px; padding: 0 4px; }}
QTabBar::tab {{ background: {t.header}; color: {t.text_dim}; padding: 4px 14px;
               border: 1px solid {t.border}; border-bottom: none;
               border-top-left-radius: 4px; border-top-right-radius: 4px; }}
QTabBar::tab:selected {{ background: {t.base}; color: {t.accent}; font-weight: 600; }}
QStatusBar {{ background: {t.surface}; color: {t.text_dim}; border-top: 1px solid {t.border}; }}
QPlainTextEdit, QTextEdit, QTextBrowser {{ background: {t.base}; color: {t.text};
                                          selection-background-color: {t.accent};
                                          selection-color: {t.accent_text}; }}
QToolTip {{ background: {t.base}; color: {t.text}; border: 1px solid {t.border}; }}
QScrollBar:vertical {{ background: {t.header}; width: 12px; }}
QScrollBar::handle:vertical {{ background: {t.border}; border-radius: 5px; min-height: 30px; }}
QScrollBar:horizontal {{ background: {t.header}; height: 12px; }}
QScrollBar::handle:horizontal {{ background: {t.border}; border-radius: 5px; min-width: 30px; }}
QDialog {{ background: {t.window}; }}
""" + t.flourish_qss


def palette(t: Theme):
    """A full QPalette so widgets render identically regardless of the OS theme."""
    from PySide6.QtGui import QColor, QPalette

    p = QPalette()
    C = QColor
    p.setColor(QPalette.Window, C(t.window))
    p.setColor(QPalette.WindowText, C(t.text))
    p.setColor(QPalette.Base, C(t.base))
    p.setColor(QPalette.AlternateBase, C(t.alt_base))
    p.setColor(QPalette.Text, C(t.text))
    p.setColor(QPalette.PlaceholderText, C(t.disabled_text))
    p.setColor(QPalette.Button, C(t.surface))
    p.setColor(QPalette.ButtonText, C(t.text))
    p.setColor(QPalette.BrightText, C("#ff5555"))
    p.setColor(QPalette.ToolTipBase, C(t.base))
    p.setColor(QPalette.ToolTipText, C(t.text))
    p.setColor(QPalette.Highlight, C(t.accent))
    p.setColor(QPalette.HighlightedText, C(t.accent_text))
    p.setColor(QPalette.Link, C(t.accent))
    for grp in (QPalette.Disabled,):
        for role in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
            p.setColor(grp, role, C(t.disabled_text))
    return p


def describe(t: Theme) -> str:
    """One line for a tooltip: the scheme followed and the three colours the
    swatch shows (window, accent, plot background)."""
    head = f"{t.scheme} · " if t.scheme else ""
    return (f"{head}window {t.window} · accent {t.accent} · plot {t.plot_bg} · "
            f"{'dark' if t.is_dark else 'light'}")


def swatch_icon(t: Theme, width: int = 30, height: int = 14):
    """A QIcon previewing a theme: three bands -- window, accent, plot
    background -- inside a hairline of the theme's border colour, so the
    entries of the Theme menu can be told apart before one is picked."""
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap

    pm = QPixmap(width, height)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    third = width // 3
    for i, colour in enumerate((t.window, t.accent, t.plot_bg)):
        w = third if i < 2 else width - 2 * third
        p.fillRect(i * third, 0, w, height, QColor(colour))
    p.setPen(QColor(t.border))
    p.drawRect(0, 0, width - 1, height - 1)
    p.end()
    return QIcon(pm)


def apply(app, name: str) -> Theme:
    """Set the active theme and apply it to the QApplication + pyqtgraph."""
    import pyqtgraph as pg

    t = set_active(name)
    pg.setConfigOptions(background=t.plot_bg, foreground=t.axis)
    if app is not None:
        app.setPalette(palette(t))
        app.setStyleSheet(qss(t))
    return t
