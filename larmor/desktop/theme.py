"""Colour themes for the LARMOR desktop app.

Ten presets the user can switch between (View ▸ Theme). Each theme is defined by
a handful of *primitive* colours; every other role (surfaces, borders, hovers,
disabled text, plot axes/curves) is **derived** by blending, so a theme is always
internally consistent. Series (site) colours and the signal marks (model /
baseline / pivot / measure) come from two shared, contrast-checked palettes —
one for light plot backgrounds, one for dark — so lines stay distinguishable and
readable whatever the theme.

`tests/test_theme.py` checks that every text-on-background and
mark-on-plot-background pair clears a WCAG contrast floor, so "keep everything
visible and readable" is enforced, not eyeballed.
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
          flourish_qss=""):
    return Theme(name=name, is_dark=is_dark, window=window, base=base, text=text,
                 text_dim=text_dim, border=border, accent=accent, plot_bg=plot_bg,
                 series=(DARK_SERIES if is_dark else LIGHT_SERIES),
                 flourish_qss=flourish_qss)


# ---- the ten presets ------------------------------------------------------
THEMES: dict[str, Theme] = {t.name: t for t in [
    _theme("Light", False, "#f0f2f0", "#ffffff", "#16202a", "#37424a",
           "#cfd6d1", "#0e7c86", "#fcfdfc"),
    _theme("Sepia (paper)", False, "#ece3d0", "#f7f1e3", "#402f1d", "#6b5c45",
           "#d6c9ac", "#9a5b1f", "#f7f1e3"),
    _theme("Solarized Light", False, "#eee8d5", "#fdf6e3", "#073642", "#4a5b60",
           "#d9d2bf", "#1e6fa8", "#fdf6e3"),
    _theme("High Contrast Light", False, "#ffffff", "#ffffff", "#000000",
           "#1a1a1a", "#3a3a3a", "#0033cc", "#ffffff"),
    _theme("Dark", True, "#262b2e", "#1e2225", "#e6ebe8", "#aab4ad",
           "#3a4247", "#2dd4bf", "#1b1f22"),
    _theme("Slate", True, "#22272e", "#1b2027", "#d0d7de", "#9aa5b1",
           "#373e47", "#539bf5", "#1b2027"),
    _theme("Nord", True, "#2e3440", "#272c36", "#e5e9f0", "#b8c0cf",
           "#3b4252", "#88c0d0", "#2b303b"),
    _theme("Ocean", True, "#0f2231", "#0b1b27", "#d9e8f4", "#9cb6c9",
           "#244055", "#39b6c4", "#0b1b27"),
    _theme("Solarized Dark", True, "#073642", "#002b36", "#eee8d5", "#93a1a1",
           "#0f4a58", "#2aa198", "#002b36"),
    _theme("High Contrast Dark", True, "#000000", "#0a0a0a", "#ffffff",
           "#e6e6e6", "#8a8a8a", "#ffd400", "#000000"),
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
_Y2K_QSS = """
QPushButton { border: 1px solid #4a5a54; border-radius: 7px;
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #34423c, stop:1 #1b2224); }
QPushButton:hover { border-color: #39ff8f; }
QPushButton:pressed { background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #1b2224, stop:1 #263029); }
QGroupBox, QTabWidget::pane { border: 1px solid #33403c; border-radius: 6px; }
"""
_VAPORWAVE_QSS = """
QPushButton { border: 1px solid #ff2fd6; border-radius: 4px;
  background: qlineargradient(x1:0,y1:0,x2:1,y2:0, stop:0 #241238, stop:1 #2e1648); }
QPushButton:hover { border-color: #4be8ff; }
QMenuBar { border-bottom: 1px solid #ff2fd6; }
"""
_SOFT_QSS = """
QPushButton { border: 1px solid #cfe8de; border-radius: 10px;
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #f8fdfb, stop:1 #e9f6f0); }
QPushButton:hover { border-color: #3f9d7a; }
"""
_DREAMCORE_QSS = """
QPushButton { border: 1px solid #e6d5e8; border-radius: 10px;
  background: qlineargradient(x1:0,y1:0,x2:0,y2:1, stop:0 #faf5f8, stop:1 #f3e8f0); }
QPushButton:hover { border-color: #9b6bb3; }
"""

AESTHETIC_THEMES: dict[str, Theme] = {t.name: t for t in [
    _theme("Y2K", True, "#14181a", "#1b2224", "#e8fff2", "#7fd9a8",
          "#33403c", "#39ff8f", "#101314", flourish_qss=_Y2K_QSS),
    _theme("Dreamcore", False, "#f3ecf0", "#faf5f8", "#4a3b52", "#8a7a92",
          "#e6d5e8", "#9b6bb3", "#fbf7fa", flourish_qss=_DREAMCORE_QSS),
    _theme("Gen X Soft Club", False, "#eef6f2", "#f8fdfb", "#2f4d43",
          "#6b9686", "#cfe8de", "#3f9d7a", "#fbfefc", flourish_qss=_SOFT_QSS),
    _theme("Vaporwave", True, "#1a0b2e", "#241238", "#f2e6ff", "#c9a0e8",
          "#3d1f57", "#ff2fd6", "#160a26", flourish_qss=_VAPORWAVE_QSS),
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


def apply(app, name: str) -> Theme:
    """Set the active theme and apply it to the QApplication + pyqtgraph."""
    import pyqtgraph as pg

    t = set_active(name)
    pg.setConfigOptions(background=t.plot_bg, foreground=t.axis)
    if app is not None:
        app.setPalette(palette(t))
        app.setStyleSheet(qss(t))
    return t
