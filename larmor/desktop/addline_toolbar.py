"""The "＋ Add line" tool on the main toolbar.

One prominent split button carries the CURRENT model ("＋ Add line: Czjzek");
its main part toggles the placing mode for that model, its arrow opens a menu
of every model grouped by use (simple lines · quadrupolar, ordered ·
disordered / glasses · CSA · other), each entry with a one-line description
and the key parameters as tooltip. Next to it sit quick buttons for the most
used models, and a "placing … click on the spectrum, Esc to stop" label that
shows while a placing mode is on.

The menu and the quick buttons hold the SAME checkable ``QAction`` objects as
``MainWindow._model_actions`` (built with the Decomposition ▸ Add line
submenu), so ``_set_add_mode``, the command palette, the self-test and every
test that checks a model action keep working unchanged; this module only
arranges them. The current model is remembered per nucleus in QSettings
(``addLine/model/<nucleus>``, with ``addLine/model`` as the fallback).

``build_addline_toolbar(win, toolbar)`` is the one call ``_build_toolbar``
makes; it returns the ``AddLineToolbar`` (kept as ``win._addline_toolbar``).
"""
from __future__ import annotations

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import QLabel, QMenu, QToolButton

from larmor import models as model_registry
from larmor.desktop import theme

#: (group title, registry names) -- the order of the split-button menu. A
#: registered model missing here is appended to "Other", never hidden.
GROUPS = (
    ("Simple lines", ("gauss_lor", "gl_norm", "voigt", "jmultiplet")),
    ("Quadrupolar (ordered)", ("quad_ct", "quad_first", "quad_csa")),
    ("Disordered / glasses", ("czjzek", "czjzek_d", "czjzek_corr", "ext_czjzek",
                              "amorphous")),
    ("CSA", ("csa_mas", "csa_czjzek")),
    ("Other", ("sidebands", "exchange2", "function")),
)
#: the quick buttons beside the split button, in order
QUICK = ("gauss_lor", "czjzek", "quad_ct", "csa_mas", "ext_czjzek")
DEFAULT_MODEL = "gauss_lor"
SETTINGS_KEY = "addLine/model"
#: the short names the tooltips list as "key parameters"
_SHORT = {"isotropic_chemical_shift_ppm": "position", "amplitude": "amplitude",
          "shift_fwhm_ppm": "width (dCS)", "line_fwhm_ppm": "lb",
          "sigma_Cq_MHz": "σ(Cq)", "Cq_MHz": "Cq", "eta": "η", "eta_q": "η(Q)",
          "eta_cs": "η(CSA)", "zeta_ppm": "ζ", "sigma_zeta_ppm": "σ(ζ)",
          "gl": "Gauss/Lorentz mix", "eps": "ε", "czjzek_d": "d",
          "shift_slope_ppm_per_MHz": "dδ/dC_Q", "Cq_fwhm_MHz": "ΔCq",
          "eta_fwhm": "Δη", "gauss_fwhm_ppm": "Gaussian width",
          "lorentz_fwhm_ppm": "Lorentzian width", "split_ppm": "Δδ A−B",
          "pop_a": "p(A)", "k_ex_hz": "k_ex", "j_hz": "J", "n_j": "n",
          "ssb_ratio": "sideband ratio", "n_ssb": "n(SSB)"}


def _first_sentence(text: str) -> str:
    text = " ".join(str(text or "").split())
    for stop in (". ", " -- ", " — "):
        if stop in text:
            return text.split(stop, 1)[0].rstrip(".") + "."
    return text


def model_tooltip(name: str) -> str:
    """'Czjzek (quad. distribution) — Czjzek distribution of … Key parameters:
    position, σ(Cq), width (dCS), amplitude (lb held).'"""
    try:
        m = model_registry.get(name)
    except ValueError:
        return name
    free = [_SHORT.get(p.name, p.key) for p in m.params if p.vary]
    held = [_SHORT.get(p.name, p.key) for p in m.params if not p.vary]
    tip = f"{m.label} — {_first_sentence(m.description)}"
    if free:
        tip += "\nKey parameters: " + ", ".join(free)
    if held:
        tip += " (" + ", ".join(held) + " held unless freed)"
    return tip


class AddLineToolbar:
    """The split button, its grouped menu, the quick buttons and the placing
    label; owns no model action, only arranges the window's."""

    def __init__(self, win, toolbar):
        self.win = win
        self.actions: dict[str, QAction] = win._model_actions
        self._current = self._remembered()
        # the split button: main part toggles the current model, arrow = menu
        self.button = QToolButton()
        self.button.setObjectName("addLineButton")
        self.button.setPopupMode(QToolButton.MenuButtonPopup)
        self.button.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self.button.setCheckable(True)
        self.button.setCursor(Qt.PointingHandCursor)
        self.menu = QMenu(self.button)
        self.menu.setToolTipsVisible(True)
        self._fill_menu()
        self.button.setMenu(self.menu)
        self.button.clicked.connect(self._toggle_current)
        toolbar.addWidget(self.button)
        # the quick buttons: the same QAction objects, so their checked state
        # is the placing mode itself
        for name in QUICK:
            act = self.actions.get(name)
            if act is not None:
                toolbar.addAction(act)
        self.placing = QLabel()
        self.placing.setObjectName("addLinePlacing")
        self.placing.setVisible(False)
        toolbar.addWidget(self.placing)
        for act in self.actions.values():
            act.toggled.connect(self._refresh)
        self.menu.aboutToShow.connect(self._refresh)
        self._refresh()

    # ------------------------------------------------------------- menu
    def _fill_menu(self) -> None:
        placed = set()
        for title, names in GROUPS:
            group = [self.actions[n] for n in names if n in self.actions]
            extra = []
            if title == "Other":
                extra = [a for n, a in self.actions.items()
                         if n not in placed and n not in names]
                extra.sort(key=lambda a: a.text().lower())
            if not group and not extra and title != "Other":
                continue
            self.menu.addSection(title)
            for act in group + extra:
                name = next(n for n, a in self.actions.items() if a is act)
                act.setToolTip(model_tooltip(name))
                self.menu.addAction(act)
                placed.add(name)
            if title == "Other":
                # the background spectrum is a dialog, not a placing mode
                a_bg = QAction("Spectrum (background)…", self.menu)
                a_bg.setToolTip(model_tooltip("spectrum")
                                + "\nOpens Edit ▸ Add background spectrum…")
                a_bg.triggered.connect(lambda *_: self.win.add_background_spectrum())
                self.menu.addAction(a_bg)
                self.action_background = a_bg

    # ------------------------------------------------------------ state
    def _nucleus(self) -> str:
        try:
            return str((self.win.recipe or {}).get("nucleus") or "")
        except Exception:
            return ""

    def _remembered(self) -> str:
        s = QSettings("LARMOR", "app")
        nuc = self._nucleus()
        for key in ((f"{SETTINGS_KEY}/{nuc}",) if nuc else ()) + (SETTINGS_KEY,):
            v = s.value(key, "")
            if isinstance(v, str) and v in self.actions:
                return v
        return DEFAULT_MODEL if DEFAULT_MODEL in self.actions else next(iter(self.actions), "")

    def _remember(self, name: str) -> None:
        s = QSettings("LARMOR", "app")
        s.setValue(SETTINGS_KEY, name)
        nuc = self._nucleus()
        if nuc:
            s.setValue(f"{SETTINGS_KEY}/{nuc}", name)

    def current_model(self) -> str:
        """The model the split button adds: the one being placed, else the
        last one chosen (remembered per nucleus)."""
        checked = [n for n, a in self.actions.items() if a.isChecked()]
        return checked[0] if checked else self._current

    def placing_model(self) -> str | None:
        checked = [n for n, a in self.actions.items() if a.isChecked()]
        return checked[0] if checked else None

    def _toggle_current(self, *_) -> None:
        """The main part of the split button: start placing the current model,
        or stop when it is being placed (QAction.trigger toggles the checkable
        action, which runs _set_add_mode through its triggered slot)."""
        name = self.current_model()
        act = self.actions.get(name)
        if act is None:
            return
        act.trigger()

    def _label(self, name: str) -> str:
        try:
            return model_registry.get(name).label
        except ValueError:
            return name

    def _refresh(self, *_) -> None:
        placing = self.placing_model()
        if placing is not None:
            self._current = placing
            self._remember(placing)
        name = self.current_model()
        t = theme.active()
        self.button.blockSignals(True)
        self.button.setChecked(placing is not None)
        self.button.blockSignals(False)
        self.button.setText(f"＋ Add line: {self._label(name)}")
        self.button.setToolTip(
            ("click on the spectrum to place lines; click again or press Esc to stop"
             if placing else f"click, then click on the spectrum to add a "
             f"{self._label(name)} line; the arrow picks another model")
            + "\n\n" + model_tooltip(name))
        self.button.setStyleSheet(
            f"QToolButton#addLineButton {{ font-weight: 600; padding: 2px 8px; "
            f"border: 1px solid {t.accent if placing else t.border}; "
            f"border-radius: 3px; }}"
            f"QToolButton#addLineButton:checked {{ background: {t.accent}; "
            f"color: {theme.best_text_on(t.accent)}; }}")
        if placing is not None:
            self.placing.setText(f"placing {self._label(placing)} — click on the "
                                 "spectrum · Esc to stop")
            self.placing.setStyleSheet(
                f"color: {t.accent}; font-weight: 600; padding: 0 8px;")
        self.placing.setVisible(placing is not None)

    def apply_theme(self) -> None:
        self._refresh()


def build_addline_toolbar(win, toolbar) -> AddLineToolbar:
    """Add the "＋ Add line" tool to ``toolbar``; the one call _build_toolbar
    makes."""
    return AddLineToolbar(win, toolbar)
