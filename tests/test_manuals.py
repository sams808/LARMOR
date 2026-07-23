"""Every user manual must ship, render (equations included), and be reachable
from the ? ▸ User manuals menu."""
import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# (file stem, menu title) — must match the ? ▸ User manuals registration.
MANUALS = [
    ("getting-started", "Getting started"),
    ("spectra-1d", "1D spectra — processing & fitting"),
    ("2d-processing", "2D processing"),
    ("mqmas", "MQMAS (2D)"),
    ("correlation-hmqc", "HMQC & correlation"),
    ("relaxation", "Relaxation (T1/T2)"),
    ("qcpmg", "QCPMG"),
    ("multi-dataset", "Multi-dataset & co-fitting"),
    ("processing-reference", "Processing reference"),
]


@pytest.mark.parametrize("stem,_title", MANUALS)
def test_manual_ships_and_is_substantial(stem, _title):
    from larmor.desktop.help_dialog import help_path

    p = help_path(stem)
    assert p is not None, f"larmor/help/{stem}.md is missing"
    text = p.read_text(encoding="utf-8")
    assert len(text) > 1500, f"{stem}.md looks too thin to be a real manual"
    assert "## References" in text or "references" in text.lower(), \
        f"{stem}.md has no references section"
    assert text.rstrip().endswith(
        "*LARMOR — Sam Soudani, McCloy group, Washington State University.*"), \
        f"{stem}.md missing the standard footer"


@pytest.mark.parametrize("stem,_title", MANUALS)
def test_manual_renders_to_html(stem, _title):
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from larmor.desktop.help_dialog import help_path
    from larmor.desktop.mdrender import render_help_html

    text = help_path(stem).read_text(encoding="utf-8")
    html = render_help_html(text)                  # must not raise
    assert "MDMATH" not in html, "an equation token was left unreplaced"
    assert "$$" not in html, "a display equation was left as raw LaTeX"
    assert "<p" in html or "<img" in html          # produced real HTML


# Every model + reference manual, not just the parametrized user manuals.
ALL_HELP = [s for s, _ in MANUALS] + ["lineshapes"]


@pytest.mark.parametrize("stem", ALL_HELP)
def test_every_equation_typesets(stem):
    """Each $…$ / $$…$$ must render — an unsupported LaTeX construct (\\tfrac,
    \\big, \\le, …) would silently pass through as raw text otherwise."""
    import re

    from larmor.desktop.help_dialog import help_path
    from larmor.desktop.mdrender import _math_png, _protect_code

    md, _ = _protect_code(help_path(stem).read_text(encoding="utf-8"))
    tex = [m.group(1).strip()
           for m in re.finditer(r"\$\$(.+?)\$\$", md, flags=re.DOTALL)]
    inline_src = re.sub(r"\$\$(.+?)\$\$", " ", md, flags=re.DOTALL)
    tex += [m.group(1).strip() for m in
            re.finditer(r"(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)", inline_src,
                        flags=re.DOTALL)]
    broken = [t for t in tex if not _math_png(t)[0]]
    assert not broken, f"{stem}.md has un-typesettable equations: {broken}"


def test_all_manuals_in_help_menu():
    from PySide6.QtWidgets import QApplication, QMenu

    os.environ["LARMOR_NO_SESSION"] = "1"
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    try:
        man = next(m for m in win.menuBar().findChildren(QMenu)
                   if m.title() == "User &manuals")
        labels = [a.text() for a in man.actions() if a.text()]
        for _stem, title in MANUALS:
            assert any(title in x for x in labels), f"{title} not in menu"
    finally:
        win.close()
