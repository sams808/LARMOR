"""Render a bundled Markdown manual in a scrollable dialog."""
from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtWidgets import (
    QDialog, QPushButton, QTextBrowser, QVBoxLayout,
)


def help_path(name: str) -> Path | None:
    """Locate larmor/help/<name>.md, whether from source or a frozen build."""
    here = Path(__file__).resolve()
    for base in (here.parent.parent / "help",                 # larmor/help
                 Path(getattr(sys, "_MEIPASS", "")) / "larmor" / "help"):
        p = base / f"{name}.md"
        if p.is_file():
            return p
    return None


def show_help(parent, name: str, title: str = "Help") -> None:
    p = help_path(name)
    text = (p.read_text(encoding="utf-8") if p
            else f"# {name}\n\nManual not found.")
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.resize(880, 720)
    v = QVBoxLayout(dlg)
    tb = QTextBrowser()
    tb.setOpenExternalLinks(True)
    tb.setStyleSheet("QTextBrowser { background: #ffffff; padding: 6px 10px; }")
    try:
        from larmor.desktop.mdrender import HELP_CSS, render_help_html

        tb.document().setDefaultStyleSheet(HELP_CSS)
        tb.setHtml(render_help_html(text))           # typeset the equations
    except Exception:
        try:
            tb.setMarkdown(text)
        except Exception:
            tb.setPlainText(text)
    v.addWidget(tb)
    btn = QPushButton("Close")
    btn.clicked.connect(dlg.accept)
    v.addWidget(btn)
    dlg.exec()
