"""Figure studio dialog: JSON spec editor with live preview and export."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox,
    QPlainTextEdit, QPushButton, QScrollArea, QVBoxLayout,
)


class FigureDialog(QDialog):
    def __init__(self, parent, source_path: str | None, recipe: dict | None):
        super().__init__(parent)
        self.setWindowTitle("Figure studio")
        self.resize(980, 640)
        self.source_path, self.recipe = source_path, recipe

        root = QHBoxLayout(self)

        left = QVBoxLayout()
        tpl_row = QHBoxLayout()
        btnTpl = QPushButton("Templates for loaded data")
        btnTpl.clicked.connect(self.load_templates)
        tpl_row.addWidget(btnTpl)
        self.tpl_box = QHBoxLayout()
        tpl_row.addLayout(self.tpl_box)
        tpl_row.addStretch(1)
        left.addLayout(tpl_row)

        self.spec = QPlainTextEdit()
        self.spec.setStyleSheet("font-family: Consolas, monospace; font-size: 11px;")
        self.spec.setPlaceholderText(
            '{"kind": "1d", "style": "article-wide", "traces": [...]}\n'
            "see docs/tutorials/03-figures.md")
        left.addWidget(self.spec, 1)

        actions = QHBoxLayout()
        btnPrev = QPushButton("Preview")
        btnPrev.setDefault(True)
        btnPrev.clicked.connect(self.preview)
        actions.addWidget(btnPrev)
        self.fmt_png = QCheckBox("png"); self.fmt_png.setChecked(True)
        self.fmt_svg = QCheckBox("svg"); self.fmt_svg.setChecked(True)
        self.fmt_pdf = QCheckBox("pdf"); self.fmt_pdf.setChecked(True)
        for w in (self.fmt_png, self.fmt_svg, self.fmt_pdf):
            actions.addWidget(w)
        btnExp = QPushButton("Export…")
        btnExp.clicked.connect(self.export)
        actions.addWidget(btnExp)
        actions.addStretch(1)
        left.addLayout(actions)
        root.addLayout(left, 1)

        right = QScrollArea()
        right.setWidgetResizable(True)
        self.preview_label = QLabel("preview appears here")
        self.preview_label.setAlignment(Qt.AlignCenter)
        right.setWidget(self.preview_label)
        root.addWidget(right, 1)

    # ----------------------------------------------------------------
    def load_templates(self):
        if not self.source_path:
            return
        from larmor.figures import STYLES  # noqa: F401  (import check)

        p = Path(self.source_path)
        templates: dict[str, dict] = {}
        if p.suffix.lower() in (".fxmla", ".fxml") or \
                (p.is_dir() and not (p / "acqu2s").exists()):
            templates["1d"] = {
                "kind": "1d", "style": "article-wide",
                "traces": [{"path": str(p), "label": "experiment",
                            "color": "black", "linewidth": 0.9}]}
        if (p / "acqu2s").exists():
            templates["2d"] = {"kind": "2d", "style": "thesis", "path": str(p),
                               "levels": {"mode": "log", "n": 12}}
        pdata1 = p / "pdata" / "1"
        if (pdata1 / "t1ints.txt").exists():
            templates["satrec"] = {"kind": "series", "mode": "satrec",
                                   "style": "article", "path": str(p),
                                   "stretched": True}
        if (pdata1 / "redor.txt").exists():
            templates["redor"] = {"kind": "series", "mode": "redor",
                                  "style": "article", "path": str(p)}
        while self.tpl_box.count():
            item = self.tpl_box.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        for name, spec in templates.items():
            b = QPushButton(name)
            b.clicked.connect(lambda _, s=spec: self.spec.setPlainText(
                json.dumps(s, indent=2)))
            self.tpl_box.addWidget(b)

    def _spec(self) -> dict | None:
        try:
            return json.loads(self.spec.toPlainText())
        except Exception as exc:
            QMessageBox.warning(self, "Invalid spec", f"not valid JSON: {exc}")
            return None

    def preview(self):
        spec = self._spec()
        if spec is None:
            return
        from larmor import figures

        try:
            png = figures.render_png_bytes(spec, dpi=110)
        except Exception as exc:
            QMessageBox.warning(self, "Figure failed", str(exc))
            return
        pix = QPixmap()
        pix.loadFromData(png)
        self.preview_label.setPixmap(pix)

    def export(self):
        spec = self._spec()
        if spec is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export figure (base name, extensions added)", "",
            "Figure base name (*)")
        if not path:
            return
        base = Path(path).with_suffix("")
        for parent in [base.parent, *base.parent.parents]:
            if (parent / "acqus").exists() or (parent / "fid").exists() \
                    or (parent / "ser").exists():
                QMessageBox.warning(self, "Refused",
                                    f"{parent} is an instrument data folder.")
                return
        formats = [f for f, cb in (("png", self.fmt_png), ("svg", self.fmt_svg),
                                   ("pdf", self.fmt_pdf)) if cb.isChecked()]
        from larmor import figures

        try:
            saved = figures.export(spec, base, formats=tuple(formats))
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
            return
        QMessageBox.information(self, "Exported", "\n".join(saved))
