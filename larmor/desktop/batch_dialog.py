"""Batch fit report — turn a set of already-made fits into a publication table,
per-fit plots and a report in a chosen folder. Front-end for larmor.batch.
"""
from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog,
    QHBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem, QProgressBar,
    QPushButton, QSpinBox, QVBoxLayout,
)

from larmor.desktop import theme


def _peek(path: str) -> tuple[str, str, list[str]]:
    """Cheap metadata read for the list (no spectrum load)."""
    p = Path(path)
    try:
        if p.suffix.lower() == ".json":
            d = json.loads(p.read_text(encoding="utf-8"))
            return (d.get("sample") or p.stem, d.get("nucleus", ""),
                    [s.get("model", "?") for s in d.get("sites", [])])
        if p.suffix.lower() in (".fxmla", ".fxml"):
            from larmor.io import fxmla
            dm = fxmla.read(p)
            dim = dm.dimensions[0]
            return (dm.comment or p.stem, dim.nucleus,
                    [ln.model_name for ln in dim.lines])
    except Exception:
        pass
    return (p.stem, "", [])


class BatchReportDialog(QDialog):
    def __init__(self, parent, start_dir: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Batch fit report — publication table, plots & report")
        self.resize(720, 640)
        self._paths: list[str] = []
        self._stop = False
        self._start_dir = start_dir
        v = QVBoxLayout(self)

        v.addWidget(QLabel(
            "Select several saved fits (LARMOR .recipe.json, dmfit .fxmla, or "
            ".larproj) — ideally the same nucleus and acquisition. LARMOR re-fits "
            "each for fresh errors, quantifies the populations, and writes a "
            "publication table, per-fit plots and a report."))
        v.itemAt(0).widget().setWordWrap(True)

        addrow = QHBoxLayout()
        b_add = QPushButton("Add fits…"); b_add.clicked.connect(self._add)
        b_del = QPushButton("Remove"); b_del.clicked.connect(self._remove)
        addrow.addWidget(b_add); addrow.addWidget(b_del); addrow.addStretch(1)
        v.addLayout(addrow)

        self.list = QListWidget()
        self.list.setSelectionMode(QListWidget.ExtendedSelection)
        v.addWidget(self.list, 1)
        self.warn = QLabel("")
        self.warn.setWordWrap(True)
        self.warn.setStyleSheet(f"color:{theme.active().text_dim};")
        v.addWidget(self.warn)

        # options
        opt = QHBoxLayout()
        opt.addWidget(QLabel("errors"))
        self.method = QComboBox()
        self.method.addItem("Covariance (fast)", "covariance")
        self.method.addItem("Monte-Carlo (bootstrap)", "montecarlo")
        self.method.currentIndexChanged.connect(
            lambda: self.mc_n.setEnabled(self.method.currentData() == "montecarlo"))
        opt.addWidget(self.method)
        opt.addWidget(QLabel("MC trials"))
        self.mc_n = QSpinBox(); self.mc_n.setRange(20, 5000); self.mc_n.setValue(200)
        self.mc_n.setEnabled(False)
        opt.addWidget(self.mc_n)
        opt.addStretch(1)
        v.addLayout(opt)

        fmt = QHBoxLayout()
        fmt.addWidget(QLabel("output:"))
        self.cCsv = QCheckBox("CSV"); self.cCsv.setChecked(True)
        self.cTex = QCheckBox("LaTeX"); self.cTex.setChecked(True)
        self.cMd = QCheckBox("Markdown report"); self.cMd.setChecked(True)
        self.cPlots = QCheckBox("per-fit plots"); self.cPlots.setChecked(True)
        for c in (self.cCsv, self.cTex, self.cMd, self.cPlots):
            fmt.addWidget(c)
        fmt.addStretch(1)
        v.addLayout(fmt)

        outrow = QHBoxLayout()
        outrow.addWidget(QLabel("folder"))
        self.out = QLineEdit()
        outrow.addWidget(self.out, 1)
        b_browse = QPushButton("Browse…"); b_browse.clicked.connect(self._browse)
        outrow.addWidget(b_browse)
        v.addLayout(outrow)

        self.prog = QProgressBar(); v.addWidget(self.prog)
        self.status = QLabel("")
        self.status.setStyleSheet(f"font-weight:600; color:{theme.active().accent};")
        v.addWidget(self.status)

        bb = QDialogButtonBox(QDialogButtonBox.Close | QDialogButtonBox.Help)
        self.btnGen = bb.addButton("Generate", QDialogButtonBox.ApplyRole)
        self.btnGen.clicked.connect(self._generate)
        self.btnStop = bb.addButton("Stop", QDialogButtonBox.ResetRole)
        self.btnStop.setEnabled(False)
        self.btnStop.clicked.connect(lambda: setattr(self, "_stop", True))
        self.btnOpen = bb.addButton("Open folder", QDialogButtonBox.ActionRole)
        self.btnOpen.setEnabled(False)
        self.btnOpen.clicked.connect(self._open_folder)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.helpRequested.connect(self._help)
        v.addWidget(bb)

    # ------------------------------------------------------------------
    def _add(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add fits", self._start_dir,
            "Fits (*.recipe.json *.json *.fxmla *.fxml *.larproj.json);;All (*)")
        for p in paths:
            if p not in self._paths:
                self._paths.append(p)
                sample, nuc, models = _peek(p)
                m = ", ".join(dict.fromkeys(models)) or "?"
                QListWidgetItem(f"{sample}   ·   {nuc or '?'}   ·   {m}", self.list)
        if self._paths and not self.out.text():
            self.out.setText(str(Path(self._paths[0]).parent / "batch_report"))
        self._check()

    def _remove(self):
        for it in self.list.selectedItems():
            r = self.list.row(it)
            self.list.takeItem(r)
            del self._paths[r]
        self._check()

    def _check(self):
        nuclei = {_peek(p)[1] for p in self._paths}
        msgs = []
        if len([n for n in nuclei if n]) > 1:
            msgs.append("⚠ mixed nuclei — the table is only comparable within a "
                        "nucleus")
        self.warn.setText("  ".join(msgs))

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Output folder",
                                             self.out.text() or self._start_dir)
        if d:
            self.out.setText(d)

    def _generate(self):
        from larmor import batch

        if not self._paths:
            self.status.setText("add at least one fit"); return
        if not self.out.text().strip():
            self.status.setText("choose an output folder"); return
        formats = [f for f, c in (("csv", self.cCsv), ("latex", self.cTex),
                                  ("markdown", self.cMd)) if c.isChecked()]
        self._stop = False
        self.btnGen.setEnabled(False); self.btnStop.setEnabled(True)
        self.prog.setValue(0)

        def prog(stage, k, n):
            self.prog.setRange(0, n); self.prog.setValue(k)
            self.status.setText(f"{stage} {k}/{n}…")
            QApplication.processEvents()

        try:
            res = batch.run_batch(
                self._paths, self.out.text().strip(),
                error_method=self.method.currentData(), n_mc=self.mc_n.value(),
                make_plots=self.cPlots.isChecked(), formats=tuple(formats),
                progress=prog, should_stop=lambda: self._stop)
        except Exception as exc:  # noqa: BLE001
            self.status.setText(f"failed: {exc}")
            self.btnGen.setEnabled(True); self.btnStop.setEnabled(False)
            return

        self.btnGen.setEnabled(True); self.btnStop.setEnabled(False)
        self.btnOpen.setEnabled(True)
        self.status.setText(res.summary + ("  ·  stopped early" if self._stop else ""))

    def _open_folder(self):
        import os
        import subprocess
        import sys

        d = self.out.text().strip()
        if not d:
            return
        try:
            if sys.platform.startswith("win"):
                os.startfile(d)  # noqa: S606
            elif sys.platform == "darwin":
                subprocess.Popen(["open", d])
            else:
                subprocess.Popen(["xdg-open", d])
        except Exception:
            pass

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "multi-dataset", "Multi-dataset & co-fitting")
