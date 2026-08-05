"""Plotting studio: a flexible publication-figure builder.

Compose a figure — 1D overlay/stack, 2D contour (contour / density / filled,
projections, iso & quadrupolar reference lines, contour values) or a series —
from any of LARMOR's data (spectra, dmfit fits, saved recipes, fit components),
customise it (title, labels, size, limits, colours, line styles, colour map,
contour levels…), preview it live, and export it at any DPI/size/format.

It is a thin GUI over :mod:`larmor.figures`, which does the actual rendering and
is fully spec-driven — so a figure here is a plain dict you can save, reload, or
build in a script.
"""
from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QColorDialog, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFileDialog, QFormLayout, QHBoxLayout, QInputDialog, QLabel, QLineEdit,
    QListWidget, QListWidgetItem, QMenu, QPushButton, QSpinBox, QVBoxLayout,
    QWidget,
)

from larmor import figures
from larmor.desktop import theme

_LINESTYLES = ["-", "--", "-.", ":"]
_CMAPS = ["viridis", "plasma", "magma", "cividis", "coolwarm", "Blues",
          "Reds", "Greens", "turbo", "gray"]


class _TraceEditor(QDialog):
    """Edit one 1D trace: label, colour, line style, offset and scale."""

    def __init__(self, parent, trace: dict):
        super().__init__(parent)
        self.setWindowTitle("Trace")
        self.trace = dict(trace)
        form = QFormLayout(self)
        self.label = QLineEdit(self.trace.get("label", ""))
        form.addRow("Label", self.label)
        self._color = self.trace.get("color")
        self.btnColor = QPushButton(self._color or "auto")
        self.btnColor.clicked.connect(self._pick_color)
        form.addRow("Colour", self.btnColor)
        self.ls = QComboBox(); self.ls.addItems(_LINESTYLES)
        self.ls.setCurrentText(self.trace.get("linestyle", "-"))
        form.addRow("Line style", self.ls)
        self.off = QDoubleSpinBox(); self.off.setRange(-1e6, 1e6)
        self.off.setValue(float(self.trace.get("offset", 0.0)))
        form.addRow("Offset", self.off)
        self.scale = QDoubleSpinBox(); self.scale.setRange(-1e6, 1e6)
        self.scale.setDecimals(3); self.scale.setValue(float(self.trace.get("scale", 1.0)))
        form.addRow("Scale", self.scale)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        form.addRow(bb)

    def _pick_color(self):
        c = QColorDialog.getColor()
        if c.isValid():
            self._color = c.name()
            self.btnColor.setText(self._color)

    def values(self) -> dict:
        t = dict(self.trace)
        t["label"] = self.label.text() or None
        t["linestyle"] = self.ls.currentText()
        t["offset"] = self.off.value()
        t["scale"] = self.scale.value()
        if self._color:
            t["color"] = self._color
        return t


class PlottingStudio(QDialog):
    def __init__(self, parent=None, spec: dict | None = None):
        super().__init__(parent)
        self.setWindowTitle("Plotting studio")
        self.resize(1080, 720)
        self._traces: list[dict] = []
        self._iso: list[dict] = []

        root = QHBoxLayout(self)

        # ---- controls ----
        ctl = QWidget(); cv = QVBoxLayout(ctl); ctl.setMaximumWidth(390)
        root.addWidget(ctl)

        krow = QHBoxLayout(); krow.addWidget(QLabel("Plot:"))
        self.kind = QComboBox()
        self.kind.addItems(["1D overlay / stack", "2D contour", "Series"])
        self.kind.currentIndexChanged.connect(self._kind_changed)
        krow.addWidget(self.kind, 1)
        cv.addLayout(krow)

        # 1D controls
        self.box1d = QWidget(); b1 = QVBoxLayout(self.box1d)
        b1.setContentsMargins(0, 0, 0, 0)
        self.traceList = QListWidget()
        self.traceList.itemDoubleClicked.connect(self._edit_trace)
        b1.addWidget(QLabel("Traces (double-click to style):"))
        b1.addWidget(self.traceList)
        trow = QHBoxLayout()
        b_add = QPushButton("Add…"); b_add.clicked.connect(self._add_trace_menu)
        b_del = QPushButton("Remove"); b_del.clicked.connect(self._remove_trace)
        trow.addWidget(b_add); trow.addWidget(b_del)
        b1.addLayout(trow)
        srow = QHBoxLayout(); srow.addWidget(QLabel("Stack offset:"))
        self.stack = QDoubleSpinBox(); self.stack.setRange(0, 1e6)
        self.stack.setToolTip("shift each trace up by this × its index (a stacked plot)")
        srow.addWidget(self.stack)
        srow.addWidget(QLabel("normalize:"))
        self.norm = QComboBox(); self.norm.addItems(["none", "max", "area", "noise"])
        self.norm.setToolTip("scale every trace to unit peak / area / noise before "
                             "stacking (for honest series comparison)")
        self.norm.currentIndexChanged.connect(self._refresh)
        srow.addWidget(self.norm)
        b1.addLayout(srow)
        self.chkDiff = QCheckBox("difference vs first trace")
        self.chkDiff.setToolTip("subtract the first trace from the others "
                                "(what changed across the series)")
        self.chkDiff.toggled.connect(self._refresh)
        b1.addWidget(self.chkDiff)
        cv.addWidget(self.box1d)

        # 2D controls
        self.box2d = QWidget(); b2 = QFormLayout(self.box2d)
        self.path2d = QLineEdit()
        p2row = QWidget(); p2l = QHBoxLayout(p2row); p2l.setContentsMargins(0, 0, 0, 0)
        b_pick2 = QPushButton("EXPNO…"); b_pick2.clicked.connect(self._pick_2d)
        p2l.addWidget(self.path2d, 1); p2l.addWidget(b_pick2)
        b2.addRow("2D path", p2row)
        self.cmode = QComboBox()
        self.cmode.addItems(["contour", "density", "filled", "both"])
        b2.addRow("Contour mode", self.cmode)
        self.cmap = QComboBox(); self.cmap.addItems(_CMAPS)
        b2.addRow("Colour map", self.cmap)
        self.nlev = QSpinBox(); self.nlev.setRange(2, 60); self.nlev.setValue(12)
        b2.addRow("Levels", self.nlev)
        self.levmode = QComboBox(); self.levmode.addItems(["log", "linear"])
        b2.addRow("Level spacing", self.levmode)
        self.chkValues = QCheckBox("show contour values")
        b2.addRow(self.chkValues)
        self.chkTop = QCheckBox("top projection"); self.chkTop.setChecked(True)
        self.chkRight = QCheckBox("side projection"); self.chkRight.setChecked(True)
        b2.addRow(self.chkTop); b2.addRow(self.chkRight)
        self.chkNeg = QCheckBox("negative contours")
        b2.addRow(self.chkNeg)
        b_iso = QPushButton("Add iso / quad line…"); b_iso.clicked.connect(self._add_iso)
        self.isoLabel = QLabel("0 lines")
        b2.addRow(b_iso, self.isoLabel)
        cv.addWidget(self.box2d)

        # series controls
        self.boxSer = QWidget(); bs = QFormLayout(self.boxSer)
        self.pathSer = QLineEdit()
        psrow = QWidget(); psl = QHBoxLayout(psrow); psl.setContentsMargins(0, 0, 0, 0)
        b_picks = QPushButton("EXPNO…"); b_picks.clicked.connect(self._pick_series)
        psl.addWidget(self.pathSer, 1); psl.addWidget(b_picks)
        bs.addRow("Series path", psrow)
        self.serMode = QComboBox(); self.serMode.addItems(["satrec", "redor"])
        bs.addRow("Mode", self.serMode)
        self.serStretch = QCheckBox("stretched (β)")
        bs.addRow(self.serStretch)
        cv.addWidget(self.boxSer)

        # common
        common = QFormLayout()
        self.title = QLineEdit(); common.addRow("Title", self.title)
        self.xlabel = QLineEdit(); common.addRow("x label", self.xlabel)
        self.ylabel = QLineEdit(); common.addRow("y label", self.ylabel)
        self.style = QComboBox(); self.style.addItems(list(figures.STYLES))
        common.addRow("Style", self.style)
        xr = QWidget(); xrl = QHBoxLayout(xr); xrl.setContentsMargins(0, 0, 0, 0)
        self.xhi = QDoubleSpinBox(); self.xhi.setRange(-1e6, 1e6)
        self.xlo = QDoubleSpinBox(); self.xlo.setRange(-1e6, 1e6)
        self.chkXlim = QCheckBox("x-range")
        xrl.addWidget(self.chkXlim); xrl.addWidget(self.xhi); xrl.addWidget(self.xlo)
        common.addRow("Limits", xr)
        wr = QWidget(); wrl = QHBoxLayout(wr); wrl.setContentsMargins(0, 0, 0, 0)
        self.wcm = QDoubleSpinBox(); self.wcm.setRange(2, 60); self.wcm.setValue(12)
        self.hcm = QDoubleSpinBox(); self.hcm.setRange(2, 60); self.hcm.setValue(9)
        wrl.addWidget(QLabel("W")); wrl.addWidget(self.wcm)
        wrl.addWidget(QLabel("H")); wrl.addWidget(self.hcm)
        common.addRow("Size (cm)", wr)
        cv.addLayout(common)
        cv.addStretch(1)

        # ---- preview + actions ----
        right = QWidget(); rv = QVBoxLayout(right); root.addWidget(right, 1)
        self.preview = QLabel("preview")
        self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setMinimumSize(560, 480)
        self.preview.setStyleSheet(f"background:{theme.active().plot_bg};")
        rv.addWidget(self.preview, 1)
        arow = QHBoxLayout()
        b_ref = QPushButton("Refresh preview"); b_ref.clicked.connect(self._refresh)
        b_exp = QPushButton("Export…"); b_exp.clicked.connect(self._export)
        b_save = QPushButton("Save spec…"); b_save.clicked.connect(self._save_spec)
        b_load = QPushButton("Load spec…"); b_load.clicked.connect(self._load_spec)
        for b in (b_ref, b_exp, b_save, b_load):
            arow.addWidget(b)
        arow.addStretch(1)
        rv.addLayout(arow)
        self.msg = QLabel(""); self.msg.setWordWrap(True)
        rv.addWidget(self.msg)

        if spec:
            self._apply_spec(spec)
        self._kind_changed()

    # ------------------------------------------------------------------ traces
    def _add_trace_menu(self):
        m = QMenu(self)
        m.addAction("Spectrum (EXPNO / dmfit)…", self._add_spectrum)
        m.addAction("Fit total (recipe)…", lambda: self._add_recipe("total"))
        m.addAction("Fit component…", lambda: self._add_recipe("site"))
        m.addAction("Residual…", lambda: self._add_recipe("residual"))
        m.exec(self.cursor().pos())

    def _add_spectrum(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Spectrum", "", "dmfit (*.fxmla *.fxml);;All (*)")
        if path:
            self._push_trace({"path": path, "label": Path(path).stem})
        else:
            folder = QFileDialog.getExistingDirectory(self, "EXPNO folder")
            if folder:
                self._push_trace({"path": folder, "label": Path(folder).name})

    def _add_recipe(self, part):
        path, _ = QFileDialog.getOpenFileName(
            self, "Recipe", "", "LARMOR recipe (*.json);;All (*)")
        if not path:
            return
        t = {"recipe": path, "part": part, "label": Path(path).stem}
        if part == "site":
            i, ok = QInputDialog.getInt(self, "Component", "Site index (0-based):", 0, 0)
            if not ok:
                return
            t["site"] = i
        self._push_trace(t)

    def _push_trace(self, t):
        self._traces.append(t)
        self.traceList.addItem(QListWidgetItem(t.get("label") or t.get("path", "trace")))
        self._refresh()

    def _remove_trace(self):
        r = self.traceList.currentRow()
        if r >= 0:
            self._traces.pop(r); self.traceList.takeItem(r); self._refresh()

    def _edit_trace(self, item):
        r = self.traceList.row(item)
        dlg = _TraceEditor(self, self._traces[r])
        if dlg.exec() == QDialog.Accepted:
            self._traces[r] = dlg.values()
            item.setText(self._traces[r].get("label") or "trace")
            self._refresh()

    def _add_iso(self):
        slope, ok = QInputDialog.getDouble(self, "Reference line",
                                           "Slope (F1 per F2):", 1.0, -20, 20, 4)
        if not ok:
            return
        inter, ok = QInputDialog.getDouble(self, "Reference line",
                                           "Intercept (ppm):", 0.0, -1e5, 1e5, 3)
        if not ok:
            return
        label, _ = QInputDialog.getText(self, "Reference line", "Label:")
        self._iso.append({"slope": slope, "intercept": inter,
                          "label": label or None, "linestyle": "--"})
        self.isoLabel.setText(f"{len(self._iso)} line(s)")
        self._refresh()

    # ------------------------------------------------------------------ pickers
    def _pick_2d(self):
        folder = QFileDialog.getExistingDirectory(self, "2D EXPNO folder")
        if folder:
            self.path2d.setText(folder); self._refresh()

    def _pick_series(self):
        folder = QFileDialog.getExistingDirectory(self, "Series EXPNO folder")
        if folder:
            self.pathSer.setText(folder); self._refresh()

    def _kind_changed(self):
        k = self.kind.currentIndex()
        self.box1d.setVisible(k == 0)
        self.box2d.setVisible(k == 1)
        self.boxSer.setVisible(k == 2)
        self._refresh()

    # ------------------------------------------------------------------ spec
    def _spec(self) -> dict:
        common = {"style": self.style.currentText(),
                  "figsize": (self.wcm.value() / 2.54, self.hcm.value() / 2.54)}
        if self.title.text():
            common["title"] = self.title.text()
        if self.xlabel.text():
            common["xlabel"] = self.xlabel.text()
        if self.ylabel.text():
            common["ylabel"] = self.ylabel.text()
        if self.chkXlim.isChecked():
            common["xlim"] = (self.xhi.value(), self.xlo.value())
        k = self.kind.currentIndex()
        if k == 0:
            traces = []
            for i, t in enumerate(self._traces):
                tt = dict(t)
                if self.stack.value():
                    tt["offset"] = float(tt.get("offset", 0.0)) + i * self.stack.value()
                traces.append(tt)
            spec = {"kind": "1d", "traces": traces, **common}
            if self.norm.currentText() != "none":
                spec["norm"] = self.norm.currentText()
            if self.chkDiff.isChecked():
                spec["difference"] = True
            return spec
        if k == 1:
            return {"kind": "2d", "path": self.path2d.text(),
                    "contour_mode": self.cmode.currentText(),
                    "cmap": self.cmap.currentText(),
                    "levels": {"n": self.nlev.value(),
                               "mode": self.levmode.currentText()},
                    "contour_values": self.chkValues.isChecked(),
                    "proj_top": self.chkTop.isChecked(),
                    "proj_right": self.chkRight.isChecked(),
                    "negative": self.chkNeg.isChecked(),
                    "iso_lines": self._iso, **common}
        return {"kind": "series", "path": self.pathSer.text(),
                "mode": self.serMode.currentText(),
                "stretched": self.serStretch.isChecked(), **common}

    def _apply_spec(self, spec: dict):
        self.style.setCurrentText(spec.get("style", "article"))
        self.title.setText(spec.get("title", ""))
        self.xlabel.setText(spec.get("xlabel", ""))
        self.ylabel.setText(spec.get("ylabel", ""))
        kind = {"1d": 0, "2d": 1, "series": 2}.get(spec.get("kind", "1d"), 0)
        self.kind.setCurrentIndex(kind)
        if kind == 0:
            for t in spec.get("traces", []):
                self._traces.append(t)
                self.traceList.addItem(t.get("label") or "trace")
        elif kind == 1:
            self.path2d.setText(spec.get("path", ""))
            self.cmode.setCurrentText(spec.get("contour_mode", "contour"))
            self.cmap.setCurrentText(spec.get("cmap", "viridis"))
            self._iso = list(spec.get("iso_lines", []))
            self.isoLabel.setText(f"{len(self._iso)} line(s)")
        elif kind == 2:
            self.pathSer.setText(spec.get("path", ""))
            self.serMode.setCurrentText(spec.get("mode", "satrec"))

    # ------------------------------------------------------------------ render
    def _refresh(self):
        try:
            spec = self._spec()
            if spec["kind"] == "1d" and not spec["traces"]:
                self.preview.setText("add a trace to preview"); return
            if spec["kind"] in ("2d", "series") and not spec.get("path"):
                self.preview.setText("pick an EXPNO folder to preview"); return
            png = figures.render_png_bytes(spec, dpi=110)
            pix = QPixmap(); pix.loadFromData(png)
            self.preview.setPixmap(pix.scaled(
                self.preview.width(), self.preview.height(),
                Qt.KeepAspectRatio, Qt.SmoothTransformation))
            self.msg.setText("")
        except Exception as exc:  # noqa: BLE001
            self.msg.setText(f"cannot render: {exc}")

    def _export(self):
        from larmor.desktop.export_dialog import export_matplotlib
        try:
            fig = figures.render(self._spec())
        except Exception as exc:  # noqa: BLE001
            self.msg.setText(f"cannot render: {exc}"); return
        path = export_matplotlib(self, fig, "figure")
        if path:
            self.msg.setText(f"exported {path}")

    def _save_spec(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save figure spec",
                                              "figure.json", "JSON (*.json)")
        if path:
            Path(path).write_text(json.dumps(self._spec(), indent=2))
            self.msg.setText(f"saved {path}")

    def _load_spec(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load figure spec", "",
                                              "JSON (*.json)")
        if path:
            self._apply_spec(json.loads(Path(path).read_text()))
            self._kind_changed()
