"""Plotting studio: a flexible publication-figure builder.

Compose a figure — 1D overlay/stack, 2D contour (contour / density / filled,
projections, iso & quadrupolar reference lines, contour values) or a series —
from any of LARMOR's data (spectra, dmfit fits, saved recipes, fit components),
customise it fully (title, labels, colours, line styles, colour map, contour
levels, x/y limits, major-tick spacing, minor ticks, tick direction, grid,
legend position & columns, font and line-width, figure size) and — for series /
parameter traces that carry them — **error bars**, preview it live, and export
it at any DPI/size/format.

It is a thin GUI over :mod:`larmor.figures`, which does the actual rendering and
is fully spec-driven — so a figure here is a plain dict you can save, reload, or
build in a script.
"""
from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import Qt, QSettings
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QColorDialog, QComboBox, QDialog,
    QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout,
    QInputDialog, QLabel, QLineEdit, QListWidget, QListWidgetItem, QMenu,
    QMessageBox, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
    QToolButton, QVBoxLayout, QWidget,
)

from larmor import figures, series_grid

_KINDS = ["1D overlay / stack", "2D contour", "Series", "Batch grid",
         "Species distribution"]

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

        has_yerr = bool((self.trace.get("data") or {}).get("yerr"))
        note = QLabel("(no error data on this trace)" if not has_yerr else "")
        note.setStyleSheet("color:#888; font-size:10px;")
        form.addRow("Error bars", note if not has_yerr else self._error_row())
        if has_yerr:
            self.errWidth = QDoubleSpinBox(); self.errWidth.setRange(0.2, 6.0)
            self.errWidth.setDecimals(2)
            self.errWidth.setValue(float(self.trace.get("err_width", 1.2)))
            form.addRow("  width", self.errWidth)
            self.errCap = QDoubleSpinBox(); self.errCap.setRange(0.0, 12.0)
            self.errCap.setDecimals(1)
            self.errCap.setValue(float(self.trace.get("err_capsize", 3.5)))
            form.addRow("  cap size", self.errCap)
            self._err_color = self.trace.get("err_color")
            self.btnErrColor = QPushButton(self._err_color or "match line")
            self.btnErrColor.clicked.connect(self._pick_err_color)
            form.addRow("  colour", self.btnErrColor)
        else:
            self.errVisible = None

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        form.addRow(bb)

    def _error_row(self):
        self.errVisible = QCheckBox("show")
        self.errVisible.setChecked(self.trace.get("err_visible", True))
        return self.errVisible

    def _pick_color(self):
        c = QColorDialog.getColor()
        if c.isValid():
            self._color = c.name()
            self.btnColor.setText(self._color)

    def _pick_err_color(self):
        c = QColorDialog.getColor()
        if c.isValid():
            self._err_color = c.name()
            self.btnErrColor.setText(self._err_color)

    def values(self) -> dict:
        t = dict(self.trace)
        t["label"] = self.label.text() or None
        t["linestyle"] = self.ls.currentText()
        t["offset"] = self.off.value()
        t["scale"] = self.scale.value()
        if self.errVisible is not None:
            t["err_visible"] = self.errVisible.isChecked()
            t["err_width"] = self.errWidth.value()
            t["err_capsize"] = self.errCap.value()
            if self._err_color:
                t["err_color"] = self._err_color
        if self._color:
            t["color"] = self._color
        return t


class _ReferenceLineDialog(QDialog):
    """A 2D reference line for an MQMAS/correlation figure: either typed by
    hand, or "Compute"d from larmor.twod's own physics (f1_cs_scale for the
    chemical-shift/diagonal axis, qis_slope for the quadrupolar-induced-shift
    axis) for a given nucleus/Larmor frequency/method — always still
    editable afterward, never silently overriding a hand-tuned line."""

    def __init__(self, parent, nucleus: str = "", larmor_MHz: float = 0.0):
        super().__init__(parent)
        self.setWindowTitle("Reference line")
        form = QFormLayout(self)
        self.kind = QComboBox()
        self.kind.addItems(["Manual", "CS axis (computed)", "QIS axis (computed)"])
        self.kind.currentIndexChanged.connect(self._kind_changed)
        form.addRow("Kind", self.kind)
        self.nucleus = QLineEdit(nucleus)
        form.addRow("Nucleus", self.nucleus)
        self.larmor = QDoubleSpinBox(); self.larmor.setRange(0, 2000)
        self.larmor.setDecimals(3); self.larmor.setValue(larmor_MHz)
        form.addRow("Larmor (MHz)", self.larmor)
        self.method = QComboBox(); self.method.addItems(["3QMAS", "5QMAS", "ST1"])
        form.addRow("Method", self.method)
        self.anchor = QDoubleSpinBox(); self.anchor.setRange(-1e4, 1e4)
        self.anchor.setDecimals(2)
        self.anchor.setToolTip(
            "the site's own isotropic shift (ppm) — the QIS axis is drawn "
            "starting from this point on the CS axis, moving away as Cq grows")
        form.addRow("Anchor δiso (ppm, QIS only)", self.anchor)
        b_compute = QPushButton("Compute")
        b_compute.setToolTip(
            "runs a reference mrsimulator simulation for this nucleus/method "
            "(cached) — a few seconds the first time")
        b_compute.clicked.connect(self._compute)
        form.addRow(b_compute)
        self.slope = QDoubleSpinBox(); self.slope.setRange(-20, 20)
        self.slope.setDecimals(4); self.slope.setValue(1.0)
        form.addRow("Slope (F1 per F2)", self.slope)
        self.intercept = QDoubleSpinBox(); self.intercept.setRange(-1e5, 1e5)
        self.intercept.setDecimals(3)
        form.addRow("Intercept (ppm)", self.intercept)
        self.labelEdit = QLineEdit()
        form.addRow("Label", self.labelEdit)
        self._kind_changed()
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        form.addRow(bb)

    def _kind_changed(self):
        manual = self.kind.currentIndex() == 0
        for w in (self.nucleus, self.larmor, self.method):
            w.setEnabled(not manual)
        self.anchor.setEnabled(self.kind.currentIndex() == 2)   # QIS only

    def _compute(self):
        nucleus = self.nucleus.text().strip()
        larmor = self.larmor.value()
        if self.kind.currentIndex() == 0:
            return
        if not nucleus or not larmor:
            QMessageBox.warning(self, "Reference line",
                                "enter a nucleus and Larmor frequency first")
            return
        method = self.method.currentText()
        try:
            from larmor.twod import f1_cs_scale, qis_slope
            cs = f1_cs_scale(nucleus, larmor, method)
            if self.kind.currentIndex() == 1:                    # CS axis
                self.slope.setValue(cs)
                self.intercept.setValue(0.0)
                if not self.labelEdit.text():
                    self.labelEdit.setText("CS axis")
            else:                                                 # QIS axis
                raw_slope = qis_slope(nucleus, larmor, method) * cs
                anchor = self.anchor.value()
                self.slope.setValue(raw_slope)
                self.intercept.setValue(cs * anchor - raw_slope * anchor)
                if not self.labelEdit.text():
                    self.labelEdit.setText("QIS axis")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Reference line", f"couldn't compute: {exc}")

    def values(self) -> dict:
        return {"slope": self.slope.value(), "intercept": self.intercept.value(),
                "label": self.labelEdit.text() or None, "linestyle": "--"}


class PlottingStudio(QDialog):
    def __init__(self, parent=None, spec: dict | None = None):
        super().__init__(parent)
        self.setWindowTitle("Plotting studio")
        self.resize(1240, 720)
        self._traces: list[dict] = []
        self._iso: list[dict] = []
        self._panels: list[dict] = []       # batch-grid: [{"path","sample","include"}]
        self._trace_defaults: dict = {}     # from the selected template, applied to new traces
        self._component_colors: dict[int, str] = {}   # batch-grid: site index -> "#hex"
        self._legend_hide: set[int] = set()            # batch-grid: site indices, legend only

        from PySide6.QtWidgets import QSplitter
        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ---- built-in file explorer (so you needn't leave LARMOR to find files) ----
        from larmor.desktop.explorer import ExplorerPanel
        self.files = ExplorerPanel()
        self.files.btnBatch.setVisible(False)
        self.files.open_requested.connect(self._add_from_explorer)
        self.files.setMinimumWidth(220)
        splitter.addWidget(self.files)

        # ---- controls ----
        ctl = QWidget(); cv = QVBoxLayout(ctl); ctl.setMaximumWidth(390)
        splitter.addWidget(ctl)

        trow0 = QHBoxLayout(); trow0.addWidget(QLabel("Template:"))
        self.template = QComboBox()
        self.template.addItem("(none — build from scratch)", None)
        for name in figures.TEMPLATES:
            self.template.addItem(name, name)
        self.template.setToolTip("A named, generic starting point (layout + "
                                 "sensible defaults) for a common NMR figure "
                                 "type — pick one, then customise anything "
                                 "below. Combines freely with any Style.")
        self.template.currentIndexChanged.connect(self._template_changed)
        trow0.addWidget(self.template, 1)
        cv.addLayout(trow0)
        self.templateDesc = QLabel("")
        self.templateDesc.setWordWrap(True)
        self.templateDesc.setStyleSheet("color:#888; font-size:10px;")
        cv.addWidget(self.templateDesc)

        krow = QHBoxLayout(); krow.addWidget(QLabel("Plot:"))
        self.kind = QComboBox()
        self.kind.addItems(_KINDS)
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
        self.nuc2d = QLineEdit()
        self.nuc2d.setPlaceholderText("e.g. 27Al — for axis labels + computed lines")
        b2.addRow("Nucleus", self.nuc2d)
        self.larmor2d = QDoubleSpinBox()
        self.larmor2d.setRange(0, 2000); self.larmor2d.setDecimals(3)
        self.larmor2d.setToolTip("¹H Larmor frequency (MHz) — needed for "
                                 "computed reference lines (Add iso/quad line…)")
        b2.addRow("Larmor (MHz)", self.larmor2d)
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
        self.fitRecipe2d = QLineEdit()
        fr2row = QWidget(); fr2l = QHBoxLayout(fr2row); fr2l.setContentsMargins(0, 0, 0, 0)
        b_pickfit2 = QPushButton("Recipe…"); b_pickfit2.clicked.connect(self._pick_fit_2d)
        fr2l.addWidget(self.fitRecipe2d, 1); fr2l.addWidget(b_pickfit2)
        b2.addRow("Fit overlay", fr2row)
        self.fitRecipe2d.setToolTip(
            "a saved 2D fit (.recipe.json from Decomposition ▸ Fit on a 2D "
            "map) — overlays it as a dashed contour on the experimental map, "
            "for a publication figure that shows both")
        self.mqmasMethod = QComboBox()
        self.mqmasMethod.addItems(["3QMAS", "5QMAS", "ST1"])
        self.mqmasMethod.setToolTip(
            "the MQMAS method the fit overlay was actually fit with — not "
            "stored on the recipe itself, so specify it here")
        b2.addRow("MQMAS method", self.mqmasMethod)
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

        # batch grid: a publication small-multiples figure from already-fitted
        # spectra — load a batch_table*.csv (auto-matched to its saved fits) or
        # a folder of .recipe.json fits, pick/reorder which panels to show
        self.boxGrid = QWidget(); bg = QVBoxLayout(self.boxGrid)
        bg.setContentsMargins(0, 0, 0, 0)
        grow = QHBoxLayout()
        b_csv = QPushButton("Load CSV…")
        b_csv.setToolTip("a batch_table*.csv from the batch-fit dialog — "
                         "auto-matches its rows to sibling .recipe.json fits "
                         "saved from the same session (Save individual fits…)")
        b_csv.clicked.connect(self._grid_load_csv)
        b_folder = QPushButton("Load folder…")
        b_folder.setToolTip("a folder of saved .recipe.json / .fxmla fits")
        b_folder.clicked.connect(self._grid_load_folder)
        b_roots = QPushButton("Data folders…")
        b_roots.setToolTip("default folder(s) the 'locate data for…' popups "
                           "start from — set this once to your raw-data "
                           "root(s) so relocating an unresolved sample is "
                           "usually one click, not a full browse")
        b_roots.clicked.connect(self._grid_set_data_roots)
        grow.addWidget(b_csv); grow.addWidget(b_folder); grow.addWidget(b_roots)
        bg.addLayout(grow)
        self.gridList = QListWidget()
        self.gridList.setToolTip("check the spectra to include; drag a row "
                                 "to reorder, or select one and use the ↑↓ "
                                 "buttons. Double-click a row needing data "
                                 "to locate it, or a resolved one to rename "
                                 "its panel title")
        self.gridList.setDragDropMode(QAbstractItemView.InternalMove)
        self.gridList.setDefaultDropAction(Qt.MoveAction)
        self.gridList.itemChanged.connect(self._grid_item_changed)
        self.gridList.itemDoubleClicked.connect(self._grid_item_double_clicked)
        self.gridList.model().rowsMoved.connect(self._grid_rows_reordered)
        bg.addWidget(QLabel("Panels (checked = included; drag to reorder):"))
        bg.addWidget(self.gridList)
        gmrow = QHBoxLayout()
        self.gridBtnUp = QToolButton(); self.gridBtnUp.setArrowType(Qt.UpArrow)
        self.gridBtnUp.setToolTip("move selected panel up")
        self.gridBtnUp.clicked.connect(lambda: self._grid_move(-1))
        self.gridBtnDown = QToolButton(); self.gridBtnDown.setArrowType(Qt.DownArrow)
        self.gridBtnDown.setToolTip("move selected panel down")
        self.gridBtnDown.clicked.connect(lambda: self._grid_move(1))
        b_rm = QPushButton("Remove")
        b_rm.clicked.connect(self._grid_remove)
        gmrow.addWidget(self.gridBtnUp); gmrow.addWidget(self.gridBtnDown)
        gmrow.addWidget(b_rm)
        gmrow.addStretch(1)
        bg.addLayout(gmrow)
        gf = QFormLayout()
        self.gridCols = QSpinBox(); self.gridCols.setRange(0, 12)
        self.gridCols.setToolTip("grid columns (0 = automatic)")
        gf.addRow("Columns", self.gridCols)
        self.gridComp = QComboBox()
        self.gridComp.addItems(["fill", "dashed", "hidden"])
        self.gridComp.setToolTip("how each site's component curve is drawn — "
                                 "fill = classic deconvolution look, dashed = "
                                 "outline only, hidden = experiment + total fit only")
        gf.addRow("Components", self.gridComp)
        self.gridShade = QLineEdit()
        self.gridShade.setPlaceholderText("empty = show all")
        self.gridShade.setToolTip("site indices to draw (0-based, comma-"
                                  "separated) — leave empty to show every "
                                  "site, or list e.g. '1' to highlight just "
                                  "one component per panel (a composition-"
                                  "series style)")
        gf.addRow("Shade only", self.gridShade)
        self.gridHide = QLineEdit()
        self.gridHide.setPlaceholderText("empty = hide none")
        self.gridHide.setToolTip("site indices to drop ENTIRELY (0-based, "
                                 "comma-separated) — no line, fill, or "
                                 "legend entry in any panel, unlike 'Shade "
                                 "only' this doesn't affect other sites")
        gf.addRow("Hide", self.gridHide)
        self.gridLabels = QComboBox()
        self.gridLabels.addItems(["none", "position", "label", "position+pct"])
        self.gridLabels.setToolTip("label each shown component with its "
                                   "position, its letter, or position + "
                                   "integrated population %")
        gf.addRow("Peak labels", self.gridLabels)
        self.gridTotal = QCheckBox("total fit"); self.gridTotal.setChecked(True)
        self.gridExp = QCheckBox("experiment"); self.gridExp.setChecked(True)
        gtrow = QHBoxLayout(); gtrow.addWidget(self.gridTotal); gtrow.addWidget(self.gridExp)
        gf.addRow("Show", gtrow)
        b_comps = QPushButton("Component colors / legend…")
        b_comps.setToolTip("per-component color override and legend "
                           "visibility (detected from the first resolved "
                           "panel's fit)")
        b_comps.clicked.connect(self._grid_edit_components)
        gf.addRow(b_comps)
        bg.addLayout(gf)
        self.gridMsg = QLabel(""); self.gridMsg.setWordWrap(True)
        self.gridMsg.setStyleSheet("color:#888; font-size:10px;")
        bg.addWidget(self.gridMsg)
        cv.addWidget(self.boxGrid)

        # species distribution: a 100%-stacked bar from a small category x
        # species table (typed in, or pivoted from a batch CSV's own values)
        self.boxBar = QWidget(); bb_ = QVBoxLayout(self.boxBar)
        bb_.setContentsMargins(0, 0, 0, 0)
        b_barcsv = QPushButton("Load from batch CSV…")
        b_barcsv.setToolTip("pivot a chosen parameter (e.g. amplitude) from a "
                            "batch_table*.csv into this table, one row per "
                            "sample — each row is normalized to 100% "
                            "automatically, so raw amplitudes work directly")
        b_barcsv.clicked.connect(self._bar_load_csv)
        bb_.addWidget(b_barcsv)
        self.barTable = QTableWidget(1, 1)
        self.barTable.setHorizontalHeaderLabels(["category"])
        self.barTable.setToolTip("one row per composition point; column 0 is "
                                 "its category label, each further column is "
                                 "one species — edit the column header by "
                                 "double-clicking it")
        self.barTable.horizontalHeader().sectionDoubleClicked.connect(
            self._bar_rename_column)
        bb_.addWidget(self.barTable)
        btrow = QHBoxLayout()
        b_barrow = QPushButton("+ row"); b_barrow.clicked.connect(self._bar_add_row)
        b_barcol = QPushButton("+ species"); b_barcol.clicked.connect(self._bar_add_col)
        b_barrmrow = QPushButton("− row"); b_barrmrow.clicked.connect(self._bar_remove_row)
        b_barrmcol = QPushButton("− species"); b_barrmcol.clicked.connect(self._bar_remove_col)
        for b in (b_barrow, b_barcol, b_barrmrow, b_barrmcol):
            btrow.addWidget(b)
        bb_.addLayout(btrow)
        self.barValueLabels = QCheckBox("value labels on bars")
        self.barValueLabels.setChecked(True)
        bb_.addWidget(self.barValueLabels)
        cv.addWidget(self.boxBar)

        # common
        common = QFormLayout()
        self.title = QLineEdit(); common.addRow("Title", self.title)
        self.xlabel = QLineEdit(); common.addRow("x label", self.xlabel)
        self.ylabel = QLineEdit(); common.addRow("y label", self.ylabel)
        self.style = QComboBox(); self.style.addItems(list(figures.STYLES))
        common.addRow("Style", self.style)
        self.chkPpm = QCheckBox("NMR ppm axis (high → low)")
        self.chkPpm.setChecked(True)
        self.chkPpm.setToolTip("on: chemical-shift axis, inverted, intensity axis "
                               "hidden. off: an ordinary x–y plot (e.g. a "
                               "parameter-vs-sample series)")
        common.addRow("x-axis", self.chkPpm)
        self.xticks = QLineEdit()
        self.xticks.setPlaceholderText("custom tick labels, comma-separated "
                                       "(e.g. 0Ca, 1Ca, 2Ca) or pos:label")
        self.xticks.setToolTip("leave empty for automatic numeric ticks; or list "
                               "labels placed at 1,2,3… or explicit pos:label pairs")
        common.addRow("x-ticks", self.xticks)
        xr = QWidget(); xrl = QHBoxLayout(xr); xrl.setContentsMargins(0, 0, 0, 0)
        self.xhi = QDoubleSpinBox(); self.xhi.setRange(-1e9, 1e9)
        self.xlo = QDoubleSpinBox(); self.xlo.setRange(-1e9, 1e9)
        self.chkXlim = QCheckBox("set")
        self.xhi.setToolTip("x max (left on a ppm axis)")
        self.xlo.setToolTip("x min (right on a ppm axis)")
        xrl.addWidget(self.chkXlim); xrl.addWidget(self.xhi); xrl.addWidget(self.xlo)
        common.addRow("x-limits", xr)
        yr = QWidget(); yrl = QHBoxLayout(yr); yrl.setContentsMargins(0, 0, 0, 0)
        self.ylo = QDoubleSpinBox(); self.ylo.setRange(-1e9, 1e9)
        self.yhi = QDoubleSpinBox(); self.yhi.setRange(-1e9, 1e9)
        self.chkYlim = QCheckBox("set")
        self.ylo.setToolTip("y min (bottom)"); self.yhi.setToolTip("y max (top)")
        yrl.addWidget(self.chkYlim); yrl.addWidget(self.ylo); yrl.addWidget(self.yhi)
        common.addRow("y-limits", yr)

        tr = QWidget(); trl = QHBoxLayout(tr); trl.setContentsMargins(0, 0, 0, 0)
        self.xstep = QDoubleSpinBox(); self.xstep.setRange(0, 1e6); self.xstep.setDecimals(3)
        self.xstep.setToolTip("major x-tick spacing (0 = automatic)")
        self.ystep = QDoubleSpinBox(); self.ystep.setRange(0, 1e6); self.ystep.setDecimals(3)
        self.ystep.setToolTip("major y-tick spacing (0 = automatic)")
        self.chkMinor = QCheckBox("minor")
        self.chkMinor.setToolTip("show minor ticks")
        trl.addWidget(QLabel("Δx")); trl.addWidget(self.xstep)
        trl.addWidget(QLabel("Δy")); trl.addWidget(self.ystep)
        trl.addWidget(self.chkMinor)
        common.addRow("Tick step", tr)
        dr = QWidget(); drl = QHBoxLayout(dr); drl.setContentsMargins(0, 0, 0, 0)
        self.tickdir = QComboBox(); self.tickdir.addItems(["in", "out", "inout"])
        self.tickdir.setToolTip("tick mark direction")
        self.chkGrid = QCheckBox("grid")
        drl.addWidget(self.tickdir, 1); drl.addWidget(self.chkGrid)
        common.addRow("Ticks / grid", dr)

        lr = QWidget(); lrl = QHBoxLayout(lr); lrl.setContentsMargins(0, 0, 0, 0)
        self.legloc = QComboBox()
        self.legloc.addItems(["best", "upper right", "upper left", "lower left",
                              "lower right", "center", "center left",
                              "center right", "upper center", "lower center",
                              "none"])
        self.legloc.setToolTip("legend position ('none' hides it)")
        self.legncol = QSpinBox(); self.legncol.setRange(1, 8); self.legncol.setValue(1)
        self.legncol.setToolTip("legend columns")
        lrl.addWidget(self.legloc, 1); lrl.addWidget(QLabel("cols")); lrl.addWidget(self.legncol)
        common.addRow("Legend", lr)

        fr = QWidget(); frl = QHBoxLayout(fr); frl.setContentsMargins(0, 0, 0, 0)
        self.fontsz = QDoubleSpinBox(); self.fontsz.setRange(0, 48); self.fontsz.setValue(0)
        self.fontsz.setToolTip("base font size (0 = the style's default)")
        self.lwd = QDoubleSpinBox(); self.lwd.setRange(0, 10); self.lwd.setDecimals(2); self.lwd.setValue(0)
        self.lwd.setToolTip("line width (0 = the style's default)")
        frl.addWidget(QLabel("font")); frl.addWidget(self.fontsz)
        frl.addWidget(QLabel("line")); frl.addWidget(self.lwd)
        common.addRow("Font / line", fr)

        wr = QWidget(); wrl = QHBoxLayout(wr); wrl.setContentsMargins(0, 0, 0, 0)
        self.wcm = QDoubleSpinBox(); self.wcm.setRange(2, 60); self.wcm.setValue(12)
        self.hcm = QDoubleSpinBox(); self.hcm.setRange(2, 60); self.hcm.setValue(9)
        wrl.addWidget(QLabel("W")); wrl.addWidget(self.wcm)
        wrl.addWidget(QLabel("H")); wrl.addWidget(self.hcm)
        common.addRow("Size (cm)", wr)
        cv.addLayout(common)
        cv.addStretch(1)

        # ---- preview (crisp matplotlib canvas) + actions ----
        right = QWidget(); rv = QVBoxLayout(right); splitter.addWidget(right)
        splitter.setSizes([230, 380, 620])
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        self._canvas = FigureCanvasQTAgg(Figure(figsize=(5, 4)))
        self._canvas.setMinimumSize(560, 460)
        rv.addWidget(self._canvas, 1)
        prow = QHBoxLayout()
        self.chkAuto = QCheckBox("Auto update")
        self.chkAuto.setToolTip(
            "re-render on every control change — off by default since a "
            "batch grid with many panels (each a full fit reconstruction + "
            "population-% integral) can be slow to redo on every tweak; use "
            "Preview to render on demand, or turn this on for cheap plots")
        self.chkAuto.toggled.connect(self._auto_toggled)
        b_preview = QPushButton("Preview")
        b_preview.setToolTip("render now, regardless of Auto update")
        b_preview.clicked.connect(self._refresh)
        prow.addWidget(self.chkAuto)
        prow.addWidget(b_preview)
        prow.addStretch(1)
        rv.addLayout(prow)

        arow = QHBoxLayout()
        b_exp = QPushButton("Export…"); b_exp.clicked.connect(self._export)
        b_save = QPushButton("Save spec…"); b_save.clicked.connect(self._save_spec)
        b_load = QPushButton("Load spec…"); b_load.clicked.connect(self._load_spec)
        for b in (b_exp, b_save, b_load):
            arow.addWidget(b)
        arow.addStretch(1)
        rv.addLayout(arow)
        self.msg = QLabel(""); self.msg.setWordWrap(True)
        rv.addWidget(self.msg)

        # auto-update: OFF by default (see chkAuto above) -- when on, a short
        # debounce follows every control change without re-rendering on each
        # keystroke; when off, _schedule() is a no-op and only Preview (or an
        # explicit action like loading a CSV) renders
        from PySide6.QtCore import QTimer
        self._debounce = QTimer(self); self._debounce.setSingleShot(True)
        self._debounce.setInterval(180)
        self._debounce.timeout.connect(self._refresh)
        for w in (self.title, self.xlabel, self.ylabel, self.xticks):
            w.textChanged.connect(self._schedule)
        for w in (self.style, self.norm, self.cmap, self.cmode, self.levmode,
                  self.tickdir, self.legloc):
            w.currentIndexChanged.connect(self._schedule)
        for w in (self.wcm, self.hcm, self.xhi, self.xlo, self.yhi, self.ylo,
                  self.nlev, self.stack, self.xstep, self.ystep, self.fontsz,
                  self.lwd, self.legncol):
            w.valueChanged.connect(self._schedule)
        for w in (self.chkPpm, self.chkXlim, self.chkYlim, self.chkValues,
                  self.chkTop, self.chkRight, self.chkNeg, self.chkMinor,
                  self.chkGrid):
            w.toggled.connect(self._schedule)
        for w in (self.gridComp, self.gridLabels):
            w.currentIndexChanged.connect(self._schedule)
        self.gridCols.valueChanged.connect(self._schedule)
        self.gridShade.textChanged.connect(self._schedule)
        self.gridHide.textChanged.connect(self._schedule)
        self.nuc2d.textChanged.connect(self._schedule)
        self.larmor2d.valueChanged.connect(self._schedule)
        self.fitRecipe2d.textChanged.connect(self._schedule)
        self.mqmasMethod.currentIndexChanged.connect(self._schedule)
        for w in (self.gridTotal, self.gridExp):
            w.toggled.connect(self._schedule)
        self.barTable.itemChanged.connect(self._schedule)
        self.barValueLabels.toggled.connect(self._schedule)

        if spec:
            self._apply_spec(spec)
        self._kind_changed()

    def _schedule(self, *a):
        if self.chkAuto.isChecked():
            self._debounce.start()

    def _auto_toggled(self, on: bool):
        if on:
            self._refresh()   # don't leave a stale preview when turning it on

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

    def _add_from_explorer(self, path: str):
        """A file double-clicked in the built-in explorer → add it sensibly for the
        current plot kind (a 1D trace, a 2D path, or a series path)."""
        from pathlib import Path as _P
        low = path.lower()
        k = self.kind.currentIndex()
        if k == 1:                                  # 2D contour
            self.path2d.setText(path); self._refresh(); return
        if k == 2:                                  # series
            self.pathSer.setText(path); self._refresh(); return
        # name the trace by the SAMPLE (from the recipe, else the sample folder),
        # not the bare file name ("1r")
        from larmor.desktop.batchfit_dialog import sample_label
        label = sample_label(path, {})
        if label.lower().endswith((".fxml", ".fxmla", ".json")):
            label = _P(label).stem                  # a fit file → its own name
        if low.endswith((".recipe.json", ".json")):
            self._push_trace({"recipe": path, "part": "total", "label": label})
        else:                                        # spectrum / dmfit fit / EXPNO
            self._push_trace({"path": path, "label": label})

    def _push_trace(self, t):
        # a selected template's trace_defaults (e.g. end_label: True for a
        # stacked-series style) fill in anything this trace doesn't already set
        merged = {**self._trace_defaults, **t}
        self._traces.append(merged)
        self.traceList.addItem(QListWidgetItem(merged.get("label") or merged.get("path", "trace")))
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
        dlg = _ReferenceLineDialog(self, self.nuc2d.text(), self.larmor2d.value())
        if dlg.exec() != QDialog.Accepted:
            return
        self._iso.append(dlg.values())
        self.isoLabel.setText(f"{len(self._iso)} line(s)")
        self._refresh()

    # ------------------------------------------------------------------ batch grid
    def _grid_load_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Batch table CSV", "", "CSV (*.csv)")
        if path:
            self._grid_load(path)

    def _grid_load_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Folder of saved fits")
        if folder:
            self._grid_load(folder)

    def _grid_load(self, source):
        panels, warnings = series_grid.load_panels(source)
        if not panels:
            self.gridMsg.setText(" ".join(warnings) or "no spectra found there")
            return
        self._panels = [{
            "path": p.path, "sample": p.sample, "nucleus": p.nucleus,
            "models": list(p.models), "has_data": p.has_data,
            "data_path": p.data_path, "needs_manual": p.needs_manual,
            "reconstructed": p.reconstructed, "title": None, "include": True,
        } for p in panels]
        self._grid_resolve_manual()
        self._grid_refresh_list()
        self.gridMsg.setText(" ".join(warnings))
        self._refresh()

    @staticmethod
    def _grid_data_roots() -> list[str]:
        raw = QSettings("LARMOR", "app").value("plottingStudio/dataRoots", "")
        return [r for r in (raw or "").split(";") if r]

    def _grid_set_data_roots(self):
        current = "; ".join(self._grid_data_roots())
        text, ok = QInputDialog.getText(
            self, "Data folders",
            "Default folder(s) for 'locate data for…' — separate multiple "
            "with ';':", text=current)
        if ok:
            QSettings("LARMOR", "app").setValue(
                "plottingStudio/dataRoots",
                "; ".join(p.strip() for p in text.split(";") if p.strip()))

    def _grid_ask_manual_path(self, sample: str) -> str:
        """The "successive popups" fallback: series_grid couldn't pair this
        sample with a saved fit or a source_path hint (an older CSV, moved
        files, or no matching .recipe.json at all) — ask directly, one dialog
        per sample, rather than dropping it from the figure. Offers a file
        first (e.g. a dmfit fit) and, if cancelled, an EXPNO/pdata folder
        (a Bruker "1r") instead. Starts browsing from the configured data
        root(s) (Data folders… button) when one is set, so relocating is
        usually one click rather than a full re-browse."""
        roots = self._grid_data_roots()
        start = next((r for r in roots if Path(r).exists()), roots[0] if roots else "")
        path, _ = QFileDialog.getOpenFileName(
            self, f"Locate data for '{sample}'", start,
            "dmfit (*.fxmla *.fxml);;All files (*)")
        if path:
            return path
        folder = QFileDialog.getExistingDirectory(
            self, f"…or pick the EXPNO/pdata folder for '{sample}' (containing 1r)",
            start)
        return folder or ""

    def _grid_resolve_manual(self):
        for p in self._panels:
            if not p.get("needs_manual"):
                continue
            got = self._grid_ask_manual_path(p["sample"])
            if got:
                p["data_path"] = got
                p["has_data"] = True
                p["needs_manual"] = False

    def _grid_refresh_list(self):
        self.gridList.blockSignals(True)
        self.gridList.clear()
        for i, p in enumerate(self._panels):
            label = p.get("title") or p["sample"]
            if p.get("needs_manual"):
                label += "  ⚠ locate data…"
            elif p.get("reconstructed"):
                label += "  (rebuilt from CSV)"
            elif not p.get("path") and p.get("data_path"):
                label += "  (data only, no fit)"
            item = QListWidgetItem(label)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if p.get("include", True) else Qt.Unchecked)
            item.setData(Qt.UserRole, i)
            self.gridList.addItem(item)
        self.gridList.blockSignals(False)

    def _grid_rows_reordered(self, *_args):
        """Native drag-and-drop reorder in the QListWidget (setDragDropMode
        InternalMove) -- fold the new visual order back into self._panels by
        each item's own stamped original index, then re-stamp for next time."""
        new_order = []
        for i in range(self.gridList.count()):
            old_idx = self.gridList.item(i).data(Qt.UserRole)
            if old_idx is not None and 0 <= old_idx < len(self._panels):
                new_order.append(self._panels[old_idx])
        if len(new_order) == len(self._panels):
            self._panels = new_order
        self._grid_refresh_list()
        self._refresh()

    def _grid_item_changed(self, item):
        r = self.gridList.row(item)
        if 0 <= r < len(self._panels):
            self._panels[r]["include"] = item.checkState() == Qt.Checked
            self._refresh()

    def _grid_item_double_clicked(self, item):
        r = self.gridList.row(item)
        if not (0 <= r < len(self._panels)):
            return
        p = self._panels[r]
        if p.get("needs_manual") or not (p.get("path") or p.get("data_path")):
            got = self._grid_ask_manual_path(p["sample"])
            if got:
                p["data_path"] = got; p["has_data"] = True; p["needs_manual"] = False
                self._grid_refresh_list(); self._refresh()
            return
        title, ok = QInputDialog.getText(self, "Panel title", "Title:",
                                         text=p.get("title") or p["sample"])
        if ok:
            p["title"] = title.strip() or None
            self._grid_refresh_list(); self._refresh()

    def _grid_move(self, delta):
        r = self.gridList.currentRow()
        nr = r + delta
        if r < 0 or not (0 <= nr < len(self._panels)):
            return
        self._panels[r], self._panels[nr] = self._panels[nr], self._panels[r]
        self._grid_refresh_list()
        self.gridList.setCurrentRow(nr)
        self._refresh()

    def _grid_remove(self):
        r = self.gridList.currentRow()
        if r >= 0:
            self._panels.pop(r)
            self._grid_refresh_list()
            self._refresh()

    def _grid_detect_sites(self) -> list[tuple[int, str]]:
        """(index, label) for the first loaded panel with a resolvable fit --
        used to populate the component colors/legend editor. Any panel works
        since the shared model's site count/labels are the same across a
        batch (an excluded/zeroed site is still listed here — it just won't
        draw or need a color, per render_batch_grid's own handling)."""
        from larmor.recipe import Recipe
        for p in self._panels:
            if not p.get("path"):
                continue
            try:
                rec = Recipe.load(p["path"])
            except Exception:
                continue
            if rec.sites:
                return [(i, s.label or f"s{i}") for i, s in enumerate(rec.sites)]
        return []

    def _grid_edit_components(self):
        sites = self._grid_detect_sites()
        if not sites:
            self.gridMsg.setText(
                "load panels with at least one resolvable fit first to edit "
                "component colors/legend")
            return
        dlg = QDialog(self)
        dlg.setWindowTitle("Component colors / legend")
        lay = QVBoxLayout(dlg)
        form = QFormLayout()
        rows: dict[int, tuple] = {}
        for i, label in sites:
            row = QWidget(); rl = QHBoxLayout(row); rl.setContentsMargins(0, 0, 0, 0)
            color = self._component_colors.get(i) or figures.site_color(i)
            btn = QPushButton(color)
            btn.setStyleSheet(f"background:{color};")

            def pick(_checked=False, ii=i, b=btn):
                c = QColorDialog.getColor()
                if c.isValid():
                    self._component_colors[ii] = c.name()
                    b.setText(c.name()); b.setStyleSheet(f"background:{c.name()};")
            btn.clicked.connect(pick)
            chk = QCheckBox("in legend")
            chk.setChecked(i not in self._legend_hide)
            rl.addWidget(btn, 1); rl.addWidget(chk)
            form.addRow(label, row)
            rows[i] = (btn, chk)
        lay.addLayout(form)
        b_reset = QPushButton("Reset colors/legend to defaults")
        b_reset.clicked.connect(lambda: self._grid_reset_components(dlg))
        lay.addWidget(b_reset)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject)
        lay.addWidget(bb)
        if dlg.exec() != QDialog.Accepted:
            return
        for i, (_btn, chk) in rows.items():
            if chk.isChecked():
                self._legend_hide.discard(i)
            else:
                self._legend_hide.add(i)
        self._refresh()

    def _grid_reset_components(self, dlg):
        self._component_colors.clear()
        self._legend_hide.clear()
        dlg.reject()
        self._refresh()

    # ------------------------------------------------------------------ species bar
    def _bar_load_csv(self):
        """Pivot one parameter (e.g. amplitude) out of a batch_table*.csv into
        the table, one row per sample/scope, one column per (site, param)
        found — the 100%-stacked normalization then happens at render time."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Batch table CSV", "", "CSV (*.csv)")
        if not path:
            return
        rows_by_scope = series_grid.csv_rows_by_scope(path)
        categories = [s for s in rows_by_scope if s.lower() != "shared"]
        if not categories:
            self.msg.setText("no per-spectrum rows found in that CSV")
            return
        available = sorted({r.get("param", "") for rows in rows_by_scope.values()
                            for r in rows} - {""})
        param, ok = QInputDialog.getItem(
            self, "Species distribution",
            "Parameter to pivot into species columns:", available, 0, False)
        if not ok or not param:
            return
        site_labels: list[tuple] = []
        seen = set()
        for cat in categories:
            for r in rows_by_scope[cat]:
                if r.get("param") != param:
                    continue
                key = (r.get("site"), r.get("label"))
                if key not in seen:
                    seen.add(key); site_labels.append(key)
        if not site_labels:
            self.msg.setText(f"no rows found for parameter '{param}'")
            return

        self.barTable.blockSignals(True)
        self.barTable.setRowCount(len(categories))
        self.barTable.setColumnCount(1 + len(site_labels))
        self.barTable.setHorizontalHeaderLabels(
            ["category"] + [(lbl or site) for site, lbl in site_labels])
        for r, cat in enumerate(categories):
            self.barTable.setItem(r, 0, QTableWidgetItem(cat))
            by_key = {(row.get("site"), row.get("label")): row.get("value")
                     for row in rows_by_scope[cat] if row.get("param") == param}
            for c, key in enumerate(site_labels):
                self.barTable.setItem(r, c + 1, QTableWidgetItem(str(by_key.get(key, "0"))))
        self.barTable.blockSignals(False)
        self.msg.setText(f"loaded {len(categories)} × {len(site_labels)} "
                         f"from {Path(path).name}")
        self._refresh()

    def _bar_add_row(self):
        r = self.barTable.rowCount()
        self.barTable.insertRow(r)
        self.barTable.setItem(r, 0, QTableWidgetItem(f"c{r + 1}"))
        self._schedule()

    def _bar_add_col(self):
        c = self.barTable.columnCount()
        self.barTable.insertColumn(c)
        self.barTable.setHorizontalHeaderItem(c, QTableWidgetItem(f"species{c}"))
        self._schedule()

    def _bar_remove_row(self):
        r = self.barTable.currentRow()
        if r >= 0 and self.barTable.rowCount() > 1:
            self.barTable.removeRow(r)
            self._schedule()

    def _bar_remove_col(self):
        c = self.barTable.currentColumn()
        if c >= 1:                                  # keep the category column
            self.barTable.removeColumn(c)
            self._schedule()

    def _bar_rename_column(self, index):
        if index == 0:
            return                                  # "category" header is fixed
        header = self.barTable.horizontalHeaderItem(index)
        current = header.text() if header else f"species{index}"
        name, ok = QInputDialog.getText(self, "Species name", "Name:", text=current)
        if ok and name.strip():
            self.barTable.setHorizontalHeaderItem(index, QTableWidgetItem(name.strip()))
            self._schedule()

    # ------------------------------------------------------------------ pickers
    def _pick_2d(self):
        folder = QFileDialog.getExistingDirectory(self, "2D EXPNO folder")
        if folder:
            self.path2d.setText(folder); self._refresh()

    def _pick_fit_2d(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "2D fit recipe", "", "LARMOR recipe (*.json);;All (*)")
        if not path:
            return
        self.fitRecipe2d.setText(path)
        try:
            from larmor.recipe import Recipe
            rec = Recipe.load(path)
            if rec.nucleus and not self.nuc2d.text():
                self.nuc2d.setText(rec.nucleus)
            if rec.larmor_frequency_MHz and not self.larmor2d.value():
                self.larmor2d.setValue(rec.larmor_frequency_MHz)
        except Exception:
            pass
        self._refresh()

    def _pick_series(self):
        folder = QFileDialog.getExistingDirectory(self, "Series EXPNO folder")
        if folder:
            self.pathSer.setText(folder); self._refresh()

    def _kind_changed(self):
        k = self.kind.currentIndex()
        self.box1d.setVisible(k == 0)
        self.box2d.setVisible(k == 1)
        self.boxSer.setVisible(k == 2)
        self.boxGrid.setVisible(k == 3)
        self.boxBar.setVisible(k == 4)
        self._refresh()

    # ------------------------------------------------------------------ templates
    def _template_changed(self):
        name = self.template.currentData()
        if name is None:
            self.templateDesc.setText("")
            self._trace_defaults = {}
            return
        tpl = figures.TEMPLATES[name]
        self.templateDesc.setText(tpl["description"])
        self._trace_defaults = dict(tpl.get("trace_defaults") or {})
        self.kind.setCurrentIndex(_KINDS.index(
            {"1d": "1D overlay / stack", "2d": "2D contour", "series": "Series",
             "batch_grid": "Batch grid",
             "species_bar": "Species distribution"}[tpl["kind"]]))
        self._apply_common(tpl["spec"])
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
        if self.chkYlim.isChecked():
            common["ylim"] = (self.ylo.value(), self.yhi.value())
        if self.xstep.value():
            common["xtick_step"] = self.xstep.value()
        if self.ystep.value():
            common["ytick_step"] = self.ystep.value()
        if self.chkMinor.isChecked():
            common["minor_ticks"] = True
        common["tick_direction"] = self.tickdir.currentText()
        common["legend_loc"] = self.legloc.currentText()
        common["legend_ncol"] = self.legncol.value()
        if self.chkGrid.isChecked():
            common["grid"] = True
        if self.fontsz.value():
            common["font_size"] = self.fontsz.value()
        if self.lwd.value():
            common["line_width"] = self.lwd.value()
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
            spec["x_is_ppm"] = self.chkPpm.isChecked()
            if not self.chkPpm.isChecked():
                spec["hide_yaxis"] = False
            ticks = self._parse_xticks()
            if ticks:
                spec["xticks"] = ticks
                spec["xtick_rotation"] = 45
            return spec
        if k == 1:
            spec = {"kind": "2d", "path": self.path2d.text(),
                    "contour_mode": self.cmode.currentText(),
                    "cmap": self.cmap.currentText(),
                    "levels": {"n": self.nlev.value(),
                               "mode": self.levmode.currentText()},
                    "contour_values": self.chkValues.isChecked(),
                    "proj_top": self.chkTop.isChecked(),
                    "proj_right": self.chkRight.isChecked(),
                    "negative": self.chkNeg.isChecked(),
                    "iso_lines": self._iso, **common}
            if self.nuc2d.text().strip():
                spec["nucleus"] = self.nuc2d.text().strip()
            if self.fitRecipe2d.text().strip():
                spec["fit_recipe"] = self.fitRecipe2d.text().strip()
                spec["mqmas_method"] = self.mqmasMethod.currentText()
            return spec
        if k == 2:
            return {"kind": "series", "path": self.pathSer.text(),
                    "mode": self.serMode.currentText(),
                    "stretched": self.serStretch.isChecked(), **common}
        # batch_grid / species_bar: the shared "legend"/"title" fields don't
        # mean the same thing to these renderers as to 1d/2d/series (a fixed
        # shared legend, a figure suptitle) -- keep style/tick/limit fields
        # from `common` but leave legend_loc/legend_ncol out unless the user
        # actually left "best" for something else, so each renderer's own
        # sensible default legend placement (species_bar's row of swatches,
        # batch_grid's one shared legend) isn't silently overridden.
        common_grid = {kk: v for kk, v in common.items()
                       if kk not in ("legend_loc", "legend_ncol")}
        if self.legloc.currentText() != "best":
            common_grid["legend_loc"] = self.legloc.currentText()
            common_grid["legend_ncol"] = self.legncol.value()
        if k == 3:
            title = common_grid.pop("title", None)
            if title:
                common_grid["suptitle"] = title
            if common_grid.get("legend_loc") == "none":
                common_grid["legend"] = False
                common_grid.pop("legend_loc", None)
            panels = []
            for p in self._panels:
                if not p.get("include", True):
                    continue
                entry = {"title": p.get("title") or p.get("sample")}
                if p.get("path"):
                    entry["recipe"] = p["path"]
                elif p.get("data_path"):
                    entry["data_path"] = p["data_path"]
                if p.get("nucleus"):
                    entry["nucleus"] = p["nucleus"]
                panels.append(entry)
            spec = {"kind": "batch_grid", "panels": panels,
                    "component_mode": self.gridComp.currentText(),
                    "peak_labels": (self.gridLabels.currentText()
                                   if self.gridLabels.currentText() != "none" else None),
                    "show_total": self.gridTotal.isChecked(),
                    "show_experiment": self.gridExp.isChecked(),
                    "x_is_ppm": self.chkPpm.isChecked(), **common_grid}
            if self.gridCols.value():
                spec["cols"] = self.gridCols.value()
            shade = []
            for tok in self.gridShade.text().split(","):
                tok = tok.strip()
                if tok:
                    try:
                        shade.append(int(tok))
                    except ValueError:
                        pass
            if shade:
                spec["shade_only"] = shade
            hide = []
            for tok in self.gridHide.text().split(","):
                tok = tok.strip()
                if tok:
                    try:
                        hide.append(int(tok))
                    except ValueError:
                        pass
            if hide:
                spec["hide_components"] = hide
            if self._legend_hide:
                spec["legend_hide"] = sorted(self._legend_hide)
            if self._component_colors:
                spec["component_colors"] = dict(self._component_colors)
            return spec
        # k == 4: species distribution
        nrows, ncols = self.barTable.rowCount(), self.barTable.columnCount()
        categories = []
        for r in range(nrows):
            it = self.barTable.item(r, 0)
            categories.append(it.text() if it else "")
        series = []
        for c in range(1, ncols):
            hdr = self.barTable.horizontalHeaderItem(c)
            values = []
            for r in range(nrows):
                it = self.barTable.item(r, c)
                try:
                    values.append(float(it.text()) if it and it.text().strip() else 0.0)
                except ValueError:
                    values.append(0.0)
            series.append({"label": hdr.text() if hdr else f"s{c - 1}",
                          "values": values})
        return {"kind": "species_bar", "categories": categories, "series": series,
                "value_labels": self.barValueLabels.isChecked(), **common_grid}

    def _parse_xticks(self):
        """The 'x-ticks' field → [[pos, label], …]. Accepts 'a, b, c' (placed at
        1,2,3…) or explicit 'pos:label' pairs."""
        txt = self.xticks.text().strip()
        if not txt:
            return None
        out = []
        for i, part in enumerate(p for p in txt.split(",") if p.strip()):
            part = part.strip()
            if ":" in part:
                pos, lab = part.split(":", 1)
                try:
                    out.append([float(pos), lab.strip()])
                except ValueError:
                    out.append([float(i + 1), part])
            else:
                out.append([float(i + 1), part])
        return out or None

    def _apply_common(self, spec: dict):
        """Apply whichever COMMON fields `spec` sets (style/labels/limits/
        ticks/legend/font) without touching kind-specific data — used both
        for a full "Load spec…" and for a template's partial pre-fill."""
        if spec.get("style"):
            self.style.setCurrentText(spec["style"])
        if "title" in spec:
            self.title.setText(spec.get("title", ""))
        if "xlabel" in spec:
            self.xlabel.setText(spec.get("xlabel", ""))
        if "ylabel" in spec:
            self.ylabel.setText(spec.get("ylabel", ""))
        if "x_is_ppm" in spec:
            self.chkPpm.setChecked(spec["x_is_ppm"])
        if spec.get("xticks"):
            self.xticks.setText(", ".join(
                f"{p:g}:{lab}" for p, lab in spec["xticks"]))
        if spec.get("xlim"):
            self.chkXlim.setChecked(True)
            self.xhi.setValue(spec["xlim"][0]); self.xlo.setValue(spec["xlim"][1])
        if spec.get("ylim"):
            self.chkYlim.setChecked(True)
            self.ylo.setValue(spec["ylim"][0]); self.yhi.setValue(spec["ylim"][1])
        if "xtick_step" in spec:
            self.xstep.setValue(float(spec.get("xtick_step", 0) or 0))
        if "ytick_step" in spec:
            self.ystep.setValue(float(spec.get("ytick_step", 0) or 0))
        if "minor_ticks" in spec:
            self.chkMinor.setChecked(bool(spec["minor_ticks"]))
        if "tick_direction" in spec:
            self.tickdir.setCurrentText(spec["tick_direction"])
        if "grid" in spec:
            self.chkGrid.setChecked(bool(spec["grid"]))
        if "legend_loc" in spec:
            self.legloc.setCurrentText(spec["legend_loc"])
        if "legend_ncol" in spec:
            self.legncol.setValue(int(spec["legend_ncol"]))
        if "font_size" in spec:
            self.fontsz.setValue(float(spec.get("font_size", 0) or 0))
        if "line_width" in spec:
            self.lwd.setValue(float(spec.get("line_width", 0) or 0))

    def _apply_spec(self, spec: dict):
        self._apply_common(spec)
        kind = {"1d": 0, "2d": 1, "series": 2, "batch_grid": 3,
               "species_bar": 4}.get(spec.get("kind", "1d"), 0)
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
            self.nuc2d.setText(spec.get("nucleus", ""))
            self.fitRecipe2d.setText(spec.get("fit_recipe", ""))
            self.mqmasMethod.setCurrentText(spec.get("mqmas_method", "3QMAS"))
        elif kind == 2:
            self.pathSer.setText(spec.get("path", ""))
            self.serMode.setCurrentText(spec.get("mode", "satrec"))
        elif kind == 3:
            self._panels = []
            for p in spec.get("panels", []):
                recipe = p.get("recipe", "")
                data_path = p.get("data_path", "")
                sample = (p.get("title") or (Path(recipe).stem if recipe else None)
                         or "?")
                self._panels.append({
                    "path": recipe, "sample": sample, "nucleus": p.get("nucleus", ""),
                    "models": [], "has_data": bool(recipe or data_path),
                    "data_path": data_path, "needs_manual": False,
                    "reconstructed": False, "title": p.get("title"), "include": True,
                })
            self._grid_refresh_list()
            self.gridCols.setValue(int(spec.get("cols") or 0))
            self.gridComp.setCurrentText(spec.get("component_mode", "fill"))
            self.gridShade.setText(",".join(str(i) for i in spec.get("shade_only", [])))
            self.gridHide.setText(",".join(str(i) for i in spec.get("hide_components", [])))
            self.gridLabels.setCurrentText(spec.get("peak_labels") or "none")
            self.gridTotal.setChecked(spec.get("show_total", True))
            self.gridExp.setChecked(spec.get("show_experiment", True))
            self._legend_hide = set(spec.get("legend_hide", []))
            self._component_colors = {int(k): v for k, v in
                                      (spec.get("component_colors") or {}).items()}
        elif kind == 4:
            cats = spec.get("categories", [])
            series = spec.get("series", [])
            self.barTable.blockSignals(True)
            self.barTable.setRowCount(max(len(cats), 1))
            self.barTable.setColumnCount(1 + len(series))
            self.barTable.setHorizontalHeaderLabels(
                ["category"] + [s.get("label", f"s{i}") for i, s in enumerate(series)])
            for r, cat in enumerate(cats):
                self.barTable.setItem(r, 0, QTableWidgetItem(str(cat)))
                for c, s in enumerate(series):
                    vals = s.get("values", [])
                    v = vals[r] if r < len(vals) else 0
                    self.barTable.setItem(r, c + 1, QTableWidgetItem(str(v)))
            self.barTable.blockSignals(False)
            self.barValueLabels.setChecked(bool(spec.get("value_labels", True)))

    # ------------------------------------------------------------------ render
    def _refresh(self):
        spec = self._spec()
        if spec["kind"] == "1d" and not spec["traces"]:
            self.msg.setText("add a trace to preview"); return
        if spec["kind"] in ("2d", "series") and not spec.get("path"):
            self.msg.setText("pick an EXPNO folder to preview"); return
        if spec["kind"] == "batch_grid" and not spec.get("panels"):
            self.msg.setText("load a batch CSV or a folder of saved fits, "
                             "then check the panels to include"); return
        if spec["kind"] == "species_bar" and not any(spec.get("categories", [])):
            self.msg.setText("add categories and species to the table "
                             "(or load them from a batch CSV)"); return
        try:
            # render straight onto the live canvas figure (crisp, no bitmap scaling)
            new_fig = figures.render(spec)
            old = self._canvas.figure
            self._canvas.figure = new_fig
            new_fig.set_canvas(self._canvas)
            self._canvas.draw_idle()
            import matplotlib.pyplot as plt
            plt.close(old)
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
