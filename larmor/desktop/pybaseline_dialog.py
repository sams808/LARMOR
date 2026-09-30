"""Baseline correction through pybaselines -- the interactive dialog.

Pick a method (grouped by family), tune its parameters and watch the
estimated baseline and the corrected spectrum update live, then Apply. The
registry, the defaults, the validation and the computation live in
``larmor.pybaseline`` (Qt-free); this module only renders them, and
``params()`` hands back ``{"method": key, **values}`` ready for the
``pybaseline`` processing op.

A spectrum longer than PREVIEW_MAX points is previewed on every n-th point
with λ and the windows rescaled (``pybaseline.preview_params``) so the
preview stays instant on a 64k-point spectrum; the status line says so, and
Apply records the full-length parameters that the op then runs on all
points.
"""
from __future__ import annotations

import math
import time

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QStandardItem, QStandardItemModel, QValidator
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFormLayout, QHBoxLayout, QLabel, QPushButton, QSlider, QSpinBox,
    QVBoxLayout, QWidget,
)

from larmor import pybaseline as pyb
from larmor.desktop import theme
from larmor.desktop.axes import plain_units

#: preview on at most this many points (a stride; Apply uses them all)
PREVIEW_MAX = 8192
#: slider steps per decade of a log-scale parameter (ticks every decade)
_PER_DECADE = 10
_STATUS_OK = "color:#555"
_STATUS_WARN = "color:#8a5a00"


class _SciSpinBox(QDoubleSpinBox):
    """A spin box for a log-scale parameter: shows four significant digits
    (``1e+05`` above 1e4), accepts typed scientific notation, and an arrow
    step doubles / halves instead of adding a fixed amount."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDecimals(10)
        self.setKeyboardTracking(False)

    def textFromValue(self, value: float) -> str:      # noqa: N802 - Qt name
        return "0" if value == 0 else f"{value:.4g}"

    def valueFromText(self, text: str) -> float:       # noqa: N802 - Qt name
        try:
            return float(text.strip().replace(",", "."))
        except ValueError:
            return self.value()

    def validate(self, text: str, pos: int):           # noqa: D102
        t = text.strip().replace(",", ".")
        if t in ("", "-", "+", ".") or t[-1] in "eE+-.":
            return (QValidator.State.Intermediate, text, pos)
        try:
            float(t)
        except ValueError:
            return (QValidator.State.Invalid, text, pos)
        return (QValidator.State.Acceptable, text, pos)

    def stepBy(self, steps: int):                      # noqa: N802 - Qt name
        v = self.value() * 2.0 ** steps
        self.setValue(min(max(v, self.minimum()), self.maximum()))


class PybaselineDialog(QDialog):
    """Preview + apply a pybaselines baseline (Erb 2022)."""

    def __init__(self, parent, ppm: np.ndarray, amp: np.ndarray,
                 method: str | None = None):
        super().__init__(parent)
        self.setWindowTitle("Baseline — pybaselines")
        self.resize(900, 720)
        self.ppm = np.asarray(ppm, float)
        self.amp = np.asarray(amp, float)
        self._px, self._py, self._factor = pyb.decimate(self.ppm, self.amp,
                                                        PREVIEW_MAX)
        self._widgets: dict[str, QWidget] = {}
        self._sliders: dict[str, QSlider] = {}
        self._form_widget: QWidget | None = None
        self._method: str | None = None

        v = QVBoxLayout(self)
        intro = QLabel(
            "Baseline algorithms of the <b>pybaselines</b> library: the "
            "Whittaker smoothers (a smooth, slowly curving baseline under "
            "peaks of any width — start with arPLS), the iterative "
            "polynomials (a gently curved baseline) and the window methods "
            "(narrow lines on a wide rolling baseline). The estimated "
            "baseline and the corrected spectrum update as you change a "
            "control; Apply records the step in the recipe.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        # ---- plots: spectrum + baseline (top), corrected (bottom) ----
        t = theme.active()
        glw = pg.GraphicsLayoutWidget()
        glw.setBackground(t.plot_bg)
        self.p_top = glw.addPlot(row=0, col=0)
        plain_units(self.p_top)
        self.p_top.showGrid(x=True, y=True, alpha=0.12)
        self.p_top.setLabel("left", "intensity")
        self.p_top.addLegend(offset=(-10, 10))
        self.p_bot = glw.addPlot(row=1, col=0)
        plain_units(self.p_bot)
        self.p_bot.showGrid(x=True, y=True, alpha=0.12)
        self.p_bot.setLabel("left", "corrected")
        self.p_bot.setLabel("bottom", "shift", units="ppm")
        # addPlot() returns a PlotItem: getAxis on it, never getPlotItem
        self.p_bot.getAxis("bottom").enableAutoSIPrefix(False)
        self.p_top.getAxis("bottom").enableAutoSIPrefix(False)
        self.p_bot.setXLink(self.p_top)
        glw.ci.layout.setRowStretchFactor(0, 3)
        glw.ci.layout.setRowStretchFactor(1, 2)
        v.addWidget(glw, 1)

        self.c_raw = self.p_top.plot(pen=pg.mkPen(t.experiment, width=1),
                                     name="spectrum", connect="finite")
        self.c_base = self.p_top.plot(pen=pg.mkPen(t.baseline, width=1.6),
                                      name="baseline", connect="finite")
        self.c_corr = self.p_bot.plot(pen=pg.mkPen(t.measure, width=1),
                                      name="corrected", connect="finite")
        if self.ppm.size:
            for p in (self.p_top, self.p_bot):
                p.setXRange(self.ppm.max(), self.ppm.min())   # NMR: decreasing

        # ---- method ----
        row = QHBoxLayout()
        lbl = QLabel("Method:")
        self.combo = QComboBox()
        self.combo.setToolTip("Grouped by family; the line under the box says "
                              "which baseline shape the method suits.")
        model = QStandardItemModel(self.combo)
        for fam, heading in pyb.FAMILIES.items():
            head = QStandardItem(heading)
            head.setFlags(Qt.ItemFlag.NoItemFlags)       # a heading, not a choice
            f = head.font()
            f.setBold(True)
            head.setFont(f)
            model.appendRow(head)
            for spec in pyb.METHODS.values():
                if spec.family != fam:
                    continue
                it = QStandardItem("    " + spec.label)
                it.setData(spec.key, Qt.ItemDataRole.UserRole)
                it.setToolTip(spec.description)
                model.appendRow(it)
        self.combo.setModel(model)
        row.addWidget(lbl)
        row.addWidget(self.combo, 1)
        v.addLayout(row)
        self.description = QLabel("")
        self.description.setWordWrap(True)
        self.description.setStyleSheet("color:#555; font-style:italic")
        v.addWidget(self.description)

        # ---- the dynamic parameter form ----
        slot_host = QWidget()
        self._form_slot = QVBoxLayout(slot_host)
        self._form_slot.setContentsMargins(0, 0, 0, 0)
        v.addWidget(slot_host)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setStyleSheet(_STATUS_OK)
        v.addWidget(self.status)

        cite = QLabel(
            f"Method library: {pyb.CITATION} — cite it, and the method's own "
            f"paper, if you publish results. Documentation: {pyb.DOCS_URL}")
        cite.setWordWrap(True)
        cite.setStyleSheet("color:#666; font-size:11px")
        v.addWidget(cite)

        # ---- buttons: ? · Reset defaults · … · Apply / Cancel ----
        self.btn_help = QPushButton("?")
        self.btn_help.setFixedWidth(28)
        self.btn_help.setToolTip("Open the Processing reference (pybaselines "
                                 "section: which method when, the λ rule of "
                                 "thumb, the citations)")
        self.btn_help.clicked.connect(self._help)
        self.btn_reset = QPushButton("Reset defaults")
        self.btn_reset.setToolTip("Back to the method's defaults for a "
                                  f"{self.amp.size}-point spectrum")
        self.btn_reset.clicked.connect(self.reset_defaults)
        self.bb = QDialogButtonBox(QDialogButtonBox.Apply | QDialogButtonBox.Cancel)
        self.btn_apply = self.bb.button(QDialogButtonBox.Apply)
        self.btn_apply.clicked.connect(self.accept)
        self.bb.rejected.connect(self.reject)
        brow = QHBoxLayout()
        brow.addWidget(self.btn_help)
        brow.addWidget(self.btn_reset)
        brow.addStretch(1)
        brow.addWidget(self.bb)
        v.addLayout(brow)

        # live preview, debounced so a slider drag stays smooth
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(150)
        self._timer.timeout.connect(self._preview)

        self.combo.currentIndexChanged.connect(self._on_method_changed)
        self.set_method(method if method in pyb.METHODS else next(iter(pyb.METHODS)))
        # the first preview runs once the dialog is on screen (the first
        # pybaselines import can take a few seconds: never a white window)
        self.status.setText("computing the preview…")
        QTimer.singleShot(50, self._preview)

    # ------------------------------------------------------------ method
    def method(self) -> str:
        """The selected method's key -- always a METHODS key: a family
        heading made current programmatically is bounced back."""
        return self._method

    def method_label(self) -> str:
        return pyb.METHODS[self.method()].label

    def _row_of(self, key: str) -> int | None:
        for i in range(self.combo.count()):
            if self.combo.itemData(i, Qt.ItemDataRole.UserRole) == key:
                return i
        return None

    def set_method(self, key: str):
        """Select ``key`` and build its form (unknown keys are ignored)."""
        row = self._row_of(key)
        if row is None:
            return
        if self.combo.currentIndex() == row:
            self._on_method_changed(row)
        else:
            self.combo.setCurrentIndex(row)

    def _on_method_changed(self, index: int):
        key = self.combo.itemData(index, Qt.ItemDataRole.UserRole)
        if key is None:                    # a heading: the popup never offers
            row = self._row_of(self._method) if self._method else None
            if row is not None and row != index:   # it, only code can land here
                self.combo.setCurrentIndex(row)
            return
        rebuild = key != self._method or not self._widgets
        self._method = key
        self.description.setText(pyb.METHODS[key].description)
        if rebuild:
            self._build_form(pyb.default_params(key, self.amp.size))
        self._timer.start()

    # -------------------------------------------------------------- form
    def _build_form(self, values: dict):
        spec = pyb.METHODS[self.method()]
        n = self.amp.size
        if self._form_widget is not None:
            self._form_slot.removeWidget(self._form_widget)
            self._form_widget.setParent(None)
            self._form_widget.deleteLater()
        self._widgets = {}
        self._sliders = {}
        host = QWidget()
        form = QFormLayout(host)
        form.setContentsMargins(0, 0, 0, 0)
        for p in spec.params:
            val = values.get(p.key, p.default)
            hi = pyb.resolved_max(p, n)
            row: QWidget
            if p.kind == "bool":
                w = QCheckBox()
                w.setChecked(bool(val))
                w.toggled.connect(self._changed)
                row = w
            elif p.kind == "choice":
                w = QComboBox()
                for c in p.choices:
                    w.addItem(str(c), c)
                idx = list(p.choices).index(val) if val in p.choices else 0
                w.setCurrentIndex(idx)
                w.currentIndexChanged.connect(self._changed)
                row = w
            elif p.kind == "int":
                w = QSpinBox()
                w.setKeyboardTracking(False)          # act on "293", not "29"
                w.setRange(int(p.min if p.min is not None else -10 ** 9),
                           int(hi if hi is not None else 10 ** 9))
                w.setSingleStep(int(p.step or 1))
                w.setValue(int(val))
                w.valueChanged.connect(self._changed)
                row = w
            elif p.log:
                w = _SciSpinBox()
                w.setRange(float(p.min), float(hi))
                w.setValue(float(val))
                sl = QSlider(Qt.Orientation.Horizontal)
                sl.setRange(round(math.log10(float(p.min)) * _PER_DECADE),
                            round(math.log10(float(hi)) * _PER_DECADE))
                sl.setTickInterval(_PER_DECADE)
                sl.setTickPosition(QSlider.TickPosition.TicksBelow)
                sl.setValue(round(math.log10(float(val)) * _PER_DECADE))
                sl.setToolTip("one tick = one decade")
                self._link_log(w, sl)
                w.valueChanged.connect(self._changed)
                box = QWidget()
                h = QHBoxLayout(box)
                h.setContentsMargins(0, 0, 0, 0)
                h.addWidget(w)
                h.addWidget(sl, 1)
                self._sliders[p.key] = sl
                row = box
            else:
                w = QDoubleSpinBox()
                w.setKeyboardTracking(False)
                step = float(p.step or 0.1)
                w.setDecimals(max(2, 1 - int(math.floor(math.log10(step)))))
                w.setRange(float(p.min if p.min is not None else -1e12),
                           float(hi if hi is not None else 1e12))
                w.setSingleStep(step)
                w.setValue(float(val))
                w.valueChanged.connect(self._changed)
                row = w
            w.setToolTip(p.tooltip)
            row.setToolTip(p.tooltip)
            lbl = QLabel(p.label + ":")
            lbl.setToolTip(p.tooltip)
            form.addRow(lbl, row)
            self._widgets[p.key] = w
        self._form_slot.addWidget(host)
        self._form_widget = host

    @staticmethod
    def _link_log(spin: QDoubleSpinBox, slider: QSlider):
        """Keep a log-scale spin box and its decade slider in step."""
        guard = {"busy": False}

        def from_slider(ticks: int):
            if guard["busy"]:
                return
            guard["busy"] = True
            try:
                spin.setValue(10.0 ** (ticks / _PER_DECADE))
            finally:
                guard["busy"] = False

        def from_spin(value: float):
            if guard["busy"] or value <= 0:
                return
            guard["busy"] = True
            try:
                slider.setValue(round(math.log10(value) * _PER_DECADE))
            finally:
                guard["busy"] = False

        slider.valueChanged.connect(from_slider)
        spin.valueChanged.connect(from_spin)

    def _changed(self, *_):
        self._timer.start()

    def values(self) -> dict:
        """The current parameter values, typed as the op records them."""
        spec = pyb.METHODS[self.method()]
        out: dict = {}
        for p in spec.params:
            w = self._widgets[p.key]
            if p.kind == "bool":
                out[p.key] = bool(w.isChecked())
            elif p.kind == "choice":
                out[p.key] = w.currentData()
            elif p.kind == "int":
                out[p.key] = int(w.value())
            else:
                out[p.key] = float(w.value())
        return out

    def params(self) -> dict:
        """``{"method": key, **values}`` -- the ``pybaseline`` op's kwargs."""
        return {"method": self.method(), **self.values()}

    def set_values(self, values: dict):
        """Put ``values`` into the widgets without a preview per widget."""
        spec = pyb.METHODS[self.method()]
        for p in spec.params:
            if p.key not in values:
                continue
            w = self._widgets[p.key]
            val = values[p.key]
            w.blockSignals(True)
            try:
                if p.kind == "bool":
                    w.setChecked(bool(val))
                elif p.kind == "choice":
                    i = w.findData(val)
                    w.setCurrentIndex(max(i, 0))
                else:
                    w.setValue(val)
                    sl = self._sliders.get(p.key)
                    if sl is not None and float(val) > 0:
                        sl.blockSignals(True)
                        sl.setValue(round(math.log10(float(val)) * _PER_DECADE))
                        sl.blockSignals(False)
            finally:
                w.blockSignals(False)
        self._timer.start()

    def reset_defaults(self):
        self.set_values(pyb.default_params(self.method(), self.amp.size))

    # ----------------------------------------------------------- preview
    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "processing-reference", "Processing reference")

    def _show_error(self, text: str):
        self.status.setText(text)
        self.status.setStyleSheet(_STATUS_WARN)
        self.c_base.setData([], [])            # never a stale baseline
        self.c_corr.setData([], [])
        self.btn_apply.setEnabled(False)

    def _preview(self):
        self._timer.stop()
        if not self._widgets:
            return
        method = self.method()
        spec = pyb.METHODS[method]
        values = self.values()
        decimated = self._factor > 1
        x, y = (self._px, self._py) if decimated else (self.ppm, self.amp)
        t0 = time.perf_counter()
        try:
            kw = pyb.preview_params(method, values, self._factor) if decimated else values
            base = pyb.compute(x, y, method, **kw)
        except (ValueError, ImportError) as exc:
            self._show_error(str(exc))
            return
        except Exception as exc:                # noqa: BLE001 - keep the dialog alive
            self._show_error(f"{spec.label} failed: {exc}")
            return
        dt = time.perf_counter() - t0
        self.c_raw.setData(x, y)
        self.c_base.setData(x, base)
        self.c_corr.setData(x, y - base)
        peak = float(np.nanmax(np.abs(base))) if base.size else 0.0
        ref = float(np.nanmax(np.abs(y))) if y.size else 0.0
        msg = (f"{spec.label} · baseline peak {peak:.3g} "
               f"({100 * peak / (ref or 1.0):.1f} % of the spectrum) · "
               f"{dt * 1000:.0f} ms")
        if decimated:
            msg += (f" · preview on every {self._factor}th point "
                    f"({y.size} of {self.amp.size}); Apply uses all points")
        self.status.setText(msg)
        self.status.setStyleSheet(_STATUS_OK)
        self.btn_apply.setEnabled(True)
