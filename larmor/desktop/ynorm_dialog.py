"""View > Y axis > Normalise to area of a region...: the ppm range whose
trapezoid area the plotting area scales to 1 (larmor.display.y_factor,
mode "region"). Typed in, or taken from the fit zones or the current view.
"""
from __future__ import annotations

from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QDoubleSpinBox,
                               QFormLayout, QHBoxLayout, QLabel, QPushButton,
                               QVBoxLayout)


class RegionDialog(QDialog):
    """Two ppm bounds with 'Use fit zones' / 'Use current view' shortcuts;
    ``region()`` is (hi, lo) whatever order was typed."""

    def __init__(self, parent=None, region=None, zones=None, view_range=None):
        super().__init__(parent)
        self.setWindowTitle("Normalise to the area of a region")
        self._zones = [list(z) for z in (zones or []) if len(z) == 2]
        self._view = tuple(view_range) if view_range else None
        v = QVBoxLayout(self)
        intro = QLabel("The trapezoid area of the active spectrum over this "
                       "range reads 1, and each compared spectrum is scaled to "
                       "its own area over it. Display only: the fit and every "
                       "export keep the raw intensities.")
        intro.setWordWrap(True)
        v.addWidget(intro)
        form = QFormLayout()
        self.hi = QDoubleSpinBox()
        self.lo = QDoubleSpinBox()
        for sb in (self.hi, self.lo):
            sb.setRange(-1e6, 1e6)
            sb.setDecimals(2)
            sb.setSuffix(" ppm")
        form.addRow("from (high ppm)", self.hi)
        form.addRow("to (low ppm)", self.lo)
        v.addLayout(form)
        row = QHBoxLayout()
        self.btnZones = QPushButton("Use fit zones")
        self.btnZones.setToolTip("the span of the dmfit-style fit zones "
                                 "(or the recipe's fit window)")
        self.btnZones.setEnabled(bool(self._zones))
        self.btnZones.clicked.connect(self.use_zones)
        self.btnView = QPushButton("Use current view")
        self.btnView.setToolTip("the ppm window shown on the plot")
        self.btnView.setEnabled(self._view is not None)
        self.btnView.clicked.connect(self.use_view)
        row.addWidget(self.btnZones)
        row.addWidget(self.btnView)
        row.addStretch(1)
        v.addLayout(row)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
        if region is not None:
            self.set_region(region)
        elif self._view is not None:
            self.use_view()

    def set_region(self, region):
        self.hi.setValue(float(max(region)))
        self.lo.setValue(float(min(region)))

    def use_zones(self):
        """The union of the fit zones (their outermost bounds)."""
        if self._zones:
            flat = [float(b) for z in self._zones for b in z]
            self.set_region((max(flat), min(flat)))

    def use_view(self):
        if self._view is not None:
            self.set_region(self._view)

    def region(self) -> tuple[float, float]:
        a, b = float(self.hi.value()), float(self.lo.value())
        return (max(a, b), min(a, b))

    def accept(self):
        if self.hi.value() == self.lo.value():      # an empty region: stay open
            self.lo.setFocus()
            self.lo.selectAll()
            return
        super().accept()
