"""Reusable graphical-export options: format, DPI and physical size in cm.

One dialog, two sinks:
  * ``export_pyqtgraph`` — save a live pyqtgraph PlotItem (PNG / TIFF / JPG raster
    at the requested px = cm·dpi, or SVG vector);
  * ``export_matplotlib`` — save a Matplotlib Figure (PNG / PDF / SVG / TIFF / EPS)
    sized in cm at the requested DPI.

Used everywhere LARMOR exports a figure so the controls are identical.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QSpinBox,
)

CM_PER_IN = 2.54

_RASTER = {"PNG": "png", "TIFF": "tiff", "JPG": "jpg"}
_VECTOR = {"SVG": "svg", "PDF": "pdf", "EPS": "eps"}


class ExportOptions(QDialog):
    """Pick a format, a DPI and a size in centimetres."""

    def __init__(self, parent, formats: list[str], *, dpi=300,
                 width_cm=12.0, height_cm=9.0):
        super().__init__(parent)
        self.setWindowTitle("Export figure")
        form = QFormLayout(self)

        self.cbFormat = QComboBox(); self.cbFormat.addItems(formats)
        form.addRow("Format", self.cbFormat)

        self.sbDpi = QSpinBox(); self.sbDpi.setRange(50, 1200); self.sbDpi.setValue(dpi)
        self.sbDpi.setSuffix(" dpi")
        form.addRow("Resolution", self.sbDpi)

        self.sbW = QDoubleSpinBox(); self.sbW.setRange(1, 100); self.sbW.setValue(width_cm)
        self.sbW.setSuffix(" cm"); self.sbW.setDecimals(1)
        form.addRow("Width", self.sbW)
        self.sbH = QDoubleSpinBox(); self.sbH.setRange(1, 100); self.sbH.setValue(height_cm)
        self.sbH.setSuffix(" cm"); self.sbH.setDecimals(1)
        form.addRow("Height", self.sbH)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        form.addRow(bb)

    def values(self) -> dict:
        return {"format": self.cbFormat.currentText(), "dpi": self.sbDpi.value(),
                "width_cm": self.sbW.value(), "height_cm": self.sbH.value()}


def choose(parent, formats, **kw) -> dict | None:
    dlg = ExportOptions(parent, formats, **kw)
    return dlg.values() if dlg.exec() == QDialog.Accepted else None


def _ask_path(parent, default_name, ext) -> str:
    path, _ = QFileDialog.getSaveFileName(
        parent, "Save figure", f"{default_name}.{ext}",
        f"{ext.upper()} (*.{ext});;All files (*)")
    return path


def export_pyqtgraph(parent, plotitem, default_name="figure") -> str | None:
    """Export a pyqtgraph PlotItem with the shared options dialog."""
    opt = choose(parent, list(_RASTER) + ["SVG"])
    if opt is None:
        return None
    ext = {**_RASTER, "SVG": "svg"}[opt["format"]]
    path = _ask_path(parent, default_name, ext)
    if not path:
        return None
    if opt["format"] == "SVG":
        from pyqtgraph.exporters import SVGExporter
        SVGExporter(plotitem).export(path)
    else:
        from pyqtgraph.exporters import ImageExporter
        ex = ImageExporter(plotitem)
        px = int(round(opt["width_cm"] / CM_PER_IN * opt["dpi"]))
        try:
            ex.parameters()["width"] = px          # height follows the aspect
        except Exception:
            pass
        ex.export(path)
    return path


def export_matplotlib(parent, fig, default_name="figure") -> str | None:
    """Export a Matplotlib Figure with the shared options dialog."""
    opt = choose(parent, list(_RASTER) + list(_VECTOR))
    if opt is None:
        return None
    ext = {**_RASTER, **_VECTOR}[opt["format"]]
    path = _ask_path(parent, default_name, ext)
    if not path:
        return None
    fig.set_size_inches(opt["width_cm"] / CM_PER_IN, opt["height_cm"] / CM_PER_IN)
    fig.savefig(path, dpi=opt["dpi"], bbox_inches="tight")
    return path
