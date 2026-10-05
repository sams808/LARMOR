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
    from larmor.desktop.paths import FIGURE_DIR_KEY, remember_dir, remembered_dir
    start = remembered_dir(FIGURE_DIR_KEY)
    seed = (str(Path(start) / f"{default_name}.{ext}") if start
            else f"{default_name}.{ext}")
    path, _ = QFileDialog.getSaveFileName(
        parent, "Save figure", seed,
        f"{ext.upper()} (*.{ext});;All files (*)")
    if path:                # the NEXT export starts where this one landed
        remember_dir(FIGURE_DIR_KEY, path)
    return path


def _svg_closepath_shim() -> None:
    """pyqtgraph 0.14.0's SVG exporter parses every token of a path's ``d``
    as ``x,y``; Qt 6.11's QSvgGenerator closes shapes with a bare ``Z`` (the
    ViewBox background: ``M0,0 L603.5,0 … L0,0 Z``), so every SVG export of
    a live plot died with "not enough values to unpack" (found by the
    distribution check). Upstream pyqtgraph master passes single tokens
    through; until that release ships, the bare Z becomes the explicit line
    back to its subpath's start -- the geometry closepath stands for --
    before pyqtgraph's parser sees it. Installed once, idempotent, and a
    no-op for an exporter that already copes (there is no Z left to expand
    that it would not have handled)."""
    import importlib

    # the module, not the class of the same name the package re-exports
    mod = importlib.import_module("pyqtgraph.exporters.SVGExporter")
    orig = mod.correctCoordinates
    if getattr(orig, "_larmor_closepath", False):
        return

    def expand_closepath(d: str) -> str:
        out, start = [], None
        for tok in d.strip().split(" "):
            if not tok:
                continue
            if tok in ("Z", "z"):
                if start is not None:
                    out.append("L" + start)
                continue
            if tok[0] in "Mm":
                start = tok[1:]
            out.append(tok)
        return " ".join(out)

    def patched(node, defs, item, options):
        for el in node.getElementsByTagName("path"):
            d = el.getAttribute("d")
            if "Z" in d or "z" in d:
                el.setAttribute("d", expand_closepath(d))
        return orig(node, defs, item, options)

    patched._larmor_closepath = True
    mod.correctCoordinates = patched


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
        _svg_closepath_shim()
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
