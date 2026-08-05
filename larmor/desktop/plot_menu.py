"""A shared right-click menu for pyqtgraph plots.

Two things every plot in LARMOR should offer:
  * **Export figure…** through LARMOR's own dialog (format / DPI / size in cm) —
    the built-in pyqtgraph exporter can freeze under PySide6, so we disable it;
  * **Send to Plotting studio** — hand the plotted curves to the studio for a
    publication figure.

The studio traces are read straight off the plot's data items, so any 1D plot
works with a single ``attach_plot_menu(widget)`` call.
"""
from __future__ import annotations

from PySide6.QtGui import QAction


def disable_native_export_globally() -> None:
    """Neuter pyqtgraph's built-in Export dialog application-wide.

    Under PySide6 that dialog leaves a dangling ``selectBox`` QGraphicsRectItem
    and crashes the whole app on close (``Internal C++ object already deleted``).
    Replacing ``showExportDialog`` with a no-op guarantees it never opens on ANY
    plot — LARMOR's own 'Export figure…' (attach_plot_menu) does exporting."""
    try:
        import pyqtgraph as pg
        pg.GraphicsScene.showExportDialog = lambda self: None
    except Exception:
        pass


def _disable_native_export(widget) -> None:
    """Remove pyqtgraph's built-in 'Export…' (its dialog can hang under PySide6)."""
    try:
        widget.getPlotItem().scene().contextMenu = []
    except Exception:
        pass


def traces_from_plot(widget) -> list[dict]:
    """Figure-spec traces (inline data) for every non-empty curve on the plot."""
    out: list[dict] = []
    try:
        items = widget.getPlotItem().listDataItems()
    except Exception:
        return out
    for it in items:
        x = getattr(it, "xData", None)
        y = getattr(it, "yData", None)
        if x is None or y is None or len(x) == 0:
            continue
        try:
            name = it.name()
        except Exception:
            name = None
        out.append({"data": {"x": [float(v) for v in x],
                             "y": [float(v) for v in y]},
                    "label": name})
    return out


def attach_plot_menu(widget, *, title: str = "figure", parent=None,
                     studio: bool = True, studio_spec=None) -> None:
    """Add Export + Send-to-studio to a pyqtgraph PlotWidget's right-click menu
    and disable the freeze-prone native exporter.

    ``studio_spec`` (optional) is a callable returning a full figure spec for
    "Send to Plotting studio" — use it when the plot isn't an NMR spectrum (e.g. a
    parameter-vs-sample series), so the right data, axis and ticks transfer. When
    omitted, the plotted curves are sent as a generic 1D overlay."""
    _disable_native_export(widget)
    try:
        menu = widget.getPlotItem().getViewBox().menu
    except Exception:
        return
    menu.addSeparator()
    act_exp = QAction("Export figure…", menu)
    act_exp.triggered.connect(lambda: _export(widget, title, parent))
    menu.addAction(act_exp)
    if studio:
        act_std = QAction("Send to Plotting studio", menu)
        act_std.triggered.connect(lambda: _to_studio(widget, parent, studio_spec))
        menu.addAction(act_std)


def _export(widget, title, parent):
    from larmor.desktop.export_dialog import export_pyqtgraph
    try:
        export_pyqtgraph(parent or widget.window(), widget.getPlotItem(), title)
    except Exception as exc:  # noqa: BLE001
        from PySide6.QtWidgets import QMessageBox
        QMessageBox.warning(widget, "Export", f"export failed: {exc}")


def _to_studio(widget, parent, studio_spec=None):
    from larmor.desktop.plotting_studio import PlottingStudio
    if studio_spec is not None:
        spec = studio_spec()
    else:
        traces = [t for t in traces_from_plot(widget) if t["data"]["x"]]
        spec = {"kind": "1d", "traces": traces}
    PlottingStudio(parent or widget.window(), spec).exec()
