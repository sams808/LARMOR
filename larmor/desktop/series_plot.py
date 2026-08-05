"""Series evolution plot: how fitted parameters change along a batch series.

Given a batch-fit result (one shared model, many spectra), plot any site's
δ_iso / width / C_Q / η / amplitude — or its amplitude-fraction population —
against the series index, with error bars where available. Export the plotted
numbers as CSV or the figure itself.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QDialog, QDialogButtonBox, QFileDialog, QHBoxLayout,
    QLabel, QListWidget, QListWidgetItem, QPushButton, QVBoxLayout, QWidget,
)

from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS
from larmor.desktop.plot import site_color


def series_options(result) -> list[dict]:
    """Every plottable (site, param) pair in the result, plus a derived
    amplitude-fraction population per site."""
    out: list[dict] = []
    master = result.recipes[0]
    for i, site in enumerate(master.sites):
        label = site.label or site.model
        for pn in site.params:
            if pn == "gl":
                continue
            out.append({"site": i, "param": pn, "kind": "param",
                        "text": f"s{i} {label}: {PARAM_LABELS.get(pn, pn)}"})
        out.append({"site": i, "param": "amplitude", "kind": "popfrac",
                    "text": f"s{i} {label}: population % (by amplitude)"})
    return out


def series_values(result, opt: dict):
    """(values, errors) of one option across every spectrum in the series."""
    vals, errs = [], []
    for rec in result.recipes:
        site = rec.sites[opt["site"]]
        p = site.params.get(opt["param"])
        v = float(p.value) if p is not None else np.nan
        e = float(p.stderr) if (p is not None and p.stderr is not None) else np.nan
        if opt["kind"] == "popfrac":
            tot = sum(abs(float(s.params["amplitude"].value))
                      for s in rec.sites if "amplitude" in s.params) or 1.0
            v = 100.0 * abs(v) / tot
            e = np.nan
        vals.append(v); errs.append(e)
    return np.array(vals, float), np.array(errs, float)


class SeriesPlotDialog(QDialog):
    def __init__(self, parent, result):
        super().__init__(parent)
        self.setWindowTitle("Series evolution")
        self.resize(820, 560)
        self._result = result
        self._opts = series_options(result)
        self._labels = list(result.labels)

        v = QVBoxLayout(self)
        v.addWidget(QLabel("Pick one or more series to plot (Ctrl-click for "
                           "several); they share the spectrum-index axis."))
        body = QHBoxLayout(); v.addLayout(body, 1)

        self.list = QListWidget()
        self.list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.list.setMaximumWidth(280)
        for o in self._opts:
            it = QListWidgetItem(o["text"]); it.setData(Qt.UserRole, o)
            self.list.addItem(it)
        if self._opts:
            self.list.item(0).setSelected(True)
        self.list.itemSelectionChanged.connect(self._plot)
        body.addWidget(self.list)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        self.plot.addLegend()
        self.plot.setLabel("bottom", "spectrum")
        from larmor.desktop.plot_menu import attach_plot_menu
        attach_plot_menu(self.plot, title="series", parent=self)
        body.addWidget(self.plot, 1)

        btns = QHBoxLayout()
        b_csv = QPushButton("Export parameters (CSV)…")
        b_csv.clicked.connect(self._export_csv)
        b_fig = QPushButton("Export figure…")
        b_fig.clicked.connect(self._export_fig)
        btns.addWidget(b_csv); btns.addWidget(b_fig); btns.addStretch(1)
        v.addLayout(btns)

        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(self.reject)
        bb.accepted.connect(self.accept)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        v.addWidget(bb)
        self._plot()

    def _selected(self) -> list[dict]:
        return [it.data(Qt.UserRole) for it in self.list.selectedItems()]

    def _plot(self):
        self.plot.clear()
        try:
            self.plot.plotItem.legend.clear()
        except Exception:
            pass
        x = np.arange(1, len(self._labels) + 1)
        ax = self.plot.getAxis("bottom")
        ax.setTicks([[(int(i), lab) for i, lab in zip(x, self._labels)]])
        for k, o in enumerate(self._selected()):
            vals, errs = series_values(self._result, o)
            col = site_color(o["site"])
            self.plot.plot(x, vals, pen=pg.mkPen(col, width=2),
                           symbol="o", symbolBrush=col, symbolSize=8,
                           name=o["text"])
            if np.isfinite(errs).any():
                err = pg.ErrorBarItem(x=x, y=vals,
                                      height=2 * np.nan_to_num(errs), pen=col)
                self.plot.addItem(err)
        self.plot.setLabel("left", self._selected()[0]["text"]
                           if len(self._selected()) == 1 else "value")

    def _export_csv(self):
        sel = self._selected()
        if not sel:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export series", "series.csv", "CSV (*.csv)")
        if not path:
            return
        import csv
        cols = {"spectrum": self._labels}
        for o in sel:
            vals, errs = series_values(self._result, o)
            cols[o["text"]] = vals
            if np.isfinite(errs).any():
                cols[o["text"] + " ±"] = errs
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(cols.keys())
            for i in range(len(self._labels)):
                w.writerow([cols[c][i] if not isinstance(cols[c][i], float)
                            else f"{cols[c][i]:.6g}" for c in cols])

    def _export_fig(self):
        from larmor.desktop.export_dialog import export_pyqtgraph
        export_pyqtgraph(self, self.plot.getPlotItem(), "series")
