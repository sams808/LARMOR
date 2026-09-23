"""Series evolution plot: how fitted parameters change along a batch series.

Pick one or more **sites** (lines) on the left; the right shows **one subplot per
parameter** (δ_iso, width, C_Q, η, amplitude, and the integral **population %**),
each tracing the chosen sites across the series. Export the numbers as CSV or the
whole panel as a figure.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog,
    QGridLayout, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton,
    QScrollArea, QVBoxLayout, QWidget,
)

from larmor import series_table as stab
from larmor.desktop import theme
from larmor.desktop.panels import PARAM_LABELS
from larmor.desktop.plot import site_color


def _param_specs(result) -> list[dict]:
    """One entry per plottable parameter (union across the master sites, minus the
    Gauss/Lorentz mix), then amplitude and the integral population %."""
    master = result.recipes[0]
    seen: list[str] = []
    for s in master.sites:
        for pn in s.params:
            if pn == "gl" or pn == "amplitude" or pn in seen:
                continue
            seen.append(pn)
    specs = [{"param": pn, "kind": "param", "label": PARAM_LABELS.get(pn, pn)}
             for pn in seen]
    specs.append({"param": "amplitude", "kind": "param", "label": "amplitude"})
    specs.append({"param": "population_pct", "kind": "pop_integral",
                  "label": "population % (integral)"})
    return specs


def population_integral(result, error_method: str | None = None
                        ) -> tuple[np.ndarray, np.ndarray]:
    """(n_spectra × n_sites) integral populations (%) and their errors, from the
    same integrate-over-the-window quantification as Report (F6).

    A site's integral is proportional to its amplitude for a fixed lineshape, so
    the population's error is the amplitude's *relative* error under the chosen
    ``error_method`` (covariance / Monte-Carlo / χ² profile), applied to the
    fraction — first-order, same approximation as quantify.py's own table
    (the other sites' amplitude errors, which also shift the total, are
    neglected)."""
    from larmor.quantify import quantify
    n = len(result.recipes)
    ns = len(result.recipes[0].sites)
    vals = np.full((n, ns), np.nan)
    errs = np.full((n, ns), np.nan)
    for k, rec in enumerate(result.recipes):
        try:
            q = quantify(rec, getattr(rec, "fit_window_ppm", None))
        except Exception:
            continue
        for i, row in enumerate(q["rows"]):
            if i >= ns:
                continue
            vals[k, i] = row["fraction_pct"]
            amp = rec.sites[i].params.get("amplitude")
            if amp is None or not amp.value:
                continue
            amp_err = _param_error(result, i, "amplitude", k, error_method)
            if np.isfinite(amp_err):
                errs[k, i] = row["fraction_pct"] * abs(amp_err / amp.value)
    return vals, errs


# ----- kept for scripting / CSV export / tests -----------------------------
def series_options(result) -> list[dict]:
    """Every (site, parameter) pair, plus amplitude-fraction and integral
    population per site."""
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
        out.append({"site": i, "param": "population_pct", "kind": "pop_integral",
                    "text": f"s{i} {label}: population % (integral)"})
    return out


def error_methods(result) -> list[str]:
    """Which error-calculation methods this batch result can display: always
    'none' and 'covariance' (from the fit), plus any that were computed
    (Monte-Carlo, χ² profile)."""
    methods = ["none", "covariance"]
    for m in getattr(result, "error_detail", {}) or {}:
        if m not in methods:
            methods.append(m)
    return methods


def _param_error(result, site: int, param: str, k: int, method: str | None):
    """The error of (site, param) for spectrum k under the chosen method:
    a specific computed method (Monte-Carlo / χ² profile / stored covariance),
    or the covariance stderr on the fit; NaN for 'none'."""
    if method in (None, "none"):
        return np.nan
    detail = (getattr(result, "error_detail", {}) or {}).get(method)
    if detail is not None and k < len(detail):
        pe = detail[k].get((site, param))
        if pe is not None and pe.stderr is not None:
            return float(pe.stderr)
        return np.nan
    if method == "covariance":                    # fall back to the fit's stderr
        p = result.recipes[k].sites[site].params.get(param)
        if p is not None and p.stderr is not None:
            return float(p.stderr)
    return np.nan


def series_values(result, opt: dict, error_method: str | None = "covariance"):
    """(values, errors) of one option across every spectrum in the series.

    ``error_method`` selects which computed error to show — 'covariance' (the
    least-squares stderr, default), 'montecarlo', 'profile', or 'none'. Both
    population kinds get an error too (first-order, from the amplitude's error
    under the chosen method — see :func:`population_integral`)."""
    if opt["kind"] == "pop_integral":
        vals, errs = population_integral(result, error_method)
        return np.asarray(vals[:, opt["site"]], float), \
               np.asarray(errs[:, opt["site"]], float)
    vals, errs = [], []
    for k, rec in enumerate(result.recipes):
        site = rec.sites[opt["site"]]
        p = site.params.get(opt["param"])
        v = float(p.value) if p is not None else np.nan
        if opt["kind"] == "popfrac":
            tot = sum(abs(float(s.params["amplitude"].value))
                      for s in rec.sites if "amplitude" in s.params) or 1.0
            frac = 100.0 * abs(v) / tot
            amp_err = _param_error(result, opt["site"], "amplitude", k, error_method)
            e = (frac * abs(amp_err / v)
                 if (np.isfinite(amp_err) and v) else np.nan)
            v = frac
        else:
            e = _param_error(result, opt["site"], opt["param"], k, error_method)
        vals.append(v); errs.append(e)
    return np.array(vals, float), np.array(errs, float)


class SeriesPlotDialog(QDialog):
    def __init__(self, parent, result, series=None):
        super().__init__(parent)
        self.setWindowTitle("Series evolution")
        self.resize(980, 660)
        self._result = result
        self._labels = list(result.labels)
        self._params = _param_specs(result)
        self._pop = None
        # the Series table (larmor.series_table): names, replicate groups and
        # numeric columns for the x axis. None (or a table that does not pair
        # 1:1 with the result) keeps the plain series-order axis.
        self._series = (series if series is not None
                        and len(series.rows) == len(self._labels) else None)
        self._can_avg = bool(self._series is not None and self._series.has_replicates())

        v = QVBoxLayout(self)
        top = QHBoxLayout()
        top.addWidget(QLabel("Pick one or more <b>lines</b> (Ctrl-click for "
                             "several) — each parameter is plotted across the "
                             "series on the right."), 1)
        top.addWidget(QLabel("Error bars:"))
        self.errSel = QComboBox()
        self.errSel.setToolTip("which computed error to draw as error bars "
                               "(and export). Run Monte-Carlo or χ² profile in "
                               "the batch dialog to add those options.")
        self._fill_error_methods()
        self.errSel.currentIndexChanged.connect(self._draw)
        top.addWidget(self.errSel)
        top.addWidget(QLabel("x axis:"))
        self.xSel = QComboBox()
        self.xSel.setToolTip("series order (sample names as ticks) or any numeric "
                             "column of the Series table — its ± becomes x error bars")
        self.xSel.addItem("series order (names)", None)
        for col in (self._series.numeric_columns() if self._series is not None else []):
            self.xSel.addItem(col.label, col.key)
        self.xSel.currentIndexChanged.connect(self._draw)
        top.addWidget(self.xSel)
        self.chkAvg = QCheckBox("average replicates")
        self.chkAvg.setToolTip("collapse each replicate group (the Series table's group "
                               "column) to its mean; the bar is the sample std (ddof = 1) "
                               "and n is exported")
        self.chkAvg.setChecked(self._can_avg)
        self.chkAvg.setVisible(self._can_avg)
        self.chkAvg.toggled.connect(self._draw)
        top.addWidget(self.chkAvg)
        self.chkOls = QCheckBox("fit line")
        self.chkOls.setToolTip("ordinary least-squares line through 3 or more points: "
                               "slope ± error, Pearson r and n in the legend")
        self.chkOls.toggled.connect(self._draw)
        top.addWidget(self.chkOls)
        v.addLayout(top)
        body = QHBoxLayout(); v.addLayout(body, 1)

        self.list = QListWidget()
        self.list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.list.setMaximumWidth(220)
        for i, s in enumerate(result.recipes[0].sites):
            it = QListWidgetItem(f"s{i}  {s.label or s.model}")
            it.setData(Qt.UserRole, i)
            it.setForeground(pg.mkColor(site_color(i)))
            self.list.addItem(it)
        if self.list.count():
            self.list.item(0).setSelected(True)
        self.list.itemSelectionChanged.connect(self._draw)
        body.addWidget(self.list)

        self._scroll = QScrollArea(); self._scroll.setWidgetResizable(True)
        self._grid_host = QWidget(); self._grid = QGridLayout(self._grid_host)
        self._scroll.setWidget(self._grid_host)
        body.addWidget(self._scroll, 1)
        self._subplots: dict = {}
        self._build_subplots()

        btns = QHBoxLayout()
        b_csv = QPushButton("Export parameters (CSV)…"); b_csv.clicked.connect(self._export_csv)
        b_fig = QPushButton("Export figure…"); b_fig.clicked.connect(self._export_fig)
        btns.addWidget(b_csv); btns.addWidget(b_fig); btns.addStretch(1)
        b_bar = QPushButton("Species bar…")
        b_bar.setToolTip("Plotting studio: a 100 %-stacked bar of every site's integral "
                         "population per spectrum (per replicate group when averaging), "
                         "categories in the current x order")
        b_bar.clicked.connect(self._species_bar)
        b_dust = QPushButton("Export DUST CSV…")
        b_dust.setToolTip("Sample + the Series table's oxide columns under DUST's canonical "
                          "names + N4_measured / N4_measured_err (the selected sites' integral "
                          "population fraction) — DUST ignores the N4 columns on import; join "
                          "them against its Results CSV")
        b_dust.clicked.connect(lambda: self._export_dust())
        btns.addWidget(b_bar); btns.addWidget(b_dust)
        v.addLayout(btns)
        self.msg = QLabel("")
        self.msg.setWordWrap(True)
        v.addWidget(self.msg)

        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        v.addWidget(bb)
        self._draw()

    # ------------------------------------------------------------------
    def _fill_error_methods(self):
        labels = {"none": "none", "covariance": "covariance (fit)",
                  "montecarlo": "Monte-Carlo", "profile": "χ² profile"}
        for m in error_methods(self._result):
            self.errSel.addItem(labels.get(m, m), m)
        want = getattr(self._result, "error_method", "covariance") or "covariance"
        idx = self.errSel.findData(want)
        self.errSel.setCurrentIndex(idx if idx >= 0 else self.errSel.findData("covariance"))

    def _error_method(self) -> str:
        return self.errSel.currentData() or "covariance"

    def _selected(self) -> list[int]:
        return [it.data(Qt.UserRole) for it in self.list.selectedItems()]

    def _build_subplots(self):
        from larmor.desktop.plot_menu import attach_plot_menu
        x = np.arange(1, len(self._labels) + 1)
        ticks = [[(int(i), lab) for i, lab in zip(x, self._labels)]]
        for idx, spec in enumerate(self._params):
            pw = pg.PlotWidget(background=theme.active().plot_bg)
            pw.setMinimumHeight(200); pw.setMinimumWidth(320)
            pw.showGrid(x=True, y=True, alpha=0.2)
            pw.addLegend(labelTextSize="8pt")
            pw.setLabel("left", spec["label"])
            pw.getAxis("bottom").setTicks(ticks)
            attach_plot_menu(pw, title=spec["param"], parent=self,
                             studio_spec=lambda s=spec: self._studio_spec_for(s))
            self._grid.addWidget(pw, idx // 2, idx % 2)
            self._subplots[spec["param"]] = pw

    # ------------------------------------------------------------------ x axis / replicates
    def _axis(self) -> dict:
        """The x axis: ``key`` None = series order (names as ticks), else a
        Series-table column (its label; its ± drawn as x error bars)."""
        key = self.xSel.currentData() if hasattr(self, "xSel") else None
        if key is not None and self._series is not None:
            col = self._series.column(key)
            return {"key": key, "label": col.label if col else str(key),
                    "categorical": False}
        return {"key": None, "label": "sample", "categorical": True}

    def _averaging(self) -> bool:
        return bool(self._can_avg and self.chkAvg.isChecked())

    def _points(self, spec: dict, i: int) -> dict:
        """One site's points for one parameter on the current axis:
        {names, x, xerr, y, yerr, n}. Replicate groups are collapsed when
        'average replicates' is on (series_table.replicate_stats: mean, sample
        std, n); on a numeric axis rows without an x are dropped and the rest
        sorted by x."""
        vals, errs = series_values(
            self._result, {"site": i, "param": spec["param"], "kind": spec["kind"]},
            self._error_method())
        ax = self._axis()
        x = xerr = None
        if ax["key"] is not None:
            x, xerr = self._series.x_values(ax["key"])
        if self._averaging():
            stt = stab.replicate_stats(self._series.groups(), vals, errs, x, xerr)
            names, y, yerr, n = stt["labels"], stt["y"], stt["yerr"], stt["n"]
            x, xerr = stt["x"], stt["xerr"]
        else:
            names = list(self._labels)
            y, yerr = np.asarray(vals, float), np.asarray(errs, float)
            n = np.ones(len(names), int)
        if ax["categorical"]:
            x = np.arange(1, len(names) + 1, dtype=float)
            xerr = None
            order = np.arange(len(names))
        else:
            order = np.array([k for k in np.argsort(x, kind="stable")
                              if np.isfinite(x[k])], int)

        def pick(a):
            return None if a is None else np.asarray(a)[order]
        return {"names": [names[k] for k in order], "x": pick(x), "xerr": pick(xerr),
                "y": pick(y), "yerr": pick(yerr), "n": pick(n)}

    def _apply_axis(self, ax: dict, names: list):
        """Ticks = the names on the series-order axis; the column label on a
        numeric one."""
        for pw in self._subplots.values():
            axis = pw.getAxis("bottom")
            if ax["categorical"]:
                axis.setTicks([[(k + 1, lab) for k, lab in enumerate(names)]])
                axis.showLabel(False)
            else:
                axis.setTicks(None)
                pw.setLabel("bottom", ax["label"])

    def _studio_spec_for(self, spec: dict) -> dict:
        """A publication-figure spec for THIS subplot: the selected sites' values
        vs the series -- the sample names as x-ticks on the series-order axis,
        or the Series-table column (data.xerr, and the dashed OLS trace when
        'fit line' is on) on a numeric one."""
        ax = self._axis()
        traces = []
        ref = None
        for i in self._selected() or [0]:
            pts = self._points(spec, i)
            ref = pts
            data = {"x": [float(v) for v in pts["x"]], "y": [float(v) for v in pts["y"]]}
            if np.isfinite(pts["yerr"]).any():
                data["yerr"] = [float(e) if np.isfinite(e) else 0.0 for e in pts["yerr"]]
            if pts["xerr"] is not None and np.isfinite(pts["xerr"]).any():
                data["xerr"] = [float(e) if np.isfinite(e) else 0.0 for e in pts["xerr"]]
            traces.append({"data": data, "label": f"s{i}", "marker": "o",
                           "color": site_color(i), "linestyle": "-"})
            if self.chkOls.isChecked():
                fit = stab.ols(pts["x"], pts["y"])
                if fit["n"] >= 3 and np.isfinite(fit["slope"]):
                    traces.append(stab.ols_trace(pts["x"], fit, f"s{i}", site_color(i)))
        out = {"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
               "xlabel": ax["label"], "ylabel": spec["label"],
               "xtick_rotation": 45, "traces": traces,
               "title": spec["label"]}
        if ax["categorical"] and ref is not None:
            out["xticks"] = [[int(xi), lab] for xi, lab in zip(ref["x"], ref["names"])]
        return out

    def _draw(self):
        sites = self._selected()
        ax = self._axis()
        ref = self._points(self._params[0], sites[0] if sites else 0)
        self._apply_axis(ax, ref["names"])
        for spec in self._params:
            pw = self._subplots[spec["param"]]
            pw.clear()
            try:
                pw.plotItem.legend.clear()
            except Exception:
                pass
            for i in sites:
                pts = self._points(spec, i)
                x, vals, errs, xerr = pts["x"], pts["y"], pts["yerr"], pts["xerr"]
                col = site_color(i)
                pw.plot(x, vals, pen=pg.mkPen(col, width=2), symbol="o",
                        symbolBrush=col, symbolSize=7, name=f"s{i}")
                bars = {}
                if np.isfinite(errs).any():
                    bars["height"] = 2 * np.nan_to_num(errs)
                if xerr is not None and np.isfinite(xerr).any():
                    bars["width"] = 2 * np.nan_to_num(xerr)
                if bars:
                    pw.addItem(pg.ErrorBarItem(x=x, y=vals, pen=col, **bars))
                if self.chkOls.isChecked():
                    fit = stab.ols(x, vals)
                    if fit["n"] >= 3 and np.isfinite(fit["slope"]):
                        xx = np.array([np.nanmin(x), np.nanmax(x)])
                        pw.plot(xx, fit["slope"] * xx + fit["intercept"],
                                pen=pg.mkPen(col, width=1.2, style=Qt.DashLine),
                                name=stab.ols_label(fit))

    # ------------------------------------------------------------------
    def _export_csv(self):
        sites = self._selected()
        if not sites:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export series", "series.csv",
                                              "CSV (*.csv)")
        if not path:
            return
        import csv
        method = self._error_method()
        ax = self._axis()
        ref = self._points(self._params[0], sites[0])
        cols = {"spectrum": list(ref["names"])}
        if not ax["categorical"]:
            cols[f"x ({ax['label']})"] = ref["x"]
            cols["x ±"] = (ref["xerr"] if ref["xerr"] is not None
                           else np.full(len(ref["x"]), np.nan))
        if not ax["categorical"] or self._averaging():
            cols["n"] = ref["n"]
        etag = f"sd; {method} when n=1" if self._averaging() else method
        for i in sites:
            for spec in self._params:
                pts = self._points(spec, i)
                cols[f"s{i}:{spec['param']}"] = pts["y"]
                if np.isfinite(pts["yerr"]).any():
                    cols[f"s{i}:{spec['param']} ±({etag})"] = pts["yerr"]
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(cols.keys())
            for r in range(len(ref["names"])):
                w.writerow([cols[c][r] if isinstance(cols[c][r], str)
                            else f"{cols[c][r]:.6g}" for c in cols])

    def _export_fig(self):
        sites = self._selected()
        if not sites:
            return
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        specs = self._params
        ncol = 2
        nrow = int(np.ceil(len(specs) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(9, 2.6 * nrow), squeeze=False)
        axis = self._axis()
        ref = self._points(specs[0], sites[0])
        for idx, spec in enumerate(specs):
            ax = axes[idx // ncol][idx % ncol]
            for i in sites:
                pts = self._points(spec, i)
                x, vals, errs, xerr = pts["x"], pts["y"], pts["yerr"], pts["xerr"]
                ax.errorbar(x, vals,
                            yerr=np.nan_to_num(errs) if np.isfinite(errs).any() else None,
                            xerr=(np.nan_to_num(xerr) if xerr is not None
                                  and np.isfinite(xerr).any() else None),
                            marker="o", capsize=2, label=f"s{i}", color=site_color(i))
                if self.chkOls.isChecked():
                    fit = stab.ols(x, vals)
                    if fit["n"] >= 3 and np.isfinite(fit["slope"]):
                        xx = np.array([np.nanmin(x), np.nanmax(x)])
                        ax.plot(xx, fit["slope"] * xx + fit["intercept"], "--", lw=1.0,
                                color=site_color(i), label=stab.ols_label(fit))
            ax.set_ylabel(spec["label"])
            if axis["categorical"]:
                ax.set_xticks(ref["x"])
                ax.set_xticklabels(ref["names"], rotation=45, ha="right", fontsize=7)
            else:
                ax.set_xlabel(axis["label"])
            if len(sites) > 1 or self.chkOls.isChecked():
                ax.legend(fontsize=7)
        for j in range(len(specs), nrow * ncol):
            axes[j // ncol][j % ncol].set_visible(False)
        fig.tight_layout()
        from larmor.desktop.export_dialog import export_matplotlib
        export_matplotlib(self, fig, "series")

    # ------------------------------------------------------------------ species bar / DUST
    def _species_bar_spec(self) -> dict:
        """A species_bar spec of EVERY site's integral population (a 100 %-
        stacked bar of a subset would renormalise among the subset), categories
        = the names (or replicate groups) in the current x order."""
        pop = next(s for s in self._params if s["kind"] == "pop_integral")
        sites = list(range(len(self._result.recipes[0].sites)))
        pts = [self._points(pop, j) for j in sites]
        values = np.column_stack([p["y"] for p in pts])
        labels = [f"s{j} {s.label or s.model}"
                  for j, s in enumerate(self._result.recipes[0].sites)]
        ax = self._axis()
        return stab.species_bar_spec(
            pts[0]["names"], labels, values,
            xlabel=None if ax["categorical"] else f"in order of {ax['label']}",
            colors=[site_color(j) for j in sites])

    def _species_bar(self):
        from larmor.desktop.plotting_studio import PlottingStudio
        PlottingStudio(self, self._species_bar_spec()).exec()

    def _dust_rows(self) -> tuple:
        """(header, rows, skipped) for the DUST CSV: N4_measured = the selected
        sites' summed integral population fraction (errors in quadrature)."""
        sites = self._selected()
        vals, errs = population_integral(self._result, self._error_method())
        v, e = vals[:, sites] / 100.0, errs[:, sites] / 100.0
        n4 = np.where(np.isfinite(v).any(axis=1), np.nansum(v, axis=1), np.nan)
        n4e = np.where(np.isfinite(e).any(axis=1), np.sqrt(np.nansum(e ** 2, axis=1)), np.nan)
        return stab.dust_rows(self._series, n4, n4e, self._averaging(), names=self._labels)

    def _export_dust(self, path: str | None = None):
        """Export DUST CSV…: Sample + oxide columns (DUST's canonical names) +
        N4_measured / N4_measured_err; the message names what was written,
        what was skipped and how DUST treats the N4 columns."""
        sites = self._selected()
        if not sites:
            self.msg.setText("select the site(s) that make up N4 (the BO4 line) in the "
                             "list first")
            return
        header, rows, skipped = self._dust_rows()
        if path is None:
            from larmor.desktop.paths import suggest_save_dir
            start = (suggest_save_dir(self._series.rows[0].source_path)
                     if self._series is not None and self._series.rows else "")
            path, _ = QFileDialog.getSaveFileName(
                self, "Export DUST CSV",
                str(Path(start) / "dust_compositions.csv") if start else "dust_compositions.csv",
                "CSV (*.csv)")
        if not path:
            return
        import csv
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
        oxides = header[1:-2]
        self.msg.setText(
            f"wrote {Path(path).name}: {len(rows)} row(s) · oxide columns: "
            + (", ".join(oxides) if oxides
               else "none — join a composition CSV in the Series table first")
            + " · N4_measured = " + " + ".join(f"s{i}" for i in sites)
            + " integral population / 100"
            + (" (replicate means; N4_measured_err = spread)" if self._averaging() else "")
            + (f" · skipped non-oxide columns: {', '.join(skipped)}" if skipped else "")
            + " · DUST's import dialog maps N4_measured to Ignore — join it against "
              "DUST's Results CSV by Sample")
