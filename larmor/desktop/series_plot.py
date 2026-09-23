"""Series evolution plot: how fitted parameters change along a batch series.

Pick one or more **sites** (lines) on the left; the right shows **one subplot per
parameter** (δ_iso, width, C_Q, η, amplitude, and the integral **population %**),
each tracing the chosen sites across the series. Export the numbers as CSV or the
whole panel as a figure.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDialog, QDialogButtonBox, QFileDialog,
    QGridLayout, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton,
    QScrollArea, QVBoxLayout, QWidget,
)

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
    # a ratio subplot only when the master recipe's families define one (N3)
    if _group_names(result)[1]:
        specs.append({"param": "ratio", "kind": "ratio", "label": "named ratio"})
    return specs


def _group_names(result) -> tuple[list[str], list[str]]:
    """(family names, defined ratio names) of the master recipe's tags
    (larmor.families) -- ([], []) when nothing is tagged or quantify fails."""
    from larmor.quantify import quantify
    master = result.recipes[0]
    try:
        q = quantify(master, getattr(master, "fit_window_ppm", None))
    except Exception:
        return [], []
    fams = [f["family"] for f in q.get("families") or []]
    ratios = [r["name"] for r in q.get("ratios") or []
              if r.get("defined") and r.get("value") is not None]
    return fams, ratios


def _stored_group_error(result, k: int, key: str, method) -> float:
    """The exact per-spectrum family/ratio error an error run stored under
    ``(-1, key)`` ('family:NAME' / 'ratio:NAME'), else NaN. Never reads the
    recipe (site -1 is not a site), unlike _param_error."""
    if method in (None, "none"):
        return np.nan
    detail = (getattr(result, "error_detail", {}) or {}).get(method)
    if detail is not None and k < len(detail):
        pe = detail[k].get((-1, key))
        if pe is not None and pe.stderr is not None:
            return float(pe.stderr)
    return np.nan


def _group_series(result, opt: dict, method) -> tuple[np.ndarray, np.ndarray]:
    """(values, errors) across the series of one family's summed % or one
    named ratio: quantify() per spectrum, the error the stored exact value
    (Monte-Carlo / covariance run) or quantify's own flagged fallback; NaN
    for 'none', for a spectrum whose every member is excluded, or when the
    ratio is undefined there."""
    from larmor.batchfit import is_zeroed_out
    from larmor.quantify import quantify
    kind = opt["kind"]
    name = opt["family"] if kind == "family" else opt["ratio"]
    vals, errs = [], []
    for k, rec in enumerate(result.recipes):
        v = e = np.nan
        try:
            q = quantify(rec, getattr(rec, "fit_window_ppm", None))
        except Exception:
            q = None
        row, fallback = None, None
        if q is not None and kind == "family":
            row = next((f for f in q["families"] if f["family"] == name), None)
            if row is not None and all(
                    is_zeroed_out(rec.sites[i].params.get("amplitude"))
                    for i in row["sites"] if i < len(rec.sites)):
                row = None                         # every member excluded here
            if row is not None:
                v, fallback = row["fraction_pct"], row["fraction_err_pct"]
        elif q is not None:
            row = next((r for r in q["ratios"] if r["name"] == name
                        and r.get("defined") and r.get("value") is not None), None)
            if row is not None:
                v, fallback = row["value"], row["err"]
        if row is not None and method not in (None, "none"):
            stored = _stored_group_error(result, k, f"{kind}:{name}", method)
            e = stored if np.isfinite(stored) else (
                float(fallback) if fallback is not None else np.nan)
        vals.append(v)
        errs.append(e)
    return np.array(vals, float), np.array(errs, float)


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
    fams, ratios = _group_names(result)
    for name in fams:
        out.append({"kind": "family", "family": name, "param": "family_pct",
                    "site": None, "text": f"Σ {name}: population % (integral)"})
    for name in ratios:
        out.append({"kind": "ratio", "ratio": name, "param": "ratio",
                    "site": None, "text": f"{name}: named ratio"})
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
    under the chosen method — see :func:`population_integral`). The 'family'
    and 'ratio' kinds (N3: ``opt['family']`` / ``opt['ratio']``) take the
    exact per-spectrum error an error run stored, else quantify's own."""
    if opt["kind"] in ("family", "ratio"):
        return _group_series(result, opt, error_method)
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
    def __init__(self, parent, result):
        super().__init__(parent)
        self.setWindowTitle("Series evolution")
        self.resize(980, 660)
        self._result = result
        self._labels = list(result.labels)
        self._params = _param_specs(result)
        self._pop = None

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
        # Σ families feed the population subplot, ratios the ratio subplot
        # (N3); UserRole 'f:NAME' / 'r:NAME' beside the sites' ints
        fams, ratios = _group_names(result)
        self._group_index = {}
        for j, name in enumerate(fams):
            it = QListWidgetItem(f"Σ {name}")
            it.setData(Qt.UserRole, f"f:{name}")
            it.setForeground(pg.mkColor(self._sel_color(f"f:{name}", j)))
            it.setToolTip("summed population % of the lines tagged "
                          f"{name} (Report F6 families)")
            self.list.addItem(it)
            self._group_index[f"f:{name}"] = j
        for j, name in enumerate(ratios):
            it = QListWidgetItem(name)
            it.setData(Qt.UserRole, f"r:{name}")
            it.setForeground(pg.mkColor(self._sel_color(f"r:{name}", j)))
            it.setToolTip(f"named ratio {name} of the tagged families")
            self.list.addItem(it)
            self._group_index[f"r:{name}"] = j
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
        v.addLayout(btns)

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

    def _selected(self) -> list:
        """Site indices (int) and group keys ('f:NAME' / 'r:NAME', str)."""
        return [it.data(Qt.UserRole) for it in self.list.selectedItems()]

    # -- one selection x one subplot -> a series option (or None) ---------
    @staticmethod
    def _opt_for(sel, spec: dict) -> dict | None:
        """Sites feed every parameter subplot but the ratio one; a Σ family
        feeds only the population subplot; a ratio only the ratio subplot."""
        if isinstance(sel, str):
            kind, name = sel.split(":", 1)
            if kind == "f" and spec["kind"] == "pop_integral":
                return {"kind": "family", "family": name, "param": "family_pct",
                        "site": None}
            if kind == "r" and spec["kind"] == "ratio":
                return {"kind": "ratio", "ratio": name, "param": "ratio",
                        "site": None}
            return None
        if spec["kind"] == "ratio":
            return None
        return {"site": sel, "param": spec["param"], "kind": spec["kind"]}

    @staticmethod
    def _sel_label(sel) -> str:
        if isinstance(sel, str):
            kind, name = sel.split(":", 1)
            return f"Σ {name}" if kind == "f" else name
        return f"s{sel}"

    def _sel_color(self, sel, j: int | None = None) -> str:
        if isinstance(sel, str):
            if j is None:
                j = getattr(self, "_group_index", {}).get(sel, 0)
            # hues away from the site palette; families warm, ratios cool
            return pg.intColor(j, hues=6, values=2, maxValue=200,
                               minValue=90, sat=180,
                               alpha=255).name() if sel.startswith("f:") else \
                pg.intColor(j + 3, hues=6, values=2, maxValue=160, minValue=60,
                            sat=200, alpha=255).name()
        return site_color(sel)

    @staticmethod
    def _sel_key(sel, spec: dict) -> str:
        """CSV column: 's0:amplitude', 'family:BO4', 'ratio:N4'."""
        if isinstance(sel, str):
            kind, name = sel.split(":", 1)
            return ("family:" if kind == "f" else "ratio:") + name
        return f"s{sel}:{spec['param']}"

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

    def _studio_spec_for(self, spec: dict) -> dict:
        """A publication-figure spec for THIS subplot: the selected sites' values
        vs the series, on an upright axis with the sample names as x-ticks."""
        x = list(range(1, len(self._labels) + 1))
        traces = []
        for sel in self._selected() or [0]:
            opt = self._opt_for(sel, spec)
            if opt is None:
                continue
            vals, errs = series_values(self._result, opt, self._error_method())
            data = {"x": x, "y": [float(v) for v in vals]}
            if np.isfinite(errs).any():
                data["yerr"] = [float(e) if np.isfinite(e) else 0.0 for e in errs]
            traces.append({"data": data, "label": self._sel_label(sel), "marker": "o",
                           "color": self._sel_color(sel), "linestyle": "-"})
        return {"kind": "1d", "x_is_ppm": False, "hide_yaxis": False,
                "xlabel": "sample", "ylabel": spec["label"],
                "xticks": [[xi, lab] for xi, lab in zip(x, self._labels)],
                "xtick_rotation": 45, "traces": traces,
                "title": spec["label"]}

    def _draw(self):
        sites = self._selected()
        x = np.arange(1, len(self._labels) + 1)
        for spec in self._params:
            pw = self._subplots[spec["param"]]
            pw.clear()
            try:
                pw.plotItem.legend.clear()
            except Exception:
                pass
            for sel in sites:
                opt = self._opt_for(sel, spec)
                if opt is None:
                    continue
                vals, errs = series_values(self._result, opt, self._error_method())
                col = self._sel_color(sel)
                pw.plot(x, vals, pen=pg.mkPen(col, width=2), symbol="o",
                        symbolBrush=col, symbolSize=7, name=self._sel_label(sel))
                if np.isfinite(errs).any():
                    pw.addItem(pg.ErrorBarItem(x=x, y=vals,
                               height=2 * np.nan_to_num(errs), pen=col))

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
        cols = {"spectrum": self._labels}
        for sel in sites:
            for spec in self._params:
                opt = self._opt_for(sel, spec)
                if opt is None:
                    continue
                vals, errs = series_values(self._result, opt, method)
                key = self._sel_key(sel, spec)
                cols[key] = vals
                if np.isfinite(errs).any():
                    cols[f"{key} ±({method})"] = errs
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(cols.keys())
            for r in range(len(self._labels)):
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
        x = np.arange(1, len(self._labels) + 1)
        for idx, spec in enumerate(specs):
            ax = axes[idx // ncol][idx % ncol]
            for sel in sites:
                opt = self._opt_for(sel, spec)
                if opt is None:
                    continue
                vals, errs = series_values(self._result, opt, self._error_method())
                ax.errorbar(x, vals, yerr=np.nan_to_num(errs) if np.isfinite(errs).any()
                            else None, marker="o", capsize=2,
                            label=self._sel_label(sel), color=self._sel_color(sel))
            ax.set_ylabel(spec["label"]); ax.set_xticks(x)
            ax.set_xticklabels(self._labels, rotation=45, ha="right", fontsize=7)
            if len(sites) > 1:
                ax.legend(fontsize=7)
        for j in range(len(specs), nrow * ncol):
            axes[j // ncol][j % ncol].set_visible(False)
        fig.tight_layout()
        from larmor.desktop.export_dialog import export_matplotlib
        export_matplotlib(self, fig, "series")
