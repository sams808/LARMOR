"""Two-field (or multi-field) infinite-field extrapolation of the isotropic
chemical shift from the central-transition centre of gravity — for QCPMG data
of half-integer quadrupolar nuclei measured at more than one magnetic field
(Sandland et al. 2004; Baasner et al. 2014). See larmor.qcpmg_fields."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from larmor.desktop import theme
from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from larmor.qcpmg_fields import (
    ERR_FLOOR_PPM, FieldPoint, field_plausibility_warning, fmt_result_lines,
    infinite_field_diso, read_field_spectrum,
)


#: window-mode combo labels and the qcpmg.measure_cg modes behind them
WINDOW_MODES = ("first minima / dragged band", "whole manifold",
                "centreband (< ν_r)")
WINDOW_MODE_KEYS = ("manual", "manifold", "centreband")


def _spin_of(nucleus: str) -> float:
    try:
        from mrsimulator.spin_system.isotope import ISOTOPE_DATA
        d = ISOTOPE_DATA.get(nucleus or "")
        if d:
            return (d["spin_multiplicity"] - 1) / 2.0
    except Exception:
        pass
    return 1.5


class QcpmgFieldsDialog(QDialog):
    def __init__(self, parent, nucleus: str = "", current=None):
        super().__init__(parent)
        self.setWindowTitle("QCPMG — infinite-field δiso (2+ fields)")
        self.resize(860, 720)
        self.setMinimumSize(560, 420)
        self.setSizeGripEnabled(True)
        # this window is meant to stay open across several processing
        # sessions, so it needs to be minimisable and never modal
        self.setWindowFlags(self.windowFlags() | Qt.Window
                            | Qt.WindowMinimizeButtonHint
                            | Qt.WindowMaximizeButtonHint)
        self.setModal(False)
        # an empty nucleus stays empty (spin defaults to 3/2, the label says
        # so) rather than being silently relabelled 35Cl
        self._nucleus = nucleus or ""
        self._current = current            # (larmor_MHz, ppm, amp) of the open spectrum
        v = QVBoxLayout(self)

        intro = QLabel(
            "Extrapolate the central-transition centre of gravity to infinite "
            "field to remove the second-order quadrupolar shift and obtain the "
            "true isotropic shift δiso and C_Q (Sandland 2004 Eq. 1 / Baasner "
            "2014 Fig. 6). Enter δcg at each field, or grab it from the open "
            "spectrum. In the large-C_Q limit CT-selective and non-selective "
            "fields can be combined (see the Lineshapes/QCPMG manual).")
        intro.setWordWrap(True)
        intro.setStyleSheet(f"color:{theme.active().text_dim};")
        v.addWidget(intro)

        top = QHBoxLayout()
        self.lblNuc = QLabel()
        self._set_nucleus(self._nucleus, spin=False)
        top.addWidget(self.lblNuc)
        top.addWidget(QLabel("  ·  spin I ="))
        self.spin = QDoubleSpinBox(); self.spin.setDecimals(1)
        self.spin.setRange(1.5, 4.5); self.spin.setSingleStep(1.0)
        self.spin.setValue(_spin_of(self._nucleus))
        top.addWidget(self.spin)
        top.addSpacing(16)
        top.addWidget(QLabel("η (assumed)"))
        self.eta = QDoubleSpinBox(); self.eta.setRange(0.0, 1.0)
        self.eta.setSingleStep(0.05)
        from larmor.qcpmg_fields import DEFAULT_ETA
        self.eta.setValue(DEFAULT_ETA)
        self.eta.setToolTip("η is not determined by two centres of gravity; "
                            "0.7 is the conventional choice (Stebbins & Du 2002)")
        top.addWidget(self.eta)
        top.addStretch(1)
        v.addLayout(top)

        # per-field table: Larmor (MHz), δcg (ppm), ±err, FWHM (ppm), CT-selective
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["Larmor ν₀ (MHz)", "δcg (ppm)", "± err (ppm)",
             "FWHM (ppm)", "CT-selective (declared)"])
        self.table.horizontalHeader().setStretchLastSection(True)
        # a QTableWidget's own minimum is ~1.5 rows, so at the default dialog
        # size the plot's fixed minimum squeezed the table nearly out of view
        self.table.setMinimumHeight(150)             # header + ~4 field rows
        v.addWidget(self.table, 1)

        row = QHBoxLayout()
        b_add = QPushButton("＋ Add field")
        b_add.clicked.connect(lambda: self._add_row())
        b_del = QPushButton("Remove selected")
        b_del.clicked.connect(self._del_row)
        b_ds = QPushButton("Add from datasets…")
        b_ds.setToolTip("pick the sum-echo datasets (.csv from Save as "
                        "dataset…) or processed 1r spectra measured at each "
                        "field: δcg ± σ and FWHM are read off automatically, "
                        "and selecting the row shows the spectrum with a "
                        "draggable band to supervise the values. A QCPMG 1r "
                        "is a spikelet comb: its window is seeded from the "
                        "envelope and the row is flagged")
        b_ds.clicked.connect(self._pick_datasets)
        self.b_cur = QPushButton("δcg from open spectrum (visible range)")
        self.b_cur.setToolTip("centre of gravity of the currently open spectrum "
                              "over the visible x-range — zoom to the CT band first")
        self.b_cur.clicked.connect(self._from_current)
        self.b_cur.setEnabled(self._current is not None)
        row.addWidget(b_add); row.addWidget(b_del); row.addWidget(b_ds)
        row.addWidget(self.b_cur)
        row.addStretch(1)
        row.addWidget(QLabel("window"))
        self.winMode = QComboBox()
        self.winMode.addItems(WINDOW_MODES)
        self.winMode.setToolTip(
            "how the selected dataset row's δcg window is defined:\n"
            "• first minima / dragged band — the default, supervise it\n"
            "• whole manifold — integrate the full axis (edge floor "
            "subtracted); exact for a distribution of sites under MAS "
            "provided the whole sideband manifold is in the spectrum\n"
            "• centreband — peak ± ν_r/2, valid only when it reproduces the "
            "whole-manifold CG (a pattern narrower than ν_r)")
        self.winMode.currentIndexChanged.connect(self._mode_changed)
        row.addWidget(self.winMode)
        v.addLayout(row)

        self.plot = pg.PlotWidget(background=theme.active().plot_bg)
        self.plot.setLabel("bottom", "1 / ν₀²", units="MHz⁻²")
        self.plot.setLabel("left", "δcg", units="ppm")
        self.plot.getPlotItem().getAxis("left").enableAutoSIPrefix(False)
        self.plot.showGrid(x=True, y=True, alpha=0.15)
        self.plot.setMinimumHeight(140)
        v.addWidget(self.plot, 1)        # table and plot share extra height

        self.result = QLabel("add at least two fields, then Compute")
        self.result.setStyleSheet(f"font-weight:600; color:{theme.active().accent};")
        self.result.setWordWrap(True)
        v.addWidget(self.result)

        self.wresult = QLabel(
            "Two-field linewidth split (Sandland Eq. 2): also fill FWHM (ppm) "
            "at both fields, then 'Split W_q / W_csd'.")
        self.wresult.setStyleSheet(f"color:{theme.active().text_dim};")
        self.wresult.setWordWrap(True)
        self.wresult.setTextFormat(Qt.RichText)
        v.addWidget(self.wresult)

        bb = QDialogButtonBox()
        b_comp = bb.addButton("Compute δiso", QDialogButtonBox.ApplyRole)
        b_comp.clicked.connect(self._compute)
        b_w = bb.addButton("Split W_q / W_csd", QDialogButtonBox.ApplyRole)
        b_w.clicked.connect(self._compute_widths)
        self.btnReport = bb.addButton("Export report…",
                                      QDialogButtonBox.ActionRole)
        self.btnReport.setToolTip("every input point, every fitted number "
                                  "with its uncertainty, and the assumptions")
        self.btnReport.clicked.connect(self._export_report)
        self.btnFig = bb.addButton("Export figure…", QDialogButtonBox.ActionRole)
        self.btnFig.setToolTip("the extrapolation as a publication figure "
                               "(.png + .svg + .pdf, 600 dpi)")
        self.btnFig.clicked.connect(self._export_figure)
        bb.addButton(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.addButton(QDialogButtonBox.Help).clicked.connect(self._help)
        v.addWidget(bb)

        # dataset supervision: rows added from a dataset keep their spectrum,
        # and selecting one shows it with a draggable band
        self._ds: dict[int, dict] = {}
        self._ds_seq = 0
        self._region = None
        self._cg_line = None
        self.table.itemSelectionChanged.connect(self._show_selected_dataset)

        # seed two rows; prefill the first from the open spectrum's field
        self._add_row(self._current[0] if self._current else 0.0)
        self._add_row()

    # ------------------------------------------------- datasets (supervised)
    def _pick_datasets(self):
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Choose the processed spectrum measured at each field",
            "", "LARMOR spectrum (*.csv *.txt);;Bruker processed (1r);;"
                "All files (*)")
        errors = []
        for p in paths:
            try:
                self.add_dataset_file(p)
            except Exception as exc:                          # noqa: BLE001
                errors.append(f"{p}: {exc}")
        if errors:
            QMessageBox.warning(self, "Some datasets were skipped",
                                "\n".join(errors))

    def add_dataset_file(self, path: str) -> int:
        """One field from a file: the sum-echo .csv written by *Save as
        dataset...* (its header carries the Larmor frequency, nucleus and
        processing mode) or a Bruker 1r (mode from procs PH_mod; a QCPMG 1r
        is a spikelet comb and is seeded from its envelope, see
        :func:`larmor.qcpmg.seed_window`)."""
        fs = read_field_spectrum(path)
        r = self.add_dataset_spectrum(
            fs["larmor"], fs["ppm"], fs["amp"], magnitude=fs["magnitude"],
            nucleus=fs["nucleus"], source=fs["source"], meta=fs["meta"],
            rotor_Hz=fs["rotor_Hz"])
        if r >= 0 and fs.get("rotor_note"):
            self.wresult.setText(
                f"<span style='color:#c0392b'>⚠ row {r + 1}: {fs['rotor_note']}</span>")
        return r

    def _set_nucleus(self, nucleus: str, spin: bool = True):
        self._nucleus = nucleus or ""
        self.lblNuc.setText(f"nucleus <b>{self._nucleus or '—'}</b>")
        if spin:
            self.spin.setValue(_spin_of(self._nucleus))

    def has_data_rows(self) -> bool:
        """True once any row holds a δcg (typed or from a dataset); a seeded
        field alone does not count."""
        for r in range(self.table.rowCount()):
            it = self.table.item(r, 1)
            if (it is not None and it.text().strip()) or self._row_ds_id(r) is not None:
                return True
        return False

    def accept_nucleus(self, nucleus: str) -> bool:
        """Reconcile an incoming spectrum's nucleus with the dialog's. An
        empty or equal nucleus is fine; a different one is adopted while the
        table is still empty (label + spin follow) and REFUSED once rows are
        present -- a 35Cl row fitted with the 27Al spin would double C_Q
        with no sign of it. Returns False on refusal (a warning is shown)."""
        nucleus = (nucleus or "").strip()
        if not nucleus or nucleus == self._nucleus:
            return True
        if not self._nucleus or not self.has_data_rows():
            self._set_nucleus(nucleus)
            return True
        self.wresult.setText(
            f"<span style='color:#c0392b'>⚠ this dialog holds {self._nucleus} "
            f"rows; a {nucleus} spectrum was not added -- remove the rows "
            "(or open a new dialog) before extrapolating another "
            "nucleus.</span>")
        return False

    def add_dataset_spectrum(self, larmor_MHz: float, ppm, amp,
                             window=None, magnitude: bool | None = False,
                             nucleus: str = "", source: str = "",
                             meta: dict | None = None,
                             rotor_Hz: float = 0.0) -> int:
        """Add one field's spectrum: δcg ± σ and FWHM are computed over the
        given window (or a seeded one) and written into a new row; the
        spectrum stays attached so selecting the row shows it for
        supervision. ``magnitude`` None = mode unknown. ``source`` names
        where the spectrum came from (1r path, sum-echo, workspace) for the
        row tooltip and the report. Returns the row, or -1 when the
        spectrum's ``nucleus`` is not this dialog's (:meth:`accept_nucleus`)."""
        from larmor import qcpmg
        if not self.accept_nucleus(nucleus):
            return -1
        ppm = np.asarray(ppm, float); amp = np.asarray(amp, float)
        seed = None
        if window is not None:
            lo, hi = float(min(window)), float(max(window))
        else:
            seed = qcpmg.seed_window(ppm, amp, meta)
            hi, lo = seed.hi_ppm, seed.lo_ppm
        if not (np.isfinite(hi) and np.isfinite(lo)) or hi <= lo:
            span = float(ppm.max() - ppm.min())
            mid = float(ppm.min()) + span / 2.0
            lo, hi = mid - span / 6.0, mid + span / 6.0
        self._ds_seq += 1
        ds_id = self._ds_seq
        self._ds[ds_id] = {"ppm": ppm, "amp": amp, "window": (lo, hi),
                           "magnitude": None if magnitude is None else bool(magnitude),
                           "source": source or "",
                           "comb": bool(seed is not None and seed.comb),
                           "seed_note": seed.note if seed is not None else "",
                           "larmor": float(larmor_MHz or 0.0),
                           "rotor_Hz": float(rotor_Hz or 0.0),
                           "mode": "manual" if window is not None else "minima"}
        self._add_row(larmor_MHz)
        r = self.table.rowCount() - 1
        self.table.item(r, 0).setData(Qt.UserRole, ds_id)
        mode = ("MAGNITUDE (mc)" if magnitude else
                "absorption" if magnitude is not None else
                "mode unknown -- taken from the workspace")
        tip = f"δcg measured on the {mode} spectrum"
        if source:
            tip += f"\nsource: {source}"
        if rotor_Hz:
            tip += f"\nMAS {float(rotor_Hz):.0f} Hz = {float(rotor_Hz) / float(larmor_MHz):.0f} ppm"
        if seed is not None and seed.comb:
            tip += f"\n⚠ {seed.note}"
            self.wresult.setText(
                f"<span style='color:#c0392b'>⚠ row {r + 1}: {seed.note}</span>")
        self._ds[ds_id]["tip"] = tip
        pp = str((meta or {}).get("pulse_program", "") or "")
        if pp:
            w = self.table.cellWidget(r, 4)
            chk = w.findChild(QCheckBox) if w is not None else None
            if chk is not None:
                chk.setToolTip(f"pulse program: {pp} -- selectivity is still "
                               "your declaration (nu_rf vs nu_Q)")
        self._apply_ds_values(r, ds_id)
        self._warn_mixed_modes()
        self.table.selectRow(r)
        return r

    def _warn_mixed_modes(self):
        """Sandland's extrapolation compares centres of gravity ACROSS
        fields; one measured in magnitude and one in absorption are not the
        same observable, so say so rather than fitting them together
        silently."""
        modes = {bool(d.get("magnitude")) for d in self._ds.values()
                 if d.get("magnitude") is not None}
        if len(modes) > 1:
            self.wresult.setText(
                "<span style='color:#c0392b'>⚠ these fields mix magnitude "
                "(mc) and absorption δcg values — they are different "
                "observables; reprocess them the same way before "
                "extrapolating.</span>")

    def _apply_ds_values(self, r: int, ds_id: int):
        from larmor import qcpmg
        d = self._ds[ds_id]
        mode = d.get("mode", "minima")
        try:
            m = qcpmg.measure_cg(
                d["ppm"], d["amp"],
                window=None if mode == "manifold" else tuple(d["window"]),
                mode="manual" if mode in ("minima", "manual") else mode,
                rotor_Hz=d.get("rotor_Hz", 0.0), larmor_MHz=d.get("larmor", 0.0),
                magnitude=d.get("magnitude"))
        except ValueError as exc:                     # centreband refused
            self.wresult.setText(f"<span style='color:#c0392b'>⚠ {exc}</span>")
            d["mode"] = "manual"
            return
        m.mode = mode
        d["meas"] = m
        d["window"] = tuple(m.window)
        cg, sigma, fw_ppm = m.cg_ppm, m.sigma_ppm, m.fwhm_ppm
        if np.isfinite(cg):
            self.table.setItem(r, 1, QTableWidgetItem(f"{cg:.2f}"))
            self.table.setItem(r, 2, QTableWidgetItem(f"{max(sigma, ERR_FLOOR_PPM):.1f}"))
        self.table.setItem(r, 3, QTableWidgetItem(f"{fw_ppm:.2f}"))
        tip = d.get("tip", "")
        tip = "\n".join(bit for bit in [tip, f"window {m.window[0]:.1f} … "
                                              f"{m.window[1]:.1f} ppm ({mode})",
                                        (m.convergence.sequence()
                                         if m.convergence is not None else ""),
                                        *m.flags] if bit)
        for c in range(4):
            it = self.table.item(r, c)
            if it is not None:
                it.setToolTip(tip)
                if c == 2:
                    it.setForeground(QBrush(QColor("#c0392b")) if m.flags
                                     else QBrush())
        if m.flags:
            self.wresult.setText(
                f"<span style='color:#c0392b'>row {r + 1}: "
                + "  ·  ".join(m.flags) + "</span>")
        if self._cg_line is not None and np.isfinite(cg):
            self._cg_line.setValue(cg)
        if self._region is not None and self._row_ds_id(self.table.currentRow()) == ds_id:
            self._region.blockSignals(True)
            self._region.setRegion(m.window)
            self._region.blockSignals(False)
        self._draw_ticks(m)

    def _draw_ticks(self, m):
        """Dotted lines at the tallest peak ± k·ν_r/ν0 on the supervision
        plot, so a window that catches one sideband is visible."""
        for ln in getattr(self, "_tick_lines", []):
            try:
                self.plot.removeItem(ln)
            except Exception:                                 # noqa: BLE001
                pass
        self._tick_lines = []
        if self._region is None or not m.ticks_ppm:
            return
        pen = pg.mkPen(theme.active().text_dim, style=Qt.DotLine)
        for t in m.ticks_ppm:
            ln = pg.InfiniteLine(pos=t, angle=90, movable=False, pen=pen)
            self.plot.addItem(ln)
            self._tick_lines.append(ln)

    def _mode_changed(self, idx: int):
        r = self.table.currentRow()
        ds_id = self._row_ds_id(r) if r >= 0 else None
        if ds_id is None or ds_id not in self._ds:
            return
        self._ds[ds_id]["mode"] = WINDOW_MODE_KEYS[idx]
        self._apply_ds_values(r, ds_id)

    def _row_ds_id(self, r: int):
        it = self.table.item(r, 0)
        return it.data(Qt.UserRole) if it is not None else None

    def _show_selected_dataset(self):
        r = self.table.currentRow()
        ds_id = self._row_ds_id(r) if r >= 0 else None
        if ds_id is None or ds_id not in self._ds:
            return
        d = self._ds[ds_id]
        self.plot.clear()
        self._region = self._cg_line = None
        self._tick_lines = []
        self.winMode.blockSignals(True)
        self.winMode.setCurrentIndex(
            WINDOW_MODE_KEYS.index(d.get("mode", "minima"))
            if d.get("mode", "minima") in WINDOW_MODE_KEYS else 0)
        self.winMode.blockSignals(False)
        self.plot.getPlotItem().invertX(True)
        self.plot.setLabel("bottom", "shift", units="ppm")
        self.plot.setLabel("left", "intensity", units="")
        self.plot.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)
        self.plot.setTitle("dataset — drag the band edges; δcg / FWHM in the "
                           "row follow", color=theme.active().text_dim,
                           size="9pt")
        self.plot.plot(d["ppm"], d["amp"],
                       pen=pg.mkPen(theme.active().experiment, width=1.2))
        lo, hi = d["window"]
        self._region = pg.LinearRegionItem(values=(lo, hi), movable=True)
        self.plot.addItem(self._region)
        self._region.sigRegionChangeFinished.connect(
            lambda *_: self._region_moved(r, ds_id))
        self._cg_line = pg.InfiniteLine(pos=0.0, angle=90, movable=False,
                                        pen=pg.mkPen(theme.active().pivot,
                                                     style=Qt.DashLine))
        self.plot.addItem(self._cg_line)
        self._apply_ds_values(r, ds_id)

    def _region_moved(self, r: int, ds_id: int):
        if self._region is None or ds_id not in self._ds:
            return
        a, b = self._region.getRegion()
        self._ds[ds_id]["window"] = (min(a, b), max(a, b))
        self._ds[ds_id]["mode"] = "manual"          # the user placed it
        self.winMode.blockSignals(True)
        self.winMode.setCurrentIndex(0)
        self.winMode.blockSignals(False)
        if r < self.table.rowCount() and self._row_ds_id(r) == ds_id:
            self._apply_ds_values(r, ds_id)

    # --------------------------------------------------------- exports
    def _result_map(self):
        """{sample: InfiniteFieldResult} from the current table, so the report
        and the figure describe exactly what Compute just showed."""
        from larmor.qcpmg_fields import fit_samples
        pts = self._points()
        if len(pts) < 2:
            return {}
        name = self._nucleus or "sample"
        return fit_samples([(name, p) for p in pts],
                           spin=self.spin.value(), eta=self.eta.value())

    def _figure_spec(self) -> dict:
        from larmor.qcpmg_fields import InfiniteFieldResult
        samples = [{"label": k,
                    "points": [[p.larmor_MHz, p.dcg_ppm, p.dcg_err_ppm]
                               for p in r.points]}
                   for k, r in self._result_map().items()
                   if isinstance(r, InfiniteFieldResult)]
        return {"kind": "infinite_field", "style": "article",
                "nucleus": self._nucleus, "spin": self.spin.value(),
                "eta": self.eta.value(), "samples": samples}

    def _export_report(self):
        from PySide6.QtWidgets import QFileDialog
        from larmor.desktop.paths import (FIGURE_DIR_KEY, remember_dir,
                                          remembered_dir)
        from larmor.qcpmg_fields import report_text, two_field_widths
        results = self._result_map()
        if not results:
            self.result.setText("need δcg at ≥ 2 fields before a report")
            return
        widths = {}
        fw = self._fields_fwhm()
        if len(fw) >= 2:
            widths = {k: two_field_widths(fw[0][0], fw[0][1], fw[1][0], fw[1][1])
                      for k in results}
        start = remembered_dir(FIGURE_DIR_KEY)
        seed = str(Path(start) / "infinite_field_report.txt") if start             else "infinite_field_report.txt"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export report", seed, "Text (*.txt);;All files (*)")
        if not path:
            return
        remember_dir(FIGURE_DIR_KEY, path)
        Path(path).write_text(
            report_text(results, self.spin.value(), self.eta.value(),
                        self._nucleus, widths), encoding="utf-8")
        self.result.setText(f"report written — {Path(path).name}")

    def _export_figure(self):
        from PySide6.QtWidgets import QFileDialog
        from larmor import figures
        from larmor.desktop.paths import (FIGURE_DIR_KEY, remember_dir,
                                          remembered_dir)
        spec = self._figure_spec()
        if not spec["samples"]:
            self.result.setText("need δcg at ≥ 2 fields before a figure")
            return
        start = remembered_dir(FIGURE_DIR_KEY)
        seed = str(Path(start) / "infinite_field") if start else "infinite_field"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export figure (base name — extensions added)", seed,
            "Figure base name (*)")
        if not path:
            return
        remember_dir(FIGURE_DIR_KEY, path)
        try:
            written = figures.export(spec, Path(path).with_suffix(""))
        except Exception as exc:                              # noqa: BLE001
            self.result.setText(f"figure export failed: {exc}")
            return
        self.result.setText(f"figure written — {Path(written[0]).stem}.png/.svg/.pdf")

    def _help(self):
        from larmor.desktop.help_dialog import show_help
        show_help(self, "qcpmg", "QCPMG")

    def _add_row(self, larmor=0.0):
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, 0, QTableWidgetItem(f"{larmor:g}" if larmor else ""))
        self.table.setItem(r, 1, QTableWidgetItem(""))
        # the +- cell starts EMPTY: a missing sigma is "no uncertainty known"
        # and the fitter then says so, rather than a default 5 ppm silently
        # setting the weights and the reported +-
        err_item = QTableWidgetItem("")
        err_item.setToolTip("δcg uncertainty (ppm); leave empty when unknown "
                            "-- the fit is then unweighted and says so")
        self.table.setItem(r, 2, err_item)
        self.table.setItem(r, 3, QTableWidgetItem(""))          # FWHM (ppm)
        # tri-state, starting UNKNOWN: selectivity is the operator's judgement
        # (nu_rf vs nu_Q), not a property of the data or the pulse program,
        # so a prefilled 'yes' would be fabricated provenance
        chk = QCheckBox(); chk.setTristate(True)
        chk.setCheckState(Qt.PartiallyChecked)
        chk.setToolTip("was this field acquired with a CT-selective pulse? "
                       "declared by you, recorded for provenance, not used by "
                       "the fit -- click through unknown / yes / no")
        w = QWidget(); lay = QHBoxLayout(w); lay.setContentsMargins(0, 0, 0, 0)
        lay.setAlignment(Qt.AlignCenter); lay.addWidget(chk)
        self.table.setCellWidget(r, 4, w)

    def _del_row(self):
        r = self.table.currentRow()
        if r >= 0:
            self.table.removeRow(r)

    def _from_current(self):
        """δcg of the open workspace spectrum over the VISIBLE x-range --
        through the same dataset route as 'Add from datasets…' (signed CG,
        data-derived σ, FWHM, a stored window and the draggable band), never
        a clipped whole-spectrum integral. The visible range is the window
        when the parent has a view; otherwise the first-minima window is
        seeded. The workspace does not know its processing mode, so the
        row says 'mode unknown' unless the recipe recorded it."""
        if self._current is None:
            return
        larmor, ppm, amp = self._current
        window = None
        parent = self.parent()
        if parent is not None and hasattr(parent, "view"):
            x0, x1 = parent.view.getPlotItem().getViewBox().viewRange()[0]
            window = (float(min(x0, x1)), float(max(x0, x1)))
        mag = None
        recipe = getattr(parent, "recipe", None)
        if isinstance(recipe, dict) and "qcpmg_magnitude" in recipe:
            mag = bool(recipe.get("qcpmg_magnitude"))
        self.add_dataset_spectrum(float(larmor), ppm, amp, window=window,
                                  magnitude=mag, source="workspace spectrum")

    def _points(self) -> list[FieldPoint]:
        pts = []
        for r in range(self.table.rowCount()):
            try:
                nu = float(self.table.item(r, 0).text())
                dcg = float(self.table.item(r, 1).text())
            except (AttributeError, ValueError):
                continue
            try:
                err = float(self.table.item(r, 2).text())
            except (AttributeError, ValueError):
                err = 0.0                     # missing: no uncertainty known
            if not (np.isfinite(err) and err > 0.0):
                err = 0.0                     # negative / NaN: not a sigma
            w = self.table.cellWidget(r, 4)
            sel = None
            if w is not None:
                state = w.findChild(QCheckBox).checkState()
                sel = (None if state == Qt.PartiallyChecked
                       else state == Qt.Checked)
            ds = self._ds.get(self._row_ds_id(r))
            if ds is not None and ds.get("meas") is not None:
                pts.append(FieldPoint.from_measurement(
                    nu, ds["meas"], magnitude=ds.get("magnitude"),
                    source=ds.get("source", ""), rotor_Hz=ds.get("rotor_Hz", 0.0),
                    ct_selective=sel, dcg_ppm=dcg, dcg_err_ppm=err))
            else:
                pts.append(FieldPoint(nu, dcg, err, sel, source="manual"))
        return pts

    def _fields_fwhm(self):
        """(larmor, fwhm_ppm) for rows that have both filled — for Eq. 2."""
        out = []
        for r in range(self.table.rowCount()):
            try:
                nu = float(self.table.item(r, 0).text())
                fw = float(self.table.item(r, 3).text())
                out.append((nu, fw))
            except (AttributeError, ValueError):
                continue
        return out

    def _compute_widths(self):
        from larmor.qcpmg_fields import two_field_widths
        fw = self._fields_fwhm()
        if len(fw) < 2:
            self.wresult.setText("enter the FWHM (ppm) at two fields")
            return
        (n1, f1), (n2, f2) = fw[0], fw[1]
        ws = two_field_widths(n1, f1, n2, f2)
        if not ws.ok:
            self.wresult.setText("⚠ " + ws.note)
            return
        self.wresult.setText(
            f"quadrupolar width W_q = {ws.wq_lo_ppm:.1f} ppm (at {min(n1,n2):.0f} "
            f"MHz) / {ws.wq_hi_ppm:.1f} ppm (at {max(n1,n2):.0f} MHz)  ·  "
            f"chemical-shift-distribution width W_csd = <b>{ws.wcsd_ppm:.1f} "
            f"ppm</b> (field-independent)")

    def _compute(self):
        pts = self._points()
        if len(pts) < 2:
            self.result.setText("need δcg at ≥ 2 fields (fill Larmor + δcg)")
            return
        try:
            res = infinite_field_diso(pts, spin=self.spin.value(),
                                      eta=self.eta.value())
        except Exception as exc:
            self.result.setText(f"cannot extrapolate: {exc}")
            return
        self.plot.clear()
        # the plot may be in dataset-supervision mode (inverted ppm axis)
        self._region = self._cg_line = None
        self.plot.getPlotItem().invertX(False)
        self.plot.setTitle(None)
        self.plot.setLabel("bottom", "1 / ν₀²", units="MHz⁻²")
        self.plot.setLabel("left", "δcg", units="ppm")
        x = np.array([1.0 / p.larmor_MHz ** 2 for p in pts])
        y = np.array([p.dcg_ppm for p in pts])
        self.plot.plot(x, y, pen=None, symbol="o", symbolBrush="#1f6feb",
                       symbolSize=9)
        xs = np.linspace(0.0, float(x.max()) * 1.05, 50)
        self.plot.plot(xs, res.line(xs),
                       pen=pg.mkPen("#c0392b", width=1.6, style=Qt.DashLine))
        self.plot.plot([0.0], [res.delta_iso_ppm], pen=None, symbol="star",
                       symbolBrush="#c0392b", symbolSize=14)
        self.result.setText(self._result_html(res, pts))

    def _result_html(self, res, pts) -> str:
        """The headline: the same wording as the report (a bound is a bound,
        a NaN is '--'), plus the notes and the field-plausibility check."""
        lines = fmt_result_lines(res)
        head = (f"<b>{lines[0].replace('delta_iso', 'δiso')}</b>  "
                f"(intercept, 1/ν₀²→0)   ·   "
                + "   ·   ".join(ln.replace("eta", "η") for ln in lines[1:3]))
        extra = []
        if res.note:
            extra.append(f"<span style='color:#c0392b'>{res.note}</span>")
        if res.warning:
            extra.append(f"<span style='color:#c0392b'>⚠ {res.warning}</span>")
        for p in pts:
            msg = field_plausibility_warning(p.larmor_MHz, self._nucleus)
            if msg:
                extra.append(f"<span style='color:#c0392b'>⚠ {msg}</span>")
        return "<br>".join([head] + extra)


#: the one shared instance — fields sent from several QCPMG processing
#: sessions accumulate here until the user computes the extrapolation
_shared: QcpmgFieldsDialog | None = None


def shared_fields_dialog(parent=None, nucleus: str = "",
                         current=None) -> QcpmgFieldsDialog:
    """Get (or create) the persistent infinite-field dialog. Non-modal by
    design: process one field's dataset, send it here, process the next,
    send it too — then Compute."""
    global _shared
    if _shared is not None:
        try:
            _shared.isVisible()               # raises once the C++ side died
        except RuntimeError:
            _shared = None
    if _shared is None:
        _shared = QcpmgFieldsDialog(parent, nucleus, current)
    else:
        # a different nucleus is adopted while the table is empty; with rows
        # present the dialog keeps its nucleus and shows a warning
        _shared.accept_nucleus(nucleus)
        if current is not None:
            _shared._current = current
            _shared.b_cur.setEnabled(True)
    return _shared
