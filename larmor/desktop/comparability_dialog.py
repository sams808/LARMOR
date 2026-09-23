"""Comparability widgets shared by the Batch fit dialog, the Sequential fit
dialog and the Datasets dock (front-end for :mod:`larmor.comparability`).

- :class:`ComparabilityBar` -- the one-line verdict under a dialog's banner:
  hidden for CSV/fxmla-only series, a dim "alike" line, or the fit-health
  strip's *check* styling with **Details…** and **Reprocess all from fid…**;
- :class:`ComparabilityDialog` -- the parameter x spectrum table (majority
  bold, deviants amber, mixed nuclei red, the TopSpin audit line as tooltip)
  with Copy as text / Save CSV… / Help;
- :class:`ReprocessDialog` -- the common-pipeline form (window, LB, GB,
  TDeff, SI, phase), pre-filled from the majority or from one spectrum;
- :func:`apply_reprocess` / :func:`revert_reprocess` -- the procedural
  helpers the batch dialog applies to its per-spectrum data rows and cells.

Amber (``_CHECK_CSS_COLOR``) is a literal like batchfit_dialog's flagged red:
a warning must read as one on every theme. Red stays reserved for *bad*
(mixed nuclei / fields).
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QEventLoop, Signal
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDialog,
    QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout,
    QLabel, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout,
    QWidget,
)

from larmor import comparability as C
from larmor.desktop import theme

__all__ = ["ComparabilityBar", "ComparabilityDialog", "ReprocessDialog",
           "apply_reprocess", "revert_reprocess", "_CHECK_CSS_COLOR"]

#: the referencing dialog's amber -- 'check': differs from the series majority
_CHECK_CSS_COLOR = "#A8570F"
#: the flagged-RMSD / mixed-nuclei red -- 'bad'
_BAD_CSS_COLOR = "#c0392b"


def _pill_css(level: str) -> str:
    t = theme.active()
    if level == "ok":
        return f"color:{t.text_dim}; font-size:11px;"
    bg = t.model if level == "bad" else t.baseline
    return (f"background:{bg}; color:{theme.best_text_on(bg)}; font-weight:600; "
            "padding:3px 8px; border-radius:3px;")


# ------------------------------------------------------------------ the bar
class ComparabilityBar(QWidget):
    """One line under a dialog's banner saying whether the series was
    acquired and processed alike, with the actions that follow from it."""

    details_requested = Signal()
    reprocess_requested = Signal()
    revert_requested = Signal()

    def __init__(self, parent=None, allow_reprocess: bool = True):
        super().__init__(parent)
        self._allow_reprocess = allow_reprocess
        self._cmp = None
        h = QHBoxLayout(self)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(6)
        self.lbl = QLabel("")
        self.lbl.setWordWrap(True)
        self.btnDetails = QPushButton("Details…")
        self.btnDetails.setToolTip("every acquisition (acqus) and processing (procs) "
                                   "parameter of every spectrum, the majority "
                                   "highlighted and the deviants in amber")
        self.btnDetails.clicked.connect(self.details_requested)
        self.btnReprocess = QPushButton("Reprocess all from fid…")
        self.btnReprocess.setToolTip(
            "rebuild every spectrum that has a raw fid with ONE common pipeline "
            "(TDeff, window, zero-fill, FT referenced by its own SF, its own "
            "TopSpin phase or Autophase); the chain is saved in each fit")
        self.btnReprocess.clicked.connect(self.reprocess_requested)
        self.btnRevert = QPushButton("Use TopSpin processing")
        self.btnRevert.setToolTip("back to the spectra as TopSpin processed them "
                                  "(the 1r files); the reprocess chain is dropped")
        self.btnRevert.clicked.connect(self.revert_requested)
        h.addWidget(self.lbl, 1)
        h.addWidget(self.btnDetails)
        h.addWidget(self.btnReprocess)
        h.addWidget(self.btnRevert)
        for b in (self.btnDetails, self.btnReprocess, self.btnRevert):
            b.setAutoDefault(False)
        self.btnRevert.setVisible(False)
        self.hide()

    def text(self) -> str:
        return self.lbl.text()

    def set_comparison(self, cmp, n_with_fid: int = 0) -> None:
        """Hidden when there is nothing to compare (level 'none'); a dim
        one-liner when alike; the check/bad pill otherwise. The reprocess
        button needs at least two members with a fid."""
        self._cmp = cmp
        if cmp is None or cmp.level == "none":
            self.hide()
            return
        lvl = cmp.level
        self.lbl.setText(cmp.summary())
        self.lbl.setStyleSheet(_pill_css(lvl))
        self.btnDetails.setVisible(True)
        self.btnReprocess.setText(
            f"Reprocess all from fid… ({n_with_fid} of {cmp.n} have a fid)")
        self.btnReprocess.setVisible(self._allow_reprocess and lvl != "ok"
                                     and n_with_fid >= 2)
        self.btnRevert.setVisible(False)
        self.show()

    def set_reprocessed(self, template: dict, n_done: int, n_total: int,
                        phase: str = "own") -> None:
        self.lbl.setText("reprocessed from fid with one pipeline: "
                         f"{C.describe_template(template, phase)} — "
                         f"{n_done}/{n_total} spectra")
        self.lbl.setStyleSheet(_pill_css("ok"))
        self.btnDetails.setVisible(True)
        self.btnReprocess.setVisible(False)
        self.btnRevert.setVisible(True)
        self.show()


# ------------------------------------------------------------------ Details
class ComparabilityDialog(QDialog):
    """Parameter x spectrum table of a Comparison (the referencing-audit
    dialog's table pattern)."""

    def __init__(self, parent, cmp, title: str | None = None):
        super().__init__(parent)
        self._cmp = cmp
        self.setWindowTitle(title or "Comparability — acquisition & processing of "
                            f"{cmp.n_bruker} spectra")
        self.resize(min(1200, 420 + 120 * max(1, cmp.n)), 640)
        v = QVBoxLayout(self)
        self.lblSummary = QLabel(cmp.summary())
        self.lblSummary.setWordWrap(True)
        self.lblSummary.setStyleSheet(_pill_css(cmp.level))
        v.addWidget(self.lblSummary)
        self.chkOnlyDiff = QCheckBox("only differences")
        self.chkOnlyDiff.setChecked(True)
        self.chkOnlyDiff.setToolTip("hide the parameters every spectrum shares")
        self.chkOnlyDiff.toggled.connect(lambda _on: self._fill())
        v.addWidget(self.chkOnlyDiff)
        t = self.table = QTableWidget(0, 0)
        t.setEditTriggers(QAbstractItemView.NoEditTriggers)
        t.setSelectionMode(QAbstractItemView.NoSelection)
        t.setAlternatingRowColors(True)
        t.verticalHeader().setVisible(False)
        v.addWidget(t, 1)
        status = QLabel("majority = the most frequent value (ties: the first spectrum); "
                        "amber = differs from it; red = a shared model is meaningless; "
                        "— = no such file. Hover a processing cell for that spectrum's "
                        "TopSpin command line.")
        status.setWordWrap(True)
        status.setStyleSheet(f"color:{theme.active().text_dim}; font-size:10px;")
        v.addWidget(status)
        bb = QDialogButtonBox(QDialogButtonBox.Close | QDialogButtonBox.Help)
        self.btnCopy = bb.addButton("Copy as text", QDialogButtonBox.ActionRole)
        self.btnCopy.clicked.connect(self._copy)
        self.btnCsv = bb.addButton("Save CSV…", QDialogButtonBox.ActionRole)
        self.btnCsv.clicked.connect(self._save_csv)
        bb.helpRequested.connect(self._help)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        v.addWidget(bb)
        for b in self.findChildren(QPushButton):
            b.setAutoDefault(False)
        self._fill()

    def _shown(self) -> list:
        only = self.chkOnlyDiff.isChecked()
        return [r for r in self._cmp.reports
                if not only or r.flagged or (r.mention and r.differs)]

    def _fill(self) -> None:
        cmp, t = self._cmp, self.table
        reports = self._shown()
        t.clear()
        ncol = 2 + cmp.n
        t.setColumnCount(ncol)
        t.setHorizontalHeaderLabels(["parameter", "majority", *cmp.labels])
        for j, p in enumerate(cmp.params):
            it = t.horizontalHeaderItem(2 + j)
            if it is not None:
                it.setToolTip(f"{p.path} · proc {p.procno}" if p is not None
                              else "not a Bruker spectrum (no acqus/procs)")
        bold = QFont()
        bold.setBold(True)
        tt = theme.active()
        if not reports:
            t.setRowCount(1)
            it = QTableWidgetItem("no differences from the series majority — untick "
                                  "“only differences” to list every parameter")
            it.setForeground(QColor(tt.text_dim))
            t.setItem(0, 0, it)
            t.setSpan(0, 0, 1, ncol)
            t.resizeColumnsToContents()
            return
        groups = [(g, [r for r in reports if r.group == g])
                  for g in ("processing", "referencing", "acquisition")]
        groups = [(g, rs) for g, rs in groups if rs]
        t.setRowCount(sum(1 + len(rs) for _g, rs in groups))
        row = 0
        for group, rs in groups:
            head = QTableWidgetItem(group.capitalize())
            head.setFont(bold)
            head.setBackground(QColor(tt.header))
            t.setItem(row, 0, head)
            t.setSpan(row, 0, 1, ncol)
            row += 1
            for r in rs:
                key = QTableWidgetItem(r.label)
                key.setToolTip(r.message or r.key)
                t.setItem(row, 0, key)
                maj = QTableWidgetItem(r.majority_display)
                maj.setFont(bold)
                if r.varies:
                    maj.setText(f"{r.majority_display} (varies)")
                    maj.setFont(QFont())
                t.setItem(row, 1, maj)
                ink = _BAD_CSS_COLOR if r.level == "bad" else _CHECK_CSS_COLOR
                for k, disp in enumerate(r.display):
                    it = QTableWidgetItem(disp)
                    if k in r.deviants:
                        it.setForeground(QColor(ink))
                        it.setFont(bold)
                    elif disp != "—" and not r.varies:
                        it.setFont(bold)
                    p = cmp.params[k]
                    if r.group == "processing" and p is not None and p.commands:
                        it.setToolTip(" · ".join(p.commands))
                    t.setItem(row, 2 + k, it)
                row += 1
        t.resizeColumnsToContents()

    def _copy(self) -> None:
        QApplication.clipboard().setText(
            self._cmp.to_text(only_differences=self.chkOnlyDiff.isChecked()))

    def _save_csv(self) -> None:
        from larmor.desktop.paths import suggest_save_dir

        first = next((p.path for p in self._cmp.params if p is not None), None)
        start = suggest_save_dir(first)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save comparability table",
            str(Path(start) / "comparability.csv") if start else "comparability.csv",
            "CSV (*.csv)")
        if path:
            self._cmp.to_csv(path)

    def _help(self) -> None:
        from larmor.desktop.help_dialog import show_help
        show_help(self, "multi-dataset", "Multi-dataset & co-fitting")


# ------------------------------------------------------------------ Reprocess
class ReprocessDialog(QDialog):
    """The common pipeline every member with a fid is rebuilt with."""

    WINDOWS = (("none", 0), ("EM — exponential (TopSpin LB)", 1),
               ("GM — Gaussian (the glass rule)", 2), ("QSINE — squared sine bell", 4))
    PHASES = (("each spectrum's own PHC0/PHC1 (TopSpin)", "own"),
              ("Autophase", "auto"), ("none", "none"))

    def __init__(self, parent, cmp, n_with_fid: int):
        super().__init__(parent)
        self._cmp = cmp
        self._fcor = 1.0
        self.setWindowTitle("Reprocess all from fid — one common pipeline")
        v = QVBoxLayout(self)
        form = QFormLayout()
        self.cbFrom = QComboBox()
        self.cbFrom.addItem("majority of the series", None)
        for k, (lab, p) in enumerate(zip(cmp.labels, cmp.params)):
            if p is not None and p.proc:
                self.cbFrom.addItem(lab, k)
        self.cbFrom.setToolTip("pre-fill the pipeline from the series majority or "
                               "from one spectrum's own processing")
        form.addRow("take the common pipeline from:", self.cbFrom)
        self.cbWindow = QComboBox()
        for lab, code in self.WINDOWS:
            self.cbWindow.addItem(lab, code)
        form.addRow("window:", self.cbWindow)
        self.spLB = QDoubleSpinBox()
        self.spLB.setRange(-2000.0, 2000.0)
        self.spLB.setDecimals(1)
        self.spLB.setSuffix(" Hz")
        self.spLB.setToolTip("EM: line broadening (Hz). GM: TopSpin's LB is "
                             "NEGATIVE (the Lorentzian it removes)")
        form.addRow("LB:", self.spLB)
        self.spGB = QDoubleSpinBox()
        self.spGB.setRange(0.001, 1.0)
        self.spGB.setDecimals(3)
        self.spGB.setSingleStep(0.05)
        self.spGB.setToolTip("GM only: the Gaussian maximum as a fraction of the "
                             "acquisition time (TopSpin GB)")
        form.addRow("GB:", self.spGB)
        self.spTDeff = QSpinBox()
        self.spTDeff.setRange(0, 4_194_304)
        self.spTDeff.setSpecialValueText("whole fid")
        self.spTDeff.setToolTip("TopSpin TDeff, in REAL points as TopSpin counts them "
                                "(0 = the whole fid)")
        form.addRow("TDeff (TopSpin points):", self.spTDeff)
        self.spSI = QSpinBox()
        self.spSI.setRange(256, 4_194_304)
        self.spSI.setToolTip("zero-fill to this many points (TopSpin SI)")
        form.addRow("SI:", self.spSI)
        self.cbPhase = QComboBox()
        for lab, code in self.PHASES:
            self.cbPhase.addItem(lab, code)
        form.addRow("phase:", self.cbPhase)
        v.addLayout(form)
        n_without = cmp.n - n_with_fid
        note = QLabel(
            "Every spectrum with a raw fid is rebuilt with this one chain, "
            "referenced by its own SF (procs), and the chain is saved in each fit "
            "(processing, processing_from_raw, the EXPNO as source) so the Plotting "
            "studio and a reopened project replay it from the instrument fid. "
            + (f"{n_without} spectrum/spectra without a fid keep their TopSpin "
               "processing and stay flagged. " if n_without else "")
            + "The glass protocol (Glass fitting §2) recommends a Gaussian window "
              "(GM) rather than EM. Baselines are cleared — Fit baseline… still "
              "applies afterwards; TopSpin's abs (ABSG) is not replayed.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{theme.active().text_dim}; font-size:10px;")
        v.addWidget(note)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
        self.cbFrom.currentIndexChanged.connect(lambda _i: self._prefill())
        self.cbWindow.currentIndexChanged.connect(lambda _i: self._sync_enabled())
        self._prefill()

    def _prefill(self) -> None:
        k = self.cbFrom.currentData()
        if k is None:
            t = C.majority_template(self._cmp)
        else:
            t = C.template_from(self._cmp.params[k])
        self._fcor = float(t.get("fcor", 1.0))
        wdw = int(t.get("wdw") or 0)
        idx = self.cbWindow.findData(wdw)
        if idx < 0:                         # SINE -> QSINE, exotic windows -> none
            idx = self.cbWindow.findData(4 if wdw == 3 else 0)
        self.cbWindow.setCurrentIndex(idx)
        self.spLB.setValue(float(t.get("lb_hz") or 0.0))
        self.spGB.setValue(float(t.get("gb") or 0.1))
        self.spTDeff.setValue(int(t.get("tdeff") or 0))
        self.spSI.setValue(int(t.get("si") or 32768))
        self._sync_enabled()

    def _sync_enabled(self) -> None:
        wdw = self.cbWindow.currentData()
        self.spLB.setEnabled(wdw in (1, 2))
        self.spGB.setEnabled(wdw == 2)

    def template(self) -> dict:
        wdw = int(self.cbWindow.currentData() or 0)
        return {"wdw": wdw, "lb_hz": float(self.spLB.value()),
                "gb": float(self.spGB.value()), "ssb": 2.0,
                "tdeff": int(self.spTDeff.value()), "si": int(self.spSI.value()),
                "fcor": self._fcor}

    def phase(self) -> str:
        return str(self.cbPhase.currentData() or "own")


# ------------------------------------------------------------------ helpers
def _set_cell(cells: list, k: int, ppm, amp) -> None:
    if k >= len(cells):
        return
    cell = cells[k]
    exp = cell.get("exp")
    if exp is not None:
        exp.setData(ppm, amp)
    plot = cell.get("plot")
    if plot is not None and len(ppm):
        plot.setXRange(float(ppm.min()), float(ppm.max()))


def apply_reprocess(data: list, cells: list, template: dict, phase: str = "own",
                    status_cb=None) -> tuple[int, list]:
    """Rebuild every row whose params carry a fid with ``template`` + its own
    phase; the row's ppm/amp/amp0, proc_ops, expno, baseline_ops (cleared)
    and S/N are replaced and its cell redrawn. Members without a fid and
    failures become notes. Returns ``(n_done, notes)``."""
    from larmor.desktop.batchfit_dialog import _snr

    n_done, notes = 0, []
    for k, d in enumerate(data):
        p = d.get("params")
        name = d.get("sample") or str(k + 1)
        if p is None:
            notes.append(f"{name}: not a Bruker spectrum — kept as is")
            continue
        if not p.has_fid:
            notes.append(f"{name}: no fid — TopSpin spectrum kept")
            continue
        if status_cb is not None:
            status_cb(f"reprocessing {k + 1}/{len(data)} — {name}…")
            QApplication.processEvents(QEventLoop.ExcludeUserInputEvents)
        try:
            ops = C.ops_for(p, template, phase)
            ppm, amp, _notes = C.reprocess(p, ops)
        except Exception as exc:                          # noqa: BLE001
            notes.append(f"{name}: could not reprocess ({exc}) — TopSpin spectrum kept")
            continue
        d["ppm"], d["amp"], d["amp0"] = ppm, amp, amp.copy()
        d["proc_ops"] = [dict(o) for o in ops]
        d["expno"] = p.expno
        d["baseline_ops"] = []
        d["snr"] = _snr(amp)
        _set_cell(cells, k, ppm, amp)
        n_done += 1
    return n_done, notes


def revert_reprocess(data: list, cells: list) -> None:
    """Back to the TopSpin arrays kept in ppm_src / amp_src; the chain and any
    baseline are dropped."""
    from larmor.desktop.batchfit_dialog import _snr

    for k, d in enumerate(data):
        src_ppm, src_amp = d.get("ppm_src"), d.get("amp_src")
        if src_ppm is None or src_amp is None:
            continue
        d["ppm"], d["amp"], d["amp0"] = src_ppm.copy(), src_amp.copy(), src_amp.copy()
        d["proc_ops"] = []
        d["expno"] = ""
        d["baseline_ops"] = []
        d["snr"] = _snr(d["amp"])
        _set_cell(cells, k, d["ppm"], d["amp"])
