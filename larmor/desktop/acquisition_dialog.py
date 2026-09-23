"""Experimental section window (Tools > Experimental section…; the
'Acquisition table…' button of Batch fit and Sequential fit): the
acquisition table of a set of spectra with every column that varies across
the series highlighted, the Experimental paragraph with ranges for those
parameters, and Copy / Save actions. Needs no fit: everything comes from
acqus / pdata/N/procs / title / uxnmr.info / the booking sidecar through
larmor.acquisition (Qt-free, tested). Non-modal, kept alive on the main
window like the referencing audit; nothing is ever written into an
instrument folder."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QFileDialog, QHBoxLayout, QHeaderView, QLabel,
    QListWidget, QPlainTextEdit, QPushButton, QSplitter, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

from larmor import acquisition as A
from larmor.desktop import theme
from larmor.desktop.paths import is_instrument_file, suggest_save_dir

#: a varying cell's background and the MAS-badge red of an unconfirmed rate
_VARY_BG = "#F4DFA3"
_MAS_RED = "#9E3324"


class AcquisitionTableDialog(QDialog):
    """``paths``: spectra (EXPNO / pdata / 1r paths) to start with;
    ``session_hint``: the month folder for the referencing evidence and the
    Add-folder dialog; ``spin_rates``: {path: (rate_Hz, uncertain)} -- the
    fitted recipes' confirmed rates, keyed by the given path or by the
    EXPNO; ``recipe``: the open fit (its audit block counts as referencing
    evidence when it is the one spectrum listed)."""

    def __init__(self, parent=None, paths=(), session_hint: str = "",
                 spin_rates: dict | None = None, recipe: dict | None = None):
        super().__init__(parent)
        self.setWindowTitle("Experimental section — acquisition parameters")
        self.setModal(False)
        self.resize(1000, 650)
        self.blocks: list[dict] = []
        self.acq_table = A.AcqTable()
        self.spin_rates: dict = dict(spin_rates or {})
        self.recipe = recipe
        self.session_hint = session_hint or ""
        t = theme.active()
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(QLabel("spectra (one row per EXPNO)"))
        self.list = QListWidget()
        self.list.setSelectionMode(QListWidget.ExtendedSelection)
        lv.addWidget(self.list, 1)
        row = QHBoxLayout()
        self.btnAdd = QPushButton("Add spectra…")
        self.btnAdd.clicked.connect(self._add_dialog)
        self.btnFolder = QPushButton("Add folder…  (sample or month)")
        self.btnFolder.clicked.connect(self._folder_dialog)
        self.btnRemove = QPushButton("Remove")
        self.btnRemove.clicked.connect(self._remove)
        for w in (self.btnAdd, self.btnFolder, self.btnRemove):
            row.addWidget(w)
        lv.addLayout(row)
        top.addWidget(left, 1)
        v.addLayout(top, 0)

        self.status = QLabel("add spectra to build the table")
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color: {t.text_dim};")
        v.addWidget(self.status)

        split = QSplitter(Qt.Vertical)
        self.table = QTableWidget(0, 0)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        split.addWidget(self.table)
        self.para = QPlainTextEdit()
        self.para.setReadOnly(True)
        self.para.setPlaceholderText("the Experimental paragraph appears here — ranges "
                                     "'(Table S1)' for every parameter that varies, "
                                     "brackets for what the files do not prove")
        split.addWidget(self.para)
        split.setSizes([380, 200])
        v.addWidget(split, 1)

        bottom = QHBoxLayout()
        self.chkCompact = QCheckBox("compact LaTeX table (constant parameters go to the caption)")
        self.chkCompact.setChecked(True)
        bottom.addWidget(self.chkCompact)
        bottom.addStretch(1)
        self.btnCopyPara = QPushButton("Copy paragraph")
        self.btnCopyPara.clicked.connect(self._copy_para)
        self.btnCopyTex = QPushButton("Copy LaTeX table")
        self.btnCopyTex.clicked.connect(self._copy_tex)
        self.btnSave = QPushButton("Save CSV + LaTeX…")
        self.btnSave.clicked.connect(self._save)
        btnHelp = QPushButton("Help")
        btnHelp.clicked.connect(self._help)
        btnClose = QPushButton("Close")
        btnClose.clicked.connect(self.close)
        for w in (self.btnCopyPara, self.btnCopyTex, self.btnSave, btnHelp, btnClose):
            bottom.addWidget(w)
        v.addLayout(bottom)
        for b in self.findChildren(QPushButton):
            b.setAutoDefault(False)
        self.add_paths(list(paths or []))

    # ------------------------------------------------------------- data
    @staticmethod
    def _expno_of(path) -> str:
        from larmor.io import bruker

        p = Path(str(path))
        try:
            return str(bruker.resolve(p).expno)
        except (ValueError, FileNotFoundError, OSError):
            return str(p) if bruker.is_expno(p) else ""

    def add_paths(self, paths):
        """Read the acquisition block of every new EXPNO (duplicates by
        resolved EXPNO are ignored; a CSV / fxmla / recipe path has no block
        and is reported); the table and paragraph rebuild once."""
        have = {b["expno_path"] for b in self.blocks}
        skipped = 0
        for k, p in enumerate(paths or []):
            expno = self._expno_of(p)
            if not expno or expno in have:
                if not expno:
                    skipped += 1
                continue
            b = A.read_block(p)
            if not b:
                skipped += 1
                continue
            if str(p) in self.spin_rates and expno not in self.spin_rates:
                self.spin_rates[expno] = self.spin_rates.pop(str(p))
            self.blocks.append(b)
            have.add(expno)
            self.list.addItem(f"{b.get('sample') or Path(expno).parent.name}/{b.get('expno')}"
                              f"  ·  {b.get('nucleus')}  ·  {expno}")
            if k % 10 == 9:
                self.status.setText(f"reading… {len(self.blocks)} spectra")
                QApplication.processEvents()
        self.rebuild()
        if skipped:
            self.status.setText(self.status.text()
                                + f"  ·  {skipped} path(s) without an acqus skipped")

    def scan_folder(self, folder):
        """A sample folder adds its EXPNOs, a month folder the whole session."""
        self.add_paths(A.folder_expnos(folder))

    def _remove(self):
        rows = sorted({i.row() for i in self.list.selectedIndexes()}, reverse=True)
        for r in rows:
            if 0 <= r < len(self.blocks):
                del self.blocks[r]
                self.list.takeItem(r)
        self.rebuild()

    def _add_dialog(self):
        start = self.session_hint or (str(Path(self.blocks[0]["expno_path"]).parent)
                                      if self.blocks else str(Path.home()))
        files, _ = QFileDialog.getOpenFileNames(self, "Add spectra (1r files, or any file "
                                                      "inside the EXPNO)", start)
        if files:
            self.add_paths(files)

    def _folder_dialog(self):
        start = self.session_hint or str(Path.home())
        d = QFileDialog.getExistingDirectory(self, "Add a sample or month folder", start)
        if d:
            self.scan_folder(d)

    # ------------------------------------------------------------- view
    def _evidence(self):
        rec = None
        if self.recipe and len(self.blocks) == 1:
            src = str(self.recipe.get("source_path") or "")
            if src and Path(src) == Path(self.blocks[0]["expno_path"]):
                rec = self.recipe
        return A.referencing_evidence(self.blocks, self.session_hint or None, rec)

    def rebuild(self):
        t = A.table(self.blocks, self.spin_rates or None)
        self.acq_table = t
        self.table.clear()
        self.table.setRowCount(len(t.rows))
        self.table.setColumnCount(len(t.columns))
        for j, (key, head, *_r) in enumerate(t.columns):
            h = QTableWidgetItem(head + ("  ▲" if key in t.varying else ""))
            if key in t.varying:
                vals = [A._cell(t, key, r.get(key)) for r in t.rows]
                nums = [r.get(key) for r in t.rows
                        if isinstance(r.get(key), (int, float)) and not isinstance(r.get(key), bool)]
                span = (f"{A._cell(t, key, min(nums))} – {A._cell(t, key, max(nums))}"
                        if nums and len(nums) == len(t.rows)
                        else " / ".join(dict.fromkeys(v for v in vals if v)))
                h.setToolTip(f"varies across the series: {span}")
            self.table.setHorizontalHeaderItem(j, h)
            for i, r in enumerate(t.rows):
                txt = A._cell(t, key, r.get(key))
                it = QTableWidgetItem(txt)
                it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                it.setData(Qt.UserRole, r.get("expno_path"))
                if key in t.varying:
                    it.setBackground(QColor(_VARY_BG))
                    it.setForeground(QColor("#222222"))
                    it.setToolTip(h.toolTip())
                if key == "spin_rate_Hz" and r.get("mas_uncertain"):
                    cands = " · ".join(
                        f"{n} {A._khz(vv)[:-4]}" for n, vv in A.mas_candidates(r))
                    it.setText(f"⚠ {txt}  ({cands})" if cands else f"⚠ {txt}")
                    it.setForeground(QColor(_MAS_RED))
                    it.setToolTip("MAS rate not confirmed — the recipe's Experiment "
                                  "parameters dialog settles it")
                self.table.setItem(i, j, it)
        self.para.setPlainText(A.paragraph(self.blocks, spin_rates=self.spin_rates or None,
                                           referencing_text=self._evidence())
                               if self.blocks else "")
        nuclei = sorted({str(r.get("nucleus") or "") for r in t.rows})
        n_unc = sum(1 for r in t.rows if r.get("mas_uncertain"))
        var = t.varying_headers()
        self.status.setText(
            f"{len(t.rows)} spectra · {len(nuclei)} {'nucleus' if len(nuclei) == 1 else 'nuclei'}"
            f" · varying: {', '.join(var) if var else 'none'} · MAS unconfirmed for {n_unc}"
            if t.rows else "add spectra to build the table")
        on = bool(t.rows)
        for w in (self.btnCopyPara, self.btnCopyTex, self.btnSave):
            w.setEnabled(on)

    # ---------------------------------------------------------- actions
    def _tex(self) -> str:
        return A.table_latex(self.acq_table, caption="Acquisition parameters.",
                             compact=self.chkCompact.isChecked())

    def _copy_para(self):
        QApplication.clipboard().setText(self.para.toPlainText())
        self.status.setText("Experimental paragraph copied to the clipboard — check every "
                            "[bracket] before pasting")

    def _copy_tex(self):
        QApplication.clipboard().setText(self._tex())
        self.status.setText("LaTeX table copied to the clipboard (ASCII/TeX only)")

    def _save(self):
        if not self.blocks:
            return
        start = suggest_save_dir(self.blocks[0]["expno_path"], str(Path.home()))
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV + LaTeX table",
                                              str(Path(start) / "acquisition.csv"),
                                              "CSV (*.csv)")
        if not path:
            return
        p = Path(path)
        if is_instrument_file(p):
            self.status.setText("refused: that name would replace an acquired file")
            return
        if p.suffix.lower() != ".csv":
            p = p.with_suffix(".csv")
        tex = p.with_suffix(".tex")
        if is_instrument_file(p) or is_instrument_file(tex):
            self.status.setText("refused: that name would replace an acquired file")
            return
        try:
            p.write_text(A.table_csv(self.acq_table), encoding="utf-8")
            tex.write_text(self._tex(), encoding="utf-8")
        except OSError as exc:
            self.status.setText(f"could not write: {exc}")
            return
        self.status.setText(f"wrote {p.name} and {tex.name} to {p.parent}")

    def _help(self):
        from larmor.desktop.help_dialog import show_help

        show_help(self, "multi-dataset", "Multi-dataset & co-fitting")
