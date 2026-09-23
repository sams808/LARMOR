"""Dialogs: experiment parameters, parameter links, and fit bounds."""
from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFormLayout, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QVBoxLayout,
)

from larmor.desktop import theme


class BoundsDialog(QDialog):
    """Constrain a fitted parameter between a min and a max.

    Either bound is optional (unchecked = unbounded on that side). The fit
    keeps the parameter inside [min, max] via lmfit box constraints.
    """

    def __init__(self, parent, label: str, p: dict):
        super().__init__(parent)
        self.setWindowTitle("Constrain parameter")
        self.result_min = p.get("min")
        self.result_max = p.get("max")
        v = QVBoxLayout(self)
        v.addWidget(QLabel(f"<b>{label}</b>  (current value: {p['value']:g})"))

        span = abs(p["value"]) or 1.0
        row_lo = QHBoxLayout()
        self.use_min = QCheckBox("minimum")
        self.use_min.setChecked(p.get("min") is not None)
        self.min = QDoubleSpinBox()
        self.min.setDecimals(4); self.min.setRange(-1e12, 1e12)
        self.min.setValue(p["min"] if p.get("min") is not None
                          else p["value"] - span)
        row_lo.addWidget(self.use_min); row_lo.addWidget(self.min, 1)
        v.addLayout(row_lo)

        row_hi = QHBoxLayout()
        self.use_max = QCheckBox("maximum")
        self.use_max.setChecked(p.get("max") is not None)
        self.max = QDoubleSpinBox()
        self.max.setDecimals(4); self.max.setRange(-1e12, 1e12)
        self.max.setValue(p["max"] if p.get("max") is not None
                          else p["value"] + span)
        row_hi.addWidget(self.use_max); row_hi.addWidget(self.max, 1)
        v.addLayout(row_hi)

        note = QLabel("The fit will not let this parameter leave the range. "
                      "A value that ends the fit exactly on a bound is flagged "
                      "in the report.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color: {theme.active().text_dim};")
        v.addWidget(note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        v.addWidget(buttons)

    def _accept(self):
        lo = self.min.value() if self.use_min.isChecked() else None
        hi = self.max.value() if self.use_max.isChecked() else None
        if lo is not None and hi is not None and lo >= hi:
            self.min.setStyleSheet("background: rgba(200,70,60,0.28);")
            return
        self.result_min, self.result_max = lo, hi
        self.accept()


class ExperimentDialog(QDialog):
    """Edit nucleus / Larmor frequency / MAS rate (nu_rot) for the recipe.

    When the recipe carries ``provenance["mas_rate"]`` (a spectrum read from
    an EXPNO), a "Where the rate comes from" group lists the present sources
    (acqus MASR, title, booking sidecar) with a [Use] button each, a
    "measured" row whose [Measure] runs ``measure()`` (the main window's
    sideband detector) and a verdict line; "Remember for this session"
    records the confirmed rate for later spectra of the same session folder,
    rotor and nucleus (``larmor.masrate``). CSV / fxmla / VOCS recipes have
    no block and keep the plain four-row dialog.
    """

    _SOURCE_NAMES = {"acqus": "acqus MASR", "title": "title",
                     "booking": "booking sidecar"}

    def __init__(self, parent, recipe: dict, *, measure=None):
        super().__init__(parent)
        self.setWindowTitle("Experiment parameters")
        self.recipe = recipe
        self._measure_fn = measure
        self._measured_hz = None
        #: {'acqus' | 'title' | 'booking' | 'measured': its Use/Measure button}
        self.source_rows: dict[str, QPushButton] = {}
        #: set by [Forget]: edit_experiment appends a reversal to the store
        self.forget_requested = False
        form = QFormLayout(self)

        self.nucleus = QLineEdit(recipe.get("nucleus", ""))
        self.nucleus.setPlaceholderText("e.g. 27Al, 23Na, 29Si")
        form.addRow("nucleus", self.nucleus)

        self.larmor = QDoubleSpinBox()
        self.larmor.setDecimals(4)
        self.larmor.setRange(0.1, 2000.0)
        self.larmor.setSuffix(" MHz")
        self.larmor.setValue(recipe.get("larmor_frequency_MHz", 100.0) or 100.0)
        form.addRow("Larmor frequency", self.larmor)

        self.mas = QDoubleSpinBox()
        self.mas.setDecimals(1)
        self.mas.setRange(0.0, 300000.0)
        self.mas.setSuffix(" Hz")
        self.mas.setToolTip("MAS spinning rate nu_rot; 0 = static")
        self.mas.setValue(recipe.get("spin_rate_Hz", 0.0) or 0.0)
        form.addRow("MAS rate (νrot)", self.mas)

        # the three-source block written by the loader (larmor.masrate); a
        # CSV / fxmla / VOCS recipe has none and the dialog stays as it was
        block = (recipe.get("provenance") or {}).get("mas_rate") or {}
        self.remember_mas = QCheckBox("Remember for this session")
        self.remember_mas.setChecked(False)
        if block:
            self._build_mas_sources(form, block)
        else:
            self.remember_mas.setVisible(False)

        sr_row = QHBoxLayout()
        self.sr = QDoubleSpinBox()
        self.sr.setDecimals(2); self.sr.setRange(-1e7, 1e7); self.sr.setSuffix(" Hz")
        self.sr.setToolTip("spectral reference SR = SF − BF1; changing it shifts "
                           "the ppm axis by SR/SFO1")
        self.sr.setValue(recipe.get("sr_hz", 0.0) or 0.0)
        sr_row.addWidget(self.sr, 1)
        btnCopy = QPushButton("Copy from spectrum…")
        btnCopy.setToolTip("read SR from another dataset and reference this "
                           "spectrum to match it")
        btnCopy.clicked.connect(self._copy_sr)
        sr_row.addWidget(btnCopy)
        form.addRow("SR (reference)", sr_row)

        note = QLabel("Changing nucleus / field / νrot re-simulates every line; "
                      "changing SR re-references the ppm axis.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color: {theme.active().text_dim};")
        form.addRow(note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    # ------------------------------------------------- where the rate comes from
    def _build_mas_sources(self, form, block: dict):
        from larmor import masrate

        dim = f"color: {theme.active().text_dim};"
        group = QGroupBox("Where the rate comes from")
        grid = QGridLayout(group)
        grid.setColumnStretch(2, 1)
        present = [k for k in ("acqus", "title", "booking")
                   if block.get(f"{k}_Hz") is not None]
        row = 0
        for key in present:
            hz = float(block[f"{key}_Hz"])
            if key == "title":
                detail = block.get("title_note") or (
                    f'"{block["title_raw"]}"' if block.get("title_raw") else "")
            else:
                detail = "experiment_addenda.xml" if key == "booking" else ""
            grid.addWidget(QLabel(self._SOURCE_NAMES[key]), row, 0)
            grid.addWidget(QLabel(masrate.format_hz(hz) + " Hz"), row, 1)
            d = QLabel(detail)
            d.setStyleSheet(dim)
            d.setWordWrap(True)
            grid.addWidget(d, row, 2)
            btn = QPushButton("Use")
            btn.setToolTip(f"set νrot to {masrate.format_hz(hz)} Hz")
            btn.clicked.connect(lambda _=False, v=hz: self.mas.setValue(v))
            grid.addWidget(btn, row, 3)
            self.source_rows[key] = btn
            row += 1

        # the spectrum itself: the ±νrot repeat of the pattern
        grid.addWidget(QLabel("measured"), row, 0)
        self.measured_label = QLabel("—")
        self.measured_label.setWordWrap(True)
        grid.addWidget(self.measured_label, row, 1, 1, 2)
        self.measured_detail = QLabel("from the ±νrot repeat of this spectrum")
        self.measured_detail.setStyleSheet(dim)
        btn = QPushButton("Measure")
        btn.setToolTip("run the spinning-sideband detector on the spectrum "
                       "on screen; the button then sets νrot to the measured rate")
        btn.setEnabled(self._measure_fn is not None)
        btn.clicked.connect(self._measure_or_use)
        grid.addWidget(btn, row, 3)
        self.source_rows["measured"] = btn
        row += 1
        grid.addWidget(self.measured_detail, row, 1, 1, 2)
        row += 1

        # the verdict, and [Forget] when a stored confirmation is in force
        verdict_row = QHBoxLayout()
        self.verdict_label = QLabel(self._verdict(block))
        self.verdict_label.setWordWrap(True)
        self.verdict_label.setStyleSheet(dim)
        verdict_row.addWidget(self.verdict_label, 1)
        if block.get("source") == "confirmed" and block.get("confirmed"):
            self.btnForget = QPushButton("Forget")
            self.btnForget.setToolTip("stop applying this confirmation to the "
                                      "session's spectra (the store keeps a "
                                      "reversal record)")
            self.btnForget.clicked.connect(self._forget)
            verdict_row.addWidget(self.btnForget)
        grid.addLayout(verdict_row, row, 0, 1, 4)
        form.addRow(group)

        session = block.get("session")
        if session:
            key = (session, block.get("rotor") or "", block.get("nucleus") or "")
            self.remember_mas.setText(
                f"Remember for this session ({masrate.describe_key(key)}): "
                "spectra with the same three source values load with this "
                "rate and no warning")
            self.remember_mas.setChecked(bool(block.get("uncertain")))
            form.addRow(self.remember_mas)
        else:
            self.remember_mas.setVisible(False)

    def _verdict(self, block: dict) -> str:
        from larmor import masrate

        names = self._SOURCE_NAMES
        present = [k for k in ("acqus", "title", "booking")
                   if block.get(f"{k}_Hz") is not None]
        src = str(block.get("source") or "")
        if src == "confirmed" and block.get("confirmed"):
            key = (block.get("session") or "", block.get("rotor") or "",
                   block.get("nucleus") or "")
            return (f"Confirmed on {str(block['confirmed'])[:10]} for "
                    f"{masrate.describe_key(key)}.")
        if src == "all":
            return "All sources agree."
        if "+" in src:
            pair = src.split("+")
            text = f"{names[pair[0]].capitalize()} and {names[pair[1]]} agree"
            odd = [k for k in present if k not in pair]
            if odd:
                v = masrate.format_hz(float(block[f"{odd[0]}_Hz"]))
                text += f"; {names[odd[0]]} ({v} Hz) is outvoted"
            return text + "."
        if src == "highest":
            return ("Sources disagree. The title is typed at acquisition, the "
                    "booking when the instrument was reserved; on this "
                    "spectrometer acqus MASR is a leftover. If the spectrum "
                    "shows spinning sidebands, Measure settles it.")
        if src == "fallback":
            return (f"No source found — {masrate.format_hz(masrate.MAS_FALLBACK_HZ)}"
                    " Hz assumed.")
        if src == "static":
            return ("Static (0 Hz)" + (" — no second source confirms it."
                                       if block.get("uncertain") else "."))
        if src in names:
            text = f"Only the {names[src]} carries a rate"
            if block.get("title_note"):
                text += f" ({block['title_note']})"
            return text + "."
        return ""

    def _measure_or_use(self):
        """[Measure] runs the detector; once a rate is measured the same
        button reads [Use] and writes it into the spinbox."""
        if self._measured_hz is not None:
            self.mas.setValue(self._measured_hz)
            return
        self._measure()

    def _measure(self):
        from larmor import sidebands

        if self._measure_fn is None:
            return None
        det = self._measure_fn()
        if det is None or not getattr(det, "ok", False):
            msg = getattr(det, "message", "") or "nothing to run on"
            self.measured_label.setText("no ±νrot repeat found — " + msg)
            return det
        self._measured_hz = float(det.nu_rot_Hz)
        self.measured_label.setText(sidebands.describe(det))
        self.measured_detail.setText("")
        btn = self.source_rows["measured"]
        btn.setText("Use")
        btn.setToolTip(f"set νrot to {sidebands.format_hz(self._measured_hz)} Hz")
        return det

    def _forget(self):
        self.forget_requested = True
        self.remember_mas.setChecked(False)
        self.btnForget.setEnabled(False)
        self.verdict_label.setText(self.verdict_label.text() +
                                   " Will be forgotten on OK.")

    def _accept(self):
        nuc = self.nucleus.text().strip()
        if nuc:
            try:
                from mrsimulator.spin_system.isotope import Isotope

                Isotope(symbol=nuc)   # validates
            except Exception:
                self.nucleus.setStyleSheet("border: 1px solid #b0442e;")
                self.nucleus.setToolTip(f"unknown isotope: {nuc!r}")
                return
        self.recipe["nucleus"] = nuc
        self.recipe["larmor_frequency_MHz"] = float(self.larmor.value())
        self.recipe["spin_rate_Hz"] = float(self.mas.value())
        self.recipe["sr_hz"] = float(self.sr.value())
        self.accept()

    def _copy_sr(self):
        from PySide6.QtWidgets import QFileDialog

        path, _ = QFileDialog.getOpenFileName(
            self, "Read SR from another spectrum", "",
            "Spectra (*.fxmla *.json 1r 2rr *.csv *.txt);;All files (*)")
        if not path:
            return
        try:
            from larmor.io import bruker

            bruker.resolve(path)        # raises on a non-Bruker path
            data = bruker.read(path)
            self.sr.setValue(float(data.meta.get("sr_hz", 0.0)))
        except Exception:
            try:
                from larmor.loader import load_any

                _, _, rec, *_ = load_any(path)
                self.sr.setValue(float(rec.get("sr_hz", 0.0)))
            except Exception as exc:
                self.sr.setToolTip(f"could not read SR: {exc}")


class ProcessingStepsDialog(QDialog):
    """The applied processing pipeline as a removable list (ssNake per-step
    undo). Delete steps; the remaining pipeline is re-applied from the raw
    source."""

    def __init__(self, parent, ops: list[dict]):
        super().__init__(parent)
        self.setWindowTitle("Processing steps")
        self.resize(420, 380)
        self.ops = [dict(o) for o in ops]
        from PySide6.QtWidgets import QListWidget, QPushButton

        v = QVBoxLayout(self)
        v.addWidget(QLabel("Applied steps (top = first). Remove any, then OK to "
                           "re-apply the rest."))
        self.list = QListWidget()
        self._fill()
        v.addWidget(self.list, 1)
        row = QHBoxLayout()
        btnDel = QPushButton("Remove selected"); btnDel.clicked.connect(self._remove)
        row.addWidget(btnDel); row.addStretch(1)
        v.addLayout(row)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        v.addWidget(bb)

    def _fill(self):
        from PySide6.QtWidgets import QListWidgetItem

        self.list.clear()
        for o in self.ops:
            kw = ", ".join(f"{k}={v}" for k, v in o.items() if k != "op")
            self.list.addItem(QListWidgetItem(f"{o['op']}"
                                              + (f"  ({kw})" if kw else "")))

    def _remove(self):
        r = self.list.currentRow()
        if 0 <= r < len(self.ops):
            del self.ops[r]; self._fill()

    def result_ops(self) -> list[dict]:
        return self.ops


class ComputingParamsDialog(QDialog):
    """Tune the Czjzek / MQMAS kernel resolution (dmfit's Computing parameters):
    accuracy vs speed. Clears the kernel caches on OK."""

    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("Computing parameters")
        from larmor import engine, twod

        self.engine, self.twod = engine, twod
        form = QFormLayout(self)
        e, m = engine.KERNEL_SETTINGS, twod.MQMAS_SETTINGS

        def spin(val, lo, hi, dec=0):
            s = QDoubleSpinBox(); s.setDecimals(dec); s.setRange(lo, hi)
            s.setValue(val); return s

        form.addRow(QLabel("<b>1D Czjzek kernel</b>"))
        self.npts = spin(e["npts"], 256, 65536); form.addRow("computed points", self.npts)
        # NO "Cq max" control in 1D: the ceiling is automatic (a 25-400 MHz
        # ladder follows each model's requested width). The old spinbox wrote
        # a setting whose only reader threw its kernel away -- a decorative
        # physics control is worse than none.
        auto = QLabel("automatic (25–400 MHz ladder follows the model width)")
        auto.setStyleSheet(f"color: {theme.active().text_dim};")
        form.addRow("Cq max (1D)", auto)
        self.ncq = spin(e["n_cq"], 10, 200); form.addRow("Cq steps", self.ncq)
        self.neta = spin(e["n_eta"], 3, 41); form.addRow("η steps", self.neta)
        form.addRow(QLabel("<b>MQMAS kernel</b>"))
        self.n2 = spin(m["n2"], 48, 512); form.addRow("F2 points", self.n2)
        self.n1 = spin(m["n1"], 32, 512); form.addRow("F1 points", self.n1)
        self.mncq = spin(m["n_cq"], 8, 120); form.addRow("Cq steps (2D)", self.mncq)
        self.mneta = spin(m["n_eta"], 3, 21); form.addRow("η steps (2D)", self.mneta)
        self.mcqmax = spin(m["cq_max_MHz"], 2, 40, 1); form.addRow("Cq max (2D, MHz)", self.mcqmax)

        info = engine.kernel_cache_info()
        cache = QLabel(f"kernel cache: {info['entries']} kernel(s), "
                       f"{info['mb']:.0f} MB held "
                       f"(budget {engine.KERNEL_CACHE_BUDGET_MB:.0f} MB, LRU)")
        cache.setStyleSheet(f"color: {theme.active().text_dim};")
        form.addRow(cache)
        note = QLabel("More steps/points = more accurate but slower; a change "
                      "rebuilds the kernels on the next fit.")
        note.setWordWrap(True); note.setStyleSheet(f"color: {theme.active().text_dim};")
        form.addRow(note)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self._accept); bb.rejected.connect(self.reject)
        form.addRow(bb)

    def _accept(self):
        self.engine.KERNEL_SETTINGS.update(
            npts=int(self.npts.value()),
            n_cq=int(self.ncq.value()), n_eta=int(self.neta.value()))
        self.twod.MQMAS_SETTINGS.update(
            n2=int(self.n2.value()), n1=int(self.n1.value()),
            n_cq=int(self.mncq.value()), n_eta=int(self.mneta.value()),
            cq_max_MHz=float(self.mcqmax.value()))
        self.engine.clear_kernel_cache(); self.twod.clear_kernel_cache()
        self.accept()


def _other_lines(recipe: dict, exclude: int) -> list[tuple[int, str]]:
    return [(j, s.get("label") or f"s{j}")
            for j, s in enumerate(recipe.get("sites", [])) if j != exclude]


class LinkPositionDialog(QDialog):
    """'This line sits at a fixed offset (ppm or Hz) from another line.'"""

    def __init__(self, parent, recipe: dict, row: int):
        super().__init__(parent)
        self.setWindowTitle("Position relative to another line")
        self.recipe, self.row = recipe, row
        self.expr: str | None = None
        v = QVBoxLayout(self)
        form = QFormLayout()

        self.ref = QComboBox()
        for j, label in _other_lines(recipe, row):
            self.ref.addItem(f"s{j} — {label}", j)
        form.addRow("reference line", self.ref)

        self.offset = QDoubleSpinBox()
        self.offset.setDecimals(4)
        self.offset.setRange(-1e7, 1e7)
        form.addRow("offset", self.offset)

        self.unit = QComboBox()
        self.unit.addItems(["ppm", "Hz"])
        form.addRow("unit", self.unit)
        v.addLayout(form)

        note = QLabel("Hz offsets are converted with the recipe's Larmor "
                      "frequency. The position follows the reference during "
                      "the fit, with error propagation.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color: {theme.active().text_dim};")
        v.addWidget(note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        v.addWidget(buttons)

    def _accept(self):
        if self.ref.currentIndex() < 0:
            self.reject()
            return
        j = self.ref.currentData()
        off = float(self.offset.value())
        if self.unit.currentText() == "Hz":
            larmor = self.recipe.get("larmor_frequency_MHz", 0.0) or 1.0
            off = off / larmor          # Hz -> ppm
        self.expr = f"s{j}.isotropic_chemical_shift_ppm + {off:.6g}"
        self.accept()


class RatioDialog(QDialog):
    """'This line's <param> is a fixed multiple of another line's.'"""

    def __init__(self, parent, recipe: dict, row: int, param: str,
                 title: str, default_ratio: float = 1.0):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.expr: str | None = None
        self.param = param
        v = QVBoxLayout(self)
        form = QFormLayout()
        self.ref = QComboBox()
        for j, label in _other_lines(recipe, row):
            self.ref.addItem(f"s{j} — {label}", j)
        form.addRow("reference line", self.ref)
        self.ratio = QDoubleSpinBox()
        self.ratio.setDecimals(6)
        self.ratio.setRange(1e-6, 1e6)
        self.ratio.setValue(default_ratio)
        form.addRow("ratio", self.ratio)
        v.addLayout(form)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        v.addWidget(buttons)

    def _accept(self):
        if self.ref.currentIndex() < 0:
            self.reject()
            return
        j = self.ref.currentData()
        r = float(self.ratio.value())
        self.expr = (f"s{j}.{self.param}" if r == 1.0
                     else f"{r:.6g} * s{j}.{self.param}")
        self.accept()
