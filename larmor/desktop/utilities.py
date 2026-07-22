"""Utilities dialogs (ssNake Utilities parity): the NMR table and the
chemical-shift / quadrupole / dipolar conversion calculators."""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDoubleSpinBox, QGridLayout, QGroupBox, QHBoxLayout,
    QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from larmor import convert as C
from larmor import nuclei as N

#: colour per nuclear spin (echoes ssNake's table)
SPIN_COLOR = {0.5: "#2b6cb0", 1.0: "#c05621", 1.5: "#2f855a", 2.5: "#c53030",
              3.5: "#8a6d1a", 4.5: "#7b341e"}


def _spin_color(spin: float) -> str:
    return SPIN_COLOR.get(spin, "#4a5568")


class NmrTableDialog(QDialog):
    """Interactive periodic table of Larmor frequencies. Set B0 (T) or your
    magnet's ¹H frequency; double-click an element for its isotopes."""

    def __init__(self, parent, h1_MHz: float = 400.0):
        super().__init__(parent)
        self.setWindowTitle("NMR table")
        self.resize(1180, 560)
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        top.addWidget(QLabel("B₀ [T]"))
        self.b0 = QDoubleSpinBox(); self.b0.setDecimals(3); self.b0.setRange(0.01, 40.0)
        self.b0.setValue(N.b0_from_1H(h1_MHz))
        top.addWidget(self.b0)
        top.addWidget(QLabel("¹H [MHz]"))
        self.h1 = QDoubleSpinBox(); self.h1.setDecimals(2); self.h1.setRange(1.0, 1700.0)
        self.h1.setValue(h1_MHz)
        top.addWidget(self.h1)
        top.addStretch(1)
        leg = " ".join(f"<span style='color:{c}'>●</span> {int(s*2)}/2" if s % 1
                       else f"<span style='color:{c}'>●</span> {int(s)}"
                       for s, c in sorted(SPIN_COLOR.items()))
        lab = QLabel("spin " + leg); top.addWidget(lab)
        v.addLayout(top)
        self._guard = False
        self.b0.valueChanged.connect(self._b0_changed)
        self.h1.valueChanged.connect(self._h1_changed)

        grid = QGridLayout(); grid.setSpacing(2)
        self._buttons = {}
        for r, row in enumerate(N.PERIODIC_ROWS):
            for c, sym in enumerate(row):
                if not sym:
                    continue
                b = QPushButton(); b.setFixedSize(58, 44)
                b.clicked.connect(lambda _=False, s=sym: self._details(s))
                grid.addWidget(b, r, c)
                self._buttons[sym] = b
        holder = QWidget(); holder.setLayout(grid)
        v.addWidget(holder)
        self._refresh()

    def _b0_changed(self, val):
        if self._guard:
            return
        self._guard = True
        self.h1.setValue(val * N.GAMMA_1H)
        self._guard = False
        self._refresh()

    def _h1_changed(self, val):
        if self._guard:
            return
        self._guard = True
        self.b0.setValue(N.b0_from_1H(val))
        self._guard = False
        self._refresh()

    def _refresh(self):
        b0 = self.b0.value()
        for sym, btn in self._buttons.items():
            iso = N.primary_isotope(sym)
            if iso is None:
                btn.setText(sym); btn.setEnabled(False)
                btn.setStyleSheet("color:#a0aec0;")
                continue
            btn.setText(f"{sym}\n{iso.larmor_MHz(b0):.1f}")
            col = _spin_color(iso.spin)
            btn.setStyleSheet(f"border: 1.5px solid {col}; border-radius: 3px; "
                              f"font-size: 9px; color: #16202a;")
            btn.setToolTip(f"{iso.symbol} · spin {iso.spin} · "
                           f"{iso.abundance:.2f}% · {iso.larmor_MHz(b0):.3f} MHz")

    def _details(self, element: str):
        IsotopeDetailsDialog(self, element, self.b0.value()).exec()


class IsotopeDetailsDialog(QDialog):
    def __init__(self, parent, element: str, b0_T: float):
        super().__init__(parent)
        self.setWindowTitle(f"{element} — {N.ELEMENT_NAME.get(element, '')}")
        self.resize(620, 300)
        v = QVBoxLayout(self)
        isos = N.isotopes_by_element().get(element, [])
        cols = ["isotope", "spin", "abund %", "γ (MHz/T)", "Q (barn)",
                f"ν @ {b0_T:.2f} T (MHz)", "rec. (¹H=1)"]
        t = QTableWidget(len(isos), len(cols))
        t.setHorizontalHeaderLabels(cols)
        for i, iso in enumerate(isos):
            vals = [iso.symbol, f"{iso.spin:g}", f"{iso.abundance:.3f}",
                    f"{iso.gamma_MHz_T:.4f}", f"{iso.quad_moment_barn:.4f}",
                    f"{iso.larmor_MHz(b0_T):.3f}", f"{iso.receptivity_1H:.2e}"]
            for j, val in enumerate(vals):
                t.setItem(i, j, QTableWidgetItem(val))
        t.resizeColumnsToContents()
        v.addWidget(t)


class ConvertDialog(QDialog):
    """Chemical-shift / quadrupole / dipolar-distance conversions."""

    def __init__(self, parent, sfo_MHz: float = 100.0):
        super().__init__(parent)
        self.setWindowTitle("Conversion tools")
        self.resize(460, 520)
        v = QVBoxLayout(self)

        # --- chemical shift ---
        g1 = QGroupBox("Chemical shift"); f1 = QGridLayout(g1)
        f1.addWidget(QLabel("SFO (MHz)"), 0, 0)
        self.sfo = QDoubleSpinBox(); self.sfo.setDecimals(4); self.sfo.setRange(0.1, 1700)
        self.sfo.setValue(sfo_MHz); f1.addWidget(self.sfo, 0, 1)
        f1.addWidget(QLabel("ppm"), 1, 0)
        self.ppm = QDoubleSpinBox(); self.ppm.setDecimals(3); self.ppm.setRange(-1e6, 1e6)
        f1.addWidget(self.ppm, 1, 1)
        f1.addWidget(QLabel("Hz"), 2, 0)
        self.hz = QDoubleSpinBox(); self.hz.setDecimals(2); self.hz.setRange(-1e9, 1e9)
        f1.addWidget(self.hz, 2, 1)
        self.ppm.valueChanged.connect(
            lambda x: self._set(self.hz, C.ppm_to_Hz(x, self.sfo.value())))
        self.hz.valueChanged.connect(
            lambda x: self._set(self.ppm, C.Hz_to_ppm(x, self.sfo.value())))
        v.addWidget(g1)

        # --- quadrupole ---
        g2 = QGroupBox("Quadrupole (spin I)"); f2 = QGridLayout(g2)
        f2.addWidget(QLabel("spin I"), 0, 0)
        self.spin = QComboBox(); self.spin.addItems(["1", "1.5", "2.5", "3.5", "4.5"])
        self.spin.setCurrentText("2.5"); f2.addWidget(self.spin, 0, 1)
        f2.addWidget(QLabel("Larmor ν₀ (MHz)"), 0, 2)
        self.nu0 = QDoubleSpinBox(); self.nu0.setDecimals(3); self.nu0.setRange(1, 1700)
        self.nu0.setValue(sfo_MHz); f2.addWidget(self.nu0, 0, 3)
        f2.addWidget(QLabel("Cq (MHz)"), 1, 0)
        self.cq = QDoubleSpinBox(); self.cq.setDecimals(4); self.cq.setRange(0, 100)
        self.cq.setValue(2.0); f2.addWidget(self.cq, 1, 1)
        f2.addWidget(QLabel("η"), 1, 2)
        self.eta = QDoubleSpinBox(); self.eta.setDecimals(3); self.eta.setRange(0, 1)
        f2.addWidget(self.eta, 1, 3)
        self.qout = QLabel(""); self.qout.setWordWrap(True)
        f2.addWidget(self.qout, 2, 0, 1, 4)
        for w in (self.cq, self.eta, self.nu0):
            w.valueChanged.connect(self._quad)
        self.spin.currentTextChanged.connect(self._quad)
        v.addWidget(g2)

        # --- dipolar ---
        g3 = QGroupBox("Dipolar coupling / distance"); f3 = QGridLayout(g3)
        nmr = sorted({i.symbol for i in N.all_isotopes() if i.spin > 0},
                     key=lambda s: (N._split_symbol(s)[1], N._split_symbol(s)[0]))
        f3.addWidget(QLabel("nucleus 1"), 0, 0)
        self.n1 = QComboBox(); self.n1.addItems(nmr); self.n1.setCurrentText("1H")
        f3.addWidget(self.n1, 0, 1)
        f3.addWidget(QLabel("nucleus 2"), 0, 2)
        self.n2 = QComboBox(); self.n2.addItems(nmr); self.n2.setCurrentText("13C")
        f3.addWidget(self.n2, 0, 3)
        f3.addWidget(QLabel("distance (Å)"), 1, 0)
        self.dist = QDoubleSpinBox(); self.dist.setDecimals(4); self.dist.setRange(0.1, 100)
        self.dist.setValue(2.0); f3.addWidget(self.dist, 1, 1)
        f3.addWidget(QLabel("D (Hz)"), 1, 2)
        self.dcoup = QDoubleSpinBox(); self.dcoup.setDecimals(2); self.dcoup.setRange(0, 1e9)
        f3.addWidget(self.dcoup, 1, 3)
        self._d_guard = False
        self.dist.valueChanged.connect(self._dip_from_r)
        self.dcoup.valueChanged.connect(self._dip_from_d)
        for w in (self.n1, self.n2):
            w.currentTextChanged.connect(self._dip_from_r)
        v.addWidget(g3)
        v.addStretch(1)

        self._quad(); self._dip_from_r()

    def _set(self, spin, val):
        spin.blockSignals(True); spin.setValue(val); spin.blockSignals(False)

    def _quad(self, *_):
        I = float(self.spin.currentText()); cq = self.cq.value(); eta = self.eta.value()
        pq = C.pq_from_cq_eta(cq, eta); nq = C.nu_q(cq, I)
        ct = C.ct_second_order_shift_ppm(pq, I, self.nu0.value())
        self.qout.setText(
            f"PQ = <b>{pq:.4f}</b> MHz   ·   νQ = <b>{nq*1e3:.1f}</b> kHz   ·   "
            f"CT 2nd-order shift = <b>{ct:.2f}</b> ppm")

    def _gamma(self, sym):
        return next(i.gamma_MHz_T for i in N.all_isotopes() if i.symbol == sym)

    def _dip_from_r(self, *_):
        if self._d_guard:
            return
        self._d_guard = True
        d = C.dipolar_Hz(self._gamma(self.n1.currentText()),
                         self._gamma(self.n2.currentText()), self.dist.value())
        self.dcoup.setValue(abs(d))
        self._d_guard = False

    def _dip_from_d(self, *_):
        if self._d_guard:
            return
        self._d_guard = True
        r = C.distance_from_dipolar(self._gamma(self.n1.currentText()),
                                    self._gamma(self.n2.currentText()),
                                    self.dcoup.value())
        if r != float("inf"):
            self.dist.setValue(r)
        self._d_guard = False
