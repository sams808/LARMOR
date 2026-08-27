"""Parameter correlation heat-map from the last fit's covariance.
Shows which fitted parameters the data cannot separate — |correlation| near 1
(deep red/blue) means the pair trades off and their individual errors are large."""
from __future__ import annotations

import numpy as np
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QHeaderView, QLabel, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)


def _corr_color(c: float) -> QColor:
    """+1 → red, 0 → white, −1 → blue (intensity = |c|)."""
    a = min(abs(c), 1.0)
    if c >= 0:
        return QColor(int(255), int(255 * (1 - a)), int(255 * (1 - a)))
    return QColor(int(255 * (1 - a)), int(255 * (1 - a)), int(255))


from larmor.desktop import theme


class CorrelationDialog(QDialog):
    def __init__(self, parent, lmfit_result):
        super().__init__(parent)
        self.setWindowTitle("Parameter correlations")
        self.resize(640, 560)
        v = QVBoxLayout(self)

        names = list(getattr(lmfit_result, "var_names", []) or [])
        covar = getattr(lmfit_result, "covar", None)
        if not names or covar is None:
            v.addWidget(QLabel(
                "No covariance available — run a fit first (a fit pinned at "
                "bounds or with fixed parameters may not report one)."))
            bb = QDialogButtonBox(QDialogButtonBox.Close)
            bb.rejected.connect(self.reject)
            bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
            v.addWidget(bb)
            return

        cov = np.asarray(covar, float)
        d = np.sqrt(np.clip(np.diag(cov), 1e-300, None))
        corr = cov / np.outer(d, d)

        v.addWidget(QLabel(
            "Correlation of the fitted parameters. <b>Deep red/blue (|r|→1)</b> "
            "means the two parameters trade off — the data can't separate them, "
            "so their individual error bars are inflated."))
        v.itemAt(v.count() - 1).widget().setWordWrap(True)

        n = len(names)
        t = QTableWidget(n, n)
        short = [nm.replace("isotropic_chemical_shift_ppm", "pos")
                   .replace("sigma_Cq_MHz", "σCq").replace("shift_fwhm_ppm", "dCS")
                   .replace("amplitude", "amp").replace("_", ".") for nm in names]
        t.setHorizontalHeaderLabels(short)
        t.setVerticalHeaderLabels(short)
        for i in range(n):
            for j in range(n):
                c = float(corr[i, j])
                it = QTableWidgetItem(f"{c:+.2f}")
                it.setBackground(_corr_color(c))
                if abs(c) > 0.6:
                    it.setForeground(QColor("white") if abs(c) > 0.8
                                     else QColor("#16202a"))
                t.setItem(i, j, it)
        t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        v.addWidget(t, 1)

        # flag the strongly-correlated pairs
        pairs = [(abs(corr[i, j]), short[i], short[j], corr[i, j])
                 for i in range(n) for j in range(i + 1, n)]
        pairs.sort(reverse=True)
        strong = [f"{a}↔{b} ({c:+.2f})" for m, a, b, c in pairs[:6] if m > 0.8]
        note = QLabel("Strongest: " + ("  ·  ".join(strong) if strong
                                       else "none above 0.8 — parameters are "
                                       "well separated"))
        note.setWordWrap(True); note.setStyleSheet(f"color:{theme.active().text_dim};")
        v.addWidget(note)

        # identifiability: pairs the data genuinely CANNOT separate (|r| ≥ 0.95)
        from larmor.identifiability import (unidentifiable_pairs,
                                            IDENTIFIABILITY_THRESHOLD)
        uni = unidentifiable_pairs(lmfit_result)

        def _short(nm):
            return (nm.replace("isotropic_chemical_shift_ppm", "pos")
                      .replace("sigma_Cq_MHz", "σCq").replace("shift_fwhm_ppm", "dCS")
                      .replace("amplitude", "amp").replace("_", "."))
        if uni:
            txt = "  ·  ".join(f"{_short(a)}↔{_short(b)} ({r:+.2f})"
                               for a, b, r in uni[:6])
            banner = QLabel(f"⚠ <b>Unidentifiable</b> (|r| ≥ "
                            f"{IDENTIFIABILITY_THRESHOLD:.2f}): {txt}. "
                            "These pairs are degenerate — fix or link one, or add "
                            "a constraint; don't trust their separate values.")
            banner.setWordWrap(True)
            banner.setStyleSheet("background:#c0392b; color:white; padding:4px; "
                                 "border-radius:3px; font-weight:600;")
            v.addWidget(banner)

        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        v.addWidget(bb)
