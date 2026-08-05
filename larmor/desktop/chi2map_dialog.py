"""χ² map over a pair of fitted parameters — shows how well the data determines
them (a round basin) or whether they trade off (a diagonal valley)."""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QPushButton,
    QSpinBox, QVBoxLayout,
)

from larmor.desktop import theme
from larmor.chi2map import varying_params, chi2_surface


class Chi2MapDialog(QDialog):
    def __init__(self, parent, recipe, exp_ppm, exp_amp, window):
        super().__init__(parent)
        self.setWindowTitle("χ² map — parameter pair")
        self.resize(640, 620)
        self._recipe, self._ppm, self._amp, self._window = \
            recipe, exp_ppm, exp_amp, window
        self._params = varying_params(recipe)
        self._fig = None

        v = QVBoxLayout(self)
        row = QHBoxLayout()
        self.cbA = QComboBox(); self.cbB = QComboBox()
        for i, pn, label in self._params:
            self.cbA.addItem(label, (i, pn)); self.cbB.addItem(label, (i, pn))
        if self.cbB.count() > 1:
            self.cbB.setCurrentIndex(1)
        row.addWidget(QLabel("X:")); row.addWidget(self.cbA, 1)
        row.addWidget(QLabel("Y:")); row.addWidget(self.cbB, 1)
        row.addWidget(QLabel("grid:"))
        self.n = QSpinBox(); self.n.setRange(7, 41); self.n.setValue(15)
        row.addWidget(self.n)
        b_go = QPushButton("Compute"); b_go.clicked.connect(self._compute)
        row.addWidget(b_go)
        v.addLayout(row)

        self.img = QLabel("pick two parameters and Compute")
        self.img.setAlignment(Qt.AlignCenter)
        self.img.setMinimumSize(520, 460)
        self.img.setStyleSheet(f"background:{theme.active().plot_bg};")
        v.addWidget(self.img, 1)

        bb = QDialogButtonBox(QDialogButtonBox.Close)
        self.btnExp = bb.addButton("Export figure…", QDialogButtonBox.ActionRole)
        self.btnExp.setEnabled(False)
        self.btnExp.clicked.connect(self._export)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)

    def _compute(self):
        a = self.cbA.currentData(); b = self.cbB.currentData()
        if a is None or b is None or a == b:
            self.img.setText("choose two different parameters"); return
        self.setCursor(Qt.WaitCursor)
        try:
            A, B, Z, (a0, b0) = chi2_surface(
                self._recipe, self._ppm, self._amp, self._window, a, b,
                n=self.n.value())
        except Exception as exc:  # noqa: BLE001
            self.unsetCursor(); self.img.setText(f"failed: {exc}"); return
        finally:
            self.unsetCursor()
        self._fig = self._render(A, B, Z, a0, b0,
                                 self.cbA.currentText(), self.cbB.currentText())
        import io
        buf = io.BytesIO(); self._fig.savefig(buf, format="png", dpi=110,
                                              bbox_inches="tight")
        pix = QPixmap(); pix.loadFromData(buf.getvalue())
        self.img.setPixmap(pix.scaled(self.img.width(), self.img.height(),
                                      Qt.KeepAspectRatio, Qt.SmoothTransformation))
        self.btnExp.setEnabled(True)

    @staticmethod
    def _render(A, B, Z, a0, b0, xlabel, ylabel):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        cf = ax.contourf(A, B, Z, levels=18, cmap="viridis")
        ax.contour(A, B, Z, levels=10, colors="white", linewidths=0.4, alpha=0.5)
        ax.plot(a0, b0, "r+", ms=12, mew=2, label="fitted")
        fig.colorbar(cf, ax=ax, label="χ² (Σ residual²)")
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        return fig

    def _export(self):
        if self._fig is None:
            return
        from larmor.desktop.export_dialog import export_matplotlib
        export_matplotlib(self, self._fig, "chi2_map")
