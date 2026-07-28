"""Show the Czjzek C_Q distribution implied by each fitted σ (idea #2).
Turns the width σ into the physical P(C_Q) it stands for; see larmor.czjzek_dist
and the Lineshapes manual (Czjzek section)."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QLabel, QVBoxLayout,
)

from larmor.czjzek_dist import marginal_cq, rms_pq, suggested_cq_axis
from larmor.desktop.plot import site_color


class CzjzekDistDialog(QDialog):
    def __init__(self, parent, recipe: dict):
        super().__init__(parent)
        self.setWindowTitle("Czjzek distribution — P(C_Q)")
        self.resize(660, 480)
        v = QVBoxLayout(self)

        sites = [(i, s) for i, s in enumerate(recipe.get("sites", []))
                 if s.get("model") in ("czjzek", "ext_czjzek")
                 and "sigma_Cq_MHz" in s.get("params", {})]
        v.addWidget(QLabel(
            "The Czjzek width σ stands for a whole distribution of quadrupolar "
            "couplings, P(C_Q) ∝ C_Q⁴·⟨η terms⟩·exp(−C_Q²/2σ²) — it peaks near "
            "C_Q = 2σ. √⟨P_Q²⟩ = √5·σ is the invariant to report for a glass."))
        v.itemAt(0).widget().setWordWrap(True)

        plot = pg.PlotWidget(background="#fcfdfc")
        plot.setLabel("bottom", "C_Q", units="MHz")
        plot.setLabel("left", "P(C_Q)")
        plot.showGrid(x=True, y=True, alpha=0.15)
        plot.addLegend(offset=(-10, 10))
        v.addWidget(plot, 1)

        lines = []
        smax = max((s["params"]["sigma_Cq_MHz"]["value"] for _, s in sites),
                   default=1.0)
        cq = suggested_cq_axis(smax, 400)
        for i, s in sites:
            sigma = float(s["params"]["sigma_Cq_MHz"]["value"])
            p = marginal_cq(sigma, cq)
            name = s.get("label") or f"site {i}"
            plot.plot(cq, p, pen=pg.mkPen(site_color(i), width=1.8), name=name)
            # mark the mode (2σ) and √⟨P_Q²⟩
            plot.addItem(pg.InfiniteLine(
                pos=2 * sigma, angle=90,
                pen=pg.mkPen(site_color(i), width=1, style=Qt.DotLine)))
            lines.append(f"<b style='color:{site_color(i)}'>{name}</b>: "
                         f"σ={sigma:.3g} MHz · mode C_Q≈{2*sigma:.3g} · "
                         f"√⟨P_Q²⟩={rms_pq(sigma):.3g} MHz")

        summary = QLabel("<br>".join(lines) if lines
                         else "no Czjzek sites in the current fit")
        summary.setTextFormat(Qt.RichText)
        v.addWidget(summary)

        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(self.reject); bb.accepted.connect(self.accept)
        bb.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        v.addWidget(bb)
