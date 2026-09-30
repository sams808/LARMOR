"""Pulse program viewer (Tools ▸ Pulse program…; Explorer ▸ right-click an
EXPNO ▸ Pulse program…): the sequence TopSpin compiled for an EXPNO as a
TopSpin-like timing diagram, as text, and as a table of the parameters it
uses with their acqus values.

The diagram is NOT to scale: every element has a fixed nominal width
(``larmor.pulseprog.timeline``), one row per rf channel in use plus the
receiver row, loop brackets underneath. Read-only and non-modal
(``show_tool_window``); nothing is written into the EXPNO. Parsing, value
resolution and the layout are Qt-free in ``larmor.pulseprog`` and tested
there; this module only paints.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, QSize, Qt
from PySide6.QtGui import (
    QBrush, QColor, QFont, QFontMetrics, QImage, QPainter, QPainterPath, QPen,
    QPolygonF,
)
from PySide6.QtWidgets import (
    QApplication, QDialog, QFileDialog, QHBoxLayout, QHeaderView, QLabel, QMenu,
    QPlainTextEdit, QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
    QTabWidget, QVBoxLayout, QWidget,
)

from larmor import pulseprog
from larmor.desktop import theme
from larmor.desktop.paths import suggest_save_dir

#: (name, title) of the manual the "?" button opens
HELP_PAGE = ("pulse-programs", "Reading a pulse program")

#: table order of the parameter categories
_CATEGORY_ORDER = ("pulse", "delay", "power", "shape", "decoupling", "loop",
                   "increment", "constant", "gradient", "phase", "other")


def _natural(sym: str):
    m = re.match(r"^([A-Za-z_]+?)(\d+)$", sym)
    return (m.group(1).lower(), int(m.group(2))) if m else (sym.lower(), -1)


class TimelineWidget(QWidget):
    """Paints a :class:`larmor.pulseprog.Timeline`: channel rows with the
    nucleus from acqus, pulses as filled boxes (name above, phase below),
    shaped pulses with an envelope, delays as labelled gaps, decoupling as
    a hatched bar, acquisition as a decaying FID, loops as brackets.
    ``image()`` renders it to a QImage and ``export()`` writes PNG or SVG."""

    UNIT = 52        # px per nominal width unit
    ROW_H = 88       # px per channel row
    LEFT = 100       # channel-label gutter
    TOP = 34         # title strip
    LOOP_H = 20      # px per loop-bracket level
    PAD = 16

    def __init__(self, parent=None):
        super().__init__(parent)
        self._tl: pulseprog.Timeline | None = None
        self._nuclei: dict[str, str] = {}
        self._title = ""
        self._notes: list[str] = []
        self._placeholder = "no pulse program to draw"
        self.setMinimumSize(420, 200)

    # ------------------------------------------------------------ data
    def set_timeline(self, tl: pulseprog.Timeline | None, nuclei: dict | None = None,
                     title: str = "", notes: list[str] | None = None):
        self._tl = tl
        self._nuclei = dict(nuclei or {})
        self._title = title
        self._notes = list(notes if notes is not None else (tl.notes if tl else []))
        w, h = self._size()
        self.setMinimumSize(w, h)
        self.updateGeometry()
        self.update()

    def set_placeholder(self, text: str):
        self._placeholder = text
        self.update()

    def timeline(self) -> pulseprog.Timeline | None:
        return self._tl

    def _size(self) -> tuple[int, int]:
        tl = self._tl
        if tl is None:
            return 420, 200
        rows = len(tl.channels)
        levels = (max(lp.level for lp in tl.loops) + 1) if tl.loops else 0
        w = self.LEFT + tl.total * self.UNIT + 2 * self.PAD
        h = (self.TOP + rows * self.ROW_H + 10 + levels * self.LOOP_H + 10
             + len(self._notes) * 15 + self.PAD)
        return int(math.ceil(w)), int(math.ceil(h))

    def sizeHint(self) -> QSize:
        return QSize(*self._size())

    # --------------------------------------------------------- painting
    def paintEvent(self, _event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setRenderHint(QPainter.TextAntialiasing)
        self.render_to(p, QRectF(self.rect()))
        p.end()

    def _fonts(self):
        base = QFont(self.font())
        ps = base.pointSizeF() if base.pointSizeF() > 0 else 9.0
        small = QFont(base)
        small.setPointSizeF(max(7.0, ps - 1))
        tiny = QFont(base)
        tiny.setPointSizeF(max(6.5, ps - 2))
        bold = QFont(base)
        bold.setBold(True)
        return base, small, tiny, bold

    @staticmethod
    def _text(p: QPainter, x: float, y: float, s: str, font: QFont, color: QColor,
              align: str = "center", max_w: float | None = None):
        if not s:
            return
        p.setFont(font)
        fm = QFontMetrics(font)
        if max_w is not None and max_w > 10:
            s = fm.elidedText(s, Qt.ElideRight, int(max_w))
        w = fm.horizontalAdvance(s)
        if align == "center":
            x = x - w / 2
        elif align == "right":
            x = x - w
        p.setPen(QPen(color))
        p.drawText(QPointF(x, y), s)

    def render_to(self, p: QPainter, rect: QRectF):
        t = theme.active()
        bg, ink, dim = QColor(t.plot_bg), QColor(t.text), QColor(t.text_dim)
        accent, border = QColor(t.accent), QColor(t.border)
        p.fillRect(rect, bg)
        base_f, small, tiny, bold = self._fonts()
        tl = self._tl
        if tl is None:
            self._text(p, rect.center().x(), rect.center().y(), self._placeholder, base_f, dim)
            return
        U, L, R = self.UNIT, self.LEFT, self.ROW_H

        def X(u: float) -> float:
            return L + u * U

        # title strip
        self._text(p, L, self.TOP - 14, self._title, bold, ink, align="left")
        self._text(p, X(tl.total), self.TOP - 14, "not to scale", tiny, dim, align="right")

        # rows and their gutter labels
        baselines: dict[str, float] = {}
        for i, ch in enumerate(tl.channels):
            base = self.TOP + i * R + R * 0.55
            baselines[ch] = base
            name = "acq" if ch == "acq" else ch
            sub = ("receiver" if ch == "acq"
                   else self._nuclei.get(ch, "gradient" if ch == "grad" else ""))
            self._text(p, 10, base - 3, name, bold, ink, align="left")
            self._text(p, 10, base + 13, sub, small, dim, align="left")
            p.setPen(QPen(border, 1))
            p.drawLine(QPointF(X(0), base), QPointF(X(tl.total), base))
            if i:
                sep = self.TOP + i * R - 2
                p.setPen(QPen(border, 0.5, Qt.DotLine))
                p.drawLine(QPointF(L - 4, sep), QPointF(X(tl.total) + 4, sep))

        first_rf = tl.channels[0] if tl.channels else "f1"
        h_pulse = R * 0.40
        for ev in tl.events:
            row = ev.channel or first_rf
            if row not in baselines:
                row = first_rf
            base = baselines[row]
            x0, x1 = X(ev.start), X(ev.end)
            cx, w = (x0 + x1) / 2, x1 - x0
            if ev.kind == "delay":
                self._text(p, cx, base - 5, ev.label, small, dim, max_w=w * 1.5)
                self._text(p, cx, base + 13, ev.value, small, dim, max_w=w * 1.5)
                self._text(p, cx, base + 26, ev.note, tiny, dim, max_w=w * 1.8)
            elif ev.kind in ("pulse", "shaped_pulse"):
                if ev.kind == "pulse":
                    r = QRectF(x0 + 3, base - h_pulse, max(w - 6, 4), h_pulse)
                    p.setPen(QPen(accent.darker(115), 1))
                    p.setBrush(QBrush(accent))
                    p.drawRect(r)
                else:
                    self._shape(p, x0 + 3, x1 - 3, base, h_pulse, ev.note, accent)
                if ev.cond:
                    self._text(p, cx, base - h_pulse - 18, f"#ifdef {ev.cond}", tiny, dim, max_w=w * 2)
                self._text(p, cx, base - h_pulse - 5, ev.label, bold, ink, max_w=w * 1.8)
                self._text(p, cx, base + 13, ev.sub, base_f, ink, max_w=w * 1.8)
                self._text(p, cx, base + 26, ev.value, small, dim, max_w=w * 1.8)
                self._text(p, cx, base + 38, ev.note, tiny, dim, max_w=w * 2.2)
            elif ev.kind == "set_power":
                p.setPen(QPen(dim, 1))
                p.drawLine(QPointF(cx, base), QPointF(cx, base - 10))
                self._text(p, cx, base - 13, ev.label, tiny, dim, max_w=w * 1.6)
                self._text(p, cx, base + 13, ev.value, tiny, dim, max_w=w * 1.6)
            elif ev.kind == "decouple":
                hb = R * 0.26
                r = QRectF(x0 + 2, base - hb, max(w - 4, 4), hb)
                hatch = QColor(accent)
                hatch.setAlpha(170)
                p.setPen(QPen(accent, 1))
                p.setBrush(QBrush(hatch, Qt.BDiagPattern))
                p.drawRect(r)
                self._text(p, cx, base - hb - 5, ev.label, small, ink, max_w=w - 4)
                self._text(p, cx, base + 13, ev.value, small, dim, max_w=w - 4)
            elif ev.kind == "acquire":
                amp = R * 0.34
                pts = []
                n = max(24, int(w))
                for k in range(n + 1):
                    s = k / n
                    y = base - amp * math.exp(-3.2 * s) * math.cos(2 * math.pi * 6.5 * s)
                    pts.append(QPointF(x0 + 4 + s * (w - 8), y))
                p.setPen(QPen(accent, 1.5))
                p.setBrush(Qt.NoBrush)
                p.drawPolyline(QPolygonF(pts))
                self._text(p, cx, base - amp - 6, ev.label, bold, ink, max_w=w)
                self._text(p, cx, base + 13, ev.sub, base_f, ink, max_w=w)
                self._text(p, cx, base + 26, ev.value, small, dim, max_w=w)
                if ev.cond:
                    self._text(p, cx, base + 38, f"#ifdef {ev.cond}", tiny, dim, max_w=w)

        # loop brackets, innermost nearest the rows
        y0 = self.TOP + len(tl.channels) * R + 8
        for lp in tl.loops:
            y = y0 + lp.level * self.LOOP_H + 4
            xa, xb = X(lp.start) + 2, X(lp.end) - 2
            p.setPen(QPen(dim, 1.2))
            p.drawLine(QPointF(xa, y), QPointF(xb, y))
            p.drawLine(QPointF(xa, y), QPointF(xa, y - 6))
            p.drawLine(QPointF(xb, y), QPointF(xb, y - 6))
            self._text(p, (xa + xb) / 2, y + 12, lp.text or f"× {lp.times}", tiny, dim,
                       max_w=max(xb - xa, 40))
        levels = (max(lp.level for lp in tl.loops) + 1) if tl.loops else 0
        yn = y0 + levels * self.LOOP_H + 16
        note_f = QFont(tiny)
        note_f.setItalic(True)
        for k, note in enumerate(self._notes):
            self._text(p, L, yn + k * 15, note, note_f, dim, align="left")

    def _shape(self, p: QPainter, x0: float, x1: float, base: float, h: float,
               note: str, accent: QColor):
        """A shaped pulse: a ramp for ramp.* shapes, a swept oscillation
        under a bell for sweeps (dfs / wurst / chirp), a Gaussian otherwise."""
        shape = note.split(" · ", 1)[0].lower() if " · " in note else ""
        w = max(x1 - x0, 6)
        n = max(30, int(w))
        path = QPainterPath(QPointF(x0, base))
        for k in range(n + 1):
            s = k / n
            if shape.startswith("ramp"):
                env = 0.35 + 0.65 * s
            else:
                env = math.exp(-((s - 0.5) / 0.22) ** 2)
            path.lineTo(QPointF(x0 + s * w, base - h * env))
        path.lineTo(QPointF(x1, base))
        path.closeSubpath()
        fill = QColor(accent)
        fill.setAlpha(110)
        p.setPen(QPen(accent, 1.3))
        p.setBrush(QBrush(fill))
        p.drawPath(path)
        if any(key in shape for key in ("dfs", "wurst", "chirp", "sweep", "hs", "tanh")):
            pts = []
            for k in range(n + 1):
                s = k / n
                env = math.exp(-((s - 0.5) / 0.22) ** 2)
                pts.append(QPointF(x0 + s * w, base - h * env * 0.5
                                   * (1 + math.cos(2 * math.pi * (2.0 * s + 6.0 * s * s)))))
            p.setPen(QPen(accent.darker(130), 0.9))
            p.setBrush(Qt.NoBrush)
            p.drawPolyline(QPolygonF(pts))

    # ----------------------------------------------------------- export
    def image(self, scale: float = 2.0) -> QImage:
        w, h = self._size()
        img = QImage(int(w * scale), int(h * scale), QImage.Format_ARGB32)
        img.fill(QColor(theme.active().plot_bg))
        p = QPainter(img)
        p.setRenderHint(QPainter.Antialiasing)
        p.setRenderHint(QPainter.TextAntialiasing)
        p.scale(scale, scale)
        self.render_to(p, QRectF(0, 0, w, h))
        p.end()
        return img

    def export(self, path) -> str:
        """Write the diagram as PNG, or as SVG for a ``.svg`` path (PNG next
        to it when QtSvg is not available). Returns the path written."""
        path = str(path)
        if path.lower().endswith(".svg"):
            try:
                from PySide6.QtSvg import QSvgGenerator
            except ImportError:
                path = path[:-4] + ".png"
            else:
                w, h = self._size()
                gen = QSvgGenerator()
                gen.setFileName(path)
                gen.setSize(QSize(w, h))
                gen.setViewBox(QRectF(0, 0, w, h))
                gen.setTitle(self._title or "pulse program")
                gen.setDescription("LARMOR pulse program diagram (not to scale)")
                p = QPainter(gen)
                self.render_to(p, QRectF(0, 0, w, h))
                p.end()
                return path
        if not self.image().save(path):
            raise OSError(f"could not write {path}")
        return path


class PulseProgramDialog(QDialog):
    """``expno_path``: an EXPNO folder (or anything under it -- a pdata
    folder, a fid, a 1r). Three tabs: Diagram, Text, Parameters; a "?"
    button opens the manual page; Export… writes PNG / SVG."""

    def __init__(self, parent=None, expno_path: str = ""):
        super().__init__(parent)
        self.expno_path = str(expno_path or "")
        self.record: pulseprog.Record | None = (
            pulseprog.load(self.expno_path) if self.expno_path else None)
        rec, prog = self.record, (self.record.program if self.record else None)
        self.name = (rec.pulprog if rec else "") or (prog.name if prog else "") or "pulse program"
        folder = rec.expno.name if rec else (Path(self.expno_path).name if self.expno_path else "")
        self.setWindowTitle(f"Pulse program — {self.name}"
                            + (f"  (EXPNO {folder})" if folder else ""))
        self.setModal(False)
        self.resize(1000, 660)
        t = theme.active()
        v = QVBoxLayout(self)

        self.header = QLabel()
        self.header.setWordWrap(True)
        self.header.setTextFormat(Qt.RichText)
        self.header.setTextInteractionFlags(Qt.TextSelectableByMouse)
        v.addWidget(self.header)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet(f"color: {t.text_dim};")
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        v.addWidget(self.summary)

        self.tabs = QTabWidget()
        self.diagram = TimelineWidget()
        self.diagram.setContextMenuPolicy(Qt.CustomContextMenu)
        self.diagram.customContextMenuRequested.connect(self._diagram_menu)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.diagram)
        self.tabs.addTab(self.scroll, "Diagram")
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        self.text.setLineWrapMode(QPlainTextEdit.NoWrap)
        f = self.text.font()
        f.setFamily("Consolas")
        f.setStyleHint(QFont.Monospace)
        self.text.setFont(f)
        self.tabs.addTab(self.text, "Text")
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["symbol", "value", "unit", "comment"])
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.Stretch)
        self.tabs.addTab(self.table, "Parameters")
        self.tabs.setTabToolTip(0, "one row per rf channel plus the receiver; widths are nominal, "
                                   "not durations")
        self.tabs.setTabToolTip(1, "the pulseprogram file as TopSpin compiled it, cpp line "
                                   "markers removed")
        self.tabs.setTabToolTip(2, "every symbol the body uses, with its value from acqus or "
                                   "from the program's own definitions")
        v.addWidget(self.tabs, 1)

        bottom = QHBoxLayout()
        self.notes = QLabel()
        self.notes.setWordWrap(True)
        self.notes.setStyleSheet(f"color: {t.text_dim};")
        bottom.addWidget(self.notes, 1)
        self.btnExport = QPushButton("Export…")
        self.btnExport.setToolTip("save the diagram as a PNG image or an SVG drawing")
        self.btnExport.clicked.connect(self._export)
        self.btnHelp = QPushButton("?")
        self.btnHelp.setFixedWidth(32)
        self.btnHelp.setToolTip("Reading a pulse program — the help page: syntax, the diagram, "
                                "worked examples")
        self.btnHelp.clicked.connect(self._help)
        btnClose = QPushButton("Close")
        btnClose.clicked.connect(self.close)
        for w in (self.btnExport, self.btnHelp, btnClose):
            bottom.addWidget(w)
        v.addLayout(bottom)
        for b in self.findChildren(QPushButton):
            b.setAutoDefault(False)
        self._fill()

    # ------------------------------------------------------------- fill
    def _fill(self):
        rec = self.record
        prog = rec.program if rec else None
        if rec is None or prog is None:
            why = "; ".join(rec.warnings) if rec and rec.warnings else "no EXPNO given"
            self.header.setText(f"<b>{self.name}</b> &nbsp;·&nbsp; {self.expno_path or '—'}")
            self.summary.setText(f"{why} — TopSpin writes the compiled sequence as "
                                 f"'pulseprogram' next to acqus at acquisition time.")
            self.diagram.set_placeholder(why)
            self.text.setPlainText(f"; {why}\n")
            self.notes.setText("")
            self.btnExport.setEnabled(False)
            return
        nuclei = " · ".join(f"{ch} {nuc}" for ch, nuc in rec.nuclei.items())
        comment = prog.meta.get("COMMENT", "").strip()
        head = [f"<b>{self.name}</b> &nbsp;·&nbsp; {rec.path}"]
        line2 = " &nbsp;·&nbsp; ".join(x for x in (nuclei, comment) if x)
        if line2:
            head.append(line2)
        self.header.setText("<br>".join(head))
        self.summary.setText(pulseprog.describe(prog, rec.values))
        tl = pulseprog.timeline(prog, rec.values)
        title = f"{self.name} — {Path(rec.expno).parent.name}/{rec.expno.name}"
        self.diagram.set_timeline(tl, rec.nuclei, title)
        self.text.setPlainText(prog.text)
        self._fill_table(prog, rec.values)
        notes = list(rec.warnings) + list(tl.notes)
        self.notes.setText("  ·  ".join(notes))

    def _fill_table(self, prog: pulseprog.Program, values: dict):
        syms = list(values.keys())
        rank = {c: i for i, c in enumerate(_CATEGORY_ORDER)}
        syms.sort(key=lambda s: (rank.get(pulseprog.category(s), len(rank)), _natural(s)))
        self.table.setRowCount(len(syms))
        for r, sym in enumerate(syms):
            val = values[sym]
            legend = prog.legend_for(sym)
            comment = " — ".join(x for x in (val.note, legend) if x)
            if not comment and sym in prog.definitions:
                comment = f"= {prog.definitions[sym]}"
            cells = [sym, val.text, "" if val.unit in ("count", "list") else val.unit, comment]
            for c, s in enumerate(cells):
                it = QTableWidgetItem(s)
                if c == 1:
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                if c == 0 and val.source:
                    it.setToolTip(f"from {val.source}")
                self.table.setItem(r, c, it)

    # ---------------------------------------------------------- actions
    def export_to(self, path) -> str:
        """Write the diagram to ``path`` (PNG, or SVG for .svg); returns the
        path written."""
        return self.diagram.export(path)

    def _export(self):
        if self.record is None or self.record.program is None:
            return
        base = re.sub(r"[^\w.-]+", "_", self.name) or "pulseprogram"
        start = suggest_save_dir(self.expno_path, "")
        seed = str(Path(start) / f"{base}_diagram.png") if start else f"{base}_diagram.png"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export the diagram", seed, "PNG image (*.png);;SVG drawing (*.svg)")
        if not path:
            return
        try:
            out = self.export_to(path)
        except OSError as exc:
            self.notes.setText(f"export failed: {exc}")
            return
        self.notes.setText(f"saved {out}")

    def _copy_image(self):
        QApplication.clipboard().setImage(self.diagram.image())
        self.notes.setText("diagram copied to the clipboard")

    def _help(self):
        from larmor.desktop.help_dialog import show_help

        show_help(self, *HELP_PAGE)

    def _diagram_menu(self, pos):
        m = QMenu(self)
        m.addAction("Copy image", self._copy_image)
        m.addAction("Export…  (PNG or SVG)", self._export)
        m.addSeparator()
        m.addAction("Reading a pulse program  (help)", self._help)
        m.exec(self.diagram.mapToGlobal(pos))
