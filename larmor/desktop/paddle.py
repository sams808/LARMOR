"""dmfit-style paddle: the on-spectrum manipulator for one line.

Structure (mirrors dmfit's 'Show Paddle'):
  - a vertical stem at the line position,
  - a square TOP handle at (position, amplitude): drag it to move the line in
    BOTH position and amplitude at once,
  - two small side handles at half height: drag either to change the width
    symmetrically.
"""
from __future__ import annotations

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPen
import pyqtgraph as pg


class _Handle(pg.GraphicsObject):
    """A small draggable square in view coordinates."""

    def __init__(self, paddle, kind: str, color: QColor, size_px: int = 9):
        super().__init__()
        self.paddle, self.kind, self.color, self.size_px = paddle, kind, color, size_px
        self.setFlag(self.GraphicsItemFlag.ItemIgnoresTransformations)
        self.setCursor(Qt.SizeAllCursor if kind == "top" else Qt.SizeHorCursor)
        self.setZValue(30)

    def boundingRect(self) -> QRectF:
        s = self.size_px
        return QRectF(-s, -s, 2 * s, 2 * s)

    def paint(self, painter, *args):
        s = self.size_px / 2 + 1.5
        painter.setPen(QPen(self.color, 1.4))
        painter.setBrush(QColor(255, 255, 255, 220))
        if self.kind == "top":
            painter.drawRect(QRectF(-s, -s, 2 * s, 2 * s))
        else:
            painter.drawEllipse(QRectF(-s, -s, 2 * s, 2 * s))

    def mouseDragEvent(self, ev):
        ev.accept()
        vb = self.paddle.getViewBox()
        if vb is None:
            return
        p = vb.mapSceneToView(ev.scenePos())
        if ev.isStart():
            self.paddle._drag_started()
        self.paddle._handle_dragged(self.kind, p)
        if ev.isFinish():
            self.paddle._drag_finished()


class Paddle(pg.GraphicsObject):
    """Manipulator for one line: emits live updates while dragged."""

    moved = Signal(int, float, float, float)     # index, pos, amp, fwhm (live)
    released = Signal(int)                       # drag finished -> snapshot

    def __init__(self, index: int, color: str,
                 pos: float, amp: float, fwhm: float, movable: bool = True):
        super().__init__()
        self.index = index
        self.color = QColor(color)
        self._pos, self._amp, self._fwhm = float(pos), float(amp), float(fwhm)
        self.movable = movable
        self.setZValue(20)

        self.h_top = _Handle(self, "top", self.color)
        self.h_left = _Handle(self, "left", self.color, 7)
        self.h_right = _Handle(self, "right", self.color, 7)
        for hnd in (self.h_top, self.h_left, self.h_right):
            hnd.setParentItem(self)
            hnd.setVisible(movable)
        self._layout_handles()

    # ------------------------------------------------------------- geometry
    def boundingRect(self) -> QRectF:
        w = max(self._fwhm, 1e-6)
        h = max(abs(self._amp), 1e-12)
        return QRectF(self._pos - w, min(0.0, self._amp) - 0.05 * h,
                      2 * w, 1.2 * h)

    def paint(self, painter, *args):
        pen = QPen(self.color, 0)
        pen.setStyle(Qt.DashLine)
        pen.setCosmetic(True)
        painter.setPen(pen)
        # stem
        painter.drawLine(QPointF(self._pos, 0.0), QPointF(self._pos, self._amp))
        # width bar at half height
        y2 = self._amp / 2.0
        painter.drawLine(QPointF(self._pos - self._fwhm / 2.0, y2),
                         QPointF(self._pos + self._fwhm / 2.0, y2))

    def _layout_handles(self):
        self.h_top.setPos(self._pos, self._amp)
        y2 = self._amp / 2.0
        self.h_left.setPos(self._pos - self._fwhm / 2.0, y2)
        self.h_right.setPos(self._pos + self._fwhm / 2.0, y2)

    # ------------------------------------------------------------- dragging
    def _drag_started(self):
        pass

    def _handle_dragged(self, kind: str, p):
        if not self.movable:
            return
        if kind == "top":
            self._pos = float(p.x())
            self._amp = max(float(p.y()), 1e-12)
        else:
            half = abs(float(p.x()) - self._pos)
            self._fwhm = max(2.0 * half, 1e-3)
        self.prepareGeometryChange()
        self._layout_handles()
        self.update()
        self.moved.emit(self.index, self._pos, self._amp, self._fwhm)

    def _drag_finished(self):
        self.released.emit(self.index)

    def set_state(self, pos: float, amp: float, fwhm: float):
        self._pos, self._amp, self._fwhm = float(pos), float(amp), float(fwhm)
        self.prepareGeometryChange()
        self._layout_handles()
        self.update()
