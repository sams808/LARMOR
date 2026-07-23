"""The popup itself: a small, self-contained QDialog styled to look
like a specimen card. Pack-agnostic -- takes any zero-argument
"give me a card dict" callable and doesn't know or care which pack
(or which of several packs, chosen at random) it came from."""

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QDialog, QFrame, QGridLayout, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget,
)

from .worker import FactWorker

INK = "#171A21"
INK_SOFT = "#4A4F58"
INK_FAINT = "#7B818C"
RULE = "#D7DAD9"
BRAND = "#1F3A5F"
BG_CARD = "#FFFFFF"
BG_PAGE = "#EDEFEE"

IMG_W, IMG_H = 380, 210


def _dot(color: str, size: int = 10) -> QLabel:
    lab = QLabel()
    lab.setFixedSize(size, size)
    lab.setStyleSheet(f"background:{color}; border-radius:{size // 2}px;")
    return lab


def _scaled_cover(data: bytes, w: int, h: int) -> QPixmap:
    pix = QPixmap()
    pix.loadFromData(data)
    if pix.isNull():
        return pix
    scaled = pix.scaled(w, h, Qt.KeepAspectRatioByExpanding, Qt.SmoothTransformation)
    x = max(0, (scaled.width() - w) // 2)
    y = max(0, (scaled.height() - h) // 2)
    return scaled.copy(x, y, w, h)


class FactDialog(QDialog):
    def __init__(self, get_card, parent=None):
        super().__init__(parent)
        self._get_card = get_card
        self.setWindowTitle("Fact")
        self.setFixedWidth(IMG_W + 40)
        self.setStyleSheet(f"QDialog {{ background:{BG_PAGE}; }}")
        self._worker = None
        self._build_skeleton()
        self.refresh()

    # ---- layout skeleton, populated by _apply_card -------------------

    def _build_skeleton(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(16, 16, 16, 16)
        outer.setSpacing(10)

        self.card_frame = QFrame()
        self.card_frame.setStyleSheet(
            f"QFrame {{ background:{BG_CARD}; border:1px solid {RULE}; border-top:3px solid {BRAND}; }}")
        card = QVBoxLayout(self.card_frame)
        card.setContentsMargins(0, 0, 0, 0)
        card.setSpacing(0)

        head = QHBoxLayout()
        head.setContentsMargins(14, 12, 14, 0)
        self.catalog_lab = QLabel()
        self.catalog_lab.setStyleSheet(f"color:{INK_FAINT}; font-family:Consolas; font-size:10px; letter-spacing:1px;")
        self.class_lab = QLabel()
        self.class_lab.setStyleSheet(
            f"color:{BRAND}; background:rgba(31,58,95,0.08); font-family:Consolas; "
            f"font-size:10px; padding:3px 8px; letter-spacing:1px;")
        head.addWidget(self.catalog_lab)
        head.addStretch(1)
        head.addWidget(self.class_lab)
        card.addLayout(head)

        self.name_lab = QLabel()
        self.name_lab.setStyleSheet(f"color:{INK}; font-family:Cambria,Georgia,serif; font-size:21px; padding:8px 14px 0;")
        self.name_lab.setWordWrap(True)
        card.addWidget(self.name_lab)

        self.subtitle_lab = QLabel()
        self.subtitle_lab.setStyleSheet(f"color:{INK_SOFT}; font-size:12px; padding:2px 14px 12px;")
        self.subtitle_lab.setWordWrap(True)
        card.addWidget(self.subtitle_lab)

        self.image_lab = QLabel()
        self.image_lab.setFixedSize(IMG_W, IMG_H)
        self.image_lab.setAlignment(Qt.AlignCenter)
        self.image_lab.setStyleSheet(f"background:{RULE}; border-top:1px solid {RULE}; border-bottom:1px solid {RULE};")
        card.addWidget(self.image_lab)

        self.caption_lab = QLabel()
        self.caption_lab.setStyleSheet(
            f"color:{INK_FAINT}; font-family:Consolas; font-size:10px; padding:6px 14px; "
            f"background:{BG_CARD}; border-bottom:1px solid {RULE};")
        self.caption_lab.setWordWrap(True)
        card.addWidget(self.caption_lab)

        headline_row = QHBoxLayout()
        headline_row.setContentsMargins(14, 12, 14, 4)
        headline_row.setSpacing(8)
        self.headline_dot = _dot(INK_FAINT)
        self.headline_lab = QLabel()
        self.headline_lab.setStyleSheet(f"color:{INK}; font-family:Consolas; font-size:12px;")
        headline_row.addWidget(self.headline_dot)
        headline_row.addWidget(self.headline_lab, 1)
        card.addLayout(headline_row)

        self.grid_widget = QWidget()
        self.grid_layout = QGridLayout(self.grid_widget)
        self.grid_layout.setContentsMargins(14, 6, 14, 6)
        self.grid_layout.setHorizontalSpacing(16)
        self.grid_layout.setVerticalSpacing(8)
        card.addWidget(self.grid_widget)

        self.foot_lab = QLabel()
        self.foot_lab.setStyleSheet(f"color:{INK_FAINT}; font-size:11px; padding:10px 14px 14px; border-top:1px solid {RULE};")
        self.foot_lab.setWordWrap(True)
        card.addWidget(self.foot_lab)

        outer.addWidget(self.card_frame)

        controls = QHBoxLayout()
        self.shuffle_btn = QPushButton("\U0001F3B2 Another one")
        self.shuffle_btn.setStyleSheet(
            f"QPushButton {{ background:{BRAND}; color:white; border:none; padding:9px 18px; font-weight:600; }}"
            f"QPushButton:disabled {{ background:{INK_FAINT}; }}")
        self.shuffle_btn.clicked.connect(self.refresh)
        controls.addWidget(self.shuffle_btn)
        controls.addStretch(1)
        outer.addLayout(controls)

    # ---- data population ----------------------------------------------

    def refresh(self):
        self.shuffle_btn.setEnabled(False)
        self.shuffle_btn.setText("tuning in…")
        self.name_lab.setText("\U0001F50E Looking for a fact…")
        self.subtitle_lab.setText("")
        self.image_lab.clear()
        self.caption_lab.setText("")
        self.headline_lab.setText("")
        self.foot_lab.setText("")
        self._clear_grid()

        self._worker = FactWorker(self._get_card, self)
        self._worker.ready.connect(self._apply_card)
        self._worker.start()

    def _clear_grid(self):
        while self.grid_layout.count():
            item = self.grid_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

    def _apply_card(self, card: dict):
        self.shuffle_btn.setEnabled(True)
        self.shuffle_btn.setText("\U0001F3B2 Another one")

        tier_suffix = "  ·  offline sample" if card.get("tier") == "offline" else ""
        self.catalog_lab.setText(card["catalog_tag"] + tier_suffix)
        self.class_lab.setText(card["class_tag"].upper())
        self.name_lab.setText(card["name"])
        self.subtitle_lab.setText(card["subtitle"])

        data = card.get("image_bytes")
        if data:
            self.image_lab.setPixmap(_scaled_cover(data, IMG_W, IMG_H))
        self.caption_lab.setText(card.get("image_caption", ""))

        h = card["headline"]
        if h.get("unknown"):
            self.headline_dot.setStyleSheet("background:transparent;")
            self.headline_lab.setStyleSheet(f"color:{INK_FAINT}; font-family:Consolas; font-size:12px; font-style:italic;")
            self.headline_lab.setText(h["unknown"])
        else:
            color = h.get("dot", INK_FAINT)
            self.headline_dot.setStyleSheet(f"background:{color}; border-radius:5px;")
            self.headline_lab.setStyleSheet(f"color:{INK}; font-family:Consolas; font-size:12px;")
            self.headline_lab.setText(f"{h.get('text', '')}   ·  {h.get('tail', '')}")

        self._clear_grid()
        for row, (label, value, unit) in enumerate(card["grid"]):
            r, c = divmod(row, 2)
            lab_w = QLabel(label.upper())
            lab_w.setStyleSheet(f"color:{INK_FAINT}; font-family:Consolas; font-size:9px; letter-spacing:0.5px;")
            if unit == "na":
                val_w = QLabel("not on file")
                val_w.setStyleSheet(f"color:{INK_FAINT}; font-size:12px; font-style:italic;")
            else:
                text = f"{value}  {unit}".strip() if unit else str(value)
                val_w = QLabel(text)
                val_w.setStyleSheet(f"color:{INK}; font-size:13px;")
            val_w.setWordWrap(True)
            cell = QVBoxLayout()
            cell.setSpacing(2)
            cell.addWidget(lab_w)
            cell.addWidget(val_w)
            holder = QWidget()
            holder.setLayout(cell)
            self.grid_layout.addWidget(holder, r, c)

        self.foot_lab.setText(card["foot"])
