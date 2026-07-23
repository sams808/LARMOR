"""Generic QThread wrapper: runs any zero-argument "give me a card"
callable off the Qt main thread, so a slow API can never freeze the
host application's window. Pack-agnostic -- doesn't know or care what
get_card actually fetches."""

from PySide6.QtCore import QThread, Signal


class FactWorker(QThread):
    ready = Signal(dict)

    def __init__(self, get_card, parent=None):
        super().__init__(parent)
        self._get_card = get_card

    def run(self):
        self.ready.emit(self._get_card())
