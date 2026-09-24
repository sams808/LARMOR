"""The help window (larmor/desktop/help_dialog.py): non-modal, one window per
open page, a text-size bar whose factor applies to every help window and
persists.

Sam: "Why can't we use LARMOR when the help page is open?" and "add buttons
to increase and decrease text size". The manuals render through
``mdrender.render_help_html(md, scale)``; ``toHtml`` writes the body font as
an absolute pt size that beats any stylesheet rule, so the zoom rewrites that
size and scales the equation images with it (see mdrender._scale_body_font).
"""
import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEvent, QPoint, QPointF, QSettings, Qt  # noqa: E402
from PySide6.QtGui import QWheelEvent  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _fresh_scale(monkeypatch):
    """Every test starts from an unread factor and leaves none behind."""
    from larmor.desktop import help_dialog

    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    monkeypatch.setattr(help_dialog, "_scale", None)
    yield
    for w in list(help_dialog._open):
        try:
            w.close()
        except RuntimeError:
            pass
    help_dialog._open.clear()
    help_dialog._scale = None


@pytest.fixture()
def win(qapp):
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    w.show()
    qapp.processEvents()
    yield w
    w.close()


def _flush(qapp):
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    qapp.processEvents()


def test_help_window_is_non_modal_and_the_main_window_stays_usable(win, qapp):
    from larmor.desktop.help_dialog import HelpWindow, show_help

    hw = show_help(win, "getting-started", "Getting started")
    qapp.processEvents()
    assert isinstance(hw, HelpWindow)
    assert not hw.isModal() and hw.isVisible()
    assert hw.windowFlags() & Qt.Window
    assert hw.windowFlags() & Qt.WindowMinMaxButtonsHint
    assert hw.windowFlags() & Qt.WindowCloseButtonHint
    assert hw.parent() is win
    assert hw.windowTitle() == "Getting started"
    assert "Getting started" in hw.tb.toPlainText()
    # a menu action of the main window while the help is open
    before = win.actResid.isChecked()
    win.actResid.trigger()
    assert win.actResid.isChecked() != before
    win.zoom_full()
    # listed in the open-windows bar like any tool window
    assert hw in win._window_bar.entries()


def test_menu_slots_go_through_show_help_and_reuse_the_open_page(win, qapp):
    from larmor.desktop.help_dialog import show_help

    win._open_manual("qcpmg", "QCPMG")
    qapp.processEvents()
    reg = win._help_windows
    assert set(reg) == {("manual", "qcpmg")}
    hw = reg[("manual", "qcpmg")]
    assert show_help(win, "qcpmg", "QCPMG") is hw               # raised, not doubled
    win._open_tutorial("01-first-fit-27Al-czjzek", "1 · A first fit")
    qapp.processEvents()
    assert set(reg) == {("manual", "qcpmg"), ("tutorial", "01-first-fit-27Al-czjzek")}
    assert len(win._window_bar.entries()) == 2
    hw.close()
    _flush(qapp)
    assert set(reg) == {("tutorial", "01-first-fit-27Al-czjzek")}
    # closed and freed: the same page opens as a new window
    hw2 = show_help(win, "qcpmg", "QCPMG")
    assert hw2 is not hw and reg[("manual", "qcpmg")] is hw2


def test_help_from_a_dialog_button_is_owned_by_that_dialog(win, qapp):
    """A tool dialog's Help button opens the page as a child of THAT window
    (so it is usable while the dialog is up and closes with it), while the
    main window keeps its own list."""
    from PySide6.QtWidgets import QDialog

    from larmor.desktop.help_dialog import show_help

    dlg = QDialog(win)
    dlg.show()
    qapp.processEvents()
    hw = show_help(dlg, "lineshapes", "Lineshapes")
    assert hw.parent() is dlg
    assert dlg._help_windows[("manual", "lineshapes")] is hw
    assert not getattr(win, "_help_windows", {})
    dlg.close()


def test_text_size_buttons_change_the_rendering_for_every_help_window(win, qapp):
    from larmor.desktop import help_dialog
    from larmor.desktop.help_dialog import ZOOM_STEPS, show_help

    hw = show_help(win, "lineshapes", "Lineshapes")   # equations inside
    qapp.processEvents()
    assert hw.scale == 1.0 and hw.lblScale.text() == "100 %"
    h0 = hw.tb.document().size().height()
    html0 = hw.html
    hw.actLarger.trigger()
    assert hw.scale == ZOOM_STEPS[ZOOM_STEPS.index(1.0) + 1]
    assert hw.lblScale.text() == f"{round(hw.scale * 100)} %"
    assert hw.tb.document().size().height() > h0             # text grew
    assert hw.html != html0
    # the typeset equations grew with it: the <img> widths are larger
    import re
    w0 = [int(x) for x in re.findall(r'<img [^>]*width="(\d+)"', html0)]
    w1 = [int(x) for x in re.findall(r'<img [^>]*width="(\d+)"', hw.html)]
    assert w0 and len(w0) == len(w1) and all(b > a for a, b in zip(w0, w1))
    # the factor is one per application: a second window opens at it
    other = show_help(win, "qcpmg", "QCPMG")
    assert other.scale == hw.scale == help_dialog.text_scale()
    # ... and follows every later change at once
    hw.actSmaller.trigger()
    hw.actSmaller.trigger()
    assert hw.scale == other.scale == ZOOM_STEPS[ZOOM_STEPS.index(1.0) - 1]
    other.actReset.trigger()
    assert hw.scale == other.scale == 1.0
    assert hw.tb.document().size().height() == pytest.approx(h0)
    # the shortcuts inside the window
    assert "Ctrl+=" in [s.toString() for s in hw.actLarger.shortcuts()]
    assert "Ctrl++" in [s.toString() for s in hw.actLarger.shortcuts()]
    assert "Ctrl+-" in [s.toString() for s in hw.actSmaller.shortcuts()]
    assert hw.actReset.shortcut().toString() == "Ctrl+0"
    # A- disabled at the smallest step, A+ at the largest
    help_dialog.set_text_scale(ZOOM_STEPS[0])
    assert not hw.actSmaller.isEnabled() and hw.actLarger.isEnabled()
    help_dialog.set_text_scale(ZOOM_STEPS[-1])
    assert hw.actSmaller.isEnabled() and not hw.actLarger.isEnabled()


def test_ctrl_wheel_zooms_through_the_same_factor(win, qapp):
    from larmor.desktop.help_dialog import show_help

    hw = show_help(win, "getting-started", "Getting started")
    qapp.processEvents()

    def wheel(dy):
        ev = QWheelEvent(QPointF(10, 10), QPointF(10, 10), QPoint(0, 0),
                         QPoint(0, dy), Qt.NoButton, Qt.ControlModifier,
                         Qt.ScrollUpdate, False)
        QApplication.sendEvent(hw.tb.viewport(), ev)

    wheel(120)
    assert hw.scale > 1.0
    wheel(-120)
    wheel(-120)
    assert hw.scale < 1.0
    # a plain wheel (no Ctrl) is left to the browser: no zoom
    s = hw.scale
    ev = QWheelEvent(QPointF(10, 10), QPointF(10, 10), QPoint(0, 0), QPoint(0, 120),
                     Qt.NoButton, Qt.NoModifier, Qt.ScrollUpdate, False)
    QApplication.sendEvent(hw.tb.viewport(), ev)
    assert hw.scale == s


def test_text_size_and_geometry_persist_in_settings(win, qapp, tmp_path, monkeypatch):
    from larmor.desktop import help_dialog
    from larmor.desktop.help_dialog import show_help

    ini = tmp_path / "larmor.ini"
    monkeypatch.setattr(help_dialog, "_settings",
                        lambda: QSettings(str(ini), QSettings.IniFormat))
    monkeypatch.delenv("LARMOR_NO_SESSION", raising=False)
    assert help_dialog.text_scale() == 1.0                    # nothing stored yet
    help_dialog.set_text_scale(1.3)
    assert float(QSettings(str(ini), QSettings.IniFormat).value("helpTextScale")) == 1.3
    help_dialog._scale = None                                 # a fresh launch
    assert help_dialog.text_scale() == 1.3
    hw = show_help(win, "qcpmg", "QCPMG")
    qapp.processEvents()
    assert hw.scale == 1.3
    hw.resize(640, 480)
    qapp.processEvents()
    hw.close()
    _flush(qapp)
    assert QSettings(str(ini), QSettings.IniFormat).value("helpWindowGeometry") is not None
    hw2 = show_help(win, "qcpmg", "QCPMG")
    qapp.processEvents()
    assert (hw2.width(), hw2.height()) == (640, 480)


def test_no_session_leaves_settings_untouched(win, qapp, tmp_path, monkeypatch):
    from larmor.desktop import help_dialog
    from larmor.desktop.help_dialog import show_help

    ini = tmp_path / "larmor.ini"
    monkeypatch.setattr(help_dialog, "_settings",
                        lambda: QSettings(str(ini), QSettings.IniFormat))
    hw = show_help(win, "qcpmg", "QCPMG")
    hw.zoom_in()
    hw.close()
    _flush(qapp)
    assert not ini.exists()


def test_render_help_html_scale_rewrites_body_font_and_image_sizes():
    from larmor.desktop.mdrender import HELP_CSS, help_css, render_help_html

    md = "# T\n\nbody text with $x^2$ inline\n\n$$E = mc^2$$\n"
    h1 = render_help_html(md, 1.0)
    h2 = render_help_html(md, 2.0)
    import re
    pt1 = float(re.search(r"<body[^>]*?font-size:([0-9.]+)pt", h1).group(1))
    pt2 = float(re.search(r"<body[^>]*?font-size:([0-9.]+)pt", h2).group(1))
    assert pt2 == pytest.approx(2 * pt1)
    w1 = [int(x) for x in re.findall(r'<img [^>]*width="(\d+)"', h1)]
    w2 = [int(x) for x in re.findall(r'<img [^>]*width="(\d+)"', h2)]
    assert len(w1) == len(w2) == 2 and all(b >= 2 * a - 1 for a, b in zip(w1, w2))
    assert help_css(1.0) == HELP_CSS
    assert "font-size: 28px" in help_css(2.0)                 # body 14px doubled
