"""The spinning-sideband offer banner on its own (a bare SpectrumView, no
MainWindow): placement, wording, comb guides, signals, dismissal."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
# the banner never touches QSettings, but MainWindow-free widget tests state
# the precondition the desktop suite relies on when run alone
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def view(qapp):
    from larmor.desktop.plot import SpectrumView

    v = SpectrumView()
    v.resize(800, 500)
    yield v
    v.close()


def _det():
    """A detection on a synthetic 1H manifold: 1 kHz at 100 MHz = 10 ppm."""
    from larmor import engine, sidebands
    from larmor.recipe import Param, Recipe, SiteModel

    r = Recipe(nucleus="1H", larmor_frequency_MHz=100.0, spin_rate_Hz=1000.0,
               sites=[SiteModel(model="sidebands", label="S", params={
                   "isotropic_chemical_shift_ppm": Param(0.0),
                   "shift_fwhm_ppm": Param(1.0), "amplitude": Param(1.0),
                   "ssb_ratio": Param(0.4), "n_ssb": Param(3.0),
                   "gl": Param(1.0)})])
    x = np.linspace(-50, 50, 4000)
    xc, y, _ = engine.simulate(r, exp_ppm=x)
    det = sidebands.detect(xc, y, 100.0, nu_rot_Hz=1000.0)
    assert det.ok
    return det


def test_banner_sits_top_centre_and_follows_resize(qapp, view):
    from larmor.desktop.sideband_offer import BANNER_TOP, SidebandBanner

    view.show()
    b = SidebandBanner(view)
    assert b.isHidden()                      # nothing to offer yet
    b.set_detection(_det(), 1000.0)
    b.show()
    qapp.processEvents()
    assert b.parent() is view
    assert abs(b.x() + b.width() / 2 - 400) <= 2
    assert b.y() == BANNER_TOP == 36
    assert b.width() < 800 and b.height() < 200      # compact, not a sheet
    view.resize(1000, 500)
    qapp.processEvents()
    assert abs(b.x() + b.width() / 2 - 500) <= 2
    assert b.y() == BANNER_TOP


def test_banner_text_guides_and_signals(qapp, view):
    import pyqtgraph as pg
    from larmor.desktop.sideband_offer import SidebandBanner

    det = _det()
    b = SidebandBanner(view)
    b.set_detection(det, recipe_rate_Hz=1100.0)
    text = b.label.text()
    assert "⚠" in text and "1 100" in text and "1 000 Hz" in text
    assert "10.0 ppm" in text and "100.00 MHz" in text
    b.set_detection(det, recipe_rate_Hz=1000.0)
    assert "⚠" not in b.label.text()
    assert b.det is det

    b.show_guides(det)
    assert len(b.guides) == 1 + len(det.matched()) == 7
    items = view.getPlotItem().items
    for g in b.guides:
        assert isinstance(g, pg.InfiniteLine) and g in items
        assert not g.movable
    labels = sorted(g.label.textItem.toPlainText() for g in b.guides)
    assert labels == sorted(["centre", "+1", "−1", "+2", "−2", "+3", "−3"])
    positions = sorted(round(float(g.value()), 3) for g in b.guides)
    assert positions == [-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0]

    got = []
    b.manifold.connect(lambda: got.append("manifold"))
    b.model_line.connect(lambda: got.append("model"))
    b.copy.connect(lambda: got.append("copy"))
    b.dismissed.connect(lambda: got.append("dismissed"))
    b.btnManifold.click()
    b.actModelLine.trigger()
    b.btnCopy.click()
    b.btnClose.click()
    assert got == ["manifold", "model", "copy", "dismissed"]
    got.clear()
    QTest.keyClick(b, Qt.Key_Return)
    QTest.keyClick(b, Qt.Key_Escape)
    assert got == ["manifold", "dismissed"]
    assert b.actModelLine in b.btnManifold.menu().actions()
    assert "sidebands" in b.actModelLine.text()

    b.dismiss()
    assert b.isHidden() and b.guides == [] and b.det is None
    leftover = [it for it in view.getPlotItem().items
                if isinstance(it, pg.InfiniteLine)
                and getattr(it, "label", None) is not None
                and it.label.textItem.toPlainText() in ("centre", "+1", "−1")]
    assert leftover == []


def test_banner_does_not_swallow_view_events(qapp, view):
    """The event filter only repositions on Resize; the view still receives
    everything (a swallowed event would break drag-and-drop and clicks)."""
    from PySide6.QtCore import QEvent
    from larmor.desktop.sideband_offer import SidebandBanner

    b = SidebandBanner(view)
    assert b.eventFilter(view, QEvent(QEvent.MouseButtonPress)) is False
    assert b.eventFilter(view, QEvent(QEvent.Resize)) is False
