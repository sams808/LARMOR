"""FitHealthStrip alone (no MainWindow): rendering, click-through signals, the
details menu, stale dimming, the size policy and colour contrast."""
import os
import types

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtGui import QAction  # noqa: E402
from PySide6.QtWidgets import QApplication, QSizePolicy  # noqa: E402

from larmor import fithealth  # noqa: E402
from larmor.desktop import theme  # noqa: E402
from larmor.recipe import Param, Recipe, SiteModel  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


X = np.linspace(-50, 150, 400)
PEAK = np.exp(-0.5 * ((X - 60) / 8) ** 2)
Y = PEAK + np.random.RandomState(0).normal(0, 0.02, 400)
WINDOW = (40.0, -10.0)


def _rec():
    p = {"isotropic_chemical_shift_ppm": Param(15.0), "shift_fwhm_ppm": Param(6.0),
         "amplitude": Param(100.0), "gl": Param(0.5, vary=False)}
    return Recipe(nucleus="11B", larmor_frequency_MHz=160.0,
                  sites=[SiteModel(model="gauss_lor", label="A", params=p)])


def _ok():
    lm = types.SimpleNamespace(var_names=["s0_amplitude"], covar=np.eye(1), redchi=1.03)
    return fithealth.assess(_rec(), Y, PEAK, lmfit_result=lm, rmsd=0.0462, window=WINDOW)


def _flag(kind, level, text, target, params=(), stale=False, cov=False, detail=None):
    return fithealth.Flag(kind=kind, level=level, text=text,
                          detail=detail or f"detail of {kind}", target=target,
                          params=list(params), stale=stale, covariance_based=cov)


def _bad():
    """physical + degenerate + nocov + noise, built directly (the widget only
    renders; assess() would never pair 'degenerate' with 'no covariance')."""
    h = _ok()
    h.flags = [
        _flag("physical", "bad", "unphysical ×1", "param", [(0, "gl")],
              detail="A: Gauss/Lorentz mix = 1.4 outside [0, 1]"),
        _flag("degenerate", "bad", "degenerate ×1: s0.σCq↔s1.amp (-0.95)",
              "correlations", cov=True),
        _flag("nocov", "check", "no error bars (no covariance)", "errors", cov=True),
        _flag("noise", "check", "residual 3.2× noise", "residual"),
    ]
    return h


def _strip():
    from larmor.desktop.fithealth_strip import FitHealthStrip
    return FitHealthStrip()


def test_strip_empty_and_ok_states(qapp):
    s = _strip()
    s.set_health(None)
    assert "no fit yet" in s.pill.text()
    assert s.chips == [] and s.health() is None
    assert not s.rmsd_chip.isVisibleTo(s)
    h = _ok()
    s.set_health(h)
    assert s.pill.text() == "✓ Fit OK"
    assert s.rmsd_chip.isVisibleTo(s) and s.rmsd_chip.text() == h.chi_text()
    assert s.chips == []                                # passing checks collapse
    assert s.pill.toolTip() == h.tooltip()
    s.close()


def test_chips_mirror_flags_and_emit_their_targets(qapp):
    s = _strip()
    h = _bad()
    s.set_health(h)
    assert [c.text() for c in s.chips] == [f.text for f in h.flags]
    assert [c.toolTip() for c in s.chips] == [f.detail for f in h.flags]
    got = []
    s.focus_param.connect(lambda i, p: got.append(("param", i, p)))
    s.open_correlations.connect(lambda: got.append(("corr",)))
    s.open_errors.connect(lambda: got.append(("errors",)))
    s.show_residual.connect(lambda: got.append(("resid",)))
    s.open_report.connect(lambda: got.append(("report",)))
    for c in s.chips:
        c.click()
    s.rmsd_chip.click()
    assert got == [("param", 0, "gl"), ("corr",), ("errors",), ("resid",), ("report",)]
    s.close()


def test_details_menu_lists_flags_and_reuses_passed_actions(qapp):
    s = _strip()
    h = _bad()
    s.set_health(h)
    extra = [QAction("X"), QAction("Y")]
    m = s.details_menu(extra)
    acts = m.actions()
    texts = [a.text() for a in acts]
    for f in h.flags:
        matches = [a for a in acts if a.text().endswith(f.text)]
        assert len(matches) == 1, f.text
        assert matches[0].isEnabled() == bool(f.target)
    assert any(a is extra[0] for a in acts) and any(a is extra[1] for a in acts)
    assert "Show fit report" in texts and "What the flags mean…" in texts
    got = []
    s.focus_param.connect(lambda i, p: got.append((i, p)))
    next(a for a in acts if a.text().endswith(h.flags[0].text)).trigger()
    assert got == [(0, "gl")]
    # an info flag without a target is listed but cannot be triggered
    h.flags.append(_flag("frozen", "info", "frozen: C", ""))
    s.set_health(h)
    m2 = s.details_menu()
    frozen = next(a for a in m2.actions() if a.text().endswith("frozen: C"))
    assert not frozen.isEnabled()
    s.close()


def test_stale_dims_covariance_chips(qapp):
    s = _strip()
    h = _ok()
    h.stale = True
    h.flags = [_flag("degenerate", "bad", "degenerate ×1: a↔b (+0.97)",
                     "correlations", stale=True, cov=True),
               _flag("noise", "check", "residual 2.1× noise", "residual")]
    s.set_health(h)
    assert s.pill.text().endswith(" · edited since fit")
    assert s.chips[0].toolTip().startswith("from the last fit — F5 to refresh")
    assert not s.chips[1].toolTip().startswith("from the last fit")
    assert s.rmsd_chip.toolTip().startswith("from the last fit")
    s.set_health(_ok())
    assert not s.pill.text().endswith(" · edited since fit")
    # the drag-time hint changes only the pill wording, and only once
    s.set_stale_hint(True)
    assert s.pill.text() == "Model OK · edited since fit"
    assert s.chips == []
    s.set_stale_hint(False)
    assert s.pill.text() == "✓ Fit OK"
    # the unjudged acquisition facts are listed in the details menu, inert
    h = _ok()
    h.unchecked = ["flip angle unknown (I = 3/2) — not judged"]
    s.set_health(h)
    assert s.chips == []                                # never a chip
    (a,) = [a for a in s.details_menu().actions()
            if a.text().startswith("· flip angle unknown")]
    assert not a.isEnabled() and a.toolTip() == h.unchecked[0]
    s.close()


def test_strip_never_forces_window_width(qapp):
    s = _strip()
    h = _ok()
    h.flags = [_flag("noise", "check", "x" * 400, "residual") for _ in range(6)]
    s.set_health(h)
    assert s.sizePolicy().horizontalPolicy() == QSizePolicy.Ignored
    assert s.minimumSizeHint().width() < 200
    assert len(s.chips) == 6 and all(len(c.text()) <= 49 for c in s.chips)
    assert all(c.toolTip().startswith("detail of noise") for c in s.chips)
    s.close()


def test_level_colours_are_readable_in_every_theme(qapp):
    s = _strip()
    s.set_health(_bad())
    prev = theme.active().name
    try:
        for name in theme.names() + theme.aesthetic_names():
            theme.set_active(name)
            s.apply_theme()
            for level in ("ok", "check", "bad"):
                bg, fg = s.colours(level)
                assert theme.contrast(bg, fg) >= 4.5, (name, level, bg, fg)
            bg, fg = s.colours("none")
            assert theme.contrast(bg, fg) >= 3.0, (name, "none", bg, fg)
            assert s.colours("bad")[0] == theme.active().model
    finally:
        theme.set_active(prev)
        s.apply_theme()
    s.close()


def test_widen_relaxation_and_flip_targets_emit_their_signals(qapp):
    from PySide6.QtCore import Qt

    s = _strip()
    h = _ok()
    h.flags = [
        _flag("tail", "check", "tail outside window: A 14 %", "widen", cov=True),
        _flag("recovery", "check", "D1 = 3.0 T1 → 95 % (90° assumed)", "relaxation"),
        _flag("excitation", "info", "flip angle unknown (I = 3/2)", "flip"),
    ]
    s.set_health(h)
    assert [c.text() for c in s.chips] == [f.text for f in h.flags]
    assert [c.toolTip() for c in s.chips] == [f.detail for f in h.flags]
    assert [c.property("flag_kind") for c in s.chips] == ["tail", "recovery", "excitation"]
    assert all(c.cursor().shape() == Qt.PointingHandCursor for c in s.chips)
    got = []
    s.widen_window.connect(lambda: got.append("widen"))
    s.open_relaxation.connect(lambda: got.append("relaxation"))
    s.enter_flip.connect(lambda: got.append("flip"))
    for c in s.chips:
        c.click()
    assert got == ["widen", "relaxation", "flip"]
    acts = s.details_menu().actions()
    for f in h.flags:
        (a,) = [a for a in acts if a.text().endswith(f.text)]
        assert a.isEnabled()
    s.close()
