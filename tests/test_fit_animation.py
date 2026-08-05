"""Live fit animation: the fitter streams intermediate model curves, and the
SpectrumView draws them with a fading ghost trail."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _model_and_data():
    from larmor.recipe import Recipe, SiteModel, Param
    from larmor.engine import make_context, simulate_site
    x = np.linspace(-20, 120, 300)
    r = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="gauss_lor", label="p", params={
            "isotropic_chemical_shift_ppm": Param(60, min=0, max=120),
            "shift_fwhm_ppm": Param(8, min=1, max=40),
            "gl": Param(0.5, min=0, max=1, vary=False),
            "amplitude": Param(1e6, min=0)})])
    ctx = make_context(r, exp_ppm=x)
    y = np.sum([simulate_site(s, ctx) for s in r.sites], axis=0)
    return r, x, y


def test_fit_streams_frames_when_animated(qapp):
    from larmor.desktop.app import FitWorker
    r, x, y = _model_and_data()
    fw = FitWorker(r.to_dict(), x, y * 1.4, (120, -20), animate=True)
    frames = []
    fw.frame.connect(lambda xx, yy, it: frames.append((it, np.asarray(yy))))
    done = {}
    fw.done.connect(lambda res, mode: done.update(res=res))
    fw.run()
    assert len(frames) >= 1                         # at least one live frame
    assert frames[0][1].shape == x.shape            # a full model curve
    assert done["res"].rmsd < 1.0                   # and the fit still converged


def test_no_frames_when_not_animated(qapp):
    from larmor.desktop.app import FitWorker
    r, x, y = _model_and_data()
    fw = FitWorker(r.to_dict(), x, y * 1.4, (120, -20), animate=False)
    frames = []
    fw.frame.connect(lambda *a: frames.append(a))
    fw.done.connect(lambda *a: None)
    fw.run()
    assert frames == []


def test_spectrumview_animation_trail(qapp):
    from larmor.desktop.plot import SpectrumView
    _, x, y = _model_and_data()
    v = SpectrumView()
    v.start_fit_animation()
    for scale in (0.5, 0.8, 1.0, 1.05):
        v.set_fit_frame(x, y * scale, iteration=3, rms=0.05)
    drawn = sum(1 for g in v._anim_ghosts
                if g.xData is not None and len(g.xData) > 0)
    assert drawn >= 1                               # a ghost trail is present
    assert len(v._anim_main.xData) == len(x)        # current frame drawn
    v.stop_fit_animation()
    assert len(v._anim_main.xData or []) == 0       # cleared for the final model
