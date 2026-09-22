"""The general-d Czjzek model (`czjzek_d`): Czjzek's dimensionality as a
parameter on the shared (C_Q, eta) kernel.

Le Caer & Brand's Gaussian Isotropic Model is the d = 5 Czjzek distribution
mrsimulator implements (LARMOR's `czjzek`); its extended GIM is `ext_czjzek`.
So the deliverable is the general-d family (dmfit CzSimple's <d>), which must
(a) reproduce the fit engine's weights exactly at d = 5, (b) carry the
closed-form chi-law invariants, and (c) be wired into every hand-maintained
table, the dmfit round trip and the P(C_Q) dialog.
"""
import os

import numpy as np
import pytest

from larmor import czjzek_dist as cd
from larmor import models
from larmor.models.base import SimContext
from larmor.recipe import Param, Recipe, SiteModel

# the context of tests/test_models_v2.py, so ONE cached kernel serves every
# slow test of the Czjzek family in a run
CTX_AL = SimContext("27Al", 195.483, 20000.0, np.linspace(-150, 200, 2048))

BASE = {"isotropic_chemical_shift_ppm": 60.0, "sigma_Cq_MHz": 1.5,
        "shift_fwhm_ppm": 8.0, "line_fwhm_ppm": 0.0, "amplitude": 1.0}


def _cog(x, y):
    return float((x * y).sum() / y.sum())


def test_rms_and_mode_follow_sqrt_d():
    """P_Q is chi-distributed with d dof and scale sigma_Cz = 2 sigma, so
    sqrt<P_Q^2> = 2 sqrt(d) sigma for every d; a lower d piles weight toward
    C_Q -> 0, so the marginal mode is ordered d = 2 < 3 < 5."""
    eta = np.linspace(0.0, 1.0, 201)
    for sigma in (1.0, 2.5):
        cq = np.linspace(1e-3, 12.0 * sigma, 3000)
        CQ, ETA = np.meshgrid(cq, eta, indexing="xy")
        pq2 = (CQ ** 2 * (1.0 + ETA ** 2 / 3.0)).ravel()
        modes = []
        for d in (2, 3, 5):
            w = cd.czjzek_weights(sigma, d, cq, eta)
            assert w.shape == (eta.size * cq.size,)
            assert w.sum() == pytest.approx(1.0)
            assert np.sqrt((w * pq2).sum()) == pytest.approx(
                2.0 * np.sqrt(d) * sigma, rel=0.01)
            marg = w.reshape(eta.size, cq.size).sum(axis=0)
            modes.append(cq[int(np.argmax(marg))])
        assert modes[0] < modes[1] < modes[2]
    # d = 1 (a flat C_Q^0 prefactor) stays finite on a grid starting above 0
    w1 = cd.czjzek_weights(1.0, 1.0, np.linspace(1e-3, 12.0, 3000), eta)
    assert np.isfinite(w1).all() and w1.sum() == pytest.approx(1.0)


@pytest.mark.slow
def test_render_d5_equals_czjzek():
    """The identity that makes the model trustworthy: d = 5 (explicit or by
    default) renders as the flagship czjzek on the same kernel."""
    y_cz = models.get("czjzek").render(dict(BASE), CTX_AL)
    y_d5 = models.get("czjzek_d").render({**BASE, "czjzek_d": 5.0}, CTX_AL)
    y_def = models.get("czjzek_d").render(dict(BASE), CTX_AL)
    assert np.isfinite(y_cz).all() and y_cz.max() > 0.9
    assert np.allclose(y_d5, y_cz, rtol=1e-5, atol=1e-6)
    assert np.allclose(y_def, y_cz, rtol=1e-5, atol=1e-6)


@pytest.mark.slow
def test_lower_d_is_narrower_and_less_shifted():
    """d = 2 puts its weight at lower C_Q than d = 5: a narrower pattern with
    a smaller second-order shift (centre of gravity at higher ppm)."""
    from larmor import estimate

    x = CTX_AL.x_ppm
    y5 = models.get("czjzek_d").render({**BASE, "czjzek_d": 5.0}, CTX_AL)
    y2 = models.get("czjzek_d").render({**BASE, "czjzek_d": 2.0}, CTX_AL)
    for y in (y2, y5):
        assert np.isfinite(y).all()
        assert y.max() == pytest.approx(1.0, abs=0.01)
    assert estimate.band_width_ppm(x, y2, frac=0.10)[1] < \
        estimate.band_width_ppm(x, y5, frac=0.10)[1]
    assert _cog(x, y2) > _cog(x, y5)


def test_registry_tables_and_columns():
    from larmor import constraints_util, engine, estimate, fit
    from larmor.desktop import table

    names = {m["name"] for m in models.describe_all()}
    assert "czjzek_d" in names
    m = models.get("czjzek_d")
    assert m.needs_quadrupolar
    pd = {p.name: p for p in m.params}
    assert (pd["czjzek_d"].default, pd["czjzek_d"].min, pd["czjzek_d"].max,
            pd["czjzek_d"].vary) == (5.0, 1.0, 5.0, False)
    assert pd["sigma_Cq_MHz"].max == 40.0
    assert "czjzek_d" in fit._SIMULATED_MODELS
    assert "czjzek_d" not in fit._ANALYTIC_MODELS
    assert "czjzek_d" in engine._GRID_RESTRICTABLE
    assert "czjzek_d" in constraints_util._NOT_PEAK_FWHM_MODELS
    assert estimate._WIDTH_KEY["czjzek_d"] == ("sigma_Cq_MHz", True)
    rec = Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="czjzek_d", label="d",
                  params={"amplitude": Param(1.0)})])
    assert engine.needs_kernel(rec)
    assert "czjzek_d" in [k for k, _ in table.PARAM_COLUMNS]
    # name-collision guard: the `function` model owns the parameter literally
    # called `d`, and nothing else may -- columns, labels and multi-field
    # sharing are keyed by parameter name across models
    owners = sorted(n for n in models.REGISTRY
                    if "d" in models.get(n).param_names)
    assert owners == ["function"]


def _czsimple_fxmla(d_tag: str) -> str:
    """A minimal dmfit .fxmla with one CzSimple line (the flat layout of
    tests/test_fxmla.py's newformat file, without the spectrum block)."""
    return (
        "<dmfit>\n<FitParameters>\n"
        "<DMFitVersion>dmfit #20230120</DMFitVersion>\n"
        "<FitModeAsc>Fit 1D</FitModeAsc>\n"
        "<Dimension>F2<nucleus>27Al</nucleus><frequency>156.28</frequency>"
        "<spinrate>35714</spinrate>\n"
        "<line><ModelName>CzSimple</ModelName><Group>"
        '<amp Fix="*">129.38</amp><pos Fix="*" Unit="ppm">64.58</pos>'
        '<dCS Fix="*" Unit="ppm">9.55</dCS><lb>-3.86</lb>'
        f'<CQ>4756.57</CQ><sCZ_CQ Fix="*">2378.29</sCZ_CQ>{d_tag}'
        "</Group></line>\n"
        "</Dimension>\n</FitParameters>\n</dmfit>\n")


def test_fxml_czsimple_d_round_trip(tmp_path):
    from larmor.io import export, fxmla

    def imported(tag, d_tag):
        p = tmp_path / f"cz_{tag}.fxmla"
        p.write_text(_czsimple_fxmla(d_tag), encoding="utf-8")
        recipe, warnings = fxmla.to_recipe(fxmla.read(p))
        return recipe.sites[0], warnings

    s3, warn3 = imported("d3", "<d>3</d>")
    assert s3.model == "czjzek_d"
    assert s3.params["czjzek_d"].value == 3.0
    assert s3.params["czjzek_d"].vary is False
    assert s3.params["sigma_Cq_MHz"].value == pytest.approx(2378.29 * 0.5 / 1000.0)
    assert any("czjzek_d" in w for w in warn3)
    # d = 5, or no <d> at all, stays the plain czjzek (existing round trips
    # and the sigma-convention test are untouched)
    assert imported("d5", "<d>5</d>")[0].model == "czjzek"
    assert imported("none", "")[0].model == "czjzek"

    # export a czjzek_d site: CzSimple carrying its <d>, re-imported as czjzek_d
    site = SiteModel(model="czjzek_d", label="AlIV", params={
        "isotropic_chemical_shift_ppm": Param(60.0),
        "sigma_Cq_MHz": Param(1.6),
        "czjzek_d": Param(3.2, vary=False, min=1.0, max=5.0),
        "shift_fwhm_ppm": Param(8.0), "line_fwhm_ppm": Param(0.0),
        "amplitude": Param(1.0)})
    r1 = Recipe(nucleus="27Al", larmor_frequency_MHz=156.28,
                spin_rate_Hz=20000.0, sites=[site])
    x = np.linspace(-50.0, 150.0, 401)
    y = np.exp(-((x - 60.0) / 15.0) ** 2)
    out = tmp_path / "czd_export.fxmla"
    text = export.export_fxmla(r1, x, y, out)
    assert "<ModelName>CzSimple</ModelName>" in text
    assert "<d>3.2</d>" in text
    r2, _ = fxmla.to_recipe(fxmla.read(out))
    assert r2.sites[0].model == "czjzek_d"
    assert r2.sites[0].params["czjzek_d"].value == pytest.approx(3.2)
    assert r2.sites[0].params["sigma_Cq_MHz"].value == pytest.approx(1.6)


def test_dialog_and_batch_are_d_aware():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication, QLabel
    QApplication.instance() or QApplication([])
    from larmor import batch
    from larmor.desktop.czjzek_dist_dialog import CzjzekDistDialog

    site = {"model": "czjzek_d", "label": "AlIV", "params": {
        "sigma_Cq_MHz": {"value": 1.8}, "czjzek_d": {"value": 3.0}}}
    dlg = CzjzekDistDialog(None, {"nucleus": "27Al", "sites": [site]})
    texts = " ".join(lb.text() for lb in dlg.findChildren(QLabel))
    assert "AlIV" in texts and "d = 3" in texts
    dlg.close()
    cols = dict((h, (v, e)) for h, v, e in batch._site_columns(site, {}))
    assert cols["√⟨P_Q²⟩ (MHz)"][0] == pytest.approx(2.0 * np.sqrt(3.0) * 1.8)
    assert cols["d (Czjzek)"][0] == 3.0
