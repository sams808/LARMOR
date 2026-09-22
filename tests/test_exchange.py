"""Two-site chemical exchange (`exchange2`): the closed-form steady-state
Bloch-McConnell lineshape (Gutowsky & Holm 1956; McConnell 1958).

Physics limits are checked against their analytic forms: k -> 0 (two
population-weighted Lorentzians), k -> inf (one Lorentzian at the weighted
mean), the equal-population coalescence threshold k_ex = sqrt(2) pi dnu, the
fast-exchange residual width 4 pi p_A p_B dnu^2 / k_ex and the k-independence
of the area. Then the registry/partition wiring, seeding, fit recovery, the
k_c read-out in the table and an add/undo smoke test in the app.
"""
import os

import numpy as np
import pytest

from larmor import models
from larmor.models.analytic import two_site_exchange
from larmor.models.base import SimContext
from larmor.recipe import Param, Recipe, SiteModel


def _lorentz(nu, nu0, fwhm):
    return np.real(1.0 / (np.pi * fwhm + 2j * np.pi * (nu - nu0)))


def _n_maxima(y):
    return int(((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])).sum())


def test_slow_limit_is_two_population_weighted_lorentzians():
    nu = np.linspace(-2000.0, 3000.0, 50001)
    y = two_site_exchange(nu, 0.0, 1000.0, 0.3, 0.0, 50.0)
    ref = 0.3 * _lorentz(nu, 0.0, 50.0) + 0.7 * _lorentz(nu, 1000.0, 50.0)
    assert np.allclose(y, ref, rtol=0.0, atol=1e-12)
    # peak heights 3 : 7 (each peak carries ~0.15 % of the other line's tail
    # at 20 half-widths, hence the tolerance)
    ia, ib = int(np.argmin(np.abs(nu))), int(np.argmin(np.abs(nu - 1000.0)))
    assert y[ia] / y[ib] == pytest.approx(0.3 / 0.7, rel=5e-3)


def test_fast_limit_is_one_lorentzian_at_the_weighted_mean():
    nu = np.linspace(-2000.0, 3000.0, 50001)            # 0.1 Hz step
    y = two_site_exchange(nu, 0.0, 1000.0, 0.3, 1e9, 50.0)
    assert np.isfinite(y).all()
    assert abs(nu[int(np.argmax(y))] - 700.0) <= 0.1 + 1e-9
    ref = _lorentz(nu, 700.0, 50.0)
    assert np.allclose(y, ref, rtol=0.0, atol=1e-4 * ref.max())


def test_coalescence_threshold_is_sqrt2_pi_dnu():
    """For k_ex = k_AB + k_BA and equal populations the two maxima merge at
    k_ex = sqrt(2) pi dnu = 4.44 dnu (NOT pi dnu / sqrt(2), which is the
    one-way rate): two maxima persist at 4.3 dnu, one remains at 4.6 dnu."""
    dnu = 1000.0
    nu = np.linspace(-1500.0, 2500.0, 40001)            # 0.1 Hz step
    for f in (3.0, 4.3):
        assert _n_maxima(two_site_exchange(nu, 0.0, dnu, 0.5, f * dnu, 1.0)) == 2, f
    for f in (4.6, 6.0):
        assert _n_maxima(two_site_exchange(nu, 0.0, dnu, 0.5, f * dnu, 1.0)) == 1, f
    assert np.sqrt(2.0) * np.pi == pytest.approx(4.443, abs=1e-3)


def test_fast_exchange_residual_width_and_area_invariance():
    dnu = 1000.0
    k = 100.0 * dnu
    nu = np.linspace(-1000.0, 2000.0, 300001)           # 0.01 Hz step
    y = two_site_exchange(nu, 0.0, dnu, 0.5, k, 2.0)
    above = nu[y > 0.5 * y.max()]
    fwhm = float(above[-1] - above[0])
    assert fwhm == pytest.approx(2.0 + 4.0 * np.pi * 0.25 * dnu ** 2 / k,
                                 rel=0.10)
    # the area does not depend on k: exchange redistributes intensity
    wide = np.linspace(-2e5, 2e5, 800001)               # 0.5 Hz step
    areas = []
    for kk in (0.0, 1e2, 1e4, 1e6, 1e9):
        yy = two_site_exchange(wide, 0.0, dnu, 0.5, kk, 20.0)
        assert yy.min() >= -1e-12                       # an absorption
        areas.append(float(np.trapezoid(yy, wide)))
    assert np.allclose(areas, areas[0], rtol=1e-3)


def test_registry_render_positions_and_partitions():
    from larmor import constraints_util, engine, estimate, fit

    ctx = SimContext("31P", 100.0, 0.0, np.linspace(-10.0, 15.0, 5001))
    v = {"isotropic_chemical_shift_ppm": 1.0, "split_ppm": 5.0, "pop_a": 0.4,
         "k_ex_hz": 1.0, "lorentz_fwhm_ppm": 0.2, "amplitude": 2.0}
    m = models.get("exchange2")
    assert not m.needs_quadrupolar
    y = m.render(dict(v), ctx)
    x = ctx.x_ppm
    peaks = x[1:-1][(y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])]
    assert peaks.size == 2
    assert peaks.max() == pytest.approx(1.0 + 0.6 * 5.0, abs=0.01)  # delta_A
    assert peaks.min() == pytest.approx(1.0 - 0.4 * 5.0, abs=0.01)  # delta_B
    assert y.max() == pytest.approx(2.0)
    assert float((x * y).sum() / y.sum()) == pytest.approx(1.0, abs=0.1)
    # no Larmor frequency: still finite (ppm units stand in for Hz)
    assert np.isfinite(m.render(dict(v), SimContext("31P", 0.0, 0.0, x))).all()
    # the four partitions, the Jacobian step and the grid rule
    assert "exchange2" in fit._ANALYTIC_MODELS
    assert "exchange2" not in fit._SIMULATED_MODELS
    assert "exchange2" in engine._GRID_RESTRICTABLE
    assert "exchange2" in constraints_util._NOT_PEAK_FWHM_MODELS
    assert estimate._WIDTH_KEY["exchange2"] == ("lorentz_fwhm_ppm", False)
    rec = Recipe(nucleus="31P", larmor_frequency_MHz=162.0, sites=[
        SiteModel(model="exchange2", label="e",
                  params={"amplitude": Param(1.0)})])
    assert fit.diff_step_for(rec) is None
    assert engine.grid_restrictable(rec) and not engine.needs_kernel(rec)
    # data-driven seeds span the observed band
    xs = np.linspace(-20.0, 20.0, 4001)
    doublet = (np.exp(-((xs + 2.0) / 0.5) ** 2)
               + np.exp(-((xs - 2.0) / 0.5) ** 2))
    sv = estimate.start_values("exchange2", xs, doublet, "31P", 162.0)
    assert sv["split_ppm"] > 0.0 and sv["lorentz_fwhm_ppm"] > 0.0


def _exchange_recipe(values, larmor=162.0):
    params = models.get("exchange2").defaults()
    for k, val in values.items():
        params[k].value = val
    return Recipe(nucleus="31P", larmor_frequency_MHz=larmor, spin_rate_Hz=0.0,
                  sites=[SiteModel(model="exchange2", label="e", params=params)])


def test_fit_recovers_k_split_and_population():
    from larmor import engine
    from larmor import fit as fitmod

    truth = {"isotropic_chemical_shift_ppm": 0.0, "split_ppm": 4.0,
             "pop_a": 0.6, "k_ex_hz": 800.0, "lorentz_fwhm_ppm": 0.8,
             "amplitude": 1.0}
    x = np.linspace(-15.0, 15.0, 3001)
    _, y, _ = engine.simulate(_exchange_recipe(truth), exp_ppm=x)
    y = y + np.random.default_rng(0).normal(0.0, 0.01, x.size)
    start = _exchange_recipe({**truth, "k_ex_hz": 200.0, "split_ppm": 3.0,
                              "pop_a": 0.5})
    res = fitmod.fit(start, x, y)
    p = res.recipe.sites[0].params
    assert p["k_ex_hz"].value == pytest.approx(800.0, rel=0.15)
    assert p["split_ppm"].value == pytest.approx(4.0, abs=0.2)
    assert p["pop_a"].value == pytest.approx(0.6, abs=0.05)
    assert not res.at_bounds


# ---------------------------------------------------------------- desktop

@pytest.fixture(scope="module")
def qapp():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication

    yield QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")   # never inherit a real session
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


def test_k_ex_cell_shows_coalescence_readout(qapp):
    from larmor.desktop.table import LinesTable

    m = models.get("exchange2")
    params = {p.name: {"value": p.default, "vary": p.vary,
                       "min": p.min, "max": p.max} for p in m.params}
    params["split_ppm"]["value"] = 5.0
    rec = {"nucleus": "31P", "larmor_frequency_MHz": 100.0, "spin_rate_Hz": 0.0,
           "sites": [{"model": "exchange2", "label": "ex", "params": params}]}
    t = LinesTable()
    t.rebuild(rec, set())
    for key in ("split_ppm", "pop_a", "k_ex_hz", "lorentz_fwhm_ppm"):
        assert key in t._used_keys, key
    col = 2 + t._used_keys.index("k_ex_hz")
    cell = t.table.cellWidget(0, col)
    # sqrt(2) pi (5 ppm x 100 MHz = 500 Hz) = 2221 s^-1, shown to 3 sig. figs
    assert "k_c" in cell.derived.text() and "2.22e+03" in cell.derived.text()
    assert "coalescence" in cell.derived.toolTip()
    assert "k_ex" in t.table.horizontalHeaderItem(col).text()


def test_add_exchange_and_correlated_sites_in_app_and_undo(qapp, win):
    x = np.linspace(-20.0, 20.0, 801)
    doublet = np.exp(-((x + 2.0) / 0.5) ** 2) + np.exp(-((x - 2.0) / 0.5) ** 2)
    win._display_1d(x, doublet, "31P", 162.0, None, "A", "sA")
    qapp.processEvents()
    win._model_actions["exchange2"].setChecked(True)
    win.add_site_at(0.0, 1.0)
    win._sim_timer.stop()             # add/undo is under test, not the live sim
    qapp.processEvents()
    sites = win.recipe["sites"]
    assert len(sites) == 1 and sites[0]["model"] == "exchange2"
    assert sites[0]["params"]["isotropic_chemical_shift_ppm"]["value"] == 0.0
    assert sites[0]["params"]["split_ppm"]["value"] > 0.0
    assert "k_ex_hz" in win.lines_table._used_keys
    assert "split_ppm" in win.lines_table._used_keys
    assert "exchange2" in win.statusBar().currentMessage()
    win.undo()
    win._sim_timer.stop()
    qapp.processEvents()
    assert len(win.recipe["sites"]) == 0
    win.redo()
    win._sim_timer.stop()
    qapp.processEvents()
    assert len(win.recipe["sites"]) == 1
    win._model_actions["exchange2"].setChecked(False)

    # a quadrupolar trace: the correlated Czjzek starts at slope 0 with a
    # data-seeded sigma and dCS (no remembered defaults for a new model)
    xa = np.linspace(-100.0, 150.0, 1001)
    ya = np.exp(-((xa - 60.0) / 15.0) ** 2)
    win._display_1d(xa, ya, "27Al", 195.483, 20000.0, "B", "sB")
    qapp.processEvents()
    win._model_actions["czjzek_corr"].setChecked(True)
    win.add_site_at(60.0, 1.0)
    win._sim_timer.stop()             # no background kernel build in a smoke test
    qapp.processEvents()
    s = win.recipe["sites"][-1]
    assert s["model"] == "czjzek_corr"
    assert s["params"]["shift_slope_ppm_per_MHz"]["value"] == 0.0
    assert s["params"]["sigma_Cq_MHz"]["value"] > 0.0
    assert s["params"]["shift_fwhm_ppm"]["value"] > 0.0
    assert "shift_slope_ppm_per_MHz" in win.lines_table._used_keys
    win._model_actions["czjzek_corr"].setChecked(False)
    w = getattr(win, "_sim_worker", None)
    if w is not None and w.isRunning():
        w.wait()
