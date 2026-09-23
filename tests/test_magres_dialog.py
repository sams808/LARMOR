"""Offscreen tests of Tools > Import DFT tensors (.magres): the two-tab
MagresDialog and MainWindow._magres_add_sites. Every file input goes through
its programmatic entry, so no QFileDialog ever opens."""
import json
import os

import numpy as np
import pytest

from conftest import (caf2_like_records, diag9, kogarkoite_like_records,
                      synthetic_magres)
from larmor import dft
from larmor.shiftcal import CalibrationPoint, ShiftCalibration

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture()
def caf2_file(tmp_path):
    p = tmp_path / "CaF2.nmr.magres"
    p.write_text(synthetic_magres(caf2_like_records(), prefix="CaF2",
                                  efg=[("F", 3, diag9(-0.1, -0.1, 0.2))]))
    return str(p)


@pytest.fixture()
def kog_file(tmp_path):
    p = tmp_path / "kog.nmr.magres"
    p.write_text(synthetic_magres(
        kogarkoite_like_records(), prefix="kog",
        pspots=("Na.pbe-spn-kjpaw_psl.1.0.0.UPF", "F.pbe-n-kjpaw_psl.1.0.0.UPF")))
    return str(p)


def _ref_file(tmp_path, name, sigma, **kw):
    p = tmp_path / f"{name}.nmr.magres"
    p.write_text(synthetic_magres([("F", 1, diag9(sigma, sigma, sigma)),
                                   ("F", 2, diag9(sigma, sigma, sigma))],
                                  prefix=name, **kw))
    return str(p)


def _dialog(qapp, recipe=None, nucleus="19F", spin_rate=35714.0):
    from larmor.desktop.magres_dialog import MagresDialog

    if recipe is None:
        recipe = {"nucleus": nucleus, "larmor_frequency_MHz": 564.27,
                  "spin_rate_Hz": spin_rate, "sites": []}
    return MagresDialog(None, nucleus=nucleus, recipe=recipe,
                        spin_rate_Hz=spin_rate, exp_max=100.0)


def _col(d, col):
    return [d.table.item(r, col).text() for r in range(d.table.rowCount())]


def _three_point_line(calc):
    pts = [CalibrationPoint("NaF", 402.0654, -225.0, 0.2, calc=calc),
           CalibrationPoint("CaF2", 232.9226, -108.82, 0.2, calc=calc),
           CalibrationPoint("X", 300.0, -155.0, 0.2, calc=calc)]
    from larmor.shiftcal import fit_calibration
    return fit_calibration(pts, "19F")


def test_dialog_filters_models_by_spin_groups_atoms_and_states_the_sign(qapp, caf2_file):
    d = _dialog(qapp)
    try:
        d._load(caf2_file)
        assert d.lbl.text() == "CaF2.nmr.magres"
        assert d.calc_lbl.text() == "QE-GIPAW 7.5 · PBE · 60/480 Ry"
        assert d.iso.currentText().startswith("19F") and "I = 1/2" in d.iso.currentText()
        models = [d.model.itemData(i) for i in range(d.model.count())]
        assert models == dft.seedable_models(0.5)
        assert d.current_model() == "gl_norm"
        assert "exact multiplicity locks" in d.model.currentText()
        # grouped: one F row of eight, one odd F; C_Q shown as em dash
        assert _col(d, 1) == ["×8", "1"] and _col(d, 0) == ["F3", "F11"]
        assert _col(d, 5) == ["—", "—"]
        assert "spin-1/2" in d.table.horizontalHeaderItem(5).toolTip()
        assert "Herzfeld" in d.table.horizontalHeaderItem(4).toolTip()
        assert "δaniso = −ζ" in d.note.text() and "δ = a·σ + b" in d.note.text()
        assert d.btnAdd.text() == "Add 2 sites to the fit"
        # the lock is on (gl_norm is an area model), the width is shareable
        assert d.chk_lock.isChecked() and "amplitude links" in d.chk_lock.text()
        assert not d.chk_width.isHidden() and d.chk_width.isChecked()
        assert "shift_fwhm_ppm" in d.chk_width.text()
        # switch to the Ca isotope: quadrupolar list, quad_ct recommended
        ca = next(i for i in range(d.iso.count()) if d.iso.itemData(i) == "43Ca")
        d.iso.setCurrentIndex(ca)
        assert "I = 7/2" in d.iso.currentText()
        models = [d.model.itemData(i) for i in range(d.model.count())]
        assert models == dft.seedable_models(3.5) and d.current_model() == "quad_ct"
        assert _col(d, 1) == ["×2"]
        assert not d.chk_lock.isChecked() and "heights" in d.chk_lock.text()
        # back to F, ungrouped: eight + one rows
        d.iso.setCurrentIndex(next(i for i in range(d.iso.count())
                                   if d.iso.itemData(i) == "19F"))
        d.chk_group.setChecked(False)
        assert d.table.rowCount() == 9 and d.btnAdd.text() == "Add 9 sites to the fit"
        # a csa_mas choice defaults the lock off and says so
        d.chk_group.setChecked(True)
        d.model.setCurrentIndex([d.model.itemData(i) for i in range(d.model.count())]
                                .index("csa_mas"))
        assert not d.chk_lock.isChecked() and "heights" in d.chk_lock.text()
        assert "ζ" in d.model.currentText()
    finally:
        d.close()


def test_add_disabled_until_conversion_spectrum_and_matching_settings(qapp, caf2_file):
    d = _dialog(qapp)
    try:
        assert not d.btnAdd.isEnabled() and "open a .magres" in d.warn.text()
        d._load(caf2_file)
        assert not d.btnAdd.isEnabled()
        assert "not chemical shifts" in d.warn.text()
        assert _col(d, 3) == ["—", "—"]                      # no δ_pred yet
        # sigma_ref typed -> enabled, δ_pred = σ_ref − σ
        d.ref.setValue(560.0)
        assert d.rb_ref.isChecked() and d.btnAdd.isEnabled()
        assert _col(d, 3)[0] == f"{560.0 - 232.9226:.1f}"
        assert d.warn.text() == ""
        # a calibration from another functional blocks, names the key
        cal = ShiftCalibration("19F", -0.7, 56.5, cov=[[1e-6, 0], [0, 0.01]], n=3,
                               dof=1, calc=dft.MagresCalc(
                                   code="QE-GIPAW", version="7.5", xc="LDA",
                                   cutoff_wfc_Ry=60.0, cutoff_rho_Ry=480.0,
                                   pspots={"F": "F.pbe-n-kjpaw_psl.1.0.0.UPF"}))
        d._use_line(cal)
        assert d.rb_cal.isChecked()
        assert not d.btnAdd.isEnabled() and "calc_xcfunctional" in d.warn.text()
        assert not d.chk_override.isHidden()
        d.chk_override.setChecked(True)
        assert d.btnAdd.isEnabled() and "recorded" in d.warn.text()
        assert _col(d, 3)[0].startswith(f"{-0.7 * 232.9226 + 56.5:.1f} ±")
        d._accept()
        assert d.result is not None
        assert d.result.provenance["calc_override"] == ["calc_xcfunctional"]
        assert d.result.provenance["calibration"]["slope"] == -0.7
    finally:
        d.close()
    # no spectrum: the dialog opens, Add stays disabled with the reason
    from larmor.desktop.magres_dialog import MagresDialog
    d2 = MagresDialog(None, nucleus="", recipe=None)
    try:
        d2._load(caf2_file)
        d2.ref.setValue(560.0)
        assert not d2.btnAdd.isEnabled()
        assert "load a spectrum" in d2.btnAdd.toolTip()
        assert d2.iso.currentData() == "19F"           # first isotope, no nucleus
    finally:
        d2.close()


def test_calibration_tab_fits_flags_dof_and_hands_the_line_to_sites(qapp, caf2_file, tmp_path):
    d = _dialog(qapp)
    try:
        d._load(caf2_file)
        calc = dft.read_calculation(caf2_file)
        c = d.cal
        c.add_point(CalibrationPoint("NaF", 402.0654, -225.0, 0.3, calc=calc))
        c.add_point(CalibrationPoint("CaF2", 232.9226, -108.82, 0.3, calc=calc))
        assert c.table.rowCount() == 2
        cal = c.fit()
        assert cal is not None and cal.dof == 0
        assert "zero residual degrees of freedom" in c.res.text()
        assert "a = -0.68" in c.res.text()
        assert abs(float(c.table.item(0, 6).text())) < 1e-6
        # a third point: real residuals, dof 1
        c.add_point(CalibrationPoint("X", 300.0, -156.0, 0.3, calc=calc))
        assert c.calibration is None                      # invalidated
        cal = c.fit()
        assert cal.dof == 1 and "±" in c.res.text()
        assert "zero residual" not in c.res.text()
        # hand the line to the Sites tab
        c._use()
        assert d.rb_cal.isChecked() and d.calibration is cal
        assert d.tabs.currentWidget() is d.sites_tab
        sig = 232.9226
        assert _col(d, 3)[0] == f"{cal.delta(sig):.1f} ± {cal.delta_err(sig):.1f}"
        assert d.btnAdd.isEnabled()
        # save / load round-trips through the tab
        p = tmp_path / "f19.shiftcal.json"
        c.save(str(p))
        assert p.exists()
        c.points.clear(); c._rebuild_table()
        back = c.load(str(p))
        assert back.slope == pytest.approx(cal.slope) and c.table.rowCount() == 3
        assert c.btnUse.isEnabled()
        # a typed δ without ± makes the fit unweighted and says so
        c.table.item(2, 4).setText("")
        assert c.points[2].delta_err is None
        cal2 = c.fit()
        assert cal2 is not None and "unweighted" in c.res.text()
        # typing δ through the cell
        c.table.item(2, 3).setText("-157.5")
        assert c.points[2].delta_exp == -157.5 and c.points[2].source == "typed"
        # mismatched references refuse to fit until overridden
        c.add_point(CalibrationPoint("Y", 350.0, -190.0, 0.3, calc=dft.MagresCalc(
            code="QE-GIPAW", version="7.5", xc="PBE", cutoff_wfc_Ry=80.0,
            cutoff_rho_Ry=480.0, pspots={"F": "F.pbe-n-kjpaw_psl.1.0.0.UPF"})))
        assert c.fit() is None and "calc_cutoffenergy" in c.mismatch_lbl.text()
        assert not c.chk_mismatch.isHidden()
        c.chk_mismatch.setChecked(True)
        cal3 = c.fit()
        assert cal3 is not None and any("override" in f for f in cal3.flags)
        c.remove_selected()                                # nothing selected: no-op
        c.table.selectRow(3)
        c.remove_selected()
        assert c.table.rowCount() == 3
    finally:
        d.close()


def test_calibration_tab_reads_references_and_shifts_from_files(qapp, tmp_path):
    from larmor.recipe import Param, Recipe, SiteModel

    d = _dialog(qapp)
    try:
        c = d.cal
        naf = _ref_file(tmp_path, "NaF", 402.0654,
                        pspots=("Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
                                "F.pbe-n-kjpaw_psl.1.0.0.UPF"))
        new = c.add_reference_magres(naf)
        assert [p.name for p in new] == ["NaF"] and new[0].n_atoms == 2
        assert new[0].sigma_iso == pytest.approx(402.0654)
        assert c.table.rowCount() == 1 and c.table.item(0, 3).text() == ""
        # δ from a saved fit: the largest-amplitude site is proposed
        rec = Recipe(nucleus="19F", larmor_frequency_MHz=564.27, sr_hz=718.67,
                     sites=[SiteModel(model="gl_norm", label="ssb", params={
                         "isotropic_chemical_shift_ppm": Param(-160.0, stderr=0.5),
                         "shift_fwhm_ppm": Param(2.0), "amplitude": Param(3.0),
                         "gl": Param(1.0, vary=False)}),
                            SiteModel(model="gl_norm", label="NaF", params={
                         "isotropic_chemical_shift_ppm": Param(-223.33, stderr=0.02),
                         "shift_fwhm_ppm": Param(1.6), "amplitude": Param(100.0),
                         "gl": Param(1.0, vary=False)})])
        rp = tmp_path / "NaF.recipe.json"
        rec.save(rp)
        c.table.selectRow(0)
        delta, err = c.delta_from_recipe(str(rp))
        assert delta == -223.33 and err == 0.02
        assert c.points[0].delta_exp == -223.33 and c.points[0].delta_err == 0.02
        assert c.points[0].source.startswith("fit: NaF.recipe.json NaF")
        assert "SR 718.7 Hz" in c.points[0].source
        assert not c.site_pick.isHidden() and c.site_pick.count() == 2
        # the combo picks the other site
        c.site_pick.setCurrentIndex(0)
        assert c.points[0].delta_exp == -160.0
        c.site_pick.setCurrentIndex(1)
        # a two-site reference file gives two rows
        two = tmp_path / "cryo.nmr.magres"
        two.write_text(synthetic_magres(
            [("F", 1, diag9(300.0, 300.0, 300.0)), ("F", 2, diag9(300.0, 300.0, 300.0)),
             ("F", 3, diag9(340.0, 340.0, 340.0))], prefix="cryolite"))
        new = c.add_reference_magres(str(two))
        assert [p.name for p in new] == ["cryolite F1", "cryolite F3"]
        assert [p.n_atoms for p in new] == [2, 1]
        # δ from the current fit (the dialog's recipe)
        d.recipe.update(rec.to_dict())
        c.recipe = d.recipe
        c.table.selectRow(1)
        assert c.delta_from_current(d.recipe)[0] == -223.33
        assert c.points[1].source.startswith("fit: current fit NaF")
        assert c.delta_from_current({"sites": []}) is None
    finally:
        d.close()


def test_accept_builds_locked_sites_with_provenance_and_note(qapp, kog_file, tmp_path):
    recipe = {"nucleus": "19F", "larmor_frequency_MHz": 564.27,
              "spin_rate_Hz": 35714.0, "sites": [
                  {"model": "gl_norm", "label": "a", "params": {}},
                  {"model": "gl_norm", "label": "b", "params": {}}]}
    d = _dialog(qapp, recipe=recipe)
    try:
        d._load(kog_file)
        cal = _three_point_line(dft.read_calculation(kog_file))
        p = tmp_path / "f19.shiftcal.json"
        cal.save(str(p))
        assert d.load_calibration(str(p)) is not None
        assert "f19.shiftcal.json" in d.cal_lbl.text()
        assert d.current_model() == "gl_norm" and d.chk_lock.isChecked()
        assert _col(d, 1) == ["×2"] * 6
        assert d.btnAdd.isEnabled() and d.warn.text() == ""
        d._accept()
        r = d.result
        assert r is not None and len(r.sites) == 6 and r.model == "gl_norm"
        assert all(set(s) == {"model", "label", "params"} for s in r.sites)
        kog = dft.read_magres(kog_file)
        dft.assign_isotopes(kog)
        groups = dft.group_equivalent(dft.sites_for_isotope(kog, "19F"))
        for s, g in zip(r.sites, groups):
            assert s["params"]["isotropic_chemical_shift_ppm"]["value"] == \
                pytest.approx(cal.delta(g.shielding()["iso_ppm"]))
        assert r.sites[0]["params"]["amplitude"]["expr"] is None
        assert all(s["params"]["amplitude"]["expr"] == "s2.amplitude"
                   for s in r.sites[1:])                  # equal multiplicities
        assert all(s["params"]["shift_fwhm_ppm"]["expr"] == "s2.shift_fwhm_ppm"
                   for s in r.sites[1:])
        prov = r.provenance
        assert prov["file"].endswith("kog.nmr.magres") and len(prov["sha256"]) == 64
        assert prov["isotope"] == "19F" and prov["model"] == "gl_norm"
        assert prov["calculation"]["code"] == "QE-GIPAW"
        assert [s["multiplicity"] for s in prov["sites"]] == [2] * 6
        assert prov["sites"][0]["members"] == ["F3", "F4"]
        assert prov["sites"][0]["delta_err_ppm"] > 0
        assert prov["calibration"]["n"] == 3 and prov["calibration"]["dof"] == 1
        assert prov["calc_override"] == [] and prov["lock_amplitude"] is True
        assert r.ratios == [2] * 6 and r.n_groups == 6
        assert "kog.nmr.magres" in r.note and "a = " in r.note
        assert "f19.shiftcal.json" in r.conversion
        json.dumps(prov)                                    # JSON-safe
    finally:
        d.close()
    # csa_mas: the lock defaults off, its label says heights, ζ/η_CS seeded
    d = _dialog(qapp, recipe=dict(recipe, sites=[]), spin_rate=0.0)
    try:
        d._load(kog_file)
        assert d.current_model() == "csa_mas"              # static spectrum
        assert not d.chk_lock.isChecked() and "heights" in d.chk_lock.text()
        d.ref.setValue(560.0)
        d._accept()
        s0 = d.result.sites[0]["params"]
        # diag(s-12, s, s+12): |zeta| = 12 and eta = 1 by construction
        assert abs(s0["zeta_ppm"]["value"]) == pytest.approx(12.0, abs=1e-3)
        assert s0["eta"]["value"] == pytest.approx(1.0, abs=1e-6)
        assert d.result.sites[1]["params"]["amplitude"]["expr"] is None
        assert d.result.sites[1]["params"]["shift_fwhm_ppm"]["expr"] == \
            "s0.shift_fwhm_ppm"                           # width still shared
        assert d.result.provenance["calibration"]["kind"] == "sigma_ref"
    finally:
        d.close()


def test_app_applies_a_magres_import_and_the_model_simulates(qapp, kog_file, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QDialog

    from larmor import engine, methods
    from larmor.desktop.app import MainWindow
    from larmor.recipe import Recipe

    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    win = MainWindow()
    try:
        # nothing open: the apply step refuses politely
        from larmor.desktop.magres_dialog import MagresImport
        empty = MagresImport([], "", {"file": "x"}, "gl_norm", 0, [], "", True,
                             True, True)
        win.recipe = None
        win._magres_add_sites(empty)
        assert "load a spectrum" in win.statusBar().currentMessage()

        cal = _three_point_line(dft.read_calculation(kog_file))
        x = np.linspace(-300.0, -50.0, 3001)
        y = np.zeros_like(x)
        kog = dft.read_magres(kog_file)
        dft.assign_isotopes(kog)
        for g in dft.group_equivalent(dft.sites_for_isotope(kog, "19F")):
            y += np.exp(-0.5 * ((x - cal.delta(g.shielding()["iso_ppm"])) / 1.0) ** 2)
        win._display_1d(x, y, "19F", 564.27, 35714.0, "kog", "x")
        d = _dialog(qapp, recipe=win.recipe)
        try:
            d._load(kog_file)
            d._use_line(cal)
            d._accept()
            res = d.result
        finally:
            d.close()
        assert res is not None
        win._magres_add_sites(res)
        assert len(win.recipe["sites"]) == 6
        assert win.recipe["provenance"]["dft_import"]["calibration"]["slope"] == \
            pytest.approx(cal.slope)
        assert win.recipe["notes"][-1] == res.note
        msg = win.statusBar().currentMessage()
        assert "6 site" in msg and "2:2:2:2:2:2" in msg and "shared linewidth" in msg
        # the master amplitude was seeded from the spectrum (area model)
        a0 = win.recipe["sites"][0]["params"]["amplitude"]["value"]
        assert a0 == pytest.approx(1.0 / 6 * 5.0 * 1.064, rel=1e-6)
        _, total, per = engine.simulate(Recipe.from_dict(win.recipe), exp_ppm=x)
        assert np.isfinite(total).all() and total.max() > 0 and len(per) == 6
        sent = methods.methods_sentence(win.recipe)
        assert "δ = a·σ + b" in sent and "3 reference compounds" in sent
        win.undo()
        assert len(win.recipe["sites"]) == 0
        assert "dft_import" not in (win.recipe.get("provenance") or {})
        # the menu slot: a rejected dialog changes nothing
        monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Rejected)
        win.open_magres()
        assert len(win.recipe["sites"]) == 0
    finally:
        win.close()
