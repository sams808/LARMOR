"""Shielding-to-shift calibration (larmor.shiftcal): the weighted line with
covariance, its degrees-of-freedom honesty, [calculation] compatibility, the
`larmor shiftcal` CLI and one real-data anchor (CaF2 / NaF)."""
import json

import numpy as np
import pytest

from conftest import (CAF2_MAGRES, F19_STD_CAF2, F19_STD_NAF, NAF_MAGRES,
                      diag9, kogarkoite_like_records, require, synthetic_magres)
from larmor import dft, shiftcal
from larmor.qcpmg_fields import FieldPoint, dcg_at_field, infinite_field_diso, weighted_line
from larmor.shiftcal import CalibrationPoint, ShiftCalibration, fit_calibration


def _calc(**kw) -> dft.MagresCalc:
    base = dict(code="QE-GIPAW", version="7.5", xc="PBE", cutoff_wfc_Ry=60.0,
                cutoff_rho_Ry=480.0,
                pspots={"Ca": "Ca.pbe-spn-kjpaw_psl.1.0.0.UPF",
                        "F": "F.pbe-n-kjpaw_psl.1.0.0.UPF"},
                kgrid=(4, 4, 4), prefix="x")
    base.update(kw)
    return dft.MagresCalc(**base)


def test_weighted_line_matches_infinite_field_diso_bit_for_bit():
    """The closed form was factored out of infinite_field_diso: same
    arithmetic order, so the QCPMG two-field numbers are unchanged."""
    I, diso, cq, eta = 1.5, -50.0, 3.0, 0.7
    sets = [
        [FieldPoint(n, dcg_at_field(diso, cq, n, I, eta)) for n in (58.726, 81.599)],
        [FieldPoint(58.726, -125.9, 5.0), FieldPoint(81.599, -89.31, 5.0)],
        [FieldPoint(58.726, -125.9, 5.0), FieldPoint(81.599, -89.31, 4.0),
         FieldPoint(104.3, -70.0, 6.0)],
    ]
    for pts in sets:
        res = infinite_field_diso(pts, spin=I, eta=eta)
        x = np.array([1.0 / p.larmor_MHz ** 2 for p in pts])
        y = np.array([p.dcg_ppm for p in pts])
        err = np.array([p.dcg_err_ppm or 1.0 for p in pts])
        slope, intercept, cov = weighted_line(x, y, err)
        assert slope == res.slope and intercept == res.intercept
        assert float(np.sqrt(cov[1, 1])) == res.delta_iso_err_ppm
        assert cov[0, 1] == cov[1, 0]
    with pytest.raises(ValueError, match="degenerate"):
        weighted_line([1.0, 1.0], [0.0, 1.0], [1.0, 1.0])


def test_fit_calibration_recovers_a_known_line_with_covariance():
    rng = np.random.default_rng(7)
    x = np.array([232.92, 402.07, 300.0, 350.0])
    err = np.full(4, 0.3)
    y = -0.70 * x + 56.0 + rng.normal(0.0, 0.3, 4)
    pts = [CalibrationPoint(f"p{i}", xi, yi, ei, calc=_calc())
           for i, (xi, yi, ei) in enumerate(zip(x, y, err))]
    cal = fit_calibration(pts, "19F")
    sa, sb = cal.slope_err, cal.intercept_err
    assert abs(cal.slope + 0.70) < 3 * sa
    assert abs(cal.intercept - 56.0) < 3 * sb
    cov = np.asarray(cal.cov)
    assert cov.shape == (2, 2) and cov[0, 1] == cov[1, 0]
    assert cov[0, 0] > 0 and cov[1, 1] > 0
    # equals numpy's unscaled weighted covariance (inflated only by chi2_red > 1)
    _, ref_cov = np.polyfit(x, y, 1, w=1.0 / err, cov="unscaled")
    scale = max(1.0, cal.chi2_red)
    assert cov == pytest.approx(ref_cov * scale, rel=1e-9)
    assert cal.n == 4 and cal.dof == 2 and len(cal.residuals) == 4
    assert cal.rmse_ppm > 0 and cal.flags == [] and cal.kind == "fit"
    assert cal.calc == _calc()
    # delta_err is smallest at the weighted centroid and grows outward
    xc = float(np.mean(x))
    e_c = cal.delta_err(xc)
    assert cal.delta_err(xc - 150) > e_c and cal.delta_err(xc + 150) > e_c
    assert "δ = " in cal.describe() and "n = 4" in cal.describe()


def test_two_points_are_exact_but_flagged_and_bad_inputs_are_refused():
    naf = CalibrationPoint("NaF", 402.0654, -225.0, 0.3, calc=_calc())
    caf2 = CalibrationPoint("CaF2", 232.9226, -108.82, 0.3, calc=_calc())
    cal = fit_calibration([naf, caf2], "19F")
    assert cal.slope == pytest.approx(-0.6869, abs=1e-3)
    assert cal.intercept == pytest.approx(51.17, abs=1e-2)
    assert cal.residuals == pytest.approx([0.0, 0.0], abs=1e-9)
    assert cal.dof == 0 and cal.n == 2
    assert shiftcal.FLAG_ZERO_DOF in cal.flags
    assert "zero residual degrees of freedom" in cal.flags[0]
    assert cal.cov is not None and np.all(np.isfinite(cal.cov))
    assert cal.delta_err(300.0) is not None and cal.delta_err(300.0) > 0
    assert cal.delta(402.0654) == pytest.approx(-225.0)
    assert cal.delta(232.9226) == pytest.approx(-108.82)
    assert "exact by construction" in cal.methods_clause()
    # the same without errors: unweighted, no covariance at all
    bare = fit_calibration([CalibrationPoint("NaF", 402.0654, -225.0),
                            CalibrationPoint("CaF2", 232.9226, -108.82)], "19F")
    assert shiftcal.FLAG_UNWEIGHTED in bare.flags
    assert bare.cov is None and bare.delta_err(300.0) is None
    assert bare.slope == pytest.approx(cal.slope)
    with pytest.raises(ValueError, match="at least two"):
        fit_calibration([naf], "19F")
    with pytest.raises(ValueError, match="at least two"):
        fit_calibration([naf, CalibrationPoint("x", 300.0, None)], "19F")
    with pytest.raises(ValueError, match="identical"):
        fit_calibration([naf, CalibrationPoint("y", 402.0654, -200.0, 0.3)], "19F")


def test_weights_mixed_errors_sigma_ref_and_roundtrip(tmp_path):
    three = [CalibrationPoint("a", 230.0, -108.0, 0.3, calc=_calc()),
             CalibrationPoint("b", 300.0, -156.0, 0.3, calc=_calc()),
             CalibrationPoint("c", 400.0, -226.0, 0.3, calc=_calc())]
    c3 = fit_calibration(three, "19F")
    # a fourth point with a 100x larger error barely moves the line
    loose = CalibrationPoint("d", 350.0, -170.0, 30.0, calc=_calc())
    c4 = fit_calibration(three + [loose], "19F")
    assert c4.slope == pytest.approx(c3.slope, abs=2e-3)
    assert c4.intercept == pytest.approx(c3.intercept, abs=0.6)
    assert c4.dof == 2 and c4.flags == []
    # one point without an error makes the whole fit unweighted
    c5 = fit_calibration(three + [CalibrationPoint("e", 350.0, -170.0)], "19F")
    assert shiftcal.FLAG_UNWEIGHTED in c5.flags and c5.cov is not None
    # sigma_ref
    sr = shiftcal.from_sigma_ref(560.0, "19F")
    assert sr.delta(500.0) == 60.0 and sr.slope == -1.0
    assert sr.kind == "sigma_ref" and sr.zeta_scale == 1.0
    assert sr.delta_err(500.0) is None
    assert "σ_ref = 560.0" in sr.describe() and "σ_ref" in sr.methods_clause()
    # round trip through JSON and through save/load
    d = json.loads(json.dumps(c4.to_dict()))
    assert d["larmor_shiftcal_version"] == 1
    back = ShiftCalibration.from_dict(d)
    assert back.slope == c4.slope and back.intercept == c4.intercept
    assert back.cov == c4.cov and back.flags == c4.flags
    assert [p.name for p in back.points] == ["a", "b", "c", "d"]
    assert back.points[3].delta_err == 30.0 and back.calc == c4.calc
    p = tmp_path / "f19.shiftcal.json"
    c4.save(p)
    loaded = ShiftCalibration.load(p)
    assert loaded.file == str(p) and loaded.describe() == c4.describe()
    assert loaded.zeta_scale == pytest.approx(-c4.slope)
    bad = tmp_path / "x.json"
    bad.write_text("{}")
    with pytest.raises(ValueError, match="not a LARMOR shift calibration"):
        ShiftCalibration.load(bad)


def test_calculation_compatibility_ignores_kgrid_and_other_elements():
    caf2 = _calc(prefix="CaF2")
    naf = _calc(prefix="NaF", pspots={"Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
                                       "F": "F.pbe-n-kjpaw_psl.1.0.0.UPF"})
    assert shiftcal.compatible(caf2, naf, "F") == []
    assert shiftcal.compatible(caf2, naf, "19F") == []
    assert shiftcal.compatible(caf2, _calc(kgrid=(2, 2, 2)), "F") == []
    assert shiftcal.compatible(caf2, _calc(xc="LDA"), "F") == ["calc_xcfunctional"]
    assert shiftcal.compatible(
        caf2, _calc(pspots={"F": "F.pbe-n-rrkjus_psl.1.0.0.UPF"}), "F") == \
        ["calc_pspot[F]"]
    assert shiftcal.compatible(caf2, _calc(cutoff_wfc_Ry=80.0), "F") == \
        ["calc_cutoffenergy"]
    castep = dft.MagresCalc(code="CASTEP")
    diff = shiftcal.compatible(caf2, castep, "F")
    assert "calc_cutoffenergy" in diff and "calc_cutoffenergy_rho" in diff
    assert "calc_code" in diff
    # nothing to compare against
    assert shiftcal.compatible(None, caf2, "F") == []
    assert shiftcal.compatible(caf2, None, "F") == []
    # mismatched references refuse to fit unless overridden (then flagged)
    pts = [CalibrationPoint("a", 230.0, -108.0, 0.3, calc=caf2),
           CalibrationPoint("b", 400.0, -226.0, 0.3, calc=_calc(cutoff_wfc_Ry=80.0))]
    with pytest.raises(ValueError, match="calc_cutoffenergy"):
        fit_calibration(pts, "19F")
    cal = fit_calibration(pts, "19F", allow_mismatch=True)
    assert any("override recorded" in f and "b: calc_cutoffenergy" in f
               for f in cal.flags)
    assert shiftcal.mismatches(pts, "F") == ["b: calc_cutoffenergy"]


def test_cli_shiftcal_builds_a_line_then_magres_uses_it(tmp_path):
    from larmor import cli

    def _ref(name, sigma, cutoff=60.0):
        p = tmp_path / f"{name}.nmr.magres"
        p.write_text(synthetic_magres(
            [("F", 1, diag9(sigma, sigma, sigma)),
             ("F", 2, diag9(sigma + 1e-4, sigma, sigma))],
            prefix=name, cutoff_Ry=cutoff))
        return str(p)

    naf, caf2 = _ref("NaF", 402.0654), _ref("CaF2", 232.9226)
    third = _ref("X", 300.0)
    out = tmp_path / "f19.shiftcal.json"
    rc = cli.main(["shiftcal", "--ref", "NaF", naf, "-225.0", "0.2",
                   "--ref", "CaF2", caf2, "-108.82", "0.2",
                   "--ref", "X", third, "-150.0", "0.2",
                   "--isotope", "19F", "-o", str(out)])
    assert rc == 0 and out.exists()
    cal = ShiftCalibration.load(out)
    assert cal.n == 3 and cal.dof == 1 and cal.isotope == "19F"
    assert [p.n_atoms for p in cal.points] == [2, 2, 2]
    # two references: the CLI says the line is exact by construction
    out2 = tmp_path / "two.shiftcal.json"
    rc = cli.main(["shiftcal", "--ref", "NaF", naf, "-225.0", "0.2",
                   "--ref", "CaF2", caf2, "-108.82", "0.2", "-o", str(out2)])
    assert rc == 0                      # --isotope inferred: only 19F present
    assert shiftcal.FLAG_ZERO_DOF in ShiftCalibration.load(out2).flags
    # a reference from another cutoff is refused unless overridden
    other = _ref("Y", 350.0, cutoff=80.0)
    assert cli.main(["shiftcal", "--ref", "NaF", naf, "-225.0",
                     "--ref", "Y", other, "-190.0", "--isotope", "19F"]) == 2
    assert cli.main(["shiftcal", "--ref", "NaF", naf, "-225.0",
                     "--ref", "Y", other, "-190.0", "--isotope", "19F",
                     "--calc-override"]) == 0
    assert cli.main(["shiftcal", "--ref", "NaF", naf, "-225.0"]) == 2
    assert cli.main(["shiftcal", "--ref", "NaF", naf,
                     "--ref", "CaF2", caf2]) == 2          # missing DELTA

    # ... and `magres` consumes the saved line
    kog = tmp_path / "kog.nmr.magres"
    kog.write_text(synthetic_magres(
        kogarkoite_like_records(), prefix="kog",
        pspots=("Na.pbe-spn-kjpaw_psl.1.0.0.UPF", "F.pbe-n-kjpaw_psl.1.0.0.UPF")))
    rec_p = tmp_path / "kog.recipe.json"
    assert cli.main(["magres", str(kog), "--isotope", "19F", "--calibration",
                     str(out), "-o", str(rec_p)]) == 0
    rec = json.loads(rec_p.read_text(encoding="utf-8"))
    sites = dft.read_magres(kog)
    dft.assign_isotopes(sites)
    groups = dft.group_equivalent(dft.sites_for_isotope(sites, "19F"))
    for s, g in zip(rec["sites"], groups):
        assert s["params"]["isotropic_chemical_shift_ppm"]["value"] == \
            pytest.approx(cal.delta(g.shielding()["iso_ppm"]))
    prov = rec["provenance"]["dft_import"]
    assert prov["calibration"]["n"] == 3
    assert prov["sites"][0]["delta_err_ppm"] > 0


def test_real_caf2_naf_magres_and_2026_03_standards_give_a_negative_slope():
    """The two fluoride standards actually computed (QE-GIPAW 7.5, PBE,
    60/480 Ry) against their 19F MAS spectra measured in the 2026-03 session:
    one crystallographic F site each, a compatible header, and a two-point
    line with a clearly negative slope and the zero-dof flag."""
    for p in (CAF2_MAGRES, NAF_MAGRES, F19_STD_CAF2, F19_STD_NAF):
        require(p)
    from larmor import fit, loader
    from larmor.recipe import Param, Recipe, SiteModel

    caf2 = dft.read_magres_file(CAF2_MAGRES)
    naf = dft.read_magres_file(NAF_MAGRES)
    dft.assign_isotopes(caf2.sites)
    dft.assign_isotopes(naf.sites)
    gc = dft.group_equivalent(caf2.sites)
    assert [(g.isotope, g.multiplicity) for g in gc] == [("43Ca", 4), ("19F", 8)]
    assert gc[1].shielding()["iso_ppm"] == pytest.approx(232.9226, abs=1e-4)
    gn = dft.group_equivalent(dft.sites_for_isotope(naf.sites, "19F"))
    assert len(gn) == 1 and gn[0].multiplicity == 4
    assert gn[0].shielding()["iso_ppm"] == pytest.approx(402.066, abs=2e-3)
    assert shiftcal.compatible(caf2.calc, naf.calc, "F") == []
    assert caf2.calc.describe() == "QE-GIPAW 7.5 · PBE · 60/480 Ry"

    measured = {}
    for name, expno, lit in (("CaF2", F19_STD_CAF2, -108.82),
                             ("NaF", F19_STD_NAF, -225.0)):
        ppm, amp, rec, _meta, _w = loader.load_any(str(expno))
        assert rec["nucleus"] == "19F"
        assert rec["larmor_frequency_MHz"] == pytest.approx(564.27, abs=0.01)
        ppm, amp = np.asarray(ppm, float), np.asarray(amp, float)
        peak = float(ppm[int(np.argmax(amp))])
        recipe = Recipe(nucleus="19F", larmor_frequency_MHz=rec["larmor_frequency_MHz"],
                        spin_rate_Hz=35714.0, sites=[SiteModel(
                            model="gl_norm", label=name, params={
                                "isotropic_chemical_shift_ppm": Param(peak),
                                "shift_fwhm_ppm": Param(2.0, min=0.1),
                                "amplitude": Param(float(amp.max() * 2.0), min=0.0),
                                "gl": Param(0.5, min=0.0, max=1.0)})])
        res = fit.fit(recipe, ppm, amp, window_ppm=(peak + 15.0, peak - 15.0),
                      compute_errorbars=True)
        site = res.recipe.sites[0]
        d = site.params["isotropic_chemical_shift_ppm"]
        # the 2026-03 session reads both standards ~1.7 ppm above the
        # literature values (CaF2 -107.14, NaF -223.33): a referencing
        # offset the measured-reference line absorbs into b, which is the
        # point of calibrating against MEASURED shifts. 2.5 ppm still
        # discriminates the right peak (the two lines are 116 ppm apart).
        assert abs(d.value - lit) < 2.5, (name, d.value)
        assert d.stderr is not None and np.isfinite(d.stderr)
        measured[name] = (d.value, max(float(d.stderr), 0.05))
    pts = [shiftcal.point_from_magres(naf, "19F", "NaF", *measured["NaF"]),
           shiftcal.point_from_magres(caf2, "19F", "CaF2", *measured["CaF2"])]
    cal = fit_calibration(pts, "19F")
    assert -1.1 < cal.slope < -0.5, cal.describe()
    assert cal.dof == 0 and shiftcal.FLAG_ZERO_DOF in cal.flags
    assert cal.delta_err(300.0) is not None
