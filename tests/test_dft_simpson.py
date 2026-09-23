"""DFT (.magres) import and SIMPSON bridge."""
import json

import numpy as np
import pytest

from conftest import (caf2_like_records, diag9, kogarkoite_like_records,
                      synthetic_magres)
from larmor import dft, simpson
from larmor.recipe import Param, Recipe, SiteModel

MAGRES = """#$magres-abinitio-v1.0
[calculation]
calc_code CASTEP
[/calculation]
[atoms]
units lattice Angstrom
units atom Angstrom
atom Al Al 1 0.0 0.0 0.0
atom O O 1 1.7 0.0 0.0
[/atoms]
[magres]
units ms ppm
ms Al 1  550.0 0.0 0.0  0.0 560.0 0.0  0.0 0.0 600.0
efg Al 1  -0.30 0.0 0.0  0.0 -0.40 0.0  0.0 0.0 0.70
ms O 1  200.0 0.0 0.0  0.0 220.0 0.0  0.0 0.0 300.0
efg O 1  -0.60 0.0 0.0  0.0 -0.80 0.0  0.0 0.0 1.40
[/magres]
"""


@pytest.fixture()
def magres_file(tmp_path):
    p = tmp_path / "test.magres"
    p.write_text(MAGRES)
    return p


def test_read_magres(magres_file):
    sites = dft.read_magres(magres_file)
    assert len(sites) == 2
    al = next(s for s in sites if s.label == "Al1")
    assert al.ms_tensor.shape == (3, 3)
    assert al.efg_tensor is not None
    assert al.ms_tensor[2, 2] == pytest.approx(600.0)


def test_read_magres_rejects_other_files(tmp_path):
    p = tmp_path / "not.magres"
    p.write_text("hello world")
    with pytest.raises(ValueError, match="no ms/efg records"):
        dft.read_magres(p)


def test_assign_isotopes(magres_file):
    sites = dft.read_magres(magres_file)
    warnings = dft.assign_isotopes(sites)
    assert warnings == []
    assert {s.isotope for s in sites} == {"27Al", "17O"}
    assert len(dft.sites_for_isotope(sites, "27Al")) == 1


def test_shielding_haeberlen(magres_file):
    al = dft.read_magres(magres_file)[0]
    sh = al.shielding()
    # isotropic = mean of the three eigenvalues
    assert sh["iso_ppm"] == pytest.approx((550 + 560 + 600) / 3)
    # zeta = largest deviation from isotropic
    assert sh["zeta_ppm"] == pytest.approx(600 - 570.0)
    assert 0.0 <= sh["eta"] <= 1.0


def test_quadrupolar_from_efg(magres_file):
    """Validate the EFG -> Cq conversion against an INDEPENDENT calculation
    from fundamental constants: Cq = e*Q*Vzz/h, with Vzz converted from
    atomic units. Agreement to a few % (tabulated Q values differ slightly
    between sources) proves the import is not silently off by a unit factor.
    """
    sites = dft.read_magres(magres_file)
    dft.assign_isotopes(sites)
    al = next(s for s in sites if s.isotope == "27Al")
    q = al.quadrupolar()

    e, h = 1.602176634e-19, 6.62607015e-34
    au_to_v_per_m2 = 9.7173618e21
    q_27al_barn = 0.1466                     # literature quadrupole moment
    vzz_au = 0.70                            # from the magres fixture
    expected_mhz = (e * (q_27al_barn * 1e-28) * vzz_au * au_to_v_per_m2
                    / h / 1e6)
    assert abs(q["Cq_MHz"]) == pytest.approx(expected_mhz, rel=0.05)

    # eta = |Vxx - Vyy| / Vzz with |Vxx| <= |Vyy| <= |Vzz|
    assert q["eta"] == pytest.approx(abs((-0.30 - -0.40) / 0.70), abs=1e-6)
    assert 0.0 <= q["eta"] <= 1.0


def test_quadrupolar_scales_with_efg(magres_file):
    """Cq is linear in Vzz: the O site has exactly 2x the Al EFG here."""
    sites = dft.read_magres(magres_file)
    dft.assign_isotopes(sites)
    al = next(s for s in sites if s.isotope == "27Al")
    o = next(s for s in sites if s.isotope == "17O")
    ratio_efg = 1.40 / 0.70
    from mrsimulator.spin_system.isotope import Isotope

    expected = ratio_efg * (Isotope(symbol="17O").efg_to_Cq
                            / Isotope(symbol="27Al").efg_to_Cq)
    assert (o.quadrupolar()["Cq_MHz"] / al.quadrupolar()["Cq_MHz"]
            == pytest.approx(expected, rel=1e-6))


def test_to_site_dict_flags_shielding_without_reference(magres_file):
    sites = dft.read_magres(magres_file)
    dft.assign_isotopes(sites)
    al = next(s for s in sites if s.isotope == "27Al")
    d = al.to_site_dict(model="quad_ct")
    assert d["model"] == "quad_ct"
    assert d["params"]["Cq_MHz"]["value"] == pytest.approx(
        al.quadrupolar()["Cq_MHz"])
    # without a reference, the shielding must NOT masquerade as a shift
    assert any("SHIELDING" in n for n in d["notes"])

    d2 = al.to_site_dict(model="quad_ct", reference_ppm=560.0)
    assert d2["params"]["isotropic_chemical_shift_ppm"]["value"] == \
        pytest.approx(560.0 - al.shielding()["iso_ppm"])
    assert not any("SHIELDING" in n for n in d2["notes"])


def test_to_site_dict_is_a_valid_recipe_site(magres_file):
    """The imported site must round-trip through the recipe and simulate."""
    from larmor.engine import make_context, simulate_site

    sites = dft.read_magres(magres_file)
    dft.assign_isotopes(sites)
    al = next(s for s in sites if s.isotope == "27Al")
    site_dict = al.to_site_dict(model="quad_ct", reference_ppm=560.0)
    recipe = Recipe.from_dict({
        "nucleus": "27Al", "larmor_frequency_MHz": 195.5,
        "spin_rate_Hz": 20000.0,
        "sites": [{k: v for k, v in site_dict.items() if k != "notes"}]})
    ctx = make_context(recipe, exp_ppm=np.linspace(-200, 200, 1024))
    y = simulate_site(recipe.sites[0], ctx)
    assert np.isfinite(y).all() and y.max() > 0


# ------------------------------------------------ N10: spin-aware seeding
@pytest.fixture()
def qe_file(tmp_path):
    """CaF2-like QE-GIPAW file: 2 Ca, 8 equivalent F, 1 odd F, one F efg."""
    p = tmp_path / "CaF2.nmr.magres"
    p.write_text(synthetic_magres(
        caf2_like_records(), prefix="CaF2",
        efg=[("F", 3, diag9(-0.1, -0.1, 0.2))]))
    return p


@pytest.fixture()
def kog_file(tmp_path):
    p = tmp_path / "kog.nmr.magres"
    p.write_text(synthetic_magres(
        kogarkoite_like_records(), prefix="kog",
        pspots=("Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "F.pbe-n-kjpaw_psl.1.0.0.UPF")))
    return p


def test_to_site_dict_routes_shielding_and_quadrupolar_eta_by_model(magres_file):
    """csa_mas's ``eta`` is the SHIELDING asymmetry; quad_ct's ``eta`` the
    quadrupolar one. Matching on the shared name once wrote the
    quadrupolar eta into csa_mas and C_Q into nothing."""
    sites = dft.read_magres(magres_file)
    dft.assign_isotopes(sites)
    al = next(s for s in sites if s.isotope == "27Al")
    sh, q = al.shielding(), al.quadrupolar()
    assert sh["eta"] != pytest.approx(q["eta"])            # a real test

    csa = al.to_site_dict("csa_mas", reference_ppm=560.0)["params"]
    assert csa["eta"]["value"] == pytest.approx(sh["eta"])
    assert "Cq_MHz" not in csa
    assert csa["zeta_ppm"]["value"] == pytest.approx(sh["zeta_ppm"])

    qc = al.to_site_dict("quad_csa", reference_ppm=560.0)["params"]
    assert qc["eta_cs"]["value"] == pytest.approx(sh["eta"])
    assert qc["eta_q"]["value"] == pytest.approx(q["eta"])
    assert qc["Cq_MHz"]["value"] == pytest.approx(q["Cq_MHz"])
    assert qc["zeta_ppm"]["value"] == pytest.approx(sh["zeta_ppm"])

    ct = al.to_site_dict("quad_ct", reference_ppm=560.0)["params"]
    assert ct["eta"]["value"] == pytest.approx(q["eta"])   # unchanged
    assert ct["Cq_MHz"]["value"] == pytest.approx(q["Cq_MHz"])

    # position-only model on a quadrupolar nucleus: the EFG is dropped, noted
    d = al.to_site_dict("czjzek", reference_ppm=560.0)
    assert any("EFG tensor ignored" in n for n in d["notes"])
    assert set(d) == {"model", "label", "params", "notes"}


def test_seeded_zeta_sign_matches_the_herzfeld_berger_convention():
    """LARMOR's zeta is the SHIELDING anisotropy (delta_aniso = -zeta). With
    sigma_ref the shift tensor is sigma_ref - sigma: same anisotropy sign
    as the shielding, so zeta_seed = zeta_sigma. With a slope a the shift
    tensor is a*sigma + b, its anisotropy a*zeta_sigma and the shielding-
    convention zeta of THAT tensor is -a*zeta_sigma. Two independent in-repo
    statements of the convention are the oracle."""
    from larmor import convert as C
    from larmor.herzfeld_berger import _shift_tensor_aniso
    from larmor.shiftcal import ShiftCalibration

    s = dft.ComputedSite("F1", "19F", 1, ms_tensor=np.diag([200.0, 210.0, 260.0]))
    p = s.to_site_dict("csa_mas", reference_ppm=560.0)["params"]
    diso, zeta, eta = (p["isotropic_chemical_shift_ppm"]["value"],
                       p["zeta_ppm"]["value"], p["eta"]["value"])
    assert diso == pytest.approx(336.6667, abs=1e-3)
    assert zeta == pytest.approx(36.6667, abs=1e-3)
    assert eta == pytest.approx(0.2727, abs=1e-3)
    assert C.csa_principal_from_haeberlen(diso, zeta, eta) == \
        pytest.approx((360.0, 350.0, 300.0), abs=1e-9)

    cal = ShiftCalibration("19F", slope=-0.7, intercept=56.5)
    p2 = s.to_site_dict("csa_mas", calibration=cal)["params"]
    d2, z2, e2 = (p2["isotropic_chemical_shift_ppm"]["value"],
                  p2["zeta_ppm"]["value"], p2["eta"]["value"])
    assert z2 == pytest.approx(0.7 * 36.6667, abs=1e-3)
    assert C.csa_principal_from_haeberlen(d2, z2, e2) == pytest.approx(
        tuple(sorted((-0.7 * v + 56.5 for v in (200.0, 210.0, 260.0)),
                     reverse=True)), abs=1e-6)
    # the anisotropic shift tensor is a * (sigma_eig - sigma_iso), reordered
    aniso = np.sort(np.diag(_shift_tensor_aniso(z2, e2)))
    expect = np.sort(-0.7 * (np.array([200.0, 210.0, 260.0]) - 223.3333333))
    assert aniso == pytest.approx(expect, abs=1e-6)


def test_spin_half_has_no_quadrupolar_and_refuses_quad_models(qe_file):
    sites = dft.read_magres(qe_file)
    dft.assign_isotopes(sites)
    f3 = next(s for s in sites if s.label == "F3")
    assert f3.efg_tensor is not None
    assert f3.quadrupolar() is None                        # spin-1/2 guard
    d = f3.to_site_dict("gl_norm", reference_ppm=560.0)
    assert "Cq_MHz" not in d["params"]
    assert any("EFG tensor ignored" in n for n in d["notes"])
    with pytest.raises(ValueError, match="spin-1/2"):
        f3.to_site_dict("quad_ct", reference_ppm=560.0)
    with pytest.raises(ValueError, match="cannot be seeded"):
        f3.to_site_dict("jmultiplet", reference_ppm=560.0)
    assert dft.spin_of("19F") == 0.5
    assert dft.spin_of("27Al") == 2.5
    assert dft.spin_of("Xx99") is None
    assert dft.element_of("19F") == "F" and dft.element_of("27Al") == "Al"


def test_seedable_models_and_default_model_follow_spin():
    from larmor import models

    half = dft.seedable_models(0.5)
    assert half[0] == "gl_norm"
    assert "csa_mas" in half and "voigt" in half
    assert not any(models.get(n).needs_quadrupolar for n in half)
    quad = dft.seedable_models(2.5)
    assert quad[0] == "quad_ct" and "quad_csa" in quad and "czjzek" in quad
    assert len(dft.seedable_models(None)) == 13
    assert dft.default_model(0.5, 35714.0) == "gl_norm"
    assert dft.default_model(0.5, 0.0) == "csa_mas"
    assert dft.default_model(0.5, None) == "gl_norm"
    assert dft.default_model(1.5, 20000.0) == "quad_ct"
    assert dft.default_model(None) == "quad_ct"
    assert "position only" in dft.model_choice_label("gl_norm") or \
        "exact multiplicity" in dft.model_choice_label("gl_norm")
    assert "ζ" in dft.model_choice_label("csa_mas")


def test_read_calculation_and_fingerprint(qe_file, tmp_path, magres_file):
    from larmor.recipe import sha256_of

    calc = dft.read_calculation(qe_file)
    assert calc.code == "QE-GIPAW" and calc.version == "7.5"
    assert calc.xc == "PBE"
    assert (calc.cutoff_wfc_Ry, calc.cutoff_rho_Ry) == (60.0, 480.0)
    assert calc.pspots["F"] == "F.pbe-n-kjpaw_psl.1.0.0.UPF"
    assert calc.kgrid == (4, 4, 4) and calc.prefix == "CaF2"
    assert calc.describe() == "QE-GIPAW 7.5 · PBE · 60/480 Ry"
    fp = calc.fingerprint("19F")
    assert set(fp) == {"calc_code", "calc_code_version", "calc_xcfunctional",
                       "calc_cutoffenergy", "calc_cutoffenergy_rho",
                       "calc_pspot[F]"}
    assert "calc_pspot[Ca]" not in fp
    mf = dft.read_magres_file(qe_file)
    assert mf.sha256 == sha256_of(qe_file) and mf.calc.code == "QE-GIPAW"
    assert len(mf.sites) == 11
    # the CASTEP-style fixture has calc_code only
    castep = dft.read_calculation(magres_file)
    assert castep.code == "CASTEP"
    assert castep.cutoff_wfc_Ry is None and castep.cutoff_rho_Ry is None
    assert castep.describe() == "CASTEP"
    # a CASTEP header with eV cutoffs and 'calc_pspot F 2|...' converts
    p = tmp_path / "c.magres"
    p.write_text("[calculation]\ncalc_code CASTEP\ncalc_cutoffenergy 816.34 eV\n"
                 "calc_pspot F 2|1.4|20|24|26|20:21(qc=8)\n[/calculation]\n"
                 "[magres]\nms F 1 1 0 0 0 1 0 0 0 1\n[/magres]\n")
    c2 = dft.read_calculation(p)
    assert c2.cutoff_wfc_Ry == pytest.approx(60.0, abs=0.01)
    assert c2.pspots["F"].startswith("2|1.4")
    # no header at all -> an empty record, not an error
    assert dft.MagresCalc().describe() == "no calculation header"
    rt = dft.MagresCalc.from_dict(json.loads(json.dumps(calc.to_dict())))
    assert rt == calc


def test_group_equivalent_collapses_symmetry_equivalent_atoms(qe_file, kog_file):
    sites = dft.read_magres(qe_file)
    dft.assign_isotopes(sites)
    f = dft.sites_for_isotope(sites, "19F")
    groups = dft.group_equivalent(f)
    assert [g.multiplicity for g in groups] == [8, 1]
    assert groups[0].label == "F3"
    assert groups[0].members == [f"F{i}" for i in range(3, 11)]
    assert any("8 equivalent" in n for n in groups[0].notes)
    assert groups[0].shielding()["iso_ppm"] == pytest.approx(232.9226, abs=1e-4)
    assert groups[1].shielding()["iso_ppm"] == pytest.approx(237.9226, abs=1e-4)
    # different isotopes never merge even at equal shielding
    ca = dft.group_equivalent(sites)
    assert [(g.isotope, g.multiplicity) for g in ca] == \
        [("43Ca", 2), ("19F", 8), ("19F", 1)]

    kog = dft.read_magres(kog_file)
    dft.assign_isotopes(kog)
    kf = dft.group_equivalent(dft.sites_for_isotope(kog, "19F"))
    assert [g.multiplicity for g in kf] == [2] * 6
    assert [g.label for g in kf] == ["F3", "F5", "F7", "F9", "F11", "F13"]
    # the representative keeps the anisotropy (no orientation smearing)
    assert kf[0].shielding()["zeta_ppm"] == pytest.approx(12.0, abs=1e-3)

    # tolerance: 0.2 ppm apart stays separate at 0.05, merges at 1.0
    a = dft.ComputedSite("F1", "19F", 1, ms_tensor=np.diag([100.0, 100.0, 100.0]))
    b = dft.ComputedSite("F2", "19F", 2, ms_tensor=np.diag([100.2, 100.2, 100.2]))
    assert len(dft.group_equivalent([a, b])) == 2
    assert len(dft.group_equivalent([a, b], tol_ppm=1.0)) == 1
    # an Al pair 0.5 MHz apart in C_Q stays separate
    al1 = dft.ComputedSite("Al1", "27Al", 1, ms_tensor=np.eye(3) * 500,
                           efg_tensor=np.diag([-0.3, -0.4, 0.7]))
    al2 = dft.ComputedSite("Al2", "27Al", 2, ms_tensor=np.eye(3) * 500,
                           efg_tensor=np.diag([-0.3, -0.4, 0.7]) * 1.2)
    assert len(dft.group_equivalent([al1, al2])) == 2


def test_sites_to_recipe_dicts_locks_by_multiplicity_and_round_trips():
    from larmor import cellparse, constraints_util, fit
    from larmor.shiftcal import ShiftCalibration

    a = dft.ComputedSite("F1", "19F", 1, ms_tensor=np.eye(3) * 300.0,
                         multiplicity=2, members=["F1", "F2"])
    b = dft.ComputedSite("F3", "19F", 3, ms_tensor=np.eye(3) * 320.0,
                         multiplicity=1, members=["F3"])
    cal = ShiftCalibration("19F", -0.7, 56.5, cov=[[1e-6, 0.0], [0.0, 0.04]])
    dicts, notes = dft.sites_to_recipe_dicts([a, b], "gl_norm", calibration=cal,
                                             first_index=3)
    assert [set(d) for d in dicts] == [{"model", "label", "params"}] * 2
    assert dicts[0]["params"]["amplitude"]["expr"] is None
    assert dicts[1]["params"]["amplitude"]["expr"] == "0.5*s3.amplitude"
    assert dicts[1]["params"]["shift_fwhm_ppm"]["expr"] == "s3.shift_fwhm_ppm"
    assert cellparse.format_link("0.5*s3.amplitude", "amplitude") == "0.5D"
    assert dicts[0]["params"]["isotropic_chemical_shift_ppm"]["value"] == \
        pytest.approx(-0.7 * 300.0 + 56.5)
    assert any("×2 (F1, F2)" in n and "δ_pred" in n and "±" in n for n in notes)
    assert any("amplitude locked to s3" in n for n in notes)

    # equal multiplicities -> the bare link
    b.multiplicity = 2
    d2, _ = dft.sites_to_recipe_dicts([a, b], "gl_norm", calibration=cal,
                                      first_index=3)
    assert d2[1]["params"]["amplitude"]["expr"] == "s3.amplitude"

    # the dicts load behind three placeholder sites and the links evaluate
    ph = {"model": "gauss_lor", "label": "p", "params": {
        "isotropic_chemical_shift_ppm": {"value": 0.0}, "shift_fwhm_ppm": {"value": 1.0},
        "amplitude": {"value": 1.0}, "gl": {"value": 1.0, "vary": False}}}
    rec = Recipe.from_dict({"nucleus": "19F", "larmor_frequency_MHz": 564.27,
                            "spin_rate_Hz": 35714.0,
                            "sites": [ph, ph, ph] + dicts})
    params = fit._make_params(rec)
    assert params["s4_amp"].value == pytest.approx(0.5 * params["s3_amp"].value)
    assert params["s4_fwhm"].value == pytest.approx(params["s3_fwhm"].value)
    sites = [ph, ph, ph] + [json.loads(json.dumps(d)) for d in dicts]
    assert constraints_util.sanitize_constraints(sites) == []

    # a quadrupolar model shares no width (its breadth key is a C_Q)
    al = dft.ComputedSite("Al1", "27Al", 1, ms_tensor=np.eye(3) * 500,
                          efg_tensor=np.diag([-0.3, -0.4, 0.7]), multiplicity=2)
    al2 = dft.ComputedSite("Al2", "27Al", 2, ms_tensor=np.eye(3) * 520,
                           efg_tensor=np.diag([-0.3, -0.4, 0.7]), multiplicity=1)
    dq, _ = dft.sites_to_recipe_dicts([al, al2], "quad_ct", reference_ppm=560.0)
    assert dq[1]["params"]["amplitude"]["expr"] == "0.5*s0.amplitude"
    assert dq[1]["params"]["Cq_MHz"]["expr"] is None
    assert dq[1]["params"]["shift_fwhm_ppm"]["expr"] is None
    # locks off -> no exprs at all
    d0, n0 = dft.sites_to_recipe_dicts([a, b], "gl_norm", calibration=cal,
                                       lock_amplitude=False, share_width=False)
    assert all(p["expr"] is None for d in d0 for p in d["params"].values())
    # no conversion at all: the shielding is flagged, as before
    _, n1 = dft.sites_to_recipe_dicts([a], "gl_norm")
    assert any("SHIELDING" in n for n in n1)


def test_cli_magres_auto_model_grouping_and_calibration(qe_file, kog_file, tmp_path):
    from larmor import cli
    from larmor.shiftcal import ShiftCalibration

    cal = ShiftCalibration("19F", -0.6869, 51.17, cov=[[6e-6, 0.0], [0.0, 0.6]],
                           n=2, dof=0, calc=dft.read_calculation(qe_file))
    calp = tmp_path / "f19.shiftcal.json"
    cal.save(calp)
    out = tmp_path / "kog.recipe.json"
    rc = cli.main(["magres", str(kog_file), "--isotope", "19F", "--calibration",
                   str(calp), "--lock-multiplicity", "-o", str(out)])
    assert rc == 0 and out.exists()
    rec = json.loads(out.read_text(encoding="utf-8"))
    assert [s["model"] for s in rec["sites"]] == ["gl_norm"] * 6
    kog = dft.read_magres(kog_file)
    dft.assign_isotopes(kog)
    groups = dft.group_equivalent(dft.sites_for_isotope(kog, "19F"))
    for s, g in zip(rec["sites"], groups):
        assert s["params"]["isotropic_chemical_shift_ppm"]["value"] == \
            pytest.approx(-0.6869 * g.shielding()["iso_ppm"] + 51.17)
    assert rec["sites"][0]["params"]["amplitude"]["expr"] is None
    assert all(s["params"]["amplitude"]["expr"] == "s0.amplitude"
               for s in rec["sites"][1:])
    assert any("×2" in n for n in rec["notes"])
    prov = rec["provenance"]["dft_import"]
    assert prov["calibration"]["slope"] == pytest.approx(-0.6869)
    assert prov["sites"][0]["multiplicity"] == 2 and prov["calc_override"] == []
    assert prov["sha256"] and prov["model"] == "gl_norm"
    Recipe.load(out)                                    # loads cleanly

    # a quadrupolar model on 19F is refused, nothing written
    out2 = tmp_path / "bad.recipe.json"
    assert cli.main(["magres", str(kog_file), "--isotope", "19F", "--model",
                     "quad_ct", "-o", str(out2)]) == 2
    assert not out2.exists()

    # a calibration from an 80 Ry run does not apply to a 60 Ry file
    p80 = tmp_path / "kog80.nmr.magres"
    p80.write_text(synthetic_magres(kogarkoite_like_records(), prefix="kog",
                                    cutoff_Ry=80.0))
    cal80 = ShiftCalibration("19F", -0.7, 56.0, calc=dft.read_calculation(p80))
    cal80p = tmp_path / "f19_80.shiftcal.json"
    cal80.save(cal80p)
    out3 = tmp_path / "o3.recipe.json"
    assert cli.main(["magres", str(kog_file), "--isotope", "19F",
                     "--calibration", str(cal80p), "-o", str(out3)]) == 2
    assert not out3.exists()
    assert cli.main(["magres", str(kog_file), "--isotope", "19F",
                     "--calibration", str(cal80p), "--calc-override",
                     "-o", str(out3)]) == 0
    prov3 = json.loads(out3.read_text(encoding="utf-8"))["provenance"]["dft_import"]
    assert prov3["calc_override"] == ["calc_cutoffenergy"]

    # --no-group writes one site per atom; the old --reference form works
    out4 = tmp_path / "o4.recipe.json"
    assert cli.main(["magres", str(kog_file), "--isotope", "19F", "--no-group",
                     "--reference", "560", "-o", str(out4)]) == 0
    rec4 = json.loads(out4.read_text(encoding="utf-8"))
    assert len(rec4["sites"]) == 12
    assert rec4["provenance"]["dft_import"]["calibration"]["kind"] == "sigma_ref"
    # listing without an isotope still works (no recipe)
    assert cli.main(["magres", str(qe_file)]) == 0


# ---------------------------------------------------------------- SIMPSON
def _recipe():
    return Recipe(nucleus="13C", larmor_frequency_MHz=100.6, spin_rate_Hz=10000,
                  sites=[SiteModel(model="quad_ct", label="C1", params={
                      "isotropic_chemical_shift_ppm": Param(30.0),
                      "Cq_MHz": Param(2.5), "eta": Param(0.3),
                      "shift_fwhm_ppm": Param(1.0), "amplitude": Param(1.0)})])


def test_spinsys_block_contains_the_physics():
    ss = simpson.spinsys_block(_recipe())
    assert ss.startswith("spinsys {") and ss.rstrip().endswith("}")
    assert "channels 13C" in ss
    assert "quadrupole 1 2 2.5e+06" in ss.replace("2500000", "2.5e+06")
    assert "shift 1 30.0p" in ss


def test_spinsys_block_with_dipolar_partner():
    ss = simpson.spinsys_block(_recipe(), partner={"isotope": "15N",
                                                   "dipolar_hz": -900.0})
    assert "15N" in ss
    assert "dipole 1 2 -900" in ss


def test_redor_input_is_complete_text():
    text = simpson.redor_input(_recipe(), "15N", -900.0, spin_rate_hz=10000.0)
    assert "spinsys {" in text and "par {" in text and "proc main" in text
    assert "spin_rate        10000" in text
    assert "variable tr      0.0001" in text


def test_run_without_simpson_says_so_clearly(monkeypatch):
    monkeypatch.setattr(simpson, "simpson_available", lambda: None)
    with pytest.raises(RuntimeError, match="not on PATH"):
        simpson.run("spinsys {}")


def test_parse_fid(tmp_path):
    p = tmp_path / "sim.fid"
    p.write_text("SIMP\nNP=3\nSW=1000\nTYPE=FID\nDATA\n1.0 0.0\n0.5 0.1\n"
                 "0.25 0.05\nEND\n")
    res = simpson.parse_fid(p)
    assert res.y.size == 3
    assert res.y[0] == pytest.approx(1.0 + 0j)
    assert res.y[1] == pytest.approx(0.5 + 0.1j)
    assert res.x[1] == pytest.approx(1 / 1000)


def test_parse_fid_rejects_junk(tmp_path):
    p = tmp_path / "junk.fid"
    p.write_text("not a simpson file")
    with pytest.raises(ValueError, match="not a SIMPSON"):
        simpson.parse_fid(p)
