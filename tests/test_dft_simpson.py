"""DFT (.magres) import and SIMPSON bridge."""
import numpy as np
import pytest

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
