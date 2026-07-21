"""Letter-based inline cell parsing: links (A+20, A+20kHz, 0.5B) and bounds."""
import pytest

from larmor import cellparse as cp


def test_letters_roundtrip():
    for i, s in [(0, "A"), (1, "B"), (25, "Z"), (26, "AA"), (27, "AB"), (51, "AZ"),
                 (52, "BA")]:
        assert cp.index_to_letter(i) == s
        assert cp.letter_to_index(s) == i


def _parse(text, param_name="isotropic_chemical_shift_ppm", unit="ppm",
           this=2, n=5, larmor=130.323):
    return cp.parse_cell(text, param_name=param_name, param_unit=unit,
                         this_index=this, n_sites=n, larmor_MHz=larmor)


def test_plain_number_sets_value_and_clears_link():
    r = _parse("62.6266")
    assert r.set_value and r.value == pytest.approx(62.6266)
    assert r.set_expr and r.expr is None
    assert r.error is None


def test_equal_link():
    r = _parse("A")
    assert r.set_expr and r.expr == "s0.isotropic_chemical_shift_ppm"


def test_offset_link_ppm():
    r = _parse("A+20")
    assert r.expr == "s0.isotropic_chemical_shift_ppm + 20"


def test_offset_link_khz_converts_to_ppm():
    """The user's exact example: A+20kHz at 130.323 MHz → +153.5 ppm."""
    r = _parse("A+20kHz")
    assert r.error is None
    # 20000 Hz / 130.323 MHz = 153.46 ppm
    assert r.expr.startswith("s0.isotropic_chemical_shift_ppm + ")
    off = float(r.expr.rsplit("+", 1)[1])
    assert off == pytest.approx(20000.0 / 130.323, rel=1e-4)


def test_negative_khz_offset():
    r = _parse("B-1.5kHz")
    assert "s1.isotropic_chemical_shift_ppm - " in r.expr
    off = float(r.expr.rsplit("-", 1)[1])
    assert off == pytest.approx(1500.0 / 130.323, rel=1e-4)


def test_ratio_link_amplitude():
    r = _parse("0.5B", param_name="amplitude", unit="")
    assert r.expr == "0.5*s1.amplitude"
    r2 = _parse("0.29*A", param_name="amplitude", unit="")
    assert r2.expr == "0.29*s0.amplitude"


def test_khz_offset_on_cq_converts_to_mhz():
    r = _parse("A+500kHz", param_name="Cq_MHz", unit="MHz")
    assert r.expr == "s0.Cq_MHz + 0.5"        # 500 kHz = 0.5 MHz


def test_bounds_only():
    r = _parse("[0..100]")
    assert r.set_min and r.min == 0.0
    assert r.set_max and r.max == 100.0
    assert not r.set_value and not r.set_expr


def test_value_with_bounds():
    r = _parse("62.6 [0..100]")
    assert r.set_value and r.value == pytest.approx(62.6)
    assert r.min == 0.0 and r.max == 100.0


def test_link_with_bounds():
    r = _parse("A+20 [50:80]")
    assert r.expr == "s0.isotropic_chemical_shift_ppm + 20"
    assert r.min == 50.0 and r.max == 80.0


def test_half_open_bounds():
    r = _parse("[0..]")
    assert r.min == 0.0 and r.max is None
    r2 = _parse("[..100]")
    assert r2.min is None and r2.max == 100.0


def test_errors():
    assert _parse("A", this=1).error is None              # B links to A: fine
    assert "itself" in _parse("C", this=2).error          # C is index 2 == self
    assert "does not exist" in _parse("Z", n=3).error
    assert "min must be < max" in _parse("[10..5]").error
    assert _parse("garbage!!").error is not None


def test_format_link_roundtrip():
    assert cp.format_link("s0.isotropic_chemical_shift_ppm",
                          "isotropic_chemical_shift_ppm") == "A"
    assert cp.format_link("s0.isotropic_chemical_shift_ppm + 20",
                          "isotropic_chemical_shift_ppm") == "A+20"
    assert cp.format_link("0.5*s1.amplitude", "amplitude") == "0.5B"
    # a link to a DIFFERENT parameter is shown raw (not a simple letter form)
    raw = "s0.amplitude + 3"
    assert cp.format_link(raw, "isotropic_chemical_shift_ppm") == raw


def test_parse_format_are_inverse_for_ppm_offset():
    r = _parse("A+20")
    assert cp.format_link(r.expr, "isotropic_chemical_shift_ppm") == "A+20"
