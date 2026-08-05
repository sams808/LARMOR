"""Parameter diff between two fits (compare vs a reference)."""
import pytest

from larmor.recipe_diff import recipe_diff


def _rec(pos, width):
    return {"sites": [{"model": "gauss_lor", "label": "A", "params": {
        "isotropic_chemical_shift_ppm": {"value": pos},
        "shift_fwhm_ppm": {"value": width},
        "amplitude": {"value": 100.0}, "gl": {"value": 1.0}}}]}


def test_diff_computes_delta_and_percent():
    rows = recipe_diff(_rec(15.5, 6.0), _rec(15.0, 6.0))
    d = next(r for r in rows if r["param"] == "isotropic_chemical_shift_ppm")
    assert d["current"] == 15.5 and d["reference"] == 15.0
    assert d["delta"] == 0.5
    assert d["delta_pct"] == pytest.approx(100 * 0.5 / 15.0)
    # gl is excluded from the diff
    assert not any(r["param"] == "gl" for r in rows)


def test_diff_handles_missing_site():
    cur = {"sites": _rec(15, 6)["sites"] + _rec(2, 3)["sites"]}
    ref = _rec(15, 6)                              # only one site
    rows = recipe_diff(cur, ref)
    # site 1's params have no reference → reference/delta are None
    site1 = [r for r in rows if r["site"] == 1]
    assert site1 and all(r["reference"] is None for r in site1)
    assert all(r["delta"] is None for r in site1)
