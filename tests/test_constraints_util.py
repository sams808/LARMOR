"""Constraint remapping on line deletion — must never leave a self-reference
(which recurses forever at fit time). The site is popped first (as the caller
does), then remap_exprs_after_delete fixes the survivors' references."""
from larmor.constraints_util import remap_exprs_after_delete, site_refs

ISO = "isotropic_chemical_shift_ppm"


def _sites(exprs):
    """exprs: list of expr-string-or-None for each site's iso-shift param."""
    return [{"params": {ISO: {"value": 0.0, "expr": e}}} for e in exprs]


def _iso(site):
    return site["params"][ISO]["expr"]


def _delete(sites, idx):
    sites = list(sites)
    sites.pop(idx)                                # caller pops first
    dropped = remap_exprs_after_delete(sites, idx)
    return sites, dropped


def test_references_above_deleted_shift_down():
    # s3 links to s2; delete s0 → survivors renumber, s3(→s2) links to s1
    sites, dropped = _delete(_sites([None, None, None, f"s2.{ISO} + 5.3"]), 0)
    assert dropped == []
    assert _iso(sites[2]) == f"s1.{ISO} + 5.3"


def test_reference_to_deleted_site_is_dropped():
    sites, dropped = _delete(_sites([None, None, f"s0.{ISO} + 1"]), 0)
    assert dropped == [f"s1.{ISO}"]
    assert _iso(sites[1]) is None


def test_self_reference_is_dropped_not_kept():
    # the reported bug: after deletion an existing/created self-reference must be
    # removed, never left to recurse
    sites, dropped = _delete(_sites([None, f"s1.{ISO} + 5.3"]), 0)   # s1 → s0 (self)
    assert dropped == [f"s0.{ISO}"]
    assert _iso(sites[0]) is None


def test_references_below_deleted_unchanged():
    sites, dropped = _delete(_sites([None, f"s0.{ISO} + 2", None, None]), 3)
    assert dropped == []
    assert _iso(sites[1]) == f"s0.{ISO} + 2"


def test_move_remaps_references_to_follow_lines():
    from larmor.constraints_util import remap_exprs_after_move
    # site 0 links to site 2; swap sites 0 and 2 (reorder for comfort)
    sites = _sites([f"s2.{ISO}", None, None])
    sites = [sites[2], sites[1], sites[0]]          # new order = [old2, old1, old0]
    remap_exprs_after_move(sites, {0: 2, 1: 1, 2: 0})
    # old site0 is now at index 2; its link still points at old site2 (now index 0)
    assert _iso(sites[2]) == f"s0.{ISO}"
    assert _iso(sites[0]) is None


def test_site_refs():
    assert site_refs("0.5 * s0.amplitude + s2.eta") == {0, 2}
    assert site_refs("") == set()
