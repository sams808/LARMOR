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


def test_references_self_only_same_site_and_param():
    from larmor.constraints_util import references_self
    assert references_self(f"s4.{ISO} + 5.3", 4, ISO) is True      # self
    assert references_self(f"s3.{ISO} + 5.3", 4, ISO) is False     # other line
    # same site, DIFFERENT parameter is legitimate (no cycle)
    assert references_self("s4.amplitude * 0.5", 4, ISO) is False


def test_sanitize_drops_self_reference():
    from larmor.constraints_util import sanitize_constraints
    sites = _sites([None, f"s1.{ISO} + 5.3"])          # s1 links to itself
    dropped = sanitize_constraints(sites)
    assert dropped == [f"s1.{ISO}"]
    assert _iso(sites[1]) is None


def test_sanitize_drops_out_of_range_reference():
    from larmor.constraints_util import sanitize_constraints
    sites = _sites([f"s9.{ISO}", None])                # no s9
    dropped = sanitize_constraints(sites)
    assert dropped == [f"s0.{ISO}"] and _iso(sites[0]) is None


def test_sanitize_breaks_a_cycle():
    from larmor.constraints_util import sanitize_constraints
    # s0 ← s1, s1 ← s0  (a two-line loop)
    sites = _sites([f"s1.{ISO}", f"s0.{ISO}"])
    dropped = sanitize_constraints(sites)
    assert len(dropped) == 1                            # break exactly one link
    remaining = [i for i in (0, 1) if _iso(sites[i]) is not None]
    assert len(remaining) == 1                          # the loop is broken


def test_sanitize_keeps_valid_links():
    from larmor.constraints_util import sanitize_constraints
    sites = _sites([None, f"s0.{ISO} + 5.3"])
    assert sanitize_constraints(sites) == []
    assert _iso(sites[1]) == f"s0.{ISO} + 5.3"


def test_site_refs():
    assert site_refs("0.5 * s0.amplitude + s2.eta") == {0, 2}
    assert site_refs("") == set()


def test_restrict_glass_protocol_applies_eden_ranges():
    """Eden 2023 sec 8.3 one-step restriction: diso confined to +-window
    around the CURRENT value; peak-FWHM floor of 4 ppm ONLY for analytic
    (spin-1/2 style) peak models -- a Czjzek site's shift_fwhm_ppm is dCS
    (a shift-DISTRIBUTION width, legitimately < 4 ppm in published fits)
    and must NOT be floored; linked/pinned params untouched; vary flags
    never changed (restricted, not fixed)."""
    from larmor.constraints_util import restrict_glass_protocol

    sites = [
        {"model": "gauss_lor", "params": {
            "isotropic_chemical_shift_ppm": {"value": -85.0, "vary": True},
            "shift_fwhm_ppm": {"value": 2.0, "vary": True},
            "amplitude": {"value": 100.0, "vary": True, "min": 0.0}}},
        {"model": "czjzek", "params": {
            "isotropic_chemical_shift_ppm": {"value": 60.0, "vary": True},
            "shift_fwhm_ppm": {"value": 2.5, "vary": True},   # dCS!
            "sigma_Cq_MHz": {"value": 1.2, "vary": True},
            "amplitude": {"value": 50.0, "vary": True, "min": 0.0}}},
        {"model": "gauss_lor", "params": {
            # linked shift + pinned width: both must be left alone
            "isotropic_chemical_shift_ppm": {"value": -80.0, "vary": True,
                                             "expr": "s0.isotropic_chemical_shift_ppm + 5"},
            "shift_fwhm_ppm": {"value": 3.0, "vary": False}}},
    ]
    notes = restrict_glass_protocol(sites, shift_window_ppm=3.0)

    p0 = sites[0]["params"]
    assert p0["isotropic_chemical_shift_ppm"]["min"] == -88.0
    assert p0["isotropic_chemical_shift_ppm"]["max"] == -82.0
    assert p0["shift_fwhm_ppm"]["min"] == 4.0
    assert p0["shift_fwhm_ppm"]["value"] == 4.0      # clamped into the floor
    assert "min" not in p0["amplitude"] or p0["amplitude"]["min"] == 0.0

    p1 = sites[1]["params"]
    assert p1["isotropic_chemical_shift_ppm"]["min"] == 57.0
    assert p1["isotropic_chemical_shift_ppm"]["max"] == 63.0
    assert "min" not in p1["shift_fwhm_ppm"]          # dCS NOT floored
    assert p1["shift_fwhm_ppm"]["value"] == 2.5
    assert "min" not in p1["sigma_Cq_MHz"]            # sigma untouched

    p2 = sites[2]["params"]
    assert "min" not in p2["isotropic_chemical_shift_ppm"]   # linked: skipped
    assert "min" not in p2["shift_fwhm_ppm"]                 # pinned: skipped

    # vary flags never changed anywhere (restricted, not fixed)
    assert p0["shift_fwhm_ppm"]["vary"] is True
    assert p2["shift_fwhm_ppm"]["vary"] is False
    # exactly what changed, nothing more: s0 shift + s0 width floor + s1 shift
    # (s1's dCS/sigma untouched; s2 fully skipped)
    assert len(notes) == 3

def test_restrict_glass_protocol_note_count_matches_what_changed():
    from larmor.constraints_util import restrict_glass_protocol
    sites = [{"model": "gauss_lor", "params": {
        "isotropic_chemical_shift_ppm": {"value": 0.0, "vary": True},
        "shift_fwhm_ppm": {"value": 5.0, "vary": True}}}]
    notes = restrict_glass_protocol(sites)
    assert len(notes) == 2                    # one shift window + one floor
    # existing WIDER floor is kept (max(), never loosened)
    sites2 = [{"model": "gauss_lor", "params": {
        "shift_fwhm_ppm": {"value": 8.0, "vary": True, "min": 6.0}}}]
    restrict_glass_protocol(sites2)
    assert sites2[0]["params"]["shift_fwhm_ppm"]["min"] == 6.0
