"""Pairing lines across spectra by label (larmor.components): the rule the
series carry, the sweep's warm start / smoothing and the Series plot share."""
from larmor.components import component_map, display_names, pair_sites, usable_labels


def _s(label="", model="gauss_lor"):
    return {"model": model, "label": label, "params": {}}


def test_pair_sites_by_label_first_then_by_index_only_where_a_label_cannot_identify():
    # labels pair wherever the line sits; a name the other side lacks pairs with nothing
    assert pair_sites([_s("A"), _s("B"), _s("C")], [_s("C"), _s("A")]) == [1, None, 0]
    # no labels at all: the old index rule
    assert pair_sites([_s(), _s(), _s()], [_s(), _s()]) == [0, 1, None]
    # two different usable names are two components, never one in disguise
    assert pair_sites([_s("A")], [_s("B")]) == [None]
    # a repeated label identifies nothing: those lines fall back to index,
    # the uniquely named one still pairs by name (and its index is then taken)
    assert usable_labels([_s("A"), _s("A"), _s("X")]) == [False, False, True]
    assert pair_sites([_s("A"), _s("A"), _s("X")],
                      [_s("X"), _s("Y"), _s("Z")]) == [None, 1, 0]
    # one side unlabelled, the other named: index, the way a recipe saved
    # before labels existed meets one saved after
    assert pair_sites([_s(), _s()], [_s("A"), _s("B")]) == [0, 1]
    # whitespace is not a label
    assert pair_sites([_s("  ")], [_s("A")]) == [0]


def test_pair_sites_accepts_site_models():
    from larmor.recipe import Param, SiteModel

    a = SiteModel(model="gauss_lor", label="A", params={"amplitude": Param(1.0)})
    b = SiteModel(model="gauss_lor", label="B", params={"amplitude": Param(1.0)})
    assert pair_sites([a, b], [b, a]) == [1, 0]


def test_component_map_follows_a_series_whose_members_differ():
    r0 = {"sites": [_s("A"), _s("B")]}
    r1 = {"sites": [_s("B"), _s("A"), _s("C")]}          # reordered, C appears
    r2 = {"sites": [_s("A"), _s("C")]}                   # B gone
    comps = component_map([r0, r1, r2])
    assert [c["name"] for c in comps] == ["A", "B", "C"]
    assert [c["index"] for c in comps] == [[0, 1, 0], [1, 0, None], [None, 2, 1]]
    assert [c["n"] for c in comps] == [3, 2, 2] and [c["first"] for c in comps] == [0, 0, 1]
    assert display_names(comps) == ["A", "B", "C"]
    # one shared model (the batch tool): the map is the site list itself
    same = component_map([r0, r0, r0])
    assert [c["index"] for c in same] == [[0, 0, 0], [1, 1, 1]]
    # unlabelled members pair by index against the master; an extra line
    # beyond the master's count opens a component of its own each time
    u0 = {"sites": [_s(), _s()]}
    u1 = {"sites": [_s(), _s(), _s()]}
    comps = component_map([u0, u1])
    assert [c["index"] for c in comps] == [[0, 0], [1, 1], [None, 2]]
    assert display_names(comps) == ["s0 gauss_lor", "s1 gauss_lor", "s2 gauss_lor"]
    # a usable name meeting an unlabelled master line at the same index: by index
    comps = component_map([{"sites": [_s(), _s("Q")]}, {"sites": [_s("P"), _s("Q")]}])
    assert [c["index"] for c in comps] == [[0, 0], [1, 1]]
    assert component_map([]) == [] and component_map([{"sites": []}]) == []
