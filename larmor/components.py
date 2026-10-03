"""Which line of one spectrum is which line of another (Qt-free).

A series fitted one spectrum at a time rarely keeps one model: a component
vanishes half-way, another appears, the user removes a line here and adds
one there. Every place that pairs lines ACROSS spectra -- the carry between
neighbours (``seriesmode.carry_into``), the sweep's warm start and its
trajectory smoothing (``seqfit``), the Series plot's trajectories -- pairs
them through the two functions here, by the line's **label** (the name in
the lines table) rather than by its position in the list:

* :func:`pair_sites` -- the lines of one recipe against another's. Lines
  pair by label when both carry one that is non-empty and unique within its
  own recipe; what is left pairs by index (the old rule) only when at least
  one of the two labels cannot identify its line, so two lines named
  differently are two components, never the same one in disguise.
* :func:`component_map` -- the components of a whole series: each one with
  the site index it has in every member (None where that member's model has
  no such line), in order of first appearance.

Labels come free: Add line names every new line ``<model>-<n>`` and a copied
model keeps its labels, so a series started from one model is consistently
labelled until the user renames (the label column) -- which is exactly how
two lines are told to be the same component or not.
"""
from __future__ import annotations

__all__ = ["label_of", "model_of", "usable_labels", "pair_sites", "component_map",
           "display_names"]


def label_of(site) -> str:
    """The line's label, stripped ('' when none); a site dict or SiteModel."""
    lab = site.get("label") if isinstance(site, dict) else getattr(site, "label", "")
    return str(lab or "").strip()


def model_of(site) -> str:
    m = site.get("model") if isinstance(site, dict) else getattr(site, "model", "")
    return str(m or "")


def sites_of(recipe) -> list:
    """The site list of a recipe dict, a Recipe, or a bare list of sites."""
    if recipe is None:
        return []
    if isinstance(recipe, dict):
        return list(recipe.get("sites") or [])
    sites = getattr(recipe, "sites", None)
    if sites is not None:
        return list(sites)
    return list(recipe)


def usable_labels(sites) -> list:
    """Per line: does its label identify it -- non-empty and unique within
    this recipe? Two lines both called 'A' identify nothing."""
    labels = [label_of(s) for s in sites]
    counts: dict = {}
    for lab in labels:
        counts[lab] = counts.get(lab, 0) + 1
    return [bool(lab) and counts[lab] == 1 for lab in labels]


def pair_sites(dst_sites, src_sites) -> list:
    """For every line of ``dst_sites``, the index of the ``src_sites`` line it
    is the counterpart of, or None.

    1. Lines pair by LABEL when both labels are usable (non-empty, unique in
       their own recipe): 'Q3' pairs with 'Q3' wherever it sits in the list.
    2. What is left pairs by INDEX -- the rule every recipe without labels
       relied on -- but only when at least one of the two labels is unusable
       and that source index is still free. Two lines with different usable
       labels never pair: they are two components.
    """
    dst = list(dst_sites)
    src = list(src_sites)
    d_lab = [label_of(s) for s in dst]
    s_lab = [label_of(s) for s in src]
    d_ok = usable_labels(dst)
    s_ok = usable_labels(src)
    out: list = [None] * len(dst)
    taken: set = set()
    by_label = {lab: j for j, (lab, ok) in enumerate(zip(s_lab, s_ok)) if ok}
    for i, (lab, ok) in enumerate(zip(d_lab, d_ok)):
        if ok and lab in by_label:
            out[i] = by_label[lab]
            taken.add(by_label[lab])
    for i in range(len(dst)):
        if out[i] is not None or i >= len(src) or i in taken:
            continue
        if not d_ok[i] or not s_ok[i]:
            out[i] = i
            taken.add(i)
    return out


def component_map(recipes) -> list:
    """The components of a series whose members need not share one model.

    One entry per component, in order of first appearance along the series:
    ``{"label", "model", "name" (label or model), "usable", "first" (the
    member it first appears in), "index": [site index in each member, None
    where that member has no such line], "n" (members that have it)}``.

    A member's lines attach to the components by the :func:`pair_sites`
    rules -- by usable label first, then by index against the FIRST member
    (the master) when one of the two labels cannot identify its line -- and
    a line that attaches to nothing opens a new component. With one shared
    model (the batch tool) the map is the site list itself.
    """
    recs = [sites_of(r) for r in recipes]
    n = len(recs)
    comps: list = []
    master_at: dict = {}
    for k, sites in enumerate(recs):
        labels = [label_of(s) for s in sites]
        ok = usable_labels(sites)
        for i, s in enumerate(sites):
            target = None
            if ok[i]:
                target = next((c for c in comps if c["usable"] and c["label"] == labels[i]
                               and c["index"][k] is None), None)
            if target is None:
                c = master_at.get(i)
                if (c is not None and c["index"][k] is None
                        and (not ok[i] or not c["usable"])):
                    target = c
            if target is None:
                target = {"label": labels[i], "model": model_of(s),
                          "name": labels[i] or model_of(s), "usable": ok[i],
                          "first": k, "index": [None] * n, "n": 0}
                comps.append(target)
                if k == 0:
                    master_at[i] = target
            target["index"][k] = i
    for c in comps:
        c["n"] = sum(1 for v in c["index"] if v is not None)
    return comps


def display_names(comps) -> list:
    """One readable name per component: its label (or model) when no other
    component shares it, else ``s<k> <name>`` so two unlabelled lines of the
    same model stay apart."""
    names = [c["name"] for c in comps]
    out = []
    for k, c in enumerate(comps):
        nm = c["name"]
        out.append(nm if names.count(nm) == 1 else f"s{k} {nm}")
    return out
