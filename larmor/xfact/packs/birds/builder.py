"""Turns an iNaturalist taxon record into a normalized card dict.
Same honesty rules as every other pack: show what the record actually
contains, mark missing fields "not on file", and separate any outside
context (this module has none to add -- iNaturalist's own data already
carries real conservation status and a real descriptive summary)."""

import re

from .curated import IUCN_STATUS

_TAG_RE = re.compile(r"<[^>]+>")


def _strip_html(s: str) -> str:
    return _TAG_RE.sub("", s or "").strip()


def _na(label):
    return (label, None, "na")


def _row(label, value, unit=""):
    return (label, value, unit)


def _ancestor_name(rec: dict, rank: str) -> str | None:
    for a in rec.get("ancestors") or []:
        if a.get("rank") == rank:
            return a.get("preferred_common_name") or a.get("name")
    return None


def _iucn_code(rec: dict) -> str | None:
    for status in rec.get("conservation_statuses") or []:
        if status.get("authority") == "IUCN Red List" and status.get("status"):
            return status["status"].upper()
    top = rec.get("conservation_status") or {}
    return (top.get("iucn_status_code") or top.get("status") or "").upper() or None


def build_bird_card(rec: dict) -> dict:
    common = rec.get("preferred_common_name") or rec.get("name") or "Unknown bird"
    sci = rec.get("name", "")
    order = _ancestor_name(rec, "order")
    family = _ancestor_name(rec, "family")

    code = _iucn_code(rec)
    obs = rec.get("observations_count")

    if code and code in IUCN_STATUS:
        label, color = IUCN_STATUS[code]
        headline = {"dot": color, "text": label, "tail": f"IUCN Red List ({code})"}
    elif obs:
        headline = {"dot": "#5B6FA6", "text": f"{obs:,} observations", "tail": "logged on iNaturalist"}
    else:
        headline = {"unknown": "no conservation status or observation count on file"}

    grid = [
        _row("Order", order, "") if order else _na("Order"),
        _row("Family", family, "") if family else _na("Family"),
        _row("Observations logged", f"{obs:,}" if obs else None, "iNaturalist" if obs else "na"),
        _row("Conservation status", IUCN_STATUS[code][0], f"IUCN {code}") if code in IUCN_STATUS else _na("Conservation status"),
    ]

    summary = _strip_html(rec.get("wikipedia_summary", ""))
    foot = summary if summary else "No Wikipedia summary on file for this record."

    return {
        "category": "bird",
        "catalog_tag": f"iNaturalist · {sci}",
        "class_tag": "Bird",
        "name": common,
        "subtitle": f"{sci} · Aves",
        "headline": headline,
        "grid": grid,
        "foot": foot,
    }
