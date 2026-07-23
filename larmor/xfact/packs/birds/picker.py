"""Orchestrates one bird card: pick a species, fetch its taxon record
live, use its own attached photo (or a cat, on the rare record with
none), and fall back to the bundled offline sample if iNaturalist
itself is unreachable.

Simpler than the space pack's picker: no CONFIRMED/CAT_ONLY split
needed since a real 40-species batch test came back 100% photo hits --
see xfact/tests/simulate_cat_rate.py.
"""

import random

from ...core import http
from . import builder, offline, sources
from .curated import SPECIES

PACK_NAME = "birds"
PACK_WEIGHT = 1.0


def get_live_card() -> dict | None:
    name = random.choice(SPECIES)
    rec = sources.fetch_bird(name)
    if not rec:
        return None
    card = builder.build_bird_card(rec)

    image_bytes, attribution = sources.fetch_photo_bytes(rec)
    is_cat = False
    if image_bytes:
        caption = f"{attribution} · via iNaturalist" if attribution else "via iNaturalist"
    else:
        image_bytes = http.fetch_cat_bytes()
        is_cat = bool(image_bytes)
        caption = "cataas.com — no picture on file for this one, so: cat"

    if not image_bytes:
        return None

    card["image_bytes"] = image_bytes
    card["image_is_cat"] = is_cat
    card["image_caption"] = caption
    card["tier"] = "live"
    return card


def get_random_card() -> dict:
    return get_live_card() or get_live_card() or offline.load_offline_card()
