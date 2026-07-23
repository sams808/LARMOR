"""Explicit registry of installed fact packs. To add a new pack: build
it under xfact/packs/<name>/ exposing get_random_card(), PACK_NAME, and
PACK_WEIGHT (see packs/birds/ for the smallest complete example), then
add one import + one entry to _PACKS below."""

import random

from ..packs.birds import picker as _birds

# LARMOR vendors only the birds pack of XFact (the birdfact easter egg).
_PACKS = [_birds]


def get_random_card_from_any_pack() -> dict:
    weights = [p.PACK_WEIGHT for p in _PACKS]
    pack = random.choices(_PACKS, weights=weights)[0]
    return pack.get_random_card()


def get_pack(name: str):
    for p in _PACKS:
        if p.PACK_NAME == name:
            return p
    raise KeyError(f"no such pack: {name!r} (installed: {[p.PACK_NAME for p in _PACKS]})")


def list_packs() -> list[str]:
    return [p.PACK_NAME for p in _PACKS]
