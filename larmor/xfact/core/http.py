"""Shared HTTP helpers every pack uses: short-timeout GET wrappers that
return None instead of raising, plus the two truly cross-pack
utilities -- fetching a specific NASA image asset by id, and the
universal "no picture found" cat fallback. Pack-specific API calls
(arcsecond, iNaturalist, ...) live in each pack's own sources.py.
"""

import random
import time

import requests

TIMEOUT = 6  # seconds; short on purpose -- a hung request shouldn't stall a card
             # for long, but generous enough that a slightly slow link still
             # reaches iNaturalist rather than dropping to the offline pool

session = requests.Session()
session.headers.update({"User-Agent": "XFact/1.0 (+github.com/sams808/XFact)"})


def get_json(url, params=None):
    try:
        r = session.get(url, params=params, timeout=TIMEOUT)
        if r.status_code == 429:
            # a single real-world app click is never concurrent enough to
            # trigger this on its own -- this only guards against sharing
            # an IP/network with other traffic hitting the same API
            time.sleep(0.5)
            r = session.get(url, params=params, timeout=TIMEOUT)
        if r.status_code != 200:
            return None
        return r.json()
    except (requests.RequestException, ValueError):
        return None


def get_bytes(url):
    try:
        r = session.get(url, timeout=TIMEOUT)
        if r.status_code != 200 or not r.content:
            return None
        return r.content
    except requests.RequestException:
        return None


def fetch_nasa_image_bytes(nasa_id: str) -> bytes | None:
    """Best available resolution (medium, else small/thumb/orig) for a
    specific, individually-verified NASA image id -- see the space
    pack's curated.py for why this is a direct id, not a live search."""
    manifest = get_json(f"https://images-api.nasa.gov/asset/{nasa_id}")
    if not manifest:
        return None
    hrefs = [it.get("href", "") for it in manifest.get("collection", {}).get("items", [])]
    for size in ("~medium.jpg", "~small.jpg", "~thumb.jpg", "~orig.jpg"):
        for href in hrefs:
            if href.endswith(size):
                data = get_bytes(href.replace("http://", "https://"))
                if data:
                    return data
    return None


def fetch_cat_bytes() -> bytes | None:
    """A cat, from cataas.com -- the universal 'no picture found'
    fallback shared by every pack."""
    return get_bytes(f"https://cataas.com/cat?_={random.randint(0, 10**9)}")
