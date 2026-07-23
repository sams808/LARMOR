"""Wrapper around iNaturalist's public taxa API -- keyless, and the
photo comes structurally attached to the species record itself
(default_photo / taxon_photos), so there's no separate image search
step and no false-positive risk the way the space pack had to guard
against.

Each pick costs two sequential requests (search, then detail) to the
same host. A single real click is nowhere near iNaturalist's real
per-minute rate limit; heavy back-to-back automated testing (as in
simulate_cat_rate.py) can trip it, which reads as a "miss" here rather
than a cat -- get_random_card()'s retry + offline fallback mean an
actual user never sees that, but don't be surprised if a fast test
loop shows a lower resolve rate than real usage would."""

from ...core import http


def fetch_bird(name: str) -> dict | None:
    """Full taxon record for a species by common name."""
    search = http.get_json("https://api.inaturalist.org/v1/taxa",
                            params={"q": name, "rank": "species", "iconic_taxa": "Aves", "per_page": 1})
    if not search or not search.get("results"):
        return None
    taxon_id = search["results"][0]["id"]

    detail = http.get_json(f"https://api.inaturalist.org/v1/taxa/{taxon_id}")
    if not detail or not detail.get("results"):
        return None
    return detail["results"][0]


def fetch_photo_bytes(rec: dict) -> tuple[bytes | None, str]:
    """The species' own photo + its required attribution string."""
    photo = rec.get("default_photo") or {}
    url = photo.get("medium_url") or photo.get("square_url") or photo.get("original_url")
    if not url:
        return None, ""
    return http.get_bytes(url), photo.get("attribution", "")
