"""Tier-3 fallback: used only when iNaturalist itself is unreachable.
Real records, fetched and verified once, bundled with the app."""

from pathlib import Path

ASSETS_DIR = Path(__file__).parent / "assets" / "offline"

OFFLINE_CARDS = [
    {
        "catalog_tag": "iNaturalist · Cardinalis cardinalis", "class_tag": "Bird",
        "name": "Northern Cardinal", "subtitle": "Cardinalis cardinalis · Aves",
        "headline": {"dot": "#3E8F63", "text": "Least Concern", "tail": "IUCN Red List (LC)"},
        "grid": [("Order", "Perching Birds", ""), ("Family", "Cardinals and Allies", ""),
                 ("Observations logged", "379,246", "iNaturalist"), ("Conservation status", "Least Concern", "IUCN LC")],
        "foot": "The northern cardinal is a North American bird in the genus Cardinalis, also known as the redbird "
                "or common cardinal. Found from southern Canada through the eastern US to Mexico.",
        "image_file": "cardinal.jpg",
        "image_caption": "(c) Laura Keene, CC BY-NC, via iNaturalist (offline sample)",
    },
    {
        "catalog_tag": "iNaturalist · Haliaeetus leucocephalus", "class_tag": "Bird",
        "name": "Bald Eagle", "subtitle": "Haliaeetus leucocephalus · Aves",
        "headline": {"dot": "#3E8F63", "text": "Least Concern", "tail": "IUCN Red List (LC)"},
        "grid": [("Order", "Hawks, Eagles, Kites, and Allies", ""), ("Family", "Hawks, Eagles, and Kites", ""),
                 ("Observations logged", "229,647", "iNaturalist"), ("Conservation status", "Least Concern", "IUCN LC")],
        "foot": "The bald eagle is a bird of prey found in North America, the national bird of the United States "
                "since 1782. Its scientific name means 'white-headed sea eagle'.",
        "image_file": "eagle.jpg",
        "image_caption": "(c) Addy, CC BY-NC, via iNaturalist (offline sample)",
    },
    {
        "catalog_tag": "iNaturalist · Aptenodytes forsteri", "class_tag": "Bird",
        "name": "Emperor Penguin", "subtitle": "Aptenodytes forsteri · Aves",
        "headline": {"dot": "#8A9E3E", "text": "Near Threatened", "tail": "IUCN Red List (NT)"},
        "grid": [("Order", "Penguins", ""), ("Family", "Penguins", ""),
                 ("Observations logged", "625", "iNaturalist"), ("Conservation status", "Near Threatened", "IUCN NT")],
        "foot": "The tallest and heaviest of all living penguin species, endemic to Antarctica. Despite its iconic "
                "status, it's Near Threatened -- not the safe 'Least Concern' many would assume.",
        "image_file": "penguin.jpg",
        "image_caption": "(c) Martha de Jong-Lantink, CC BY-NC-ND, via iNaturalist (offline sample)",
    },
]


def load_offline_card() -> dict:
    import random
    card = dict(random.choice(OFFLINE_CARDS))
    card["category"] = "bird"
    card["tier"] = "offline"
    card["image_is_cat"] = False
    image_path = ASSETS_DIR / card.pop("image_file")
    card["image_bytes"] = image_path.read_bytes() if image_path.exists() else None
    return card
