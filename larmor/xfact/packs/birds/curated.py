"""Curated bird species list.

Unlike the space pack, this list didn't need a CONFIRMED/CAT_ONLY split:
iNaturalist ties a photo directly to the taxon record it returns (no
separate search step, so no false-positive risk), and a real batch
test across all 40 of these came back 40/40 with a usable photo. The
cat is kept only as a genuine safety net for the rare taxon record
that has no photo at all, not as an engineered ratio.

Mix of common backyard birds, a few icons, and some genuinely unusual
ones for variety.
"""

SPECIES = [
    "Northern Cardinal", "Bald Eagle", "Blue Jay", "American Robin", "Peregrine Falcon",
    "Great Horned Owl", "Ruby-throated Hummingbird", "Emperor Penguin", "Wandering Albatross",
    "Common Ostrich", "Resplendent Quetzal", "Atlantic Puffin", "Greater Flamingo", "Snowy Owl",
    "Bar-headed Goose", "Arctic Tern", "Kakapo", "Shoebill", "Secretarybird", "Andean Condor",
    "Superb Lyrebird", "Rainbow Lorikeet", "European Robin", "House Sparrow", "American Crow",
    "Mallard", "Canada Goose", "Great Blue Heron", "Barn Owl", "Red-tailed Hawk",
    "Golden Eagle", "California Condor", "Wild Turkey", "Mute Swan", "Common Kingfisher",
    "Toco Toucan", "Scarlet Macaw", "Southern Cassowary", "Greater Roadrunner", "Anna's Hummingbird",
]

# IUCN Red List code -> (plain label, severity color). Order matters
# for nothing here, just readability.
IUCN_STATUS = {
    "LC": ("Least Concern", "#3E8F63"),
    "NT": ("Near Threatened", "#8A9E3E"),
    "VU": ("Vulnerable", "#C97F2E"),
    "EN": ("Endangered", "#B8452F"),
    "CR": ("Critically Endangered", "#7A2734"),
    "EW": ("Extinct in the Wild", "#4A2E4A"),
    "EX": ("Extinct", "#171A21"),
    "DD": ("Data Deficient", "#7B818C"),
    "NE": ("Not Evaluated", "#7B818C"),
}
