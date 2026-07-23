"""Tier-3 fallback: used only when iNaturalist itself is unreachable.

Real iNaturalist taxon records + their own photos, fetched and verified once and
bundled with the app, so the offline experience still has genuine variety instead
of just two or three birds. Each card carries its photographer's CC attribution.
"""

from pathlib import Path

ASSETS_DIR = Path(__file__).parent / "assets" / "offline"

OFFLINE_CARDS = [{'catalog_tag': 'iNaturalist · Cardinalis cardinalis',
  'class_tag': 'Bird',
  'name': 'Northern Cardinal',
  'subtitle': 'Cardinalis cardinalis · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Perching Birds', ''],
           ['Family', 'Cardinals and Allies', ''],
           ['Observations logged', '379,267', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The northern cardinal (Cardinalis cardinalis) is a North American bird in the genus Cardinalis; it is '
          'also known colloquially as the redbird or common cardinal. It can be found in southern Canada, through '
          'the eastern United States from Maine to Texas and south through Mexico. Its habitat includes woodlands, '
          'gardens, shrublands, wetlands.',
  'image_file': 'northern_cardinal.jpg',
  'image_caption': '(c) Laura Keene, some rights reserved (CC BY-NC), uploaded by Laura Keene · via iNaturalist '
                   '(offline sample)'},
 {'catalog_tag': 'iNaturalist · Haliaeetus leucocephalus',
  'class_tag': 'Bird',
  'name': 'Bald Eagle',
  'subtitle': 'Haliaeetus leucocephalus · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Hawks, Eagles, Kites, and Allies', ''],
           ['Family', 'Hawks, Eagles, and Kites', ''],
           ['Observations logged', '229,647', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The bald eagle (Haliaeetus leucocephalus, from Greek hali "sea", aiētos "eagle", leuco "white", '
          'cephalos "head") is a bird of prey found in North America. A sea eagle, it has two known subspecies and '
          'forms a species pair with the white-tailed eagle (Haliaeetus albicilla). Its range includes most of '
          'Canada and Alaska, all of the contiguous United States, and northern Mexico. It is found near large '
          'bodies of open water with an abundant food supply and old-growth...',
  'image_file': 'bald_eagle.jpg',
  'image_caption': '(c) Addy, some rights reserved (CC BY-NC), uploaded by Addy · via iNaturalist (offline '
                   'sample)'},
 {'catalog_tag': 'iNaturalist · Aptenodytes forsteri',
  'class_tag': 'Bird',
  'name': 'Emperor Penguin',
  'subtitle': 'Aptenodytes forsteri · Aves',
  'headline': {'dot': '#8A9E3E', 'text': 'Near Threatened', 'tail': 'IUCN Red List (NT)'},
  'grid': [['Order', 'Penguins', ''],
           ['Family', 'Penguins', ''],
           ['Observations logged', '625', 'iNaturalist'],
           ['Conservation status', 'Near Threatened', 'IUCN NT']],
  'foot': 'The emperor penguin (Aptenodytes forsteri) is the tallest and heaviest of all living penguin species '
          'and is endemic to Antarctica. The male and female are similar in plumage and size, reaching 122\xa0cm '
          '(48\xa0in) in height and weighing from 22 to 45\xa0kg (49 to 99\xa0lb). The dorsal side and head are '
          'black and sharply delineated from the white belly, pale-yellow breast and bright-yellow ear patches. '
          'Like all penguins it is flightless, with a streamlined body, and wings...',
  'image_file': 'emperor_penguin.jpg',
  'image_caption': '(c) Martha de Jong-Lantink, some rights reserved (CC BY-NC-ND) · via iNaturalist (offline '
                   'sample)'},
 {'catalog_tag': 'iNaturalist · Cyanocitta cristata',
  'class_tag': 'Bird',
  'name': 'Blue Jay',
  'subtitle': 'Cyanocitta cristata · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Perching Birds', ''],
           ['Family', 'Crows, Jays, and Magpies', ''],
           ['Observations logged', '199,964', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The blue jay (Cyanocitta cristata) is a passerine bird in the family Corvidae, native to North America. '
          'It is resident through most of eastern and central United States, although western populations may be '
          'migratory. Resident populations are also found in Newfoundland, Canada, while breeding populations can '
          'be found in southern Canada. It breeds in both deciduous and coniferous forests, and is common near and '
          'in residential areas. It is predominantly blue with a white chest and...',
  'image_file': 'blue_jay.jpg',
  'image_caption': '(c) Judy Gallagher, some rights reserved (CC BY) · via iNaturalist (offline sample)'},
 {'catalog_tag': 'iNaturalist · Turdus migratorius',
  'class_tag': 'Bird',
  'name': 'American Robin',
  'subtitle': 'Turdus migratorius · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Perching Birds', ''],
           ['Family', 'Thrushes', ''],
           ['Observations logged', '461,200', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The American robin (Turdus migratorius) is a migratory songbird of the true thrush genus and Turdidae, '
          'the wider thrush family. It is named after the European robin because of its reddish-orange breast, '
          'though the two species are not closely related, with the European robin belonging to the Old World '
          'flycatcher family. The American robin is widely distributed throughout North America, wintering from '
          'southern Canada to central Mexico and along the Pacific Coast. It is the state...',
  'image_file': 'american_robin.jpg',
  'image_caption': '(c) John D Reynolds, some rights reserved (CC BY-NC), uploaded by John D Reynolds · via '
                   'iNaturalist (offline sample)'},
 {'catalog_tag': 'iNaturalist · Fratercula arctica',
  'class_tag': 'Bird',
  'name': 'Atlantic Puffin',
  'subtitle': 'Fratercula arctica · Aves',
  'headline': {'dot': '#C97F2E', 'text': 'Vulnerable', 'tail': 'IUCN Red List (VU)'},
  'grid': [['Order', 'Shorebirds and Allies', ''],
           ['Family', 'Auks, Murres, and Puffins', ''],
           ['Observations logged', '16,287', 'iNaturalist'],
           ['Conservation status', 'Vulnerable', 'IUCN VU']],
  'foot': 'The Atlantic puffin (Fratercula arctica), also known as the common puffin, is a species of seabird in '
          'the auk family. It is the only puffin native to the Atlantic Ocean; two related species, the tufted '
          'puffin and the horned puffin, are found in the northeastern Pacific. The Atlantic puffin breeds in '
          'Iceland, Norway, Greenland, Newfoundland and many North Atlantic islands, and as far south as Maine in '
          'the west and the west coast of Ireland and...',
  'image_file': 'atlantic_puffin.jpg',
  'image_caption': '(c) Paul Steeves, some rights reserved (CC BY-NC), uploaded by Paul Steeves · via iNaturalist '
                   '(offline sample)'},
 {'catalog_tag': 'iNaturalist · Bubo virginianus',
  'class_tag': 'Bird',
  'name': 'Great Horned Owl',
  'subtitle': 'Bubo virginianus · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Owls', ''],
           ['Family', 'Typical Owls', ''],
           ['Observations logged', '89,824', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The great horned owl (Bubo virginianus), also known as the tiger owl (originally derived from early '
          'naturalists\' description as the "winged tiger" or "tiger of the air") or the hoot owl, is a large owl '
          'native to the Americas. It is an extremely adaptable bird with a vast range and is the most widely '
          'distributed true owl in the Americas. Its primary diet is rabbits and hares, rats and mice and voles, '
          'although it freely hunts...',
  'image_file': 'great_horned_owl.jpg',
  'image_caption': '(c) Paul G. Johnson, some rights reserved (CC BY-NC-SA), uploaded by Paul G. Johnson · via '
                   'iNaturalist (offline sample)'},
 {'catalog_tag': 'iNaturalist · Bubo scandiacus',
  'class_tag': 'Bird',
  'name': 'Snowy Owl',
  'subtitle': 'Bubo scandiacus · Aves',
  'headline': {'dot': '#C97F2E', 'text': 'Vulnerable', 'tail': 'IUCN Red List (VU)'},
  'grid': [['Order', 'Owls', ''],
           ['Family', 'Typical Owls', ''],
           ['Observations logged', '12,485', 'iNaturalist'],
           ['Conservation status', 'Vulnerable', 'IUCN VU']],
  'foot': 'The snowy owl (Bubo scandiacus) is a large, white owl of the typical owl family. Snowy owls are native '
          'to Arctic regions in North America and Eurasia. Males are almost all white, while females have more '
          'flecks of black plumage. Juvenile snowy owls have black feathers until they turn white. The snowy owl '
          'is a ground nester that primarily hunts rodents and waterfowl, and opportunistically eats carrion. Most '
          'owls sleep during the day and hunt at...',
  'image_file': 'snowy_owl.jpg',
  'image_caption': '(c) Patrick Randall, some rights reserved (CC BY-NC-SA) · via iNaturalist (offline sample)'},
 {'catalog_tag': 'iNaturalist · Ramphastos toco',
  'class_tag': 'Bird',
  'name': 'Toco Toucan',
  'subtitle': 'Ramphastos toco · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Woodpeckers and Allies', ''],
           ['Family', 'Toucans', ''],
           ['Observations logged', '8,325', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The toco toucan (Ramphastos toco), also known as the common toucan or giant toucan, is the largest and '
          'probably the best known species in the toucan family. It is found in semi-open habitats throughout a '
          'large part of central and eastern South America. It is a common attraction in zoos.',
  'image_file': 'toco_toucan.jpg',
  'image_caption': '(c) Paul Steeves, some rights reserved (CC BY-NC), uploaded by Paul Steeves · via iNaturalist '
                   '(offline sample)'},
 {'catalog_tag': 'iNaturalist · Ara macao',
  'class_tag': 'Bird',
  'name': 'Scarlet Macaw',
  'subtitle': 'Ara macao · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Parrots', ''],
           ['Family', 'New World and African Parrots', ''],
           ['Observations logged', '14,078', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The scarlet macaw (Ara macao) is a large red, yellow, and blue South American parrot, a member of a '
          'large group of Neotropical parrots called macaws. It is native to humid evergreen forests of tropical '
          'South America. Range extends from south-eastern Mexico to the Peruvian Amazon, Colombia, Bolivia, '
          'Venezuela and Brazil in lowlands up to 500\xa0m (1,640\xa0ft) (at least formerly) up to 1,000\xa0m '
          '(3,281\xa0ft). It has suffered from local extinction through habitat destruction and capture for...',
  'image_file': 'scarlet_macaw.jpg',
  'image_caption': '(c) Matthew Patchett, some rights reserved (CC BY-NC), uploaded by Matthew Patchett · via '
                   'iNaturalist (offline sample)'},
 {'catalog_tag': 'iNaturalist · Tyto furcata',
  'class_tag': 'Bird',
  'name': 'American Barn Owl',
  'subtitle': 'Tyto furcata · Aves',
  'headline': {'dot': '#5B6FA6', 'text': '18,085 observations', 'tail': 'logged on iNaturalist'},
  'grid': [['Order', 'Owls', ''],
           ['Family', 'Barn Owls', ''],
           ['Observations logged', '18,085', 'iNaturalist'],
           ['Conservation status', None, 'na']],
  'foot': 'The American barn owl (Tyto furcata) is usually considered a subspecies group and together with the '
          'western barn owl group, the eastern barn owl group, and sometimes the Andaman masked owl, make up the '
          'barn owl, cosmopolitan in range. The barn owl is recognized by most taxonomic authorities. A few '
          "(including the International Ornithologists' Union) separate them into distinct species, as is done "
          'here. The American barn owl is native to North and South America, and',
  'image_file': 'american_barn_owl.jpg',
  'image_caption': '(c) Andria E. Hernandez, some rights reserved (CC BY-NC), uploaded by Andria E. Hernandez · '
                   'via iNaturalist (offline sample)'},
 {'catalog_tag': 'iNaturalist · Phoenicopterus roseus',
  'class_tag': 'Bird',
  'name': 'Greater Flamingo',
  'subtitle': 'Phoenicopterus roseus · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Flamingos', ''],
           ['Family', 'Flamingos', ''],
           ['Observations logged', '37,793', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': 'The greater flamingo (Phoenicopterus roseus) is the most widespread species of the flamingo family. It '
          'is found in Africa, on the Indian subcontinent, in the Middle East and southern Europe.',
  'image_file': 'greater_flamingo.jpg',
  'image_caption': '(c) cog2022, some rights reserved (CC BY-NC) · via iNaturalist (offline sample)'},
 {'catalog_tag': 'iNaturalist · Balaeniceps rex',
  'class_tag': 'Bird',
  'name': 'Shoebill',
  'subtitle': 'Balaeniceps rex · Aves',
  'headline': {'dot': '#C97F2E', 'text': 'Vulnerable', 'tail': 'IUCN Red List (VU)'},
  'grid': [['Order', 'Pelicans, Herons, Ibises, and Allies', ''],
           ['Family', 'Shoebills', ''],
           ['Observations logged', '583', 'iNaturalist'],
           ['Conservation status', 'Vulnerable', 'IUCN VU']],
  'foot': 'The shoebill (Balaeniceps rex) also known as whalehead or shoe-billed stork, is a very large stork-like '
          'bird. It derives its name from its massive shoe-shaped bill. Although it has a somewhat stork-like '
          'overall form and has previously been classified in the order Ciconiiformes, its true affiliations with '
          'other living birds is ambiguous. Some authorities now reclassify it with the Pelecaniformes. The adult '
          'is mainly grey while the juveniles are browner. It lives in tropical east Africa...',
  'image_file': 'shoebill.jpg',
  'image_caption': '(c) Nik Borrow, some rights reserved (CC BY-NC), uploaded by Nik Borrow · via iNaturalist '
                   '(offline sample)'},
 {'catalog_tag': 'iNaturalist · Calypte anna',
  'class_tag': 'Bird',
  'name': "Anna's Hummingbird",
  'subtitle': 'Calypte anna · Aves',
  'headline': {'dot': '#3E8F63', 'text': 'Least Concern', 'tail': 'IUCN Red List (LC)'},
  'grid': [['Order', 'Swifts and Hummingbirds', ''],
           ['Family', 'Hummingbirds', ''],
           ['Observations logged', '126,055', 'iNaturalist'],
           ['Conservation status', 'Least Concern', 'IUCN LC']],
  'foot': "Anna's hummingbird (Calypte anna), a medium-sized hummingbird native to the west coast of North "
          "America, was named after Anna Masséna, Duchess of Rivoli. In the early 20th century, Anna's "
          'hummingbirds bred only in northern Baja California and southern California. The transplanting of exotic '
          'ornamental plants in residential areas throughout the Pacific coast and inland deserts provided '
          'expanded nectar and nesting sites, allowing the species to expand its breeding range.',
  'image_file': 'anna_s_hummingbird.jpg',
  'image_caption': '(c) Nancy Christensen, all rights reserved, uploaded by Nancy Christensen · via iNaturalist '
                   '(offline sample)'}]


def load_offline_card() -> dict:
    import random
    card = dict(random.choice(OFFLINE_CARDS))
    card["category"] = "bird"
    card["tier"] = "offline"
    card["image_is_cat"] = False
    card["grid"] = [tuple(row) for row in card["grid"]]
    image_path = ASSETS_DIR / card.pop("image_file")
    card["image_bytes"] = image_path.read_bytes() if image_path.exists() else None
    return card
