# Dunn et al. 2014
# Klieve et al. 2015 (mammalia)
# Lee at al. 2013 (arthropoda)
#
# 'name':[NCBI taxon id, number of species in NCBI protein database, common name, number of protein sequences for that species in the NCBI nr protein database, {VEGF-A: number of hits, etc.}]
# cyclostomata = agnatha
# phoroniformea = phoronida
# tunicata = urochordata
# If there is not common name for the phylum, but a well-known species, the common name of the species is given.
# If the phylum consists of several (better known) subphyla/classes, the common names of the subphyla/classes are given
#
# Sarcopterygii includes all tetrapods!!!!
#
# Seach term for downloading gi lists:
#
# sacropterygii: txid8287[Organism] NOT txid32523[Organism]
# The number of species is off, because it counts still all tetrapods
# One could substract tetrapods, which would leave 8 species => replaced
# by actinista (!= lobe-fined fishes, because Homo is highly derived lobe finned fish)
#
# The arthropod gi list for local analysis contains bacterial contamination (3064
# sequences from Morganella morganii)
#
# The mollusca gi list for local analysis contains bacterial contamination (2047
# sequences from Coxiella burnetii)
#
# placentalia = eutheria (eutheria comprise a extinct sister clade)
#
# xenocarida was removed since NCBI does not know about them, crustacea & hexapoda are shown in the cladogram as
# sister clades, which is not correct!
#
# Removed because of too big gi id filesize (github accepts max. 100 MB):
#
#'tetrapoda': [32523, 27063, ''],
#'arthropoda': [6656, 128142, '']
#'mammalia': [40674, 5022, '']
#'craniata': [115366, 46293, '?']
#
# echinodermata = sea stars/urchins/cucumbers/lilies
# tunicata = ascidians/sea squirts/tunicates/sea pork/sea livers/sea tulips
#
{'myriapoda': [61985, 949, 'millipeds', 7435, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 1],
 'xiphosura': [6845, 5, 'horseshoe crabs', 39494, {'VEGF-A165': 18, 'VEGF-B186': 18, 'VEGF-A206': 18, 'VEGF-A121': 18, 'VEGF-B167': 18, 'PlGF-1': 27, 'VEGF-D': 2, 'VEGF-C': 11}, 1],
 'tardigrada': [42241, 168, '', 45754, {'VEGF-A165': 6, 'VEGF-B186': 4, 'VEGF-A206': 3, 'VEGF-A121': 4, 'VEGF-B167': 4, 'PlGF-1': 2, 'VEGF-D': 4, 'VEGF-C': 4}, 2],
 'coelacanthimorpha': [118072, 2, 'lobe-finned fishes', 34973, {'VEGF-A165': 14, 'VEGF-B186': 18, 'VEGF-A206': 14, 'VEGF-A121': 15, 'VEGF-B167': 17, 'PlGF-1': 14, 'VEGF-D': 17, 'VEGF-C': 17}, 1],
 'annelida': [6340, 3170, 'segmented worms', 128675, {'VEGF-A165': 3, 'VEGF-B186': 2, 'VEGF-A206': 3, 'VEGF-A121': 3, 'VEGF-B167': 3, 'PlGF-1': 3, 'VEGF-D': 2, 'VEGF-C': 3}, 5],
 'cnidaria': [6073, 3586, 'medusae/polyps', 114763, {'VEGF-A165': 21, 'VEGF-B186': 6, 'VEGF-A206': 15, 'VEGF-A121': 13, 'VEGF-B167': 79, 'PlGF-1': 3, 'VEGF-D': 26, 'VEGF-C': 55}, 18],
 'lepidosauria': [8504, 6803, 'lizards/snakes', 563435, {'VEGF-A165': 149, 'VEGF-B186': 183, 'VEGF-A206': 151, 'VEGF-A121': 170, 'VEGF-B167': 168, 'PlGF-1': 162, 'VEGF-D': 185, 'VEGF-C': 178}, 21],
 'nematomorpha': [33310, 28, 'horsehair worms', 368, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'amphibia': [8292, 5520, '', 478348, {'VEGF-A165': 58, 'VEGF-B186': 77, 'VEGF-A206': 61, 'VEGF-A121': 66, 'VEGF-B167': 73, 'PlGF-1': 67, 'VEGF-D': 79, 'VEGF-C': 77}, 6],
 'eutheria': [9347, 4692, 'placentals', 8003008, {'VEGF-A165': 1308, 'VEGF-B186': 1763, 'VEGF-A206': 1303, 'VEGF-A121': 1340, 'VEGF-B167': 1680, 'PlGF-1': 1532, 'VEGF-D': 1774, 'VEGF-C': 1586}, 181],
 'placozoa': [10226, 1, '?', 35645, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 2],
 'brachiopoda': [7568, 100, 'lamp shells', 42113, {'VEGF-A165': 1, 'VEGF-B186': 1, 'VEGF-A206': 1, 'VEGF-A121': 1, 'VEGF-B167': 2, 'PlGF-1': 0, 'VEGF-D': 1, 'VEGF-C': 3}, 2],
 'micrognathozoa': [195505, 1, '?', 2, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'monotremata': [9255, 5, '', 25687, {'VEGF-A165': 7, 'VEGF-B186': 8, 'VEGF-A206': 7, 'VEGF-A121': 7, 'VEGF-B167': 8, 'PlGF-1': 6, 'VEGF-D': 8, 'VEGF-C': 8}, 1],
 'phoroniformea': [120557, 14, 'horseshoe worms', 165, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 1],
 'dicyemida': [10215, 24, '?', 150, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'ctenophora': [1003038, 57, 'comb jellies', 6, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'echinodermata': [7586, 1724, '?', 135526, {'VEGF-A165': 17, 'VEGF-B186': 18, 'VEGF-A206': 12, 'VEGF-A121': 12, 'VEGF-B167': 19, 'PlGF-1': 12, 'VEGF-D': 17, 'VEGF-C': 22}, 11],
 'arachnida': [6854, 9644, 'spiders', 646367, {'VEGF-A165': 39, 'VEGF-B186': 35, 'VEGF-A206': 34, 'VEGF-A121': 41, 'VEGF-B167': 38, 'PlGF-1': 65, 'VEGF-D': 17, 'VEGF-C': 33}, 27],
 'dipnoi': [7878, 6, 'lungfishes', 1321, {'VEGF-A165': 4, 'VEGF-B186': 4, 'VEGF-A206': 0, 'VEGF-A121': 4, 'VEGF-B167': 4, 'PlGF-1': 4, 'VEGF-D': 4, 'VEGF-C': 4}, 0],
 'crocodylia': [1294634, 24, 'crocodiles', 179029, {'VEGF-A165': 52, 'VEGF-B186': 63, 'VEGF-A206': 46, 'VEGF-A121': 57, 'VEGF-B167': 57, 'PlGF-1': 53, 'VEGF-D': 63, 'VEGF-C': 63}, 4],
 'mollusca': [6447, 13499, '', 742099, {'VEGF-A165': 13, 'VEGF-B186': 10, 'VEGF-A206': 9, 'VEGF-A121': 13, 'VEGF-B167': 15, 'PlGF-1': 3, 'VEGF-D': 11, 'VEGF-C': 16}, 26],
 'aves': [8782, 9339, 'birds', 2582703, {'VEGF-A165': 475, 'VEGF-B186': 772, 'VEGF-A206': 475, 'VEGF-A121': 539, 'VEGF-B167': 680, 'PlGF-1': 576, 'VEGF-D': 807, 'VEGF-C': 838}, 132],
 'kinorhyncha': [51516, 56, 'mud dragons', 436, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'testudines': [8459, 353, 'turtles', 184246, {'VEGF-A165': 86, 'VEGF-B186': 98, 'VEGF-A206': 67, 'VEGF-A121': 79, 'VEGF-B167': 74, 'PlGF-1': 72, 'VEGF-D': 98, 'VEGF-C': 99}, 10],
 'bryozoa': [10205, 308, 'moss animals', 2724, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'nematoda': [6231, 3356, 'roundworms', 1550206, {'VEGF-A165': 0, 'VEGF-B186': 21, 'VEGF-A206': 0, 'VEGF-A121': 0, 'VEGF-B167': 23, 'PlGF-1': 0, 'VEGF-D': 2, 'VEGF-C': 4}, 103],
 'cyclostomata': [1476529, 103, 'hagfish/lamprey', 8221, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 3],
 'crustacea': [6657, 10373, '', 946524, {'VEGF-A165': 21, 'VEGF-B186': 26, 'VEGF-A206': 16, 'VEGF-A121': 20, 'VEGF-B167': 29, 'PlGF-1': 12, 'VEGF-D': 15, 'VEGF-C': 22}, 25],
 'pycnogonida': [57294, 185, 'sea spiders', 2057, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'gastrotricha': [33313, 118, 'hairybacks', 389, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'orthonectida': [33209, 2, '?', 8724, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 1],
 'cycliophora': [69815, 2, 'symbion', 278, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'entoprocta': [43120, 23, '?', 155, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'actinopterygii': [7898, 18329, 'ray-finned fishes', 2380000, {'VEGF-A165': 1029, 'VEGF-B186': 1392, 'VEGF-A206': 1026, 'VEGF-A121': 1446, 'VEGF-B167': 1254, 'PlGF-1': 1143, 'VEGF-D': 1250, 'VEGF-C': 1350}, 186],
 'xenacoelomorpha': [1312402, 150, '?', 925, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'nemertea': [6217, 255, 'ribbon worms', 4618, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 1],
 'hexapoda': [6960, 106986, 'insects', 7186866, {'VEGF-A165': 318, 'VEGF-B186': 255, 'VEGF-A206': 213, 'VEGF-A121': 256, 'VEGF-B167': 292, 'PlGF-1': 196, 'VEGF-D': 93, 'VEGF-C': 132}, 339],
 'hemichordata': [10219, 39, 'acorn wormws', 23424, {'VEGF-A165': 2, 'VEGF-B186': 4, 'VEGF-A206': 2, 'VEGF-A121': 3, 'VEGF-B167': 3, 'PlGF-1': 1, 'VEGF-D': 4, 'VEGF-C': 4}, 2],
 'porifera': [6040, 1318, 'sponges', 34074, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 2],
 'chondrichthyes': [7777, 812, 'cartilaginous fishes', 114763, {'VEGF-A165': 25, 'VEGF-B186': 30, 'VEGF-A206': 24, 'VEGF-A121': 25, 'VEGF-B167': 29, 'PlGF-1': 25, 'VEGF-D': 28, 'VEGF-C': 29}, 6],
 'platyhelminthes': [6157, 4171, 'flatworms', 560567, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 26],
 'chaetognatha': [10229, 56, 'arrow worms', 2167, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'rotifera': [10190, 233, 'wheel animals', 64483, {'VEGF-A165': 1, 'VEGF-B186': 2, 'VEGF-A206': 0, 'VEGF-A121': 0, 'VEGF-B167': 1, 'PlGF-1': 1, 'VEGF-D': 1, 'VEGF-C': 1}, 6],
 'gnathostomulida': [66780, 21, 'jaw worms', 79, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'loricifera': [310840, 1, '?', 1, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'cephalochordata': [7735, 11, 'lancelets', 95424, {'VEGF-A165': 7, 'VEGF-B186': 7, 'VEGF-A206': 7, 'VEGF-A121': 7, 'VEGF-B167': 8, 'PlGF-1': 5, 'VEGF-D': 7, 'VEGF-C': 13}, 4],
 'onychophora': [2074142, 89, 'velvet worms', 1713, {'VEGF-A121': 0, 'VEGF-A165': 0, 'VEGF-B186': 0, 'VEGF-A206': 0, 'VEGF-B167': 0, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 0}, 0],
 'priapulida': [33467, 7, 'penis worms', 20843, {'VEGF-A165': 1, 'VEGF-B186': 1, 'VEGF-A206': 0, 'VEGF-A121': 0, 'VEGF-B167': 1, 'PlGF-1': 0, 'VEGF-D': 0, 'VEGF-C': 1}, 1],
 'tunicata': [7712, 353, '?', 63726, {'VEGF-A165': 1, 'VEGF-B186': 1, 'VEGF-A206': 1, 'VEGF-A121': 1, 'VEGF-B167': 1, 'PlGF-1': 2, 'VEGF-D': 1, 'VEGF-C': 1}, 6],
 'metatheria': [9263, 324, 'marsupials', 142149, {'VEGF-A165': 26, 'VEGF-B186': 40, 'VEGF-A206': 26, 'VEGF-A121': 32, 'VEGF-B167': 36, 'PlGF-1': 29, 'VEGF-D': 39, 'VEGF-C': 39}, 5]}
