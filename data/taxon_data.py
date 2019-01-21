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
{'testudines': [8459, 353, 'turtles', 184246, {'VEGF-B167': 74, 'VEGF-A206': 67, 'VEGF-C': 99, 'VEGF-A165': 86, 'PlGF-1': 72, 'VEGF-A121': 79, 'VEGF-D': 98, 'VEGF-B186': 98}, 0],
 'amphibia': [8292, 5520, '', 478348, {'VEGF-B167': 73, 'VEGF-A206': 61, 'VEGF-C': 77, 'VEGF-A165': 58, 'PlGF-1': 67, 'VEGF-A121': 66, 'VEGF-D': 79, 'VEGF-B186': 77}, 6],
 'gnathostomulida': [66780, 21, 'jaw worms', 79, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'cyclostomata': [1476529, 103, 'hagfish/lamprey', 8221, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'aves': [8782, 9339, 'birds', 2582703, {'VEGF-B167': 680, 'VEGF-A206': 475, 'VEGF-C': 838, 'VEGF-A165': 475, 'PlGF-1': 576, 'VEGF-A121': 539, 'VEGF-D': 807, 'VEGF-B186': 772}, 132],
 'lepidosauria': [8504, 6803, 'lizards/snakes', 563435, {'VEGF-B167': 168, 'VEGF-A206': 151, 'VEGF-C': 178, 'VEGF-A165': 149, 'PlGF-1': 162, 'VEGF-A121': 170, 'VEGF-D': 185, 'VEGF-B186': 183}, 0],
 'xenacoelomorpha': [1312402, 150, '?', 925, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'actinopterygii': [7898, 18329, 'ray-finned fishes', 2380000, {'VEGF-B167': 1254, 'VEGF-A206': 1026, 'VEGF-C': 1350, 'VEGF-A165': 1029, 'PlGF-1': 1143, 'VEGF-A121': 1446, 'VEGF-D': 1250, 'VEGF-B186': 1392}, 0],
 'porifera': [6040, 1318, 'sponges', 34074, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'placozoa': [10226, 1, '?', 35645, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'monotremata': [9255, 5, '', 25687, {'VEGF-B167': 8, 'VEGF-A206': 7, 'VEGF-C': 8, 'VEGF-A165': 7, 'PlGF-1': 6, 'VEGF-A121': 7, 'VEGF-D': 8, 'VEGF-B186': 8}, 1],
 'cephalochordata': [7735, 11, 'lancelets/amphioxus', 95424, {'VEGF-B167': 8, 'VEGF-A206': 7, 'VEGF-C': 13, 'VEGF-A165': 7, 'PlGF-1': 5, 'VEGF-A121': 7, 'VEGF-D': 7, 'VEGF-B186': 7}, 0],
 'arachnida': [6854, 9644, 'spiders', 646367, {'VEGF-B167': 38, 'VEGF-A206': 34, 'VEGF-C': 33, 'VEGF-A165': 39, 'PlGF-1': 65, 'VEGF-A121': 41, 'VEGF-D': 17, 'VEGF-B186': 35}, 0],
 'onychophora': [2074142, 89, 'velvet worms', 1713, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'chondrichthyes': [7777, 812, 'cartilaginous fishes', 114763, {'VEGF-B167': 29, 'VEGF-A206': 24, 'VEGF-C': 29, 'VEGF-A165': 25, 'PlGF-1': 25, 'VEGF-A121': 25, 'VEGF-D': 28, 'VEGF-B186': 30}, 0],
 'echinodermata': [7586, 1724, 'sea stars/urchins/cucumbers/lilies', 135526, {'VEGF-B167': 19, 'VEGF-A206': 12, 'VEGF-C': 22, 'VEGF-A165': 17, 'PlGF-1': 12, 'VEGF-A121': 12, 'VEGF-D': 17, 'VEGF-B186': 18}, 0],
 'annelida': [6340, 3170, 'segmented worms', 128675, {'VEGF-B167': 3, 'VEGF-A206': 3, 'VEGF-C': 3, 'VEGF-A165': 3, 'PlGF-1': 3, 'VEGF-A121': 3, 'VEGF-D': 2, 'VEGF-B186': 2}, 0],
 'orthonectida': [33209, 2, '?', 8724, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'tunicata': [7712, 353, 'ascidians/sea squirts/tunicates/sea pork/sea livers/sea tulips', 63726, {'VEGF-B167': 1, 'VEGF-A206': 1, 'VEGF-C': 1, 'VEGF-A165': 1, 'PlGF-1': 2, 'VEGF-A121': 1, 'VEGF-D': 1, 'VEGF-B186': 1}, 0],
 'chaetognatha': [10229, 56, 'arrow worms', 2167, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'nematomorpha': [33310, 28, 'horsehair worms', 368, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'dipnoi': [7878, 6, 'lungfishes', 1321, {'VEGF-B167': 4, 'VEGF-A206': 0, 'VEGF-C': 4, 'VEGF-A165': 4, 'PlGF-1': 4, 'VEGF-A121': 4, 'VEGF-D': 4, 'VEGF-B186': 4}, 0],
 'ctenophora': [1003038, 57, 'comb jellies', 6, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'dicyemida': [10215, 24, '?', 150, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'entoprocta': [43120, 23, '?', 155, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'metatheria': [9263, 324, 'marsupials', 142149, {'VEGF-B167': 36, 'VEGF-A206': 26, 'VEGF-C': 39, 'VEGF-A165': 26, 'PlGF-1': 29, 'VEGF-A121': 32, 'VEGF-D': 39, 'VEGF-B186': 40}, 0],
 'xiphosura': [6845, 5, 'horseshoe crabs', 39494, {'VEGF-B167': 18, 'VEGF-A206': 18, 'VEGF-C': 11, 'VEGF-A165': 18, 'PlGF-1': 27, 'VEGF-A121': 18, 'VEGF-D': 2, 'VEGF-B186': 18}, 0],
 'loricifera': [310840, 1, '?', 1, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'nemertea': [6217, 255, 'ribbon worms', 4618, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'cycliophora': [69815, 2, 'symbion', 278, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'tardigrada': [42241, 168, '', 45754, {'VEGF-B167': 4, 'VEGF-A206': 3, 'VEGF-C': 4, 'VEGF-A165': 6, 'PlGF-1': 2, 'VEGF-A121': 4, 'VEGF-D': 4, 'VEGF-B186': 4}, 0],
 'coelacanthimorpha': [118072, 2, 'lobe-finned fishes', 34973, {'VEGF-B167': 17, 'VEGF-A206': 14, 'VEGF-C': 17, 'VEGF-A165': 14, 'PlGF-1': 14, 'VEGF-A121': 15, 'VEGF-D': 17, 'VEGF-B186': 18}, 0],
 'rotifera': [10190, 233, 'wheel animals', 64483, {'VEGF-B167': 1, 'VEGF-A206': 0, 'VEGF-C': 1, 'VEGF-A165': 1, 'PlGF-1': 1, 'VEGF-A121': 0, 'VEGF-D': 1, 'VEGF-B186': 2}, 0],
 'crocodylia': [1294634, 24, 'crocodiles', 179029, {'VEGF-B167': 57, 'VEGF-A206': 46, 'VEGF-C': 63, 'VEGF-A165': 52, 'PlGF-1': 53, 'VEGF-A121': 57, 'VEGF-D': 63, 'VEGF-B186': 63}, 0],
 'priapulida': [33467, 7, 'penis worms', 20843, {'VEGF-B167': 1, 'VEGF-A206': 0, 'VEGF-C': 1, 'VEGF-A165': 1, 'PlGF-1': 0, 'VEGF-A121': 0, 'VEGF-D': 0, 'VEGF-B186': 1}, 0],
 'myriapoda': [61985, 949, 'millipeds', 7435, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'platyhelminthes': [6157, 4171, 'flatworms', 560567, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 27],
 'kinorhyncha': [51516, 56, 'mud dragons', 436, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'micrognathozoa': [195505, 1, '?', 2, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'eutheria': [9347, 4692, 'placentals', 8003008, {'VEGF-B167': 1680, 'VEGF-A206': 1303, 'VEGF-C': 1586, 'VEGF-A165': 1308, 'PlGF-1': 1532, 'VEGF-A121': 1340, 'VEGF-D': 1774, 'VEGF-B186': 1763}, 0],
 'crustacea': [6657, 10373, '', 946524, {'VEGF-B167': 29, 'VEGF-A206': 16, 'VEGF-C': 22, 'VEGF-A165': 21, 'PlGF-1': 12, 'VEGF-A121': 20, 'VEGF-D': 15, 'VEGF-B186': 26}, 0],
 'hemichordata': [10219, 39, 'acorn wormws', 23424, {'VEGF-B167': 3, 'VEGF-A206': 2, 'VEGF-C': 4, 'VEGF-A165': 2, 'PlGF-1': 1, 'VEGF-A121': 3, 'VEGF-D': 4, 'VEGF-B186': 4}, 0],
 'nematoda': [6231, 3356, 'roundworms', 1550206, {'VEGF-B167': 23, 'VEGF-A206': 0, 'VEGF-C': 4, 'VEGF-A165': 0, 'PlGF-1': 0, 'VEGF-A121': 0, 'VEGF-D': 2, 'VEGF-B186': 21}, 0],
 'cnidaria': [6073, 3586, 'medusae/polyps', 114763, {'VEGF-B167': 79, 'VEGF-A206': 15, 'VEGF-C': 55, 'VEGF-A165': 21, 'PlGF-1': 3, 'VEGF-A121': 13, 'VEGF-D': 26, 'VEGF-B186': 6}, 0],
 'brachiopoda': [7568, 100, 'lamp shells', 42113, {'VEGF-B167': 2, 'VEGF-A206': 1, 'VEGF-C': 3, 'VEGF-A165': 1, 'PlGF-1': 0, 'VEGF-A121': 1, 'VEGF-D': 1, 'VEGF-B186': 1}, 0],
 'pycnogonida': [57294, 185, 'sea spiders', 2057, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'hexapoda': [6960, 106986, 'insects', 7186866, {'VEGF-B167': 292, 'VEGF-A206': 213, 'VEGF-C': 132, 'VEGF-A165': 318, 'PlGF-1': 196, 'VEGF-A121': 256, 'VEGF-D': 93, 'VEGF-B186': 255}, 339],
 'phoroniformea': [120557, 14, 'horseshoe worms', 165, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'bryozoa': [10205, 308, 'moss animals', 2724, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 1}, 0],
 'gastrotricha': [33313, 118, 'hairybacks', 389, {'PlGF-1': 0, 'VEGF-A206': 0, 'VEGF-C': 0, 'VEGF-A121': 0, 'VEGF-B167': 0, 'VEGF-A165': 0, 'VEGF-D': 0, 'VEGF-B186': 0}, 0],
 'mollusca': [6447, 13499, '', 742099, {'VEGF-B167': 15, 'VEGF-A206': 9, 'VEGF-C': 16, 'VEGF-A165': 13, 'PlGF-1': 3, 'VEGF-A121': 13, 'VEGF-D': 11, 'VEGF-B186': 10}, 0]}
