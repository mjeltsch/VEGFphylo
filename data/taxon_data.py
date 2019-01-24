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
{'myriapoda': [61985, 949, 'millipeds', 7435, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 1],
 'xiphosura': [6845, 5, 'horseshoe crabs', 39494, {'VEGF-A165': [18, None], 'VEGF-B186': [18, None], 'VEGF-A206': [18, None], 'VEGF-A121': [18, None], 'VEGF-B167': [18, None], 'PlGF-1': [27, None], 'VEGF-D': [2, None], 'VEGF-C': [11, None]}, 1],
 'tardigrada': [42241, 168, '', 45754, {'VEGF-A165': [6, None], 'VEGF-B186': [4, None], 'VEGF-A206': [3, None], 'VEGF-A121': [4, None], 'VEGF-B167': [4, None], 'PlGF-1': [2, None], 'VEGF-D': [4, None], 'VEGF-C': [4, None]}, 2],
 'coelacanthimorpha': [118072, 2, 'lobe-finned fishes', 34973, {'VEGF-A165': [14, None], 'VEGF-B186': [18, None], 'VEGF-A206': [14, None], 'VEGF-A121': [15, None], 'VEGF-B167': [17, None], 'PlGF-1': [14, None], 'VEGF-D': [17, None], 'VEGF-C': [17, None]}, 1],
 'annelida': [6340, 3170, 'segmented worms', 128675, {'VEGF-A165': [3, None], 'VEGF-B186': [2, None], 'VEGF-A206': [3, None], 'VEGF-A121': [3, None], 'VEGF-B167': [3, None], 'PlGF-1': [3, None], 'VEGF-D': [2, None], 'VEGF-C': [3, None]}, 5],
 'cnidaria': [6073, 3586, 'medusae/polyps', 114763, {'VEGF-A165': [21, None], 'VEGF-B186': [6, None], 'VEGF-A206': [15, None], 'VEGF-A121': [13, None], 'VEGF-B167': [79, None], 'PlGF-1': [3, None], 'VEGF-D': [26, None], 'VEGF-C': [55, None]}, 18],
 'lepidosauria': [8504, 6803, 'lizards/snakes', 563435, {'VEGF-A165': [149, None], 'VEGF-B186': [183, None], 'VEGF-A206': [151, None], 'VEGF-A121': [170, None], 'VEGF-B167': [168, None], 'PlGF-1': [162, None], 'VEGF-D': [185, None], 'VEGF-C': [178, None]}, 21],
 'nematomorpha': [33310, 28, 'horsehair worms', 368, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'amphibia': [8292, 5520, '', 478348, {'VEGF-A165': [58, None], 'VEGF-B186': [77, None], 'VEGF-A206': [61, None], 'VEGF-A121': [66, None], 'VEGF-B167': [73, None], 'PlGF-1': [67, None], 'VEGF-D': [79, None], 'VEGF-C': [77, None]}, 6],
 'eutheria': [9347, 4692, 'placentals', 8003008, {'VEGF-A165': [1308, None], 'VEGF-B186': [1763, None], 'VEGF-A206': [1303, None], 'VEGF-A121': [1340, None], 'VEGF-B167': [1680, None], 'PlGF-1': [1532, None], 'VEGF-D': [1774, None], 'VEGF-C': [1586, None]}, 181],
 'placozoa': [10226, 1, '?', 35645, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 2],
 'brachiopoda': [7568, 100, 'lamp shells', 42113, {'VEGF-A165': [1, None], 'VEGF-B186': [1, None], 'VEGF-A206': [1, None], 'VEGF-A121': [1, None], 'VEGF-B167': [2, None], 'PlGF-1': [0, None], 'VEGF-D': [1, None], 'VEGF-C': [3, None]}, 2],
 'micrognathozoa': [195505, 1, '?', 2, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'monotremata': [9255, 5, '', 25687, {'VEGF-A165': [7, None], 'VEGF-B186': [8, None], 'VEGF-A206': [7, None], 'VEGF-A121': [7, None], 'VEGF-B167': [8, None], 'PlGF-1': [6, None], 'VEGF-D': [8, None], 'VEGF-C': [8, None]}, 1],
 'phoroniformea': [120557, 14, 'horseshoe worms', 165, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 1],
 'dicyemida': [10215, 24, '?', 150, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'ctenophora': [1003038, 57, 'comb jellies', 6, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'echinodermata': [7586, 1724, '?', 135526, {'VEGF-A165': [17, None], 'VEGF-B186': [18, None], 'VEGF-A206': [12, None], 'VEGF-A121': [12, None], 'VEGF-B167': [19, None], 'PlGF-1': [12, None], 'VEGF-D': [17, None], 'VEGF-C': [22, None]}, 11],
 'arachnida': [6854, 9644, 'spiders', 646367, {'VEGF-A165': [39, None], 'VEGF-B186': [35, None], 'VEGF-A206': [34, None], 'VEGF-A121': [41, None], 'VEGF-B167': [38, None], 'PlGF-1': [65, None], 'VEGF-D': [17, None], 'VEGF-C': [33, None]}, 27],
 'dipnoi': [7878, 6, 'lungfishes', 1321, {'VEGF-A165': [4, None], 'VEGF-B186': [4, None], 'VEGF-A206': [0, None], 'VEGF-A121': [4, None], 'VEGF-B167': [4, None], 'PlGF-1': [4, None], 'VEGF-D': [4, None], 'VEGF-C': [4, None]}, 0],
 'crocodylia': [1294634, 24, 'crocodiles', 179029, {'VEGF-A165': [52, None], 'VEGF-B186': [63, None], 'VEGF-A206': [46, None], 'VEGF-A121': [57, None], 'VEGF-B167': [57, None], 'PlGF-1': [53, None], 'VEGF-D': [63, None], 'VEGF-C': [63, None]}, 4],
 'mollusca': [6447, 13499, '', 742099, {'VEGF-A165': [13, None], 'VEGF-B186': [10, None], 'VEGF-A206': [9, None], 'VEGF-A121': [13, None], 'VEGF-B167': [15, None], 'PlGF-1': [3, None], 'VEGF-D': [11, None], 'VEGF-C': [16, None]}, 26],
 'aves': [8782, 9339, 'birds', 2582703, {'VEGF-A165': [475, None], 'VEGF-B186': [772, None], 'VEGF-A206': [475, None], 'VEGF-A121': [539, None], 'VEGF-B167': [680, None], 'PlGF-1': [576, None], 'VEGF-D': [807, None], 'VEGF-C': [838, None]}, 132],
 'kinorhyncha': [51516, 56, 'mud dragons', 436, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'testudines': [8459, 353, 'turtles', 184246, {'VEGF-A165': [86, None], 'VEGF-B186': [98, None], 'VEGF-A206': [67, None], 'VEGF-A121': [79, None], 'VEGF-B167': [74, None], 'PlGF-1': [72, None], 'VEGF-D': [98, None], 'VEGF-C': [99, None]}, 10],
 'bryozoa': [10205, 308, 'moss animals', 2724, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'nematoda': [6231, 3356, 'roundworms', 1550206, {'VEGF-A165': [0, None], 'VEGF-B186': [21, None], 'VEGF-A206': [0, None], 'VEGF-A121': [0, None], 'VEGF-B167': [23, None], 'PlGF-1': [0, None], 'VEGF-D': [2, None], 'VEGF-C': [4, None]}, 103],
 'cyclostomata': [1476529, 103, 'hagfish/lamprey', 8221, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 3],
 'crustacea': [6657, 10373, '', 946524, {'VEGF-A165': [21, None], 'VEGF-B186': [26, None], 'VEGF-A206': [16, None], 'VEGF-A121': [20, None], 'VEGF-B167': [29, None], 'PlGF-1': [12, None], 'VEGF-D': [15, None], 'VEGF-C': [22, None]}, 25],
 'pycnogonida': [57294, 185, 'sea spiders', 2057, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'gastrotricha': [33313, 118, 'hairybacks', 389, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'orthonectida': [33209, 2, '?', 8724, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 1],
 'cycliophora': [69815, 2, 'symbion', 278, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'entoprocta': [43120, 23, '?', 155, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'actinopterygii': [7898, 18329, 'ray-finned fishes', 2380000, {'VEGF-A165': [1029, None], 'VEGF-B186': [1392, None], 'VEGF-A206': [1026, None], 'VEGF-A121': [1446, None], 'VEGF-B167': [1254, None], 'PlGF-1': [1143, None], 'VEGF-D': [1250, None], 'VEGF-C': [1350, None]}, 186],
 'xenacoelomorpha': [1312402, 150, '?', 925, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'nemertea': [6217, 255, 'ribbon worms', 4618, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 1],
 'hexapoda': [6960, 106986, 'insects', 7186866, {'VEGF-A165': [318, None], 'VEGF-B186': [255, None], 'VEGF-A206': [213, None], 'VEGF-A121': [256, None], 'VEGF-B167': [292, None], 'PlGF-1': [196, None], 'VEGF-D': [93, None], 'VEGF-C': [132, None]}, 339],
 'hemichordata': [10219, 39, 'acorn wormws', 23424, {'VEGF-A165': [2, None], 'VEGF-B186': [4, None], 'VEGF-A206': [2, None], 'VEGF-A121': [3, None], 'VEGF-B167': [3, None], 'PlGF-1': [1, None], 'VEGF-D': [4, None], 'VEGF-C': [4, None]}, 2],
 'porifera': [6040, 1318, 'sponges', 34074, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 2],
 'chondrichthyes': [7777, 812, 'cartilaginous fishes', 114763, {'VEGF-A165': [25, None], 'VEGF-B186': [30, None], 'VEGF-A206': [24, None], 'VEGF-A121': [25, None], 'VEGF-B167': [29, None], 'PlGF-1': [25, None], 'VEGF-D': [28, None], 'VEGF-C': [29, None]}, 6],
 'platyhelminthes': [6157, 4171, 'flatworms', 560567, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 26],
 'chaetognatha': [10229, 56, 'arrow worms', 2167, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'rotifera': [10190, 233, 'wheel animals', 64483, {'VEGF-A165': [1, None], 'VEGF-B186': [2, None], 'VEGF-A206': [0, None], 'VEGF-A121': [0, None], 'VEGF-B167': [1, None], 'PlGF-1': [1, None], 'VEGF-D': [1, None], 'VEGF-C': [1, None]}, 6],
 'gnathostomulida': [66780, 21, 'jaw worms', 79, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'loricifera': [310840, 1, '?', 1, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'cephalochordata': [7735, 11, 'lancelets', 95424, {'VEGF-A165': [7, None], 'VEGF-B186': [7, None], 'VEGF-A206': [7, None], 'VEGF-A121': [7, None], 'VEGF-B167': [8, None], 'PlGF-1': [5, None], 'VEGF-D': [7, None], 'VEGF-C': [13, None]}, 4],
 'onychophora': [2074142, 89, 'velvet worms', 1713, {'VEGF-A121': [0, None], 'VEGF-A165': [0, None], 'VEGF-B186': [0, None], 'VEGF-A206': [0, None], 'VEGF-B167': [0, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [0, None]}, 0],
 'priapulida': [33467, 7, 'penis worms', 20843, {'VEGF-A165': [1, None], 'VEGF-B186': [1, None], 'VEGF-A206': [0, None], 'VEGF-A121': [0, None], 'VEGF-B167': [1, None], 'PlGF-1': [0, None], 'VEGF-D': [0, None], 'VEGF-C': [1, None]}, 1],
 'tunicata': [7712, 353, '?', 63726, {'VEGF-A165': [1, None], 'VEGF-B186': [1, None], 'VEGF-A206': [1, None], 'VEGF-A121': [1, None], 'VEGF-B167': [1, None], 'PlGF-1': [2, None], 'VEGF-D': [1, None], 'VEGF-C': [1, None]}, 6],
 'metatheria': [9263, 324, 'marsupials', 142149, {'VEGF-A165': [26, None], 'VEGF-B186': [40, None], 'VEGF-A206': [26, None], 'VEGF-A121': [32, None], 'VEGF-B167': [36, None], 'PlGF-1': [29, None], 'VEGF-D': [39, None], 'VEGF-C': [39, None]}, 5]}
