# Dunn et al. 2014
# Klieve et al. 2015 (mammalia)
# Lee at al. 2013 (arthropoda)
#
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
# phylum: [ncbi_taxon_id, number of animnal species in this phylum in the NCBI database, common name, number of protein sequences for that taxon in the NCBI nr protein database, {}, number of fully sequenced genomes in the phylum] 

{'rotifera': [10190, 233, 'wheel animals', 64483, {}, None],
 'onychophora': [27563, 89, 'velvet worms', 1713, {}, None],
 'platyhelminthes': [6157, 4171, 'flatworms', 560567, {}, None],
 'gastrotricha': [33313, 118, 'hairybacks', 389, {}, None],
 'priapulida': [33467, 7, 'penis worms', 20843, {}, None],
 'arachnida': [6854, 9644, 'spiders', 646367, {}, None],
 'kinorhyncha': [51516, 56, 'mud dragons', 436, {}, None],
 'dipnoi': [7878, 6, 'lungfishes', 1321, {}, None],
 'aves': [8782, 9339, 'birds', 2582703, {}, None],
 'hemichordata': [10219, 39, 'acorn wormws', 23424, {}, None],
 'coelacanthimorpha': [118072, 2, 'lobe-finned fishes', 34973, {}, None],
 'crocodylia': [1294634, 24, 'crocodiles', 179029, {}, None],
 'orthonectida': [33209, 2, '?', 8724, {}, None],
 'monotremata': [9255, 5, '', 25687, {}, None],
 'nematoda': [6231, 3356, 'roundworms', 1550206, {}, None],
 'cnidaria': [6073, 3586, 'medusae/polyps', 114763, {}, None],
 'amphibia': [8292, 5520, '', 478348, {}, None],
 'micrognathozoa': [195505, 1, '?', 2, {}, None],
 'mollusca': [6447, 13499, '', 742099, {}, None],
 'gnathostomulida': [66780, 21, 'jaw worms', 79, {}, None],
 'lepidosauria': [8504, 6803, 'lizards/snakes', 563435, {}, None],
 'cephalochordata': [7735, 11, 'lancelets', 95424, {}, None],
 'chondrichthyes': [7777, 812, 'cartilaginous fishes', 114763, {}, None],
 'cyclostomata': [1476529, 103, 'hagfish/lamprey', 8221, {}, None],
 'chaetognatha': [10229, 56, 'arrow worms', 2167, {}, None],
 'tardigrada': [42241, 168, '', 45754, {}, None],
 'tunicata': [7712, 353, '?', 63726, {}, None],
 'actinopterygii': [7898, 18329, 'ray-finned fishes', 2380000, {}, None],
 'xiphosura': [6845, 5, 'horseshoe crabs', 39494, {}, None],
 'ctenophora': [10197, 57, 'comb jellies', 6, {}, None],
 'bryozoa': [10205, 308, 'moss animals', 2724, {}, None],
 'pycnogonida': [57294, 185, 'sea spiders', 2057, {}, None],
 'brachiopoda': [7568, 100, 'lamp shells', 42113, {}, None],
 'hexapoda': [6960, 106986, 'insects', 7186866, {}, None],
 'nemertea': [6217, 255, 'ribbon worms', 4618, {}, None],
 'crustacea': [6657, 10373, '', 946524, {}, None],
 'placozoa': [10226, 1, '?', 35645, {}, None],
 'nematomorpha': [33310, 28, 'horsehair worms', 368, {}, None],
 'dicyemida': [10215, 24, '?', 150, {}, None],
 'xenacoelomorpha': [1312402, 150, '?', 925, {}, None],
 'myriapoda': [61985, 949, 'millipeds', 7435, {}, None],
 'metatheria': [9263, 324, 'marsupials', 142149, {}, None],
 'porifera': [6040, 1318, 'sponges', 34074, {}, None],
 'annelida': [6340, 3170, 'segmented worms', 128675, {}, None],
 'loricifera': [310840, 1, '?', 1, {}, None],
 'eutheria': [9347, 4692, 'placentals', 8003008, {}, None],
 'testudines': [8459, 353, 'turtles', 184246, {}, None],
 'phoroniformea': [120557, 14, 'horseshoe worms', 165, {}, None],
 'entoprocta': [43120, 23, '?', 155, {}, None],
 'cycliophora': [69815, 2, 'symbion', 278, {}, None],
 'echinodermata': [7586, 1724, '?', 135526, {}, None]}
