# Dunn et al. 2014
# Klieve et al. 2015 (mammalia)
# Lee at al. 2013 (arthropoda)
#
# 'name':[NCBI taxon id, number of species in NCBI protein database, common name, number of protein sequences for that species in the NCBI nr protein database]
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
{'aves': [8782, 9339, 'birds', 2582703], 'chaetognatha': [10229, 56, 'arrow worms', 2167], 'pycnogonida': [57294, 185, 'sea spiders', 2057], 'actinopterygii': [7898, 18329, 'ray-finned fishes', 2380000], 'arachnida': [6854, 9644, 'spiders', 646367], 'hexapoda': [6960, 106986, 'insects', 7186866], 'amphibia': [8292, 5520, '', 478348], 'cyclostomata': [1476529, 103, 'hagfish/lamprey', 8221], 'hemichordata': [10219, 39, 'acorn wormws', 23424], 'chondrichthyes': [7777, 812, 'cartilaginous fishes', 114763], 'micrognathozoa': [195505, 1, '?', 2], 'entoprocta': [43120, 23, '?', 155], 'crustacea': [6657, 10373, '', 946524], 'priapulida': [33467, 7, 'penis worms', 20843], 'monotremata': [9255, 5, '', 25687], 'gnathostomulida': [66780, 21, 'jaw worms', 79], 'orthonectida': [33209, 2, '?', 8724], 'placozoa': [10226, 1, '?', 35645], 'xiphosura': [6845, 5, 'horseshoe crabs', 39494], 'echinodermata': [7586, 1724, 'sea stars/urchins/cucumbers/lilies', 135526], 'kinorhyncha': [51516, 56, 'mud dragons', 436], 'eutheria': [9347, 4693, 'placentals', 8003008], 'myriapoda': [61985, 949, 'millipeds', 7435], 'nemertea': [6217, 255, 'ribbon worms', 4618], 'loricifera': [310840, 1, '?', 1], 'mollusca': [6447, 13501, '', 742099], 'gastrotricha': [33313, 118, 'hairybacks', 389], 'onychophora': [2074142, 89, 'velvet worms', 1713], 'testudines': [8459, 353, 'turtles', 184246], 'nematoda': [6231, 3356, 'roundworms', 1550206], 'dicyemida': [10215, 24, '?', 150], 'annelida': [6340, 3170, 'segmented worms', 128675], 'cephalochordata': [7735, 11, 'lancelets/amphioxus', 95424], 'ctenophora': [1003038, 57, 'comb jellies', 6], 'brachiopoda': [7568, 100, 'lamp shells', 42113], 'nematomorpha': [33310, 28, 'horsehair worms', 368], 'platyhelminthes': [6157, 4171, 'flatworms', 560567], 'xenacoelomorpha': [1312402, 150, '?', 925], 'dipnoi': [7878, 6, 'lungfishes', 1321], 'rotifera': [10190, 233, 'wheel animals', 64483], 'metatheria': [9263, 324, 'marsupials', 142149], 'crocodylia': [1294634, 24, 'crocodiles', 179029], 'lepidosauria': [8504, 6803, 'lizards/snakes', 563435], 'tunicata': [7712, 353, 'ascidians/sea squirts/tunicates/sea pork/sea livers/sea tulips', 63726], 'cnidaria': [6073, 3586, 'medusae/polyps', 114763], 'phoroniformea': [120557, 14, 'horseshoe worms', 165], 'tardigrada': [42241, 168, '', 45754], 'bryozoa': [10205, 308, 'moss animals', 2724], 'porifera': [6040, 1318, 'sponges', 34074], 'cycliophora': [69815, 2, 'symbion', 278], 'coelacanthimorpha': [118072, 2, 'lobe-finned fishes', 34973]}
