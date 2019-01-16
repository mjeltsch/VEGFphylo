# Dunn et al. 2014
# Klieve et al. 2015 (mammalia)
# Lee at al. 2013 (arthropoda)
#
# 'name':[NCBI taxon id, number of species in NCBI protein database]
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
# One could substract tetrapods, which would leave 8 species 
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
{'mollusca': [6447, 13501, ''], 'hexapoda': [6960, 106986, 'insects'], 'platyhelminthes': [6157, 4171, 'flatworms'], 'priapulida': [33467, 7, 'penis worms'], 'crustacea': [6657, 10373, ''], 'aves': [8782, 9339, 'birds'], 'arachnida': [6854, 9644, 'spiders'], 'lepidosauria': [8504, 6803, 'lizards/snakes'], 'mammalia': [40674, 5022, ''], 'dicyemida': [10215, 24, '?'], 'orthonectida': [33209, 2, '?'], 'onychophora': [2074142, 89, 'velvet worms'], 'amphibia': [8292, 5520, ''], 'arthropoda': [6656, 128142, ''], 'eutheria': [9347, 4693, 'placentals'], 'echinodermata': [7586, 1724, 'sea stars/urchins/cucumbers/lilies'], 'myriapoda': [61985, 949, 'millipeds'], 'crocodylia': [1294634, 24, 'crocodiles'], 'micrognathozoa': [195505, 1, '?'], 'actinopterygii': [7898, 18329, 'ray-finned fishes'], 'metatheria': [9263, 324, 'marsupials'], 'annelida': [6340, 3170, 'segmented worms'], 'phoroniformea': [120557, 14, 'horseshoe worms'], 'rotifera': [10190, 233, 'wheel animals'], 'ctenophora': [1003038, 57, 'comb jellies'], 'chondrichthyes': [7777, 812, 'cartilaginous fishes'], 'gastrotricha': [33313, 118, 'hairybacks'], 'gnathostomulida': [66780, 21, 'jaw worms'], 'sarcopterygii': [8287, 27071, 'lobe-finned fishes'], 'kinorhyncha': [51516, 56, 'mud dragons'], 'nematomorpha': [33310, 28, 'horsehair worms'], 'hemichordata': [10219, 39, 'acorn wormws'], 'testudines': [8459, 353, 'turtles'], 'xenacoelomorpha': [1312402, 150, '?'], 'nemertea': [6217, 255, 'ribbon worms'], 'craniata': [115366, 46293, '?'], 'nematoda': [6231, 3356, 'roundworms'], 'monotremata': [9255, 5, ''], 'cycliophora': [69815, 2, 'symbion'], 'placozoa': [10226, 1, '?'], 'loricifera': [310840, 1, '?'], 'pycnogonida': [57294, 185, 'sea spiders'], 'brachiopoda': [7568, 100, 'lamp shells'], 'cephalochordata': [7735, 11, 'lancelets/amphioxus'], 'cnidaria': [6073, 3586, 'medusae/polyps'], 'entoprocta': [43120, 23, '?'], 'tardigrada': [42241, 168, ''], 'chaetognatha': [10229, 56, 'arrow worms'], 'bryozoa': [10205, 308, 'moss animals'], 'tunicata': [7712, 353, 'ascidians/sea squirts/tunicates/sea pork/sea livers/sea tulips'], 'porifera': [6040, 1318, 'sponges'], 'xiphosura': [6845, 5, 'horseshoe crabs'], 'cyclostomata': [1476529, 103, 'hagfish/lamprey']}