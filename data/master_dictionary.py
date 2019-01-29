# The master dictionary holds the starting sequences that are used for the individual blastp searches
# The plan is to cover all members of the VEGF family:
#
# VEGF-A, PlGF, VEGF-B, VEGF-C, VEGF-D, VEGF-E, VEGF-F
#
# The keys in the master_dictionary is the unique id of the sequence (used for retrieval)
# The values in the master dictionary are lists:
#
# 0 image (name of the svg file)
# 1 common species name
# 2 scientific species name
# 3 animal "class"
# 4 protein_id [fasta, gene_id]
# 5 outgroup [protein_id, gene_id, protein, scientific species name]
# 6 alignment-trimming [start, end]
# 7 optional local amino acid sequence (if the sequences are not available online)
# 8 number of blasts hist to return (default 100)
# 9 list of alternative names
# 10 list of close relatives that can results in blast hits, that should be excluded
#
# For some entries, one needs perhaps to find a better outgroup
# It might also be better to blast some of the proteins with an incomplete sequence: E.g. the VEGF-B blast should be performed only with the isoform-specific SEQUENCES
#
# The gi number (the second id after the accesion number.version) is only important for the outgroup and locally given sequences (in which case any arbitrary string can be given)
#
{'VEGF-C':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_005420.1", "4885653"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [55, 227], None, 20000, ['VEGF-C', 'vascular endothelial growth factor C'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']], # trimming: end of N-terminal propeptide + VHD, outgroup: starfish VEGF-C
 'VEGF-D':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_004460.1", "4758378"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [55, 205], None, 20000, ['VEGF-D', 'vascular endothelial growth factor D'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']], # trimming: end of N-terminal propeptide + VHD, outgroup: starfish VEGF-C
 'VEGF-A206':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_001165094.1", "284172459"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [184, 232], None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']], # trimming: Heparin binding domain, outgroup:
 'VEGF-A165':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_001165097.1", "284172465"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [143, 191], None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']], # trimming: Heparin binding domain
 'VEGF-A121':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_001165099.1", "284172469"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [50, 132], None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']], # trimming: VHD
 'VEGF-B186':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_003368.1", "4507887"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [137, 207], None, 20000, ['VEGF-B', 'vascular endothelial growth factor B'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']], # trimming: 186-specific sequences
 'VEGF-B167':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_001230662.1", "344179103"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [137, 188], None, 20000, ['VEGF-B', 'vascular endothelial growth factor B'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']], # trimming: 167-specific sequences
 'PlGF-1':["Homo_sapiens.svg", "Human", "Homo_sapiens", "eutheria", ["NP_001280572.1", "20149543"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [50, 132], None, 20000, ['PlGF', 'placenta growth factor'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'uncharacterized protein', 'hypothetical protein']]} # trimming: VHD
 #'VEGF-F':["lepidosauria.svg", "Western sand viper", "Vipera ammodytes ammodytes", "lepidosauria", ["APB93447.1", None], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [25, 145], None, 20000, ['VEGF-F', 'toxin'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', ['uncharacterized protein', 'hypothetical protein']]]} # trimming: none
 #'PlGF-3':["Homo_sapiens.svg", "Human", "Homo_sapiens", "Mammalia", ["ENST00000405431.2", "local"], ["XP_022098898.1", "1229171207", "VEGF-C", "Acanthaster planci"], [132, 221], None, 50]} # trimming: Heparin binding domain
