# The master dictionary holds the starting sequences that are used for the individual blastp searches
# The plan is to cover all members of the VEGF family:
#
# VEGF-A, PlGF, VEGF-B, VEGF-C, VEGF-D, VEGF-E, VEGF-F, PDGF-A, PDGF-B, PDGF-C, PDGF-D
#
# The value in the master dictionary is a list:
#
# 0 protein_id [fasta]
# 1 subrange
# 2 number of blasts hist to return (default is 50 if None is given)
# 3 list of alternative names/synonyms
# 4 list of close relatives that can result in blast hits
#
{'VEGF-C': ["NP_005420.1", None, 20000, ['VEGF-C', 'vascular endothelial growth factor C'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor B', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-D': ["NP_004460.1", None, 20000, ['VEGF-D', 'vascular endothelial growth factor D'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'placenta growth factor']],
 'VEGF-A206': ["NP_001165094.1", None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-A165': ["NP_001165097.1", None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-A121': ["NP_001165099.1", None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-B186': ["NP_003368.1", None, 20000, ['VEGF-B', 'vascular endothelial growth factor B'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-B167': ["NP_001230662.1", None, 20000, ['VEGF-B', 'vascular endothelial growth factor B'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'PlGF-1': ["NP_001280572.1", None, 20000, ['PlGF', 'placenta growth factor'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D']]}

