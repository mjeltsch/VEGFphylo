# The master dictionary holds the starting sequences that are used for the individual blastp searches
# The plan is to cover all members of the VEGF family:
#
# VEGF-A, PlGF, VEGF-B, VEGF-C, VEGF-D, VEGF-E, VEGF-F
#
# The value in the master dictionary is a list of five elements:
#
# 0 protein_id [fasta]
# 1 subrange (e.g. "50-200"; not used at the moment)
# 2 number of blasts hist to return (default is 50 if None is given)
# 3 list of alternative names (the first name is going to be the canonical ortholog group name)
# 4 list of close relatives that can results in blast hits, that should be excluded
#
{'VEGF-C':["NP_005420.1", None, 20000, ['VEGF-C', 'vascular endothelial growth factor C'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor B', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-D':["NP_004460.1", None, 20000, ['VEGF-D', 'vascular endothelial growth factor D'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'placenta growth factor']],
 'VEGF-A206':["NP_001165094.1", None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-A165':["NP_001165097.1", None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-A121':["NP_001165099.1", None, 20000, ['VEGF-A', 'vascular endothelial growth factor A'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-B186':["NP_003368.1", None, 20000, ['VEGF-B', 'vascular endothelial growth factor B'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'VEGF-B167':["NP_001230662.1", None, 20000, ['VEGF-B', 'vascular endothelial growth factor B'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'PlGF-1':["NP_001280572.1", None, 20000, ['PlGF', 'placenta growth factor'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D']],
 'VEGF-F':["APB93447.1", None, 20000, ['VEGF-F', 'toxin'], ['platelet-derived growth factor', 'PDGF and VEGF related factor', 'vascular endothelial growth factor A', 'vascular endothelial growth factor B', 'vascular endothelial growth factor C', 'vascular endothelial growth factor D', 'placenta growth factor']],
 'PDGF-A': ["NP_002598.4", None, 20000, ['PDGF-A', 'platelet-derived growth factor A'], ['platelet-derived growth factor B', 'platelet-derived growth factor C', 'platelet-derived growth factor D', 'PDGF and VEGF related factor', 'vascular endothelial growth factor', 'placenta growth factor']],
 'PDGF-B': ["NP_002599.1", None, 20000, ['PDGF-B', 'platelet-derived growth factor B'], ['platelet-derived growth factor A', 'platelet-derived growth factor C', 'platelet-derived growth factor D', 'PDGF and VEGF related factor', 'vascular endothelial growth factor', 'placenta growth factor']],
 'PDGF-C': ["NP_057289.1", None, 20000, ['PDGF-C', 'platelet-derived growth factor C'], ['platelet-derived growth factor A', 'platelet-derived growth factor B', 'platelet-derived growth factor D', 'PDGF and VEGF related factor', 'vascular endothelial growth factor', 'placenta growth factor']],
 'PDGF-D': ["NP_079484.1", None, 20000, ['PDGF-D', 'platelet-derived growth factor D'], ['platelet-derived growth factor A', 'platelet-derived growth factor B', 'platelet-derived growth factor C', 'PDGF and VEGF related factor', 'vascular endothelial growth factor', 'placenta growth factor']]}
