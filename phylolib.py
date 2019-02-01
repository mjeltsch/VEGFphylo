def load_blacklist():
    # All the phyla in the list below have been found in the first screen not to have any VEGF-like molecules
    #blacklist = ['ctenophora', 'porifera', 'placozoa', 'xenacoelomorpha', 'cyclostomata', 'onychophora', 'pycnogonida', 'myriapoda', 'nematomorpha', 'loricifera', 'kinorhyncha', 'chaetognatha', 'bryozoa', 'entoprocta', 'cycliophora', 'nemertea', 'phoroniformea', 'gastrotricha', 'platyhelminthes', 'gnathostomulida', 'micrognathozoa', 'orthonectida', 'dicyemida']
    # All hits for the phyla below were manually checked and changed in taxon_data.py => therefore, we do not want to chenge them anymore
    blacklist = ['porifera', 'xenacoelomorpha', 'myriapoda', 'bryozoa', 'entoprocta', 'platyhelminthes', 'orthonectida']
    #blacklist = []
    return blacklist
