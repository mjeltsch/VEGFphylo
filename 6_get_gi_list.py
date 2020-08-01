#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script downloads a list of gi protein identifiers 
# for a specific species
# Info used:
# https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
# https://www.ncbi.nlm.nih.gov/books/NBK25500/

import Bio, csv, os, requests, time, subprocess, re
from phylolib import sanitize_species_name, execute_subprocess

def run(SPECIES_NAME):
    SPECIES_SANITIZED = sanitize_species_name(SPECIES_NAME)
    # Determine directories of script (in order to load & save the data files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    RESULT_DIR = '{0}/data/species_id_lists/'.format(APPLICATION_PATH)
    FILENAME = '{0}{1}.gi'.format(RESULT_DIR, SPECIES_SANITIZED)
    # Delete old gi list if it exists
    if os.path.isfile(FILENAME):
        os.remove(FILENAME)
    # First reequest is only to get the number of hits
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={0}[orgn]&&retstart={1}&retmax=0'.format(SPECIES_NAME, 0)
    print('Retrieving total number of protein sequences for {0}...'.format(SPECIES_NAME))
    r = requests.get(URL)
    # Extract from the second line the number between "<Count>" and "</Count><RetMax>" (= total hits)
    print('URL: {}'.format(URL))
    #print('Line 2: {0}'.format(r.content.decode('utf-8').splitlines()[2]))
    max = int(re.search('<Count>(.*)</Count><RetMax>', r.content.decode('utf-8').splitlines()[2]).group(1))
    print('{0} protein sequences.'.format(max))
    # Download species-specific id lists
    start_position = 0
    for start_position in range(0, max, 100000):
        URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={0}[orgn]&&retstart={1}&retmax=100000'.format(SPECIES_NAME, start_position)
        print('Retrieving {0} batch of the protein gi list for {1}...'.format(start_position/100000, SPECIES_NAME))
        # Write ids to a GI list file (analogous to downlaoding the pubmed results as a file with the GI list format)
        r = requests.get(URL)
        with open(FILENAME, 'a+') as file:
            print('Writing results to {0}'.format(FILENAME))
            for line in r.content.decode('utf-8').splitlines():
                if line[0:5] == '\t<Id>':
                    #print(line[5:-5])
                    file.write(line[5:-5]+'\n')
        #print('Waiting before initiating next download...')
        time.sleep(3)

    # Creating blast database aliases
    DATABASE = '{0}/../blast_database/nr'.format(APPLICATION_PATH)
    DATABASE_SPECIES = '{0}/../blast_database/nr_{1}'.format(APPLICATION_PATH, SPECIES_SANITIZED)
    DATABASE_TITLE = 'nr_{0}'.format(SPECIES_SANITIZED)
    comment = 'Creating a database subset for {0}:'.format(SPECIES_NAME)
    bash_command = 'blastdb_aliastool -gilist {0} -db {1} -out {2} -title {3}'.format(FILENAME, DATABASE, DATABASE_SPECIES, DATABASE_TITLE)
    execute_subprocess(comment, bash_command)
    print('\n')

if __name__ == '__main__':
    run('Drosophila subobscura')
