#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script downloads a CSV file of all fully sequenced genomes
# and creates species-specific blast databases (aliases) in order
# to allow for local species-specific blasting.
# I used information form this thread to make it work: https://www.biostars.org/p/6528/.

import Bio, csv, os, requests, time, subprocess

def execute_subprocess(comment, bash_command, working_directory='.'):
    print("\n" + comment, bash_command)
    process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=working_directory)
    output, error = process.communicate()
    process_status = process.wait()
    if output.decode('utf-8') != '':
        print("Output: " + str(output))
    if error.decode('UTF-8') != '':
        print("Error: " + str(error))

def getCSVfile(CSV_FILE):
    # Open/download the list of fully sequenced genomes
    try:
        input_csv_file = csv.DictReader(open(CSV_FILE))
    except FileNotFoundError as err:
        print(err)
        # Download list of all fully sequences animal genomes:
        URL = 'https://www.ncbi.nlm.nih.gov/genomes/solr2txt.cgi?q=%5Bdisplay()%5D.from(GenomeBrowser).usingschema(%2Fschema%2FGenomeAssemblies).matching(group%3D%3D%5B%22Animals%22%5D)&fields=organism%7COrganism%20Name%2Clineage%7COrganism%20Groups%2Csize%7CSize(Mb)%2Cchromosomes%7CChromosomes%2Corganelles%7COrganelles%2Cplasmids%7CPlasmids%2Cassemblies%7CAssemblies&filename=genomes.csv&nolimit=on'
        r = requests.get(URL)
        with open('CSV_FILE', 'wb') as file:
            file.write(r.content)
        input_csv_file = csv.DictReader(open(CSV_FILE))
    return input_csv_file

def run():
    # Determine directories of script (in order to load & save the data files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    CSV_FILE = '{0}/data/genomes.csv'.format(APPLICATION_PATH)
    RESULT_DIR = '{0}/data/species_id_lists/'.format(APPLICATION_PATH)
    DATABASE = '{0}/../blast_database/nr'.format(APPLICATION_PATH)
    input_csv_file = getCSVfile(CSV_FILE)
    for line in input_csv_file:
        species_name = line['#Organism Name']
        species_name_ = species_name.replace(' ', '_')
        # Download species-specific id lists
        URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={0}[orgn]&retmax=100000000'.format(species_name)
        print('Retrieving protein id list for {0}...'.format(species_name))
        r = requests.get(URL)
        FILENAME = '{0}{1}.gi'.format(RESULT_DIR, species_name_)
        # Write ids to a GI list file (analogous to downlaoding the pubmed results as a file with the GI list format)
        with open(FILENAME, 'w') as file:
            print('Writing results to {0}'.format(FILENAME))
            for line in r.content.decode('utf-8').splitlines():
                if line[0:5] == '\t<Id>':
                    #print(line[5:-5])
                    file.write(line[5:-5]+'\n')
        print('Waiting...')
        time.sleep(1)

        # Creating blast database aliases
        DATABASE_SPECIES = '{0}/../blast_database/nr_{1}'.format(APPLICATION_PATH, species_name_)
        DATABASE_TITLE = 'nr_{0}'.format(species_name_)
        comment = 'Creating a database subset for {0}'.format(species_name)
        bash_command = 'blastdb_aliastool -gilist {0} -db {1} -out {2} -title {3}'.format(FILENAME, DATABASE, DATABASE_SPECIES, DATABASE_TITLE)
        execute_subprocess(comment, bash_command)
        print('\n')

if __name__ == '__main__':
    run()
