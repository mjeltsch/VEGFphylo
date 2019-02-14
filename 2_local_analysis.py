#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import argparse, subprocess, Bio, os, sys, shutil, re, time, datetime, socket, random, requests, xmltodict, json, csv, sqlite3
from Bio import SeqIO
from Bio import Phylo
from Bio import SearchIO
from Bio.Blast import NCBIXML
from os.path import basename, dirname, splitext, split
# To print some terminal output in color
import xml.etree.ElementTree as ET
from phylolib import load_blacklist, load_dictionary, insert_line_breaks, write_dict_to_file, read_file_to_dict

# Puropose: To programmatically retrieve the species numbers in the non-redundant NCBI protein Database
# using e utilities: https://www.ncbi.nlm.nih.gov/books/NBK25500/#chapter1.Searching_a_Database
#
# Database:
# https://www.ncbi.nlm.nih.gov/taxonomy/
# Sample query via the web form: porifera [subtree] AND specified [prop] NOT subspecies [rank]

# Parse command line
parser = argparse.ArgumentParser()
parser.add_argument("directory", help = "Specify the subdirectory where the data is!", default = '')
args = parser.parse_args()

def get_species_number(taxon):
    print('Retrieving species number for {0}... -> '.format(taxon), end='')
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[subtree]+AND+specified[prop]+NOT+subspecies[rank]'.format(taxon)
    r = requests.get(URL)
    text = r.text
    #print(text)
    return int(text.split("Count>")[1][:-2])

def get_taxon_id_and_phylum(species_name):
    # Check whether the species is in the local sqlite database
    try:
        print('Getting id for \'{0}\', checking local sqlite database first.'.format(species_name), end = '')
        data = db_retrieve_species(conn, species_name, True)
        if data != 0:
            TAXON_ID = data[1]
            phylum = data[2]
    # If the species is not in the local dictionary, get the data from NCBI
    except Exception as ex:
        print("Species not in local sqlite database. Getting id from NCBI...")
        while True:
            try:
                print('Retrieving taxon id and phylum for {0} '.format(species_name), end='')
                URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[Scientific Name]'.format(species_name)
                r = requests.get(URL, timeout=20)
                # Extract taxon_id
                list = r.text.split("Id>")
                TAXON_ID = list[1][:-2]
                # Sleep for online requests (otherwise you'll be blocked after a while!)
                time.sleep(2)
            except Exception as ex:
                print('Unexpected server response:\n{0}\n\nError: {1}\n\nTrying again...\n'.format(text, str(ex)))
                # Wait one minute if server returns an error (mostly because you are blocked or the network is down)
                time.sleep(60)
                continue
            # If there were no exceptions (= we got the taxon_id), break out of the while loop
            else:
                break
        print('=> {0}'.format(TAXON_ID))
        # Identify to which phylum the species belongs. For this, we need to download the full taxonomy tree for the species (which
        # is impossible via the API). Hence, we retrieve the full web page for the taxon and look whether any of the phyla from our
        # TAXON_DICTIONARY_FILE are part of the full taxonomy tree.
        while True:
            try:
                URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={0}&retmode=xml&rettype=full'.format(TAXON_ID)
                r = requests.get(URL)
                text = r.text
                print(text)
                i = 0
                for taxon, taxon_data in taxon_dictionary.items():
                    i += 1
                    print('Looking for TAXON_ID {0}'.format(taxon_data[0]))
                    if '<TaxId>{0}</TaxId>'.format(taxon_data[0]) in text:
                        print('Adding phylum info "{0}" for {1} to local dictionary'.format(taxon, species_name))
                        phylum = taxon
                phylum
                # The following line will be only executed if calling 'phylum' does not give a NameError
                print('{0} phyla (of 51) checked before hit was found.'.format(i))
            except NameError as err:
                print("Sleeping due to server error. Will retry in 30 seconds.")
                time.sleep(60)
                continue
            # If there were no exceptions (= we got the phylum), break out of the while loop
            else:
                if db_insert_species(conn, species_name, TAXON_ID, phylum, False) != 0:
                    print('New species {0} inserted into database.'.format(species_name))
                break
    return TAXON_ID, phylum

# How many sequences are in the local database for this taxon? The gi files were manually downloaded from
# NCBI. There should be smarted way to get these!
def get_sequence_number(taxon):
    print('Retrieving number of sequences for taxon {0}... -> '.format(taxon), end='')
    SEQUENCE_GI_LIST = '{0}/data/gi_lists/{1}.gi'.format(APPLICATION_PATH, taxon)
    with open(SEQUENCE_GI_LIST) as file:
        for i, l in enumerate(file):
            pass
    return i + 1

def get_protein_data(taxon):
    if args.directory == '':
        return []
    else:
        new_protein_data = {}
        print('\n---------------')
        print('Analysing TAXON {0}'.format(taxon))
        print('---------------\n')
        for protein, protein_data in protein_dictionary.items():
            # Here we analyze the blasts hits for each TAXON/PROTEIN pair:
            # a) number of blast hits
            # b) number of false positive according to the description (categorized in a list)
            #
            # Open the individual results file for a protein and a phylum
            XML_RESULTS_FILE = '{0}{1}/blast_results_{2}.xml'.format(args.directory, protein, taxon)
            # Create database connection
            conn = create_connection(DATABASE_FILE)

            # Get extended information about the false-positive hits
            negative_dict = {}
            for related_protein in protein_data[10]:
                negative_dict[related_protein] = 0
            # Get extended information about the true-positive hits
            positive_dict = {}
            for synonym in protein_data[9]:
                positive_dict[synonym] = 0

            print('Analysing {0}:'.format(protein))
            blastp_result = SearchIO.read(XML_RESULTS_FILE, "blast-xml")
            number = len(blastp_result)
            for hit in blastp_result:
                # get gid and protein_id
                hitident = hit.id.split('|')
                print('hit: {0} - {1}'.format(hitident[1], hitident[3]))
                # Extract species name
                species = re.search(r'\[(.*?)\]',hit.description).group(1)
                # BLASTHIT = [gid, protein_id, species, fasta_description]
                BLASTHIT = [hitident[1], hitident[3], species, hit.description]
                print('Checking false- or true-positivity.')
                for related_protein in protein_data[10]:
                    if related_protein in str(hit.description):
                        negative_dict[related_protein] += 1
                        print('Potential false-positive found: {0}'.format(related_protein))
                for synonym in protein_data[9]:
                    if synonym in str(hit.description):
                        positive_dict[synonym] += 1
                        print('True positive found: {0}'.format(synonym))
                with conn:
                    # If the gid is already in the database, the insertion fails
                    print(db_insert_protein(conn, BLASTHIT, False))
            print('Analysis for {0}/{1} completed.'.format(taxon, protein))
            print('Number of hits: {0}'.format(number))
            print('Number of true-/false-positives: {0}/{1}\n'.format(negative_dict, positive_dict))
            new_protein_data[protein] = [number, negative_dict, positive_dict]
        print('Analysis for taxon {0} completed.'.format(taxon))
        print('new_protein_data: {0}'.format(new_protein_data))
        # Returns a dictionary with related proteins and the number of hits for them (according to fasta description)
        return new_protein_data

def get_fully_sequenced_genomes(CSV_FILE):
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
    # Initialize dictionary with zeros
    full_genome_dictionary = {}
    for taxon, taxon_data in taxon_dictionary.items():
        full_genome_dictionary[taxon] = 0
    for line in input_csv_file:
        species_name = line['#Organism Name']
        print('species_name = {0}'.format(species_name))
        TAXON_ID, phylum = get_taxon_id_and_phylum(species_name)
        full_genome_dictionary[phylum] += 1
        print('Adding +1 to fully sequenced {0}'.format(phylum))
    print(str(full_genome_dictionary))
    return full_genome_dictionary

def create_connection(DATABASE_FILE):
    global conn
    try:
        conn = sqlite3.connect(DATABASE_FILE)
        return conn
    except Error as err:
        print(err)
    return None

def db_insert_protein(CONNECTION, BLASTHIT, VERBOSE):
    if VERBOSE: print('Trying to execute SQL insertion ({0})...'.format(BLASTHIT[1]))
    # In some fasta descriptions one can really find this char: ' !!!!
    BLASTHIT[2] = BLASTHIT[2].replace('\'', '')
    query = 'INSERT INTO protein (gid, protein_id, fasta_description, species) VALUES ({0}, \'{1}\', \'{2}\', \'{3}\')'.format(BLASTHIT[0], BLASTHIT[1], BLASTHIT[2], BLASTHIT[3])
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.lastrowid
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    return result

# Give a list of scientific_name, taxon_id and phylum to insert into SQLite database
# There should never be the need to update any entry from this table, only to insert
def db_insert_species(CONNECTION, SPECIES_NAME, TAXON_ID, PHYLUM, VERBOSE):
    if VERBOSE: print('Trying to execute SQL insertion ({0})...'.format(SPECIES_NAME))
    query = 'INSERT INTO species (scientific_name, taxon_id, phylum) VALUES (\'{0}\', {1}, \'{2}\')'.format(SPECIES_NAME, TAXON_ID, PHYLUM)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.lastrowid
        if VERBOSE: print('SQL insertion successful.')
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    return result

# Return a list of scientific_name, taxon_id and phylum when queried by scientific_name
def db_retrieve_species(CONNECTION, SPECIES_NAME, VERBOSE):
    if VERBOSE: print('Trying to retrieve SQLite data for {0}'.format(SPECIES_NAME))
    query = 'SELECT * FROM species WHERE scientific_name = \'{0}\''.format(SPECIES_NAME)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.fetchone()
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    # There should be only one row!
    return result

def run():
    global APPLICATION_PATH, taxon_dictionary, protein_dictionary, blacklist, DATABASE_FILE
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    #print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_data.py'.format(APPLICATION_PATH)
    PROTEIN_DICTIONARY_FILE = '{0}/data/master_dictionary.py'.format(APPLICATION_PATH)
    CSV_FILE = '{0}/data/genomes.csv'.format(APPLICATION_PATH)
    DATABASE_FILE = '{0}/data/database.sql'.format(APPLICATION_PATH)

    preamble1, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)
    #print('taxon_dictionary: {0}\n'.format(taxon_dictionary))
    #print("Populating new taxon dictionary")
    new_taxon_dictionary = taxon_dictionary

    preamble2, protein_dictionary = load_dictionary(PROTEIN_DICTIONARY_FILE)
    #print('protein_dictionary: {0}\n'.format(protein_dictionary))

    blacklist = load_blacklist()
    for taxon, taxon_data in taxon_dictionary.items():
        #print("Enter recursion")
        if taxon not in blacklist:
            #print("Passed blacklist loading")
            # Replace old data if new data is available
            # PROTEIN DATA (how many hits, false-positive list)
            try:
                new_protein_data = get_protein_data(taxon)
            except Exception as err:
                print('Error was: {0}'.format(err))
            else:
                new_taxon_dictionary[taxon][4] = new_protein_data
            # TAXON_ID
            # try:
            #     new_taxon_id, phylum = get_taxon_id_and_phylum(taxon)
            # except Exception as err:
            #     print('Could not get new taxon id and phylum. Error was: {0}',value(err))
            # # execute if no exceptions have occured
            # else:
            #     new_taxon_dictionary[taxon][0] = new_taxon_id
            # NUMBER OF SPECIES IN NCBI DATABASE
            # try:
            #     new_species_number = get_species_number(taxon)
            # except Exception as err:
            #     print('Could not get new species number. Error was: {0}',value(err))
            # else:
            #     new_taxon_dictionary[taxon][1] = new_species_number
            # NUMBER OF SEQUENCES IN NCBI DATABASE
            # try:
            #     new_sequence_number = get_sequence_number(taxon)
            # except Exception as err:
            #     print('Could not get new sequence number. Error was: {0}',value(err))
            # else:
            #     new_taxon_dictionary[taxon][3] = new_sequence_number
            # Wait in order not to overload the server
            time.sleep(1)
    # Add the number of fully sequenced genomes to the new taxon dictionary
    full_genome_dictionary = get_fully_sequenced_genomes(CSV_FILE)
    for taxon, number in full_genome_dictionary.items():
        new_taxon_dictionary[taxon][5] = number
    conn.commit()
    #print('\nWriting new_taxon_dictionary:\n{0}'.format(str(new_taxon_dictionary)))
    # Make backup file before overwriting
    os.rename(TAXON_DICTIONARY_FILE, TAXON_DICTIONARY_FILE+'~')
    # Write new taxon data to file
    write_dict_to_file(preamble1, new_taxon_dictionary, TAXON_DICTIONARY_FILE)
    if insert_line_breaks(TAXON_DICTIONARY_FILE) == True:
        print('Successfully formatted the taxon data file.')
    else:
        print('Formating the taxon data file failed.')

if __name__ == '__main__':
    run()
