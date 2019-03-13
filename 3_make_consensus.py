#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import os, sqlite3, ntpath
from Bio import SeqIO
from Bio import Entrez
from pathlib import Path
from phylolib import execute_subprocess

def create_connection(DATABASE_FILE):
    global conn
    try:
        conn = sqlite3.connect(DATABASE_FILE)
        return conn
    except Error as err:
        print(err)
    return None

# Return a list of phyla
def db_retrieve_phyla(CONNECTION, VERBOSE=True):
    if VERBOSE: print('Trying to retrieve all phyla from SQLite database.'.format())
    query = 'SELECT DISTINCT phylum FROM species'.format()
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.fetchall()
        #print('\nRESULT: {0}\n'.format(result))
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    else:
        phylum_list = []
        for item in result:
            #print('item: {0}'.format(item[0]))
            phylum_list.append(item[0])
        print('phylum_list:\n{0}'.format(phylum_list))
    return phylum_list

# Return a list of the ortholog groups
def db_retrieve_ortholog_groups(CONNECTION, VERBOSE=True):
    if VERBOSE: print('Trying to retrieve all ortholog groups from SQLite database.'.format())
    query = 'SELECT DISTINCT ortholog_group FROM protein WHERE ortholog_group != \'None\''.format()
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.fetchall()
        #print('\nRESULT: {0}\n'.format(result))
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    else:
        ortholog_group_list = []
        for item in result:
            #print('item: {0}'.format(item[0]))
            ortholog_group_list.append(item[0])
        if VERBOSE: print('ortholog_group_list:\n{0}'.format(ortholog_group_list))
    return ortholog_group_list

# Return a list of proteins (which should be later used to construct a consensus sequence)
def db_retrieve_proteins(CONNECTION, PHYLUM, ORTHOLOG_GROUP, VERBOSE=True):
    if VERBOSE: print('Trying to retrieve all {0} for phylum {1} from SQLite database.'.format(ORTHOLOG_GROUP, PHYLUM))
    query = 'SELECT gid, species, ortholog_group, fasta_description FROM protein JOIN species on protein.species = species.scientific_name WHERE species.phylum = \'{0}\' AND protein.ortholog_group = \'{1}\' AND protein.fasta_description NOT LIKE \'%partial%\''.format(PHYLUM, ORTHOLOG_GROUP)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.fetchall()
        #print('\nRESULT: {0}\n'.format(result))
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    else:
        protein_list = []
        for item in result:
            #print('item: {0}'.format(item[0]))
            protein_list.append(item)
        if VERBOSE: print('protein_list:\n{0}'.format(protein_list))
    return protein_list

def download_and_append_to_fasta_file(PHYLUM, ORTHOLOG, sequence_dictionary):
    # Get all fasta sequences from the SEQUENCE_LIST via Entrez
    Entrez.email = "michael@jeltsch.org"
    Entrez.tool = "local_script_under_development"
    # Make an empty file
    Path('{0}/{1}-{2}.fasta'.format(MULTISEQUENCE_FASTA_DIR, ORTHOLOG, PHYLUM)).touch()
    with open('{0}/{1}-{2}.fasta'.format(MULTISEQUENCE_FASTA_DIR, ORTHOLOG, PHYLUM),"a") as multiple_fasta_file:
        for item in sequence_dictionary:
            # If the sequence exists locally
            single_fasta_file = SEQUENCE_FILES+'{0}.fasta'.format(item[0])
            if os.path.isfile(single_fasta_file):
                print('Retrieving {0} locally.'.format(item))
                seq_record = SeqIO.read(single_fasta_file, "fasta")
            # get it remotely
            else:
                file = open(single_fasta_file,"w")
                print('Retrieving {0} from Entrez.'.format(item))
                try:
                    with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=item[0]) as handle:
                        seq_record = SeqIO.read(handle, "fasta")
                        file.write(seq_record.format("fasta"))
                except Exception as ex:
                    print("Problem contacting Blast server. Skipping " + str(item) + ". Error: " + str(ex))
                file.close()
            multiple_fasta_file.write(seq_record.format("fasta"))
    multiple_fasta_file.close()

def limit_to_most_informative_sequences(NUMBER, FASTA_FILE_IN, FASTA_FILE_OUT):
    execute_subprocess(
        'Reducing the complexity to the most informative {} sequences:'.format(NUMBER),
        't_coffee -other_pg seq_reformat -in {0} -action +trim _seq_n{1} -output fasta_seq > {2}'.format(FASTA_FILE_IN, NUMBER, FASTA_FILE_OUT))

def run():
    global SEQUENCE_FILES, APPLICATION_PATH, MULTISEQUENCE_FASTA_DIR
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    MULTISEQUENCE_FASTA_DIR = '{0}/data/all_multifasta_files/'.format(APPLICATION_PATH)
    INFORMATIVE_MULTISEQUENCE_FASTA_DIR = '{0}/data/informative_multifasta_files/'.format(APPLICATION_PATH)
    DATABASE_FILE = '{0}/data/database.sqlite3'.format(APPLICATION_PATH)
    SEQUENCE_FILES = '{0}/data/protein_sequences/'.format(APPLICATION_PATH)
    CONNECTION = create_connection(DATABASE_FILE)
    with CONNECTION:
        phylum_list = db_retrieve_phyla(CONNECTION)
        ortholog_group_list = db_retrieve_ortholog_groups(CONNECTION)
        counter_phyla = 0
        counter_orthologs = 0
        counter_proteins = 0
        for PHYLUM in phylum_list:
            counter_phyla += 1
            for ORTHOLOG in ortholog_group_list:
                counter_orthologs += 1
                protein_list = db_retrieve_proteins(CONNECTION, PHYLUM, ORTHOLOG, VERBOSE=False)
                if len(protein_list) > 0:
                    print('\n\n{0} - {1}:'.format(PHYLUM, ORTHOLOG))
                    print(protein_list)
                    counter_proteins += len(protein_list)
                    download_and_append_to_fasta_file(PHYLUM, ORTHOLOG, protein_list)
    print('\nphyla: {0}'.format(counter_phyla))
    print('ortholog groups: {0}'.format(round(counter_orthologs/counter_phyla)))
    print('proteins: {0}'.format(counter_proteins))
    # Reduce the number of sequences for each phylum/ortholog group combinbation to NUMBER
    os.chdir(MULTISEQUENCE_FASTA_DIR)
    NUMBER = 4
    for filename_in in os.listdir(MULTISEQUENCE_FASTA_DIR):
        #filename_in = '{0}{1}'.format(MULTISEQUENCE_FASTA_DIR, filename_in)
        print('filename_in: {0}'.format(filename_in))
        filename_out = '../informative_multifasta_files/{0}'.format(filename_in)
        print('filename_out: {0}'.format(filename_out))
        limit_to_most_informative_sequences(NUMBER, filename_in, filename_out)

if __name__ == '__main__':
    run()
