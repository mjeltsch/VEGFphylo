#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import os, sqlite3
from Bio import SeqIO
from Bio import Entrez

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

def download(sequence_dictionary):
    # Get all fasta sequences from the SEQUENCE_LIST via Entrez
    Entrez.email = "michael@jeltsch.org"
    Entrez.tool = "local_script_under_development"
    file = open(FASTA,"w")
    print("Retrieving sequences from Entrez:\n")
    for item in sequence_dictionary:
        # If the sequence is given locally
        if len(sequence_dictionary[item]) > 4:
            print(sequence_dictionary[item][4])
            file.write(sequence_dictionary[item][4])
        else:
            print("Retrieving " + item)
            try:
                with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=item) as handle:
                    seq_record = SeqIO.read(handle, "fasta")
                    file.write(seq_record.format("fasta"))
            except Exception as ex:
                print("Problem contacting Blast server. Skipping " + item + ". Error: " + str(ex))
    file.close()


def run():
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    DATABASE_FILE = '{0}/data/database.sqlite3'.format(APPLICATION_PATH)
    SEQUENCE_FILES = '{0}/data/alignments/'.format(APPLICATION_PATH)
    CONNECTION = create_connection(DATABASE_FILE)
    with CONNECTION:
        phylum_list = db_retrieve_phyla(CONNECTION)
        ortholog_group_list = db_retrieve_ortholog_groups(CONNECTION)
        for PHYLUM in phylum_list:
            for ORTHOLOG_GROUP in ortholog_group_list:
                protein_list = db_retrieve_proteins(CONNECTION, PHYLUM, ORTHOLOG_GROUP, VERBOSE=False)
                if len(protein_list) > 0:
                    print('\n\n{0} - {1}:'.format(PHYLUM, ORTHOLOG_GROUP))
                    print(protein_list)

if __name__ == '__main__':
    run()
