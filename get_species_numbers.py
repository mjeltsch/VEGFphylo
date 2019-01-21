#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import argparse, subprocess, Bio, os, sys, shutil, re, time, datetime, socket, random, requests, xmltodict, json
from Bio import Entrez
from Bio import SeqIO
from Bio import Phylo
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, SequenceFace, ImgFace, SVGFace, faces, add_face_to_node
from os.path import basename, dirname, splitext, split
# To print some terminal output in color
from termcolor import colored, cprint
import xml.etree.ElementTree as ET

# Puropose: To programmatically retrieve the species numbers in the non-redundant NCBI protein Database
# using e utilities: https://www.ncbi.nlm.nih.gov/books/NBK25500/#chapter1.Searching_a_Database
#
# Database:
# https://www.ncbi.nlm.nih.gov/taxonomy/
# Sample query via the web form: porifera [subtree] AND specified [prop] NOT subspecies [rank]

def read_file_to_dict(file_name):
    try:
        with open(file_name, "r") as file:
            string = file.read()
            file.close()
            #print('File content:\n' + string)
            dictionary = eval(string)
            # Store the first three lines of the dictionary
            buf = string.split('\n')
            preamble = ''
            for i in range(0,len(buf)):
                if buf[i][:1] == '#':
                    preamble += buf[i] + '\n'
                else:
                    break
            return preamble, dictionary
    except Exception as ex:
        print('Could not read dictionary from file {0}. Error: '.format(file_name) + str(ex))
        return '', {}

def write_dict_to_file(preamble, dictionary, file_name):
    try:
        with open(file_name, "w") as file:
            file.write(preamble)
            file.write(str(dictionary))
            file.close()
            return True
    except Exception as ex:
        print("Could not write dictionary to file. Error: " + str(ex))
        return False

def load_dictionary(FILENAME):
    #print("Enter subroutine")
    # Load a sequence_dictionary if it exists
    if os.path.isfile(FILENAME):
        try:
            preamble, dictionary = read_file_to_dict(FILENAME)
            print("\nReading in the dictionary " + FILENAME + ":\n")
            for key, value in dictionary.items():
                print(key, value)
            print("Done reading dictionary.")
        except Exception as e:
            dictionary = {}
            print('Could not read taxon dictionary {0}'.format(FILENAME))
    else:
        dictionary = {}
    return preamble, dictionary

def get_species_number(taxon):
    print('Retrieving species number for {0}... -> '.format(taxon), end='')
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[subtree]+AND+specified[prop]+NOT+subspecies[rank]'.format(taxon)
    r = requests.get(URL)
    text = r.text
    #print(text)
    return int(text.split("Count>")[1][:-2])

def get_taxon_id(taxon):
    print('Retrieving taxon id for {0}... -> '.format(taxon), end='')
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[Scientific Name]'.format(taxon)
    r = requests.get(URL)
    text = r.text
    #print(text)
    list = text.split("Id>")
    return int(list[1][:-2])

def get_sequence_number_local(taxon):
    print('Retrieving number of sequences for taxon {0}... -> '.format(taxon), end='')
    SEQUENCE_GI_LIST = '{0}/data/gi_lists/{1}.gi'.format(APPLICATION_PATH, taxon)
    with open(SEQUENCE_GI_LIST) as file:
        for i, l in enumerate(file):
            pass
    return i + 1

def get_protein_data(protein, taxon):
    return [1, 2, 4, 2222]

def run():
    global APPLICATION_PATH
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_ids.py'.format(APPLICATION_PATH)
    PROTEIN_DICTIONARY_FILE = '{0}/data/master_dictionary.py'.format(APPLICATION_PATH)
    preamble, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)
    preamble, protein_dictionary = load_dictionary(PROTEIN_DICTIONARY_FILE)
    print("Populating new taxon dictionary")
    new_taxon_dictionary = taxon_dictionary
    print('new_taxon_dictionary: {0}'.format(str(new_taxon_dictionary)))
    for taxon, taxon_data in taxon_dictionary.items():
        #print("Enter recursion")
        #new_species_number = get_species_number(taxon)
        #print(str(new_species_number))
        #new_taxon_id = get_taxon_id(taxon)
        #print(str(new_taxon_id))
        #new_sequence_number = get_sequence_number_local(taxon)
        #print(str(new_sequence_number))
        new_protein_data = {}
        for protein, proteindata in protein_dictionary.items():
            new_protein_data[protein] = get_protein_data(protein, taxon)
            print('new_protein_data for taxon {0}: {1} => {2}'.format(taxon, protein, new_protein_data[protein]))
        # Replace old data if new data is available
        try:
            new_taxon_id
        except NameError:
            pass
        else:
            new_taxon_dictionary[taxon][0] = new_taxon_id
        try:
            new_species_number
        except NameError:
            pass
        else:
            new_taxon_dictionary[taxon][1] = new_species_number
        try:
            new_sequence_number
        except NameError:
            pass
        else:
            new_taxon_dictionary[taxon][3] = new_sequence_number
        try:
            new_protein_data
        except NameError:
            pass
        else:
            print('taxon: {0}'.format(taxon))
            new_taxon_dictionary[taxon][4] = new_protein_data
        # Wait in order not to overload the server
        time.sleep(1)
    print("\Writing new_taxon_dictionary:\n" + str(new_taxon_dictionary))
    # Make backup file before overwriting
    os.rename(TAXON_DICTIONARY_FILE, TAXON_DICTIONARY_FILE+'~')
    # Write new data to file
    write_dict_to_file(preamble, new_taxon_dictionary, TAXON_DICTIONARY_FILE)

if __name__ == '__main__':
    run()
