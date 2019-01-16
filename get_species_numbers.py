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
    print("Enter subroutine")
    # Load a sequence_dictionary if it exists
    if os.path.isfile(FILENAME):
        try:
            preamble, dictionary = read_file_to_dict(FILENAME)
            print("\nReading in the dictionary " + FILENAME + ":\n")
            for key, value in dictionary.items():
                print(key, value)
        except Exception as e:
            dictionary = {}
            print('Could not read taxon dictionary {0}'.format(FILENAME))
    else:
        dictionary = {}
    return preamble, dictionary

def run():
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_ids.py'.format(APPLICATION_PATH)
    preamble, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)
    new_taxon_dictionary = {}
    for taxon, taxon_data in taxon_dictionary.items():
        print('Retrieving species number for {0}... -> '.format(taxon), end='')
        URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[subtree]+AND+specified[prop]+NOT+subspecies[rank]'.format(taxon)
        r = requests.get(URL)
        text = r.text
        print(text)
        new_species_number = int(text.split("Count>")[1][:-2])
        print(str(new_species_number))
        print('Retrieving taxon id for {0}... -> '.format(taxon), end='')
        URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[Scientific Name]'.format(taxon)
        r = requests.get(URL)
        text = r.text
        #print(text)
        list = text.split("Id>")
        new_taxon_id = int(list[1][:-2])
        print(str(new_taxon_id))
        new_taxon_dictionary[taxon] = [new_taxon_id, new_species_number, taxon_data[2]]
        time.sleep(1)
    print(new_taxon_dictionary)
    write_dict_to_file(preamble, new_taxon_dictionary, TAXON_DICTIONARY_FILE)

if __name__ == '__main__':
    run()
