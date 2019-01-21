#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import argparse, subprocess, Bio, os, sys, shutil, re, time, datetime, socket, random, requests, xmltodict, json, csv
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
from phylolib import load_blacklist

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

def insert_line_breaks(file_name):
    try:
        with open(file_name, "r") as file:
            content = file.read()
            file.close()
        content = content.replace("], '","],\n '")
        with open(file_name, "w") as file:
            file.write(content)
            file.close()
            return True
    except Exception as ex:
        print('Could not insert line breaks into file {0} Error: {1}'.format(filename, str(ex)))
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
    if args.directory == '':
        return []
    else:
        XML_RESULTS_FILE = '{0}/{1}/blast_results_{2}.xml'.format(args.directory, protein, taxon)
        with open(XML_RESULTS_FILE) as blastresults:
            print("Parsing...")
            blast_records = NCBIXML.parse(blastresults)
            blast_records = list(blast_records)
            for blast_record in blast_records:
                number = len(blast_record.descriptions)
            print(" completed.\n")
        print('number of hits: {0}'.format(str(number)))
        return number

def get_fully_sequenced_genomes(CSV_FILE):

    # Download list of all fully sequences animal genomes:
    # https://www.ncbi.nlm.nih.gov/genomes/solr2txt.cgi?q=%5Bdisplay()%5D.from(GenomeBrowser).usingschema(%2Fschema%2FGenomeAssemblies).matching(group%3D%3D%5B%22Animals%22%5D)&fields=organism%7COrganism%20Name%2Clineage%7COrganism%20Groups%2Csize%7CSize(Mb)%2Cchromosomes%7CChromosomes%2Corganelles%7COrganelles%2Cplasmids%7CPlasmids%2Cassemblies%7CAssemblies&filename=genomes.csv&nolimit=on
    full_genome_dictionary = {}
    for taxon, taxon_data in taxon_dictionary.items():
        full_genome_dictionary[taxon] = 0
    with open(CSV_FILE, newline='') as csvfile:
        # csv.DictReader stopped working with python 3.6 since it switched from reading into a dictionary
        # to reading into an ordered dictionary (which is not accepted as input for "eval")
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        # Skip the first (header) line of the CSV file
        next(reader, None)
        for line in reader:
            #print('line: {0}'.format(line))
            # str removed from eval(line)
            list = eval(str(line))
            #print('list: {0}'.format(list))
            species_name = list[0]
            print('Retrieving taxon id for {0}... -> '.format(species_name), end='')
            URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[Scientific Name]'.format(species_name)
            r = requests.get(URL)
            text = r.text
            list = text.split("Id>")
            TAXON_ID = list[1][:-2]
            URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={0}&retmode=xml&rettype=full'.format(TAXON_ID)
            r = requests.get(URL)
            text = r.text
            #print(text)
            for taxon, taxon_data in taxon_dictionary.items():
                if '<TaxId>{0}</TaxId>'.format(taxon_data[0]) in text:
                    print('Adding {0} to {1}'.format(species_name, taxon))
                    full_genome_dictionary[taxon] += 1
            time.sleep(1)
    print(str(full_genome_dictionary))
    return full_genome_dictionary

def run():
    global APPLICATION_PATH, taxon_dictionary
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_data.py'.format(APPLICATION_PATH)
    PROTEIN_DICTIONARY_FILE = '{0}/data/master_dictionary.py'.format(APPLICATION_PATH)
    CSV_FILE = '{0}/data/genomes.csv'.format(APPLICATION_PATH)
    preamble1, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)
    preamble2, protein_dictionary = load_dictionary(PROTEIN_DICTIONARY_FILE)
    print("Populating new taxon dictionary")
    new_taxon_dictionary = taxon_dictionary
    print('new_taxon_dictionary: {0}'.format(str(new_taxon_dictionary)))
    blacklist = load_blacklist()
    for taxon, taxon_data in taxon_dictionary.items():
        if taxon not in blacklist:
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
            try:
                new_protein_data
            except NameError:
                pass
            else:
                print('taxon: {0}'.format(taxon))
                new_taxon_dictionary[taxon][4] = new_protein_data
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
        # Wait in order not to overload the server
        time.sleep(1)
    # Add the number of fully sequenced genomes to the new taxon dictionary
    print('Getting list of fully sequenced genomes...')
    full_genome_dictionary = get_fully_sequenced_genomes(CSV_FILE)
    for taxon, number in full_genome_dictionary.items():
        new_taxon_dictionary[taxon][5] = number
    print('\nWriting new_taxon_dictionary:\n{0}'.format(str(new_taxon_dictionary)))
    # Make backup file before overwriting
    os.rename(TAXON_DICTIONARY_FILE, TAXON_DICTIONARY_FILE+'~')
    # Write new data to file
    write_dict_to_file(preamble1, new_taxon_dictionary, TAXON_DICTIONARY_FILE)
    if insert_line_breaks(TAXON_DICTIONARY_FILE) == True:
        print('Successfully formatted the taxon data file.')
    else:
        print('Formating the taxon data file failed.')

if __name__ == '__main__':
    run()
