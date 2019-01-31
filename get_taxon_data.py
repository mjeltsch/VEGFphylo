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
        content = content.replace("}], '","}],\n     '")
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
            preamble = '#'
            dictionary = {}
            print('Could not read taxon dictionary {0}'.format(FILENAME))
    else:
        preamble = '#'
        dictionary = {}
    return preamble, dictionary

def get_species_number(taxon):
    print('Retrieving species number for {0}... -> '.format(taxon), end='')
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[subtree]+AND+specified[prop]+NOT+subspecies[rank]'.format(taxon)
    r = requests.get(URL)
    text = r.text
    #print(text)
    return int(text.split("Count>")[1][:-2])

def get_taxon_id_online(taxon):
    print('Retrieving taxon id for {0}... -> '.format(taxon), end='')
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[Scientific Name]'.format(taxon)
    r = requests.get(URL)
    text = r.text
    #print(text)
    list = text.split("Id>")
    return int(list[1][:-2])

def get_taxon_id(species_name):
    while True:
        # Check whether the species is in the local dictionary
        try:
            print('Getting id for \'{0}\', checking local taxon_id_dictionary first.'.format(species_name), end = '')
            #print('taxon_id_dict: {0}'.format(taxon_id_dict))
            #print('taxon_id_dict[species_name]: {0}'.format(taxon_id_dict[species_name]))
            # print(taxon_id_dict[species_name][0])
            TAXON_ID = taxon_id_dict[species_name][0]
            break
        # If the species is not in the local dictionary, get the data from NCBI
        except Exception as ex:
            print("Species not in local dictionary. Getting id from NCBI...")
            try:
                print('Retrieving taxon id for {0} '.format(species_name), end='')
                URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[Scientific Name]'.format(species_name)
                r = requests.get(URL, timeout=20)
                # Extract taxon_id
                list = r.text.split("Id>")
                TAXON_ID = list[1][:-2]
                # Add new taxon_id to dictionary
                taxon_id_dict[species_name][0] = int(TAXON_ID)
                # Sleep for online requests (otherwise you'll be blocked after a while!)
                time.sleep(2)
            except Exception as ex:
                print('Unexpected server response:\n{0}\n\nError: {1}\n\nTrying again...\n'.format(text, str(ex)))
                # Wait one minute if server returns an error (mostly because you are blocked or the network is down)
                time.sleep(60)
                continue
            # If there were no exceptions, break out of the while loop
            else:
                break
    print('=> {0}'.format(TAXON_ID))
    return TAXON_ID

def get_phylum(species_name, TAXON_ID):
    # Identify to which phylum the species belongs. For this, we need to download the full taxonomy tree for the species (which
    # is impossible via the API). Hence, we retrieve the full web page for the taxon and look whether any of the phyla from our
    # TAXON_DICTIONARY_FILE are part of the full taxonomy tree.
    try:
        phylum = taxon_id_dict[species_name][1]
        if phylum == None:
            raise Exception('Phylum for {0} was None'.format(species_name))
    except Exception as err:
        print('Phylum info for {0} not available from local dictionary. Getting from NCBI...'.format(species_name))
        URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={0}&retmode=xml&rettype=full'.format(TAXON_ID)
        while True:
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
                    taxon_id_dict[species_name][1] = phylum
            try:
                phylum
                # The following line will be only executed if calling 'phylum' does not give a NameError
                print('{0} phyla (of 51) checked before hit was found.'.format(i))
                # break applies to both while and for and therefore cannot be inside the for loop
                break
            except NameError as err:
                print("Sleeping due to server error. Will retry in 30 seconds.")
                time.sleep(30)
                continue
    return phylum

def get_sequence_number_local(taxon):
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
            with open(XML_RESULTS_FILE) as blastresults:
                print('Parsing {0}'.format(XML_RESULTS_FILE))
                blast_records = NCBIXML.parse(blastresults)
                blast_records = list(blast_records)
                print('Analysing {0}:'.format(protein))
                negative_dict = {}
                # Get extended information about the matching sequences
                for related_protein in protein_data[10]:
                    i = 0
                    for blast_record in blast_records:
                        number = len(blast_record.descriptions)
                        print('Checking for group of false positives \"{0}\"'.format(related_protein))
                        for description in blast_record.descriptions:
                            #print('description: {0}'.format(description))
                            if related_protein in str(description):
                                i += 1
                                print('Potential false-positive found: {0}'.format(related_protein))
                    negative_dict[related_protein] = i
                print('Analysis for {0}/{1} completed.'.format(taxon, protein))
                print('Number of hits: {0}'.format(number))
                print('Number of false-positives: {0}\n'.format(negative_dict))
            new_protein_data[protein] = [number, negative_dict]
        print('Analysis for taxon {0} completed.'.format(taxon))
        print('new_protein_data: {0}'.format(new_protein_data))
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
        TAXON_ID = get_taxon_id(species_name)
        phylum = get_phylum(species_name, TAXON_ID)
        full_genome_dictionary[phylum] += 1
        print('Adding +1 to fully sequenced {0}'.format(phylum))
    print(str(full_genome_dictionary))
    return full_genome_dictionary

def run():
    global APPLICATION_PATH, taxon_dictionary, TAXON_ID_DICTIONARY_FILE, protein_dictionary, taxon_id_dict, blacklist
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    #print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_data.py'.format(APPLICATION_PATH)
    TAXON_ID_DICTIONARY_FILE = '{0}/data/taxon_ids.py'.format(APPLICATION_PATH)
    PROTEIN_DICTIONARY_FILE = '{0}/data/master_dictionary.py'.format(APPLICATION_PATH)
    CSV_FILE = '{0}/data/genomes.csv'.format(APPLICATION_PATH)

    preamble1, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)
    #print('taxon_dictionary: {0}\n'.format(taxon_dictionary))
    #print("Populating new taxon dictionary")
    new_taxon_dictionary = taxon_dictionary

    preamble2, protein_dictionary = load_dictionary(PROTEIN_DICTIONARY_FILE)
    #print('protein_dictionary: {0}\n'.format(protein_dictionary))

    preamble3, taxon_id_dict = load_dictionary(TAXON_ID_DICTIONARY_FILE)
    #print('taxon_id_dict: {0}\n'.format(taxon_id_dict))

    blacklist = load_blacklist()
    for taxon, taxon_data in taxon_dictionary.items():
        #print("Enter recursion")
        if taxon not in blacklist:
            #print("Passed blacklist loading")
            # Replace old data if new data is available
            # PROTEIN DATA (how many hits, false-positive list)
            try:
                new_protein_data = get_protein_data(taxon)
            except:
                pass
            else:
                new_taxon_dictionary[taxon][4] = new_protein_data
            # TAXON_ID
            try:
                new_taxon_id
            except NameError:
                pass
            else:
                new_taxon_dictionary[taxon][0] = new_taxon_id
            # NUMBER OF SPECIES IN NCBI DATABASE
            try:
                new_species_number
            except NameError:
                pass
            else:
                new_taxon_dictionary[taxon][1] = new_species_number
            # NUMBER OF SEQUENCES IN NCBI DATABASE
            try:
                new_sequence_number
            except NameError:
                pass
            else:
                new_taxon_dictionary[taxon][3] = new_sequence_number
            # Wait in order not to overload the server
            time.sleep(1)
    # Add the number of fully sequenced genomes to the new taxon dictionary
    full_genome_dictionary = get_fully_sequenced_genomes(CSV_FILE)
    for taxon, number in full_genome_dictionary.items():
        new_taxon_dictionary[taxon][5] = number
    #print('\nWriting new_taxon_dictionary:\n{0}'.format(str(new_taxon_dictionary)))
    # Make backup file before overwriting
    os.rename(TAXON_DICTIONARY_FILE, TAXON_DICTIONARY_FILE+'~')
    os.rename(TAXON_ID_DICTIONARY_FILE, TAXON_ID_DICTIONARY_FILE+'~')
    # Write new taxon data to file
    write_dict_to_file(preamble1, new_taxon_dictionary, TAXON_DICTIONARY_FILE)
    write_dict_to_file(preamble3, taxon_id_dict, TAXON_ID_DICTIONARY_FILE)
    if insert_line_breaks(TAXON_DICTIONARY_FILE) == True:
        print('Successfully formatted the taxon data file.')
    else:
        print('Formating the taxon data file failed.')

if __name__ == '__main__':
    run()
