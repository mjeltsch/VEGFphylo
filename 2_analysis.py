#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import argparse, subprocess, Bio, os, sys, shutil, re, time, datetime, socket, random, requests, xmltodict, json, csv, sqlite3, types
from Bio import SeqIO
from Bio import Phylo
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from os.path import basename, dirname, splitext, split
# To print some terminal output in color
import xml.etree.ElementTree as ET
from phylolib import load_blacklist, load_dictionary, insert_line_breaks, write_dict_to_file, read_file_to_dict, execution_time_str, make_synonym_dictionary, create_sqlite_file, expand_complex_taxa, download_proteins, execute_subprocess
from ete3 import NCBITaxa

# Puropose: To programmatically retrieve the species numbers in the non-redundant NCBI protein Database
# using e utilities: https://www.ncbi.nlm.nih.gov/books/NBK25500/#chapter1.Searching_a_Database
#
# Database:
# https://www.ncbi.nlm.nih.gov/taxonomy/
# Sample query via the web form: porifera [subtree] AND specified [prop] NOT subspecies [rank]


# This script goes through all primary blast result files and does the following:
#
# It looks at the protein description and determines which ortholog group it belongs to.
# If this cannot be figured out from the curated protein description using the list of
# aliases from master_dictionary.py (e.g. "HYPOTHETICAL PROTEIN"), then it initiates a
# blastp search ("backcheck blast") with that protein and categorizes the results
# according to the majority description of the 50 best blastp results (THRESHOLD = 0.5)
# into one of three groups: "synonym", "related protein" or "unknown". The unknown
# proteins are written to the file "check_manually.html".
#
# - Determines the number of species for each taxon in the NCBI database
# - Determines the number of fully sequenced genomes of the taxons specified in the taxon_dictionary.py file
# - Updates the taxon_dictionary.py file with the results


# The results of all tasks that require to connect to a remote service or database (NCBI, Entrez) already
# cached in order for subsequent runs to be faster. Only after deleting the local cache (manually by deleting
# the files), new remote requests are initiated.

class category_found(Exception):
    pass

class species_number_error(Exception):
    pass


# Parse command line
parser = argparse.ArgumentParser()
parser.add_argument("--directory", nargs='?', help = "Specify the subdirectory where the data is!", default = 'data/primary_blast_results/')
args = parser.parse_args()

# This will conenct to NCBI taxonomy database and retrieve the number of species for a given phylum
# If the NCBI database is down, this will loop forever!
def get_species_number_from_ncbi(taxon, taxon_data):
    print('Retrieving species number for {0}... -> '.format(taxon), end='')
    # Determine, whether a taxon is to be excluded from the request
    TAXA_LIST = expand_complex_taxa(taxon_data[0])
    while True:
        try:
            number = 0
            for item in TAXA_LIST:
                sign = item[0]
                URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=txid{0}[subtree]+AND+specified[prop]+NOT+subspecies[rank]'.format(item[1:])
                print('URL: {0}'.format(URL))
                r = requests.get(URL)
                text = r.text
                print('species_number_text:\n{0}'.format(text))
                if "<ERROR>" in text:
                    raise species_number_error
                # Add only the first number and substract all following numbers
                if sign == '+':
                    number += int(text.split("Count>")[1][:-2])
                if sign == '-':
                    number -= int(text.split("Count>")[1][:-2])
        except species_number_error:
            print('Remote database error. Waiting 30 seconds and repeating request.')
            time.sleep(30)
        except Exception as err:
            print('Failed to fetch species number. Error: {0}'.format(err))
        else:
            # the "else" is only executed if everything under the "try" has completed without error
            break
    return number

def get_taxon_id_and_phylum(species_name):
    # Check whether the species is in the local sqlite database
    print('Getting id & phylum for \'{0}\', checking local sqlite database first.\n'.format(species_name), end = '')
    if species_name in sqlite_species_dict:
        TAXON_ID = sqlite_species_dict[species_name][1]
        PHYLUM = sqlite_species_dict[species_name][2]
    else:
        print("Species not in local sqlite database. Getting data from NCBI... ")
        TAXON_ID = get_taxon_id_from_NCBI(species_name)
        print('=> {0}'.format(TAXON_ID))
        # Identify to which phylum the species belongs. For this, we need to download the full taxonomy tree for the species (which
        # is impossible via the API). Hence, we retrieve the full web page for the taxon and look whether any of the phyla from our
        # TAXON_DICTIONARY_FILE are part of the full taxonomy tree.
        PHYLUM = get_phylum_from_NCBI(TAXON_ID)
        print('=> {0}'.format(PHYLUM))
        if db_insert_species(species_name, TAXON_ID, PHYLUM) != 0:
            print('New species {0} inserted into database.'.format(species_name))
    return TAXON_ID, PHYLUM

# How many sequences are in the local database for this taxon? The gi files were manually downloaded from
# NCBI. E.g. https://www.ncbi.nlm.nih.gov/protein/?term=txid7777%5Borganism%5D and then "Send to file" -> "Format: GI List" -> "Create File"
# There should be smarter way to get these!
def get_sequence_number(taxon, taxon_data):
    print('Retrieving number of sequences for taxon {0}... -> '.format(taxon), end='')
    TAXA_LIST = expand_complex_taxa(taxon_data[0])
    i = 0
    for item in TAXA_LIST:
        SEQUENCE_GI_LIST = '{0}/data/gi_lists/txid{1}.gi'.format(APPLICATION_PATH, item[1:])
        with open(SEQUENCE_GI_LIST) as file:
            for j, m in enumerate(file):
                pass
        if item[0] == '+':
            i += j
        elif item[0] == '-':
            i -= j
    return i + 1

# This should be replaced by the corresponding function from ETE3
def get_taxon_id_from_NCBI(species_name, VERBOSE=True):
    ncbi = NCBITaxa()
    taxon_id_list = ncbi.get_name_translator([species_name])
    taxon_id = taxon_id_list[species_name][0]
    # Check that something sensinble was returned
    if isinstance(taxon_id, int) and taxon_id > 0:
        print('taxon_id for {0}: {1}'.format(species_name, taxon_id))
    else:
        taxon_id = None
        print('No taxon_id for {0}. NCBI returned: {1}'.format(species_name, taxon_id_list))
    #return TAXON_ID
    return taxon_id

# This sends a species id to NCBI and returns to which of the phyla within taxon_dictionary.py
# this species belongs to.
def get_phylum_from_NCBI(TAXON_ID, VERBOSE=True):
    # Counter to increase waiting period upon incresing number of unsuccesful trials
    i = 0
    while True:
        try:
            URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={0}&retmode=xml&rettype=full'.format(TAXON_ID)
            r = requests.get(URL)
            text = r.text
            if VERBOSE: print(text)
            i = 0
            for taxon, taxon_data in taxon_dictionary.items():
                # Determine, whether a taxon is complex (with some branches to exclude)
                try:
                    TAXON = taxon_data[0].split('-')[0]
                    SUBTRACT_TAXON = taxon_data[0].split('-')[1]
                except Exception as err:
                    SUBTRACT_TAXON = '-'
                if VERBOSE: print('Looking for TAXON_ID {0}, exluding TAXON_ID {1}'.format(TAXON, SUBTRACT_TAXON))
                if '<TaxId>{0}</TaxId>'.format(TAXON) in text and '<TaxId>{0}</TaxId>'.format(SUBTRACT_TAXON) not in text:
                    if VERBOSE: print('Adding phylum info {0} to species {1}'.format(taxon, TAXON_ID))
                    phylum = taxon
                    print(phylum)
                    break
                i += 1
            # The following line will be only executed if calling 'phylum' does not give a NameError
            # When i = 51, none of the 51 phyla was in the text that was received from NCBI.
            if i == 51:
                if VERBOSE: print('Species {0} does not belong to any phylum in the phylum list'.format(TAXON_ID))
                phylum = 'unknown'
            elif i < 51:
                i += 1
                if VERBOSE: print('{0} phyla (of 51) checked before hit was found.'.format(i))
            else:
                if VERBOSE: print('Unknown error when trying to identify phylum for species {0}.'.format(TAXON_ID))
                phylum = 'unknown'
        except NameError as err:
            if VERBOSE: print('Sleeping due to server error. Will retry in 60 seconds. Error: {0}'.format(err))
            time.sleep(60)
            continue
        # If there were no exceptions (= we got the phylum), break out of the while loop
        else:
            break
    return phylum

def write_to_html_scrutinize_file(protein, taxon, list_to_scrutinize):
    global TOTAL_UNKNOWN_COUNT_HTML_CHECK
    # This sorts the list of lists according to the third element of each list (= type of hit; synonym, related protein, unknown)
    list_to_scrutinize.sort(key=lambda x: x[2])
    preamble = '''<html>
<head>
<title>Results</title>
<style>
body { font-family: "Open Sans", Arial; }
</style>
<script type="text/javascript">
<!--
    function toggle_visibility(id) {
       var e = document.getElementById(id);
       if(e.style.display == 'block')
          e.style.display = 'none';
       else
          e.style.display = 'block';
    }
//-->
</script>
</head>
<body>'''
    #for item in sorted_list:
    new_type = ''
    i = 1
    j = 1
    html_text = ''
    TOTAL_UNKNOWN_COUNT_HTML_CHECK += len(list_to_scrutinize)
    print('list_to_scrutinize:\n{0}'.format(list_to_scrutinize))
    for item in list_to_scrutinize:
        # Make a heading row if the type of hit changes and color unknown stuff red
        if item[2] != new_type:
            if item[2] in ['hypothetical protein', 'uncharacterized protein', 'unknown']:
                red_color = ' style="color:Red;"'
            else:
                red_color = ''
            if item[2] != '':
                # End of toggle visibility section
                html_text += '</table>\n</div>\n\n'
            html_text += '<p{0}><a href="#"{0} onclick="toggle_visibility(\'{1}-{2}-{3}\');">{4}. {1}-{2}-{3}</a></p>\n'.format(red_color, item[0], item[1], item[2], i)
            # Start of toggle visibility section
            html_text += '<div id="{0}-{1}-{2}" style="display:none">\n'.format(item[0], item[1], item[2])
            # New table
            #lineitem += '<table border="1">\n<tr><td><h5>{0}</h5></td><td><h5>{1}</h5></td><td colspan="5"><h2>{2}</h2></td></tr>\n'.format(item[0], item[1], item[2])
            html_text += '<table border="1">\n'.format()
        link_to_protein = '<a href="https://www.ncbi.nlm.nih.gov/protein/{0}/" target="_blank">{1}</a>'.format(item[3], item[5])
        link_to_blast = '<a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?LAYOUT=OneWindow&PROGRAM=blastp&PAGE=Proteins&CMD=Web&DATABASE=nr&FORMAT_TYPE=HTML&NCBI_GI=on&SHOW_OVERVIEW=yes&QUERY={0}" target="_blank">blastp</a>'.format(item[3])
        link_to_local_msa = '<a href="../../protein_results/{0}.html#{1}" target="_blank">local MSA</a>'.format(taxon, item[3])
        #                                                                                                                                 taxon    protein  type_of  id       accesion_no
        #                                                                                                                                                   _hit
        html_text += '<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td><td>{4}</td><td>{5}</td><td>{6}</td><td>{7}</td></tr>\n'.format(item[0], item[1], item[2], item[3], item[4], link_to_protein, link_to_blast, link_to_local_msa)
        new_type = item[2]
        i += 1
    #lineitem += '\n</body>\n</html>\n'
    HTML_DIR = '{0}/data/analysis_results/{1}'.format(APPLICATION_PATH, protein)
    HTML_FILE = '{0}/{1}.html'.format(HTML_DIR, taxon)
    if not os.path.isdir(HTML_DIR):
        os.mkdir(HTML_DIR)
    if not os.path.isfile(HTML_FILE):
        with open(HTML_FILE, 'w') as file:
            file.write(preamble)
            print('\nCreated file {}.'.format(HTML_FILE))
    with open(HTML_FILE, 'a') as file:
        file.write(html_text)
        print('\nWrote HTML file {0}.'.format(HTML_FILE))

def add_to_protein_hitdict(taxon, protein, hit_accession_no, hit_description, ortholog_group):
    global protein_hitdict
    # Do not write the same protein twice, but add the protein to the string
    print('Adding to protein_hitdict...')
    if hit_accession_no in protein_hitdict:
        protein_hitdict[hit_accession_no][3] += ' '+ortholog_group
        print('{0} added to existing protein.'.format(hit_accession_no))
    else:
        protein_hitdict[hit_accession_no] = [taxon, protein, hit_description, ortholog_group]
        print('{0} added as new protein.'.format(hit_accession_no))

def most_frequent(List):
    return max(set(List), key = List.count)

def write_protein_hitdict_to_file(protein_hitdict):
    print('Writing protein files and doing needle alignments')
    for key, value in protein_hitdict.items():
        print('protein_hitdict:\n{0} -> {1}'.format(key, value))
        # Evaluate which protein was mostly identified as a homolog (in order to include it in the MSA)
        closest_homolog_list = value[3].split()
        # Delete all "None" elements from the list
        closest_homolog_list = list(filter(lambda x: x!= 'None', closest_homolog_list))
        if len(closest_homolog_list) > 0:
            closest_homolog = most_frequent(closest_homolog_list)
        else:
            closest_homolog = ''
        print('Closest homolog identified as "{0}"'.format(closest_homolog))
        if closest_homolog == '':
            CLOSEST_HOMOLOG_FILE = ''
        elif closest_homolog != value[1]:
            CLOSEST_HOMOLOG_FILE = '{0}/data/proteins/{1}.fasta'.format(APPLICATION_PATH, closest_homolog)
        else:
            CLOSEST_HOMOLOG_FILE = ''
        # VEGF-A.fasta needed to be added manually to data/proteins/ as the specific isoforms are stored in value[0] and
        # this would lead to a file not found error...
        PROT_FILE = '{0}/data/protein_results/{1}.html'.format(APPLICATION_PATH, value[0])
        QUERY_FILE1 = '{0}/data/proteins/{1}.fasta'.format(APPLICATION_PATH, key)
        if not os.path.isfile(QUERY_FILE1):
            download_proteins('{0}/data/proteins'.format(APPLICATION_PATH), {key: [key]})
            time.sleep(1)
        QUERY_FILE2 = '{0}/data/proteins/{1}.fasta'.format(APPLICATION_PATH, value[1])
        VEGF_signature = '{0}/data/proteins/VEGF_signature.fasta'.format(APPLICATION_PATH)
        # Simple one-to-one comparison
        #bash_command = 'needle {0} {1} -gapopen 10 -gapextend 0.5 stdout'.format(QUERY_FILE1, QUERY_FILE2)
        # MSA including additonally the VEGF signature uding edialign
        #bash_command = 'cat {0} {1} {2} | edialign -filter'.format(QUERY_FILE1, QUERY_FILE2, VEGF_signature)
        # MSA using emma (clustalx)
        bash_command = 'cat {0} {1} {2} {3}| emma -filter -osformat2 msf -dendoutfile /dev/zero'.format(QUERY_FILE1, QUERY_FILE2, VEGF_signature, CLOSEST_HOMOLOG_FILE)
        comment = 'Making alignment of {0} with {1} and {2}:\n'.format(value[1], key, closest_homolog)
        alignment = execute_subprocess(comment, bash_command)
        # Trim header from msf file (necessary for emma)
        alignment = alignment.split('//')[1]
        # Trim comments from the alignment text blob (necessary for needle)
        line_list = alignment.split('\n')
        alignment = ''
        for line in line_list:
            if not line.startswith('#'):
                alignment += line+'\n'
        if os.path.exists(PROT_FILE):
            append_or_write = 'a' # append if PROT_FILE exists
        else:
            append_or_write = 'w' # make a new file if PROT_FILE does not exist
            with open(PROT_FILE, append_or_write) as handle:
                preamble = '''<html>
            <head>
            <title>Individual alignments</title>
            </head>
            <body>'''
                handle.write(preamble)
        with open(PROT_FILE, append_or_write) as handle:
            handle.write('<a id={0}>{0} -> {1}</a>\n'.format(key, value))
            handle.write('<pre>\n{0}\n</pre>'.format(alignment))
            print('Writing {0} -> {1} to {2}'.format(key, value, PROT_FILE))

def get_protein_data(taxon):
    global TOTAL_COUNT, TOTAL_UNKNOWN_COUNT, conn
    if args.directory == '':
        return []
    else:
        new_protein_data = {}
        print('\n---------------')
        print('Analyzing taxon {0}'.format(taxon).upper())
        print('---------------\n')
        for protein, protein_data in master_dictionary.items():
            list_to_scrutinize_ortholog = []
            list_to_scrutinize_paralog = []
            list_to_scrutinize_unknown = []
            #print('\n\nSTART\nPROTEIN: {0}\nTAXON: {1}\n\n'.format(protein, taxon))
            #print('master_dictionary:\n{0}\n\n'.format(master_dictionary))
            #time.sleep(2)
            # Here we analyze the blasts hits for each TAXON/PROTEIN pair:
            # a) number of blast hits
            # b) number of false positive according to the description (categorized in a list)
            #
            # Open the individual results file for a protein and a phylum
            XML_RESULTS_FILE = '{0}{1}/blast_results_{2}.xml'.format(args.directory, protein, taxon)
            # Create database connection
            #conn = create_connection(DATABASE_FILE)

            # Get extended information about the false-positive hits
            negative_dict = {}
            for related_protein in protein_data[4]:
                negative_dict[related_protein] = 0
            # Get extended information about the true-positive hits
            positive_dict = {}
            for synonym in protein_data[3]:
                positive_dict[synonym] = 0

            print('\nAnalyzing {0} ({1}):'.format(protein, taxon).upper())
            blastp_result = SearchIO.read(XML_RESULTS_FILE, "blast-xml")
            number = len(blastp_result)
            print('Number of blastp results: {0}'.format(number))
            # counter for hits
            i = 0
            # counter for unknown proteins
            u = 0
            for hit in blastp_result:
                found = False
                # get id and accession number
                # Blast has changed the data it returns in the Hit_id field of the xml file.
                # It used to report both id and accession number, but now returns only the
                # the accession number appended with the version (e.g. "ref|XP_015200204.1|")
                # the Hit_accession field of the xml file contains the accession number without
                # the version extension
                hit_id = hit.id.split('|')[1]
                hit_accession_no = hit.accession
                print('Analyzing #{0} (from {1} {2} {3} homologs): hit_id:{4} - hit_accession_no:{5}'.format(i+1, number, taxon, protein_data[3][0], hit_id, hit_accession_no))
                # Extract species name; some entries are not according to the specs
                # and we need to handle that...
                try:
                    species = re.search(r'\[(.*?)\]',hit.description).group(1)
                except Exception as err:
                    print('Non-standard protein description ({0}). Continuing without species information...'.format(hit.description))
                    species = 'unknown'
                #print('Checking false- or true-positivity.')
                # If one protein category matches, we need to abort searching to avoid items being counted twice or more often if
                # their description contains a repetition like "VEGF-A (Vacular endothelial growth factor-A)".
                # The Try statement is simply for escaping both for-loops
                ortholog_group = 'None'
                try:
                    # Try synonyms first (more likely to find a category?)
                    for synonym in protein_data[3]:
                        if synonym.lower() in str(hit.description).lower():
                            positive_dict[synonym] += 1
                            print('True positive found: {0}'.format(synonym))
                            #list_to_scrutinize.append([taxon, protein, synonym, hit_id, hit_accession_no, hit.description])
                            ortholog_group = protein_data[3][0]
                            list_to_scrutinize_ortholog.append([taxon, protein, 'ortholog', hit_id, hit_accession_no, hit.description])
                            raise category_found
                    for related_protein in protein_data[4]:
                        if related_protein.lower() in str(hit.description).lower():
                            negative_dict[related_protein] += 1
                            print('Potential false-positive found: {0}'.format(related_protein))
                            #list_to_scrutinize.append([taxon, protein, related_protein, hit_id, hit_accession_no, hit.description])
                            # related_protein should always give the same name!!!!
                            # This is why we need the synonym dictionary
                            try:
                                ortholog_group = synonym_dictionary[related_protein]
                                list_to_scrutinize_paralog.append([taxon, protein, 'paralog', hit_id, hit_accession_no, hit.description])
                            except Exception as err:
                                ortholog_group = 'None'
                                list_to_scrutinize_paralog.append([taxon, protein, 'paralog', hit_id, hit_accession_no, hit.description])
                            raise category_found
                    # This part of the Try block gets only executed when both for-loops are finshing without
                    # raising a category_found exception
                    what_kind_of_protein = run_backcheck_blast(hit_accession_no, protein_data)
                    if what_kind_of_protein == 'synonym':
                        positive_dict[synonym] += 1
                        TOTAL_COUNT += 1
                        print('True positive found after backcheck blast: {0}'.format(protein_data[3][0]))
                        ortholog_group = protein_data[3][0]
                        list_to_scrutinize_ortholog.append([taxon, protein, 'ortholog', hit_id, hit_accession_no, hit.description])
                    elif what_kind_of_protein == 'related_protein':
                        negative_dict[related_protein] += 1
                        TOTAL_COUNT += 1
                        #print('Potential false-positive found: {0}'.format(related_protein))
                        print('Potential false-positive found.')
                        ortholog_group = 'None'
                        list_to_scrutinize_paralog.append([taxon, protein, 'paralog', hit_id, hit_accession_no, hit.description])
                    else:
                        # These are all the unknown proteins, that even a backcheck blast cannot identify
                        list_to_scrutinize_unknown.append([taxon, protein, 'unknown', hit_id, hit_accession_no, hit.description])
                        ortholog_group = 'None'
                        u += 1
                        TOTAL_COUNT += 1
                        TOTAL_UNKNOWN_COUNT += 1
                except category_found:
                    TOTAL_COUNT += 1
                # BLASTHIT = [id, accession_no, species, fasta_description, ortholog group (only if determined with high probability)]
                # id is a string (used to be an integer)!
                BLASTHIT = [hit_id, hit_accession_no, species, hit.description, ortholog_group]
                add_to_protein_hitdict(taxon, protein, hit_id, hit.description, ortholog_group)
                with conn:
                    # If the id is already in the database, the insertion fails and 0 is returned (otherwise the id is returned)
                    print('\n\nINSERT PROTEIN:\n{0}\n\n'.format(BLASTHIT))
                    print(str(db_insert_protein(BLASTHIT, True))+'\n')
                i += 1
            # Write ortholog, paralog and unknown lists individually for each protein/taxon
            if len(list_to_scrutinize_ortholog) > 0:
                write_to_html_scrutinize_file(protein, taxon, list_to_scrutinize_ortholog)
            else:
                print('No orthologs to write to HTML file.')
            if len(list_to_scrutinize_paralog) > 0:
                write_to_html_scrutinize_file(protein, taxon, list_to_scrutinize_paralog)
            else:
                print('No paralogs to write to HTML file.')
            if len(list_to_scrutinize_unknown) > 0:
                write_to_html_scrutinize_file(protein, taxon, list_to_scrutinize_unknown)
            else:
                print('No unknown proteins to write to HTML file.')

            print('Analysis for {0}/{1} completed.\n'.format(taxon, protein))
            print('Number of true-positives:'.format())
            control_counter = 0
            for key, value in positive_dict.items():
                print('                         '+key, str(value))
                control_counter += value
            print('Number of false-positives:'.format())
            for key, value in negative_dict.items():
                print('                         '+key, str(value))
                control_counter += value
            print('Number of unidentified homologous proteins: {0}\n'.format(u))
            control_counter += u
            print('{0} out of {1} hits were analyzed and {2} could be assigned into categories.\n\n'.format(i, number, control_counter))
            new_protein_data[protein] = [number, negative_dict, positive_dict]
        print('Analysis for all proteins in taxon {0} completed.\n'.format(taxon))
        print('new_protein_data: {0}'.format(new_protein_data))

        # Returns a dictionary with related proteins and the number of hits for them (according to fasta description)
        print('\nRETURNED PROTEIN_DATA FOR TAXON {0}:\n'.format(taxon))
        for protein, value in new_protein_data.items():
            print('{0}:\n{1}'.format(protein, value))
        return new_protein_data

# This function should return the most common string in the hit description list
def run_backcheck_blast(ID, proteindata):
    # We need to define that this function should use the global variable LAST_BLAST_REPLY_TIME,
    # because we are changing LAST_BLAST_REPLY_TIME inside the function!
    global LAST_BLAST_REPLY_TIME
    try:
        blastp_result = SearchIO.read('data/backcheck/{0}.xml'.format(ID), 'blast-xml')
    except Exception as ex:
        print('No local backcheck blast results for accession no {0}. Need to run remote blast.'.format(ID))
        success = False
        try:
            # Make maximally 5 tries for any protein with increasing waiting time
            i = 0
            while success == False and i < 5:
                # Sleep if the last blast request was answered less than BLAST_WAITING_TIME (default = 60) seconds ago
                # This is necessary in order not to overload the server and become blocked.
                SECONDS_SINCE_LAST_BLAST = time.time()-LAST_BLAST_REPLY_TIME
                if SECONDS_SINCE_LAST_BLAST < BLAST_WAITING_TIME*(i+1):
                    print('Waiting for {0} seconds before initiating new blast request...'.format(round(BLAST_WAITING_TIME*(i+1)-SECONDS_SINCE_LAST_BLAST)))
                    time.sleep(BLAST_WAITING_TIME*(i+1)-SECONDS_SINCE_LAST_BLAST)
                else:
                    print('{0} seconds since last blast result was received. Continue without waiting...'.format(round(SECONDS_SINCE_LAST_BLAST)))
                    # remove this after checking
                    time.sleep(2)
                blast_start_time = time.time()
                print('Executing backcheck blast for accession no {0}. Waiting for results...'.format(ID))
                result_handle = NCBIWWW.qblast("blastp", "nr", ID , hitlist_size = 50, expect = 0.01, format_type = "XML")
                i += 1
                with open('data/backcheck/{0}.xml'.format(ID), 'w') as out_handle:
                    out_handle.write(result_handle.read())
                result_handle.close()
                # Read in the data again to check that they contain really blast results
                with open('data/backcheck/{0}.xml'.format(ID), 'r') as xml_file:
                    xml_string = xml_file.read()
                    # This print command will print lots of lines...
                    #print('data:\n{0}'.format(data))
                    # The database was not searched. Sometimes there is no other error message,
                    # but all query results contain this string if there is no result
                    if '<Statistics_db-num>0</Statistics_db-num>' not in xml_string:
                        success = True
                        print('The blast search Was successful. Continuing...')
                    else:
                        print('No sequences were returned from blast server, repeating query...')
                xml_file.close()
                LAST_BLAST_REPLY_TIME = time.time()
                print('The remote blasting took {0}.'.format(execution_time_str(LAST_BLAST_REPLY_TIME-blast_start_time)))
        except Exception as ex:
            # Sequence id not found error
            if 'Message ID#24 Error' in str(ex):
                print("The blast server did not know the accession number. Error: " + str(ex))
                return 'unknown'
            else:
                print("Something went wrong with the backcheck blast. Perhaps the Blast server did not respond? More about the error: " + str(ex))
    else:
        print('Using local backcheck blast results (data/backcheck/{0}.xml).'.format(ID))
    blastp_result = SearchIO.read('data/backcheck/{0}.xml'.format(ID), 'blast-xml')
    how_many = len(blastp_result)
    print('blastp result should equal 50: {0}'.format(how_many))
    i = 0 # count synonyms
    j = 0 # count related proteins
    k = 1 # count 50 rounds
    for hit in blastp_result:
        if k == 1:
            print('Checking blast result {0}'.format(k), end=' ')
        else:
            print('{0}'.format(k), end = ' ')
        for synonym in proteindata[3]:
            if synonym.lower() in str(hit.description).lower():
                i += 1
                break
        for related_protein in proteindata[4]:
            if related_protein.lower() in str(hit.description).lower():
                j += 1
                break
        k += 1
    print('\n{0} from {1} backblast results confirm that {2} is {3}'.format(i, how_many, ID, proteindata[3][0]))
    print('{0} from {1} backblast results confirm that {2} is a related protein'.format(j, how_many, ID))
    # THRESHOLD needs to be adjusted after manually scrutinizing the results- Maybe 0.3 is sufficient!
    THRESHOLD = 0.5
    if how_many > 0:
        if i/how_many >= THRESHOLD:
            return 'synonym'
        elif j/how_many >= THRESHOLD:
            return 'related_protein'
        else:
            return 'unknown'
    else:
        return 'unknown'

def get_fully_sequenced_genomes(CSV_FILE):
    print('------------------------------------------------ ---------')
    print('--------\nGET NUMBERS OF FULLY SEQUENCED GENOMES\n--------')
    print('------------------------------------------------- --------')
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
    print('taxon_dictionary_print:\n{0}'.format(taxon_dictionary))
    for taxon, taxon_data in taxon_dictionary.items():
        full_genome_dictionary[taxon] = 0
        print('Adding {0} to full_genome_dictionary.'.format(taxon))
    for line in input_csv_file:
        species_name = line['#Organism Name']
        print('species_name = {0}'.format(species_name))
        TAXON_ID, phylum = get_taxon_id_and_phylum(species_name)
        print('TAXON_ID: {0}, phylum: {1}'.format(TAXON_ID, phylum))
        # If phylum cannot be resolved via NBI, this would give an error when trying to increment the
        # number of fully sequenced genomes!
        if phylum != 'unknown':
            full_genome_dictionary[phylum] += 1
        print('Adding +1 to fully sequenced {0}'.format(phylum))
    print(str(full_genome_dictionary))
    return full_genome_dictionary

def create_connection(DATABASE_FILE):
    global conn
    print('Trying to connect to database file {0}...'.format(DATABASE_FILE))
    if os.path.isfile(DATABASE_FILE) and os.stat(DATABASE_FILE).st_size != 0:
        try:
            conn = sqlite3.connect(DATABASE_FILE)
            return conn
        except Exception as err:
            print('Database file {0} not found.'.format(DATABASE_FILE))
            print('Error: {0}'.format(err))
    else:
        print('Database file {0} does not exist or is empty. Trying to create a new one.'.format(DATABASE_FILE))
        if create_sqlite_file(DATABASE_FILE):
            conn = sqlite3.connect(DATABASE_FILE)
            return conn
        else:
            print('Database file {0} not found and unable to create a new one.'.format(DATABASE_FILE))

def db_insert_protein(BLASTHIT, VERBOSE = True):
    global conn
    # Check whether an entry exists already
    if VERBOSE: print('Checking whether an entry with id {0} exists already in the SQLite database... '.format(BLASTHIT[0]), end='')
    #for key in sqlite_protein_dict:
    #    print(key, end=' ')
    if BLASTHIT[0] in sqlite_protein_dict.keys():
        print('Yes')
        # Check whether the entry has already the ortholog group set
        print('ortholog group: {0}'.format(sqlite_protein_dict[BLASTHIT[0]]))
        # Only try to update if the new ortholog_group is not 'None' and the old one is 'None'
        if BLASTHIT[4] != 'None' and sqlite_protein_dict[BLASTHIT[0]][4] == 'None':
            # Update database entry to include ortholog group
            if VERBOSE: print('Entry with id = {0} has no ortholog group information. Trying to update with \"{1}\".'.format(BLASTHIT[0], BLASTHIT[4]))
            query = 'UPDATE protein SET ortholog_group = \'{0}\' WHERE id = \'{1}\''.format(BLASTHIT[4], BLASTHIT[0])
            if VERBOSE: print('query: {0}'.format(query))
            try:
                cur = conn.cursor()
                cur.execute(query)
                result = cur.lastrowid
                conn.commit()
            except sqlite3.Error as err:
                if VERBOSE: print('Database error during update: {0}'.format(err))
                result = 0
            except Exception as err:
                if VERBOSE: print('Unknown error during update: {0}'.format(err))
                result = 0
            else:
                print('Ortholog group of id {0} updated to {1}.'.format(BLASTHIT[0], BLASTHIT[4]))
        else:
            print('Ortholog group is already set ({0}).'.format(BLASTHIT[4]))
            result = 0
    else:
        print('No')
        if VERBOSE: print('No entry with id = {0}. Trying to insert.'.format(BLASTHIT[0]))
        # This is a full database insertion
        # In some fasta descriptions one can really find this char: ' !!!!
        BLASTHIT[2] = BLASTHIT[2].replace('\'', '')
        query = 'INSERT INTO protein (id, accession_no, species, fasta_description, ortholog_group) VALUES (\'{0}\', \'{1}\', \'{2}\', \'{3}\', \'{4}\')'.format(BLASTHIT[0], BLASTHIT[1], BLASTHIT[2], BLASTHIT[3], BLASTHIT[4])
        if VERBOSE: print('query: {0}'.format(query))
        try:
            cur = conn.cursor()
            cur.execute(query)
            result = cur.lastrowid
            conn.commit()
        except sqlite3.Error as err:
            if VERBOSE: print('Database error during insertion: {0}'.format(err))
            result = 0
        except Exception as err:
            if VERBOSE: print('Unknown error during insertion: {0}'.format(err))
            result = 0
    return result

# Give a list of scientific_name, taxon_id and phylum to insert into SQLite database
# There should never be the need to update any entry from this table, only to insert
def db_insert_species(SPECIES_NAME, TAXON_ID, PHYLUM, VERBOSE=True):
    global conn
    if VERBOSE: print('Trying to execute SQL insertion ({0})...'.format(SPECIES_NAME))
    query = 'INSERT INTO species (scientific_name, taxon_id, phylum) VALUES (\'{0}\', {1}, \'{2}\')'.format(SPECIES_NAME, TAXON_ID, PHYLUM)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = conn.cursor()
        cur.execute(query)
        result = cur.lastrowid
        if VERBOSE: print('SQL insertion successful.')
        conn.commit()
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    return result

# Return a list of scientific_name, taxon_id and phylum when queried by scientific_name
def db_retrieve_species(SPECIES_NAME, VERBOSE=True):
    global conn
    if VERBOSE: print('Trying to retrieve SQLite data for {0}...'.format(SPECIES_NAME))
    query = 'SELECT * FROM species WHERE scientific_name = \'{0}\''.format(SPECIES_NAME)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = conn.cursor()
        cur.execute(query)
        result = cur.fetchone()
        #print('\nRESULT: {0}\n'.format(result))
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    # There should be only one row!
    return result

def load_sqlite_table_to_dict(TABLENAME, VERBOSE = True):
    global conn
    if VERBOSE: print('\n\nLoading table {0} from SQLite database...'.format(TABLENAME))
    query = 'SELECT * FROM {0}'.format(TABLENAME)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = conn.cursor()
        cur.execute(query)
        results = cur.fetchall()
    except sqlite3.Error as err:
        # No or empty database file
        if VERBOSE: print('Database error: {0}'.format(err))
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
    else:
        sqlite_dict = {}
        i = 0
        for row in results:
            sqlite_dict[row[0]] = row
            i += 1
        print('Loading of table {0} from SQLite database was successful ({1} entries):\n'.format(TABLENAME, i))
        print('table {0}:\n{1}\n\n'.format(TABLENAME, sqlite_dict))
    return sqlite_dict

def run():
    global APPLICATION_PATH, taxon_dictionary, master_dictionary, synonym_dictionary, blacklist, DATABASE_FILE, LAST_BLAST_REPLY_TIME, BLAST_WAITING_TIME, TOTAL_COUNT, TOTAL_UNKNOWN_COUNT, TOTAL_UNKNOWN_COUNT_HTML_CHECK, conn, sqlite_protein_dict, sqlite_species_dict, protein_hitdict
    TOTAL_COUNT = 0
    TOTAL_UNKNOWN_COUNT = 0
    TOTAL_UNKNOWN_COUNT_HTML_CHECK = 0
    # You can adjust this down until you see that the blat server starts blocking your requests!
    # 60 (seconds) is a very conservative (but slow) estimate, that does not result in blocking
    BLAST_WAITING_TIME = 30
    LAST_BLAST_REPLY_TIME = time.time()-BLAST_WAITING_TIME
    print(LAST_BLAST_REPLY_TIME)
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    #print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_data.py'.format(APPLICATION_PATH)
    CSV_FILE = '{0}/data/genomes.csv'.format(APPLICATION_PATH)
    DATABASE_FILE = '{0}/data/database.sqlite3'.format(APPLICATION_PATH)
    LOGFILE = '{0}/data/logfile.txt'.format(APPLICATION_PATH)
    ANALYSIS_RESULTS_DIR = '{0}/data/analysis_results'.format(APPLICATION_PATH)
    PROTEIN_RESULTS_DIR = '{0}/data/protein_results'.format(APPLICATION_PATH)
    # Remove all old analysis and protein results and recreate the directories
    if os.path.isdir(ANALYSIS_RESULTS_DIR):
        shutil.rmtree(ANALYSIS_RESULTS_DIR)
        os.mkdir(ANALYSIS_RESULTS_DIR)
    else:
        os.mkdir(ANALYSIS_RESULTS_DIR)
    if os.path.isdir(PROTEIN_RESULTS_DIR):
        shutil.rmtree(PROTEIN_RESULTS_DIR)
        os.mkdir(PROTEIN_RESULTS_DIR)
    else:
        os.mkdir(PROTEIN_RESULTS_DIR)
    preamble1, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)
    #print('taxon_dictionary: {0}\n'.format(taxon_dictionary))
    #print("Populating new taxon dictionary")
    new_taxon_dictionary = taxon_dictionary
    preamble2, master_dictionary = load_dictionary('{0}/data/master_dictionary.py'.format(APPLICATION_PATH))
    # Downlaod reference proteins and rename according to the key (human-readable name)
    download_proteins('{0}/data/proteins'.format(APPLICATION_PATH), master_dictionary, True, False)
    synonym_dictionary = make_synonym_dictionary(master_dictionary)
    print('synonym_dictionary:\n{0}'.format(synonym_dictionary))
    blacklist = load_blacklist()
    protein_hitdict = {}
    conn = create_connection(DATABASE_FILE)
    print('{0}'.format(conn))
    sqlite_protein_dict = load_sqlite_table_to_dict('protein', VERBOSE = True)
    sqlite_species_dict = load_sqlite_table_to_dict('species', VERBOSE = True)
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
            # Get taxon_id if a new taxon was added manually to taxon_data.py, but no taxon_id was given
            # This is not necessary if the taxon_id is always manually assigned when adding a new taxon
            # and it also fails when adding a complex taxon
            #try:
            #    new_taxon_id, phylum = get_taxon_id_and_phylum(taxon)
            #except Exception as err:
            #    print('Could not get new taxon id and phylum. Error was: {0}',value(err))
            # execute if no exceptions have occured
            #else:
            #    new_taxon_dictionary[taxon][0] = new_taxon_id
            #
            # NUMBER OF SPECIES IN NCBI DATABASE
            try:
                new_species_number = get_species_number_from_ncbi(taxon, taxon_data)
            except Exception as err:
                print('Could not get new species number. Error was: {0}',value(err))
            else:
                new_taxon_dictionary[taxon][1] = new_species_number
            # NUMBER OF SEQUENCES IN NCBI DATABASE
            try:
                new_sequence_number = get_sequence_number(taxon, taxon_data)
            except Exception as err:
                print('Could not get new sequence number. Error was: {0}',value(err))
            else:
                new_taxon_dictionary[taxon][3] = new_sequence_number
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
    write_protein_hitdict_to_file(protein_hitdict)
    if insert_line_breaks(TAXON_DICTIONARY_FILE) == True:
        print('Successfully formatted the taxon data file.')
    else:
        print('Formating the taxon data file failed.')
    save_stats = '\nAnalyzed sequences (hits resulting from {0} blast searches, {1} animal groups x {2} query sequences, not uniques!): {3} (out of which unclassified: {4})\nUnclassified sequences written to {5}: {6}'.format(len(taxon_dictionary)*len(master_dictionary), len(taxon_dictionary), len(master_dictionary), TOTAL_COUNT, TOTAL_UNKNOWN_COUNT, '/data/primary_blast_results/unknown/multiple_scrutinize.html', TOTAL_UNKNOWN_COUNT_HTML_CHECK)
    with open(LOGFILE, 'a') as log_file:
        log_file.write(save_stats)
    print(save_stats)

if __name__ == '__main__':
    run()
