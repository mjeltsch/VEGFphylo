#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
#
# UBUNTU PACKAGES
#
# This script has been developed to run under Linux and was specifically tested ONLY under Ubuntu 16.04 and 18.04.
# It requires the installation of the following ubuntu packages:
#
# python3-setuptools
# python3-termcolor
# clustalw
# mafft
# dialign-tx
# poa
# probcos
# muscle
# kalign
# amap-align
# proda
# prank
# t-coffee
# phyml
#
# sudo apt -y --show-progress install inkscape python3-setuptools python3-pyqt4 python3-pyqt4.qtopengl python3-pip autoconf t-coffee clustalw mafft dialign-tx poa probcons muscle kalign amap-align proda prank phyml t-coffee imagemagick build-essential libblas-dev liblapack-dev zlib1g-dev libcairo2-dev libcurl4-openssl-dev python3-numpy python3-lxml python3-six
#
#
# For some reason the Ubuntu package names some of the alignment executables differently and some links need to be
# created in order for t-coffee to find them:
#
# sudo ln -s /usr/bin/dialign-tx /bin/dialign-t
# sudo ln -s /usr/bin/clustalw /bin/clustalw2
#
#
# MANUAL INSTALLATIONS
#
# In addition, the pcma executable need to be manually downloaded from http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz
# and compiled because it is not available from the Ubuntu repository:
#
# wget http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz
# tar -xvzf pcma.tar.gz
# cd pcma
# make
# sudo cp pcma /bin
#
# Newick Utilities from http://cegg.unige.ch/newick_utils:
#
# wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
# cd newick-utils-1.6
# sudo cp src/nw_* /usr/local/bin
# cd _blast_results
# ./test_binaries.sh
#
# OPTIONAL INSTALLS
#
# The multicore-enabled version of phyml (phyml-mpi) is not available as a precompiled Ubuntu package and needs
# to be installed manually, but they single-core version works as well (is just slower). The command to execute
# the multicore version is: "mpirun -n 4 phyml-mpi -i " + PHYLIP_ALIGNED_TRIMMED_CODED +  " -d aa -b -1"
# If this script is run on a (headless) server, the xvfb package is required since the ete3 package requires the presence of x.org.
# In this case the script needs to be invoked as follows:
#
# xvfb-run ./do_analysis
#
# PYTHON MODULES
#
# The following Python modules need to be installed:
#
# biopython
# ete3
#
# sudo pip3 install biopython ete3
#
#
# Needs also NCBI blast executables >2.6 from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# for remote access to the blast server and overriding of Ubunut's old blast executables (2.2) in the user's
# ~/.profile file:
#
# # set PATH so it includes user's private blast executables if they exist
#if [ -d "$HOME/bin/blast/bin" ] ; then
#    PATH="$HOME/bin/blast/bin:$PATH"
#    fi
#
#
# Example commands:
# (xvfb-run) ./do_analysis.py blastp download gettree drawtree
# (xvfb-run) ./do_analysis.py psiblast download gettree drawtree
#
# THIS WILL LAST AN ETERNITY (since it does alignment, etc.):
# (xvfb-run) ./do_analysis.py blastp
# (xvfb-run) ./do_analysis.py psiblast
#
#
# For large datasets, the PASTA algorithm is used for alignment (https://github.com/smirarab/pasta)
# and needs to be installed according to the instructions on github.
# Install on ubuntu python-setuptools
# Also git clone https://github.com/smirarab/sate-tools-linux.git
# Both should be git-cloned into ~/bin/
#
#
# WHAT DOES THIS SCRIPT DO?
#
# This script takes all sequence ids from the master_dictionary.py and performs pblasts requesting both html and xml
# files seperately for each taxon that is mentioned in the taxon_data.py file. The results are saved under the
# primary_blast_results directory (in protein-specific subdirectories).
# The master_dictionary.py file contains also a list of all aliases of the protein names and a list of all false friends.
# The taxon_dictionary.py file also stores the statistical data associated with each taxon that is later generated with
# the 2_analysis.py script.

import argparse, Bio, os, sys, shutil, re, time, csv, operator
from functools import wraps
import errno
import signal
import os
from Bio import SearchIO
from Bio import Entrez
from Bio import SeqIO
from Bio import Phylo
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, SequenceFace, ImgFace, SVGFace, faces, add_face_to_node
from os.path import basename, dirname, splitext, split
from phylolib import load_dictionary, get_phylum_from_NCBI, get_script_style_section

parser = argparse.ArgumentParser()
parser.add_argument("local_remote", help = "local (uses local blast), remote (Use NCBI's servers)")
args = parser.parse_args()
if args.local_remote == 'remote':
    REMOTE = True
    print('This run uses NCBI\'s blast servers. This will take a long time!')
else:
    REMOTE = False
    print('This run uses the locally installed blast server.')

def getCSVfile(CSV_FILE):
    # Open/download the list of fully sequenced genomes
    try:
        input_csv_file = csv.DictReader(open(CSV_FILE))
    except FileNotFoundError as err:
        print(err)
        # Download list of all fully sequences animal genomes:
        URL = 'https://www.ncbi.nlm.nih.gov/genomes/solr2txt.cgi?q=%5Bdisplay()%5D.from(GenomeBrowser).usingschema(%2Fschema%2FGenomeAssemblies).matching(group%3D%3D%5B%22Animals%22%5D)&fields=organism%7COrganism%20Name%2Clineage%7COrganism%20Groups%2Csize%7CSize(Mb)%2Cchromosomes%7CChromosomes%2Corganelles%7COrganelles%2Cplasmids%7CPlasmids%2Cassemblies%7CAssemblies&filename=genomes.csv&nolimit=on'
        r = requests.get(URL)
        with open(CSV_FILE, 'wb') as file:
            file.write(r.content)
        input_csv_file = csv.DictReader(open(CSV_FILE))
    return input_csv_file

def write_to_file(FILENAME, content):
    with open(FILENAME, 'a+') as out_handle:
        print('Writing to {0}:\n{1}'.format(FILENAME, content))
        out_handle.write(content)

def convert_to_int(x):
    return int(x[4])

def sort_csv(CSV_FILE):
    data = []
    with open(CSV_FILE, 'r') as file:
        data2 = csv.reader(file, delimiter='\t')
        for row in data2:
            data.append(tuple(row))
            print('ROW:\n{0}'.format(row))
        # Use the 5th column
        print('DATA:\n{0}'.format(data))
        # This gets the 5th element of the list and converts it into an integer for sorting
        # This line has given me headaches as the list was never sorted (or never sorted correctly)
        #sortedlist = sorted(data, key = int(operator.itemgetter(4)))
        sortedlist = sorted(data, key = convert_to_int, reverse = True)
    print('SORTEDLIST:\n{0}'.format(sortedlist))
    # Write the sorted list into a new CSV file
    with open(FINAL_RESULTS_SORTED, 'w') as file:
        fileWriter = csv.writer(file, delimiter = '\t')
        for row in sortedlist:
            print('ROW: {0}'.format(row))
            fileWriter.writerow(row)

def convert_to_html(CSV_FILE):
    FINAL_RESULTS_HTML = '{0}/full_genomes.html'.format(APPLICATION_PATH)
    with open(FINAL_RESULTS_HTML, 'w+') as outfile:
        outfile.write('<html>\n<head>\n<title>Blast hist counts for individual species</title>\n')
        outfile.write(get_script_style_section())
        outfile.write('</head>\n<body>\n<table>')
        outfile.write('<thead><th>Phylum</th><th>Species</th><th>Blast hits</th><th>Total number of protein sequences</th><th>List of hit ids</th></thead>\n<tbody>\n')
        with open(CSV_FILE, 'r') as infile:
            line = infile.readline()
            while line:
                fields = line.split('\t')
                # Write first 5 columns to table row to html file
                outfile.write('<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td><td>'.format(fields[0], fields[1], fields[2], fields[3]))
                # Write last column (which contains zero to many protein_ids), first split into individual protein_ids
                id_list = fields[5].split(',')
                #print('id_list: {0}'.format(id_list))
                # If the last column does not contain any protein_ids, write nothing into the last table cell
                if id_list == ['[]\n']:
                    outfile.write('&nbsp;')
                else:
                    linklist = ''
                    for item in id_list:
                        isolated_id = item.strip(' \'|[]\n').split('|')[-1]
                        #print('item:\n{0}\nstripped_id:\n{1}isolated_id:\n{2}'.format(item, stripped_id, isolated_id))
                        print('add_to_linklist: <a href="https://www.ncbi.nlm.nih.gov/protein/{0}" target="_blank">{0}</a><br>'.format(isolated_id))
                        linklist += '<a href="https://www.ncbi.nlm.nih.gov/protein/{0}">{0}</a><br>'.format(isolated_id)
                    outfile.write('<button onclick="myFunction(\'align_{0}\')">Show/Hide hits</button>'.format(fields[1].replace(' ', '_')))
                    outfile.write('<div id="align_{0}" style="display: none;">{1}</div>'.format(fields[1].replace(' ', '_'), linklist))
                outfile.write('</td></tr>\n')
                line = infile.readline()
            outfile.write('</tbody>\n</table>\n</body>\n<html>\n')

def blastp(SPECIES):
    global final_result_string
    # Get phylum for the species
    phylum = get_phylum_from_NCBI(SPECIES, VERBOSE=True)
    # Make a subfolder for each species if it does not exist
    SPECIES_DIR = '{0}{1}'.format(RESULT_DIR, SPECIES.replace(' ', '_'))
    if not os.path.isdir(SPECIES_DIR):
        os.mkdir(SPECIES_DIR)
    # Database
    BLAST_DATABASE = 'nr'
    EVALUE = 0.1
    OPTIONAL_BLAST_NO = 100
    temp_protein_hit_list = []
    total_num_sequences = 0
    for protein, data in master_dictionary.items():
        PROTEIN_ID = data[0]
        # This was temporarily used to use a specific species
        #SPECIES = 'txid6239'
        if REMOTE:
            print('\nRunning remote blastp against {0}, limited to organism {1}, requesting {2} results in XML format, with the following query sequence: {3}.'.format(BLAST_DATABASE, SPECIES, OPTIONAL_BLAST_NO, PROTEIN_ID))
            ENTREZ_QUERY = '({0}[organism])'.format(SPECIES)
            try:
                print('Start remote blast job...')
                result_handle = NCBIWWW.qblast('blastp', BLAST_DATABASE, PROTEIN_ID, hitlist_size = OPTIONAL_BLAST_NO, expect = EVALUE, entrez_query = ENTREZ_QUERY, format_type = 'XML')
                # Write blast result to xml or html file
                outfile = '{0}{1}{2}.xml'.format(RESULT_DIR, SPECIES.replace(' ', '_'), PROTEIN_ID)
                print('Opening file {0} for writing'.format(outfile))
                with open(outfile, "w") as out_handle:
                    out_handle.write(result_handle.read())
                    print('{0} written to disk.'.format(outfile))
                result_handle.close()
            except Exception as ex:
                print('Remote blasting failed. The error was: '.format(ex))
            else:
                print('Remote blasting completed successfully.')
        else:
            print('\nRunning local blastp against {0}, limited to organism {1}, requesting {2} results in XML format, with the following query sequence: {3} ({4}).'.format(BLAST_DATABASE, SPECIES, OPTIONAL_BLAST_NO, protein, PROTEIN_ID))
            ENTREZ_QUERY = '({0}[organism])'.format(SPECIES)
            try:
                stdout = ''
                stderr = ''
                print('Start local blast job...')
                outfile = '{0}{1}/{2}.xml'.format(RESULT_DIR, SPECIES.replace(' ', '_'), PROTEIN_ID)
                infile = '{0}/data/reference_proteins/{1}.fasta'.format(APPLICATION_PATH, PROTEIN_ID)
                blastx_cline = NcbiblastpCommandline(query = infile, db = '/home/cloud-user/blast_database/nr_{}'.format(SPECIES.replace(' ','_')), evalue = EVALUE, outfmt=5, out = outfile)
                # Uncomment the following line if you want to ignore previous results which are locally cached
                if not os.path.isfile(outfile):
                    stdout, stderr = blastx_cline()
            except Exception as ex:
                print('Local blasting failed. The error was:\n{0}\nstdout:\n{1}\stderr:\n{2}'.format(ex, stdout, stderr))
            else:
                print('Local blasting completed successfully.')
        result_handle = open(outfile)
        try:
            number_of_results = '-'
            # Use two different methods to access the XML file (is there any other way to do this?)
            # Method 1 (to get the individual ids of the hits)
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                if blast_record.num_sequences_in_database > total_num_sequences:
                    total_num_sequences = blast_record.num_sequences_in_database
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        print(alignment.hit_id)
                        # Add hit to temporary list if it is unique
                        if alignment.hit_id not in temp_protein_hit_list:
                            isolate_id = alignment.hit_id.split('|')[1]
                            exclude_file = '{0}/data/proteins-exclude/{1}/{2}.fasta'.format(APPLICATION_PATH, phylum, isolate_id)
                            print('Checking for file {0}... '.format(exclude_file), end='')
                            if not os.path.isfile(exclude_file):
                                temp_protein_hit_list.append(alignment.hit_id)
                                print('Not found, adding to alignment hit list.')
                            else:
                                print('Found! Excluding from alignment hit list.')
            # Method 2 (to get the number of hits)
            blastp_result = SearchIO.read(outfile, 'blast-xml')
            number_of_results = len(blastp_result)
            print('{0} {1} homologs found in a total of {2} protein sequences for {3} ({4}).'.format(number_of_results, protein, total_num_sequences, SPECIES, phylum))
        except Exception as ex:
            print('Could not parse blast-xml file {0}. File might not exist or is empty.'.format(outfile))
    if total_num_sequences >= MINIMUM:
        #print('\nAdding to CSV file:\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(phylum, SPECIES, len(temp_protein_hit_list), total_num_sequences, int(total_num_sequences// (len(temp_protein_hit_list)+0.00001)), temp_protein_hit_list))
        write_to_file(FINAL_RESULTS, '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(phylum, SPECIES, len(temp_protein_hit_list), total_num_sequences, int(total_num_sequences//(len(temp_protein_hit_list)+0.00001)), temp_protein_hit_list))

def run():
    global master_dictionary, RESULT_DIR, REMOTE, APPLICATION_PATH, final_result_string, FINAL_RESULTS, FINAL_RESULTS_SORTED, MINIMUM

    # Only record if 2000 or more protein sequences exist for a species
    MINIMUM = 2000
    # Determine directories of script (in order to load & save the data & log files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    # List of fully sequenced genomes
    CSV_FILE = '{0}/data/genomes.csv'.format(APPLICATION_PATH)
    RESULT_DIR = '{0}/data/individual_species/'.format(APPLICATION_PATH)
    FINAL_RESULTS = '{0}/data/full_genome_data.csv'.format(APPLICATION_PATH)
    if os.path.isfile(FINAL_RESULTS):
        os.remove(FINAL_RESULTS)
    FINAL_RESULTS_SORTED = '{0}/data/full_genome_data_sorted.csv'.format(APPLICATION_PATH)
    preamble1, master_dictionary = load_dictionary('{0}/data/master_dictionary.py'.format(APPLICATION_PATH), VERBOSE = False)
    input_csv_file = getCSVfile(CSV_FILE)
    for line in input_csv_file:
        species_name = line['#Organism Name']
        blastp(species_name)
    print('Writing of {0} complete.'.format(FINAL_RESULTS))
    sort_csv(FINAL_RESULTS)
    print('Sorting of {0} complete.'.format(FINAL_RESULTS))
    convert_to_html(FINAL_RESULTS_SORTED)

if __name__ == '__main__':
    run()
