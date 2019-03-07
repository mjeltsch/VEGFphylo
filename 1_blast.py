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

import argparse, subprocess, Bio, os, sys, shutil, re, time, datetime, socket, random
from functools import wraps
import errno
import signal
import os
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
from phylolib import execution_time_str, load_dictionary, load_blacklist

parser = argparse.ArgumentParser()
parser.add_argument("local_remote", help = "local (tries to use local - previously received - data and only does remote blast searches if necessary), remote (do all new blast searches)")
args = parser.parse_args()

class TimeoutError(Exception):
    print("No data received within 30 seconds.")

def timeout(error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(SECONDS)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator

if args.local_remote == 'local':
    REMOTE = 'local'
    print('This run uses - if available - blastp search results from previous runs to save time.')
elif args.local_remote == 'remote':
    REMOTE = 'remote'
    print('This run does execute a new remote blastp search for every protein. This will take a long time!')
else:
    pass

def create_subdirectory(PROTEIN):
    if not os.path.exists(PROTEIN):
        os.mkdir(PROTEIN)

# Animal SVG silouette images from http://phylopic.org (public domain) or own creations
# It would be nice to have a lungfish VEGF-C, but the genomes have not been made publicly available!
# Protopterus annectens
# Neoceratodus forsteri
# Lepidosiren paradoxa

# This just puts breaks in the terminal output for better navigation of the stdout file
def page_break(section_name, protein):
    print("\n\n-------------------------------------\n-------------------------------------\n-------------------------------------\n\n")
    print("ENTERING ANALYSIS STAGE \"" + section_name + "\" FOR " + protein)
    print("\n\n-------------------------------------\n-------------------------------------\n-------------------------------------\n\n")

# The remove_extension function removes the extension from a filepath, e.g.
# /home/user/bioinformatics/data/sequence.fasta -> /home/user/bioinformatics/data/sequence
def remove_extension(path):
    path, filename_with_extension = os.path.split(path)
    #filename_with_extension = os.path.basename(path)
    filename_without_extension, extension = os.path.splitext(filename_with_extension)
    return os.path.join(path, filename_without_extension)

def execute_subprocess(comment, bash_command, working_directory='.'):
    print("\n" + comment, bash_command)
    process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=working_directory)
    output, error = process.communicate()
    process_status = process.wait()
    if output.decode('utf-8') != '':
        print("Output: " + str(output))
    if error.decode('UTF-8') != '':
        print("Error: " + str(error))

def write_dict_to_file(dictionary, file_name):
    try:
        with open(file_name, "w") as file:
            file.write(str(dictionary))
            file.close()
            return True
    except Exception as ex:
        print("Could not write dictionary to file. Error: " + str(ex))
        return False

def read_file_to_dict(file_name):
    try:
        with open(file_name, "r") as file:
            string = file.read()
            file.close()
            #print('File content:\n' + string)
            dictionary = eval(string)
            return dictionary
    except Exception as ex:
        print('Could not read dictionary from file {0}. Error: '.format(file_name) + str(ex))
        return {}

def parse_blast_result(PROTEIN, PROTEIN_DATA, TAXON):
    PROTEIN_ID, SUBRANGE, OPTIONAL_BLAST_NO, SYNONYMS, RELATED_PROTEINS = PROTEIN_DATA
    # Parse Blastp result and write results to SEQUENCE_LIST
    with open(SUMMARY_FILE, "a") as results_summary:
        #print(os.getcwd() + " " + BLAST_XMLFILE)
        with open(BLAST_XMLFILE) as blastresults:
            blast_records = NCBIXML.parse(blastresults)
            print("Parsing XML file...", end = '')
            blast_records = list(blast_records)
            for blast_record in blast_records:
                how_many_hits = len(blast_record.descriptions)
                if blast_record.num_sequences_in_database > 0:
                    quotient = str(round(blast_record.num_letters_in_database/blast_record.num_sequences_in_database, 1))
                else:
                    quotient = 'N/A'
                results_summary.write(PROTEIN + "\t" + TAXON + "\t" + str(blast_record.num_sequences_in_database) + "\t" + str(blast_record.num_letters_in_database) + "\t" + quotient + "\t" + str(how_many_hits) + "\n")
                for description in blast_record.descriptions:
                    results_summary.write("\t\t\t\t\t\t" + str(description.num_alignments) + "\t" + str(description.score) + "\t" + str(description.e) + "\t" + description.title + "\n")
                    #print("\ndescription: " + description.title)
                    #print("\ndescription: " + str(description.score))
                    #print("\ndescription: " + str(description.e))
                    #print("\ndescription: " + str(description.num_alignments))
                # Probably not necessay to list the high scoring pairs
                #for alignment in blast_record.alignments:
                #    how_many_hsps = len(alignment.hsps)
                #    for hsp in alignment.hsps:
                #        results_summary.write(str(how_many_hsps) + "\t\t\t\t\t\t" + alignment.title + "\n")
                # To find out about the record structure
                #attrs = dir(blast_record)
                #print("\nPrinting attributes:\n")
                #for i in range(0, len(attrs)):
                #    attribut = getattr(blast_record, attrs[i])
                #    print("\n" + str(i) + ".\n " + str(attribut) + "\n")
        print(" completed.\n")

def blastp(PROTEIN, PROTEIN_DATA, TAXON, TAXON_DATA, FORMAT):
    # Database
    BLAST_DATABASE = 'nr'
    EVALUE = 0.1
    PROTEIN_ID, SUBRANGE, OPTIONAL_BLAST_NO, SYNONYMS, RELATED_PROTEINS = PROTEIN_DATA
    WHATFORMAT = FORMAT[0][0]
    try:
        WHATFORMAT += '+ '+FORMAT[1][0]
    except NameError:
        pass
    # Memorize start time to figure out how long the whole script execution takes
    start_time = time.time()
    print('\nRunning remote blastp against {0}, subsection {1} ({2}), requesting {3} results in {4} format, with the following query sequence: {5}. Blast job started at {6}'.format(BLAST_DATABASE, TAXON, TAXON_DATA[2], OPTIONAL_BLAST_NO, WHAT_FORMAT, PROTEIN, str(datetime.datetime.now())[:-7]))
    # The hitlist_size seems to be ignored by the wrapper/service????
    # Reference: http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
    # The default is to use the sequence_id for the blasting. However, if no sequence_id is available, we can provide the amino acid sequence locally
    # Check here whether the sequence_id is and entry in the local sequence_disctionary and get the sequence for blasting
    #result_handle = NCBIWWW.qblast("blastp", BLAST_DATABASE, PROTEIN_ID[0], hitlist_size = str(OPTIONAL_BLAST_NO), expect = EVALUE, query_from = ALIGNMENT_TRIMMING[0], query_to = ALIGNMENT_TRIMMING[1], entrez_query='NOT '+EXCLUDE+'[organism]')
    # Alternatively to the taxon name, the taxon id can be specified like e.g. "txid9606"[organism]
    # When selecting the format ("Format_type"), please remember that HTML cannot be easily parsed. Use XML instead!
    ENTREZ_QUERY = '\"' + TAXON +'\"[organism]'
    if SUBRANGE == None:
        STARTPOS = ''
        ENDPOS = ''
    else:
        STARTPOS = SUBRANGE.split('-')[0]
        ENDPOS = SUBRANGE.split('-')[1]
    try:
        # Repeat for HTML output (how else to do this easily???)
        # Converting XML to html?
        # xsltproc --novalid blast2html.xsl blast.xml
        for return_format in FORMAT:
            result_handle = NCBIWWW.qblast('blastp', BLAST_DATABASE, PROTEIN_ID, hitlist_size = OPTIONAL_BLAST_NO, expect = EVALUE, entrez_query = ENTREZ_QUERY, format_type = return_format[0], query_from = STARTPOS, query_to = ENDPOS)
            # Write blast result to xml or html file
            with open(return_format[1], "w") as out_handle:
                out_handle.write(result_handle.read())
                print('{0}/{1} written to disk.'.format(PROTEIN, return_format[1]))
            result_handle.close()
    except Exception as ex:
        how_long = execution_time_str(time.time()-start_time)
        print('Blasting failed. It took {0}. The error was: '.format(how_long, ex))
        return False
    else:
        how_long = execution_time_str(time.time()-start_time)
        print('Blasting completed successfully in {0}.\n'.format(how_long))
        return True

# Wrapping blastp in a timeout function
blastp = timeout()(blastp)

def run():
    blacklist = load_blacklist()
    global DATA_DIR, BLAST_XMLFILE, BLAST_HTMLFILE, APPLICATION_PATH, SUMMARY_FILE, REMOTE, SECONDS

    # This enables simultaneous output to the terminal and a logfile
    class logfile(object):
        def __init__(self, filepath = 'logfile.txt'):
            self.terminal = sys.stdout
            self.log = open(filepath, "a")
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
        def flush(self):
            #For python 3 compatibility
            pass

    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    FIRST_BLAST = True
    print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    # Loading data file
    preamble1, master_dictionary = load_dictionary('{0}/data/master_dictionary.py'.format(APPLICATION_PATH))
    preamble2, taxon_dictionary = load_dictionary('{0}/data/taxon_data.py'.format(APPLICATION_PATH))
    LOGFILE = 'logfile.txt'
    DATA_DIR = 'data/primary_blast_results'
    if REMOTE == 'local':
        print('Using .xml blast result files from directory {0}/{1} (skipping remote blasts whenever possible).'.format(os.getcwd(), DATA_DIR))
    elif REMOTE == 'remote':
        # Rename old results and create a new result directory
        OLD_DATA_DIR = 'data/'+datetime.datetime.now().strftime("%y%m%d_%H%M%S")
        os.rename(DATA_DIR, OLD_DATA_DIR)
        print('\nCreating new directory {0} to store primary blast results.\n'.format(DATA_DIR))
    # Create directory for the downloaded and generated data
    if not os.path.exists(DATA_DIR):
        os.mkdir(DATA_DIR)
    # Direct all terminal output additionally into a logfile
    sys.stdout = logfile(DATA_DIR + "/" + LOGFILE)
    # Change the current working directory
    os.chdir(DATA_DIR)
    SUMMARY_FILE = '{0}/{1}/summary.csv'.format(APPLICATION_PATH, DATA_DIR)
    # Create summary file
    with open(SUMMARY_FILE, 'w') as results_summary:
        results_summary.write("PROTEIN\tTAXON\tsequences in db\tletters in db\tavrg length\thits\talignments\tscore\te-value\tdescription\n")
    for protein, protein_data in master_dictionary.items():
        # Default number to return from blast search
        if protein_data[2] == None:
            protein_data[2] = 50
        # Create a subdirectory named according to the protein and change cwd into it
        create_subdirectory(protein)
        os.chdir(protein)
        for taxon, taxon_data in taxon_dictionary.items():
            if taxon not in blacklist:
                print('Analyzing {0}:'.format(taxon))
                BLAST_XMLFILE = 'blast_results_{0}.xml'.format(taxon)
                BLAST_HTMLFILE = 'blast_results_{0}.html'.format(taxon)
                # Check whether we have both XMl and HTML version of a previous blast result file for the protein
                if not os.path.isfile(BLAST_XMLFILE) or not os.path.isfile(BLAST_HTMLFILE):
                    # Do only request the format that is not locally available
                    if not os.path.isfile(BLAST_XMLFILE) and not os.path.isfile(BLAST_HTMLFILE):
                        format = [['XML', BLAST_XMLFILE], ['HTML', BLAST_HTMLFILE]]
                    else:
                        if not os.path.isfile(BLAST_HTMLFILE):
                            format = [['HTML', BLAST_HTMLFILE]]
                        if not os.path.isfile(BLAST_XMLFILE):
                            format = [['XML', BLAST_XMLFILE]]
                    # If no local results are found, loop as long until the blasting succeeds or for maximally MAX_ITER tries
                    iteration = 0
                    MAX_ITER = 10
                    while True and iteration < MAX_ITER:
                        try:
                            # Customize timeout period before starting a new blastp request
                            # taxon_data[3] is the number of protein sequences for that taxon in the nr protein database
                            SECONDS = int(round(taxon_data[3]**(1/9)*1000, 0))
                            print('If no results are received within {0} seconds, the blast will be cancelled'.format(SECONDS))
                            print('Protein: {0}\nProtein_data: {1}'.format(protein, protein_data))
                            print('Taxon: {0}\nTaxon_data: {1}'.format(taxon, taxon_data))
                            if blastp(protein, protein_data, taxon, taxon_data, format) == True:
                                parse_blast_result(protein, protein_data, taxon)
                                break
                        except Exception: # Replace Exception with something more specific.
                            iteration += 1
                            print('{0}. blast attempt failed. Trying again after a break of {1} seconds.'.format(iteration, iteration**1.5*30))
                            time.sleep(iteration**1.5*30)
                            continue
                else:
                    print('Checking for {0}/{1} -> ok, checking for {0}{2} -> ok.'.format(protein, BLAST_XMLFILE, BLAST_HTMLFILE))
        # Go back one directory (otherwise every new analysis will be a subdirectory in the previous diretory)
        os.chdir('..')

if __name__ == '__main__':
    run()
