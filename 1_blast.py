#!/usr/bin/python3.6
# -*- coding: utf-8 -*-
#
#
# UBUNTU PACKAGES
#
# This script has been developed to run under Linux and was specifically tested ONLY under Ubuntu 16.04 and 18.04.
# It requires the installation of the following ubuntu packages:
#
# New list for ubuntu packages: emacs, python3-pip, 
# New list for python packages via pip3: biopython, ete3 
#
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
from phylolib import execution_time_str, load_dictionary, load_blacklist, execute_subprocess, get_taxon_id_from_NCBI, blast_formatter

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

if args.local_remote == 'local_recycle':
    REMOTE = 'local_recycle'
    print('This run uses - if available - blastp search results from previous runs to save time. New searchers are done locally.')
elif args.local_remote == 'remote_recycle':
    REMOTE = 'remote_recycle'
    print('This run uses - if available - blastp search results from previous runs to save time. New searches are done remotely.')
elif args.local_remote == 'local':
    REMOTE = 'local'
    print('This run does execute a new local blastp search for every protein. This will take a long time!')
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
    print('Parsing...')
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

def blastp(PROTEIN, PROTEIN_DATA, TAXON, TAXON_DATA):
    global BLAST_XMLFILE, BLAST_HTMLFILE
    # Database
    BLAST_DATABASE = 'nr'
    EVALUE = 0.1
    PROTEIN_ID, SUBRANGE, OPTIONAL_BLAST_NO, SYNONYMS, RELATED_PROTEINS = PROTEIN_DATA
    # Memorize start time to figure out how long the whole script execution takes
    start_time = time.time()
    #print('Types:\nTAXON: {0}\nTAXON_DATA: {1}\nPROTEIN: {2}\nTAXON_DATA[2]: {3}\n'.format(type(TAXON), type(TAXON_DATA), type(PROTEIN), type(TAXON_DATA[2])))
    print('\nRunning remote blastp against {0}, subsection {1}, requesting {2} results in html format, with the following query sequence: {3} ({4}).'.format(BLAST_DATABASE, TAXON, OPTIONAL_BLAST_NO, PROTEIN, PROTEIN_ID))
    # The hitlist_size seems to be ignored by the wrapper/service????
    # Reference: http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
    # The default is to use the sequence_id for the blasting. However, if no sequence_id is available, we can provide the amino acid sequence locally
    # Check here whether the sequence_id is and entry in the local sequence_disctionary and get the sequence for blasting
    #result_handle = NCBIWWW.qblast("blastp", BLAST_DATABASE, PROTEIN_ID[0], hitlist_size = str(OPTIONAL_BLAST_NO), expect = EVALUE, query_from = ALIGNMENT_TRIMMING[0], query_to = ALIGNMENT_TRIMMING[1], entrez_query='NOT '+EXCLUDE+'[organism]')
    # Alternatively to the taxon name, the taxon id can be specified like e.g. "txid9606"[organism]
    # When selecting the format ("Format_type"), please remember that HTML cannot be easily parsed. Use XML instead!
    ENTREZ_QUERY = 'txid{0}[organism]'.format(TAXON_DATA[0])
    #print('ENTREZ_QUERY = {0}'.format(ENTREZ_QUERY))
    if SUBRANGE == None:
        STARTPOS = ''
        ENDPOS = ''
    else:
        STARTPOS = SUBRANGE.split('-')[0]
        ENDPOS = SUBRANGE.split('-')[1]
    #try:
    if REMOTE == 'remote' or REMOTE == 'remote_recycle':
        print('Starting html request at {0} ... '.format(str(datetime.datetime.now())[:-7]), end='')
        result_handle = NCBIWWW.qblast('blastp', BLAST_DATABASE, PROTEIN_ID, hitlist_size = OPTIONAL_BLAST_NO, expect = EVALUE, entrez_query = ENTREZ_QUERY, format_type = 'HTML', query_from = STARTPOS, query_to = ENDPOS)
        print('Finished html request.')
        # Write blast results to an html file
        with open(BLAST_HTMLFILE, "w") as out_handle:
            html_source = result_handle.read()
            out_handle.write(html_source)
            print('{0} ({1}) written to file {2}.'.format(PROTEIN, PROTEIN_ID, BLAST_HTMLFILE))
        result_handle.close()
    elif REMOTE == 'local' or REMOTE == 'local_recycle':


        comment = 'Executing local blastp search.'
        bash_command = 'blastp -db {0} -taxids {1} -query {2} -evalue {3} -html > {4}'.format(BLAST_DATABASE, TAXID, infile, EVALUE, outfile)
        execute_subprocess(comment, bash_command)


    # except Exception as ex:
    #     print('Blasting failed. The error was: {0}, {1}'.format(ex, error))
    #     how_long = execution_time_str(time.time()-start_time)
    #     print('Blasting took {0}'.format(how_long))
    #     return False
    # else:
    #     print('Blasting succeeded. The blast output was: {0}'.format(output))
    #     how_long = execution_time_str(time.time()-start_time)
    #     print('Blasting took {0}'.format(how_long))
    #     return True
    return True

def extract_RID_from_html_file(FILE):
    # Get the RID from the html file
    with open(FILE, 'r') as handle:
        html_source = handle.read()
        RID = re.search('</span> <h1>results for RID-(.*)</h1></span>', html_source).group(1)
        print('The extracted RID is: {0}'.format(RID))
    return RID

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

    # Determine directories of script (in order to load & save the data & log files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    DATA_DIR = '{0}/user_data/primary_blast_results'.format(APPLICATION_PATH)
    LOGFILE = '{0}/prog_data/temp/logfile.txt'.format(APPLICATION_PATH)
    SUMMARY_FILE = '{0}/summary.csv'.format(DATA_DIR)
    print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    create_subdirectory(DATA_DIR)
    # Loading data file
    preamble1, master_dictionary = load_dictionary('{0}/prog_data/master_dictionary.py'.format(APPLICATION_PATH), VERBOSE = False)
    preamble2, taxon_dictionary = load_dictionary('{0}/prog_data/taxon_data.py'.format(APPLICATION_PATH), VERBOSE = False)
    # If commandline argument 'remote' is given, rename old results and execute new pblast requests
    # instead of recycling existing ones
    if REMOTE == 'local_recycle' or REMOTE == 'remote_recycle':
        print('Using .xml blast result files from directory {0} (skipping remote blasts whenever possible).'.format(DATA_DIR))
    elif REMOTE == 'remote' or REMOTE == 'local':
        # Rename old results and create a new result directory
        OLD_DATA_DIR = '{0}_{1}'.format(DATA_DIR, datetime.datetime.now().strftime("%y%m%d_%H%M%S"))
        os.rename(DATA_DIR, OLD_DATA_DIR)
        print('\nMoving old data to {0}.\n'.format(OLD_DATA_DIR))
    # Create directory for the downloaded and generated data
    if not os.path.exists(DATA_DIR):
        os.mkdir(DATA_DIR)
    # Direct all terminal output additionally into a logfile
    sys.stdout = logfile(LOGFILE)
    # Change the current working directory
    os.chdir(DATA_DIR)
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
                print('Analyzing protein {0} in taxon {1}:'.format(protein, taxon))
                BLAST_HTMLFILE = 'blast_results_{0}.html'.format(taxon)
                BLAST_XMLFILE = 'blast_results_{0}.xml'.format(taxon)
                # If no local results are found, loop as long until the blasting succeeds or for maximally MAX_ITER tries
                iteration = 0
                MAX_ITER = 10
                while True and iteration < MAX_ITER:
                    try:
                        # Check whether we have HTML file and it is > 0
                        if os.path.isfile(BLAST_HTMLFILE) and os.path.getsize(BLAST_HTMLFILE) > 0:
                            print('{0} exists and is > 0'.format(BLAST_HTMLFILE))
                            # Check whether we have an XML file and it is > 0
                            if os.path.isfile(BLAST_XMLFILE) and os.path.getsize(BLAST_XMLFILE) > 0:
                                print('{0} exists and is > 0'.format(BLAST_XMLFILE))
                                parse_blast_result(protein, protein_data, taxon)
                                break
                            else:
                                # Check whether the HMTL file is younger than 24 hours
                                age = round(time.time() - os.path.getmtime(BLAST_HTMLFILE))
                                if age < 86400:
                                    print('{0} is younger than one day ({1} seconds old)'.format(BLAST_HTMLFILE, age))
                                    RID = extract_RID_from_html_file(BLAST_HTMLFILE)
                                    if RID != '':
                                        print('RID: {0}'.format(RID))
                                        result = blast_formatter(RID, '5', BLAST_XMLFILE, False)
                                        print('blast_fomatter result:\n{0}'.format(result))
                                        if os.path.getsize(BLAST_XMLFILE) > 0:
                                            print('XML:\n{0}'.format(result))
                                            parse_blast_result(protein, protein_data, taxon)
                                            break
                                        else: raise Exception('BLAST_XMLFILE file has size 0.')
                                    else: raise Exception('RID could not be extracted from BLAST_HTMLFILE')
                                else:
                                    if os.path.isfile(BLAST_HTMLFILE):
                                        os.remove(BLAST_HTMLFILE)
                                    raise Exception('Need to get a new BLAST_HTMLFILE...')
                        else:
                            TIME_MULTIPLICATOR = 2
                            SECONDS = int(round(taxon_data[3]**(1/9)*1000*TIME_MULTIPLICATOR, 0))
                            print('{0} does not exist or has size 0. Initiating blast. If no results are received within {1} seconds, the blast will be cancelled'.format(BLAST_HTMLFILE, SECONDS))
                            if blastp(protein, protein_data, taxon, taxon_data) == True:
                                print('Blast was successful, trying to extract RID...')
                                RID = extract_RID_from_html_file(BLAST_HTMLFILE)
                                if RID != '':
                                    print('RID: {0}'.format(RID))
                                    result = blast_formatter(RID, '5', BLAST_XMLFILE, True)
                                    print('blast_fomatter result:\n{0}'.format(result))
                                    #if result != '' and result != 'Network error: Error fetching sequence data from BLAST databases at NCBI, please try again later':
                                    if os.path.getsize(BLAST_XMLFILE) > 0:
                                        #print('XML:\n{0}'.format(result))
                                        parse_blast_result(protein, protein_data, taxon)
                                        break
                                    else:
                                        raise Exception('XML file has size 0!')
                                else: raise Exception('RID could not be extracted from {0}.'.format(BLAST_HTMLFILE))
                            else:
                                raise Exception('Blast was not successful.')
                    except Exception as error: # Replace Exception with something more specific.
                        iteration += 1
                        print('{0}. round failed: {1} \nTrying again after a break of {2} seconds.'.format(iteration, repr(error), iteration**1.5*30))
                        time.sleep(iteration**1.5*30)
                        continue
        # Go back one directory (otherwise every new analysis will be a subdirectory in the previous diretory)
        os.chdir('..')

if __name__ == '__main__':
    run()
