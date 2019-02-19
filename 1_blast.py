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

parser = argparse.ArgumentParser()
parser.add_argument("local_remote", help = "local, remote or test")
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
    print("Remote blast:")
elif args.local_remote == 'remote':
    REMOTE = 'remote'
    print("Local blast:")
else:
    REMOTE = 'test'
    print("Simulated blast:")

# Converts the number of seconds into a days/hours/minutes/seconds string
def execution_time_str(elapsed_time_seconds):
    min, sec = divmod(elapsed_time_seconds, 60)
    hours, min = divmod(min, 60)
    days, hours = divmod(hours, 24)

    return (str(days) + " days " if days != 0 else '') + (str(hours) + " hours " if hours != 0 else '') + str(min)[:-2] + " min " + str(round(sec, 1)) + " sec"

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

def load_dictionary(FILENAME):
    # Load a sequence_dictionary if it exists
    if os.path.isfile(FILENAME):
        try:
            dictionary = read_file_to_dict(FILENAME)
            print("\nReading in the dictionary " + FILENAME + ":\n")
            for key, value in dictionary.items():
                print(key, value)
        except Exception as e:
            dictionary = {}
    else:
        dictionary = {}
    return dictionary

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

def parse_blast_result(PROTEIN, LIST, TAXON):
    PROTEIN_ID, SUBRANGE, OPTIONAL_BLAST_NO, ALIASES, FALSE_FRIENDS = LIST
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

def blastp(PROTEIN, LIST, TAXON, TAXON_DATA):
    # Database
    BLAST_DATABASE = 'nr'
    EVALUE = 0.1
    PROTEIN_ID, SUBRANGE, OPTIONAL_BLAST_NO, ALIASES, FALSE_FRIENDS = LIST
    # Memorize start time to figure out how long the whole script execution takes
    start_time = time.time()
    if REMOTE == 'remote':
        print('\nRunning remote blastp against {0} (subsection {1} ({2}), requesting {3} results) with the following query sequence: {5}. Blast job started at {6}'.format(BLAST_DATABASE, TAXON, TAXON_DATA[2], str(OPTIONAL_BLAST_NO), PROTEIN, str(datetime.datetime.now())[:-7]))
        try:
            # The hitlist_size seems to be ignored by the wrapper/service????
            # Reference: http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
            # The default is to use the sequence_id for the blasting. However, if no sequence_id is available, we can provide the amino acid sequence locally
            # Check here whether the sequence_id is and entry in the local sequence_disctionary and get the sequence for blasting
            #result_handle = NCBIWWW.qblast("blastp", BLAST_DATABASE, PROTEIN_ID[0], hitlist_size = str(OPTIONAL_BLAST_NO), expect = EVALUE, query_from = ALIGNMENT_TRIMMING[0], query_to = ALIGNMENT_TRIMMING[1], entrez_query='NOT '+EXCLUDE+'[organism]')
            # Alternatively to the taxon name, the taxon id can be specified like e.g. "txid9606"[organism]
            # When selecting the format ("Format_type"), please remember that HTML cannot be easily parsed. Use XML instead!
            ENTREZ_QUERY = '\"' + TAXON +'\"[organism]'
            if SUBRANGE != None:
                SUBRANGE = ', query_loc = ' + SUBRANGE
            else:
                SUBRANGE = ''
            result_handle = NCBIWWW.qblast("blastp", BLAST_DATABASE, PROTEIN_ID, hitlist_size = str(OPTIONAL_BLAST_NO), expect = EVALUE, entrez_query = ENTREZ_QUERY, format_type = "XML"SUBRANGE)
            # Write blast result to xml file
            with open(BLAST_XMLFILE, "w") as out_handle:
                out_handle.write(result_handle.read())
            # Repeat for HTML output (how else to do this easily???)
            # Converting XML to html?
            # xsltproc --novalid blast2html.xsl blast.xml
            result_handle = NCBIWWW.qblast("blastp", BLAST_DATABASE, PROTEIN_ID, hitlist_size = str(OPTIONAL_BLAST_NO), expect = EVALUE, entrez_query = ENTREZ_QUERY, format_type = "HTML"SUBRANGE)
            # Write blast result to HTML file
            with open(BLAST_HTMLFILE, "w") as out_handle:
                out_handle.write(result_handle.read())
                return True
            result_handle.close()
        except Exception as ex:
            print("Something went wrong with the blasting. Perhaps the Blast server did not respond? More about the error: " + str(ex))
            return False
    elif REMOTE == 'local':
        # The local blast branch is not kept up-to-date with the remote blast branch!!!! The protein sequences in data/proteins are only needed for the local blast search. 
        START_SEQUENCE = '{0}/data/proteins/{1}.fasta'.format(APPLICATION_PATH, PROTEIN)
        PHYLUMFILE = '{0}/data/gi_lists/{1}.gi'.format(APPLICATION_PATH, TAXON)
        blastp_cline = NcbiblastpCommandline(cmd='blastp', query=START_SEQUENCE, db=BLAST_DATABASE, gilist=PHYLUMFILE, evalue=EVALUE, remote=False, outfmt=5, max_target_seqs=OPTIONAL_BLAST_NO, out=BLAST_XMLFILE, num_threads=4)
        print("\nRunning local blastp with the following query:")
        print(blastp_cline)
        print("\nPlease be patient. Retrieving XML file. ", end ='')
        try:
            stdout, stderr = blastp_cline()
        except Exception as ex:
            print("\nSomething went wrong with the local blastp. Most likely the local standalone blast server did not respond. More about the error:\n" + str(ex))
        # Repeat for HTML output (how else to do this easily???)
        blastp_cline = NcbiblastpCommandline(cmd='blastp', query=START_SEQUENCE, db=BLAST_DATABASE, gilist=PHYLUMFILE, evalue=EVALUE, remote=False, html=True, num_descriptions=OPTIONAL_BLAST_NO, num_alignments=OPTIONAL_BLAST_NO, out=BLAST_HTMLFILE, num_threads=4)
        print("Retrieving HTML file.")
        try:
            stdout, stderr = blastp_cline()
            return True
        except Exception as ex:
            print("\nSomething went wrong with the local blastp. Most likely the local standalone blast server did not respond. More about the error:\n" + str(ex))
            return False
    else:
        # Use blast result files from the test directory
        pass
    how_long = execution_time_str(time.time()-start_time)
    print('Blasting completed in {0}.\n'.format(how_long))

# Wrapping blastp in a timeout function
blastp = timeout()(blastp)

def run():
    blacklist = ['ctenophora', 'porifera', 'placozoa', 'xenacoelomorpha', 'cyclostomata', 'onychophora', 'pycnogonida', 'myriapoda', 'nematomorpha', 'loricifera', 'kinorhyncha', 'chaetognatha', 'bryozoa', 'entoprocta', 'cycliophora', 'nemertea', 'phoroniformea', 'gastrotricha', 'platyhelminthes', 'gnathostomulida', 'micrognathozoa', 'orthonectida', 'dicyemida']
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
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    # Loading data file
    master_dictionary = load_dictionary('{0}/data/master_dictionary.py'.format(APPLICATION_PATH))
    taxon_dictionary = load_dictionary('{0}/data/taxon_data.py'.format(APPLICATION_PATH))
    LOGFILE = 'logfile.txt'
    if REMOTE == 'test':
        DATA_DIR = 'test'
        print('Using data from directory {0}/{1} (skipping local/remote blast).'.format(os.getcwd(), DATA_DIR))
    else:
        DATA_DIR = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
        print('\nCreating new directory {0} to store results.\n'.format(DATA_DIR))
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
        if protein_data[8] == None:
            protein_data[8] = 50
        # Create a subdirectory named according to the protein and change cwd into it
        create_subdirectory(protein)
        os.chdir(protein)
        for taxon, taxon_data in taxon_dictionary.items():
            if taxon not in blacklist:
                print('Analyzing {0}:'.format(taxon), end='')
                BLAST_XMLFILE = 'blast_results_{0}.xml'.format(taxon)
                BLAST_HTMLFILE = 'blast_results_{0}.html'.format(taxon)
                # Loop as long as the blasting succeeds
                while True:
                    try:
                        SECONDS = int(round(taxon_data[3]**(1/9)*600, 0))
                        #SECONDS = 30
                        print(' Timeout = {0} seconds.'.format(SECONDS))
                        print('protein: {0}\nprotein_data: {1}\ntaxon: {2}\ntaxon_data: {3}'.format(protein, protein_data, taxon, taxon_data))
                        if blastp(protein, protein_data, taxon, taxon_data) == True:
                            print('ready')
                            parse_blast_result(protein, protein_data, taxon)
                            print('Blasting succeeded.')
                            break
                    except Exception: # Replace Exception with something more specific.
                        print('Blasting failed. Trying again after a break...')
                        time.sleep(random.randint(120, 240))
                        continue
                # Wait between taxa not to upset the server
                if REMOTE == 'test':
                    pass
                else:
                    time.sleep(random.randint(5, 15))
        # Go back one directory (otherwise every new analysis will be a subdirectory in the previous diretory)
        os.chdir('..')

if __name__ == '__main__':
    run()
