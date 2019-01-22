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

import argparse, subprocess, Bio, os, sys, shutil, re, time, datetime, socket, math
#from Bio.Blast import NCBIWWW
#from Bio.Blast import NCBIXML
#from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio import Entrez
from Bio import SeqIO
from Bio import Phylo
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, SequenceFace, ImgFace, SVGFace, faces, add_face_to_node
from os.path import basename, dirname, splitext, split
# To print some terminal output in color
from termcolor import colored, cprint

def execute_subprocess(comment, bash_command, working_directory='.'):
    print("\n" + comment, bash_command)
    process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=working_directory)
    output, error = process.communicate()
    process_status = process.wait()
    if output.decode('utf-8') != '':
        print("Output: " + str(output))
    if error.decode('UTF-8') != '':
        print("Error: " + str(error))

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

def drawtree(TREEFILE):
    color_dict = {  '-':"White",
                    '0':"#FF0000",
                    '1':"#FF1111",
                    '2':"#FF4444",
                    '3':"#FF8888",
                    '4':"#FFBBBB",
                    '5':"#FFDDDD",
                    '6':"#FFEEEE",
                    '7':"#FFFFFF",
                    '8':"#FFFFFF",
                    '9':"#FFFFFF",
                    '10':"#FFFFFF"}
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_data.py'.format(APPLICATION_PATH)
    preamble, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)

    print('\nLoading tree file {0}:\n'.format(TREEFILE))
    # Display tree file
    with open(TREEFILE, "r") as file:
        treestring = file.read()
        print(treestring)

    t = Tree(TREEFILE)
    ts = TreeStyle()
    ts.show_leaf_name = False
    # Zoom in x-axis direction
    ts.scale = 40
    # This makes all branches the same length!!!!!!!
    ts.force_topology = True
    #Tree.render(t, "final_tree_decoded.svg")
    #t.set_outgroup(t & OUTGROUP[0])
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.draw_guiding_lines = True
    ts.branch_vertical_margin = 10 # 10 pixels between adjacent branches
    print(t)

    # Define node styles for different animal classes
    # Mammals
    mammalia = NodeStyle()
    mammalia["bgcolor"] = "Chocolate"
    # Reptiles
    reptilia = NodeStyle()
    reptilia["bgcolor"] = "Olive"
    # Cartilaginous fish
    chondrichthyes = NodeStyle()
    chondrichthyes["bgcolor"] = "SteelBlue"
    # Ray-fenned fish
    actinopterygii = NodeStyle()
    actinopterygii["bgcolor"] = "CornflowerBlue"
    # Lobe-finned fish
    sarcopterygii = NodeStyle()
    sarcopterygii["bgcolor"] = "DarkCyan"
    # Birds
    aves = NodeStyle()
    aves["bgcolor"] = "DarkSalmon"
    # Amphibia
    amphibia = NodeStyle()
    amphibia["bgcolor"] = "DarkSeaGreen"
    # Myxini
    myxini = NodeStyle()
    myxini["bgcolor"] = "LightBlue"

    general_leaf_style = NodeStyle()
    # size of the blue ball
    general_leaf_style["size"] = 15

    # Draws nodes as small red spheres of diameter equal to 10 pixels
    nstyle = NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 15
    nstyle["fgcolor"] = "darkred"
    # Gray dashed branch lines
    #nstyle["hz_line_type"] = 1
    #nstyle["hz_line_color"] = "#cccccc"

    # Applies the same static style to all nodes in the tree if they are not leaves.
    # Note that if "nstyle" is modified, changes will affect to all nodes
    # Apply a separate style to all leaf nodes
    for node in t.traverse():
        if node.is_leaf():
            print("Setting leaf style for node " + node.name)
            animal_class_name = ''
            print(animal_class_name)
            if animal_class_name == 'Mammalia':
                node.set_style(mammalia)
            elif animal_class_name == 'Reptilia':
                node.set_style(reptilia)
            elif animal_class_name == 'Chondrichthyes':
                node.set_style(chondrichthyes)
            elif animal_class_name == 'Actinopterygii':
                node.set_style(actinopterygii)
            elif animal_class_name == 'Sarcopterygii':
                node.set_style(sarcopterygii)
            elif animal_class_name == 'Aves':
                node.set_style(aves)
            elif animal_class_name == 'Amphibia':
                node.set_style(amphibia)
            elif animal_class_name == 'Myxini':
                node.set_style(myxini)
            else:
                node.set_style(nstyle)
            # Set general leaf attributes
            node.set_style(general_leaf_style)
        else:
            node.set_style(nstyle)

    # ADD TEXT
    # for key, value in sequence_dictionary.items():
    #     # For multiple sequences of a single species an error is generated
    #     # Check first if the key is present in the tree file!
    #     if key in treestring:
    #         if key == outgroup_name:
    #             print(key, value[0], value[1], value[2], value[3], value[4])
    #             textFace = TextFace(value[1] + " (" + value[2] + ", " + key +  ", " + value[4] + ")   ", fsize = 16)
    #             (t & key).add_face(textFace, 2, "aligned")
    #             textFace.margin_left = 10
    #             print(key, value[0], value[1], value[2], value[3])
    #             textFace = TextFace(" ", fsize = 16)
    #             (t & key).add_face(textFace, 2, "aligned")
    #             textFace.margin_left = 10
    #             print(key, value[0], value[1], value[2], value[3])
    #             textFace = TextFace("Consensus score", fsize = 16)
    #             (t & key).add_face(textFace, 2, "aligned")
    #             textFace.margin_left = 10
    #         else:
    #             print(key, value[0], value[1], value[2], value[3], value[4], end = '')
    #             textFace = TextFace(value[1] + " (" + value[2] + ", " + key +  ", " + value[4] + ")   ", fsize = 16)
    #             (t & key).add_face(textFace, 2, "aligned")
    #             textFace.margin_left = 10
    #         print(" Adding...")
    #     else:
    #         print("key " + key + " is not present in the tree file. Skipping...")

    # ADD SVG IMAGES
    for node in t.traverse():
        if node.is_leaf():
            animal_class_name = node.name
            # Get the common name for this animal class if it exists
            #print('taxon_dictionary: {0}'.format(str(taxon_dictionary)))
            try:
                number_of_fully_sequenced_genomes = taxon_dictionary[animal_class_name][5]
                number_of_animal_species = taxon_dictionary[animal_class_name][1]
                number_of_protein_sequences = taxon_dictionary[animal_class_name][3]

                # Heuristic number for reliability
                reliability = number_of_fully_sequenced_genomes**2/number_of_animal_species*number_of_protein_sequences+1
                reliability = math.log10(reliability)*1.25
                reliability = str(int(round(reliability, 0)))
                #print('reliability for {0} is {1}'.format(animal_class_name, reliability))

                if number_of_protein_sequences < 1000:
                    number_of_protein_sequences = ' ' + str(number_of_protein_sequences)
                elif number_of_protein_sequences < 1000000:
                    number_of_protein_sequences = ' ' + str(int(round(number_of_protein_sequences/1000, 0))) + 'k'
                elif number_of_protein_sequences < 1000000000:
                    number_of_protein_sequences = ' ' + str(int(round(number_of_protein_sequences/1000000, 0))) + 'M'
                animal_class_name_common = taxon_dictionary[animal_class_name][2]

                # Add a "fake" header (is actually part of the first phylum row)
                if animal_class_name == 'ctenophora':
                    textFace = TextFace('# animal\n species', fsize = 16)
                    (t & animal_class_name).add_face(textFace, 0, "aligned")
                    textFace = TextFace(' # se-\n quences', fsize = 16)
                    (t & animal_class_name).add_face(textFace, 1, "aligned")
                    textFace = TextFace(' # compl.\n genomes', fsize = 16)
                    (t & animal_class_name).add_face(textFace, 2, "aligned")
                    textFace = TextFace(' ', fsize = 16)
                    (t & animal_class_name).add_face(textFace, 3, "aligned")
                    textFace = TextFace(' ', fsize = 16)
                    (t & animal_class_name).add_face(textFace, 4, "aligned")
                    protein_dict = taxon_dictionary[animal_class_name][4]
                    #(t & animal_class_name).add_face(textFace, 13, "aligned")
                    #textFace = TextFace(' reliability\n (1-10)', fsize = 16)
                    i = 5
                    for protein, value in protein_dict.items():
                        textFace = TextFace(' '+protein, fsize = 16)
                        (t & animal_class_name).add_face(textFace, i, "aligned")
                        i += 1

                # Number of animal species in this class in the NCBI protein sequence database
                textFace = TextFace(' ' + str(number_of_animal_species), fsize = 16)
                (t & animal_class_name).add_face(textFace, 0, "aligned")
                textFace = TextFace(' ' + str(number_of_protein_sequences), fsize = 16)
                (t & animal_class_name).add_face(textFace, 1, "aligned")
                textFace = TextFace(' ' + str(number_of_fully_sequenced_genomes), fsize = 16)
                (t & animal_class_name).add_face(textFace, 2, "aligned")

                # TEXT
                print('Adding text for {0}'.format(animal_class_name))
                # Do not display anything if there is no common animal class name
                if animal_class_name_common == '?' or animal_class_name_common == '':
                    description = '{0}'.format(animal_class_name)
                else:
                    description = '{0} ({1})'.format(animal_class_name, animal_class_name_common)
                textFace = TextFace(description, fsize = 16)
                (t & animal_class_name).add_face(textFace, 4, "aligned")
                # Reliability measure
                #textFace = TextFace(reliability, fsize = 16)
                #(t & animal_class_name).add_face(textFace, 13, "aligned")
                #textFace.margin_left = 10

                # PROTEIN HITS
                protein_dict = taxon_dictionary[animal_class_name][4]
                i = 5
                for protein, value in protein_dict.items():
                    textFace = TextFace(' '+str(value), fsize = 16)
                    (t & animal_class_name).add_face(textFace, i, "aligned")
                    textFace.background.color = color_dict[reliability]
                    i += 1

                # IMAGE
                print("Adding image")
                svgFace = SVGFace('{0}{1}.svg'.format(IMG_BASENAME, animal_class_name), height = 40)
                (t & animal_class_name).add_face(svgFace, 3, "aligned")
                svgFace.margin_right = 10
                svgFace.margin_left = 10
                svgFace.hzalign = 2
                print('Adding svg image for {0}.'.format(animal_class_name))
            except Exception as ex:
                print('Could not get common animal class name for {0} dictionary from file {1}. Error: '.format(animal_class_name, TAXON_DICTIONARY_FILE) + str(ex))

    # ROTATING SOME NODES CAN BE DONE HERE:
    # MOVE ACTINOPTERYGII NEXT TO THE OTHER FISH

    # spotted gar and coelacant
    #n1 = t.get_common_ancestor("XP_006632034.2", "XP_006006690.1")
    #n1.swap_children()

    # Add description to treefile
    description_text = "Analysis performed " + datetime.datetime.now().strftime("%y%m%d_%H%M%S") + "\n"
    # Add stuff
    description_text += " ";

    # Add legend
    textFaceLegend = TextFace(description_text, fsize = 11)
    ts.legend.add_face(textFaceLegend, column = 0)
    t.render(SVGFILE, tree_style = ts, units = "mm", w = 300)
    print('Drawing tree completed.')

def run():
    global APPLICATION_PATH, TREEFILE, SVGFILE, IMG_BASENAME, LOGFILE

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
    TREEFILE = '{0}/data/animalia.newick'.format(APPLICATION_PATH)
    SVGFILE = '{0}/animalia.svg'.format(APPLICATION_PATH)
    IMG_BASENAME = '{0}/images/'.format(APPLICATION_PATH)
    LOGFILE = 'logfile.txt'

    if os.path.isfile(TREEFILE):
        drawtree(TREEFILE)
        command = 'inkscape --export-pdf {0} {1}'.format('animalia.pdf', SVGFILE)
        execute_subprocess("Generating PDF file with the following command:\n", command)
    else:
        print('Not drawing tree since tree file {0} does not exist.'.format(TREEFILE))

if __name__ == '__main__':
    run()
