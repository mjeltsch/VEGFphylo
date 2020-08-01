#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
#
# UBUNTU PACKAGES
#
# This script has been developed to run under Linux and was specifically tested ONLY under Ubuntu 18.04.
# It requires the installation of the following ubuntu packages:
# -
#
# PYTHON MODULES
#
# The following Python modules need to be installed:
# sudo pip3 install biopython
#
#
# NCBI BLAST EXECUTABLES
#
# Needs also NCBI blast executables >2.6 from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# for remote access to the blast server and overriding of Ubuntu's old blast executables (2.2) in the user's
# ~/.profile file:
#
# # set PATH so it includes user's private blast executables if they exist
# if [ -d "$HOME/bin/blast/bin" ] ; then
#    PATH="$HOME/bin/blast/bin:$PATH"
# fi
#
#
# WHAT DOES THIS SCRIPT DO?
#
# This script takes the downloaded list of species with sequenced genomes from the NCBI.
# Then it runs multiple blast searches (one for each PDGF/VEGF reference protein) against
# the non-redundant protein database. From the results
# false-positives are eliminated (using a manually assembled list).
# Then the script generates an html file, which contains a sorted table with all relevant
# sequences, annotations, and links. 

import argparse, Bio, os, csv, operator, sys, re, time, requests
from Bio import SearchIO
from Bio import SeqIO
from Bio import Entrez
from time import sleep
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from os.path import dirname
from phylolib import load_dictionary, get_phylum_from_NCBI, get_script_style_section, download_fasta_file, sanitize_species_name, get_taxon_id_from_NCBI, execute_subprocess
csv.field_size_limit(sys.maxsize)
Entrez.email = 'michael@jeltsch.org'

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
    # Open/download the list of fully sequenced genomes and return it as a dictionary
    try:
        input_csv_file = csv.DictReader(open(CSV_FILE))
    except FileNotFoundError as err:
        print(err)
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
        outfile.write('<html>\n<head>\n<title>PDGF/VEGF-like Blast hit counts for individual species</title>\n')
        outfile.write(get_script_style_section())
        outfile.write('</head>\n<body>\n')
        outfile.write('<h1>PDGF/VEGF-like Blast hist counts for individual species</h1> ')
        outfile.write('''<p>This table contains all species that fulfil the following criteria:
<ul>
<li>A BioProject (formely known as "Genome Project") with NCBI was registered with the data type "Genome sequencing" or "Genome sequencing and assembly" for a species (<a href="''')
        outfile.write(URL)
        outfile.write('''">link</a>); sequencing projects of whole genera are excluded.</li>
<li>Results have been published for the project.</li>
<li>More than 2000 protein sequences exist for the species in the non-redundant NCBI protein database or at least one PDGF/VEGF-like seuquence has been identified among the protein sequences.</li>
</ul>
</p>''')
        outfile.write('<table>')
        outfile.write('<thead><th>Phylum</th><th>Species</th><th>Blast<br>hits</th><th>protein<br>sequences</th><th>List of hit ids (link to NCBI entry)/<br>description (link to the fasta file)</th></thead>\n<tbody>\n')
        with open(CSV_FILE, 'r') as infile:
            line = infile.readline()
            while line:
                fields = line.split('\t')
                # Write first 5 columns to table row to html file
                outfile.write('<tr><td>{0}</td><td><a href="https://en.wikipedia.org/wiki/{1}" target="_blank">{1}</a></td><td>{2}</td><td>{3}</td><td>'.format(fields[0], fields[1], fields[2], fields[3]))
                # eval chokes on some long dictionaries and evaluates them as strings
                # Need to implement an own parser
                id_list = fields[5].strip("\{\}'").split("', '")
                if id_list == ['\n']:
                    outfile.write(' </td></td> ')
                    print('id_list is empty: {0}'.format(id_list))
                else:
                    print('id_list: {0}'.format(id_list))
                    linklist = ''
                    for item in id_list:
                        isolated_accno = item.split('|')[3]
                        isolated_gi = item.split('|')[1]
                        description = item.split('|')[-1].strip(": '\n}")
                        description2 = re.sub(r'\[[^)]*\]', '', description)
                        #print('add_to_linklist: <a href="https://www.ncbi.nlm.nih.gov/protein/{0}" target="_blank">{0}</a><br>'.format(isolated_accno))
                        linklist += '<span class="a1"><a href="https://www.ncbi.nlm.nih.gov/protein/{0}" target="_blank">{0}/gi{1}</a></span> <span class="b1"><a href="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id={1}" style="color:black";>{2}</a></span><br>'.format(isolated_accno, isolated_gi, description2)
                    outfile.write('<button onclick="myFunction(\'align_{0}\')">Show/Hide hits</button>'.format(fields[1].replace(' ', '_')))
                    outfile.write('<div id="align_{0}" style="display: none;">{1}</div>'.format(fields[1].replace(' ', '_'), linklist))
                outfile.write('</td></tr>\n')
                line = infile.readline()
            outfile.write('</tbody>\n</table>\n</body>\n<html>\n')

def blastp(SPECIES):
    global final_result_string
    SPECIES_SANITIZED = sanitize_species_name(SPECIES)
    # Get phylum for the species
    phylum = get_phylum_from_NCBI(SPECIES, VERBOSE=False)
    # Get the taxid for the species
    # Replace maybe later with the new get_species_taxids.sh script
    # from the blast tools
    TAXID = get_taxon_id_from_NCBI(SPECIES, VERBOSE=False)
    # Make a subfolder for each species if it does not exist
    SPECIES_DIR = '{0}{1}'.format(RESULT_DIR, SPECIES_SANITIZED)
    if not os.path.isdir(SPECIES_DIR):
        os.mkdir(SPECIES_DIR)
    # Database
    BLAST_DATABASE = '{0}/../blastdb/nr'.format(APPLICATION_PATH)
    EVALUE = 0.1
    OPTIONAL_BLAST_NO = 20000
    temp_protein_hit_dict = {}
    total_num_sequences = 0
    for protein, data in master_dictionary.items():
        PROTEIN_ID = data[0]
        outfile = '{0}{1}/{2}.xml'.format(RESULT_DIR, SPECIES_SANITIZED, PROTEIN_ID)
        # This was temporarily used to use a specific species
        #SPECIES = 'txid6239'
        if REMOTE:
            print('\nRunning remote blastp against {0}, limited to organism {1}, requesting {2} results in XML format, with the following query sequence: {3} ({4}).'.format(BLAST_DATABASE, SPECIES, OPTIONAL_BLAST_NO, protein, PROTEIN_ID))
            ENTREZ_QUERY = 'txid{0}[organism]'.format(TAXID)
            try:
                print('Start remote blast job...')
                result_handle = NCBIWWW.qblast('blastp', BLAST_DATABASE, PROTEIN_ID, hitlist_size = OPTIONAL_BLAST_NO, expect = EVALUE, entrez_query = ENTREZ_QUERY, format_type = 'XML')
                # Write blast result to xml or html file
                print('Opening file {0} for writing'.format(outfile))
                with open(outfile, "w") as out_handle:
                    out_handle.write(result_handle.read())
                    print('{0} written to disk.'.format(outfile))
                result_handle.close()
            except Exception as ex:
                print('Remote blasting failed. The error was:\n{0}'.format(ex))
            else:
                print('Remote blasting completed successfully.')
        else:
            print('Running local blast...')
            try:
                infile = '{0}/prog_data/reference_proteins/{1}.fasta'.format(APPLICATION_PATH, PROTEIN_ID)
                # The Biopython wrapper for blastp is not yet taxon-aware! Hence we need to execute a subprocess.
                #blastx_cline = 'NcbiblastpCommandline(query = {0}, db = {1}, evalue = {2}, outfmt=5, out = {3})'.format(infile, BLAST_DATABASE, EVALUE, outfile)
                comment = 'Executing local blastp search.'
                bash_command = 'blastp -db {0} -taxids {1} -query {2} -evalue {3} -out {4} -outfmt 5'.format(BLAST_DATABASE, TAXID, infile, EVALUE, outfile)
                # Uncomment the following line if you want to ignore previous results which are locally cached
                if not os.path.isfile(outfile):
                    execute_subprocess(comment, bash_command)
            except Exception as ex:
                print('Local blasting failed. The error was: \n{0}'.format(ex))
            else:
                print('Local blasting completed successfully.')
        result_handle = open(outfile, "r")
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
                        # Add hit to temporary list only if it is unique (not yet in temp_protein_hit_dict)
                        if alignment.hit_id in temp_protein_hit_dict:
                            print('This protein ({0}) was already added to the temp_protein_hit_dict.\n'.format(alignment.hit_id))
                        else:
                            isolate_id = alignment.hit_id.split('|')[1]
                            print('isolate_id: {0}'.format(isolate_id))
                            exclude_file = '{0}/prog_data/proteins-exclude/{1}/{2}.fasta'.format(APPLICATION_PATH, phylum, isolate_id)
                            print('Checking for file {0}... '.format(exclude_file), end='')
                            if not os.path.isfile(exclude_file):
                                temp_protein_hit_dict[alignment.hit_id] = alignment.hit_def
                                print('Not found, adding {0} to alignment hit list.\n'.format(alignment.hit_id))
                                FASTA_FILE = '{0}{1}/{2}.fasta'.format(RESULT_DIR, SPECIES_SANITIZED, isolate_id)
                                print('Checking for locally cached file ({0})...'.format(FASTA_FILE), end='')
                                if os.path.isfile(FASTA_FILE):
                                    print('Success! Using locally cached file.')
                                else:
                                    print('Nope. Need to fetch it from NCBI.')
                                    if download_fasta_file(isolate_id, FASTA_FILE):
                                        print('Successfully downloaded {0} from Entrez.'.format(isolate_id))
                                    else:
                                        print('Downloading {0} from Entrez failed.'.format(isolate_id))
                                    time.sleep(2)
                            else:
                                print('Found! Excluding {0} from alignment hit list.\n'.format(alignment.hit_id))
            # Method 2 (to get the number of hits)
            blastp_result = SearchIO.read(outfile, 'blast-xml')
            number_of_results = len(blastp_result)
            print('{0} {1} homologs found in a total of {2} protein sequences for {3} ({4}).'.format(number_of_results, protein, total_num_sequences, SPECIES, phylum))
        except Exception as ex:
            print('Could not parse blast-xml file {0}. File might not exist or is empty.'.format(outfile))
        result_handle.close()
    if len(temp_protein_hit_dict) > 0 or total_num_sequences >= MINIMUM:
        #print('\nAdding to CSV file:\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(phylum, SPECIES, len(temp_protein_hit_dict), total_num_sequences, int(total_num_sequences// (len(temp_protein_hit_dict)+0.00001)), temp_protein_hit_dict))
        write_to_file(FINAL_RESULTS, '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(phylum, SPECIES, len(temp_protein_hit_dict), total_num_sequences, int(total_num_sequences//(len(temp_protein_hit_dict)+0.00001)), temp_protein_hit_dict))

def run():
    global master_dictionary, RESULT_DIR, REMOTE, APPLICATION_PATH, final_result_string, FINAL_RESULTS, FINAL_RESULTS_SORTED, MINIMUM, URL
    #URL = 'https://www.ncbi.nlm.nih.gov/genomes/solr2txt.cgi?q=%5Bdisplay()%5D.from(GenomeBrowser).usingschema(%2Fschema%2FGenomeAssemblies).matching(group%3D%3D%5B%22Animals%22%5D)&fields=organism%7COrganism%20Name%2Clineage%7COrganism%20Groups%2Csize%7CSize(Mb)%2Cchromosomes%7CChromosomes%2Corganelles%7COrganelles%2Cplasmids%7CPlasmids%2Cassemblies%7CAssemblies&filename=genomes.csv&nolimit=on'
    # This is probably the better resource since it lists the number of CDS for each genome
    URL = 'https://www.ncbi.nlm.nih.gov/genomes/solr2txt.cgi?q=%5Bdisplay()%5D.from(GenomeAssemblies).usingschema(%2Fschema%2FGenomeAssemblies).matching(tab%3D%3D%5B%22Eukaryotes%22%5D%20and%20group%3D%3D%5B%22Animals%22%5D)&fields=organism%7COrganism%20Name%2Clineage%7COrganism%20Groups%2Cstrain%7CStrain%2Cbiosample%7CBioSample%2Cbioproject%7CBioProject%2Cassembly%7CAssembly%2Clevel%7CLevel%2Csize%7CSize(Mb)%2Cgc_content%7CGC%25%2Creplicons%7CReplicons%2Cwgs%7CWGS%2Cscaffolds%7CScaffolds%2Cproteins%7CCDS%2Crelease_date%7CRelease%20Date%2Cftp_path_genbank%7CGenBank%20FTP%2Cftp_path_refseq%7CRefSeq%20FTP&filename=eukaryotes.csv&nolimit=on'
    # Only record species with 0 hits if MINUMUM or more protein sequences exist for a species
    # (Species with 1 or more hists are always recorded.
    MINIMUM = 0
    # Determine directories of script (in order to load & save the data & log files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    # List of fully sequenced genomes
    CSV_FILE = '{0}/prog_data/genomes.csv'.format(APPLICATION_PATH)
    RESULT_DIR = '{0}/user_data/individual_species/'.format(APPLICATION_PATH)
    FINAL_RESULTS = '{0}/user_data/full_genome_data.csv'.format(APPLICATION_PATH)
    if os.path.isfile(FINAL_RESULTS):
        os.remove(FINAL_RESULTS)
    FINAL_RESULTS_SORTED = '{0}/user_data/full_genome_data_sorted.csv'.format(APPLICATION_PATH)
    preamble1, master_dictionary = load_dictionary('{0}/prog_data/master_dictionary.py'.format(APPLICATION_PATH), VERBOSE = False)
    input_csv_file = getCSVfile(CSV_FILE)
    # For some species there are multiple projects, deduplication is needed
    temp_species_list = []
    for line in input_csv_file:
        print('Iterating through {0}...'.format(CSV_FILE))
        species_name = line['#Organism Name']
        CDS = line['CDS']
        # Only deal with sequenceing projects that have identified or predicted coding sequences (CDS), i.e. proteins
        if int(CDS) > 0:
            # Prevent from analyzing a species more than once (if there are more than one genome projects)
            if species_name not in temp_species_list:
                temp_species_list.append(species_name)
                # Check that it is really a species and not only the genus (i.e. that the sceintific name consists of at least 2 whitespace-separated parts)
                if len(species_name.strip().split()) > 1:
                    print("Go to blast...")
                    blastp(species_name)
    print('Writing of {0} complete.'.format(FINAL_RESULTS))
    sort_csv(FINAL_RESULTS)
    print('Sorting of {0} complete.'.format(FINAL_RESULTS))
    convert_to_html(FINAL_RESULTS_SORTED)

if __name__ == '__main__':
    run()
