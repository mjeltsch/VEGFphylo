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
from phylolib import load_blacklist, load_dictionary, insert_line_breaks, write_dict_to_file, read_file_to_dict, execution_time_str

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

def get_species_number(taxon):
    print('Retrieving species number for {0}... -> '.format(taxon), end='')
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[subtree]+AND+specified[prop]+NOT+subspecies[rank]'.format(taxon)
    r = requests.get(URL)
    text = r.text
    #print(text)
    return int(text.split("Count>")[1][:-2])

def get_taxon_id_and_phylum(species_name):
    # Check whether the species is in the local sqlite database
    try:
        print('Getting id & phylum for \'{0}\', checking local sqlite database first.\n'.format(species_name), end = '')
        data = db_retrieve_species(conn, species_name, True)
        print('\nDATABASE RESULT: {0}\n'.format(data))
        # Check that the database answer is a tuple with 3 elements (fetchone gives a tuple, fetchall a list of tuples)
        if isinstance(data, tuple) and len(data) == 3:
            TAXON_ID = data[1]
            PHYLUM = data[2]
        else:
            raise Exception('{0} not in local sqlite database.'.format(species_name))
    # If the species is not in the local dictionary, get the data from NCBI
    except Exception as ex:
        print("Species not in local sqlite database. Getting data from NCBI... ")
        TAXON_ID = get_taxon_id_from_NCBI(species_name)
        print('=> {0}'.format(TAXON_ID))
        # Identify to which phylum the species belongs. For this, we need to download the full taxonomy tree for the species (which
        # is impossible via the API). Hence, we retrieve the full web page for the taxon and look whether any of the phyla from our
        # TAXON_DICTIONARY_FILE are part of the full taxonomy tree.
        PHYLUM = get_phylum_from_NCBI(TAXON_ID)
        print('=> {0}'.format(PHYLUM))
        if db_insert_species(conn, species_name, TAXON_ID, PHYLUM) != 0:
            print('New species {0} inserted into database.'.format(species_name))
    return TAXON_ID, PHYLUM

# How many sequences are in the local database for this taxon? The gi files were manually downloaded from
# NCBI. There should be smarted way to get these!
def get_sequence_number(taxon):
    print('Retrieving number of sequences for taxon {0}... -> '.format(taxon), end='')
    SEQUENCE_GI_LIST = '{0}/data/gi_lists/{1}.gi'.format(APPLICATION_PATH, taxon)
    with open(SEQUENCE_GI_LIST) as file:
        for i, l in enumerate(file):
            pass
    return i + 1

def get_taxon_id_from_NCBI(species_name, VERBOSE=True):
    # Counter to increase waiting period upon incresing number of unsuccesful trials
    i = 0
    while True:
        try:
            if VERBOSE: print('Retrieving taxon id for {0} '.format(species_name), end='')
            URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={0}[Scientific Name]'.format(species_name)
            r = requests.get(URL, timeout=20)
            # Extract taxon_id
            list = r.text.split("Id>")
            TAXON_ID = list[1][:-2]
            # Sleep for online requests (otherwise you'll be blocked after a while!)
            time.sleep(2)
        except Exception as ex:
            i += 6
            if VERBOSE: print('Unexpected server response:\n{0}\n\nError: {1}\n\nTrying again...\n'.format(text, str(ex)))
            # Wait if server returns an error (mostly because you are blocked or the network is down)
            time.sleep(i**2)
            continue
        # If there were no exceptions (= we got the taxon_id), break out of the while loop
        else:
            break
    return TAXON_ID

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
                if VERBOSE: print('Looking for TAXON_ID {0}'.format(taxon_data[0]))
                if '<TaxId>{0}</TaxId>'.format(taxon_data[0]) in text:
                    if VERBOSE: print('Adding phylum info {0} to taxon_id {1}'.format(taxon, TAXON_ID))
                    phylum = taxon
                    break
                i += 1
            print(phylum)
            # The following line will be only executed if calling 'phylum' does not give a NameError
            # When i = 51, none of the 51 phyla was in the text that was received from NCBI.
            if i == 51:
                if VERBOSE: print('{0} does not belong to any phylum in the phylum list'.format(species_name))
            elif i < 51:
                i += 1
                if VERBOSE: print('{0} phyla (of 51) checked before hit was found.'.format(i))
            else:
                if VERBOSE: print('Unknown error when trying to identify phylum for {0}.'.format(species_name))
        except NameError as err:
            if VERBOSE: print("Sleeping due to server error. Will retry in 60 seconds.")
            time.sleep(60)
            continue
        # If there were no exceptions (= we got the phylum), break out of the while loop
        else:
            break
    return phylum

def write_to_html_scrutinize_file(list_to_scrutinize):
    # This sorts the list of lists according to the third element of each list (= type of hit; synonym, related protein, unknown)
    list_to_scrutinize.sort(key=lambda x: x[2])
    preamble = '''<html>
<head>
<title>To scrutinize</title>
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
    print('Writing to HTML file...')
    #for item in sorted_list:
    new_type = ''
    i = 1
    j = 1
    html_text = ''
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
        #                                                                                                                    taxon    protein  type_of  gid      prot_id
        #                                                                                                                                      _hit
        html_text += '<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td><td>{4}</td><td>{5}</td><td>{6}</td></tr>\n'.format(item[0], item[1], item[2], item[3], item[4], link_to_protein, link_to_blast)
        new_type = item[2]
        i += 1
    #lineitem += '\n</body>\n</html>\n'
    if not os.path.isfile(HTML_SCRUTINIZE_FILE):
        with open(HTML_SCRUTINIZE_FILE, 'w') as file:
            file.write(preamble)
    with open(HTML_SCRUTINIZE_FILE, 'a') as file:
        file.write(html_text)

def get_protein_data(taxon):
    global TOTAL_COUNT, TOTAL_UNKNOWN_COUNT
    if args.directory == '':
        return []
    else:
        new_protein_data = {}
        list_to_scrutinize = []
        print('\n---------------')
        print('Analysing taxon {0}'.format(taxon).upper())
        print('---------------\n')
        for protein, protein_data in master_dictionary.items():
            # Here we analyze the blasts hits for each TAXON/PROTEIN pair:
            # a) number of blast hits
            # b) number of false positive according to the description (categorized in a list)
            #
            # Open the individual results file for a protein and a phylum
            XML_RESULTS_FILE = '{0}{1}/blast_results_{2}.xml'.format(args.directory, protein, taxon)
            # Create database connection
            conn = create_connection(DATABASE_FILE)

            # Get extended information about the false-positive hits
            negative_dict = {}
            for related_protein in protein_data[4]:
                negative_dict[related_protein] = 0
            # Get extended information about the true-positive hits
            positive_dict = {}
            for synonym in protein_data[3]:
                positive_dict[synonym] = 0

            print('\nAnalysing {0} ({1}):'.format(protein, taxon).upper())
            blastp_result = SearchIO.read(XML_RESULTS_FILE, "blast-xml")
            number = len(blastp_result)
            # counter for hits
            i = 0
            # counter for unknown proteins
            u = 0
            for hit in blastp_result:
                found = False
                # get gid and protein_id
                hitident = hit.id.split('|')
                print('Analyzing #{0} (from {1} {2} {3} homologs): {4} - {5}'.format(i+1, number, taxon, protein_data[3][0], hitident[1], hitident[3]))
                # Extract species name
                species = re.search(r'\[(.*?)\]',hit.description).group(1)
                # BLASTHIT = [gid, protein_id, species, fasta_description]
                BLASTHIT = [hitident[1], hitident[3], species, hit.description]
                #print('Checking false- or true-positivity.')
                for related_protein in protein_data[4]:
                    if related_protein.lower() in str(hit.description).lower():
                        negative_dict[related_protein] += 1
                        print('Potential false-positive found: {0}'.format(related_protein))
                        found = True
                        TOTAL_COUNT += 1
                        #list_to_scrutinize.append([taxon, protein, related_protein, hitident[1], hitident[3], hit.description])
                        #break
                for synonym in protein_data[3]:
                    if synonym.lower() in str(hit.description).lower():
                        positive_dict[synonym] += 1
                        print('True positive found: {0}'.format(synonym))
                        #list_to_scrutinize.append([taxon, protein, synonym, hitident[1], hitident[3], hit.description])
                        found = True
                        TOTAL_COUNT += 1
                        #break
                if found == False:
                    what_kind_of_protein = run_backcheck_blast(hitident[1], protein_data)
                    if what_kind_of_protein == 'synonym':
                        positive_dict[synonym] += 1
                        TOTAL_COUNT += 1
                        print('True positive found after backcheck blast: {0}'.format(protein_data[3][0]))
                    elif what_kind_of_protein == 'related_protein':
                        negative_dict[related_protein] += 1
                        TOTAL_COUNT += 1
                        #print('Potential false-positive found: {0}'.format(related_protein))
                        print('Potential false-positive found.')
                    else:
                        # These are all the unknown proteins, that even a backcheck blast cannot identify
                        list_to_scrutinize.append([taxon, protein, 'unknown', hitident[1], hitident[3], hit.description])
                        u += 1
                        TOTAL_COUNT += 1
                        TOTAL_UNKNOWN_COUNT += 1
                with conn:
                    # If the gid is already in the database, the insertion fails and 0 is returned (otherwise the GID)
                    print(str(db_insert_protein(conn, BLASTHIT, False))+'\n')
                i += 1
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
        print('Analysis for taxon {0} completed.\n'.format(taxon))
        print('new_protein_data: {0}'.format(new_protein_data))
        write_to_html_scrutinize_file(list_to_scrutinize)
        # Returns a dictionary with related proteins and the number of hits for them (according to fasta description)
        return new_protein_data

# This function should return the most common string in the hit description list
def run_backcheck_blast(GID, proteindata):
    # We need to define that this function should use the global variable LAST_BLAST_REPLY_TIME,
    # because we are changing LAST_BLAST_REPLY_TIME inside the function!
    global LAST_BLAST_REPLY_TIME
    try:
        blastp_result = SearchIO.read('data/backcheck/{0}.xml'.format(GID), 'blast-xml')
    except Exception as ex:
        print('No local backcheck blast results for GID {0}. Need to run remote blast.'.format(GID))
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
                print('Executing backcheck blast for gid {0}. Waiting for results...'.format(GID))
                result_handle = NCBIWWW.qblast("blastp", "nr", GID , hitlist_size = 50, expect = 0.01, format_type = "XML")
                i += 1
                with open('data/backcheck/{0}.xml'.format(GID), 'w') as out_handle:
                    out_handle.write(result_handle.read())
                result_handle.close()
                # Read in the data again to check that they contain really blast results
                with open('data/backcheck/{0}.xml'.format(GID), 'r') as xml_file:
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
            print("Something went wrong with the backcheck blast. Perhaps the Blast server did not respond? More about the error: " + str(ex))
    else:
        print('Using local backcheck blast results (data/backcheck/{0}.xml).'.format(GID))
    blastp_result = SearchIO.read('data/backcheck/{0}.xml'.format(GID), 'blast-xml')
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
    print('\n{0} from {1} backblast results confirm that {2} is {3}'.format(i, how_many, GID, proteindata[3][0]))
    print('{0} from {1} backblast results confirm that {2} is a related protein'.format(j, how_many, GID))
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
        TAXON_ID, phylum = get_taxon_id_and_phylum(species_name)
        full_genome_dictionary[phylum] += 1
        print('Adding +1 to fully sequenced {0}'.format(phylum))
    print(str(full_genome_dictionary))
    return full_genome_dictionary

def create_connection(DATABASE_FILE):
    global conn
    try:
        conn = sqlite3.connect(DATABASE_FILE)
        return conn
    except Error as err:
        print(err)
    return None

def db_insert_protein(CONNECTION, BLASTHIT, VERBOSE=True):
    if VERBOSE: print('Trying to execute SQL insertion ({0})...'.format(BLASTHIT[1]))
    # In some fasta descriptions one can really find this char: ' !!!!
    BLASTHIT[2] = BLASTHIT[2].replace('\'', '')
    query = 'INSERT INTO protein (gid, protein_id, species, fasta_description) VALUES ({0}, \'{1}\', \'{2}\', \'{3}\')'.format(BLASTHIT[0], BLASTHIT[1], BLASTHIT[2], BLASTHIT[3])
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.lastrowid
        CONNECTION.commit()
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    return result

# Give a list of scientific_name, taxon_id and phylum to insert into SQLite database
# There should never be the need to update any entry from this table, only to insert
def db_insert_species(CONNECTION, SPECIES_NAME, TAXON_ID, PHYLUM, VERBOSE=True):
    if VERBOSE: print('Trying to execute SQL insertion ({0})...'.format(SPECIES_NAME))
    query = 'INSERT INTO species (scientific_name, taxon_id, phylum) VALUES (\'{0}\', {1}, \'{2}\')'.format(SPECIES_NAME, TAXON_ID, PHYLUM)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
        cur.execute(query)
        result = cur.lastrowid
        if VERBOSE: print('SQL insertion successful.')
        CONNECTION.commit()
    except sqlite3.Error as err:
        if VERBOSE: print('Database error: {0}'.format(err))
        result = 0
    except Exception as err:
        if VERBOSE: print('Unknown error: {0}'.format(err))
        result = 0
    return result

# Return a list of scientific_name, taxon_id and phylum when queried by scientific_name
def db_retrieve_species(CONNECTION, SPECIES_NAME, VERBOSE=True):
    if VERBOSE: print('Trying to retrieve SQLite data for {0}...'.format(SPECIES_NAME))
    query = 'SELECT * FROM species WHERE scientific_name = \'{0}\''.format(SPECIES_NAME)
    if VERBOSE: print('query: {0}'.format(query))
    try:
        cur = CONNECTION.cursor()
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

def run():
    global APPLICATION_PATH, taxon_dictionary, master_dictionary, blacklist, DATABASE_FILE, HTML_SCRUTINIZE_FILE, LAST_BLAST_REPLY_TIME, BLAST_WAITING_TIME, TOTAL_COUNT, TOTAL_UNKNOWN_COUNT
    TOTAL_COUNT = 0
    TOTAL_UNKNOWN_COUNT = 0
    BLAST_WAITING_TIME = 60
    LAST_BLAST_REPLY_TIME = time.time()-BLAST_WAITING_TIME
    print(LAST_BLAST_REPLY_TIME)
    # Determine directory of script (in order to load the data files)
    APPLICATION_PATH =  os.path.abspath(os.path.dirname(__file__))
    #print('\nThe script is located in {0}'.format(APPLICATION_PATH))
    TAXON_DICTIONARY_FILE = '{0}/data/taxon_data.py'.format(APPLICATION_PATH)
    CSV_FILE = '{0}/data/genomes.csv'.format(APPLICATION_PATH)
    DATABASE_FILE = '{0}/data/database.sqlite3'.format(APPLICATION_PATH)
    HTML_SCRUTINIZE_FILE = '{0}/data/check_manually.html'.format(APPLICATION_PATH)
    # Delete the old HTML file
    if os.path.isfile(HTML_SCRUTINIZE_FILE):
        os.remove(HTML_SCRUTINIZE_FILE)
    preamble1, taxon_dictionary = load_dictionary(TAXON_DICTIONARY_FILE)
    #print('taxon_dictionary: {0}\n'.format(taxon_dictionary))
    #print("Populating new taxon dictionary")
    new_taxon_dictionary = taxon_dictionary
    preamble2, master_dictionary = load_dictionary('{0}/data/master_dictionary.py'.format(APPLICATION_PATH))
    blacklist = load_blacklist()
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
            # try:
            #     new_taxon_id, phylum = get_taxon_id_and_phylum(taxon)
            # except Exception as err:
            #     print('Could not get new taxon id and phylum. Error was: {0}',value(err))
            # # execute if no exceptions have occured
            # else:
            #     new_taxon_dictionary[taxon][0] = new_taxon_id
            # NUMBER OF SPECIES IN NCBI DATABASE
            # try:
            #     new_species_number = get_species_number(taxon)
            # except Exception as err:
            #     print('Could not get new species number. Error was: {0}',value(err))
            # else:
            #     new_taxon_dictionary[taxon][1] = new_species_number
            # NUMBER OF SEQUENCES IN NCBI DATABASE
            # try:
            #     new_sequence_number = get_sequence_number(taxon)
            # except Exception as err:
            #     print('Could not get new sequence number. Error was: {0}',value(err))
            # else:
            #     new_taxon_dictionary[taxon][3] = new_sequence_number
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
    if insert_line_breaks(TAXON_DICTIONARY_FILE) == True:
        print('Successfully formatted the taxon data file.')
    else:
        print('Formating the taxon data file failed.')
    print('From a total of {0} analyzed sequences, {1} were not classified and need to be manually checked.'.format(TOTAL_COUNT, TOTAL_UNKNOWN_COUNT))

if __name__ == '__main__':
    run()
