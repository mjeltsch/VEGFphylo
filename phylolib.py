#!/usr/bin/python3
# -*- coding: utf-8 -*-
#

import os, sqlite3

def load_blacklist():
    # All the phyla in the list below have been found in the first screen not to have any VEGF-like molecules
    #blacklist = ['ctenophora', 'porifera', 'placozoa', 'xenacoelomorpha', 'cyclostomata', 'onychophora', 'pycnogonida', 'myriapoda', 'nematomorpha', 'loricifera', 'kinorhyncha', 'chaetognatha', 'bryozoa', 'entoprocta', 'cycliophora', 'nemertea', 'phoroniformea', 'gastrotricha', 'platyhelminthes', 'gnathostomulida', 'micrognathozoa', 'orthonectida', 'dicyemida']
    # All hits for the phyla below were manually checked and changed in taxon_data.py => therefore, we do not want to chenge them anymore
    #blacklist = ['porifera', 'xenacoelomorpha', 'myriapoda', 'bryozoa', 'entoprocta', 'platyhelminthes', 'orthonectida']
    #blacklist = ['kinorhyncha']
    blacklist = []
    return blacklist

def make_synonym_dictionary(master_dictionary):
    #print('master_dictionary:\n{0}'.format(master_dictionary))
    synonym_dictionary = {}
    for key, value in master_dictionary.items():
        #print('key:\n{0}'.format(key))
        #print('value:\n{0}'.format(value))
        #print('value[3][0]:\n{0}'.format(value[3][0]))
        canonical_ortholog_group_name = value[3][0]
        print('canonical_ortholog_group_name: {0}'.format(canonical_ortholog_group_name))
        for item in value[3]:
            synonym_dictionary[item] = canonical_ortholog_group_name
    return synonym_dictionary

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

def create_sqlite_file(FILE_NAME):
    try:
        conn = sqlite3.connect(FILE_NAME)
        print(sqlite3.version)
        commands = '''CREATE TABLE IF NOT EXISTS `species` (
        	`scientific_name`	TEXT NOT NULL,
        	`taxon_id`	INTEGER NOT NULL,
        	`phylum`	TEXT NOT NULL,
        	PRIMARY KEY(`taxon_id`));
        CREATE TABLE IF NOT EXISTS `protein` (
        	`gid`	INTEGER NOT NULL,
        	`protein_id`	TEXT,
        	`fasta_description`	TEXT NOT NULL,
        	`species`	TEXT NOT NULL,
        	`ortholog_group`	TEXT,
        	`curated_manually_by`	TEXT DEFAULT ('NULL'),
        	PRIMARY KEY(`gid`));
        CREATE TABLE IF NOT EXISTS `ortholog_groups` (
        	`ortholog_group`	TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS `curator` (
        	`curator`	TEXT NOT NULL);'''
        cur = conn.cursor()
        cur.execute(commands)
        conn.commit()
    except Error as err:
        print(err)
        return False
    else:
        return True
    finally:
        conn.close()

# Converts the number of seconds into a days/hours/minutes/seconds string
def execution_time_str(elapsed_time_seconds):
    min, sec = divmod(elapsed_time_seconds, 60)
    hours, min = divmod(min, 60)
    days, hours = divmod(hours, 24)
    return (str(days) + " days " if days != 0 else '') + (str(hours) + " hours " if hours != 0 else '') + str(min)[:-2] + " min " + str(round(sec, 1)) + " sec"
