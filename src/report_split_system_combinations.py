# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report split system combinations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_split_system_combinations.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_split_system_combinations.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
REPORT_LABEL = 'split_system_combinations'
REPORT_TITLE = 'FlyBase Split System Combination Report'
HEADER_LIST = [
    'FB_id',
    'Symbol',
    'Component_Alleles',
    'Stocks',
    'Synonyms',
    'References',
]

# Proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
DATABASE = set_up_dict['database']
OUTPUT_FILENAME = set_up_dict['output_filename']
log = set_up_dict['log']
CONN = set_up_dict['conn']

# Process more input parameters (-c and -v handled by set_up_db_reading() function above).
parser = argparse.ArgumentParser(description='inputs')
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    split_system_dict = get_initial_split_system_info()
    get_split_system_synonyms(split_system_dict)
    get_split_system_references(split_system_dict)

    #######################################
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = process_database_info(split_system_dict)
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_initial_split_system_info():
    """Retrieve split system combinations."""
    global CONN

    log.info('Retrieve split system combinations.')
    fb_ss_query = """
        SELECT DISTINCT f.feature_id, f.uniquename, f.name
        FROM feature f
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBco[0-9]{7}$';
    """
    ret_ss_info = connect(fb_ss_query, 'no_query', CONN)
    DB_ID = 0
    UNAME = 1
    NAME = 2
    split_system_dict = {}
    counter = 0
    for row in ret_ss_info:
        ss_result = {
            'db_id': row[DB_ID],
            'FB_id': row[UNAME],
            'Symbol': row[NAME],
            'Component_Alleles': [],
            'Stocks': [],
            'Synonyms': [],
            'References': [],
        }
        split_system_dict[row[DB_ID]] = ss_result
        counter += 1
    log.info(f'Found {counter} split system combinations in chado.')
    return split_system_dict


def get_split_system_synonyms(split_system_dict):
    """Retrieve split system combination synonyms."""
    global CONN
    log.info('Retrieve split system combination synonyms.')
    fb_ss_syno_query = """
        SELECT DISTINCT f.feature_id, s.name
        FROM feature f
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBco[0-9]{7}$'
          AND s.name != f.name;
    """
    ret_ss_syno_info = connect(fb_ss_syno_query, 'no_query', CONN)
    DB_ID = 0
    ALIAS = 1
    counter = 0
    for row in ret_ss_syno_info:
        split_system_dict[row[DB_ID]]['Synonyms'].append(row[ALIAS])
        counter += 1
    log.info(f'Found {counter} split system combination synonyms in chado.')
    return


def get_split_system_references(split_system_dict):
    """Retrieve split system combination references."""
    global CONN
    log.info('Retrieve split system combination references.')
    fb_ss_pub_query = """
        SELECT DISTINCT f.feature_id, p.uniquename
        FROM feature f
        JOIN feature_pub fp ON fp.feature_id = f.feature_id
        JOIN pub p ON p.pub_id = fp.pub_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBco[0-9]{7}$'
          AND p.is_obsolete IS FALSE
          AND p.uniquename ~ '^FBrf[0-9]{7}$';
    """
    ret_ss_pub_info = connect(fb_ss_pub_query, 'no_query', CONN)
    DB_ID = 0
    PUB_ID = 1
    counter = 0
    for row in ret_ss_pub_info:
        split_system_dict[row[DB_ID]]['References'] = row[PUB_ID]
        counter += 1
    log.info(f'Found {counter} split system combinations in chado.')
    return


def process_database_info(input_data):
    """Convert the split system combination dict to a list of data elements."""
    log.info('Convert the split system combination dict to a list of data elements.')
    data_list = []
    counter = 0
    for i in input_data.values():
        for key, value in i.items():
            if type(value) is list:
                if not value:
                    i[key] = ''
                else:
                    value = list(set(value))
                    value.sort()
                    i[key] = '|'.join(value)
        data_list.append(i)
        counter += 1
    log.info(f'Sending {counter} split system combination entries to the export file.')
    return data_list


if __name__ == "__main__":
    main()
