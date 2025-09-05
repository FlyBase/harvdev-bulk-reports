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
from harvdev_utils.char_conversions import sub_sup_sgml_to_plain_text
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
    get_component_info(split_system_dict)
    get_stock_info(split_system_dict)
    get_split_system_synonyms(split_system_dict)
    get_split_system_references(split_system_dict)
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


def get_component_info(split_system_dict):
    """Retrieve split system combination components."""
    global CONN
    log.info('Retrieve split system combination components.')
    fb_ss_component_query = """
        SELECT DISTINCT s.feature_id, o.uniquename||';'||o.name
        FROM feature s
        JOIN feature_relationship fr ON fr.subject_id = s.feature_id
        JOIN feature o ON o.feature_id = fr.object_id
        JOIN cvterm c ON c.cvterm_id = fr.type_id
        WHERE s.is_obsolete IS FALSE
          AND s.uniquename ~ '^FBco[0-9]{7}$'
          AND o.is_obsolete IS FALSE
          AND o.uniquename ~ '^FBal[0-9]{7}$'
          AND c.name = 'partially_produced_by';
    """
    ret_ss_component_info = connect(fb_ss_component_query, 'no_query', CONN)
    DB_ID = 0
    COMPONENT = 1
    counter = 0
    for row in ret_ss_component_info:
        split_system_dict[row[DB_ID]]['Component_Alleles'].append(row[COMPONENT])
        counter += 1
    log.info(f'Found {counter} split system combination components in chado.')
    return


def get_stock_info(split_system_dict):
    """Retrieve split system combination stocks."""
    global CONN
    log.info('Retrieve split system combination stocks.')
    fb_ss_stock_query = """
        SELECT DISTINCT f.feature_id, fp.value
        FROM feature f
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm c ON c.cvterm_id = fp.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBco[0-9]{7}$'
          AND c.name ~ '^derived_stock_';
    """
    ret_ss_stock_info = connect(fb_ss_stock_query, 'no_query', CONN)
    DB_ID = 0
    STOCK_PROP = 1
    counter = 0
    for row in ret_ss_stock_info:
        raw_stock_prop_list = row[STOCK_PROP].split('\n')
        cleaned_stock_prop_list = [sub_sup_sgml_to_plain_text(i.split('\t')[1].replace('@', '')) for i in raw_stock_prop_list]
        # log.debug(f'Found this stock prop: {cleaned_stock_prop_list}')
        split_system_dict[row[DB_ID]]['Stocks'].extend(cleaned_stock_prop_list)
        counter += len(cleaned_stock_prop_list)
    log.info(f'Found {counter} split system combination stocks in chado.')
    return


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
        split_system_dict[row[DB_ID]]['References'].append(row[PUB_ID])
        counter += 1
    log.info(f'Found {counter} split system combination references in chado.')
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
