# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report human disease model data.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_hdm.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_hdm.py -v -c /path/to/config.cfg

"""

import argparse
# import configparser
# import csv
# import datetime
# import logging
# import os
# import pickle
# import psycopg2
# import re
# import sys
from harvdev_utils.char_conversions import clean_free_text
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect    # other useful functions: add_unique_info, add_list_info, add_unique_dict_info
)
# from harvdev_utils.psycopg_functions.sql_queries import (
#     current_feat_symbol_sgmls, current_feat_fullname_sgmls, feat_symbol_synonyms, feat_fullname_synonyms,
#     feat_secondary_fbids, orgid_abbr, orgid_genus, indirect_rel_features, rel_features, rel_dmel_features,
#     featureprops, feat_cvterm_cvtprop    # current_features, feat_id_symbol_sgml
# )

# Global variables for the output file. Header order will match list order below.
report_label = 'human_disease_models'
report_title = 'FlyBase Human Disease Models Report'
header_list = [

]

# Proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
assembly = set_up_dict['assembly']
annotation_release = set_up_dict['annotation_release']
database = set_up_dict['database']
database_release = set_up_dict['database_release']
alliance_schema = set_up_dict['alliance_schema']
output_dir = set_up_dict['output_dir']
output_filename = set_up_dict['output_filename']
log = set_up_dict['log']
CONN = set_up_dict['conn']
the_time = set_up_dict['the_time']

# Process more input parameters (-c and -v handled by set_up_db_reading() function above).
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-i', '--input_filename', help='Input TSV file.', required=False)
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
input_filename = args.input_filename


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    database_info = get_database_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(database_info)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_database_info():
    """Retrieve a list of gene-gene/TE relationships"""
    global CONN
    log.info('Querying database for paralog info.')
    fb_gene_rel_query = """
        SELECT s.uniquename, s.name, t.name, o.uniquename, o.name
        FROM feature s
        JOIN feature_relationship fr ON fr.subject_id = s.feature_id
        JOIN feature o ON o.feature_id = fr.object_id
        JOIN cvterm t ON t.cvterm_id = fr.type_id
        JOIN organism org ON org.organism_id = s.organism_id
        WHERE s.is_obsolete IS FALSE
          AND s.uniquename ~ '^FBgn[0-9]{7}$'
          AND org.abbreviation = 'Dmel'
          AND o.is_obsolete IS FALSE
          AND o.uniquename ~ '^FB(gn|te)[0-9]{7}$'
          AND t.name IN ('has_component_gene', 'encoded_by', 'member_gene_of')
    ;"""
    ret_gene_rel_info = connect(fb_gene_rel_query, 'no_query', CONN)
    log.info(f'Found {len(ret_gene_rel_info)} gene-gene/TE relationships for Dmel gene subjects in chado.')
    return ret_gene_rel_info

def process_database_info(input_data):
    """Convert the list of SQL results into dicts for TSV output.

    Args:
        input_data (list): A list of tuples representing SQL query output.

    Returns:
        A list of dictionaries representing, in this case, paralog information.
    """
    log.info('Starting to process gene relationship info retrieved from database.')
    data_list = []
    SBJ_ID = 0
    SBJ_SYMB = 1
    REL_TYPE = 2
    OBJ_ID = 3
    OBJ_SYMB = 4
    counter = 0
    reciprocal_counter = 0
    for row in input_data:
        # For this relation type, fix chado reversed roles.
        if row[REL_TYPE] == 'has_component_gene':
            processed_rel_type = 'encoded_by'
        else:
            processed_rel_type = row[REL_TYPE]
        row_dict = {
        'Subject_FBgn_ID': row[SBJ_ID],
        'Subject_Symbol': row[SBJ_SYMB],
        'Relationship': processed_rel_type,
        'Object_FB_ID': row[OBJ_ID],
        'Object_Symbol': row[OBJ_SYMB],
        }
        data_list.append(row_dict)
        counter += 1
        # Mimic web reporting of reciprocal relationships for "member_gene_of" cases.
        if processed_rel_type == 'member_gene_of':
            reciprocal_row_dict = {
                'Subject_FBgn_ID': row[OBJ_ID],
                'Subject_Symbol': row[OBJ_SYMB],
                'Relationship': 'has_component_gene',
                'Object_FB_ID': row[SBJ_ID],
                'Object_Symbol': row[SBJ_SYMB],
            }
            data_list.append(reciprocal_row_dict)
            reciprocal_counter += 1
    log.info(f'Processed {counter} gene relationships.')
    log.info(f'Inferred an additional {reciprocal_counter} reciprocal gene relationships.')
    return data_list


if __name__ == "__main__":
    main()
