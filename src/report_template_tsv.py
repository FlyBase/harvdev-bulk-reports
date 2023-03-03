# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Template script for data retrieval of chado info and printing to TSV file.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_tsv_template.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_tsv_template.py -v -c /path/to/config.cfg

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
# from harvdev_utils.char_conversions import *

# Global variables for the output file. Header order will match list order below.
report_label = 'this_report_label'
report_title = 'FlyBase FILL_IN_THE_BLANKS report'
header_list = [
    'FBgn_ID',
    'Gene_Symbol',
    'Arm/Scaffold',
    'Location',
    'Strand',
    'Paralog_FBgn_ID',
    'Paralog_GeneSymbol',
    'Paralog_Arm/Scaffold',
    'Paralog_Location',
    'Paralog_Strand',
    'DIOPT_score'
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
conn = set_up_dict['conn']
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
    database_info = get_database_info(conn)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(database_info)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of specific data types.
def get_database_info(db_connection):
    """Retrieve a list of tuples containing paralog pairs and related info.

    Args:
        arg1 (psycopg2.extensions.connection)

    Returns:
        A list of tuples representing SQL query output.

    """
    log.info('Querying database for paralog info.')
    fb_paralog_query = """
        SELECT f1.uniquename, f1.name, src1.uniquename, fl1.fmin+1||'..'||fl1.fmax, fl1.strand,
        f2.uniquename, f2.name, src2.uniquename, fl2.fmin+1||'..'||fl2.fmax, fl2.strand, fr.value
        FROM feature f1
        JOIN feature_relationship fr ON (fr.subject_id = f1.feature_id)
        JOIN feature_relationshipprop frp ON (frp.feature_relationship_id = fr.feature_relationship_id)
        JOIN organism o1 ON (o1.organism_id = f1.organism_id)
        JOIN featureloc fl1 ON (fl1.feature_id = f1.feature_id)
        JOIN feature src1 ON (src1.feature_id = fl1.srcfeature_id)
        JOIN cvterm cvtfr ON (cvtfr.cvterm_id = fr.type_id)
        JOIN feature f2 ON (f2.feature_id = fr.object_id)
        JOIN organism o2 ON (o2.organism_id = f2.organism_id)
        JOIN featureloc fl2 ON (fl2.feature_id = f2.feature_id)
        JOIN feature src2 ON (src2.feature_id = fl2.srcfeature_id)
        WHERE f1.is_obsolete = false and f1.uniquename ~ '^FBgn[0-9]{7}$' and o1.abbreviation = 'Dmel'
          and f2.is_obsolete = false and f2.uniquename ~ '^FBgn[0-9]{7}$' and o2.abbreviation = 'Dmel'
          and cvtfr.name = 'paralogous_to' and frp.value = 'DIOPT'
          and src1.is_obsolete = false and src1.organism_id = f1.organism_id and src1.type_id = 553
          and src2.is_obsolete = false and src2.organism_id = f2.organism_id and src2.type_id = 553;"""
    ret_paralog_info = connect(fb_paralog_query, 'no_query', db_connection)
    log.info('Found {} paralogous relationships among Dmel genes.'.format(len(ret_paralog_info)))

    return ret_paralog_info


def process_database_info(input_data):
    """Convert a list of SQL results into a list of dictionaries for TSV output.

    Args:
        arg1 (list): A list of tuples representing SQL query output.

    Returns:
        A list of dictionaries representing, in this case, paralog information.
    """
    log.info('Starting to process paralog info retrieved from database.')
    data_list = [{
        'FBgn_ID': i[0],
        'GeneSymbol': i[1],
        'Arm/Scaffold': i[2],
        'Location': i[3],
        'Strand': i[4],
        'Paralog_FBgn_ID': i[5],
        'Paralog_GeneSymbol': i[6],
        'Paralog_Arm/Scaffold': i[7],
        'Paralog_Location': i[8],
        'Paralog_Strand': i[9],
        'DIOPT_score': i[10].count(',') + 1    # Does this count commas correctly?
        }
        for i in input_data]
    log.info('Done processing paralog info into a list of dictionaries.')

    return data_list


if __name__ == "__main__":
    main()
