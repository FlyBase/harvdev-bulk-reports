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
# from harvdev_utils.char_conversions import clean_free_text
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
    'FB_id',
    'name',
    'name_synonyms',
    'sub-datatype',
    'category',
    'parent_disease_FBhh',
    'parent_disease_name',
    'related_FBhh',
    'parent_entity_children_FBhh',
    'group_entity_children_FBhh',
    'OMIM_disease_MIM',
    'OMIM_disease_name',
    'OMIM_gene_MIM',
    'OMIM_gene_name',
    'HGNC_gene',
    'HGNC_name',
    'DO_ID',
    'DO_name',
    'external_links',
    'related_specific_diseases',
    'implicated_human_gene',
    'implicated_Dmel_gene',
    'implicated_other_gene',
    'description_overview',
    'description_symptoms',
    'description_genetics',
    'description_cellular',
    'description_molecular',
    'BDSC_link',
]

# Proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
assembly = set_up_dict['assembly']
database = set_up_dict['database']
database_release = set_up_dict['database_release']
output_dir = set_up_dict['output_dir']
output_filename = set_up_dict['output_filename']
log = set_up_dict['log']
CONN = set_up_dict['conn']

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
    hdm_dict = get_initial_hdm_info()
    get_hdm_synonyms(hdm_dict)
    get_hdm_subtypes(hdm_dict)
    get_hdm_categories(hdm_dict)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(hdm_dict)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_initial_hdm_info():
    """Retrieve human health disease models."""
    global CONN
    log.info('Retrieve human health disease models.')
    fb_hdm_query = """
        SELECT DISTINCT hh.humanhealth_id, hh.uniquename, hh.name
        FROM humanhealth hh
        WHERE hh.is_obsolete IS FALSE
          AND hh.uniquename ~ '^FBhh[0-9]{7}$';
    """
    ret_hdm_info = connect(fb_hdm_query, 'no_query', CONN)
    DB_ID = 0
    UNAME = 1
    NAME = 2
    hdm_dict = {}
    counter = 0
    for row in ret_hdm_info:
        hdm_result = {
            'db_id': row[DB_ID],
            'FB_id': row[UNAME],
            'name': row[NAME],
            'name_synonyms_list': [],
            'name_synonyms': None,
            'sub-datatype': 'disease',    # The default. Update the rare exceptions.
            'category': None,
            'parent_disease_FBhh': None,
            'parent_disease_name': None,
            'related_FBhh': None,
            'parent_entity_children_FBhh': None,
            'group_entity_children_FBhh': None,
            'OMIM_disease_MIM': None,
            'OMIM_disease_name': None,
            'OMIM_gene_MIM': None,
            'OMIM_gene_name': None,
            'HGNC_gene': None,
            'HGNC_name': None,
            'DO_ID': None,
            'DO_name': None,
            'external_links': None,
            'related_specific_diseases': None,
            'implicated_human_gene': None,
            'implicated_Dmel_gene': None,
            'implicated_other_gene': None,
            'description_overview': None,
            'description_symptoms': None,
            'description_genetics': None,
            'description_cellular': None,
            'description_molecular': None,
            'BDSC_link': None,
        }
        hdm_dict[row[DB_ID]] = hdm_result
        counter += 1
    log.info(f'Found {counter} human health disease models in chado.')
    return hdm_dict


def get_hdm_synonyms(hdm_dict):
    """Retrieve human health disease model synonyms."""
    global CONN
    log.info('Retrieve human health disease model synonyms.')
    fb_hdm_syno_query = """
        SELECT DISTINCT hh.humanhealth_id, s.name
        FROM humanhealth hh
        JOIN humanhealth_synonym hhs ON hhs.humanhealth_id = hh.humanhealth_id
        JOIN synonym s ON s.synonym_id = hhs.synonym_id
        WHERE hh.is_obsolete IS FALSE
          AND hh.uniquename ~ '^FBhh[0-9]{7}$'
          AND s.name != hh.name;
    """
    ret_hdm_syno_info = connect(fb_hdm_syno_query, 'no_query', CONN)
    DB_ID = 0
    ALIAS = 1
    counter = 0
    for row in ret_hdm_syno_info:
        hdm_dict[row[DB_ID]]['name_synonyms_list'].append(row[ALIAS])
        counter += 1
    for hdm in hdm_dict.values():
        hdm['name_synonyms_list'].sort()
        hdm['name_synonyms'] = '|'.join(hdm['name_synonyms_list'])
    log.info(f'Found {counter} human health disease model synonyms in chado.')
    return


def get_hdm_subtypes(hdm_dict):
    """Retrieve human health disease model subtypes."""
    global CONN
    log.info('Retrieve human health disease model subtypes.')
    # All HDMs are "disease" subtype at the time of composing this script.
    # So, use a default "disease" value, and only update exceptions when they happen.
    fb_hdm_subtype_query = """
        SELECT DISTINCT hh.humanhealth_id, hhp.value
        FROM humanhealth hh
        JOIN humanhealthprop hhp ON hhp.humanhealth_id = hh.humanhealth_id
        JOIN cvterm t ON t.cvterm_id = hhp.type_id
        WHERE hh.is_obsolete IS FALSE
          AND hh.uniquename ~ '^FBhh[0-9]{7}$'
          AND t.name = 'sub_datatype'
          AND hhp.value != 'disease';
    """
    ret_hdm_subtype_info = connect(fb_hdm_subtype_query, 'no_query', CONN)
    DB_ID = 0
    SUBTYPE = 1
    counter = 0
    for row in ret_hdm_subtype_info:
        hdm_dict[row[DB_ID]]['sub-datatype'] = row[SUBTYPE]
        counter += 1
    log.info(f'Found {counter} non-"disease" subtype annotations for human health disease models in chado.')
    return


def get_hdm_categories(hdm_dict):
    """Retrieve human health disease model categories."""
    global CONN
    log.info('Retrieve human health disease model categories.')
    fb_hdm__category_query = """
        SELECT DISTINCT hh.humanhealth_id, hhp.value
        FROM humanhealth hh
        JOIN humanhealthprop hhp ON hhp.humanhealth_id = hh.humanhealth_id
        JOIN cvterm t ON t.cvterm_id = hhp.type_id
        WHERE hh.is_obsolete IS FALSE
          AND hh.uniquename ~ '^FBhh[0-9]{7}$'
          AND t.name = 'category';
    """
    ret_hdm_category_info = connect(fb_hdm__category_query, 'no_query', CONN)
    DB_ID = 0
    CATEGORY = 1
    counter = 0
    for row in ret_hdm_category_info:
        hdm_dict[row[DB_ID]]['category'] = row[CATEGORY]
        counter += 1
    log.info(f'Found {counter} category annotations for human health disease models in chado.')
    return


def process_database_info(input_data):
    """Convert the HDM dict to a list of data elements."""
    log.info('Convert the HDM dict to a list of data elements.')
    data_list = []
    counter = 0
    for i in input_data.values():
        data_list.append(i)
        counter += 1
    log.info(f'Sending {counter} HDM entries to the export file.')
    return data_list


if __name__ == "__main__":
    main()
