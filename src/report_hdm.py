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
    'related_disease_FBhh',
    'related_disease_name',
    'child_disease_FBhh',
    'child_disease_name',
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
    get_parent_hdms(hdm_dict)
    get_related_hdms(hdm_dict)
    get_child_hdms(hdm_dict)
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
            'sub-datatype': 'disease',
            'category': None,
            'parent_disease_list': [],
            'parent_disease_FBhh': None,
            'parent_disease_name': None,
            'related_disease_list': [],
            'related_disease_FBhh': None,
            'related_disease_name': None,
            'child_disease_list': [],
            'child_disease_FBhh': None,
            'child_disease_name': None,
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
    fb_hdm_category_query = """
        SELECT DISTINCT hh.humanhealth_id, hhp.value
        FROM humanhealth hh
        JOIN humanhealthprop hhp ON hhp.humanhealth_id = hh.humanhealth_id
        JOIN cvterm t ON t.cvterm_id = hhp.type_id
        WHERE hh.is_obsolete IS FALSE
          AND hh.uniquename ~ '^FBhh[0-9]{7}$'
          AND t.name = 'category';
    """
    ret_hdm_category_info = connect(fb_hdm_category_query, 'no_query', CONN)
    DB_ID = 0
    CATEGORY = 1
    counter = 0
    for row in ret_hdm_category_info:
        hdm_dict[row[DB_ID]]['category'] = row[CATEGORY]
        counter += 1
    log.info(f'Found {counter} category annotations for human health disease models in chado.')
    return


def get_parent_hdms(hdm_dict):
    """Retrieve parent human health disease models."""
    global CONN
    log.info('Retrieve parent human health disease models.')
    fb_parent_hdm_query = """
        SELECT DISTINCT s.humanhealth_id, o.uniquename, o.name
        FROM humanhealth s
        JOIN humanhealth_relationship hhr ON hhr.subject_id = s.humanhealth_id
        JOIN humanhealth o ON o.humanhealth_id = hhr.object_id
        JOIN cvterm cvthhr ON cvthhr.cvterm_id = hhr.type_id
        WHERE s.is_obsolete IS FALSE
          AND s.uniquename ~ '^FBhh[0-9]{7}$'
          AND o.is_obsolete IS FALSE
          AND o.uniquename ~ '^FBhh[0-9]{7}$'
          AND cvthhr.name = 'belongs_to'
        ORDER BY o.name;
    """
    ret_parent_hdm_info = connect(fb_parent_hdm_query, 'no_query', CONN)
    DB_ID = 0
    PARENT_UNAME = 1
    PARENT_NAME = 2
    counter = 0
    for row in ret_parent_hdm_info:
        parent_tuple = (row[DB_ID], row[PARENT_UNAME], row[PARENT_NAME])
        hdm_dict[row[DB_ID]]['parent_disease_list'].append(parent_tuple)
        counter += 1
    for hdm in hdm_dict.values():
        hdm['parent_disease_FBhh'] = '|'.join([i[PARENT_UNAME] for i in hdm['parent_disease_list']])
        hdm['parent_disease_name'] = '|'.join([i[PARENT_NAME] for i in hdm['parent_disease_list']])
    log.info(f'Found {counter} parent human health disease models in chado.')
    return


def get_related_hdms(hdm_dict):
    """Retrieve related human health disease models."""
    global CONN
    log.info('Retrieve related human health disease models.')
    fb_related_hdm_query = """
        SELECT DISTINCT s.humanhealth_id, s.uniquename, s.name, o.humanhealth_id, o.uniquename, o.name
        FROM humanhealth s
        JOIN humanhealth_relationship hhr ON hhr.subject_id = s.humanhealth_id
        JOIN humanhealth o ON o.humanhealth_id = hhr.object_id
        JOIN cvterm cvthhr ON cvthhr.cvterm_id = hhr.type_id
        WHERE s.is_obsolete IS FALSE
          AND s.uniquename ~ '^FBhh[0-9]{7}$'
          AND o.is_obsolete IS FALSE
          AND o.uniquename ~ '^FBhh[0-9]{7}$'
          AND cvthhr.name = 'associated_with';
    """
    ret_related_hdm_info = connect(fb_related_hdm_query, 'no_query', CONN)
    SBJ_ID = 0
    SBJ_UNAME = 1
    SBJ_NAME = 2
    OBJ_ID = 3
    OBJ_UNAME = 4
    OBJ_NAME = 5
    counter = 0
    for row in ret_related_hdm_info:
        counter += 1
        sbj_tuple = (row[SBJ_ID], row[SBJ_UNAME], row[SBJ_NAME])
        obj_tuple = (row[OBJ_ID], row[OBJ_UNAME], row[OBJ_NAME])
        hdm_dict[row[SBJ_ID]]['related_disease_list'].append(obj_tuple)
        hdm_dict[row[OBJ_ID]]['related_disease_list'].append(sbj_tuple)
    REL_UNAME = 1
    REL_NAME = 2
    for hdm in hdm_dict.values():
        rel_hdm_dict = {}
        sorted_rel_hdms = []
        for rel_tuple in hdm['related_disease_list']:
            rel_hdm_dict[rel_tuple[REL_NAME]] = rel_tuple
        sorted_rel_hdm_names = sorted(list(rel_hdm_dict.keys()))
        for rel_hdm_name in sorted_rel_hdm_names:
            sorted_rel_hdms.append(rel_hdm_dict[rel_hdm_name])
        hdm['related_disease_FBhh'] = '|'.join([i[REL_UNAME] for i in sorted_rel_hdms])
        hdm['related_disease_name'] = '|'.join([i[REL_NAME] for i in sorted_rel_hdms])
    log.info(f'Found {counter} related human health disease models in chado.')
    return


def get_child_hdms(hdm_dict):
    """Retrieve child human health disease models."""
    global CONN
    log.info('Retrieve child human health disease models.')
    fb_child_hdm_query = """
        SELECT DISTINCT o.humanhealth_id, s.uniquename, s.name
        FROM humanhealth s
        JOIN humanhealth_relationship hhr ON hhr.subject_id = s.humanhealth_id
        JOIN humanhealth o ON o.humanhealth_id = hhr.object_id
        JOIN cvterm cvthhr ON cvthhr.cvterm_id = hhr.type_id
        WHERE s.is_obsolete IS FALSE
          AND s.uniquename ~ '^FBhh[0-9]{7}$'
          AND o.is_obsolete IS FALSE
          AND o.uniquename ~ '^FBhh[0-9]{7}$'
          AND cvthhr.name = 'belongs_to'
        ORDER BY s.name;
    """
    ret_child_hdm_info = connect(fb_child_hdm_query, 'no_query', CONN)
    DB_ID = 0
    CHILD_UNAME = 1
    CHILD_NAME = 2
    counter = 0
    for row in ret_child_hdm_info:
        child_tuple = (row[DB_ID], row[CHILD_UNAME], row[CHILD_NAME])
        hdm_dict[row[DB_ID]]['child_disease_list'].append(child_tuple)
        counter += 1
    for hdm in hdm_dict.values():
        hdm['child_disease_FBhh'] = '|'.join([i[CHILD_UNAME] for i in hdm['child_disease_list']])
        hdm['child_disease_name'] = '|'.join([i[CHILD_NAME] for i in hdm['child_disease_list']])
    log.info(f'Found {counter} child human health disease models in chado.')
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
