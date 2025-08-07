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
import re
from harvdev_utils.char_conversions import clean_free_text
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
REPORT_LABEL = 'human_disease_models'
REPORT_TITLE = 'FlyBase Human Disease Models Report'
HEADER_LIST = [
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
    'OMIM_disease_ID',
    'OMIM_disease_name',
    'OMIM_gene_ID',
    'OMIM_gene_name',
    'HGNC_gene_ID',
    'HGNC_gene_name',
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
    hdm_dict = get_initial_hdm_info()
    get_hdm_synonyms(hdm_dict)
    get_hdm_subtypes(hdm_dict)
    get_hdm_categories(hdm_dict)
    get_parent_hdms(hdm_dict)
    get_related_hdms(hdm_dict)
    get_child_hdms(hdm_dict)
    get_hdm_omim_pheno_series(hdm_dict)
    get_hdm_omim_pheno_xrefs(hdm_dict)
    get_hdm_omim_table_xrefs(hdm_dict)
    get_hdm_omim_table_prop(hdm_dict)    # BOB
    hdm_relevant_gene_dict = build_hdm_gene_dict()
    get_hdm_genes(hdm_dict, hdm_relevant_gene_dict)
    get_hdm_do_terms(hdm_dict)
    get_external_links(hdm_dict)
    get_hdm_props(hdm_dict)
    get_hdm_omim_bdsc_links(hdm_dict)
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = process_database_info(hdm_dict)
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
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
            'OMIM_pheno_series': [],
            'OMIM_disease_xrefs': [],
            'OMIM_disease_ID': None,
            'OMIM_disease_name': None,
            'implicated_gene_feature_ids': [],
            'OMIM_gene_xrefs': [],
            'OMIM_gene_ID': None,
            'OMIM_gene_name': None,
            'HGNC_gene_xrefs': [],
            'HGNC_gene_ID': None,
            'HGNC_gene_name': None,
            'DO_cvterms': [],
            'DO_ID': None,
            'DO_name': None,
            'xref_urls': [],
            'external_links': None,
            'OMIM_pheno_table_xrefs': [],
            'related_specific_diseases': None,
            'implicated_human_gene': None,
            'implicated_Dmel_gene': None,
            'implicated_other_gene': None,
            'description_overview': [],
            'description_symptoms': [],
            'description_genetics': [],
            'description_cellular': [],
            'description_molecular': [],
            'bdsc_links': [],
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


def get_hdm_omim_pheno_series(hdm_dict):
    """Retrieve human disease model OMIM SERIES xrefs."""
    global CONN
    log.info('Retrieve human disease model OMIM SERIES xrefs.')
    fb_hdm_omim_series_query = """
        SELECT DISTINCT hh.humanhealth_id, dbx.accession, dbx.description
        FROM humanhealth hh
        JOIN humanhealth_dbxref hhdbx ON hhdbx.humanhealth_id = hh.humanhealth_id
        JOIN dbxref dbx ON dbx.dbxref_id = hhdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE hh.is_obsolete IS FALSE
          AND hhdbx.is_current IS TRUE
          AND db.name = 'OMIM_series';
    """
    ret_hdm_omim_series_info = connect(fb_hdm_omim_series_query, 'no_query', CONN)
    DB_ID = 0
    DBX_ACC = 1
    DBX_DESC = 2
    counter = 0
    for row in ret_hdm_omim_series_info:
        omim_series_tuple = (row[DB_ID], row[DBX_ACC], row[DBX_DESC])
        hdm_dict[row[DB_ID]]['OMIM_pheno_series'].append(omim_series_tuple)
        counter += 1
    log.info(f'Found {counter} human disease model OMIM SERIES xrefs in chado.')
    return


def get_hdm_omim_pheno_xrefs(hdm_dict):
    """Retrieve human disease model OMIM PHENOTYPE xrefs."""
    global CONN
    log.info('Retrieve human disease model OMIM PHENOTYPE xrefs.')
    fb_hdm_omim_pheno_query = """
        SELECT DISTINCT hh.humanhealth_id, dbx.accession, dbx.description
        FROM humanhealth hh
        JOIN humanhealth_dbxref hhdbx ON hhdbx.humanhealth_id = hh.humanhealth_id
        JOIN humanhealth_dbxrefprop hhdbxp ON hhdbxp.humanhealth_dbxref_id = hhdbx.humanhealth_dbxref_id
        JOIN cvterm cvt ON cvt.cvterm_id = hhdbxp.type_id
        JOIN dbxref dbx ON dbx.dbxref_id = hhdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE hh.is_obsolete IS FALSE
          AND hhdbx.is_current IS TRUE
          AND db.name = 'OMIM_PHENOTYPE'
          AND cvt.name = 'hh2c_link';
    """
    ret_hdm_omim_pheno_info = connect(fb_hdm_omim_pheno_query, 'no_query', CONN)
    DB_ID = 0
    DBX_ACC = 1
    DBX_DESC = 2
    counter = 0
    for row in ret_hdm_omim_pheno_info:
        omim_pheno_tuple = (row[DB_ID], row[DBX_ACC], row[DBX_DESC])
        hdm_dict[row[DB_ID]]['OMIM_disease_xrefs'].append(omim_pheno_tuple)
        counter += 1
    log.info(f'Found {counter} human disease model OMIM PHENOTYPE xrefs in chado.')
    uniq_counter = 0
    for hdm in hdm_dict.values():
        if len(hdm['OMIM_disease_xrefs']) == 1:
            uniq_tuple = hdm['OMIM_disease_xrefs'][0]
            hdm['OMIM_disease_ID'] = f'MIM:{uniq_tuple[DBX_ACC]}'
            hdm['OMIM_disease_name'] = uniq_tuple[DBX_DESC]
            uniq_counter += 1
        elif len(hdm['OMIM_disease_xrefs']) > 1:
            log.warning(f'{hdm["FB_id"]} has MANY "hh2c" OMIM_PHENOTYPE xrefs')
    log.info(f'Found {uniq_counter} human disease models having a single OMIM PHENOTYPE xref in chado.')
    return


def get_hdm_omim_table_xrefs(hdm_dict):
    """Retrieve human disease model OMIM table xrefs."""
    global CONN
    log.info('Retrieve human disease model OMIM table xrefs.')
    fb_hdm_omim_pheno_table_query = """
        SELECT DISTINCT hh.humanhealth_id, dbx.accession, dbx.description
        FROM humanhealth hh
        JOIN humanhealth_dbxref hhdbx ON hhdbx.humanhealth_id = hh.humanhealth_id
        JOIN humanhealth_dbxrefprop hhdbxp ON hhdbxp.humanhealth_dbxref_id = hhdbx.humanhealth_dbxref_id
        JOIN cvterm cvt ON cvt.cvterm_id = hhdbxp.type_id
        JOIN dbxref dbx ON dbx.dbxref_id = hhdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE hh.is_obsolete IS FALSE
          AND hhdbx.is_current IS TRUE
          AND db.name = 'OMIM_PHENOTYPE'
          AND cvt.name = 'OMIM_pheno_table';
    """
    ret_hdm_omim_pheno_table_info = connect(fb_hdm_omim_pheno_table_query, 'no_query', CONN)
    DB_ID = 0
    DBX_ACC = 1
    DBX_DESC = 2
    counter = 0
    for row in ret_hdm_omim_pheno_table_info:
        omim_pheno_tuple = (row[DB_ID], row[DBX_ACC], row[DBX_DESC])
        hdm_dict[row[DB_ID]]['OMIM_pheno_table_xrefs'].append(omim_pheno_tuple)
        counter += 1
    log.info(f'Found {counter} human disease model OMIM table xrefs in chado.')
    # This returns all MIM IDs, not just those that have a related gene.
    # for hdm in hdm_dict.values():
    #     if hdm['OMIM_pheno_table_xrefs']:
    #         mim_ids = [f'MIM:{i[DBX_ACC]}' for i in hdm['OMIM_pheno_table_xrefs']]
    #         mim_ids.sort()
    #         hdm['related_specific_diseases'] = '|'.join(mim_ids)
    return


# Alternative, unused method for getting related_specific_diseases.
def get_hdm_omim_table_prop(hdm_dict):
    """Retrieve human disease model derived OMIM series table."""
    global CONN
    log.info('Retrieve human disease model derived OMIM series table.')
    fb_hdm_omim_series_query = """
        SELECT DISTINCT hh.humanhealth_id, hh.uniquename, hh.name, hhp.value
        FROM humanhealth hh
        JOIN humanhealthprop hhp ON hh.humanhealth_id = hhp.humanhealth_id
        JOIN cvterm cvt ON cvt.cvterm_id = hhp.type_id
        WHERE hh.is_obsolete IS FALSE
          AND cvt.name = 'derived_disease_tbl';
    """
    ret_hdm_omim_pheno_table_info = connect(fb_hdm_omim_series_query, 'no_query', CONN)
    DB_ID = 0
    # UNAME = 1
    # NAME = 2
    TABLE_TEXT = 3
    counter = 0
    # Example: "[FASPS1](https://omim.org/entry/604348) [PER2](https://omim.org/entry/603426) ..."
    table_line_rgx = r'^\[(.*?)\]\(https://omim.org/entry/(\d+)\)\s*\[(.*?)\]\(https://omim.org/entry/(\d+)\)'
    for row in ret_hdm_omim_pheno_table_info:
        # log.debug(f'Processing row: {row[DB_ID]}, {row[UNAME]}, {row[NAME]}')
        omim_table_lines = row[TABLE_TEXT].split('\n')
        omim_disease_ids = []
        for line in omim_table_lines:
            # log.debug(f'Processing line: {line}')
            try:
                # omim_disease_symbol = re.match(table_rgx, line).group(1)
                omim_disease_id = f'MIM:{re.match(table_line_rgx, line).group(2)}'
                # log.debug(f'Found OMIM disease symbol: {omim_disease_symbol}, ID: {omim_disease_id}')
                omim_disease_ids.append(omim_disease_id)
            except AttributeError:
                log.debug(f'No match for line: {line}')
                pass
        omim_disease_ids.sort()
        hdm_dict[row[DB_ID]]['related_specific_diseases'] = '|'.join(omim_disease_ids)
        counter += 1
    log.info(f'Found {counter} human disease model derived OMIM series tables chado.')
    return


def build_hdm_gene_dict():
    """Create a lookup of HDM-relevant genes."""
    global CONN
    log.info('Create a lookup of HDM-relevant genes.')
    # 1. Build initial gene list.
    fb_hdm_relevant_gene_query = """
        SELECT DISTINCT f.feature_id, f.uniquename, f.name, o.abbreviation
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN humanhealth_feature hhf ON hhf.feature_id = f.feature_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$';
    """
    ret_hdm_relevant_genes = connect(fb_hdm_relevant_gene_query, 'no_query', CONN)
    DB_ID = 0
    UNAME = 1
    NAME = 2
    ORG_ABBR = 3
    counter = 0
    hdm_relevant_gene_dict = {}
    for row in ret_hdm_relevant_genes:
        gene_dict = {
            'feature_id': row[DB_ID],
            'uniquename': row[UNAME],
            'name': row[NAME],
            'org_abbr': row[ORG_ABBR],
            'omim_xref': None,
            'hgnc_xref': None,
        }
        hdm_relevant_gene_dict[row[DB_ID]] = gene_dict
        counter += 1
    log.info(f'Found {counter} HDM-relevant genes in chado.')
    # 2. Fold in OMIM xrefs.
    fb_gene_omim_xref_query = """
        SELECT DISTINCT f.feature_id, dbx.accession, dbx.description
        FROM feature f
        JOIN feature_dbxref fdbx ON fdbx.feature_id = f.feature_id
        JOIN dbxref dbx ON dbx.dbxref_ID = fdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND fdbx.is_current IS TRUE
          AND db.name = 'OMIM_GENE';
    """
    ret_omim_gene_xrefs = connect(fb_gene_omim_xref_query, 'no_query', CONN)
    DB_ID = 0
    ACC = 1
    DESC = 2
    omim_counter = 0
    for row in ret_omim_gene_xrefs:
        if row[DB_ID] in hdm_relevant_gene_dict.keys():
            if row[DESC]:
                omim_xref_tuple = (row[DB_ID], row[ACC], row[DESC])
            else:
                omim_xref_tuple = (row[DB_ID], row[ACC], '')
            hdm_relevant_gene_dict[row[DB_ID]]['omim_xref'] = omim_xref_tuple
            omim_counter += 1
    log.info(f'Added {omim_counter} OMIM_GENE xrefs from chado.')
    # 3. Fold in HGNC xrefs.
    fb_gene_hgnc_xref_query = """
        SELECT DISTINCT f.feature_id, dbx.accession, dbx.description
        FROM feature f
        JOIN feature_dbxref fdbx ON fdbx.feature_id = f.feature_id
        JOIN dbxref dbx ON dbx.dbxref_ID = fdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND fdbx.is_current IS TRUE
          AND db.name = 'HGNC';
    """
    ret_hgnc_gene_xrefs = connect(fb_gene_hgnc_xref_query, 'no_query', CONN)
    DB_ID = 0
    ACC = 1
    DESC = 2
    hgnc_counter = 0
    for row in ret_hgnc_gene_xrefs:
        if row[DB_ID] in hdm_relevant_gene_dict.keys():
            if row[DESC]:
                hgnc_xref_tuple = (row[DB_ID], row[ACC], row[DESC])
            else:
                hgnc_xref_tuple = (row[DB_ID], row[ACC], '')
            hdm_relevant_gene_dict[row[DB_ID]]['hgnc_xref'] = hgnc_xref_tuple
            hgnc_counter += 1
    log.info(f'Added {hgnc_counter} hgnc_GENE xrefs from chado.')
    return hdm_relevant_gene_dict


def get_hdm_genes(hdm_dict, hdm_relevant_gene_dict):
    """Retrieve human disease model genes."""
    global CONN
    log.info('Retrieve human disease model genes.')
    fb_hdm_gene_query = """
        SELECT DISTINCT hh.humanhealth_id, f.feature_id
        FROM humanhealth hh
        JOIN humanhealth_feature hhf ON hhf.humanhealth_id = hh.humanhealth_id
        JOIN humanhealth_featureprop hhfp ON hhfp.humanhealth_feature_id = hhf.humanhealth_feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = hhfp.type_id
        JOIN feature f ON f.feature_id = hhf.feature_id
        JOIN organism o ON o.organism_id = f.organism_id
        WHERE hh.is_obsolete IS FALSE
          AND f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$';
    """
    ret_hdm_gene_info = connect(fb_hdm_gene_query, 'no_query', CONN)
    HDM_DB_ID = 0
    GENE_DB_ID = 1
    counter = 0
    for row in ret_hdm_gene_info:
        hdm_dict[row[HDM_DB_ID]]['implicated_gene_feature_ids'].append(row[GENE_DB_ID])
        counter += 1
    log.info(f'Found {counter} HDM-human gene associations in chado.')
    ACC = 1
    GENE_NAME = 2
    for hdm in hdm_dict.values():
        gene_by_org_dict = {
            'Hsap': [],
            'Dmel': [],
            'other': [],
        }
        # Sort genes by organism, and gather xrefs.
        for feature_id in hdm['implicated_gene_feature_ids']:
            gene = hdm_relevant_gene_dict[feature_id]
            try:
                gene_by_org_dict[gene['org_abbr']].append(gene)
            except KeyError:
                gene_by_org_dict['other'].append(gene)
            if gene['org_abbr'] == 'Hsap' and gene['omim_xref']:
                hdm['OMIM_gene_xrefs'].append(gene['omim_xref'])
            if gene['org_abbr'] == 'Hsap' and gene['hgnc_xref']:
                hdm['HGNC_gene_xrefs'].append(gene['hgnc_xref'])
        # Format OMIM and HGNC xrefs.
        if hdm['OMIM_gene_xrefs']:
            hdm['OMIM_gene_ID'] = '|'.join([f'MIM:{i[ACC]}' for i in hdm['OMIM_gene_xrefs']])
            hdm['OMIM_gene_name'] = ' | '.join([i[GENE_NAME] for i in hdm['OMIM_gene_xrefs']])
        if hdm['HGNC_gene_xrefs']:
            hdm['HGNC_gene_ID'] = '|'.join([f'HGNC:{i[ACC]}' for i in hdm['HGNC_gene_xrefs']])
            hdm['HGNC_gene_name'] = ' | '.join([i[GENE_NAME] for i in hdm['HGNC_gene_xrefs']])
        # Format implicated genes by organism.
        slots_by_org = {
            'Hsap': 'implicated_human_gene',
            'Dmel': 'implicated_Dmel_gene',
            'other': 'implicated_other_gene',
        }
        for org_abbr, slot in slots_by_org.items():
            if gene_by_org_dict[org_abbr]:
                hdm[slot] = ' | '.join([f'{i["uniquename"]};{i["name"]}' for i in gene_by_org_dict[org_abbr]])
            else:
                hdm[slot] = None
    return


def get_hdm_do_terms(hdm_dict):
    """Retrieve human disease model DO terms."""
    global CONN
    log.info('Retrieve human disease model DO terms.')
    fb_hdm_doid_query = """
        SELECT DISTINCT hh.humanhealth_id, 'DOID:'||dbx.accession, cvt.name
        FROM humanhealth hh
        JOIN humanhealth_cvterm hhcvt ON hhcvt.humanhealth_id = hh.humanhealth_id
        JOIN cvterm cvt ON cvt.cvterm_id = hhcvt.cvterm_id
        JOIN cv ON cv.cv_id = cvt.cv_id
        JOIN dbxref dbx ON dbx.dbxref_id = cvt.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE hh.is_obsolete IS FALSE
          AND cvt.is_obsolete = 0
          AND cv.name = 'disease_ontology'
          AND db.name = 'DOID'
        ORDER BY cvt.name;
    """
    ret_hdm_doid_info = connect(fb_hdm_doid_query, 'no_query', CONN)
    DB_ID = 0
    DOID = 1
    CVTERM = 2
    counter = 0
    for row in ret_hdm_doid_info:
        doid_tuple = (row[DB_ID], row[DOID], row[CVTERM])
        hdm_dict[row[DB_ID]]['DO_cvterms'].append(doid_tuple)
        counter += 1
    log.info(f'Found {counter} human disease model DOID annotations in chado.')
    for hdm in hdm_dict.values():
        if hdm['DO_cvterms']:
            hdm['DO_ID'] = '|'.join([i[DOID] for i in hdm['DO_cvterms']])
            hdm['DO_name'] = '|'.join([i[CVTERM] for i in hdm['DO_cvterms']])
    return


def get_external_links(hdm_dict):
    """Retrieve HDM xrefs."""
    global CONN
    log.info('Retrieve HDM xrefs.')
    fb_hdm_xref_query = """
        SELECT DISTINCT hh.humanhealth_id, db.urlprefix||dbx.accession
        FROM humanhealth hh
        JOIN humanhealth_dbxref hhdbx ON hhdbx.humanhealth_id = hh.humanhealth_id
        JOIN dbxref dbx ON dbx.dbxref_id = hhdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE hh.is_obsolete IS FALSE
          AND hhdbx.is_current IS TRUE
          AND db.name NOT IN ('OMIM_GENE', 'OMIM_PHENOTYPE', 'OMIM_series', 'FlyBase', 'HGNC', 'BDSC_HD')
        ORDER BY db.urlprefix||dbx.accession;
    """
    ret_hdm_xref_info = connect(fb_hdm_xref_query, 'no_query', CONN)
    DB_ID = 0
    URL = 1
    counter = 0
    for row in ret_hdm_xref_info:
        hdm_dict[row[DB_ID]]['xref_urls'].append(row[URL])
        counter += 1
    log.info(f'Found {counter} HDM xrefs in chado.')
    for hdm in hdm_dict.values():
        hdm['external_links'] = '|'.join(hdm['xref_urls'])
    return


def get_hdm_props(hdm_dict):
    """Retrieve human health disease model props (various types)."""
    global CONN
    log.info('Retrieve human health disease model props (various types).')
    hhprops = {
        'description_overview': 'dros_model_overview',
        'description_symptoms': 'phenotype_description',
        'description_genetics': 'genetics_description',
        'description_cellular': 'cellular_description',
        'description_molecular': 'molecular_description',
    }
    for slot, prop_type in hhprops.items():
        fb_hdm_prop_query = f"""
            SELECT DISTINCT hh.humanhealth_id, hhp.value
            FROM humanhealth hh
            JOIN humanhealthprop hhp ON hhp.humanhealth_id = hh.humanhealth_id
            JOIN cvterm t ON t.cvterm_id = hhp.type_id
            WHERE hh.is_obsolete IS FALSE
              AND t.name = '{prop_type}';
        """
        ret_hdm_prop_info = connect(fb_hdm_prop_query, 'no_query', CONN)
        DB_ID = 0
        PROP_VALUE = 1
        counter = 0
        for row in ret_hdm_prop_info:
            cleaned_text = clean_free_text(row[PROP_VALUE])
            hdm_dict[row[DB_ID]][slot].append(cleaned_text)
            counter += 1
        log.info(f'Found {counter} {slot} annotations for human health disease models in chado.')
    for hdm in hdm_dict.values():
        for slot in hhprops.keys():
            if hdm[slot]:
                hdm[slot].sort()
                hdm[slot] = '|'.join(hdm[slot])
            else:
                hdm[slot] = None
    return


def get_hdm_omim_bdsc_links(hdm_dict):
    """Retrieve human disease model BDSC links."""
    global CONN
    log.info('Retrieve human disease model BDSC links.')
    fb_hdm_bdsc_link_query = """
        SELECT DISTINCT hh.humanhealth_id, db.urlprefix||dbx.description
        FROM humanhealth hh
        JOIN humanhealth_dbxref hhdbx ON hhdbx.humanhealth_id = hh.humanhealth_id
        JOIN humanhealth_dbxrefprop hhdbxp ON hhdbxp.humanhealth_dbxref_id = hhdbx.humanhealth_dbxref_id
        JOIN cvterm cvt ON cvt.cvterm_id = hhdbxp.type_id
        JOIN dbxref dbx ON dbx.dbxref_id = hhdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE hh.is_obsolete IS FALSE
          AND hhdbx.is_current IS TRUE
          AND db.name = 'BDSC_HD'
          AND cvt.name = 'data_link_bdsc';
    """
    ret_hdm_bdsc_link_info = connect(fb_hdm_bdsc_link_query, 'no_query', CONN)
    DB_ID = 0
    URL = 1
    counter = 0
    for row in ret_hdm_bdsc_link_info:
        hdm_dict[row[DB_ID]]['bdsc_links'].append(row[URL])
        counter += 1
    log.info(f'Found {counter} HDM BDSC HD links in chado.')
    for hdm in hdm_dict.values():
        hdm['BDSC_link'] = '|'.join(hdm['bdsc_links'])
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
