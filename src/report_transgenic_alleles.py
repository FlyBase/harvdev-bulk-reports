# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report transgenic alleles info.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_transgenic_alleles.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_transgenic_alleles.py -v -c /path/to/config.cfg

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
report_label = 'transgenic_construct_descriptions'
report_title = 'FlyBase D. melanogaster transgenic construct descriptions report'
header_list = [
    'Component Allele (symbol)',
    'Component Allele (id)',
    'Transgenic Construct (symbol)',
    'Transgenic Construct (id)',
    'Transgenic Product class (term)',
    'Transgenic Product class (id)',
    'Regulatory region (symbol)',
    'Regulatory region (id)',
    'Encoded product/tool (symbol)',
    'Encoded product/tool (id)',
    'Tagged with (symbol)',
    'Tagged with (id)',
    'Also carries (symbol)',
    'Also carries (id)',
    'Description (text)',
    'Description (supporting reference)',
    'Stocks (number)',
    # 'Stock (list)'
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
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))


def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    database_info = get_database_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(database_info)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def get_alleles():
    """Query chado for all alleles.

    Returns:
        allele_dict (dict): A feature_id-keyed dict of alleles (represented by dicts).

    """
    global conn
    log.info('Query chado for all alleles.')
    fb_alleles_query = """
        SELECT DISTINCT f.feature_id, f.name, f.uniquename
        FROM feature f
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBal[0-9]{7}$'
        ORDER BY f.uniquename;
    """
    ret_fb_alleles = connect(fb_alleles_query, 'no_query', conn)
    FEAT_ID = 0
    SYMBOL = 1
    FB_CURIE = 2
    fb_allele_dict = {}
    for result in ret_fb_alleles:
        fb_allele_dict[result[FEAT_ID]] = {
            'Component Allele (symbol)': result[SYMBOL],
            'Component Allele (id)': result[FB_CURIE],
            'is_transgenic': False,
            'Gene (symbol)': None,
            'Gene (id)': None,
            'Transgenic Construct (symbol)': [],
            'Transgenic Construct (id)': [],
            'Transgenic Product class (term)': [],
            'Transgenic Product class (id)': [],
            'Regulatory region (symbol)': [],
            'Regulatory region (id)': [],
            'Encoded product/tool (symbol)': [],
            'Encoded product/tool (id)': [],
            'Tagged with (symbol)': [],
            'Tagged with (id)': [],
            'Also carries (symbol)': [],
            'Also carries (id)': [],
            'Description (text)': [],
            'Description (supporting reference)': [],
            'Stocks (number)': 0,
            'Stock (list)': [],
        }
    log.info(f'Found {len(ret_fb_alleles)} current alleles in chado.')
    return fb_allele_dict


def get_transgenic_constructs(fb_allele_dict):
    """Get transgenic constructs."""
    global conn
    log.info('Get transgenic constructs.')
    fb_transgenic_constructs_query = """
        SELECT DISTINCT a.feature_id, c.name, c.uniquename
        FROM feature a
        JOIN feature_relationship fr ON fr.subject_id = a.feature_id
        JOIN feature c ON c.feature_id = fr.object_id
        JOIN cvterm t ON t.cvterm_id = fr.type_id
        WHERE a.is_obsolete IS FALSE
          AND a.uniquename ~ '^FBal[0-9]{7}$'
          AND c.is_obsolete IS FALSE
          AND c.uniquename ~ '^FBtp[0-9]{7}$'
          AND t.name = 'associated_with'
        ORDER BY c.name;
    """
    ret_fb_transgenic_constructs = connect(fb_transgenic_constructs_query, 'no_query', conn)
    FEAT_ID = 0
    CONS_NAME = 1
    CONS_CURIE = 2
    counter = 0
    for result in ret_fb_transgenic_constructs:
        fb_allele_dict[result[FEAT_ID]]['is_transgenic'] = True
        fb_allele_dict[result[FEAT_ID]]['Transgenic Construct (symbol)'].append(result[CONS_NAME])
        fb_allele_dict[result[FEAT_ID]]['Transgenic Construct (id)'].append(result[CONS_CURIE])
        counter += 1
    log.info(f'Found {counter} current constructs for current alleles.')
    return


def get_allele_genes(fb_allele_dict):
    """Get allele genes."""
    global conn
    log.info('Get allele genes.')
    fb_allele_genes_query = """
        SELECT DISTINCT a.feature_id, g.name, g.uniquename
        FROM feature a
        JOIN feature_relationship fr ON fr.subject_id = a.feature_id
        JOIN feature g ON g.feature_id = fr.object_id
        JOIN cvterm t ON t.cvterm_id = fr.type_id
        WHERE a.is_obsolete IS FALSE
          AND a.uniquename ~ '^FBal[0-9]{7}$'
          AND g.is_obsolete IS FALSE
          AND g.uniquename ~ '^FBgn[0-9]{7}$'
          AND t.name = 'alleleof';
    """
    ret_fb_allele_genes = connect(fb_allele_genes_query, 'no_query', conn)
    FEAT_ID = 0
    GENE_NAME = 1
    GENE_CURIE = 2
    counter = 0
    for result in ret_fb_allele_genes:
        if fb_allele_dict[result[FEAT_ID]]['is_transgenic'] is False:
            continue
        fb_allele_dict[result[FEAT_ID]]['Gene (symbol)'] = result[GENE_NAME]
        fb_allele_dict[result[FEAT_ID]]['Gene (id)'] = result[GENE_CURIE]
        counter += 1
    log.info(f'Found {counter} parental genes for current transgenic alleles in chado.')
    return


def get_allele_transgenic_product_classes(fb_allele_dict):
    """Get allele classes."""
    global conn
    log.info('Get allele classes.')
    fb_allele_classes_query = """
        SELECT DISTINCT f.feature_id, cvt.name, db.name||':'||dbx.accession
        FROM feature f
        JOIN feature_cvterm fcvt ON fcvt.feature_id = f.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
        JOIN dbxref dbx ON dbx.dbxref_id = cvt.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        JOIN feature_cvtermprop fcvtp ON fcvtp.feature_cvterm_id = fcvt.feature_cvterm_id
        JOIN cvterm t ON t.cvterm_id = fcvtp.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBal[0-9]{7}$'
          AND fcvt.is_not IS FALSE
          AND cvt.is_obsolete = 0
          AND t.name = 'transgenic_product_class'
          AND db.name = 'SO'
        ORDER BY cvt.name;
    """
    ret_fb_allele_classes = connect(fb_allele_classes_query, 'no_query', conn)
    FEAT_ID = 0
    TERM_NAME = 1
    TERM_CURIE = 2
    counter = 0
    for result in ret_fb_allele_classes:
        if fb_allele_dict[result[FEAT_ID]]['is_transgenic'] is False:
            continue
        fb_allele_dict[result[FEAT_ID]]['Transgenic Product class (term)'].append(result[TERM_NAME])
        fb_allele_dict[result[FEAT_ID]]['Transgenic Product class (id)'].append(result[TERM_CURIE])
        counter += 1
    log.info(f'Found {counter} allele class annotations for current transgenic alleles in chado.')
    return


def get_allele_descriptions(fb_allele_dict):
    """Get allele descriptions."""
    global conn
    log.info('Get allele descriptions.')
    fb_allele_descriptions_query = """
        SELECT DISTINCT f.feature_id, fp.value, STRING_AGG(p.uniquename, '+')
        FROM feature f
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fp.type_id
        JOIN featureprop_pub fpp ON fpp.featureprop_id = fp.featureprop_id
        JOIN pub p ON p.pub_id = fpp.pub_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBal[0-9]{7}$'
          AND cvt.name IN ('aminoacid_rep', 'molecular_info', 'nucleotide_sub')
        GROUP BY f.feature_id, fp.value;
    """
    ret_fb_allele_descriptions = connect(fb_allele_descriptions_query, 'no_query', conn)
    FEAT_ID = 0
    DESC_TEXT = 1
    PUB_ID = 2
    counter = 0
    for result in ret_fb_allele_descriptions:
        if fb_allele_dict[result[FEAT_ID]]['is_transgenic'] is False:
            continue
        desc_text = clean_free_text(result[DESC_TEXT])
        fb_allele_dict[result[FEAT_ID]]['Description (text)'].append(desc_text)
        fb_allele_dict[result[FEAT_ID]]['Description (supporting reference)'].append(result[PUB_ID])
        counter += 1
    log.info(f'Found {counter} allele descriptions for current transgenic alleles in chado.')
    return


def get_allele_stock_info(fb_allele_dict):
    """Get allele stock info."""
    global conn
    log.info('Get allele stock info.')
    fb_allele_stocks_query = """
        SELECT DISTINCT f.feature_id, fp.value
        FROM feature f
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fp.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBal[0-9]{7}$'
          AND cvt.name ~ '^derived_stock';
    """
    ret_fb_allele_stocks = connect(fb_allele_stocks_query, 'no_query', conn)
    FEAT_ID = 0
    STOCK_TEXT = 1
    stock_id_rgx = r'FBst[0-9]{7}'
    for result in ret_fb_allele_stocks:
        if fb_allele_dict[result[FEAT_ID]]['is_transgenic'] is False:
            continue
        stock_ids = re.findall(stock_id_rgx, result[STOCK_TEXT])
        fb_allele_dict[result[FEAT_ID]]['Stock (list)'].extend(stock_ids)
    counter = 0
    for allele in fb_allele_dict.values():
        allele['Stock (list)'] = list(set(allele['Stock (list)']))
        allele['Stocks (number)'] = len(allele['Stock (list)'])
        counter += allele['Stocks (number)']
    log.info(f'Found {counter} stocks for current transgenic alleles in chado.')
    return


def get_direct_component_info(fb_allele_dict):
    """Get allele component info (direct)."""
    global conn
    log.info('Get allele component info (direct).')
    component_associations = {
        'has_reg_region': 'Regulatory region',
        'encodes_tool': 'Encoded product/tool',
        'tagged_with': 'Tagged with',
        'carries_tool': 'Also carries',
    }
    for asso_type, slot_name in component_associations.items():
        log.info(f'Get direct "{asso_type}" info.')
        fb_allele_component_query = f"""
            SELECT DISTINCT a.feature_id, component.name, component.uniquename
            FROM feature a
            --LEFT OUTER JOIN featureprop fp ON fp.feature_id = a.feature_id
            --  AND fp.type_id IN (SELECT cvterm_id FROM cvterm WHERE name = 'propagate_transgenic_uses')
            JOIN feature_relationship fr ON fr.subject_id = a.feature_id
              AND fr.type_id IN (SELECT cvterm_id FROM cvterm WHERE name = '{asso_type}')
            JOIN feature component ON component.feature_id = fr.object_id
            WHERE a.is_obsolete IS FALSE
              AND a.uniquename ~ '^FBal[0-9]{{7}}$'
            --  AND fp.value IS NULL
              AND component.is_obsolete IS FALSE
              AND component.uniquename ~ '^FB[a-z]{{2}}[0-9]{{7,10}}$'
            ORDER BY component.name;
        """
        ret_fb_allele_components = connect(fb_allele_component_query, 'no_query', conn)
        FEAT_ID = 0
        COMPONENT_NAME = 1
        COMPONENT_CURIE = 2
        counter = 0
        for result in ret_fb_allele_components:
            if fb_allele_dict[result[FEAT_ID]]['is_transgenic'] is False:
                continue
            symbol_slot_name = f'{slot_name} (symbol)'
            id_slot_name = f'{slot_name} (id)'
            fb_allele_dict[result[FEAT_ID]][symbol_slot_name].append(result[COMPONENT_NAME])
            fb_allele_dict[result[FEAT_ID]][id_slot_name].append(result[COMPONENT_CURIE])
            counter += 1
        log.info(f'Found {counter} current {asso_type} direct component associations for current non-transgenic alleles.')
    return


def get_database_info():
    """Run chado db queries in sequence."""
    global conn
    log.info('Query database.')
    allele_dict = get_alleles()
    get_transgenic_constructs(allele_dict)
    get_allele_genes(allele_dict)
    get_allele_transgenic_product_classes(allele_dict)
    get_allele_descriptions(allele_dict)
    get_allele_stock_info(allele_dict)
    get_direct_component_info(allele_dict)
    return allele_dict.values()


def process_database_info(input_data):
    """Process a list of data elements for TSV output.

    Args:
        input_data (list): A list of dicts representing current FB alleles.

    Returns:
        A list of dictionaries representing classical/insertion alleles.
    """
    log.info('Starting to process current alleles info retrieved from database.')
    data_list = []
    keep_counter = 0
    skip_counter = 0
    for i in input_data:
        if i['is_transgenic'] is False:
            skip_counter += 1
            continue
        keep_counter += 1
        # Minor modification if no encode product/tool is found.
        if not i['Encoded product/tool (id)']:
            i['Encoded product/tool (symbol)'] = i['Gene (symbol)']
            i['Encoded product/tool (id)'] = i['Gene (id)']
        for k, v in i.items():
            if type(v) is list:
                i[k] = '|'.join(v)
        data_list.append(i)
    log.info('Done processing current transgenic alleles into a list of dictionaries.')
    log.info(f'Exporting {keep_counter} current transgenic alleles.')
    log.info(f'Skipping {skip_counter} current classical/insertion alleles.')
    return data_list


if __name__ == "__main__":
    main()
