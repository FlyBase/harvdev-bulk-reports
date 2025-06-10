# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script for reporting FlyBase disease implicated diseases and related information.

Author(s):
    Ian longden ilongden@morgan.harvard.edu

Usage:
    report_div.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_div.py -v -c /path/to/config.cfg

Notes:

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
report_label = 'divs'
report_title = 'FlyBase disease implicated variants report'
header_list = [
    'FB_name',
    'FBrf',
    'Clinvar_id',
    'UP_SP_variant_id',
    'FB_synonyms',
    'HumanHealth_id',
    'Allele_id',
]

# Proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
database = set_up_dict['database']
database_release = set_up_dict['database_release']
output_dir = set_up_dict['output_dir']
output_filename = set_up_dict['output_filename']
log = set_up_dict['log']
conn = set_up_dict['conn']
the_time = set_up_dict['the_time']

parser = argparse.ArgumentParser(description='inputs')
args, extra_args = parser.parse_known_args()
log.info(f'Parsing args specific to this script; ignoring these: {extra_args}')


def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    # Get chem data and sort it by FBch ID.
    # div type is 142324
    type_id = 142324
    fb_divs = run_div_queries(conn, type_id)
    id_sorted_fb_chems = [fb_divs[i] for i in sorted(fb_divs.keys())]
    # Export the data to file.
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = id_sorted_fb_chems
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def get_fb_divs(db_connection, type_id):
    """Generate an div dict,.

    Args:
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    Returns:
        An din uniquename -keyed dict of div dicts.

    """
    log.info('Querying database for current FlyBase divs.')

    fb_chem_query = f"""
        select f.feature_id, f.uniquename, p.uniquename
        from feature f ,feature_pub fp, pub p
        where f.feature_id = fp.feature_id
              and fp.pub_id = p.pub_id
              and f.type_id = {type_id}
              and f.is_obsolete IS FALSE
    """
    ret_chem_info = connect(fb_chem_query, 'no_query', db_connection)
    log.info(f'Found {len(ret_chem_info)} chems in chado.')
    ID = 0
    NAME = 1
    FBRF = 2
    fb_div_dict = {}
    for result in ret_chem_info:
        result_chem_dict = {
            'FB_name': result[NAME],
            'FBrf': result[FBRF],
            'Clinvar_id': None,
            'UP_SP_variant_id': None,
            'FB_synonyms': [],
            'HumanHealth_id': [],
            'Allele_id': [],
        }
        fb_div_dict[result[NAME]] = result_chem_dict
    return fb_div_dict


def get_alleles(fb_div_dict, db_connection, type_id):
    """ Add human health and alle data"""
    fb_relation_query = f"""
        select object.uniquename, subject.uniquename
        from feature_relationship fr, feature subject, feature object, cvterm cvt
        where fr.subject_id = subject.feature_id
         and  fr.object_id  = object.feature_id
         and cvt.cvterm_id = subject.type_id
         and object.type_id = {type_id}
         and object.is_obsolete IS FALSE
         and cvt.name = 'allele'
    """
    ID = 0
    SUB_ID = 1
    ret_rela_ids = connect(fb_relation_query, 'no_query', db_connection)
    for result in ret_rela_ids:
        fb_div_dict[result[ID]]['Allele_id'].append(result[SUB_ID])


def get_humanhealth(fb_div_dict, db_connection, type_id):
    """
    Get humanhealth relationships.
    """
    hh_query = f"""
        select f.uniquename, hh.uniquename
        from humanhealth_feature hhf, humanhealth hh, feature f
        where hhf.humanhealth_id = hh.humanhealth_id
        and hhf.feature_id = f.feature_id
        and f.type_id = {type_id}
        and f.is_obsolete IS FALSE
        and hh.is_obsolete = 'f' 
    """
    ID = 0
    SUB_ID = 1
    ret_rela_ids = connect(hh_query, 'no_query', db_connection)
    for result in ret_rela_ids:
        fb_div_dict[result[ID]]['HumanHealth_id'].append(result[SUB_ID])
        print(f"{result[ID]}, {result[SUB_ID]}")


def get_external_ids(fb_chem_dict, db_connection, type_id):
    """Add external div IDs/names to the div ID-keyed dict of FlyBase chemical info.

    Args:
        fb_chem_dict (dict): An div ID-keyed dict of div dicts.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    """
    log.info('Get dbxrefs.')
    fb_external_query = f"""
        SELECT DISTINCT f.uniquename, db.name, dbx.accession, dbx.description
        FROM feature f
        JOIN feature_dbxref fdbx ON fdbx.feature_id = f.feature_id
        JOIN dbxref dbx ON dbx.dbxref_id = fdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE f.is_obsolete IS FALSE
          AND f.type_id = {type_id}
          AND fdbx.is_current IS TRUE;
    """
    ret_ext_ids = connect(fb_external_query, 'no_query', db_connection)
    ID = 0
    DB = 1
    ACC = 2
    NAME = 3
    for result in ret_ext_ids:
        if result[DB] == 'ClinVar':
            fb_chem_dict[result[ID]]['Clinvar_id'] = f'{result[DB]}:{result[ACC]}'
        elif result[DB] == 'UP_SP_variant':
            fb_chem_dict[result[ID]]['UP_SP_variant_id'] = f'{result[DB]}:{result[ACC]}'
    return


def get_div_synonyms(fb_div_dict, db_connection, type_id):
    """Add FlyBase-only synonyms to the div  ID-keyed dict of FlyBase divinfo.

    Args:
        fb_div_dict (dict): An FBch ID-keyed dict of div dicts having a "FB_synonyms" key.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.
        type_id (int): The type ID.
    """
    log.info('Get FlyBase-only div synonyms.')
    fb_div_synonym_query = f"""
        SELECT DISTINCT f.uniquename, s.name
        FROM feature f
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN pub p ON p.pub_id = fs.pub_id
        WHERE f.is_obsolete IS FALSE
          AND f.type_id = {type_id}
          AND p.is_obsolete IS FALSE
          AND s.name != f.name
          AND NOT s.name LIKE '%|%';
    """
    ret_chem_synonym_info = connect(fb_div_synonym_query, 'no_query', db_connection)
    log.info(f'Found {len(ret_chem_synonym_info)} synonyms in chado.')
    ID = 0
    SYNONYM_TEXT = 1
    for result in ret_chem_synonym_info:
        fb_div_dict[result[ID]]['FB_synonyms'].append(result[SYNONYM_TEXT])


def process_div_dict(fb_div_dict):
    """Process dict of div into final output desired."""
    log.info('Process chemical info for print output.')
    for chem in fb_div_dict.values():
        # convert lists to a string with pipes separating distinct entries.
        for chem_key, chem_attribute in chem.items():
            if type(chem_attribute) is list:
                chem[chem_key] = ' | '.join(sorted(set(chem_attribute)))
    return


def run_div_queries(db_connection, type_id):
    """Generate a fully populated div ID-keyed dict of FlyBase divinfo.

    Args:
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    Returns:
        fb_div_dict (dict): The div-keyed dict of div info.

    """
    log.info('Generate full disease implicated disease dict.')
    fb_div_dict = get_fb_divs(db_connection, type_id)
    get_external_ids(fb_div_dict, db_connection, type_id)
    get_div_synonyms(fb_div_dict, db_connection, type_id)
    get_alleles(fb_div_dict, db_connection, type_id)
    get_humanhealth(fb_div_dict, db_connection, type_id)
    process_div_dict(fb_div_dict)
    return fb_div_dict


if __name__ == "__main__":
    main()
