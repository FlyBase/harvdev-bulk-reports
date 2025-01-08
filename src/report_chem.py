# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script for reporting FlyBase chemical entities and related information.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_chem.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_chem.py -v -c /path/to/config.cfg

Notes:
    For 10 chems, 22 synonyms contain a pipe ("|") character; these synonyms
    are excluded since we want to preserve the pipe character as a separator.
    As per FB bulk file convention, only ASCII-names/aliases will be printed.

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
report_label = 'chemicals'
report_title = 'FlyBase chemicals report'
header_list = [
    'FB_id',
    'FB_name',
    'FB_synonyms',
    'InChIKey',
    'PubChem_id',
    'PubChem_name',
    'PubChem_synonyms',
    'ChEBI_id',
    'ChEBI_name',
    'ChEBI_synonyms',
    'ChEBI_definition',
    'ChEBI_roles',
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
    fb_chems = run_chem_queries(conn)
    id_sorted_fb_chems = [fb_chems[i] for i in sorted(fb_chems.keys())]
    # Export the data to file.
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = id_sorted_fb_chems
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def get_fb_chems(db_connection):
    """Generate an FBch ID-keyed dict of chem dicts.

    Args:
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    Returns:
        An FBch ID-keyed dict of chem dicts.

    """
    log.info('Querying database for current FlyBase chemicals.')
    fb_chem_query = """
        SELECT DISTINCT f.uniquename, f.name
        FROM feature f
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBch[0-9]{7}$';
    """
    ret_chem_info = connect(fb_chem_query, 'no_query', db_connection)
    log.info(f'Found {len(ret_chem_info)} chems in chado.')
    ID = 0
    NAME = 1
    fb_chem_dict = {}
    for result in ret_chem_info:
        result_chem_dict = {
            'FB_id': result[ID],
            'FB_name': result[NAME],
            'FB_synonyms': [],
            'InChIKey': [],
            'PubChem_id': None,
            'PubChem_name': None,
            'PubChem_synonyms': [],
            'ChEBI_id': None,
            'ChEBI_name': None,
            'ChEBI_synonyms': [],
            'ChEBI_definition': [],
            'ChEBI_roles': [],
            'internal_list_chebi_pubchem_synonyms': [],
        }
        fb_chem_dict[result[ID]] = result_chem_dict
    return fb_chem_dict


def get_external_ids(fb_chem_dict, db_connection):
    """Add external PubChem and CheBI IDs/names to the FBch ID-keyed dict of FlyBase chemical info.

    Args:
        fb_chem_dict (dict): An FBch ID-keyed dict of chemical dicts.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    """
    log.info('Get PubChem and ChEBI IDs and names.')
    fb_external_query = """
        SELECT DISTINCT f.uniquename, db.name, dbx.accession, dbx.description
        FROM feature f
        JOIN feature_dbxref fdbx ON fdbx.feature_id = f.feature_id
        JOIN dbxref dbx ON dbx.dbxref_id = fdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBch[0-9]{7}$'
          AND fdbx.is_current IS TRUE
          AND db.name IN ('CHEBI', 'PubChem');
    """
    ret_ext_ids = connect(fb_external_query, 'no_query', db_connection)
    chebi_counter = 0
    pubchem_counter = 0
    ID = 0
    DB = 1
    ACC = 2
    NAME = 3
    for result in ret_ext_ids:
        if result[DB] == 'PubChem':
            fb_chem_dict[result[ID]]['PubChem_id'] = f'PubChem:{result[ACC]}'
            fb_chem_dict[result[ID]]['PubChem_name'] = result[NAME]
            pubchem_counter += 1
        elif result[DB] == 'CHEBI':
            fb_chem_dict[result[ID]]['ChEBI_id'] = f'CHEBI:{result[ACC]}'
            fb_chem_dict[result[ID]]['ChEBI_name'] = result[NAME]
            chebi_counter += 1
    log.info(f'Found {pubchem_counter} PubChem IDs for chemicals.')
    log.info(f'Found {chebi_counter} ChEBI IDs for chemicals.')
    return


def get_chebi_pubchem_synonyms(fb_chem_dict, db_connection):
    """Add ChEBI and PubChem-attributed synonyms to the FBch ID-keyed dict of FlyBase chemical info.

    Args:
        fb_chem_dict (dict): An FBch ID-keyed dict of chemical dicts having a "FB_synonyms" key.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    """
    log.info('Get ChEBI and PubChem chemical synonyms.')
    fb_chem_synonym_query = """
        SELECT DISTINCT f.uniquename, s.name, p.uniquename
        FROM feature f
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN pub p ON p.pub_id = fs.pub_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBch[0-9]{7}$'
          AND p.is_obsolete IS FALSE
          AND s.name != f.name
          AND NOT s.name LIKE '%|%'
          AND NOT s.name ILIKE 'pubchem:%'
          AND p.uniquename IN ('FBrf0243181', 'FBrf0243186');
    """
    ret_chem_synonym_info = connect(fb_chem_synonym_query, 'no_query', db_connection)
    log.info(f'Found {len(ret_chem_synonym_info)} ChEBI/PubChem synonyms in chado.')
    ID = 0
    SYNONYM_TEXT = 1
    PUB = 2
    for result in ret_chem_synonym_info:
        fb_chem_dict[result[ID]]['internal_list_chebi_pubchem_synonyms'].append(result[SYNONYM_TEXT])
        if result[PUB] == 'FBrf0243186':
            fb_chem_dict[result[ID]]['PubChem_synonyms'].append(result[SYNONYM_TEXT])
        elif result[PUB] == 'FBrf0243181':
            fb_chem_dict[result[ID]]['ChEBI_synonyms'].append(result[SYNONYM_TEXT])
    return


def get_fb_chem_synonyms(fb_chem_dict, db_connection):
    """Add FlyBase-only synonyms to the FBch ID-keyed dict of FlyBase chemical info.

    Args:
        fb_chem_dict (dict): An FBch ID-keyed dict of chemical dicts having a "FB_synonyms" key.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    """
    log.info('Get FlyBase-only chemical synonyms.')
    fb_chem_synonym_query = """
        SELECT DISTINCT f.uniquename, s.name
        FROM feature f
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN pub p ON p.pub_id = fs.pub_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBch[0-9]{7}$'
          AND p.is_obsolete IS FALSE
          AND s.name != f.name
          AND NOT s.name LIKE '%|%'
          AND NOT s.name ILIKE 'pubchem:%'
          AND NOT p.uniquename IN ('FBrf0243181', 'FBrf0243186');
    """
    ret_chem_synonym_info = connect(fb_chem_synonym_query, 'no_query', db_connection)
    log.info(f'Found {len(ret_chem_synonym_info)} synonyms in chado.')
    ID = 0
    SYNONYM_TEXT = 1
    for result in ret_chem_synonym_info:
        # Exclude any synonyms also attributed to the ChEBI or PubChem references.
        if result[SYNONYM_TEXT] in fb_chem_dict[result[ID]]['internal_list_chebi_pubchem_synonyms']:
            continue
        elif result[SYNONYM_TEXT] == fb_chem_dict[result[ID]]['PubChem_name']:
            continue
        elif result[SYNONYM_TEXT] == fb_chem_dict[result[ID]]['ChEBI_name']:
            continue
        fb_chem_dict[result[ID]]['FB_synonyms'].append(result[SYNONYM_TEXT])
    return


def get_chebi_definitions(fb_chem_dict, db_connection):
    """Add ChEBI definitions to the FBch ID-keyed dict of FlyBase chemical info.

    Args:
        fb_chem_dict (dict): An FBch ID-keyed dict of chemical dicts.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    """
    # log.info('Get ChEBI definitions.')
    # fb_chebi_definition_query = """
    #     SELECT DISTINCT f.uniquename, fp.value
    #     FROM feature f
    #     JOIN featureprop fp ON fp.feature_id = f.feature_id
    #     JOIN cvterm cvt ON cvt.cvterm_id = fp.type_id
    #     WHERE f.is_obsolete IS FALSE
    #       AND f.uniquename ~ '^FBch[0-9]{7}$'
    #       AND cvt.name = 'description'
    #       AND fp.value ~ '^ChEBI: ';
    # """
    # ret_chebi_definitions = connect(fb_chebi_definition_query, 'no_query', db_connection)
    # log.info(f'Found {len(ret_chebi_definitions)} ChEBI definitions.')
    # ID = 0
    # DEFINITION = 1
    # for result in ret_chebi_definitions:
    #     fb_chem_dict[result[ID]]['ChEBI_definition'].append(result[DEFINITION])
    return


def get_inchikeys(fb_chem_dict, db_connection):
    """Add InChiKey IDs to the FBch ID-keyed dict of FlyBase chemical info.

    Args:
        fb_chem_dict (dict): An FBch ID-keyed dict of chemical dicts.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    """
    # log.info('Get InChiKey IDs.')
    # fb_inchikey_query = """
    #     SELECT DISTINCT f.uniquename, fp.value
    #     FROM feature f
    #     JOIN featureprop fp ON fp.feature_id = f.feature_id
    #     JOIN cvterm cvt ON cvt.cvterm_id = fp.type_id
    #     WHERE f.is_obsolete IS FALSE
    #       AND f.uniquename ~ '^FBch[0-9]{7}$'
    #       AND cvt.name = 'inchikey';
    # """
    # ret_inchikeys = connect(fb_inchikey_query, 'no_query', db_connection)
    # log.info(f'Found {len(ret_inchikeys)} InChiKeys.')
    # ID = 0
    # INCHIKEY = 1
    # for result in ret_inchikeys:
    #     fb_chem_dict[result[ID]]['InChIKey'].append(result[INCHIKEY])
    return


def get_chem_featureprops(fb_chem_dict, db_connection):
    """Add featureprops to the FBch ID-keyed dict of FlyBase chemical info: InChiKey, ChEBI definition and roles.

    Args:
        fb_chem_dict (dict): An FBch ID-keyed dict of chemical dicts.
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    """
    log.info('Get chemical featureprops: InChiKey, ChEBI definition, ChEBI roles.')
    fprop_types = {
        'InChIKey': " = 'inchikey' ",
        'ChEBI_definition': " = 'description' ",
        'ChEBI_roles': " IN ('biological_role', 'application_role') ",
    }
    for chem_attribute, prop_type_filter in fprop_types.items():
        log.debug(f'Get info for this attribute: {chem_attribute}')
        fprop_query = f"""
            SELECT DISTINCT f.uniquename, fp.value
            FROM feature f
            JOIN featureprop fp ON fp.feature_id = f.feature_id
            JOIN cvterm cvt ON cvt.cvterm_id = fp.type_id
            WHERE f.is_obsolete IS FALSE
            AND f.uniquename ~ '^FBch[0-9]+$'
            AND fp.value IS NOT NULL
            AND fp.value != ''
            AND cvt.name {prop_type_filter};
        """
        fprop_counter = 0
        ret_fprops = connect(fprop_query, 'no_query', db_connection)
        log.info(f'Found {len(ret_fprops)} initial {chem_attribute} featureprops.')
        ID = 0
        FPROP_VALUE = 1
        for result in ret_fprops:
            if chem_attribute == 'ChEBI_definition' and result[FPROP_VALUE] and not result[FPROP_VALUE].startswith('ChEBI: '):
                continue
            fb_chem_dict[result[ID]][chem_attribute].append(result[FPROP_VALUE])
            fprop_counter += 1
        log.info(f'For attribute {chem_attribute}, added {fprop_counter} featureprop entries.')
    return


def process_chem_dict(fb_chem_dict):
    """Process dict of chemical into final output desired."""
    log.info('Process chemical info for print output.')
    for chem in fb_chem_dict.values():
        # First, remove internal fields.
        del chem['internal_list_chebi_pubchem_synonyms']
        # Second, convert lists to a string with pipes separating distinct entries.
        for chem_key, chem_attribute in chem.items():
            if type(chem_attribute) is list:
                chem[chem_key] = ' | '.join(sorted(set(chem_attribute)))
    return


def run_chem_queries(db_connection):
    """Generate a fully populated FBch ID-keyed dict of FlyBase chemical info.

    Args:
        db_connection (psycopg2.extensions.connection): The object used to interact with the database.

    Returns:
        fb_chem_dict (dict): The FBch-keyed dict of chemical info.

    """
    log.info('Generate full chemical dict.')
    fb_chem_dict = get_fb_chems(db_connection)
    get_external_ids(fb_chem_dict, db_connection)
    get_chebi_pubchem_synonyms(fb_chem_dict, db_connection)
    get_fb_chem_synonyms(fb_chem_dict, db_connection)
    get_chem_featureprops(fb_chem_dict, db_connection)
    process_chem_dict(fb_chem_dict)
    return fb_chem_dict


if __name__ == "__main__":
    main()
