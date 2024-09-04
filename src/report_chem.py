# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script for reporting FlyBase chemical entities and related information.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_chem.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_chem.py -v -c /path/to/config.cfg

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
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))


def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    fb_chems = run_chem_queries(conn)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = fb_chems.values()
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def get_fb_chems():
    """Generate an FBch ID-keyed dict of chem dicts.

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
    ret_chem_info = connect(fb_chem_query, 'no_query', conn)
    log.info(f'Found {len(ret_chem_info)} chems in chado.'.format())
    ID = 0
    NAME = 1
    fb_chem_dict = {}
    for result in ret_chem_info:
        result_chem_dict = {
            'FB_id': result[ID],
            'FB_name': result[NAME],
            'FB_synonyms': [],
            'InChIKey': None,
            'PubChem_id': None,
            'PubChem_name': None,
            'PubChem_synonyms': {},
            'ChEBI_id': None,
            'ChEBI_name': None,
            'ChEBI_synonyms': [],
            'ChEBI_definition': None,
            'ChEBI_roles': [],
        }
        fb_chem_dict[result[ID]] = result_chem_dict
    return fb_chem_dict


def run_chem_queries(db_connection):
    """Generate a fully populated FBch ID-keyed dict of FlyBase chemical info.

    Args:
        db_connection (): The psycopg2 db connection entity used to interact with the database.

    Returns:
        fb_chem_dict (dict): The FBch-keyed dict of chemical info.

    """
    log.info('Generate full chemical dict.')
    log.debug(f'Have this conn type: {type(db_connection)}.')
    fb_chem_dict = get_fb_chems(db_connection)
    return fb_chem_dict


if __name__ == "__main__":
    main()
