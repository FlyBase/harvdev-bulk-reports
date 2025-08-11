# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report chemical synonyms.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_chem_synonyms.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_chem_synonyms.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
REPORT_LABEL = 'chem_synonyms'
REPORT_TITLE = 'FlyBase Chemical Synonyms Report'
HEADER_LIST = [
    'Publication_ID',
    'FB_Chemical_ID',
    'FB_Chemical_Name',
    'Author Synonym',
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
    chem_synonym_list = get_chem_synonym_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = chem_synonym_list
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_chem_synonym_info():
    """Retrieve chemical synonyms."""
    global CONN
    log.info('Retrieve chemical synonyms.')
    fb_chem_synonym_query = """
        SELECT DISTINCT f.uniquename, f.name, p.uniquename, s.name
        FROM feature f
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN pub p ON p.pub_id = fs.pub_id
        JOIN cvterm cvt ON cvt.cvterm_id = p.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBch[0-9]{7}$'
          AND fs.is_current IS FALSE
          AND p.is_obsolete IS FALSE
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
          AND cvt.name != 'database';
    """
    ret_chem_synonym_info = connect(fb_chem_synonym_query, 'no_query', CONN)
    CHEM_ID = 0
    CHEM_NAME = 1
    PUB_ID = 2
    ALIAS = 3
    chem_synonym_list = []
    for row in ret_chem_synonym_info:
        chem_syno_result = {
            'Publication_ID': row[PUB_ID],
            'FB_Chemical_ID': row[CHEM_ID],
            'FB_Chemical_Name': row[CHEM_NAME],
            'Author Synonym': row[ALIAS],
        }
        chem_synonym_list.append(chem_syno_result)
    log.info(f'Found {len(ret_chem_synonym_info)} chemical synonyms in chado.')
    return chem_synonym_list


if __name__ == "__main__":
    main()
