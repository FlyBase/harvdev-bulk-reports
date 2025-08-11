# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report InterPro xrefs.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_interpro_xrefs.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_interpro_xrefs.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
REPORT_LABEL = 'gene_interpro_xrefs'
REPORT_TITLE = 'FlyBase Gene InterPro Signatures Report'
HEADER_LIST = [
    'FBgn_ID',
    'FBgn_Symbol',
    'InterPro_Signature',
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
    gene_interpro_list = get_interpro_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = gene_interpro_list
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_interpro_info():
    """Retrieve gene-InterPro xrefs."""
    global CONN
    log.info('Retrieve gene-InterPro xrefs.')
    fb_gene_interpro_query = """
        SELECT DISTINCT f.uniquename, f.name, dbx.accession, dbx.description
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN feature_dbxref fdbx ON fdbx.feature_id = f.feature_id
        JOIN dbxref dbx ON dbx.dbxref_id = fdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND o.abbreviation = 'Dmel'
          AND fdbx.is_current IS TRUE
          AND db.name = 'INTERPRO'
        ORDER BY f.uniquename, dbx.accession;
    """
    ret_gene_interpro_info = connect(fb_gene_interpro_query, 'no_query', CONN)
    GENE_UNAME = 0
    GENE_NAME = 1
    ACC = 2
    DESC = 3
    gene_interpro_list = []
    for row in ret_gene_interpro_info:
        gene_interpro_result = {
            'FBgn_ID': row[GENE_UNAME],
            'FBgn_Symbol': row[GENE_NAME],
            'InterPro_Signature': f'{row[ACC]}|{row[DESC]}',
        }
        gene_interpro_list.append(gene_interpro_result)
    log.info(f'Found {len(ret_gene_interpro_info)} gene-InterPro xrefs in chado.')
    return gene_interpro_list


if __name__ == "__main__":
    main()
