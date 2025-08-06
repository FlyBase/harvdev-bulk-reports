# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report gene model annotation comments.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_gene_model_annotation_comments.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_gene_model_annotation_comments.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.char_conversions import clean_free_text
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
REPORT_LABEL = 'gene_model_annotation_comments'
REPORT_TITLE = 'FlyBase Gene Model Annotation Comments Report'
HEADER_LIST = [
    'FB_id',
    'Symbol',
    'Annotation_Status',
    'Annotation_Comment',
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
    gene_comment_list = get_gene_annotation_comments()
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = gene_comment_list
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_gene_annotation_comments():
    """Retrieve gene model annotation comments."""
    global CONN
    log.info('Retrieve gene model annotation comments.')
    fb_gene_comment_query = """
        SELECT DISTINCT f.uniquename, f.name, fp1.value, fp2.value
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN featureprop fp1 ON fp1.feature_id = f.feature_id
        JOIN cvterm c1 ON c1.cvterm_id = fp1.type_id AND c1.name = 'derived_gene_model_status'
        LEFT OUTER JOIN featureprop fp2 ON fp2.feature_id = f.feature_id
        LEFT OUTER JOIN cvterm c2 ON c2.cvterm_id = fp2.type_id AND c2.name = 'comment'
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND o.abbreviation = 'Dmel'
        ORDER BY f.uniquename, fp2.value;
    """
    ret_gene_comment_info = connect(fb_gene_comment_query, 'no_query', CONN)
    UNAME = 0
    NAME = 1
    STATUS = 2
    COMMENT = 3
    gene_comment_list = []
    counter = 0
    for row in ret_gene_comment_info:
        if row[COMMENT]:
            try:
                cleaned_comment = clean_free_text(row[COMMENT])
            except ValueError as e:
                cleaned_comment = row[COMMENT]
                log.error(f'For {row[UNAME]}, could not clean "{row[COMMENT]}": {e}')
        else:
            cleaned_comment = ''
        result = {
            'FB_id': row[UNAME],
            'Symbol': row[NAME],
            'Annotation_Status': row[STATUS],
            'Annotation_Comment': cleaned_comment,
        }
        gene_comment_list.append(result)
        counter += 1
    log.info(f'Found {counter} gene model annotation comments in chado.')
    return gene_comment_list


if __name__ == "__main__":
    main()
