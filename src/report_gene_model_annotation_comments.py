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
    'Annotation_Comments',
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
    gene_dict = get_genes()
    get_gene_annotation_comments(gene_dict)
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = process_database_info(gene_dict)
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_genes():
    """Retrieve Dmel genes from chado."""
    global CONN
    log.info('Retrieve Dmel genes from chado.')
    fb_gene_query = """
        SELECT f.uniquename, f.name, fp.value
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm c ON c.cvterm_id = fp.type_id AND c.name = 'derived_gene_model_status'
        WHERE f.is_obsolete IS FALSE
          AND o.abbreviation = 'Dmel'
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
        ORDER BY f.uniquename;
    """
    ret_gene_info = connect(fb_gene_query, 'no_query', CONN)
    UNAME = 0
    NAME = 1
    MODEL_STATUS = 2
    gene_dict = {}
    for row in ret_gene_info:
        this_gene = {
            'FB_id': row[UNAME],
            'Symbol': row[NAME],
            'Annotation_Status': row[MODEL_STATUS],
            'Annotation_Comments': [],
        }
        gene_dict[row[UNAME]] = this_gene
    log.info(f'Found {len(gene_dict)} genes in chado.')
    return gene_dict


def get_gene_annotation_comments(gene_dict):
    """Retrieve gene model annotation comments."""
    global CONN
    log.info('Retrieve gene model annotation comments.')
    fb_gene_comment_query = """
        SELECT DISTINCT f.uniquename, fp.value
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm c ON c.cvterm_id = fp.type_id AND c.name = 'comment'
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND o.abbreviation = 'Dmel';
    """
    ret_gene_comment_info = connect(fb_gene_comment_query, 'no_query', CONN)
    UNAME = 0
    COMMENT = 1
    for row in ret_gene_comment_info:
        if row[COMMENT]:
            try:
                cleaned_comment = clean_free_text(row[COMMENT])
            except KeyError as e:
                cleaned_comment = row[COMMENT]
                log.error(f'For {row[UNAME]}, could not clean "{row[COMMENT]}": {e}')
        else:
            cleaned_comment = ''
        gene_dict[row[UNAME]]['Annotation_Comments'].append(cleaned_comment)
    log.info(f'Found {len(ret_gene_comment_info)} gene model annotation comments in chado.')
    return


def process_database_info(input_data):
    """Convert the gene dict into a list of data elements."""
    log.info('Convert the gene dict to a list of data elements.')
    data_list = []
    counter = 0
    for i in input_data.values():
        for key, value in i.items():
            if type(value) is list:
                if not value:
                    i[key] = ''
                else:
                    value = list(set(value))
                    value.sort()
                    i[key] = '|'.join(value)
        data_list.append(i)
        counter += 1
    log.info(f'Sending {counter} split system combination entries to the export file.')
    return data_list


if __name__ == "__main__":
    main()
