# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Retrieve SO terms for localized Dmel genes and print data to TSV file.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu, Julie Agapite jagapite@morgan.harvard.edu

Usage:
    report_gene_so_terms.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_gene_so_terms.py -v -c /path/to/config.cfg

"""

from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
report_label = 'dmel_gene_sequence_ontology_annotations'
report_title = 'FlyBase Sequence Ontology annotations for localized D. melanogaster genes'
header_list = [
    'gene_primary_id',
    'gene_symbol',
    'so_term_name',
    'so_term_id'
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


def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')

    database_info = get_database_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(database_info)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()

    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of specific data types.
def get_database_info():
    """Retrieve a list of localized dmel genes and their SO annotations.

    Args:
        None.

    Returns:
        A list of tuples representing SQL query output.

    """
    log.info('Querying database for dmel localized genes and their SO annotations.')

    # Driver query.
    fb_gene_so_annotation_query = """
        SELECT DISTINCT f.uniquename,
               s.name,
               cvt.name,
               db.name||':'||dbx.accession
        FROM feature f
        JOIN featureloc fl ON fl.feature_id = f.feature_id
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN cvterm cvts ON cvts.cvterm_id = s.type_id
        JOIN feature_cvterm fcvt ON fcvt.feature_id = f.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
        JOIN cv ON cv.cv_id = cvt.cv_id
        JOIN dbxref dbx ON dbx.dbxref_id = cvt.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE f.is_obsolete = false and
              f.is_analysis = false and
              f.uniquename ~ '^FBgn[0-9]{7}$' and
              fs.is_current = true and
              fs.is_internal = false and
              cvts.name = 'symbol' and
              cv.name = 'SO' and
              db.name = 'SO' and
              cvt.name ~ 'gene';"""
    ret_so_annotation_info = connect(fb_gene_so_annotation_query, 'no_query', conn)
    log.info('Found {} gene-SO annotations for localized Dmel genes.'.format(len(ret_so_annotation_info)))

    return ret_so_annotation_info


def process_database_info(input_data):
    """Convert a list of SQL results into a list of dictionaries for TSV output.

    Args:
        arg1 (list): A list of tuples representing SQL query output.

    Returns:
        A list of dictionaries representing gene-SO annotations.
    """
    log.info('Starting to process gene-SO annotations retrieved from database.')
    UNIQUENAME = 0
    GENE_SYMBOL = 1
    CVTERM_NAME = 2
    CVTERM_ID = 3
    data_list = []
    for result in input_data:
        annotation_dict = {
            'gene_primary_id': result[UNIQUENAME],
            'gene_symbol': result[GENE_SYMBOL].replace('<up>', '[').replace('</up>', ']'),    # For su(w[a]) gene.
            'so_term_name': result[CVTERM_NAME],
            'so_term_id': result[CVTERM_ID]
        }
        data_list.append(annotation_dict)
    log.info('Done processing gene-SO annotation data into a list of dictionaries.')

    return data_list


if __name__ == "__main__":
    main()
