# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report human orthologs.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_human_orthologs.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_human_orthologs.py -v -c /path/to/config.cfg

"""

import argparse
# import configparser
# import datetime
# import re
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)
from harvdev_utils.char_conversions import *

# Global variables for the output file. Header order will match list order below.
report_label = 'dmel_human_orthologs_disease'
report_title = 'FlyBase D. melanogaster-Human Orthologs and associated diseases report'
header_list = [
    'Dmel_gene_ID',
    'Dmel_gene_symbol',
    'Human_gene_HGNC_ID',
    'Human_gene_OMIM_ID',
    'Human_gene_symbol',
    'DIOPT_score',
    'OMIM_Phenotype_IDs',
    'OMIM_Phenotype_IDs[name]'
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
    database_info = get_database_info(conn)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(database_info)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def get_database_info(db_connection):
    """Retrieve a list of tuples containing Dmel-Hsap ortholog pairs.

    Args:
        arg1 (psycopg2.extensions.connection)

    Returns:
        A list of tuples representing SQL query output.

    """
    log.info('Querying database for ortholog info.')
    fb_paralog_query = """
        
    ;"""
    ret_paralog_info = connect(fb_paralog_query, 'no_query', db_connection)
    log.info('Found {} paralogous relationships among Dmel genes.'.format(len(ret_paralog_info)))

    return ret_paralog_info


def process_database_info(input_data):
    """Convert a list of SQL results into a list of dictionaries for TSV output.

    Args:
        arg1 (list): A list of tuples representing SQL query output.

    Returns:
        A list of dictionaries representing, in this case, paralog information.
    """
    log.info('Starting to process paralog info retrieved from database.')
    data_list = [{
        'FBgn_ID': i[0],
        'GeneSymbol': i[1],
        'Arm/Scaffold': i[2],
        'Location': i[3],
        'Strand': i[4],
        'Paralog_FBgn_ID': i[5],
        'Paralog_GeneSymbol': i[6],
        'Paralog_Arm/Scaffold': i[7],
        'Paralog_Location': i[8],
        'Paralog_Strand': i[9],
        'DIOPT_score': i[10].count(',') + 1    # Does this count commas correctly?
        }
        for i in input_data]
    log.info('Done processing paralog info into a list of dictionaries.')

    return data_list


if __name__ == "__main__":
    main()
