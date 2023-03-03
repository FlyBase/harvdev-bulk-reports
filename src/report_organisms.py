# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script that retrieves organism chado info and prints to TSV file.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_organisms.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_organisms.py -v -c /path/to/config.cfg

"""

import argparse
# import configparser
# import csv
# import datetime
# import logging
# import os
# import pickle
# import psycopg2
# import re
# import sys
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    add_unique_info, connect, set_up_db_reading    # other useful functions: add_list_info, add_unique_dict_info
)
# from harvdev_utils.psycopg_functions.sql_queries import (
#     current_feat_symbol_sgmls, current_feat_fullname_sgmls, feat_symbol_synonyms, feat_fullname_synonyms,
#     feat_secondary_fbids, orgid_abbr, orgid_genus, indirect_rel_features, rel_features, rel_dmel_features,
#     featureprops, feat_cvterm_cvtprop    # current_features, feat_id_symbol_sgml
# )
# from harvdev_utils.char_conversions import *

# Global variables for the output file. Header order will match list order below.
report_label = 'organism_list'
report_title = 'FlyBase organisms report'
header_list = [
        'genus',
        'species',
        'abbreviation',
        'common_name',
        'NCBI_taxon_ID',
        'drosophilid?'
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
parser.add_argument('-i', '--input_filename', help='Input TSV file.', required=False)
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
input_filename = args.input_filename


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    organism_info = get_organism_info(conn)
    organism_info = get_taxon_ids(organism_info, conn)
    organism_info = get_drosophilid_info(organism_info, conn)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = list(organism_info.values())
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of specific data types.
def get_organism_info(db_connection):
    """Generate an organism_id-keyed dict of chado organism table entries.

    Args:
        arg1 (psycopg2.extensions.connection)

    Returns:
        An organism_id-keyed dict of chado organisms.

    """
    log.info('Querying database for organism info.')
    # Get initial list of organisms.
    fb_organism_query = """
        SELECT DISTINCT organism_id,
                        abbreviation,
                        genus,
                        species,
                        common_name
        FROM organism
        ORDER BY genus, species;
        """
    ret_organism_info = connect(fb_organism_query, 'no_query', db_connection)
    log.info('Found {} organisms.'.format(len(ret_organism_info)))
    organism_dict = {}
    for row in ret_organism_info:
        organism_dict[row[0]] = {
                                    'abbreviation': row[1],
                                    'genus': row[2],
                                    'species': row[3],
                                    'common_name': row[4]
                                }

    return organism_dict


def get_taxon_ids(organism_dict, db_connection):
    """Add NCBI taxon IDs to an organism_id-keyed dict of organisms.

    Args:
        arg1 (dict): An organism_id-keyed dict of organisms.
        arg2 (psycopg2.extensions.connection)

    Returns:
        An organism_id-keyed dict of chado organisms with NCBI taxon IDs added.

    """
    log.info('Querying database for NCBI taxon ID info.')
    # Get NCBI taxon IDs.
    fb_taxon_id_query = """
        SELECT DISTINCT o.organism_id,
                        dbx.accession
        FROM organism o
        JOIN organism_dbxref odbx ON odbx.organism_id = o.organism_id
        JOIN dbxref dbx ON dbx.dbxref_id = odbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE odbx.is_current = true and
              db.name = 'NCBITaxon';
        """
    organism_dict = add_unique_info(organism_dict, 'NCBI_taxon_ID', db_connection, fb_taxon_id_query)

    return organism_dict


def get_drosophilid_info(organism_dict, db_connection):
    """Add boolean Drosophilid info to an organism_id-keyed dict of organisms.

    Args:
        arg1 (dict): An organism_id-keyed dict of organisms.
        arg2 (psycopg2.extensions.connection)

    Returns:
        An organism_id-keyed dict of chado organisms with boolean Drosophilid information added.

    """
    log.info('Querying database for Drosophilid info.')
    # Get Drosophilid info.
    fb_drosophilid_query = """
        SELECT DISTINCT o.organism_id,
                        'y'
        FROM organism o
        JOIN organismprop op ON op.organism_id = o.organism_id
        JOIN cvterm cvt ON cvt.cvterm_id = op.type_id
        WHERE cvt.name = 'taxgroup' and
              op.value = 'drosophilid';
        """
    organism_dict = add_unique_info(organism_dict, 'drosophilid?', db_connection, fb_drosophilid_query)

    return organism_dict


if __name__ == "__main__":
    main()
