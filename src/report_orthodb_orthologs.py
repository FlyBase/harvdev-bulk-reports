#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Retrieve OrthoDB Drosophila orthologs for Dmel genes.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_orthodb_orthologs.py [-h] [-v VERBOSE] [-l LOCAL]

Example:
    python report_orthodb_orthologs.py -v -l

"""

import argparse
import configparser
import csv
import datetime
import logging
import os
import psycopg2
import sys
from harvdev_utils.psycopg_functions import (
    connect
)

report_name = 'dmel_orthologs_in_drosophila_species'
report_title = 'FlyBase OrthoDB Drosophila ortholog report'

parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-v', '--verbose', action='store_true', help='DEBUG-level logging.', required=False)
parser.add_argument('-l', '--local', action='store_true', help='Use local credentials.', required=False)
args = parser.parse_args()

local = args.local
if local is True:
    # Run in a local python virtual environment.
    config = configparser.ConfigParser()
    config.read('/data/credentials/reporting/config.cfg')
    database_host = config['default']['Server']
    database = config['default']['Database']
    username = config['default']['User']
    password = config['default']['PGPassword']
    annotation_release = config['default']['AnnotationRelease']
    database_release = config['default']['Release']
    output_filename = './' + report_name + '_' + database + '.tsv'
    log_filename = './' + report_name + '_' + database + '.log'
else:
    # Run within docker in GoCD pipeline.
    database_host = os.environ['SERVER']
    database = os.environ['DATABASE']
    username = os.environ['USER']
    password = os.environ['PGPASSWORD']
    annotation_release = os.environ['ANNOTATIONRELEASE']
    database_release = os.environ['RELEASE']
    output_filename = '/src/output/' + report_name + '_' + database + '.tsv'
    log_filename = '/src/logs/' + report_name + '_' + database + '.log'

verbose = args.verbose
if verbose is True:
    logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.DEBUG)
else:
    logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.INFO)
log = logging.getLogger(__name__)
sys.stdout = open(log_filename, 'a')

# Establish the db connection.
conn_string = "host=%s dbname=%s user=%s password='%s'" % (database_host, database, username, password)
conn = psycopg2.connect(conn_string)
log.info('Successful connection to database {} on host {}.'.format(database, database_host))


def main():
    """Retrieve and print out Drosophila OrthoDB orthologs for Dmel genes."""
    log.info('TIME: {}. Started main function.'.format(now()))
    orthodb_info = get_orthodb_info()
    to_export_as_tsv = create_tsv_data_structure()
    to_export_as_tsv['data'] = process_orthodb_info(orthodb_info)
    tsv_dump(to_export_as_tsv, output_filename)
    conn.close()
    log.info('TIME: {}. Ended main function.'.format(now()))


def now():
    """Create a simple human-readable timestamp for logging."""
    now = datetime.datetime.now().strftime("%H:%M:%S")
    return now


def get_orthodb_info():
    """Return a list of tuples containing Dmel gene, Dmel featureloc and Dros orthologous gene.

    Returns:
        OrthoDB info as a list of tuples.
    """
    log.info('TIME: {}. Querying database for orthodb info.'.format(now()))
    fb_orthodb_query = """
        SELECT f1.uniquename, f1.name, src1.uniquename, fl1.fmin+1||'..'||fl1.fmax, fl1.strand,
        f2.uniquename, f2.name, f2.feature_id, fr.value
        FROM feature f1
        JOIN feature_relationship fr ON fr.object_id = f1.feature_id
        JOIN feature_relationshipprop frp ON frp.feature_relationship_id = fr.feature_relationship_id
        JOIN organism o1 ON o1.organism_id = f1.organism_id
        JOIN featureloc fl1 ON fl1.feature_id = f1.feature_id
        JOIN feature src1 ON src1.feature_id = fl1.srcfeature_id
        JOIN cvterm cvtsrc ON cvtsrc.cvterm_id = src1.type_id
        JOIN cvterm cvtfr ON cvtfr.cvterm_id = fr.type_id
        JOIN feature f2 ON f2.feature_id = fr.subject_id
        JOIN organism o2 ON o2.organism_id = f2.organism_id
        WHERE f1.is_obsolete = false and
              f1.uniquename ~ '^FBgn[0-9]{7}$' and
              o1.abbreviation = 'Dmel' and
              src1.is_obsolete = false and
              src1.organism_id = f1.organism_id and
              cvtsrc.name ~ '^golden_path' and
              cvtfr.name = 'orthologous_to' and
              frp.value = 'ORTHODB' and
              f2.is_obsolete = false and
              f2.uniquename ~ '^FB(gn|og)[0-9]{7,10}$' and
              o2.abbreviation in ('Dsim', 'Dsec', 'Dyak', 'Dere', 'Dana', 'Dper', 'Dpse', 'Dwil', 'Dvir', 'Dmoj', 'Dgri', 'Dsuz');
        """
    ret_orthodb_info = connect(fb_orthodb_query, 'no_query', conn)
    log.info('TIME: {}. Found {} orthodb relationships for Dmel genes.'.format(now(), len(ret_orthodb_info)))
    return ret_orthodb_info


def get_featureloc_info(feature_id):
    """Return scaffold uniquename, fmin, fmax and strand for a given feature_id.

    Args:
        arg1 (integer): the Dmel gene feature_id.

    Returns:
        The featureloc for the gene.
    """
    fb_featureloc_query = """
        SELECT DISTINCT src.uniquename, fl.fmin, fl.fmax, fl.strand
        FROM featureloc fl
        JOIN feature src ON src.feature_id = fl.srcfeature_id
        WHERE src.is_obsolete = false and src.type_id in (204, 553)
          and fl.feature_id = %s;"""
    ret_featureloc_query = connect(fb_featureloc_query, (feature_id,), conn)
    return ret_featureloc_query


def process_orthodb_info(input_data):
    """Take SQL results and return a list of dictionaries for tsv output.

    Args:
        arg1 (list): a list of tuples representing chado query results for orthologs.

    Returns:
        A list of dictionaries representing OrthoDB ortholog info.
    """
    log.info('TIME: {}. Starting to process orthodb info retrieved from database.'.format(now()))
    data_list = []
    FBGN_ID = 0
    GENE_SYMBOL = 1
    SCAFFOLD = 2
    LOCATION = 3
    STRAND = 4
    ORTHO_FBGN_ID = 5
    ORTHO_GENE_SYMBOL = 6
    ORTHO_FEAT_ID = 7
    ORTHODB_GROUP_ID = 8
    for i in input_data:
        orthodb_item = {
            'FBgn_ID': i[FBGN_ID],
            'GeneSymbol': i[GENE_SYMBOL],
            'Arm/Scaffold': i[SCAFFOLD],
            'Location': i[LOCATION],
            'Strand': i[STRAND],
            'Ortholog_FBgn_ID': i[ORTHO_FBGN_ID],
            'Ortholog_GeneSymbol': i[ORTHO_GENE_SYMBOL],
            'Ortholog_feature_id': i[ORTHO_FEAT_ID],
            'OrthoDB_Group_ID': i[ORTHODB_GROUP_ID]
        }
        ortholog_featureloc_info = get_featureloc_info(orthodb_item['Ortholog_feature_id'])
        if ortholog_featureloc_info:
            orthodb_item['Ortholog_Arm/Scaffold'] = ortholog_featureloc_info[0][0]
            fmin = ortholog_featureloc_info[0][1] + 1
            fmax = ortholog_featureloc_info[0][2]
            orthodb_item['Ortholog_Location'] = '{}..{}'.format(fmin, fmax)
            orthodb_item['Ortholog_Strand'] = ortholog_featureloc_info[0][3]
        else:
            orthodb_item['Ortholog_Arm/Scaffold'] = None
            orthodb_item['Ortholog_Location'] = None
            orthodb_item['Ortholog_Strand'] = None
        data_list.append(orthodb_item)
    log.info('TIME: {}. Done processing ortholog info into a list of dictionaries.'.format(now()))
    return data_list


def create_tsv_data_structure():
    """Build a generic dictionary for tsv export."""
    to_export_as_tsv = {}
    to_export_as_tsv['metaData'] = {}
    to_export_as_tsv['metaData']['title'] = report_title
    to_export_as_tsv['metaData']['dateProduced'] = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    to_export_as_tsv['metaData']['database'] = database
    to_export_as_tsv['data'] = []
    return to_export_as_tsv


def tsv_dump(data_object, output_filename):
    """Send Data object to specified output file.

    Args:
        arg1 (list): a list of dictionaries to print to file.
        arg2 (string): filename for the output file.

    Returns:
        Nothing, per se.
    """
    log.info('TIME: {}. Writing data to output tsv file.'.format(now()))
    output_file = open(output_filename, 'w')

    output_file.write('## {}\n'.format(data_object['metaData']['title']))
    output_file.write('## Generated: {}\n'.format(data_object['metaData']['dateProduced']))
    output_file.write('## Using datasource: {}\n'.format(data_object['metaData']['database']))
    output_file.write('##\n## ')
    headers = [
        'FBgn_ID',
        'GeneSymbol',
        'Arm/Scaffold',
        'Location',
        'Strand',
        'Ortholog_FBgn_ID',
        'Ortholog_GeneSymbol',
        'Ortholog_Arm/Scaffold',
        'Ortholog_Location',
        'Ortholog_Strand',
        'OrthoDB_Group_ID'
        ]
    csv_writer = csv.DictWriter(output_file, fieldnames=headers, delimiter='\t', extrasaction='ignore')
    csv_writer.writeheader()

    for data_item in data_object['data']:
        csv_writer.writerow(data_item)

    output_file.write('## Finished {}.'.format(data_object['metaData']['title']))
    log.info('TIME: {}. Done writing data to output file.'.format(now()))
    return


if __name__ == "__main__":
    main()
