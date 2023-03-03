#!/usr/bin/env python3

# report_paralogs.py
# author: Gil dos Santos
# usage: report_paralogs.py [-h] [-v VERBOSE] [-l LOCAL]
################################################################################

import argparse
import configparser
import csv
import datetime
import logging
import os
# import pickle
import psycopg2
# import re
import sys
# from harvdev_utils.char_conversions import *

# Global variables for the output file. Header order will match list order below.
report_name = 'dmel_paralogs'
report_title = 'FlyBase DIOPT Drosophila melanogaster paralog report'
header_list = [
    'FBgn_ID',
    'GeneSymbol',
    'Arm/Scaffold',
    'Location',
    'Strand',
    'Paralog_FBgn_ID',
    'Paralog_GeneSymbol',
    'Paralog_Arm/Scaffold',
    'Paralog_Location',
    'Paralog_Strand',
    'DIOPT_score'
    ]

# Process input parameters.
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-v', '--verbose', action='store_true', help='DEBUG-level logging.', required=False)
parser.add_argument('-l', '--local', action='store_true', help='Use local credentials.', required=False)
args = parser.parse_args()

# Environment and output filenames/locations.
local = args.local
if local is True:
    # Run in a local python virtual environment.
    config = configparser.ConfigParser()
    config.read('/data/credentials/reporting/config.cfg')
    database_host = config['connection']['DatabaseHost']
    database = config['connection']['Database']
    username = config['connection']['Username']
    password = config['connection']['Password']
    annotation_release = config['connection']['AnnotationRelease']
    database_release = config['connection']['DatabaseRelease']
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

# Specify verbosity of log file.
verbose = args.verbose
if verbose is True:
    logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.DEBUG)
else:
    logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.INFO)
log = logging.getLogger(__name__)
sys.stdout = open(log_filename, 'a')


# Basic process of the script.
def main():
    log.info('TIME: {}. Started main function.'.format(now()))
    conn = establish_db_connection()
    paralog_info = get_database_info(conn)
    data_to_export_as_tsv = create_tsv_data_structure()
    data_to_export_as_tsv['data'] = process_database_info(paralog_info)
    tsv_dump(data_to_export_as_tsv, headers=header_list)
    conn.close()
    log.info('TIME: {}. Ended main function.'.format(now()))


# BELOW: Functions specific retrieval and processing of specific data types.
def get_database_info(db_connection):
    """Returns a list of tuples containing paralog pairs and related info.
       Required libraries: psycopg2, datetime.
       Required functions: connect(), now().
       Required global variables: none.
       Input: a database connection (db_connection).
       Output: paralog info (list of tuples)."""
    log.info('TIME: {}. Querying database for paralog info.'.format(now()))
    fb_paralog_query = """
        SELECT f1.uniquename, f1.name, src1.uniquename, fl1.fmin+1||'..'||fl1.fmax, fl1.strand,
        f2.uniquename, f2.name, src2.uniquename, fl2.fmin+1||'..'||fl2.fmax, fl2.strand, fr.value
        FROM feature f1
        JOIN feature_relationship fr ON (fr.subject_id = f1.feature_id)
        JOIN feature_relationshipprop frp ON (frp.feature_relationship_id = fr.feature_relationship_id)
        JOIN organism o1 ON (o1.organism_id = f1.organism_id)
        JOIN featureloc fl1 ON (fl1.feature_id = f1.feature_id)
        JOIN feature src1 ON (src1.feature_id = fl1.srcfeature_id)
        JOIN cvterm cvtfr ON (cvtfr.cvterm_id = fr.type_id)
        JOIN feature f2 ON (f2.feature_id = fr.object_id)
        JOIN organism o2 ON (o2.organism_id = f2.organism_id)
        JOIN featureloc fl2 ON (fl2.feature_id = f2.feature_id)
        JOIN feature src2 ON (src2.feature_id = fl2.srcfeature_id)
        WHERE f1.is_obsolete = false and f1.uniquename ~ '^FBgn[0-9]{7}$' and o1.abbreviation = 'Dmel'
          and f2.is_obsolete = false and f2.uniquename ~ '^FBgn[0-9]{7}$' and o2.abbreviation = 'Dmel'
          and cvtfr.name = 'paralogous_to' and frp.value = 'DIOPT'
          and src1.is_obsolete = false and src1.organism_id = f1.organism_id and src1.type_id = 553
          and src2.is_obsolete = false and src2.organism_id = f2.organism_id and src2.type_id = 553;"""
    ret_paralog_info = connect(fb_paralog_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} paralogous relationships among Dmel genes.'.format(now(), len(ret_paralog_info)))

    return ret_paralog_info


def process_database_info(input_data):
    """Takes SQL results and returns a list of dictionaries for tsv output.
       Required libraries: datetime.
       Required functions: now().
       Required global variables: none.
       Input: a list of tuples representing paralog info.
       Output: a list of dictionaries representing paralog info."""
    log.info('TIME: {}. Starting to process paralog info retrieved from database.'.format(now()))
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
    log.info('TIME: {}. Done processing paralog info into a list of dictionaries.'.format(now()))

    return data_list


# BELOW: All generic functions for any bulk file script: database reading, file writing, etc.
def now():
    """Creates a simple human-readable timestamp for logging.
       Required libraries: datetime.
       Required functions: none.
       Required global variables: none.
       Input: none.
       Output: human-readable timestamp (string)."""
    now = datetime.datetime.now().strftime("%H:%M:%S")

    return now


def establish_db_connection():
    """Establishes the database connection.
       Required libraries: psycopg2, datetime.
       Required functions: now().
       Required global variables: database_host, database, username, password.
       Input: none.
       Output: a database connection."""
    conn_string = "host=%s dbname=%s user=%s password='%s'" % (database_host, database, username, password)
    db_connection = psycopg2.connect(conn_string)
    log.info('TIME: {}. Successful connection to database {} on host {}.'.format(now(), database, database_host))

    return db_connection


def connect(sql, query, db_connection):
    """Retrieves information from the database, with optional argument for the query.
       Required libraries: psycopg2.
       Required functions: none.
       Required global variables: none.
       Input: an sql query (sql), query variable (either 'no_query' or a tuple), and a db connection (db_connection).
       Output: db results (a list of tuples)."""
    cursor = db_connection.cursor()
    if query == 'no_query':         # If SQL query lacks a variable.
        cursor.execute(sql)
    else:
        cursor.execute(sql, query)  # If SQL query has a variable.
    records = cursor.fetchall()
    cursor.close()

    return records


def create_tsv_data_structure():
    """Builds a generic dictionary for tsv export.
       Required libraries: datetime.
       Required functions: now().
       Required global variables: the_time, database_release, report_title.
       Input: none.
       Output: generic Alliance JSON data structure (dictionary)."""
    to_export_as_tsv = {}
    to_export_as_tsv['metaData'] = {}
    to_export_as_tsv['metaData']['title'] = report_title
    to_export_as_tsv['metaData']['dateProduced'] = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    to_export_as_tsv['metaData']['database'] = database
    to_export_as_tsv['data'] = []

    return to_export_as_tsv


def tsv_dump(data_object, **kwargs):
    """Sends a data object to specified output file (global variable).
       Required libraries: datetime, csv.
       Required functions: now().
       Required global variables: output_filename.
       Input: a data dictionary with 'metadata' key and a list of data dictionaries under 'data' key.
       Output: writes a flat tsv file, returns nothing itself."""
    log.info('TIME: {}. Writing data to output tsv file.'.format(now()))
    output_file = open(output_filename, 'w')

    # Supply a list of headers to specify the print order.
    # Make sure each element in header list matches a data dictionary key.
    if 'headers' in kwargs.keys():
        headers = kwargs['headers']
    # Otherwise, it just takes the dictionary keys from the first data object.
    else:
        headers = data_object['data'][0].keys()

    output_file.write('## {}\n'.format(data_object['metaData']['title']))
    output_file.write('## Generated: {}\n'.format(data_object['metaData']['dateProduced']))
    output_file.write('## Using datasource: {}\n'.format(data_object['metaData']['database']))
    output_file.write('##\n## ')

    csv_writer = csv.DictWriter(output_file, fieldnames=headers, delimiter='\t', extrasaction='ignore')
    csv_writer.writeheader()

    for data_item in data_object['data']:
        csv_writer.writerow(data_item)

    output_file.write('## Finished {}.'.format(data_object['metaData']['title']))
    log.info('TIME: {}. Done writing data to output file.'.format(now()))

    return


if __name__ == "__main__":
    main()
