#!/usr/bin/env python3

# report_json_template.py
# author
# usage: report_json_template.py [-h] [-v VERBOSE]
################################################################################

import argparse
import configparser
import datetime
import json
import logging
import os
import pickle
import psycopg2
import re
import strict_rfc3339
import sys
# from pprint import pformat
# from harvdev_utils.char_conversions import *

report_name = 'this_report'

parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-v', '--verbose', action='store_true', help='DEBUG-level logging.', required=False)
args = parser.parse_args()

# For running script locally.
# config = configparser.ConfigParser()
# config.read('/data/credentials/reporting/config.cfg')
# database_host = config['connection']['DatabaseHost']
# database = config['connection']['Database']
# username = config['connection']['Username']
# password = config['connection']['Password']
# annotation_release = config['connection']['AnnotationRelease']
# database_release = config['connection']['DatabaseRelease']
# output_filename = '/data/build-public-release/' + database + '/bulk_reports/' + report_name + '_' + database + '.json'
# log_filename = '/data/build-public-release/' + database + '/loading_output/' + report_name + '_' + database + '.log'

# For running script within GoCD pipeline.
database_host = os.environ['SERVER']
database = os.environ['DATABASE']
username = os.environ['USER']
password = os.environ['PGPASSWORD']
annotation_release = os.environ['ANNOTATIONRELEASE']
database_release = os.environ['RELEASE']
output_filename = './output/' + report_name + '_' + database + '.json'
log_filename = './logs/' + report_name + '_' + database + '.log'

verbose = args.verbose
if verbose is True:
    logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.DEBUG)
else:
    logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.INFO)
log = logging.getLogger(__name__)
sys.stdout = open(log_filename, 'a')


def main():

    # Log start time and arguments.
    log.info('TIME: {}. Started main function.'.format(now()))

    # Set this as the "official" time the file was generated.
    the_time = strict_rfc3339.now_to_rfc3339_localoffset()

    # Database connection.
    conn_string = "host=%s dbname=%s user=%s password='%s'" % (database_host, database, username, password)
    conn = psycopg2.connect(conn_string)
    log.info('TIME: {}. Successful connection to database {} on host {}.'.format(now(), database, database_host))

    # Create the JSON data structure and fill in the file metadata.
    to_export_as_json = {}
    to_export_as_json['metaData'] = {}
    to_export_as_json['metaData']['dataProvider'] = 'FlyBase'
    to_export_as_json['metaData']['reportName'] = report_name
    to_export_as_json['metaData']['dateProduced'] = the_time
    to_export_as_json['metaData']['databaseRelease'] = database_release
    to_export_as_json['metaData']['AnnotationRelease'] = annotation_release
    to_export_as_json['data'] = []

    # Sample query.
    fb_sample_query = """
    SELECT f.name, cvt.name
    FROM feature f JOIN cvterm cvt ON (f.type_id = cvt.cvterm_id and f.is_obsolete = false)
    WHERE f.name = 'wg' and f.uniquename ~ '^FBgn[0-9]{7}$'
    ;"""

    # Retrieve data from the database.
    log.info('TIME: {}. Querying database for info.'.format(now()))
    fb_sample_info = connect(fb_sample_query, 'no_query', conn)

    conn.close()

    # Build the list of data elements for export.
    to_export_as_json['data'] = [{
        'gene_name': i[0],
        'feature_type': i[1]}
        for i in fb_sample_info]

    if len(to_export_as_json['data']) == 0:
        log.error('No data to report.')
        raise

    # Generate the output file.
    log.info('TIME: {}. Writing data to output file.'.format(now()))

    with open(output_filename, 'w') as output_file:
        json.dump(to_export_as_json, output_file, sort_keys=True, indent=2, separators=(',', ': '))
        output_file.close()

    log.info('TIME: {}. Done writing data to output file.'.format(now()))

    # Log end of main function.
    log.info('TIME: {}. Ended main function.'.format(now()))


# Compact function for getting time stamps.
def now():
    now = datetime.datetime.now().strftime("%H:%M:%S")
    return now


# Function for db queries.
def connect(sql, query, conn):
    # Return the cursor and use it to perform queries.
    cursor = conn.cursor()
    # Execute the query, with or without some variable.
    if query == 'no_query':
        cursor.execute(sql)
    else:
        cursor.execute(sql, query)
    # Grab the results.
    records = cursor.fetchall()
    # Close the cursor
    cursor.close()
    # Return a list of tuples
    return records


if __name__ == "__main__":
    main()
