#!/usr/bin/env python3

# report_pheno_for_drsc.py
# author: Gil dos Santos
# modified from "report_allele_phenotype.pl" script
# usage: report_pheno_for_drsc.py [-h] [-v VERBOSE] [-l LOCAL]
################################################################################

import argparse
import configparser
import csv
import datetime
import logging
import os
# import pickle
import psycopg2
import re
import sys
from harvdev_utils.char_conversions import sgml_to_plain_text

# Global variables for the output file. Header order will match list order below.
report_name = 'pheno_data_for_drsc'
report_title = 'FlyBase Drosophila melanogaster RNAi and CRISPR phenotype for DRSC report'
header_list = [
    'FBgn',
    'gene_symbol',
    'FBal',
    'allele_symbol',
    'reagent_id',
    'reagent_type',
    'reagent_source',
    'driver',
    'phenotype',
    'FBrf',
    'PMID'
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
    pheno_info = get_pheno_info(conn)
    data_to_export_as_tsv = create_tsv_data_structure()
    data_to_export_as_tsv['data'] = process_pheno_info(pheno_info, conn)
    tsv_dump(data_to_export_as_tsv, headers=header_list)
    conn.close()
    log.info('TIME: {}. Ended main function.'.format(now()))


# BELOW: Functions for retrieval and processing of specific data types.
def get_pheno_info(db_connection):
    """Returns a list of tuples containing phenotype info.
       Required libraries: psycopg2, datetime.
       Required functions: connect(), now().
       Required global variables: none.
       Input: a database connection (db_connection).
       Output: phenotype info (list of tuples)."""
    log.info('TIME: {}. Querying database for phenotype info.'.format(now()))

    fb_pheno_query = """
        SELECT DISTINCT f.uniquename, f.name, fp.value, p.uniquename
        FROM feature f
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fp.type_id
        JOIN featureprop_pub fpp ON fpp.featureprop_id = fp.featureprop_id
        JOIN pub p ON p.pub_id = fpp.pub_id
        WHERE f.is_obsolete = false and f.is_analysis = false and f.uniquename ~ '^FBal[0-9]{7}$'
          and cvt.name in ('derived_pheno_class','derived_pheno_manifest');"""
    ret_pheno_info = connect(fb_pheno_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} derived phenotype featureprops for alleles.'.format(now(), len(ret_pheno_info)))

    return ret_pheno_info


def get_rnai_info(db_connection):
    """Returns all RNAi alleles and their related collections, if any, from the database.
       Required libraries: psycopg2, datetime.
       Required functions: connect(), now().
       Required global variables: none.
       Input: a database connection (db_connection).
       Output: an {FBal ID: collection} set of FBal IDs."""
    log.info('TIME: {}. Querying database for RNAi alleles.'.format(now()))

    # Get all RNAi alleles
    rnai_fbal_query = """
        SELECT DISTINCT f.uniquename, f.name
        FROM feature f
        JOIN feature_cvterm fcvt ON f.feature_id = fcvt.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
        JOIN cv ON cv.cv_id = cvt.cv_id
        WHERE f.is_obsolete = false and f.is_analysis = false
          and f.uniquename ~ '^FBal[0-9]{7}$'
          and cv.name = 'SO'
          and cvt.name = 'RNAi_reagent';"""
    ret_rnai_fbal = connect(rnai_fbal_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} RNAi alleles.'.format(now(), len(ret_rnai_fbal)))

    # Get FBal-FBsf dict.
    fbal_fbsf_query = """
        SELECT DISTINCT a.uniquename, s.uniquename
        FROM feature a
        JOIN feature_relationship ac ON ac.subject_id = a.feature_id
        JOIN feature t ON t.feature_id = ac.object_id
        JOIN feature_relationship sc ON sc.object_id = t.feature_id
        JOIN feature s ON s.feature_id = sc.subject_id
        WHERE a.is_obsolete is false
          AND a.is_analysis is false
          AND a.uniquename ~ '^FBal[0-9]{7}$'
          AND t.is_obsolete is false
          AND t.uniquename ~ '^FBtp[0-9]{7}$'
          AND s.is_obsolete is false
          AND s.is_analysis is false
          AND s.uniquename ~ '^FBsf[0-9]{10}$';"""
    ret_fbal_fbsf = connect(fbal_fbsf_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} FBal-FBtp-FBsf associations.'.format(now(), len(ret_fbal_fbsf)))
    FBAL_ID = 0
    FBSF_ID = 1
    fbal_fbsf_dict = {i[FBAL_ID]: i[FBSF_ID] for i in ret_fbal_fbsf}

    # Get RNAi FBsf-FBlc dict.
    rnai_fbsf_coll_query = """
        SELECT DISTINCT s.uniquename, l.name
        FROM feature s
        JOIN library_feature lf ON lf.feature_id = s.feature_id
        JOIN library l ON l.library_id = lf.library_id
        WHERE s.is_obsolete is false
          AND s.is_analysis is false
          AND s.uniquename ~ '^FBsf[0-9]{10}$'
          AND l.is_obsolete is false
          AND l.name ~ '(^TRiP-[0-9]{1}$|^VDRC-(GD|KK|SH)$)';"""
    ret_rnai_fbsf_coll  = connect(rnai_fbsf_coll_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} RNAi FBsf features.'.format(now(), len(ret_rnai_fbsf_coll)))
    FBSF_ID = 0
    COLL_NAME = 1
    rnai_fbsf_coll_dict = {i[FBSF_ID]: i[COLL_NAME] for i in ret_rnai_fbsf_coll}

    rnai_dict = {}
    FBAL_ID = 0
    FBAL_NAME = 1
    for row in ret_rnai_fbal:
        if row[FBAL_ID] in fbal_fbsf_dict.keys():
            related_fbsf = fbal_fbsf_dict[row[FBAL_ID]]
            try:
                rnai_dict[row[FBAL_ID]] = rnai_fbsf_coll_dict[related_fbsf]
            except KeyError:
                if '[NIG' in row[FBAL_NAME]:
                    rnai_dict[row[FBAL_ID]] = 'NIG'
                else:
                    rnai_dict[row[FBAL_ID]] = 'lab'
        else:
            if '[NIG' in row[FBAL_NAME]:
                rnai_dict[row[FBAL_ID]] = 'NIG'
            else:
                rnai_dict[row[FBAL_ID]] = 'lab'

    log.info('TIME: {}. Final RNAi dict has {} alleles.'.format(now(), len(rnai_dict.keys())))

    return rnai_dict


def get_crispr_info(db_connection):
    """Returns all crispr alleles and their related collections from the database.
       Required libraries: psycopg2, datetime.
       Required functions: connect(), now().
       Required global variables: none.
       Input: a database connection (db_connection).
       Output: an {FBal ID: collection} set of FBal IDs."""
    log.info('TIME: {}. Querying database for CRISPR alleles.'.format(now()))

    # Get FBal-FBsf dict.
    fbal_fbsf_query = """
        SELECT DISTINCT a.uniquename, s.uniquename
        FROM feature a
        JOIN feature_relationship ac ON ac.subject_id = a.feature_id
        JOIN feature t ON t.feature_id = ac.object_id
        JOIN feature_relationship sc ON sc.object_id = t.feature_id
        JOIN feature s ON s.feature_id = sc.subject_id
        WHERE a.is_obsolete is false
          AND a.is_analysis is false
          AND a.uniquename ~ '^FBal[0-9]{7}$'
          AND t.is_obsolete is false
          AND t.uniquename ~ '^FBtp[0-9]{7}$'
          AND s.is_obsolete is false
          AND s.is_analysis is false
          AND s.uniquename ~ '^FBsf[0-9]{10}$';"""
    ret_fbal_fbsf = connect(fbal_fbsf_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} FBal-FBtp-FBsf associations.'.format(now(), len(ret_fbal_fbsf)))

    # Get RNAi FBsf-FBlc dict.
    crispr_fbsf_coll_query = """
        SELECT DISTINCT s.uniquename, l.name
        FROM feature s
        JOIN library_feature lf ON lf.feature_id = s.feature_id
        JOIN library l ON l.library_id = lf.library_id
        WHERE s.is_obsolete is false
          AND s.is_analysis is false
          AND s.uniquename ~ '^FBsf[0-9]{10}$'
          AND l.is_obsolete is false
          AND l.name ~ '(^Weizmann-KO$|^TRiP-KO$|^TRiP-OE-)';"""
    ret_crispr_fbsf_coll  = connect(crispr_fbsf_coll_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} CRISPR FBsf features.'.format(now(), len(ret_crispr_fbsf_coll)))
    FBSF_ID = 0
    COLL_NAME = 1
    crispr_fbsf_coll_dict = {i[FBSF_ID]: i[COLL_NAME] for i in ret_crispr_fbsf_coll}
    log.info('TIME: {}. Have {} CRISPR FBsf-FBlc pairs in the dict.'.format(now(), len(crispr_fbsf_coll_dict.keys())))

    crispr_dict = {}
    FBAL_ID = 0
    FBSF_ID = 1
    for row in ret_fbal_fbsf:
        my_fbal = row[FBAL_ID]
        my_fbsf = row[FBSF_ID]
        log.debug('Found allele {} and seqfeat {}'.format(my_fbal, my_fbsf))
        try:
            crispr_dict[row[FBAL_ID]] = crispr_fbsf_coll_dict[row[FBSF_ID]]
        except (IndexError, KeyError):
            pass
    log.info('TIME: {}. Final CRISPR dict has {} alleles.'.format(now(), len(crispr_dict.keys())))

    return crispr_dict


def convert_pheno_string(input_string):
    """Removes FB/GO IDs, "@" signs and FB-sgml from some string.
       Required libraries: re, harvdev-utils
       Required functions: none.
       Required global variables: none.
       Input: a string.
       Output: a string with FB/GO IDs and "@" signs removed."""

    # Step 1 is to filter out IDs and odd FB chars/strings.
    substitution_dict = {
        'FB[a-z]{2}[0-9]{7,8}:': '',
        'GO[0-9]{8}:': '',
        '@': '',
        '<up>': '[',
        '</up>': ']',
        '<down>': '[[',
        '</down>': ']]',
    }
    for k, v in substitution_dict.items():
        input_string = re.sub(k, v, input_string)

    # Step 2 uses harvdev_utils/char_conversions/sgml_to_plain_text.py to filter Greek sgmls.
    output_string = sgml_to_plain_text(input_string)

    return output_string


def get_drivers(input_string):
    """Finds all GAL4 drivers mentioned in a phenotype string.
       Required libraries: re, psycopg2.
       Required functions: connect().
       Required global variables: none.
       Input: a string.
       Output: a list of GAL4s (just the allele designator)."""

    driver_list = []
    # driver_regex = r'(?<=GAL4\<up\>).+(?=</up\>)'             # Option 1 - if many alleles, returns concetanates.
    # driver_id_regex = r'FBal[0-9]{7}(?=\:Scer\\GAL4<up\>)'    # Option 2 - unambiguous, but many db queries slow.
    # Option 3 - break up string where driver string starts.
    # This will not work if ever an allele designator has superscripts within. None right now though.
    string_parts = input_string.split('GAL4<up>')
    # Skip over the first element of the split list (won't have GAL4 info, ever).
    for i in range(1, len(string_parts)):
        # The driver's allele designator is before the first "</up>" string, always.
        driver = string_parts[i].split('</up>')[0]
        driver_list.append(driver)

    return driver_list


def make_allele_to_gene_dict(db_connection):
    """Returns a dictionary of alleles (keys) and genes (values).
       Required libraries: psycopg2, datetime.
       Required functions: connect(), now().
       Required global variables: none.
       Input: a database connection (db_connection).
       Output: an {allele_id: (gene_id, gene_name)} dictionary."""
    log.info('TIME: {}. Querying database for allele-gene associations.'.format(now()))

    fb_allele_gene_query = """
        SELECT DISTINCT a.uniquename, g.uniquename, g.name
        FROM feature a
        JOIN feature_relationship fr ON fr.subject_id = a.feature_id
        JOIN feature g ON g.feature_id = fr.object_id
        JOIN cvterm cvt ON cvt.cvterm_id = fr.type_id
        WHERE cvt.name = 'alleleof'
          and a.is_obsolete = false and a.is_analysis = false and a.uniquename ~ '^FBal[0-9]{7}$'
          and g.is_obsolete = false and g.is_analysis = false and g.uniquename ~ '^FBgn[0-9]{7}$';"""
    ret_allele_gene_info = connect(fb_allele_gene_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} allele-gene pairs.'.format(now(), len(ret_allele_gene_info)))
    allele_to_gene_dict = {i[0]: (i[1], i[2]) for i in ret_allele_gene_info}

    return allele_to_gene_dict


def make_fbrf_to_pmid_dict(db_connection):
    """Returns a dictionary of pubs with FBrf IDs (keys) and PMIDs (values).
       Required libraries: psycopg2, datetime.
       Required functions: connect(), now().
       Required global variables: none.
       Input: a database connection (db_connection).
       Output: an {FBrf ID: PMID} dictionary."""
    log.info('TIME: {}. Querying database for FBrf-PMID associations.'.format(now()))

    fb_pmid_query = """
        SELECT DISTINCT p.uniquename, dbx.accession
        FROM pub p
        JOIN pub_dbxref pdbx ON pdbx.pub_id = p.pub_id
        JOIN dbxref dbx ON dbx.dbxref_id = pdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE p.is_obsolete = false
          and pdbx.is_current = true
          and db.name = 'pubmed';"""

    ret_pmid_info = connect(fb_pmid_query, 'no_query', db_connection)
    log.info('TIME: {}. Found {} current pubs.'.format(now(), len(ret_pmid_info)))
    fbrf_to_pmid_dict = {i[0]: i[1] for i in ret_pmid_info}

    return fbrf_to_pmid_dict


def process_pheno_info(input_data, db_connection):
    """Takes SQL results and returns a list of dictionaries for tsv output.
       Required libraries: datetime, re, harvdev-utils.
       Required functions: now(), replace_ids(), get_drivers(), make_allele_to_gene_dict(), make_fbrf_to_pmid_dict().
       Required global variables: none.
       Input: a list of tuples representing phenotype info.
       Output: a list of dictionaries representing phenotype info."""
    log.info('TIME: {}. Starting to process phenotype info retrieved from database.'.format(now()))

    log.info('TIME: {}. Getting additional information from the database.'.format(now()))
    fbal2fbgn = make_allele_to_gene_dict(db_connection)
    fbrf2pmid = make_fbrf_to_pmid_dict(db_connection)
    rnai_dict = get_rnai_info(db_connection)
    crispr_dict = get_crispr_info(db_connection)

    data_list = []
    for i in input_data:
        allele_id = i[0]
        allele_symbol = i[1]
        gene_id = fbal2fbgn[i[0]][0]
        gene_symbol = fbal2fbgn[i[0]][1]
        if allele_id in rnai_dict.keys():
            reagent_type = 'RNAi'
            reagent_source = rnai_dict[allele_id]
        elif allele_id in crispr_dict.keys():
            reagent_type = 'CRISPR'
            reagent_source = crispr_dict[allele_id]
        else:
            continue
        pheno = {}

        pheno['FBgn'] = gene_id
        pheno['gene_symbol'] = gene_symbol
        pheno['FBal'] = allele_id
        pheno['allele_symbol'] = allele_symbol
        pheno['reagent_type'] = reagent_type
        pheno['reagent_source'] = reagent_source
        # Chose to split off gene symbol since genes like 'su(w[a])' make regex difficult.
        pheno['reagent_id'] = allele_symbol.split(gene_symbol)[1].lstrip('[').rstrip(']')
        pheno['driver'] = ', '.join(get_drivers(i[2]))
        pheno['phenotype'] = convert_pheno_string(i[2])
        pheno['FBrf'] = i[3]
        if i[3] in fbrf2pmid.keys():
            pheno['PMID'] = fbrf2pmid[i[3]]
        else:
            pheno['PMID'] = None
        data_list.append(pheno)

    log.info('TIME: {}. Done processing phenotype. Sending {} annotations to bulk file.'.format(now(), len(data_list)))

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
       Input: a data dictionary with 'metadata' key with a list of data dictionaries under 'data' key.
              Optional "headers" kwarg takes a list of headers (should match dict keys).
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
