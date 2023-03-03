#!/usr/bin/env python3

# report_more_current_gene_product_ids.py
# author: Gil dos Santos
# usage: report_more_current_gene_product_ids.py [-h] [-v VERBOSE]
################################################################################

import argparse
import csv
import datetime
import logging
import os
import psycopg2
import sys

report_name = 'fbgn_fbtr_fbpp_expanded'

parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-v', '--verbose', action='store_true', help='DEBUG-level logging.', required=False)
args = parser.parse_args()

# For running script within GoCD pipeline.
database_host = os.environ['SERVER']
database = os.environ['DATABASE']
username = os.environ['USER']
password = os.environ['PGPASSWORD']
annotation_release = os.environ['ANNOTATIONRELEASE']
database_release = os.environ['RELEASE']
output_filename = './output/' + report_name + '_' + database + '.tsv'
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

    # Database connection.
    conn_string = "host=%s dbname=%s user=%s password='%s'" % (database_host, database, username, password)
    conn = psycopg2.connect(conn_string)
    log.info('TIME: {}. Successful connection to database {} on host {}.'.format(now(), database, database_host))

    # Create the TSV data structure and fill in the file metadata.
    to_export_as_tsv = {}
    to_export_as_tsv['metaData'] = {}
    to_export_as_tsv['metaData']['title'] = 'FlyBase expanded current gene product report'
    to_export_as_tsv['metaData']['dateProduced'] = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    to_export_as_tsv['metaData']['database'] = database
    to_export_as_tsv['data'] = []

    # Queries to retrieve different types of gene products and proper synonyms. Test: mir-11, ban, roX1, alphaTub84B.
    fb_coding_gene_query = """
    SELECT DISTINCT o.abbreviation, g.uniquename, g.feature_id, g.name, t.feature_id,
                    cvtrna.name, t.uniquename, t.name, p.feature_id, p.uniquename, p.name
    FROM feature g
    JOIN featureloc flg ON (flg.feature_id = g.feature_id)
    JOIN feature_relationship tg ON (tg.object_id = g.feature_id)
    JOIN feature t ON (t.feature_id = tg.subject_id)
    JOIN featureloc flt ON (flt.feature_id = t.feature_id)
    JOIN cvterm cvtrna ON (cvtrna.cvterm_id = t.type_id)
    JOIN feature_relationship tp ON (tp.object_id = t.feature_id)
    JOIN feature p ON (p.feature_id = tp.subject_id)
    JOIN featureloc flp ON (flp.feature_id = p.feature_id)
    JOIN organism o ON (o.organism_id = g.organism_id)
    WHERE g.is_obsolete = false and g.type_id = 219 and g.uniquename ~ '^FBgn[0-9]{7}$'
      and p.is_obsolete = false and p.uniquename ~ '^FBpp[0-9]{7}$' and NOT p.name LIKE '%-XP'
      and t.is_obsolete = false and t.uniquename ~ '^FBtr[0-9]{7}$' and NOT t.name LIKE '%-XR'
      and cvtrna.name = 'mRNA' and tp.type_id = 27;"""

    fb_ncgene_query = """
    SELECT DISTINCT o.abbreviation, g.uniquename, g.feature_id, g.name,
                    t.feature_id, cvtrna.name, t.uniquename, t.name
    FROM feature g
    JOIN featureloc flg ON (flg.feature_id = g.feature_id)
    JOIN feature_relationship tg ON (tg.object_id = g.feature_id)
    JOIN feature t ON (t.feature_id = tg.subject_id)
    JOIN featureloc flt ON (flt.feature_id = t.feature_id)
    JOIN cvterm cvtrna ON (cvtrna.cvterm_id = t.type_id)
    JOIN organism o ON (o.organism_id = g.organism_id)
    WHERE g.is_obsolete = false and g.type_id = 219 and g.uniquename ~ '^FBgn[0-9]{7}$'
      and t.is_obsolete = false and t.uniquename ~ '^FBtr[0-9]{7}$' and NOT t.name LIKE '%-XR'
      and cvtrna.name in ('pre_miRNA','tRNA','ncRNA','rRNA','snRNA','snoRNA','pseudogene');"""

    fb_mirna_query = """
    SELECT DISTINCT o.abbreviation, g.uniquename, g.feature_id, g.name,
                    mi.feature_id, cvtmirna.name, mi.uniquename, mi.name
    FROM feature g
    JOIN featureloc flg ON (flg.feature_id = g.feature_id)
    JOIN feature_relationship tg ON (tg.object_id = g.feature_id)
    JOIN feature t ON (t.feature_id = tg.subject_id)
    JOIN featureloc flt ON (flt.feature_id = t.feature_id)
    JOIN cvterm cvtrna ON (cvtrna.cvterm_id = t.type_id)
    JOIN feature_relationship tmi ON (tmi.object_id = t.feature_id)
    JOIN feature mi ON (mi.feature_id = tmi.subject_id)
    JOIN cvterm cvtmirna ON (cvtmirna.cvterm_id = mi.type_id)
    JOIN featureloc flmi ON (flmi.feature_id = mi.feature_id)
    JOIN organism o ON (o.organism_id = g.organism_id)
    WHERE g.is_obsolete = false and g.type_id = 219 and g.uniquename ~ '^FBgn[0-9]{7}$'
      and t.is_obsolete = false and t.uniquename ~ '^FBtr[0-9]{7}$' and NOT t.name LIKE '%-XR'
      and cvtrna.name in ('pre_miRNA') and tmi.type_id = 27 and cvtmirna.name = 'miRNA'
      and mi.is_obsolete = false and mi.uniquename ~ '^FBtr[0-9]{7}$' and NOT mi.name LIKE '%-XR';"""

    fb_symbol_query = """
    SELECT DISTINCT f.feature_id, s.name
    FROM feature f
    JOIN feature_synonym fs ON (fs.feature_id = f.feature_id)
    JOIN synonym s ON (s.synonym_id = fs.synonym_id)
    JOIN cvterm cvts ON (cvts.cvterm_id = s.type_id)
    JOIN featureloc fl ON (fl.feature_id = f.feature_id)
    WHERE fs.is_current = true and f.is_obsolete = false and f.uniquename ~ '^FB(gn|tr|pp)[0-9]{7}$'
      and cvts.name = 'symbol';"""

    fb_fullname_query = """
    SELECT DISTINCT f.feature_id, s.name
    FROM feature f
    JOIN feature_synonym fs ON (fs.feature_id = f.feature_id)
    JOIN synonym s ON (s.synonym_id = fs.synonym_id)
    JOIN cvterm cvts ON (cvts.cvterm_id = s.type_id)
    JOIN featureloc fl ON (fl.feature_id = f.feature_id)
    WHERE fs.is_current = true and f.is_obsolete = false and f.uniquename ~ '^FBgn[0-9]{7}$'
      and cvts.name = 'fullname';"""

    fb_annotation_query = """
    SELECT DISTINCT f.feature_id, dbx.accession
    FROM feature f
    JOIN feature_dbxref fdbx ON (fdbx.feature_id = f.feature_id)
    JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
    JOIN db ON (db.db_id = dbx.db_id)
    WHERE fdbx.is_current = true and f.is_obsolete = false and f.uniquename ~ '^FBgn[0-9]{7}$'
      and dbx.accession ~ '(C|G)[A-Z][0-9]{4,5}$'
      and db.name = 'FlyBase Annotation IDs';"""

    fb_gene_type_query = """
    SELECT DISTINCT f.feature_id, fp.value
    FROM feature f
    JOIN featureprop fp ON (fp.feature_id = f.feature_id)
    JOIN cvterm cvt ON (cvt.cvterm_id = fp.type_id)
    WHERE f.is_obsolete = false and f.uniquename ~ '^FBgn[0-9]{7}$'
      and cvt.name = 'promoted_gene_type';"""

    # Retrieve gene products from the database.
    log.info('TIME: {}. Querying database for info.'.format(now()))

    coding_product_info = connect(fb_coding_gene_query, 'no_query', conn)
    log.info('The "fb_coding_gene_query" returned {} results.'.format(len(coding_product_info)))
    coding_product_info.sort()

    ncgene_product_info = connect(fb_ncgene_query, 'no_query', conn)
    log.info('The "fb_ncgene_query" returned {} results.'.format(len(ncgene_product_info)))
    ncgene_product_info.sort()

    mirna_product_info = connect(fb_mirna_query, 'no_query', conn)
    log.info('The "fb_mirna_query" returned an additional {} results.'.format(len(mirna_product_info)))
    mirna_product_info.sort()

    gene_product_info = []
    gene_product_info.extend(coding_product_info)
    gene_product_info.extend(ncgene_product_info)
    gene_product_info.extend(mirna_product_info)
    gene_product_info.sort()

    log.info('Overall, there are {} results.'.format(len(gene_product_info)))

    # Build the list of data elements for export.
    for line in gene_product_info:
        this_dict = {}
        this_dict['organism'] = line[0]
        this_dict['gene_ID'] = line[1]
        this_dict['gene_feature_id'] = line[2]
        this_dict['gene_name'] = line[3]
        this_dict['transcript_feature_id'] = line[4]
        this_dict['transcript_type'] = line[5]
        this_dict['transcript_ID'] = line[6]
        this_dict['transcript_name'] = line[7]
        if len(line) == 11:
            this_dict['polypeptide_feature_id'] = line[8]
            this_dict['polypeptide_ID'] = line[9]
            this_dict['polypeptide_name'] = line[10]
        else:
            this_dict['polypeptide_feature_id'] = None
            this_dict['polypeptide_ID'] = None
            this_dict['polypeptide_name'] = None
        to_export_as_tsv['data'].append(this_dict)

    # Retrieve fullname and symbol synonyms from the database.
    log.info('TIME: {}. Querying database for symbol and fullname synonyms.'.format(now()))

    symbol_info = connect(fb_symbol_query, 'no_query', conn)
    symbol_dict = {i[0]: i[1] for i in symbol_info}

    fullname_info = connect(fb_fullname_query, 'no_query', conn)
    fullname_dict = {i[0]: i[1] for i in fullname_info}

    annotation_info = connect(fb_annotation_query, 'no_query', conn)
    annotation_dict = {i[0]: i[1] for i in annotation_info}

    gene_type_info = connect(fb_gene_type_query, 'no_query', conn)
    gene_type_dict = {i[0]: i[1][11:-1] for i in gene_type_info}

    conn.close()

    # Fold in synonym info into the data list.
    log.info('TIME: {}. Adding synonyms to gene product list.'.format(now()))

    for item in to_export_as_tsv['data']:
        log.debug('Adding synonyms to this item:\n')
        log.debug('\t{}'.format(item))
        try:
            item['gene_type'] = gene_type_dict[item['gene_feature_id']]
        except KeyError:
            item['gene_type'] = 'gene'
        item['gene_symbol'] = symbol_dict[item['gene_feature_id']]
        if item['gene_feature_id'] in annotation_dict.keys():
            item['annotation_ID'] = annotation_dict[item['gene_feature_id']]
        else:
            item['annotation_ID'] = 'NOANNOID'
        if item['gene_feature_id'] in fullname_dict.keys():
            item['gene_fullname'] = fullname_dict[item['gene_feature_id']]
        else:
            item['gene_fullname'] = None
        item['transcript_symbol'] = symbol_dict[item['transcript_feature_id']]
        if item['polypeptide_feature_id']:
            item['polypeptide_symbol'] = symbol_dict[item['polypeptide_feature_id']]
        else:
            item['polypeptide_symbol'] = None

    if len(to_export_as_tsv['data']) == 0:
        log.error('No data to report.')
        raise

    # Generate the output file.
    log.info('TIME: {}. Writing data to output file.'.format(now()))

    output_file = open(output_filename, 'w')

    # Write the headers.
    output_file.write('## {}\n'.format(to_export_as_tsv['metaData']['title']))
    output_file.write('## Generated: {}\n'.format(to_export_as_tsv['metaData']['dateProduced']))
    output_file.write('## Using datasource: {}\n'.format(to_export_as_tsv['metaData']['database']))
    output_file.write('##\n## ')

    headers = [
        'organism',
        'gene_type',
        'gene_ID',
        'gene_symbol',
        'gene_fullname',
        'annotation_ID',
        'transcript_type',
        'transcript_ID',
        'transcript_symbol',
        'polypeptide_ID',
        'polypeptide_symbol'
        ]

    csv_writer = csv.DictWriter(output_file, fieldnames=headers, delimiter='\t', extrasaction='ignore')
    csv_writer.writeheader()

    # Write out the data.
    for data_item in to_export_as_tsv['data']:
        csv_writer.writerow(data_item)

    output_file.write('## Finished {}.'.format(to_export_as_tsv['metaData']['title']))
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
