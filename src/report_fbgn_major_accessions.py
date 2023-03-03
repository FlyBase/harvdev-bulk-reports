#!/usr/bin/env python3

# report_fbgn_major_accessions
# author: Gil dos Santos (modeled after original perl script by Dave Emmert)
# usage: report_fbgn_major_accessions.py [-h] [-v VERBOSE]
################################################################################

import argparse
# import configparser
import csv
import datetime
import logging
import os
import psycopg2
import re
import sys

report_name = 'fbgn_NAseq_Uniprot'

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
    to_export_as_tsv['metaData']['title'] = 'FlyBase FBgn-Major Accessions Table'
    to_export_as_tsv['metaData']['dateProduced'] = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    to_export_as_tsv['metaData']['database'] = database
    to_export_as_tsv['data'] = []

    # Query for supporting accession tables (featureprop)
    # Will need to split the multi-line prop, and look for one vs two accessions.
    fb_supporting_acc_query = """
    SELECT f.name, o.abbreviation, f.uniquename, fp.value
    FROM feature f
    JOIN organism o ON (o.organism_id = f.organism_id)
    JOIN cvterm cvtf ON (cvtf.cvterm_id = f.type_id)
    JOIN featureprop fp ON (fp.feature_id = f.feature_id)
    JOIN cvterm cvtfp ON (cvtfp.cvterm_id = fp.type_id)
    WHERE cvtf.name = 'gene' and f.is_obsolete = false and f.is_analysis = false and f.uniquename ~ '^FBgn[0-9]{7}$'
      and cvtfp.name = 'derived_supporting_accessions';
    """

    # Query for EntrezGene accessions (feature_dbxref).
    fb_uniprot_query = """
    SELECT f.name, o.abbreviation, f.uniquename, dbx.accession
    FROM feature f
    JOIN cvterm cvtf ON (cvtf.cvterm_id = f.type_id)
    JOIN organism o ON (o.organism_id = f.organism_id)
    JOIN feature_dbxref fdbx ON (f.feature_id = fdbx.feature_id)
    JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
    JOIN db ON (db.db_id = dbx.db_id)
    WHERE cvtf.name = 'gene' and f.is_obsolete = false and f.is_analysis = false and f.uniquename ~ '^FBgn[0-9]{7}$'
      and fdbx.is_current = true and db.name in ('UniProt/Swiss-Prot','UniProt/TrEMBL');
    """

    # Query for UniprotKB/Swiss-Prot/TreMBL accessions (feature_dbxref).
    fb_entrez_query = """
    SELECT f.name, o.abbreviation, f.uniquename, dbx.accession
    FROM feature f
    JOIN cvterm cvtf ON (cvtf.cvterm_id = f.type_id)
    JOIN organism o ON (o.organism_id = f.organism_id)
    JOIN feature_dbxref fdbx ON (f.feature_id = fdbx.feature_id)
    JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
    JOIN db ON (db.db_id = dbx.db_id)
    WHERE cvtf.name = 'gene' and f.is_obsolete = false and f.is_analysis = false and f.uniquename ~ '^FBgn[0-9]{7}$'
      and fdbx.is_current = true and db.name = 'EntrezGene';
    """

    # Query for all non-coding gene REFSEQ accessions (feature_dbxref).
    fb_refseq_ncrna_query = """
    SELECT g.name, o.abbreviation, g.uniquename, dbx.accession
    FROM feature g
    JOIN cvterm cvtg ON (cvtg.cvterm_id = g.type_id)
    JOIN organism o ON (o.organism_id = g.organism_id)
    JOIN feature_relationship fr ON (fr.object_id = g.feature_id)
    JOIN cvterm cvtfr ON (cvtfr.cvterm_id = fr.type_id)
    JOIN feature t ON (t.feature_id = fr.subject_id)
    JOIN cvterm cvtt ON (cvtt.cvterm_id = t.type_id)
    JOIN feature_dbxref fdbx ON (t.feature_id = fdbx.feature_id)
    JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
    JOIN db ON (db.db_id = dbx.db_id)
    WHERE cvtg.name = 'gene' and g.is_obsolete = false and g.is_analysis = false and g.uniquename ~ '^FBgn[0-9]{7}$'
      and t.is_obsolete = false and t.is_analysis = false and t.uniquename ~ '^FBtr[0-9]{7}$'
      and not cvtt.name in ('mRNA','pre_miRNA') and cvtfr.name = 'partof'
      and fdbx.is_current = true and db.name = 'REFSEQ';
    """

    # Query for transcript/protein REFSEQ accession pairs (feature_dbxref).
    fb_refseq_coding_query = """
    SELECT g.name, o.abbreviation, g.uniquename, dbxt.accession, dbxp.accession
    FROM feature g
    JOIN cvterm cvtg ON (cvtg.cvterm_id = g.type_id)
    JOIN organism o ON (o.organism_id = g.organism_id)
    JOIN feature_relationship tg ON (tg.object_id = g.feature_id)
    JOIN cvterm cvttg ON (cvttg.cvterm_id = tg.type_id)
    JOIN feature t ON (t.feature_id = tg.subject_id)
    JOIN cvterm cvtt ON (cvtt.cvterm_id = t.type_id)
    JOIN feature_relationship pt ON (pt.object_id = t.feature_id)
    JOIN feature p ON (p.feature_id = pt.subject_id)
    JOIN cvterm cvtpt ON (cvtpt.cvterm_id = pt.type_id)
    JOIN cvterm cvtp ON (cvtp.cvterm_id = p.type_id)
    JOIN feature_dbxref fdbxt ON (t.feature_id = fdbxt.feature_id)
    JOIN dbxref dbxt ON (dbxt.dbxref_id = fdbxt.dbxref_id)
    JOIN feature_dbxref fdbxp ON (p.feature_id = fdbxp.feature_id)
    JOIN dbxref dbxp ON (dbxp.dbxref_id = fdbxp.dbxref_id)
    JOIN db ON (db.db_id = dbxt.db_id)
    WHERE db.db_id = dbxp.db_id
      and cvtg.name = 'gene' and g.is_obsolete = false and g.is_analysis = false and g.uniquename ~ '^FBgn[0-9]{7}$'
      and cvtt.name = 'mRNA' and t.is_obsolete = false and t.is_analysis = false and t.uniquename ~ '^FBtr[0-9]{7}$'
      and cvtp.name = 'polypeptide' and p.is_obsolete = false and p.is_analysis = false and p.uniquename ~ '^FBpp[0-9]{7}$'
      and cvttg.name = 'partof' and cvtpt.name = 'producedby'
      and fdbxt.is_current = true and fdbxp.is_current = true and db.name = 'REFSEQ';
    """

    # Recursive query for pre_miRNA and miRNA REFSEQ accessions (feature_dbxref).
    fb_refseq_mirna_query = """
    WITH RECURSIVE rel_object AS
    (SELECT g.name as gname, o.abbreviation as spp, g.uniquename as guname, fr.subject_id as key_id, dbx.accession
    FROM feature g
    JOIN cvterm cvtg ON (cvtg.cvterm_id = g.type_id)
    JOIN feature_relationship fr ON (fr.object_id = g.feature_id)
    JOIN feature t ON (t.feature_id = fr.subject_id)
    JOIN cvterm cvtfr ON (cvtfr.cvterm_id = fr.type_id)
    JOIN cvterm cvtt ON (cvtt.cvterm_id = t.type_id)
    JOIN organism o ON (o.organism_id = g.organism_id)
    JOIN feature_dbxref fdbx ON (fdbx.feature_id = t.feature_id)
    JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
    JOIN db ON (db.db_id = dbx.db_id)
    WHERE g.is_obsolete = false and g.is_analysis = false and g.uniquename ~ '^FBgn[0-9]{7}$'
    and t.is_obsolete = false and t.is_analysis = false and t.uniquename ~ '^FBtr[0-9]{7}$'
    and cvtfr.name = 'partof' and cvtg.name = 'gene' and cvtt.name = 'pre_miRNA'
    and t.organism_id = o.organism_id
    and fdbx.is_current = true
    and db.name = 'REFSEQ'
    UNION
    SELECT ro.gname, ro.spp, ro.guname, fr.subject_id, dbx.accession
    FROM feature_relationship fr
    JOIN rel_object ro ON (ro.key_id = fr.object_id)
    JOIN feature t ON (t.feature_id = fr.subject_id)
    JOIN cvterm cvtfr ON (cvtfr.cvterm_id = fr.type_id)
    JOIN cvterm cvtt ON (cvtt.cvterm_id = t.type_id)
    JOIN feature_dbxref fdbx ON (fdbx.feature_id = t.feature_id)
    JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
    JOIN db ON (db.db_id = dbx.db_id)
    WHERE t.is_obsolete = false and t.is_analysis = false and t.uniquename ~ '^FBtr[0-9]{7}$'
    and cvtfr.name = 'producedby' and cvtt.name = 'miRNA'
    and fdbx.is_current = true and db.name = 'REFSEQ'
    )
    SELECT * FROM rel_object ORDER BY guname;
    """

    # Retrieve data from the database.
    log.info('TIME: {}. Querying database for supporting accession props.'.format(now()))
    fb_supporting_acc_info = connect(fb_supporting_acc_query, 'no_query', conn)
    for row in fb_supporting_acc_info:
        prop_list = row[3].split('\n')
        for prop in prop_list:
            if re.search(r'\t', prop):
                accession = {
                    'gene_symbol': row[0],
                    'organism_abbreviation': row[1],
                    'primary_FBgn#': row[2],
                    'nucleotide_accession': prop.split('\t')[0].strip(),
                    'na_based_protein_accession': prop.split('\t')[1].strip()
                }
                to_export_as_tsv['data'].append(accession)
            else:
                accession = {
                    'gene_symbol': row[0],
                    'organism_abbreviation': row[1],
                    'primary_FBgn#': row[2],
                    'nucleotide_accession': prop.strip()
                }
                to_export_as_tsv['data'].append(accession)

    log.info('TIME: {}. Querying database for UniProt accessions.'.format(now()))
    fb_uniprot_info = connect(fb_uniprot_query, 'no_query', conn)
    uniprot_acc = [{
        'gene_symbol': i[0],
        'organism_abbreviation': i[1],
        'primary_FBgn#': i[2],
        'UniprotKB/Swiss-Prot/TrEMBL_accession': i[3]
        } for i in fb_uniprot_info]
    to_export_as_tsv['data'].extend(uniprot_acc)

    log.info('TIME: {}. Querying database for entrez accessions.'.format(now()))
    fb_entrez_info = connect(fb_entrez_query, 'no_query', conn)
    entrez_acc = [{
        'gene_symbol': i[0],
        'organism_abbreviation': i[1],
        'primary_FBgn#': i[2],
        'EntrezGene_ID': i[3]
        } for i in fb_entrez_info]
    to_export_as_tsv['data'].extend(entrez_acc)

    log.info('TIME: {}. Querying database for non-coding REFSEQ accessions.'.format(now()))
    fb_refseq_ncrna_info = connect(fb_refseq_ncrna_query, 'no_query', conn)
    refseq_nc_acc = [{
        'gene_symbol': i[0],
        'organism_abbreviation': i[1],
        'primary_FBgn#': i[2],
        'RefSeq_transcripts': i[3]
        } for i in fb_refseq_ncrna_info]
    to_export_as_tsv['data'].extend(refseq_nc_acc)

    log.info('TIME: {}. Querying database for coding REFSEQ accessions.'.format(now()))
    fb_refseq_coding_info = connect(fb_refseq_coding_query, 'no_query', conn)
    refseq_coding_acc = [{
        'gene_symbol': i[0],
        'organism_abbreviation': i[1],
        'primary_FBgn#': i[2],
        'RefSeq_transcripts': i[3],
        'RefSeq_proteins': i[4]
        } for i in fb_refseq_coding_info]
    to_export_as_tsv['data'].extend(refseq_coding_acc)

    log.info('TIME: {}. Querying database for miRNA REFSEQ accessions.'.format(now()))
    fb_refseq_mirna_info = connect(fb_refseq_mirna_query, 'no_query', conn)
    refseq_mirna_acc = [{
        'gene_symbol': i[0],
        'organism_abbreviation': i[1],
        'primary_FBgn#': i[2],
        'RefSeq_transcripts': i[4]
        } for i in fb_refseq_mirna_info]
    to_export_as_tsv['data'].extend(refseq_mirna_acc)

    conn.close()

    # Sort the data elements.
    log.info('TIME: {}. Sorting accessions.'.format(now()))
    to_export_as_tsv['data'] = sorted(to_export_as_tsv['data'], key=lambda i: (i['primary_FBgn#']))

    # Instead of printing multiple tabs to get various xrfs into the right column ...
    # Use csv_writer to print to specific field.
    # csv_writer will skip over a field if a dict lacks that key (nice!).

    # Generate the output file.
    log.info('TIME: {}. Writing data to output file.'.format(now()))

    output_file = open(output_filename, 'w')

    # Write the headers.
    output_file.write('## {}\n'.format(to_export_as_tsv['metaData']['title']))
    output_file.write('## Generated: {}\n'.format(to_export_as_tsv['metaData']['dateProduced']))
    output_file.write('## Using datasource: {}\n'.format(to_export_as_tsv['metaData']['database']))
    output_file.write('##\n## ')

    headers = [
        'gene_symbol',
        'organism_abbreviation',
        'primary_FBgn#',
        'nucleotide_accession',
        'na_based_protein_accession',
        'UniprotKB/Swiss-Prot/TrEMBL_accession',
        'EntrezGene_ID',
        'RefSeq_transcripts',
        'RefSeq_proteins'
    ]

    csv_writer = csv.DictWriter(output_file, fieldnames=headers, delimiter='\t', extrasaction='ignore')
    csv_writer.writeheader()

    # Write the data.
    for data_item in to_export_as_tsv['data']:
        csv_writer.writerow(data_item)

    output_file.write('## Finished {} report.'.format(to_export_as_tsv['metaData']['title']))
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
