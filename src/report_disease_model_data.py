#!/usr/bin/env python3

# report_disease_model_data.py
# author: Gil dos Santos
# usage: report_disease_model_data.py [-h] [-v VERBOSE]
################################################################################

import argparse
import csv
import datetime
import logging
import os
import psycopg2
import re
import sys

report_name = 'disease_model_annotations'

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
    to_export_as_tsv['metaData']['title'] = 'FlyBase disease model annotation report'
    to_export_as_tsv['metaData']['dateProduced'] = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    to_export_as_tsv['metaData']['database'] = database
    to_export_as_tsv['data'] = []

    # Query for phenotype-based annotations.
    fb_pheno_query = """
        SELECT DISTINCT g.uniquename, g.name, org.abbreviation, qual.value, db.name||':'||dbx.accession,
                        cvt.name, a.uniquename, a.name, evid.value, p.uniquename
        FROM feature g
        JOIN organism org ON (org.organism_id = g.organism_id)
        JOIN feature_relationship fr ON (fr.object_id = g.feature_id)
        JOIN feature a ON (a.feature_id = fr.subject_id)
        JOIN cvterm cvtfr ON (cvtfr.cvterm_id = fr.type_id)
        JOIN feature_cvterm fcvt ON (fcvt.feature_id = a.feature_id)
        JOIN cvterm cvt ON (cvt.cvterm_id = fcvt.cvterm_id)
        JOIN cv ON (cv.cv_id = cvt.cv_id)
        JOIN dbxref dbx ON (dbx.dbxref_id = cvt.dbxref_id)
        JOIN db ON (db.db_id = dbx.db_id)
        JOIN pub p ON (p.pub_id = fcvt.pub_id)
        JOIN feature_cvtermprop qual ON (qual.feature_cvterm_id = fcvt.feature_cvterm_id and qual.type_id = 60522)
        JOIN feature_cvtermprop evid ON (evid.feature_cvterm_id = fcvt.feature_cvterm_id and evid.type_id = 60483)
        WHERE g.is_obsolete = false and g.uniquename ~ '^FBgn[0-9]{7}$'
        and cvtfr.name = 'alleleof'
        and a.is_obsolete = false and a.uniquename ~ '^FBal[0-9]{7}$'
        and cv.name = 'disease_ontology' and cvt.is_obsolete = 0
        and db.name = 'DOID' and p.is_obsolete = false and p.uniquename != 'FBrf0241599'
        and qual.rank = evid.rank
        ORDER BY g.uniquename, a.uniquename, cvt.name """

    # Query for orthology-based annotations.
    # I get the exact same fly-human-DOterm results as Chris' script (see Jira DB-549).
    fb_ortho_query = """
        SELECT DISTINCT fly.uniquename, fly.name, qual.value, db.name||':'||dbx.accession,
                        cvt.name, human.uniquename, human.name, evid.value, p.uniquename
        FROM feature fly
        JOIN organism org1 ON (org1.organism_id = fly.organism_id)
        JOIN feature_relationship fr ON (fr.subject_id = fly.feature_id)
        JOIN feature human ON (human.feature_id = fr.object_id)
        JOIN organism org2 ON (org2.organism_id = human.organism_id)
        JOIN cvterm cvtfr ON (cvtfr.cvterm_id = fr.type_id)
        JOIN feature_cvterm fcvt ON (fcvt.feature_id = human.feature_id)
        JOIN cvterm cvt ON (cvt.cvterm_id = fcvt.cvterm_id)
        JOIN cv ON (cv.cv_id = cvt.cv_id)
        JOIN dbxref dbx ON (dbx.dbxref_id = cvt.dbxref_id)
        JOIN db ON (db.db_id = dbx.db_id)
        JOIN pub p ON (p.pub_id = fcvt.pub_id)
        JOIN feature_cvtermprop qual ON (qual.feature_cvterm_id = fcvt.feature_cvterm_id and qual.type_id = 60522)
        JOIN feature_cvtermprop evid ON (evid.feature_cvterm_id = fcvt.feature_cvterm_id and evid.type_id = 60483)
        WHERE fly.is_obsolete = false and fly.uniquename ~ '^FBgn[0-9]{7}$' and org1.abbreviation = 'Dmel'
        and cvtfr.name = 'orthologous_to'
        and ((fr.value LIKE '%,%,%' and fr.value LIKE '%Yes%No%') or
            (fr.value LIKE '%,%,%' and fr.value LIKE '%No%Yes%') or
            (fr.value LIKE '%,%' and fr.value LIKE '%Yes%Yes%'))
        and human.is_obsolete = false and human.uniquename ~ '^FBog[0-9]{10}$' and org2.abbreviation = 'Hsap'
        and cv.name = 'disease_ontology' and cvt.is_obsolete = 0
        and db.name = 'DOID' and p.is_obsolete = false and p.uniquename = 'FBrf0241599'
        and qual.rank = evid.rank
        ORDER BY fly.uniquename, human.uniquename, cvt.name """

    # Query to get symbol for a feature (gene or allele).
    # Will need to trim leading "Hsap\" for some columns but not others.
    fb_synonym_sgml_query = """
        SELECT s.name
        FROM feature f
        JOIN feature_synonym fs ON (fs.feature_id = f.feature_id)
        JOIN synonym s ON (s.synonym_id = fs.synonym_id)
        JOIN cvterm cvts ON (cvts.cvterm_id = s.type_id)
        WHERE fs.is_current = true and fs.is_internal = false and cvts.name = 'symbol'
        and f.is_obsolete = false and f.uniquename ~ '^FB(og|gn|al)([0-9]{7}|[0-9]{10})$'
        and f.uniquename = %s """

    # Query to get HGNC ID if FBgn is human.
    fb_hgnc_id_query = """
    SELECT db.name||':'||dbx.accession
    FROM feature f
    JOIN feature_dbxref fdbx ON (fdbx.feature_id = f.feature_id)
    JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
    JOIN db ON (db.db_id = dbx.db_id)
    WHERE fdbx.is_current = true
      and db.name = 'HGNC'
      and f.uniquename ~ '^FB(og|gn)([0-9]{10}|[0-9]{7})$'
      and f.is_obsolete = false and f.uniquename = %s """

    # Retrieve pheno annotations from the database.
    log.info('TIME: {}. Querying database for phenotype-based annotations.'.format(now()))
    ret_fb_pheno = connect(fb_pheno_query, 'no_query', conn)
    log.info('TIME: {}. Retrieved phenotype-based annotations.'.format(now()))
    for row in ret_fb_pheno:
        disease_annotation = {
            'gene_uname': row[0],
            'gene_name': row[1],    # internal use
            'gene_sgml': None,      # add later
            'gene_org': row[2],      # internal use
            'hgnc_id': None,        # add later for Hsap gene
            'do_qual': row[3],
            'do_id': row[4],
            'do_term': row[5],
            'allele_uname': row[6],
            'allele_name': row[7],  # internal use
            'allele_sgml': None,    # add later
            'ortho_hgnc_id': None,  # not applicable
            'ortho_uname': None,    # not applicable
            'ortho_name': None,     # not applicable
            'ortho_sgml': None,     # not applicable
            'evidence': row[8],     # internal use
            'evidence_sgml': None,   # add later
            'pub': row[9]
        }
        gene_uname = (disease_annotation['gene_uname'],)
        allele_uname = (disease_annotation['allele_uname'],)
        disease_annotation['gene_sgml'] = connect(fb_synonym_sgml_query, gene_uname, conn)[0][0]
        log.debug('This is my gene_sgml: {}'.format(disease_annotation['gene_sgml']))
        if disease_annotation['gene_org'] == 'Hsap':
            ret_hgnc_id = connect(fb_hgnc_id_query, gene_uname, conn)
            if ret_hgnc_id:
                disease_annotation['hgnc_id'] = ret_hgnc_id[0][0]
            else:
                log.info('No HGNC ID found for this gene: {}'.format(gene_uname))
        log.debug('This is my HGNC ID: {}'.format(disease_annotation['hgnc_id']))
        disease_annotation['allele_sgml'] = connect(fb_synonym_sgml_query, allele_uname, conn)[0][0]
        disease_annotation['allele_sgml'] = disease_annotation['allele_sgml'].replace('<up>', '[').replace('</up>', ']')
        log.debug('This is my allele_sgml: {}'.format(disease_annotation['allele_sgml']))
        disease_annotation['evidence_sgml'] = disease_annotation['evidence'].replace('<up>', '[').replace('</up>', ']')
        log.debug('This is my evidence_sgml: {}'.format(disease_annotation['evidence_sgml']))
        to_export_as_tsv['data'].append(disease_annotation)

    # Retrieve ortho annotations from the database.
    log.info('TIME: {}. Querying database for orthology-based annotations.'.format(now()))
    ret_fb_ortho = connect(fb_ortho_query, 'no_query', conn)
    log.info('TIME: {}. Retrieved orthology-based annotations.'.format(now()))
    for row in ret_fb_ortho:
        disease_annotation = {
            'gene_uname': row[0],
            'gene_name': row[1],    # internal use
            'gene_sgml': None,      # add later
            'gene_org': 'Dmel',     # internal use
            'hgnc_id': None,        # not applicable
            'do_qual': row[2],
            'do_id': row[3],
            'do_term': row[4],
            'allele_uname': None,   # not applicable
            'allele_name': None,  # not applicable
            'allele_sgml': None,    # not applicable
            'ortho_hgnc_id': None,  # add later
            'ortho_uname': row[5],    # internal use
            'ortho_name': row[6],     # internal use
            'ortho_sgml': None,     # add later
            'evidence': row[7],     # internal use
            'evidence_sgml': None,   # add later
            'pub': row[8]
        }
        log.debug('Assessing this annotation:\n\t{}'.format(disease_annotation))
        gene_uname = (disease_annotation['gene_uname'],)
        ortho_uname = (disease_annotation['ortho_uname'],)
        disease_annotation['gene_sgml'] = connect(fb_synonym_sgml_query, gene_uname, conn)[0][0]
        log.debug('This is my gene_sgml: {}'.format(disease_annotation['gene_sgml']))
        ret_hgnc_id = connect(fb_hgnc_id_query, ortho_uname, conn)
        if ret_hgnc_id:
            disease_annotation['ortho_hgnc_id'] = ret_hgnc_id[0][0]
        else:
            log.info('No HGNC ID found for this ortholog: {}'.format(ortho_uname))
        log.debug('This is my HGNC ID: {}'.format(disease_annotation['ortho_hgnc_id']))
        try:
            ortho_sgml = connect(fb_synonym_sgml_query, ortho_uname, conn)[0][0]
        except IndexError:
            ortho_sgml = disease_annotation['ortho_name']
            log.warning('Could not find current symbol synonym for {} ({})'.format(ortho_sgml, ortho_uname))
        disease_annotation['ortho_sgml'] = re.search(r'(?<=Hsap\\).*', ortho_sgml).group(0)
        log.debug('This is my ortho_sgml: {}'.format(disease_annotation['ortho_sgml']))
        disease_annotation['evidence_sgml'] = disease_annotation['evidence'].replace('<up>', '[').replace('</up>', ']')
        log.debug('This is my evidence_sgml: {}'.format(disease_annotation['evidence_sgml']))
        to_export_as_tsv['data'].append(disease_annotation)

    conn.close()

    # Generate the output file.
    log.info('TIME: {}. Writing data to output file.'.format(now()))

    output_file = open(output_filename, 'w')

    # Write the headers.
    output_file.write('## {}\n'.format(to_export_as_tsv['metaData']['title']))
    output_file.write('## Generated: {}\n'.format(to_export_as_tsv['metaData']['dateProduced']))
    output_file.write('## Using datasource: {}\n'.format(to_export_as_tsv['metaData']['database']))
    output_file.write('##\n## ')

    headers = [
        'FBgn ID',
        'Gene symbol',
        'HGNC ID',
        'DO qualifier',
        'DO ID',
        'DO term',
        'Allele used in model (FBal ID)',
        'Allele used in model (symbol)',
        'Based on orthology with (HGNC ID)',
        'Based on orthology with (symbol)',
        'Evidence/interacting alleles',
        'Reference (FBrf ID)'
    ]

    csv_writer = csv.DictWriter(output_file, fieldnames=headers, delimiter='\t', extrasaction='ignore')
    csv_writer.writeheader()

    data_keys = [
        'gene_uname',
        'gene_sgml',
        'hgnc_id',
        'do_qual',
        'do_id',
        'do_term',
        'allele_uname',
        'allele_sgml',
        'ortho_hgnc_id',
        'ortho_sgml',
        'evidence_sgml',
        'pub'
    ]
    csv_writer = csv.DictWriter(output_file, fieldnames=data_keys, delimiter='\t', extrasaction='ignore')

    # Write the data.
    # Sort the list of data elements by gene, allele, then CV term.
    to_export_as_tsv['data'] = sorted(to_export_as_tsv['data'], key=lambda i: (i['gene_uname'], i['do_term'], i['pub']))

    for annotation in to_export_as_tsv['data']:
        csv_writer.writerow(annotation)

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
