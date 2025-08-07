# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report representative publications.

- ORDER PUBS BY DESC SCORE, DESC fbRF id
- CUT-OFF AT 100.
- Use (score, FBrf ID) tuples as keys to sort FBrf IDs.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_representative_publications.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_representative_publications.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
REPORT_LABEL = 'representative_publications'
REPORT_TITLE = 'FlyBase Representative Publications Report'
HEADER_LIST = [
    'FBgn_ID',
    'Symbol',
    'References',
]

# Proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
DATABASE = set_up_dict['database']
OUTPUT_FILENAME = set_up_dict['output_filename']
log = set_up_dict['log']
CONN = set_up_dict['conn']

# Process more input parameters (-c and -v handled by set_up_db_reading() function above).
parser = argparse.ArgumentParser(description='inputs')
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    gene_dict = get_dmel_genes()
    pmid_dict = get_pmids()
    get_ranked_pubs(gene_dict, pmid_dict)
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = process_database_info(gene_dict)
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of data from chado.
def get_dmel_genes():
    """Retrieve Dmel genes."""
    global CONN
    log.info('Retrieve Dmel genes.')
    fb_dmel_gene_query = """
        SELECT DISTINCT f.feature_id, f.uniquename, f.name
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND o.abbreviation = 'Dmel'
        ORDER BY f.uniquename;
    """
    ret_dmel_gene_info = connect(fb_dmel_gene_query, 'no_query', CONN)
    DB_ID = 0
    UNAME = 1
    NAME = 2
    gene_dict = {}
    counter = 0
    for row in ret_dmel_gene_info:
        gene_result = {
            'db_id': row[DB_ID],
            'FBgn_ID': row[UNAME],
            'Symbol': row[NAME],
            'All_References': [],
            'References': [],
        }
        gene_dict[row[DB_ID]] = gene_result
        counter += 1
    log.info(f'Found {counter} Dmel genes in chado.')
    return gene_dict


def get_pmids():
    """Retrieve PubMed IDs."""
    global CONN
    log.info('Retrieve PubMed IDs.')
    fb_pmid_query = """
        SELECT DISTINCT p.uniquename, 'PMID:'||dbx.accession
        FROM pub p
        JOIN pub_dbxref pdbx ON pdbx.pub_id = p.pub_id
        JOIN dbxref dbx ON dbx.dbxref_id = pdbx.dbxref_id
        JOIN db d ON d.db_id = dbx.db_id AND d.name = 'pubmed'
        WHERE p.is_obsolete IS FALSE
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
          AND pdbx.is_current IS TRUE;
    """
    ret_pmid_info = connect(fb_pmid_query, 'no_query', CONN)
    FBRF_ID = 0
    PMID = 1
    pmid_dict = {}
    for row in ret_pmid_info:
        pmid_dict[row[FBRF_ID]] = row[PMID]
    log.info(f'Found {len(pmid_dict)} PubMed IDs in chado.')
    return pmid_dict


def get_ranked_pubs(gene_dict, pmid_dict):
    """Retrieve per gene ranked pub scores."""
    global CONN
    log.info('Retrieve per gene ranked pub scores.')
    fb_ranked_pub_query = """
        SELECT DISTINCT f.feature_id, p.uniquename, fpp.value
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN feature_pub fp ON fp.feature_id = f.feature_id
        JOIN pub p ON p.pub_id = fp.pub_id
        JOIN feature_pubprop fpp ON fpp.feature_pub_id = fp.feature_pub_id
        JOIN cvterm cvt ON cvt.cvterm_id = fpp.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND o.abbreviation = 'Dmel'
          AND p.is_obsolete IS FALSE
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
          AND cvt.name = 'computed_gene_pub_score';
    """
    ret_gene_pub_info = connect(fb_ranked_pub_query, 'no_query', CONN)
    GENE_ID = 0
    PUB_ID = 1
    SCORE = 2
    counter = 0
    for row in ret_gene_pub_info:
        score = float(row[SCORE])
        if row[PUB_ID] in pmid_dict.keys():
            pub_id = f'{row[PUB_ID]}|{pmid_dict[row[PUB_ID]]}'
        else:
            pub_id = f'{row[PUB_ID]}|-'
        ranked_pub = (score, pub_id)
        try:
            gene_dict[row[GENE_ID]]['All_References'].append(ranked_pub)
        except KeyError:
            log.warning(f'Gene ID {row[GENE_ID]} not found in gene_dict.')
        counter += 1
    log.info(f'Found {counter} gene-pub scores in chado.')
    return


def process_database_info(input_data):
    """Convert the gene dict with ranked pub info to a list of data elements."""
    log.info('Convert the gene dict with ranked pub info to a list of data elements.')
    data_list = []
    counter = 0
    for gene in input_data.values():
        sorted_references = list(set(gene['All_References']))
        sorted_references.sort(reverse=True)  # Sort by descending score, then by descending FBrf ID.
        rep_pub_counter = 0
        for pub in sorted_references:
            if rep_pub_counter >= 100:
                break
            gene['References'].append(pub[1])  # Append FBrf ID.
            rep_pub_counter += 1
        if not gene['References']:
            gene['References'] = ''
        else:
            gene['References'] = ','.join(gene['References'])
        data_list.append(gene)
        counter += 1
    log.info(f'Sending {counter} genes with representative pubs to the export file.')
    return data_list


if __name__ == "__main__":
    main()
