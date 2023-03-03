# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report simple gene-publication associations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_simple_gene_publication_associations.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_simple_gene_publication_associations.py -v -c /path/to/config.cfg

"""

from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
report_label = 'gene_paper'
report_title = 'FlyBase Gene Paper Associations report'
header_list = [
    'entity_id',
    'entity_name',
    'entity_display_name',
    'FlyBase_publication_id',
    'PubMed_id'
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


def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    db_results = get_database_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['metaData']['note'] = 'Only direct associations between current D. melanogaster genes and "paper" publications are reported here.'
    data_to_export_as_tsv['data'] = process_db_results(db_results)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def make_pub_dict():
    """Make an FBrf-keyed dict of PMIDs.

    Returns:
        An FBrf-keyed dict of PubMed IDs.
        Value will be None if the FlyBase publication has no corresponding PMID.

    """
    log.info('Generating FBrf-to-PMID dict.')

    all_fb_pubs_query = """
        SELECT DISTINCT uniquename
        FROM pub
        WHERE is_obsolete is false
          AND uniquename ~ '^FBrf[0-9]{7}$'
        ;"""

    fbrf_pmid_query = """
        SELECT DISTINCT p.uniquename,
                        dbx.accession
        FROM pub p
        JOIN pub_dbxref pdbx ON pdbx.pub_id = p.pub_id
        JOIN dbxref dbx ON dbx.dbxref_id = pdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE p.is_obsolete is false
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
          AND pdbx.is_current is true
          AND db.name = 'pubmed'
          AND dbx.accession ~ '^[0-9]{1,}$'
    ;"""

    ret_all_fb_pubs = connect(all_fb_pubs_query, 'no_query', conn)
    log.info('Found {} current FB pubs.'.format(len(ret_all_fb_pubs)))
    pub_dict = {i[0]: None for i in ret_all_fb_pubs}

    ret_fbrf_pmid = connect(fbrf_pmid_query, 'no_query', conn)
    log.info('Found {} current FB pubs with PubMed ID.'.format(len(ret_fbrf_pmid)))
    FBRF_ID = 0
    PMID = 1
    for i in ret_fbrf_pmid:
        pub_dict[i[FBRF_ID]] = i[PMID]

    return pub_dict


def make_gene_symbol_dict():
    """Get current gene symbols.

    Returns:
        An FBgn-keyed dict of symbols (synonym.synonym_sgml).

    """
    log.info('Querying for current gene symbols.')
    gene_symbol_query = """
        SELECT DISTINCT f.uniquename,
                        s.synonym_sgml
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN cvterm cvtf ON cvtf.cvterm_id = f.type_id
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN cvterm cvts ON cvts.cvterm_id = s.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND cvtf.name = 'gene'
          AND o.abbreviation = 'Dmel'
          AND fs.is_current IS TRUE
          AND cvts.name = 'symbol'
    """
    ret_gene_symbols = connect(gene_symbol_query, 'no_query', conn)
    log.info('Found {} current Dmel gene symbols.'.format(len(ret_gene_symbols)))
    gene_symbol_dict = {}
    FBGN_ID = 0
    SYMBOL = 1
    for row in ret_gene_symbols:
        gene_symbol_dict[row[FBGN_ID]] = row[SYMBOL]
    return gene_symbol_dict


def get_database_info():
    """Retrieve gene-to-paper associations.

    Returns:
        A list of (entity_id, entity_name, pub_fbrf_id) tuples.

    """
    log.info('Querying database for gene-paper associations.')
    gene_paper_query = """
        SELECT DISTINCT f.uniquename,
                        f.name,
                        p.uniquename
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id AND o.abbreviation = 'Dmel'
        JOIN cvterm cvtf ON cvtf.cvterm_id = f.type_id AND cvtf.name = 'gene'
        JOIN feature_pub fp ON fp.feature_id = f.feature_id
        JOIN pub p ON p.pub_id = fp.pub_id
        JOIN cvterm cvtp ON cvtp.cvterm_id = p.type_id AND cvtp.name = 'paper'
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}$'
          AND p.is_obsolete IS FALSE
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
        ORDER BY f.uniquename,
                 p.uniquename
        ;"""

    ret_gene_papers = connect(gene_paper_query, 'no_query', conn)
    log.info('Found {} current Dmel gene-to-paper associations.'.format(len(ret_gene_papers)))
    return ret_gene_papers


def process_db_results(db_results):
    """Process db query tuples into header-keyed dicts for printing to bulk report.

    Args:
        arg1 (db_results): A list of (entity_id, pub_fbrf_id, entity_name) tuples.

    Returns:
        A list of dicts having keys that match the global "header_list" elements.

    """
    log.info('Processing db results into exportable dicts.')

    fbrf_pmid_dict = make_pub_dict()
    gene_symbol_dict = make_gene_symbol_dict()

    log.info('Sorting {} db results.'.format(len(db_results)))
    db_results = sorted(list(set(db_results)))

    log.info('Converting {} unique db result tuples to dicts.'.format(len(db_results)))
    formatted_result_list = []
    ENTITY_UNIQUENAME = 0
    ENTITY_NAME = 1
    FBRF_ID = 2
    for row in db_results:
        result_dict = {
            'entity_id': row[ENTITY_UNIQUENAME],
            'entity_name': row[ENTITY_NAME],
            'entity_display_name': gene_symbol_dict[row[ENTITY_UNIQUENAME]],
            'FlyBase_publication_id': row[FBRF_ID]
        }
        try:
            result_dict['PubMed_id'] = fbrf_pmid_dict[row[FBRF_ID]]
        except KeyError:
            continue
        formatted_result_list.append(result_dict)
    log.info('The db results have been converted to dicts.')

    return formatted_result_list


if __name__ == "__main__":
    main()
