# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Retrieve RPKM for a given project (set of samples) and return a gene by sample RPKM matrix.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_tsv_template.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_tsv_template.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect, get_features_by_uname_regex, add_unique_info, add_unique_dict_info
)
from harvdev_utils.psycopg_functions.sql_queries import (
    current_feat_symbol_sgmls, current_feat_fullname_sgmls, orgid_abbr
)

# Global variables for the output file. Header order will match list order below.
report_label = 'gene_rpkm_matrix'
report_title = 'FlyBase gene RPKM matrix report'
header_list = [
    'gene_primary_id',
    'gene_symbol',
    'gene_fullname',
    'gene_type'
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
parser.add_argument('-p', '--prefix', help='Prefix to remove from start of sample names.', required=False)
parser.add_argument('-s', '--suffix', help='Suffix to remove from end of sample names.', required=False)

# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
dataset_name_prefix = args.prefix
dataset_name_suffix = args.suffix


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    gene_info = get_fb_gene_info()
    gene_info = get_rpkm_data(gene_info)

    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(gene_info)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of specific data types.
def get_fb_gene_info():
    """Retrieve FlyBase genes and return as an FBgn-keyed dict of Gene class objects.

    Returns:
        An FBgn-keyed dict of Gene class objects.

    """
    log.info('Getting list of FlyBase genes.')

    # Query for gene biotype.
    fb_gene_type_query = """
        SELECT DISTINCT f.uniquename,
                        fp.value
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fp.type_id
        WHERE f.is_obsolete = false and
              f.is_analysis = false and
              f.uniquename ~ '{}' and
              o.abbreviation = 'Dmel' and
              cvt.name = 'promoted_gene_type';
        """

    # Get all genes.
    fbgnwc = '^FBgn[0-9]{7}$'
    feat_dict = get_features_by_uname_regex(conn, fbgnwc)
    feat_dict = add_unique_dict_info(feat_dict, 'organism_id', 'org_abbr', conn, orgid_abbr)

    # Filter for Dmel genes.
    dmel_feat_dict = {}
    for gene in feat_dict.values():
        if gene.org_abbr == 'Dmel':
            dmel_feat_dict[gene.uniquename] = gene

    # Get basic gene info.
    dmel_feat_dict = add_unique_info(dmel_feat_dict, 'gene_symbol', conn, current_feat_symbol_sgmls, fbgnwc)
    dmel_feat_dict = add_unique_info(dmel_feat_dict, 'gene_fullname', conn, current_feat_fullname_sgmls, fbgnwc)
    dmel_feat_dict = add_unique_info(dmel_feat_dict, 'gene_type', conn, fb_gene_type_query, fbgnwc)
    for gene in dmel_feat_dict.values():
        gene.gene_symbol = gene.gene_symbol.replace('<up>', '[').replace('</up>', ']')    # Just for su(w[a]) gene.
        gene.gene_primary_id = gene.uniquename    # Need to re-assign this attribute to one that matches the header label.
        if hasattr(gene, 'gene_type'):
            gene.gene_type = gene.gene_type.split(':')[1].replace('@', '')    # Strip out SO ID, etc.

    return dmel_feat_dict


def get_rpkm_data(feat_dict):
    """Get RPKM data for a set of genes.

    Args:
        arg1 (dict): An FBgn-keyed dict of Gene objects.

    Returns:
        A dict of Gene objects with RPKM added, each as an attribute using dataset sample name.

    """
    log.info('Getting RPKM data for genes.')
    fbgnwc = '^FBgn[0-9]{7}$'
    fblcwc = '^FBlc[0-9]{7}$'

    # First get an ordered list of the relevant datasets.
    fb_rpkm_dataset_query = """
        SELECT DISTINCT p.uniquename, p.name, l.uniquename, l.name
        FROM library l
        JOIN organism o ON o.organism_id = l.organism_id
        JOIN cvterm cvtl ON cvtl.cvterm_id = l.type_id
        JOIN library_relationship lr ON lr.subject_id = l.library_id
        JOIN library p ON p.library_id = lr.object_id
        JOIN cvterm cvtp ON cvtp.cvterm_id = p.type_id
        JOIN cvterm cvtlr ON cvtlr.cvterm_id = lr.type_id
        JOIN library_feature lf ON lf.library_id = l.library_id
        JOIN library_featureprop lfp ON lfp.library_feature_id = lf.library_feature_id
        JOIN cvterm cvtlfp ON cvtlfp.cvterm_id = lfp.type_id
        WHERE l.is_obsolete = false and
              l.uniquename ~ '{}' and
              cvtl.name = 'result' and
              o.abbreviation = 'Dmel' and
              p.is_obsolete = false and
              p.uniquename ~ '{}' and
              cvtp.name = 'project' and
              cvtlr.name = 'belongs_to' and
              cvtlfp.name = 'RPKM'
        ORDER BY p.uniquename, l.uniquename;
        """
    ret_fb_rpkm_dataset = connect(fb_rpkm_dataset_query.format(fblcwc, fblcwc), 'no_query', conn)
    DATASET_ID = 2
    DATASET_NAME = 3
    dataset_dict = {i[DATASET_ID]: '{}_({})'.format(i[DATASET_NAME], i[DATASET_ID]) for i in ret_fb_rpkm_dataset}
    header_list.extend(dataset_dict.values())

    # Now get RPKM for all genes for each dataset.
    fb_gene_rpkm_by_dataset_query = """
        SELECT DISTINCT f.uniquename,
                        lfp.value
        FROM feature f
        JOIN library_feature lf ON lf.feature_id = f.feature_id
        JOIN library l ON l.library_id = lf.library_id
        JOIN library_featureprop lfp ON lfp.library_feature_id = lf.library_feature_id
        JOIN cvterm cvtlfp ON cvtlfp.cvterm_id = lfp.type_id
        WHERE f.is_obsolete = false and
              f.is_analysis = false and
              f.uniquename ~ '{}' and
              l.uniquename = '{}' and
              cvtlfp.name = 'RPKM';
        """
    # for dataset in dataset_list:
    for fblc_id in dataset_dict.keys():
        dataset = dataset_dict[fblc_id]
        log.info('Getting RPKM data for {}.'.format(dataset))
        feat_dict = add_unique_info(feat_dict, dataset, conn, fb_gene_rpkm_by_dataset_query, fbgnwc, fblc_id)

    # Filter for genes having at least one RPKM value (no sense reporting gene with nothing but "NA" values).
    filtered_feat_dict = {}
    for gene in feat_dict.values():
        rpkm_counter = 0
        for dataset in dataset_dict.values():
            if hasattr(gene, dataset):
                rpkm_counter += 1
        if rpkm_counter > 0:
            filtered_feat_dict[gene.uniquename] = gene

    return filtered_feat_dict


def process_database_info(feat_dict):
    """Convert dict of Gene objects into list of gene dicts for printing.

    Args:
        arg1 (dict): FBgn-keyed dict of Gene objects having RPKM data.

    Returns:
        A list of gene dicts with keys matching the "header list".
        Header list includes basic gene info, as well as RNA-seq sample names.

    """
    log.info('Processing dict of Genes into data list for printing.')

    data_list = []
    for gene in feat_dict.values():
        gene_dict = {}
        for attribute in header_list:
            try:
                gene_dict[attribute] = getattr(gene, attribute)
            except AttributeError:
                if attribute == 'gene_fullname':
                    gene_dict[attribute] = ""
                else:
                    gene_dict[attribute] = "NA"
        data_list.append(gene_dict)

    return data_list


if __name__ == "__main__":
    main()
