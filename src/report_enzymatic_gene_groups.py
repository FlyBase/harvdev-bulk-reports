# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create TSV report on ENZYMATIC gene groups and their associated genes.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_enzymatic_gene_groups.py.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_enzymatic_gene_groups.py.py -v -c /path/to/config.cfg

"""

from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect    # other useful functions: add_unique_info, add_list_info, add_unique_dict_info
)

# Global variables for the output file. Header order will match list order below.
report_label = 'Dmel_enzyme_data'
report_title = 'FlyBase D. melanogaster enzyme data for genes and gene groups'
header_list = [
    'gene_group_id',
    'gene_group_name',
    'gene_group_GO_id(s)',
    'gene_group_GO_name(s)',
    'gene_group_EC_number(s)',
    'gene_group_EC_name(s)',
    'gene_id',
    'gene_symbol',
    'gene_name',
    'gene_EC_number(s)',
    'gene_EC_name(s)'
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


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    fbgg_fbgn_list = get_gene_groups()
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(fbgg_fbgn_list)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def get_gene_groups():
    """Get list of all enzymatic gene groups and their associated genes.

    Returns:
        A list of tuples having this form.
        ('FBgg1234567', 'Random enzymatic gene group', 'FBgn0000001')

    """
    log.info('Querying database for enzymatic gene group info.')
    fbgg_fbgn_query = """
        SELECT DISTINCT
            s.uniquename,
            title.name,
            f.uniquename
        FROM grp s
        JOIN grp_synonym gs ON (gs.grp_id = s.grp_id AND gs.is_current is true)
        JOIN synonym title ON title.synonym_id = gs.synonym_id
        JOIN cvterm cvts ON (cvts.cvterm_id = title.type_id AND cvts.name = 'fullname')
        JOIN grp_relationship gr ON gr.subject_id = s.grp_id
        JOIN grp o ON o.grp_id = gr.object_id
        JOIN cvterm cvt ON cvt.cvterm_id = gr.type_id
        JOIN grpmember gm ON gm.grp_id = s.grp_id
        JOIN feature_grpmember fgm ON fgm.grpmember_id = gm.grpmember_id
        JOIN feature f ON f.feature_id = fgm.feature_id
        LEFT OUTER JOIN grpmemberprop gmp ON gmp.grpmember_id = gm.grpmember_id
        WHERE s.is_obsolete is false
          AND f.is_obsolete is false
          AND f.uniquename ~ '^FBgn[0-9]{7}'
          AND cvt.name = 'component_grp'
          AND o.uniquename = 'FBgg0001715'
          AND o.name = 'ENZ'
          AND gmp.value IS NULL
        ORDER BY s.uniquename, f.uniquename;
    """
    fbgg_fbgn_results = connect(fbgg_fbgn_query, 'no_query', conn)
    log.info('Found {} enzymatic gene_group-gene relationships .'.format(len(fbgg_fbgn_results)))
    return fbgg_fbgn_results


def get_fbgg_go_mf_data():
    """Get GO MF terms associated with gene groups.

    Returns:
        A FBgg-keyed dict of lists of GO MF cvterm_ids associated with the gene group.

    """
    log.info('Looking up GO MF terms associated with gene groups.')
    fbgg_go_mf_query = """
        SELECT DISTINCT
            grp.uniquename,
            cvt.cvterm_id
        FROM grp
        JOIN grp_cvterm gcvt ON gcvt.grp_id = grp.grp_id
        JOIN cvterm cvt ON cvt.cvterm_id = gcvt.cvterm_id
        JOIN cv ON (cv.cv_id = cvt.cv_id AND cv.name = 'molecular_function')
        WHERE grp.is_obsolete is false
          AND gcvt.is_not is false
          AND cvt.is_obsolete = 0;
    """
    fbgg_go_mf_results = connect(fbgg_go_mf_query, 'no_query', conn)
    log.info('Found {} FBgg-GO MF term associations.'.format(len(fbgg_go_mf_results)))
    fbgg_go_mf_dict = {}
    FBGG_ID = 0
    CVTERM_ID = 1
    for i in fbgg_go_mf_results:
        try:
            fbgg_go_mf_dict[i[FBGG_ID]].append(i[CVTERM_ID])
        except KeyError:
            fbgg_go_mf_dict[i[FBGG_ID]] = [i[CVTERM_ID]]
    return fbgg_go_mf_dict


def get_fbgn_go_mf_data():
    """Get GO MF terms associated with genes.

    Returns:
        A FBgn-keyed dict of lists of GO MF cvterm_ids associated with the gene.

    """
    log.info('Looking up GO MF terms associated with genes.')
    fbgn_go_mf_query = """
        SELECT DISTINCT
            g.uniquename,
            cvt.cvterm_id
        FROM feature g
        JOIN feature_cvterm fcvt ON fcvt.feature_id = g.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
        JOIN cv ON (cv.cv_id = cvt.cv_id AND cv.name = 'molecular_function')
        WHERE g.is_obsolete is false
          AND g.uniquename ~ '^FBgn[0-9]{7}$'
          AND fcvt.is_not is false
          AND cvt.is_obsolete = 0;
    """
    fbgn_go_mf_results = connect(fbgn_go_mf_query, 'no_query', conn)
    log.info('Found {} FBgn-GO MF term associations.'.format(len(fbgn_go_mf_results)))
    fbgn_go_mf_dict = {}
    FBGN_ID = 0
    CVTERM_ID = 1
    for i in fbgn_go_mf_results:
        try:
            fbgn_go_mf_dict[i[FBGN_ID]].append(i[CVTERM_ID])
        except KeyError:
            fbgn_go_mf_dict[i[FBGN_ID]] = [i[CVTERM_ID]]
    return fbgn_go_mf_dict


def get_negative_fbgn_go_mf_data():
    """Get negative GO MF gene annotations.

    Returns:
        A FBgn-keyed dict of lists of GO MF cvterm_ids negatively associated with the gene.

    """
    log.info('Looking up negative GO MF gene annotations.')
    neg_fbgn_go_mf_query = """
        SELECT DISTINCT
            g.uniquename,
            cvt.cvterm_id
        FROM feature g
        JOIN feature_cvterm fcvt ON fcvt.feature_id = g.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
        JOIN cv ON (cv.cv_id = cvt.cv_id AND cv.name = 'molecular_function')
        WHERE g.is_obsolete is false
          AND g.uniquename ~ '^FBgn[0-9]{7}$'
          AND fcvt.is_not is true
          AND cvt.is_obsolete = 0;
    """
    neg_fbgn_go_mf_results = connect(neg_fbgn_go_mf_query, 'no_query', conn)
    log.info('Found {} negative FBgn-GO MF term associations.'.format(len(neg_fbgn_go_mf_results)))
    neg_fbgn_go_mf_dict = {}
    FBGN_ID = 0
    CVTERM_ID = 1
    for i in neg_fbgn_go_mf_results:
        try:
            neg_fbgn_go_mf_dict[i[FBGN_ID]].append(i[CVTERM_ID])
        except KeyError:
            neg_fbgn_go_mf_dict[i[FBGN_ID]] = [i[CVTERM_ID]]
    return neg_fbgn_go_mf_dict


def get_go_mf_term_info():
    """Get GO MF info.

    Returns:
        A cvterm_id-keyed dict of a list of GO MF info. List elements are dicts of this form:
        {
            'go_id': 'GO:12345678',
            'go_name': 'randomase activity'
        }

    """
    log.info('Getting GO MF information.')
    go_mf_query = """
        SELECT DISTINCT
            cvt.cvterm_id,
            cvt.name,
            db.name||':'||dbx.accession
        FROM cvterm cvt
        JOIN cv ON cv.cv_id = cvt.cv_id
        JOIN dbxref dbx ON dbx.dbxref_id = cvt.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE cvt.is_obsolete = 0
          AND cv.name = 'molecular_function'
          AND db.name = 'GO';
    """
    go_mf_results = connect(go_mf_query, 'no_query', conn)
    log.info('Found {} GO MF terms.'.format(len(go_mf_results)))
    go_mf_dict = {}
    CVTERM_ID = 0
    GO_NAME = 1
    GO_ID = 2
    for i in go_mf_results:
        go_mf_dict[i[CVTERM_ID]] = {
            'go_id': i[GO_ID],
            'go_name': i[GO_NAME]
        }
    return go_mf_dict


def get_go_ec_info():
    """Get EC data for all GO MF terms.

    Returns:
        A cvterm_id-keyed dict of a list of related EC tuples (ec_number, ec_name).

    """
    log.info('Getting GO-EC information.')
    go_ec_query = """
        SELECT DISTINCT
            cvt.cvterm_id,
            dbx.accession,
            dbxp.value
        FROM cvterm cvt
        JOIN cv ON cv.cv_id = cvt.cv_id
        JOIN cvterm_dbxref cvtdbx ON cvtdbx.cvterm_id = cvt.cvterm_id
        JOIN dbxref dbx ON dbx.dbxref_id = cvtdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        LEFT OUTER JOIN dbxrefprop dbxp ON (dbxp.dbxref_id = dbx.dbxref_id AND dbxp.type_id = (SELECT cvterm_id FROM cvterm WHERE name = 'ec_description'))
        WHERE cvt.is_obsolete = 0
          AND cv.name = 'molecular_function'
          AND db.name = 'EC';
    """
    go_ec_results = connect(go_ec_query, 'no_query', conn)
    log.info('Found {} GO-EC relationships.'.format(len(go_ec_results)))
    go_ec_dict = {}
    CVTERM_ID = 0
    EC_NUMBER = 1
    EC_NAME = 2
    for i in go_ec_results:
        ec_number = i[EC_NUMBER]
        ec_name = i[EC_NAME]
        if ec_name is not None:
            if ec_name.endswith('.'):
                ec_name = ec_name[:-1]
        ec_tuple = (ec_number, ec_name)
        try:
            go_ec_dict[i[CVTERM_ID]].append(ec_tuple)
        except KeyError:
            go_ec_dict[i[CVTERM_ID]] = [ec_tuple]
    return go_ec_dict


def get_gene_symbol_info():
    """Generate an FBgn-keyed dict of gene symbols."""
    log.info('Getting gene symbol info.')
    fbgn_symbol_query = """
        SELECT DISTINCT f.uniquename,
                        s.name
        FROM feature f
        JOIN feature_synonym fs ON (fs.feature_id = f.feature_id AND fs.is_current is true)
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN cvterm cvt ON (cvt.cvterm_id = s.type_id AND cvt.name = 'symbol')
    """
    fbgn_symbol_results = connect(fbgn_symbol_query, 'no_query', conn)
    log.info('Found {} gene symbols.'.format(len(fbgn_symbol_results)))
    fbgn_symbol_dict = {}
    FBGN_ID = 0
    SYMBOL = 1
    for i in fbgn_symbol_results:
        fbgn_symbol_dict[i[FBGN_ID]] = i[SYMBOL]
    return fbgn_symbol_dict


def get_gene_fullname_info():
    """Generate an FBgn-keyed dict of gene full names."""
    log.info('Getting gene symbol info.')
    fbgn_symbol_query = """
        SELECT DISTINCT f.uniquename,
                        s.name
        FROM feature f
        JOIN feature_synonym fs ON (fs.feature_id = f.feature_id AND fs.is_current is true)
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN cvterm cvt ON (cvt.cvterm_id = s.type_id AND cvt.name = 'fullname')
    """
    fbgn_symbol_results = connect(fbgn_symbol_query, 'no_query', conn)
    log.info('Found {} gene symbols.'.format(len(fbgn_symbol_results)))
    fbgn_fullname_dict = {}
    FBGN_ID = 0
    FULLNAME = 1
    for i in fbgn_symbol_results:
        fbgn_fullname_dict[i[FBGN_ID]] = i[FULLNAME]
    return fbgn_fullname_dict


def process_database_info(input_data):
    """Convert a list of input dicts into a list of dicts for TSV output.

    Args:
        arg1 (list): A list of dicts representing FBgg-FBgn associations.

    Returns:
        A list of dictionaries having FBgg-FBgn associations and all relevant GO MF and EC info.
    """
    log.info('Starting to process FBgg-FBgn info retrieved from database.')
    # Get needed dicts.
    fbgg_go_mf_dict = get_fbgg_go_mf_data()
    fbgn_go_mf_dict = get_fbgn_go_mf_data()
    neg_fbgn_go_mf_dict = get_negative_fbgn_go_mf_data()
    go_mf_term_dict = get_go_mf_term_info()
    go_ec_dict = get_go_ec_info()
    fbgn_symbol_dict = get_gene_symbol_info()
    fbgn_fullname_dict = get_gene_fullname_info()

    # Start data list.
    data_list = []
    FBGG_ID = 0
    FBGG_NAME = 1
    FBGN_ID = 2
    for i in input_data:
        data_dict = {
            'gene_group_id': i[FBGG_ID],
            'gene_group_name': i[FBGG_NAME],
            'gene_id': i[FBGN_ID],
            'gene_symbol': fbgn_symbol_dict[i[FBGN_ID]],
            'gene_name': None,
            'gene_group_GO_id(s)': None,
            'gene_group_GO_name(s)': None,
            'gene_group_EC_number(s)': None,
            'gene_group_EC_name(s)': None,
            'gene_EC_number(s)': None,
            'gene_EC_name(s)': None
        }
        try:
            data_dict['gene_name'] = fbgn_fullname_dict[i[FBGN_ID]]
        except KeyError:
            pass
        # Constants for handling EC info.
        EC_NUMBER = 0
        EC_NAME = 1
        # Get gene group GO MF and EC info.
        if i[FBGG_ID] in fbgg_go_mf_dict.keys():
            gene_group_go_mf_ids = [go_mf_term_dict[cvterm_id]['go_id'] for cvterm_id in fbgg_go_mf_dict[i[FBGG_ID]]]
            data_dict['gene_group_GO_id(s)'] = '|'.join(gene_group_go_mf_ids)
            gene_group_go_mf_names = [go_mf_term_dict[cvterm_id]['go_name'] for cvterm_id in fbgg_go_mf_dict[i[FBGG_ID]]]
            data_dict['gene_group_GO_name(s)'] = '|'.join(gene_group_go_mf_names)
            # Get EC for gene groups.
            gene_group_ec_list = []
            for cvterm_id in fbgg_go_mf_dict[i[FBGG_ID]]:
                try:
                    gene_group_ec_list.extend(go_ec_dict[cvterm_id])
                except KeyError:
                    pass
            if gene_group_ec_list:
                gene_group_ec_set = set(gene_group_ec_list)
                gene_group_ec_numbers = [ec[EC_NUMBER] for ec in gene_group_ec_set]
                data_dict['gene_group_EC_number(s)'] = '|'.join(gene_group_ec_numbers)
                gene_group_ec_names = [ec[EC_NAME] for ec in gene_group_ec_set if ec[EC_NAME] is not None]
                data_dict['gene_group_EC_name(s)'] = '|'.join(gene_group_ec_names)
        # Now get EC info for genes.
        if i[FBGN_ID] in fbgn_go_mf_dict.keys():
            gene_ec_list = []
            all_positive_annotations = set(fbgn_go_mf_dict[i[FBGN_ID]])
            try:
                filtered_annotations = all_positive_annotations - set(neg_fbgn_go_mf_dict[i[FBGN_ID]])
            except KeyError:
                filtered_annotations = all_positive_annotations
            for cvterm_id in filtered_annotations:
                try:
                    gene_ec_list.extend(go_ec_dict[cvterm_id])
                except KeyError:
                    pass
            if gene_ec_list:
                gene_ec_set = set(gene_ec_list)
                gene_ec_numbers = []
                gene_ec_names = []
                for ec in gene_ec_set:
                    # Filter out generic EC numbers (ending with "-") and specific EC numbers with no name.
                    if ec[EC_NAME]:
                        gene_ec_numbers.append(ec[EC_NUMBER])
                        gene_ec_names.append(ec[EC_NAME])
                data_dict['gene_EC_number(s)'] = '|'.join(gene_ec_numbers)
                data_dict['gene_EC_name(s)'] = '|'.join(gene_ec_names)
        data_list.append(data_dict)
    log.info('Done processing FBgg-FBgn info.')
    return data_list


if __name__ == "__main__":
    main()
