# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report entity-publication associations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_entity_publication_associations.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_entity_publication_associations.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
report_label = 'entity_publication'
report_title = 'FlyBase Entity Publication Associations report'
header_list = [
    'entity_id',
    'entity_name',
    'FlyBase_publication_id',
    'PubMed_id'
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
parser.add_argument('-i', '--input_filename', help='Input TSV file.', required=False)
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
input_filename = args.input_filename


def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    db_results = get_database_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['metaData']['note'] = 'The associations in this file may be more expansive than those reported on FlyBase Reference web reports.'
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


def get_database_info():
    """Retrieve entity-to-publication associations.

    Returns:
        A list of (entity_id, pub_fbrf_id, entity_name) tuples.

    """
    log.info('Querying database for entity-publication associations.')

    entity_regex = '^FB[a-z]{2}[0-9]{1,10}$'
    excluded_entity_regex = '^FB(og|pp|tr)'    # Features that are excluded, or have special queries.
    pub_regex = '^FBrf[0-9]{7}$'

    # Various queries used, some specific, some generic.
    cell_line_pub_query = """
        SELECT DISTINCT cl.uniquename,
                        p.uniquename,
                        cl.name
        FROM cell_line cl
        JOIN cell_line_pub clp ON clp.cell_line_id = cl.cell_line_id
        JOIN pub p ON p.pub_id = clp.pub_id
        WHERE cl.uniquename ~ '^FBtc[0-9]{7}$'
          AND p.is_obsolete is false
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
        ORDER BY cl.uniquename,
                 p.uniquename
        ;"""

    gene_product_pub_query = """
        SELECT DISTINCT f.uniquename,
                        p.uniquename,
                        f.name
        FROM feature f
        JOIN featureloc fl ON fl.feature_id = f.feature_id
        JOIN feature_pub fp ON fp.feature_id = f.feature_id
        JOIN pub p ON p.pub_id = fp.pub_id
        WHERE f.is_obsolete is false
          AND f.uniquename ~ '^FB(tr|pp)[0-9]{7}'
          AND p.is_obsolete is false
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
        ORDER BY f.uniquename,
                 p.uniquename
    ;"""

    interaction_group_pub_query = """
        SELECT DISTINCT ig.uniquename,
                        p.uniquename,
                        ig.description
        FROM interaction i
        JOIN interaction_pub ip ON ip.interaction_id = i.interaction_id
        JOIN pub p ON p.pub_id = ip.pub_id
        JOIN feature_interaction fi ON fi.interaction_id = i.interaction_id
        JOIN interaction_group_feature_interaction igfi ON igfi.feature_interaction_id = fi.feature_interaction_id
        JOIN interaction_group ig ON ig.interaction_group_id = igfi.interaction_group_id
        WHERE i.is_obsolete is false
          AND ig.is_obsolete is false
          AND ig.uniquename ~ '^FBig[0-9]{10}$'
          AND p.is_obsolete is false
          AND p.uniquename ~ '^FBrf[0-9]{7}$'
        ORDER BY ig.uniquename,
                 p.uniquename
        ;"""

    # Expression is annotated to internal features - this pulls out the related public features, usually FBgn/FBal.
    # Reporter expression will pull in appropriate genes via "attributed_as_expression_of" feature_relationships.
    feature_expression_pub_query = """
        SELECT DISTINCT o.uniquename,
                        p.uniquename,
                        o.name
        FROM feature s
        JOIN feature_expression fe ON fe.feature_id = s.feature_id
        JOIN pub p ON p.pub_id = fe.pub_id
        JOIN feature_relationship fr ON fr.subject_id = s.feature_id
        JOIN feature o ON o.feature_id = fr.object_id
        WHERE s.is_obsolete is false
          AND o.is_obsolete is false
          AND o.uniquename ~ '{entity_regex}'
          AND NOT o.uniquename ~ '{excluded_entity_regex}'
          AND p.is_obsolete is false
          AND p.uniquename ~ '{pub_regex}'
        ;""".format(entity_regex=entity_regex, excluded_entity_regex=excluded_entity_regex, pub_regex=pub_regex)

    generic_entity_pub_query = """
        SELECT DISTINCT {entity}.uniquename,
                        p.uniquename,
                        {entity}.name
        FROM {entity}
        JOIN {entity}_{attribute_type} ON {entity}_{attribute_type}.{entity}_id = {entity}.{entity}_id
        JOIN pub p ON p.pub_id = {entity}_{attribute_type}.pub_id
        WHERE {entity}.is_obsolete is false
          AND {entity}.uniquename ~ '{entity_regex}'
          AND NOT {entity}.uniquename ~ '{excluded_entity_regex}'
          AND p.is_obsolete is false
          AND p.uniquename ~ '{pub_regex}'
        ORDER BY {entity}.uniquename,
                 p.uniquename
    ;"""

    generic_entityprop_pub_query = """
        SELECT DISTINCT {entity}.uniquename,
                        p.uniquename,
                        {entity}.name
        FROM {entity}
        JOIN {entity}prop ON {entity}prop.{entity}_id = {entity}.{entity}_id
        JOIN {entity}prop_pub ON {entity}prop_pub.{entity}prop_id = {entity}prop.{entity}prop_id
        JOIN pub p ON p.pub_id = {entity}prop_pub.pub_id
        WHERE {entity}.is_obsolete is false
          AND {entity}.uniquename ~ '{entity_regex}'
          AND NOT {entity}.uniquename ~ '{excluded_entity_regex}'
          AND p.is_obsolete is false
          AND p.uniquename ~ '{pub_regex}'
        ORDER BY {entity}.uniquename,
                 p.uniquename
        ;"""

    generic_rel_subject_query = """
        SELECT DISTINCT s.uniquename,
                        p.uniquename,
                        s.name
        FROM {entity} s
        JOIN {entity}_relationship rel ON rel.subject_id = s.{entity}_id
        JOIN {entity} o ON o.{entity}_id = rel.object_id
        JOIN {entity}_relationship_pub relpub ON relpub.{entity}_relationship_id = rel.{entity}_relationship_id
        JOIN pub p ON p.pub_id = relpub.pub_id
        WHERE s.is_obsolete is false
          AND s.uniquename ~ '{entity_regex}'
          AND NOT s.uniquename ~ '{excluded_entity_regex}'
          AND o.is_obsolete is false
          AND p.is_obsolete is false
          AND p.uniquename ~ '{pub_regex}'
        ;"""

    generic_rel_object_query = """
        SELECT DISTINCT o.uniquename,
                        p.uniquename,
                        o.name
        FROM {entity} s
        JOIN {entity}_relationship rel ON rel.subject_id = s.{entity}_id
        JOIN {entity} o ON o.{entity}_id = rel.object_id
        JOIN {entity}_relationship_pub relpub ON relpub.{entity}_relationship_id = rel.{entity}_relationship_id
        JOIN pub p ON p.pub_id = relpub.pub_id
        WHERE s.is_obsolete is false
          AND o.is_obsolete is false
          AND o.uniquename ~ '{entity_regex}'
          AND NOT o.uniquename ~ '{excluded_entity_regex}'
          AND p.is_obsolete is false
          AND p.uniquename ~ '{pub_regex}'
        ;"""

    # Comments below indicate low-priority/infrequent data types that are probably missed.
    # Primarily these are database entities referenced in the proforma of some other entity type.
    # For example, cell_line_library is curated in DATASET proforma.
    #     So, while library gets a direct library_pub, the cell_line does not.
    #     Hence, for cell_line_library data, library is probably covered by generic queries, but cell_line not.
    entity_dict = {
        'cell_line': cell_line_pub_query,       # misses: cell_line_library, interaction_cell_line
        'gene_product': gene_product_pub_query,
        'interaction_group': interaction_group_pub_query,
        'feature_expression': feature_expression_pub_query,
        'relationship': 'generic_relationship',
        'feature': generic_entity_pub_query,    # misses: cell_line_feature, humanhealth_feature, some library_feature
        'grp': generic_entity_pub_query,
        'humanhealth': generic_entity_pub_query,
        'library': generic_entity_pub_query,    # misses: library_interaction
        'strain': generic_entity_pub_query      # misses: cell_line_strain, library_strain
    }
    attribute_types = ['pub', 'synonym', 'cvterm']
    relationship_types = ['feature', 'grp', 'humanhealth', 'library', 'strain']    # cell_line_relationship has no pub info.

    all_query_results = []
    for entity, entity_query in entity_dict.items():
        log.info('Querying {} data.'.format(entity))
        if entity_query == generic_entity_pub_query:
            for attribute_type in attribute_types:
                log.info('Querying {} entity, {} attribute.'.format(entity, attribute_type))
                formatted_entity_query = entity_query.format(entity=entity, attribute_type=attribute_type, entity_regex=entity_regex,
                                                             excluded_entity_regex=excluded_entity_regex, pub_regex=pub_regex)
                entity_query_results = connect(formatted_entity_query, 'no_query', conn)
                log.info('Found {} entity-pub rows for {}_{} query.'.format(len(entity_query_results), entity,
                                                                            attribute_type))
                all_query_results.extend(entity_query_results)
            # And for generic scenario, use also generic prop query.
            log.info('Querying {}prop.'.format(entity))
            formatted_entity_query = generic_entityprop_pub_query.format(entity=entity, entity_regex=entity_regex,
                                                                         excluded_entity_regex=excluded_entity_regex, pub_regex=pub_regex)
            entity_query_results = connect(formatted_entity_query, 'no_query', conn)
            all_query_results.extend(entity_query_results)
        elif entity_query == 'generic_relationship':
            for entity in relationship_types:
                log.info('Querying {}_relationship.'.format(entity))
                formatted_entity_query = generic_rel_subject_query.format(entity=entity, entity_regex=entity_regex,
                                                                          excluded_entity_regex=excluded_entity_regex,
                                                                          pub_regex=pub_regex)
                rel_subject_query_results = connect(formatted_entity_query, 'no_query', conn)
                log.info('Found {} subject-pub rows for {} relationships.'.format(len(rel_subject_query_results), entity))
                all_query_results.extend(rel_subject_query_results)
                formatted_entity_query = generic_rel_object_query.format(entity=entity, entity_regex=entity_regex,
                                                                         excluded_entity_regex=excluded_entity_regex,
                                                                         pub_regex=pub_regex)
                rel_object_query_results = connect(formatted_entity_query, 'no_query', conn)
                log.info('Found {} object-pub rows for {} relationships.'.format(len(rel_object_query_results), entity))
                all_query_results.extend(rel_object_query_results)
        else:
            log.info('Querying {} entity.'.format(entity))
            formatted_entity_query = entity_query
            entity_query_results = connect(formatted_entity_query, 'no_query', conn)
            log.info('Found {} entity-pub rows for {} query.'.format(len(entity_query_results), entity))
            all_query_results.extend(entity_query_results)

    return all_query_results


def process_db_results(db_results):
    """Process db query tuples into header-keyed dicts for printing to bulk report.

    Args:
        arg1 (db_results): A list of (entity_id, pub_fbrf_id, entity_name) tuples.

    Returns:
        A list of dicts having keys that match the global "header_list" elements.

    """
    log.info('Processing db results into exportable dicts.')

    pub_dict = make_pub_dict()

    log.info('Sorting {} db results.'.format(len(db_results)))
    db_results = sorted(list(set(db_results)))

    log.info('Converting {} unique db result tuples to dicts.'.format(len(db_results)))
    formatted_result_list = []
    ENTITY_UNIQUENAME = 0
    FBRF_ID = 1
    ENTITY_NAME = 2
    for row in db_results:
        result_dict = {
            'entity_id': row[ENTITY_UNIQUENAME],
            'entity_name': row[ENTITY_NAME],
            'FlyBase_publication_id': row[FBRF_ID],
            'PubMed_id': pub_dict[row[FBRF_ID]]
        }
        formatted_result_list.append(result_dict)
    log.info('The db results have been converted to dicts.')

    return formatted_result_list


if __name__ == "__main__":
    main()
