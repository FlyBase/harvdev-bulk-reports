# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report experimental tools info.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_experimental_tools.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_experimental_tools.py -v -c /path/to/config.cfg

"""

import argparse
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
report_label = 'experimental_tool_data'
report_title = 'FlyBase experimental tool report'
header_list = [
    'Symbol',
    'FlyBase ID',
    'Name',
    'Uses (term)',
    'Uses (id)',
    'Description',
    'Compatible tools (symbol)',
    'Compatible tools (id)',
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
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))


def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    database_info = get_database_info()
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = process_database_info(database_info)
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    conn.close()
    log.info('Ended main function.')


def get_exp_tools():
    """Query chado for all experimental tools.

    Returns:
        tool_dict (dict): A feature_id-keyed dict of tools (represented by dicts).

    """
    global conn
    log.info('Query chado for all experimental tools.')
    fb_tools_query = """
        SELECT DISTINCT f.feature_id, f.name, f.uniquename
        FROM feature f
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBto[0-9]{7}$';
    """
    ret_fb_tools = connect(fb_tools_query, 'no_query', conn)
    FEAT_ID = 0
    SYMBOL = 1
    FB_CURIE = 2
    fb_tool_dict = {}
    for result in ret_fb_tools:
        fb_tool_dict[result[FEAT_ID]] = {
            'Symbol': result[SYMBOL],
            'FlyBase ID': result[FB_CURIE],
            'Name': None,
            'Uses (term)': [],
            'Uses (id)': [],
            'Description': [],
            'compatible_tool_tuples': [],
            'Compatible tools (symbol)': [],
            'Compatible tools (id)': [],
        }
    log.info(f'Found {len(ret_fb_tools)} current experimental tools in chado.')
    return fb_tool_dict


def get_tool_fullnames(fb_tool_dict):
    """Get tool full names."""
    global conn
    log.info('Get tool full names.')
    fb_tool_fullnames_query = """
        SELECT DISTINCT f.feature_id, s.name
        FROM feature f
        JOIN feature_synonym fs ON fs.feature_id = f.feature_id
        JOIN synonym s ON s.synonym_id = fs.synonym_id
        JOIN cvterm t ON t.cvterm_id = s.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBto[0-9]{7}$'
          AND fs.is_current IS TRUE
          AND t.name = 'fullname';
    """
    ret_fb_tool_fullnames = connect(fb_tool_fullnames_query, 'no_query', conn)
    FEAT_ID = 0
    FULLNAME = 1
    counter = 0
    for result in ret_fb_tool_fullnames:
        fb_tool_dict[result[FEAT_ID]]['Name'] = result[FULLNAME]
        counter += 1
    log.info(f'Found {counter} current full names for experimental tools in chado.')
    return


def get_tool_uses(fb_tool_dict):
    """Get tool uses."""
    global conn
    log.info('Get tool uses.')
    fb_tool_uses_query = """
        SELECT DISTINCT f.feature_id, cvt.name, db.name||':'||dbx.accession
        FROM feature f
        JOIN feature_cvterm fcvt ON fcvt.feature_id = f.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
        JOIN feature_cvtermprop fcvtp ON fcvtp.feature_cvterm_id = fcvt.feature_cvterm_id
        JOIN cvterm t ON t.cvterm_id = fcvtp.type_id
        JOIN dbxref dbx ON dbx.dbxref_id = cvt.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBto[0-9]{7}$'
          AND cvt.is_obsolete = 0
          AND t.name = 'tool_uses'
    """
    ret_fb_tool_uses = connect(fb_tool_uses_query, 'no_query', conn)
    FEAT_ID = 0
    TERM_NAME = 1
    TERM_CURIE = 2
    counter = 0
    for result in ret_fb_tool_uses:
        fb_tool_dict[result[FEAT_ID]]['Uses (term)'].append(result[TERM_NAME])
        fb_tool_dict[result[FEAT_ID]]['Uses (id)'].append(result[TERM_CURIE])
        counter += 1
    log.info(f'Found {counter} "tool_uses" annotations for experimental tools in chado.')
    return


def get_tool_descriptions(fb_tool_dict):
    """Get tool descriptions."""
    global conn
    log.info('Get tool descriptions.')
    fb_tool_desc_query = """
        SELECT DISTINCT f.feature_id, fp.value
        FROM feature f
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm t ON t.cvterm_id = fp.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.uniquename ~ '^FBto[0-9]{7}$'
          AND t.name = 'description';
    """
    ret_fb_tool_desc = connect(fb_tool_desc_query, 'no_query', conn)
    FEAT_ID = 0
    DESCRIPTION = 1
    counter = 0
    for result in ret_fb_tool_desc:
        fb_tool_dict[result[FEAT_ID]]['Description'].append(result[DESCRIPTION])
        counter += 1
    log.info(f'Found {counter} descriptions for experimental tools in chado.')
    return


def get_compatible_tools(fb_tool_dict):
    """Get compatible tools."""
    global conn
    log.info('Get compatible tools.')
    fb_compatible_tools_query = """
        SELECT DISTINCT s.feature_id, s.name, s.uniquename, o.feature_id, o.name, o.uniquename
        FROM feature s
        JOIN feature_relationship fr ON fr.subject_id = s.feature_id
        JOIN feature o ON o.feature_id = fr.object_id
        JOIN cvterm t ON t.cvterm_id = fr.type_id
        WHERE s.is_obsolete IS FALSE
          AND s.uniquename ~ '^FBto[0-9]{7}$'
          AND o.is_obsolete IS FALSE
          AND o.uniquename ~ '^FBto[0-9]{7}$'
          AND t.name = 'compatible_tool';
    """
    ret_fb_compatible_tools = connect(fb_compatible_tools_query, 'no_query', conn)
    SBJ_ID = 0
    SBJ_NAME = 1
    SBJ_CURIE = 2
    OBJ_ID = 3
    OBJ_NAME = 4
    OBJ_CURIE = 5
    counter = 0
    for result in ret_fb_compatible_tools:
        sbj_tuple = (result[SBJ_NAME], result[SBJ_CURIE])
        obj_tuple = (result[OBJ_NAME], result[OBJ_CURIE])
        fb_tool_dict[result[SBJ_ID]]['compatible_tool_tuples'].append(obj_tuple)
        fb_tool_dict[result[OBJ_ID]]['compatible_tool_tuples'].append(sbj_tuple)
        counter += 1
    for fb_tool in fb_tool_dict.values():
        fb_tool['compatible_tool_tuples'] = set(fb_tool['compatible_tool_tuples'])
        for compatible_tool_tuple in fb_tool['compatible_tool_tuples']:
            fb_tool['Compatible tools (symbol)'].append(compatible_tool_tuple[0])
            fb_tool['Compatible tools (id)'].append(compatible_tool_tuple[1])
    log.info(f'Found {counter} "compatible_tool" relationships in chado.')
    return


def get_database_info():
    """Run chado db queries in sequence."""
    global conn
    log.info('Query database.')
    tool_dict = get_exp_tools()
    get_tool_fullnames(tool_dict)
    get_tool_uses(tool_dict)
    get_tool_descriptions(tool_dict)
    get_compatible_tools(tool_dict)
    return tool_dict.values()


def process_database_info(input_data):
    """Process a list of data elements for TSV output.

    Args:
        input_data (list): A list of dicts representing FB experimental tools.

    Returns:
        A list of dictionaries representing, in this case, paralog information.
    """
    log.info('Starting to process paralog info retrieved from database.')
    data_list = []
    for i in input_data:
        i['Uses (term)'] = '|'.join(i['Uses (term)'])
        i['Uses (id)'] = '|'.join(i['Uses (id)'])
        i['Description'] = '|'.join(i['Description'])
        i['Compatible tools (symbol)'] = '|'.join(i['Compatible tools (symbol)'])
        i['Compatible tools (id)'] = '|'.join(i['Compatible tools (id)'])
    data_list.append(i)
    log.info('Done processing paralog info into a list of dictionaries.')
    return data_list


if __name__ == "__main__":
    main()
