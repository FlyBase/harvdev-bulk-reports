# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate a bulk report with the best available gene summary for each gene.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_best_gene_summary.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_best_gene_summary.py -v -c /foo/bar/config.cfg

"""

from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
from harvdev_utils.char_conversions import (
    sgml_to_plain_text, sub_sup_sgml_to_plain_text
)
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import set_up_db_reading
from harvdev_utils.reporting import (
    Cvterm, Db, Dbxref, Dbxrefprop, Feature, FeatureDbxref, Featureprop, Organism
)


# Important label for output files.
report_label = 'best_gene_summary'
report_title = 'FlyBase Best Gene Summary report'
header_list = [
    'FBgn_ID',
    'Gene_Symbol',
    'Summary_Source',
    'Summary'
]

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
input_dir = set_up_dict['input_dir']
output_dir = set_up_dict['output_dir']
output_filename = set_up_dict['output_filename']
log = set_up_dict['log']

# Create SQL Alchemy engines from environmental variables.
log.info('Open db connection to {} {}'.format(server, database))
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
insp = inspect(engine)


# The main process.
def main():    # noqa E501
    """Get all gene summaries and report the best one."""
    log.info('Started main function.')

    # Create the summary handler instance, get db info, pick the best summary.
    summary_handler = SummaryHandler()
    db_query_transaction(summary_handler)
    summary_handler.pick_best_summary()

    # Print out the best info.
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    notes = []
    notes.append('The single best available gene summary is reported for each D. melanogaster gene.')
    notes.append('Gene summaries are taken from the following sources, in order of decreasing rank:')
    notes.append('FlyBase gene snapshots, UniProtKB functional descriptions, InteractiveFly summaries, Alliance of Genome Resources automated descriptions, FlyBase automatically generated summaries.')
    notes.append('For other non-D. melanogaster genes, please see FlyBase\'s "automated_gene_summaries.tsv.gz" file.')
    data_to_export_as_tsv['metaData']['note'] = notes
    data_to_export_as_tsv['data'] = summary_handler.export_data()
    tsv_report_dump(data_to_export_as_tsv, output_filename, print_footer=False, headers=header_list)

    log.info('Ended main function.\n')


class SummaryHandler(object):
    """Object that gets gene summaries and picks the best one to write as a new Featureprop."""
    def __init__(self):
        """Create the object that picks the best gene summary."""
        # Create a FBgn_ID-keyed dict of summary dicts.
        # Each summary dict has gene ID, name, lists of summaries (by type), and keys for the best summary.
        self.summary_dict = {}

    # List gene summary types, in order of preference.
    rank_ordered_summary_types = [
        'FlyBase Gene Snapshot',
        'UniProtKB',
        'Interactive Fly',
        'Alliance',
        'FlyBase Auto Summary'
    ]

    def get_dmel_genes(self, session):
        """Get all Dmel genes."""
        log.info('Get Dmel genes.')
        fbgn_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(fbgn_regex),
            Organism.abbreviation == 'Dmel'
        )
        results = session.query(Feature).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            self.summary_dict[result.uniquename] = {
                'uniquename': result.uniquename,
                'name': result.name,
                'best_summary_type': None,
                'best_summary_text': None
            }
        log.info('Found {} Dmel genes.'.format(counter))
        return

    def get_gene_snapshots(self, session):
        """Get FlyBase Gene Snapshot summaries."""
        summary_type = 'FlyBase Gene Snapshot'
        log.info('Get "{}" gene summaries.'.format(summary_type))
        fbgn_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(fbgn_regex),
            Organism.abbreviation == 'Dmel',
            Cvterm.name == 'gene_summary_text'
        )
        results = session.query(Feature, Featureprop).\
            join(Feature, (Feature.feature_id == Featureprop.feature_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            summary_string = result.Featureprop.value.replace('@', '')
            summary_string = sgml_to_plain_text(summary_string)
            summary_string = sub_sup_sgml_to_plain_text(summary_string)
            try:
                self.summary_dict[result.Feature.uniquename][summary_type].append(summary_string)
                log.warning('Have multiple {} summaries for {}'.format(summary_type, result.Feature.uniquename))
            except KeyError:
                self.summary_dict[result.Feature.uniquename][summary_type] = [summary_string]
        log.info('Found {} "{}" summaries.'.format(counter, summary_type))
        return

    def get_uniprot_function_comments(self, session):
        """Get UniProt/GCRP "UniProt_Function_comment" summaries for genes."""
        summary_type = 'UniProtKB'
        log.info('Get "{}" gene summaries.'.format(summary_type))
        fbgn_regex = r'^FBgn[0-9]{7}$'
        gcrp = aliased(Db, name='gcrp')
        gcrp_xref = aliased(Dbxref, name='gcrp_xref')
        swissprot = aliased(Db, name='swissprot')
        swissprot_xref = aliased(Dbxref, name='swissprot_xref')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(fbgn_regex),
            Organism.abbreviation == 'Dmel',
            FeatureDbxref.is_current.is_(True),
            swissprot.name == 'UniProt/Swiss-Prot',
            gcrp.name == 'UniProt/GCRP',
            Cvterm.name == 'UniProt_Function_comment'
        )
        results = session.query(Feature, gcrp_xref, Dbxrefprop).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(gcrp_xref, (gcrp_xref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(gcrp, (gcrp.db_id == gcrp_xref.db_id)).\
            join(swissprot_xref, (swissprot_xref.accession == gcrp_xref.accession)).\
            join(swissprot, (swissprot.db_id == swissprot_xref.db_id)).\
            join(Dbxrefprop, (Dbxrefprop.dbxref_id == swissprot_xref.dbxref_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Dbxrefprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            cleaned_prop_value = result.Dbxrefprop.value.replace('\n', ' ').replace('  ', ' ')
            summary_string = '{} (UniProtKB:{})'.format(cleaned_prop_value, result.gcrp_xref.accession)
            try:
                self.summary_dict[result.Feature.uniquename][summary_type].append(summary_string)
                log.warning('Have multiple {} summaries for {}'.format(summary_type, result.Feature.uniquename))
            except KeyError:
                self.summary_dict[result.Feature.uniquename][summary_type] = [summary_string]
        log.info('Found {} "{}" summaries.'.format(counter, summary_type))
        return

    def get_interactive_fly_summaries(self, session):
        """Get Interactive Fly summaries."""
        summary_type = 'Interactive Fly'
        log.info('Get "{}" gene summaries.'.format(summary_type))
        fbgn_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(fbgn_regex),
            Organism.abbreviation == 'Dmel',
            FeatureDbxref.is_current.is_(True),
            Db.name == 'INTERACTIVEFLY',
            Cvterm.name == 'if_summary'
        )
        results = session.query(Feature, Dbxrefprop).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            join(Dbxrefprop, (Dbxrefprop.dbxref_id == Dbxref.dbxref_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Dbxrefprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            try:
                self.summary_dict[result.Feature.uniquename][summary_type].append(result.Dbxrefprop.value)
                log.warning('Have multiple {} summaries for {}'.format(summary_type, result.Feature.uniquename))
            except KeyError:
                self.summary_dict[result.Feature.uniquename][summary_type] = [result.Dbxrefprop.value]
        log.info('Found {} "{}" summaries.'.format(counter, summary_type))
        return

    def get_alliance_summaries(self):
        """Get Alliance gene summaries from download file."""
        summary_type = 'Alliance'
        log.info('Get "{}" gene summaries.'.format(summary_type))
        input_filename = '{}alliance_gene_descriptions_fb_{}.tsv'.format(input_dir, database_release)
        counter = 0
        try:
            input_file = open(input_filename, 'r')
        except (TypeError, FileNotFoundError):
            log.error('Could not get Alliance summary file: {}'.format(input_filename))
            raise Exception
        FBGN_ID = 0
        GENE_SYMBOL = 1
        SUMMARY_TEXT = 2
        for line in input_file:
            if line.startswith('#'):
                continue
            if line == '\n':
                continue
            line_parts = line.split('\t')
            fbgn_id = line_parts[FBGN_ID].replace('FB:', '')
            symbol = line_parts[GENE_SYMBOL]
            summary_text = line_parts[SUMMARY_TEXT].strip()
            if summary_text.startswith('No description available'):
                continue
            if fbgn_id not in self.summary_dict.keys():
                log.warning('The gene "{}" ({}) does not correspond to a current FBgn ID.'.format(symbol, fbgn_id))
                continue
            try:
                self.summary_dict[fbgn_id][summary_type].append(summary_text)
                log.warning('Have multiple {} summaries for {}'.format(summary_type, fbgn_id))
            except KeyError:
                self.summary_dict[fbgn_id][summary_type] = [summary_text]
            counter += 1
        log.info('Found {} "{}" summaries.'.format(counter, summary_type))
        return

    def query_chado(self, session):
        """Query for each type of gene summary."""
        self.get_dmel_genes(session)
        self.get_gene_snapshots(session)
        self.get_uniprot_function_comments(session)
        self.get_interactive_fly_summaries(session)
        self.get_alliance_summaries()
        return

    def pick_best_summary(self):
        """Pick the best summary for each gene."""
        log.info('Pick the best gene summary for each gene.')
        for gene in self.summary_dict.values():
            log.debug('Assessing gene "{}" ({})'.format(gene['name'], gene['uniquename']))
            for summary_type in self.rank_ordered_summary_types:
                if summary_type in gene.keys():
                    gene['best_summary_type'] = summary_type
                    break
            if gene['best_summary_type']:
                gene['best_summary_text'] = ' | '.join(gene[gene['best_summary_type']])
        return

    def export_data(self):
        """Export the data as simple dicts for TSV output."""
        log.info('Exporting best summaries to output file.')
        input_counter = 0
        export_counter = 0
        export_data_list = []
        for gene in self.summary_dict.values():
            input_counter += 1
            if not gene['best_summary_type']:
                continue
            simple_summary_dict = {
                'FBgn_ID': gene['uniquename'],
                'Gene_Symbol': gene['name'],
                'Summary_Source': gene['best_summary_type'],
                'Summary': gene['best_summary_text']
            }
            export_data_list.append(simple_summary_dict)
            export_counter += 1
        log.info('Exported summaries for {} of {} Dmel genes.'.format(export_counter, input_counter))
        return export_data_list


def db_query_transaction(object_to_execute):
    """Query the chado database given an object that has a "query_chado()" method.

    Function assumes a global "engine" variable for SQLAlchemy processes.

    Args:
        arg1 (object_to_execute): (object) Some object that has a "query_chado()" method.

    Returns:
        None.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        object_to_execute.query_chado(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during chado query; rolling back and exiting.')
        raise
    return


if __name__ == "__main__":
    main()
