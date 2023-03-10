# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report scRNA-seq data.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_scrna_seq_data.py [-h] [-c CONFIG] [-v VERBOSE]

Example:
    python report_scrna_seq_data.py -v -t -c /foo/bar/config.cfg

Notes:
    This script reports mean expression and spread for scRNA-seq clusters.

"""

import argparse
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
# from harvdev_utils.chado_functions import get_or_create
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, ExpressionCvterm, Feature, Library, LibraryCvterm,
    LibraryCvtermprop, LibraryExpression, LibraryFeature, LibraryFeatureprop,
    LibraryPub, LibraryRelationship, Pub
)
from harvdev_utils.psycopg_functions import set_up_db_reading


# Global variables for the output file. Header order will match list order below.
report_label = 'scRNA-seq'
report_title = 'FlyBase scRNA-seq gene expression'
header_list = [
    'Pub_ID',
    'Pub_miniref',
    'Clustering_Analysis_ID',
    'Clustering_Analysis_Name',
    'Source_Tissue_Sex',
    'Source_Tissue_Stage',
    'Source_Tissue_Anatomy',
    'Cluster_ID',
    'Cluster_Name',
    'Cluster_Cell_Type_ID',
    'Cluster_Cell_Type_Name',
    'Gene_ID',
    'Gene_Symbol',
    'Mean_Expression',
    'Spread'
]

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
output_dir = set_up_dict['output_dir']
output_filename = set_up_dict['output_filename']
log = set_up_dict['log']

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
insp = inspect(engine)

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')


# The main process.
def main():
    """Run the steps for updating gene associations to sequence targeting reagents."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')

    # Instantiate the handler and run its "write_chado()" method.
    data_reporter = SingleCellRNASeqReporter()
    db_query_transaction(data_reporter)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    notes = []
    notes.append('Note: Mean_Expression is the average level of expression of the gene across all cells of the cluster in which the gene is detected at all.')
    notes.append('Note: Spread is the proportion of cells in the cluster in which the gene is detected.')
    notes.append('Note: In "Source_Tissue_*" columns, "mixed" is shown when there are many applicable terms - please see the dataset report for details.')
    data_to_export_as_tsv['metaData']['note'] = notes
    data_to_export_as_tsv['data'] = data_reporter.data_to_export
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)

    # Close
    log.info('Ended main function.\n')


class ClusteringAnalysis(object):
    """A clustering analysis along with its child clusters."""

    def __init__(self, library):
        """Create a ClusteringAnalysis object from a Library object.

        Args:
            feature: (Library) An SQLAlchemy object for the library.

        Returns:
            An object of the ClusteringAnalysis class.

        """
        # Data from chado.
        self.library = library                 # The Library object.
        self.id = library.uniquename           # The FB ID.
        self.child_clusters = []               # Will be list of child Library objects.
        self.papers = []                       # Will be a list of Pub objects for the clustering analysis.
        self.source_tissue_stage = []          # Will be source tissue stage terms.
        self.source_tissue_sex = []            # Will be source tissue sex terms.
        self.source_tissue_anatomy = []        # Will be source tissue anatomy terms.
        self.source_tissue_stage_str = ''      # Will be source tissue stage terms.
        self.source_tissue_sex_str = ''        # Will be source tissue sex terms.
        self.source_tissue_anatomy_str = ''    # Will be source tissue anatomy terms.

        # Bins for note, warning, error and action messages.
        self.warnings = []                    # Issues that may complicate analysis of the seqfeat.
        self.notes = []                       # Other notes about the seqfeat.

    def __str__(self):
        """Informative string for this sequence targeting reagent."""
        return f'{self.library.name} ({self.library.uniquename})'


class SingleCellRNASeqReporter(object):
    """An object that gets scRNA-seq data and exports it to file."""

    def __init__(self):
        """Create the SingleCellRNASeqReporter object."""

    # Data dicts used for data processing.
    cluster_dict = {}              # library_id-keyed dict of ClusteringAnalysis objects.
    cluster_cell_type_dict = {}    # library_id-keyed dict of cluster cell type (Cvterm object).
    mean_expr_spread_dict = {}     # library_id-keyed dict of list of (Feature, mean_expr, spread) tuples.
    data_to_export = []
    # Uniquename regexes.
    lib_regex = r'^FBlc[0-9]{7}$'
    gene_regex = r'^FBgn[0-9]{7}$'
    pub_regex = r'^FBrf[0-9]{7}$'

    def get_clustering_analyses(self, session):
        """Get datasets for clusters and parent clustering analyses."""
        log.info('Get datasets for clusters and parent clustering analyses.')
        analysis = aliased(Library, name='analysis')
        cluster = aliased(Library, name='cluster')
        lib_type = aliased(Cvterm, name='lib_type')
        lib_rel_type = aliased(Cvterm, name='lib_rel_type')
        lcvt1 = aliased(LibraryCvterm, name='lcvt1')
        lcvt2 = aliased(LibraryCvterm, name='lcvt2')
        analysis_type = aliased(Cvterm, name='analysis_type')
        cluster_type = aliased(Cvterm, name='cluster_type')
        filters = (
            analysis.is_obsolete.is_(False),
            analysis.uniquename.op('~')(self.lib_regex),
            lib_type.name == 'result',
            analysis_type.name == 'cell clustering analysis',
            cluster.is_obsolete.is_(False),
            cluster.uniquename.op('~')(self.lib_regex),
            cluster.type_id == analysis.type_id,
            cluster_type.name == 'transcriptional cell cluster',
            analysis.uniquename == 'FBlc0003731',    # BOB: DEV
            lib_rel_type.name == 'belongs_to'
        )
        results = session.query(analysis, cluster).\
            join(lib_type, (lib_type.cvterm_id == analysis.type_id)).\
            join(lcvt1, (lcvt1.library_id == analysis.library_id)).\
            join(analysis_type, (analysis_type.cvterm_id == lcvt1.cvterm_id)).\
            join(LibraryRelationship, (LibraryRelationship.object_id == analysis.library_id)).\
            join(lib_rel_type, (lib_rel_type.cvterm_id == LibraryRelationship.type_id)).\
            join(cluster, (cluster.library_id == LibraryRelationship.subject_id)).\
            join(lcvt2, (lcvt2.library_id == cluster.library_id)).\
            join(cluster_type, (cluster_type.cvterm_id == lcvt2.cvterm_id)).\
            filter(*filters).\
            distinct()
        parent_counter = 0
        cluster_counter = 0
        for result in results:
            try:
                self.cluster_dict[result.analysis.library_id].child_clusters.append(result.cluster)
                cluster_counter += 1
            except KeyError:
                self.cluster_dict[result.analysis.library_id] = ClusteringAnalysis(result.analysis)
                self.cluster_dict[result.analysis.library_id].child_clusters.append(result.cluster)
                parent_counter += 1
                cluster_counter += 1
        log.info(f'Found {parent_counter} groups of {cluster_counter} clusters.')
        return

    def get_cluster_pubs(self, session):
        """Get publication for scRNA-seq clustering analysis."""
        log.info('Get publication for scRNA-seq clustering analysis.')
        filters = (
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.lib_regex),
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(self.pub_regex),
            Cvterm.name == 'paper'
        )
        results = session.query(Library, Pub).\
            join(LibraryPub, (LibraryPub.library_id == Library.library_id)).\
            join(Pub, (Pub.pub_id == LibraryPub.pub_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Pub.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            try:
                self.cluster_dict[result.Library.library_id].papers.append(result.Pub)
                counter += 1
            except KeyError:
                pass
        log.info(f'Found {counter} papers for clustering analyses.')
        return

    def get_source_tissue_sex_and_stage(self, session):
        """Get source tissue stage for clustering analyses."""
        log.info('Get source tissue stage for clustering analyses.')
        lcvt1 = aliased(LibraryCvterm, name='lcvt1')
        lcvt2 = aliased(LibraryCvterm, name='lcvt2')
        lib_type = aliased(Cvterm, name='lib_type')
        stage_term = aliased(Cvterm, name='stage_term')
        term_association_type = aliased(Cvterm, name='term_association_type')
        filters = (
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.lib_regex),
            lib_type.name == 'cell clustering analysis',
            term_association_type.name == 'derived_stage'
        )
        results = session.query(Library, Db, stage_term).\
            join(lcvt1, (lcvt1.library_id == Library.library_id)).\
            join(lib_type, (lib_type.cvterm_id == lcvt1.cvterm_id)).\
            join(lcvt2, (lcvt2.library_id == Library.library_id)).\
            join(stage_term, (stage_term.cvterm_id == lcvt2.cvterm_id)).\
            join(LibraryCvtermprop, (LibraryCvtermprop.library_cvterm_id == lcvt2.library_cvterm_id)).\
            join(term_association_type, (term_association_type.cvterm_id == LibraryCvtermprop.type_id)).\
            join(Dbxref, (Dbxref.dbxref_id == stage_term.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        stage_counter = 0
        sex_counter = 0
        for result in results:
            if result.Db.name == 'FBdv':
                try:
                    self.cluster_dict[result.Library.library_id].source_tissue_stage.append(result.stage_term.name)
                    stage_counter += 1
                except KeyError:
                    log.warning('Failure to add FBdv term.')
                    pass
            elif result.Db.name == 'FBcv':
                try:
                    self.cluster_dict[result.Library.library_id].source_tissue_sex.append(result.stage_term.name)
                    sex_counter += 1
                except KeyError:
                    log.warning('Failure to add FBcv term.')
                    pass
            else:
                log.warning(f'Cannot handle term={result.stage_term.name}, db={result.Db.name}')

        log.info(f'Found source tissue stage info for {stage_counter} clustering analyses.')
        log.info(f'Found source tissue sex info for {sex_counter} clustering analyses.')
        return

    def get_source_tissue_anatomy(self, session):
        """Get source tissue anatomy for clustering analyses."""
        log.info('Get source tissue anatomy for clustering analyses.')
        lib_type = aliased(Cvterm, name='lib_type')
        ec_type = aliased(Cvterm, name='ec_type')
        anatomy = aliased(Cvterm, name='anatomy')
        filters = (
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.lib_regex),
            lib_type.name == 'cell clustering analysis',
            ec_type.name == 'anatomy'
        )
        results = session.query(Library, anatomy).\
            join(LibraryCvterm, (LibraryCvterm.library_id == Library.library_id)).\
            join(lib_type, (lib_type.cvterm_id == LibraryCvterm.cvterm_id)).\
            join(LibraryExpression, (LibraryExpression.library_id == Library.library_id)).\
            join(ExpressionCvterm, (ExpressionCvterm.expression_id == LibraryExpression.expression_id)).\
            join(ec_type, (ec_type.cvterm_id == ExpressionCvterm.cvterm_type_id)).\
            join(anatomy, (anatomy.cvterm_id == ExpressionCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            try:
                self.cluster_dict[result.Library.library_id].source_tissue_anatomy.append(result.anatomy.name)
            except KeyError:
                pass
        log.info(f'Found source tissue anatomy info for {counter} clustering analyses.')
        return

    def process_source_tissue_info(self, session):
        """Process source tissue info."""
        log.info('Process source tissue info.')
        for cluster_analysis in self.cluster_dict.values():
            # Stage
            if len(cluster_analysis.source_tissue_stage) == 1:
                cluster_analysis.source_tissue_stage_str = cluster_analysis.source_tissue_stage[0]
            elif len(cluster_analysis.source_tissue_stage) > 1:
                cluster_analysis.source_tissue_stage_str = 'mixed'
            # Anatomy
            if len(cluster_analysis.source_tissue_anatomy) == 1:
                cluster_analysis.source_tissue_anatomy_str = cluster_analysis.source_tissue_anatomy[0]
            elif len(cluster_analysis.source_tissue_anatomy) > 1:
                cluster_analysis.source_tissue_anatomy_str = 'mixed'
            # Sex
            male = False
            female = False
            for sex_term in cluster_analysis.source_tissue_sex:
                if 'female' in sex_term:
                    female = True
                elif 'male' in sex_term:
                    male = True
            if male and female:
                cluster_analysis.source_tissue_sex_str = 'mixed'
            elif female:
                cluster_analysis.source_tissue_sex_str = 'female'
            elif male:
                cluster_analysis.source_tissue_sex_str = 'male'
        return

    def get_cluster_cell_types(self, session):
        """Get cell types for expression clusters."""
        log.info('Get cell types for expression clusters.')
        cvterm_type = aliased(Cvterm, name='cvterm_type')
        cell_type = aliased(Cvterm, name='cell_type')
        lib_type = aliased(Cvterm, name='lib_type')
        filters = (
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.lib_regex),
            lib_type.name == 'transcriptional cell cluster',
            cvterm_type.name == 'anatomy',
            Db.name == 'FBbt'
        )
        counter = 0
        results = session.query(Library, cell_type).\
            join(LibraryCvterm, (LibraryCvterm.library_id == Library.library_id)).\
            join(lib_type, (lib_type.cvterm_id == LibraryCvterm.cvterm_id)).\
            join(LibraryExpression, (LibraryExpression.library_id == Library.library_id)).\
            join(ExpressionCvterm, (ExpressionCvterm.expression_id == LibraryExpression.expression_id)).\
            join(cvterm_type, (cvterm_type.cvterm_id == ExpressionCvterm.cvterm_type_id)).\
            join(cell_type, (cell_type.cvterm_id == ExpressionCvterm.cvterm_id)).\
            join(Dbxref, (Dbxref.dbxref_id == cell_type.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        for result in results:
            self.cluster_cell_type_dict[result.Library.library_id] = result.cell_type
            counter += 1
        log.info(f'Found {counter} cluster-to-cell_type term associations.')
        return

    def get_mean_expr_spread_values(self, session):
        """Get mean_expr and spread values for scRNA-seq data."""
        log.info('Get mean_expr and spread values for scRNA-seq data.')
        mean_expr = aliased(LibraryFeatureprop, name='mean_expr')
        spread = aliased(LibraryFeatureprop, name='spread')
        mean_expr_type = aliased(Cvterm, name='mean_expr_type')
        spread_type = aliased(Cvterm, name='spread_type')
        lib_counter = 0
        data_counter = 0
        for analysis in self.cluster_dict.values():
            for cluster in analysis.child_clusters:
                lib_counter += 1
                if lib_counter % 100 == 0:
                    log.info(f'Getting data for cluster #{lib_counter}.')
                log.debug(f'Getting mean_expr and spread data for {cluster.name} ({cluster.uniquename}).')
                self.mean_expr_spread_dict[cluster.library_id] = []
                filters = (
                    LibraryFeature.library_id == cluster.library_id,
                    Feature.is_obsolete.is_(False),
                    Feature.uniquename.op('~')(self.gene_regex),
                    mean_expr_type.name == 'mean_expr',
                    spread_type.name == 'spread'
                )
                results = session.query(Feature, mean_expr, spread).\
                    join(LibraryFeature, (LibraryFeature.feature_id == Feature.feature_id)).\
                    join(mean_expr, (mean_expr.library_feature_id == LibraryFeature.library_feature_id)).\
                    join(mean_expr_type, (mean_expr_type.cvterm_id == mean_expr.type_id)).\
                    join(spread, (spread.library_feature_id == LibraryFeature.library_feature_id)).\
                    join(spread_type, (spread_type.cvterm_id == spread.type_id)).\
                    filter(*filters).\
                    distinct()
                for result in results:
                    datum = {
                        'id': result.Feature.uniquename,
                        'name': result.Feature.name,
                        'mean_expr': result.mean_expr.value,
                        'spread': result.spread.value
                    }
                    self.mean_expr_spread_dict[cluster.library_id].append(datum)
                    data_counter += 1
        log.info(f'Found {data_counter} scRNA-seq "spread" data points.')
        return

    def process_database_info(self):
        """Print out scRNA-seq data."""
        log.info('Print out scRNA-seq data.')
        for analysis in self.cluster_dict.values():
            for cluster in analysis.child_clusters:
                data_key = cluster.library_id
                log.debug(f'Export data for {cluster.name} ({cluster.uniquename})')
                for datum in self.mean_expr_spread_dict[data_key]:
                    data_dict = {
                        'Pub_ID': f'{analysis.papers[0].uniquename}',
                        'Pub_miniref': f'{analysis.papers[0].miniref}',
                        'Clustering_Analysis_ID': f'{analysis.library.uniquename}',
                        'Clustering_Analysis_Name': f'{analysis.library.name}',
                        'Source_Tissue_Sex': f'{analysis.source_tissue_sex_str}',
                        'Source_Tissue_Stage': f'{analysis.source_tissue_stage_str}',
                        'Source_Tissue_Anatomy': f'{analysis.source_tissue_anatomy_str}',
                        'Cluster_ID': f'{cluster.uniquename}',
                        'Cluster_Name': f'{cluster.name}',
                        'Cluster_Cell_Type_ID': f'FBbt:{self.cluster_cell_type_dict[cluster.library_id].dbxref.accession}',
                        'Cluster_Cell_Type_Name': f'{self.cluster_cell_type_dict[cluster.library_id].name}',
                        'Gene_ID': f'{datum["id"]}',
                        'Gene_Symbol': f'{datum["name"]}',
                        'Mean_Expression': f'{datum["mean_expr"]}',
                        'Spread': f'{datum["spread"]}'
                    }
                    self.data_to_export.append(data_dict)
        return

    def query_chado(self, session):
        """Run write methods."""
        log.info('Starting "write_to_chado" method.')
        self.get_clustering_analyses(session)
        self.get_cluster_pubs(session)
        self.get_source_tissue_sex_and_stage(session)
        self.get_source_tissue_anatomy(session)
        self.process_source_tissue_info(session)
        self.get_cluster_cell_types(session)
        self.get_mean_expr_spread_values(session)
        self.process_database_info()
        log.info('Method "write_to_chado" is done.')
        return


def db_query_transaction(object_to_execute):
    """Query the chado database given an object that has a "query_chado()" method.

    Function assumes a global "engine" variable for SQLAlchemy processes.

    Args:
        arg1 (class): some object that has a "query_chado()" method.

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
