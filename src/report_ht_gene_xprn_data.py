# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report high-throughput gene expression data.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_ht_gene_xprn_data.py [-h] [-c CONFIG] [-v VERBOSE]

Example:
    python report_ht_gene_xprn_data.py -v -c /foo/bar/config.cfg

Notes:
    This script reports gene expression values (with units) from various high-
    -throughput expression datasets featured in gene report bar graphs. It is
    limited to cases where there is a single value type per sample, which
    excludes scRNA-seq (for which we have a separate file).

"""

import argparse
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
# from harvdev_utils.chado_functions import get_or_create
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.production import (
    Cvterm, Feature, Library, LibraryFeature, LibraryFeatureprop, LibraryRelationship
)
from harvdev_utils.psycopg_functions import set_up_db_reading


# Global variables for the output file. Header order will match list order below.
report_label = 'high-throughput_gene_expression'
report_title = 'FlyBase high-throughput gene expression'
header_list = [
    'High_Throughput_Expression_Section',
    'Dataset_ID',
    'Dataset_Name',
    'Sample_ID',
    'Sample_Name',
    'Gene_ID',
    'Gene_Symbol',
    'Expression_Unit',
    'Expression_Value'
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
    data_reporter = HTXprnReporter()
    db_query_transaction(data_reporter)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = data_reporter.data_to_export
    notes = ['This file reports high-throughput gene expression, as reported in the "High-Throughput Expression Data" section of FlyBase gene reports.']
    notes.append('This file does not include scRNA-Seq data, which is structured differently and available in other download files.')
    data_to_export_as_tsv['metaData']['note'] = notes
    data_to_export_as_tsv['data'] = data_reporter.data_to_export
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)

    # Close
    log.info('Ended main function.\n')


class HTXprnReporter(object):
    """An object that gets high-throughput expression data and exports it to file."""

    def __init__(self):
        """Create the HTXprnReporter object."""
        # Data bins.
        self.data_to_export = []    # Will be the data, sorted by key in self.all_data_dict above.

    # Uniquename regexes.
    lib_regex = r'^FBlc[0-9]{7}$'
    gene_regex = r'^FBgn[0-9]{7}$'
    # Parental datasets to report.
    datasets_to_report = {
        'FlyAtlas2': 'FlyAtlas2 Anatomy RNA-Seq',    # FlyAtlas2 spans two graphs: RNA-seq & miRNA-seq (TPM).
        'modENCODE_mRNA-Seq_tissues': 'modENCODE Anatomy RNA-Seq',
        'modENCODE_mRNA-Seq_development': 'modENCODE Development RNA-Seq',
        'modENCODE_mRNA-Seq_cell.B': 'modENCODE Cell Lines RNA-Seq',
        'modENCODE_mRNA-Seq_treatments': 'modENCODE Treatments RNA-Seq',
        'Knoblich_Neural_Cell_RNA-Seq': 'Knoblich Neural Cell RNA-Seq',
        'Lai_miRNA_RPMM_expression_development': 'Anatomy miRNA RNA-Seq',
        'Lai_miRNA_RPMM_expression_tissues': 'Development miRNA RNA-Seq',
        'Lai_miRNA_RPMM_expression_cells': 'Cell Lines miRNA RNA-Seq',
        'Casas-Vila_proteome_life_cycle': 'Developmental Proteome: Life Cycle',
        'Casas-Vila_proteome_embryogenesis': 'Developmental Proteome: Embryogenesis',
    }
    # HT datasets lacking a parental dataset.
    samples_to_report = {
        'testis_specificity_index_2021_Vedelek': 'Testis Specificity Index',
    }
    xprn_section_order = {
        'FlyAtlas2 Anatomy RNA-Seq': 0,
        'FlyAtlas2 Anatomy miRNA RNA-Seq': 1,
        'modENCODE Anatomy RNA-Seq': 2,
        'modENCODE Development RNA-Seq': 3,
        'modENCODE Cell Lines RNA-Seq': 4,
        'modENCODE Treatments RNA-Seq': 5,
        'Knoblich Neural Cell RNA-Seq': 6,
        'Anatomy miRNA RNA-Seq': 7,
        'Development miRNA RNA-Seq': 8,
        'Cell Lines miRNA RNA-Seq': 9,
        'Developmental Proteome: Life Cycle': 10,
        'Developmental Proteome: Embryogenesis': 11,
        'Testis Specificity Index': 12,
    }
    # Data types to report.
    xprn_types_to_report = [
        'RPKM',
        'RPMM',
        'TPM',
        'LFQ_geom_mean_intensity',
        'testis_specificity_index_score',
    ]

    def get_ht_project_data(self, session):
        """Get high-throughput data for datasets having many samples."""
        log.info('Get high-throughput data for datasets having many samples.')
        dataset = aliased(Library, name='dataset')
        sample = aliased(Library, name='sample')
        gene = aliased(Feature, name='gene')
        value = aliased(LibraryFeatureprop, name='value')
        unit = aliased(Cvterm, name='unit')
        lib_rel_type = aliased(Cvterm, name='lib_rel_type')
        counter = 0
        for dataset_name, xprn_section in self.datasets_to_report.items():
            log.info(f'Get expression data for {dataset_name}.')
            xprn_section_rank = self.xprn_section_order[xprn_section]
            # Get the data.
            filters = (
                gene.is_obsolete.is_(False),
                gene.uniquename.op('~')(self.gene_regex),
                sample.is_obsolete.is_(False),
                sample.uniquename.op('~')(self.lib_regex),
                dataset.is_obsolete.is_(False),
                dataset.uniquename.op('~')(self.lib_regex),
                dataset.name == dataset_name,
                unit.name.in_((self.xprn_types_to_report)),
                lib_rel_type.name == 'belongs_to'
            )
            results = session.query(dataset, sample, gene, unit, value).\
                select_from(gene).\
                join(LibraryFeature, (LibraryFeature.feature_id == gene.feature_id)).\
                join(sample, (sample.library_id == LibraryFeature.library_id)).\
                join(value, (value.library_feature_id == LibraryFeature.library_feature_id)).\
                join(unit, (unit.cvterm_id == value.type_id)).\
                join(LibraryRelationship, (LibraryRelationship.subject_id == sample.library_id)).\
                join(lib_rel_type, (lib_rel_type.cvterm_id == LibraryRelationship.type_id)).\
                join(dataset, (dataset.library_id == LibraryRelationship.object_id)).\
                filter(*filters).\
                distinct()
            # Process the data into data dicts for export.
            this_counter = 0
            this_data_dict = {}
            for result in results:
                # Make an adjustment for FlyAtlas2 bar graphs.
                if dataset_name == 'FlyAtlas2' and result.sample.name.startswith('microRNA'):
                    xprn_section_to_use = 'FlyAtlas2 Anatomy miRNA RNA-Seq'
                    xprn_section_rank_to_use = self.xprn_section_order['FlyAtlas2 Anatomy miRNA RNA-Seq']
                else:
                    xprn_section_to_use = xprn_section
                    xprn_section_rank_to_use = xprn_section_rank
                # Make an adjustment for FlyAtlas2 FPKM data (temporarily in chado as RPKM so as to not break web).
                if dataset_name == 'FlyAtlas2' and result.unit.name == 'RPKM':
                    unit_to_use = 'FPKM'
                else:
                    unit_to_use = result.unit.name
                # Record the xprn_section, sample id and gene id as the data dict key for sorting.
                data_dict_key = (xprn_section_rank_to_use, result.sample.uniquename, result.gene.uniquename)
                # Build the dict itself.
                data_dict = {
                    'High_Throughput_Expression_Section': xprn_section_to_use,
                    'Dataset_ID': result.dataset.uniquename,
                    'Dataset_Name': result.dataset.name,
                    'Sample_ID': result.sample.uniquename,
                    'Sample_Name': result.sample.name,
                    'Gene_ID': result.gene.uniquename,
                    'Gene_Symbol': result.gene.name,
                    'Expression_Unit': unit_to_use,
                    'Expression_Value': result.value.value
                }
                this_data_dict[data_dict_key] = data_dict
                this_counter += 1
            # Sort all data before sending it to the export list.
            data_keys = list(this_data_dict.keys())
            data_keys.sort()
            for i in data_keys:
                self.data_to_export.append(this_data_dict[i])
            counter += this_counter
            log.info(f'Found {this_counter} expression values for {dataset_name}.')
        log.info(f'Found {counter} expression values for dataset projects overall.')
        return

    def get_ht_sample_data(self, session):
        """Get high-throughput data for individual samples/analyses."""
        log.info('Get high-throughput data for individual samples/analyses.')
        sample = aliased(Library, name='sample')
        gene = aliased(Feature, name='gene')
        value = aliased(LibraryFeatureprop, name='value')
        unit = aliased(Cvterm, name='unit')
        counter = 0
        for sample_name, xprn_section in self.samples_to_report.items():
            log.info(f'Get expression data for {sample_name}.')
            xprn_section_rank = self.xprn_section_order[xprn_section]
            # Get the data.
            filters = (
                gene.is_obsolete.is_(False),
                gene.uniquename.op('~')(self.gene_regex),
                sample.is_obsolete.is_(False),
                sample.uniquename.op('~')(self.lib_regex),
                unit.name.in_((self.xprn_types_to_report)),
            )
            results = session.query(sample, gene, unit, value).\
                select_from(gene).\
                join(LibraryFeature, (LibraryFeature.feature_id == gene.feature_id)).\
                join(sample, (sample.library_id == LibraryFeature.library_id)).\
                join(value, (value.library_feature_id == LibraryFeature.library_feature_id)).\
                join(unit, (unit.cvterm_id == value.type_id)).\
                filter(*filters).\
                distinct()
            # Process the data into data dicts for export.
            this_counter = 0
            this_data_dict = {}
            for result in results:
                xprn_section_to_use = xprn_section
                xprn_section_rank_to_use = xprn_section_rank
                unit_to_use = result.unit.name
                # Record the xprn_section, sample id and gene id as the data dict key for sorting.
                data_dict_key = (xprn_section_rank_to_use, result.sample.uniquename, result.gene.uniquename)
                # Build the dict itself.
                data_dict = {
                    'High_Throughput_Expression_Section': xprn_section_to_use,
                    'Dataset_ID': None,
                    'Dataset_Name': None,
                    'Sample_ID': result.sample.uniquename,
                    'Sample_Name': result.sample.name,
                    'Gene_ID': result.gene.uniquename,
                    'Gene_Symbol': result.gene.name,
                    'Expression_Unit': unit_to_use,
                    'Expression_Value': result.value.value
                }
                this_data_dict[data_dict_key] = data_dict
                this_counter += 1
            # Sort all data before sending it to the export list.
            data_keys = list(this_data_dict.keys())
            data_keys.sort()
            for i in data_keys:
                self.data_to_export.append(this_data_dict[i])
            counter += this_counter
            log.info(f'Found {this_counter} expression values for {sample_name}.')
        log.info(f'Found {counter} expression values for individual datasets overall.')
        return

    def query_chado(self, session):
        """Run query methods."""
        log.info('Starting "query_chado" method.')
        self.get_ht_project_data(session)
        self.get_ht_sample_data(session)
        log.info('Method "query_chado" is done.')
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
