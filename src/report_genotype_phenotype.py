# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report genotype-level phenotypes. Jira DB-813.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_genotype_phenotype.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_genotype_phenotype.py -v -c /path/to/config.cfg

"""

import argparse
import re
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker    # could add aliased here
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Feature, FeatureGenotype, Genotype, Phenotype, PhenotypeCvterm, Phenstatement, Pub
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Global variables for the output file. Header order will match list order below.
report_label = 'genotype_phenotype_data'
report_title = 'FlyBase genotype-phenotype report'
header_list = [
    'genotype_symbols',
    'genotype_FBids',
    'phenotype_name',
    'phenotype_id',
    'qualifier_names',
    'qualifier_ids',
    'reference'
]

# Proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
output_dir = set_up_dict['output_dir']
output_filename = set_up_dict['output_filename']
# report_filename = '{}{}_{}.report'.format(output_dir, report_label, database)
log = set_up_dict['log']
the_time = set_up_dict['the_time']

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
insp = inspect(engine)

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')

# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    phenotype_reporter = PhenotypeReporter()
    db_query_transaction(phenotype_reporter)
    data_to_export_as_tsv = generic_FB_tsv_dict(report_title, database)
    data_to_export_as_tsv['data'] = phenotype_reporter.process_database_info()
    tsv_report_dump(data_to_export_as_tsv, output_filename, headers=header_list)
    log.info('Ended main function.')


class PhenStatementAnnotation(object):
    """A PhenStatementAnnotation and related attributes."""

    def __init__(self, phenstatement):
        """Create a base PhenStatementAnnotation from a chado Phenotype object.

        Args:
            arg1 (phenstatement): (PhenStatementAnnotation) The Phenotype object.

        Returns:
            An object of the PhenStatementAnnotation class.

        """
        # Organism attributes, some filled in by downstream steps.
        self.phenstatement = phenstatement               # The PhenStatementAnnotation object.
        self.for_export = True                           # Change to False if any errors detected.
        self.reference = None                            # Propagated from related pub.
        self.genotype_symbols = None                     # Propagated from related genotype.
        self.genotype_FBids = None                       # Propagated from related genotype.
        self.phenotype_name = None                       # Propagated from related phenotype.
        self.phenotype_id = None                         # Propagated from related phenotype.
        self.qualifier_names = None                      # Propagated from related phenotype.
        self.qualifier_ids = None                        # Propagated from related phenotype.
        self.desc = None                                 # Will be derived from geno, pheno and ref info.
        self.errors = []                                 # Errors preventing data export.


class GenotypeAnnotation(object):
    """A GenotypeAnnotation and related attributes."""

    def __init__(self, genotype):
        """Create a base GenotypeAnnotation from a chado Genotype object.

        Args:
            arg1 (genotype): (Genotype) The Genotype object.

        Returns:
            An object of the GenotypeAnnotation class.

        """
        # Organism attributes, some filled in by downstream steps.
        self.genotype = genotype              # The Genotype object.
        self.components = []                  # Will be list of FB IDs for component features (excludes bogus).
        self.genotype_FBids = None            # Will be the IDs for genotype components, concatenated.
        self.genotype_symbols = None          # Will be the genotype name displayed in the export file.
        self.errors = []                      # Errors preventing data export.


class PhenotypeAnnotation(object):
    """A PhenotypeAnnotation and related attributes."""

    def __init__(self, phenotype):
        """Create a base PhenotypeAnnotation from a chado Phenotype object.

        Args:
            arg1 (phenotype): (Phenotype) The Phenotype object.

        Returns:
            An object of the PhenotypeAnnotation class.

        """
        # Organism attributes, some filled in by downstream steps.
        self.phenotype = phenotype    # The Phenotype object.
        self.pheno_cvterms = []       # Will be list of PhenotypeCvterm objects.
        self.phenotype_name = None    # Will be name for phenotype CV term.
        self.phenotype_id = None      # Will be the ID for the phenotype CV term.
        self.qualifier_names = ''     # Will be a string of qualifier names separated by pipes.
        self.qualifier_ids = ''       # Will be a string of qualifier term IDs separated by pipes.
        self.errors = []              # Errors preventing data export.


class PhenotypeReporter(object):
    """Create the PhenotypeReporter object."""

    def __init__(self):
        """Create the PhenotypeReporter object."""

    phenstmt_dict = {}     # phenstatement_id-keyed dict of PhenstatementAnnotation objects.
    pub_dict = {}          # pub_id-keyed dict of Pub objects
    geno_dict = {}         # genotype_id-keyed dict of GenotypeAnnotation objects.
    pheno_dict = {}        # phenotype_id-keyed dict of PhenotypeAnnotation objects.
    feat_name_dict = {}    # Feature_id-keyed dict of feature names for output file.
    cvterm_dict = {}       # Cvterm_id-keyed CV term dicts: {60858: {'name': 'lethal', 'id': 'FBcv:0000351'}, ...}

    def __get_phenstatements(self, session):
        """Get phenstatements."""
        log.info('Get phenstatements.')
        results = session.query(Phenstatement).distinct()
        counter = 0
        for result in results:
            self.phenstmt_dict[result.phenstatement_id] = PhenStatementAnnotation(result)
            counter += 1
        log.info(f'Found {counter} phenstatements.')
        return

    def __get_pubs(self, session):
        """Get pubs related to phenstatements."""
        log.info('Get pubs related to phenstatements.')
        results = session.query(Pub).\
            join(Phenstatement, (Phenstatement.pub_id == Pub.pub_id)).\
            distinct()
        counter = 0
        for result in results:
            self.pub_dict[result.pub_id] = result
            counter += 1
        log.info(f'Found {counter} pubs related to phenstatements.')
        return

    def __get_genotypes(self, session):
        """Get genotypes related to phenstatements."""
        log.info('Get genotypes related to phenstatements.')
        results = session.query(Genotype).\
            join(Phenstatement, (Phenstatement.genotype_id == Genotype.genotype_id)).\
            distinct()
        counter = 0
        for result in results:
            self.geno_dict[result.genotype_id] = GenotypeAnnotation(result)
            counter += 1
        log.info(f'Found {counter} genotypes related to phenstatements.')
        return

    def __get_features(self, session):
        """Get features related to genotypes."""
        log.info('Get features related to genotypes.')
        feat_rgx = r'^FB[a-z]{2}[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(feat_rgx)
        )
        results = session.query(Feature).\
            join(FeatureGenotype, (FeatureGenotype.feature_id == Feature.feature_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.feat_name_dict[result.uniquename] = result.name
            counter += 1
        log.info(f'Found {counter} current features related to genotypes.')
        return

    def __get_phenotypes(self, session):
        """Get phenotypes related to phenstatements."""
        log.info('Get phenotypes related to phenstatements.')
        # Start with phenotype table.
        results = session.query(Phenotype).\
            join(Phenstatement, (Phenstatement.phenotype_id == Phenotype.phenotype_id)).\
            distinct()
        counter = 0
        for result in results:
            self.pheno_dict[result.phenotype_id] = PhenotypeAnnotation(result)
            counter += 1
        log.info(f'Found {counter} phenotypes related to phenstatements.')
        # Then get phenotype_cvterm entries.
        results = session.query(PhenotypeCvterm).\
            join(Phenstatement, (Phenstatement.phenotype_id == PhenotypeCvterm.phenotype_id)).\
            distinct()
        counter = 0
        for result in results:
            self.pheno_dict[result.phenotype_id].pheno_cvterms.append(result)
            counter += 1
        log.info(f'Found {counter} phenotype_cvterm entries.')
        # Check for phenotypes with many qualifiers.
        test_cnt = 0
        for phenotype in self.pheno_dict.values():
            if len(phenotype.pheno_cvterms) > 1:
                test_cnt += 1
        log.debug(f'Found {test_cnt} phenotypes with many qualifiers.')
        return

    def __get_cvterms(self, session):
        """Get current CV terms and their IDs."""
        log.info('Get current CV terms and their IDs.')
        cvs_allowed = ['FBcv', 'FBbt', 'FBdv', 'GO', 'SO']
        filters = (
            Cvterm.is_obsolete == 0,
            Db.name.in_((cvs_allowed))
        )
        results = session.query(Cvterm, Db, Dbxref).\
            join(Dbxref, (Dbxref.dbxref_id == Cvterm.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            cvt = {
                'name': result.Cvterm.name,
                'id': f'{result.Db.name}:{result.Dbxref.accession}'
            }
            self.cvterm_dict[result.Cvterm.cvterm_id] = cvt
            counter += 1
        log.info(f'Found {counter} CV terms in these CVs: {cvs_allowed}.')
        # Check CV terms.
        for k, v in self.cvterm_dict.items():
            log.debug(f'CHECK CV TERMS: cvterm_id={k}, name={v["name"]}, id={v["id"]}.')
        return

    def __process_genotypes(self):
        """Process genotypes to get name and component IDs for display."""
        log.info('Process genotypes to get name and component IDs for display.')
        for genotype in self.geno_dict.values():
            if genotype.genotype.is_obsolete is True:
                genotype.errors.append('obsolete')
                continue
            if genotype.genotype.uniquename.startswith('PROBLEMATIC'):
                genotype.errors.append('problematic')
                continue
            if not genotype.genotype.description:
                genotype.errors.append('no description')
                continue
            # Convert description (component ID string) for export.
            # Start with IDs for each cgroup and convert each as needed (when we find bogus "[+]" or "[-]" symbols).
            edited_cgroups = []
            cgroups = genotype.genotype.description.split('_')
            for cgroup in cgroups:
                cgroup_parts = cgroup.split('|')
                edited_cgroup_parts = []
                for cgroup_part in cgroup_parts:
                    if cgroup_part.endswith('[+]'):
                        edited_cgroup_parts.append('+')
                    elif cgroup_part.endswith('[-]') and cgroup_part != 'Tn10\\tetR[-]':    # Temporary exception.
                        edited_cgroup_parts.append('-')
                    else:
                        edited_cgroup_parts.append(cgroup_part)
                # Different rules for cgroups with one or two parts.
                if len(edited_cgroup_parts) == 1:
                    edited_cgroup_id_str = f'{edited_cgroup_parts[0]}'
                else:
                    # Reverse components if we have generic "+" or "-" symbol.
                    if edited_cgroup_parts[0] == '-' or edited_cgroup_parts[0] == '+':
                        edited_cgroup_id_str = f'{edited_cgroup_parts[1]}/{edited_cgroup_parts[0]}'
                    # Otherwise, preserve order.
                    else:
                        edited_cgroup_id_str = f'{edited_cgroup_parts[0]}/{edited_cgroup_parts[1]}'
                edited_cgroups.append(edited_cgroup_id_str)
            # Reassemble the edited cgroups into a new component ID string.
            edited_cgroups.sort()
            genotype.genotype_FBids = ' '.join(edited_cgroups)
            # Now generate genotype label by replacing IDs in component ID string with feature names.
            if not genotype.genotype_FBids:
                continue
            genotype.genotype_symbols = genotype.genotype_FBids
            feat_id_rgx = 'FB[a-z]{2}[0-9]{7}'
            genotype.components = list(set(re.findall(feat_id_rgx, genotype.genotype_FBids)))
            for feat_id in genotype.components:
                genotype.genotype_symbols = genotype.genotype_symbols.replace(feat_id, self.feat_name_dict[feat_id])
            # Check test component IDs conversion and genotype label
            log.debug(f'CHECK IDs: desc={genotype.genotype.description}, ids={genotype.genotype_FBids}')
            log.debug(f'CHECK GENO NAME: uname={genotype.genotype.uniquename}, symbol={genotype.genotype_symbols}')
        return

    def __process_phenotypes(self):
        """Process phenotypes for export."""
        log.info('Process phenotypes for export.')
        for pheno in self.pheno_dict.values():
            # Start with main phenotype CV term.
            primary_cvterm = None
            # Phenotype will have a CV term in either observable_id or cvalue_id column (not both).
            # First check observable_id.
            try:
                primary_cvterm = self.cvterm_dict[pheno.phenotype.observable_id]
            except KeyError:
                if pheno.phenotype.observable_id != 60468:    # i.e., ignore "unspecified"
                    pheno.errors.append(f'internal or obsolete cvterm_id={pheno.phenotype.observable_id}')
            # Then check cvalue_id.
            try:
                primary_cvterm = self.cvterm_dict[pheno.phenotype.cvalue_id]
            except KeyError:
                if pheno.phenotype.cvalue_id != 60468:    # i.e., ignore "unspecified"
                    pheno.errors.append(f'internal or obsolete cvterm_id={pheno.phenotype.cvalue_id}')
            if primary_cvterm:
                pheno.phenotype_name = primary_cvterm['name']
                pheno.phenotype_id = primary_cvterm['id']
            else:
                pheno.errors.append('primary cvterm is internal or obsolete')
            # Now get qualifiers.
            cvterm_list = []
            bad_cvterm_counter = 0
            for pheno_cvterm in pheno.pheno_cvterms:
                try:
                    cvterm_list.append(self.cvterm_dict[pheno_cvterm.cvterm_id])
                except KeyError:
                    bad_cvterm_counter += 1
                    pheno.errors.append(f'internal or obsolete cvterm_id={pheno_cvterm.cvterm_id}')
            if cvterm_list and not bad_cvterm_counter:
                pheno.qualifier_names = '|'.join([i['name'] for i in cvterm_list])
                pheno.qualifier_ids = '|'.join([i['id'] for i in cvterm_list])
            elif bad_cvterm_counter:
                pheno.errors.append('qualifier cvterm is internal or obsolete')
        return

    def __process_phenstatements(self):
        """Process phenstatements for export."""
        log.info('Process phenstatements for export.')
        for phenstmt in self.phenstmt_dict.values():
            # Look up pub and propagate info.
            pub = self.pub_dict[phenstmt.phenstatement.pub_id]
            phenstmt.reference = pub.uniquename
            if pub.is_obsolete is True:
                phenstmt.errors.append('pub: obsolete')
            # Look up genotype and propagate info.
            genotype = self.geno_dict[phenstmt.phenstatement.genotype_id]
            phenstmt.genotype_symbols = genotype.genotype_symbols
            phenstmt.genotype_FBids = genotype.genotype_FBids
            if genotype.errors:
                phenstmt.errors.append(f"genotype: {', '.join(genotype.errors)}")
            # Look up phenotype and propagate info.
            phenotype = self.pheno_dict[phenstmt.phenstatement.phenotype_id]
            phenstmt.phenotype_name = phenotype.phenotype_name
            phenstmt.phenotype_id = phenotype.phenotype_id
            phenstmt.qualifier_names = phenotype.qualifier_names
            phenstmt.qualifier_ids = phenotype.qualifier_ids
            if phenotype.errors:
                phenstmt.errors.append(f"phenotype: {', '.join(phenotype.errors)}")
            # Create description.
            id = phenstmt.phenstatement.phenstatement_id
            geno_desc = genotype.genotype_symbols
            pheno_desc = phenotype.phenotype.uniquename
            ref_desc = phenstmt.reference
            phenstmt.desc = f'phenstatement_id={id}, genotype={geno_desc}, phenotype={pheno_desc}, ref={ref_desc}'
            # Export filter.
            if phenstmt.errors:
                phenstmt.for_export = False
                log.warning(f'REJECT: {phenstmt.desc}. REASONS: {phenstmt.errors}')
        return

    def query_chado(self, session):
        """Run wrapper method for querying and processing info."""
        self.__get_phenstatements(session)
        self.__get_pubs(session)
        self.__get_genotypes(session)
        self.__get_features(session)
        self.__get_phenotypes(session)
        self.__get_cvterms(session)
        self.__process_genotypes()
        self.__process_phenotypes()
        self.__process_phenstatements()
        return

    def process_database_info(self):
        """Process PhenStatementAnnotations into export dicts."""
        log.info('Processing phenstatements into export dicts.')
        # First create genotype-grouped lists of phenstatements.
        counter = 0
        export_counter = 0
        reject_counter = 0
        grouped_phenstmts = {}
        for phenstmt in self.phenstmt_dict.values():
            counter += 1
            if phenstmt.for_export is False:
                reject_counter += 1
                continue
            try:
                grouped_phenstmts[phenstmt.genotype_symbols].append(phenstmt)
            except KeyError:
                grouped_phenstmts[phenstmt.genotype_symbols] = [phenstmt]
            export_counter += 1
        # Now export in order of genotype.
        ordered_genotypes = list(grouped_phenstmts.keys())
        ordered_genotypes.sort()
        output_data = []
        for genotype in ordered_genotypes:
            for phenstmt in grouped_phenstmts[genotype]:
                phenstmt_output_dict = {
                    'genotype_symbols': phenstmt.genotype_symbols,
                    'genotype_FBids': phenstmt.genotype_FBids,
                    'phenotype_name': phenstmt.phenotype_name,
                    'phenotype_id': phenstmt.phenotype_id,
                    'qualifier_names': phenstmt.qualifier_names,
                    'qualifier_ids': phenstmt.qualifier_ids,
                    'reference': phenstmt.reference
                }
                output_data.append(phenstmt_output_dict)
        log.info(f'Of {counter} phenstatements, exported {export_counter} and rejected {reject_counter}.')
        return output_data


def db_query_transaction(object_to_execute):
    """Run an object's "query_chado()" method.

    Args:
        arg1 (object_to_execute): An object that has an SQL ORM "query_chado()" method.

    Returns:
        None.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    log.info('Querying chado.')
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        object_to_execute.query_chado(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred, rolling back and exiting.')
        raise
    return


if __name__ == "__main__":
    main()
