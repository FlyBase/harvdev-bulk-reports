# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate FlyCyc files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    generate_flycyc_files.py [-h] [-v VERBOSE] [-f FASTA] [-c CONFIG]

Example:
    python generate_flycyc_files.py -v -f -c /foo/bar/config.cfg

Notes:
    The script reports relevant data for genes on major chromosome scaffolds.
    If "-f" option is specified, it also generates FASTA files for each of the
    major chromosome scaffolds; ~10m/chr to generate FASTA files.

"""

import argparse
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
from textwrap import TextWrapper
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.production import (
    Cv, Cvterm, CvtermDbxref, Db, Dbxref, Dbxrefprop, Feature, FeatureCvterm,
    FeatureCvtermprop, FeatureDbxref, FeatureRelationship, FeatureSynonym,
    Featureloc, Organism, Pub, PubDbxref, Synonym
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Important label for output files.
report_label = 'flycyc'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
annotation_release = set_up_dict['annotation_release']
database_release = set_up_dict['database_release']
output_dir = set_up_dict['output_dir']
log = set_up_dict['log']

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
insp = inspect(engine)

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-f', '--fasta', action='store_true', help='Write out fasta files.', required=False)
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
chr_fasta = args.fasta


# The main process.
def main():
    """Run the steps for generating FlyCyc files."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')

    flycyc_generator = FlyCycGenerator()
    db_query_transaction(flycyc_generator)
    flycyc_generator.print_files()

    log.info('Ended main function.\n')


class FlyCycGenerator(object):
    """This object gets data from chado and generates FlyCyc files."""
    def __init__(self):
        """Create the FlyCycGenerator object."""
        self.chr_dict = {}              # Will be a uniquename-keyed dict of Feature objects representing each major chr scaffold.
        self.chr_id_dict = {}           # Will be a feature_id-uniquename dict for major chr scaffolds.
        self.chr_gene_dict = {}         # Will be a chr-uniquename-keyed list of gene dicts.
        self.fbrf_to_pmid = {}          # Will be an FBrf-to-PubMed ID dict.
        self.gene_product_type = {}     # Will be an FBgn-keyed dict of product type (e.g., "P", "TRNA", etc.)
        self.trpt_exon_locs = {}        # Will be FBtr-keyed dict of exon location lists.
        self.gene_xref_dict = {}        # Will be FBgn-keyed dict of dbxref lists for each gene.
        self.gene_fullname_dict = {}    # Will be FBgn-keyed dict of gene fullname.
        self.gene_synonym_dict = {}     # Will be FBgn-keyed dict of non-current synonym lists for each gene.
        self.gene_go_annos = {}         # Will be FBgn-keyed dict of simple GO annotations dicts.
        self.gene_neg_go_annos = {}     # Will be FBgn-keyed dict of lists GO term IDs from negative GO annotations for each gene.
        self.go_ec_dict = {}            # Will be a GO ID-keyed dict of a list of related EC numbers (no "GO:" prefix).
        self.go_metacyc_dict = {}       # Will be a GO ID-keyed dict of a list of related METACYC xrefs (no "GO:" prefix).

    # This list defines the "major" chr scaffolds for which reports are generated.
    chr_scaffolds_to_report = [
        'X',
        'Y',
        '2L',
        '2R',
        '3L',
        '3R',
        '4',
        'mitochondrion_genome',
        'rDNA',
        'Unmapped_Scaffold_8_D1580_D1567'
    ]
    # chr_scaffolds_to_report = ['mitochondrion_genome', 'rDNA']    # For faster testing, limit to sample chr test set.

    def query_chr(self, session):
        """Get Dmel chr sequences."""
        log.info('Retrieving chromosome scaffold sequences for FASTA output.')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.in_((self.chr_scaffolds_to_report)),
            Organism.abbreviation == 'Dmel',
            Cvterm.name == 'golden_path'
        )
        chr_results = session.query(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            filter(*filters).\
            distinct()
        for result in chr_results:
            self.chr_dict[result.uniquename] = result
            self.chr_id_dict[result.feature_id] = result.uniquename
        return

    def generate_fbrf_pubmed_dict(self, session):
        """Create an FBrf-to-PMID dict."""
        pub_regex = r'FBrf[0-9]{7}$'
        filters = (
            Pub.uniquename.op('~')(pub_regex),
            Pub.is_obsolete.is_(False),
            PubDbxref.is_current.is_(True),
            Db.name == 'pubmed'
        )
        pubmed_xrefs = session.query(Pub, PubDbxref, Dbxref).\
            join(PubDbxref, (PubDbxref.pub_id == Pub.pub_id)).\
            join(Dbxref, (Dbxref.dbxref_id == PubDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        for xref in pubmed_xrefs:
            self.fbrf_to_pmid[xref.Pub.uniquename] = xref.Dbxref.accession
        return

    def query_gene_negative_go_annotations(self, session):
        """Get negative GO annotations for genes."""
        log.info('Getting negative GO annotations for genes.')
        gene_uniquename_regex = r'^FBgn[0-9]{7}$'
        go_cv_list = ['biological_process', 'cellular_component', 'molecular_function']
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(gene_uniquename_regex),
            Organism.abbreviation == 'Dmel',
            Cvterm.is_obsolete == 0,
            FeatureCvterm.is_not.is_(True),
            Cv.name.in_((go_cv_list))
        )
        neg_go_anno_results = session.query(Feature, Cvterm).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
            filter(*filters).\
            distinct()
        for result in neg_go_anno_results:
            try:
                self.gene_neg_go_annos[result.Feature.uniquename].append(result.Cvterm.dbxref.accession)
            except KeyError:
                self.gene_neg_go_annos[result.Feature.uniquename] = [result.Cvterm.dbxref.accession]
        return

    def query_gene_go_annotations(self, session):
        """Get GO annotations for genes."""
        log.info('Getting GO annotations for genes.')
        # Define the evidence_code key:
        evidence_code_dict = {
            'inferred from mutant phenotype': 'IMP',
            'inferred from genetic interaction': 'IGI',
            'inferred from physical interaction': 'IPI',
            'inferred from sequence or structural similarity': 'ISS',
            'inferred from sequence model': 'ISM',
            'inferred from sequence alignment': 'ISA',
            'inferred from sequence orthology': 'ISO',
            'inferred from experiment': 'EXP',
            'inferred from direct assay': 'IDA',
            'inferred from electronic annotation': 'IEA',
            'inferred from expression pattern': 'IEP',
            'inferred from reviewed computational analysis': 'RCA',
            'traceable author statement': 'TAS',
            'non-traceable author statement': 'NAS',
            'inferred by curator': 'IC',
            'inferred from genomic context': 'IGC',
            'no biological data available': 'ND',
            'inferred from biological aspect of ancestor': 'IBA',
            'inferred from biological aspect of descendant': 'IBD',
            'inferred from key residues': 'IKR',
            'inferred from rapid divergence': 'IRD',
            'inferred from high throughput experiment': 'HTP',
            'inferred from high throughput direct assay': 'HDA',
            'inferred from high throughput expression pattern': 'HEP',
            'inferred from high throughput genetic interaction': 'HGI',
            'inferred from high throughput mutant phenotype': 'HMP'
        }
        # Now get the GO annotations.
        gene_uniquename_regex = r'^FBgn[0-9]{7}$'
        go_cv_list = ['biological_process', 'cellular_component', 'molecular_function']
        cvterm = aliased(Cvterm, name='cvterm')
        qualifier_type = aliased(Cvterm, name='qualifier_type')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(gene_uniquename_regex),
            Organism.abbreviation == 'Dmel',
            FeatureCvterm.is_not.is_(False),
            cvterm.is_obsolete == 0,
            Cv.name.in_((go_cv_list)),
            qualifier_type.name == 'evidence_code',
            Pub.is_obsolete.is_(False),
            Pub.uniquename != 'unattributed'
        )
        go_anno_results = session.query(Feature, cvterm, Pub, FeatureCvtermprop).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(cvterm, (cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(qualifier_type, (qualifier_type.cvterm_id == FeatureCvtermprop.type_id)).\
            join(Pub, (Pub.pub_id == FeatureCvterm.pub_id)).\
            filter(*filters).\
            distinct()
        for result in go_anno_results:
            go_anno_dict = {
                'name': result.cvterm.name,
                'go_id': result.cvterm.dbxref.accession,
                'pub': result.Pub.uniquename
            }
            # Filter out any gene-GO term annotation that has some negative annotations.
            if result.Feature.uniquename in self.gene_neg_go_annos.keys():
                if result.cvterm.dbxref.accession in self.gene_neg_go_annos[result.Feature.uniquename]:
                    log.debug('For {} ({}): negative annotation to "{}" (GO:{}) will be filtered out.'.
                              format(result.Feature.uniquename, result.Feature.name, result.cvterm.name, result.cvterm.dbxref.accession))
                    continue
            # Tease out the evidence.
            evidence_code = result.FeatureCvtermprop.value
            for evidence_type in evidence_code_dict.keys():
                evidence_code = evidence_code.replace(evidence_type, evidence_code_dict[evidence_type])
            go_anno_dict['evi_code'] = evidence_code.split(' ')[0]
            # Get the PMID, if available.
            try:
                go_anno_dict['pub'] = self.fbrf_to_pmid[go_anno_dict['pub']]
            except KeyError:
                pass
            # Derive the string.
            go_anno_dict['string'] = '{}|{}|{}|{}'.format(go_anno_dict['name'], go_anno_dict['go_id'], go_anno_dict['pub'], go_anno_dict['evi_code'])
            try:
                self.gene_go_annos[result.Feature.uniquename].append(go_anno_dict)
            except KeyError:
                self.gene_go_annos[result.Feature.uniquename] = [go_anno_dict]
        return

    def query_go_ec_numbers(self, session):
        """Get EC numbers related to each GO term."""
        log.info('Getting GO-EC associations.')
        go_cv_list = ['biological_process', 'cellular_component', 'molecular_function']
        cvterm = aliased(Cvterm, name='cvterm')
        proptype = aliased(Cvterm, name='proptype')
        filters = (
            cvterm.is_obsolete == 0,
            Cv.name.in_((go_cv_list)),
            Db.name == 'EC',
            proptype.name == 'ec_description'
        )
        go_ec_results = session.query(cvterm, Dbxref).\
            join(Cv, (Cv.cv_id == cvterm.cv_id)).\
            join(CvtermDbxref, (CvtermDbxref.cvterm_id == cvterm.cvterm_id)).\
            join(Dbxref, (Dbxref.dbxref_id == CvtermDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            join(Dbxrefprop, (Dbxrefprop.dbxref_id == Dbxref.dbxref_id)).\
            join(proptype, (proptype.cvterm_id == Dbxrefprop.type_id)).\
            filter(*filters).\
            distinct()
        for result in go_ec_results:
            try:
                self.go_ec_dict[result.cvterm.dbxref.accession].append(result.Dbxref.accession)
            except KeyError:
                self.go_ec_dict[result.cvterm.dbxref.accession] = [result.Dbxref.accession]
        return

    def query_go_metacyc(self, session):
        """Get METACYC xrefs related to each GO term."""
        log.info('Getting GO-METACYC associations.')
        filters = (
            Cvterm.is_obsolete == 0,
            Cv.name == 'molecular_function',
            Db.name == 'MetaCyc',
        )
        go_metacyc_results = session.query(Cvterm, Dbxref).\
            join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
            join(CvtermDbxref, (CvtermDbxref.cvterm_id == Cvterm.cvterm_id)).\
            join(Dbxref, (Dbxref.dbxref_id == CvtermDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        for result in go_metacyc_results:
            try:
                self.go_metacyc_dict[result.Cvterm.dbxref.accession].append(result.Dbxref.accession)
            except KeyError:
                self.go_metacyc_dict[result.Cvterm.dbxref.accession] = [result.Dbxref.accession]
        return

    def query_gene_xrefs(self, session):
        """Get dbxrefs for genes."""
        dbxrefs_to_get = {
            'EntrezGene': 'GeneID',
            'UniProt/GCRP': 'UNIPROT',
            'RNAcentral': 'RNAcentral'
        }
        log.info('Getting gene dbxrefs for these databases: {}'.format(dbxrefs_to_get.keys()))
        gene_uniquename_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(gene_uniquename_regex),
            Cvterm.name == 'gene',
            Organism.abbreviation == 'Dmel',
            FeatureDbxref.is_current.is_(True),
            Db.name.in_((dbxrefs_to_get.keys()))
        )
        dbxref_results = session.query(Feature, Db, Dbxref).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        for result in dbxref_results:
            xref_name = '{}:{}'.format(dbxrefs_to_get[result.Db.name], result.Dbxref.accession)
            try:
                self.gene_xref_dict[result.Feature.uniquename].append(xref_name)
            except KeyError:
                self.gene_xref_dict[result.Feature.uniquename] = [xref_name]
        return

    def query_transcript_exon_locations_bad(self, session):
        """Get gene exon locations."""
        log.info('Getting gene exon locations.')
        transcript_uniquename_regex = r'^FBtr[0-9]{7}$'
        transcript = aliased(Feature, name='transcript')
        transcript_part = aliased(Feature, name='transcript_part')
        chr = aliased(Feature, name='chr')
        rel_type = aliased(Cvterm, name='rel_type')
        part_type = aliased(Cvterm, name='part_type')
        chr_type = aliased(Cvterm, name='chr_type')
        filters = (
            transcript.is_obsolete.is_(False),
            transcript.uniquename.op('~')(transcript_uniquename_regex),
            Organism.abbreviation == 'Dmel',
            transcript_part.is_obsolete.is_(False),
            chr.is_obsolete.is_(False),
            rel_type.name == 'partof',
            part_type.name == 'exon',
            chr_type == 'golden_path',
        )
        exon_locs = session.query(transcript, Featureloc).\
            select_from(transcript).\
            join(Organism, (Organism.organism_id == transcript.organism_id)).\
            join(FeatureRelationship, (FeatureRelationship.object_id == transcript.feature_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            join(transcript_part, (transcript_part.feature_id == FeatureRelationship.subject_id)).\
            join(part_type, (part_type.cvterm_id == transcript_part.type_id)).\
            join(Featureloc, (Featureloc.feature_id == transcript_part.feature_id)).\
            join(chr, (chr.feature_id == Featureloc.srcfeature_id)).\
            join(chr_type, (chr_type.cvterm_id == chr.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for exon_loc in exon_locs:
            trpt_name = f'{exon_loc.transcript.name} ({exon_loc.transcript.uniquename})'
            if exon_loc.Featureloc.strand == -1:
                exon_string = f'{exon_loc.Featureloc.fmax}--{exon_loc.Featureloc.fmin + 1}'
                log.debug(f'BOB MINUS: {trpt_name} exon at {exon_string}')
            else:
                exon_string = f'{exon_loc.Featureloc.fmin + 1}--{exon_loc.Featureloc.fmax}'
                log.debug(f'BOB PLUS: {trpt_name} exon at {exon_string}')
            try:
                self.trpt_exon_locs[exon_loc.transcript.uniquename].append(exon_string)
            except KeyError:
                self.trpt_exon_locs[exon_loc.transcript.uniquename] = [exon_string]
            counter += 1
        log.info(f'Found {counter} exons for {len(self.trpt_exon_locs.keys())} current Dmel transcripts.')
        return

    def query_transcript_exon_locations(self, session):
        """Get gene exon locations."""
        log.info('Getting gene exon locations.')
        transcript_uniquename_regex = r'^FBtr[0-9]{7}$'
        transcript = aliased(Feature, name='transcript')
        transcript_part = aliased(Feature, name='transcript_part')
        chr = aliased(Feature, name='chr')
        rel_type = aliased(Cvterm, name='rel_type')
        part_type = aliased(Cvterm, name='part_type')
        chr_type = aliased(Cvterm, name='chr_type')
        filters = (
            transcript.is_obsolete.is_(False),
            transcript.uniquename.op('~')(transcript_uniquename_regex),
            Organism.abbreviation == 'Dmel',
            transcript_part.is_obsolete.is_(False),
            rel_type.name == 'partof',
            part_type.name == 'exon',
        )
        exon_locs = session.query(transcript, transcript_part).\
            select_from(transcript).\
            join(Organism, (Organism.organism_id == transcript.organism_id)).\
            join(FeatureRelationship, (FeatureRelationship.object_id == transcript.feature_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            join(transcript_part, (transcript_part.feature_id == FeatureRelationship.subject_id)).\
            join(part_type, (part_type.cvterm_id == transcript_part.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for exon_loc in exon_locs:
            counter += 1
        log.info(f'Found {counter} exon results.')
        exit()
        return

    def query_gene_fullnames(self, session):
        """Get gene current full names."""
        log.info('Getting gene full names.')
        feature_type = aliased(Cvterm, name='feature_type')
        synonym_type = aliased(Cvterm, name='synonym_type')
        gene_uniquename_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(gene_uniquename_regex),
            feature_type.name == 'gene',
            Organism.abbreviation == 'Dmel',
            FeatureSynonym.is_current.is_(True),
            synonym_type.name == 'fullname'
        )
        fullname_results = session.query(Feature, Synonym).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureSynonym, (FeatureSynonym.feature_id == Feature.feature_id)).\
            join(Synonym, (Synonym.synonym_id == FeatureSynonym.synonym_id)).\
            join(synonym_type, (synonym_type.cvterm_id == Synonym.type_id)).\
            filter(*filters).\
            distinct()
        for result in fullname_results:
            self.gene_fullname_dict[result.Feature.uniquename] = result.Synonym.name
        return

    def query_gene_synonyms(self, session):
        """Get gene synonyms."""
        log.info('Getting gene synonyms.')
        gene_uniquename_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(gene_uniquename_regex),
            Cvterm.name == 'gene',
            Organism.abbreviation == 'Dmel',
            FeatureSynonym.is_current.is_(False)
        )
        synonym_results = session.query(Feature, Synonym).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureSynonym, (FeatureSynonym.feature_id == Feature.feature_id)).\
            join(Synonym, (Synonym.synonym_id == FeatureSynonym.synonym_id)).\
            filter(*filters).\
            distinct()
        for result in synonym_results:
            try:
                self.gene_synonym_dict[result.Feature.uniquename].append(result.Synonym.name)
            except KeyError:
                self.gene_synonym_dict[result.Feature.uniquename] = [result.Synonym.name]
        return

    def query_gene_products(self, session):
        """Query genes and related products."""
        log.info('Querying for gene product types.')
        product_type_key = {
            'mRNA': 'P',
            'pseudogene': 'PSEUDO',
            'tRNA': 'TRNA',
            'rRNA': 'RRNA',
            'miRNA': 'MISC-RNA',
            'ncRNA': 'MISC-RNA',
            'pre_miRNA': 'MISC-RNA',
            'snRNA': 'MISC-RNA',
            'snoRNA': 'MISC-RNA'
        }
        gene_feature = aliased(Feature, name='gene_feature')
        transcript_feature = aliased(Feature, name='transcript_feature')
        gene_type = aliased(Cvterm, name='gene_type')
        transcript_type = aliased(Cvterm, name='transcript_type')
        gene_uniquename_regex = r'^FBgn[0-9]{7}$'
        transcript_uniquename_regex = r'FBtr[0-9]{7}$'
        filters = (
            gene_feature.is_obsolete.is_(False),
            gene_feature.uniquename.op('~')(gene_uniquename_regex),
            gene_type.name == 'gene',
            transcript_feature.is_obsolete.is_(False),
            transcript_feature.uniquename.op('~')(transcript_uniquename_regex),
            Organism.abbreviation == 'Dmel'
        )
        gene_product_results = session.query(gene_feature, transcript_type).\
            join(gene_type, (gene_type.cvterm_id == gene_feature.type_id)).\
            join(Organism, (Organism.organism_id == gene_feature.organism_id)).\
            join(FeatureRelationship, (FeatureRelationship.object_id == gene_feature.feature_id)).\
            join(transcript_feature, (transcript_feature.feature_id == FeatureRelationship.subject_id)).\
            join(Featureloc, (Featureloc.feature_id == transcript_feature.feature_id)).\
            join(transcript_type, (transcript_type.cvterm_id == transcript_feature.type_id)).\
            filter(*filters).\
            distinct()
        for result in gene_product_results:
            self.gene_product_type[result.gene_feature.uniquename] = product_type_key[result.transcript_type.name]
        return

    def query_genes(self, session):
        """Get Dmel genes and create a chr dict of gene lists."""
        log.info('Retrieving localized Dmel genes.')
        # First create chr-keyed dict of genes.
        for chr_uniquename in self.chr_scaffolds_to_report:
            self.chr_gene_dict[chr_uniquename] = []
        gene_uniquename_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(gene_uniquename_regex),
            Cvterm.name == 'gene',
            Organism.abbreviation == 'Dmel'
        )
        # Note: This query does not pull in chr Feature objects corresponding to Featureloc.srcfeature_id.
        #     That's because doing so causes the query to hang.
        #     I imagine it's taxing to pull in chr features with 30M bases for each gene result.
        #     Instead, to get chr info, I use self.chr_id_dict.
        gene_results = session.query(Feature, Featureloc).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(Featureloc, (Featureloc.feature_id == Feature.feature_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(FeatureRelationship, (FeatureRelationship.object_id == Feature.feature_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in gene_results:
            counter += 1
            if counter % 1000 == 0:
                log.debug('Processing gene number {}'.format(counter))
            # Filter for genes on major chr scaffolds.
            try:
                chr_uniquename = self.chr_id_dict[result.Featureloc.srcfeature_id]
            except KeyError:
                log.warning('Skipping gene on minor chr scaffold: {} ({})'.
                            format(result.Feature.name, result.Feature.uniquename))
                continue
            # log.debug('Processing gene on major chr scaffold: {} ({})'.
            #           format(result.Feature.name, result.Feature.uniquename))
            gene_dict = {
                'ID': result.Feature.uniquename,
                'NAME': result.Feature.name,
                'DBLINK': [
                    'FlyBase:{}'.format(result.Feature.uniquename),
                    'Alliance:FB:{}'.format(result.Feature.uniquename)
                ],
                'CODING-SEGMENT': [],
                'SYNONYM': [],
                'GO': [],
                'EC': [],
                'METACYC': []
            }
            # Add coordinates.
            if result.Featureloc.strand == 1:
                gene_dict['STARTBASE'] = str(result.Featureloc.fmin + 1)
                gene_dict['ENDBASE'] = str(result.Featureloc.fmax)
            else:
                gene_dict['STARTBASE'] = str(result.Featureloc.fmax)
                gene_dict['ENDBASE'] = str(result.Featureloc.fmin + 1)
            # Add product_type. Filter out if there is no transcript: e.g., mt:ori FBgn0013687).
            try:
                gene_dict['PRODUCT-TYPE'] = self.gene_product_type[result.Feature.uniquename]
            except KeyError:
                log.warning('Skipping gene lacking transcript: {} ({})'.format(result.Feature.name, result.Feature.uniquename))
                continue
            # Add fullname.
            try:
                gene_dict['FUNCTION'] = self.gene_fullname_dict[result.Feature.uniquename]
            except KeyError:
                gene_dict['FUNCTION'] = result.Feature.name
            # Add dbxrefs.
            try:
                gene_dict['DBLINK'].extend(self.gene_xref_dict[result.Feature.uniquename])
            except KeyError:
                pass
            # Add synonyms.
            try:
                gene_dict['SYNONYM'] = self.gene_synonym_dict[result.Feature.uniquename]
            except KeyError:
                pass
            # Add GO terms.
            try:
                gene_dict['GO'] = [i['string'] for i in self.gene_go_annos[result.Feature.uniquename]]
            except KeyError:
                pass
            # Add EC numbers.
            if result.Feature.uniquename in self.gene_go_annos.keys():
                for go_id in [i['go_id'] for i in self.gene_go_annos[result.Feature.uniquename]]:
                    try:
                        gene_dict['EC'].extend(self.go_ec_dict[go_id])
                    except KeyError:
                        pass
                gene_dict['EC'] = set(gene_dict['EC'])
            # Add METACYC xrefs.
            if result.Feature.uniquename in self.gene_go_annos.keys():
                for go_id in [i['go_id'] for i in self.gene_go_annos[result.Feature.uniquename]]:
                    try:
                        gene_dict['METACYC'].extend(self.go_metacyc_dict[go_id])
                    except KeyError:
                        pass
                gene_dict['METACYC'] = set(gene_dict['METACYC'])
            # Append gene info to the appropriate chr.
            self.chr_gene_dict[chr_uniquename].append(gene_dict)
        return

    def query_chado(self, session):
        """Wrapper query method."""
        # There are some dependencies on the order of these methods.
        self.query_transcript_exon_locations(session)
        self.query_chr(session)
        self.generate_fbrf_pubmed_dict(session)
        self.query_gene_negative_go_annotations(session)
        # The query_gene_go_annotations() method should be run after:
        # 1. generate_fbrf_pubmed_dict()
        # 2. query_gene_negative_go_annotations()
        self.query_gene_go_annotations(session)
        self.query_go_ec_numbers(session)
        self.query_go_metacyc(session)
        self.query_gene_xrefs(session)
        self.query_gene_fullnames(session)
        self.query_gene_synonyms(session)
        self.query_gene_products(session)
        # The query_genes() method should be run last (needs info from methods above).
        self.query_genes(session)
        return

    def print_chr_files(self):
        """Print out FASTA file for each chromosome scaffold."""
        if not chr_fasta:
            log.info('Skipping FASTA output of chr scaffold sequences.')
            return
        log.info('Printing chromosome scaffold FASTA files.')
        for chr_uniquename in self.chr_scaffolds_to_report:
            chr = self.chr_dict[chr_uniquename]
            log.info('Processing FASTA sequence for chr scaffold {}'.format(chr_uniquename))
            output_filename = 'chr_{}.fsa'.format(chr_uniquename)
            output_file = open(output_filename, 'w')
            header = '>D. melanogaster Release 6 chromosome {}\n'.format(chr_uniquename)
            output_file.write(header)
            text_wrapper = TextWrapper(width=80)
            wrapped_sequence = '\n'.join(text_wrapper.wrap(chr.residues))
            output_file.write(wrapped_sequence)
        log.info('Done printing FASTA files.')
        return

    def print_gene_files(self):
        """Print out FlyCyc gene files for each chromosome."""
        log.info('Printing FlyCyc gene files.')
        for chr_uniquename in self.chr_scaffolds_to_report:
            log.info('Printing report for genes on chr scaffold {}'.format(chr_uniquename))
            output_filename = '{}chr_{}.pf'.format(output_dir, chr_uniquename)
            output_file = open(output_filename, 'w')
            header = ';;\n;; Pathologic format annotation file for D. melanogaster Release {}/FB{} chromosome {}\n;;\n;;\n'.\
                format(annotation_release.replace('R', ''), database_release, chr_uniquename)
            output_file.write(header)
            field_list = [
                'ID',
                'NAME',
                'FUNCTION',
                'SYNONYM',
                'PRODUCT-TYPE',
                'STARTBASE',
                'ENDBASE',
                'CODING-SEGMENT',
                'GO',
                'EC',
                'METACYC',
                'DBLINK'
            ]
            for gene in self.chr_gene_dict[chr_uniquename]:
                for field in field_list:
                    if type(gene[field]) in [set, list]:
                        for element in gene[field]:
                            output_file.write('{}\t{}\n'.format(field, element))
                    else:
                        output_file.write('{}\t{}\n'.format(field, gene[field]))
                output_file.write('//\n\n')
        return

    def print_files(self):
        """Wrapper write method."""
        self.print_gene_files()
        self.print_chr_files()
        return


def db_query_transaction(object_to_execute):
    """Query the chado database given an object that has a "query_chado()" method.

    Function assumes a global "engine" variable for SQLAlchemy processes.

    Args:
        arg1 (object_to_execute): Some object that has an SQL ORM "query_chado()" method.

    Returns:
        An sqlalchemy session that can be used later on for updating the queried objects.

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
