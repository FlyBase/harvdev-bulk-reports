#!/usr/bin/env python3

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate JSON files of non-coding RNAs for RNAcentral submission.

Author(s):
    Dave Emmert
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_rnacentral_json.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_rnacentral_json.py -v -c /foo/bar/config.cfg

"""

import calendar
from datetime import datetime, timezone
from harvdev_utils.psycopg_functions import (
    connect, set_up_db_reading
)
import json
import re
import time

# Important label for output files.
report_label = 'ncRNA_genes'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
database_host = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
annotation_release = set_up_dict['annotation_release']
database_release = set_up_dict['database_release']
output_dir = set_up_dict['output_dir']
output_filename = '{}{}_{}.json'.format(output_dir, report_label, database)
log = set_up_dict['log']
conn = set_up_dict['conn']


def main():
    """Retrieve ncRNA info and export as RNAcentral-compliant JSON file."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')

    # Instantiate main export dict and add metadata.
    ncrna_dict = {}
    ncrna_dict['metaData'] = {}
    ncrna_dict['metaData']['dataProvider'] = 'FlyBase'
    ncrna_dict['metaData']['publications'] = ['PMID:35266522']
    ncrna_dict['metaData']['schemaVersion'] = '0.4.0'
    ncrna_dict['metaData']['release'] = 'fb_' + database_release
    ncrna_dict['metaData']['genomicCoordinateSystem'] = '1-start, fully-closed'
    iso_time = datetime.fromtimestamp(calendar.timegm(time.gmtime()), tz=timezone.utc).isoformat()
    ncrna_dict['metaData']['dateProduced'] = iso_time

    # Instantiate data list and fill it up.
    ncrna_dict['data'] = []
    get_ncrna_json(database, ncrna_dict)
    with open(output_filename, 'w') as outfile:
        json.dump(ncrna_dict, outfile, indent=5, separators=(', ', ': '), ensure_ascii=False)

    log.info('Done writing data to JSON output file.')
    log.info('Ended main function.\n')


def pop_assembly_dict(assembly_dict):
    # Hard-coded dict of assembly info not in chado (doubleblech) (release number and GB GCA acc#)
    # See comments from DB-479
    assembly_dict['Dmel'] = {}
    assembly_dict['Dmel']['rel'] = 'R6'
    assembly_dict['Dmel']['gbacc'] = 'GCA_000001215.4'

    assembly_dict['Dana'] = {}
    assembly_dict['Dana']['rel'] = 'R1'
    assembly_dict['Dana']['gbacc'] = 'GCA_000005115.1'

    assembly_dict['Dere'] = {}
    assembly_dict['Dere']['rel'] = 'R1'
    assembly_dict['Dere']['gbacc'] = 'GCA_000005135.1'

    assembly_dict['Dgri'] = {}
    assembly_dict['Dgri']['rel'] = 'R1'
    assembly_dict['Dgri']['gbacc'] = 'GCA_000005155.1'

    assembly_dict['Dmoj'] = {}
    assembly_dict['Dmoj']['rel'] = 'R1'
    assembly_dict['Dmoj']['gbacc'] = 'GCA_000005175.1'

    assembly_dict['Dper'] = {}
    assembly_dict['Dper']['rel'] = 'R1'
    assembly_dict['Dper']['gbacc'] = 'GCA_000005195.1'

    assembly_dict['Dpse'] = {}
    assembly_dict['Dpse']['rel'] = 'R3'
    assembly_dict['Dpse']['gbacc'] = 'GCA_000001765.2'

    assembly_dict['Dsec'] = {}
    assembly_dict['Dsec']['rel'] = 'R1'
    assembly_dict['Dsec']['gbacc'] = 'GCA_000005215.1'

    assembly_dict['Dsim'] = {}
    assembly_dict['Dsim']['rel'] = 'R2'
    assembly_dict['Dsim']['gbacc'] = 'GCA_000754195.3'

    assembly_dict['Dvir'] = {}
    assembly_dict['Dvir']['rel'] = 'R1'
    assembly_dict['Dvir']['gbacc'] = 'GCA_000005245.1'

    assembly_dict['Dwil'] = {}
    assembly_dict['Dwil']['rel'] = 'R1'
    assembly_dict['Dwil']['gbacc'] = 'GCA_000005925.1'

    assembly_dict['Dyak'] = {}
    assembly_dict['Dyak']['rel'] = 'R1'
    assembly_dict['Dyak']['gbacc'] = 'GCA_000005975.1'


def get_synonyms(fid, stype):
    """Get various types of synonym for a given feature_id."""
    snout = list()

    if (stype == 'synonym'):
        get_synonym = ('SELECT distinct(s.name) '
                       'FROM feature_synonym fs, synonym s, cvterm cvt '
                       'WHERE fs.feature_id = %s '
                       'AND fs.is_internal = \'f\' '
                       'AND fs.is_current = \'f\' '
                       'AND fs.synonym_id = s.synonym_id '
                       'AND s.type_id = cvt.cvterm_id '
                       'AND cvt.name = \'symbol\' ')

    elif (stype == 'fullname'):
        get_synonym = ('SELECT distinct(s.name) '
                       'FROM feature_synonym fs, synonym s, cvterm cvt '
                       'WHERE fs.feature_id = %s '
                       'AND fs.is_internal = \'f\' '
                       'AND fs.is_current = \'t\' '
                       'AND fs.synonym_id = s.synonym_id '
                       'AND s.type_id = cvt.cvterm_id '
                       'AND cvt.name = \'fullname\' ')

    synonyms = connect(get_synonym, fid, conn)
    for syn in synonyms:
        log.debug('Synonym ({}): {}'.format(stype, syn[0]))
        snout.extend((syn[0], ))

    if snout:
        return snout


def get_taxonid(oid):
    # Get NCBI TaxonID for an organism, given organism_id
    get_txid = ('SELECT accession '
                'FROM organism_dbxref od, dbxref dx, db '
                'WHERE od.organism_id = %s '
                'AND od.dbxref_id = dx.dbxref_id '
                'AND dx.db_id = db.db_id '
                'AND db.name = \'NCBITaxon\' ')
    txids = connect(get_txid, oid, conn)
    for txid in txids:
        log.debug('taxonID: {}'.format(txid[0]))
        return txid[0]


def get_soid(cvid, fbtr, sofix_dict, tsymb):
    # Get SO ID given various inputs  (See JIRA DB-479 for gory details.)
    # log.debug("\t\t\tIn get_soid: cvid is:\t", cvid, "\t(", type(cvid), ")\ttx symbol is:\t", tsymb[0])

    # Handle scaRNA (we detect these using Tx symbol -- ugh!)
    s = re.search('scaRNA', tsymb[0])
    if s:
        log.debug('\t\t\tSetting soTermID (scaRNA): SO:0002095')
        return 'SO:0002095'

    # Handle RNase_P_RNA and RNase_MRP_RNA (we detect these using Tx symbol -- ugh again!)
    p = re.search('RNaseP', tsymb[0])
    if p:
        log.debug('\t\t\tSetting soTermID (RNaseP): SO:0000386')
        return 'SO:0000386'
    m = re.search('RNaseM', tsymb[0])
    if m:
        log.debug('\t\t\tSetting soTermID (RNaseM): SO:0000385')
        return 'SO:0000385'

    # No longer need to hardcode this since pre_miRNA is now a proper SO term.
    # # Handle pre_miRNA cases (cvterm 'pre_miRNA' is missing SO term in chado, evidently)
    # if cvid[0] == 97396:
    #     log.debug('\t\t\t\tsetting soTermId (pre_miRNA):\tSO:0001244')
    #     return('SO:0001244')

    # Handle cases that get hard-fixed using data in sofix_dict
    if fbtr[0] in sofix_dict:
        log.debug('Setting soTermId (from sofix_dict): {}'.format(sofix_dict[fbtr[0]]))
        return sofix_dict[fbtr[0]]

    # The common case -- simple SO dbxref lookup for cvterm linked from feature.type_id
    get_sid = ('SELECT db.name||\':\'||dbx.accession '
               'FROM cvterm cvt '
               'JOIN dbxref dbx ON dbx.dbxref_id = cvt.dbxref_id '
               'JOIN db ON db.db_id = dbx.db_id '
               'WHERE cvt.cvterm_id = %s and '
               'db.name = \'SO\' ')
    sids = connect(get_sid, cvid, conn)
    if sids:
        if sids[0][0] == 'SO:0000655':
            # Handle SRP_RNA_gene, lncRNA and antisense_lncRNA cases (See JIRA DB-479 for explanation)
            log.debug("\t\t\tChecking SO:0000655...")
            get_alr = ('SELECT cvt.name '
                       'FROM feature_relationship fr, feature t, feature f, cvterm gt, feature_cvterm fc, cvterm cvt '
                       'WHERE t.uniquename = %s '
                       'AND t.feature_id = fr.subject_id '
                       'AND fr.object_id = f.feature_id '
                       'AND f.is_obsolete = \'f\' '
                       'AND f.type_id = gt.cvterm_id '
                       'AND gt.name = \'gene\' '
                       'AND f.feature_id = fc.feature_id '
                       'AND fc.cvterm_id = cvt.cvterm_id '
                       'AND cvt.name in (\'antisense_lncRNA_gene\', \'SRP_RNA_gene\', \'lncRNA_gene\') ')
            alrs = connect(get_alr, fbtr, conn)
            if alrs:
                for alr in alrs:
                    if alr[0] == 'lncRNA_gene':
                        log.debug('\t\t\t\tSetting soTermId (lncRNA):\tSO:0001877')
                        return 'SO:0001877'
                    elif alr[0] == 'antisense_lncRNA_gene':
                        log.debug('\t\t\t\tSetting soTermId (antisense_lncRNA):\tSO:0001904')
                        return 'SO:0001904'
                    else:
                        log.debug('\t\t\t\tSetting soTermId (SRP_RNA):\tSO:0000590')
                        return 'SO:0000590'
            else:
                log.debug('\t\t\t\tSetting soTermId (ncRNA):\tSO:0000655')
                return 'SO:0000655'
        else:
            log.debug('\t\t\tsoTermId: {}'.format(sids[0][0]))
            return sids[0][0]
    else:
        log.debug("\t\t\tsoTermID: WARNING - NONE FOUND")


def pop_sofix_dict(sofix_dict):
    # Hard-coded dict fixing cases where wrong SOIDs are associated w/ Tx in chado (blech)
    sofix_dict['FBtr0346873'] = 'SO:0000209'
    sofix_dict['FBtr0346877'] = 'SO:0000209'
    sofix_dict['FBtr0346881'] = 'SO:0000209'
    sofix_dict['FBtr0346899'] = 'SO:0000209'


def get_xrefs(tid):
    # Get REFSEQ & MIR dbxrefs given a feature_id
    xout = []

    get_dbxrefs = ('SELECT db.name, accession, version '
                   'FROM feature_dbxref fd, dbxref dx, db '
                   'WHERE fd.feature_id = %s '
                   'AND fd.dbxref_id = dx.dbxref_id '
                   'AND fd.is_current = \'t\' '
                   'AND dx.db_id = db.db_id '
                   'AND db.name in (\'REFSEQ\', \'MIR\')')
    dxrefs = connect(get_dbxrefs, tid, conn)
    for dxref in dxrefs:
        log.debug('\t\t\tcrossReferenceIDs: {}:{}.{}'.format(dxref[0], dxref[1], dxref[2]))

        if dxref[2]:
            accno = dxref[1] + '.' + dxref[2]
        else:
            accno = dxref[1]

        if dxref[0] == 'MIR':    # 'MIRBASE' is specified in the RNAcentral json schema
            accout = 'MIRBASE:' + accno
            xout.append(accout)
        else:
            accout = 'REFSEQ:' + accno
            xout.append(accout)

    return xout


def get_relseqs(tid):
    # Get related sequences: for now, set up only for related precursor and matureProduct of miRNAs
    routs = []

    # Get precursor miRNA
    get_precseqs = ('SELECT ft.name, f.uniquename, f.name '
                    'FROM feature f, feature_relationship fr, cvterm frt, cvterm ft '
                    'WHERE fr.type_id = frt.cvterm_id '
                    'AND frt.name = \'producedby\' '
                    'AND object_id = f.feature_id '
                    'AND f.is_obsolete = \'f\' '
                    'AND f.type_id = ft.cvterm_id '
                    'AND fr.subject_id = %s')
    prseqs = connect(get_precseqs, tid, conn)
    for prseq in prseqs:
        rrec = {}
        log.debug('\t\t\trelated sequence (precursor):\t{}\t{}\t{}'.format(prseq[0], prseq[1], prseq[2]))
        rrec['sequenceId'] = 'FLYBASE:' + prseq[1]
        rrec['relationship'] = 'precursor'
        routs.append(rrec)

    # Get matureProduct miRNA(s)
    get_prodseqs = ('SELECT ft.name, f.uniquename, f.name '
                    'FROM feature f, feature_relationship fr, cvterm frt, cvterm ft '
                    'WHERE fr.type_id = frt.cvterm_id '
                    'AND frt.name = \'producedby\' '
                    'AND subject_id = f.feature_id '
                    'AND f.is_obsolete = \'f\' '
                    'AND f.type_id = ft.cvterm_id '
                    'AND fr.object_id = %s')
    prdseqs = connect(get_prodseqs, tid, conn)
    for prdseq in prdseqs:
        rrec = {}
        log.debug('\t\t\trelated sequence (matureProduct):\t{}\t{}\t{}'.format(prdseq[0], prdseq[1], prdseq[2]))
        rrec['sequenceId'] = 'FLYBASE:' + prdseq[1]
        rrec['relationship'] = 'matureProduct'
        routs.append(rrec)

    return routs


def get_rep_pubs(gene_uniquename):
    """Get representative publications for a transcript's parent gene."""
    log.debug('Get represenative pubs for gene {}'.format(gene_uniquename))
    rep_pub_query = """
        SELECT DISTINCT dbx.accession
        FROM feature f
        JOIN feature_pub fp ON fp.feature_id = f.feature_id
        JOIN feature_pubprop fpp ON fpp.feature_pub_id = fp.feature_pub_id
        JOIN cvterm cvt ON cvt.cvterm_id = fpp.type_id
        JOIN pub_dbxref pdbx ON pdbx.pub_id = fp.pub_id
        JOIN dbxref dbx ON dbx.dbxref_id = pdbx.dbxref_id
        JOIN db ON db.db_id = dbx.db_id
        WHERE db.name = 'pubmed'
          AND cvt.name = 'computed_gene_pub_score'
          AND f.uniquename = %s;
    """
    rep_pub_results = connect(rep_pub_query, (gene_uniquename, ), conn)
    if rep_pub_results:
        rep_pubs = ['PMID:{}'.format(i[0]) for i in rep_pub_results]
        log.debug('Found {} represenative pubs for gene {}'.format(len(rep_pubs), gene_uniquename))
    else:
        rep_pubs = None
        log.debug('NO PUBS for {}'.format(gene_uniquename))
    return rep_pubs


def get_annoid(gid):
    # Get annotation ID for a gene, given feature_id
    get_annid = ('SELECT db.name, accession, is_current '
                 'FROM feature_dbxref fd, dbxref dx, db '
                 'WHERE fd.is_current = \'t\' '
                 'AND fd.dbxref_id = dx.dbxref_id '
                 'AND dx.db_id = db.db_id '
                 'AND db.name = \'FlyBase Annotation IDs\' '
                 'AND fd.feature_id = %s')
    annids = connect(get_annid, gid, conn)
    for annid in annids:
        log.debug('\t\t\tannotation ID: {}'.format(annid[1]))
        return annid[1]


def get_mirnatx(tid):
    # Get info on miRNA transcript(s) (associated w/ -RM transcripts)
    mirecs = []

    get_mitxdata = ('SELECT DISTINCT cvt.name, t.feature_id, t.uniquename, '
                    '       t.name, t.organism_id, cvt.cvterm_id, t.residues '
                    'FROM feature_relationship fr, feature t, cvterm cvt, cv '
                    'WHERE fr.object_id = %s '
                    'AND fr.subject_id = t.feature_id '
                    'AND t.type_id = cvt.cvterm_id '
                    'AND cvt.name = \'miRNA\' '
                    'AND cvt.cv_id = cv.cv_id '
                    'AND cv.name = \'SO\'')
    mitxs = connect(get_mitxdata, tid, conn)
    for mitx in mitxs:
        log.debug('\t\tmiRNA:\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.
                  format(mitx[0], mitx[1], mitx[2], mitx[3], mitx[4], mitx[5], mitx[6]))
        mirecs.append(mitx)

    return mirecs


def get_glocinfo(tid):
    # Get location information
    ex_dict = {}
    ex_list = []

    str_dict = {1: '+', -1: '-'}

    get_glocex = ('SELECT f.feature_id, f.uniquename, f.name, fmin, fmax, '
                  '       strand, srcfeature_id, s.name, accession, version '
                  'FROM feature f, feature_relationship fr, feature s, featureloc fl, '
                  '     feature_dbxref fd, dbxref dx, db, cvterm cvt, cvterm cvt2 '
                  'WHERE fr.object_id = %s '
                  'AND fr.subject_id = f.feature_id '
                  'AND f.is_obsolete = \'f\' '
                  'AND f.type_id = cvt2.cvterm_id '
                  'AND cvt2.name = \'exon\' '
                  'AND fr.type_id = cvt.cvterm_id '
                  'AND cvt.name = \'partof\' '
                  'AND f.feature_id = fl.feature_id '
                  'AND fl.srcfeature_id = s.feature_id '
                  'AND s.feature_id = fd.feature_id '
                  'AND fd.is_current = \'t\' '
                  'AND fd.dbxref_id = dx.dbxref_id '
                  'AND dx.db_id = db.db_id '
                  'AND db.name = \'GB\' ')
    UNAME = 1
    NAME = 2
    FMIN = 3
    FMAX = 4
    STRAND = 5
    SRC_NAME = 7
    SRC_ACC = 8
    VERSION = 9
    exes = connect(get_glocex, tid, conn)
    if exes:    # The case for all except miRNAs
        for ex in exes:
            log.debug('\t\t\texon:\t{}\t{}\tARM:\t{}\t{}:{}..{} ({})'.
                      format(ex[NAME], ex[UNAME], ex[SRC_ACC], ex[SRC_NAME], ex[FMIN], ex[FMAX], ex[STRAND]))
            ex_dict[ex[UNAME]] = {}
            ex_dict[ex[UNAME]]['INSDC_accession'] = ex[SRC_ACC] + '.' + ex[VERSION]
            ex_dict[ex[UNAME]]['chromosome'] = ex[SRC_NAME]
            ex_dict[ex[UNAME]]['strand'] = str_dict[ex[STRAND]]
            ex_dict[ex[UNAME]]['startPosition'] = ex[FMIN] + 1    # Convert to 1-base
            ex_dict[ex[UNAME]]['endPosition'] = ex[FMAX]

            ex_list.append(ex_dict[ex[UNAME]])
    else:    # For miRNAs (that dont have exons like others)
        get_mirex = ('SELECT f.feature_id, f.uniquename, f.name, fmin, '
                     '       fmax, strand, srcfeature_id, s.name, accession, version '
                     'FROM feature f, feature s, featureloc fl, feature_dbxref fd, dbxref dx, db '
                     'WHERE f.feature_id = %s '
                     'AND f.is_obsolete = \'f\' '
                     'AND f.feature_id = fl.feature_id '
                     'AND fl.srcfeature_id = s.feature_id '
                     'AND s.feature_id = fd.feature_id '
                     'AND fd.is_current = \'t\' '
                     'AND fd.dbxref_id = dx.dbxref_id '
                     'AND dx.db_id = db.db_id '
                     'AND db.name = \'GB\' ')
        mixes = connect(get_mirex, tid, conn)
        for ex in mixes:
            log.debug('\t\t\texon:\t{}\t{}\tARM:\t{}\t{}:{}..{} ({})'.
                      format(ex[NAME], ex[UNAME], ex[SRC_ACC], ex[SRC_NAME], ex[FMIN], ex[FMAX], ex[STRAND]))
            ex_dict[ex[UNAME]] = {}
            ex_dict[ex[UNAME]]['INSDC_accession'] = ex[SRC_ACC] + '_' + ex[VERSION]
            ex_dict[ex[UNAME]]['chromosome'] = ex[SRC_NAME]
            ex_dict[ex[UNAME]]['strand'] = str_dict[ex[STRAND]]
            ex_dict[ex[UNAME]]['startPosition'] = ex[FMIN] + 1    # Convert to 1-base
            ex_dict[ex[UNAME]]['endPosition'] = ex[FMAX]

            ex_list.append(ex_dict[ex[UNAME]])

    return ex_list


def get_ncrna_json(database, ncrna_dict):
    # Build JSON records for each ncRNA to be reported
    # Instantiate assembly_dict, containing GB acc# and release# for species genome assemblies (not in chado)
    assembly_dict = {}
    pop_assembly_dict(assembly_dict)

    # Instantiate sofix_dict, containing hard fix of SO terms to use for some FBtr (yeah ick)
    sofix_dict = {}
    pop_sofix_dict(sofix_dict)

    # Main driver query for ncRNA genes
    get_ncgenes = """
        SELECT DISTINCT g.feature_id,
                        g.uniquename,
                        g.name,
                        srcfeature_id,
                        o.abbreviation,
                        cvtt.name,
                        t.feature_id,
                        t.uniquename,
                        t.name,
                        t.organism_id,
                        cvtt.cvterm_id,
                        t.residues
        FROM feature t
        JOIN feature_relationship fr ON t.feature_id = fr.subject_id
        JOIN feature g ON fr.object_id = g.feature_id
        JOIN cvterm cvtt ON (t.type_id = cvtt.cvterm_id AND NOT cvtt.name in ('mRNA', 'pseudogene'))
        JOIN cvterm cvtfr ON fr.type_id = cvtfr.cvterm_id
        JOIN featureloc fl ON t.feature_id = fl.feature_id
        JOIN organism o ON g.organism_id = o.organism_id
        WHERE o.abbreviation = 'Dmel'
          AND cvtfr.name = 'partof'
          AND t.is_obsolete IS FALSE
          AND t.uniquename ~ '^FBtr[0-9]{7}$'
          AND g.is_obsolete IS FALSE
          AND g.uniquename ~ '^FBgn[0-9]{7}$';
    """
    ncgenes = connect(get_ncgenes, 'no_query', conn)
    log.info('Found {} ncRNA genes.'.format(len(ncgenes)))
    GENE_UNAME = 1
    GENE_NAME = 2
    TX_TYPE = 5
    TX_FEAT_ID = 6
    TX_UNAME = 7
    TX_NAME = 8
    for ncgene in ncgenes:
        log.debug('PROCESSING: GENE: {} ({})\tTX_TYPE: {}\tTX: {} ({})'.
                  format(ncgene[GENE_NAME], ncgene[GENE_UNAME], ncgene[TX_TYPE], ncgene[TX_NAME], ncgene[TX_UNAME]))
        record_dict = {}
        ncrg = []

        # If we got an -RM transcript, this is a possibly a pre_miRNA.
        # We look for an associated transcript which is the one we want for the report.
        if ncgene[TX_NAME].endswith('-RM'):
            mixl = []
            log.debug('STEP1. Before line 520 query.')
            mixl = get_mirnatx((ncgene[TX_FEAT_ID], ))
            log.debug('STEP2. After line 520 query.')
            # Produce report for the mature miRNAs
            if mixl:
                for mid in mixl:
                    ncrg = []
                    record_dict = {}
                    for i in range(5):
                        ncrg.append(ncgene[i])
                    ncrg.extend(mid)
                    log.debug('Calling pop_json_record: {}'.format(ncrg[TX_UNAME]))
                    pop_json_record(database, record_dict, ncrg, assembly_dict, sofix_dict)
                    log.debug('STEP3. After pop_json_record call.')
                    ncrna_dict['data'].append(record_dict)
            # Produce report for the precursor miRNA
            ncrg = ncgene
            record_dict = {}
            log.debug('Calling pop_json_record (RM): {}'.format(ncrg[TX_UNAME]))
            pop_json_record(database, record_dict, ncrg, assembly_dict, sofix_dict)
            ncrna_dict['data'].append(record_dict)
        else:
            ncrg = ncgene
            log.debug('Calling pop_json_record: {}'.format(ncrg[TX_UNAME]))
            pop_json_record(database, record_dict, ncrg, assembly_dict, sofix_dict)
            ncrna_dict['data'].append(record_dict)


def pop_json_record(database, record_dict, ncrg, assembly_dict, sofix_dict):
    """Initialize & start populating record_dict (JSON object for each ncRNA) w/ data at hand."""
    GENE_FEAT_ID = 0
    GENE_UNAME = 1
    GENE_NAME = 2
    GENE_ORG = 4
    TX_FEAT_ID = 6
    TX_UNAME = 7
    TX_NAME = 8
    TX_ORG = 9
    TX_TYPE_ID = 10
    TX_SEQ = 11
    record_dict['primaryId'] = 'FLYBASE:' + ncrg[TX_UNAME]
    record_dict['symbol'] = ncrg[TX_NAME]
    record_dict['sequence'] = ncrg[TX_SEQ].upper()
    record_dict['url'] = 'http://flybase.org/reports/' + ncrg[TX_UNAME] + '.html'

    # Get synonym(s) associated w/ FBtr
    tsyns = []
    tsyns = get_synonyms((ncrg[TX_FEAT_ID], ), 'synonym')
    if tsyns:
        record_dict['symbolSynonyms'] = tsyns

    # Get taxonID associated w/ FBtr
    ntxid = get_taxonid((ncrg[TX_ORG], ))
    record_dict['taxonId'] = 'NCBITaxon:' + ntxid

    # Get soTermId associated w/ FBtr
    stid = get_soid((ncrg[TX_TYPE_ID], ), (ncrg[TX_UNAME], ), sofix_dict, (ncrg[TX_NAME], ))
    record_dict['soTermId'] = stid

    # Get crossReferenceIds associated w/ FBtr
    xrefs = []
    xrefs = get_xrefs((ncrg[TX_FEAT_ID], ))
    if xrefs:
        record_dict['crossReferenceIds'] = xrefs

        # Get related sequences associated w/ FBtr
        relseqs = []
        relseqs = get_relseqs((ncrg[TX_FEAT_ID], ))
        if relseqs:
            record_dict['relatedSequences'] = relseqs

    # Get genomeLocation assembly info associated w/ FBtr
    record_dict['genomeLocations'] = []
    gloc = {}
    gloc['assembly'] = assembly_dict[ncrg[GENE_ORG]]['rel']
    gloc['gca_accession'] = assembly_dict[ncrg[GENE_ORG]]['gbacc']

    # Get genomeLocation location info associated w/ FBtr
    gloc['exons'] = []
    gloc['exons'] = get_glocinfo((ncrg[TX_FEAT_ID], ))
    record_dict['genomeLocations'].append(gloc)

    # Get info associated with parent gene (data at hand)
    record_dict['gene'] = {}
    record_dict['gene']['geneId'] = 'FLYBASE:' + ncrg[GENE_UNAME]
    record_dict['gene']['symbol'] = ncrg[GENE_NAME]
    record_dict['gene']['url'] = 'http://flybase.org/reports/' + ncrg[GENE_UNAME] + '.html'

    # Get representative pubs associated with the parent gene.
    rep_pubs = get_rep_pubs(ncrg[GENE_UNAME])
    if rep_pubs:
        record_dict['publications'] = rep_pubs

    # Get the annotation ID associated w/ the parent gene and combine w/ abbreviation for locusTag
    # See JIRA DB-357 for rationale of Dsim -> Dsimw501 fudge
    anoid = get_annoid((ncrg[GENE_FEAT_ID], ))
    if ncrg[GENE_ORG] == 'Dsim':
        record_dict['gene']['locusTag'] = "Dsimw501_" + anoid
    else:
        record_dict['gene']['locusTag'] = ncrg[GENE_ORG] + "_" + anoid

    # Get synonym(s) associated w/ parent gene
    gsyns = []
    gsyns = get_synonyms((ncrg[GENE_FEAT_ID], ), 'synonym')
    if gsyns:
        record_dict['gene']['synonyms'] = gsyns

    # Get fullname associated w/ parent gene (if any)
    fname = []
    fname = get_synonyms((ncrg[GENE_FEAT_ID], ), 'fullname')
    if fname:
        record_dict['gene']['name'] = fname[0]


if __name__ == "__main__":
    main()
