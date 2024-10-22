#!/usr/bin/env python3

# report_fu_gal4_table.py
# author: David Emmert, modified by Gil dos Santos
# usage: report_tsv_template.py [-h] [-v VERBOSE]
# notes: gets info on frequently used GAL4 drivers and writes JSON for tabular output at IU
################################################################################

import argparse
import logging
import strict_rfc3339
import psycopg2
import re
import json
import sys
import os


# Function for db queries.
def connect(sql, query, conn):
    # Return the cursor and use it to perform queries.
    cursor = conn.cursor()
    # Execute the query, with or without some variable.
    if query == 'no_query':
        cursor.execute(sql)
    else:
        cursor.execute(sql, query)
    # Grab the results.
    records = cursor.fetchall()
    # Close the cursor
    cursor.close()
    # Return a list of tuples
    return records


def get_synonym_sgml(gid, conn):
    # Get sgml-ized synonym for a given feature_id
    get_sgml_synonym = ('SELECT synonym_sgml '
                        'FROM feature_synonym fs, synonym s, cvterm cvt '
                        'WHERE fs.feature_id = %s '
                        'AND fs.is_internal = \'f\' '
                        'AND fs.is_current = \'t\' '
                        'AND fs.synonym_id = s.synonym_id '
                        'AND s.type_id = cvt.cvterm_id '
                        'AND cvt.name = \'symbol\' ')
    sgml_synonyms = connect(get_sgml_synonym,gid,conn)
    for ssyn in sgml_synonyms:
        logging.debug('\t\tSGML_Synonym:\t%s' % (ssyn[0]))
    return(ssyn[0])


# Get image metadata from local file.
def get_image_metadata():
    image_file = open('./harvdev-reports/image_metadata.txt','r')
    image_dict = {}
    line_counter = 0
    accepted_lines_counter = 0
    ignored_lines_counter = 0
    rejected_lines_counter = 0
    files_processed = []
    for line in image_file:
        line_counter += 1
        if re.search(r'^#',line):
            logging.debug('IMAGES: Ignoring comment line: %s' % (line.rstrip()))
            ignored_lines_counter += 1
            continue
        else:
            n_fields = len(line.split('\t'))
            if n_fields != 5:
                logging.warning('IMAGES: Line %d in the input file has %d fields, not the 5 expected. Line was not parsed.' % (line_counter, n_fields))
                logging.warning(line)
                rejected_lines_counter += 1
                continue
            else:
                logging.debug('IMAGES: Line %d has 5 fields and so it was parsed.' % (line_counter)) 
                key_checks_failed_counter = 0
                this_image = {}
                this_image['imageDescription'] = line.split('\t')[0]
                this_filename = line.split('\t')[1]
                this_allele = line.split('\t')[1].split('_')[0]
                this_image['publicationId'] = line.split('\t')[2]
                this_image['pubFigure'] = line.split('\t')[3]
                this_image['permission'] = line.split('\t')[4].rstrip()
                # REJECT cases where...
                # 1. Reject if imageFilename has unexpected file extension. Note that *.tif/*.tiff files are not allowed.
                if not re.search(r'^fbal[0-9]{7}(_|-)[0-9]{1,2}(\.jpg|\.jpeg|\.gif|\.png)$', this_filename.lower()):
                    logging.warning(f'IMAGES: Line {line_counter} "imageFileName" has unexpected file extension: {this_filename}.')
                    key_checks_failed_counter += 1
                # 2. Reject if publicationId is specified (not empty string) but is not an FBrf ID.
                if this_image['publicationId'] and not re.search(r'^FBrf[0-9]{7}$', this_image['publicationId']):
                    logging.warning(f'IMAGES: Line {line_counter} "publicationId" has unexpected pub ID: {this_image["publicationId"]}')
                    key_checks_failed_counter += 1
                # 3. Reject if permission not granted.
                permission_granted = False
                permission_keywords = ['open', 'grant', 'provid', 'give', 'vfb', 'virtual fly brain', 'virtualflybrain']
                for keyword in permission_keywords:
                    if keyword in this_image['permission'].lower():
                        permission_granted = True
                if permission_granted is False:
                    logging.warning(f'IMAGES: Line {line_counter} "permission" is not correctly indicated: {this_image["permission"]}')
                    key_checks_failed_counter += 1
                # FLAG (but keep) cases where ...
                # 1. Flag if imageDescription appears many times.
                if this_image['imageDescription'] in files_processed:
                    logging.warning(f'IMAGES: Line {line_counter} lists "imageFileName" {this_image["imageDescription"]} again.')
                else:
                    files_processed.append(this_image['imageDescription'])
                # 2. Flag if pubFigure does not contain "Fig" string.
                if 'fig' not in this_image['pubFigure'].lower() and this_image['pubFigure'] != 'Virtual Fly Brain':
                    logging.warning(f'IMAGES: Line {line_counter} "pubFigure" does not mention any "fig": {this_image["pubFigure"]}')
                # FINAL assessment of image metadata.
                if key_checks_failed_counter > 0:
                    logging.warning(f'IMAGES: Line {line_counter} was NOT processed: failed {key_checks_failed_counter} checks.')
                    rejected_lines_counter += 1
                else:
                    logging.debug(f'IMAGES: Line {line_counter} was processed.')
                    logging.debug(this_image)
                    accepted_lines_counter += 1
                    if this_allele in image_dict.keys():
                        image_dict[this_allele][this_filename] = this_image
                    else:
                        image_dict[this_allele] = {}
                        image_dict[this_allele][this_filename] = this_image
    logging.info('IMAGES: SUMMARY: Accepted metadata for %d images.' % (accepted_lines_counter))
    logging.info('IMAGES: SUMMARY: Rejected metadata for %d images.' % (rejected_lines_counter))
    logging.info('IMAGES: SUMMARY: Ignored %d lines.' % (ignored_lines_counter))
    for element in image_dict:
        logging.debug(element)
    return image_dict


def fb_repchar(input_string, input_format, output_format):

    valid_input_format = ('proforma', 'chado_sgml', 'html', 'text_file')
    valid_output_format = ('ascii', 'proforma', 'chado_sgml', 'html', 'text_file')

    char_dicts = {
        'alpha': {'ascii': 'alpha', 'proforma': '&agr;', 'chado_sgml': '\u03B1', 'html': '\u03B1', 'text_file': '\u03B1'},
        'Alpha': {'ascii': 'Alpha', 'proforma': '&Agr;', 'chado_sgml': '\u0391', 'html': '\u0391', 'text_file': '\u0391'},
        'beta': {'ascii': 'beta', 'proforma': '&bgr;', 'chado_sgml': '\u03B2', 'html': '\u03B2', 'text_file': '\u03B2'},
        'Beta': {'ascii': 'Beta', 'proforma': '&Bgr;', 'chado_sgml': '\u0392', 'html': '\u0392', 'text_file': '\u0392'},
        'gamma': {'ascii': 'gamma', 'proforma': '&ggr;', 'chado_sgml': '\u03B3', 'html': '\u03B3', 'text_file': '\u03B3'},
        'Gamma': {'ascii': 'Gamma', 'proforma': '&Ggr;', 'chado_sgml': '\u0393', 'html': '\u0393', 'text_file': '\u0393'},
        'delta': {'ascii': 'delta', 'proforma': '&dgr;', 'chado_sgml': '\u03B4', 'html': '\u03B4', 'text_file': '\u03B4'},
        'Delta': {'ascii': 'Delta', 'proforma': '&Dgr;', 'chado_sgml': '\u0394', 'html': '\u0394', 'text_file': '\u0394'},
        'epsilon': {'ascii': 'epsilon', 'proforma': '&egr;', 'chado_sgml': '\u03B5', 'html': '\u03B5', 'text_file': '\u03B5'},
        'Epsilon': {'ascii': 'Epsilon', 'proforma': '&Egr;', 'chado_sgml': '\u0395', 'html': '\u0395', 'text_file': '\u0395'},
        'zeta': {'ascii': 'zeta', 'proforma': '&zgr;', 'chado_sgml': '\u03B6', 'html': '\u03B6', 'text_file': '\u03B6'},
        'Zeta': {'ascii': 'Zeta', 'proforma': '&Zgr;', 'chado_sgml': '\u0396', 'html': '\u0396', 'text_file': '\u0396'},
        'eta': {'ascii': 'eta', 'proforma': '&eegr;', 'chado_sgml': '\u03B7', 'html': '\u03B7', 'text_file': '\u03B7'},
        'Eta': {'ascii': 'Eta', 'proforma': '&EEgr;', 'chado_sgml': '\u0397', 'html': '\u0397', 'text_file': '\u0397'},
        'theta': {'ascii': 'theta', 'proforma': '&thgr;', 'chado_sgml': '\u03B8', 'html': '\u03B8', 'text_file': '\u03B8'},
        'Theta': {'ascii': 'Theta', 'proforma': '&THgr;', 'chado_sgml': '\u0398', 'html': '\u0398', 'text_file': '\u0398'},
        'iota': {'ascii': 'iota', 'proforma': '&igr;', 'chado_sgml': '\u03B9', 'html': '\u03B9', 'text_file': '\u03B9'},
        'Iota': {'ascii': 'Iota', 'proforma': '&Igr;', 'chado_sgml': '\u0399', 'html': '\u0399', 'text_file': '\u0399'},
        'kappa': {'ascii': 'kappa', 'proforma': '&kgr;', 'chado_sgml': '\u03BA', 'html': '\u03BA', 'text_file': '\u03BA'},
        'Kappa': {'ascii': 'Kappa', 'proforma': '&Kgr;', 'chado_sgml': '\u039A', 'html': '\u039A', 'text_file': '\u039A'},
        'lambda': {'ascii': 'lambda', 'proforma': '&lgr;', 'chado_sgml': '\u03BB', 'html': '\u03BB', 'text_file': '\u03BB'},
        'Lambda': {'ascii': 'Lambda', 'proforma': '&Lgr;', 'chado_sgml': '\u039B', 'html': '\u039B', 'text_file': '\u039B'},
        'mu': {'ascii': 'mu', 'proforma': '&mgr;', 'chado_sgml': '\u03BC', 'html': '\u03BC', 'text_file': '\u03BC'},
        'Mu': {'ascii': 'Mu', 'proforma': '&Mgr;', 'chado_sgml': '\u039C', 'html': '\u039C', 'text_file': '\u039C'},
        'nu': {'ascii': 'nu', 'proforma': '&ngr;', 'chado_sgml': '\u03BD', 'html': '\u03BD', 'text_file': '\u03BD'},
        'Nu': {'ascii': 'Nu', 'proforma': '&Ngr;', 'chado_sgml': '\u039D', 'html': '\u039D', 'text_file': '\u039D'},
        'xi': {'ascii': 'xi', 'proforma': '&xgr;', 'chado_sgml': '\u03BE', 'html': '\u03BE', 'text_file': '\u03BE'},
        'Xi': {'ascii': 'Xi', 'proforma': '&Xgr;', 'chado_sgml': '\u039E', 'html': '\u039E', 'text_file': '\u039E'},
        'omicron': {'ascii': 'omicron', 'proforma': '&ogr;', 'chado_sgml': '\u03BF', 'html': '\u03BF', 'text_file': '\u03BF'},
        'Omicron': {'ascii': 'Omicron', 'proforma': '&Ogr;', 'chado_sgml': '\u039F', 'html': '\u039F', 'text_file': '\u039F'},
        'pi': {'ascii': 'pi', 'proforma': '&pgr;', 'chado_sgml': '\u03C0', 'html': '\u03C0', 'text_file': '\u03C0'},
        'Pi': {'ascii': 'Pi', 'proforma': '&Pgr;', 'chado_sgml': '\u03A0', 'html': '\u03A0', 'text_file': '\u03A0'},
        'rho': {'ascii': 'rho', 'proforma': '&rgr;', 'chado_sgml': '\u03C1', 'html': '\u03C1', 'text_file': '\u03C1'},
        'Rho': {'ascii': 'Rho', 'proforma': '&Rgr;', 'chado_sgml': '\u03A1', 'html': '\u03A1', 'text_file': '\u03A1'},
        'sigma': {'ascii': 'sigma', 'proforma': '&sgr;', 'chado_sgml': '\u03C3', 'html': '\u03C3', 'text_file': '\u03C3'},
        'Sigma': {'ascii': 'Sigma', 'proforma': '&Sgr;', 'chado_sgml': '\u03A3', 'html': '\u03A3', 'text_file': '\u03A3'},
        'tau': {'ascii': 'tau', 'proforma': '&tgr;', 'chado_sgml': '\u03C4', 'html': '\u03C4', 'text_file': '\u03C4'},
        'Tau': {'ascii': 'Tau', 'proforma': '&Tgr;', 'chado_sgml': '\u03A4', 'html': '\u03A4', 'text_file': '\u03A4'},
        'upsilon': {'ascii': 'upsilon', 'proforma': '&ugr;', 'chado_sgml': '\u03C5', 'html': '\u03C5', 'text_file': '\u03C5'},
        'Upsilon': {'ascii': 'Upsilon', 'proforma': '&Ugr;', 'chado_sgml': '\u03A5', 'html': '\u03A5', 'text_file': '\u03A5'},
        'phi': {'ascii': 'phi', 'proforma': '&phgr;', 'chado_sgml': '\u03C6', 'html': '\u03C6', 'text_file': '\u03C6'},
        'Phi': {'ascii': 'Phi', 'proforma': '&PHgr;', 'chado_sgml': '\u03A6', 'html': '\u03A6', 'text_file': '\u03A6'},
        'chi': {'ascii': 'chi', 'proforma': '&khgr;', 'chado_sgml': '\u03C7', 'html': '\u03C7', 'text_file': '\u03C7'},
        'Chi': {'ascii': 'Chi', 'proforma': '&KHgr;', 'chado_sgml': '\u03A7', 'html': '\u03A7', 'text_file': '\u03A7'},
        'psi': {'ascii': 'psi', 'proforma': '&psgr;', 'chado_sgml': '\u03C8', 'html': '\u03C8', 'text_file': '\u03C8'},
        'Psi': {'ascii': 'Psi', 'proforma': '&PSgr;', 'chado_sgml': '\u03A8', 'html': '\u03A8', 'text_file': '\u03A8'},
        'omega': {'ascii': 'omega', 'proforma': '&ohgr;', 'chado_sgml': '\u03C9', 'html': '\u03C9', 'text_file': '\u03C9'},
        'Omega': {'ascii': 'Omega', 'proforma': '&OHgr;', 'chado_sgml': '\u03A9', 'html': '\u03A9', 'text_file': '\u03A9'},
        # 'subscript_left': {'ascii': '[[', 'proforma': '[[', 'chado_sgml': '<down>', 'html': '<sub>', 'text_file': '[['},
        # 'subscript_right': {'ascii': ']]', 'proforma': ']]', 'chado_sgml': '</down>', 'html': '</sub>', 'text_file': ']]'},
        # 'superscript_left': {'ascii': '[', 'proforma': '[', 'chado_sgml': '<up>', 'html': '<sup>', 'text_file': '['},
        # 'superscript_right': {'ascii': ']', 'proforma': ']', 'chado_sgml': '</up>', 'html': '</sup>', 'text_file': ']'}
    }

    if input_format == 'ascii':
        logging.error('REPCHAR: invalid input_format specified: "'+str(input_format)+'". Not supported by the "fb_repchar" function.')
        logging.error('REPCHAR: Spellings of Greek characters like "pi" are too easily confused with other strings.')
        logging.error('REPCHAR: valid input formats are: "proforma", "chado_sgml", "html", "text_file".')
        raise
    elif input_format not in valid_input_format:
        logging.error('REPCHAR: invalid input_format specified: "'+str(input_format)+""'. Valid input formats are: "proforma", "chado_sgml", "html", "text_file".')
        raise
    elif output_format not in valid_output_format:
        logging.error('REPCHAR: invalid output_format specified: "'+str(output_format)+'". Valid output formats are: "ascii", "proforma", "chado_sgml", "html", "text_file".')
        raise

    for char_dict in char_dicts:
        # print(type(char_rep_dict))
        # print(char_rep_dict[input_format])
        # print(char_rep_dict[output_format])
        input_string = input_string.replace(char_dicts[char_dict][input_format], char_dicts[char_dict][output_format])
    return(input_string)


# The core function. Get data, write into a JSON and dump it.
def get_fu_gal4_json(database_host, database, username, password, annotation_release, database_release, report_name, output_filename):

    # Set this as the "official" time the file was generated.
    the_time = strict_rfc3339.now_to_rfc3339_localoffset()

    # Database connection.
    conn_string = "host=%s dbname=%s user=%s password='%s'" % (database_host, database, username, password)
    conn = psycopg2.connect(conn_string)
    logging.info('Successful database connection to database %s on host %s. Time: %s' % (database, database_host, strict_rfc3339.now_to_rfc3339_localoffset()))

    # Gather image metadata.
    image_dict = get_image_metadata()
    logging.info('Gathered image metadata. Time: %s' % (strict_rfc3339.now_to_rfc3339_localoffset()))

    ## Setup JSON
    fug4_json = {}
    fug4_json['metaData'] = {}
    fug4_json['metaData']['dataProvider'] = 'FlyBase'
    fug4_json['metaData']['title'] = 'Frequently Used GAL4 Table'
    fug4_json['metaData']['dateProduced'] = the_time
    fug4_json['metaData']['databaseRelease'] =  database_release
    fug4_json['metaData']['annotationRelease'] =  annotation_release
    fug4_json['data'] = []

    ## This gets table column A (GAL4 Driver) & drives subsequent queries
    ## Query for gene features linked to tx features where tx feature_pub = FBrf0237128
    get_fu_gal4_genes = ('SELECT g.feature_id, g.uniquename, g.name, t.feature_id, t.uniquename, t.name '
                        'FROM feature g, feature t, feature_relationship fr, cvterm cvt, feature_pub fp, pub p '
                        'WHERE p.uniquename = \'FBrf0237128\' '
                        'AND p.pub_id = fp.pub_id '
                        'AND fp.feature_id = t.feature_id '
                        'AND t.is_obsolete = \'f\' '
                        'AND t.feature_id = fr.subject_id '
                        'AND fr.type_id = cvt.cvterm_id and cvt.name = \'associated_with\' '
                        'AND fr.object_id = g.feature_id '
                        'AND g.is_obsolete = \'f\'')
    fu_genes = connect(get_fu_gal4_genes, 'no_query', conn)
    allele_ids_already_handled = []
    for gene in fu_genes:
        logging.debug('\nProcessing gene: %s\t%s\t%s\t%s\t%s\t%s' % (gene))
        
        qid = (gene[3],)
        gid = (gene[0],)
        allele_uniquename = gene[1]
        allele_name = gene[2]
        xprn_uniquename = gene[4]
        xprn_name = gene[5]

        if gid in allele_ids_already_handled:
            logging.warning(f'Skip duplicate xprn feature found for allele {allele_name} ({allele_uniquename}): {xprn_name} ({xprn_uniquename})')
            continue
        else:
            allele_ids_already_handled.append(gid)

        gsyn = get_synonym_sgml(gid,conn)

        record_dict = {}
        record_dict['driver'] = {}
        record_dict['driver']['name'] = gsyn
        record_dict['driver']['fbid'] = gene[1]
        logging.debug("Added this gene to the data list: %s" % (gene[1]))

        if gene[1] in image_dict.keys():
            record_dict['driver']['images'] = image_dict[gene[1]]
            logging.debug("Found image metadata for this gene: %s" % (gene[1]))
        else:
            record_dict['driver']['images'] = None
            logging.debug("No image metadata found for this gene: %s" % (gene[1]))

        ## Start list of ids for querying stocks
        sids = list()
        sids.extend(gid)

## Excluded from the table per Josh (See JIRA)
#        ## Get table column B (Synonyms)
#        ## MAKE SURE TO GET SGML-IZED SYNONYMS
#        ## Query using Josh's algorithm (using g.feature_id for now)
#        record_dict['driver']['synonyms'] = list()
#        get_synonyms = ('SELECT DISTINCT s.name '
#                   'FROM feature_synonym fs, synonym s, cvterm cvt '
#                   'WHERE fs.is_current = \'f\' ' 
#                   'AND fs.is_internal = \'f\' '
#                   'AND fs.synonym_id = s.synonym_id '
#                   'AND s.type_id = cvt.cvterm_id '
#                   'AND cvt.name = \'symbol\' '
#                   'AND fs.feature_id = %s')
#        synonyms = connect(get_synonyms,gid,conn)
#        for synonym in synonyms:
#            print("\tSynonym:\t",synonym[0])
#            record_dict['driver']['synonyms'].append(synonym[0])

        ## Get table column K (Pubs)
        ## Query using g.feature_id
        record_dict['driver']['pubs'] = {}
        get_pubs = ('SELECT DISTINCT p.uniquename, p.miniref '
                   'FROM feature_pub fp, pub p, cvterm cvt '
                   'WHERE fp.pub_id = p.pub_id '
                   'AND p.is_obsolete = \'f\' '
                   'AND p.type_id = cvt.cvterm_id '
                   'AND cvt.name = \'paper\' '
                   'AND fp.feature_id = %s')
        pubs = connect(get_pubs,gid,conn)
        logging.debug('\tPubs found:\t%d' % (len(pubs)))
        for pub in pubs:
#            print("\tPub:\t",pub[0],"\t",pub[1])
            record_dict['driver']['pubs'][pub[0]] = pub[1]
        
        ## Get table column C (Reflects Expression of Gene)
        ## Query using t.feature_id & f_r.type = attributed_as_expression_of
        get_rex_genes = ('SELECT DISTINCT r.feature_id, r.uniquename, r.name '
                        'FROM feature r, feature_relationship fr, cvterm cvt '
                        'WHERE fr.subject_id = %s '
                        'AND fr.type_id = cvt.cvterm_id '
                        'AND cvt.name = \'attributed_as_expression_of\' '
                        'AND fr.object_id = r.feature_id '
                        'AND r.is_obsolete = \'f\' ')
        rex_genes = connect(get_rex_genes,qid,conn)
        for rgene in rex_genes:
            logging.debug('\tREX GENE:\t%s\t%s\t%s' % (rgene[0], rgene[1], rgene[2]))
            rsyn = get_synonym_sgml((rgene[0],),conn)
            record_dict['driver']['rex_gene'] = {}
            record_dict['driver']['rex_gene'][rgene[1]] = rsyn

        ## Get table column D (Common terms used to describe expression pattern)
        ## Query using t.feature_id & feature_expressionprop.type = GAL4_table_note
        get_cterms = ('SELECT value '
                     'FROM feature_expression fe, feature_expressionprop fep, pub p, cvterm cvt '
                     'WHERE fe.feature_expression_id = fep.feature_expression_id '
                     'AND fep.type_id = cvt.cvterm_id '
                     'AND cvt.name = \'GAL4_table_note\' '
                     'AND fe.pub_id = p.pub_id '
                     'AND p.uniquename = \'FBrf0237128\' '
                     'AND fe.feature_id = %s')
        cterms = connect(get_cterms,qid,conn)
        for cterm in cterms:
            logging.debug('\tCommon terms\t%s' % (cterm[0]))
            record_dict['driver']['common_terms'] = cterm[0]

        ## Get table columns E (Major Tissue FBbt) & F (Major Stage FBdv)
        ## Query using t.feature_id & expression_cvterm (cv = FlyBase anatomy CV?) and feature_expression.pub = FBrf0237128
        get_cvterms = ('SELECT cv.name, cvt.name, db.name || accession '
                      'FROM feature_expression fe, expression_cvterm ec, cvterm cvt, cv, dbxref dx, db, pub p '
                      'WHERE fe.feature_id = %s '
                      'AND fe.pub_id = p.pub_id '
                      'AND p.uniquename = \'FBrf0237128\' '
                      'AND fe.expression_id = ec.expression_id '
                      'AND ec.cvterm_id = cvt.cvterm_id '
                      'AND cvt.cv_id = cv.cv_id '
                      'AND cvt.dbxref_id = dx.dbxref_id '
                      'AND dx.db_id = db.db_id')
        cvterms = connect(get_cvterms,qid,conn)
        record_dict['driver']['major_stages'] = {}
        record_dict['driver']['major_tissues'] = {}
        record_dict['driver']['transposons'] = {}          
        fbcv_dict = {}
        fbbt_dict = {}
        ## Handle cvterms depending on cv (developmental/FBdv, anatomy/FBbt, or flybase_controlled/FBcv)
        for cvterm in cvterms:
            logging.debug('%s\t%s\t%s' % (cvterm[0], cvterm[1], cvterm[2]))
            if cvterm[0] == 'FlyBase development CV':
                record_dict['driver']['major_stages'][cvterm[2]] = cvterm[1]
            elif cvterm[0] == 'FlyBase anatomy CV':
                fbbt_dict[cvterm[1]] = cvterm[2]
            else:
                fbcv_dict[cvterm[1]] = cvterm[2]
        ## For FBbt, construct compount statement w/ FBcv value whenever we see 'organism'  (e.g., "organism | dorsal")
        for akey in fbbt_dict:
            if akey == 'organism':
                for qkey in fbcv_dict:
                            stmt_val = "organism | " + qkey
                            stmt_ids = fbbt_dict['organism'] + " | " + fbcv_dict[qkey]
                            logging.debug('\t\tFBbt (compound):\t%s\t%s' % (stmt_val,stmt_ids))
                            record_dict['driver']['major_tissues'][stmt_ids] = stmt_val
            else:
                record_dict['driver']['major_tissues'][fbbt_dict[akey]] = akey
                logging.debug('\t\tFBbt:\t%s\t%s' % (akey, fbbt_dict[akey]))     

        ## Get table column G (Text description of GAL4 expression patterns)
        ## Query using t.feature_id & fp.type = bodypart_expression_text and fp_pub = FBrf0237128
        get_dexps = ('SELECT value '
                    'FROM featureprop fp, featureprop_pub fpp, pub p, cvterm cvt '
                    'WHERE fp.feature_id = %s '
                    'AND fp.type_id = cvt.cvterm_id '
                    'AND cvt.name = \'bodypart_expression_text\' '
                    'AND fp.featureprop_id = fpp.featureprop_id '
                    'AND fpp.pub_id = p.pub_id '
                    'AND p.uniquename = \'FBrf0237128\'')
        dexps = connect(get_dexps,qid,conn)
        for dexp in dexps:
            polished_dexp = fb_repchar(dexp[0], "proforma", "text_file")
            logging.debug('\tExp Pattern description:\t%s' % (polished_dexp))
            record_dict['driver']['expression_desc_text'] = polished_dexp

        ## Get table column I (Insertion)
        ## Query using feature_relationship.type "associated_with" from FBgn to FBti (type "transposable_element_insertion_site")
        get_tis = ('SELECT i.feature_id, i.uniquename, i.name '
                  'FROM feature_relationship fr, feature i, cvterm cvt1, cvterm cvt2 '
                  'WHERE subject_id = %s '
                  'AND fr.type_id = cvt1.cvterm_id '
                  'AND cvt1.name = \'associated_with\' '
                  'AND object_id = i.feature_id '
                  'AND i.type_id = cvt2.cvterm_id '
                  'AND i.is_obsolete = \'f\' '
                  'AND cvt2.name = \'transposable_element_insertion_site\'')
        tis = connect(get_tis,gid,conn)
        ## If theres a TI, get associated TP; If not, try to find TP via gene
        if 0 < len(tis):
            record_dict['driver']['insertions'] = {}            
            for ti in tis:
                logging.debug('\tInsertion:\t%s\t%s\t%s' % (ti[0], ti[1], ti[2]))
                tid = (ti[0],)
                sids.extend(tid)
                tisyn = get_synonym_sgml(tid,conn)
                record_dict['driver']['insertions'][ti[1]] = tisyn


        ## Get table column J (Construct) via TI
        ## Query using feature_relationship.type = producedby from FBti to FBtp (type = transgenic_transposon)
        ## Suppress reporting of base GAL4 and GawB constructs, since this is already explicit in FBti symbol.
                get_tps = ('SELECT t.feature_id, t.uniquename, t.name '
                          'FROM feature_relationship fr, feature t, cvterm cvt1, cvterm cvt2 '
                          'WHERE subject_id = %s '
                          'AND fr.type_id = cvt1.cvterm_id '
                          'AND cvt1.name = \'producedby\' '
                          'AND object_id = t.feature_id '
                          'AND t.type_id = cvt2.cvterm_id '
                          'AND t.is_obsolete = \'f\' '
                          'AND cvt2.name = \'transgenic_transposable_element\''
                          'AND NOT t.uniquename in (\'FBtp0000352\',\'FBtp0001433\')' )
                tps = connect(get_tps,tid,conn)

                for tp in tps:
                    logging.debug('\tConstruct:\t%s\t%s\t%s' % (tp[0],tp[1],tp[2]))
                    qtpid = (tp[0],)
                    tpsyn = get_synonym_sgml(qtpid,conn)
                    record_dict['driver']['transposons'][tp[1]] = tpsyn
        else: 
        ## Get table column J (Construct) via Gene
        ## Query using feature_relationship.type = associated_with from FBal to FBtp (type = transgenic_transposon)
#            print("\t\tChecking for TP via gene...")
            get_gtps = ('SELECT t.feature_id, t.uniquename, t.name '
                        'FROM feature_relationship fr, feature t, cvterm cvt1, cvterm cvt2 '
                        'WHERE subject_id = %s '
                        'AND fr.type_id = cvt1.cvterm_id '
                        'AND cvt1.name = \'associated_with\' '
                        'AND object_id = t.feature_id '
                        'AND t.type_id = cvt2.cvterm_id '
                        'AND t.is_obsolete = \'f\' '                        
                        'AND cvt2.name = \'transgenic_transposable_element\'')
            gtps = connect(get_gtps,gid,conn)
            for gtp in gtps:
                logging.debug('\tConstruct (via FBal):\t%s\t%s\t%s' % (gtp[0], gtp[1], gtp[2]))
                qtpid = (gtp[0],)
                tpsyn = get_synonym_sgml(qtpid,conn)
                record_dict['driver']['transposons'][gtp[1]] = tpsyn

        ## Get table column H (Stocks) AFTER I & J (because we use stocks linked to A+I+J for this field)
        ## We use two queries to make sure we get 'em all
        ## Query using feature_genotype, stock_genotype, stock
        stock_dict = {}
        for sid in sids:
#            print("\t\tChecking for stock associated with feature (via stock_genotype):\t",sid)
            sqid = (sid,)
            get_stocks = ('SELECT s.stock_id, s.uniquename, s.name '
                         'FROM genotype g, feature_genotype fg, stock_genotype sg, stock s '
                         'WHERE fg.feature_id = %s '
                         'AND g.is_obsolete IS FALSE '
                         'AND NOT g.uniquename ~ \'^PROBLEMATIC\' '
                         'AND g.genotype_id = fg.genotype_id '
                         'AND fg.genotype_id = sg.genotype_id '
                         'AND sg.stock_id = s.stock_id ')
            stocks = connect(get_stocks,sqid,conn)
            for stock in stocks:
                stock_dict[stock[2]] = stock[1]
#                print("\tStock (stock_genotype):\t",stock[0],"\t",stock[1],"\t",stock[2])
        for sid in sids:
#            print("\t\tChecking for stock associated with feature (via featureprop type derived_stock_*):\t",sid)
            sqid = (sid,)
            get_dstocks = ('SELECT value '
                           'FROM featureprop fp, cvterm cvt '
                           'WHERE fp.feature_id = %s '
                           'AND fp.type_id = cvt.cvterm_id '
                           'AND cvt.name in (\'derived_stock_Bloomington\',\'derived_stock_FlyORF\',\'derived_stock_Harvard\',\'derived_stock_Kyoto\',\'derived_stock_SD\',\'derived_stock_Szeged\',\'derived_stock_VDRC\') ')
            dstocks = connect(get_dstocks,sqid,conn)
            for dstock in dstocks:
#                print("\t\tStock (fprop):\t",dstock[0])
                for line in dstock[0].split('\n'):
#                    print("\t\tsplit line:\t",line)
                    stockid = re.match('^[a-zA-Z]{0,1}[0-9]+', line).group(0)   # This now accounts for the third of stock IDs beginning with a single letter.
                    fbstid = re.search("FBst[0-9]{7}",line).group(0)
#                    print("\tStock (derived_stock_*)\t",stockid,"\t",fbstid)
                    stock_dict[stockid] = fbstid
        if 0 < len(stock_dict):
            record_dict['driver']['stocks'] = {}
        for fstock in stock_dict:
                logging.debug('\tStock (stock_dict):\t%s\t%s' % (fstock,stock_dict[fstock]))
                record_dict['driver']['stocks'][stock_dict[fstock]]  = fstock

        fug4_json['data'].append(record_dict)

    with open(output_filename,'w') as outfile:
        json.dump(fug4_json, outfile, indent=2, separators=(',', ': '), ensure_ascii=False)

        outfile.close()

    conn.close()


def main():

    report_name = 'fu_gal4_table'

    parser = argparse.ArgumentParser(description='inputs')
    parser.add_argument('-v', '--verbose', action='store_true', help='Provide verbose "DEBUG"-level logging.', required=False)
    args = parser.parse_args()

    database_host = os.environ['SERVER']
    database = os.environ['DATABASE']
    username = os.environ['USER']
    password = os.environ['PGPASSWORD']
    annotation_release = os.environ['ANNOTATIONRELEASE']
    database_release = os.environ['RELEASE']

    output_filename = './output/' + report_name + '_' + database + '.json'
    log_filename = './logs/' + report_name + '_' + database + '.log'

    verbose = args.verbose
    if verbose == True:
        logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', filename=log_filename, level=logging.INFO)
    sys.stdout = open(log_filename, 'a')
    sys.stderr = open(log_filename, 'a')
    
    logging.info('Started main function. Time: %s' % (strict_rfc3339.now_to_rfc3339_localoffset()))
    get_fu_gal4_json(database_host, database, username, password, annotation_release, database_release, report_name, output_filename)
    logging.info('Ended main function. Time: %s' % (strict_rfc3339.now_to_rfc3339_localoffset()))


if __name__ == "__main__":
    main()
