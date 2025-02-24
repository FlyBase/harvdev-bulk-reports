#!/usr/local/bin/perl -w

# Script to build the GA file by querying GO data in chado producing gene_association.fb
# Usage: perl GA_file_builder host:dbname outputfile optional_log

use strict;
use DBI;
use LWP::UserAgent;
use LWP::Protocol::https;
use Data::Dumper;
use lib '../perl_modules';
use Utils;

if ( @ARGV < 5 ) {
    print "\n USAGE: perl GA_file_builder server db user password outfile (optional: log_file)\n";
    exit();
}
my $server  = shift @ARGV;
my $dbname  = shift @ARGV;
my $user    = shift @ARGV;
my $pass    = shift @ARGV;
my $outfile = shift @ARGV;
if (@ARGV) {
    my $logfile = shift @ARGV;
    open( STDOUT, ">>$logfile" ) or die print_log("Can't redirect STDOUT to $logfile");
    open( STDERR, ">>$logfile" ) or die print_log("Can't redirect STDERR to $logfile");
}
my $start = localtime();
print_log("INFO: Start GA_file_builder.");

############################## DB CONNECTION #############################
my $chado = "dbi:Pg:dbname=$dbname; host=$server;port=5432";
my $dbh   = DBI->connect( $chado, $user, $pass ) or die print_log("ERROR: Can't connect to $chado");
print_log("INFO: Connected to $chado");
##########################################################################

# Generate hash to hold the FBrf to GO ref lookup.
my %gorefs;
%gorefs = %{ fetch_and_parse_gorefs() };

# Set up output files.
open( OUT, ">>$outfile" )
  or die print_log("Can't open $outfile to write output.");
open( OUT2, ">$outfile.lines2review" )
  or die print_log("Can't open file to write report of lines to review.");
open( OUT3, ">$outfile.pubs_wo_xrefs" )
  or die print_log("Can't open file to report of pubs without dbxrefs.");

# Header info.
my $RELEASE = "$dbname";
$RELEASE = $1 if $dbname =~ /(fb_20[0-9]{2}_[0-9]{2}).*/;
$RELEASE =~ s/fb_/FB/g;
my $DATE = sprintf(
    "%04d-%02d-%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ) }
      ->(localtime)
);
my $header = '';
$header .= "!gaf-version: 2.2\n";
$header .= "!generated-by: FlyBase\n";
$header .= "!date-generated: $DATE\n";
$header .= "!saved-by: FlyBase GOcur gocur\@morgan.harvard.edu\n";
$header .= "!FlyBase release: $RELEASE\n";
print OUT $header;

# Create an FBgn-keyed hash to non-current synonyms.
# NOTE - If GO terms get moved to proteins/transcripts this will require rejiggering.
print_log("INFO: Querying for synonyms");
my %synonyms;
my $syn_query = $dbh->prepare(
    sprintf("
        SELECT DISTINCT f.uniquename, s.name
        FROM synonym s
        JOIN feature_synonym fs ON fs.synonym_id = s.synonym_id
        JOIN feature f ON f.feature_id = fs.feature_id
        WHERE f.is_obsolete IS FALSE
          AND f.is_analysis IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}\$'
          AND fs.is_current IS FALSE
    ")
);
$syn_query->execute or die print_log("ERROR: Can't fetch synonyms.");
while ( ( my $fbid, my $syn ) = $syn_query->fetchrow_array() ) {
    push @{ $synonyms{$fbid} }, $syn;
}

# Create a FBgn-keyed hash of current fullnames.
print_log("INFO: Querying for fullnames");
my %fullnames;
my $fn_query = $dbh->prepare(
    sprintf("
        SELECT DISTINCT f.uniquename, s.name
        FROM synonym s
        JOIN feature_synonym fs ON fs.synonym_id = s.synonym_id
        JOIN feature f ON f.feature_id = fs.feature_id
        JOIN cvterm cvt ON cvt.cvterm_id = s.type_id
        WHERE f.is_obsolete IS FALSE
        AND f.is_analysis IS FALSE
        AND f.uniquename ~ '^FBgn[0-9]{7}\$'
        AND fs.is_current IS TRUE
        AND cvt.name = 'fullname'
    ")
);
$fn_query->execute or die print_log("Can't fetch fullnames");
while ( ( my $fbid, my $fn ) = $fn_query->fetchrow_array() ) {
    $fullnames{$fbid} = $fn;
}

# Create a pub types hash.
my %pubsbytype;
my $stmt =
"SELECT cvterm_id, cvterm.name FROM cvterm, cv WHERE cvterm.cv_id = cv.cv_id and cv.name = 'pub type'";
my $query = $dbh->prepare($stmt);
$query->execute or die print_log("Can't do pub type query!");
while ( ( my $cid, my $name ) = $query->fetchrow_array() ) {
    $pubsbytype{$cid} = $name;
}

# This is global so it can be used in the parse_evidence_bits subroutine.
our %EVC = (
    'inferred from mutant phenotype'                    => 'IMP',
    'inferred from genetic interaction'                 => 'IGI',
    'inferred from physical interaction'                => 'IPI',
    'inferred from sequence or structural similarity'   => 'ISS',
    'inferred from sequence model'                      => 'ISM',
    'inferred from sequence alignment'                  => 'ISA',
    'inferred from sequence orthology'                  => 'ISO',
    'inferred from experiment'                          => 'EXP',
    'inferred from direct assay'                        => 'IDA',
    'inferred from electronic annotation'               => 'IEA',
    'inferred from expression pattern'                  => 'IEP',
    'inferred from reviewed computational analysis'     => 'RCA',
    'traceable author statement'                        => 'TAS',
    'non-traceable author statement'                    => 'NAS',
    'inferred by curator'                               => 'IC',
    'inferred from genomic context'                     => 'IGC',
    'no biological data available'                      => 'ND',
    'inferred from biological aspect of ancestor'       => 'IBA',
    'inferred from biological aspect of descendant'     => 'IBD',
    'inferred from key residues'                        => 'IKR',
    'inferred from rapid divergence'                    => 'IRD',
    'inferred from high throughput experiment'          => 'HTP',
    'inferred from high throughput direct assay'        => 'HDA',
    'inferred from high throughput expression pattern'  => 'HEP',
    'inferred from high throughput genetic interaction' => 'HGI',
    'inferred from high throughput mutant phenotype'    => 'HMP',
);

my %ASP = (
    'biological_process' => 'P',
    'cellular_component' => 'C',
    'molecular_function' => 'F',
);

my %TAX = ( Dmel => '7227' );

# Prepare the main gene GO annotation query.
my $ga_query = $dbh->prepare(
    sprintf("
        SELECT DISTINCT gene.feature_id,
                        gene.uniquename AS fbid,
                        symb.name AS symbol,
                        fcvt.feature_cvterm_id,
                        cv.name AS aspect,
                        godb.accession AS GO_id,
                        ppub.uniquename AS ppub,
                        o.abbreviation AS species,
                        evc.value AS evidence_code,
                        prv.value AS provenance,
                        date.value AS date,
                        fcvt.is_not AS is_not,
                        evc.rank AS evidence_code_rank
        FROM cvterm goterm
        JOIN cv ON goterm.cv_id = cv.cv_id
        JOIN dbxref godb ON goterm.dbxref_id = godb.dbxref_id
        JOIN feature_cvterm fcvt ON goterm.cvterm_id = fcvt.cvterm_id
        JOIN feature gene ON fcvt.feature_id = gene.feature_id
        JOIN organism o ON gene.organism_id = o.organism_id
        JOIN feature_synonym fs1 ON gene.feature_id = fs1.feature_id AND fs1.is_current IS TRUE
        JOIN synonym symb ON fs1.synonym_id = symb.synonym_id
        JOIN cvterm stype on symb.type_id = stype.cvterm_id and stype.name = 'symbol'
        JOIN pub ppub ON fcvt.pub_id = ppub.pub_id
        JOIN feature_cvtermprop prv ON fcvt.feature_cvterm_id = prv.feature_cvterm_id
        JOIN cvterm prvname ON prv.type_id = prvname.cvterm_id and prvname.name = 'provenance'
        JOIN feature_cvtermprop evc ON  fcvt.feature_cvterm_id = evc.feature_cvterm_id
        JOIN cvterm evcname ON evc.type_id = evcname.cvterm_id and evcname.name = 'evidence_code'
        JOIN feature_cvtermprop date ON fcvt.feature_cvterm_id = date.feature_cvterm_id
        JOIN cvterm dtname ON date.type_id = dtname.cvterm_id and dtname.name = 'date'
        WHERE gene.is_obsolete IS FALSE
          AND gene.is_analysis IS FALSE
          AND gene.uniquename ~ '^FBgn[0-9]{7}\$'
          AND o.abbreviation = 'Dmel'
          AND goterm.is_obsolete = 0
          AND cv.name in ('cellular_component','molecular_function','biological_process')
    ")
);

# A a separate query to get qualifier information to stick in post-big-query.
print_log("INFO: Getting the qualifier information.");
my %quals;
my $qual_query = $dbh->prepare(
    sprintf("
        SELECT fcvtp.feature_cvterm_id, qual.name
        FROM feature_cvtermprop fcvtp
        JOIN cvterm qual ON fcvtp.type_id = qual.cvterm_id
        WHERE qual.name IN (
            'enables',
            'contributes_to',
            'involved_in',
            'acts_upstream_of',
            'acts_upstream_of_positive_effect',
            'acts_upstream_of_negative_effect',
            'located_in',
            'part_of',
            'is_active_in',
            'colocalizes_with'
        )
    ")
);
$qual_query->execute or die print_log("Can't query for qualifiers");
while ( my ( $fcvtid, $qual ) = $qual_query->fetchrow_array() ) {
    $quals{$fcvtid} = $qual;
}

# DB-893: Hash to store GO extensions: keys are feature_cvterm_id plus rank with intervening underscore char.
print_log("INFO: Getting the GO extension information.");
my %go_xtns;
my $go_xtn_query = $dbh->prepare(
    sprintf("
        SELECT fcvtp.feature_cvterm_id, fcvtp.rank, fcvtp.value
        FROM feature_cvtermprop fcvtp
        JOIN cvterm cvt ON cvt.cvterm_id = fcvtp.type_id
        WHERE cvt.name = 'go_annotation_extension'
    ")
);
my $go_xtn_counter = 0;
$go_xtn_query->execute or die print_log("Can't query for GO extensions");
while ( my ( $fcvtid, $rank, $go_xtn_text ) = $go_xtn_query->fetchrow_array() )
{
    $go_xtns{ $fcvtid . '_' . $rank } = $go_xtn_text;
    $go_xtn_counter++;
}
print_log("INFO: Found $go_xtn_counter GO annotation extensions.");

# Build a lookup has to exclude genes annotated with 'transposable_element_gene' term.
my %te_genes;
my $te_query = $dbh->prepare(
    sprintf(
        "SELECT f.uniquename
     FROM   feature f, feature_cvterm fc, cvterm c, cv
     WHERE  f.feature_id = fc.feature_id and fc.cvterm_id = c.cvterm_id
     and    c.cv_id = cv.cv_id and cv.name = 'SO' and c.name = 'transposable_element_gene'"
    )
);
$te_query->execute or die print_log("Can't query for te genes");
while ( ( my $teg_uname ) = $te_query->fetchrow_array() ) {
    $te_genes{$teg_uname} = 1;
}

# we can prepare this  next query once and then bind the parameters of the FBrf
# in the place of the ? in the query as we go through the loop of returned results from the $ga_query
# this is in theory more efficient than both preparing and then executing the statement within the loop
# and I think in general we can do many fetches within the loop rather than setting up a lookup first
# as we did for synonyms because in this case only a minority subset of all references will be checked
# but it may turn out more efficient to do similar to synonyms (remains to be seen)
my $pmid_query = $dbh->prepare(
    sprintf(
        "SELECT CASE WHEN pt.name = 'supplementary material' 
            THEN (SELECT d.accession
                  FROM   pub_relationship pr, pub sm, pub p, pub_dbxref pd, dbxref d, db, cvterm c
                  WHERE  sm.uniquename = ? 
                    and ((sm.pub_id = pr.subject_id and pr.object_id = p.pub_id) or (sm.pub_id = pr.object_id and pr.subject_id = p.pub_id))  
                    and  pr.type_id = c.cvterm_id and c.name = 'related_to'
                    and  p.pub_id = pd.pub_id and pd.dbxref_id = d.dbxref_id and pd.is_current = true
                    and  d.db_id = db.db_id and db.name = 'pubmed')
            ELSE (SELECT d.accession
                  FROM   pub p, pub_dbxref pd, dbxref d, db
                  WHERE  p.uniquename = ? and p.pub_id = pd.pub_id
                    and  pd.dbxref_id = d.dbxref_id and pd.is_current = true
                    and  d.db_id = db.db_id and db.name = 'pubmed')
            END as PMID
      FROM pub, cvterm pt where pub.uniquename = ? and pub.type_id = pt.cvterm_id"
    )
);

# DB-943: Improve gene product types listed for ncRNA genes.
# Build, in layers, a gene feature_id-keyed hash of gene product types.
print_log("INFO: Get gene product types.");
# 1. Create feature_id-keyed hash of types based on annotated transcripts.
#    This assumes one transcript type per gene, the current convention.
my %gene_product_types;
my $transcript_type_query = $dbh->prepare(
    sprintf("
        SELECT DISTINCT gene.feature_id, cvt.name
        FROM feature trpt
        JOIN cvterm cvt ON cvt.cvterm_id = trpt.type_id
        JOIN featureloc fl ON fl.feature_id = trpt.feature_id
        JOIN feature_relationship fr ON fr.subject_id = trpt.feature_id
        JOIN cvterm cvtfr ON cvtfr.cvterm_id = fr.type_id AND cvtfr.name = 'partof'
        JOIN feature gene ON gene.feature_id = fr.object_id
        JOIN organism o ON o.organism_id = gene.organism_id
        WHERE gene.is_obsolete IS FALSE
          AND gene.is_analysis IS FALSE
          AND gene.uniquename ~ '^FBgn[0-9]{7}\$'
          AND o.abbreviation = 'Dmel'
          AND trpt.is_obsolete IS FALSE
          AND trpt.is_analysis IS FALSE
          AND trpt.uniquename ~ '^FBtr[0-9]{7}\$'
          AND trpt.name !~ '-XR\$'
    ")
);
my $trpt_type_counter = 0;
$transcript_type_query->execute or die print_log("Can't query for gene transcript types.");
while ( my ( $fid, $trpt_type ) = $transcript_type_query->fetchrow_array() )
{
    if ( $trpt_type eq 'mRNA' ) {
        $gene_product_types{$fid} = 'protein';
    }
    elsif ( $trpt_type eq 'pre_miRNA' ) {
        $gene_product_types{$fid} = 'miRNA';
    }
    else {
        $gene_product_types{$fid} = $trpt_type;
    }
    $trpt_type_counter++;
}
print_log("INFO: Found $trpt_type_counter gene product types for current localized Dmel genes.");

# 2. Get more detailed gene product types for some ncRNA genes.
#    This only works on reporting builds where 'promoted_gene_type' is available.
#    For a production db, ncRNA genes simply have 'ncRNA' gene product.
my %ncrna_gene_class_mapping = (
    '@SO0001269:SRP_RNA_gene@' => 'SRP_RNA',
    '@SO0001640:RNase_MRP_RNA_gene@' => 'RNase_MRP_RNA',
    '@SO0001639:RNase_P_RNA_gene@' => 'RNase_P_RNA',
    '@SO0002127:lncRNA_gene@' => 'lncRNA',
    '@SO0002182:antisense_lncRNA_gene@' => 'lncRNA',
    '@SO0002353:sbRNA_gene@' => 'sbRNA',
);
my $promoted_ncrna_gene_type_query = $dbh->prepare(
    sprintf("
        SELECT DISTINCT f.feature_id, fp.value
        FROM feature f
        JOIN organism o ON o.organism_id = f.organism_id
        JOIN featureprop fp ON fp.feature_id = f.feature_id
        JOIN cvterm t ON t.cvterm_id = fp.type_id
        WHERE f.is_obsolete IS FALSE
          AND f.is_analysis IS FALSE
          AND f.uniquename ~ '^FBgn[0-9]{7}\$'
          AND o.abbreviation = 'Dmel'
          AND t.name = 'promoted_gene_type'
          AND fp.value ~ 'RNA_gene'
    ")
);
my $result_counter = 0;
my $update_gp_type_counter = 0;
$promoted_ncrna_gene_type_query->execute or die print_log("Can't query for ncRNA gene promoted gene types.");
while ( my ( $fid, $pgt_value ) = $promoted_ncrna_gene_type_query->fetchrow_array() )
{
    if ( $ncrna_gene_class_mapping{$pgt_value} ) {
        $gene_product_types{$fid} = $ncrna_gene_class_mapping{$pgt_value};
        $update_gp_type_counter++;
    }
    $result_counter++;
}
print_log("INFO: Found $result_counter ncRNA promoted_gene_type values for current Dmel genes.");
print_log("INFO: Updated gene product type for $update_gp_type_counter ncRNA genes.");

# Query for gene model status.
( my $promoted_prop_tyid ) = $dbh->selectrow_array(
    sprintf(
"SELECT cvterm_id FROM cvterm c WHERE c.name = 'derived_gene_model_status'"
    )
);
my $tyq;
if ( defined $promoted_prop_tyid ) {
    $tyq = $dbh->prepare(
        sprintf(
            "SELECT value FROM featureprop WHERE type_id = $promoted_prop_tyid
        and feature_id = ?"
        )
    );
}

# Execute the big query
print_log("INFO: Querying for ga info");
$ga_query->execute or die print_log("Can't get ga info");

# Fetch the results

# Do we want to declare some hash or array to hold the info for sorting here? see below.
my %GA_results;
my %pseudos_withdrawn;

my %seen_pubs
  ; # tracking hash so we don't need to query for so many pubs with complex query
my %seen_gns;    # tracking hash with type as value
my $rows;

# Main loop: assess ~150,000 GO annotations (~3h).
# note I am explicitly declaring all variables in results
# we could also just fetch an array and refer to array index numbers which might be slightly
# more efficient but harder to keep track of
print_log("INFO: Processing results of GA query");
while (
    my (
        $fid, $fbid, $symb, $fcvtid, $asp, $goid, $pub,
        $orgn, $evid, $src, $date, $is_not, $ev_rank
    )
    = $ga_query->fetchrow_array()
  )
{
    $rows++;
    # print_log("DEBUG: 1. Start processing row #$rows: feature_cvterm_id=$fcvtid");
    my $pseudoflag;
    my $extratabsflag;

    # skip te genes
    next if $te_genes{$fbid};

    # here is where some decisions need to be made on formatting
    # do you want the lines of the file sorted in any particular order for example by gene symbol
    # or go id
    # if so we will need to build up some type of data structure to be able to sort the values
    # it can be as simple as an array of lines but then the sorting function would be somewhat more
    # complicated or it could be some combination of hashes allowing sequential sorting of the keys

    # basically you will want to set up the line for the GA file as you like it and then either print
    # it to your output or put it in the data structure for later sorting and printing
    my $line = "FB\t";    # here is the variable to store the GA file line with column one info

    # this bit here will do the query for the pub med id for the pub returned
    my $pmid;

    if ( exists $seen_pubs{$pub} ) {
        $pmid = $seen_pubs{$pub} if defined $seen_pubs{$pub};
        # print_log("DEBUG: 2a. pub already seen.");
    }
    else {
        $pmid_query->bind_param( 1, $pub );
        $pmid_query->bind_param( 2, $pub );
        $pmid_query->bind_param( 3, $pub );
        $pmid_query->execute or die print_log("Can't fetch pub med id for $pub");
        # print_log("DEBUG: 2b1. Queried for pub info.");

        # because any pub or supplementary info of a pub should only return one pubmed id we will do
        # a single fetch - but it is still possible that we won't find a PMID so check
        ($pmid) = $pmid_query->fetchrow_array();

        # add to seen hash if found
        if ($pmid) {
            $pmid = "PMID:$pmid";
            $seen_pubs{$pub} = $pmid;
        }
        else {
           # do other checks and associations based on JIRA DB-233 specification
            $pmid = check4doi( $dbh, $pub );
            if ($pmid) {
                $pmid = "DOI:$pmid";
                $seen_pubs{$pub} = $pmid;
            }
            else {
                ( my $ptyid ) = $dbh->selectrow_array(
                    "SELECT type_id FROM pub WHERE uniquename = '$pub'");
                if (   $pubsbytype{$ptyid} eq 'review'
                    or $pubsbytype{$ptyid} eq 'paper'
                    or $pub eq 'FBrf0045941' )
                {
                    $pmid = 'GO_REF:0000095';
                    print OUT3 "$pub|$pmid\n";
                }
                elsif ( $gorefs{$pub} ) {
                    $pmid = $gorefs{$pub};
                }
                elsif (
                    $pubsbytype{$ptyid} eq 'personal communication to FlyBase' )
                {
                    $pmid = 'GO_REF:0000097';
                }
                elsif ( $pubsbytype{$ptyid} eq 'abstract' ) {
                    $pmid = 'GO_REF:0000098';
                }
                elsif ( $pubsbytype{$ptyid} eq 'DNA/RNA sequence record' ) {
                    $pmid = 'GO_REF:0000099';
                }
                elsif ( $pubsbytype{$ptyid} eq 'protein sequence record' ) {
                    $pmid = 'GO_REF:0000106';
                }
                else {
                    $pmid = undef;
                    print OUT3 "$pub\n";
                }
                $seen_pubs{$pub} = $pmid;
            }
        # print_log("ERROR: No PMID for $pub");
        }
        # print_log("DEBUG: 2b2. Assessed query results for pub info.");

    }

    # going to need to parse the $evid field (i.e., feature_cvtermprop.value of type evidence) looking for with's and ands
    # note that there can be multiple evidence statements in the same property value
    # and one or more of them may have info for col 8
    # lets set up a sub to do the parsing and then figure out the best way to deal with
    # multiple lines
    my @cols_7_8 = parse_evidence_bits( $fcvtid, $evid, $dbh );
    # print_log("DEBUG: 3. Parsed evidence bits.");

    # start building the line
    # cols 2 and 3
    $line .= "$fbid\t$symb\t";
    # print_log("DEBUG: 4. Built line col1-3.");

    # handle negation (part of col 4).
    # $line .= "IS_NOT VALUE: $is_not |||";
    if ($is_not) {
        $line .= "NOT|";
    }
    # print_log("DEBUG: 5a. Filled in col4: NOT.");

    # handle mandatory gp2term qualifiers (part of col 4).
    if ( exists( $quals{$fcvtid} ) ) {
        $line .= "$quals{$fcvtid}\t";
    }
    else {
        print_log("ERROR: Missing gp2term qualifier for this annotation: fcvt_id = $fcvtid, gene = $symb ($fbid), cvterm = GO:$goid, pub = $pub, is_not = $is_not.");
        die print_log("ERROR: FAILING due to data issue: some annotations are missing gp2term qualifers.");
    }
    # print_log("DEBUG: 5b. Filled in col4: gp2term.");

    # cols 5 and 6
    # $line .= "GO:$goid\tFB:$pub";    # DB-1002
    $line .= "GO:$goid\t";
    # print_log("DEBUG: 6. Filled in cols 5-6: GO & FB pub IDs.");

    # check for PMID
    $line .= "$pmid" if $pmid;
    $line .= "\t";
    # print_log("DEBUG: 7. Filled in col 7: PubMed ID.");

    # evidence code col 7 and optional with col 8
    # NOTE: as these can have values that might need to be split over multiple lines
    # at this point we just put in a place holder to substitute in the values
    $line .= "PUT_COLS_8_9_HERE\t";
    # print_log("DEBUG: 8. Filled in cols 7-8 with evidence placeholed.");

    # aspect col 9
    $line .= "$ASP{$asp}\t";
    # print_log("DEBUG: 9. Filled in col 9: aspect.");

    # optional fullname col 10
    if ( $fullnames{$fbid} ) {
        my $fn = $fullnames{$fbid};
        if ( $fn =~ /\t/ ) {
            $extratabsflag = $fn;
            $fn =~ s/\t/ /g;
        }
        $line .= $fn;
        $fullnames{$fbid} = $fn;
    }
    $line .= "\t";
    # print_log("DEBUG: 10. Filled in col 10: fullname.");

    # synonyms col 11
    if ( $synonyms{$fbid} ) {
        my @syns = sort @{ $synonyms{$fbid} };
        foreach my $s (@syns) {
            next if ( $s eq 'unnamed' );
            if ( $s =~ /\t/ ) {
                $extratabsflag = $s;
                $s =~ s/\t/ /g;
            }
            next
              if ( ( $s eq $symb )
                or ( $fullnames{$fbid} and ( $s eq $fullnames{$fbid} ) ) );
            $line .= "$s|";
        }
        # this little bit just checks to make sure there are any synonyms left
        if ( $line =~ /\|$/ ) {
            $line =~ s/\|$/\t/;
        }
        else {
            $line .= "\t";
        }
    }
    else {
        $line .= "\t";
    }
    # print_log("DEBUG: 11. Filled in col 11: synonyms.");

    # col 12 - determine the gene product type.
    # DB-943: Use more detailed gene product descriptors.
    # Note: In old code, the gene product look up was the rate limiting step.
    #       When a new gene was encountered in the list of GO annotations, 
    #       it would take about 700ms to get the gene product value. For
    #       ~14,763 genes in FB2024_01, this adds up to ~2h52m, essentially
    #       the entire run time for this script. New method builds the lookup
    #       at the start of the script to avoid repeated slow db queries.
    my $gp_type = 'gene_product';
    if ( $gene_product_types{$fid} ) {
        if ( $gene_product_types{$fid} ne 'pseudogene' ) {
            $gp_type = $gene_product_types{$fid};
        }
        else {
            $pseudoflag = 1;
        }
    }
    $line .= "$gp_type\t";
    # print_log("DEBUG: 12. Filled in col 12: gene product type.");

    # col 13
    $line .= "taxon:$TAX{$orgn}\t";
    # print_log("DEBUG: 13. Filled in col 13: taxon.");

    # col 14
    $line .= "$date\t";
    # print_log("DEBUG: 14. Filled in col 14: date.");

    # col 15
    $line .= "$src\t";
    # print_log("DEBUG: 15. Filled in col15: source.");

    # col 16
    # handle GO annotation extensions.
    if ( exists( $go_xtns{ $fcvtid . '_' . $ev_rank } ) ) {
        $line .= "$go_xtns{$fcvtid . '_' . $ev_rank}\n";
    }
    else {
        $line .= "\n";
    }
    # print_log("DEBUG: 16. Filled in col16: GO annotation extension.");

    # add the evidence and with cols to as many lines as needed
    print_log("ERROR: Missing evidence for:\nt$line") and next unless @cols_7_8;
    my $line_w_ev     = '';
    my $mismatch_line = '';
    foreach my $c (@cols_7_8) {
        my $evc_gn_mismatch;
        if ( $c =~ /PROBLEM/ ) {
            print_log("ERROR: $c IS A PROBLEM FOR $line");
            next;
        }
        if ( $c =~ /MISMATCH/ ) {
            $c =~ s/MISMATCH:(.*)$//;
            $evc_gn_mismatch = $1;
        }
        my $line_copy = $line;
        $line_copy =~ s/PUT_COLS_8_9_HERE/$c/;
        if ($evc_gn_mismatch) {
            my $mline = $line_copy;
            chomp $mline;
            $mismatch_line .=
              $mline . "evc_has_fbgn_symbol_mismatch ($evc_gn_mismatch)\n";
        }
        $line_w_ev .= "$line_copy";
    }
    $line = $line_w_ev;
    # print_log("DEBUG: 17. Went back and filled in cols 7-8: evidence.");

    # and here's where to decide what to do with the line
    #  print $line;
    # and we'd like to sort if first by gene symbol and then by pub so build a HOH
    push @{ $GA_results{$symb}{$pub} }, $line;
    # print_log("DEBUG: 18. Add line to hash.");

    # generate lines for problem reports
    if ($pseudoflag) {
        my $rep_line = $line;
        $rep_line =~ s/\n/pseudogene\n/;
        push @{ $pseudos_withdrawn{$symb}{$pub} }, $rep_line;
    }

    if ($tyq) {
        $tyq->bind_param( 1, $fid );
        $tyq->execute or die print_log("Can't do gene type query");
        ( my $status ) = $tyq->fetchrow_array();
        if ( $status eq 'Withdrawn' ) {
            my $rep_line = $line;
            $rep_line =~ s/\n/withdrawn\n/;
            push @{ $pseudos_withdrawn{$symb}{$pub} }, $rep_line;
        }
    }

    if ($mismatch_line) {
        push @{ $pseudos_withdrawn{$symb}{$pub} }, $mismatch_line;
    }

    if ($extratabsflag) {
        my $rep_line = $line;
        $rep_line =~ s/\n/TAB IN: $extratabsflag\n/;
        push @{ $pseudos_withdrawn{$symb}{$pub} }, $rep_line;
    }
    # print_log("DEBUG: 19. Done with additional lines checks.");
}
# End of main loop assessing GO annotations.

# and here we output the sorted lines sorted first by symbol and then by FBrf number
print_log("INFO: Producing output file");
my $lcnt = 0;
foreach my $s ( sort keys %GA_results ) {
    foreach my $r ( sort keys %{ $GA_results{$s} } ) {
        foreach my $l ( @{ $GA_results{$s}->{$r} } ) {
            print OUT $l;
            $lcnt++;

            #      print $l;
        }
    }
}
print_log("INFO: Printed out $lcnt lines for GA file.");

my $rlcnt = 0;
foreach my $s ( sort keys %pseudos_withdrawn ) {
    foreach my $r ( sort keys %{ $pseudos_withdrawn{$s} } ) {
        foreach my $l ( @{ $pseudos_withdrawn{$s}->{$r} } ) {
            print OUT2 $l;
            $rlcnt++;

            #      print $l;
        }
    }
}
print_log("INFO: Printed out $rlcnt lines for the lines to review report.");

my $end = localtime();
print_log("INFO: STARTED: $start\tENDED:   $end");

#$dbh->finish();
#$dbh->disconnect();

# checks to see if pub has a DOI number
sub check4doi {
    my $dbh  = shift;
    my $fbrf = shift;
    my $stmt =
"SELECT accession FROM pub p, pub_dbxref pd, dbxref d, db WHERE p.pub_id = pd.pub_id and pd.dbxref_id = d.dbxref_id and pd.is_current = true and d.db_id = db.db_id and db.name = 'DOI' and p.uniquename = ?";
    my $query = $dbh->prepare($stmt);
    $query->bind_param( 1, $fbrf );
    $query->execute or print_log("WARNING: Can't do doi query");
    ( my $doi ) = $query->fetchrow_array();    #expecting only one current one
    return $doi;
}

# For GO annotations, the "evidence_code" feature_cvtermprop.value starts with
# the actual evidence code (fully written out); this must be converted to the
# appropriate abbreviation for the GAF "evidence" column.
# The evidence code in the feature_cvtermprop.value may be followed by the word
# "with" and a list of cross-references to other objects (e.g., FB genes,
# UniProt IDs, etc); these xrefs will be reported in the GAF "with" column.
# The xref list can be a bit dirty - various separators, many instances of the
# the word "with", extra prefixes (e.g., FLYBASE) or suffices (e.g., the name
# of the FB gene) that should not be reported). All this extra stuff must be
# stripped out to report a comma-separated list of db:accession strings.
sub parse_evidence_bits {
    my $fcvtid = shift;    # The feature_cvterm_id for the annotation. 
    my $evbit = shift;      # The evidence_code feature_cvtermprop.value.
    my $dbh = shift;
    my @evlines;

# First, substitute any spelled out evidence_code(s) with the corresponding abbreviation(s).
# Any instances of the word "with" in the written out ev code are removed.
    foreach my $code ( sort keys %EVC ) {
        $evbit =~ s/$code/$EVC{$code}/g;
    }

    # Second, split up the value in case there are many evidence codes in there.
    # It was once the case that the feature_cvtermprop.value could contain
    # many evidence code bits separated by the word "AND".
    # This no longer seems to be true, but we're keeping this approach just to
    # be safe.
    my @evidence = trim( split / AND /, $evbit );

    # Third, in each evidence code string, split the evidence code abbreviation
    # from the list of cross-references (the word "with" separates them).
    foreach my $e (@evidence) {
        my $line;
        ( my $evc, my $dbxrefs ) = trim( split / with | from /, $e );

        # Check that there is a recognized evidence code abbreviation.
        unless ( grep $evc eq $_, values %EVC ) {
            print_log("WARNING: Unrecognized evidence code: $evc");
            push @evlines, "PROBLEM: $e";
            next;
        }

        # Fourth, get the list of xrefs using the get_dbxrefs() subroutine.
        ( my $col8, my $mismatch ) = get_dbxrefs( $fcvtid, $dbxrefs, $dbh, $evc )
          if $dbxrefs;

        # Fifth, combine the evidence code abbreviation and the xrefs as a
        # col7\tcol8 string and push to the array.
        if ($col8) {
            $line = "$evc\t$col8";
            $line = "${line}MISMATCH:$mismatch" if ($mismatch);
        }
        else {
            $line = "$evc\t";
        }
        push @evlines, $line;
    }
    return @evlines;
}

# DB-919: Ensure that all IDs are returned (not just the first), and, that any
# extraneous bits are stripped from the db:dbxref.accession string.
# For example, for "FLYBASE:rod; FB:FBgn0003268,FLYBASE:Zw10; FB:FBgn0004643",
# we should return "FB:FBgn0003268,FB:FBgn0004643",
# rather than "FB:FBgn0003268,FLYBASE:Zw10" (currently the case).
sub get_dbxrefs {
    my $fcvtid = shift;    # The feature_cvterm_id for the annotation.
    my $inline  = shift;
    my $dbh     = shift;
    my $evc     = shift;
    my $outline = '';
    my $nomatch;

    # Split out xrefs into an array, trimming flanking white space.
    # NOTE - comma may or may not have a trailing space.
    my @dbxrefs = trim( split /,/, $inline );

    # Prepare a query that confirms that FB IDs are still current.
    my $stmt =
"SELECT feature_id FROM feature WHERE is_obsolete = false and name = ? and uniquename = ?";
    my $query = $dbh->prepare($stmt);

# Clean up each xref in the list.
# Special handling for FB xrefs - restricted to genes, having this pattern:
# FLYBASE:feature.name; FB:feature.uniquename
# e.g., "FLYBASE:rod; FB:FBgn0003268"
# In evidence codes, the semi-colon is only used in this way; it's always followed by a space.
    foreach my $d (@dbxrefs) {
        my @parts = split /; /, $d;
        if ( $parts[1] ) {
            my $xref = trim( $parts[1] );
            $outline .= "$xref|";
        }
        else {
            my $xref = $parts[0];
            $outline .= "$xref|";
        }

        # For FB xrefs, check that the name and uniquename match.
        if ( $parts[0] =~ /FLYBASE:(.+)/ ) {
            my $symb = $1;
            $symb = decon($symb);
            print_log("WARNING: feature_cvterm_id=$fcvtid: missing symbol in evidence code line: $inline\n")
              unless ($symb);
            if ( $parts[1] =~ /FB:(FBgn[0-9]{7})/ ) {
                my $fbgn = $1;

                # print "CHECKING FOR MATCH BETWEEN $symb and $fbgn");
                $query->bind_param( 1, $symb );
                $query->bind_param( 2, $fbgn );
                $query->execute
                  or
                  print_log("WARNING: feature_cvterm_id=$fcvtid: can't execute $stmt FOR $symb:$fbgn");
                unless ( $query->rows() > 0 ) {
                    $nomatch = $fbgn . ':' . $symb;
                    print_log("WARNING: feature_cvterm_id=$fcvtid: WE HAVE A MISMATCH for $symb:$fbgn");
                }
            }
            else {
                print_log("WARNING: feature_cvterm_id=$fcvtid: no FBgn provided for $symb in evidence code line: $inline");
            }
        }
    }

    # Trim off the xref divider at end of line.
    $outline =~ s/\|$//;

    # For certain evidence codes, use commas instead of pipes.
    if ( $evc =~ /HGI|IBA|IC|IGC|IGI|IMP|IPI|ISA|ISS/ ) {
        $outline =~ s/\|/,/g;
    }

    return ( $outline, $nomatch );
}

# attempts to fetch GO.references file from geneontology.org site and parse to lookup hash
# keyed by FBrf with GO ref as value
sub fetch_and_parse_gorefs {
    my $ua = LWP::UserAgent->new;
    my $req =
      HTTP::Request->new( GET =>
'https://raw.githubusercontent.com/geneontology/go-site/master/metadata/gorefs.yaml'
      );
    my $response = $ua->request($req);
    my $page     = $response->content;
    unless ( $page =~ /^- id: GO_REF:/ ) {
        print_log("ERROR: Cannot open the gorefs.yaml document.");
        return;
    }
    my @lines = split /\n/, $page;
    my %fbrf2goref;
    my $goid = '';
    my $fbrf = '';
    my $current_go_ref = 1;
    foreach my $l (@lines) {
        if ( $l =~ /^-\sid:\s(GO_REF:[0-9]+)/ ) {
            # Before parsing next stanza, check if previous one represented a
            # current GO_REF-to-FBrf xref and add to hash if it does.
            if ( $goid && $fbrf && $current_go_ref ) {
                $fbrf2goref{$fbrf} = $goid;
            }
            # Reset stanza info once the previous one has been recorded.
            $goid = $1;
            $fbrf = '';
            $current_go_ref = 1;
        }
        elsif ( $l =~ /-\sFB:(FBrf[0-9]{7})/ ) {
            $fbrf = $1 if ( $l =~ /-\sFB:(FBrf[0-9]{7})/ );
        }
        elsif ( $l =~ /is_obsolete:\strue/ ) {
            $current_go_ref = 0;
        }
    }
    # Just in case the last stanza had a GO-FB xref, add it.
    # Because in loop above, hash additions happen when a new stanza is found.
    # There is no marker to the end of a GO_REF other than the start of the next
    # one, or, the end of the file.
    if ( $goid && $fbrf && $current_go_ref ) {
        $fbrf2goref{$fbrf} = $goid;
    }
    # Suppressing hard-coded associations as these seem to be in the gorefs.yaml file.
    # $fbrf2goref{'FBrf0253064'} = 'GO_REF:0000115';    # DB-767
    # # $fbrf2goref{'FBrf0253063'} = 'GO_REF:0000024';    # DB-823
    # $fbrf2goref{'FBrf0255270'} = 'GO_REF:0000024';    # DB-823
    # $fbrf2goref{'FBrf0254415'} = 'GO_REF:0000047';    # DB-811
    # $fbrf2goref{'FBrf0258542'} = 'GO_REF:0000033';    # DB-928
    print_log("INFO: Constructed FBrf -> GO_REF Mapping:");
    print Dumper( \%fbrf2goref );
    return \%fbrf2goref;
}

# Converts SGML-formatted symbols to 'symbol_plain' format (modified from conv_greeks)
sub decon {
    my $string = shift;
    $string =~ s/&agr\;/alpha/g;
    $string =~ s/&Agr\;/Alpha/g;
    $string =~ s/&bgr\;/beta/g;
    $string =~ s/&Bgr\;/Beta/g;
    $string =~ s/&ggr\;/gamma/g;
    $string =~ s/&Ggr\;/Gamma/g;
    $string =~ s/&dgr\;/delta/g;
    $string =~ s/&Dgr\;/Delta/g;
    $string =~ s/&egr\;/epsilon/g;
    $string =~ s/&Egr\;/Epsilon/g;
    $string =~ s/&zgr\;/zeta/g;
    $string =~ s/&Zgr\;/Zeta/g;
    $string =~ s/&eegr\;/eta/g;
    $string =~ s/&EEgr\;/Eta/g;
    $string =~ s/&thgr\;/theta/g;
    $string =~ s/&THgr\;/Theta/g;
    $string =~ s/&igr\;/iota/g;
    $string =~ s/&Igr\;/Iota/g;
    $string =~ s/&kgr\;/kappa/g;
    $string =~ s/&Kgr\;/Kappa/g;
    $string =~ s/&lgr\;/lambda/g;
    $string =~ s/&Lgr\;/Lambda/g;
    $string =~ s/&mgr\;/mu/g;
    $string =~ s/&Mgr\;/Mu/g;
    $string =~ s/&ngr\;/nu/g;
    $string =~ s/&Ngr\;/Nu/g;
    $string =~ s/&xgr\;/xi/g;
    $string =~ s/&Xgr\;/Xi/g;
    $string =~ s/&ogr\;/omicron/g;
    $string =~ s/&Ogr\;/Omicron/g;
    $string =~ s/&pgr\;/pi/g;
    $string =~ s/&Pgr\;/Pi/g;
    $string =~ s/&rgr\;/rho/g;
    $string =~ s/&Rgr\;/Rho/g;
    $string =~ s/&sgr\;/sigma/g;
    $string =~ s/&Sgr\;/Sigma/g;
    $string =~ s/&tgr\;/tau/g;
    $string =~ s/&Tgr\;/Tau/g;
    $string =~ s/&ugr\;/upsilon/g;
    $string =~ s/&Ugr\;/Upsilon/g;
    $string =~ s/&phgr\;/phi/g;
    $string =~ s/&PHgr\;/Phi/g;
    $string =~ s/&khgr\;/chi/g;
    $string =~ s/&KHgr\;/Chi/g;
    $string =~ s/&psgr\;/psi/g;
    $string =~ s/&PSgr\;/Psi/g;
    $string =~ s/&ohgr\;/omega/g;
    $string =~ s/&OHgr\;/Omega/g;
    $string =~ s/\<\/down\>/\]\]/g;
    $string =~ s/\<down\>/\[\[/g;
    $string =~ s/\<up\>/\[/g;
    $string =~ s/\<\/up\>/\]/g;
    return $string;
}
