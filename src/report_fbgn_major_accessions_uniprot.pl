#!/usr/local/bin/perl
# report_fbgn_major_accessions_uniprot
#
#   Perl script generates tab-delimited file of gene, annotation IDs and UniprotKB/Swiss-Prot/TrEMBL xrefs.
#   Intended for use by UniProt. Contact is Dushi (djyothi@ebi.ac.uk, uniprot-prod@ebi.ac.uk)
#   Modified version of Dave's "fbgn_major_accessions" script.
#   As per DB-644, this script restricts reporting to genes passing the "coding" criteria.
#   There is an optional input for a SO term list that defines "coding".
#   If not given, script limits xrefs to genes with "protein_coding_gene" term.
#   In command line, list of SO terms should be enclosed in quotes, terms separated by spaces.

use DBI;
use Getopt::Long;
require "conversions.pm";
# require "/users/emmert/work/perl-sub/conversions.pm";

if (@ARGV < 5) {
    print "\n USAGE: report_fbgn_major_accessions pg_server db_name pg_username pg_password output_filename -s \"SO term list to filter genes\" (optional)\n\n";
    print "\toutput_filename is the output file for std and error output.\n\n";
    exit;
}
# Look for optional SO term list.
my $so_filter;
GetOptions ("so_filter=s" => \$so_filter);
if ($so_filter) {
  $so_filter = "('" . join("', '", split(' ', $so_filter)) . "')";
} else {
  $so_filter = "('protein_coding_gene')";
}

if (@ARGV != 5) {
    print "\n USAGE: report_fbgn_major_accessions pg_server db_name pg_username pg_password output_filename -s \"SO term list to filter genes\" (optional)\n\n";
    print "\toutput_filename is the output file for std and error output.\n\n";
    exit;
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);

$OUTFILE = shift(@ARGV);
open(STDERR,">>$OUTFILE");
open(STDOUT,">>$OUTFILE");

#
##  DB Connections
#
## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
## Data source (g4)
$dbh1 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";

#
##  General setup
#
## Setup file header
$jetzt = scalar localtime;
print "## FlyBase gene Uniprot accessions\n";
print "## Using chado datasource: $dsource\n";
print "##primary_FBgn#\tgene_symbol\tCG#\tUniprotKB/Swiss-Prot/TrEMBL_accession\n";

#
## Main method
#
## Get genes...
my $fbgnwc = 'FBgn%';
## Main driver
my $gq = $dbh1->prepare("SELECT DISTINCT f.uniquename, f.name, f.feature_id, cvtf.name as ftype
                        FROM feature f
                        JOIN cvterm cvtf ON (cvtf.cvterm_id = f.type_id)
                        JOIN featureprop fp ON fp.feature_id = f.feature_id
                        JOIN cvterm cvtfp ON cvtfp.cvterm_id = fp.type_id
                        WHERE cvtf.name = 'gene' and
                              f.uniquename ~ '^FBgn[0-9]{7}\$' and
                              f.is_obsolete = false and
                              f.is_analysis = false and
                              cvtfp.name = 'derived_gene_model_status' and
                              fp.value != 'Withdrawn'
                        ORDER BY f.uniquename");
$gq->execute or die "WARNING: ERROR: Unable to execute gene query...\n";

## Get associated information for each gene.
while (my %gr = %{$gq->fetchrow_hashref}) {
    ## For each gene, confirm that it is considered coding. This is based on associated SO terms.
    ## Default is "protein_coding_gene", but can pass list of terms in command line.
    ## Current criteria are association with either "protein_coding_gene" or "transposable_element_gene" terms.
    ## Per SM: JIRA DB-418, DB-644
    my $soq = $dbh2->prepare(sprintf("SELECT DISTINCT *
                                      FROM feature_cvterm fc
                                      JOIN cvterm cvt ON (cvt.cvterm_id = fc.cvterm_id)
                                      JOIN cv ON (cv.cv_id = cvt.cv_id)
                                      WHERE cvt.name in %s
                                        and cv.name = 'SO'
                                        and fc.feature_id = %d", $so_filter, $gr{feature_id}));
    $soq->execute or die "WARNING: ERROR: Unable to execute SO query\n";
    my $soq_cnt = $soq->rows;
    next if ($soq_cnt < 1);

    ## Get annotation id for the gene (if there is one)
    my $aq = $dbh3->prepare(sprintf("SELECT DISTINCT *
                                     FROM feature_dbxref fdbx
                                     JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
                                     JOIN db ON (db.db_id = dbx.db_id)
                                     WHERE fdbx.is_current = true
                                       and db.name = 'FlyBase Annotation IDs'
                                       and fdbx.feature_id = %d",$gr{feature_id}));
    $aq->execute or die "WARNING: ERROR: Unable to execute annotation ID query\n";
    my $aq_cnt = $aq->rows;
    if ($aq_cnt > 0) {
        while (my %ar = %{$aq->fetchrow_hashref}) {
            $gr{annoid} = $ar{accession};
        }
    } else {
        $gr{annoid} = '';
    }

    ## Get UniprotKB/Swiss-Prot/TreMBL acc# (via dbxref)
    my $xrefq = $dbh4->prepare(sprintf("SELECT DISTINCT *
                                     FROM feature_dbxref fdbx
                                     JOIN dbxref dbx ON (dbx.dbxref_id = fdbx.dbxref_id)
                                     JOIN db ON (db.db_id = dbx.db_id)
                                     WHERE fdbx.is_current = true
                                       and db.name in ('UniProt/Swiss-Prot','UniProt/TrEMBL') 
                                       and fdbx.feature_id = %d",$gr{feature_id}));
    $xrefq->execute or die "WARNING: ERROR: Unable to execute dbxref query\n";
    my $xrefq_cnt = $xrefq->rows;
    if ($xrefq_cnt > 0) {
	    while (my %xrefr = %{$xrefq->fetchrow_hashref}) {
	        print "$gr{uniquename}\t$gr{name}\t$gr{annoid}\t$xrefr{accession}\n";
	    }
    }
}