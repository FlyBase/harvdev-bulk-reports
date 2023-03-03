#!/usr/local/bin/perl
# report_gene_snapshots
#
#       Perl script generates tab-delimited report_gene_snapshots table with the 
#       following fields: 
#	
#	FBgn_ID
#	gene_symbol
#	gene_name
#	datestamp (date of gene summary curation)
#	gene_snapshot_text
#
#-----------------------------------------------------------------------------#
#
#	NOTES:
#	          See JIRA DB-337 for details
#	
#
#-----------------------------------------------------------------------------#
use DBI;
use POSIX;
require "conversions.pm";
# require "/users/emmert/work/perl-sub/conversions.pm";

if (@ARGV != 5) {
    print "\n USAGE: report_gene_snapshots pg_server db_name pg_username pg_password report_output_filename\n\n";
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
#my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh5 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


#
##  General setup
#


## Setup file header
$jetzt = scalar localtime;
print  "## FlyBase Gene Snapshot Report\n## Generated: $jetzt\n";
print  "## Using datasource: $dsource...\n";
print  "## Genes that currently do not have a gene snapshot are included in this file and marked with \"Contributions welcome.\" in the gene_snapshot_text column.\n\n";
print  "##FBgn_ID\tGeneSymbol\tGeneName\tdatestamp\tgene_snapshot_text\n";


#
## Main methods
#


## Get localized Dmel genes having featureprop of type 'promoted_gene_type' with featureprop.value 'protein_coding_gene'
my $fid;
my $funame;
my $fname;
my $gq = $dbh->prepare("SELECT f.feature_id, f.uniquename, f.name from feature f, featureloc fl, featureprop fp, organism o, cvterm cvt where f.is_obsolete = 'f' and f.feature_id = fl.feature_id and f.organism_id = o.organism_id and o.abbreviation = 'Dmel' and f.feature_id = fp.feature_id and fp.type_id = cvt.cvterm_id and cvt.name = 'promoted_gene_type' and fp.value = E'\@SO0001217\:protein\_coding\_gene\@'");
$gq->execute or die "WARNING: ERROR: Unable to execute gene query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
#    print "\nProcessing gene:\t$gr{feature_id}\t$gr{uniquename}\t$gr{name}\n";
    $fid = $gr{feature_id};
    $funame = $gr{uniquename};
    $fname = $gr{name};
## Get gene fullname or use "-"
    my $fullname;
    my $nq = $dbh2->prepare(sprintf("SELECT distinct(s.name) from feature_synonym fs, synonym s, cvterm cvt where fs.is_current = 't' and fs.synonym_id = s.synonym_id and s.type_id = cvt.cvterm_id and cvt.name = 'fullname' and fs.feature_id = %d",$fid));
    $nq->execute or die "WARNING: ERROR: Unable to execute fullname query\n";
    my $nq_cnt = $nq->rows;
    if ($nq_cnt == 1) {
	while (my %nr = %{$nq->fetchrow_hashref}) {
	    $fullname = $nr{name};
	}
    }
    else {
	$fullname = '-';
    }
## Get gene_summary_date
    my $datestamp;
    my $dq = $dbh5->prepare(sprintf("SELECT * from featureprop fp, cvterm cvt where feature_id = %d and fp.type_id = cvt.cvterm_id and cvt.name = 'gene_summary_date'",$fid));
    $dq->execute or die "WARNING: ERROR: Unable to execute datestamp query\n";
    my $dq_cnt = $dq->rows;
    if ($dq_cnt > 0) {
	while (my %dr = %{$dq->fetchrow_hashref}) {
	    $datestamp = $dr{value};
	}
    }

## Check for featureprop of type 'gene_summary_text' (the snapshot)
    my $sq = $dbh3->prepare(sprintf("SELECT * from featureprop fp, cvterm cvt where feature_id = %d and fp.type_id = cvt.cvterm_id and cvt.name = 'gene_summary_text'",$fid));
    $sq->execute or die "WARNING: ERROR: Unable to execute gene_summary_text query\n";
    my $sq_cnt = $sq->rows;
    if ($sq_cnt > 0) {
	while (my %sr = %{$sq->fetchrow_hashref}) {
	    my $plainval = &decon($sr{value});
	    $plainval =~ s/\@//g;
	    $plainval =~ s/\@//g;
	    print "$funame\t$fname\t$fullname\t$datestamp\t$plainval\n";
	}
    }
    else {
	    $datestamp = '-';
	    print "$funame\t$fname\t$fullname\t$datestamp\tContributions welcome.\n";
	}
}


$jetzt = scalar localtime;
print "## Finished report_gene_snapshots: $jetzt\n\n";

