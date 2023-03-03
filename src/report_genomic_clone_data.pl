#!/usr/local/bin/perl
# report_genomic_clone_data
#
#       Perl script generates tab-delimited report of genomic clone data 
#       containing the following fields: FBcl#, clone_name, acc#
#
#-----------------------------------------------------------------------------#
#
#	NOTES
#	acc# are linked via feature_dbxref to features type 
#	'BAC_cloned_genomic_insert' which are linked to genomic_clone features
#
#-----------------------------------------------------------------------------#
use DBI;
require "conversions.pm";

if (@ARGV < 5) {
    print "\n USAGE: report_genomic_clone_data pg_server db_name pg_username pg_password output_filename\n\n";
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
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";


## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
# $dbh = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh2 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh3 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";


#
##  General setup
#
## Setup file header
$jetzt = scalar localtime;
print "## FlyBase genomic clone table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print "##FBcl#\torganism_abbreviation\tclone_name\taccession\n";

#
## Main method
#
## Get clones...

my $fbclwc = 'FBcl%';
## Main query for features type genomic_clone
my $cq = $dbh->prepare(sprintf("SELECT f.uniquename, f.name, f.feature_id, cvt.name as ctype, abbreviation from feature f, cvterm cvt, organism o where f.organism_id = o.organism_id and f.type_id = cvt.cvterm_id and cvt.name = 'genomic_clone' and uniquename like '%s' and f.is_obsolete = 'f' and is_analysis = 'f'",$fbclwc));
$cq->execute or die "WARNING: ERROR: Unable to execute clone query\n";
while (my %cr = %{$cq->fetchrow_hashref}) {
#    print "\nProcessing clone: $cr{feature_id}\t$cr{uniquename}\t$cr{name}\t$cr{ctype}\n";
##  Get f_r associated BAC_cloned_genomic_insert features & acc#
    my $accno;
    my $bq = $dbh2->prepare(sprintf("SELECT f.name, f.uniquename, accession, version from feature_relationship, feature f, cvterm cvt, feature_dbxref fd, dbxref dx, db where object_id = %d and subject_id = f.feature_id and f.type_id = cvt.cvterm_id and cvt.name = 'BAC_cloned_genomic_insert' and f.feature_id = fd.feature_id and fd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name = 'GB'",$cr{feature_id}));
    $bq->execute or die "WARNINIG: ERROR: Unable to execute clone acc# query\n";
    my $bq_cnt = $bq->rows;
    if ($bq_cnt > 0) {
	while (my %br = %{$bq->fetchrow_hashref}) {
#	    print "\tinsert found: $br{accession}.$br{version}\n";
	    $accno = $br{accession};
	}
    }
## OUTPUT: FBcl#, clone name, acc#
    print "$cr{uniquename}\t$cr{abbreviation}\t$cr{name}\t$accno\n";
}
$jetzt = scalar localtime;
print "## Finished FlyBase genomic clone table: $jetzt\n";
