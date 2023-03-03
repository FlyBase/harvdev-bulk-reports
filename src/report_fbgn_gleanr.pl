#!/usr/local/bin/perl
# fbgn_gleanr
#
#       Perl script generates tab-delimited fbgn_gleanr table with the 
#       following fields: gene_symbol, FBgn#, GLEANR_ID
#
#-----------------------------------------------------------------------------#
#
#	NOTES
#	
#	
#
#-----------------------------------------------------------------------------#
use DBI;
require "conversions.pm";

if (@ARGV < 5) {
    print "\n USAGE: fbgn_gleanr pg_server db_name pg_username pg_password output_filename\n\n";
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
$jetzt = scalar localtime;
print "## FlyBase FBgn-GLEANR ID Correspondence Table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print "##organism_abbreviation\tgene_symbol\tprimary_FBgn#\tGLEANR_ID\n";


#
## Main method
#
## Get genes...
my $fbgnwc = 'FBgn%';
my $gleanrwc = '%GLEANR%';
## Main driver
my $gq = $dbh->prepare(sprintf("SELECT f.uniquename, f.name, s.name as sname, abbreviation from feature f, organism o, synonym s, feature_synonym fs, cvterm cvt where f.organism_id = o.organism_id and s.type_id = cvt.cvterm_id and cvt.name = 'symbol' and f.is_obsolete = 'f' and f.is_analysis = 'f' and f.uniquename like '%s' and f.feature_id = fs.feature_id and fs.synonym_id = s.synonym_id and fs.is_current = 'f' and fs.is_internal = 'f' and s.name like '%s'",$fbgnwc,$gleanrwc));
$gq->execute or die "WARNING: ERROR: Unable to execute gleanr query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
## OUTPUT: FBgn#, GLEANR_ID
    print(sprintf("%s\t%s\t%s\t%s\n",$gr{abbreviation},$gr{name},$gr{uniquename},$gr{sname}));
}
