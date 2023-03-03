#!/usr/bin/perl
# report_funcomps
#
#       Perl script reports Dmel genes and genes that functionally complement 
#       (using derived functionally_complements f_r)
#       
#       See JIRA HDIS-11 for details
#       
#    Fields are:   
#       
#       1. Dmel gene (symbol)
#       2. Dmel gene (FBgn)
#       3. Functionally complementing ortholog (symbol)
#       4. Functionally complementing ortholog (FBgn#)
#       5. Supporting_FBrf
#       
#       
#-----------------------------------------------------------------------------#
#
#	
#	
#-----------------------------------------------------------------------------#
use DBI;

if (@ARGV != 5) {
    print "\n USAGE: report_funcomps pg_server db_name pg_username pg_password output_filename\n\n";
    exit();
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
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


#
##  General setup
#

## Setup output file header
$jetzt = scalar localtime;
print "## FlyBase Gene Functional Complementation Report\n## Using data-source: $db\n## $jetzt\n##\n";
print "## Dmel_gene_symbol\tDmel_gene_FBgn_ID\tortholog_gene_symbol\tortholog_gene_FBgn_ID\treference_FBrf_ID\n";


#
## Main driver
#

## Get set of orthologous organisms and then get their specific organism databases
my $gq = $dbh->prepare("SELECT f.uniquename as funame, f.name as fname, o.uniquename as ouname, o.name as oname, p.uniquename as puname from feature f, feature o, feature_relationship fr, cvterm cvt, feature_relationship_pub frp, pub p where f.feature_id = fr.object_id and o.feature_id = fr.subject_id and fr.type_id = cvt.cvterm_id and cvt.name = 'functionally_complements' and fr.feature_relationship_id = frp.feature_relationship_id and frp.pub_id = p.pub_id order by fname, oname, puname");
$gq->execute or die "WARNING: ERROR: Unable to execute funcomp query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
    print "$gr{fname}\t$gr{funame}\t$gr{oname}\t$gr{ouname}\t$gr{puname}\n";
}

$jetzt = scalar localtime;
print "\n\n## Finished report_funcomps: $jetzt\n";
exit();


