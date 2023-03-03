#!/usr/bin/perl
# report_gp2protein
#
#       Perl script produces gp2protein.fb report for GO Consortium
#	
#       See JIRA DB-231 for details
#	
#	This report gets included in with other ontology files produced with
#	each release
#
#-----------------------------------------------------------------------------#
#
#	NOTES
#	
#	perl report_gp2protein flysql9 fb_2015_01_reporting pguser pgpwd log_file output_fb_file
#	
#	
#	
#-----------------------------------------------------------------------------#
use DBI;

if (@ARGV != 6) {
    print "\n USAGE: report_gp2protein pg_server db_name pg_username pg_password logging_output_filename report_output_filename\n\n";
    exit();
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);

$OUTFILE = shift(@ARGV);
$OUTFILE = '/dev/null';  ## <=== Comment this line out and run script to see log output 
open(STDERR,">>$OUTFILE");
open(STDOUT,">>$OUTFILE");
 
$REPFILE = shift(@ARGV);
open(REPOUT,">>$REPFILE");


#
##  DB Connections
#

## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
# $dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
# $dbh5 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
# $dbh6 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


#
##  General setup
#

## Setup output file header
$jetzt = scalar localtime;
# print "## Starting report_gp2protein\t$jetzt\n## Using data-source: $dsource\n\n";
# print REPOUT "## FlyBase gp2protein.fb report\n";
# print REPOUT "## Generated: $jetzt\n";
# print REPOUT "## Using datasource: $dsource\n\n";


#
# Main methods
#

## Get localized Dmel genes that have links to accessions from 'UniProtKB/TrEMBL' or 'UniProtKB/Swiss-Prot' dbs
print "\nGetting genes & uniprot ids...\n";
my %gh;
my $fbgnwc = 'FBgn%';
my $gq = $dbh->prepare(sprintf("select f.is_obsolete, f.feature_id, f.uniquename, f.name, fd.is_current, dx.db_id, db.name as dbname, dx.dbxref_id, dx.accession from organism o, feature f natural join featureloc, feature_dbxref fd, dbxref dx, db where f.organism_id = o.organism_id and o.abbreviation = 'Dmel' and f.is_obsolete = 'f' and f.feature_id = fd.feature_id and fd.is_current = 't' and fd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name in ('UniProt/TrEMBL','UniProt/Swiss-Prot') and f.uniquename like '%s'",$fbgnwc));
$gq->execute or die "WARNING: ERROR: Unable to execute gene/dbxref query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
     print "Processing gene/dbxref:\t$gr{feature_id}\t$gr{uniquename}\t$gr{name}\t$gr{is_current}\t$gr{dbname}\t$gr{accession}\n";
## Check SO terms to eliminate pseudogene_attribute and non_protein_coding_gene genes
    my $skippit = 0;
    my $cq = $dbh2->prepare(sprintf("SELECT cvt.cvterm_id, cvt.name from feature_cvterm fc, cvterm cvt, cv where fc.feature_id = %d and fc.cvterm_id = cvt.cvterm_id and cvt.cv_id = cv.cv_id and cv.name = 'SO'",$gr{feature_id}));
    $cq->execute or die "WARNING: ERROR: Unable to execute SO query\n";
    my $cq_cnt = $cq->rows;
    if ($cq_cnt > 0) {
	while (my %cr = %{$cq->fetchrow_hashref}) {
	    if (($cr{name} =~ /pseudogene/) || ($cr{name} =~ /non\_protein\_coding/)) {
		$skippit++;
		print "\tSkipping gene ($cr{name})\n";
	    }
	}
    }
    else {
	print "\tWARNING: ERROR: No SO cvterm returned for gene: $gr{feature_id}\t$gr{uniquename}\t$gr{name}\n";
    }
    if ($skippit == 0) {
## Push acc# into hash for report output 
	push(@{$gh{$gr{uniquename}}{$gr{dbname}}},$gr{accession});
    }
}

## Write report -- if theres a SwissProt acc#, we report that alone, if theres not, we report 
## all TrEMBL acc#
print "\nWriting report...\n";
foreach my $gene (sort(keys(%gh))) {
#    my @accnos = uniq(@{$gh{$gene}});
#    print(sprintf("REPORT\t%s\t%s\n",$gene,join(';',@accnos)));
    if ($gh{$gene}{'UniProt/Swiss-Prot'}) {
	my @accnos = uniq(@{$gh{$gene}{'UniProt/Swiss-Prot'}});
	print REPOUT (sprintf("FB:%s\tUniProtKB:%s\n",$gene,join(';UniProtKB:',@accnos)));
    }
    else {
	my @accnos = uniq(@{$gh{$gene}{'UniProt/TrEMBL'}});
	print REPOUT (sprintf("FB:%s\tUniProtKB:%s\n",$gene,join(';UniProtKB:',@accnos)));
    }
}


sub uniq {
my %seen;
return grep { !$seen{$_}++ } @_;
}


$jetzt = scalar localtime;
print "\n\n## Finished report_gp2protein: $jetzt\n";
# print REPOUT "\n\n## Finished report_gp2protein: $jetzt\n";
exit();

