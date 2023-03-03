#!/usr/local/bin/perl
# report_pmid_fbgn_uniprot
#
#       Perl script generates tab-delimited report_pmid_fbgn_uniprot table 
#       listing references with a PMID that are associated with a gene which 
#       is associated with a UniProt (UniProt/Swiss-Prot, UniProt/TrEMBL)
#       accession.  Report has the following fields: 
#       FBrf#, PMID, FBgn#, UniProt_acc#
#
#-----------------------------------------------------------------------------#
#
#	NOTES
#	See email thread w/ Kathy Matthews et al. June 29, 2009
#	Subject: Re: FlyBase bibliography
#
#
#	The method for getting these is a bit strange, owing to performance 
#	considerations.  The obvious method for this report, namely getting 
#	FB pubs w/ PMIDs, then getting genes linked to UniProt IDs, was really
#	slow -- it was taking > 18hours.  The alternate method used here is to
#	first populate a hash with UniProt ids and associated features, and 
#	then using the list of features, querying by feature for features linked 
#	to papers with PMIDs.  This executes in approx. 5 minutes (on flysql)
#
#-----------------------------------------------------------------------------#
use DBI;
require "conversions.pm";

if (@ARGV < 5) {
    print "\n USAGE: report_pmid_fbgn_uniprot pg_server db_name pg_username pg_password output_filename\n\n";
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

## Setup header
$jetzt = scalar localtime;
print "## FlyBase PMID_FBgn_UniProt Table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print "##FBrf_id\tPMID\tFBgn_id\tUniProt_database\tUniProt_id\n";


## Get hash of dbxref_id & accession for UniProt/TrEMBL and UniProt/Swiss-Prot records, to
## speed processing
$jetzt = scalar localtime;
# print "\nGetting Uniprot ids: $jetzt\n";
my %uh;
my $uq = $dbh->prepare("SELECT fd.dbxref_id, accession, db.name, feature_id from feature_dbxref fd, dbxref dx, db where fd.is_current = 't' and fd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name in ('UniProt/TrEMBL','UniProt/Swiss-Prot')");
$uq->execute or die "WARNING: ERROR: Unable to execute UniProt query\n";
while (my %ur = %{$uq->fetchrow_hashref}) {
    $uh{$ur{dbxref_id}}{acc} = $ur{accession};
    $uh{$ur{dbxref_id}}{db} = $ur{name};
    push(@{$uh{$ur{dbxref_id}}{fids}},$ur{feature_id});
}
$jetzt = scalar localtime;
# print "\nFinish getting Uniprot ids: $jetzt\n";

## Get db_id for UniProt/TrEMBL and UniProt/Swiss-Prot (this speeds subsequent query vastly)
my %dbids;
my $uq = $dbh->prepare("SELECT db_id, name from db where name in ('UniProt/TrEMBL','UniProt/Swiss-Prot')");
$uq->execute or die "WARNING: ERROR: Unable to execute db query\n";
while (my %ur = %{$uq->fetchrow_hashref}) {
  $dbids{$ur{db_id}} = $ur{name};
}


#
## Main method
#
## Get genes...
my $fbgnwc = 'FBgn%';
my $fbrfwc = 'FBrf%';
## Main driver


## Now using the feature_id(s) associated with each UniProt ID, we query for associated genes
## with pubmed-linked pubs attached
foreach my $dx (keys(%uh)) {
# print "\nChecking $dx ($uh{$dx}{acc})\n";
  foreach my $feat (@{$uh{$dx}{fids}}) {
    my $fq = $dbh->prepare(sprintf("SELECT f.uniquename as funame, p.uniquename, accession from feature f, feature_pub fp, pub p, pub_dbxref pd, dbxref dx, db where p.uniquename like '%s' and p.pub_id = pd.pub_id and pd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name = 'pubmed' and f.feature_id = fp.feature_id and fp.pub_id = p.pub_id and f.uniquename like '%s' and f.feature_id = %d",$fbrfwc,$fbgnwc,$feat));
    $fq->execute or die "WARNING: ERROR: Unable to execute feature/PMID query\n";
    my $fq_cnt = $fq->rows;
    if ($fq_cnt > 0) {
      while (my %fr = %{$fq->fetchrow_hashref}) {
#	print "\tOUTPUT: $fr{uniquename}\t$fr{accession}\t$fr{funame}\t$uh{$dx}{db}\t$uh{$dx}{acc}\n";
	print "$fr{uniquename}\t$fr{accession}\t$fr{funame}\t$uh{$dx}{db}\t$uh{$dx}{acc}\n";
      }
    }
    else {
#      print "\t\tNo pubmed pub found for this dx/feat:\t$dx/$feat\n";
    }
  }
}


