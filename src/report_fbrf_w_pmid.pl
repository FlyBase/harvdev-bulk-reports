#!/usr/bin/perl -w
# perl report_fbrf_w_pmid server db user pass output

# generates a report of FBrfs and their PubMed IDs
# report is used in textpresso paper retrieval pipeline (and provided to users)
# format tabbed separated
# FBrf    PMID    pubtype    miniref	release_added
#use strict;
use DBI;

usage() unless @ARGV == 5;

my $start = localtime();

print "STARTED: $start\n";

# db connection
my $server = shift @ARGV;
my $dbname = shift @ARGV;
my $user = shift @ARGV;
my $pass = shift @ARGV;
my $outfile = shift @ARGV;

############################### db connection ############################
my $chado="dbi:Pg:dbname=$dbname; host=$server;port=5432";
my $dbh = DBI->connect($chado,$user,$pass) or die "cannot connect to $chado";
my $dbh2 = DBI->connect($chado,$user,$pass) or die "cannot connect to $chado";

print "Connected to db $chado\n";
#######################################################################

open(my $out, ">", $outfile);

(my $release = $dbname) =~ s/_reporting//; 

my $stmt = "SELECT p.pub_id, p.uniquename , d.accession, c.name, p.miniref FROM pub p, pub_dbxref pd, dbxref d, db, cvterm c WHERE p.pub_id = pd.pub_id and pd.dbxref_id = d.dbxref_id and d.db_id = db.db_id and db.name = 'pubmed' and p.type_id = c.cvterm_id and  p.is_obsolete = false";

my $ppstmt = "select max(value) FROM pubprop pp, cvterm c WHERE pp.type_id = c.cvterm_id and c.name = 'pmid_added' and pp.pub_id = ?";
my $ppq = $dbh->prepare($ppstmt);

my $pubq = $dbh->prepare($stmt);
$pubq->execute or die "Can't do $stmt\n";

my %pubs;
while ((my $pid, my $fbrf, my $pmid, my $ptype, my $miniref) = $pubq->fetchrow_array()) {
  # check for pmid added value if none found log and add current release
  my $added;
  $ppq->bind_param(1, $pid);
  $ppq->execute or die "Can't do $ppstmt for $pid\n";
  ($added) = $ppq->fetchrow_array();
  unless ($added) {
    $added = $release;
    print STDERR "WARNING $fbrf with PMID=$pmid IS MISSING THE 'pmid_added' PUB PROP IN $dbname\n";
  }
  $pubs{$added}{$fbrf} = {
			  pmid => $pmid,
			  type => $ptype,
			  miniref => $miniref,
			  pmid_added => $added,
			 };
}

######## Add PMCID and DOI to this report ###############
##  DE 9-Oct-2014 (per JIRA DB-225)

foreach my $dad (keys(%pubs)) {
    foreach my $ref (keys(%{$pubs{$dad}})) {
#	print(sprintf("\tSELECT db.name, accession from pub p, pub_dbxref pd, dbxref dx, db where p.pub_id = pd.pub_id and pd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name in ('PMCID','DOI') and p.uniquename = '%s'\n",$ref));
	my $pq2 = $dbh2->prepare(sprintf("SELECT db.name, accession from pub p, pub_dbxref pd, dbxref dx, db where p.pub_id = pd.pub_id and pd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name in ('PMCID','DOI') and p.uniquename = '%s'",$ref));
	$pq2->execute or die "WARNING: ERROR: Unable to execute PMCID/DOI query for ref:\t$ref\n";
	while (my %pr2 = %{$pq2->fetchrow_hashref}) {
	    if ($pr2{name} eq 'PMCID') {
		$pubs{$dad}{$ref}{pmcid} = $pr2{accession};
	    }
	    elsif ($pr2{name} eq 'DOI') {
		$pubs{$dad}{$ref}{doi} = $pr2{accession};
	    }
	    else {
		print "\tWARNING: ERROR: non PMCID/DOI record returned for PMCID/DOI query on: .$ref.\n";
	    }
	}
    }
}

  


# print out the results
print $out "#PUBLICATIONS IN $release WITH  PUB MED IDs\n";
print $out "#PLEASE NOTE: this report first generated for fb_2012_01 release.  All publications associated with a Pub Med ID prior to this release have pmid_added = fb_2011_10\n";
print $out "#FBrf\tPMID\tPMCID\tDOI\tpub_type\tminiref\tpmid_added\n";
print $out "#---------------------------------------------------------------------------------------------------------------------------------------\n";
foreach my $r (sort {$b cmp $a} keys %pubs) {
  foreach my $p (sort keys %{$pubs{$r}}) {
    print  $out "$p\t$pubs{$r}{$p}->{pmid}\t$pubs{$r}{$p}->{pmcid}\t$pubs{$r}{$p}->{doi}\t$pubs{$r}{$p}->{type}\t";
    print $out $pubs{$r}{$p}->{miniref} if $pubs{$r}{$p}->{miniref};
    print $out "\t$pubs{$r}{$p}->{pmid_added}\n";
  }
}

my $end = localtime();
print $out "#START: $start\tEND:$end\n";


sub usage {
  print STDOUT "Usage - perl $0 pgserver db user password outfile\n";
  die;
}
