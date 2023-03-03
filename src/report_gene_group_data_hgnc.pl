#!/usr/bin/perl
# report_gene_group_data_hgnc
#
#       Perl script reports FlyBase FBgg gene groups and associated HGNC Gene 
#       Family IDs (if any) for report to HGNC.
#       
#       See below for reported fields.
#
#-----------------------------------------------------------------------------#
#
#       NOTES
#       
#       See JIRA DB-266 for more details
#       
#       This script is a ripoff of report_gene_group_data 
#
#-----------------------------------------------------------------------------#
use DBI;

if (@ARGV != 5) {
  print "\n USAGE: report_gene_group_data_hgnc pg_server db_name pg_username pg_password report_output_filename\n\n";
  exit();
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);

$OUTFILE = shift(@ARGV);
open(STDERR,">$OUTFILE");
open(STDOUT,">$OUTFILE");


#
##  DB Connections
#

# ## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


#
##  General setup
#

## Print header
$jetzt = scalar localtime;
print "## FlyBase Gene Groups to HGNC Gene Families correspondence table\n## Generated: $jetzt\n## Using chado datasource: $dsource\n\n";


print "## The absence of an HGNC_family_ID entry indicates there is no equivalent HGNC gene family for that FlyBase gene group.\n";
print "## Because of different sub-group structures (etc), a single HGNC family may be associated with multiple FlyBase gene groups.\n";
print "## Similarly, a single FlyBase gene group may be associated with multiple HGNC gene families - these are shown on separate lines.\n\n";

print "## FB_group_id\tFB_group_symbol\tFB_group_name\tHGNC_family_ID\n";


#
## Main method
#

## Populate hash %gnh of all grp grp_id, uniquename, name, and gene(s) if any
my %gnh;
my $grq = $dbh3->prepare("SELECT grp_id, name, uniquename from grp where is_obsolete = 'f'");
$grq->execute or die "WARNING: ERROR: Unable to execute grp name/uniquename query\n";
while (my %grr = %{$grq->fetchrow_hashref}) {
#    print "Adding grp to gnh:\t$grr{grp_id}\t$grr{name}\t$grr{uniquename}\n";
    $gnh{$grr{grp_id}}{name} = $grr{name};
    $gnh{$grr{grp_id}}{uniquename} = $grr{uniquename};
}

## Add fullname for each grp record to hash %gnh 
foreach $rec (keys(%gnh)) {
    my $sq = $dbh4->prepare(sprintf("SELECT s.name as sname from grp_synonym gs, synonym s, cvterm cvt where gs.synonym_id = s.synonym_id and gs.is_current = 't' and gs.is_internal = 'f' and s.type_id = cvt.cvterm_id and cvt.name = 'fullname' and gs.grp_id = %d",$rec));
    $sq->execute or die "WARNING: ERROR: Unable to execute grp_synonym query\n";
    my $sq_cnt = $sq->rows;
    if ($sq_cnt > 0) {
	if ($sq_cnt > 1) {
	    print "\tWARNING: ERROR: Multiple fullname synonyms found for grp: $rec\t$gnh{$rec}{uniquename}\t$gnh{$rec}{name}\n";
	}
	while (my %sr = %{$sq->fetchrow_hashref}) {
#	    print "\tSetting fullname for grp $rec\t$gnh{$rec}{uniquename}\t$gnh{$rec}{name}:\t$sr{sname}\n";
	    $gnh{$rec}{fullname} = $sr{sname};
	}
    }
}

## Add associated HGNC accession(s) (if any) to grp records in %gnh
foreach my $rec (keys(%gnh)) {
    my $hq = $dbh2->prepare(sprintf("SELECT db.name as dbname, dx.accession from grp_dbxref gdx, dbxref dx, db where gdx.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name = 'HGNC-GG1' and gdx.grp_id = %d",$rec));
    $hq->execute or die "WARNING: ERROR: Unable to execute HGNC query\n";
    my $hq_cnt = $hq->rows;
    if ($hq_cnt > 0) {
	while (my %hr = %{$hq->fetchrow_hashref}) {
	    $gnh{$rec}{hgnc}{$hr{accession}} = 't';
#	    print "\t\tAdding hgnc for grp $rec\t$gnr{grpuname}\t$gnr{grpname}:\t$hr{accession}\n";
	}
    }
}


#
## Output report
#

## For each group
foreach my $p (keys(%gnh)) {   
#    print "\tGetting hgnc for grp: .$p.\n";
    if ($gnh{$p}{hgnc}) { ## And report associated HGNC
	foreach my $g (keys(%{$gnh{$p}{hgnc}})) {
#	    print "REPORT\t$p\t$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t$g\n";
	    print "$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t$g\n";
	}
    }
    else {
#	print "REPORT_NOHGNC\t$p\t$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t\n";
	print "$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t\n";
    }
}


$jetzt = scalar localtime;
print "\n\n## Finished HGNC Gene Groups report:\t$jetzt\n";
exit();
