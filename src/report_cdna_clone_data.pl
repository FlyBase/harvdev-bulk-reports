#!/usr/local/bin/perl
# report_cdna_clone_data
#
#       Perl script generates tab-delimited report of cDNA clone data 
#       containing the following fields: FBcl#, clone_name, library_name, 
#       cDNA_acc#, EST_acc#
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
    print "\n USAGE: report_cdna_clone_data pg_server db_name pg_username pg_password output_filename\n\n";
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
#  $dbh = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh2 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh3 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";


#
##  General setup
#
## Setup file header
$jetzt = scalar localtime;
print "## FlyBase cDNA clone data table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print "##FBcl#\torganism_abbreviation\tclone_name\tdataset_metadata_name\tcDNA_accession(s)\tEST_accession(s)\n";

#
## Main method
#
## Get clones...

my $fbclwc = 'FBcl%';
## Main query
my $cq = $dbh->prepare(sprintf("SELECT f.uniquename, f.name, f.feature_id, cvt.name as ctype, abbreviation from feature f, organism o, cvterm cvt where f.organism_id = o.organism_id and f.type_id = cvt.cvterm_id and cvt.name = 'cDNA_clone' and uniquename like '%s' and f.is_obsolete = 'f' and is_analysis = 'f'",$fbclwc));
$cq->execute or die "WARNING: ERROR: Unable to execute clone query\n";
while (my %cr = %{$cq->fetchrow_hashref}) {
#    print "\nProcessing clone: $cr{feature_id}\t$cr{uniquename}\t$cr{name}\t$cr{ctype}\n";
## Get clone library (if exists; there are a few clones w/out libraries...)
    my $lib;
    my $lq = $dbh2->prepare(sprintf("SELECT * from library l, library_feature lf where l.library_id = lf.library_id and lf.feature_id = %d",$cr{feature_id}));
    $lq->execute or die "WARNING: ERROR: Unable to execute library query\n";
    my $lq_cnt = $lq->rows;
    if ($lq_cnt > 0) {
	while (my %lr = %{$lq->fetchrow_hashref}) {
#	    print "\tlibrary: $lr{uniquename}\n";
	    $lib = $lr{name};
	}
    }
## Get associated cDNA(s)
    my @cdnas;
    my @ests;
    my $aq = $dbh3->prepare(sprintf("SELECT f.feature_id, f.uniquename, f.name, f.type_id, cvt.name as ctype from feature f, feature_relationship fr, cvterm cvt where fr.subject_id = f.feature_id and f.is_obsolete = 'f' and f.type_id = cvt.cvterm_id and cvt.name in ('EST','cDNA') and fr.object_id = %d",$cr{feature_id}));
    $aq->execute or die "WARNING: ERROR: Unable to execute cDNA query\n";
    my $aq_cnt = $aq->rows;
    if ($aq_cnt > 0) {
	while (my %ar = %{$aq->fetchrow_hashref}) {
#	    print "\tcDNA found: $ar{feature_id}, $ar{uniquename}\t$ar{name}\t$ar{ctype}\n";
	    if (($ar{ctype} eq 'cDNA') && ($ar{uniquename} !~ /\:contig\d*$/)) {
		push(@cdnas,$ar{uniquename});
	    }
	    elsif (($ar{ctype} eq 'EST')  && ($ar{uniquename} !~ /\:contig\d*$/)) {
		push(@ests,$ar{uniquename});
	    }
	}
    }
    else {
#	print "\tNo cDNAs associated with this clone?\n";
    }
## output   
    print(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n",$cr{uniquename},$cr{abbreviation},$cr{name},$lib,join(',',@cdnas),join(',',@ests)));
}

$jetzt = scalar localtime;
print "## Finished cDNA clone data table: $jetzt\n";
