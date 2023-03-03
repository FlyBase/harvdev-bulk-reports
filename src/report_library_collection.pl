#!/usr/local/bin/perl
# report_library_collection
#
#       Perl script generates tab-delimited library_collection table with the 
#       following fields: Library/Collection_ID, Library/Collection_Name, 
#       Item_ID, Item_Name
#       
#       See Roberts email: Subject: Library/Collection Pre-computed Files
#       Date: Fri, 26 Jun 2009 13:37:08
#
#-----------------------------------------------------------------------------#
#
#	NOTES
#	
#	
#
#-----------------------------------------------------------------------------#
use DBI;

if (@ARGV < 5) {
    print "\n USAGE: report_library_collection pg_server db_name pg_username pg_password output_filename\n\n";
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

## Build a hash of feature types and their ids
my $ftq = $dbh->prepare("SELECT cvt.name as termname, cvterm_id from cvterm cvt, cv c where c.cv_id = cvt.cv_id and c.name = 'SO'");
$ftq->execute or die "WARNING: ERROR: Cannot get feature types\n";
while (my %ftp = %{$ftq->fetchrow_hashref}) {
  $ftype{$ftp{cvterm_id}} = $ftp{termname};
  $ftype2{$ftp{termname}} = $ftp{cvterm_id};
}


## Setup file header
$jetzt = scalar localtime;
print "## FlyBase Dataset Metadata Report\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print "##Dataset_Metadata_ID, Dataset_Metadata_Name, Item_ID, Item_Name\n\n";


#
## Main methods
#

## SQL query (provided by Robert)

## Main driver
## my $lq = $dbh->prepare("SELECT lc.uniquename as luname, lc.name as lname, f.uniquename as funame, f.name as fname, f.feature_id as fid, cvt.name as ctype from feature f, library_feature lf, library lc, cvterm cvt where lc.type_id = cvt.cvterm_id and lc.is_obsolete = 'f' and lf.library_id=lc.library_id and lf.feature_id = f.feature_id and f.is_obsolete = 'f' and f.is_analysis = 'f' order by f.uniquename, f.name desc");
my $lq = $dbh->prepare("SELECT lc.uniquename as luname, lc.name as lname, f.uniquename as funame, f.name as fname, f.feature_id as fid, cvt.name as ctype, cvt2.name as ftype from feature f, library_feature lf, library lc, cvterm cvt, cvterm cvt2 where lc.type_id = cvt.cvterm_id and lc.is_obsolete = 'f' and lf.library_id=lc.library_id and lf.feature_id = f.feature_id and f.is_obsolete = 'f' and f.type_id = cvt2.cvterm_id order by f.uniquename, f.name desc");
$lq->execute or die "WARNING: ERROR: Unable to execute library/collection query\n";
while (my %lr = %{$lq->fetchrow_hashref}) {
    next if ($lr{ctype} eq 'predicted_gene_model_set');
    next if ($lr{ftype} eq 'match');
    if ($lr{ctype} eq 'gene expression cluster') { 
    my $gq = $dbh2->prepare(sprintf("SELECT g.uniquename, g.name from cvterm cvt, cvterm cvt2, feature g, feature_relationship fr where fr.subject_id = %d and fr.object_id = g.feature_id and g.type_id = cvt2.cvterm_id and cvt2.name = 'gene' and g.is_obsolete = 'f' and fr.type_id = cvt.cvterm_id and cvt.name = 'associated_with'",$lr{fid}));
    $gq->execute or die "WARNING: ERROR: Unable to execute XR gene query on: $lr{fid}\t$lr{funame}\t$lr{fname}\n";
    my $gq_cnt = $gq->rows;
    if ($gq_cnt > 0) {
      while (my %gr = %{$gq->fetchrow_hashref}) {
	  print "$lr{luname}\t$lr{lname}\t$gr{uniquename}\t$gr{name}\n";
      }
    }
    else {
      print "\tWARNING: ERROR: Gene query on XR transcript returns no rows: $lr{fid}\t$lr{funame}\t$lr{fname}\n";
    }
  }
  else {
    print "$lr{luname}\t$lr{lname}\t$lr{funame}\t$lr{fname}\n";
  }
}


$jetzt = scalar localtime;
print "\n\n## Finished report_library_collection: $jetzt\n";
