#!/usr/local/bin/perl

# report_allele_phenotype
#
#       Perl script generates allele phenotype table 
#       with the following fields:
#       allele_symbol, allele_FBid, phenotype, FBrf#
#       
#
#-----------------------------------------------------------------------------#
#
#	NOTES
#	
#	Values are taken from the following derived fields:
#	   derived_pheno_class
#	   derived_pheno_manifest
#
#-----------------------------------------------------------------------------#
use DBI;
require "conversions.pm";

if (@ARGV < 5) {
    print "\n USAGE: report_allele_phenotype pg_server db_name pg_username pg_password output_filename\n\n";
    print "\toutput_filename is the output file for std and error output.\n\n";
    exit;
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);

$OUTFILE = shift(@ARGV);
# open(STDERR,">>$OUTFILE");
open(STDOUT,">>$OUTFILE");

#
## DB connection
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

## Array of derived featureprops for allele phenotype data
my @dfprops = ('derived_pheno_class','derived_pheno_manifest');

## allele FBal wildcard (for querying)
my $alwc = 'FBal%';

## hash for data to output
my %ah;

## Setup header
$jetzt = scalar localtime;
print "\n## FlyBase Allele Phenotype Table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource\n\n";
print "##allele_symbol\tallele_FBal#\tphenotype\tFBrf#\n";


#
# Main method
#
## Main driver
foreach my $fprop (@dfprops) {
##  Main driver
    my $aq = $dbh->prepare(sprintf("SELECT f.uniquename as fbal, f.name as alsym, fp.value, p.uniquename as fbrf from feature f, featureprop fp, cvterm cvt, featureprop_pub fpp, pub p where f.feature_id = fp.feature_id and fp.type_id = cvt.cvterm_id and cvt.name = '%s' and fp.featureprop_id = fpp.featureprop_id and fpp.pub_id = p.pub_id and f.uniquename like '%s' and f.is_obsolete = 'f'",$fprop,$alwc));
## Test driver (Limited set)
#     my $aq = $dbh->prepare(sprintf("SELECT f.uniquename as fbal, f.name as alsym, fp.value, p.uniquename as fbrf from feature f, featureprop fp, cvterm cvt, featureprop_pub fpp, pub p where f.feature_id = fp.feature_id and fp.type_id = cvt.cvterm_id and cvt.name = '%s' and fp.featureprop_id = fpp.featureprop_id and fpp.pub_id = p.pub_id and f.uniquename like '%s' and f.is_obsolete = 'f' and f.uniquename < 'FBal0001000'",$fprop,$alwc));
    $aq->execute or die "WARNING: ERROR: Unable to execute allele interaction query ($fprop)\n";
    my $aq_cnt = $aq->rows;
    if ($aq_cnt > 0) {
	while (my %ar = %{$aq->fetchrow_hashref}) {
	    $ar{value} =~ s/\@//g;
	    $ar{value} =~ s/FB\D{2}\d{7}\://g;
	    $ar{value} =~ s/FB\D{2}\d{8}\://g;
	    $ar{value} = &decon($ar{value});
#	    print "$ar{alsym}\t$ar{fbal}\t$ar{value}\t$ar{fbrf}\n";
	    my $t1 = sprintf("$ar{alsym}\t$ar{fbal}");
	    my $t2 = sprintf("$ar{value}\t$ar{fbrf}");
	    push(@{$ah{$t1}},$t2);
	}
    }
}

## Print OUTPUT
foreach my $al (sort(keys(%ah))) {
    my $last_in;
    foreach my $in (sort(@{$ah{$al}})) {
	next if ($in eq $last_in);
	print "$al\t$in\n";
	$last_in = $in;
    }
}

$jetzt = scalar localtime;
print "\n## Finished report_allele_phenotype: $jetzt\n\n";
