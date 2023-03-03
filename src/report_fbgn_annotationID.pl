#!/usr/local/bin/perl
# cg_fbgn
#
#       Perl script generates tab-delimited cg_fbgn table with the following 
#       fields: valid_symbol, FBgn#, secondary_FBgn(s), CG#, secondar_CG#(s)
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
# require "../perl_modules/conversions.pm";
if (@ARGV < 5) {
    print "\n USAGE: cg_fbgn pg_server db_name pg_username pg_password output_filename\n\n";
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
print "## FlyBase FBgn-Annotation ID Correspondence Table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print "##gene_symbol\torganism_abbreviation\tprimary_FBgn#\tsecondary_FBgn#(s)\tannotation_ID\tsecondary_annotation_ID(s)\n";


#
## Main method
#
## Get genes...
my $fbgnwc = 'FBgn%';
my $fbalwc = 'FBal%';
## Main driver
my $gq = $dbh->prepare(sprintf("SELECT f.uniquename, f.name, f.feature_id, cvt.name as ftype, abbreviation from feature f, cvterm cvt, organism o where f.organism_id = o.organism_id and f.type_id = cvt.cvterm_id and cvt.name = 'gene' and f.is_obsolete = 'f' and f.is_analysis = 'f' and f.uniquename like '%s' and o.abbreviation = 'Dmel'",$fbgnwc));
$gq->execute or die "WARNING: ERROR: Unable to execute gene query...\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
#    print "\nProcessing gene: $gr{feature_id}\t$gr{uniquename} / $gr{name}...\n";
## Get CG#(s)
    my $primary_cg;
    my @secondary_cgs;
    my @secondary_fbgns;
    my @primcgs; ## An array we shouldnt need, consider it a workaround for now...
    my $cgcnt = 0;
    my $cq = $dbh2->prepare(sprintf("SELECT is_current, accession, db.name from feature_dbxref fd, dbxref d, db where fd.feature_id = %d and fd.dbxref_id = d.dbxref_id and d.db_id = db.db_id and db.name in ('FlyBase Annotation IDs','FlyBase')",$gr{feature_id}));
    $cq->execute or die "WARNING: ERROR: Unable to execute CG query\n";
    my $cq_cnt = $cq->rows;
    if ($cq_cnt > 0) {
	while (my %cr = %{$cq->fetchrow_hashref}) {
	    next if ($cr{accession} =~ /\-/);
## Handle CG#
	    if ($cr{name} eq 'FlyBase Annotation IDs') {
		if ($cr{is_current} == 0) {
		    push(@secondary_cgs,$cr{accession});
#		    print "\tpushing secondary: $cr{accession}\n";
		}
		else {
		    push(@primcgs,$cr{accession});
		    $primary_cg = $cr{accession};
#		    print "\tprimary: $cr{accession}\n";
		    $cgcnt++;
		}
	    }
## Handle FBgn#
	    elsif (($cr{name} eq 'FlyBase') && ($cr{accession} =~ /^FBgn/)) {
#		print "\t\tSecondary fbgn?: $cr{accession} ($cr{is_current})\n";
		if (($cr{is_current} == 0) || ($cr{accession} ne $gr{uniquename})) {
#		    print "\t\t\tPushing secondary FBgn: $cr{accession}\n";
		    push(@secondary_fbgns,$cr{accession});
		}
	    }
	}

## Check for any CG# among synonyms...
	my $sq = $dbh2->prepare(sprintf("SELECT is_current, s.name from feature_synonym fs, synonym s, cvterm st where fs.feature_id = %d and fs.is_current = 'f' and fs.synonym_id = s.synonym_id and s.type_id = st.cvterm_id and st.name = 'symbol'",$gr{feature_id}));
	$sq->execute or die "WARNING: ERROR: Unable to execute synonym query...\n";
	my $sq_cnt = $sq-rows;
	if ($sq_cnt > 0) {
	    while (my %sr = %{$sq->fetchrow_hashref}) {
		if ($sr{name} =~ /^CG\d+$/) {
#		    print "\t\tCG# found in synonyms: $sr{name}\n";
		    my $not_new;
		    foreach my $cg (@secondary_cgs) {
			if ($sr{name} eq $cg) {
			    $not_new++;
			}
		    }
		    if ((!$not_new) && ($sr{name} ne $primary_cg)) {
#			print "\t\t\tPushing this cg#: $sr{name}\n";
			push(@secondary_cgs,$sr{name});
		    }
		}
	    }
	}

	if ($cgcnt > 1) {
## THIS SHOULDNT BE HAPPENNING!  BUT WE RESOLVE WHAT WE CAN HERE...
	    print(sprintf("\tWARNING: ERROR: Multiple current CG# linked to this gene: %s\t%s\t%s\n",$gr{feature_id},$gr{uniquename},$gr{name}));
	    foreach my $pcg (@primcgs) {
		if ($pcg eq $gr{name}) {
		    $primary_cg = $pcg;
		}
		elsif ($pcg ne $primary_cg) {
		    push(@secondary_cgs,$pcg);
		}
	    }
	}
    }
## OUTPUT: gene_symbol, FBgn#, secondary_FBgn#(s), CG#, secondary_CG#(s)    
    if ($primary_cg) {
	print(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n",$gr{name},$gr{abbreviation},$gr{uniquename},join(",",@secondary_fbgns),$primary_cg,join(",",@secondary_cgs)));
    }
}
