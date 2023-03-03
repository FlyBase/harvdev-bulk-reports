#!/usr/local/bin/perl
# report_unique_proteins
#
#       Perl script generates unique protein records for new and updated
#       genes
#       
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
    print "\n USAGE: report_unique_proteins organism_abbreviation pg_server db_name pg_username pg_password output_filename\n\n";
    print "\toutput_filename is the output file for std and error output.\n\n";
    exit;
}

my $insp = shift(@ARGV);
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
print "## FlyBase unique proteins by gene table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n";
## print "## Script: /users/emmert/work/gmod/schema/chado/unique_proteinset/report_unique_proteins\n\n";
print "##FBgn\tFB_gene_symbol\trepresentative_protein\tidentical_protein(s)\n\n";


#
## Main methods
#
## Start with query for all FBgn and for each of these, get all associated FBtr & FBpp
my $fbgnwc = 'FBgn%';
## Main driver
my $gq = $dbh->prepare(sprintf("SELECT uniquename, f.feature_id, genus, species, f.name from feature f, cvterm cvt, organism o where f.type_id = cvt.cvterm_id and cvt.name = 'gene' and f.organism_id = o.organism_id and o.abbreviation = '%s' and uniquename like '%s' and f.is_obsolete = 'f'",$insp,$fbgnwc));
## Test driver (limited query)
# my $gq = $dbh->prepare(sprintf("SELECT uniquename, f.feature_id, genus, species, f.name from feature f, cvterm cvt, organism o where f.type_id = cvt.cvterm_id and cvt.name = 'gene' and f.organism_id = o.organism_id and o.abbreviation = '%s' and uniquename like '%s' and f.is_obsolete = 'f' limit 500",$insp,$fbgnwc));
$gq->execute or die "WARNING: ERROR: Unable to execute evidence query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
    my %gprots;
#    print "\nProcessing gene: $gr{uniquename}\t$gr{name}\t$gr{genus}\t$gr{species}\n";
## Query for transcript(s)
    my $tq = $dbh2->prepare(sprintf("SELECT frtp.name, t.uniquename as fbtr, t.name, t.feature_id, ttp.name as ttype from feature t, feature_relationship frt, cvterm ttp, cvterm frtp, featureloc fl where frt.object_id = %d and frt.type_id = frtp.cvterm_id and frtp.name = 'partof' and frt.subject_id = t.feature_id and t.type_id = ttp.cvterm_id and t.is_obsolete = 'f' and t.feature_id = fl.feature_id order by t.name",$gr{feature_id}));
    $tq->execute or die "WARNING: ERROR: Unable to execute trancsript query\n";
    my $tq_cnt = $tq->rows;
    if ($tq_cnt > 0) {
	while (my %tr = %{$tq->fetchrow_hashref}) {
	    if ($tr{ttype} =~ /mRNA|tRNA|miRNA|ncRNA|pre\-miRNA|rRNA|snRNA|snoRNA|pseudogene/) {
## Query for protein associated with each transcript...
#		print "\tProcessing transcript: $tr{fbtr}\t$tr{name}\t$tr{ttype}\n";
		my $isamatch = 0;
		my $pq = $dbh3->prepare(sprintf("SELECT p.uniquename as fbpp, p.name, residues from feature p, feature_relationship frp, cvterm ptp, featureloc fl where frp.object_id = %d and frp.subject_id = p.feature_id and p.type_id = ptp.cvterm_id and ptp.name = 'polypeptide' and p.is_obsolete = 'f' and p.is_analysis = 'f' and p.feature_id = fl.feature_id",$tr{feature_id}));
		$pq->execute or die "WARNING: ERROR: Unable to execute protein query\n";
		my $pq_cnt = $pq->rows;
		if ($pq_cnt > 0) {
		    while (my %pr = %{$pq->fetchrow_hashref}) {
## When a protein is found, compare it to other proteins from this gene.  Proteins which are identical 
## to proteins found earlier will be listed as "identical_protein(s)" to the earlier protein
#			print "\t\tprotein found: $pr{fbpp}\t$pr{name}\n";
			foreach my $prot (keys(%gprots)) {
#			  print "\t\t\tComparing...\n\t\t\t$pr{residues}\n\t\t\t$gprots{$prot}{residues}\n";
			    if ($pr{residues} eq $gprots{$prot}{residues}) {
				push(@{$gprots{$prot}{matches}},$pr{name});
#				print "\t\t\t$pr{fbpp} matches $prot (pushing)...\n";
				$gprots{$pg{name}}{unique}--;
				$isamatch++;
				last;
			    }
			}
			if (!$isamatch) {
			    $gprots{$pr{name}}{residues} = $pr{residues};
			    $gprots{$pr{name}}{unique}++;
#			    print "\t\t\t$pr{fbpp} is unique (maybe)...\n";
			}
		    }
		}
	    }
	}
    }
    foreach my $prot (sort(keys(%gprots))) {
	if ((@{$gprots{$prot}{matches}}) || ($gprots{$prot}{unique} > 0)) {
	    print(sprintf("%s\t%s\t%s\t%s\n",$gr{uniquename},$gr{name},$prot,join(',',@{$gprots{$prot}{matches}})));
	}
    }
}

$jetzt = scalar localtime;
print "\n\n## Finished report_unique_proteins: $jetzt\n";
