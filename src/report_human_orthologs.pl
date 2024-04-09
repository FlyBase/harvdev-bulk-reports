#!/usr/local/bin/perl
# report_human_orthologs
#
#       Perl script generates tab-delimited report_human_orthologs table with the 
#       following fields: 
#	
#	Dmel_gene_ID
#	Dmel_gene_symbol
#	Human_gene_HGNC_ID
#	Human_gene_OMIM_ID
#	Human_gene_symbol
#	DIOPT_score
#	OMIM_Phenotype_ID  (name [,name,name,...])
#
#-----------------------------------------------------------------------------#
#
#	NOTES:
#	          See JIRA DB-271 for details
#	
#
#-----------------------------------------------------------------------------#
use DBI;
use POSIX;
require "conversions.pm";

if (@ARGV != 5) {
    print "\n USAGE: report_human_orthologs pg_server db_name pg_username pg_password report_output_filename\n\n";
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
## Data source (g4)
#my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
#$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
#$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
#$dbh5 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


#
##  General setup
#


## Setup file header
$jetzt = scalar localtime;
print  "## FlyBase D. melanogaster-Human Orthologs and associated diseases report\n## Generated: $jetzt\n";
print  "## Using datasource: $dsource...\n\n";
print  "##Dmel_gene_ID\tDmel_gene_symbol\tHuman_gene_HGNC_ID\tHuman_gene_OMIM_ID\tHuman_gene_symbol\tDIOPT_score\tOMIM_Phenotype_IDs\tOMIM_Phenotype_IDs[name]\n";


my $dioptwc = 'diopt%';

#
## Main methods
#


my @pouts;
## Get Dmel genes w/ Hsap orthologs
my $oq = $dbh->prepare("SELECT f.feature_id as fid, f.uniquename as funame, f.name as fname, h.feature_id as hid, h.uniquename as huname, h.name as hname, fr.value from feature f, feature h, feature_relationship fr, cvterm cvt, organism ho, organism fo where f.is_obsolete = 'f' and f.organism_id = fo.organism_id and fo.abbreviation = 'Dmel' and f.feature_id = fr.subject_id and fr.object_id = h.feature_id and h.is_obsolete = 'f' and h.organism_id = ho.organism_id and ho.abbreviation = 'Hsap' and fr.type_id = cvt.cvterm_id and cvt.name = 'orthologous_to' UNION SELECT f.feature_id as fid, f.uniquename as funame, f.name as fname, h.feature_id as hid, h.uniquename as huname, h.name as hname, fr.value from feature f, feature h, feature_relationship fr, cvterm cvt, organism ho, organism fo where f.is_obsolete = 'f' and f.organism_id = fo.organism_id and fo.abbreviation = 'Dmel' and h.feature_id = fr.subject_id and fr.object_id = f.feature_id and h.is_obsolete = 'f' and h.organism_id = ho.organism_id and ho.abbreviation = 'Hsap' and fr.type_id = cvt.cvterm_id and cvt.name = 'orthologous_to'");
$oq->execute or die "WARNING: ERROR: Unable to execute fly-hsap query\n";
while (my %or = %{$oq->fetchrow_hashref}) {
#    print "\nProcessing ortholog pair:\t$or{fid}\t$or{funame}\t$or{fname}\t$or{hid}\t$or{huname}\t$or{hname}\t";
    my %vhash;
    my @vlines = split(/\n/,$or{value});
    foreach my $vline (@vlines) {
	if ($vline =~ /(^.+)\:\s(.*)/) {
	    $vhash{$1} = $2;
	}
    }
#    print "\tmethods:\t.$vhash{'diopt methods'}.\n";
    my $dscore = (split(",",$vhash{'diopt methods'}));
#    print "\tdscore:\t$dscore\n";
    if ($dscore > 2) {  # skip orthologs w/ diopt score < 3
	my %ohash;
## Get accno associated with this Hsap ortholog here...	
	my $xq = $dbh2->prepare(sprintf("SELECT db.name, dx.* from feature_dbxref fd, dbxref dx, db where fd.feature_id = %d and fd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and accession not like '%s'",$or{hid},$dioptwc));
	$xq->execute or die "WARNING: ERROR: Unable to execute xref query\n";
	my $xq_cnt = $xq->rows;
	if ($xq_cnt > 0) {
	    while (my %xr = %{$xq->fetchrow_hashref}) {
#		print "\tFound xref:\t$xr{name}\t$xr{accession}\t$xr{description}\n";
		if ($xr{name} eq 'HGNC') {
		    $ohash{$or{hname}}{HGNC} = $xr{accession};
		}
		if ($xr{name} eq 'OMIM_GENE') {
		    $ohash{$or{hname}}{OMIM_GENE} = $xr{accession};
		}
		if ($xr{name} eq 'OMIM_PHENOTYPE') {
		    my $ophef = sprintf("%d[%s]",$xr{accession},$xr{description});
		    push(@{$ohash{$or{hname}}{OMIM_PHENO_FULL}},$ophef);
		    push(@{$ohash{$or{hname}}{OMIM_PHENO_IDS}},$xr{accession});
		}
	    }
## Now print out report...
	    if ($ohash{$or{hname}}{OMIM_GENE}) {  # Sometimes theres not an OMIM_GENE record, e.g.: FBog0000403309 / Hsap\PROP1
		my $hsymb = $or{hname};
		$hsymb =~ s/Hsap\\//g;
		push(@pouts,sprintf("%s\t%s\tHGNC:%d\tMIM:%d\t%s\t%d\t%s\t%s\n",$or{funame},$or{fname},$ohash{$or{hname}}{HGNC},$ohash{$or{hname}}{OMIM_GENE},$hsymb,$dscore,join(',',sort(@{$ohash{$or{hname}}{OMIM_PHENO_IDS}})),join(',',sort(@{$ohash{$or{hname}}{OMIM_PHENO_FULL}}))));
	    }
	    else {
		my $hsymb = $or{hname};
		$hsymb =~ s/Hsap\\//g;
		push(@pouts,sprintf("%s\t%s\tHGNC:%d\t\t%s\t%d\t%s\t%s\n",$or{funame},$or{fname},$ohash{$or{hname}}{HGNC},$hsymb,$dscore,join(',',sort(@{$ohash{$or{hname}}{OMIM_PHENO_IDS}})),join(',',sort(@{$ohash{$or{hname}}{OMIM_PHENO_FULL}}))));
	    }
	}
	else {
	    print "\tWARNING: ERROR(?): No xrefs found for Hsap gene:\t$or{hid}\t$or{huname}\t$or{hname}\n";
	}
    }
}

## Uniquify & print report
foreach my $pout (&uniq(@pouts)) {
    print "$pout";
}



$jetzt = scalar localtime;
print "## Finished report_human_orthologs: $jetzt\n\n";

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

