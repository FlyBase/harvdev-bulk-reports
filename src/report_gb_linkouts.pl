#!/usr/local/bin/perl

# report_gb_linkouts
#
#       Perl script produces XML reports for GenBank LinkOuts:
#	
#	
#       - nt-flybase.xml is LinkOut from Entrez Nucleotide records to FB Gene 
#         reports    
#
#       - pr-flybase.xml is LinkOut from Entrez Protein records to FB Gene
#         reports
#
#       - pm-flybase.xml is LinkOut from PubMed records to FB Publication
#         reports
#	
#	
#	
#-----------------------------------------------------------------------------#
#
#	NOTES:
#
#   -Before putting to GB, validate XML at:
#    http://www.ncbi.nlm.nih.gov/entrez/linkout/doc/validate.shtml
#
#   -FTP files nt-flybase.xml, pr-flybase.xml, and pm-flybase.xml to:
#    Host:  ftp-private.ncbi.nih.gov
#    User: flybase
#    Password: NgJ<l2's
#    directory: holdings/
#	
#	
#	
#	
#	
#	
#	
#
#-----------------------------------------------------------------------------#
use DBI;


if (@ARGV < 6) {
    print "\n USAGE: report_gb_linkouts pg_server db_name pg_username pg_password linkout_file_output_directory debug_output_filename\n\n";
    exit;
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);
my $outdir = shift(@ARGV);

$OUTFILE = shift(@ARGV);
open(STDERR,">>$OUTFILE");
open(STDOUT,">>$OUTFILE");

#
## Set up DB connections
#
## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";

## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
# $dbh = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";

$jetzt = scalar localtime;
print "\nStarting report_gb_linkouts $jetzt\n";
print "Using datasource: $dsource\n";


#
##  General setup
#

## Set up output XML files
## pr-flybase.xml: LinkOut from Entrez Protein records to FB Gene reports
my $prout = $outdir . 'pr-flybase.xml';
#open(PR,"> pr-flybase.xml") or die "Cannot open PR output file!\n";
open(PR,"> $prout") or die "Cannot open PR output file: $prout!\n";
print(PR "<!DOCTYPE LinkSet PUBLIC \"-//NLM//DTD LinkOut //EN\" \"LinkOut.dtd\"\n");
print(PR "\[\n");
# print(PR "<!ENTITY icon.url \"http://www.ncbi.nlm.nih.gov/PMGifs/Toolbar/flybase.gif\">\n");
# print(PR "<!ENTITY icon.url \"http://flybase.org/images/fly_logo.png\">\n");
print(PR "<!ENTITY base.url \"http://flybase.org/cgi-bin/uniq.html?db=fbgn&amp;field=cross_references&amp;context=\">\n");
print(PR "\]>\n");
print(PR "<LinkSet>\n\n");
print(PR "\t<Link>\n");
print(PR "\t\t<LinkId>1</LinkId>\n");
print(PR "\t\t<ProviderId>2006</ProviderId>\n");
print(PR "\t\t<IconUrl>&icon.url;</IconUrl>\n");
print(PR "\t\t<ObjectSelector>\n");
print(PR "\t\t\t<Database>Protein</Database>\n");
print(PR "\t\t\t<ObjectList>\n");
$gimls = "\t\t\t\t<Query>";  # will go at the start of PR XML objects printed
# $gimle = " </Query>";  # will go at the end of PR XML objects printed
$gimle = " [pacc]</Query>";  # will go at the end of PR XML objects printed

## nt-flybase.xml: LinkOut from Entrez Nucleotide records to FB Gene reports
my $ntout = $outdir . 'nt-flybase.xml';
#open(NT,"> nt-flybase.xml") or die "Cannot open NT output file!\n";
open(NT,"> $ntout") or die "Cannot open NT output file: $ntout!\n";
print(NT "<!DOCTYPE LinkSet PUBLIC \"-//NLM//DTD LinkOut //EN\" \"LinkOut.dtd\"\n");
print(NT "\[\n");
# print(NT "<!ENTITY icon.url \"http://www.ncbi.nlm.nih.gov/PMGifs/Toolbar/flybase.gif\">\n");
# print(NT "<!ENTITY icon.url \"http://flybase.org/images/fly_logo.png\">\n");
print(NT "<!ENTITY base.url \"http://flybase.org/cgi-bin/uniq.html?db=fbgn&amp;field=cross_references&amp;context=\">\n");
print(NT "\]>\n");
print(NT "<LinkSet>\n\n");
print(NT "\t<Link>\n");
print(NT "\t\t<LinkId>1</LinkId>\n");
print(NT "\t\t<ProviderId>2006</ProviderId>\n");
print(NT "\t\t<IconUrl>&icon.url;</IconUrl>\n");
print(NT "\t\t<ObjectSelector>\n");
print(NT "\t\t\t<Database>Nucleotide</Database>\n");
print(NT "\t\t\t<ObjectList>\n");
$acmls = "\t\t\t\t<Query>";  # will go at the start of NT XML objects printed
# $acmle = " </Query>";  # will go at the end of NT XML objects printed
$acmle = " [pacc]</Query>";  # will go at the end of NT XML objects printed

## pm-flybase.xml: LinkOut from PubMed records to FB Publication reports
my $pmout = $outdir . 'pm-flybase.xml';
#open(PM,"> pm-flybase.xml") or die "Cannot open output file!\n";
open(PM,"> $pmout") or die "Cannot open PM output file: $pmout!\n";
print(PM "<!DOCTYPE LinkSet PUBLIC \"-//NLM//DTD LinkOut //EN\" \"LinkOut.dtd\"\n");
print(PM "\[\n");
# print(PM "<!ENTITY icon.url \"http://www.ncbi.nlm.nih.gov/PMGifs/Toolbar/flybase.gif\">\n");
# print(PM "<!ENTITY icon.url \"http://flybase.org/images/fly_logo.png\">\n");
print(PM "<!ENTITY base.url \"http://flybase.org/cgi-bin/uniq.html?db=fbrf&amp;field=pubmed_id&amp;context=\">\n");
print(PM "\]>\n");
print(PM "<LinkSet>\n\n");
print(PM "\t<Link>\n");
print(PM "\t\t<LinkId>1</LinkId>\n");
print(PM "\t\t<ProviderId>2006</ProviderId>\n");
print(PM "\t\t<IconUrl>&icon.url;</IconUrl>\n");
print(PM "\t\t<ObjectSelector>\n");
print(PM "\t\t\t<Database>PubMed</Database>\n");
print(PM "\t\t\t<ObjectList>\n");
$pmmls = "\t\t\t\t<ObjId>";  # will go at the start of PM XML objects printed
$pmmle = "</ObjId>";         # will go at the end of PM XML objects printed

#
## Main methods
#

#
##  Get PubMed IDs...
#
print "\nGenerating PubMed LinkOuts...\n";
## Main driver for pub query...
my $pq = $dbh->prepare("SELECT * from pub p, pub_dbxref pd, dbxref dx, db  where db.name = 'pubmed' and db.db_id = dx.db_id and dx.dbxref_id = pd.dbxref_id and pd.pub_id = p.pub_id and p.is_obsolete = 'f'");
## Test driver (limited query)
## my $pq = $dbh->prepare("SELECT * from pub p, pub_dbxref pd, dbxref dx, db  where db.name = 'pubmed' and db.db_id = dx.db_id and dx.dbxref_id = pd.dbxref_id and pd.pub_id = p.pub_id and p.is_obsolete = 'f' limit 500");
$pq->execute or die "WARNING: ERROR: Unable to execute pub/pubmed query\n";
while (my %pr = %{$pq->fetchrow_hashref}) {
    print PM sprintf("%s%s%s\n",$pmmls,$pr{accession},$pmmle);
}


#
##  Get Nucleotide & Protein IDs...
#
print "\nGenerating Nucleotide and Protein LinkOuts...\n";
my @acs;
my @prs;
my $fbgnwc = 'FBgn%';
my $gq = $dbh->prepare(sprintf("SELECT accession, db.name from feature f, feature_dbxref fd, dbxref dx, db where f.feature_id = fd.feature_id and fd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name in ('GB','GB_protein') and f.uniquename like '%s'",$fbgnwc));
## Test driver (limited query)
## my $gq = $dbh->prepare(sprintf("SELECT accession, db.name from feature f, feature_dbxref fd, dbxref dx, db where f.feature_id = fd.feature_id and fd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name in ('GB','GB_protein') and f.uniquename like '%s' limit 500",$fbgnwc));
$gq->execute or die "WARNING: ERROR: Unable to execute GB accession query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
    my $gracc = $gr{accession};
    if ($gr{accession} =~ /\.\d+/) {
	$gracc =~ s/\.\d+//;  ## strip out .versions (shouldnt be in dx.accession anyway!)
    }
    if ($gr{name} eq 'GB') {
#	push(@acs,$gr{accession});
	push(@acs,$gracc);
    }
    else {
#	push(@prs,$gr{accession});
	push(@prs,$gracc);
    }
}
## Make unique & output
my $last_ac;
foreach my $ac (sort(@acs)) {
    if ($ac ne $last_ac) {
	print NT sprintf("%s%s%s\n",$acmls,$ac,$acmle);
    }
    $last_ac = $ac;
}

my $last_pr;
foreach my $pr (sort(@prs)) {
    if ($pr ne $last_pr) {
	print PR sprintf("%s%s%s\n",$gimls,$pr,$gimle);
    }
}


#
## Finish off input files...
#
# PM LinkOut
print(PM "\t\t\t</ObjectList>\n");
print(PM "\t\t</ObjectSelector>\n");
print(PM "\t\t<ObjectUrl>\n");
print(PM "\t\t\t<Base>&base.url;</Base>\n");
print(PM "\t\t\t<Rule>&lo.id;</Rule>\n");
print(PM "\t\t</ObjectUrl>\n");
print(PM "\t</Link>\n");
print(PM "</LinkSet>\n");

# PR LinkOut
print(PR "\t\t\t</ObjectList>\n");
print(PR "\t\t</ObjectSelector>\n");
print(PR "\t\t<ObjectUrl>\n");
print(PR "\t\t\t<Base>&base.url;</Base>\n");
print(PR "\t\t\t<Rule>&lo.pacc;</Rule>\n");
##print(PR "\t\t\t<Rule>&lo.id;</Rule>\n");
print(PR "\t\t</ObjectUrl>\n");
print(PR "\t</Link>\n");
print(PR "</LinkSet>\n");

# NT Linkout
print(NT "\t\t\t</ObjectList>\n");
print(NT "\t\t</ObjectSelector>\n");
print(NT "\t\t<ObjectUrl>\n");
print(NT "\t\t\t<Base>&base.url;</Base>\n");
print(NT "\t\t\t<Rule>&lo.pacc;</Rule>\n");
##print(NT "\t\t\t<Rule>&lo.id;</Rule>\n");
print(NT "\t\t</ObjectUrl>\n");
print(NT "\t</Link>\n");
print(NT "</LinkSet>\n");


$jetzt = scalar localtime;
print "\nFinished report_gb_linkouts $jetzt\n\n";
