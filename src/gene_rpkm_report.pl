#!/usr/bin/perl
# gene_rpkm_report
#
#       Perl script produces gene - RPKM report
#       
#    Fields are:   
#       
#       GENE RPKM REPORT:
#       01. FB_ReleaseID
#       01. Gene_FBgnID
#       02. Gene_Symbol
#       03. Parent_library_FBlc#
#       04. Parent_library_name
#       05. RNASource_FBlc# (Library_FBlc#)
#       06. RNASource_name (Library_name)
#       07. RPKM_value
#       08. Bin_value
#       09. Unique_exon_base_count
#       10. Total_exon_base_count
#       11. Count_used
#       
#-----------------------------------------------------------------------------#
#
#	
#	
#-----------------------------------------------------------------------------#
use DBI;

if (@ARGV != 9) {
    print "\n USAGE: gene_rpkm_report pg_server db_name pg_username pg_password organism_abbrivation release_designator input_unstranded_transcriptome_gff input_stranded_transcriptome_gff output_filename\n\n";
    exit();
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);
my $org = shift(@ARGV);
my $rel = shift(@ARGV);

my $unstranded_gffin = shift(@ARGV);
my $stranded_gffin = shift(@ARGV);


$OUTFILE = shift(@ARGV);
open(STDERR,">>$OUTFILE");
open(STDOUT,">>$OUTFILE");


my $relid = $org . "_" . $rel;

#
##  DB Connections
#
## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh5 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh6 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";



#
##  General setup
#


## Setup output file header
$jetzt = scalar localtime;
print "## gene_rpkm_report report:\t$jetzt\n";
print "## using datasource: $dsource\n";
print "## using stranded transcriptome file: $stranded_gffin\n";
print "## using unstranded transcriptome file: $unstranded_gffin\n";
print "##\n";


#
## Main methods
#

## Get list of genes (polycistronic, etc) where we dont exclude shared regions from rpkm analysis
my $ciswc = '%cistronic%';
my $fbgnwc = 'FBgn%';
my %excluded_genes;
my $exgq = $dbh->prepare(sprintf("SELECT f.feature_id, f.uniquename, f.name as gname, cvt.name as gene_type, accession FROM feature f, feature_cvterm fc, cvterm cvt, dbxref dx, db, organism o, cvterm cvt2 WHERE f.is_obsolete = false and f.type_id = cvt2.cvterm_id and cvt2.name = 'gene' and f.uniquename like '%s' and f.organism_id = o.organism_id and o.abbreviation = '%s' and f.feature_id = fc.feature_id and fc.cvterm_id = cvt.cvterm_id and cvt.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and db.name = 'SO' and cvt.name like '%s'",$fbgnwc,$org,$ciswc));
$exgq->execute or die "WARNING: ERROR: Unable to execute excluded genes query\n";
while (my %exgr = %{$exgq->fetchrow_hashref}) {
  $excluded_genes{$exgr{uniquename}}{name} = $exgr{gname};
  $excluded_genes{$exgr{uniquename}}{type} = $exgr{gene_type};
#  print "Setting excluded gene:\t$exgr{uniquename}\t$exgr{gname}\t$exgr{gene_type}\n";
}

## Read stranded & unstranded transcriptome gff files into memory
my %stranded_txome;
my %unstranded_txome;
open(SIN,$stranded_gffin) or die "WARNING: Unable to open stranded tx gff input: .$stranded_gffin.\n";
open(UNSIN,$unstranded_gffin) or die "WARNING: Unable to open unstranded tx gff input: .$unstranded_gffin.\n";
#print "\nParsing stranded transcriptome gff file: $stranded_gffin...\n";
while (<SIN>) {
  next if ($_ =~ /^\#\#/);
  chomp;
#  print "\nProcessing: .$_.\n";
  my ($arm,$blah1,$blah2,$fmin,$fmax,$blah3,$strand,$blah4,$attrib) = split(/\t/,$_);
  my ($id,$rent_gene,$rent_tx,$score) = split(/\;/,$attrib);
  $score =~ s/score_text\=//;
  my $slen = $fmax - $fmin;
#  print "\t$fmin..$fmax\t($slen)\t$rent_gene\t$rent_tx\t$score\n";
  if ($rent_gene =~ /Parent\_gene\=(.*)/) {
    my @rents = split(/\,/,$1);
    foreach my $rent (@rents) {
#      print "\t\trent: $rent\n";
      $stranded_txome{$rent}{$score} += $slen;
    }
  }
  else {
    print "\tWARNING: ERROR: Unable to parser parent_gene attribute: .$rent_gene.\n";
    exit();
  }
}
#print "\nParsing unstranded transcriptome gff file: $unstranded_gffin...\n";
while (<UNSIN>) {
  next if ($_ =~ /^\#\#/);
  chomp;
#  print "\nProcessing: .$_.\n";
  my ($arm,$blah1,$blah2,$fmin,$fmax,$blah3,$strand,$blah4,$attrib) = split(/\t/,$_);
  my ($id,$rent_gene,$rent_tx,$score) = split(/\;/,$attrib);
  $score =~ s/score_text\=//;
  my $slen = $fmax - $fmin;
#  print "\t$fmin..$fmax\t($slen)\t$rent_gene\t$rent_tx\t$score\n";
  if ($rent_gene =~ /Parent\_gene\=(.*)/) {
    my @rents = split(/\,/,$1);
    foreach my $rent (@rents) {
#      print "\t\trent: $rent\n";
      $unstranded_txome{$rent}{$score} += $slen;
    }
  }
  else {
    print "\tWARNING: ERROR: Unable to parser parent_gene attribute: .$rent_gene.\n";
    exit();
  }
}

## For debugging...
# print "\n\n";
# foreach my $gene (sort keys(%unstranded_txome)) {
#   next if ($gene =~ /^\s*$/);
#   print "REPORT\t", $gene, "\t", $unstranded_txome{$gene}{unique} || "0",  "\t", $unstranded_txome{$gene}{shared} || "0", "\n";
# }
# print "\n\n";

print "##Release_ID\tFBgn#\tGeneSymbol\tParent_library_FBlc#\tParent_library_name\tRNASource_FBlc#\tRNASource_name\tRPKM_value\tBin_value\tUnique_exon_base_count\tTotal_exon_base_count\tCount_used\n";

## Main driver gets gene & rpkm values from chado
## Driver
my $fbgnwc = 'FBgn%';
my $icwc = 'Local wiggle file location%';
## Main driver (Modified for new metadata implementation - 20171002)
my $rq = $dbh2->prepare(sprintf("SELECT f.uniquename as funame, f.name as fname, p.uniquename as puname, p.name as pname, l.library_id as libid, l.uniquename as luname, l.name as lname, lfp.value from organism o, library l, feature f, library_feature lf, cvterm cvt, library_featureprop lfp, library p, library_relationship lr , cvterm cvt2, cvterm cvt3 WHERE f.is_obsolete = false and f.uniquename LIKE 'FBgn%' and f.organism_id = o.organism_id and o.abbreviation = '%s' and f.feature_id = lf.feature_id and lf.library_id = l.library_id and lf.library_feature_id = lfp.library_feature_id and lfp.type_id = cvt2.cvterm_id and cvt2.name in ('RPKM') and l.type_id = cvt.cvterm_id and cvt.name = 'result' and l.library_id = lr.subject_id and lr.object_id = p.library_id and p.type_id = cvt3.cvterm_id and cvt3.name = 'project' order by funame, puname, luname",$org));
$rq->execute or die "WARNING: ERROR: Unable to execute rpkm query\n";
while (my %rr = %{$rq->fetchrow_hashref}) {
## Set count_status -- this is set depending on whether a gene is one of those excluded from RPKM calculation only on unique bases
  my $count_status = 'Unique';   ## If the gene is one of our excluded genes, we uses 'total' bases in rpkm calculation; if not, we used 'unique'
  if ($excluded_genes{$rr{funame}}) {
    $count_status = 'Total';
  }
## Set bin value
  my $bin = &rpkmbin($rr{value});
## Get internalnotes that indicate whether this library is stranded or non-stranded and set $strandedness
  my $strandedness;
  my $iq = $dbh3->prepare(sprintf("SELECT * from libraryprop lp, cvterm cvt where lp.library_id = %d and lp.type_id = cvt.cvterm_id and cvt.name = 'internalnotes' and value like '%s'",$rr{libid},$icwc));
  $iq->execute or die "WARNING: ERROR: Unable to execute library internalnote (strandedness) query for library: $rr{libid}\t$rr{luname}\t$rr{lname}\n";
  my $iq_cnt = $iq->rows;
  if ($iq_cnt > 0) {
    while (my %ir = %{$iq->fetchrow_hashref}) {
      if ($ir{value} !~ /\(strand: 0/) {
	$strandedness = 'stranded'
      }
      else {
	$strandedness = 'unstranded'
      }
    }
    if ($strandedness eq 'stranded') {
	my $stotal = $stranded_txome{$rr{funame}}{unique} + $stranded_txome{$rr{funame}}{shared};
#	print "REPORT_STRANDED\t", $rr{funame}, "\t", $rr{fname}, "\t", $rr{puname}, "\t", $rr{pname}, "\t", $rr{luname}, "\t", $rr{lname}, "\t", $rr{value}, "\t", $bin, "\t", $stranded_txome{$rr{funame}}{unique}, "\t", $stotal, "\t", $count_status, "\n";
	print $relid, "\t", $rr{funame}, "\t", $rr{fname}, "\t", $rr{puname}, "\t", $rr{pname}, "\t", $rr{luname}, "\t", $rr{lname}, "\t", $rr{value}, "\t", $bin, "\t", $stranded_txome{$rr{funame}}{unique}, "\t", $stotal, "\t", $count_status, "\n";
    }
    else {
	my $unstotal = $unstranded_txome{$rr{funame}}{unique} + $unstranded_txome{$rr{funame}}{shared};
#	print "REPORT_UNSTRANDED\t", $rr{funame}, "\t", $rr{fname}, "\t", $rr{puname}, "\t", $rr{pname}, "\t", $rr{luname}, "\t", $rr{lname}, "\t", $rr{value}, "\t", $bin, "\t", $unstranded_txome{$rr{funame}}{unique}, "\t", $unstotal, "\t", $count_status, "\n";
	print $relid, "\t", $rr{funame}, "\t", $rr{fname}, "\t", $rr{puname}, "\t", $rr{pname}, "\t", $rr{luname}, "\t", $rr{lname}, "\t", $rr{value}, "\t", $bin, "\t", $unstranded_txome{$rr{funame}}{unique}, "\t", $unstotal, "\t", $count_status, "\n";
    }
  }
  else {
    print "\tWARNING: ERROR: No strandedness internalnote info found for library: $rr{libid}\t$rr{luname}\t$rr{lname}\n";
  }
}

$jetzt = scalar localtime;
print "\n\n## Finished gene_rpkm_report2: $jetzt\n";
exit();


sub rpkmbin {
## Get RPKM expression level bin# for a given RPKM

  my $inrpkm = shift;

  my $outbin;

  if ($inrpkm < 1) {
    $outbin = 0;
  }
  elsif ($inrpkm < 4) {
    $outbin = 1;
  }
  elsif ($inrpkm < 11) {
    $outbin = 2;
  }
  elsif ($inrpkm < 26) {
    $outbin = 3;
  }
  elsif ($inrpkm < 51) {
    $outbin = 4;
  }
  elsif ($inrpkm < 101) {
    $outbin = 5;
  }
  elsif ($inrpkm < 1001) {
    $outbin = 6;
  }
  else {
    $outbin = 7;
  }

  return($outbin);
}
