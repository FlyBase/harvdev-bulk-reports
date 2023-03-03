#!/usr/bin/perl
# flatten_transcriptome_unstranded
#
#       Perl script produces UNSTRANDED 'flattened' transcriptome
#       gff3 for all arms of given organism genome
#       
#       The gff output file is used in generating RPKM reports
#	
#-----------------------------------------------------------------------------#
#
#	NOTES
#	The algorithm used in generating the flattened transcriptome is the 
#	same as used in coverage2genes scripts
#	
#
#-----------------------------------------------------------------------------#
use DBI;
use Bio::Seq;
use Bio::SeqIO;

if (@ARGV != 7) {
    print "\n USAGE: flatten_transcriptome_unstranded pg_server db_name pg_username pg_password organism_abbreviation log_output_filename gff_output_filename\n\n";
    exit();
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);
my $org = shift(@ARGV);

$OUTFILE = shift(@ARGV);
open(STDERR,">>$OUTFILE");
open(STDOUT,">>$OUTFILE");

$GFFOUT = shift(@ARGV);
open(GFFOUT,">>$GFFOUT");

#
##  DB Connections
#
## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
# $dbh5 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
# $dbh6 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";


#
##  General setup
#

## Setup output file header
$jetzt = scalar localtime;
print "\n\n## flatten_transcriptome_unstranded report\n## Using database: $dsource\n## Starting: $jetzt\n\n";

## Setup gff output header
print GFFOUT "## gff-version 3\n";

## Build a hash of feature_id for chromosome_arm in chado
# print "\nGetting arm/scaffold for organism $org...\n";
my %armh;
my %armh2;
my %armres;
my $armscaff;
$armscaff = sprintf("'golden_path','ultra_scaffold'");

# Simple feature uniquename patterns.
my $fbgnwc = 'FBgn%';
my $fbtrwc = 'FBtr%';

## Get arm/scaffold feature_ids, uniquenames & residues (only those having genes).
my $caq = $dbh->prepare(sprintf("
    SELECT DISTINCT src.feature_id,
                    src.uniquename
    FROM feature f
    JOIN cvterm cvtf ON (cvtf.cvterm_id = f.type_id and cvtf.name = 'gene')
    JOIN featureloc fl ON fl.feature_id = f.feature_id
    JOIN feature src ON src.feature_id = fl.srcfeature_id
    JOIN cvterm cvts ON cvts.cvterm_id = src.type_id
    JOIN organism o ON o.organism_id = src.organism_id
    WHERE f.is_obsolete is false
      AND f.is_analysis is false
      AND f.uniquename LIKE '%s'
      AND src.is_obsolete is false
      AND src.is_analysis is false
      AND cvts.name in (%s)
      AND o.abbreviation = '%s'
	",$fbgnwc, $armscaff, $org));
$caq->execute or die "WARNING: ERROR: Cannot get feature_id / chromosome_arm name hash\n";
while (%cq = %{$caq->fetchrow_hashref}) {
#    print "\tPushing arm/scaffold: $cq{feature_id}\t$cq{uniquename}...\n";
    $armh{$cq{feature_id}}{uname} = $cq{uniquename};
    $armh2{$cq{uniquename}} = $cq{feature_id};
#    $armres{$cq{feature_id}} = Bio::Seq->new (-seq => $cq{residues}, -alphabet => 'dna');
}

 ## Get genus and species for input organism abbreviation
 my $genus;
 my $species;
 my $orq = $dbh->prepare(sprintf("SELECT * from organism where abbreviation = '%s'",$org));
 $orq->execute or die "WARNING: ERROR: Unable to execute organism query\n";
 my $orq_cnt = $orq->rows;
 if ($orq_cnt == 1) {
   while (my %orr = %{$orq->fetchrow_hashref}) {
     $genus = $orr{genus};
     $species = $orr{species};
   }
 }
 else {
   print "WARNING: ERROR: organism record not found for input abbreviation: .$org.\n";
   exit();
 }


#
## Main methods
#

## For each arm we get all exon edges in %edge object
my $segid = 1;
foreach my $arm (keys(%armh)) {
  my %gnameh;
  my %txgeneh;
  my %edges;   ## Hash of edge objects
  my %ghash;  ## Hash of gene info (for report)
#  next if ($armh{$arm}{uname} ne '3R');
  print "\nFetching exon edges on arm: $armh{$arm}{uname}...\n";
## Query below copied from flatten_transcriptome3
  my $gq = $dbh->prepare(sprintf("
    SELECT g.uniquename AS guname,
           g.name AS gname,
           gl.fmin AS gfmin,
    	   gl.fmax AS gfmax,
           t.feature_id,
           t.uniquename AS tuname,
           t.name AS tname,
           cvt.name AS ttype,
           fl.fmin,
           fl.fmax,
           fl.strand
    FROM feature g
    JOIN featureloc gl ON gl.feature_id = g.feature_id
    JOIN feature_relationship fr ON fr.object_id = g.feature_id
    JOIN feature t ON t.feature_id = fr.subject_id
    JOIN featureloc fl ON fl.feature_id = t.feature_id
    JOIN organism o ON o.organism_id = g.organism_id
    JOIN cvterm cvt ON (cvt.cvterm_id = t.type_id AND cvt.name in ('mRNA','tRNA','ncRNA','pre_miRNA','rRNA','snRNA','snoRNA','pseudogene'))
    JOIN cvterm cvt2 ON (cvt2.cvterm_id = g.type_id AND cvt2.name = 'gene')
    JOIN cvterm cvt3 ON (cvt3.cvterm_id = fr.type_id AND cvt3.name = 'partof')
    WHERE g.is_obsolete is false
      AND g.is_analysis is false
      AND t.is_obsolete is false
      AND t.is_analysis is false
      AND o.abbreviation = '%s'
      AND fl.srcfeature_id = %d
      AND g.uniquename LIKE '%s'
      AND t.uniquename LIKE '%s'
    ", $org, $arm, $fbgnwc, $fbtrwc));
  $gq->execute or die "WARNING: ERROR: Unable to execute gn/tr query for arm: .$armh{$arm}{uname}\n";
  while (my %gr = %{$gq->fetchrow_hashref}) {
    $gnameh{$gr{guname}} = $gr{gname};
    $txgeneh{$gr{tuname}} = $gr{guname};
    print "\tProcessing tx: $gr{guname}\t$gr{gname}\t$gr{feature_id}\t$gr{tuname}\t$gr{tname}\n";
    my $eq = $dbh3->prepare(sprintf("SELECT e.uniquename, fmin, fmax, strand from feature e, featureloc fl, feature_relationship fr, cvterm cvt1, cvterm cvt2 where e.feature_id = subject_id and e.type_id = cvt1.cvterm_id and cvt1.name in ('exon','miRNA') and fr.type_id = cvt2.cvterm_id and cvt2.name in ('partof','producedby') and e.feature_id = fl.feature_id and object_id = %d",$gr{feature_id}));
    $eq->execute or die "WARNING: ERROR: Unable to execute exon query\n";
    my $eq_cnt = $eq->rows;
    if ($eq_cnt > 0) {
      while (my %er = %{$eq->fetchrow_hashref}) {
#	print "\t\tExon found: $er{uniquename}\t$er{fmin}\t$er{fmax}\n";
## For each exon from a tx, add info to edges hash...

## Add FBtr	to their array (there wont be duplicates... will there?)  
	push(@{$edges{$er{fmin}}{fmin}{transcripts}},$gr{tuname});
	push(@{$edges{$er{fmax}}{fmax}{transcripts}},$gr{tuname});


## Only add (fmin) FBgn which are not already in the array...
	my %fming;
	map { $fming{$_} = 1 } @{$edges{$er{fmin}}{fmin}{genes}};
	push(@{$edges{$er{fmin}}{fmin}{genes}}, grep { !exists $fming{$_} } $gr{guname});

## Only add (fmax) FBgn which are not already in the array...
	my %fmaxg;
	map { $fmaxg{$_} = 1 } @{$edges{$er{fmax}}{fmax}{genes}};
	push(@{$edges{$er{fmax}}{fmax}{genes}}, grep { !exists $fmaxg{$_} } $gr{guname});
      }
    }
    else {
## If the tx doesnt have exons (e.g., miRNAs), use tx location 
      print "\t\tWARNING: No exon found for tx; using tx location... $gr{feature_id}\t$gr{tuname}\t$gr{tname}\n";
      exit();
    }
  }
#
## Now process edges
## We step through the exon edges (on each strand), building segments as we go, one %active_seg at a time...
#
  my %active_seg;
  print "\tNow processing edges on strand $strand...\n";
  foreach my $edge (sort { $a <=> $b } keys(%edges)) {

#
## Both FMAX & FMIN edges
#
    if ( ($edges{$edge}{fmin}) && ($edges{$edge}{fmax}) ) {
      print(sprintf("Edge (BOTH fmin & fmax): %d\t(%s)\n",$edge,join(", ",@{$edges{$edge}{fmin}{genes}})));

      if (%active_seg) {
	  
	print "\t\tEXISTING SEG (new fmin)...\n";
	print "\t\t\tWRITE SEG NOW...\n";
	print "\t\t\t\twrite loc: $active_seg{fmin}..$edge\n";
	my $seqlen = $edge - $active_seg{fmin};
	print(sprintf("\t\t\t\tseqlen is: %d\n",$seqlen));
	print(sprintf("\t\t\t\twrite tx: %s\n",join(", ",sort(@{$active_seg{transcripts}}))));
	print(sprintf("\t\t\t\twrite gn: %s\n",join(", ",sort(@{$active_seg{genes}}))));
## Print GFF -- We use "0" for strand (since this is unstranded)
	my $gn_cnt = scalar @{$active_seg{genes}};
	my $etype = 'unique';
	if ($gn_cnt > 1) {
	  $etype = 'shared';
	}
#	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",join("_",sort(@{$active_seg{transcripts}})),join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",$segid,join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	print GFFOUT sprintf("%s\tFlyBase\texon_region\t%d\t%d\t.\t0\t.\t%s\n",$armh{$arm}{uname},$active_seg{fmin},$edge,$ID);
	  $segid++;

## Pop gene info hash (for report)
	my $gn_cnt = scalar @{$active_seg{genes}};
	if ($gn_cnt > 1) {
	  foreach my $gene (@{$active_seg{genes}}) {
	    $ghash{$gene}{shared} += $seqlen;
	  }
	}
	else {
	  foreach my $gene (@{$active_seg{genes}}) {
	    $ghash{$gene}{unique} += $seqlen;
	  }
	}

## Reset active_seg
	$active_seg{fmin} = $edge;

## Push any genes who are not already in active_seg
	my %g;
	map { $g{$_} = 1 } @{$active_seg{genes}};
	push(@{$active_seg{genes}}, grep { !exists $g{$_} } @{$edges{$edge}{fmin}{genes}});
## Push any transcripts who are not already in active_seg
	my %t;
	map { $t{$_} = 1 } @{$active_seg{transcripts}};
	push(@{$active_seg{transcripts}}, grep { !exists $t{$_} } @{$edges{$edge}{fmin}{transcripts}});

	print "\t\t\tactive_seg fmin:\t$active_seg{fmin}\n";
	print(sprintf("\t\t\tactive_seg tx: %s\n",join(", ",@{$active_seg{transcripts}})));
	print(sprintf("\t\t\tactive_seg gn: %s\n",join(", ",@{$active_seg{genes}})));

## Remove FBtr that are not continued	 in next active_seg
	my @remaining_tx = grep { not $_ ~~ @{$edges{$edge}{fmax}{transcripts}} } @{$active_seg{transcripts}};
	print(sprintf("\t\t\tremaining_tx is: %s\n",join(", ",@remaining_tx)));
	@{$active_seg{transcripts}} = @remaining_tx;

## Remove FBgn that are not continued in next active_seg
#	my @remaining_gn = grep { not $_ ~~ @{$edges{$edge}{fmax}{genes}} } @{$active_seg{genes}};
#	print(sprintf("\t\t\tremaining_gn is: %s\n",join(", ",@remaining_gn)));
#	my $rgn_cnt = scalar @remaining_gn;
#	@{$active_seg{genes}} = @remaining_gn if ($rgn_cnt > 0);
## Set @remaining_gn to contain gene for each tx
	  my @remaining_gn;
	  foreach my $gtx (@remaining_tx) {
	    my %tg;
	    map { $tg{$_} = 1 } @remaining_gn; 
	    push(@remaining_gn, grep {!exists $tg{$_}} $txgeneh{$gtx});
	  }
 	  print(sprintf("\t\t\tremaining_gn is: %s\n",join(", ",@remaining_gn)));
 	  my $rgn_cnt = scalar @remaining_gn;
 	  @{$active_seg{genes}} = @remaining_gn if ($rgn_cnt > 0);

## Undef active_seg of all transcripts have ended
	my $tx_cnt = scalar @remaining_tx;
	if (1 > $tx_cnt) {
	  undef(%active_seg);
	}

      }
      else {
	print "\t\t\tWARNING: ERROR: fmax found (shared fmin/fmax edge) when no active_seg!\n";
	exit;
      }

      }

#
## FMIN edge
#
    elsif ($edges{$edge}{fmin}) {

      print(sprintf("Edge (fmin): %d\t(%s)\n",$edge,join(", ",@{$edges{$edge}{fmin}{genes}})));
	
## For debugging...
      print(sprintf("  genes are: %s\n",join(", ",@{$edges{$edge}{fmin}{genes}})));
      print(sprintf("  transcripts are: %s\n",join(", ",@{$edges{$edge}{fmin}{transcripts}})));

      if (%active_seg) {
	print "\t\tEXISTING SEG (new fmin)...\n";
	print "\t\t\tWRITE SEG NOW...\n";
	print "\t\t\t\twrite loc: $active_seg{fmin}..$edge\n";
	my $seqlen = $edge - $active_seg{fmin};
	print(sprintf("\t\t\t\tseqlen is: %d\n",$seqlen));
	print(sprintf("\t\t\t\twrite tx: %s\n",join(", ",sort(@{$active_seg{transcripts}}))));
	print(sprintf("\t\t\t\twrite gn: %s\n",join(", ",sort(@{$active_seg{genes}}))));
## Print GFF -- We use "0" for strand (since this is unstranded)
	my $gn_cnt = scalar @{$active_seg{genes}};
	my $etype = 'unique';
	if ($gn_cnt > 1) {
	  $etype = 'shared';
	}
#	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",join("_",sort(@{$active_seg{transcripts}})),join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",$segid,join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	print GFFOUT sprintf("%s\tFlyBase\texon_region\t%d\t%d\t.\t0\t.\t%s\n",$armh{$arm}{uname},$active_seg{fmin},$edge,$ID);
	  $segid++;


## Pop gene info hash (for report)
	my $gn_cnt = scalar @{$active_seg{genes}};
	if ($gn_cnt > 1) {
	  foreach my $gene (@{$active_seg{genes}}) {
	    $ghash{$gene}{shared} += $seqlen;
	  }
	}
	else {
	  foreach my $gene (@{$active_seg{genes}}) {
	    $ghash{$gene}{unique} += $seqlen;
	  }
	}

## Reset active_seg
	$active_seg{fmin} = $edge;

## Push any genes who are not already in active_seg
	my %g;
	map { $g{$_} = 1 } @{$active_seg{genes}};
	push(@{$active_seg{genes}}, grep { !exists $g{$_} } @{$edges{$edge}{fmin}{genes}});
## Push any transcripts who are not already in active_seg
	my %t;
	map { $t{$_} = 1 } @{$active_seg{transcripts}};
	push(@{$active_seg{transcripts}}, grep { !exists $t{$_} } @{$edges{$edge}{fmin}{transcripts}});

	print "\t\t\tactive_seg fmin:\t$active_seg{fmin}\n";
	print(sprintf("\t\t\tactive_seg tx: %s\n",join(", ",@{$active_seg{transcripts}})));
	print(sprintf("\t\t\tactive_seg gn: %s\n",join(", ",@{$active_seg{genes}})));

      }
      else {
	print "\t\tSTARTING NEW SEG...\n";

## Set active_seg
	$active_seg{fmin} = $edge;
	$active_seg{genes} = $edges{$edge}{fmin}{genes};
	$active_seg{transcripts} = $edges{$edge}{fmin}{transcripts};

	print "\t\t\tactive_seg fmin:\t$active_seg{fmin}\n";
	print(sprintf("\t\t\tactive_seg tx: %s\n",join(", ",@{$active_seg{transcripts}})));
	print(sprintf("\t\t\tactive_seg gn: %s\n",join(", ",@{$active_seg{genes}})));
      }
    }


#
## FMAX edge
#
    elsif ($edges{$edge}{fmax}) {
      print(sprintf(" Edge (fmax): %d\t(%s)\n",$edge,join(", ",@{$edges{$edge}{fmax}{genes}})));

## For debugging...
      print(sprintf("  genes are: %s\n",join(", ",@{$edges{$edge}{fmax}{genes}})));
      print(sprintf("  transcripts are: %s\n",join(", ",@{$edges{$edge}{fmax}{transcripts}})));

      if (%active_seg) {
	print "\t\tEXISTING SEG (new fmax)...\n";
	
	$active_seg{fmax} = $edge;
	print "\t\t\tWRITE SEG NOW...\n";
	print "\t\t\t\twrite loc: $active_seg{fmin}..$active_seg{fmax}\n";
	my $seqlen = $active_seg{fmax} - $active_seg{fmin};
	print(sprintf("\t\t\t\tseqlen is: %d\n",$seqlen));
	print(sprintf("\t\t\t\twrite tx: %s\n",join(", ",sort(@{$active_seg{transcripts}}))));
	print(sprintf("\t\t\t\twrite gn: %s\n",join(", ",sort(@{$active_seg{genes}}))));
## Print GFF -- We use "0" for strand (since this is unstranded)
	my $gn_cnt = scalar @{$active_seg{genes}};
	my $etype = 'unique';
	if ($gn_cnt > 1) {
	  $etype = 'shared';
	}
#	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",join("_",sort(@{$active_seg{transcripts}})),join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",$segid,join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	print GFFOUT sprintf("%s\tFlyBase\texon_region\t%d\t%d\t.\t0\t.\t%s\n",$armh{$arm}{uname},$active_seg{fmin},$active_seg{fmax},$ID);
	  $segid++;

## Pop gene info hash (for report)
	my $gn_cnt = scalar @{$active_seg{genes}};
	if ($gn_cnt > 1) {
	  foreach my $gene (@{$active_seg{genes}}) {
	    $ghash{$gene}{shared} += $seqlen;
	  }
	}
	else {
	  foreach my $gene (@{$active_seg{genes}}) {
	    $ghash{$gene}{unique} += $seqlen;
	  }
	}

## Evaluate and reset active_seg...

	$active_seg{fmin} = $edge;


## Remove FBtr that are not continued	 in next active_seg
	my @remaining_tx = grep { not $_ ~~ @{$edges{$edge}{fmax}{transcripts}} } @{$active_seg{transcripts}};
	print(sprintf("\t\t\tremaining_tx is: %s\n",join(", ",@remaining_tx)));
	@{$active_seg{transcripts}} = @remaining_tx;

## Remove FBgn that are not continued in next active_seg
#	my @remaining_gn = grep { not $_ ~~ @{$edges{$edge}{fmax}{genes}} } @{$active_seg{genes}};
#	print(sprintf("\t\t\tremaining_gn is: %s\n",join(", ",@remaining_gn)));
#	my $rgn_cnt = scalar @remaining_gn;
#	@{$active_seg{genes}} = @remaining_gn if ($rgn_cnt > 0);
## Set @remaining_gn to contain gene for each tx
	my @remaining_gn;
	foreach my $gtx (@remaining_tx) {
	  my %tg;
	  map { $tg{$_} = 1 } @remaining_gn; 
	  push(@remaining_gn, grep {!exists $tg{$_}} $txgeneh{$gtx});
	}
	print(sprintf("\t\t\tremaining_gn is: %s\n",join(", ",@remaining_gn)));
	my $rgn_cnt = scalar @remaining_gn;
	@{$active_seg{genes}} = @remaining_gn if ($rgn_cnt > 0);


## Undef active_seg of all transcripts have ended
	my $tx_cnt = scalar @remaining_tx;
	if (1 > $tx_cnt) {
	  undef(%active_seg);
	}

## For debugging...
	print(sprintf("\t active genes are: %s\n",join(@{$active_seg{genes}})));
	print(sprintf("\t active transcripts are: %s\n",join(@{$active_seg{transcripts}})));
      }
      else {
	print "\t\t\tWARNING: ERROR: fmax found when no active_seg!\n";
	exit;
      }
    }
  }
}
foreach my $gene (sort(keys(%ghash))) {
  print(sprintf("REPORT\t%s\t%s\t%s\t%d\t%d\n",$armh{$arm}{uname},$gene,$gnameh{$gene},$ghash{$gene}{unique},$ghash{$gene}{shared}));
}

$jetzt = scalar localtime;
print "\n\n## Finished flatten_transcriptome_unstranded: $jetzt\n";
exit();
