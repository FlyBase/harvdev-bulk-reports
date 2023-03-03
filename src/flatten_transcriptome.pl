#!/usr/bin/perl
# flatten_transcriptome
#
#       Perl script produces 'flattened' transcriptome gff3 for all arms
#       of given organism genome
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
    print "\n USAGE: flatten_transcriptome pg_server db_name pg_username pg_password organism_abbreviation log_output_filename gff_output_filename\n\n";
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

my %strandhash;
$strandhash{1} = '+';
$strandhash{-1} = '-';

## Setup output file header
$jetzt = scalar localtime;
print "## Starting flatten_transcriptome report\n## Using database: $dsource\n## Starting: $jetzt\n\n";

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

## This is the hash of completely overlapped additional genes from coverage2gene... hang onto it for now (for testing)
# my %noverlap = ('FBgn0035889','mkg-p','FBgn0261794','kcc','FBgn0261790','SmE','FBgn0261679','CG42726','FBgn0261680','CG42727','FBgn0261683','CG42730','FBgn0261697','tectonic','FBgn0261862','whd','FBgn0261710','nocte','FBgn0082983','snoRNA:Psi28S-2626','FBgn0261714','Cpn','FBgn0261791','SmG','FBgn0261933','SmD1','FBgn0261799','dsx-c73A','FBgn0083952','CG34116','FBgn0051210','CG31210','FBgn0261823','Asx','FBgn0261844','CG42776','FBgn0261786','mi','FBgn0261850','Xpd','FBgn0261709','CR42746','FBgn0261787','bru','FBgn0040995','CG12508','FBgn0033614','Obp47b','FBgn0038522','CG18139','FBgn0261922','CG10260','FBgn0261800','LanB1','FBgn0261686','CG42733','FBgn0040623','Spase12','FBgn0053251','CG33251','FBgn0261881','l(2)35Be','FBgn0033848','CG13330','FBgn0083942','CG34106','FBgn0032089','CG9573','FBgn0261793','Trf2','FBgn0261885','osa','FBgn0086068','snoRNA:Or-CD12','FBgn0261699','CG42735','FBgn0065059','snoRNA:660','FBgn0086056','snoRNA:Me18S-A425','FBgn0038482','CG4053','FBgn0030733','CG3560','FBgn0261789','SmD2','FBgn0261808','cu','FBgn0261925','CG42792','FBgn0261704','CG42740','FBgn0261858','CG42787','FBgn0053052','CG33052','FBgn0065084','snmRNA:763','FBgn0051086','CG31086','FBgn0065080','snoRNA:Me18S-A1374','FBgn0261882','l(2)35Bc','FBgn0038157','CG12538','FBgn0261860','CG42789','FBgn0053725','CG33725','FBgn0261845','CG42777','FBgn0261792','snRNP-U1-C','FBgn0261818','CG42760','FBgn0040291','Roc1b','FBgn0031715','tomb','FBgn0261708','CR42745','FBgn0261859','CG42788','FBgn0034169','CG9013','FBgn0053199','CG33199','FBgn0261688','Gyc76C','FBgn0261698','CG42732','FBgn0261703','gce','FBgn0261705','CG42741','FBgn0261722','flower','FBgn0261723','Dbx','FBgn0261788','Ank2','FBgn0261797','Dhc64C','FBgn0261801','CG42747','FBgn0261802','CG42748','FBgn0261803','CG42749','FBgn0261804','CG42750','FBgn0261811','pico','FBgn0261813','CG42755','FBgn0261822','Bsg','FBgn0261836','Msp-300','FBgn0261837','CG42769','FBgn0261838','CG42770','FBgn0261839','CG42771','FBgn0261840','CG42772','FBgn0261841','CG42773','FBgn0261842','CG42774','FBgn0261843','CG42775','FBgn0261854','aPKC','FBgn0261855','CG42784','FBgn0261871','dpr2','FBgn0261872','scaf6','FBgn0261873','sdt','FBgn0261914','Glut1','FBgn0261928','CG42795','FBgn0261929','CG42796','FBgn0261930','vnd','FBgn0261931','CG42797','FBgn0263757','Rpb4','FBgn0085292','CG34263','FBgn0002863','Acp95EF','FBgn0040687','CG14645','FBgn0032595','CG17996');


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
	push(@{$edges{$er{strand}}{$er{fmin}}{fmin}{transcripts}},$gr{tuname});
	push(@{$edges{$er{strand}}{$er{fmax}}{fmax}{transcripts}},$gr{tuname});


## Only add (fmin) FBgn which are not already in the array...
	my %fming;
	map { $fming{$_} = 1 } @{$edges{$er{strand}}{$er{fmin}}{fmin}{genes}};
	push(@{$edges{$er{strand}}{$er{fmin}}{fmin}{genes}}, grep { !exists $fming{$_} } $gr{guname});

## Only add (fmax) FBgn which are not already in the array...
	my %fmaxg;
	map { $fmaxg{$_} = 1 } @{$edges{$er{strand}}{$er{fmax}}{fmax}{genes}};
	push(@{$edges{$er{strand}}{$er{fmax}}{fmax}{genes}}, grep { !exists $fmaxg{$_} } $gr{guname});
      }
    }
    else {
## If the tx doesnt have exons (e.g., miRNAs), use tx location 
      print "\t\tWARNING: No exon found for tx; using tx location... $gr{feature_id}\t$gr{tuname}\t$gr{tname}\n";
  }
  }

#
## Now process edges
## We step through the exon edges (on each strand), building segments as we go, one %active_seg at a time...
#
  foreach my $strand (keys(%edges)) {
    my %active_seg;
    print "\tNow processing edges on strand $strand...\n";
    foreach my $edge (sort { $a <=> $b } keys(%{$edges{$strand}})) {

#
## Both FMAX & FMIN edges
#
      if ( ($edges{$strand}{$edge}{fmin}) && ($edges{$strand}{$edge}{fmax}) ) {
	print(sprintf("Edge (BOTH fmin & fmax): %d\t(%s)\n",$edge,join(", ",@{$edges{$strand}{$edge}{fmin}{genes}})));

	if (%active_seg) {
	  
	  print "\t\tEXISTING SEG (new fmin)...\n";
	  print "\t\t\tWRITE SEG NOW...\n";
	  print "\t\t\t\twrite loc: $active_seg{fmin}..$edge\n";
	  my $seqlen = $edge - $active_seg{fmin};
	  print(sprintf("\t\t\t\tseqlen is: %d\n",$seqlen));
	  print(sprintf("\t\t\t\twrite tx: %s\n",join(", ",sort(@{$active_seg{transcripts}}))));
	  print(sprintf("\t\t\t\twrite gn: %s\n",join(", ",sort(@{$active_seg{genes}}))));
## Print GFF
	my $gn_cnt = scalar @{$active_seg{genes}};
	my $etype = 'unique';
	if ($gn_cnt > 1) {
	  $etype = 'shared';
	}
#	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",join("_",sort(@{$active_seg{transcripts}})),join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	  my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",$segid,join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	  print GFFOUT sprintf("%s\tFlyBase\texon_region\t%d\t%d\t.\t%s\t.\t%s\n",$armh{$arm}{uname},$active_seg{fmin},$edge,$strandhash{$strand},$ID);
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
	  push(@{$active_seg{genes}}, grep { !exists $g{$_} } @{$edges{$strand}{$edge}{fmin}{genes}});
## Push any transcripts who are not already in active_seg
	  my %t;
	  map { $t{$_} = 1 } @{$active_seg{transcripts}};
	  push(@{$active_seg{transcripts}}, grep { !exists $t{$_} } @{$edges{$strand}{$edge}{fmin}{transcripts}});

	  print "\t\t\tactive_seg fmin:\t$active_seg{fmin}\n";
	  print(sprintf("\t\t\tactive_seg tx: %s\n",join(", ",@{$active_seg{transcripts}})));
	  print(sprintf("\t\t\tactive_seg gn: %s\n",join(", ",@{$active_seg{genes}})));

## Remove FBtr that are not continued	 in next active_seg
	  my @remaining_tx = grep { not $_ ~~ @{$edges{$strand}{$edge}{fmax}{transcripts}} } @{$active_seg{transcripts}};
	  print(sprintf("\t\t\tremaining_tx is: %s\n",join(", ",@remaining_tx)));
	  @{$active_seg{transcripts}} = @remaining_tx;

## Remove FBgn that are not continued in next active_seg
#	  my @remaining_gn = grep { not $_ ~~ @{$edges{$strand}{$edge}{fmax}{genes}} } @{$active_seg{genes}};
#	  print(sprintf("\t\t\tremaining_gn is: %s\n",join(", ",@remaining_gn)));
#	  my $rgn_cnt = scalar @remaining_gn;
#	  @{$active_seg{genes}} = @remaining_gn if ($rgn_cnt > 0);

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

## Undef active_seg if all transcripts have ended
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
      elsif ($edges{$strand}{$edge}{fmin}) {

	print(sprintf("Edge (fmin): %d\t(%s)\n",$edge,join(", ",@{$edges{$strand}{$edge}{fmin}{genes}})));
	
## For debugging...
	print(sprintf("  genes are: %s\n",join(", ",@{$edges{$strand}{$edge}{fmin}{genes}})));
	print(sprintf("  transcripts are: %s\n",join(", ",@{$edges{$strand}{$edge}{fmin}{transcripts}})));

	if (%active_seg) {
	  print "\t\tEXISTING SEG (new fmin)...\n";
	  print "\t\t\tWRITE SEG NOW...\n";
	  print "\t\t\t\twrite loc: $active_seg{fmin}..$edge\n";
	  my $seqlen = $edge - $active_seg{fmin};
	  print(sprintf("\t\t\t\tseqlen is: %d\n",$seqlen));
	  print(sprintf("\t\t\t\twrite tx: %s\n",join(", ",sort(@{$active_seg{transcripts}}))));
	  print(sprintf("\t\t\t\twrite gn: %s\n",join(", ",sort(@{$active_seg{genes}}))));
## Print GFF
	my $gn_cnt = scalar @{$active_seg{genes}};
	my $etype = 'unique';
	if ($gn_cnt > 1) {
	  $etype = 'shared';
	}
#	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",join("_",sort(@{$active_seg{transcripts}})),join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",$segid,join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	print GFFOUT sprintf("%s\tFlyBase\texon_region\t%d\t%d\t.\t%s\t.\t%s\n",$armh{$arm}{uname},$active_seg{fmin},$edge,$strandhash{$strand},$ID);
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
	  push(@{$active_seg{genes}}, grep { !exists $g{$_} } @{$edges{$strand}{$edge}{fmin}{genes}});
## Push any transcripts who are not already in active_seg
	  my %t;
	  map { $t{$_} = 1 } @{$active_seg{transcripts}};
	  push(@{$active_seg{transcripts}}, grep { !exists $t{$_} } @{$edges{$strand}{$edge}{fmin}{transcripts}});

	  print "\t\t\tactive_seg fmin:\t$active_seg{fmin}\n";
	  print(sprintf("\t\t\tactive_seg tx: %s\n",join(", ",@{$active_seg{transcripts}})));
	  print(sprintf("\t\t\tactive_seg gn: %s\n",join(", ",@{$active_seg{genes}})));

	}
	else {
	  print "\t\tSTARTING NEW SEG...\n";

## Set active_seg
	  $active_seg{fmin} = $edge;
	  $active_seg{genes} = $edges{$strand}{$edge}{fmin}{genes};
	  $active_seg{transcripts} = $edges{$strand}{$edge}{fmin}{transcripts};

	  print "\t\t\tactive_seg fmin:\t$active_seg{fmin}\n";
	  print(sprintf("\t\t\tactive_seg tx: %s\n",join(", ",@{$active_seg{transcripts}})));
	  print(sprintf("\t\t\tactive_seg gn: %s\n",join(", ",@{$active_seg{genes}})));
	}
      }


#
## FMAX edge
#
      elsif ($edges{$strand}{$edge}{fmax}) {
	print(sprintf(" Edge (fmax): %d\t(%s)\n",$edge,join(", ",@{$edges{$strand}{$edge}{fmax}{genes}})));

## For debugging...
	print(sprintf("  genes are: %s\n",join(", ",@{$edges{$strand}{$edge}{fmax}{genes}})));
	print(sprintf("  transcripts are: %s\n",join(", ",@{$edges{$strand}{$edge}{fmax}{transcripts}})));

	if (%active_seg) {
	  print "\t\tEXISTING SEG (new fmax)...\n";
	
	  $active_seg{fmax} = $edge;
	  print "\t\t\tWRITE SEG NOW...\n";
	  print "\t\t\t\twrite loc: $active_seg{fmin}..$active_seg{fmax}\n";
	  my $seqlen = $active_seg{fmax} - $active_seg{fmin};
	  print(sprintf("\t\t\t\tseqlen is: %d\n",$seqlen));
	  print(sprintf("\t\t\t\twrite tx: %s\n",join(", ",sort(@{$active_seg{transcripts}}))));
	  print(sprintf("\t\t\t\twrite gn: %s\n",join(", ",sort(@{$active_seg{genes}}))));
## Print GFF
	my $gn_cnt = scalar @{$active_seg{genes}};
	my $etype = 'unique';
	if ($gn_cnt > 1) {
	  $etype = 'shared';
	}
#	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",join("_",sort(@{$active_seg{transcripts}})),join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	my $ID = sprintf("ID=%s;Parent_gene=%s;Parent=%s;score_text=%s;",$segid,join(",",sort(@{$active_seg{genes}})),join(",",sort(@{$active_seg{transcripts}})),$etype);
	print GFFOUT sprintf("%s\tFlyBase\texon_region\t%d\t%d\t.\t%s\t.\t%s\n",$armh{$arm}{uname},$active_seg{fmin},$active_seg{fmax},$strandhash{$strand},$ID);
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
	  my @remaining_tx = grep { not $_ ~~ @{$edges{$strand}{$edge}{fmax}{transcripts}} } @{$active_seg{transcripts}};
	  print(sprintf("\t\t\tremaining_tx is: %s\n",join(", ",@remaining_tx)));
	  @{$active_seg{transcripts}} = @remaining_tx;

## Remove FBgn that are not continued in next active_seg
# 	  my @remaining_gn = grep { not $_ ~~ @{$edges{$strand}{$edge}{fmax}{genes}} } @{$active_seg{genes}};
# 	  print(sprintf("\t\t\tremaining_gn is: %s\n",join(", ",@remaining_gn)));
# 	  my $rgn_cnt = scalar @remaining_gn;
# 	  @{$active_seg{genes}} = @remaining_gn if ($rgn_cnt > 0);

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



## Undef active_seg if all transcripts have ended
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
}

$jetzt = scalar localtime;
print "\n\n## Finished flatten_transcriptome: $jetzt\n";
exit();
