#!/usr/local/bin/perl
# report_gene_assessment_new
#
#       Perl script generates tab-delimited report_gene_assessment table with the 
#       following fields: gene_symbol, FBgn#, arm, fmin, fmax, strand, #_associated_alleles, 
#       #_assoc_alleles_w_phenotype_data, #_primary_fbrfs, #_reviews, 
#       #_papers_since_1999, #_papers_since_2003, #_papers_since_2006, 
#       #_papers_since_2006_harv_triaged_outstanding, #_papers_since_2006_harv_triaged_outstanding
#       #_papers_since_2009, #_papers_since_2009_harv_triaged_outstanding, 
#       #_papers_since_2009_harv_triaged_outstanding
#       #_genetic_interactions, #_phenotype_statements, 
#       #_tr_epression_tap_statements, #_pr_epression_tap_statements, #_interpro_acc#,
#       #_go_biological_process_terms, #_go_molecular_function_terms,
#       #_go_cellular_component_terms
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

if (@ARGV < 6) {
    print "\n USAGE: report_gene_assessment pg_server db_name pg_username pg_password organism_abbreviation output_filename\n\n";
    print "\toutput_filename is the output file for std and error output.\n\n";
    exit;
}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);
my $org = shift(@ARGV);

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
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh5 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh6 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh7 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh8 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh9 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


#
##  General setup
#

## Build a hash of feature_id for chromosome_arm in chado
# print "\nGetting arm/scaffold for organism $org...\n";
my %armh;
my %armh2;
my %armres;
my $armscaff;
$armscaff = sprintf("'golden_path_region','ultra_scaffold'");

## Get arm/scaffold feature_ids, uniquenames & residues
my $caq = $dbh->prepare(sprintf("SELECT feature_id, uniquename from feature f, cvterm cvt, organism o where cvt.name in (%s) and cvterm_id = f.type_id and f.is_obsolete = 'f' and f.is_analysis = 'f' and f.organism_id = o.organism_id and o.abbreviation = '%s' and uniquename != 'dmel_mitochondrion_genome'",$armscaff,$org));
$caq->execute or die "WARNING: ERROR: Cannot get feature_id / chromosome_arm name hash\n";
while (%cq = %{$caq->fetchrow_hashref}) {
#  print "\tPushing arm/scaffold: $cq{feature_id}\t$cq{uniquename}...\n";
  $armh{$cq{feature_id}}{uname} = $cq{uniquename};
  $armh2{$cq{uniquename}} = $cq{feature_id};
}

## Get uniquename & residues for dmel mitochondrion specially
if ($org eq 'Dmel') {
  my $caq = $dbh->prepare(sprintf("SELECT feature_id, uniquename from feature f, cvterm cvt, organism o where cvt.name = 'chromosome' and cvterm_id = f.type_id
and f.is_obsolete = 'f' and f.is_analysis = 'f' and f.organism_id = o.organism_id and o.abbreviation = '%s' and uniquename = 'dmel_mitochondrion_genome'",$org));
  $caq->execute or die "WARNING: ERROR: Cannot get feature_id / chromosome_arm name hash\n";
  while (my %cq = %{$caq->fetchrow_hashref}) {
#    print "\tPushing arm/scaffold: $cq{feature_id}\t$cq{uniquename}...\n";
    $armh{$cq{feature_id}}{uname} = $cq{uniquename};
    $armh2{$cq{uniquename}} = $cq{feature_id};
  }
}

## Setup file header
$jetzt = scalar localtime;
print "## FlyBase Gene assessment report ($org)\n## Generated: $jetzt\n";
print "## Using script: (fb_cvs)FB/scripts/reports/report_gene_assessment\n";
print "## Using datasource: $dsource...\n\n";
print "##primary_FBgn#\tgene_symbol\tarm\tfmin\tfmax\tstrand\t#_associated_alleles\t#_assoc_alleles_w_phenotype_data\t#_primary_FBrfs\t#_review_FBrfs\t#_FBrfs_since_1999\t#_FBrfs_since_2003\t#_FBrfs_since_2006\t#_since_2006_harv_triaged_outstanding\t$_since_2006_camb_triaged_outstanding\t#_FBrfs_since_2009\t#_since_2009_harv_triaged_outstanding\t#_since_2009_camb_triaged_outstanding\t#_genetic_interactions\t#_phenotype_statements\t#_tr_expression_tap_statements\t#_pr_expression_tap_statements\t#_reporter_constructs_w_expression_data\t#_enhancer_traps_w_expression_data\t#_interpro_acc\t#_go_biological_process_terms\t#_go_molecular_function_terms\t#_go_subcellular_component_terms\n";


#
## Main methods
#

#
## Get genes...
#
my $fbgnwc = 'FBgn%';  ## Wild type FBgn for query
my $fbalwc = 'FBal%'; ## Wild type FBal for query
## my %gh;  ## Hash of gene information  -- we report from this
## Main driver
my $gq = $dbh->prepare(sprintf("SELECT g.feature_id as gfid, g.uniquename as fbgn, g.name as gsym, fmin, fmax, strand, srcfeature_id from feature g natural join featureloc, organism o, cvterm cvt where g.type_id = cvt.cvterm_id and cvt.name = 'gene' and g.organism_id = o.organism_id and o.abbreviation = '%s' and g.uniquename like '%s' and g.is_obsolete = 'f' and g.is_analysis = 'f'",$org,$fbgnwc));
## Test driver (limited query)
#my $gq = $dbh->prepare(sprintf("SELECT g.feature_id as gfid, g.uniquename as fbgn, g.name as gsym, fmin, fmax, strand, srcfeature_id from feature g natural join featureloc, organism o, cvterm cvt where g.type_id = cvt.cvterm_id and cvt.name = 'gene' and g.organism_id = o.organism_id and o.abbreviation = '%s' and g.uniquename like '%s' and g.is_obsolete = 'f' and g.is_analysis = 'f' limit 400",$org,$fbgnwc));
$gq->execute or die "WARNING: ERROR: Unable to execute gene query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
#     print "\nProcessing: $gr{fbgn}\t$gr{gsym}\n";

#
## Get transcript-associated expression/tap statements (count)
#
    my $tr_expression = 0;
    my $xrwc = '%XR';
    my $rq = $dbh7->prepare(sprintf("SELECT * from cvterm cvt, feature_relationship fr, feature f, feature_expression fe, expression e where f.feature_id = fr.subject_id and fr.type_id = cvt.cvterm_id and cvt.name = 'associated_with' and f.feature_id = fe.feature_id and fe.expression_id = e.expression_id and f.name like '%s' and fr.object_id = %d",$xrwc,$gr{gfid}));
    $rq->execute or die "WARNING: ERROR: Unable to execute XR transcript expression query\n";
    $tr_expression = $rq->rows;

#
## Get transcript-associated expression/tap statements (count)
#
    my $pr_expression = 0;
    my $xpwc = '%XP';
    my $rq = $dbh7->prepare(sprintf("SELECT * from cvterm cvt, feature_relationship fr, feature f, feature_expression fe, expression e where f.feature_id = fr.subject_id and fr.type_id = cvt.cvterm_id and cvt.name = 'associated_with' and f.feature_id = fe.feature_id and fe.expression_id = e.expression_id and f.name like '%s' and fr.object_id = %d",$xpwc,$gr{gfid}));
    $rq->execute or die "WARNING: ERROR: Unable to execute XP protein expression query\n";
    $pr_expression = $rq->rows;

#
## Get associated expressed reporter constructs & enhancer-traps  (counts)
#
    my $rc_count = 0;
    my $et_count = 0;
## Use query: /Users/emmert/work/psql/gene_transgenic_expression.sql for this; group by fr.type_id; get input from LC on appropriate breakdown
    my $eq = $dbh7->prepare(sprintf("SELECT subject_id, cvt2.name as frtype, object_id, f.feature_id as fid, f.uniquename as funame, f.name as fname, frp.value from feature f, feature_relationship fr, cvterm cvt1, cvterm cvt2, feature_relationshipprop frp, cvterm cvt3 where f.type_id = cvt1.cvterm_id and cvt1.name in ('transgenic_transposon','transposable_element_insertion_site') and f.feature_id = object_id and fr.feature_relationship_id = frp.feature_relationship_id and frp.type_id = cvt3.cvterm_id and cvt3.name = 'has_expression_data' and frp.value = 'Yes' and fr.type_id = cvt2.cvterm_id and cvt2.name in ('derived_assoc_vital_reporter_construct','derived_assoc_reporter_construct','derived_assoc_insertion_of_enhancer_trap','derived_assoc_insertion_of_enhancer_trap_binary_system') and subject_id = %d",$gr{gfid}));
    $eq->execute or die "WARNING: ERROR: Unable to execute reporter construct / enhancer-trap query\n";
    my $eq_cnt = $eq->rows;
    if ($eq_cnt > 0) {
	while (my %er = %{$eq->fetchrow_hashref}) {
	    if ($er{frtype} =~ /reporter\_construct/) {
		$rc_count++;
	    }
	    else {
		$et_count++;
	    }
	}
    }

#
## Get associated GO (molecular_function, biological_process, cellular_component) terms (counts)
#
    my %go;
##     print(sprintf("\tSELECT cv.name as cvname, cvt.name as goterm from feature_cvterm fc, cvterm cvt, cv where fc.feature_id = %d and fc.cvterm_id = cvt.cvterm_id and cvt.cv_id = cv.cv_id and cv.name in ('biological_process','molecular_function')\n",$gr{gfid}));
    my $goq = $dbh8->prepare(sprintf("SELECT cv.name as cvname, cvt.name as goterm from feature_cvterm fc, cvterm cvt, cv where fc.feature_id = %d and fc.cvterm_id = cvt.cvterm_id and cvt.cv_id = cv.cv_id and cv.name in ('biological_process','molecular_function','cellular_component')",$gr{gfid}));
    $goq->execute or die "WARNING: ERROR: Unable to execute GO ($t) query\n";
    while (my %gor = %{$goq->fetchrow_hashref}) {
	$go{$gor{cvname}}++;
    }
 

#
## Get allele information (count of alleles)
#
    my @alleles;
    my $allele_phenotypes = 0;
    my $aq = $dbh2->prepare(sprintf("SELECT a.feature_id as afid, a.uniquename as fbal, a.name as asym from feature a, organism o, cvterm cvt, cvterm cvt2, feature_relationship fr where a.type_id = cvt.cvterm_id and cvt.name = 'gene' and a.organism_id = o.organism_id and o.abbreviation = '%s' and a.is_obsolete = 'f' and a.is_analysis = 'f' and subject_id = a.feature_id and object_id = %d and fr.type_id = cvt2.cvterm_id and cvt2.name = 'alleleof'",$org,$gr{gfid}));
    $aq->execute or die "WARNING: ERROR: Unable to execute allele query\n";
    my $aq_cnt = $aq->rows;
    if ($aq_cnt > 0) {
	while (my %ar = %{$aq->fetchrow_hashref}) {
#	    print "\tAllele found: $ar{fbal}\t$ar{asym}\n";
	    push(@alleles,$ar{afid});
## For each allele, check whether theres any phenotype data
	    my $phq = $dbh6->prepare(sprintf("SELECT * from featureprop fp, cvterm cvt where fp.feature_id = %d and fp.type_id = cvt.cvterm_id and cvt.name in ('derived_allele_pheno_description','derived_pheno_class')",$ar{afid}));
	    $phq->execute or die "WARNING: ERROR: Unable to execute allele phenotype query\n";
	    my $phq_cnt = $phq->rows;
	    if ($phq_cnt > 0) {
		$allele_phenotypes++;
	    }
	}
    }
    else {
#	print "\tNo allele associated with this gene: $gr{fbgn}\t$gr{gsym}\n";
    }

#
## Get publication information (counts of primary refs and reviews, counts of post-1999 and post-2003 refs)
#
    my @papers;
    my @reviews;
    my $year_gt_1999 = 0;
    my $year_gt_2003 = 0;
    my $year_gt_2006 = 0;
    my $harv_flag_outstanding_2006 = 0;
    my $camb_flag_outstanding_2006 = 0;
    my $year_gt_2009 = 0;
    my $harv_flag_outstanding_2009 = 0;
    my $camb_flag_outstanding_2009 = 0;
#     print(sprintf("\tSELECT * from feature_pub fp, pub p, cvterm cvt where fp.feature_id = %d and fp.pub_id = p.pub_id and p.type_id = cvt.cvterm_id\n",$gr{gfid}));
    my $rq = $dbh3->prepare(sprintf("SELECT * from feature_pub fp, pub p, cvterm cvt where fp.feature_id = %d and fp.pub_id = p.pub_id and p.type_id = cvt.cvterm_id",$gr{gfid}));
    $rq->execute or die "WARNING: ERROR: Unable to execute pub query\n";
    my $rq_cnt = $rq->rows;
    if ($rq_cnt > 0) {
	while (my %rr = %{$rq->fetchrow_hashref}) {
#	    print "\tPub found: $rr{name}\t$rr{uniquename}\n";
	    if ($rr{name} eq 'paper') {
		push(@papers,$rr{uniquename});
		if ($rr{pyear} > 1999) {
		    $year_gt_1999++;
		}
		if ($rr{pyear} > 2003) {
		    $year_gt_2003++;
		}
		if ($rr{pyear} > 2006) {
		    $year_gt_2006++;
		    $harv_flag_outstanding_2006 += &check_triaged_nc('harv_flag',$rr{pub_id});
		    $camb_flag_outstanding_2006 += &check_triaged_nc('cam_flag',$rr{pub_id});
		}
		if ($rr{pyear} > 2009) {
		    $year_gt_2009++;
		    $harv_flag_outstanding_2009 += &check_triaged_nc('harv_flag',$rr{pub_id});
		    $camb_flag_outstanding_2009 += &check_triaged_nc('cam_flag',$rr{pub_id});
		}
	    }
	    elsif ($rr{name} eq 'review') {
		push(@reviews,$rr{uniquename});
	    }
	}
    }
    else {
#	print "\tNo pub associated with this gene: $gr{fbgn}\t$gr{gsym}\n";
    }

#
## Get genetic interaction count
#
    my $interactions = 0;
    foreach my $al (@alleles) {
##	my $iq = $dbh4->prepare(sprintf("SELECT * from featureprop fp, cvterm cvt where feature_id = '%d' and fp.type_id = cvt.cvterm_id and cvt.name in ('derived_suppressor_manifest','derived_non-enhanceable_manifest','derived_non-suppressible_manifest','derived_enhanceable_manifest','derived_enhancer_manifest','derived_suppressible_manifest','derived_non-suppressor_manifest','derived_other_manifest','derived_non-enhancer_manifest','derived_enhancer_class','derived_suppressible_class','derived_non-suppressor_class','derived_non-enhancer_class','derived_non-suppressible_class','derived_other_class','derived_non-enhanceable_class','derived_enhanceable_class','derived_suppressor_class','Derived_allele_xeno_interaction_comment','derived_allele_interaction_comment')",$al));
## Limit query to exclude "non_' interactions (per LC) -- these are the "not" ones bill wants excluded...?
	my $iq = $dbh4->prepare(sprintf("SELECT * from featureprop fp, cvterm cvt where feature_id = '%d' and fp.type_id = cvt.cvterm_id and cvt.name in ('derived_suppressor_manifest','derived_enhanceable_manifest','derived_enhancer_manifest','derived_suppressible_manifest','derived_other_manifest','derived_enhancer_class','derived_suppressible_class','derived_other_class','derived_enhanceable_class','derived_suppressor_class')",$al));
	$iq->execute or die "WARNING: ERROR: Unable to execute interaction query\n";
	$interactions = $interactions + $iq->rows;
    }
#    print "\tinteractions count is: $interactions\n";

#
## Get phenotype count
#
    my $phenotypes = 0;
    foreach my $al (@alleles) {
	my $pq = $dbh5->prepare(sprintf("SELECT distinct(value) from featureprop fp, cvterm cvt where feature_id = %d and fp.type_id = cvt.cvterm_id and cvt.name in ('derived_allele_pheno_description','derived_pheno_class','derived_pheno_manifest') and value != '\@FBcv0000349\:viable\@' and value != '\@FBcv0000374\:fertile\@'",$al));
	$pq->execute or die "WARNING: ERROR: Unable to execute phenotype query\n";
	$phenotypes = $phenotypes + $pq->rows;
    }

#
## Get the interpro count
#
    my $ipro = 0;
    my $ipq = $dbh6->prepare(sprintf("SELECT accession from feature_dbxref fd, dbxref dx, db where fd.feature_id = %d and fd.dbxref_id = dx.dbxref_id and dx.db_id = db.db_id and fd.is_current = 't' and db.name = 'INTERPRO'",$gr{gfid}));
    $ipq->execute or die "WARNING: ERROR: Unable to execute interpro query\n";
    my $ipq_cnt = $ipq->rows;
    if ($ipq_cnt > 0) {
      while (my %ipr = %{$ipq->fetchrow_hashref}) {
	$ipro++;
      }
    }

#
## Now print out the mess...
#
#    print(sprintf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$gr{fbgn},$gr{gsym},$armh{$gr{srcfeature_id}}{uname},$gr{fmin},$gr{fmax},$gr{strand},scalar(@alleles),$allele_phenotypes,scalar(@papers),scalar(@reviews),$year_gt_1999,$year_gt_2003,$year_gt_2006,$triaged_not_curated_2006,$year_gt_2009,$triaged_not_curated_2009,$interactions,$phenotypes,$tr_expression,$pr_expression,$rc_count,$et_count,$ipro,$go{biological_process},$go{molecular_function},$go{cellular_component}));
    print(sprintf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$gr{fbgn},$gr{gsym},$armh{$gr{srcfeature_id}}{uname},$gr{fmin},$gr{fmax},$gr{strand},scalar(@alleles),$allele_phenotypes,scalar(@papers),scalar(@reviews),$year_gt_1999,$year_gt_2003,$year_gt_2006,$harv_flag_outstanding_2006,$camb_flag_outstanding_2006,$year_gt_2009,$harv_flag_outstanding_2009,$camb_flag_outstanding_2009,$interactions,$phenotypes,$tr_expression,$pr_expression,$rc_count,$et_count,$ipro,$go{biological_process},$go{molecular_function},$go{cellular_component}));
#    print(sprintf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",$gr{fbgn},$gr{gsym},$year_gt_2006,$harv_flag_outstanding_2006,$camb_flag_outstanding_2006,$year_gt_2009,$harv_flag_outstanding_2009,$camb_flag_outstanding_2009));

}

$jetzt = scalar localtime;
print "## Finished report_gene_assessment: $jetzt\n\n";

sub check_triaged_nc {
## Given a pub_id, checks whether the pub is one which has been triaged but not 
## fully curated.   Returns 1 if the paper has been triaged but not fully 
## curated, otherwise returns 0

    my $flag = @_[0];  ##  One of: harv_flag, cam_flag
    my $pub_id = @_[1];


    my $has_flag = 0;
    my $is_curated = 0;

## For cam_flags, we simply count whether or not there is a non-excluded curated_by flag, which we assume will be a 
## camcur flag, and that it indicates that the entire paper was curated (from Cam's POV, at least)
    if ($flag eq 'cam_flag') {
	my $pq = $dbh2->prepare(sprintf("SELECT * from pubprop pp, cvterm cvt where pp.type_id = cvt.cvterm_id and cvt.name in ('cam_flag','curated_by') and pp.pub_id = %d",$pub_id));
	$pq->execute or die "WARNING: ERROR: Unable to execute pubprop query\n";
	my $pq_cnt = $pq->rows;
	if ($pq_cnt > 0) {
	    while (my %pr = %{$pq->fetchrow_hashref}) {
#		print "\t$pr{name}\t$pr{value}\n";
		if ($pr{name} =~ /harv|cam/) {
		    if ($pr{value} !~ /nocur/) {
			$has_flag++;
		    }
		}
		if ($pr{name} =~ /curated/) {
		    if ($pr{value} !~ /Leyland|bibl|skim|Author\sSubmission/) {
			$is_curated++;
		    }
		}
	    }
	}
	else {
	    return(0);
	}

    }
## For harv_flags, we check each harv_flag value to see if it has a DONE (since Harvcur doesnt curate whole paper, but rather dataset-by-dataset).
## If there are more harv_flags than DONE harv_flags, we count this as having outstanding harv-flags
    elsif ($flag eq 'harv_flag') {
	my $pq = $dbh2->prepare(sprintf("SELECT * from pubprop pp, cvterm cvt where pp.type_id = cvt.cvterm_id and cvt.name = 'harv_flag' and pp.pub_id = %d",$pub_id));
	$pq->execute or die "WARNING: ERROR: Unable to execute pubprop query\n";
	my $pq_cnt = $pq->rows;
	if ($pq_cnt > 0) {
	    while (my %pr = %{$pq->fetchrow_hashref}) {
		if (($pr{value} !~ /no\_flag/) && ($pr{value} !~ /DONE/)) {
		    $has_flag++;
		}
		else {
		    $is_curated++;
		}
	    }
	    if ($has_flag > $is_curated) {
		return(1);
	    }
	}
	else {
	    return(0);
	}
    }
    else {
	print "WARNING: ERROR: Invalid value in check_triaged_ac: $flag\n";
	exit();
    }

    if (($has_flag > 0) && ($is_curated == 0)) {
	return(1);
    }
    return(0);

}
