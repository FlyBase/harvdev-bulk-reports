#!/usr/local/bin/perl

# report_gene_mapping
#
#       Perl script generates tab-delimited gene mapping table with the 
#       following fields:
#       gene_symbol, FBid, chromosome, recombination_loc, cytogenetic_loc, 
#       squence_loc
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

if (@ARGV < 5) {
    print "\n USAGE: report_gene_mapping pg_server db_name pg_username pg_password output_filename\n\n";
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
## DB connection
#
## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";

## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
# $dbh = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh2 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh3 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";


#
##  General setup
#

## Setup header
$jetzt = scalar localtime;
print "\n## FlyBase Gene Mapping Table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print "##organism_abbreviation\tcurrent_symbol\tprimary_FBid\trecombination_loc\tcytogenetic_loc\tsequence_loc\n";

## Array of drosophilids (taken from mol5_may_2006 where taxgroup = 'drosophilid')
## Once FBsp.obo is implemented we should get this data by querying the organism data in chado
## Added additional abbreviations from fb_2013_05 (October, 2013)
my %drosophilids_hash = ('Amag','blah','Apic','blah','Camo','blah','Cbic','blah','Ccos','blah','Cpro','blah','Deni','blah','Dabu','blah','Daca','blah','Dacu','blah','Dada','blah','Dadi','blah','Dadu','blah','Daff','blah','Dafs','blah','Dagl','blah','Dagu','blah','Dalb','blah','Dald','blah','Dalg','blah','Dali','blah','Dalo','blah','Dalp','blah','Damb','blah','Damx','blah','Dame','blah','Dana','blah','Danc','blah','Dand','blah','Dang','blah','Danm','blah','Dann','blah','Dano','blah','Dant','blah','Dara','blah','Darc','blah','Dari','blah','Darr','blah','Darx','blah','Darw','blah','Dark','blah','Darz','blah','Dasa','blah','Dass','blah','Dast','blah','Data','blah','Dath','blah','Dato','blah','Datr','blah','Datt','blah','Daur','blah','Dawd','blah','Dazt','blah','Dbax','blah','Dbai','blah','Dbak','blah','Dban','blah','Dbar','blah','Dbdu','blah','Dbia','blah','Dbic','blah','Dbif','blah','Dbii','blah','Dbim','blah','Dbio','blah','Dbip','blah','Dbir','blah','Dbis','blah','Dbiu','blah','Dbiz','blah','Dboa','blah','Dboc','blah','Dbol','blah','Dbok','blah','Dbor','blah','Dbos','blah','Dbpo','blah','Dbra','blah','Dbrb','blah','Dbrn','blah','Dbro','blah','Dbry','blah','Dbun','blah','Dbur','blah','Dbus','blah','Dbuz','blah','Dcad','blah','Dcai','blah','Dcal','blah','Dcam','blah','Dcan','blah','Dcap','blah','Dcar','blah','Dcau','blah','Dces','blah','Dcha','blah','Dchy','blah','Dcil','blah','Dcit','blah','Dcla','blah','Dcli','blah','Dcly','blah','Dcna','blah','Dcnd','blah','Dcni','blah','Dcnp','blah','Dcns','blah','Dcog','blah','Dcol','blah','Dcon','blah','Dcor','blah','Dcra','blah','Dcru','blah','Dcua','blah','Dcub','blah','Dcur','blah','Dcuv','blah','Dcvi','blah','Dcyr','blah','Ddac','blah','Ddar','blah','Ddas','blah','Ddav','blah','Dddx','blah','Dddt','blah','Dddu','blah','Ddef','blah','Ddes','blah','Ddfl','blah','Ddic','blah','Ddif','blah','Ddim','blah','Ddip','blah','Ddis','blah','Ddol','blah','Ddos','blah','Ddun','blah','Dele','blah','Dell','blah','Delp','blah','Dema','blah','Deng','blah','Deoh','blah','Dequ','blah','Derc','blah','Dere','blah','Derm','blah','Desk','blah','Deug','blah','Deur','blah','Dexi','blah','Dezo','blah','Dfal','blah','Dfas','blah','Dfic','blah','Dfim','blah','Dfin','blah','Dfla','blah','Dflg','blah','Dfle','blah','Dflp','blah','Dflv','blah','Dflx','blah','Dfor','blah','Dfra','blah','Dfrm','blah','Dfug','blah','Dful','blah','Dfum','blah','Dfun','blah','Dfuv','blah','Dfuy','blah','Dgan','blah','Dgas','blah','Dgau','blah','Dgib','blah','Dgnu','blah','Dgou','blah','Dgre','blah','Dgri','blah','Dgrs','blah','Dgua','blah','Dgum','blah','Dgun','blah','Dgur','blah','Dgut','blah','Dguy','blah','Dgym','blah','Dhaf','blah','Dhal','blah','Dham','blah','Dhan','blah','Dhaw','blah','Dhee','blah','Dhel','blah','Dhem','blah','Dhet','blah','Dhex','blah','Dhis','blah','Dhid','blah','Dhua','blah','Dhub','blah','Dhuc','blah','Dhui','blah','Dhyd','blah','Dhye','blah','Dhyp','blah','Diki','blah','Dima','blah','Dimm','blah','Dimp','blah','Dinc','blah','Dinf','blah','Ding','blah','Dini','blah','Dinn','blah','Dino','blah','Dins','blah','Diri','blah','Djam','blah','Dkal','blah','Dkan','blah','Dkep','blah','Dkha','blah','Dkik','blah','Dkit','blah','Dkoe','blah','Dkoh','blah','Dkun','blah','Dlac','blah','Dlae','blah','Dlah','blah','Dlai','blah','Dlar','blah','Dlas','blah','Dlat','blah','Dleb','blah','Dleo','blah','Dlet','blah','Dlev','blah','Dlie','blah','Dlii','blah','Dlix','blah','Dlim','blah','Dlin','blah','Dlit','blah','Dliu','blah','Dlng','blah','Dlon','blah','Dlog','blah','Dlow','blah','Dluc','blah','Dlue','blah','Dlum','blah','Dlus','blah','Dlut','blah','Dmac','blah','Dmad','blah','Dmaf','blah','Dmai','blah','Dmal','blah','Dmar','blah','Dmax','blah','Dmas','blah','Dmat','blah','Dmmm','blah','Dmpl','blah','Dmau','blah','Dmay','blah','Dmcr','blah','Dmdf','blah','Dmdi','blah','Dmdk','blah','Dmda','blah','Dmdl','blah','Dmdp','blah','Dmds','blah','Dmea','blah','Dmed','blah','Dmei','blah','Dmel','blah','Dmen','blah','Dmep','blah','Dmex','blah','Dmer','blah','Dmes','blah','Dmet','blah','Dmez','blah','Dmic','blah','Dmie','blah','Dmil','blah','Dmim','blah','Dmir','blah','Dmju','blah','Dmla','blah','Dmli','blah','Dmll','blah','Dmln','blah','Dmlp','blah','Dmls','blah','Dmlt','blah','Dmng','blah','Dmnl','blah','Dmob','blah','Dmoi','blah','Dmoj','blah','Dmon','blah','Dmor','blah','Dmpa','blah','Dmri','blah','Dmro','blah','Dmrp','blah','Dmsp','blah','Dmrt','blah','Dmul','blah','Dmun','blah','Dmut','blah','Dmys','blah','Dnaf','blah','Dnag','blah','Dnai','blah','Dnaj','blah','Dnan','blah','Dnar','blah','Dnas','blah','Dnst','blah','Dnav','blah','Dneb','blah','Dnec','blah','Dneh','blah','Dnem','blah','Dnen','blah','Dneo','blah','Dnep','blah','Dnek','blah','Dner','blah','Dnet','blah','Dnhy','blah','Dngo','blah','Dngr','blah','Dnic','blah','Dnid','blah','Dnig','blah','Dnih','blah','Dnik','blah','Dnil','blah','Dnim','blah','Dnin','blah','Dnis','blah','Dniv','blah','Dnom','blah','Dnov','blah','Dnpe','blah','Dnpi','blah','Doah','blah','Dobc','blah','Dobs','blah','Dobu','blah','Docc','blah','Docg','blah','Doch','blah','Docr','blah','Dogu','blah','Dohn','blah','Doka','blah','Dora','blah','Dore','blah','Dorf','blah','Dori','blah','Dorn','blah','Doro','blah','Dorp','blah','Dors','blah','Dort','blah','Dpaa','blah','Dpac','blah','Dpae','blah','Dpai','blah','Dpal','blah','Dpam','blah','Dpan','blah','Dpap','blah','Dpar','blah','Dpas','blah','Dpat','blah','Dpau','blah','Dpav','blah','Dpbp','blah','Dpch','blah','Dpcr','blah','Dpct','blah','Dpec','blah','Dpee','blah','Dpeg','blah','Dpen','blah','Dper','blah','Dpet','blah','Dpgt','blah','Dpha','blah','Dphe','blah','Dpic','blah','Dpil','blah','Dpin','blah','Dpiv','blah','Dpla','blah','Dpli','blah','Dpme','blah','Dpmo','blah','Dpol','blah','Dpos','blah','Dpra','blah','Dprc','blah','Dprg','blah','Dpri','blah','Dprl','blah','Dprm','blah','Dprn','blah','Dpro','blah','Dprp','blah','Dprs','blah','Dprt','blah','Dpru','blah','Dprv','blah','Dpsa','blah','Dpsp','blah','Dpsb','blah','Dpse','blah','Dpsi','blah','Dpst','blah','Dpsu','blah','Dpua','blah','Dpuc','blah','Dpug','blah','Dpul','blah','Dpun','blah','Dput','blah','Dpuu','blah','Dpvu','blah','Dqas','blah','Dqua','blah','Dqud','blah','Dqui','blah','Dqur','blah','Dqus','blah','Drac','blah','Draj','blah','Drct','blah','Drec','blah','Dred','blah','Drel','blah','Drep','blah','Dric','blah','Drit','blah','Drob','blah','Drua','blah','Drbf','blah','Drub','blah','Druf','blah','Drum','blah','Drur','blah','Dsal','blah','Dsan','blah','Dsbd','blah','Dsbs','blah','Dsch','blah','Dsci','blah','Dsec','blah','Dseg','blah','Dsei','blah','Dser','blah','Dset','blah','Dsex','blah','Dsfu','blah','Dsia','blah','Dsig','blah','Dsix','blah','Dsii','blah','Dsil','blah','Dsim','blah','Dsin','blah','Dsiv','blah','Dsla','blah','Dslb','blah','Dslf','blah','Dsoo','blah','Dsor','blah','Dsp','blah','Dspa','blah','Dspe','blah','Dspn','blah','Dspr','blah','Dsqu','blah','Dsri','blah','Dssa','blah','Dsta','blah','Dste','blah','Dsti','blah','Dstm','blah','Dsto','blah','Dstr','blah','Dsts','blah','Dstu','blah','Dsua','blah','Dsub','blah','Dsuc','blah','Dsue','blah','Dsul','blah','Dsuo','blah','Dsus','blah','Dsut','blah','Dsuz','blah','Dtak','blah','Dtai','blah','Dtal','blah','Dtan','blah','Dtei','blah','Dtes','blah','Dtet','blah','Dtex','blah','Dtol','blah','Dtra','blah','Dtre','blah','Dtri','blah','Dtrl','blah','Dtrn','blah','Dtro','blah','Dtrp','blah','Dtrt','blah','Dtru','blah','Dtrv','blah','Dtrz','blah','Dtsa','blah','Dtsi','blah','Dtsu','blah','Dunm','blah','Duni','blah','Dunp','blah','Duns','blah','Dust','blah','Dval','blah','Dvar','blah','Dvir','blah','Dvra','blah','Dvnz','blah','Dvou','blah','Dvul','blah','Dwad','blah','Dwas','blah','Dwat','blah','Dwhe','blah','Dwil','blah','Dwts','blah','Dyak','blah','Dyoo','blah','Dyun','blah','Dzot','blah','Exsp','blah','Gbiv','blah','Lacu','blah','Land','blah','Laer','blah','Lcla','blah','Lcol','blah','Lfen','blah','Lkur','blah','Lmac','blah','Lmag','blah','Lmik','blah','Lori','blah','Lpse','blah','Lsta','blah','Lten','blah','Msp1','blah','Msp2','blah','Mpoe','blah','Mshi','blah','Mtak','blah','Psp1','blah','Robe','blah','Sbro','blah','Schi','blah','Sadu','blah','Salb','blah','Sano','blah','Scy3','blah','Selm','blah','Sexi','blah','Sfla','blah','Sgra','blah','Sleo','blah','Soka','blah','Spal','blah','Spam','blah','Spat','blah','Sph1','blah','Zbad','blah','Zbog','blah','Zcap','blah','Zcer','blah','Zghe','blah','Zind','blah','Zine','blah','Zkol','blah','Zlin','blah','Zmas','blah','Zorn','blah','Zsep','blah','Zspi','blah','Ztar','blah','Ztub','blah','Zver','blah','Zvit','blah','Dptn','blah','Zdav','blah','Ztsa','blah','Zlac','blah','Zcam','blah','Dnml','blah','Deyp','blah','Dmme','blah','Dfma','blah','Dson','blah','Cpar','blah','Dflr','blah','Drio','blah','Dmos','blah','Dmow','blah','Dpon','blah','Drho','blah','Dspc','blah','Zngl','blah');


#
# Main method
#
## Get genes...                                                                 
my $fbgnwc = 'FBgn%';
## Main driver                                                                  
my $gq = $dbh->prepare(sprintf("SELECT f.uniquename, f.name, f.feature_id, cvt.name as ftype, abbreviation from feature f, organism o, cvterm cvt where f.organism_id = o.organism_id and f.type_id = cvt.cvterm_id and cvt.name = 'gene' and f.is_obsolete = 'f' and f.is_analysis = 'f' and f.uniquename like '%s'",$fbgnwc));
$gq->execute or die "WARNING: ERROR: Unable to execute gene query...\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
#    print "\nProcessing gene: $gr{feature_id}\t$gr{uniquename} / $gr{name}...\n";
    next if (!$gr{name});  ## #Exclude unnamed gene records
    next if (!$drosophilids_hash{$gr{abbreviation}});  ## Exclude non-drosophilid genes
## Exclude genes that are functional tags or engineered fusion genes
    my $isengineered = 'F';
    my $eq = $dbh4->prepare(sprintf("SELECT cvt.name from feature_cvterm fcv, feature_cvtermprop fcvp, cvterm cvt, cv, cvterm cvt2 where fcv.feature_id = %d and fcv.cvterm_id = cvt.cvterm_id and cvt.cv_id = cv.cv_id and cv.name = 'SO' and fcv.feature_cvterm_id = fcvp.feature_cvterm_id and fcvp.type_id = cvt2.cvterm_id and cvt2.name = 'gene_class'",$gr{feature_id}));
    $eq->execute or die "WARNING: ERROR: Unable to execute feature_cvterm query for SO terms\n";
    my $eq_cnt = $eq->rows;
    if ($eq_cnt > 0) {
	while (my %er = %{$eq->fetchrow_hashref}) {
#	    print "\tgene_class tag found: .$er{name}.\n";
	    if (($er{name} eq 'engineered_tag') || ($er{name} eq 'engineered_region') || ($er{name} eq 'engineered_fusion_gene')) {
#		print "\t\tSKIPPING ENGINEERED GENE...\n";
		$isengineered = 'T';
	    }
	}
    }
    next if ($isengineered eq 'T');

## Get chrom (not all genes will have this...)
    my $chrom;
    my $seqloc;
    my $fmin;
    my $fmax;
    my $strand;
    my $cq = $dbh2->prepare(sprintf("SELECT fmin, fmax, strand, c.uniquename from featureloc fl, feature c where fl.srcfeature_id = c.feature_id and fl.feature_id = %d",$gr{feature_id}));
    $cq->execute or die "WARNING: ERROR: Unable to execute genome location query\n";
    my $cq_cnt = $cq->rows;
    if ($cq_cnt > 0) {
	if ($cq_cnt > 1) {
	    print "\tWARNING: Multiple flocs found for this gene: $gr{feature_id}\t$gr{uniquename} / $gr{name}\n";
	}
	while (my %cr = %{$cq->fetchrow_hashref}) {
	    $chrom = $cr{uniquename};
	    $seqloc = sprintf("%s:%d..%d(%d)",$chrom,$cr{fmin}+1,$cr{fmax},$cr{strand});
	    $fmin = $cr{fmin};
	    $fmax = $cr{fmax};
	    $strand = $cr{strand};
	}
    }
##   else {
## Get the chromosome location from elsewhere...?
##    }

## Get recombination map location 
    my $recom;
    my $cyto;
    my $dcyto;
    my $rq = $dbh3->prepare(sprintf("SELECT value, cvt.name from featureprop fp, cvterm cvt where fp.feature_id = %d and fp.type_id = cvt.cvterm_id and cvt.name in ('promoted_genetic_location','inferred_cyto','promoted_gene_type','cyto_range','derived_computed_cyto')",$gr{feature_id}));
    $rq->execute or die "WARNING: ERROR: Unable to execute recomb map query\n";
    my $rq_cnt = $rq->rows;
    if ($rq_cnt > 0) {
	while (my %rr = %{$rq->fetchrow_hashref}) {
	    if ($rr{name} eq 'promoted_genetic_location') {
		$recom = $rr{value};
#		print "\trecom is: $recom\n";
	    }
	    elsif ($rr{name} eq 'inferred_cyto') {
		$cyto = $rr{value};
#		print "\tcyto is: $cyto\n";
	    }
	    elsif ($rr{name} eq 'cyto_range') {
		if (!$cyto) {
		    $cyto = $rr{value};
		}
	    }
	    elsif ($rr{name} eq 'derived_computed_cyto') {
	      $dcyto = $rr{value};
	      $dcyto =~ s/\;.*$//;
	    }
## Do we need to exclude certain types of genes?  If so...
	    elsif ($rr{name} eq 'promoted_gene_type') {  
		if ($rr{value} eq 'transposable_element') {
#		    print "\tThis is a transposable_element...\n";
		}
	    }
	  }
	if ((!$cyto) && ($dcyto)) {
	  $cyto = $dcyto;
	}
    }
## OUTPUT
    print(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n",$gr{abbreviation},$gr{name},$gr{uniquename},$recom,$cyto,$seqloc));
}

$jetzt = scalar localtime;
print "\n# Finished report_gene_mapping $jetzt\n\n";
