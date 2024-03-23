#!/usr/local/bin/perl -w

# script to build the GA file by querying GO data in chado producing gene_association.fb
# Usage: perl GA_file_builder host:dbname outputfile optional_error_log

use strict;
use DBI;
use LWP::UserAgent;
use LWP::Protocol::https;
use Data::Dumper;

if (@ARGV < 5) {
  print "\n USAGE: perl GA_file_builder server db user password outfile\n";
  exit();
}

my $start = localtime();
my $now = localtime();

my $server = shift @ARGV;
my $dbname = shift @ARGV;
my $user = shift @ARGV;
my $pass = shift @ARGV;
my $outfile = shift @ARGV;


# connect to db
################################ db connection ############################
$now = localtime();
print STDOUT "$now: DEBUG: About to connect to the database\n";
my $chado="dbi:Pg:dbname=$dbname; host=$server;port=5432";
my $dbh = DBI->connect($chado,$user,$pass) or die "cannot connect to $chado";
$now = localtime();
print STDOUT "$now: INFO: Connected to $chado\n";
##########################################################################

# hash to hold the FBrf to GO ref lookup
my %gorefs;
%gorefs = %{fetch_and_parse_gorefs()};

open(OUT, ">>$outfile") or die "Can't open $outfile to write output\n";

open(OUT2, ">$outfile.lines2review") or die "Can't open file to write report of lines to review\n";

open(OUT3, ">$outfile.pubs_wo_xrefs") or die "Can't open file to report of pubs without dbxrefs\n";

if (@ARGV) {
  my $logfile = shift @ARGV;
  open(STDERR, ">>$logfile") or die "Can't redirect STDERR to $logfile\n";
}

# deal with file header info
my $RELEASE = "$dbname";
$RELEASE = $1 if $dbname =~ /(fb_20[0-9]{2}_[0-9]{2}).*/;
$RELEASE =~ s/fb_/FB/g;

my $DATE = sprintf("%04d-%02d-%02d", sub {($_[5] + 1900, $_[4]+1, $_[3])}->(localtime));
my $header = '';
$header .= "!gaf-version: 2.2\n";
$header .= "!generated-by: FlyBase\n";
$header .= "!date-generated: $DATE\n";
$header .= "!saved-by: FlyBase GOcur gocur\@morgan.harvard.edu\n";
$header .= "!FlyBase release: $RELEASE\n";

print OUT $header;

# $dbh is a database handle that is used to talk to the db 
# the prepare method does all the formatting and such to get the query ready to execute

# I think the best way to deal with the synonyms at the moment is to set up a hash keyed by FBgn number with
# value a reference to an array of all the non_current synonyms for every gene 
# * note that if the go terms get moved to proteins and/or transcripts then the synonym bit will
# need to be rejiggered
my $syn_query = $dbh->prepare
  (sprintf
   ("SELECT distinct f.uniquename, s.name
     FROM   synonym s, feature_synonym fs, feature f
     WHERE  f.is_obsolete = false and f.is_analysis = false and f.uniquename like 'FBgn_______'
     and    f.feature_id = fs.feature_id and fs.synonym_id = s.synonym_id
     and    fs.is_current = false"
   )
  );

my %synonyms; # hash keyed by fbid with array of synonyms as value

# execute the query
$now = localtime();
print STDOUT "$now: INFO: Querying for synonyms\n";
$syn_query->execute or die "Can't fetch synonyms\n";

# retrieve the results and build the hash
while ((my $fbid, my $syn) = $syn_query->fetchrow_array()) {
  push @{$synonyms{$fbid}}, $syn;
}

# create a hash for fullname lookups so we don't need to do this in GA_query
# keyed by uniquename with current fullname as value
my %fullnames;

my $fn_query = $dbh->prepare
  (sprintf
   ("SELECT distinct f.uniquename, s.name
     FROM   synonym s, feature_synonym fs, feature f, cvterm c
     WHERE  f.is_obsolete = false and f.is_analysis = false
     and    f.feature_id = fs.feature_id and fs.synonym_id = s.synonym_id 
     and    fs.is_current = true
     and    s.type_id = c.cvterm_id and c.name = 'fullname'"
   )
  );

# execute the query
$now = localtime();
print STDOUT "$now: INFO: Querying for fullnames\n";
$fn_query->execute or die "Can't fetch fullnames\n";

while ((my $fbid, my $fn) = $fn_query->fetchrow_array()) {
  $fullnames{$fbid} = $fn;
}

# set up various mapping hashes

# pub types by type_id
my %pubsbytype;
my $stmt = "SELECT cvterm_id, cvterm.name FROM cvterm, cv WHERE cvterm.cv_id = cv.cv_id and cv.name = 'pub type'";
my $query = $dbh->prepare($stmt);
$query->execute or die "Can't do pub type query!\n";
while ((my $cid, my $name) = $query->fetchrow_array()) {
    $pubsbytype{$cid} = $name;
}

# this is global so can use it in the parse evidence code sub
our %EVC = ('inferred from mutant phenotype' => 'IMP',
	    'inferred from genetic interaction' => 'IGI',
	    'inferred from physical interaction' => 'IPI',
	    'inferred from sequence or structural similarity' => 'ISS',
	    'inferred from sequence model' => 'ISM',
	    'inferred from sequence alignment' => 'ISA',
	    'inferred from sequence orthology' => 'ISO',
	    'inferred from experiment' => 'EXP',
	    'inferred from direct assay' => 'IDA',
	    'inferred from electronic annotation' => 'IEA',
	    'inferred from expression pattern' => 'IEP',
	    'inferred from reviewed computational analysis' => 'RCA',
	    'traceable author statement' => 'TAS',
	    'non-traceable author statement' => 'NAS',
	    'inferred by curator'  => 'IC',
	    'inferred from genomic context' => 'IGC',
	    'no biological data available'  => 'ND',
	    'inferred from biological aspect of ancestor' => 'IBA',
	    'inferred from biological aspect of descendant' => 'IBD',
	    'inferred from key residues' => 'IKR',
	    'inferred from rapid divergence' => 'IRD',
	    'inferred from high throughput experiment' => 'HTP',
	    'inferred from high throughput direct assay' => 'HDA',
	    'inferred from high throughput expression pattern' => 'HEP',
	    'inferred from high throughput genetic interaction' => 'HGI',
	    'inferred from high throughput mutant phenotype' => 'HMP',
	  );

my %ASP = ('biological_process' => 'P',
	   'cellular_component' => 'C',
	   'molecular_function' => 'F',
	  );

my %TAX = (Dmel => '7227');


# prepare the biggie
my $ga_query = $dbh->prepare
  (sprintf
   ("SELECT DISTINCT gene.feature_id, gene.uniquename as fbid, symb.name as symbol,
       fcvt.feature_cvterm_id, cv.name as aspect,
       godb.accession as GO_id, ppub.uniquename as ppub, o.abbreviation as species,
       evc.value as evidence_code, prv.value as provenance, date.value as date, fcvt.is_not as is_not, evc.rank as evidence_code_rank
       FROM
       cvterm goterm
       JOIN   cv
       ON   (goterm.cv_id = cv.cv_id and goterm.is_obsolete = 0 and 
             cv.name in ('cellular_component','molecular_function','biological_process'))
       JOIN   dbxref godb
       ON   (goterm.dbxref_id = godb.dbxref_id)
       JOIN   feature_cvterm fcvt
       ON   (goterm.cvterm_id = fcvt.cvterm_id)
       JOIN   feature gene
       ON   (fcvt.feature_id = gene.feature_id and gene.is_obsolete = false and gene.is_analysis = false
            and gene.uniquename like 'FBgn_______')
       JOIN   organism o
       -- to remove organism restriction comment out 'and o.abbreviation = 'Dmel''
       -- but remember to close paren
       ON   (gene.organism_id = o.organism_id and o.abbreviation = 'Dmel')
       JOIN   feature_synonym fs1
       ON   (gene.feature_id = fs1.feature_id and fs1.is_current = true)
       JOIN   synonym symb
       ON   (fs1.synonym_id = symb.synonym_id)
       JOIN   cvterm stype
       ON   (symb.type_id = stype.cvterm_id and stype.name = 'symbol')
       JOIN   pub ppub
       ON   (fcvt.pub_id = ppub.pub_id)
       JOIN   feature_cvtermprop prv
       ON   (fcvt.feature_cvterm_id = prv.feature_cvterm_id)
       JOIN   cvterm prvname
       ON   (prv.type_id = prvname.cvterm_id and prvname.name = 'provenance')
       JOIN   feature_cvtermprop evc
       ON   (fcvt.feature_cvterm_id = evc.feature_cvterm_id)
       JOIN   cvterm evcname
       ON   (evc.type_id = evcname.cvterm_id and evcname.name = 'evidence_code')
       JOIN   feature_cvtermprop date
       ON   (fcvt.feature_cvterm_id = date.feature_cvterm_id)
       JOIN   cvterm dtname
       ON   (date.type_id = dtname.cvterm_id and dtname.name = 'date')"
   )
  );

# this is a separate query we have set up to retrieve qualifier information to stick
# in post big query to avoid row duplication problem (due to LEFT JOIN misbehavior) in previous version

# hash to store the feature_cvterms with qualifiers
my %quals;
# set up the query
my $qual_query = $dbh->prepare
  (sprintf
   ("SELECT fcvtp.feature_cvterm_id, qual.name
     FROM   feature_cvtermprop fcvtp, cvterm qual
     WHERE  fcvtp.type_id = qual.cvterm_id and
            qual.name IN ('enables', 'contributes_to', 'involved_in', 'acts_upstream_of',
                          'acts_upstream_of_positive_effect', 'acts_upstream_of_negative_effect',
                         'located_in', 'part_of', 'is_active_in', 'colocalizes_with')"
   )
  );
# build the quals hash
$now = localtime();
print STDOUT "$now: INFO: Getting the qualifier information.\n";
$qual_query->execute or die "Can't query for qualifiers\n";
while ( my ($fcvtid, $qual) = $qual_query->fetchrow_array()) {
  $quals{$fcvtid} = $qual;
}

# DB-893: Hash to store GO extensions: keys are feature_cvterm_id plus rank with intervening underscore char.
my %go_xtns;
my $go_xtn_query = $dbh->prepare
  (sprintf
    ("SELECT fcvtp.feature_cvterm_id, fcvtp.rank, fcvtp.value
      FROM feature_cvtermprop fcvtp
      JOIN cvterm cvt ON cvt.cvterm_id = fcvtp.type_id
      WHERE cvt.name = 'go_annotation_extension'"
    )
  );
$now = localtime();
print STDOUT "$now: INFO: Getting the GO extension information.\n";
my $go_xtn_counter = 0;
$go_xtn_query->execute or die "Can't query for GO extensions\n";
while ( my ( $fcvtid, $rank, $go_xtn_text ) = $go_xtn_query->fetchrow_array() ) {
  $go_xtns{$fcvtid . '_' . $rank} = $go_xtn_text;
  $go_xtn_counter += 1;
}
$now = localtime();
print STDOUT "$now: INFO: Found $go_xtn_counter GO annotation extensions.\n";

# here is a query to build a lookup hash to exclude genes annotated with 'transposable_element_gene' term
my %te_genes;

# here's the query
my $te_query = $dbh->prepare
  (sprintf
   ("SELECT f.uniquename
     FROM   feature f, feature_cvterm fc, cvterm c, cv
     WHERE  f.feature_id = fc.feature_id and fc.cvterm_id = c.cvterm_id
     and    c.cv_id = cv.cv_id and cv.name = 'SO' and c.name = 'transposable_element_gene'"));

$te_query->execute or die "Can't query for te genes\n";

while ((my $teg_uname) = $te_query->fetchrow_array()) {
  $te_genes{$teg_uname} = 1;
}

# we can prepare this  next query once and then bind the parameters of the FBrf
# in the place of the ? in the query as we go through the loop of returned results from the $ga_query
# this is in theory more efficient than both preparing and then executing the statement within the loop
# and I think in general we can do many fetches within the loop rather than setting up a lookup first
# as we did for synonyms because in this case only a minority subset of all references will be checked
# but it may turn out more efficient to do similar to synonyms (remains to be seen)
my $pmid_query = $dbh->prepare
  (sprintf
   ("SELECT CASE WHEN pt.name = 'supplementary material' 
            THEN (SELECT d.accession
                  FROM   pub_relationship pr, pub sm, pub p, pub_dbxref pd, dbxref d, db, cvterm c
                  WHERE  sm.uniquename = ? 
                    and ((sm.pub_id = pr.subject_id and pr.object_id = p.pub_id) or (sm.pub_id = pr.object_id and pr.subject_id = p.pub_id))  
                    and  pr.type_id = c.cvterm_id and c.name = 'related_to'
                    and  p.pub_id = pd.pub_id and pd.dbxref_id = d.dbxref_id and pd.is_current = true
                    and  d.db_id = db.db_id and db.name = 'pubmed')
            ELSE (SELECT d.accession
                  FROM   pub p, pub_dbxref pd, dbxref d, db
                  WHERE  p.uniquename = ? and p.pub_id = pd.pub_id
                    and  pd.dbxref_id = d.dbxref_id and pd.is_current = true
                    and  d.db_id = db.db_id and db.name = 'pubmed')
            END as PMID
      FROM pub, cvterm pt where pub.uniquename = ? and pub.type_id = pt.cvterm_id"
   )
  );

# query for gene type using 'promoted_gene_type' featureprop
# and map to value in col 12
#my %gene_types = ('protein_coding_gene' => 'protein',
#		  'miRNA_gene' => 'miRNA',
#		  'tRNA_gene' => 'tRNA',
#		  'rRNA_gene' => 'rRNA',
#		  'snoRNA_gene' => 'snoRNA',
#		  'snRNA_gene' => 'snRNA',
#		  'non_protein_coding_gene' => 'ncRNA',
#		 );

my %gene_types = ('mRNA' => 'protein',
		  'miRNA' => 'miRNA',
		  'pre_miRNA' => 'miRNA',
		  'tRNA' => 'tRNA',
		  'rRNA' => 'rRNA',
		  'snoRNA' => 'snoRNA',
		  'snRNA' => 'snRNA',
		  'ncRNA' => 'ncRNA',
		 );

(my $promoted_prop_tyid) = $dbh->selectrow_array
  (sprintf
   ("SELECT cvterm_id FROM cvterm c WHERE c.name = 'derived_gene_model_status'"));

my $tyq;
if (defined $promoted_prop_tyid) {
  $tyq = $dbh->prepare
  (sprintf
   ("SELECT value FROM featureprop WHERE type_id = $promoted_prop_tyid
        and feature_id = ?"));
}

my $partof = $dbh->selectrow_array
  (sprintf
   ("SELECT cvterm_id FROM cvterm c WHERE c.name = 'partof'"));

my $regex = '%-XR';

my $partyq = $dbh->prepare
  (sprintf
   ("SELECT distinct c.name
       FROM feature t, feature_relationship fr, cvterm c
      WHERE fr.type_id = $partof and fr.subject_id = t.feature_id
        and t.uniquename like 'FBtr_______' and t.name not like '%s' and t.is_obsolete = false
        and t.type_id = c.cvterm_id and fr.object_id = ?", $regex));

# execute the big query
$now = localtime();
print STDOUT "$now: INFO: Querying for ga info\n";
$ga_query->execute or die "Can't get ga info\n";

# fetch the results

# do we want to declare some hash or array to hold the info for sorting here? see below.
my %GA_results;
my %pseudos_withdrawn;

my %seen_pubs; # tracking hash so we don't need to query for so many pubs with complex query
my %seen_gns; # tracking hash with type as value
my $rows;

# note I am explicitly declaring all variables in results
# we could also just fetch an array and refer to array index numbers which might be slightly
# more efficient but harder to keep track of

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
# BOB: This "while" loop currently takes ~ 3h to process.
$now = localtime();
print STDOUT "$now: INFO: Processing results of GA query\n";
while ( my ($fid, $fbid, $symb, $fcvtid, $asp, $goid, $pub, $orgn, $evid, $src, $date, $is_not, $ev_rank)
	= $ga_query->fetchrow_array()) {
  $rows++;
  print STDOUT "$now: DEBUG: 1. Start processing row #$rows: feature_cvterm_id=$fcvtid\n";
  my $pseudoflag;
  my $extratabsflag;

  # skip te genes
  next if $te_genes{$fbid};

  # here is where some decisions need to be made on formatting
  # do you want the lines of the file sorted in any particular order for example by gene symbol
  # or go id
  # if so we will need to build up some type of data structure to be able to sort the values
  # it can be as simple as an array of lines but then the sorting function would be somewhat more
  # complicated or it could be some combination of hashes allowing sequential sorting of the keys

  # basically you will want to set up the line for the GA file as you like it and then either print 
  # it to your output or put it in the data structure for later sorting and printing
  my $line = "FB\t"; # here is the variable to store the GA file line with column one info

  # this bit here will do the query for the pub med id for the pub returned
  my $pmid;

  if (exists $seen_pubs{$pub}) {
    $pmid = $seen_pubs{$pub} if defined $seen_pubs{$pub};
    print STDOUT "$now: DEBUG: 2a. pub already seen.\n";
  } else {
    $pmid_query->bind_param(1, $pub);
    $pmid_query->bind_param(2, $pub);
    $pmid_query->bind_param(3, $pub);
    $pmid_query->execute or die "Can't fetch pub med id for $pub\n";
    print STDOUT "$now: DEBUG: 2b1. Queried for pub info.\n";
    # because any pub or supplementary info of a pub should only return one pubmed id we will do
    # a single fetch - but it is still possible that we won't find a PMID so check
    ($pmid) = $pmid_query->fetchrow_array();

    # add to seen hash if found
    if ($pmid) {
	$pmid = "PMID:$pmid";
	$seen_pubs{$pub} = $pmid;
    } else { 
	# do other checks and associations based on JIRA DB-233 specification
	$pmid = check4doi($dbh, $pub);
	if ($pmid) {
	    $pmid = "DOI:$pmid";
	    $seen_pubs{$pub} = $pmid;
	} else {
	    (my $ptyid) = $dbh->selectrow_array("SELECT type_id FROM pub WHERE uniquename = '$pub'");
	    if ($pubsbytype{$ptyid} eq 'review' or $pubsbytype{$ptyid} eq 'paper' or $pub eq 'FBrf0045941') {
		$pmid = 'GO_REF:0000095';
		print OUT3 "$pub|$pmid\n";
	    } elsif ($gorefs{$pub}) {
		$pmid = $gorefs{$pub};
	    } elsif ($pubsbytype{$ptyid} eq 'personal communication to FlyBase') {
		$pmid = 'GO_REF:0000097';
	    } elsif ($pubsbytype{$ptyid} eq 'abstract') {
		$pmid = 'GO_REF:0000098';
	    } elsif ($pubsbytype{$ptyid} eq 'DNA/RNA sequence record') {
		$pmid = 'GO_REF:0000099';
	    } elsif ($pubsbytype{$ptyid} eq 'protein sequence record') {
		$pmid = 'GO_REF:0000106';
	    } else {
		$pmid = undef;
		print OUT3 "$pub\n";
	    }
	    $seen_pubs{$pub} = $pmid;
	}
      # we can print out an optional warning if no PMID - comment out if no warning wanted
  #      print STDERR "NO PMID for $pub\n";
    }
    print STDOUT "$now: DEBUG: 2b2. Assessed query results for pub info.\n";

  }

  # going to need to parse the $evid field (i.e., feature_cvtermprop.value of type evidence) looking for with's and ands
  # note that there can be multiple evidence statements in the same property value
  # and one or more of them may have info for col 8
  # lets set up a sub to do the parsing and then figure out the best way to deal with
  # multiple lines
  my @cols_7_8 = parse_evidence_bits($evid, $dbh);


  # start building the line
  # cols 2 and 3
  $line .= "$fbid\t$symb\t";

  # handle negation (part of col 4).
  # $line .= "IS_NOT VALUE: $is_not |||";
  if ($is_not) {
      $line .= "NOT|";
  }
  # handle mandatory gp2term qualifiers (part of col 4).
  if (exists($quals{$fcvtid})) {
      $line .= "$quals{$fcvtid}\t";
  } else {
      print STDERR "CRITICAL ERROR: Missing gp2term qualifier for this annotation: fcvt_id = $fcvtid, gene = $symb ($fbid), cvterm = GO:$goid, pub = $pub, is_not = $is_not.\n";
      die "FAILING due to data issue: some annotations are missing gp2term qualifers.";
  }

  # cols 5 and 6
  $line .= "GO:$goid\tFB:$pub";

  # check for PMID
  $line .= "|$pmid" if $pmid;
  $line .= "\t";

  # evidence code col 7
  # and optional with col 8
  # NOTE: as these can have values that might need to be split over multiple lines 
  # at this point we just put in a place holder to substitute in the values
  $line .= "PUT_COLS_8_9_HERE\t";

  # aspect col 9
  $line .= "$ASP{$asp}\t";

  # optional fullname col 10
  if ($fullnames{$fbid}) {
    my $fn = $fullnames{$fbid};
    if ($fn =~ /\t/) {
      $extratabsflag = $fn;
      $fn =~ s/\t/ /g;
    }
    $line .= $fn;
    $fullnames{$fbid} = $fn;
  }
  $line .= "\t";

  # synonyms col 11
  if ($synonyms{$fbid}) {
    my @syns = sort @{$synonyms{$fbid}};
    foreach my $s (@syns) {
      next if ($s eq 'unnamed');
      if ($s =~ /\t/) {
	$extratabsflag = $s;
	$s =~ s/\t/ /g;
      }
      next if (($s eq $symb) or ($fullnames{$fbid} and ($s eq $fullnames{$fbid})));
      $line .= "$s|";
    }
    # this little bit just checks to make sure there are any synonyms left
    if ($line =~ /\|$/) {
      $line =~ s/\|$/\t/;
    } else {
      $line .= "\t";
    }
  } else {
    $line .= "\t";
  }

  # col 12 - determine type of transcript - assume only one but with miRNA exception
  # current GAF1 default = 'gene'
  #  $line .= "gene\t";

  # query for GAF2 format getting product type from promoted gene type
  # for product type - do lookup
  my $type = 'gene_product';
  if ($seen_gns{$fid}) {
    $type = $seen_gns{$fid};
    $type = 'gene_product' if $type eq 'pseudogene';

  } else {
    $partyq->bind_param(1, $fid);
    $partyq->execute or die "Can't do gene type query\n";

    if ($partyq->rows > 1) {
      print STDERR "WARNING - more than one type of transcript associated with $fbid\n";
    }
    (my $gtype) = $partyq->fetchrow_array();
  #    print "$gtype\n";
    if ($gtype) {
      if ($gtype eq 'pseudogene') {
	$pseudoflag = 1;
	$seen_gns{fid} = $gtype;
      } else {
	$type = $gene_types{$gtype} if $gene_types{$gtype};
	$seen_gns{$fid} = $type;
      }
    }
  }

  $line .= "$type\t";

  # col 13
  $line .= "taxon:$TAX{$orgn}\t";

  # col 14
  $line .= "$date\t";

  # col 15
  $line .= "$src\t";

  # col 16
  # handle GO annotation extensions.
  if (exists($go_xtns{$fcvtid . '_' . $ev_rank})) {
    $line .= "$go_xtns{$fcvtid . '_' . $ev_rank}\n";
  } else {
    $line .= "\n";
  }

  # add the evidence and with cols to as many lines as needed
  print STDERR "MISSING EVIDENCE for:\nt$line" and next unless @cols_7_8;
  my $line_w_ev = '';
  my $mismatch_line = '';
  foreach my $c (@cols_7_8) {
    my $evc_gn_mismatch;
    if ($c =~ /PROBLEM/) {
      print STDERR "$c IS A PROBLEM FOR $line\n";
      next;
    }
    if ($c =~ /MISMATCH/) {
      $c =~ s/MISMATCH:(.*)$//;
      $evc_gn_mismatch = $1;
    }
    my $line_copy = $line;
    $line_copy =~ s/PUT_COLS_8_9_HERE/$c/;
    if ($evc_gn_mismatch) {
      my $mline = $line_copy;
      chomp $mline;
      $mismatch_line .= $mline . "evc_has_fbgn_symbol_mismatch ($evc_gn_mismatch)\n";
    }      
    $line_w_ev .= "$line_copy";
  }
  $line = $line_w_ev;

  # and here's where to decide what to do with the line
  #  print $line;
  # and we'd like to sort if first by gene symbol and then by pub so build a HOH
  push @{$GA_results{$symb}{$pub}}, $line;

  # generate lines for problem reports
  if ($pseudoflag) {
    my $rep_line = $line;
    $rep_line =~ s/\n/pseudogene\n/;
    push @{$pseudos_withdrawn{$symb}{$pub}}, $rep_line;
  } 

  if ($tyq) {
    $tyq->bind_param(1, $fid);
    $tyq->execute or die "Can't do gene type query\n";
    (my $status) = $tyq->fetchrow_array();
    if ($status eq 'Withdrawn') {
      my $rep_line = $line;
      $rep_line =~ s/\n/withdrawn\n/;
      push @{$pseudos_withdrawn{$symb}{$pub}}, $rep_line;
    }
  }
 
  if ($mismatch_line) {
    push @{$pseudos_withdrawn{$symb}{$pub}}, $mismatch_line;
  }

  if ($extratabsflag) {
    my $rep_line = $line;
    $rep_line =~ s/\n/TAB IN: $extratabsflag\n/;
    push @{$pseudos_withdrawn{$symb}{$pub}}, $rep_line;
  }    
}
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

# and here we output the sorted lines sorted first by symbol and then by FBrf number
$now = localtime();
print STDOUT "$now: INFO: Producing output file\n";
my $lcnt;
foreach my $s (sort keys %GA_results) {
  foreach my $r (sort keys %{$GA_results{$s}}) {
    foreach my $l (@{$GA_results{$s}->{$r}}) {
      print OUT $l;
      $lcnt++;
#      print $l;
    }
  }
}
print "PRINTED OUT $lcnt LINES FOR GA FILE\n";

my $rlcnt;
foreach my $s (sort keys %pseudos_withdrawn) {
  foreach my $r (sort keys %{$pseudos_withdrawn{$s}}) {
    foreach my $l (@{$pseudos_withdrawn{$s}->{$r}}) {
      print OUT2 $l;
      $rlcnt++;
#      print $l;
    }
  }
}
print "PRINTED OUT $rlcnt LINES FOR LINES TO REVIEW REPORT\n";

my $end = localtime();
print "STARTED: $start\nENDED:   $end\n";

#$dbh->finish();
#$dbh->disconnect();

# checks to see if pub has a DOI number
sub check4doi {
    my $dbh = shift;
    my $fbrf = shift;
    my $stmt = "SELECT accession FROM pub p, pub_dbxref pd, dbxref d, db WHERE p.pub_id = pd.pub_id and pd.dbxref_id = d.dbxref_id and pd.is_current = true and d.db_id = db.db_id and db.name = 'DOI' and p.uniquename = ?";
    my $query = $dbh->prepare($stmt);
    $query->bind_param(1, $fbrf);
    $query->execute or warn "Can't do doi query\n";
    (my $doi) = $query->fetchrow_array(); #expecting only one current one
    return $doi;
}

# For GO annotations, the "evidence_code" feature_cvtermprop.value starts with
# the actual evidence code (fully written out); this must be converted to the
# appropriate abbreviation for the GAF "evidence" column.
# The evidence code in the feature_cvtermprop.value may be followed by the word
# "with" and a list of cross-references to other objects (e.g., FB genes,
# UniProt IDs, etc); these xrefs will be reported in the GAF "with" column.
# The xref list can be a bit dirty - various separators, many instances of the
# the word "with", extra prefixes (e.g., FLYBASE) or suffices (e.g., the name
# of the FB gene) that should not be reported). All this extra stuff must be
# stripped out to report a comma-separated list of db:accession strings.
  sub parse_evidence_bits {
  my $evbit = shift;    # The $evbit corresponds to the evidence_code feature_cvtermprop.value passed to the sub.
  my $dbh = shift;
  my @evlines;
  
  # First, substitute any spelled out evidence_code(s) with the corresponding abbreviation(s).
  # Any instances of the word "with" in the written out ev code are removed.
  foreach my $code (sort keys %EVC) {
    $evbit =~ s/$code/$EVC{$code}/g;
  }
  
  # Second, split up the value in case there are many evidence codes in there.
  # It was once the case that the feature_cvtermprop.value could contain
  # many evidence code bits separated by the word "AND".
  # This no longer seems to be true, but we're keeping this approach just to
  # be safe.
  my @evidence = trim(split / AND /, $evbit);

  # Third, in each evidence code string, split the evidence code abbreviation
  # from the list of cross-references (the word "with" separates them).
  foreach my $e (@evidence) {
    my $line;
    (my $evc, my $dbxrefs) = trim(split / with | from /, $e);
    # Check that there is a recognized evidence code abbreviation.
    unless (grep $evc eq $_, values %EVC) {
      print STDERR "WARNING - unrecognized evidence code: $evc\n";
      push @evlines, "PROBLEM: $e";
      next;
    }

    # Fourth, get the list of xrefs using the get_dbxrefs() subroutine.
    (my $col8, my $mismatch) = get_dbxrefs($dbxrefs, $dbh, $evc) if $dbxrefs;

    # Fifth, combine the evidence code abbreviation and the xrefs as a
    # col7\tcol8 string and push to the array.
    if ($col8) {
      $line = "$evc\t$col8";
      $line = "${line}MISMATCH:$mismatch" if ($mismatch);
    } else {
      $line= "$evc\t";
    }
    push @evlines, $line;
  }
  return @evlines;
}

# DB-919: Ensure that all IDs are returned (not just the first), and, that any
# extraneous bits are stripped from the db:dbxref.accession string.
# For example, for "FLYBASE:rod; FB:FBgn0003268,FLYBASE:Zw10; FB:FBgn0004643",
# we should return "FB:FBgn0003268,FB:FBgn0004643",
# rather than "FB:FBgn0003268,FLYBASE:Zw10" (currently the case).
sub get_dbxrefs {
  my $inline = shift;
  my $dbh = shift;
  my $evc = shift;
  my $outline = '';
  my $nomatch;

  # Split out xrefs into an array, trimming flanking white space.
  # NOTE - comma may or may not have a trailing space.
  my @dbxrefs = trim(split /,/, $inline);

  # Prepare a query that confirms that FB IDs are still current.
  my $stmt = "SELECT feature_id FROM feature WHERE is_obsolete = false and name = ? and uniquename = ?";
  my $query = $dbh->prepare($stmt);

  # Clean up each xref in the list.
  # Special handling for FB xrefs - restricted to genes, having this pattern:
  # FLYBASE:feature.name; FB:feature.uniquename
  # e.g., "FLYBASE:rod; FB:FBgn0003268"
  # In evidence codes, the semi-colon is only used in this way; it's always followed by a space.
  foreach my $d (@dbxrefs) {
    my @parts = split /; /, $d;
    if ($parts[1]) {
      my $xref = trim($parts[1]);
      $outline .= "$xref|";
    } else {
      my $xref = $parts[0];
      $outline .= "$xref|";
    }

    # For FB xrefs, check that the name and uniquename match.
    if ($parts[0] =~ /FLYBASE:(.+)/) {
      my $symb = $1;
      $symb = decon($symb);
      print STDERR "WARNING - missing Symbol in evidence code line $inline\n" unless ($symb);
      if ($parts[1] =~ /FB:(FBgn[0-9]{7})/) {
	      my $fbgn = $1;
        # print "CHECKING FOR MATCH BETWEEN $symb and $fbgn\n";
        $query->bind_param(1, $symb);
        $query->bind_param(2, $fbgn);
        $query->execute or warn "Can't execute $stmt FOR $symb:$fbgn\n";
        unless ($query->rows() > 0) {
          $nomatch = $fbgn.':'.$symb;
          $now = localtime();
          print STDOUT "$now: WARNING: WE HAVE A MISMATCH for $symb:$fbgn\n";
        }
      } else {
        print STDERR "WARNING - No FBgn provided for $symb in evidence code line\n";
      }
    }
  }

  # Trim off the xref divider at end of line.
  $outline =~ s/\|$//;

  # For certain evidence codes, use commas instead of pipes.
  if ($evc =~ /HGI|IBA|IC|IGC|IGI|IMP|IPI|ISA|ISS/) {
    $outline =~ s/\|/,/g;
  }

  return( $outline, $nomatch);
}

# attempts to fetch GO.references file from geneontology.org site and parse to lookup hash
# keyed by FBrf with GO ref as value
sub fetch_and_parse_gorefs {
  my $ua = LWP::UserAgent->new;
  my $req = HTTP::Request->new( GET => 'https://raw.githubusercontent.com/geneontology/go-site/master/metadata/gorefs/README.md');
  my $response = $ua->request($req);
  my $page = $response->content;
  unless ($page =~ /^# GO REFs/) {
	  return;
  }
  my @lines = split /\n/, $page;
  my %fbrf2goref;
  my $goid = '';
  foreach my $l (@lines) {
    next if $l =~ /^#/; # skip comments
	  next if $l =~ /^\s*$/; # and blank lines
	$goid = $1 if ($l =~ /^\s\*\sid:\s\[(GO_REF:[0-9]+)\]/);
  	next if ($goid eq "GO_REF:0000033");
	  if ($l =~ /^\s\*\sext\sxref:\sFB:(FBrf[0-9]{7})/) {
	    my $fbrf = $1;
	    $fbrf2goref{$fbrf} = $goid;
	  }
  }
  $fbrf2goref{'FBrf0253064'} = 'GO_REF:0000115';      # DB-767
  # $fbrf2goref{'FBrf0253063'} = 'GO_REF:0000024';    # DB-823
  $fbrf2goref{'FBrf0255270'} = 'GO_REF:0000024';      # DB-823
  $fbrf2goref{'FBrf0254415'} = 'GO_REF:0000047';      # DB-811
  $fbrf2goref{'FBrf0258542'} = 'GO_REF:0000033';      # DB-928

  $now = localtime();
  print "$now: INFO: Constructed FBrf -> GO_REF Mapping:\n";
  print Dumper(\%fbrf2goref);
  return \%fbrf2goref;
}


# trims any leading or trailing white space
# can provide a string or an array of strings
sub trim {
    my @s = @_;
    for (@s) {s/^\s+//; s/\s+$//;}
    return wantarray ? @s : $s[0];
}

sub decon {
# Converts SGML-formatted symbols to 'symbol_plain' format (modified from conv_greeks)

    my $string = shift ;

    $string =~ s/&agr\;/alpha/g;
    $string =~ s/&Agr\;/Alpha/g;
    $string =~ s/&bgr\;/beta/g;
    $string =~ s/&Bgr\;/Beta/g;
    $string =~ s/&ggr\;/gamma/g;
    $string =~ s/&Ggr\;/Gamma/g;
    $string =~ s/&dgr\;/delta/g;
    $string =~ s/&Dgr\;/Delta/g;
    $string =~ s/&egr\;/epsilon/g;
    $string =~ s/&Egr\;/Epsilon/g;
    $string =~ s/&zgr\;/zeta/g;
    $string =~ s/&Zgr\;/Zeta/g;
    $string =~ s/&eegr\;/eta/g;
    $string =~ s/&EEgr\;/Eta/g;
    $string =~ s/&thgr\;/theta/g;
    $string =~ s/&THgr\;/Theta/g;
    $string =~ s/&igr\;/iota/g;
    $string =~ s/&Igr\;/Iota/g;
    $string =~ s/&kgr\;/kappa/g;
    $string =~ s/&Kgr\;/Kappa/g;
    $string =~ s/&lgr\;/lambda/g;
    $string =~ s/&Lgr\;/Lambda/g;
    $string =~ s/&mgr\;/mu/g;
    $string =~ s/&Mgr\;/Mu/g;
    $string =~ s/&ngr\;/nu/g;
    $string =~ s/&Ngr\;/Nu/g;
    $string =~ s/&xgr\;/xi/g;
    $string =~ s/&Xgr\;/Xi/g;
    $string =~ s/&ogr\;/omicron/g;
    $string =~ s/&Ogr\;/Omicron/g;
    $string =~ s/&pgr\;/pi/g;
    $string =~ s/&Pgr\;/Pi/g;
    $string =~ s/&rgr\;/rho/g;
    $string =~ s/&Rgr\;/Rho/g;
    $string =~ s/&sgr\;/sigma/g;
    $string =~ s/&Sgr\;/Sigma/g;
    $string =~ s/&tgr\;/tau/g;
    $string =~ s/&Tgr\;/Tau/g;
    $string =~ s/&ugr\;/upsilon/g;
    $string =~ s/&Ugr\;/Upsilon/g;
    $string =~ s/&phgr\;/phi/g;
    $string =~ s/&PHgr\;/Phi/g;
    $string =~ s/&khgr\;/chi/g;
    $string =~ s/&KHgr\;/Chi/g;
    $string =~ s/&psgr\;/psi/g;
    $string =~ s/&PSgr\;/Psi/g;
    $string =~ s/&ohgr\;/omega/g;
    $string =~ s/&OHgr\;/Omega/g;
    $string =~ s/\<\/down\>/\]\]/g;
    $string =~ s/\<down\>/\[\[/g;
    $string =~ s/\<up\>/\[/g;
    $string =~ s/\<\/up\>/\]/g;

    return $string;
}
