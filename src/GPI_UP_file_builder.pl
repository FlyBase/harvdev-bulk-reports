#!/usr/local/bin/perl -w

# script to build the GPI file for submission to protein2go group that contains UniProt
# IDs in the dbxref field - built for only annotated protein coding genes
# specification for this file can be found in comments in JIRA DC-331
# Usage: perl GPI_UP_file_builder server dbname user pass outfile

use strict;
use DBI;

# constant strings
my $FB = 'FB';
my $TAXID = 'taxon:7227';

if (@ARGV < 5) {
  print "\n USAGE: perl GPI_UP_file_builder server db user password outfile\n";
  exit();
}

my $start = localtime();

my $server = shift @ARGV;
my $dbname = shift @ARGV;
my $user = shift @ARGV;
my $pass = shift @ARGV;
my $outfile = shift @ARGV;

# connect to db
################################ db connection ############################
my $chado="dbi:Pg:dbname=$dbname; host=$server;port=5432";
my $dbh = DBI->connect($chado,$user,$pass) or die "cannot connect to $chado";
print STDOUT "Connected to $chado\n";
##########################################################################

open(OUT, ">>$outfile") or die "Can't open $outfile to write output\n";

if (@ARGV) {
  my $logfile = shift @ARGV;
  open(STDERR, ">>$logfile") or die "Can't redirect STDERR to $logfile\n";
}

# deal with file header info
my $RELEASE = 'unknown';
$RELEASE = $1 if $dbname =~ /(fb_20[0-9]{2}_[0-9]{2}).*/;

my $DATE = sprintf("%02d\/%02d\/%04d", sub {($_[4]+1, $_[3], $_[5] + 1900)}->(localtime));
my $header = '';
$header .= "!gpi-version: 1.2\n";
$header .= "!date: $DATE \$\n!from: FlyBase\n!saved-by: Helen Attrill hla28\@gen.cam.ac.uk\n!Data from FlyBase version $RELEASE\n";

print OUT $header;

## set up some queries to get various bits of info given a feature_id


# gene fullname
my $fullnameq = $dbh->prepare(
    ("SELECT DISTINCT s.name FROM feature_synonym fs, synonym s, cvterm c
       WHERE fs.synonym_id = s.synonym_id and fs.is_current = true
         and s.type_id = c.cvterm_id and c.name = 'fullname'
         and fs.feature_id = ?"));

# annotation id
my $annidq = $dbh->prepare(
    ("SELECT d.accession FROM feature_dbxref fd, dbxref d, db
       WHERE fd.dbxref_id = d.dbxref_id and fd.is_current = true
         and d.db_id = db.db_id and db.name = 'FlyBase Annotation IDs'
         and fd.feature_id = ?"));



# query for annotated protein-coding genes - results fetched in main loop
my $pcgq = $dbh->prepare(
    ("
    SELECT DISTINCT g.feature_id, g.uniquename, tty.name
    FROM feature g
    JOIN organism o ON o.organism_id = g.organism_id
    JOIN feature_relationship fr ON fr.object_id = g.feature_id
    JOIN feature t ON t.feature_id = fr.subject_id
    JOIN featureloc fl ON fl.feature_id = t.feature_id
    JOIN cvterm rty ON rty.cvterm_id = fr.type_id
    JOIN cvterm gty ON gty.cvterm_id = g.type_id
    JOIN cvterm tty ON tty.cvterm_id = t.type_id
    WHERE g.is_obsolete = false
      AND gty.name = 'gene'
      AND g.uniquename LIKE 'FBgn%'
      AND o.abbreviation = 'Dmel'
      AND t.is_obsolete = false 
      AND t.uniquename LIKE 'FBtr%'
      AND tty.name IN ('mRNA', 'rRNA', 'ncRNA', 'pre_miRNA', 'tRNA', 'snoRNA', 'snRNA')
      AND rty.name = 'partof'
    UNION
    SELECT DISTINCT g.feature_id, g.uniquename, 'gene_product'
    FROM feature g
    JOIN cvterm cvtg ON cvtg.cvterm_id = g.type_id
    LEFT OUTER JOIN featureloc fl ON fl.feature_id = g.feature_id
    JOIN organism o ON o.organism_id = g.organism_id
    JOIN feature_cvterm fcvt ON fcvt.feature_id = g.feature_id
    JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
    JOIN dbxref dbx ON dbx.dbxref_id = cvt.dbxref_id
    JOIN db ON db.db_id = dbx.db_id
    WHERE g.is_obsolete = false
      AND cvtg.name = 'gene'
      AND g.uniquename LIKE 'FBgn%'
      AND o.abbreviation = 'Dmel'
      AND fl.featureloc_id IS NULL
      AND cvt.is_obsolete = 0
      AND db.name = 'GO'"
    ));

# fetch the results
my %gpi_lines;
my $rows;

$pcgq->execute or die "Can't do PCG query\n";

print "Processing results of PCG query\n";
while ( my ($fid, $uniquename, $transcript_type) = $pcgq->fetchrow_array()) {
  $rows++;

  # get info from fid specific queries use direct select for those with only one 
  # anticipated value
  # current symbol
  (my $symb) = $dbh->selectrow_array(
    ("SELECT DISTINCT s.name FROM feature_synonym fs, synonym s, cvterm c
       WHERE fs.synonym_id = s.synonym_id and fs.is_current = true
         and s.type_id = c.cvterm_id and c.name = 'symbol'
         and fs.feature_id = $fid"));

  # annotation id
  (my $cgid) = $dbh->selectrow_array(
    ("SELECT d.accession FROM feature_dbxref fd, dbxref d, db
       WHERE fd.dbxref_id = d.dbxref_id and fd.is_current = true
         and d.db_id = db.db_id and db.name = 'FlyBase Annotation IDs'
         and fd.feature_id = $fid"));

  # gene fullname - may not be one
  my $fullname;
  ($fullname) = $dbh->selectrow_array(
    ("SELECT DISTINCT s.name FROM feature_synonym fs, synonym s, cvterm c
       WHERE fs.synonym_id = s.synonym_id and fs.is_current = true
         and s.type_id = c.cvterm_id and c.name = 'fullname'
         and fs.feature_id = $fid"));

  # sub to do logic to get wanted UniProts as it could change
  my $uline = '';
  if ( $transcript_type eq 'mRNA') {
    my @upids = get_uniprots($dbh, $fid);
    $uline .= "UniProtKB:$_|" for @upids;
    $uline =~ s/\|$//;
  }
  elsif ( $transcript_type eq 'gene_product') {
    my @upids = get_uniprots($dbh, $fid);
    $uline .= "UniProtKB:$_|" for @upids;
    $uline =~ s/\|$//;
  }


  else {
    my @rnacids = get_rnacentral_xrefs($dbh, $fid);
    $uline .= "RNAcentral:$_|" for @rnacids;
    $uline =~ s/\|$//;
  }

  # print the line to file
  print OUT "$FB\t$uniquename\t$symb\t";
  print OUT $fullname if $fullname;
  if ( $transcript_type eq 'mRNA') {
    print OUT "\t$cgid\tprotein\t$TAXID\t\t$uline\t\n";
  }
  elsif ( $transcript_type eq 'pre_miRNA') {
    print OUT "\t$cgid\tmiRNA\t$TAXID\t\t$uline\t\n";
  }
  else {
    print OUT "\t$cgid\t$transcript_type\t$TAXID\t\t$uline\t\n";
  }

}
my $end = localtime();
print "STARTED: $start\tENDED: $end\n";

# returns uniprot ids based on specification provided 
# may evolve 
# initially checks for SwissProt IDs associated with gene and returns if found
# if not then return all TrEMBL IDs
sub get_uniprots {
    my $dbh = shift;
    my $fid = shift;
    my $SP = 'UniProt/Swiss-Prot';
    my $TR = 'UniProt/TrEMBL';

    # Get Swiss-Prot first.
    my @sp = get_up_acc4id_by_db($dbh, $SP, $fid);

    # Get TrEMBL next.
    my @tr = get_up_acc4id_by_db($dbh, $TR, $fid);

    # Combine the two lists (Jira DB-752).
    push(@sp, @tr);
    return @sp;
}

# returns uniprot ids based on specification provided 
# may evolve 
# initially checks for SwissProt IDs associated with gene and returns if found
# if not then return all TrEMBL IDs
sub get_rnacentral_xrefs {
    my $dbh = shift;
    my $fid = shift;
    my $DB = 'RNAcentral';

    # Get RNAcentral IDs.
    my @rc = get_up_acc4id_by_db($dbh, $DB, $fid);

    return @rc;
}

sub get_up_acc4id_by_db {
    my $dbh = shift;
    my $db = shift;
    my $fid = shift;
    my $upq = $dbh->prepare
	("SELECT d.accession FROM dbxref d, feature_dbxref fd, db
          WHERE  fd.dbxref_id = d.dbxref_id and d.db_id = db.db_id 
            and  db.name = '$db' and fd.feature_id = $fid");
    $upq->execute or die "Can't do Uniprot query!\n";
    my @ups;
    while ((my $upid) = $upq->fetchrow_array()) {
	push @ups, $upid;
    }
    return @ups;
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
