#!/usr/local/bin/perl
# report_insertion_mapping
#
#       Perl script generates tab-delimited insertion mapping table with the 
#       following fields: valid_symbol, FBti#, genomic_loc, range, orientation, 
#       estimated_cytogenetic_loc, observed_cytogenetic_loc
#
#-----------------------------------------------------------------------------#
#
#    NOTES
#    
#    
#
#-----------------------------------------------------------------------------#
use DBI;
require "conversions.pm";
# require "/users/emmert/work/perl-sub/conversions.pm";    # For local testing.

if (@ARGV < 5) {
    print "\n USAGE: report_insertion_mapping pg_server db_name pg_username pg_password output_filename\n\n";
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
##  DB Connections
#
## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $data_source\n";


## Data source (g4)
my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
# $dbh = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh2 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh3 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";


#
##  General setup
#
## Setup file header
$jetzt = scalar localtime;
print "## FlyBase Insertion Mapping Table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource\n\n";
print "##insertion_symbol\tFBti#\tgenomic_location\trange\torientation\tpublications\testimated_cytogenetic_location\tobserved_cytogenetic_location\n";

#
##
#

my $fbtiwc = 'FBti%';
## Main driver
my $iq = $dbh->prepare(sprintf("SELECT i.uniquename, i.name, i.feature_id, cvt.name as ftype from feature i, cvterm cvt where i.type_id = cvt.cvterm_id and i.is_obsolete = 'f' and is_analysis = 'f' and cvt.name = 'transposable_element_insertion_site' and uniquename like '%s' ORDER BY i.uniquename",$fbtiwc));
$iq->execute or die "WARNING: ERROR: Unable to execute TI query\n";
while (my %ir = %{$iq->fetchrow_hashref}) {
    ## Get estimated & observed cytogenetic locations
    my $est_cloc;
    my $cur_cloc;
    my $cq = $dbh2->prepare(sprintf("SELECT * from featureprop fp, cvterm cvt where fp.feature_id = %d and fp.type_id = cvt.cvterm_id and cvt.name in ('derived_genomic_cyto_location','curated_cytological_location')",$ir{feature_id}));
    $cq->execute or die "WARNING: ERROR: Unable to execute cyto fprop query\n";
    my $cq_rows = $cq->rows;
    if ($cq_rows > 0) {
        while (my %cr = %{$cq->fetchrow_hashref}) {
            if ($cr{name} eq 'derived_genomic_cyto_location') {
                $est_cloc = $cr{value};
                # print "\test_cloc is: .$est_cloc.\n";
            }
            elsif ($cr{name} eq 'curated_cytological_location') {
                $cur_cloc = $cr{value};
                # print "\tcur_cloc is: .$cur_cloc.\n";
            }
        }
    }
    ## Get genomic location (if any), orientation and range.
    my $gq = $dbh2->prepare(sprintf("SELECT DISTINCT featureloc_id, fmin, fmax, strand, is_fmin_partial, is_fmax_partial, g.uniquename as arm from featureloc fl, feature g  where fl.feature_id = %d and fl.srcfeature_id = g.feature_id",$ir{feature_id}));
    $gq->execute or die "WARNING: ERROR: Unable to execute genomic location query\n";
    my $gq_rows = $gq->rows;
    if ($gq_rows > 0) {
        while (my %gr = %{$gq->fetchrow_hashref}) {
            ## Process the location data.
            my $gloc;
            my $orien;
            my $range;
            ## Adjust for interbase
            if ($gr{fmin} != $gr{fmax}) { $gr{fmin}++; }
            $gloc = sprintf("%s:%s..%s",$gr{arm},$gr{fmin},$gr{fmax});
            # print "\tGenomic location: .$gloc.\n";
            ## Set range (t/f)
            if ($gr{fmin} != $gr{fmax}) {
                $range = 't';
            }
            else {
                $range = 'f';
            }
            ## Set strand to 1 when it isnt explicitly set.
            if (defined($gr{strand})) {  
                $orien = $gr{strand};
            }
            else {
                $orien = 0;
            }
            # QC check.
            if ((defined($ir{is_fmin_partial})) || (defined($ir{is_fmax_partial}))) {
                print "\tWARNING: unhandled partial floc not handled for: $ir{feature_id}\t$ir{uniquename}\t$ir{name}\n";
                exit();
            }
            ## Get featureloc attribution.
            my @pub_id_list = ();
            my $pub_id_string;
            my $pub_query = $dbh3->prepare(sprintf("SELECT DISTINCT uniquename FROM featureloc_pub flp JOIN pub p ON p.pub_id = flp.pub_id WHERE p.is_obsolete IS FALSE AND p.uniquename != 'unattributed' AND flp.featureloc_id = %d", $gr{featureloc_id}));
            $pub_query->execute or die "WARNING: ERROR: Unable to execute featureloc_pub query\n";
            while (my $row = $pub_query->fetchrow_hashref) {
                push @pub_id_list, $row->{uniquename};
            }
            if (@pub_id_list) {
                $pub_id_string = join('|', @pub_id_list);
            }
            print(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$ir{name},$ir{uniquename},$gloc,$range,$orien,$pub_id_string,$est_cloc,$cur_cloc));
        }
    }
    else {
        print(sprintf("%s\t%s\t\t\t\t\t%s\t%s\n",$ir{name},$ir{uniquename},$est_cloc,$cur_cloc));
    }
}

$jetzt = scalar localtime;
print "## finished report_insertion_mapping: $jetzt\n";

