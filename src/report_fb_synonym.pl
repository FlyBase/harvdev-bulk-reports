#!/usr/local/bin/perl

# fb_synonym
#
#       Perl script generates tab-delimited synonym table with the following fields:
#       FBid, current_symbol, current_fullname, synonym(s)
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
# require "/users/emmert/work/perl-sub/conversions.pm";

if ( @ARGV < 5 ) {
    print
"\n USAGE: fb_synonym pg_server db_name pg_username pg_password output_filename\n\n";
    print "\toutput_filename is the output file for std and error output.\n\n";
    exit;
}

my $server = shift(@ARGV);
my $db     = shift(@ARGV);
my $user   = shift(@ARGV);
my $pwd    = shift(@ARGV);

$OUTFILE = shift(@ARGV);
open( STDERR, ">>$OUTFILE" );
open( STDOUT, ">>$OUTFILE" );

#
## DB Connections
#
## Data source (system)
my $dsource = sprintf( "dbi:Pg:dbname=%s;host=%s;port=5432", $db, $server );
$dbh = DBI->connect( $dsource, $user, $pwd )
  or die "cannot connect to $data_source\n";
$dbh2 = DBI->connect( $dsource, $user, $pwd )
  or die "cannot connect to $data_source\n";
$dbh3 = DBI->connect( $dsource, $user, $pwd )
  or die "cannot connect to $data_source\n";

## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
# $dbh = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh2 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";
# $dbh3 = DBI->connect($dsource,$user) or die "cannot connect to $data_source\n";

#
##  General setup
#
$jetzt = scalar localtime;
print
  "\n## FlyBase Symbol-Synonym Correspondence Table\n## Generated: $jetzt\n";
print "## Using datasource: $dsource...\n\n";
print
"##primary_FBid\torganism_abbreviation\tcurrent_symbol\tcurrent_fullname\tfullname_synonym(s)\tsymbol_synonym(s)\n";

## Hash of fbids and associated feature.type_id(s)
push( @{ $fbids{'FBgn'} }, "'gene'" );
push( @{ $fbids{'FBtr'} }, "'mRNA'" );
push( @{ $fbids{'FBpp'} }, "'polypeptide'" );
push( @{ $fbids{'FBal'} }, "'allele'" );
push( @{ $fbids{'FBab'} }, "'chromosome_structure_variation'" );
push( @{ $fbids{'FBba'} }, "'chromosome_structure_variation'" );
push( @{ $fbids{'FBti'} },
    "'transposable_element'", "'transposable_element_insertion_site'" );
push(
    @{ $fbids{'FBtp'} },
    "'transgenic_transposable_element'",
    "'natural_transposon'"
);

## Jira DB-630: suppress reporting of unlocalized expression features (e.g., wg-XR, hh-PP).
## Array of FB-ID prefixes that require filtering for localized features only.
my @loc_feat_type = ( 'FBtr', 'FBpp' );

#
# Main method
#
## For each FBid, we get basic feature information...
foreach my $id ( keys(%fbids) ) {
    my $idwc = $id . '%';

    #    print(sprintf("\nProcessing %s (%s)\n",$id,join(",",@{$fbids{$id}})));
    if ( $id ~~ @loc_feat_type ) {
        our $fq = $dbh->prepare(
            sprintf(
"SELECT f.feature_id, f.uniquename, f.name, abbreviation from feature f, featureloc fl, organism o, cvterm cvt where f.organism_id = o.organism_id and f.type_id = cvt.cvterm_id and cvt.name in (%s) and f.is_obsolete = 'f' and f.is_analysis = 'f' and f.feature_id = fl.feature_id and uniquename like '%s'",
                join( ",", @{ $fbids{$id} } ), $idwc
            )
        );

    }
    else {
        our $fq = $dbh->prepare(
            sprintf(
"SELECT f.feature_id, f.uniquename, f.name, abbreviation from feature f, organism o, cvterm cvt where f.organism_id = o.organism_id and f.type_id = cvt.cvterm_id and cvt.name in (%s) and f.is_obsolete = 'f' and f.is_analysis = 'f' and uniquename like '%s'",
                join( ",", @{ $fbids{$id} } ), $idwc
            )
        );
    }
    $fq->execute
      or die "WARNING: ERROR: Unable to execute feature query for $idwc\n";
    while ( my %fr = %{ $fq->fetchrow_hashref } ) {

        #	print "\tProcessing: $fr{feature_id}\t$fr{uniquename}\t$fr{name}\n";
## Get fullname (if any) and synonym(s) (if any) for each record
        my @syns;
        my @fullnames;
        my $curr_symbol;
        my $curr_fname;
        my $sq = $dbh2->prepare(
            sprintf(
"SELECT DISTINCT s.name as sname, s.synonym_sgml as synonym_sgml, st.name as stype, fs.is_current from feature_synonym fs, synonym s, cvterm st where fs.feature_id = %d and fs.synonym_id = s.synonym_id and s.type_id = st.cvterm_id and st.name in ('symbol','fullname')",
                $fr{feature_id} )
        );
        $sq->execute or die "WARNING: ERROR: Unable to execute synonym query\n";
        my $sq_cnt = $sq->rows;

        if ( $sq_cnt > 0 ) {
            while ( my %sr = %{ $sq->fetchrow_hashref } ) {
                next if ( $sr{sname} eq 'unnamed' );
##		print "\t\tsynonym: $sr{sname}\t$sr{stype} ($sr{is_current})\n";
                if ( ( $sr{is_current} == 1 ) && ( $sr{stype} eq 'symbol' ) ) {
                    $curr_symbol = $sr{sname};
                    $curr_symbol =~ s/<up>/[/g;
                    $curr_symbol =~ s/<\/up>/]/g;
                    $curr_fname  =~ s/<down>/[[/g;
                    $curr_fname  =~ s/<\/down>/]]/g;

                    #		    print "\t\tCURRENT symbol: $sr{sname}\t$sr{stype}\n";
                }
                elsif (( $sr{is_current} == 1 )
                    && ( $sr{stype} eq 'fullname' ) )
                {
                    $curr_fname = $sr{sname};
                    $curr_fname =~ s/<up>/[/g;
                    $curr_fname =~ s/<\/up>/]/g;
                    $curr_fname =~ s/<down>/[[/g;
                    $curr_fname =~ s/<\/down>/]]/g;

                  #		    print "\t\tCURRENT fullname: $sr{sname}\t$sr{stype}\n";
                }
                elsif (( $sr{is_current} == 0 )
                    && ( $sr{stype} eq 'fullname' ) )
                {
                    my $synonym_string = trim_quote_chars($sr{sname});
                    push( @fullnames, $synonym_string );
                    if ( $sr{synonym_sgml} ne $sr{sname} ) {
                        $sr{synonym_sgml} =~ s/<up>/[/g;
                        $sr{synonym_sgml} =~ s/<\/up>/]/g;
                        $sr{synonym_sgml} =~ s/<down>/[[/g;
                        $sr{synonym_sgml} =~ s/<\/down>/]]/g;
                        my $synonym_sgml_string = trim_quote_chars($sr{synonym_sgml});
                        push( @fullnames, $synonym_sgml_string );
                    }
                }
                elsif ( ( $sr{is_current} == 0 ) && ( $sr{stype} eq 'symbol' ) )
                {
                    my $synonym_string = trim_quote_chars($sr{sname});
                    push( @syns, $synonym_string );
                    if ( $sr{synonym_sgml} ne $sr{sname} ) {
                        $sr{synonym_sgml} =~ s/<up>/[/g;
                        $sr{synonym_sgml} =~ s/<\/up>/]/g;
                        $sr{synonym_sgml} =~ s/<down>/[[/g;
                        $sr{synonym_sgml} =~ s/<\/down>/]]/g;
                        my $synonym_sgml_string = trim_quote_chars($sr{synonym_sgml});
                        push( @syns, $synonym_sgml_string );
                    }

                   #		    print "\t\tpushing synonym: $sr{sname}\t$sr{stype}\n";
                }
            }
        }
        else {
            print "\t\tno synonyms linked to this gene...\n";
        }
        print(
            sprintf(
                "%s\t%s\t%s\t%s\t%s\t%s\n",
                $fr{uniquename},         $fr{abbreviation},
                $curr_symbol,            $curr_fname,
                join( "|", @fullnames ), join( "|", @syns )
            )
        );
    }
}

# This trims off flanking double-quote characters that are usually curated unintentionally.
sub trim_quote_chars{
    my $input_string = shift;
    print "Have this input: $input_string\n";
    if ( $input_string =~ /(^"|"$)/ ) {
        print "Trim flanking double-quote chars!\n";
        $input_string =~ s/(^"|"$)//g;
    }
    return $input_string;
}

## print "\nFinished fb_synonym $jetzt\n\n";
