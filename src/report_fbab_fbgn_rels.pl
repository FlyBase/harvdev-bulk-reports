#!/usr/local/bin/perl
use strict;
use warnings;
use DBI;
use Encode;
binmode(STDOUT, ":utf8");

=head1

report_insertion_mapping

Perl script that generates tab-delimited aberration_experimental_gene_del_dup_data bulk report file

=cut


unless (@ARGV == 5) {

    print "\n USAGE: $0 pg_server db_name pg_username pg_password output_filename\n\n";
    exit;


}

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);

my $output_file_name = shift(@ARGV);


##  DB Connections
#
## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
my $dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


my $type_mapping = {

	'deletes' => 'completely deleted/disrupted (complementation)',
	'part_deletes' => 'partially deleted/disrupted (complementation)',
	'nondeletes' => 'not deleted/disrupted (complementation)',
	'molec_deletes' => 'completely deleted (molecular)',
	'molec_partdeletes' => 'partially deleted (molecular)',
	'molec_nondeletes' => 'not deleted/disrupted (molecular)',


	'duplicates' => 'completely duplicated (complementation)',
	'part_duplicates' => 'partially duplicated (complementation)',
	'nonduplicates' => 'not duplicated (complementation)',
	'molec_dups' => 'completely duplicated (molecular)',
	'molec_partdups' => 'partially duplicated (molecular)',
	'molec_nondups' => 'not duplicated (molecular)',



};

open my $output, ">", "$output_file_name"
	or die "Can't open $output_file_name ($!)\n";
binmode($output, ":utf8");

my @output_rows;
# get data
foreach my $fr_type (keys %{$type_mapping}) {

	my $data = &get_feature_relationship_details_by_frtype($dbh, $fr_type, 'FBgn', "object", 'FBab');

	foreach my $FBgn (keys %{$data}) {

		foreach my $FBab (keys %{$data->{$FBgn}->{f_r}}) {

			my $reference_list = join '|', sort keys %{$data->{$FBgn}->{f_r}->{$FBab}->{FBrf}};

			# remove unattributed
			$reference_list =~ s/unattributed$//;
			$reference_list =~ s/|$//;

			push @output_rows, "$FBgn\t$data->{$FBgn}->{symbol}\t$type_mapping->{$fr_type}\t$FBab\t$data->{$FBgn}->{f_r}->{$FBab}->{symbol}\t$reference_list\n";

		}
	}

}

# print output

# header info
print $output "## FlyBase experimental aberration gene deletion/duplication information file.\n";
print $output '## Generated: ' . scalar localtime . "\n";
print $output '## Using datasource: ' . $db . "\n";
print $output "##\n";
print $output "## Column descriptions:\n";
print $output "## gene_id: Current FlyBase identifier (FBgn#) of the gene.\n";
print $output "## gene_symbol: Current symbol of the gene.\n";
print $output "## type: Description of how the gene is affected by the aberration listed in the aberration_id/aberration_symbol column pair, determined experimentally (either by genetic complementation analysis or molecular mapping).\n";
print $output "## aberration_id: Current FlyBase identifier (FBab#) of the aberration.\n";
print $output "## aberration_symbol: Current symbol of the aberration.\n";
print $output "## references: FlyBase identifier(s) (FBrf#) of the source reference(s).\n";
print $output "##\n";
print $output "## KEY to type column:\n";
print $output "## completely deleted/disrupted (complementation): the gene is reported to be completely deleted/disrupted by the aberration, as determined by genetic complementation analysis.\n";
print $output "## partially deleted/disrupted (complementation): the gene is reported to be partially deleted/disrupted by the aberration, as determined by genetic complementation analysis.\n";
print $output "## not deleted/disrupted (complementation): the gene is reported not to be removed or broken by the aberration, as determined by genetic complementation analysis.\n";
print $output "## completely deleted (molecular): the gene is reported to be completely deleted by the aberration, as determined by molecular mapping.\n";
print $output "## partially deleted (molecular): the gene is reported to be partially deleted by the aberration, as determined by molecular mapping.\n";
print $output "## not deleted/disrupted (molecular): the gene is reported not to be removed or broken by the aberration, as determined by molecular mapping.\n";
print $output "## completely duplicated (complementation): the gene is reported to be fully duplicated within the aberration, as determined by genetic complementation analysis.\n";
print $output "## partially duplicated (complementation): the gene is reported to be partially duplicated within the aberration, as determined by genetic complementation analysis.\n";
print $output "## not duplicated (complementation): the gene is reported not to be duplicated within the aberration, as determined by genetic complementation analysis.\n";
print $output "## completely duplicated (molecular): the gene is reported to be fully duplicated within the aberration, as determined by molecular mapping.\n";
print $output "## partially duplicated (molecular): the gene is reported to be partially duplicated within the aberration, as determined by molecular mapping.\n";
print $output "## not duplicated (molecular): the gene is reported not to be duplicated within the aberration, as determined by molecular mapping.\n";
print $output "##\n";


print $output "#gene_id\tgene_symbol\ttype\taberration_id\taberration_symbol\treferences\n";

# the data
foreach my $line (sort @output_rows) {

	print $output $line;
}


close $output;



=head1 SUBROUTINE:
=cut

=head1

	Title:    get_feature_relationship_details_by_frtype
	Function: The get_feature_relationship_details_by_frtype subroutine gets all feature_relationships of a particular feature_relationship type, where the user specifies the expected FBid type of both the subject and object.  Attribution details are also stored.
	Example:  my $data = &get_feature_relationship_details_by_frtype($dbh, "deletes", 'FBab', "subject", 'FBgn');

Arguments:

- $dbh = database handle

- $fr_type = 'type' of feature_relationship

- $feature_type = FBid type of the features in one half of the feature_relationship (the one that you want to sort the output by). If a single type is required/expected, fill in the FBid 'type' e.g. FBal or FBti. If multiple types of feature are required/expected fill in 'FB' (the sql query uses 'like' so FB will get all types of feature whose uniquename starts with FB).

- $feature_position = 'position' in the feature_relationship of the list of the features given in 3rd argument (i.e. are they the 'subject' or 'object')

- $other_feature_type = FBid 'type' of the other feature that is expected in the other half of the feature_relationship. If a single type is required/expected, fill in the FBid 'type' e.g. FBal or FBti. If multiple types of feature are required/expected fill in 'FB' (the sql query uses 'like' so FB will get all types of feature whose uniquename starts with FB).


The $data reference structure that is returned stores the following information:

$data->{$feature_id}->{f_r}->{$related_object_id}->{FBrf}->{$FBrf}++;
$data->{$feature_id}->{f_r}->{$related_object_id}->{symbol} = $related_object_symbol;
$data->{$feature_id}->{symbol} = $feature_symbol;

=cut



sub get_feature_relationship_details_by_frtype {

	unless (@_ == 5) {

		die "Wrong number of parameters passed to the get_feature_relationship_details_by_frtype subroutine, something is wrong with the script\n";
	}

	my ($dbh, $fr_type, $feature_type, $feature_position, $other_feature_type) = @_;

	my $data = {};

	my $related_object_position;

	$feature_type = $feature_type . '%';
	$other_feature_type = $other_feature_type . '%';

	if ($feature_position eq 'subject') {

		$related_object_position = 'object';


	} elsif ($feature_position eq 'object') {

		$related_object_position = 'subject';

	} else {

		warn "The \$feature_position variable in get_feature_relationship_details_by_frtype subroutine must be one of 'subject' or 'object', but the value $feature_position was passed to the subroutine. This needs fixing before the script will work\n";
		return ($data);

	}



# 

	my $sql_query = sprintf("SELECT f1.uniquename, f1.name, f2.uniquename, f2.name, pub.uniquename from feature f1, feature_relationship fr, cvterm cvt, feature f2, feature_relationship_pub frp, pub where fr.type_id = cvt.cvterm_id and cvt.name = '%s' and fr.%s_id = f1.feature_id and fr.%s_id = f2.feature_id and f1.is_obsolete = 'f' and f2.is_obsolete = 'f' and f1.uniquename like '%s' and f2.uniquename like '%s' and frp.feature_relationship_id = fr.feature_relationship_id and frp.pub_id = pub.pub_id and pub.is_obsolete = 'f'", $fr_type, $feature_position, $related_object_position, $feature_type, $other_feature_type);


	my $db_query = $dbh->prepare($sql_query);
	$db_query->execute or warn "WARNING: ERROR: Unable to execute get_feature_relationship_details_by_frtype query ($!)\n";

	my $row_count = $db_query->rows;

	if ($row_count > 0) {
		while (my ($feature_id, $feature_symbol, $related_object_id, $related_object_symbol, $FBrf) = $db_query->fetchrow_array) {


			$data->{$feature_id}->{f_r}->{$related_object_id}->{FBrf}->{$FBrf}++;
			$data->{$feature_id}->{f_r}->{$related_object_id}->{symbol} = $related_object_symbol;
			$data->{$feature_id}->{symbol} = $feature_symbol;

		}
	}

	$db_query->finish();

	return ($data);


}