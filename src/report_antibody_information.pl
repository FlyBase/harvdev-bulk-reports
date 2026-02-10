#!/usr/local/bin/perl
use strict;
use warnings;
use DBI;
use Encode;
binmode(STDOUT, ":utf8");

=head1

antibody_information

Perl script that generates tab-delimited antibody_information bulk report file

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


open my $output, ">", "$output_file_name"
	or die "Can't open $output_file_name ($!)\n";
binmode($output, ":utf8");


# 1. get commercial antibody info

my $supplier_list = ['DSHB', 'Cell Signaling Technology'];

my $antibody_info = &get_commercial_antibody_info($dbh,$supplier_list);


# 2. add lab-generated antibody info


$antibody_info = &get_lab_antibody_info($dbh,$antibody_info);


# print output

# header info
print $output "## FlyBase antibody information file.\n";
print $output '## Generated: ' . scalar localtime . "\n";
print $output '## Using datasource: ' . $db . "\n";
print $output "##\n";
print $output "#gene_id\tgene_symbol\tantibody_source\tantibody_clonality\tpub_id\tcitation\tsupplier\tproduct_number\n";


foreach my $row (sort @{$antibody_info}) {

	print $output "$row\n";

}


close $output;


sub get_commercial_antibody_info {

=head1 SUBROUTINE:
=cut

=head1

	Title:    get_commercial_antibody_info
	Function: The get_commercial_antibody_info subroutine gets commercial antibody information for the list of suppliers provided in the second argument. This data has been provided by suppliers as linkout information which will have been added into the reporting chado instancec in feature_dbxref, dbxref and dbxrefprop.
	Example:  my $commercial_antibody_data = &get_commercial_antibody_info($dbh,$supplier_list);

Arguments:

- $dbh = database handle

- $supplier_list = list of suppliers (each entry corresponds to a db entry in chado)


=cut

	unless (@_ == 2) {

		die "Wrong number of parameters passed to the get_commercial_antibody_info subroutine, something is wrong with the script\n";
	}

	my ($dbh, $supplier_list) = @_;

	my $data;

	foreach my $supplier (@{$supplier_list}) {
	

	my $sql_query =sprintf("SELECT distinct f.uniquename, f.name, dbxp.value, db.name, dbx.accession FROM feature f, feature_dbxref fdbx, dbxref dbx, dbxrefprop dbxp, db WHERE db.name = '%s' and fdbx.is_current='true' and f.is_obsolete = 'f' and f.uniquename like '%s' and f.feature_id = fdbx.feature_id and fdbx.dbxref_id = dbx.dbxref_id and dbx.db_id=db.db_id and dbxp.dbxref_id = dbx.dbxref_id", $supplier, 'FBgn%');

	my $db_query = $dbh->prepare($sql_query);
	$db_query->execute or warn "WARNING: ERROR: Unable to execute get_commercial_antibody_info query ($!)\n";

		while (my ($FBid, $symbol, $clonality, $supplier, $product_number) = $db_query->fetchrow_array) {


			$clonality =~ s/\n/ /g;
			push @{$data}, "$FBid\t$symbol\tcommercial\t$clonality\t\t\t$supplier\t$product_number";


		}
	}

	return $data;

}



sub get_lab_antibody_info {

=head1 SUBROUTINE:
=cut

=head1

	Title:    get_lab_antibody_info
	Function: The get_lab_antibody_info subroutine gets lab-generated antibody information. This data has been curated by FlyBase from references.
	Example:  my $antibody_info = &get_lab_antibody_info($dbh,$antibody_info);

Arguments:

- $dbh = database handle

- $data = reference to an array into which the results are to be pushed (so can easily combine with commercial antibody info)


=cut



	unless (@_ == 2) {

		die "Wrong number of parameters passed to the get_lab_antibody_info subroutine, something is wrong with the script\n";
	}


	my ($dbh, $data) = @_;

	my $sql_query =sprintf("SELECT distinct f.uniquename, f.name, fp.value, pub.uniquename, pub.miniref from feature f, featureprop fp, cvterm cvt, featureprop_pub fpp, pub where f.uniquename like '%s' and f.feature_id = fp.feature_id and fp.type_id = cvt.cvterm_id and cvt.name = 'reported_antibod_gen' and f.is_obsolete = 'f' and fpp.featureprop_id = fp.featureprop_id and fpp.pub_id = pub.pub_id and pub.is_obsolete = 'f'", 'FBgn%');


	my $db_query = $dbh->prepare($sql_query);
	$db_query->execute or warn "WARNING: ERROR: Unable to execute get_lab_antibody_info query ($!)\n";

	while (my ($FBid, $symbol, $clonality, $FBrf, $miniref) = $db_query->fetchrow_array) {

		$clonality =~ s/\n/ /g;
		push @{$data}, "$FBid\t$symbol\tlab generated\t$clonality\t$FBrf\t$miniref";


	}

	return $data;

}