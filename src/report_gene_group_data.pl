#!/usr/bin/perl
# report_gene_group_data
#
#       Perl script reports FlyBase FBgg gene group or pathway data for public
#       bulk report. The "-t" option restricts output to a given set of groups;
#       if the -t option is not used, groups in defined sets will be excluded
#       instead so that only non-pathway gene groups are reported.
#
#-----------------------------------------------------------------------------#
#
#       NOTES
#       
#       See JIRA DB-261, DB-794, DB-992 for more details
#       
#
#-----------------------------------------------------------------------------#
use DBI;
use Getopt::Long;

if (@ARGV <5) {
  print "\n USAGE: report_gene_group_data pg_server db_name pg_username pg_password report_output_filename (-t (signaling|metabolic))\n\n";
  exit();
}
# Look for pathway option.
my $grp_type = '';
GetOptions ("t=s" => \$grp_type);
my %grp_type_to_term = (
    'signaling' => 'signaling pathway group',
    'metabolic' => 'metabolic pathway group',
);
print "Running the " . __FILE__ . " script with the -t option set to grp_type=" . $grp_type . ".\n";
my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);

$OUTFILE = shift(@ARGV);
open(STDERR,">$OUTFILE");
open(STDOUT,">$OUTFILE");


#
##  DB Connections
#

# ## Data source (system)
my $dsource = sprintf("dbi:Pg:dbname=%s;host=%s;port=5432",$db,$server);
## Data source (g4)
# my $dsource = sprintf("dbi:Pg:dbname=%s;",$db);
$dbh = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh2 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh3 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh4 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";
$dbh5 = DBI->connect($dsource,$user,$pwd) or die "cannot connect to $dsource\n";


#
##  General setup
#

## Print header
$jetzt = scalar localtime;
print "## FlyBase $ggtype{$pathway} report\n## Generated: $jetzt\n## Using chado datasource: $dsource\n\n";
print "## Where groups are arranged into hierarchies, note that:\n";
print "## i) the member genes are only associated with the terminal subgroups\n";
print "## ii) the immediate parent of any subgroup is identified in the 'Parent_FB_group_id' and 'Parent_FB_group_symbol' columns\n\n";
print "# FB_group_id\tFB_group_symbol\tFB_group_name\tParent_FB_group_id\tParent_FB_group_symbol\tGroup_member_FB_gene_id\tGroup_member_FB_gene_symbol\n";

#
## Main method
#

## Populate hash %gnh of all grp grp_id, uniquename, name, and gene(s) if any
# Modify query depending on the $grp_type argument.
my $cvterm_query_bit;
if ($grp_type) {
    $cvterm_query_bit = "AND cvt.name = " . $dbh3->quote($grp_type_to_term{$grp_type});
} else {
    my @cvterms_to_exclude = values %grp_type_to_term;
    my $cvterms_to_exclude_str = join("', '", map { $dbh3->quote($_) } @cvterms_to_exclude);
    $cvterm_query_bit = "AND NOT cvt.name IN ('$cvterms_to_exclude_str')";
}

my %gnh;
my $grq = $dbh3->prepare(sprintf("
    SELECT g.grp_id,
           g.name,
           g.uniquename
    FROM grp g
    JOIN grp_cvterm gcvt ON gcvt.grp_id = g.grp_id
    JOIN cvterm cvt ON cvt.cvterm_id = gcvt.cvterm_id
    JOIN cv ON cv.cv_id = cvt.cv_id
    WHERE g.is_obsolete IS FALSE
      AND cv.name = 'FlyBase miscellaneous CV'
      %s
", $cvterm_query_bit));
$grq->execute or die "WARNING: ERROR: Unable to execute grp name/uniquename query\n";
while (my %grr = %{$grq->fetchrow_hashref}) {
#    print "Adding grp to gnh:\t$grr{grp_id}\t$grr{name}\t$grr{uniquename}\n";
    $gnh{$grr{grp_id}}{name} = $grr{name};
    $gnh{$grr{grp_id}}{uniquename} = $grr{uniquename};
}

## Add fullname for each grp record to hash %gnh 
foreach $rec (keys(%gnh)) {
    my $sq = $dbh4->prepare(sprintf("SELECT s.name as sname from grp_synonym gs, synonym s, cvterm cvt where gs.synonym_id = s.synonym_id and gs.is_current = 't' and gs.is_internal = 'f' and s.type_id = cvt.cvterm_id and cvt.name = 'fullname' and gs.grp_id = %d",$rec));
    $sq->execute or die "WARNING: ERROR: Unable to execute grp_synonym query\n";
    my $sq_cnt = $sq->rows;
    if ($sq_cnt > 0) {
	if ($sq_cnt > 1) {
	    print "\tWARNING: ERROR: Multiple fullname synonyms found for grp: $rec\t$gnh{$rec}{uniquename}\t$gnh{$rec}{name}\n";
	}
	while (my %sr = %{$sq->fetchrow_hashref}) {
#	    print "\tSetting fullname for grp $rec\t$gnh{$rec}{uniquename}\t$gnh{$rec}{name}:\t$sr{sname}\n";
	    $gnh{$rec}{fullname} = $sr{sname};
	}
    }
}

## Add associated gene(s) to grp records in %gnh (if any)
foreach my $rec (keys(%gnh)) {
    my $gnq = $dbh2->prepare(sprintf("SELECT gm.type_id, g.uniquename as grpuname, g.name as grpname, f.uniquename as guname, f.name as gname from grp g, grpmember gm, feature_grpmember fgm, feature f, cvterm cvt where g.grp_id = gm.grp_id and gm.type_id = cvt.cvterm_id and cvt.name = 'grpmember_feature' and g.grp_id = %d and gm.grpmember_id = fgm.grpmember_id and fgm.feature_id = f.feature_id and f.is_obsolete = false",$rec));
    $gnq->execute or die "WARNING: ERROR: Unable to execute gene query\n";
    $gnq_cnt = $gnq->rows;
    if ($gnq_cnt > 0) {
	while (my %gnr = %{$gnq->fetchrow_hashref}) {
	    $gnh{$rec}{genes}{$gnr{guname}} = $gnr{gname};
#	    print "Adding gene for grp $rec\t$gnr{grpuname}\t$gnr{grpname}:\t$gnr{guname}\t$gnr{gname}\n";
	}
    }
}

## Query all records from grp_relationship, making hash %grph keying subject_id to array of all object_id(s)
my %grph;
my $gq = $dbh->prepare("
    SELECT gr.*
    FROM grp_relationship gr
    JOIN cvterm cvt ON gr.type_id = cvt.cvterm_id
    JOIN grp g1 ON g1.grp_id = gr.subject_id
    JOIN grp g2 ON g2.grp_id = gr.object_id
    WHERE cvt.name = 'parent_grp'
      AND g1.is_obsolete IS FALSE
      AND g2.is_obsolete IS FALSE
");
$gq->execute or die "WARNING: ERROR: Unable to execute grp_relationship query\n";
while (my %gr = %{$gq->fetchrow_hashref}) {
#    print "Loading grp_relationship record: $gr{grp_relationship_id}:\t$gr{subject_id}\t$gr{object_id}\n";
    push(@{$grph{$gr{subject_id}}{top_grs}},$gr{object_id});
}

## Add any grp records that didnt have a record in grp_relationship.  These wont have a parent.
foreach my $rec (keys(%gnh)) {
    if (!$grph{$rec}) {
#	print "\tPushing non-grp_relationship record: $rec\n";
	push(@{$grph{$rec}{top_grs}},'');
    }
}

## For testing -- comment out later
# print "\n";
# foreach my $p (keys(%grph)) {
#     print(sprintf("Consolidated_parent (subject->object)\t%d\t%s\n",$p,join(",",sort(@{$grph{$p}{top_grs}}))));
# }


#
## Output report
#

## For each group
foreach my $p (keys(%grph)) {   
#     print(sprintf("Final parent(s) (subject->object)\t%d\t%s\n",$p,join(",",sort(@{$grph{$p}{top_grs}}))));
    foreach my $c (sort(@{$grph{$p}{top_grs}})) {  ## Report parent
# 	print "\tGetting genes for child: .$p.\n";
	if ($gnh{$p}{genes}) { ## And report associated gene
	    foreach my $g (keys(%{$gnh{$p}{genes}})) {
#		print "REPORT\t$p\t$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t$gnh{$c}{uniquename}\t$gnh{$c}{name}\t$g\t$gnh{$p}{genes}{$g}\n";
		print "$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t$gnh{$c}{uniquename}\t$gnh{$c}{name}\t$g\t$gnh{$p}{genes}{$g}\n";
	    }
	}
	elsif ($gnh{$p}{uniquename}) {
#	    print "REPORT_NOGENE(w/child)\t$p\t$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t$gnh{$c}{uniquename}\t$gnh{$c}{name}\t\t\n";
	    print "$gnh{$p}{uniquename}\t$gnh{$p}{name}\t$gnh{$p}{fullname}\t$gnh{$c}{uniquename}\t$gnh{$c}{name}\t\t\n";
	}
    }
}


$jetzt = scalar localtime;
print "\n\n## Finished $ggtype{$pathway} report:\t$jetzt\n";
exit();
