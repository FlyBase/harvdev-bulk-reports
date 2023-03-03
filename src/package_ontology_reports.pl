#!/usr/local/bin/perl

# package_ontology_reports
#
#	Perl script makes tar archive of files created in 
#	bulk-processing and reporting build for Susan 
#	Tweedy to use as reference in doing ontology 
#	annotation
#
#-----------------------------------------------------------------------------#
#
#       NOTES:
#
#
#
#-----------------------------------------------------------------------------#

if (@ARGV != 1) {
    print "\n USAGE: package_ontology_reports db_name tarfile_filename\n\n";
    print "\ttarfile_filename is the name for the .tar archive of reports\n\n";
    exit;
}

# my $db = shift(@ARGV);
my $tfile = shift(@ARGV);




#
## General setup
#

# For docker, mount "'/data/build-public-release/' . $db . '/bulk_reports/'" to ./output/
my $repdir = '/src/output/';
print "\tMy bulk_reports directory is: $repdir\n";
my $tarfile = $repdir . $tfile;

## For docker, mount the bulk_processing directory to ./input/
my $proddir = '/src/input/';
print "\tMy bulk_processing directory is: $proddir\n";

#
## Main method
#

## Setup the first part of the system call: tar cvf 
my $tarcall = sprintf("tar cvf %s -C %s", $tarfile, $repdir);


print "Start looking here: $proddir\n";
opendir(my $pd, $proddir);
my @prodfiles = grep { /unknown\_links/ || /gene2go/  || /GO\_PubMed/ } readdir($pd);
closedir($pd);
foreach my $prodfile (@prodfiles) {
  print "\tGetting production file: $prodfile\n";
  my $mvcall = "cp $proddir/$prodfile $repdir";
  print "\tThe mvcall is: $mvcall\n";
  system($mvcall);
  $tarcall = $tarcall . " $prodfile";
}
print "\nThe tarcall so far is: $tarcall\n\n";

print "Next, look here: $repdir\n\n";
opendir(my $rd, $repdir);
my @gafiles = grep { /gene\_association.*\.fb/ || /gp_information\.fb/ } readdir ($rd);
close($rd);
foreach my $gafile (@gafiles) {
  print "\tGetting report file: $gafile\n";
  $tarcall = $tarcall . " $gafile";
}

print "\nThe final tarcall is: $tarcall\n\n";

system($tarcall);
