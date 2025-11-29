#!/usr/local/bin/perl

# run_bulk_report_scripts
#
#       Perl script calls bulk report scripts for public release generation
#       
#-----------------------------------------------------------------------------#
#
#       NOTES:
#
#
#
#-----------------------------------------------------------------------------#

if (@ARGV < 6) {
    print "\n USAGE: run_bulk_report_scripts pg_server db_name Dmel_release_designator pg_username pg_password output_filename\n\n";
    print "\trelease_designator is the release designator for Dmel release being built, e.g., 'R5.50'\n\n";
    print "\toutput_filename is the output file for std and error output.\n\n";
    exit;
}

print "\n\nDID THE SCRIPT FOR GENERATING APOLLO CHXML FOR IU GET INCORPORATED INTO THIS SCRIPT YET?  \n\n\n";

my $server = shift(@ARGV);
my $db = shift(@ARGV);
my $relno = shift(@ARGV);
my $user = shift(@ARGV);
my $pwd = shift(@ARGV);

$ROUTFILE = shift(@ARGV);
open(STDERR,">>$ROUTFILE");
open(STDOUT,">>$ROUTFILE");

$jetzt = scalar localtime;
print "Starting run_bulk_report_scripts $jetzt...\n";
print "Using datasource: $server / $db\n\n";

#
## General setup
#

my $repdir = '/data/build-public-release/' . $db . '/bulk_reports/';
my $metdir = '/data/build-public-release/' . $db . '/metrics/';
my $bpdir =  '/data/build-public-release/' . $db . '/bulk_processing/';
#
## Call scripts...
#

# ## SCRIPT: report_fb_synonym pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_fb_synonym $jetzt\n";
my $outfile = $repdir . 'fb_synonym_' . $db . '.tsv';
my $call = sprintf("perl report_fb_synonym %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_fb_synonym $jetzt\n\n";

## RETIRED: discontinued following Jim Thurmond email of 14-Mar-2014
## Subject: obsolete bulk data file
## SCRIPT: report_gene_genetic_interactions pg_server db_name pg_username pg_password output_filename.
# $jetzt = scalar localtime;
# print "Calling report_gene_genetic_interactions $jetzt\n";
# my $outfile = $repdir . 'gene_genetic_interactions_' . $db . '.tsv';
# my $call = sprintf("perl report_gene_genetic_interactions %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
# print "\t.$call.\n";
# system($call);
# $jetzt = scalar localtime;
# print "Finished report_gene_genetic_interactions $jetzt\n\n";

## RETIRED: replaced by "mitab.py" in "harvdev-reports" repo. See Jira DB-559.
## SCRIPT: report_physical_interactions pg_server db_name pg_username pg_password output_filename.
# $jetzt = scalar localtime;
# print "Calling report_physical_interactions $jetzt\n";
# my $outfile = $repdir . 'physical_interactions_' . $db . '.tsv';
# my $call = sprintf("perl report_physical_interactions %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
# print "\t.$call.\n";
# system($call);
# $jetzt = scalar localtime;
# print "Finished report_physical_interactions $jetzt\n\n";

## SCRIPT: report_fbrf_w_pmid pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_fbrf_w_pmid $jetzt\n";
my $outfile = $repdir . 'fbrf_pmid_pmcid_doi_' . $db . '.tsv';
my $call = sprintf("perl report_fbrf_w_pmid %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_report_fbrf_w_pmid $jetzt\n\n";


## SCRIPT: report_fbgn_major_accessions pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_fbgn_major_accessions $jetzt\n";
my $outfile = $repdir . 'fbgn_NAseq_Uniprot_' . $db . '.tsv';
my $call = sprintf("perl report_fbgn_major_accessions %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_fbgn_major_accessions $jetzt\n\n";


## SCRIPT: report_fbgn_major_accessions_uniprot pg_server db_name pg_username pg_password output_filename.
# Note - output file omits release number so we can provide an invariant URL to UniProt for file access.
$jetzt = scalar localtime;
print "Calling report_fbgn_major_accessions_uniprot $jetzt\n";
my $outfile = $repdir . 'fbgn_uniprot_' . $db . '.tsv';
my $call = sprintf("perl report_fbgn_major_accessions_uniprot %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_fbgn_major_accessions_uniprot $jetzt\n\n";


## SCRIPT: report_fbgn_annotationID pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_fbgn_annotationID $jetzt\n";
my $outfile = $repdir . 'fbgn_annotation_ID_' . $db . '.tsv';
my $call = sprintf("perl report_fbgn_annotationID %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_fbgn_annotationID $jetzt\n\n";


## SCRIPT: report_gene_mapping pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_gene_mapping $jetzt\n";
my $outfile = $repdir . 'gene_map_table_' . $db . '.tsv';
my $call = sprintf("perl report_gene_mapping %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_gene_mapping $jetzt\n\n";


## SCRIPT: report_insertion_mapping pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_insertion_mapping $jetzt\n";
my $outfile = $repdir . 'insertion_mapping_' . $db . '.tsv';
my $call = sprintf("perl report_insertion_mapping %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_insertion_mapping $jetzt\n\n";


## SCRIPT: report_cdna_clone_data pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_cdna_clone_data $jetzt\n";
my $outfile = $repdir . 'cDNA_clone_data_' . $db . '.tsv';
my $call = sprintf("perl report_cdna_clone_data %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_cdna_clone_data $jetzt\n\n";


## SCRIPT: report_genomic_clone_data pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_genomic_clone_data $jetzt\n";
my $outfile = $repdir . 'genomic_clone_data_' . $db . '.tsv';
my $call = sprintf("perl report_genomic_clone_data %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_genomic_clone_data $jetzt\n\n";

## SCRIPT: report_allele_genetic_interactions pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_allele_genetic_interactions $jetzt\n";
my $outfile = $repdir . 'allele_genetic_interactions_' . $db . '.tsv';
my $call = sprintf("perl report_allele_genetic_interactions %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_allele_genetic_interactions $jetzt\n\n";

## SCRIPT: GA_file_builder pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling GA_file_builder $jetzt\n";
my $outfile = $repdir . 'gene_association.fb';
my $call = sprintf("perl GA_file_builder %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished GA_file_builder $jetzt\n\n";

## SCRIPT: GPI_UP_file_builder pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling GPI_UP_file_builder $jetzt\n";
my $outfile = $repdir . 'gp_information.fb';
my $call = sprintf("perl GPI_UP_file_builder %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished GPI_UP_file_builder $jetzt\n\n";


## SCRIPT: report_allele_phenotype pg_server db_name pg_username pg_password output_filename.
$jetzt = scalar localtime;
print "Calling report_allele_phenotype $jetzt\n";
my $outfile = $repdir . 'allele_phenotypic_data_' . $db . '.tsv';
my $call = sprintf("perl report_allele_phenotype %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_allele_phenotype $jetzt\n\n";

## SCRIPT: report_fbgn_gleanr pg_server db_name pg_username pg_password output_filename
$jetzt = scalar localtime;
print "Calling report_fbgn_gleanr $jetzt\n";
my $outfile = $repdir . 'fbgn_gleanr_' . $db . '.tsv';
my $call = sprintf("perl report_fbgn_gleanr %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_fbgn_gleanr $jetzt\n\n";

## RETIRED: replaced by "report_orthodb_orthologs.py" in "harvdev-reports" repo. See Jira DB-530.
## SCRIPT: report_gene_orthologs pg_server db_name pg_username pg_password output_filename
# $jetzt = scalar localtime;
# print "Calling report_gene_orthologs $jetzt\n";
# my $outfile = $repdir . 'dmel_orthologs_in_drosophila_species_' . $db . '.tsv';
# my $call = sprintf("perl report_gene_orthologs %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
# print "\t.$call.\n";
# system($call);
# $jetzt = scalar localtime;
# print "Finished report_gene_orthologs $jetzt\n\n";


## SCRIPT: report_pmid_fbgn_uniprot pg_server db_name pg_username pg_password output_filename
$jetzt = scalar localtime;
print "Calling report_pmid_fbgn_uniprot $jetzt\n";
my $outfile = $repdir . 'pmid_fbgn_uniprot_' . $db . '.tsv';
my $call = sprintf("perl report_pmid_fbgn_uniprot %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_pmid_fbgn_uniprot $jetzt\n\n";

## RETIRED: replaced by "report_current_gene_product_ids.py" in "harvdev-reports" repo. See Jira DB-535.
## SCRIPT: report_current_gene_product_ids pg_server db_name pg_username pg_password output_filename
# $jetzt = scalar localtime;
# print "Calling report_current_gene_product_ids $jetzt\n";
# my $outfile = $repdir . 'fbgn_fbtr_fbpp_' . $db . '.tsv';
# my $call = sprintf("perl report_current_gene_product_ids %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
# print "\t.$call.\n";
# system($call);
# $jetzt = scalar localtime;
# print "Finished report_current_gene_product_ids $jetzt\n\n";


## SCRIPT: report_library_collection pg_server db_name pg_username pg_password output_filename
$jetzt = scalar localtime;
print "Calling report_library_collection $jetzt\n";
my $outfile = $repdir . 'dataset_metadata_' . $db . '.tsv';
my $call = sprintf("perl report_library_collection %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_library_collection $jetzt\n\n";

## RETIRED
## SCRIPT: report_harvcur_flags pg_server db_name pg_username pg_password output_filename
#$jetzt = scalar localtime;
#print "Calling report_harvcur_flags $jetzt\n";
#my $config = '/users/ftpuser/CUR/curation_guidelines/configs/report_harvcur_flags.conf';
#my $call = sprintf("perl report_harvcur_flags %s %s %s %s %s",$server,$db,$user,$pwd,$config);
#print "\t.$call.\n";
#system($call);
#$jetzt = scalar localtime;
#print "Finished report_harvcur_flags $jetzt\n\n";


## SCRIPT: report_unique_proteins organism_abbreviation pg_server db_name pg_username pg_password output_filename
$jetzt = scalar localtime;
print "Calling report_unique_proteins $jetzt\n";
my $outfile = $repdir . 'dmel_unique_protein_isoforms_' . $db . '.tsv';
my $call = sprintf("perl report_unique_proteins %s %s %s %s %s %s",'Dmel',$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_unique_proteins $jetzt\n\n";

## SCRIPT: gogenes_pubmed_nc.pl non-protein coding genes with GO annotation & pub_med status
$jetzt = scalar localtime;
print "Calling  gogenes_pubmed_nc.pl $jetzt\n";
my $outfile = $bpdir . 'ggnc_' . $db . '.out';
my $call = sprintf("perl gogenes_pubmed_nc.pl Dmel %s %s %s %s %s >& %s",$server,$db,$user,$pwd,$odir,$outfile);
print "\tcalling: $call\n";
system($call);
$jetzt = scalar localtime;
print "Finished gogenes_pubmed_nc.pl\n\n";

## SCRIPT: report_gp2protein pg_server db_name pg_username pg_password logging_output_filename report_output_filename (Dave)
# Note that currently $outfile is not produced.   In script report_gp2protein is set to /dev/null
$jetzt = scalar localtime;
print "Calling report_gp2protein $jetzt\n";
my $outfile = $repdir . 'rgp2p_' . $db . '.out';
my $repfile = $repdir . 'gp2protein.fb';
my $call = sprintf("perl report_gp2protein %s %s %s %s %s %s",$server,$db,$user,$pwd,$outfile,$repfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_gp2protein $jetzt\n\n";

## SCRIPT: package_ontology_reports db_name (Dave)
$jetzt = scalar localtime;
print "Calling package_ontology_reports $jetzt\n";
my $call = sprintf("perl package_ontology_reports %s ontology_reports.tar",$db);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished package_ontology_reports $jetzt\n";

## SCRIPT: report_gene_assessment (Dave)
##USAGE: report_gene_assessment pg_server db_name pg_username pg_password organism_abbreviation output_filename
$jetzt = scalar localtime;
print "Calling report_gene_assessment $jetzt\n";
my $outfile = $metdir . 'gene_assessment_' . $db . '.tsv';
my $call = sprintf("perl report_gene_assessment %s %s %s %s Dmel %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_gene_assessment $jetzt\n";

## SCRIPT: report_gb_linkouts (Dave)
##USAGE: report_gb_linkouts pg_server db_name pg_username pg_password debug_output_filename
$jetzt = scalar localtime;
print "Calling report_gb_linkouts $jetzt\n";
my $outfile = $repdir . 'report_gb_linkouts' . $db . '.out';
my $call = sprintf("perl report_gb_linkouts %s %s %s %s %s %s",$server,$db,$user,$pwd,$repdir,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_gb_linkouts $jetzt\n";

## SCRIPT: flatten_transcriptome (Dave)
##USAGE: flatten_transcriptome pg_server db_name pg_username pg_password organism_abbreviation log_output_filename gff_output_filename
##NOTE: This script produces output used in RPKM report (below)
$jetzt = scalar localtime;
print "Calling flatten_transcriptome $jetzt\n";
my $outfile = $repdir . 'flat_transcriptome_' . $db . '.out';
my $gffoutfile = $repdir . 'flat_transcriptome_' . $db . '.gff';
my $call = sprintf("perl flatten_transcriptome %s %s %s %s Dmel %s %s",$server,$db,$user,$pwd,$outfile,$gffoutfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished flatten_transcriptome $jetzt\n";


## SCRIPT: flatten_transcriptome_unstranded (Dave)
##USAGE: flatten_transcriptome_unstranded pg_server db_name pg_username pg_password organism_abbreviation log_output_filename gff_output_filename
##NOTE: This script produces output used in RPKM report (below)
$jetzt = scalar localtime;
print "Calling flatten_transcriptome_unstranded $jetzt\n";
my $outfile = $repdir . 'flat_transcriptome_unstranded_' . $db . '.out';
my $gffoutfile_u = $repdir . 'flat_transcriptome_unstranded_' . $db . '.gff';
my $call = sprintf("perl flatten_transcriptome_unstranded %s %s %s %s Dmel %s %s",$server,$db,$user,$pwd,$outfile,$gffoutfile_u);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished flatten_transcriptome $jetzt\n";


## SCRIPT: gene_rpkm_report (Dave)
##USAGE: gene_rpkm_report pg_server db_name pg_username pg_password organism_abbrivation release_designator input_unstranded_transcriptome_gff input_stranded_transcriptome_gff output_filename
## NOTE THAT THIS SCRIPT USES FILES SPECIFIED ABOVE (input_unstranded_transcriptome_gff & input_stranded_transcriptome_gff)
$jetzt = scalar localtime;
print "Calling gene_rpkm_report $jetzt\n";
my $outfile = $repdir . 'gene_rpkm_report_' . $db . '.tsv';
my $call = sprintf("perl gene_rpkm_report %s %s %s %s Dmel %s %s %s %s",$server,$db,$user,$pwd,$relno,$gffoutfile_u,$gffoutfile,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished gene_rpkm_report $jetzt\n";

# RETIRED: replaced by "report_disease_model_data.py" in "harvdev-reports" repo. See Jira DB-560.
# ## SCRIPT: report_allele_doid (Dave)
# ## USAGE: report_allele_doid pg_server db_name pg_username pg_password output_filename
# $jetzt = scalar localtime;
# print "Calling report_allele_doid $jetzt\n";
# my $outfile = $repdir . 'allele_human_disease_model_data_' . $db . '.tsv';
# my $call = sprintf("perl report_allele_doid %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
# print "\t.$call.\n";
# system($call);
# $jetzt = scalar localtime;
# print "Finished report_allele_doid $jetzt\n";

## SCRIPT: report_gene_group_data (Dave)
## Reports FB Gene Groups, their parent groups (if any), and asociated genes (if any)
## USAGE: report_gene_group_data pg_server db_name, pg_username pg_password report_output_filename
$jetzt = scalar localtime;
print "Calling report_gene_group_data $jetzt\n";
my $outfile = $repdir . 'gene_group_data_' . $db . '.tsv';
my $call = sprintf("perl report_gene_group_data %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_gene_group_data $jetzt\n";

## SCRIPT: report_gene_group_data_hgnc (Dave)
## Reports FB Gene Group -> HGNC Family ID correspondence
## USAGE: report_gene_group_data_hgnc pg_server db_name, pg_username pg_password report_output_filename
$jetzt = scalar localtime;
print "Calling report_gene_group_data_hgnc $jetzt\n";
my $outfile = $repdir . 'gene_groups_HGNC_' . $db . '.tsv';
my $call = sprintf("perl report_gene_group_data_hgnc %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_gene_group_data_hgnc $jetzt\n";

## SCRIPT: report_gene_snapshots (Dave)
## Reports gene snapshots curated for localized Dmel genes 
## USAGE: report_gene_snapshots pg_server db_name pg_username pg_password report_output_filename
$jetzt = scalar localtime;
print "Calling report_gene_snapshots $jetzt\n";
my $outfile = $repdir . 'gene_snapshots_' . $db . '.tsv';
my $call = sprintf("perl report_gene_snapshots %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_gene_snapshots $jetzt\n";

## RETIRED: replaced by "report_rnacentral_json.py" in "harvdev-reports" repo. See Jira DB-533.
## SCRIPT: report_rnacentral (Dave)
## Reports ncRNA genes for ENA / RNAcentral
##  USAGE: report_rnacentral pg_server db_name pg_username pg_password output_filename
# $jetzt = scalar localtime;
# print "Calling report_rnacentral $jetzt\n";
# my $outfile = $repdir . 'ncRNA_genes_' . $db . '.tsv';
# my $call = sprintf("perl report_rnacentral %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
# print "\t.$call.\n";
# system($call);
# $jetzt = scalar localtime;
# print "Finished report_rnacentral $jetzt\n";


## SCRIPT: report_human_orthologs (Dave)
## Reports human disease genes orthologous to fly genes
##   USAGE: report_human_orthologs pg_server db_name pg_username pg_password report_output_filename
$jetzt = scalar localtime;
print "Calling report_human_orthologs $jetzt\n";
my $outfile = $repdir . 'dmel_human_orthologs_disease_' . $db . '.tsv';
my $call = sprintf("perl report_human_orthologs %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_human_orthologs $jetzt\n";

## SCRIPT: report_funcomps (Dave)
## Reports fly genes and the ortholog genes that functionally complement them
## USAGE: report_funcomps pg_server db_name pg_username pg_password output_filename
$jetzt = scalar localtime;
print "Calling report_funcomps $jetzt\n";
my $outfile = $repdir . 'gene_functional_complementation_' . $db . '.tsv';
my $call = sprintf("perl report_funcomps %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_funcomps $jetzt\n";


## REPORT IS NOT READY YET -- SEE JIRA DB-239    DE 25-Aug-2015
# ## SCRIPT: report_gene_expression_annotations (Dave)
# ## USAGE: report_gene_expression_annotations pg_server db_name pg_username pg_password organism_abbreviation logging_output_filename report_output_filename  (Note: log output is set to null in the script at the mo)
# $jetzt = scalar localtime;
# print "Calling report_gene_expression_annotations $jetzt\n";
# my $outfile = $repdir . 'gene_expression_data_' . $db . '.tsv';
# my $logfile = $repdir . 'gene_expression_data_' . $db . '.out';
# my $call = sprintf("perl report_gene_expression_annotations %s %s %s %s Dmel %s %s",$server,$db,$user,$pwd,$logfile,$outfile);
# print "\t.$call.\n";
# system($call);
# $jetzt = scalar localtime;
# print "Finished report_gene_expression_annotations $jetzt\n";

## SCRIPT: report_fbab_fbgn_rels
## Reports experimental aberration gene deletion/duplication information
## USAGE: report_fbab_fbgn_rels pg_server db_name pg_username pg_password output_filename
$jetzt = scalar localtime;
print "Calling report_fbab_fbgn_rels.pl $jetzt\n";
my $outfile = $repdir . 'aberration_experimental_gene_del_dup_data.' . $db . '.tsv';
my $call = sprintf("perl report_fbab_fbgn_rels.pl %s %s %s %s %s",$server,$db,$user,$pwd,$outfile);
print "\t.$call.\n";
system($call);
$jetzt = scalar localtime;
print "Finished report_fbab_fbgn_rels.pl $jetzt\n";


$jetzt = scalar localtime;
print "Finished run_bulk_report_scripts $jetzt...\n";

