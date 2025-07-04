# harvdev-bulk-reports

<!-- toc -->

- [Overview](#Overview)
- [BulkReportsSOP](#BulkReportsSOP)
  * [EnvironmentVariables](#EnvironmentVariables)
  * [PipelineSummary](#PipelineSummary)
  * [NextSteps](#NextSteps)
- [Gal4DriverImages](#Gal4DriverImages)
- [TroubleShooting](#TroubleShooting)

## Overview
This repo contains scripts that are used to generate most FlyBase bulk reports, as well as reports for FlyBase curators and external collaborators. This is a new public repository that replaces the private, retired `harvdev-reports` (renamed to `harvdev-reports-old`) repository. This repo contains a mix of newer python scripts and older perl (taken from the `fb_cvs/FB/scripts/reports` repo). Bulk reports for FlyBase users are shipped off to IUDev and incorporated into the public FB site on the [Downloads page](http://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release) and [FTP site](ftp://ftp.flybase.net/releases/current). Other reports are posted to internal and external FTP sites for use.  
This code is intended to run in docker using the `Bulk_Reports` GoCD pipeline (`Reporting_Build` pipeline group).  
This `Bulk_Reports` GoCD pipeline is automatically triggered by the completion of the related, upstream `Reporting_Build` GoCD pipeline.  

## BulkReportsSOP

### EnvironmentVariables
Variables are set at the outset of the entire reporting release build for the `Reporting_Build` pipeline group environment that runs all pipelines.   
The variables relevant to the generation of bulk reports are as follows:  

**Here is the list of plain text variables that change each release:**  
`ANNOTATIONRELEASE` - the Dmel annotation release version for the reporting build - check the FB public release notes and increment that by one.  
`RELEASE` - the current release under construction: e.g., `2024_01`.  
`PREV_RELEASE` - the release number for the previous reporting build: e.g., if `RELEASE=2024_01`, then `PREV_RELEASE=2023_06`.  

**Here is the list of plain text variables that rarely change from release to release:**  
`BUILD_DIR` - the directory in which files for the db build are stored: i.e., `/data/build-public-release`.  
`USER` - the postgres use name, which does not change: i.e., `go`.  
`DOCKER_USER` - the docker username for the FlyBase docker account: i.e., `harvdevgocd`.  
`REPORTING_SERVER` - the server on which the reporting dbs are kept: i.e., `flysql25`.  
`REPORTING_DATABASE` - the name of the reporting build that is in progress: e.g., `fb_2024_01_reporting`.  
`ASSEMBLY` - the name of the current Dmel reference genome assembly as it is known at the Alliance: i.e., `R6`.  

### PipelineSummary
Download files are generated by the [Bulk_Reports](http://flysql22:8153/go/admin/pipelines/Bulk_Reports/general) in `Reporting_Build` pipeline group.  
The pipeline automates these steps:  
1. Gets HarvDev docker container and builds the appropriate docker image using this repo.  
2. Saves output bulk files to the `/data/build-reporting/fb_${RELEASE}_reporting/bulk_reports` directory.  
3. Checks file sizes relative to a reference release (usually the previous release) to detect any missing files, or files that are <99% expected size.  
- For FB2024_01 and earlier, file sizes were manually checked.
- See the [Reporting Builds](https://drive.google.com/drive/folders/1lHjCrX-ee7pSaThbo4UuMJ3LWGjngKja) Google Drive directory for examples.  
4. Notifies HarvDev by email that the files have been generated.  

### DetailedSOP
1. Review the `Environment variables` for the `Reporting_Build` GoCD pipeline group.  
- They should've been updated in earlier release build steps, but double check them.  
2. Run the pipeline.  
3. When the files have been generated, upon receiving the email, check files sizes.  

### NextSteps
1. If both this `Bulk_Reports` GoCD pipeline, and the upstream `Reporting_Build` GoCD pipeline, seem to have completed without issue and all file sizes are normal, then manually start the `Upload_Reporting_Build` GoCD pipeline, which will upload all files related to the release build for various users.  

## Gal4DriverImages
Occasionally (up to once per FlyBase public release), Sian provides a flash drive with updated images for the Gal4 table. These new images need to be added to the existing set (rarely, there will be replacements). The updated image set needs to be uploaded to the FB FTP for IUDev as a tarball. Also, the metadata needs to be updated so that image names are incorporated into the JSON file that feeds the table.

### Images
1. Make a new directory for the images (where YYYY-MM-DD in the command below represents today's date).
cd /data/harvcur/gal4images/
mkdir -m "FUG4_YYYY-MM-DD"
2. Move all images from the most recent directory to the new one.
mv FUG4_YYYY-MM-DD1/** FUG4_YYYY-MM-DD2
3. Copy new images from the flash drive to the new directory (GUI is easiest).
4. Make a new tarball.
cd /data/harvcur/gal4images/
tar -zcvf FUG4_YYYY-MM-DD.tar.gz FUG4_YYYY-MM-DD
5. Check the tarball (want to see that the directory holding the images is included).
tar -tzf FUG4_YYYY-MM-DD
6. Copy the tarball to the internal FlyBase FTP site (foriu/gal4images directory).

## Metadata
1. Make a new branch.
2. Copy image metadata for new images from this google sheet ("TEMPLATE" tab):
https://docs.google.com/spreadsheets/d/1JRdVJezUPZPxlPJhHSE0bJ73aA5Yw8FlBR5ASkNmq6M/edit?pli=1&gid=1067073457#gid=1067073457
3. Add the new lines to the `src/image_metadata.txt` file in this repo, save, and create a pull request.


## TroubleShooting
The [Reporting Build SOP](https://github.com/FlyBase/harvdev-docs/blob/master/reporting_build/reporting_build_sop.md#TroubleShooting) discusses various troubleshooting scenarios for dealing with failed scripts and GoCD pipelines.
