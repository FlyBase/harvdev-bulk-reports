# harvdev-reports

<!-- toc -->

- [Overview](#Overview)
- [GoCDPipeline](#GoCDPipeline)

## Overview
This is a new public repository that replaces the private, retired `harvdev-reports` (renamed to `harvdev-reports-old`) repository.

FlyBase in-house scripts that repackage data from reporting chado db into FlyBase bulk reports (flat files).  
These bulk reports are shipped off to IUDev and incorporated into the public FB site.  
This repo contains a mix of newer python scripts and older perl (taken from the `fb_cvs/FB/scripts/reports` repo).     
Files available to the public on the "Downloads > [Current Release](http://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release) web page, 
as well as the [FTP site](ftp://ftp.flybase.net/releases/current).  
Modifications should be made in new branches, then merged back into the "master" branch.  
The master branch is used to spawn "releases", which are used to generate the data for export.  
A GoCD pipeline runs these scripts and sends output off to IUDev. See [GoCDPipeline](#GoCDPipeline) section.  
Some scripts (not all yet) have a "-l" option so that script can be run on local machine for development.  

## GoCDPipeline
Most of this process is now automated in the a GoCD pipeline.  
See [Bulk_Reports](http://flysql22:8153/go/admin/pipelines/Bulk_Reports/general) in "Reporting Build" pipeline group.  
See also general [GoCD documentation](https://github.com/FlyBase/harvdev-docs/blob/master/gocd/gocd.md).  

The pipeline automates these steps:  
1. Gets HarvDev docker container.  
2. Builds docker image with FB `alliance-flybase` and Alliance `agr_schemas` repos.  
3. Fetches data using `harvdev-reports` scripts from a specified branch or release tag.  
4. Saves output bulk files to the `/data/build-reporting/fb_<release>_reporting/bulk_reports` folder.  
5. Packages up files into various tarballs and uploads them to the appropriate ftp sites.  
6. Sends e-mails to relevant personnel upon successful tarball upload.  

Some manual steps are still required:  
1. Generation of an appropriate `harvdev-reports` release tag.  
2. Updating `Environment variables` for the GoCD pipeline.  
  - SERVER - flysql machine where the reporting db is located: e.g., `flysql20`  
  - DATABASE - the name of the reporting db to use: e.g., `fb_2019_03_reporting`  
  - USER - for GoCD, use `go`  
  - ANNOTATIONRELEASE - the Dmel annotation release for the db: e.g., `R6.28`  
  - RELEASE - YYYY_NN: e.g., `2019_03`  
  - TAG - the release tag for `harvdev-reports`` Git repo that is to be used: e.g., `2019_04`.  
3. Always good to check file sizes generated against previous releases.  

