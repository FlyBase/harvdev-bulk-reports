#!/bin/bash

# This shell script uploads one file to the EBI EPMC FTP site.
# Three arguments are required: 1) the FB release number: e.g., 2024_01; 2) the EBI FTP password; 3) the EBI directory.

# FTP details for labslink.ebi.ac.uk
HOST="193.62.193.162"
USER="elinks"
cd /data/build-public-release/fb_$1_reporting/bulk_reports

# FTP login and upload.
ftp -inv $HOST <<EOF
user $USER $2
cd $3
put epmc-flybase.xml.gz
bye
EOF


