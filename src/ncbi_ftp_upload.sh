#!/bin/bash

# This shell script uploads three files to the NCBI FTP site.
# Two arguments are required: 1) the FB release number: e.g., 2024_01; 2) the NCBI FTP password.

# FTP details for ftp-private.ncbi.nih.gov
HOST="130.14.250.6"
USER="flybase"
cd /data/build-public-release/fb_$1_reporting/bulk_reports
pwd

# FTP login and upload.
ftp -inv $HOST <<EOF
user $USER $2
cd holdings
put nt-flybase.xml
put pm-flybase.xml
put pr-flybase.xml
bye
EOF


