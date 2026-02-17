#!/bin/bash
#
# copy_scrna_metadata.sh
#
# Copies scRNAseq metadata files forward from the previous release,
# renaming them for the current release. These files don't change
# between releases (until a new scRNAseq curator is hired at Cambridge).
#
# Expects GoCD environment variables:
#   RELEASE          - current release (e.g., 2025_06)
#   PREV_RELEASE     - previous release (e.g., 2025_05)
#   BUILD_DIR        - build directory (e.g., /data/build-public-release)
#
# Usage:
#   copy_scrna_metadata.sh
#

set -euo pipefail

# Validate required environment variables
for var in RELEASE PREV_RELEASE BUILD_DIR; do
    if [ -z "${!var:-}" ]; then
        echo "ERROR: Required environment variable $var is not set." >&2
        exit 1
    fi
done

PREV_DIR="${BUILD_DIR}/fb_${PREV_RELEASE}_reporting/bulk_reports"
CURR_DIR="${BUILD_DIR}/fb_${RELEASE}_reporting/bulk_reports"

FILES=(
    "scRNAseq_metadata_fb_YYYY_MM_reporting.json"
    "scRNAseq_metadata_schema_fb_YYYY_MM_reporting.yaml"
)

echo "Copying scRNAseq metadata from ${PREV_RELEASE} to ${RELEASE}..."
echo "  Source: ${PREV_DIR}"
echo "  Dest:   ${CURR_DIR}"

# Verify source directory exists
if [ ! -d "${PREV_DIR}" ]; then
    echo "ERROR: Previous release directory not found: ${PREV_DIR}" >&2
    exit 1
fi

# Verify destination directory exists
if [ ! -d "${CURR_DIR}" ]; then
    echo "ERROR: Current release directory not found: ${CURR_DIR}" >&2
    exit 1
fi

for template in "${FILES[@]}"; do
    src_file="${PREV_DIR}/${template//YYYY_MM/${PREV_RELEASE}}"
    dst_file="${CURR_DIR}/${template//YYYY_MM/${RELEASE}}"

    if [ ! -f "${src_file}" ]; then
        echo "ERROR: Source file not found: ${src_file}" >&2
        exit 1
    fi

    cp "${src_file}" "${dst_file}"
    echo "  OK: $(basename "${dst_file}") ($(stat -c%s "${dst_file}") bytes)"
done

echo "Done. scRNAseq metadata files copied successfully."
