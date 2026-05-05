#!/bin/bash

set -eu

if [[ ! -f dac ]]; then echo "I'm assuming the DAC data validation and submission tool available here as 'dac', but it seems inaccessible for me."; exit 1; fi
if [ "$#" -ne 1 ]; then echo "need to provide a donor ID, in the format of 0XX"; exit 1; fi

donor_id=$1
METADATA_WORKBOOK="${donor_id}/output/all_metadata.xlsx"
if [[ ! -f "${METADATA_WORKBOOK}" ]]; then echo "${METADATA_WORKBOOK} does not exist" ; exit 1; fi
# donor_id=$(echo "${METADATA_WORKBOOK}" | grep -Eo "SMHT[0-9]{3}" | sed 's/SMHT//')
# if [[ -z "${donor_id}" ]]; then
#     donor_id=$(echo "${METADATA_WORKBOOK}" | awk -F '/' '{print $1}')
# fi
echo "=================================================="
echo "Validating metadata for ${donor_id} with DAC submitr tool..."
echo "=================================================="
./dac \
submit-metadata-bundle \
    "${METADATA_WORKBOOK}" \
    --env data \
    --validate \
| tee \
    "${donor_id}/output/"metadata_validation_report.log
echo "!!!DONE!!!"
