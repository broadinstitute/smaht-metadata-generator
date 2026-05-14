#!/bin/bash

set -eu

if [[ ! -f dac ]]; then echo "I'm assuming the DAC data validation and submission tool available here as 'dac', but it seems inaccessible for me."; exit 1; fi
if [ "$#" -ne 1 ]; then echo "Need to provide either donor ID, in the format of 0XX, or an explicit path to an existing metadata workbook"; exit 1; fi

input=$1
if [[ ${input} =~ ^[0-9]+$ ]];
then 
    donor_id=${input}
    METADATA_WORKBOOK="${donor_id}/output/all_metadata.xlsx"
else
    METADATA_WORKBOOK="${input}"
    donor_id=$(echo "${METADATA_WORKBOOK}" | grep -Eo "SMHT[0-9]{3}" | sed 's/SMHT//')
    if [[ -z "${donor_id}" ]]; then donor_id=$(echo "${METADATA_WORKBOOK}" | awk -F '/' '{print $1}'); fi
fi
if [[ ! -f "${METADATA_WORKBOOK}" ]]; then echo "${METADATA_WORKBOOK} does not exist" ; exit 1; fi
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
