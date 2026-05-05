#!/bin/bash

set -eu

# script_name=$(basename "$0")
# donor_id=$(echo "${script_name}" | grep -Eo "SMHT[0-9]{3}")
#
if [ $# -ne 1 ]; then echo "Please provide 3-digit donor_id" && exit 1; fi
donor_id=$1
metadata_g_sheet_downloaded_validated=$(find . -type f -name "Broad_*${donor_id}*.xlsx")

while true; do
    read -p "Do you want to update the metadata submitted? (y/n): " answer
    case $answer in
        [Yy]* )
            echo "User opted to update submitted metadata..."
            flag_opted='--submit-update'
            break
            ;;
        [Nn]* )
            echo "User opted to not update submitted metadata..."
            flag_opted='--submit'
            break
            ;;
        * )
            echo "Please answer y or n."
            ;;
    esac
done
echo "${metadata_g_sheet_downloaded_validated}"
echo "${flag_opted}"

./dac \
submit-metadata \
    "${metadata_g_sheet_downloaded_validated}" \
    --env data \
    "${flag_opted}"
