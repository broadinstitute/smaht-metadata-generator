#!/bin/bash

set -eu

if [ "$#" -ne 1 ]; then echo "need to provide a donor ID, in the format of 0XX"; exit 1; fi

donor_id=$1

bash 1.generate_metadata.sh "${donor_id}"
bash 2.validate_metadata.sh "${donor_id}"/output/all_metadata.xlsx
