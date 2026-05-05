#!/bin/bash
set -eu


if [ "$#" -ne 1 ]; then
    echo "need to provide a donor ID, in the format of 0XX"; exit 1;
else
    donor_id=$1
fi


mkdir -p "${donor_id}/output"
cp templates/0XX.lrDNA.terra.tsv \
   "${donor_id}/${donor_id}.lrDNA.terra.tsv"
cp templates/0XX.srDNA.terra.tsv \
   "${donor_id}/${donor_id}.srDNA.terra.tsv"
cp templates/0XX.srRNA.terra.tsv \
   "${donor_id}/${donor_id}.srRNA.terra.tsv"
sed "s/SMHTXXX/SMHT${donor_id}/" \
   templates/SMHTXXX.donor_information.template.tsv \
   > "${donor_id}/SMHT${donor_id}.donor_information.tsv"
