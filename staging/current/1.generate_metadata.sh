#!/bin/bash
set -eu

if [ "$#" -ne 1 ]; then echo "need to provide a donor ID, in the format of 0XX"; exit 1; fi
#########################################################################
#                   CHANGE THIS INPUT FOR EACH DONOR                    #
#########################################################################
donor_id=$1

#########################################################################
# if you follow a naming convention, no need to change the following:
srDNA_TERRA_TSV="${donor_id}/${donor_id}.srDNA.terra.tsv"
srRNA_TERRA_TSV="${donor_id}/${donor_id}.srRNA.terra.tsv"
lrDNA_TERRA_TSV="${donor_id}/${donor_id}.lrDNA.terra.tsv"
srDNA_TRIMMED_TSV="${donor_id}/SMHT${donor_id}_Short_Read_DNA_trimmed.tsv"
srRNA_TRIMMED_TSV="${donor_id}/SMHT${donor_id}_Short_Read_RNA_trimmed.tsv"
lrDNA_TRIMMED_TSV="${donor_id}/SMHT${donor_id}_Long_Read_DNA_trimmed.tsv"
DONOR_INFORMATION_TSV="${donor_id}/SMHT${donor_id}.donor_information.tsv"

echo "=================================================="
echo "         Working on donor SMHT${donor_id}         "
echo "=================================================="

#########################################################################
# reformat the Terra=generated TSV
echo "=================================================="
echo "Reformatting Terra-derived TSV to generator......." 
echo "=================================================="
# allows for different modalities to be ready at different times
par1=""
arg1=""
if [[ -f "${srDNA_TERRA_TSV}" ]]; then
    python3 ../../utils/0A_strip_cram_and_bam_paths.py \
        --input  "${srDNA_TERRA_TSV}" \
        --output "${srDNA_TRIMMED_TSV}"
    par1="--sr-dna"
    arg1="${srDNA_TRIMMED_TSV}"
fi 
par2=""
arg2=""
if [[ -f "${srRNA_TERRA_TSV}" ]]; then
    python3 ../../utils/0A_strip_cram_and_bam_paths.py \
        --input  "${srRNA_TERRA_TSV}" \
        --output "${srRNA_TRIMMED_TSV}"
    par2="--rna"
    arg2="${srRNA_TRIMMED_TSV}"
fi
par3=""
arg3=""
if [[ -f "${lrDNA_TERRA_TSV}" ]]; then
    python3 ../../utils/0A_strip_cram_and_bam_paths.py \
        --input  "${lrDNA_TERRA_TSV}" \
        --output "${lrDNA_TRIMMED_TSV}"
    par3="--lr-dna"
    arg3="${lrDNA_TRIMMED_TSV}"
fi

modalities=()
[ -n "${par1}" ] && { modalities+=("$par1") ; modalities+=("${arg1}"); }
[ -n "${par2}" ] && { modalities+=("$par2") ; modalities+=("${arg2}"); }
[ -n "${par3}" ] && { modalities+=("$par3") ; modalities+=("${arg3}"); }

#########################################################################
# usually, no need the change the following
mkdir -p "${donor_id}/output"

echo "=================================================="
echo "Generating metadata sheet with in-house scripts..."
echo "=================================================="
python3 \
../../scripts/generate_metadata.py \
     "${modalities[@]}" \
     --donor-info  "${DONOR_INFORMATION_TSV}" \
     --uberon ../../input_examples/tissue_uberon_identifiers.tsv \
     --submitter-prefix BROAD \
     --out-donor-tsv                      "${donor_id}/output/"donor.tsv \
     --out-donor-xlsx                     "${donor_id}/output/"donor.xlsx \
     --out-tissue-tsv                     "${donor_id}/output/"tissue.tsv \
     --out-tissue-xlsx                    "${donor_id}/output/"tissue.xlsx \
     --out-tissuesample-tsv               "${donor_id}/output/"tissuesample.tsv \
     --out-tissuesample-xlsx              "${donor_id}/output/"tissuesample.xlsx \
     --out-analyte-tsv                    "${donor_id}/output/"analyte.tsv \
     --out-analyte-xlsx                   "${donor_id}/output/"analyte.xlsx \
     --out-library-tsv                    "${donor_id}/output/"library.tsv \
     --out-library-xlsx                   "${donor_id}/output/"library.xlsx \
     --out-librarypreparation-tsv         "${donor_id}/output/"librarypreparation.tsv \
     --out-librarypreparation-xlsx        "${donor_id}/output/"librarypreparation.xlsx \
     --out-sequencing-tsv                 "${donor_id}/output/"sequencing.tsv \
     --out-sequencing-xlsx                "${donor_id}/output/"sequencing.xlsx \
     --out-fileset-tsv                    "${donor_id}/output/"fileset.tsv \
     --out-fileset-xlsx                   "${donor_id}/output/"fileset.xlsx \
     --out-unalignedreads-tsv             "${donor_id}/output/"unalignedreads.tsv \
     --out-unalignedreads-xlsx            "${donor_id}/output/"unalignedreads.xlsx \
     --out-alignedreads-tsv               "${donor_id}/output/"alignedreads.tsv \
     --out-alignedreads-xlsx              "${donor_id}/output/"alignedreads.xlsx \
     --out-software-tsv                   "${donor_id}/output/"software.tsv \
     --out-software-xlsx                  "${donor_id}/output/"software.xlsx \
     --out-combined-xlsx       "${donor_id}/output/"all_metadata.xlsx
echo "!!!DONE!!!"
