# smaht-metadata-generator
Pipeline to generate standardized metadata sheets for SMaHT sequencing projects

# Metadata Pipeline

This repository provides a **pipeline** for generating standardized metadata sheets required for genomic sequencing submissions, covering Donor, Tissue, TissueSample, Analyte, Library, FileSet, UnalignedReads, and AlignedReads.

---

## Overview

This pipeline reads in: donor metadata, sample tables for short read, long read, and RNA inputs, plus a UBERON mapping file.  It then produces eight metadata sheets with properly formatted `submitted_id` values:

- **Donor**: one row per donor
- **Tissue**: one row per core tissue code
- **TissueSample**: one row per individual sample
- **Analyte**: DNA/RNA preparations
- **Library**: sequencing libraries
- **FileSet**: grouping of libraries
- **UnalignedReads**: raw read files (BAM)
- **AlignedReads**: alignment outputs (CRAM or BAM)

These are output both as TSV and XLSX.

**For CODEC sequencing workflows**, this repository also includes a specialized script (`generate_codec_metadata.py`) that produces additional metadata sheets specific to CODEC protocols, including LibraryPreparation, Sequencing, VariantCalls, and Software sheets.

---

## Optional Preprocessing Utilities

To assist with metadata preparation, we provide an optional helper script in the [`utils/`](./utils) directory. This is **not required**, but it can help standardize and automate input formatting prior to running the metadata generation pipeline.

---

### `0A_strip_cram_and_bam_paths.py`

This script strips full file paths in your sample tables down to just the base file name (e.g., `gs://bucket/data/sample.bam` → `sample.bam`).  

**It works for both regular sequencing and CODEC workflows**, handling columns such as:

**Regular sequencing:**
- `cram_path`
- `input_bam`

**CODEC sequencing:**
- `MolConsensusBAM`
- `RAW_BAM`
- `vcf`

**Input:**  
A tab-delimited `.tsv` file with file path columns.

**Output:**  
A `.tsv` with the same structure, but with only the file names in the path fields.

**Example:**

**Input:**

| sample_id | RAW_BAM                         | MolConsensusBAM                   |
|-----------|-----------------------------------|-----------------------------------|
| SMHT001   | gs://bucket/data/SMHT001.raw.aligned.bam      | gs://bucket/data/SMHT001.mol_consensus.aligned.bam  |

**Output:**

| sample_id | input_bam     | MolConsensusBAM   |
|-----------|---------------|-------------------|
| SMHT001   | SMHT001.raw.aligned.bam  | SMHT001.mol_consensus.aligned.bam   |

**Usage:**
```bash
python utils/0A_strip_cram_and_bam_paths.py --input <input_tsv> --output <output_tsv>
```

---

## Supported Sequencing Technologies

- **WGS_ILLUMINA** (short read)  
- **WGS_PACBIO**  (long read)  
- **RNA_WATCHMAKER** (RNA total)  
- **RNA_TRUSEQ**     (RNA mRNA)  
- **CODEC** (specialized sequencing with one protocol: DDBTP)
-   NOTE: On August 15, 2025, notified by the CODEC group at the Broad that for 
    production samples, only DDBTP protocol is provided (originally CODEC had 
    three protocols: DDBTP, DRV1, DRV2).

The pipeline infers these from file names or analyte IDs.

---
## Tissue Label Normalization

The pipeline automatically normalizes specific brain tissue labels for consistency:
- `TEMPORAL_LOBE` → `TEMPORAL-LOBE`
- `FRONTAL_LOBE` → `FRONTAL-LOBE`
- `L_HIPPOCAMPUS` → `L-HIPPOCAMPUS`
- `R_HIPPOCAMPUS` → `R-HIPPOCAMPUS`

All other tissue names retain their original formatting (spaces converted to hyphens, uppercase).

---

## Input Files

### Donor Info TSV

- Provided via `--donor-info`  
- Must include at least columns:  
  - `donor` (unique donor code, e.g. `SMHT001`)  
  - `gender` (`Male`/`Female`)  
  - `age` (integer)

### Sample TSVs (Short Read DNA/Long Read DNA/RNA)

- **Short-read DNA**: Passed via `--sr-dna` (accepts multiple files)
- **Long-read DNA**: Passed via `--lr-dna` (accepts multiple files)  
- **RNA**: Passed via `--rna` (accepts multiple files)

At least one of these options must be provided. Filenames should reflect their type:
- Short read tables: filename should contain `short_read` or similar
- Long read tables: filename should contain `long_read` or `long`
- RNA tables: filename should contain `watchmaker` or `truseq`
Each must include at least:

| Column                 | Description                                  |
|------------------------|----------------------------------------------|
| `collaborator_sample_id` | `{donor}-{tissue}-{core}` (e.g. `SMHT001-3S-001A1`) |
| `donor_id`               | Matching donor code (e.g. `SMHT001`)           |

Additional columns (e.g. `input_bam`, `cram_path`) are used for unaligned/aligned sheets.

### CODEC Sample TSVs

For **CODEC workflows**, use the specialized `generate_codec_metadata.py` script with:

- **CODEC files**: Passed via `--codec` (accepts multiple files)

CODEC files must include at least:

| Column                 | Description                                  |
|------------------------|----------------------------------------------|
| `sample_id`            | Unique sample identifier                     |
| `collaborator_sample_id` | `{donor}-{tissue}-{core}` (e.g. `SMHT001-3S-001A1`) |
| `donor_id`             | Matching donor code (e.g. `SMHT001`)        |
| `MolConsensusBAM`      | Path to processed duplex BAM file           |
| `RAW_BAM`              | Path to raw alignment BAM file              |
| `vcf`                  | Path to variant calls VCF file              |

### UBERON Mapping TSV

- Provided via `--uberon`  
- Must include:
  - `tissue_identifier_code` (e.g. `3S`)  
  - `corresponding_tissue`   (e.g. `Heart`)

The `input_examples/` folder includes sample TSVs that demonstrate the expected format for:

- Donor info files
- Short-read, long-read, and RNA sample tables
- CODEC sample tables
- UBERON tissue mapping files

These examples are intended to help you test the pipeline or structure your own metadata inputs.
#### Updating Tissue Identifiers
The metadata generation system relies on `tissue_uberon_identifiers.tsv` to map tissue identifier codes to standardized tissue names and UBERON IDs. If you encounter samples with tissue codes that are not included in this mapping file (e.g., code `3Y`), you'll need to update the file with the correct information.

**How to Update Missing Tissue Identifiers**

1. **Identify the missing tissue code** from your sample data
2. **Look up the correct tissue information** using the SMaHT database links below
3. **Update** `tissue_uberon_identifiers.tsv` with the proper format

***Required Information Format***

For each tissue identifier, you need:
- **Tissue identifier code** (e.g., `3Y`)
- **Corresponding tissue name** (standardized format)
- **UBERON ID** (ontology identifier)

#### SMaHT Database Resources

Use these example links to find the correct tissue information:

#### Sample-Specific Searches
- [SMHT022 Donor Samples](https://data.smaht.org/search/?donor.display_title=SMHT022&submission_centers.display_title=NDRI+TPC&type=SampleSource)
- [SMHT005 Donor Samples](https://data.smaht.org/search/?donor.display_title=SMHT005&submission_centers.display_title=NDRI+TPC&type=SampleSource)
**Search Strategy:** If you can't find your tissue of interest from these donors, systematically search other donors by replacing SMHT022 with other donor IDs (e.g., `SMHT001`, `SMHT003`, etc.) until you locate the tissue identifier you need.

#### Ontology Reference
- [UBERON Ontology Terms](https://data.smaht.org/search/?type=OntologyTerm) - Browse all available tissue ontology terms

**Example Update Process**

1. **Find missing code**: Sample has tissue code `3Y` but it's not in `tissue_uberon_identifiers.tsv`
2. **Search SMaHT database**: Use the donor-specific links above to find what tissue `3Y` represents
3. **Add to mapping file**: Update `tissue_uberon_identifiers.tsv` with the new row:
   ```
    tissue_identifier_code	corresponding_tissue	uberon_id
    3Y	L_OVARY	UBERON:0002119
    1D	LUNG	UBERON:0008952
   ```

#### Notes
- Always verify the tissue name formatting matches existing entries in the file
- Ensure UBERON IDs are accurate and follow the `UBERON:xxxxxxx` format
- Test your updated mapping file before running the full metadata generation

---

## Generated Sheets

Each sheet is written as both TSV and XLSX. The CLI flags control their paths.

| Sheet Name      | CLI Flag               | Description                             |
|-----------------|------------------------|-----------------------------------------|
| Donor           | `--out-donor-*`        | One row per donor                       |
| Tissue          | `--out-tissue-*`       | One row per core tissue code            |
| TissueSample    | `--out-tissuesample-*` | One row per sample suffix               |
| Analyte         | `--out-analyte-*`      | DNA/RNA preparations                    |
| Library         | `--out-library-*`      | Sequencing libraries                    |
| FileSet         | `--out-fileset-*`      | Grouped libraries                       |
| UnalignedReads  | `--out-unalignedreads-*` | Raw BAM files                         |
| AlignedReads    | `--out-alignedreads-*`  | Sorted CRAM files                      |
| Combined Excel  | `--out-combined-xlsx`  | Multi sheet workbook of all above      |

### Additional CODEC-Specific Sheets

When using `generate_codec_metadata.py`, additional sheets are generated:

| Sheet Name           | CLI Flag                      | Description                           |
|---------------------|-------------------------------|---------------------------------------|
| LibraryPreparation  | `--out-library-preparation-*` | CODEC protocol details (DDBTP) |
| Sequencing          | `--out-sequencing-*`          | Sequencing parameters                 |
| VariantCalls        | `--out-variantcalls-*`        | VCF files and variant metadata       |
| Software            | `--out-software-*`            | Analysis pipeline information         |

---

## File Naming & Column Requirements

- **Submitted IDs** must follow patterns like `PREFIX_<SHEET>_<donor>_<tissue>_<type>_<suffix>`.  See individual scripts for regex patterns.
- Input TSVs must have the columns listed above; additional columns are optional but required ones cannot be renamed.

---

## Installation

Requires Python 3.8+ and:

```
pip install pandas openpyxl
```

Optionally, create a `requirements.txt`:

```text
pandas>=1.3
openpyxl>=3.0
```

and install via:

```
pip install -r requirements.txt
```

---

## Usage

### Regular Sequencing (Short Read DNA, Long Read DNA, Short Read RNA)

```bash
python generate_metadata.py \
  --donor-info donor_info.tsv \
  --sr-dna short_read1.tsv short_read2.tsv \
  --lr-dna long_read.tsv \
  --rna rna_watchmaker.tsv \
  --uberon tissue_uberon_identifiers.tsv \
  --submitter-prefix BROAD \
  --out-donor-tsv donor.tsv \
  --out-donor-xlsx donor.xlsx \
  --out-tissue-tsv tissue.tsv \
  --out-tissue-xlsx tissue.xlsx \
  --out-tissuesample-tsv ts.tsv \
  --out-tissuesample-xlsx ts.xlsx \
  --out-analyte-tsv analyte.tsv \
  --out-analyte-xlsx analyte.xlsx \
  --out-library-tsv library.tsv \
  --out-library-xlsx library.xlsx \
  --out-sequencing-tsv sequencing.tsv \
  --out-sequencing-xlsx sequencing.xlsx \
  --out-fileset-tsv fileset.tsv \
  --out-fileset-xlsx fileset.xlsx \
  --out-unalignedreads-tsv unaligned.tsv \
  --out-unalignedreads-xlsx unaligned.xlsx \
  --out-alignedreads-tsv aligned.tsv \
  --out-alignedreads-xlsx aligned.xlsx \
  --out-combined-xlsx all_metadata.xlsx
```

### CODEC Sequencing

For CODEC workflows, use the specialized script:

```bash
python generate_codec_metadata.py \
    --donor-info donor_information.tsv \
    --codec trimmed_codec.tsv \
    --uberon tissue_uberon_identifiers.tsv \
    --submitter-prefix BROAD \
    --out-donor-tsv output/donor.tsv \
    --out-donor-xlsx output/donor.xlsx \
    --out-tissue-tsv output/tissue.tsv \
    --out-tissue-xlsx output/tissue.xlsx \
    --out-tissuesample-tsv output/tissuesample.tsv \
    --out-tissuesample-xlsx output/tissuesample.xlsx \
    --out-analyte-tsv output/analyte.tsv \
    --out-analyte-xlsx output/analyte.xlsx \
    --out-library-tsv output/library.tsv \
    --out-library-xlsx output/library.xlsx \
    --out-library-preparation-tsv output/library_preparation.tsv \
    --out-library-preparation-xlsx output/library_preparation.xlsx \
    --out-sequencing-tsv output/sequencing.tsv \
    --out-sequencing-xlsx output/sequencing.xlsx \
    --out-fileset-tsv output/fileset.tsv \
    --out-fileset-xlsx output/fileset.xlsx \
    --out-unalignedreads-tsv output/unalignedreads.tsv \
    --out-unalignedreads-xlsx output/unalignedreads.xlsx \
    --out-alignedreads-tsv output/alignedreads.tsv \
    --out-alignedreads-xlsx output/alignedreads.xlsx \
    --out-variantcalls-tsv output/variantcalls.tsv \
    --out-variantcalls-xlsx output/variantcalls.xlsx \
    --out-software-tsv output/software.tsv \
    --out-software-xlsx output/software.xlsx \
    --out-combined-xlsx output/codec_metadata_combined.xlsx
```

**Recommended:** Create an output directory for organization:
```bash
mkdir output
```

---

## Examples

### Regular Sequencing Examples

1. Generate all sheets for a project with multiple sequencing types:
   ```bash
   python generate_metadata.py \
     --donor-info donor_info.tsv \
     --sr-dna short_read_batch1.tsv short_read_batch2.tsv \
     --lr-dna filtered_long.tsv \
     --rna filtered_rna_watchmaker.tsv \
     --uberon tissue_uberon_identifiers.tsv \
     [all output flags…] \
     --out-combined-xlsx project_metadata.xlsx
   ```
2. Process only RNA data:
 ```bash
python generate_metadata.py \
  --donor-info donor_info.tsv \
  --rna rna_truseq.tsv rna_watchmaker.tsv \
  --uberon tissue_uberon_identifiers.tsv \
  [all output flags…] \
  --out-combined-xlsx rna_only_metadata.xlsx
```

3. Inspect only the Library sheet:
   ```bash
   python generate_metadata.py … \
     --out-library-xlsx library_only.xlsx
   ```

### CODEC Examples

1. Complete CODEC workflow with preprocessing:
   ```bash
   # Step 1: Preprocess to trim file paths
   python utils/0A_strip_cram_and_bam_paths.py \
     --input raw_codec_data.tsv \
     --output trimmed_codec.tsv
   
   # Step 2: Generate CODEC metadata
   mkdir output
   python generate_codec_metadata.py \
     --donor-info donor_information.tsv \
     --codec trimmed_codec.tsv \
     --uberon tissue_uberon_identifiers.tsv \
     --submitter-prefix BROAD \
     [all CODEC output flags as shown above] \
     --out-combined-xlsx output/codec_metadata_combined.xlsx
   ```

2. Multiple CODEC input files:
   ```bash
   python generate_codec_metadata.py \
     --donor-info donor_info.tsv \
     --codec batch1_codec.tsv batch2_codec.tsv \
     --uberon tissue_uberon_identifiers.tsv \
     [remaining flags...]
   ```

---

## Future Development

- **Quality metrics**: Incorporate additional QC metrics (e.g. coverage, duplication rates) into the library and file-level sheets.  
- **Automation & CI**: Integrate with GitHub Actions or a scheduler to auto-run whenever new raw data files appear.  

---

## Repository Structure

```
metadata-pipeline/
├── scripts/                    # Python scripts for generating metadata sheets
│   ├— generate_metadata.py          # Regular sequencing workflows
│   ├— generate_codec_metadata.py    # CODEC sequencing workflows
│   ├— generate_fileset_sheet.py
│   ├— generate_unaligned_reads_sheet.py
│   └— generate_aligned_reads_sheet.py
├── utils/                      # Optional helper scripts for input preprocessing
│   └— 0A_strip_cram_and_bam_paths.py  # For both regular and CODEC files
├── input_examples/             # Sample input TSVs to test or guide formatting
│   ├— donor_info.tsv
│   ├— short_read_example.tsv
│   ├— long_read_example.tsv
│   ├— rna_watchmaker_example.tsv
│   ├— codec_example.tsv             # CODEC sample format
│   └— tissue_uberon_identifiers.tsv
├── files/                      # Reference and input data files
│   └— tissue_uberon_identifiers.tsv
├── requirements.txt            # Required Python packages
└— README.md                   # This documentation
```
# Other Branches:
## `full_sheet` Branch

The `full_sheet` branch includes support for generating **all required metadata sheets**, including empty placeholder sheets for those not yet populated. This ensures the final combined Excel workbook aligns fully with the SMaHT metadata submission template, even if some sheets are left blank for now.

You can switch to this branch using:

```bash
git checkout full_sheet
```
This branch introduces two additional required inputs:

- overview-path: Path to the overview guidelines Excel file

- template-path: Path to the submission template with all expected sheet names and order

## Example Usage
1. Generate all sheets for a project with multiple sequencing types:

```bash
python scripts/generate_metadata.py \
  --overview-path files/overview_guidelines.xlsx \
  --template-path files/gcc_automated_submission_example_v1.4.0.xlsx \
  --donor-info files/donor_info.tsv \
  --sr-dna short_read_batch1.tsv short_read_batch2.tsv \
  --lr-dna filtered_long.tsv \
  --rna filtered_rna_watchmaker.tsv \
  --uberon files/tissue_uberon_identifiers.tsv \
  --submitter-prefix BROAD \
  --out-donor-tsv output/donor_sheet.tsv \
  --out-donor-xlsx output/donor_sheet.xlsx \
  --out-tissue-tsv output/tissue_sheet.tsv \
  --out-tissue-xlsx output/tissue_sheet.xlsx \
  --out-tissuesample-tsv output/tissuesample_sheet.tsv \
  --out-tissuesample-xlsx output/tissuesample_sheet.xlsx \
  --out-analyte-tsv output/analyte_sheet.tsv \
  --out-analyte-xlsx output/analyte_sheet.xlsx \
  --out-library-tsv output/library_sheet.tsv \
  --out-library-xlsx output/library_sheet.xlsx \
  --out-fileset-tsv output/fileset_sheet.tsv \
  --out-fileset-xlsx output/fileset_sheet.xlsx \
  --out-unalignedreads-tsv output/unalignedreads_sheet.tsv \
  --out-unalignedreads-xlsx output/unalignedreads_sheet.xlsx \
  --out-alignedreads-tsv output/alignedreads_sheet.tsv \
  --out-alignedreads-xlsx output/alignedreads_sheet.xlsx \
  --out-combined-xlsx output/combined_metadata.xlsx
```

2. Process only RNA data:

```bash
 python scripts/generate_metadata.py \
  --donor-info files/donor_info.tsv \
  --uberon files/tissue_uberon_identifiers.tsv \
  --rna files/filtered_rna_watchmaker.tsv \
  --submitter-prefix BROAD \
  --out-donor-tsv output/donor.tsv \
  --out-donor-xlsx output/donor.xlsx \
  --out-tissue-tsv output/tissue.tsv \
  --out-tissue-xlsx output/tissue.xlsx \
  --out-tissuesample-tsv output/tissuesample.tsv \
  --out-tissuesample-xlsx output/tissuesample.xlsx \
  --out-analyte-tsv output/analyte.tsv \
  --out-analyte-xlsx output/analyte.xlsx \
  --out-library-tsv output/library.tsv \
  --out-library-xlsx output/library.xlsx \
  --out-fileset-tsv output/fileset.tsv \
  --out-fileset-xlsx output/fileset.xlsx \
  --out-unalignedreads-tsv output/unalignedreads.tsv \
  --out-unalignedreads-xlsx output/unalignedreads.xlsx \
  --out-alignedreads-tsv output/alignedreads.tsv \
  --out-alignedreads-xlsx output/alignedreads.xlsx \
  --out-combined-xlsx output/all_metadata.xlsx
```

---

## License

**License:** Broad Institute of MIT and Harvard
