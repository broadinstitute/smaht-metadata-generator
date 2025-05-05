# smaht-metadata-generator
Pipeline to generate standardized metadata sheets for SMaHT sequencing projects
# Metadata Pipeline

This repository provides a **pipeline** for generating standardized metadata sheets required for genomic sequencing submissions, covering Donor, Tissue, TissueSample, Analyte, Library, FileSet, UnalignedReads, and AlignedReads.

---

## Overview

This pipeline reads in: donor metadata, sample tables for short‑read, long‑read, and RNA inputs, plus a UBERON mapping file.  It then produces eight metadata sheets with properly formatted `submitted_id` values:

- **Donor**: one row per donor
- **Tissue**: one row per core tissue code
- **TissueSample**: one row per individual sample
- **Analyte**: DNA/RNA preparations
- **Library**: sequencing libraries
- **FileSet**: grouping of libraries
- **UnalignedReads**: raw read files (BAM)
- **AlignedReads**: alignment outputs (CRAM)

These are output both as TSV and XLSX.
---

## Optional Preprocessing Utilities

To assist with metadata preparation, we provide optional helper scripts in the [`utils/`](./utils) directory. These are **not required**, but they can help standardize and automate input formatting prior to running the metadata generation pipeline.

---

### `0B_strip_cram_and_bam_paths.py`

This script strips full file paths in your sample tables down to just the base file name (e.g., `gs://bucket/data/sample.bam` → `sample.bam`).  
It works for columns commonly used in Unaligned and Aligned reads inputs, such as:

- `bam_path`
- `bam_index_path`
- `cram_path`
- `cram_index_path`

**Input:**  
A tab-delimited `.tsv` file with at least the following columns:

| Column            | Description                            |
|-------------------|----------------------------------------|
| `sample_id`       | Unique ID for the sample (any format)  |
| `bam_path`        | Full BAM path (e.g. `gs://...`)        |
| `bam_index_path`  | Full BAM index path (optional)         |
| `cram_path`       | Full CRAM path (e.g. `gs://...`)       |
| `cram_index_path` | Full CRAM index path (optional)        |

**Output:**  
A `.tsv` with the same structure, but with only the file names in the above fields.

**Example:**

**Input:**

| sample_id | bam_path                          | cram_path                          |
|-----------|-----------------------------------|------------------------------------|
| SMHT001   | gs://bucket/data/SMHT001.bam      | gs://bucket/data/SMHT001.cram      |

**Output:**

| sample_id | bam_path     | cram_path     |
|-----------|--------------|---------------|
| SMHT001   | SMHT001.bam  | SMHT001.cram  |

**Usage:**
```bash
python utils/0B_strip_cram_and_bam_paths.py <input_tsv> <output_tsv>

---

## Supported Sequencing Technologies

- **WGS_ILLUMINA** (short‑read)  
- **WGS_PACBIO**  (long‑read)  
- **RNA_WATCHMAKER** (RNA total)  
- **RNA_TRUSEQ**     (RNA mRNA)  

The pipeline infers these from file names or analyte IDs.

---

## Input Files

### Donor Info TSV

- Provided via `--donor-info`  
- Must include at least columns:  
  - `donor` (unique donor code, e.g. `SMHT001`)  
  - `gender` (`M`/`F`)  
  - `age` (integer)

### Sample TSVs (Short/Long/RNA)

- Passed via `--inputs` (one or more) and `--rna` (exactly one).  Filenames should reflect their type:
  - Short‑read tables: filename contains `short_read` or similar
  - Long‑read tables: filename contains `long_read` or `long`
  - RNA tables: filename contains `watchmaker` or `truseq`

Each must include at least:

| Column                 | Description                                  |
|------------------------|----------------------------------------------|
| `collaborator_sample_id` | `{donor}-{core}-{suffix}` (e.g. `SMHT001-3S-001A1`) |
| `donor_id`               | Matching donor code (e.g. `SMHT001`)           |

Additional columns (e.g. `input_bam`, `cram_path`) are used for unaligned/aligned sheets.

### UBERON Mapping TSV

- Provided via `--uberon`  
- Must include:
  - `tissue_identifier_code` (e.g. `3S`)  
  - `corresponding_tissue`   (e.g. `Heart`)
  
The input_examples/ folder includes sample TSVs that demonstrate the expected format for:

Donor info files

Short-read, long-read, and RNA sample tables

UBERON tissue mapping files

These examples are intended to help you test the pipeline or structure your own metadata inputs.
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
| Combined Excel  | `--out-combined-xlsx`  | Multi‑sheet workbook of all above      |

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

```bash
python generate_metadata.py \
  --donor-info donor_info.tsv \
  --inputs short.tsv long.tsv \
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
  --out-fileset-tsv fileset.tsv \
  --out-fileset-xlsx fileset.xlsx \
  --out-unalignedreads-tsv unaligned.tsv \
  --out-unalignedreads-xlsx unaligned.xlsx \
  --out-alignedreads-tsv aligned.tsv \
  --out-alignedreads-xlsx aligned.xlsx \
  --out-combined-xlsx all_metadata.xlsx
```

---

## Examples

1. Generate all sheets for a project and combine:
   ```bash
   python generate_metadata.py \
     --donor-info donor_info.tsv \
     --inputs filtered_long.tsv short_read_stripped.tsv \
     --rna filtered_rna_watchmaker.tsv \
     --uberon tissue_uberon_identifiers.tsv \
     [all output flags…] \
     --out-combined-xlsx project_metadata.xlsx
   ```

2. Inspect only the Library sheet:
   ```bash
   python generate_metadata.py … \
     --out-library-xlsx library_only.xlsx
   ```

---

## Future Development

- **CODEC support**: Add support for CODEC sequencing workflows, including custom sheet templates for CODEC-specific libraries and metadata fields.  
- **VariantFiles sheets**: Generate VARIANTCALLS metadata sheets (`VariantFiles`) with fields like `submitted_id`, `data_category`, `data_type`, `filename`, `submitted_md5sum`,	`file_format`, and	`software`.
- **Quality metrics**: Incorporate additional QC metrics (e.g. coverage, duplication rates) into the library and file-level sheets.  
- **Automation & CI**: Integrate with GitHub Actions or a scheduler to auto-run whenever new raw data files appear.  

---
---

## Repository Structure

```
metadata-pipeline/
├── scripts/                    # Python scripts for generating metadata sheets
│   ├── generate_metadata.py
│   ├── generate_fileset_sheet.py
│   ├── generate_unaligned_reads_sheet.py
│   └── generate_aligned_reads_sheet.py
├── files/                      # Reference and input data files
│   └── tissue_uberon_identifiers.tsv
├── input_examples/                    # Example of input files
│   ├── Short-read, long-read, and RNA sample tables
│   ├── Donor info files
│   └── UBERON tissue mapping files
├── README.md                   # This documentation
```

---

## License

**License:** Broad Institute of MIT and Harvard

