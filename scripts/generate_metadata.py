#!/usr/bin/env python3
"""
generate_metadata.py (Modified)

Author: Shadi Zaheri
Date:   2025-04-25
Description:
    Orchestration script for generating Donor, Tissue, TissueSample,
    Analyte, Library, FileSet, UnalignedReads and AlignedReads sheets.
    
Modified to have explicit optional arguments for DNA sequencing types.
"""
__author__ = "Shadi Zaheri"
__version__ = "1.1"

import pandas as pd
import argparse

from generate_sequencing_sheet import generate_sequencing_sheet
from generate_fileset_sheet import generate_fileset_sheet
from generate_unaligned_reads_sheet import generate_unalignedreads_sheet
from generate_aligned_reads_sheet import generate_alignedreads_sheet


def normalize_tissue_label(tissue_label):
    """Normalize specific tissue labels by replacing underscores with hyphens."""
    if tissue_label is None:
        return tissue_label
    
    # Define specific tissues to normalize
    tissue_replacements = {
        'TEMPORAL_LOBE': 'TEMPORAL-LOBE',
        'FRONTAL_LOBE': 'FRONTAL-LOBE',
        'L_HIPPOCAMPUS': 'L-HIPPOCAMPUS',
        'R_HIPPOCAMPUS': 'R-HIPPOCAMPUS'
    }
    
    # First apply standard formatting
    formatted = tissue_label.replace(" ", "-").upper()
    
    # Then apply specific replacements
    for old_name, new_name in tissue_replacements.items():
        if formatted == old_name:
            return new_name
    
    return formatted


def format_tissue_submitted_id(donor_id, core_code, tissue_label):
    normalized_tissue = normalize_tissue_label(tissue_label)
    return f"NDRI_TISSUE_{core_code.upper()}-{normalized_tissue}"

def format_tissuesample_submitted_id(prefix, donor_id, tissue_label, collaborator_sample_id, category):
    sample_suffix = collaborator_sample_id.split("-")[-1]
    normalized_tissue = normalize_tissue_label(tissue_label)
    return f"{prefix}_TISSUE-SAMPLE_{donor_id}_{normalized_tissue}_{sample_suffix}_{category}".upper()


def format_analyte_id(prefix, label, sample_id):
    return f"{prefix}_ANALYTE_{label}_{sample_id}".replace(" ", "-").upper()

def infer_category(tissue_label):
    if isinstance(tissue_label, str):
        return "Liquid" if "BLOOD" in tissue_label.upper() else "Core"
    return "Core"

def generate_metadata_sheets(donor_info_tsv, sr_dna_files, lr_dna_files, rna_tsv, uberon_tsv,
                              out_donor_tsv, out_donor_xlsx,
                              out_tissue_tsv, out_tissue_xlsx,
                              out_tissuesample_tsv, out_tissuesample_xlsx,
                              submitter_prefix, out_analyte_tsv, out_analyte_xlsx, 
                              out_library_tsv, out_library_xlsx,
                              out_sequencing_tsv, out_sequencing_xlsx,
                              out_fileset_tsv, out_fileset_xlsx,
                              out_unalignedreads_tsv, out_unalignedreads_xlsx,
                              out_alignedreads_xlsx, out_alignedreads_tsv):
    
    uberon_df   = pd.read_csv(uberon_tsv, sep="\t")
    tissue_map  = dict(
        zip(
            uberon_df["tissue_identifier_code"],
            uberon_df["corresponding_tissue"]
        )
    )

    # Collect all input files with their types
    input_files = []
    if sr_dna_files:
        for f in sr_dna_files:
            input_files.append((f, "short_read"))
    if lr_dna_files:
        for f in lr_dna_files:
            input_files.append((f, "long_read"))

    # --- DONOR SHEET ---
    donor_df = pd.read_csv(donor_info_tsv, sep="\t")
    donor_df["submitted_id"] = donor_df["donor"].apply(lambda x: f"NDRI_DONOR_{x}")
    donor_df["external_id"] = donor_df["donor"]
    donor_df["sex"] = donor_df["gender"]
    donor_df["tpc_submitted"] = "True"
    donor_df["eligibility"] = ""
    donor_df["hardy_scale"] = ""
    donor_sheet = donor_df[[
        "submitted_id", "age", "external_id", "sex", "tpc_submitted", "eligibility", "hardy_scale"
    ]]
    donor_sheet.to_csv(out_donor_tsv, sep="\t", index=False)
    donor_sheet.to_excel(out_donor_xlsx, index=False)
    print(f"✅ Donor sheet saved to: {out_donor_tsv} and {out_donor_xlsx}")

    # --- TISSUE SHEET ---
    uberon_df = pd.read_csv(uberon_tsv, sep="\t")
    df_list = []
    
    # Process DNA files
    for f, file_type in input_files:
        df = pd.read_csv(f, sep="\t")
        df = df.dropna(subset=["collaborator_sample_id", "donor_id"])
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)
        df["donor_id"] = df["donor_id"].astype(str)
        df_list.append(df)
    
    # Process RNA files if provided
    if rna_tsv:
        for rna_file in rna_tsv:
            df_rna = pd.read_csv(rna_file, sep="\t")
            df_rna = df_rna.dropna(subset=["collaborator_sample_id", "donor_id"])
            df_rna["collaborator_sample_id"] = df_rna["collaborator_sample_id"].astype(str)
            df_rna["donor_id"] = df_rna["donor_id"].astype(str)
            df_list.append(df_rna)

    # Skip tissue processing if no input files
    if not df_list:
        print("❌ No input files provided. At least one of --sr-dna, --lr-dna, or --rna is required.")
        return

    df_all = pd.concat(df_list, ignore_index=True)
    df_all["core_id"] = df_all["collaborator_sample_id"].apply(lambda x: x.split("-")[1])
    df_all["external_id"] = df_all["collaborator_sample_id"].apply(lambda x: "-".join(x.split("-")[:2]))

    tissue_df = df_all[["collaborator_sample_id", "external_id", "core_id", "donor_id"]].drop_duplicates()
    tissue_df = tissue_df.merge(uberon_df, how="left", left_on="core_id", right_on="tissue_identifier_code")
    tissue_df["donor"] = tissue_df["donor_id"].apply(lambda x: f"NDRI_DONOR_{x}")
    tissue_df["submitted_id"] = tissue_df.apply(
        lambda row: format_tissue_submitted_id(row["donor_id"], row["external_id"], row["corresponding_tissue"])
        if pd.notnull(row["corresponding_tissue"]) else f"NDRI_TISSUE_{row['external_id']}",
        axis=1
    )

    tissue_sheet = pd.DataFrame(columns=[
        "submitted_id", "external_id", "sample_count", "preservation_type", "preservation_medium",
        "ischemic_time", "anatomical_location", "pathology_notes", "ph", "prosector_notes",
        "size", "size_unit", "volume", "weight", "donor", "uberon_id"
    ])
    tissue_sheet["submitted_id"] = tissue_df["submitted_id"]
    tissue_sheet["external_id"] = tissue_df["external_id"]
    tissue_sheet["donor"] = tissue_df["donor"]
    tissue_sheet["uberon_id"] = tissue_df["uberon_id"]
    tissue_sheet = tissue_sheet.drop_duplicates(subset=["submitted_id"])
    tissue_sheet.to_csv(out_tissue_tsv, sep="\t", index=False)
    tissue_sheet.to_excel(out_tissue_xlsx, index=False)
    print(f"✅ Tissue sheet saved to: {out_tissue_tsv} and {out_tissue_xlsx}")

    # --- TISSUE SAMPLE SHEET ---
    ext_to_sample_source = dict(zip(tissue_sheet["external_id"], tissue_sheet["submitted_id"]))
    
    # FIXED: Get tissue labels directly from tissue_map to avoid extraction issues
    ext_to_label = {}
    for _, row in df_all[["external_id", "core_id"]].drop_duplicates().iterrows():
        core_id = row["core_id"]
        external_id = row["external_id"]
        # Get tissue label directly from tissue_map and normalize it
        raw_tissue = tissue_map.get(core_id, core_id)
        tissue_label = normalize_tissue_label(raw_tissue)
        ext_to_label[external_id] = tissue_label

    ext_to_donor = {
        row["external_id"]: row["submitted_id"].split("_TISSUE_")[1].split("-")[0]
        for _, row in tissue_sheet.iterrows()
    }

    df_all["sample_sources"] = df_all["external_id"].map(ext_to_sample_source)
    df_all["tissue_label"] = df_all["external_id"].map(ext_to_label)
    df_all["donor"] = df_all["external_id"].map(ext_to_donor)
    df_all["category"] = df_all["tissue_label"].apply(infer_category)
    df_all["submitted_id"] = df_all.apply(
        lambda row: format_tissuesample_submitted_id(submitter_prefix, row["donor"], row["tissue_label"], row["collaborator_sample_id"], row["category"]),
        axis=1
    )

    tissuesample_sheet = pd.DataFrame(columns=[
        "submitted_id", "category", "external_id", "preservation_type", "preservation_medium", "description",
        "core_size", "processing_date", "processing_notes", "weight", "sample_sources", "parent_samples"
    ])
    tissuesample_sheet["submitted_id"] = df_all["submitted_id"]
    tissuesample_sheet["category"] = df_all["category"]
    tissuesample_sheet["external_id"] = df_all["collaborator_sample_id"]
    tissuesample_sheet["sample_sources"] = df_all["sample_sources"]
    tissuesample_sheet = tissuesample_sheet.drop_duplicates(subset=["submitted_id"])
    tissuesample_sheet.to_csv(out_tissuesample_tsv, sep="\t", index=False)
    tissuesample_sheet.to_excel(out_tissuesample_xlsx, index=False)
    print(f"✅ TissueSample sheet saved to: {out_tissuesample_tsv} and {out_tissuesample_xlsx}")

    # --- ANALYTE SHEET ---
    file_labels = {}
    
    # Label DNA files based on explicit type
    if sr_dna_files:
        for f in sr_dna_files:
            file_labels[f] = "GDNA_SHORTREAD"
    if lr_dna_files:
        for f in lr_dna_files:
            file_labels[f] = "GDNA_LONGREAD"
    
    # Label RNA files
    if rna_tsv:
        for rna_file in rna_tsv:
            if "truseq" in rna_file.lower():
                file_labels[rna_file] = "RNA_TRUSEQ"
            else:
                file_labels[rna_file] = "RNA_WATCHMAKER"

    analyte_records = []
    for f, label in file_labels.items():
        df = pd.read_csv(f, sep="\t").dropna(subset=["collaborator_sample_id"])
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)

        for sample_id in df["collaborator_sample_id"]:
            ts_match = df_all[df_all["collaborator_sample_id"] == sample_id]
            if ts_match.empty:
                continue

            if label == "RNA_WATCHMAKER":
                molecule, molecule_detail = "RNA", "Total RNA"
            elif label == "RNA_TRUSEQ":
                molecule, molecule_detail = "RNA", "mRNA"
            else:
                molecule, molecule_detail = "DNA", "Total DNA"

            donor_part = sample_id.split("-")[0]
            core = sample_id.split("-")[1]
            tissue_label = normalize_tissue_label(tissue_map.get(core, core))

            analyte_id = (
                f"{submitter_prefix}_ANALYTE_"
                f"{donor_part}_{tissue_label}_{label}"
            ).upper()

            analyte_records.append({
                "submitted_id":   analyte_id,
                "molecule":       molecule,
                "molecule_detail":molecule_detail,
                "external_id":    sample_id,
                "samples":        ts_match["submitted_id"].values[0]
            })

    analyte_df = pd.DataFrame(analyte_records)
    analyte_columns = [
        "submitted_id", "molecule", "molecule_detail", "external_id", "description", "a260_a280_ratio",
        "average_fragment_size", "concentration", "concentration_unit", "dna_integrity_number",
        "dna_integrity_number_instrument", "dna_quality_number", "dna_quality_number_instrument",
        "dna_quality_size_threshold", "genomic_quality_number", "genomic_quality_number_instrument",
        "genomic_quality_size_threshold", "quantitation_method", "ribosomal_rna_ratio",
        "rna_integrity_number", "rna_integrity_number_instrument", "sample_quantity", "sample_quantity_unit",
        "volume", "volume_unit", "total_yield", "yield_per_unit", "yield_unit", "samples", "analyte_preparation"
    ]
    analyte_sheet = pd.DataFrame(columns=analyte_columns)
    analyte_sheet["submitted_id"] = analyte_df["submitted_id"]
    analyte_sheet["molecule"] = analyte_df["molecule"]
    analyte_sheet["molecule_detail"] = analyte_df["molecule_detail"]
    analyte_sheet["external_id"] = analyte_df["external_id"]
    analyte_sheet["samples"] = analyte_df["samples"]
    analyte_sheet = analyte_sheet.drop_duplicates(subset=["submitted_id"])
    analyte_sheet.to_csv(out_analyte_tsv, sep="\t", index=False)
    analyte_sheet.to_excel(out_analyte_xlsx, index=False)
    print(f"✅ Analyte sheet saved to: {out_analyte_tsv} and {out_analyte_xlsx}")

    # --- LIBRARY SHEET ---
    analyte_df = analyte_df.copy()

    def get_library_type(analyt_id):
        if "SHORTREAD" in analyt_id:
            return "WGS_ILLUMINA"
        elif "LONGREAD" in analyt_id:
            return "WGS_PACBIO"
        elif "WATCHMAKER" in analyt_id:
            return "RNA_WATCHMAKER"
        elif "TRUSEQ" in analyt_id:
            return "RNA_TRUSEQ"
        else:
            return "UNKNOWN"

    def infer_assay(library_type):
        if library_type in {"RNA_WATCHMAKER", "RNA_TRUSEQ"}:
            return "bulk_rna_seq"
        elif library_type in {"WGS_ILLUMINA", "WGS_PACBIO"}:
            return "bulk_wgs_pcr_free"
        else:
            return "other"

    lib_counter = {}
    library_records = []

    for _, row in analyte_df.iterrows():
        analyte_id = row["submitted_id"]
        ext_id = row["external_id"]

        donor = ext_id.split("-")[0]
        core = ext_id.split("-")[1]  # Extract core from collaborator_sample_id
        tissue_label = normalize_tissue_label(tissue_map.get(core, core))
        library_type = get_library_type(analyte_id)
        assay = infer_assay(library_type)

        lib_key = (donor, tissue_label, library_type)
        lib_counter[lib_key] = lib_counter.get(lib_key, 0) + 1
        lib_suffix = f"_{lib_counter[lib_key]}" 

        lib_id = f"{submitter_prefix}_LIBRARY_{donor}_{tissue_label}_{library_type}{lib_suffix}".replace(" ", "-").upper()

        library_records.append({
            "submitted_id": lib_id,
            "analytes": analyte_id,
            "assay": assay
        })

    library_df = pd.DataFrame(library_records)

    library_columns_order = [
        "submitted_id", "comments", "external_id", "description", "a260_a280_ratio", "adapter_name",
        "adapter_sequence", "amplification_cycles", "analyte_weight", "antibody", "barcode_sequences",
        "concatenated_reads", "fragment_mean_length", "guide_sequence", "insert_coefficient_of_variation",
        "insert_maximum_length", "insert_mean_length", "insert_minimum_length", "preparation_date",
        "dna_target", "target_fragment_size", "target_insert_maximum_length", "target_insert_minimum_length",
        "target_monomer_size", "analytes", "assay", "library_preparation"
    ]

    for col in library_columns_order:
        if col not in library_df.columns:
            library_df[col] = ""

    library_df = library_df[library_columns_order]
    library_df.to_csv(out_library_tsv, sep="\t", index=False)
    library_df.to_excel(out_library_xlsx, index=False)
    print(f"✅ Library sheet saved to: {out_library_tsv} and {out_library_xlsx}")
    
     # --- SEQUENCING SHEET --- (ADD THIS SECTION)
    # Create temporary args object for sequencing sheet generation
    seq_args = argparse.Namespace()
    seq_args.sr_dna = sr_dna_files
    seq_args.lr_dna = lr_dna_files  
    seq_args.rna = rna_tsv
    seq_args.out_sequencing_tsv = out_sequencing_tsv
    seq_args.out_sequencing_xlsx = out_sequencing_xlsx
    
    generate_sequencing_sheet(seq_args)

def main():
    parser = argparse.ArgumentParser(
        description="Generate metadata sheets including Donor, Tissue, TissueSample, "
                    "Analyte, Library, FileSet, UnalignedReads, and AlignedReads."
    )
    parser.add_argument(
        "--donor-info", required=True,
        help="Path to the donor information TSV (with columns: donor, gender, age, etc.)"
    )
    parser.add_argument(
        "--sr-dna", nargs="*", default=None,
        help="Short-read DNA sequencing TSV files (optional)"
    )
    parser.add_argument(
        "--lr-dna", nargs="*", default=None,
        help="Long-read DNA sequencing TSV files (optional)"
    )
    parser.add_argument(
        "--rna", nargs="*", default=None,
        help="RNA sequencing TSV files (filtered_rna_watchmaker or filtered_rna_truseq, optional)"
    )
    parser.add_argument(
        "--uberon", required=True,
        help="Path to the tissue_uberon_identifiers TSV mapping core codes to tissue names"
    )
    parser.add_argument(
        "--submitter-prefix", default="BROAD",
        help="Prefix to use in all submitted_id fields (default: BROAD)"
    )
    parser.add_argument(
        "--out-donor-tsv", required=True,
        help="Where to write the Donor sheet TSV"
    )
    parser.add_argument(
        "--out-donor-xlsx", required=True,
        help="Where to write the Donor sheet XLSX"
    )
    parser.add_argument(
        "--out-tissue-tsv", required=True,
        help="Where to write the Tissue sheet TSV"
    )
    parser.add_argument(
        "--out-tissue-xlsx", required=True,
        help="Where to write the Tissue sheet XLSX"
    )
    parser.add_argument(
        "--out-tissuesample-tsv", required=True,
        help="Where to write the TissueSample sheet TSV"
    )
    parser.add_argument(
        "--out-tissuesample-xlsx", required=True,
        help="Where to write the TissueSample sheet XLSX"
    )
    parser.add_argument(
        "--out-analyte-tsv", required=True,
        help="Where to write the Analyte sheet TSV"
    )
    parser.add_argument(
        "--out-analyte-xlsx", required=True,
        help="Where to write the Analyte sheet XLSX"
    )
    parser.add_argument(
        "--out-library-tsv", required=True,
        help="Where to write the Library sheet TSV"
    )
    parser.add_argument(
        "--out-library-xlsx", required=True,
        help="Where to write the Library sheet XLSX"
    )
    parser.add_argument(
        "--out-sequencing-tsv", required=True,
        help="Where to write the Sequencing sheet TSV"
    )
    parser.add_argument(
        "--out-sequencing-xlsx", required=True,
        help="Where to write the Sequencing sheet XLSX"
    )
    parser.add_argument(
        "--out-fileset-tsv", required=True,
        help="Where to write the FileSet sheet TSV"
    )
    parser.add_argument(
        "--out-fileset-xlsx", required=True,
        help="Where to write the FileSet sheet XLSX"
    )
    parser.add_argument(
        "--out-unalignedreads-tsv", required=True,
        help="Where to write the UnalignedReads sheet TSV"
    )
    parser.add_argument(
        "--out-unalignedreads-xlsx", required=True,
        help="Where to write the UnalignedReads sheet XLSX"
    )
    parser.add_argument(
        "--out-alignedreads-tsv", required=True,
        help="Where to write the AlignedReads sheet TSV"
    )
    parser.add_argument(
        "--out-alignedreads-xlsx", required=True,
        help="Where to write the AlignedReads sheet XLSX"
    )
    parser.add_argument(
        "--out-combined-xlsx", required=True,
        help="Where to write the combined Excel workbook containing all sheets"
    )

    args = parser.parse_args()
    
    # Check that at least one sequencing type is provided
    if not any([args.sr_dna, args.lr_dna, args.rna]):
        parser.error("At least one of --sr-dna, --lr-dna, or --rna must be provided")
    
    generate_metadata_sheets(
        args.donor_info, args.sr_dna, args.lr_dna, args.rna, args.uberon,
        args.out_donor_tsv, args.out_donor_xlsx,
        args.out_tissue_tsv, args.out_tissue_xlsx,
        args.out_tissuesample_tsv, args.out_tissuesample_xlsx,
        args.submitter_prefix,
        args.out_analyte_tsv, args.out_analyte_xlsx, 
        args.out_library_tsv, args.out_library_xlsx,
        args.out_sequencing_tsv, args.out_sequencing_xlsx,
        args.out_fileset_tsv, args.out_fileset_xlsx,
        args.out_unalignedreads_tsv, args.out_unalignedreads_xlsx,
        args.out_alignedreads_tsv, args.out_alignedreads_xlsx
    )
    
    # Create a temporary args object for the helper functions
    # They expect an 'inputs' attribute, so we'll create it from our new structure
    helper_args = argparse.Namespace()
    for attr, value in vars(args).items():
        setattr(helper_args, attr, value)
    
    # Create the inputs list for helper functions
    helper_args.inputs = []
    if args.sr_dna:
        helper_args.inputs.extend(args.sr_dna)
    if args.lr_dna:
        helper_args.inputs.extend(args.lr_dna)
    
    # For backward compatibility with helper functions that expect single RNA file,
    # we'll pass the first RNA file if multiple are provided
    if args.rna:
        helper_args.rna = args.rna[0]  # Use first RNA file for helpers
        # DO NOT add RNA files to inputs - they are processed separately
    
    generate_fileset_sheet(helper_args)
    generate_unalignedreads_sheet(helper_args)
    generate_alignedreads_sheet(helper_args)
    
    with pd.ExcelWriter(args.out_combined_xlsx, engine="openpyxl") as writer:
        df_donor = pd.read_excel(args.out_donor_xlsx, dtype={"tpc_submitted": str})
        df_donor["tpc_submitted"] = df_donor["tpc_submitted"].str.capitalize()
        df_donor.to_excel(writer, sheet_name="(Donor)", index=False)
        pd.read_excel(args.out_tissue_xlsx)       .to_excel(writer, sheet_name="(Tissue)",          index=False)
        pd.read_excel(args.out_tissuesample_xlsx) .to_excel(writer, sheet_name="TissueSample",    index=False)
        pd.read_excel(args.out_analyte_xlsx)      .to_excel(writer, sheet_name="Analyte",         index=False)
        pd.read_excel(args.out_library_xlsx)      .to_excel(writer, sheet_name="Library",         index=False)
        pd.read_excel(args.out_sequencing_xlsx)   .to_excel(writer, sheet_name="Sequencing",      index=False)
        pd.read_excel(args.out_fileset_xlsx)      .to_excel(writer, sheet_name="FileSet",         index=False)
        pd.read_excel(args.out_unalignedreads_xlsx).to_excel(writer, sheet_name="UnalignedReads",   index=False)
        pd.read_excel(args.out_alignedreads_xlsx) .to_excel(writer, sheet_name="AlignedReads",     index=False)
    print(f"✅ The combined metadata Excel is saved to: {args.out_combined_xlsx}")

if __name__ == "__main__":
    main()
