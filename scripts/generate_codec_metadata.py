#!/usr/bin/env python3
"""
generate_codec_metadata.py

Author: Shadi Zaheri
Date: 2025-08-15
Description:
    Specialized metadata generator for CODEC samples.
    Generates Donor, Tissue, TissueSample, Analyte, Library, FileSet,
    AlignedReads, VariantCalls, and Software sheets specifically for CODEC data.
    
    FIXES APPLIED:
    - Correct tissue label extraction preserving compound names (L-HIPPOCAMPUS, etc.)
    - Include sample suffixes in IDs to ensure uniqueness
    - Maintain 1:1:1:1 hierarchical relationships
    - Preserve selective tissue normalization (only 4 special tissues get hyphens)
"""

import pandas as pd
import argparse


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


def infer_category(tissue_label):
    if isinstance(tissue_label, str):
        return "Liquid" if "BLOOD" in tissue_label.upper() else "Core"
    return "Core"


def generate_codec_metadata_sheets(donor_info_tsv, codec_files, uberon_tsv,
                                   out_donor_tsv, out_donor_xlsx,
                                   out_tissue_tsv, out_tissue_xlsx,
                                   out_tissuesample_tsv, out_tissuesample_xlsx,
                                   submitter_prefix, out_analyte_tsv, out_analyte_xlsx,
                                   out_library_tsv, out_library_xlsx,
                                   out_fileset_tsv, out_fileset_xlsx,
                                   out_unalignedreads_tsv, out_unalignedreads_xlsx,
                                   out_alignedreads_tsv, out_alignedreads_xlsx,
                                   out_variantcalls_tsv, out_variantcalls_xlsx,
                                   out_software_tsv, out_software_xlsx,
                                   out_sequencing_tsv, out_sequencing_xlsx,
                                   out_library_preparation_tsv, out_library_preparation_xlsx):
    
    # Load tissue mapping
    uberon_df = pd.read_csv(uberon_tsv, sep="\t")
    tissue_map = dict(zip(uberon_df["tissue_identifier_code"], 
                         uberon_df["corresponding_tissue"]))

    # --- DONOR SHEET (same as regular) ---
    donor_df = pd.read_csv(donor_info_tsv, sep="\t")
    donor_df["submitted_id"] = donor_df["donor"].apply(lambda x: f"NDRI_DONOR_{x}")
    donor_df["external_id"] = donor_df["donor"]
    donor_df["sex"] = donor_df["gender"]
    donor_df["tpc_submitted"] = "True"
    donor_df["eligibility"] = ""
    donor_df["hardy_scale"] = ""
    donor_sheet = donor_df[["submitted_id", "age", "external_id", "sex", "tpc_submitted", "eligibility", "hardy_scale"]]
    donor_sheet.to_csv(out_donor_tsv, sep="\t", index=False)
    donor_sheet.to_excel(out_donor_xlsx, index=False)
    print(f"✅ Donor sheet saved to: {out_donor_tsv} and {out_donor_xlsx}")

    # Load CODEC data
    df_list = []
    for codec_file in codec_files:
        df = pd.read_csv(codec_file, sep="\t")
        df = df.dropna(subset=["collaborator_sample_id", "donor_id"])
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)
        df["donor_id"] = df["donor_id"].astype(str)
        df_list.append(df)

    if not df_list:
        print("❌ No CODEC input files provided.")
        return

    df_all = pd.concat(df_list, ignore_index=True)
    df_all["core_id"] = df_all["collaborator_sample_id"].apply(lambda x: x.split("-")[1])
    df_all["external_id"] = df_all["collaborator_sample_id"].apply(lambda x: "-".join(x.split("-")[:2]))

    # --- TISSUE SHEET (FIXED: consistent normalization) ---
    tissue_df = df_all[["collaborator_sample_id", "external_id", "core_id", "donor_id"]].drop_duplicates()
    tissue_df = tissue_df.merge(uberon_df, how="left", left_on="core_id", right_on="tissue_identifier_code")
    tissue_df["donor"] = tissue_df["donor_id"].apply(lambda x: f"NDRI_DONOR_{x}")
    
    # FIXED: Apply normalization consistently
    tissue_df["submitted_id"] = tissue_df.apply(
        lambda row: f"NDRI_TISSUE_{row['external_id']}-{normalize_tissue_label(row['corresponding_tissue'])}"
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

    # --- TISSUE SAMPLE SHEET (FIXED: correct tissue label extraction) ---
    ext_to_sample_source = dict(zip(tissue_sheet["external_id"], tissue_sheet["submitted_id"]))
    
    # FIXED: Correct tissue label extraction that preserves compound names
    ext_to_label = {}
    for _, row in tissue_sheet.iterrows():
        tissue_id = row["submitted_id"]  # e.g., "NDRI_TISSUE_SMHT001-001A1-L-HIPPOCAMPUS"
        if tissue_id.startswith("NDRI_TISSUE_"):
            # Remove prefix to get "SMHT001-001A1-L-HIPPOCAMPUS"
            remainder = tissue_id.replace("NDRI_TISSUE_", "")
            # Split on "-" and take everything after the second part
            parts = remainder.split("-")
            if len(parts) >= 3:  # donor-core-tissue format
                tissue_label = "-".join(parts[2:])  # Join remaining parts for compound tissues
            else:
                tissue_label = parts[-1]  # Fallback to last part
            ext_to_label[row["external_id"]] = tissue_label
        else:
            # Fallback for unexpected formats
            ext_to_label[row["external_id"]] = tissue_id.split("-")[-1]

    ext_to_donor = {row["external_id"]: row["submitted_id"].split("_TISSUE_")[1].split("-")[0]
                   for _, row in tissue_sheet.iterrows()}

    df_all["sample_sources"] = df_all["external_id"].map(ext_to_sample_source)
    df_all["tissue_label"] = df_all["external_id"].map(ext_to_label)
    df_all["donor"] = df_all["external_id"].map(ext_to_donor)
    df_all["category"] = df_all["tissue_label"].apply(infer_category)
    df_all["submitted_id"] = df_all.apply(
        lambda row: format_tissuesample_submitted_id(submitter_prefix, row["donor"], 
                                                   row["tissue_label"], row["collaborator_sample_id"], 
                                                   row["category"]), axis=1
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

    # --- ANALYTE SHEET (FIXED: include sample suffix for uniqueness) ---
    analyte_records = []

    for _, ts_row in tissuesample_sheet.iterrows():
        tissue_sample_id = ts_row["submitted_id"]
        sample_id = ts_row["external_id"]
        
        donor_part = sample_id.split("-")[0]
        core = sample_id.split("-")[1] 
        sample_suffix = sample_id.split("-")[-1]  # e.g., "001A1"
        tissue_label = normalize_tissue_label(tissue_map.get(core, core))

        # FIXED: Include sample suffix to ensure uniqueness
        analyte_id = f"{submitter_prefix}_ANALYTE_{donor_part}_{tissue_label}_{sample_suffix}_GDNA_SHORTREAD".upper()

        analyte_records.append({
            "submitted_id": analyte_id,
            "molecule": "DNA", 
            "molecule_detail": "Total DNA", 
            "external_id": sample_id,
            "samples": tissue_sample_id
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

    # --- LIBRARY SHEET (FIXED: use sample suffixes, remove counter logic) ---
    codec_protocols = ["DDBTP"]
    library_records = []

    for _, analyte_row in analyte_df.iterrows():
        analyte_id = analyte_row["submitted_id"]
        sample_id = analyte_row["external_id"]
        
        donor = sample_id.split("-")[0]
        core = sample_id.split("-")[1]
        sample_suffix = sample_id.split("-")[-1]
        tissue_label = normalize_tissue_label(tissue_map.get(core, core))
        
        for protocol in codec_protocols:
            # FIXED: Include sample suffix for uniqueness - no counter needed
            lib_id = f"{submitter_prefix}_LIBRARY_CODEC_{donor}_{tissue_label}-{protocol}_{sample_suffix}".upper()
            lib_prep_id = f"{submitter_prefix}_LIBRARY-PREPARATION_CODEC_{protocol}".upper()

            library_records.append({
                "submitted_id": lib_id,
                "analytes": analyte_id,
                "assay": "codec", 
                "library_preparation": lib_prep_id
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

    # --- FILESET SHEET (FIXED: align with new library IDs) ---
    fileset_records = []
    
    for _, analyte_row in analyte_df.iterrows():
        sample_id = analyte_row["external_id"]
        donor = sample_id.split("-")[0]
        core = sample_id.split("-")[1] 
        sample_suffix = sample_id.split("-")[-1]
        tissue_label = normalize_tissue_label(tissue_map.get(core, core))
        
        for protocol in codec_protocols:
            # FIXED: Use sample suffix to match library IDs
            fileset_id = f"{submitter_prefix}_FILE-SET_CODEC_{donor}_{tissue_label}-{protocol}_{sample_suffix}".upper()
            lib_id = f"{submitter_prefix}_LIBRARY_CODEC_{donor}_{tissue_label}-{protocol}_{sample_suffix}".upper()
            
            fileset_records.append({
                "submitted_id": fileset_id,
                "description": "",
                "submitter_comments": "",
                "libraries": lib_id,
                "sequencing": f"{submitter_prefix}_SEQUENCING_CODEC".upper(),
                "samples": ""
            })

    # Create DataFrame with correct column order
    fileset_columns = ["submitted_id", "description", "submitter_comments", "libraries", "sequencing", "samples"]
    fileset_df = pd.DataFrame(fileset_records)
    fileset_df = fileset_df[fileset_columns]
    fileset_df.to_csv(out_fileset_tsv, sep="\t", index=False)
    fileset_df.to_excel(out_fileset_xlsx, index=False)
    print(f"✅ FileSet sheet saved to: {out_fileset_tsv} and {out_fileset_xlsx}")

    # --- UNALIGNED READS SHEET (Empty for CODEC) ---
    unaligned_cols = [
        "submitted_id", "filename", "submitted_md5sum", "data_category", "data_type", "n50",
        "flow_cell_barcode", "flow_cell_lane", "description", "read_pair_number", "file_format",
        "file_sets", "derived_from", "external_quality_metrics", "paired_with", "software"
    ]
    empty_unaligned_df = pd.DataFrame(columns=unaligned_cols)
    empty_unaligned_df.to_csv(out_unalignedreads_tsv, sep="\t", index=False)
    empty_unaligned_df.to_excel(out_unalignedreads_xlsx, index=False)
    print(f"✅ Empty UnalignedReads sheet saved to: {out_unalignedreads_tsv} and {out_unalignedreads_xlsx}")

    # --- ALIGNED READS SHEET (FIXED: use sample suffixes for consistent linking) ---
    aligned_records = []
    
    for codec_file in codec_files:
        df = pd.read_csv(codec_file, sep="\t")
        df = df.dropna(subset=["collaborator_sample_id"])
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)
        
        for _, row in df.iterrows():
            sample_id = row["collaborator_sample_id"]
            donor = sample_id.split("-")[0]
            core = sample_id.split("-")[1]
            sample_suffix = sample_id.split("-")[-1]
            tissue_label = normalize_tissue_label(tissue_map.get(core, core))
            
            for protocol in codec_protocols:
                # FIXED: Use sample suffix to match fileset IDs
                fileset_id = f"{submitter_prefix}_FILE-SET_CODEC_{donor}_{tissue_label}-{protocol}_{sample_suffix}".upper()
                
                # RAW BAM (from RAW_BAM column)
                if pd.notna(row.get("RAW_BAM")):
                    raw_id = f"{submitter_prefix}_ALIGNED-READS_RAW_CODEC-{sample_id}-{protocol}_{sample_suffix}".upper()
                    aligned_records.append({
                        "submitted_id": raw_id,
                        "filename": row["RAW_BAM"],
                        "submitted_md5sum": "",
                        "data_category": "Sequencing Reads",
                        "data_type": "",
                        "n50": "",
                        "flow_cell_barcode": "",
                        "flow_cell_lane": "",
                        "description": f"Raw bam CODEC_{sample_id}-{protocol.lower()}",
                        "alignment_details": "Sorted",
                        "file_format": "bam",
                        "file_sets": fileset_id,
                        "reference_genome": "broad_grch38",
                        "derived_from": "",
                        "external_quality_metrics": "",
                        "software": f"{submitter_prefix}_SOFTWARE_CODEC_ALIGNMENT_PIPELINE_V1".upper()
                    })
                
                # PROCESSED BAM (from MolConsensusBAM column)
                if pd.notna(row.get("MolConsensusBAM")):
                    processed_id = f"{submitter_prefix}_ALIGNED-READS_PROCESSED_CODEC-{sample_id}-{protocol}_{sample_suffix}".upper()
                    aligned_records.append({
                        "submitted_id": processed_id,
                        "filename": row["MolConsensusBAM"],
                        "submitted_md5sum": "",
                        "data_category": "Consensus Reads",
                        "data_type": "",
                        "n50": "",
                        "flow_cell_barcode": "",
                        "flow_cell_lane": "",
                        "description": f"Processed duplex bam CODEC_{sample_id}-{protocol}",
                        "alignment_details": "Sorted",
                        "file_format": "bam",
                        "file_sets": fileset_id,
                        "reference_genome": "broad_grch38",
                        "derived_from": raw_id if pd.notna(row.get("RAW_BAM")) else "",
                        "external_quality_metrics": "",
                        "software": f"{submitter_prefix}_SOFTWARE_CODEC_BAM_PROCESSING_PIPELINE_V1".upper()
                    })

    aligned_cols = [
        "submitted_id", "filename", "submitted_md5sum", "data_category", "data_type", "n50",
        "flow_cell_barcode", "flow_cell_lane", "description", "alignment_details", "file_format",
        "file_sets", "reference_genome", "derived_from", "external_quality_metrics", "software"
    ]
    
    aligned_df = pd.DataFrame(aligned_records)
    if not aligned_df.empty:
        aligned_df = aligned_df[aligned_cols]
    else:
        aligned_df = pd.DataFrame(columns=aligned_cols)
    
    aligned_df.to_csv(out_alignedreads_tsv, sep="\t", index=False)
    aligned_df.to_excel(out_alignedreads_xlsx, index=False)
    print(f"✅ AlignedReads sheet saved to: {out_alignedreads_tsv} and {out_alignedreads_xlsx}")

    # --- VARIANT CALLS SHEET (FIXED: use sample suffixes for consistent linking) ---
    variantcalls_records = []
    
    for codec_file in codec_files:
        df = pd.read_csv(codec_file, sep="\t")
        df = df.dropna(subset=["collaborator_sample_id"])
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)
        
        for _, row in df.iterrows():
            sample_id = row["collaborator_sample_id"]
            donor = sample_id.split("-")[0]
            core = sample_id.split("-")[1]
            sample_suffix = sample_id.split("-")[-1]
            tissue_label = normalize_tissue_label(tissue_map.get(core, core))
            
            for protocol in codec_protocols:
                if pd.notna(row.get("vcf")):
                    # FIXED: Use sample suffix to match other IDs
                    fileset_id = f"{submitter_prefix}_FILE-SET_CODEC_{donor}_{tissue_label}-{protocol}_{sample_suffix}".upper()
                    processed_id = f"{submitter_prefix}_ALIGNED-READS_PROCESSED_CODEC-{sample_id}-{protocol}_{sample_suffix}".upper()
                    variant_id = f"{submitter_prefix}_VARIANT-CALLS_CODEC-{sample_id}-{protocol}_{sample_suffix}".upper()
                    
                    variantcalls_records.append({
                        "submitted_id": variant_id,
                        "data_category": "Somatic Variant Calls",
                        "data_type": "SNV | Indel",
                        "filename": row["vcf"],
                        "submitted_md5sum": "",
                        "description": f"CODEC_{sample_id}-{protocol}",
                        "comparators": "",
                        "external_databases": "",
                        "filtering_methods": "",
                        "mode": "",
                        "file_format": "vcf",
                        "file_sets": fileset_id,
                        "reference_genome": "broad_grch38",
                        "derived_from": processed_id,
                        "external_quality_metrics": "",
                        "software": f"{submitter_prefix}_SOFTWARE_CODEC_VARIANT_CALL_PIPELINE_V1".upper()
                    })

    variantcalls_cols = [
        "submitted_id", "data_category", "data_type", "filename", "submitted_md5sum", "description",
        "comparators", "external_databases", "filtering_methods", "mode", "file_format", "file_sets",
        "reference_genome", "derived_from", "external_quality_metrics", "software"
    ]
    
    variantcalls_df = pd.DataFrame(variantcalls_records)
    if not variantcalls_df.empty:
        variantcalls_df = variantcalls_df[variantcalls_cols]
    else:
        variantcalls_df = pd.DataFrame(columns=variantcalls_cols)
    
    variantcalls_df.to_csv(out_variantcalls_tsv, sep="\t", index=False)
    variantcalls_df.to_excel(out_variantcalls_xlsx, index=False)
    print(f"✅ VariantCalls sheet saved to: {out_variantcalls_tsv} and {out_variantcalls_xlsx}")

    # --- SOFTWARE SHEET (CODEC-specific) ---
    software_records = [
        {
            "submitted_id": f"{submitter_prefix}_SOFTWARE_CODEC_ALIGNMENT_PIPELINE_V1".upper(),
            "category": "Alignment",
            "title": "CODEC alignment pipeline",
            "version": "1",
            "description": "",
            "binary_url": "",
            "commit": "",
            "source_url": "",
            "gpu": "",
            "model": "",
            "modification_tags": ""
        },
        {
            "submitted_id": f"{submitter_prefix}_SOFTWARE_CODEC_BAM_PROCESSING_PIPELINE_V1".upper(),
            "category": "Read Manipulation",
            "title": "CODEC bam process pipeline",
            "version": "1",
            "description": "",
            "binary_url": "",
            "commit": "",
            "source_url": "",
            "gpu": "",
            "model": "",
            "modification_tags": ""
        },
        {
            "submitted_id": f"{submitter_prefix}_SOFTWARE_CODEC_VARIANT_CALL_PIPELINE_V1".upper(),
            "category": "Variant Calling",
            "title": "CODEC variant calling pipeline",
            "version": "1",
            "description": "",
            "binary_url": "",
            "commit": "",
            "source_url": "",
            "gpu": "",
            "model": "",
            "modification_tags": ""
        }
    ]

    software_cols = [
        "submitted_id", "category", "title", "version", "description", "binary_url",
        "commit", "source_url", "gpu", "model", "modification_tags"
    ]
    
    software_df = pd.DataFrame(software_records)
    software_df = software_df[software_cols]
    software_df.to_csv(out_software_tsv, sep="\t", index=False)
    software_df.to_excel(out_software_xlsx, index=False)
    print(f"✅ Software sheet saved to: {out_software_tsv} and {out_software_xlsx}")

    # --- SEQUENCING SHEET (CODEC-specific) ---
    sequencing_records = [
        {
            "submitted_id": f"{submitter_prefix}_SEQUENCING_CODEC".upper(),
            "read_type": "Paired-end",
            "target_read_length": 150,
            "flow_cell": "",
            "movie_length": "",
            "on_target_rate": "",
            "target_coverage": 2,
            "target_read_count": "",
            "target_monomer_length": "",
            "sequencer": "illumina_novaseq_x_plus",
            "preparation_kits": ""
        }
    ]
    
    sequencing_cols = [
        "submitted_id", "read_type", "target_read_length", "flow_cell", "movie_length",
        "on_target_rate", "target_coverage", "target_read_count", "target_monomer_length",
        "sequencer", "preparation_kits"
    ]
    
    sequencing_df = pd.DataFrame(sequencing_records)
    sequencing_df = sequencing_df[sequencing_cols]
    
    # Fill empty strings for missing columns
    for col in sequencing_cols:
        if col not in sequencing_df.columns:
            sequencing_df[col] = ""
    
    sequencing_df.to_csv(out_sequencing_tsv, sep="\t", index=False)
    sequencing_df.to_excel(out_sequencing_xlsx, index=False)
    print(f"✅ Sequencing sheet saved to: {out_sequencing_tsv} and {out_sequencing_xlsx}")

    # --- LIBRARY PREPARATION SHEET (CODEC-specific: DDBTP only) ---
    library_preparation_records = [
        {
            "submitted_id": f"{submitter_prefix}_LIBRARY-PREPARATION_CODEC_DDBTP".upper(),
            "description": "ddBTP End Repair is an end repair strategy that utilizes ddBTPs to optimize DNA fragment ends for ligation while preventing unintended polymerase extension and exonuclease degradation.",
            "adapter_inclusion_method": "",
            "amplification_method": "",
            "fragmentation_method": "",
            "insert_selection_method": "",
            "enzymes": "",
            "rna_seq_protocol": "",
            "size_selection_method": "",
            "strand": "",
            "trim_adapter_sequence": "",
            "preparation_kits": "",
            "treatments": ""
        }
    ]
    
    library_preparation_cols = [
        "submitted_id", "description", "adapter_inclusion_method", "amplification_method",
        "fragmentation_method", "insert_selection_method", "enzymes", "rna_seq_protocol",
        "size_selection_method", "strand", "trim_adapter_sequence", "preparation_kits", "treatments"
    ]
    
    library_preparation_df = pd.DataFrame(library_preparation_records)
    library_preparation_df = library_preparation_df[library_preparation_cols]
    
    library_preparation_df.to_csv(out_library_preparation_tsv, sep="\t", index=False)
    library_preparation_df.to_excel(out_library_preparation_xlsx, index=False)
    print(f"✅ LibraryPreparation sheet saved to: {out_library_preparation_tsv} and {out_library_preparation_xlsx}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate metadata sheets for CODEC samples including Donor, Tissue, TissueSample, "
                    "Analyte, Library, FileSet, AlignedReads, VariantCalls, and Software."
    )
    parser.add_argument(
        "--donor-info", required=True,
        help="Path to the donor information TSV (with columns: donor, gender, age, etc.)"
    )
    parser.add_argument(
        "--codec", nargs="+", required=True,
        help="CODEC TSV files (with columns: sample_id, donor_id, collaborator_sample_id, MolConsensusBAM, RAW_BAM, vcf)"
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
        "--out-variantcalls-tsv", required=True,
        help="Where to write the VariantCalls sheet TSV"
    )
    parser.add_argument(
        "--out-variantcalls-xlsx", required=True,
        help="Where to write the VariantCalls sheet XLSX"
    )
    parser.add_argument(
        "--out-software-tsv", required=True,
        help="Where to write the Software sheet TSV"
    )
    parser.add_argument(
        "--out-software-xlsx", required=True,
        help="Where to write the Software sheet XLSX"
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
        "--out-library-preparation-tsv", required=True,
        help="Where to write the LibraryPreparation sheet TSV"
    )
    parser.add_argument(
        "--out-library-preparation-xlsx", required=True,
        help="Where to write the LibraryPreparation sheet XLSX"
    )
    parser.add_argument(
        "--out-combined-xlsx", required=True,
        help="Where to write the combined Excel workbook containing all sheets"
    )

    args = parser.parse_args()
    
    generate_codec_metadata_sheets(
        args.donor_info, args.codec, args.uberon,
        args.out_donor_tsv, args.out_donor_xlsx,
        args.out_tissue_tsv, args.out_tissue_xlsx,
        args.out_tissuesample_tsv, args.out_tissuesample_xlsx,
        args.submitter_prefix,
        args.out_analyte_tsv, args.out_analyte_xlsx,
        args.out_library_tsv, args.out_library_xlsx,
        args.out_fileset_tsv, args.out_fileset_xlsx,
        args.out_unalignedreads_tsv, args.out_unalignedreads_xlsx,
        args.out_alignedreads_tsv, args.out_alignedreads_xlsx,
        args.out_variantcalls_tsv, args.out_variantcalls_xlsx,
        args.out_software_tsv, args.out_software_xlsx,
        args.out_sequencing_tsv, args.out_sequencing_xlsx,
        args.out_library_preparation_tsv, args.out_library_preparation_xlsx
    )
    
    # Create combined Excel workbook
    with pd.ExcelWriter(args.out_combined_xlsx, engine="openpyxl") as writer:
        df_donor = pd.read_excel(args.out_donor_xlsx, dtype={"tpc_submitted": str})
        df_donor["tpc_submitted"] = df_donor["tpc_submitted"].str.capitalize()
        df_donor.to_excel(writer, sheet_name="(Donor)", index=False)
        pd.read_excel(args.out_tissue_xlsx).to_excel(writer, sheet_name="(Tissue)", index=False)
        pd.read_excel(args.out_tissuesample_xlsx).to_excel(writer, sheet_name="TissueSample", index=False)
        pd.read_excel(args.out_analyte_xlsx).to_excel(writer, sheet_name="Analyte", index=False)
        pd.read_excel(args.out_library_xlsx).to_excel(writer, sheet_name="Library", index=False)
        pd.read_excel(args.out_library_preparation_xlsx).to_excel(writer, sheet_name="LibraryPreparation", index=False)
        pd.read_excel(args.out_sequencing_xlsx).to_excel(writer, sheet_name="Sequencing", index=False)
        pd.read_excel(args.out_fileset_xlsx).to_excel(writer, sheet_name="FileSet", index=False)
        pd.read_excel(args.out_unalignedreads_xlsx).to_excel(writer, sheet_name="UnalignedReads", index=False)
        pd.read_excel(args.out_alignedreads_xlsx).to_excel(writer, sheet_name="AlignedReads", index=False)
        pd.read_excel(args.out_variantcalls_xlsx).to_excel(writer, sheet_name="VariantCalls", index=False)
        pd.read_excel(args.out_software_xlsx).to_excel(writer, sheet_name="Software", index=False)
    print(f"✅ The combined CODEC metadata Excel is saved to: {args.out_combined_xlsx}")


if __name__ == "__main__":
    main()
