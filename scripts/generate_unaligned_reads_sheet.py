import argparse
import pandas as pd

# --- HELPER FUNCTIONS ---
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


def get_library_type(analyte_id):
    if "SHORTREAD" in analyte_id:
        return "WGS_ILLUMINA"
    elif "LONGREAD" in analyte_id:
        return "WGS_PACBIO"
    elif "WATCHMAKER" in analyte_id:
        return "RNA_WATCHMAKER"
    elif "TRUSEQ" in analyte_id:
        return "RNA_TRUSEQ"
    return "UNKNOWN"

def infer_sequencing(library_type):
    if "PACBIO" in library_type:
        return "BROAD_SEQUENCING_PACBIO_20000BP_12X"
    elif "WATCHMAKER" in library_type:
        return "BROAD_SEQUENCING_NOVASEQX_10B_150BP_100M"
    elif "TRUSEQ" in library_type:
        return "BROAD_SEQUENCING_NOVASEQ6000_150BP_75M"
    else:
        return "BROAD_SEQUENCING_NOVASEQX_25B_150BP_160X"

def get_unalignedread_type():
    return "WGS_PACBIO_SMRT"

def generate_unalignedreads_sheet(args):
    # Define required columns first
    required_cols = [
        "submitted_id", "filename", "submitted_md5sum", "data_category", "data_type", "n50",
        "flow_cell_barcode", "flow_cell_lane", "description", "read_pair_number", "file_format",
        "file_sets", "derived_from", "external_quality_metrics", "paired_with", "software"
    ]
    
    longread_file = next((f for f in args.inputs if "long" in f.lower()), None)
    if not longread_file:
        print("❌ No long read file found in inputs - creating empty UnalignedReads sheet")
        # CREATE EMPTY FILES
        empty_df = pd.DataFrame(columns=required_cols)
        empty_df.to_csv(args.out_unalignedreads_tsv, sep="\t", index=False)
        empty_df.to_excel(args.out_unalignedreads_xlsx, index=False)
        print(f"✅ Empty UnalignedReads sheet saved to: {args.out_unalignedreads_tsv} and {args.out_unalignedreads_xlsx}")
        return

    long_df = pd.read_csv(longread_file, sep="\t")
    long_df = long_df.dropna(subset=["collaborator_sample_id", "input_bam"])

    long_df["collaborator_sample_id"] = long_df["collaborator_sample_id"].astype(str)
    long_df["input_bam"] = long_df["input_bam"].astype(str)

    # UBERON info
    uberon_df = pd.read_csv(args.uberon, sep="\t")
    tissue_map = dict(zip(uberon_df["tissue_identifier_code"], uberon_df["corresponding_tissue"]))

    sample_seen = {}
    unaligned_records = []

    for _, row in long_df.iterrows():
        sample_id = row["collaborator_sample_id"]
        donor = sample_id.split("-")[0]
        core_id = sample_id.split("-")[1]
        tissue_label = normalize_tissue_label(tissue_map.get(core_id, core_id))

        sample_seen[sample_id] = sample_seen.get(sample_id, 0) + 1
        smrt_suffix = f"SMRT{sample_seen[sample_id]}"

        submitted_id = f"{args.submitter_prefix}_UNALIGNED-READS_{donor}_{tissue_label}_{get_unalignedread_type()}_{smrt_suffix}".upper()

        file_set_id = f"{args.submitter_prefix}_FILE-SET_{donor}_{tissue_label}_WGS_PACBIO_{sample_seen[sample_id]}".upper()

        unaligned_records.append({
            "submitted_id": submitted_id,
            "filename": row["input_bam"],
            "description": sample_id,
            "file_format": "bam",
            "file_sets": file_set_id
        })

    unaligned_df = pd.DataFrame(unaligned_records)
    
    # Add missing columns
    for col in required_cols:
        if col not in unaligned_df.columns:
            unaligned_df[col] = ""

    unaligned_df = unaligned_df[required_cols]
    unaligned_df.to_csv(args.out_unalignedreads_tsv, sep="\t", index=False)
    unaligned_df.to_excel(args.out_unalignedreads_xlsx, index=False)
    print(f"✅ UnalignedReads sheet saved to: {args.out_unalignedreads_tsv} and {args.out_unalignedreads_xlsx}")