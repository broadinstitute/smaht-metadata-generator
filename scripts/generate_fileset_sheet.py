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


def get_library_type(identifier):
    """Return standardized library type based on source_type or analyte_id patterns."""
    if identifier in {"WGS_ILLUMINA", "WGS_PACBIO", "RNA_WATCHMAKER", "RNA_TRUSEQ"}:
        return identifier
    # fallback: infer from analyte_id-like string
    id_str = identifier.upper()
    if "SHORTREAD" in id_str:
        return "WGS_ILLUMINA"
    elif "LONGREAD" in id_str:
        return "WGS_PACBIO"
    elif "WATCHMAKER" in id_str:
        return "RNA_WATCHMAKER"
    elif "TRUSEQ" in id_str:
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


def generate_fileset_sheet(args):
    # Load UBERON tissue map
    uberon_df = pd.read_csv(args.uberon, sep="\t")
    tissue_map = dict(zip(uberon_df["tissue_identifier_code"], uberon_df["corresponding_tissue"]))

    # Build lookup for source types from ALL input files
    all_input = []
    
    # Process DNA files (both short-read and long-read)
    for f in args.inputs:
        print(f"📁 Processing DNA file: {f}")
        df = pd.read_csv(f, sep="\t")
        df["donor_id"] = df["donor_id"].astype(str)
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)
        df["sample_id"] = df.get("sample_id", df["collaborator_sample_id"])
        
        # IMPROVED: Better file type detection
        f_lower = f.lower()
        if "short" in f_lower or "short_read" in f_lower:
            df["source_type"] = "WGS_ILLUMINA"
        elif "long" in f_lower or "long_read" in f_lower:
            df["source_type"] = "WGS_PACBIO"
        else:
            # Default to short-read if not specified
            df["source_type"] = "WGS_ILLUMINA"
            print(f"⚠️  File type not clear from filename '{f}', defaulting to WGS_ILLUMINA")
        
        all_input.append(df[["donor_id", "collaborator_sample_id", "sample_id", "source_type"]])

    # Process RNA files if provided
    if getattr(args, "rna", None):
        print(f"📁 Processing RNA file: {args.rna}")
        rna_df = pd.read_csv(args.rna, sep="\t")
        rna_df["donor_id"] = rna_df["donor_id"].astype(str)
        rna_df["collaborator_sample_id"] = rna_df["collaborator_sample_id"].astype(str)
        rna_df["sample_id"] = rna_df.get("sample_id", rna_df["collaborator_sample_id"])
        
        rna_lower = args.rna.lower()
        if "watchmaker" in rna_lower:
            rna_df["source_type"] = "RNA_WATCHMAKER"
        elif "truseq" in rna_lower:
            rna_df["source_type"] = "RNA_TRUSEQ"
        else:
            rna_df["source_type"] = "RNA_WATCHMAKER"  # Default to WATCHMAKER
            print(f"⚠️  RNA type not clear from filename '{args.rna}', defaulting to RNA_WATCHMAKER")
        
        all_input.append(rna_df[["donor_id", "collaborator_sample_id", "sample_id", "source_type"]])

    if not all_input:
        print("❌ No input files to process for FileSet sheet")
        return

    lookup = pd.concat(all_input, ignore_index=True)
    print(f"📊 Combined lookup table has {len(lookup)} entries from {len(all_input)} files")

    # Generate FileSet entries directly from sample lookup
    sample_seen = {}
    records = []
    
    # Loop over each sample from long-read, short-read, and RNA tables
    for _, sample_row in lookup.iterrows():
        donor = sample_row["donor_id"]
        collaborator_id = sample_row["collaborator_sample_id"]
        description = sample_row["sample_id"]
        source_type = sample_row["source_type"]

        # derive core and tissue_label just like UnalignedReads
        core = collaborator_id.split("-")[1]
        tissue_label = normalize_tissue_label(tissue_map.get(core, core))

        # determine library type (WGS_ILLUMINA, WGS_PACBIO, etc.)
        lib_type = get_library_type(source_type)

        # increment counter per donor-tissue-library combo
        key = f"{donor}_{tissue_label}_{lib_type}"
        sample_seen[key] = sample_seen.get(key, 0) + 1

        # build submitted_id
        submitted_id = (
            f"{args.submitter_prefix}_FILE-SET_{donor}_{tissue_label}_{lib_type}_{sample_seen[key]}"
            .upper()
        )
        sequencing = infer_sequencing(lib_type)

        records.append({
            "submitted_id": submitted_id,
            "description": description,
            "libraries": "",
            "sequencing": sequencing
        })

    print(f"📊 Generated {len(records)} FileSet records")

    # assemble and write outputs
    out_df = pd.DataFrame(records)
    library_df = pd.read_csv(args.out_library_tsv, sep="\t")
    library_ids = set(library_df["submitted_id"])

    def match_libraries(fileset_id):
        parts = fileset_id.split("_")
        donor = parts[2]
        tissue_label = "_".join(parts[3:-2])  # handles multi-word tissue names
        library_type = parts[-2]
        replicate_suffix = parts[-1]

        lib_id = f"{args.submitter_prefix}_LIBRARY_{donor}_{tissue_label}_{library_type}_{replicate_suffix}".upper()
        return lib_id if lib_id in library_ids else ""

    out_df["libraries"] = out_df["submitted_id"].apply(match_libraries)
    
    # Ensure all required columns exist
    for c in ["submitted_id", "description", "libraries", "sequencing", "samples"]:
        if c not in out_df.columns:
            out_df[c] = ""
    out_df = out_df[["submitted_id", "description", "libraries", "sequencing", "samples"]]
    
    out_df.to_csv(args.out_fileset_tsv, sep="\t", index=False)
    out_df.to_excel(args.out_fileset_xlsx, index=False)
    print(f"✅ FileSet sheet saved to: {args.out_fileset_tsv} and {args.out_fileset_xlsx}")