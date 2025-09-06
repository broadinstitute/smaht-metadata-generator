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


def generate_alignedreads_sheet(args):
    import pandas as pd

    # Define required columns first
    cols = [
        "submitted_id","filename","submitted_md5sum","data_category","data_type","n50",
        "flow_cell_barcode","flow_cell_lane","description","alignment_details",
        "file_format","file_sets","reference_genome","derived_from",
        "external_quality_metrics","software"
    ]

    uberon_df = pd.read_csv(args.uberon, sep="\t")
    tissue_map = dict(zip(uberon_df["tissue_identifier_code"], uberon_df["corresponding_tissue"]))

    file_labels = {
        f: ("WGS_ILLUMINA" if "short" in f.lower() else "WGS_PACBIO")
        for f in args.inputs
        if "long" not in f.lower()
    }
    
    if getattr(args, "rna", None):
        file_labels[args.rna] = (
            "RNA_TRUSEQ" if "truseq" in args.rna.lower() else "RNA_WATCHMAKER"
        )

    # Check if we have any files with CRAM data
    files_with_cram = []
    for f in file_labels.keys():
        try:
            df = pd.read_csv(f, sep="\t")
            if "cram_path" in df.columns and not df["cram_path"].dropna().empty:
                files_with_cram.append(f)
        except:
            continue

    if not files_with_cram:
        print("⚠ No files with CRAM data found - creating empty AlignedReads sheet")
        # CREATE EMPTY FILES
        empty_df = pd.DataFrame(columns=cols)
        empty_df.to_csv(args.out_alignedreads_tsv, sep="\t", index=False)
        empty_df.to_excel(args.out_alignedreads_xlsx, index=False)
        print(f"✅ Empty AlignedReads sheet saved to: {args.out_alignedreads_tsv} and {args.out_alignedreads_xlsx}")
        return

    fileset_df = pd.read_csv(args.out_fileset_tsv, sep="\t")

    aligned_records = []
    # CHANGE: Track counts per (donor, tissue_label, lib_type) instead of (sample_id, lib_type)
    align_seen = {}
    fileset_ids = {}

    for f, lib_type in file_labels.items():
        if f not in files_with_cram:
            continue
            
        df = pd.read_csv(f, sep="\t").dropna(subset=["collaborator_sample_id", "cram_path"])
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)

        for _, row in df.iterrows():
            sample_id = row["collaborator_sample_id"]
            donor, core = sample_id.split("-")[:2]
            tissue_label = normalize_tissue_label(tissue_map.get(core, core))

            # CHANGE: count repeats per donor-tissue-library combination
            key = (donor, tissue_label, lib_type)
            align_seen[key] = align_seen.get(key, 0) + 1
            count = align_seen[key]

            # build seg part
            if lib_type.startswith("RNA"):
                seg = f"{lib_type}-CRAM"
            else:
                seg = f"{lib_type}-CRAM" if count == 1 else f"{lib_type}-{count}-CRAM"

            submitted_id = (
                f"{args.submitter_prefix}_ALIGNED-READS_{donor}-{tissue_label}-{seg}"
            ).upper()

            # filename and file_format
            filename = row["cram_path"].split("/")[-1]
            file_format = filename.split('.')[-1]

            fs_key = (donor, tissue_label, lib_type, count)
            if fs_key not in fileset_ids:
                fs_id = f"{args.submitter_prefix}_FILE-SET_{donor}_{tissue_label}_{lib_type}_{count}".upper()
                match = fileset_df[fileset_df["submitted_id"] == fs_id]
                fileset_ids[fs_key] = fs_id if not match.empty else ""
            file_set_id = fileset_ids[fs_key]

            aligned_records.append({
                "submitted_id": submitted_id,
                "filename": filename,
                "submitted_md5sum": "",
                "data_category": "",
                "data_type": "",
                "n50": "",
                "flow_cell_barcode": "",
                "flow_cell_lane": "",
                "description": sample_id,
                "alignment_details": "Sorted",
                "file_format": file_format,
                "file_sets": file_set_id,
                "reference_genome": "broad_grch38",
                "derived_from": "",
                "external_quality_metrics": "",
                "software": "BROAD_SOFTWARE_DRAGEN_07.021.604.3.7.8"
            })

    aligned_df = pd.DataFrame(aligned_records)
    aligned_df = aligned_df[cols]

    aligned_df.to_csv(args.out_alignedreads_tsv, sep="\t", index=False)
    aligned_df.to_excel(args.out_alignedreads_xlsx, index=False)
    print(f"✅ AlignedReads sheet saved to: {args.out_alignedreads_tsv} and {args.out_alignedreads_xlsx}")