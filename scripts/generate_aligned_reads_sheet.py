def generate_alignedreads_sheet(args):
    import pandas as pd

    uberon_df = pd.read_csv(args.uberon, sep="\t")
    tissue_map = dict(zip(uberon_df["tissue_identifier_code"], uberon_df["corresponding_tissue"]))

    file_labels = {
        f: ("WGS_ILLUMINA" if "short" in f.lower() else "WGS_PACBIO")
        for f in args.inputs
        if "long" not in f.lower()
    }
    file_labels[args.rna] = ("RNA_TRUSEQ" if "truseq" in args.rna.lower() else "RNA_WATCHMAKER")

    fileset_df = pd.read_csv(args.out_fileset_tsv, sep="\t")

    aligned_records = []
    align_seen = {}
    fileset_ids = {}

    for f, lib_type in file_labels.items():
        df = pd.read_csv(f, sep="\t").dropna(subset=["collaborator_sample_id", "cram_path"])
        df["collaborator_sample_id"] = df["collaborator_sample_id"].astype(str)

        for _, row in df.iterrows():
            sample_id = row["collaborator_sample_id"]
            donor, core = sample_id.split("-")[:2]
            tissue_label = tissue_map.get(core, core).replace(" ", "-").upper()

            # count repeats
            key = (sample_id, lib_type)
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

            # lookup or cache FileSet ID
            fs_key = (donor, tissue_label, lib_type)
            if fs_key not in fileset_ids:
                # find matching FileSet
                mask = fileset_df["submitted_id"].str.contains(
                    f"{donor}_{tissue_label}_{lib_type}", na=False
                )
                fileset_ids[fs_key] = fileset_df.loc[mask, "submitted_id"].iat[0] if mask.any() else ""
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
                "software": "STAR"
            })

    aligned_df = pd.DataFrame(aligned_records)
    cols = [
        "submitted_id","filename","submitted_md5sum","data_category","data_type","n50",
        "flow_cell_barcode","flow_cell_lane","description","alignment_details",
        "file_format","file_sets","reference_genome","derived_from",
        "external_quality_metrics","software"
    ]
    aligned_df = aligned_df[cols]

    aligned_df.to_csv(args.out_alignedreads_tsv, sep="\t", index=False)
    aligned_df.to_excel(args.out_alignedreads_xlsx, index=False)
    print(f"✅ AlignedReads sheet saved to: {args.out_alignedreads_tsv} and {args.out_alignedreads_xlsx}")

