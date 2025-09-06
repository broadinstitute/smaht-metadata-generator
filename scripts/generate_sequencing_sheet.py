import pandas as pd

def generate_sequencing_sheet(args):
    """Generate Sequencing sheet based on input file types provided."""
    
    # Define all possible columns
    sequencing_columns = [
        "submitted_id", "read_type", "target_read_length", "flow_cell", 
        "movie_length", "on_target_rate", "target_coverage", "target_read_count", 
        "target_monomer_length", "sequencer", "preparation_kits"
    ]
    
    sequencing_records = []
    seen_sequencing_ids = set()
    
    # Check for Long Read DNA files
    if getattr(args, 'lr_dna', None) and args.lr_dna:
        sequencing_id = "BROAD_SEQUENCING_PACBIO_15000BP_20X"
        if sequencing_id not in seen_sequencing_ids:
            sequencing_records.append({
                "submitted_id": sequencing_id,
                "read_type": "Single-end",
                "target_read_length": 15000,
                "flow_cell": "",
                "movie_length": "",
                "on_target_rate": "",
                "target_coverage": 20,
                "target_read_count": "",
                "target_monomer_length": "",
                "sequencer": "pacbio_revio_hifi",
                "preparation_kits": ""
            })
            seen_sequencing_ids.add(sequencing_id)
    
    # Check for Short Read DNA files  
    if getattr(args, 'sr_dna', None) and args.sr_dna:
        sequencing_id = "BROAD_SEQUENCING_NOVASEQXPLUS_25B_150BP_160X"
        if sequencing_id not in seen_sequencing_ids:
            sequencing_records.append({
                "submitted_id": sequencing_id,
                "read_type": "Paired-end",
                "target_read_length": 150,
                "flow_cell": "25B",
                "movie_length": "",
                "on_target_rate": "",
                "target_coverage": 160,
                "target_read_count": "",
                "target_monomer_length": "",
                "sequencer": "illumina_novaseq_x_plus",
                "preparation_kits": ""
            })
            seen_sequencing_ids.add(sequencing_id)
    
    # Check for RNA Watchmaker files
    if getattr(args, 'rna', None) and args.rna:
        for rna_file in args.rna:
            if "watchmaker" in rna_file.lower():
                sequencing_id = "BROAD_SEQUENCING_NOVASEQX_25B_145BP_100M"
                if sequencing_id not in seen_sequencing_ids:
                    sequencing_records.append({
                        "submitted_id": sequencing_id,
                        "read_type": "Paired-end", 
                        "target_read_length": 145,
                        "flow_cell": "25B",
                        "movie_length": "",
                        "on_target_rate": "",
                        "target_coverage": "",
                        "target_read_count": 100000000,
                        "target_monomer_length": "",
                        "sequencer": "illumina_novaseq_x_plus",
                        "preparation_kits": ""
                    })
                    seen_sequencing_ids.add(sequencing_id)
    
    # Create DataFrame
    if sequencing_records:
        sequencing_df = pd.DataFrame(sequencing_records)
        print(f"📊 Generated {len(sequencing_records)} Sequencing records")
    else:
        print("⚠ No sequencing records to generate - creating empty Sequencing sheet")
        sequencing_df = pd.DataFrame(columns=sequencing_columns)
    
    # Ensure all columns are present and in correct order
    sequencing_df = sequencing_df.reindex(columns=sequencing_columns, fill_value="")
    
    # Save to files
    sequencing_df.to_csv(args.out_sequencing_tsv, sep="\t", index=False)
    sequencing_df.to_excel(args.out_sequencing_xlsx, index=False) 
    print(f"✅ Sequencing sheet saved to: {args.out_sequencing_tsv} and {args.out_sequencing_xlsx}")
