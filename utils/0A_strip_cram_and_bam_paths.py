#!/usr/bin/env python3

import csv
import argparse

def strip_paths(input_file, output_file):
    """
    Strip full paths from file columns and keep only filenames.
    
    Handles columns for both regular sequencing files and CODEC files:
    - Regular: cram_path, input_bam
    - CODEC: MolConsensusBAM, RAW_BAM, vcf
    """
    
    # Define columns that may contain file paths
    path_columns = [
        "cram_path",        # Regular sequencing
        "input_bam",        # Regular sequencing  
        "MolConsensusBAM",  # CODEC
        "RAW_BAM",          # CODEC
        "vcf"               # CODEC
    ]
    
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")

        writer.writeheader()
        
        files_processed = 0
        columns_found = []
        
        for row in reader:
            files_processed += 1
            
            # Process each potential path column
            for col in path_columns:
                if col in row and row[col]:
                    # Track which columns we found for reporting
                    if col not in columns_found:
                        columns_found.append(col)
                    
                    # Strip path - keep only filename
                    # Handle both Unix/Linux paths (/) and Windows paths (\)
                    if "/" in row[col]:
                        row[col] = row[col].split("/")[-1]
                    elif "\\" in row[col]:
                        row[col] = row[col].split("\\")[-1]
            
            writer.writerow(row)

    print(f"Processed {files_processed} rows")
    print(f"Found and processed these path columns: {', '.join(columns_found) if columns_found else 'None'}")
    print(f"Output written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Strip full paths from file columns and keep only filenames. "
                    "Works with both regular sequencing files (cram_path, input_bam) "
                    "and CODEC files (MolConsensusBAM, RAW_BAM, vcf)."
    )
    parser.add_argument("--input", required=True, help="Path to the input TSV file")
    parser.add_argument("--output", required=True, help="Path to save the output TSV with stripped paths")

    args = parser.parse_args()
    strip_paths(args.input, args.output)

if __name__ == "__main__":
    main()