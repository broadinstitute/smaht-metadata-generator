
#!/usr/bin/env python3

import csv
import argparse

def strip_paths(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")

        writer.writeheader()
        for row in reader:
            if "cram_path" in row and row["cram_path"]:
                row["cram_path"] = row["cram_path"].split("/")[-1]
            if "input_bam" in row and row["input_bam"]:
                row["input_bam"] = row["input_bam"].split("/")[-1]
            writer.writerow(row)

    print(f"✅ Output written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Strip full paths in 'cram_path' and/or 'input_bam' columns and keep only filenames."
    )
    parser.add_argument("--input", required=True, help="Path to the input TSV file")
    parser.add_argument("--output", required=True, help="Path to save the output TSV with stripped paths")

    args = parser.parse_args()
    strip_paths(args.input, args.output)

if __name__ == "__main__":
    main()
