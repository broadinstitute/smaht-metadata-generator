# python
import click
import sys
from pathlib import Path
import pandas as pd

# Add the scripts directory to the Python path
scripts_dir = Path(__file__).parents[2] / "scripts"
sys.path.append(str(scripts_dir))

from .generate_metadata import generate_metadata_sheets
from .generate_fileset_sheet import generate_fileset_sheet
from .generate_unaligned_reads_sheet import generate_unalignedreads_sheet
from .generate_aligned_reads_sheet import generate_alignedreads_sheet

@click.group()
def cli():
    """metagen - Tool for generating SMaHT metadata sheets."""
    pass

@cli.command()
@click.option('--donor-info', required=True, type=click.Path(exists=True), help='Path to donor info TSV')
@click.option('--inputs', required=True, multiple=True, type=click.Path(exists=True), help='Paths to input sample TSVs')
@click.option('--rna', required=True, type=click.Path(exists=True), help='Path to RNA input TSV')
@click.option('--uberon', required=True, type=click.Path(exists=True), help='Path to UBERON mapping TSV')
@click.option('--submitter-prefix', required=True, help='Prefix for submitted IDs')
@click.option('--output-dir', required=True, type=click.Path(), help='Directory for output files')
@click.option('--combined-filename', default='all_metadata.xlsx', help='Filename for combined Excel workbook')
def metadata(donor_info, inputs, rna, uberon, submitter_prefix, output_dir, combined_filename):
    """Generate all metadata sheets for SMaHT."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Define output paths
    out_paths = {
        'donor_tsv': str(output_dir / "donor.tsv"),
        'donor_xlsx': str(output_dir / "donor.xlsx"),
        'tissue_tsv': str(output_dir / "tissue.tsv"),
        'tissue_xlsx': str(output_dir / "tissue.xlsx"),
        'tissuesample_tsv': str(output_dir / "tissuesample.tsv"),
        'tissuesample_xlsx': str(output_dir / "tissuesample.xlsx"),
        'analyte_tsv': str(output_dir / "analyte.tsv"),
        'analyte_xlsx': str(output_dir / "analyte.xlsx"),
        'library_tsv': str(output_dir / "library.tsv"),
        'library_xlsx': str(output_dir / "library.xlsx"),
        'fileset_tsv': str(output_dir / "fileset.tsv"),
        'fileset_xlsx': str(output_dir / "fileset.xlsx"),
        'unalignedreads_tsv': str(output_dir / "unalignedreads.tsv"),
        'unalignedreads_xlsx': str(output_dir / "unalignedreads.xlsx"),
        'alignedreads_tsv': str(output_dir / "alignedreads.tsv"),
        'alignedreads_xlsx': str(output_dir / "alignedreads.xlsx"),
        'combined_xlsx': str(output_dir / combined_filename)
    }

    click.echo("Generating metadata sheets...")

    # Call the main metadata generation function
    generate_metadata_sheets(
        donor_info, list(inputs), rna, uberon,
        out_paths['donor_tsv'], out_paths['donor_xlsx'],
        out_paths['tissue_tsv'], out_paths['tissue_xlsx'],
        out_paths['tissuesample_tsv'], out_paths['tissuesample_xlsx'],
        submitter_prefix,
        out_paths['analyte_tsv'], out_paths['analyte_xlsx'],
        out_paths['library_tsv'], out_paths['library_xlsx'],
        out_paths['fileset_tsv'], out_paths['fileset_xlsx'],
        out_paths['unalignedreads_tsv'], out_paths['unalignedreads_xlsx'],
        out_paths['alignedreads_tsv'], out_paths['alignedreads_xlsx']
    )

    # Generate additional sheets
    generate_fileset_sheet(vars(type("Args", (), {
        "inputs": list(inputs),
        "rna": rna,
        "uberon": uberon,
        "submitter_prefix": submitter_prefix,
        "out_fileset_tsv": out_paths['fileset_tsv'],
        "out_fileset_xlsx": out_paths['fileset_xlsx']
    })()))

    generate_unalignedreads_sheet(vars(type("Args", (), {
        "donor_info": donor_info,
        "inputs": list(inputs),
        "rna": rna,
        "uberon": uberon,
        "submitter_prefix": submitter_prefix,
        "out_unalignedreads_tsv": out_paths['unalignedreads_tsv'],
        "out_unalignedreads_xlsx": out_paths['unalignedreads_xlsx']
    })()))

    generate_alignedreads_sheet(vars(type("Args", (), {
        "donor_info": donor_info,
        "inputs": list(inputs),
        "rna": rna,
        "uberon": uberon,
        "submitter_prefix": submitter_prefix,
        "out_alignedreads_tsv": out_paths['alignedreads_tsv'],
        "out_alignedreads_xlsx": out_paths['alignedreads_xlsx']
    })()))

    # Create combined Excel workbook
    with pd.ExcelWriter(out_paths['combined_xlsx'], engine="openpyxl") as writer:
        pd.read_excel(out_paths['donor_xlsx']).to_excel(writer, sheet_name="Donor", index=False)
        pd.read_excel(out_paths['tissue_xlsx']).to_excel(writer, sheet_name="Tissue", index=False)
        pd.read_excel(out_paths['tissuesample_xlsx']).to_excel(writer, sheet_name="TissueSample", index=False)
        pd.read_excel(out_paths['analyte_xlsx']).to_excel(writer, sheet_name="Analyte", index=False)
        pd.read_excel(out_paths['library_xlsx']).to_excel(writer, sheet_name="Library", index=False)
        pd.read_excel(out_paths['fileset_xlsx']).to_excel(writer, sheet_name="FileSet", index=False)
        pd.read_excel(out_paths['unalignedreads_xlsx']).to_excel(writer, sheet_name="UnalignedReads", index=False)
        pd.read_excel(out_paths['alignedreads_xlsx']).to_excel(writer, sheet_name="AlignedReads", index=False)

    click.echo(f"All metadata sheets generated. Combined workbook saved to {out_paths['combined_xlsx']}.")


# python
@cli.command()
@click.option('--inputs', required=True, multiple=True, type=click.Path(exists=True), help='Paths to input sample TSVs')
@click.option('--rna', required=True, type=click.Path(exists=True), help='Path to RNA input TSV')
@click.option('--uberon', required=True, type=click.Path(exists=True), help='Path to UBERON mapping TSV')
@click.option('--submitter-prefix', required=True, help='Prefix for submitted IDs')
@click.option('--output-dir', required=True, type=click.Path(), help='Directory for output files')
def fileset(inputs, rna, uberon, submitter_prefix, output_dir):
    """
    Generate only the FileSet metadata sheet.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    out_fileset_tsv = str(output_dir / "fileset.tsv")
    out_fileset_xlsx = str(output_dir / "fileset.xlsx")

    click.echo("Generating FileSet metadata sheet...")

    fs_args = type("Args", (), {
        "inputs": list(inputs),
        "rna": rna,
        "uberon": uberon,
        "submitter_prefix": submitter_prefix,
        "out_fileset_tsv": out_fileset_tsv,
        "out_fileset_xlsx": out_fileset_xlsx
    })()

    generate_fileset_sheet(fs_args)

    click.echo(f"FileSet sheet saved to {out_fileset_tsv} and {out_fileset_xlsx}")

# python
@cli.command()
@click.option('--inputs', required=True, multiple=True, type=click.Path(exists=True), help='Paths to input sample TSVs')
@click.option('--uberon', required=True, type=click.Path(exists=True), help='Path to UBERON mapping TSV')
@click.option('--submitter-prefix', required=True, help='Prefix for submitted IDs')
@click.option('--output-dir', required=True, type=click.Path(), help='Directory for output files')
def unalignedreads(inputs, uberon, submitter_prefix, output_dir):
    """
    Generate only the UnalignedReads metadata sheet.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    out_unalignedreads_tsv = str(output_dir / "unalignedreads.tsv")
    out_unalignedreads_xlsx = str(output_dir / "unalignedreads.xlsx")

    click.echo("Generating UnalignedReads metadata sheet...")

    ur_args = type("Args", (), {
        "inputs": list(inputs),
        "uberon": uberon,
        "submitter_prefix": submitter_prefix,
        "out_unalignedreads_tsv": out_unalignedreads_tsv,
        "out_unalignedreads_xlsx": out_unalignedreads_xlsx
    })()

    generate_unalignedreads_sheet(ur_args)

    click.echo(f"UnalignedReads sheet saved to {out_unalignedreads_tsv} and {out_unalignedreads_xlsx}")

@cli.command()
@click.option('--inputs', required=True, multiple=True, type=click.Path(exists=True), help='Paths to input sample TSVs')
@click.option('--rna', required=True, type=click.Path(exists=True), help='Path to RNA input TSV')
@click.option('--uberon', required=True, type=click.Path(exists=True), help='Path to UBERON mapping TSV')
@click.option('--submitter-prefix', required=True, help='Prefix for submitted IDs')
@click.option('--output-dir', required=True, type=click.Path(), help='Directory for output files')
def alignedreads(inputs, rna, uberon, submitter_prefix, output_dir):
    """
    Generate only the AlignedReads metadata sheet.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    out_alignedreads_tsv = str(output_dir / "alignedreads.tsv")
    out_alignedreads_xlsx = str(output_dir / "alignedreads.xlsx")

    click.echo("Generating AlignedReads metadata sheet...")

    ar_args = type("Args", (), {
        "inputs": list(inputs),
        "rna": rna,
        "uberon": uberon,
        "submitter_prefix": submitter_prefix,
        "out_alignedreads_tsv": out_alignedreads_tsv,
        "out_alignedreads_xlsx": out_alignedreads_xlsx
    })()

    generate_alignedreads_sheet(ar_args)

    click.echo(f"AlignedReads sheet saved to {out_alignedreads_tsv} and {out_alignedreads_xlsx}")

if __name__ == '__main__':
    cli()