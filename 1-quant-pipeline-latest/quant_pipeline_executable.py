#!/usr/bin/env python
# coding: utf-8
"""
Run the four quantification steps in order, passing each step's output to the next.

Equivalent to running the four scripts one after another:

    Map_and_count_peptides.py -> Relative_counts.py
                              -> Generate_PTM_matrices.py
                              -> Combine_PTM_matrices.py

Example
-------
    python quant_pipeline_executable.py \\
        ./v0113-2026/unprocessed/20260220_Peptidoforms_IDs_v0113.csv.gz \\
        ./v0113-2026/unprocessed/20260220_Peptidoforms_counts_v0113_breast_cancer.csv.gz \\
        ./Human_2026_01_canonical.fasta.gz \\
        -p Validation-20260812 \\
        -m ./PTMs-of-interest-submission.csv \\
        -c ./MQcontaminants_2023_11_14.fasta.gz
"""
import argparse, os
from datetime import date

from Map_and_count_peptides import get_peptidoforms_counts_mapped
from Relative_counts import get_relative_ptms
from Generate_PTM_matrices import generate_ptm_submatrices
from Combine_PTM_matrices import combine_ptm_submatrices


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("Example")[0].strip(),
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('peptidoform_ids', type=existing_file, help="Path to 'Peptidoforms_IDs' file")
    p.add_argument('peptidoform_counts', type=existing_file, help="Path to 'peptidoform_counts' file")
    p.add_argument('fasta', type=existing_file, help='Path to FASTA file. Used both to map peptidoforms to proteins and to prefilter the PTM matrices')
    p.add_argument("-p", "--prefix", type=str, default=date.today().isoformat(), help="Prefix for output files (Default: today's date)")
    p.add_argument("-t", "--threads", type=int, default=-1, help="Number of threads to use (Default: all available)")
    p.add_argument('-m', '--myptms', dest='myptms', type=existing_file, default='./PTMs-of-interest.csv',
                   help="Path to list of PTMs to analyze. (default: ./PTMs-of-interest.csv)")
    p.add_argument('-c', '--contaminants', dest='contaminants', type=existing_file, default=None,
                   help="Path to contaminants FASTA file, used to drop contaminant proteins (default: none)")
    p.add_argument('-s', '--std_filter', dest='std_filter', type=float, default=.05,
                   help="Standard deviation filtering cutoff (default: 0.05)")
    p.add_argument('--thresh', dest='thresh', type=int, default=5,
                   help="Minimum number of observations required in rows and columns (default: 5)")
    return p.parse_args()

def existing_file(path: str) -> str:
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found! --> {path}")
    else:
        return path


def main():
    args = parse_cli()

    # Step 1: map peptidoforms to proteins and count PSMs per PTM.
    PATH = get_peptidoforms_counts_mapped(
        args.peptidoform_ids,
        args.peptidoform_counts,
        args.fasta,
        p=args.prefix,
        t=args.threads
    )

    # Step 2: absolute counts to relative PTM counts.
    PATH = get_relative_ptms(PATH)

    # Step 3: one PTM-by-experiment submatrix per PTM of interest.
    PATH = generate_ptm_submatrices(
        PATH,
        target_fasta_path=args.fasta,
        contam_fasta_path=args.contaminants,
        myptms=args.myptms
    )

    # Step 4: combine the submatrices into the MoDPA matrix.
    combine_ptm_submatrices(
        PATH,
        prefix=args.prefix,
        myptms=args.myptms,
        std_filter=args.std_filter,
        thresh=args.thresh
    )


if __name__ == "__main__":
    main()
