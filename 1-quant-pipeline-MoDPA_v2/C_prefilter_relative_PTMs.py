#!/usr/bin/env python
# coding: utf-8
import polars as pl
import pandas as pd
import numpy as np
import re, argparse, os
from PTMmap import Fasta, PTMs_remapping
from string import ascii_uppercase

# In[ ]:
def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    # p.add_argument("data_folder", type=existing_folder, help="Path to processed folder with the data.")
    # p.add_argument("date", type=str, help="Date in YYYY-MM-DD format.")
    p.add_argument("input_path", type=str, help="Path to the input CSV file (***_PTMs_counts_relative.csv.gz).")
    p.add_argument('--targets', dest='target_fasta_path', type=existing_file, default=None, help="Path to FASTA file with target proteins to keep.")
    p.add_argument('--contaminants', dest='contam_fasta_path', type=existing_file, default=None, help="Path to FASTA file with contaminant proteins to remove.")
    return p.parse_args()

def existing_file(path: str) -> str:
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found! --> {path}")
    else:
        return path

def existing_folder(path: str) -> str:
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Folder not found! --> {path}")
    else:
        return path

def get_protein_set_from_fasta(fasta_path):
    fasta = {}
    Fasta.getFasta(fasta_path, fasta)
    return set(fasta.keys())

# In[ ]:
def main():
    args = parse_cli()

    if args.target_fasta_path:
        prot_targets = get_protein_set_from_fasta(args.target_fasta_path)
        print('# Target proteins =',len(prot_targets))
    if args.contam_fasta_path:
        prot_contaminants = get_protein_set_from_fasta(args.contam_fasta_path)
        print('# Contaminants =',len(prot_contaminants))

    in_path  = args.input_path
    out_path = in_path.replace('.csv.gz', '_prefiltered.csv.gz')

    relcounts = pl.read_csv(in_path)
    print(relcounts.shape)

    relcounts = relcounts.lazy()
    relcounts = relcounts.filter(
        (pl.col('relative_psm_counts') > 0)&
        (pl.col('relative_psm_counts') < 1)    
    )
    if args.target_fasta_path:
        relcounts = relcounts.filter(
            pl.col('UniAcc').is_in(prot_targets)
        )
    if args.contam_fasta_path:
        relcounts = relcounts.filter(
            ~pl.col('UniAcc').is_in(prot_contaminants)
        )
    relcounts = relcounts.collect()
    print(relcounts.shape)
    print("PTMs in dataset:")
    print(relcounts.select('ptm_name').unique())

    relcounts.to_pandas().to_csv(out_path, index=False, compression='gzip', encoding='utf-8')
    print("Prefiltered data saved to:", out_path)

# In[ ]:
if __name__ == "__main__":
    main()