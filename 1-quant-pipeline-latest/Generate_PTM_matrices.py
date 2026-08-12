#!/usr/bin/env python
# coding: utf-8
import polars as pl
import pandas as pd
import numpy as np
import re, argparse, os
from datetime import datetime
from PTMmap import Fasta

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("input_path", type=existing_file, help="Path to the input CSV file (***_PTMs_counts_relative.csv.gz).")
    p.add_argument('-t', '--targets', dest='target_fasta_path', type=existing_file, default=None,
                   help="Path to FASTA file with target proteins to keep.")
    p.add_argument('-c', '--contaminants', dest='contam_fasta_path', type=existing_file, default=None,
                   help="Path to FASTA file with contaminant proteins to remove.")
    p.add_argument('-m', '--myptms', dest='myptms', type=existing_file, default='./PTMs-of-interest.csv',
                   help="Path to list of PTMs to analyze. (default: ./PTMs-of-interest.csv)")
    return p.parse_args()

def existing_file(path: str) -> str:
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found! --> {path}")
    return path

# def existing_folder(path: str) -> str:
#     if not os.path.isdir(path):
#         raise argparse.ArgumentTypeError(f"Folder not found! --> {path}")
#     return path

def get_protein_set_from_fasta(fasta_path):
    fasta = {}
    Fasta.getFasta(fasta_path, fasta)
    return set(fasta.keys())


def get_MoDPA_matrix_newer(relcounts_, myptms_):
    RES, ID, MOD = myptms_
    mod_of_interest = f"{RES}|[{ID}]{MOD}"

    relcounts_ = relcounts_[relcounts_.PTM_type == mod_of_interest].copy(deep=True)
    relcounts_.sort_values('ptm_name', inplace=True)

    try:
        return relcounts_.pivot(columns='file_name', index='PTM_ID', values='relative_psm_counts')
    except ValueError:
        pass

    relcounts_matrix = pd.DataFrame(relcounts_.groupby(['PTM_ID', 'file_name'])['relative_psm_counts'].max())
    ptms = list(relcounts_matrix.index.get_level_values(0).unique())

    tmp = pd.DataFrame()
    for i, j in enumerate(ptms):
        tmp2 = relcounts_matrix.loc[j]
        if len(tmp2) >= 10:
            tmp2.columns = [j]
            tmp = pd.concat((tmp, tmp2.T))

    tmp.index.names = ['PTM_ID']
    tmp = tmp.reset_index().sort_values('PTM_ID').set_index('PTM_ID').T
    tmp = tmp.reset_index().sort_values('file_name').set_index('file_name').T
    return tmp


def generate_ptm_submatrices(
    input_path,
    target_fasta_path=None,
    contam_fasta_path=None,
    myptms='./PTMs-of-interest.csv'
    ):
    # --- Prefiltering (C) ---
    if target_fasta_path:
        prot_targets = get_protein_set_from_fasta(target_fasta_path)
        print('# Target proteins =', len(prot_targets))
    if contam_fasta_path:
        prot_contaminants = get_protein_set_from_fasta(contam_fasta_path)
        print('# Contaminants =', len(prot_contaminants))

    relcounts = pl.scan_csv(
        input_path,
        schema_overrides={'relative_psm_counts': pl.Float64}
    ).select(
        ['file_name', 'UniAcc', 'ptm_loc', 'ptm_res', 'ptm_name', 'classification', 'relative_psm_counts']
    )

    relcounts = relcounts.filter(
        (pl.col('relative_psm_counts') > 0) &
        (pl.col('relative_psm_counts') < 1)
    )
    if target_fasta_path:
        relcounts = relcounts.filter(
            pl.col('UniAcc').is_in(prot_targets)
        )
    if contam_fasta_path:
        relcounts = relcounts.filter(
            ~pl.col('UniAcc').is_in(prot_contaminants)
        )
    relcounts = relcounts.collect(engine="streaming")
    print('Shape after filtering:', relcounts.shape)  # shape only available after collect
    print("PTMs in dataset:")
    print(relcounts.select('ptm_name').unique())

    # --- MoDPA matrix computation (D) ---
    relcounts = relcounts.with_columns(
        pl.struct(['UniAcc', 'ptm_loc', 'ptm_res', 'ptm_name']).map_elements(
            lambda row: f"{row['UniAcc']}|{row['ptm_loc']}|{row['ptm_res']}|{row['ptm_name']}",
            return_dtype=pl.String
        ).alias('PTM_ID'),
        pl.col('ptm_name').map_elements(
            lambda x: int(re.match(r'\[(\d+)\]', x).groups()[0]),
            return_dtype=pl.Int32
        ).alias('unimod_id'),
        (pl.col('ptm_res') + '|' + pl.col('ptm_name')).alias('PTM_type')
    )
    print('#Raw files =', relcounts.select('file_name').n_unique())
    print('#Modifications =', relcounts.select('PTM_ID').n_unique())
    relcounts = relcounts.to_pandas()
    print(relcounts.ptm_name.value_counts())

    myptms = pl.read_csv(myptms)
    myptms = myptms.rows()
    print(myptms)

    parent_dir, input_filename = os.path.split(input_path)
    prefix = input_filename.replace('_PTMs_counts_relative.csv.gz', '')
    out_dir = os.path.join(
        parent_dir, 
        f'MoDPA_matrices_{prefix}'
    )
    os.makedirs(out_dir, exist_ok=True)
    for n, mod in enumerate(myptms):
        savepath = f'{out_dir}/MoDPA_Rel_[{mod[1]}]{mod[2]}_{mod[0]}.pkl.gz'
        relcounts_matrix = get_MoDPA_matrix_newer(relcounts, mod)
        if len(relcounts_matrix) > 0:
            relcounts_matrix.T.to_pickle(savepath, compression='gzip')
            print(n + 1, '/', len(myptms), '\t',
                  savepath, relcounts_matrix.size, relcounts_matrix.notna().sum().sum())
        else:
            print(n + 1, '/', len(myptms), '\t',
                  mod, relcounts_matrix.size, relcounts_matrix.notna().sum().sum())
    return out_dir


def main():
    START = datetime.now()
    args = parse_cli()
    generate_ptm_submatrices(
        input_path=args.input_path,
        target_fasta_path=args.target_fasta_path,
        contam_fasta_path=args.contam_fasta_path,
        myptms=args.myptms
    )
    END = datetime.now()
    print("Done!!")
    print("Started: ", START.isoformat())
    print("Finished:", END.isoformat(), '\n')

if __name__ == "__main__":
    main()


