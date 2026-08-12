#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import polars as pl
import numpy as np
from MoDPA import MoDPA
from datetime import date, datetime
import os, argparse, re

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("data_folder", type=existing_folder, help="Path to processed folder with the data.")
    p.add_argument("-p", "--prefix",  type=str, default=date.today().isoformat() , help="Prefix for output files (Default: today's date)")
    p.add_argument('-m', '--myptms', dest='myptms', type=existing_file, default='./PTMs-of-interest.csv', 
                   help="Path to list of PTMs to analyze. (default: ./PTMs-of-interest.csv)")
    p.add_argument('-s', '--std_filter', dest='std_filter', type=float, default=.05,
                   help="Standard deviation filtering cutoff (default: 0.05)")
    p.add_argument('-t', '--thresh', dest='thresh', type=int, default=5, help="Minimum number of nobservations required in rows and columns (default: 5)")
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

def percent_nonzero(x):
    y = x > 0
    nonzero = y.sum().sum() / (y.shape[0]*y.shape[1])
    print(f"% nonzero cells = {nonzero:.2%}")
    # print(y.sum(axis=1).describe())


def combine_ptm_submatrices(
    data_folder,
    prefix,
    myptms="./PTMs-of-interest.csv",
    std_filter=0.05,
    thresh=5
):
    myptms = pl.read_csv(myptms) # csv with aminoacid, unimod id, and name
    myptms = myptms.rows()
    
    # In[3]:
    data = []
    theor_tot_rows, theor_tot_cols = 0, 0
    for n,(i,j,k) in enumerate(myptms):
        print(n+1,'/',len(myptms), (i,j,k))
        filename = os.path.join(data_folder, f'MoDPA_Rel_[{j}]{k}_{i}.pkl.gz')
        try:
            tmp = pd.read_pickle(filename).T
        except Exception as e:
            print(e)
            print("skipping...")
            continue
        # tmp.dropna(thresh=10, inplace=True)
        tmp = tmp[tmp.apply(np.nanstd, axis=1)>=std_filter].copy(deep=True)
        print(tmp.shape)
        if not tmp.empty:
            data.append(tmp)
        theor_tot_rows += tmp.shape[0]
        theor_tot_cols += tmp.shape[1]
        del tmp
    if not data:
        raise ValueError("No PTM submatrices were available for combination.")
    data = pd.concat(data)
    data.index.name, data.columns.name = None, None
    print('sanity check:', data.shape, '|', (theor_tot_rows, theor_tot_cols))
    
    
    # In[4]:
    print("Dataset size =", data.shape)
    print('# Proteins =', len(set([_.split('|')[0] for _ in data.index])))
    print('# PTMs     =', len(set([_ for _ in data.index])))
    print('value range =', [np.nanmin(data.values), np.nanmax(data.values)])
    
    
    # In[5]:
    os.makedirs(data_folder, exist_ok=True)
    
    r, c = thresh, thresh
    analyzed_data = MoDPA(
        f"{prefix}-PTMs-thresh{r}r{c}c-std{std_filter}",
        data,
        thresh=(r,c)
    )
    
    outpath = os.path.join(data_folder, analyzed_data.name+'.pkl.gz')
    
    tmp = pd.DataFrame(analyzed_data.modpa_matrix)
    percent_nonzero(tmp)
    
    MoDPA_matrix = pd.DataFrame(analyzed_data.modpa_matrix, 
                                index=analyzed_data._ptms, 
                                columns=analyzed_data._exps)
    MoDPA_matrix = MoDPA_matrix.apply(lambda col: col/np.nanmax(col), axis=0)
    MoDPA_matrix.to_pickle(outpath, compression='gzip')
    
    analyzed_ptms = MoDPA_matrix.index.tolist()
    analyzed_ptms_df = pd.DataFrame(
        [ [re.sub(r'\[(\d+)\].+', r'\1', _)] + _.split('|') for _ in analyzed_ptms ],
        columns=['PTM_ID','Gene','POS','RES','MOD']
    )
    analyzed_ptms_df.to_csv(
        os.path.join(data_folder, analyzed_data.name+'-analyzed-ptms.csv'), 
        index=False, 
        encoding='utf-8'
    )
    print("\nPTM counts summary:")
    print(pd.DataFrame(analyzed_ptms_df.MOD.value_counts()).reset_index())

    return outpath


def main():
    START = datetime.now()
    args = parse_cli()
    combine_ptm_submatrices(
        data_folder=args.data_folder,
        prefix=args.prefix,
        myptms=args.myptms,
        std_filter=args.std_filter,
        thresh=args.thresh
    )
    END = datetime.now()
    print("Started: ", START.isoformat())
    print("Finished:", END.isoformat(), '\n')

if __name__ == "__main__":
    main()