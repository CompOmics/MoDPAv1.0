#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import polars as pl
import numpy as np
from MoDPA import MoDPA
from datetime import date, datetime
import os, argparse

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("data_folder", type=str, help="Path to processed folder with the data.")
    p.add_argument("date", type=str, help="Date in YYYY-MM-DD format.")
    p.add_argument('--myptms', dest='myptms', type=str, default='./PTMs-of-interest.csv', 
                   help="Path to list of PTMs to analyze. (default: ./PTMs-of-interest.csv)")
    p.add_argument('--std_filter', dest='std_filter', type=float, default=.05,
                   help="Standard deviation filtering cutoff (default: 0.05)")
    return p.parse_args()

def percent_nonzero(x):
    y = x > 0
    nonzero = y.sum().sum() / (y.shape[0]*y.shape[1])
    print(f"% nonzero cells = {nonzero:.2%}")
    print(y.sum(axis=1).describe())


# In[2]:
args = parse_cli()

myptms = pl.read_csv(args.myptms) # csv with aminoacid, unimod id, and name
myptms = myptms.rows()

# In[3]:
data = []
theor_tot_rows, theor_tot_cols = 0, 0
for n,(i,j,k) in enumerate(myptms):
    print(n+1,'/',len(myptms), (i,j,k))
    filename = f"{args.data_folder}/MoDPA_Rel_[{j}]{k}_{i}.pkl.gz"
    tmp = pd.read_pickle(filename).T
    # tmp.dropna(thresh=10, inplace=True)
    tmp = tmp[tmp.apply(np.nanstd, axis=1)>=args.std_filter].copy(deep=True)
    print(tmp.shape)
    data.append(tmp)
    theor_tot_rows += tmp.shape[0]
    theor_tot_cols += tmp.shape[1]
    del tmp
data = pd.concat(data)
data.index.name, data.columns.name = None, None
print('sanity check:', data.shape, '|', (theor_tot_rows, theor_tot_cols))


# In[4]:
print("Dataset size =", data.shape)
print('# Proteins =', len(set([_.split('|')[0] for _ in data.index])))
print('# PTMs     =', len(set([_ for _ in data.index])))
print('value range =', [np.nanmin(data.values), np.nanmax(data.values)])

PTMs = pd.DataFrame([ [_] + _.split("|") for _ in data.index], columns=['PTM_ID','Gene','POS','RES','MOD'])
PTMs['res_mod'] = PTMs.RES + PTMs.MOD
pd.DataFrame(PTMs.MOD.value_counts()).reset_index()


# In[5]:
out_subfld = 'modpa-matrices'
os.makedirs(os.path.join(args.data_folder,out_subfld), exist_ok=True)

r,c = (5,5)
analyzed_data = MoDPA(f"{args.date}-PTMs-thresh{r}r{c}c-std{args.std_filter}", data.fillna(0), thresh=(r,c))

outpath = os.path.join(args.data_folder,out_subfld, analyzed_data.name+'.pkl.gz')
print(outpath)

tmp = pd.DataFrame(analyzed_data.modpa_matrix)
percent_nonzero(tmp)

MoDPA_matrix = pd.DataFrame(analyzed_data.modpa_matrix, 
                            index=analyzed_data._ptms, 
                            columns=analyzed_data._exps)
print(MoDPA_matrix.shape)
MoDPA_matrix = MoDPA_matrix.apply(lambda col: col/np.nanmax(col), axis=0)
print(np.nanmin(MoDPA_matrix),'-',np.nanmax(MoDPA_matrix))
MoDPA_matrix.to_pickle(outpath, compression='gzip')
outpath

PTMs = pd.DataFrame([ [_] + _.split("|") for _ in analyzed_data._ptms], columns=['PTM_ID','Gene','POS','RES','MOD'])
PTMs.to_csv(os.path.join(args.data_folder,out_subfld,f'{args.date}-analyzed-ptms.csv'), index=False)
PTMs['res_mod'] = PTMs.RES + PTMs.MOD
print(pd.DataFrame(PTMs.MOD.value_counts()).reset_index())
